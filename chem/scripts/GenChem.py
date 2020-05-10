#!/usr/bin/env python3
"""
GenChem.py - a Python script to read in equations and output production and
loss terms for Fortran programs.

@author Alan Briolat, SEIY  2013,
          converted from Dave Simpson perl script GenChem.pl
        Hannah Imhof, Chalmers,  2016
        John Johansson, Chalmers, 2017-2018
        Dave Simpson, 2013-2019

"""

# READING classes:
# - ShorthandMap
# - SpeciesReader
# - ReactionsReader (types of species, e.g. explicit catalyst in 
#   ._read_reaction_term
# WRITING classes:
# - DimensionsWriter (new: CM_ChemDims_mod.f90)
# - SpeciesWriter (CM_ChemSpecs_mod.f90, CM_DryDep.inc, CM_WetDep.inc)
# - GroupsWriter (CM_ChemGroups_mod.f90)
# - ReactionsWriter (CM_Reaction1/2.inc, CM_ChemRates_mod.f90, CM_EmisFile.inc,
#    CM_EmisSpecs.inc, femis.defaults, CM_emislist.csv)
# OTHER classes:
# - ChemicalScheme (-> species, groups, emissions, reactions)
# - Species (-> name, formula, dry/wet deposition)
# - Reaction (-> rate, left/right hand side)
# - CodeGenerator (writing header to output files, (un)indenting)
# - PrettyStreamHandler (colours in log output)


import logging
import csv
import re
import collections
import itertools
from textwrap import dedent
import sys

# how many continuation lines for production and loss terms
n_continuation_lines = 20

schemename = ""

def indent(data, amount=1, string='  '):
    return ''.join(string * amount + line for line in data.splitlines(True))


def split(s, sep=None, maxsplit=-1):
    """Call ``s.split(...)`` and strip whitespace.

    Returns an empty list if *s* is empty.  This differs from ``str.split(...)``
    behaviour which would return ``['']``."""
    if s:
        return [x.strip() for x in s.split(sep, maxsplit)]
    else:
        return []


def ichunk(iterable, size):
    i = iter(iterable)
    while True:
        chunk = list(itertools.islice(i, size))
        if len(chunk) == 0:
            break
        else:
            yield chunk


def element_remainder_pairs(elements):
    """'ABC' -> [('A', 'BC'), ('B', 'AC'), ('C', 'AB')]"""
    for i, x in enumerate(elements):
        yield (x, elements[:i] + elements[i+1:])


def expression_wrap(expr, maxlen, split):
    """Wrap an expression to a maximum length.

    After splitting *expr* after any run of the characters in *split*, return a
    list of expression fragments that are ideally no longer than *maxlen*. Some
    may be longer if the expression couldn't be split any smaller. If *maxlen*
    is None, then there will be a new fragment for every split.

    >>> expression_wrap('abc+def+ghi', 8, '+')
    ['abc+def+', 'ghi+jkl']
    >>> expression_wrap('abc+def+ghi', None, '+')
    ['abc+', 'def+', 'ghi+', 'jkl']
    """
    regex = re.compile(r'([{}]+)'.format(re.escape(split)))
    parts = regex.split(expr)

    # Join separators with the string they follow
    parts = [a + b for a, b in itertools.zip_longest(parts[0::2], parts[1::2],fillvalue='')]

    # If we're not recombining up to maxlen, return the elements now
    if maxlen is None:
        return parts

    # Combine parts into chunks that are still no larger than maxlen
    chunks = [parts[0]]
    for p in parts[1:]:
        if len(chunks[-1]) + len(p) > maxlen:
            chunks.append(p)
        else:
            chunks[-1] += p
    # Strip extraneous whitespace
    return [c.strip() for c in chunks]


class SimpleOrderedSet(collections.OrderedDict):
    """A simple replacement class for OrderedSet from external package ordered-set
    
    Thin wrapper around OrderedDict that only implements the functionality needed here.
    This was included to eliminate external dependencies.
    """
    def add(self, elem):
        self[elem] = None
        return list(self.keys()).index(elem)
    
    def __ior__(self, other):
        for item in other:
            self[item] = None
        return self

class DefaultListOrderedDict(collections.OrderedDict):
    """An ordered dictionary that also acts like ``defaultdict(list)``.

    We have a lot of key -> list mappings in this script, it's nice to keep
    insertion order, so this is a special case that does what we want.

    (this is just like a collections.OrderedDict, but if a key is given that
    is not in the dict yet, now an empty list will be added.)
    """
    def __getitem__(self, key):
        try:
            return collections.OrderedDict.__getitem__(self, key)
        except KeyError:
            self[key] = new_list = list()
            return new_list


class IndentingLogger(logging.LoggerAdapter):
    """Stateful logger adapter that indents messages.

    Provides :meth:`indent` and :meth:`outdent` to increase and decrease the
    indent level.  All messages have the indentation prepended to them when
    they pass through the adapter.

    >>> log = IndentingLogger(logging.getLogger('foo'))
    >>> log.debug('hello world')
    >>> log.indent()
    >>> log.debug('I am indented')
    >>> log.outdent()
    >>> log.debug('and now I am not')
    """
    def __init__(self, log):
        super(IndentingLogger, self).__init__(log, None)
        self._indent_level = 0

    def indent(self):
        self._indent_level += 1

    def outdent(self):
        self._indent_level = max(0, self._indent_level - 1)

    def process(self, msg, kwargs):
        return (indent(msg, self._indent_level), kwargs)


LOG = IndentingLogger(logging.getLogger(''))


class IndentingStreamWriter(object):
    """Stateful stream wrapper that indents written strings.

    Provides :meth:`indent` and :meth:`outdent` to increase and decrease the
    indent level.  All messages have the indentation prepended to them when
    they pass through :meth:`write` unless ``indent=False`` is present.
    """
    def __init__(self, stream):
        self._stream = stream
        self._indent_level = 0

    def indent(self):
        self._indent_level += 1

    def outdent(self):
        self._indent_level = max(0, self._indent_level - 1)

    def write(self, data, skip_indent=False):
        if not skip_indent:
            data = indent(data, self._indent_level)
        self._stream.write(data)


def is_numeric(n):
    """Check if *n* is numeric.

    This replaces ``is_integer`` and ``is_float`` from ``GenChem.pl``.  Valid
    floats are a superset of valid integers, so just uses ``float()`` to check.
    """
    try:
        x = float(n)
        return True
    except (TypeError, ValueError):
        return False


class ShorthandMap(object):
    """Read shorthands from *stream*.

    Given a file with a format like::

        * Ignore this line
        XT           temp(iq)
        H2O          H2O(iq)        ! A comment
        FH2O         (1.+1.4e-21*h2o*exp(2200./XT))

    creates a mapping from shorthand to expansions, which can be applied to
    other equations with the :meth:`expand` method.  Shorthand that appears in
    other expansions is expanded during reading, so in the above example::

        self.mapping['FH2O'] = '(1.+1.4e-21*H2O(IQ)*EXP(2200./TEMP(IQ)))'
    """
    def __init__(self, stream):
        self.mapping = collections.OrderedDict()

        dtxt='SHMap:'
        LOG.info(dtxt+'Processing shorthands...')
        LOG.indent()
        for line in stream:
            line = line.strip().upper()
            # Skip empty lines and comments
            if line and not line.startswith('*'):
                # Split into pattern, expansion, comments
                parts = line.split(None, 3)

                pattern = parts[0]
                # Expand any shorthand that appears in the expansion
                #  add condition, so that not too many unnecessary
                #    brackets are being added
                OPS = list('+-*/')
                if any(s in parts[1] for s in OPS):
                    expanded = '(' + self.expand(parts[1]) + ')'
                else:
                    expanded = self.expand(parts[1])
                LOG.debug(dtxt+'DSPARTS %s %s : %s'% ( parts[0], parts[1], expanded ))

                self.mapping[pattern] = expanded
                LOG.debug(dtxt+'%-12s  =>  %s', pattern, expanded)

        LOG.outdent()
        LOG.info(dtxt+'END %s shorthands processed.', len(self.mapping))

    def expand(self, eqn):
        """Expand shorthand in *eqn*."""
        """ This method dailed when the key had parentheses, e.g. cbphot(1). Haven't
            figured out whyv yet. Used very crude ds stuff below """
        dtxt='SHMap:expand:'
        dbg= len(self.mapping) > 138
        if dbg: LOG.debug('\n'+dtxt+'DSEXPANDING n=%d ' % len(self.mapping) + eqn )
        if len(self.mapping) == 0:
            LOG.debug(dtxt+'DSRETURN QUERY?')
            return eqn

        dsret=eqn
        dsfix=False
        for key in self.mapping:
         #LOG.debug(dtxt+'DStry %s map %s'% ( key, self.mapping[key] )) # => DSKEY RCPHOT(1)
         if 'CBPHOT' in key: # just for CB6 now. Expand/fix later
           if key in dsret:
              new   = self.mapping[key]
              dsret = dsret.replace(key,new)
              # code won't work with factors yet if () kept
              if dsret.startswith('('): dsret = dsret[1:-2] # Argh, 
              LOG.debug('DSKEY %s map %s %s'% ( 
                  key, self.mapping[key], dsret )) # => DSKEY RCPHOT(1)
              dsfix=True
        pattern = r'\b(' + '|'.join(self.mapping) + r')\b' # gives e.g |KRO2NO|KNO3|rcphot(1)...

        def replace(matchobj):
            return self.mapping[matchobj.group(0)]

        LOG.debug(dtxt+'DSRETURNS %s '% re.sub(pattern, replace, eqn))
        if dsfix:
          return dsret
        else:
          return re.sub(pattern, replace, eqn)


#: Atomic weights of atoms we care about
ATOMS = {
    'C': 12,
    'H': 1,
    'N': 14,
    'O': 16,
    'S': 32,
}


def count_atoms(formula):
    """formula -> (counts, molwt)

    Process a formula to find out how many of each atom in *ATOMS* it contains
    and the molecular weight from only those atoms.
    """
    counts = collections.Counter()
    # Count atoms in the formula
    for atom, n in re.findall('([A-Z][a-z]?)([+-]?(?:[0-9]*[.])?[0-9]*)', formula):
        n = float(n) if n else 1
        if atom in ATOMS:
            counts.update({atom: n})
    # Sum the molecular weight
    molwt = float(sum(ATOMS[atom] * n for atom, n in counts.items()))
    return (counts, molwt)


class ChemicalScheme(object):
    def __init__(self):
        self.species = collections.OrderedDict()
        self._species_nocase = {}
        self.groups = DefaultListOrderedDict()
        self.groups_factors = collections.OrderedDict()
        self.groups_maps = collections.OrderedDict()
        self.emis_files = SimpleOrderedSet()
        self.reactions = []

    def add_species(self, spec):
        """Add *spec* to species, updating groups as necessary."""
        if spec.name in self.species:
            raise ValueError("SPEC %s already defined!" % spec.name)

        self.species[spec.name] = spec
        self._species_nocase[spec.name.upper()] = spec

        for group in spec.groups:
            self.groups[group].append(spec.name)
        # HI also for new groups. But here I now have a OrderedDict in a 
        #  OrderedDict to keep the name of the new group, the species and the
        #  value
        for group_factor in list(spec.groups_factors.keys()):
            if group_factor in list(self.groups_factors.keys()):
                self.groups_factors[group_factor].update({spec.name: \
                    spec.groups_factors[group_factor]})
            else:
                self.groups_factors[group_factor] = collections.OrderedDict()
                self.groups_factors[group_factor].update({spec.name: \
                    spec.groups_factors[group_factor]})
        for group_map in list(spec.groups_maps.keys()):
            if group_map in list(self.groups_maps.keys()):
                self.groups_maps[group_map].update({spec.name: \
                    spec.groups_maps[group_map]})
            else:
                self.groups_maps[group_map] = collections.OrderedDict()
                self.groups_maps[group_map].update({spec.name: \
                    spec.groups_maps[group_map]})

    def find_species(self, name):
        """Find a species by case-insensitive lookup."""
        return self._species_nocase[name.upper()]

    def add_emis_files(self, files):
        """Add *files* to emis_files."""
        self.emis_files |= files
        LOG.info('EMISFILES added %s, now: %s', files, self.emis_files)

    def add_reaction(self, reaction):
        """Add *reaction* to reactions."""
        self.reactions.append(reaction)

    def get_species_list(self):
        """Get a list of species sorted by type."""
        return sorted(list(self.species.values()), key=lambda spec: spec.type)

    def get_group_list(self):
        """Get a list of groups."""
        return list(self.groups.items())

    # HI same for new groups, but not items because I need to preserve
    #  the dict in dict kind
    def get_group_factor_list(self):
        """Get a list of factor groups."""
        return self.groups_factors

    def get_group_map_list(self):
        """Get a list of map groups."""
        return self.groups_maps

    def get_dry_deposition_map(self):
        """Get mapping of species to dry deposition surrogates."""
        mapping = collections.OrderedDict()
        for spec in self.species.values():
            if spec.dry is not None:
                mapping[spec.name] = spec.dry
        return mapping

    def get_wet_deposition_map(self):
        """Get mapping of species to wet deposition surrogates."""
        mapping = collections.OrderedDict()
        for spec in self.species.values():
            if spec.wet is not None:
                mapping[spec.name] = spec.wet
        return mapping


class Species(object):
    # Species type constants
    SHORT_LIVED = 0
    ADVECTED = 1
    SEMIVOL = 2
    SLOW = 3

    def __init__(self, name, type, formula, \
                 wet=None, dry=None, molwt=None, groups=None, \
                 groups_factors=None, groups_maps=None):
        # Arguments simply copied to attributes
        self.name = name
        self.type = type
        self.formula = formula
        self.wet = wet              # Wet deposition surrogate species
        self.dry = dry              # Dry deposition surrogate species

        # Calculate atom counts and molecular weight from formula
        self.counts, self.molwt = count_atoms(self.formula or '')
        # Is non-methane hydrocarbon?
        self.NMHC = set(self.counts) == {'C', 'H'} and self.counts['C'] >= 2

        LOG.debug('processed formula: %s  =>  MOLWT=%s, NMHC=%s, COUNTS=%r',
                  self.formula, self.molwt, self.NMHC, dict(self.counts))

        # Override molecular weight?
        if molwt is not None:
            self.molwt = float(molwt)
            LOG.debug('override MOLWT: %s', self.molwt)

        # Groups, default to empty list
        self.groups = groups or []

        #  also define this for new groups. Default values are also set to
        #  the same in SpeciesReader() ~ l. 420
        # groups_factors, default to empty ordered dict
        self.groups_factors = groups_factors or collections.OrderedDict()

        # groups_maps, default to empty ordered dict
        self.groups_maps = groups_maps or collections.OrderedDict()

    def __getitem__(self, key):
        """Redirect dict-like access to attribute access.

        Particularly useful for the old-style string formatting still used by
        the logging module, e.g. ``"%(name)s" % spec``.
        """
        return getattr(self, key)

    def is_advected(self):
        """Is this an advected species?"""
        return self.type > 0


class SpeciesReader(object):
    FIELDS = ('Spec', 'type', 'formula', 'in_rmm', 'dry', 'wet', 'groups', 
                'comment')
    def __init__(self, scheme):
        self.scheme = scheme

    def read(self, stream):
        dtxt='SpeciesReader:'
        """Read species from *stream*."""
        LOG.info(dtxt+'Processing species...')
        LOG.indent()

        reader = csv.DictReader(stream, self.FIELDS)
        for row in reader:
            # The #SLOW line is the old way of doing adv=3, and is inflexible.
            # We don't allow it here.
            if row['Spec'].startswith('#SLOW'):
                raise ValueError('#SLOW not allowed, use adv=3')

            # Skip empty/comment rows
            if not row['Spec'] or row['Spec'] == 'Spec' or row['Spec'][0] in {'*', '#'}:
                continue

            # Replace "blank" fields with None
            for k, v in row.items():
                if v == 'xx' or v == '':
                    row[k] = None

            LOG.debug('SPEC %s', row['Spec'])
            LOG.indent()

            # Process groups
            groups = []
            groups_factors = collections.OrderedDict()
            groups_maps = collections.OrderedDict()
            if row['groups'] is not None:

                for g in row['groups'].split(';'):
                    #  check for ':' in groups, then recognize type of 
                    #  parameter and store it
                    if ':' in g:
                        (special_group, value) = g.split(':')
                        if is_numeric(value): # convert to upper case
                            groups_factors[special_group.upper()] = value
                        else:
                            groups_maps[special_group.upper()] = value
                    # that was what happened in the older version anytime in
                    #  the loop. But now I don't want the new groups in the
                    #  normal list of groups
                    else:
                        #  convert to upper case here as well, though it's
                        #  not really necessary (so far)
                        groups.append(g.upper())

                    if row['wet'] is not None:
                        #  we don't want to add the new groups here
                        if ':' not in g:
                            groups.append('WDEP_' + g.upper())
                    if row['dry'] is not None:
                        # same here
                        if ':' not in g:
                            groups.append('DDEP_' + g.upper() )
                LOG.debug('In groups: ' + ', '.join(groups))
                LOG.debug('In factor groups: ' + ', '.join(groups_factors))
                LOG.debug('In map groups: ' + ', '.join(groups_maps))

            spec = Species(name=row['Spec'],
                           type=int(row['type']),
                           formula=row['formula'],
                           wet=row['wet'],
                           dry=row['dry'],
                           molwt=row['in_rmm'],
                           groups=groups,
                           groups_factors=groups_factors,
                           groups_maps=groups_maps)

            self.scheme.add_species(spec)

            LOG.outdent()

        LOG.outdent()


Term = collections.namedtuple('Term', ['species', 'factor', 'type'])


class Reaction(object):
    def __init__(self, rate, LHS, RHS):
        self.rate = rate
        # Make sure we don't have factors on the LHS
        assert all(_.factor == '1.0' for _ in LHS), 'factors not allowed in LHS'
        self.LHS = LHS
        # Make sure we don't have [] terms on the RHS
        assert all(_.type != 'passive' for _ in RHS), '[] terms not allowed in RHS'
        self.RHS = RHS

    def __repr__(self):
        return 'Reaction(rate=%r, LHS=%r, RHS=%r)' % (self.rate, self.LHS, self.RHS)

    def __str__(self):
        return self.__repr__()

    @property
    def reactants(self):
        """Reaction reactants (from LHS)."""
        return [_ for _ in self.LHS if _.type is None]

    @property
    def passives(self):
        """Reaction passives (from LHS)."""
        return [_ for _ in self.LHS if _.type == 'passive']

    @property
    def explicit_catalysts(self):
        """Reaction explicit_catalysts (from LHS)."""
        return [_ for _ in self.LHS if _.type == 'explicit_catalyst']

    @property
    def products(self):
        """Reaction products (from RHS)."""
        return [_ for _ in self.RHS if _.type is None]

    @property
    def dummy_catalysts(self):
        """Reaction dummy_catalysts (from RHS)."""
        return [_ for _ in self.RHS if _.type == 'dummy_catalyst']

    def get_full_rate(self):
        """Full rate, including rates contributed by LHS passives and \
        explicit_catalysts."""
        rate = self.rate[:]
        for term in self.passives:
            rate.extend(['*', ('amount', term.species)])
        return rate

    def get_loss_rates(self):
        """Iterate over ``(species, rate)`` pairs for this reaction's losses."""
        base_rate = self.get_full_rate()
        for term, others in element_remainder_pairs(self.reactants):
            # Multiply rate by amounts of all other reactants
            term_rate = []
            for t in others:
                term_rate.extend(['*', ('amount', t.species)])
            yield (term.species, base_rate + term_rate)

    def get_prod_rates(self):
        """Iterate over ``(species, rate)`` pairs for this reaction's production."""
        base_rate = self.get_full_rate()
        # Multiply rate by amounts of all reactants
        for term in self.reactants:
            base_rate.extend(['*', ('amount', term.species)])
        for term in self.products:
            # Include product factor if present
            if term.factor == '1.0':
                term_rate = []
            else:
                term_rate = [term.factor + '*']
            yield (term.species, term_rate + base_rate)


class ReactionsReader(object):
    def __init__(self, scheme, shorthand):
        self.scheme = scheme
        self.shorthand = shorthand
        # ?? add ct_coeff to keep track of number of coefficients for the
        #  reactions log file. Better idea?
        #  and also rcts to keep track of which rates have already been used
        self.ct_coeff = 0
        self.rcts = []

    def read(self, stream, log_reac):
        """Read reactions from *stream*."""
        LOG.info('Processing reactions...')
        LOG.indent()
        from_previous_lines = ''
        for linenum, line in enumerate(stream, 1):
            # Strip whitespace, convert to uppercase
            line = line.strip().upper()
            # Skip empty lines and comments
            if not line or line.startswith('*'):
                continue
            
            # emisfiles line now treated separately
            if line.startswith('EMISFILES:'):
                LOG.info('Line %5d: %s', linenum, line)
                LOG.indent()
                files = split(line, ':', 1)[1].lower()
                self.scheme.add_emis_files(split(files, ','))
                LOG.outdent()
                from_previous_lines = ''
                continue
            
            # Code added to append consecutive lines until ; is found
            if ';' in line:
                line = from_previous_lines + line
                from_previous_lines = ''
            else:
                from_previous_lines = from_previous_lines +  line
                continue
            # Strip comment from end of line
            line = line.partition(';')[0].strip()

            LOG.info('Line %5d: %s', linenum, line)
            LOG.indent()

            self.scheme.add_reaction(self._read_reaction(line, log_reac))

            LOG.outdent()

        LOG.outdent()

    #  changed format of reactions, added ':' between rate and reaction
    #  and added reaction log file out_log
    def _read_reaction(self, reaction, out_log):
        """Turn a "<rate> : <LHS> = <RHS>" definition into a Reaction object."""
        assert '=' in reaction
        # Split up "<rate> : <reactants> = <products>"
        # now split on ':', also check for amount of ':'
        if reaction.count(':') != 1:
            raise EnvironmentError('more than one or no colon at all in this reaction!')
        rate, terms = reaction.split(':', 1)
        rate = rate.strip().replace(' ', '')
        terms = terms.strip()
        lhs, rhs = split(terms, '=', 1)
        # Parse LHS and RHS terms
        lhs = [self._read_reaction_term(_) for _ in split(lhs, '+')]
        rhs = [self._read_reaction_term(_) for _ in split(rhs, '+')]

        #  add explicit_catalysts here so these reaction rates will be
        #  written into CM_ChemRates_mod
        #  !! assuming that concentration can simply be multiplied to rate
        #     coefficient!!
        for lt in lhs:
            if lt[2] == 'explicit_catalyst':
                rate = '(' + rate + ')*' + lt[0]

        # Process rate and create reaction
        rate = self._read_reaction_rate(rate)
        self._log_reaction(out_log, rate, reaction)
        reaction = Reaction(rate, lhs, rhs)
        LOG.debug('%s', reaction)

        # Calculate atom sums for both sides
        lhs_sum = self._sum_atoms(reaction.LHS)
        rhs_sum = self._sum_atoms(reaction.RHS)
        LOG.debug('LHS atoms: %r', lhs_sum)
        LOG.debug('RHS atoms: %r', rhs_sum)
        # Warn if the two sides don't match
        if any(abs(lhs_sum[a] - rhs_sum[a]) > 0.001 for a in ATOMS):
            LOG.warning('Reaction unbalanced! LHS=%r, RHS=%r', lhs_sum, rhs_sum)

        LOG.debug('full rate including passives: %r', reaction.get_full_rate())

        return reaction

    def _sum_atoms(self, terms):
        """Create a sum of the atoms in all of *terms*."""
        # TODO: check atom counts work properly, match GenChem.pl
        atoms = collections.Counter()
        for term in terms:
            if term.species in self.scheme.species:
                spec = self.scheme.species[term.species]
                spec_atoms = spec.counts
            else:
                spec_atoms, _ = count_atoms(term.species)
            # If factor is not a constant (i.e. it's an expression),
            # it cannot be evaluated here, so instead a factor of 1.0 is used.
            try:
                factor = float(term.factor)
            except ValueError:
                factor = 1.0
            # Accumulate atom counts, scaled by the term's factor
            atoms.update({a: n * factor for a, n in spec_atoms.items()})
        return atoms

    def _read_reaction_term(self, term):
        """reactant or product -> (type, factor, species)"""
        # Split "<factor> <species>" pair
        factor, _, species = term.rpartition(' ')
        # Default factor if not supplied
        factor = factor or '1.0'
        
        # Remove space and pipes (pipes are intended for expressions)
        factor = factor.strip(" |")

        # Figure out type of term
        if species.startswith('[') and species.endswith(']'):
            type = 'passive'
            species = species[1:-1]
            # Check species exists and use its canonical name
            species = self.scheme.find_species(species).name
        elif species.startswith('{') and species.endswith('}'):
            type = 'dummy_catalyst'
            species = species[1:-1]
        elif species.startswith('<') and species.endswith('>'):
            type = 'explicit_catalyst'
            species = species[1:-1]
        else:
            type = None
            # Check species exists and use its canonical name
            species = self.scheme.find_species(species).name

        return Term(species, factor, type)

    def _read_reaction_rate(self, rate):
        """Process *rate* into a list of rate parts.

        The rate first has shorthand expansions and some basic cleanup applied
        to it.  It's then split by "top-level" operators (operators that don't
        appear inside parentheses) and each "part" is checked for special forms
        that need special handling instead of being included as-is in the rate.
        The special forms are replaced with ``(kind, arg)`` tuples, and
        everything else is left as literal strings.

        For example, ``rcphot(IDNO2)*0.222`` becomes ``[('photol', 'IDNO2'), '*0.222']``.
        The responsibility for turning this into a Fortran expression describing
        the rate (e.g. ``rcphot(IDNO2,k)*0.222``) belongs to the code generator.
        (special forms end up in CM_Reactionsx.inc)

        Non-trivial rates that do not cause any special handling become
        ``[('coeff', rate)]``, which gets used in the code generator to re-use
        rate calculations, and end up in CM_ChemRates_mod

        ``KDIM`` in rates that will later be written to CM_Reactionsx
        (special forms; numerics, aqrck or _func_ ones) will be converted to ``k``,
        e.g. ``rcphot(IDNO2(KDIM))`` becomes ``[('photol', 'IDNO2(k)')]``;
        KDIM in coeff-rates (will be written to CM_ChemRates_mod) will be
        converted to ``:``, e.g. ``4.0e-12*Fgas(ASOC_ug1,KDIM)`` becomes 
        ``[('coeff', '4.0e-12*Fgas(ASOC_ug1,:)')]``
        """
        dtxt='ReadR:'
        # Expand shorthands in the rate
        LOG.debug(dtxt+'DSTMP0 %s' % rate)
        expanded_rate = self.shorthand.expand(rate)
        if expanded_rate != rate:
            LOG.debug(dtxt+'expanded rate: %s', expanded_rate)
        rate = expanded_rate
        LOG.debug(dtxt+'DSTMPR %s' % rate)

        # Clean up the rate a bit
        # Lowercase EXP() (TODO: do we need to?)
        #rate = re.sub(r'\bEXP\b', 'exp', rate)
        # Normalise e-notation
        rate = re.sub(r'([\d\.])[EdD]([\+\-]?\d)', r'\1e\2', rate)
        # Replace 1. -> 1.0
        rate = re.sub(r'\.(?=\D)', r'.0', rate)
        LOG.debug(dtxt+'cleaned rate: %s', rate)

        # Short-circuit some cases where we just want the literal rate
        # will end up in CM_Reactions
        # replacement of KDIM, put rcemis here and don't treat it
        #         specially any more
        if is_numeric(rate) or 'AQRCK' in rate or rate.startswith('RCEMIS'):
            return [rate.replace('KDIM', 'k')]
        elif rate.startswith('DJPHOT'):
            return [rate.replace('KDIM', 'k')]
        elif rate.startswith('_FUNC_'):
            return [rate[6:].replace('KDIM', 'k')]

        LOG.debug(dtxt+'handling parts: %s', rate)
        parts = []
        def handle_part(part, op=None):
            # Recognise and replace special forms with (kind, arg) tuples
            # HI took out emissions (rcemis) here!
            # DS if part.startswith('RCPHOT('):
            if 'RCPHOT(' in part.upper():
                # RCPHOT(foo) -> ('photol', 'foo'), HI new: replacement of KDIM
                part = ('photol', part[7:-1].replace('KDIM', 'k'))
                LOG.debug(dtxt+'handled RCPHOT part')

            # If this is the first part, a special case, or follows a special
            # case, save as a new part, otherwise append it to the previous one.
            # (Keeps it down to the most concise form that still represents the
            # rate properly.)
            if isinstance(part, tuple) or len(parts) == 0 or isinstance(parts[-1], tuple):
                parts.append(part)
            else:
                parts[-1] += part

            # If the part has an associated operator (rather than terminated by
            # the end of the string), add that to the most recent part if it
            # was a string, or start a new part if it was a special case.
            if op is not None:
                if isinstance(parts[-1], tuple):
                    parts.append(op)
                else:
                    parts[-1] += op

        # Process the rate, splitting on operators, but never breaking inside
        # parentheses.  Use the handle_part function above to handle each part.
        OPS = set('+-/*')
        part_start = 0
        paren_depth = 0
        for i, c in enumerate(rate):
            if c == '(':
                paren_depth += 1
            elif c == ')':
                paren_depth -= 1
            elif paren_depth == 0 and c in OPS:
                # Found an operator outside of paretheses, end the current part
                # and process it
                handle_part(rate[part_start:i], rate[i])
                # Start a new part
                part_start = i + 1
            else:
                # Proceed to next character, accumulating the current part
                pass
        assert paren_depth == 0, 'mismatched parentheses'
        # Handle final part
        print(dtxt+'doing parts:',rate, part_start, rate[part_start:], type(rate) )
        handle_part(rate[part_start:])

        # If we just have just a single literal expression by this point, it
        # should probably be a rate coefficient (see the short-circuit section
        # at the top of this function for cases where this won't happen).  The
        # code generator can build a database of rate coefficients and reuse
        # calculations.
        # replacement of KDIM
        if parts == [rate]:
            parts = [('coeff', rate.replace('KDIM', ':'))]

        print(dtxt+'done parts:',rate, parts)
#        LOG.debug(dtxt+'done parts: %s', ';'.join(parts))
        return parts

    def _log_reaction(self, out_file, rate, reaction):
        """Logfile containing the reaction rate as noted in CM-files and
        the reaction itself.
        """
        reaction_print = reaction.split(':')[1]
        if isinstance(rate[0], tuple):
            if rate[0][0] == 'coeff':
                if rate[0][1] in self.rcts:
                    i_out = next(j+1 for j, ra in enumerate(self.rcts) 
                        if ra == rate[0][1])
                else:
                    self.ct_coeff += 1
                    i_out = self.ct_coeff
                    self.rcts.append(rate[0][1])
                outlog = 'rct, {:>13d}: {:}\n'.format(i_out, 
                    reaction_print)
            else:
                outlog = 'rc{:4}, {:>10}: {:}\n'.format(rate[0][0][:4],
                    rate[0][1], reaction_print)
        else:
            outlog = '{:18}: {}\n'.format('k, {:>15}'.format(rate[0]), 
                reaction_print)
        out_file.write(outlog)


class CodeGenerator(object):
    def write_module_header(self, stream, module, use=None):
        """Write module header to *stream*, leaving *stream* indented."""
        stream.write(dedent("""\
        ! Generated by GenChem.py - DO NOT EDIT
        ! scheme(s) %s
        """ % schemename))
        stream.write('module {}\n\n'.format(module))
        stream.indent()
        if use is not None:
            for line in use:
                stream.write(line + '\n')
            stream.write('\n')
        stream.write('implicit none\n')
        # Duplicate behaviour of GenChem.pl (TODO: can we *alway* add this?)
        if use is not None:
            stream.write('private\n')
        # DS 2017-03-21 add scheme names as public variables
        modroot=module.replace('_mod','') #DS
        stream.write('character(len=*),parameter, public :: CM_schemes_{} = "{}"\n\n'.format(modroot,schemename))

    def write_contains(self, stream):
        stream.outdent()
        stream.write('\ncontains\n')
        stream.indent()

    def write_module_footer(self, stream, module):
        stream.outdent()
        stream.write('\nend module {}\n'.format(module))


class DimensionsWriter(CodeGenerator):
    def write_beginning(self, stream):
        LOG.info('Writing dimensions')

        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)

        self.write_module_header(stream, 'ChemDims_mod')

    def write_end(self, stream):
        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)
        self.write_module_footer(stream, 'ChemDims_mod')


class SpeciesWriter(CodeGenerator):
    # Groups of species indices to write out
    INDEX_GROUPS = [
        {
            'filter': None,
            'tag': 'TOT',
            'desc': 'All reacting species',
            'ixlab': '',
        },
        {
            'filter': lambda s: s.is_advected(),
            'tag': 'ADV',
            'desc': 'Advected species',
            'ixlab': 'IXADV_',
        },
        {
            'filter': lambda s: s.type == Species.SHORT_LIVED,
            'tag': 'SHL',
            'desc': 'Short-lived (non-advected) species',
            'ixlab': 'IXSHL_',
        },
        {
            'filter': lambda s: s.type == Species.SEMIVOL,
            'tag': 'SEMIVOL',
            'desc': 'Semi-volatile organic aerosols',
            'ixlab': 'IXSOA_',
        },
    ]

    INDICES_HEADER = ("""
    ! NSPEC for {tag} : {desc}
    integer, public, parameter :: NSPEC_{tag}={count}
    """)
    INDICES_HEADER_WITH_EXTENT = dedent("""
    !+ Defines indices for {tag} : {desc}
    integer, public, parameter :: FIRST_{tag}={first}, &
                                   LAST_{tag}={last}
    """)

# HI took out CiStar and DeltaH
    DECLS = dedent("""
    !/--   Characteristics of species:
    !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)

    public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

    type, public :: Chemical
         character(len=20) :: name
         real              :: molwt
         integer           :: nmhc      ! nmhc (1) or not(0)
         real              :: carbons   ! Carbon-number
         real              :: nitrogens ! Nitrogen-number
         real              :: sulphurs  ! Sulphur-number
    endtype Chemical
    type(Chemical), public, dimension(NSPEC_TOT), target :: species

    ! Pointers to parts of species (e.g. short-lived, advected)
    """)

    DEFINE_HEADER = 'subroutine define_chemicals()\n'

    DEFINE_FOOTER = 'end subroutine define_chemicals\n'

    def __init__(self, scheme):
        self.scheme = scheme

    # added second stream for ChemDims
    def write(self, stream_specs, stream_dims):
        LOG.info('Writing species')
        LOG.indent()

        all_species = self.scheme.get_species_list()

        # Wrap stream in indenting writer, HI added second one
        stream_specs = IndentingStreamWriter(stream_specs)
        stream_dims = IndentingStreamWriter(stream_dims)

        self.write_module_header(stream_specs, 'ChemSpecs_mod',
            # necessary with using CM_ChemDims
            ['use ChemDims_mod      ! => NSPEC_TOT, NCHEMRATES, ....'])

        # Write each defined group of indices
        for info in self.INDEX_GROUPS:
            if info['filter'] is None:
                species = all_species
                stream_dims.write(self.INDICES_HEADER.format(count=len(species), **info))
            else:
                species, first, last = self._find_extent(all_species, info['filter'])
                # separate INDICES_HEADER and INDICES_HEADER_WITH_EXTENT
                stream_dims.write(self.INDICES_HEADER.format(count=len(species), **info))
                stream_specs.write(self.INDICES_HEADER_WITH_EXTENT.format(
                    first=first, last=last, **info))

            LOG.info('PROCESS %s NSPEC %s', info['tag'], len(species))
            self._write_indices(stream_specs, species, info['ixlab'])

        stream_specs.write(self.DECLS)

        for info in self.INDEX_GROUPS:
            if info['filter'] is not None:
                stream_specs.write('type(Chemical), public, dimension(:), pointer :: species_{tag}=>null()\n'
                             .format(tag=info['tag'].lower()))

        self.write_contains(stream_specs)

        stream_specs.write(self.DEFINE_HEADER)
        stream_specs.indent()

        stream_specs.write(dedent("""\
        integer :: istart ! For NAG compliance
        !+
        ! Pointers to parts of species (e.g. short-lived, advected), only assigned if
        ! non-empty.
        !
        """))
        for info in self.INDEX_GROUPS:
            if info['filter'] is not None:
                stream_specs.write(dedent("""\
                istart = 1
                if (NSPEC_{tag_uc} > 0) then
                  if( FIRST_{tag_uc} > 0 ) istart = FIRST_{tag_uc}
                  species_{tag_lc} => species(istart:LAST_{tag_uc})
                end if
                """).format(tag_uc=info['tag'], tag_lc=info['tag'].lower()))

        stream_specs.write(dedent("""
        !+
        ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
        ! array, using indices from total list of species (advected + short-lived).
        !                                                      MW  NM   C    N   S
        """))
        for spec in all_species:
            self._write_spec(stream_specs, spec)

        stream_specs.outdent()
        stream_specs.write(self.DEFINE_FOOTER)

        self.write_module_footer(stream_specs, 'ChemSpecs_mod')

        LOG.outdent()

    def _write_indices(self, stream, species, prefix):
        INDEX_HEADER = '\ninteger, public, parameter :: &\n  '
        INDEX = '{prefix}{spec.name:<12}={i:4d}'
        indices = (INDEX.format(prefix=prefix, i=i, spec=spec) for i, spec in enumerate(species, 1))
        for chunk in ichunk(indices, 10):
            stream.write(INDEX_HEADER)
            stream.write('  &\n  , '.join(chunk) + '\n')

    def _write_spec(self, stream, spec):
        SPEC_FORMAT = ('species({s.name:<12}) = Chemical("{s.name:<12}",{s.molwt:9.4f},'
                       '{nmhc:3d},{s.counts[C]:9.4f},{s.counts[N]:9.4f},{s.counts[S]:9.4f} )\n')
# HI took out cstar and DeltaH variables
        stream.write(SPEC_FORMAT.format(s=spec, nmhc=(1 if spec.NMHC else 0)))

    def _find_extent(self, all_species, predicate):
        """Find ``(species, first_ix, last_ix)`` for species matching *predicate*.

        Assumes *all_species* is ordered in such a way that everything matching
        *predicate* is in one contiguous group.  Returns ``([], -999, -999)``
        if nothing is matched.
        """
        offset = 0
        for check, species in itertools.groupby(all_species, predicate):
            if check:
                species = list(species)
                return (species, offset + 1, offset + len(species))
            else:
                offset += len(list(species))
        else:
            return ([], -999, -999)

    # separated into what is written to CM_ChemDims and CM_XXXDep.inc
    DEPMAP_DECLS_DIMS = dedent("""
    ! No. {kind} deposition species
    integer, public, parameter :: N{kind}DEP_ADV = {count}
    """)
    #!        species     surrogate      adv   idef  setrate
    DEPMAP = dedent("""\
    type(dep_t), public, dimension(N{kind}DEP_ADV), save :: CM_{kind:.1}DepMap = (/ &
    !        species     surrogate      setrate
      {map} &
    /)
    """)

    DEPMAP_SPEC = 'dep_t("{spec:10}", "{other:10}", -999.0 )'

    # add second stream for ChemDims
    def _write_depmap(self, stream, kind, depmap, stream_dims):
        """Write a deposition surrogate map for *depmap* to *stream*."""
        # add CM_ChemDims
        stream_dims = IndentingStreamWriter(stream_dims)
        stream_dims.indent()
        stream_dims.indent()
        stream_dims.write(self.DEPMAP_DECLS_DIMS.format(kind=kind, count=len(depmap)))
        stream.write(self.DEPMAP.format(kind=kind, 
            map=' &\n, '.join((self.DEPMAP_SPEC.format(kind=kind, spec=k, other=v)
                               for k, v in depmap.items()))))

    def write_depmaps(self, dry, wet, stream_dims):
        """Write maps for dry/wet deposition surrogates."""
        self._write_depmap(dry, 'DRY', self.scheme.get_dry_deposition_map(), 
            stream_dims)
        self._write_depmap(wet, 'WET', self.scheme.get_wet_deposition_map(), 
            stream_dims)



class GroupsWriter(CodeGenerator):
    DECLS = dedent("""
    ! Assignment of groups from GenIn_Species.csv
    public :: Init_ChemGroups

    type(typ_sp), dimension({ngroups}), public, save :: chemgroups
    type(typ_factors), dimension({ngroups_factors}), public, save :: chemgroups_factors
    type(typ_maps), dimension({ngroups_maps}), public, save :: chemgroups_maps
    """)

    DECLARE_GROUP = dedent("""
    integer, public, target, save, dimension ({size}) :: &
      {group}_GROUP = (/ {specs} /)
    """)

    # new function for group_factor
    DECLARE_GROUP_FACTOR = dedent("""
    real, public, target, save, dimension ({size}) :: &
      {group}_GROUP_FACTORS = (/ {factors} /)
    """)

    # function for group_map
    DECLARE_GROUP_MAP = dedent("""
    character(len=TXTLEN_SHORT), public, target, save, dimension ({size}) :: &
      {group}_GROUP_MAPBACK = [ character(len=TXTLEN_SHORT) :: &
      {maps} ]
    """)

    DEFINE_GROUP = dedent("""
    chemgroups({i})%name="{group}"
    chemgroups({i})%specs=>{group}_GROUP
    """)

    # function for group_factor
    DEFINE_GROUP_FACTOR = dedent("""
    chemgroups_factors({i})%name="{group}"
    chemgroups_factors({i})%species=>{group}_GROUP
    chemgroups_factors({i})%factors=>{group}_GROUP_FACTORS
    """)

    # function for group_map
    DEFINE_GROUP_MAP = dedent("""
    chemgroups_maps({i})%name="{group}"
    chemgroups_maps({i})%species=>{group}_GROUP
    chemgroups_maps({i})%maps=>{group}_GROUP_MAPBACK
    """)

    def __init__(self, scheme):
        self.scheme = scheme

    def write(self, stream):
        LOG.info('Writing groups')
        LOG.indent()

        groups = self.scheme.get_group_list()
        # HI new groups
        groups_factors = self.scheme.get_group_factor_list()
        # HI if there are no factors groups, define default
        if len(list(groups_factors.keys())) == 0:
            groups_factors['no_factor'] = {'-100': '0.'}
        groups_maps = self.scheme.get_group_map_list()

        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)

        # add new types for maps and factors groups
        self.write_module_header(stream, 'ChemGroups_mod',
                                 ['use ChemSpecs_mod        ! => species indices',
                                  ('use OwnDataTypes_mod  ! => typ_sp, ' +
                                    'typ_factors, typ_maps')])

        stream.write(self.DECLS.format(ngroups=len(groups), \
            ngroups_factors=len(list(groups_factors.keys())), \
            ngroups_maps=len(list(groups_maps.keys()))))

        for g, specs in groups:
            LOG.info('PROCESS %-12s N=%2d  =>  %s', g, len(specs), ', '.join(specs))
            speclist = ', '.join(specs)
            wrapped_speclist = expression_wrap(speclist, 60, ',')
            if len(wrapped_speclist) == 1:
                final_speclist = wrapped_speclist[0]
            else:
                final_speclist = ' &\n    ' + '  &\n    '.join(wrapped_speclist) + '  &\n '
            stream.write(self.DECLARE_GROUP.format(size=len(specs), group=g,
                                                   specs=final_speclist))

        # factor group, now same as maps again
        for g in groups_factors:
            specs = list(groups_factors[g].keys())
            factors = list(groups_factors[g].values())
            LOG.info('PROCESS %-12s N=%2d  =>  %s', g, len(specs), ', '.join(specs))
            speclist = ', '.join(specs)
            wrapped_speclist = expression_wrap(speclist, 60, ',')
            if len(wrapped_speclist) == 1:
                final_speclist = wrapped_speclist[0]
            else:
                final_speclist = (' &\n    ' + '  &\n    '.join(
                    wrapped_speclist) + '  &\n ')
            stream.write(self.DECLARE_GROUP.format(size=len(specs), group=g,
                                                   specs=final_speclist))
            # same for list of factors
            factorlist = ', '.join(factors)
            wrapped_factorlist = expression_wrap(factorlist, 60, ',')
            if len(wrapped_factorlist) == 1:
                final_factorlist = wrapped_factorlist[0]
            else:
                final_factorlist = (' &\n    ' + '  &\n    '.join(
                    wrapped_factorlist) + '  &\n ')
            stream.write(self.DECLARE_GROUP_FACTOR.format(size=len(factors),
                group=g, factors=final_factorlist))

        # map group
        for g in groups_maps:
            specs = list(groups_maps[g].keys())
            maps = list(groups_maps[g].values())
            LOG.info('PROCESS %-12s N=%2d  =>  %s', g, len(specs), ', '.join(specs))
            speclist = ', '.join(specs)
            wrapped_speclist = expression_wrap(speclist, 60, ',')
            if len(wrapped_speclist) == 1:
                final_speclist = wrapped_speclist[0]
            else:
                final_speclist = ' &\n    ' + '  &\n    '.join(wrapped_speclist) + '  &\n '
            stream.write(self.DECLARE_GROUP.format(size=len(specs), group=g,
                                                   specs=final_speclist))
            # HI same for list of maps, now want them to stay strings in
            #  the Fortran code as well!
            maplist = '"' + '", "'.join(maps) + '"'
            wrapped_maplist = expression_wrap(maplist, 60, ',')
            if len(wrapped_maplist) == 1:
                final_maplist = wrapped_maplist[0]
            else:
                final_maplist = ('  &\n    '.join(wrapped_maplist) +'  &\n ')
            stream.write(self.DECLARE_GROUP_MAP.format(size=len(maps),
                group=g, maps=final_maplist))

        self.write_contains(stream)

        stream.write('\nsubroutine Init_ChemGroups()\n')
        stream.indent()
        for i, (g, specs) in enumerate(groups, 1):
            stream.write(self.DEFINE_GROUP.format(i=i, group=g))
        # same for new groups
        for i, g in enumerate(groups_factors, 1):
            stream.write(self.DEFINE_GROUP_FACTOR.format(i=i, group=g))
        for i, g in enumerate(groups_maps, 1):
            stream.write(self.DEFINE_GROUP_MAP.format(i=i, group=g))
        stream.outdent()
        stream.write('\nend subroutine Init_ChemGroups\n')

        self.write_module_footer(stream, 'ChemGroups_mod')

        LOG.outdent()


class ReactionsWriter(CodeGenerator):
    # HI took out emissions from special handling (apart from remembering emission species)
    PHOTOL_REF = 'rcphot({spec},k)'
    COEFF_REF = 'rct({i},k)'
    COEFF_ALL = 'rct({i},:)'
    AMOUNT_REF = 'xnew({spec})'

    CHEMEQN_FULL = 'xnew({spec}) = (xold({spec}) + dt2 * P) / (1.0 + dt2 * L)\n'
    CHEMEQN_NOPROD = 'xnew({spec}) = xold({spec}) / (1.0 + dt2 * L)\n'
    CHEMEQN_NOLOSS = 'xnew({spec}) = xold({spec}) + dt2 * P\n'
    CHEMEQN_NONE = '! Nothing to do for {spec}! xnew({spec}) = max(0.0, xold({spec}))\n'
    # Map (has_prod, has_loss) to an equation string
    CHEMEQN_LOOKUP = {
        (True, True): CHEMEQN_FULL,
        (False, True): CHEMEQN_NOPROD,
        (True, False): CHEMEQN_NOLOSS,
        (False, False): CHEMEQN_NONE,
    }
    # Special case for RO2 pool
    CHEMEQN_RO2POOL = 'xnew(RO2POOL) = sum(xnew(RO2_GROUP))\n'

    # add ChemDims
    CHEMRATES_USE = [
        'use AeroConstants_mod  ! => AERO%PM etc, ...',
        'use AeroFunctions_mod  ! => UptakeRates, ...',
        'use ChemDims_mod       ! => NSPEC_TOT, NCHEMRATES, ....',
        'use ChemFunctions_mod  ! => IUPAC_troe, RiemerN2O5, ....',
        'use ChemSpecs_mod      ! => PINALD, .... for FgasJ08',
        'use DefPhotolysis_mod  ! => IDNO2 etc.',
        'use ZchemData_mod      ! => rct',
    ]

    CHEMRATES_DECLS = dedent("""
    public :: setChemRates
    public :: setPhotolUsed

    ! Photolysis rates
    integer, save, public, dimension(NPHOTOLRATES) :: photol_used
    """)

    CHEMRATES_DECLS_DIMS = dedent("""
    ! No. rate coefficients
    integer, parameter, public :: NCHEMRATES = {coeff_count}

    ! No. photolysis rates used
    integer, parameter, public :: NPHOTOLRATES = {photol_count}
    """)

    def __init__(self, scheme):
        self.scheme = scheme
        self.coefficients = SimpleOrderedSet()

        self.prod = DefaultListOrderedDict()
        self.loss = DefaultListOrderedDict()
        self.emis_specs = SimpleOrderedSet()
        self.photol_specs = SimpleOrderedSet()
        self.coefficients = SimpleOrderedSet()

        LOG.info('Extracting prod/loss from reactions...')
        LOG.indent()

        def process_rate_part(part):
            if isinstance(part, tuple):
                kind, arg = part
                if kind == 'photol':
                    if arg in self.photol_specs:
                        LOG.debug('PHOTOL found: %s', arg)
                    else:
                        LOG.debug('PHOTOL new: %s', arg)
                        self.photol_specs.add(arg)
                    return self.PHOTOL_REF.format(spec=arg)
                elif kind == 'coeff':
                    if arg not in self.coefficients:
                        LOG.debug('RCT NEW %2d: %s', len(self.coefficients), arg)
                    index = self.coefficients.add(arg) + 1
                    return self.COEFF_REF.format(i=index)
                elif kind == 'amount':
                    return self.AMOUNT_REF.format(spec=arg)
                else:
                    raise ValueError('unhandled rate part', part)
            else:
                #  check for emission species here
                if 'RCEMIS' in part:                     # e.g. RCEMIS(NO,k)
                    arg = part[part.find('(')+1:part.find(',')] # e.g. NO
                    #if 'BRO' in part: LOG.debug('DS RCEMIS IN: %s %s'% ( part, arg) )
                    if arg in self.emis_specs:
                        LOG.debug('RCEMIS found: %s', arg)
                    else:
                        self.emis_specs.add(arg)
                        LOG.debug('RCEMIS new: %s', arg)
                    part = part.replace('RCEMIS', 'rcemis')
                return part
                #if 'djphot' in part:                     # e.g. DJPHOT(NO,k)
                #    sys.exit('BRODJ')
                #if 'DJPHOT' in part:                     # e.g. DJPHOT(NO,k)
                #    arg = part[part.find('(')+1:part.find(',')] # e.g. NO
                #    if 'BRO' in part: LOG.debug('DS DJPHOT IN: %s %s'% ( part, arg) )
                #    if arg in self.djphotol_specs:
                #        LOG.debug('DSDJPHOT found: %s', arg)
                #    else:
                #        self.djphotol_specs.add(arg)
                #        LOG.debug('DSDJPHOT new: %s', arg)
                #    part = part.replace('DJPHOT', 'djphot')
                #return part

        for reaction in self.scheme.reactions:
            LOG.info('PROCESS %r', reaction)
            LOG.indent()
            for spec, rate in reaction.get_prod_rates():
                LOG.info('PROD  %-12s: %s', spec, rate)
                LOG.indent()
                processed_rate = ' '.join(process_rate_part(r) for r in rate)
                # add brackets around sum of photolysis rates
                OPS = ['+', '-']
                if 'rcphot(' in processed_rate.lower() and any(s in 
                    processed_rate.split(' * xnew')[0] for s in OPS):
                    where = processed_rate.find(' * xnew')
                    processed_rate = ('(' + processed_rate[:where] + ')' +
                        processed_rate[where:])
                self.prod[spec].append(processed_rate)
                LOG.debug('RATE: %s', processed_rate)
                LOG.outdent()
            for spec, rate in reaction.get_loss_rates():
                LOG.info('LOSS  %-12s: %s', spec, rate)
                LOG.indent()
                processed_rate = ' '.join(process_rate_part(r) for r in rate)
                # add brackets around sum of photolysis rates
                OPS = ['+', '-']
                if 'rcphot(' in processed_rate.lower() and any(s in 
                    processed_rate.split(' * xnew')[0] for s in OPS):
                    where = processed_rate.find(' * xnew')
                    processed_rate = ('(' + processed_rate[:where] + ')' +
                        processed_rate[where:])
                self.loss[spec].append(processed_rate)
                LOG.debug('RATE: %s', processed_rate)
                LOG.outdent()
            LOG.outdent()

        LOG.outdent()

    # taken from _write_indices, for breaking the production/loss terms
    #     every n_continuation_lines'th line
    def _break_prod_loss(self, stream, species, what, prodloss):
        PRODLOSS_HEADER = '\n!-> {species} cont.\n  {what} = {what} + '
        for chunk in ichunk(prodloss, n_continuation_lines):
            stream.write(PRODLOSS_HEADER.format(species=species, what=what))
            stream.write('  &\n  + '.join(chunk) + '\n')

    def write_prod_loss(self, normal, slow):
        LOG.info('Writing prod/loss files...')
        LOG.indent()

        # Wrap streams in indenting writers
        normal = IndentingStreamWriter(normal)
        slow = IndentingStreamWriter(slow)

        for spec in self.scheme.get_species_list():
            if spec.name == "RO2POOL":
                # Special case for RO2 pool pseudo-species
                prod = []
                loss = []
                LOG.info('RO2POOL FOUND')
                eqn = self.CHEMEQN_RO2POOL
            else:
                # All other species
                prod = self.prod[spec.name]
                loss = self.loss[spec.name]
                LOG.info('SPEC %-12s: nprod=%2d, nloss=%2d', spec.name, len(prod), len(loss))
                eqn = self.CHEMEQN_LOOKUP[(bool(prod), bool(loss))].format(spec=spec.name)

            # Choose appropriate output stream
            stream = slow if ( spec.type == Species.SLOW or spec.type == Species.SEMIVOL) else normal

            stream.write('!-> {spec}\n\n'.format(spec=spec.name))
            stream.indent()

            if prod:
                # break production after every n_continuation_lines'th line
                if len(prod) <= n_continuation_lines:
                    stream.write('P = ' + '  &\n  + '.join(prod) + '\n')
                else:
                    stream.write('P = ' + '  &\n  + '.join(
                    	prod[:n_continuation_lines]) + '\n')
                    self._break_prod_loss(stream, spec.name, 'P', 
                    	prod[n_continuation_lines:])
            else:
                stream.write('! P = 0.0\n')
            stream.write('\n')

            if loss:
                # break losses after every n_continuation_lines'th line
                if len(loss) <= n_continuation_lines:
                    stream.write('L = ' + '  &\n  + '.join(loss) + '\n')
                else:
                    stream.write('L = ' + '  &\n  + '.join(
                    	loss[:n_continuation_lines]) + '\n')
                    self._break_prod_loss(stream, spec.name, 'L', 
                    	loss[n_continuation_lines:])
            else:
                stream.write('! L = 0.0\n')
            stream.write('\n')

            stream.write(eqn)

            stream.outdent()
            stream.write('\n\n')

        LOG.outdent()

    # add second stream for ChemDims
    def write_rates(self, stream_reacs, stream_dims):
        LOG.info("Writing rates...")
        LOG.indent()

        stream_reacs = IndentingStreamWriter(stream_reacs)
        stream_dims = IndentingStreamWriter(stream_dims)

        self.write_module_header(stream_reacs, 'ChemRates_mod', self.CHEMRATES_USE)
        stream_reacs.write(self.CHEMRATES_DECLS)
        stream_dims.indent()
        stream_dims.indent()
        stream_dims.write(self.CHEMRATES_DECLS_DIMS.format(
            coeff_count=len(self.coefficients),photol_count=len(self.photol_specs)))
        self.write_contains(stream_reacs)

        stream_reacs.write('\nsubroutine setPhotolUsed()\n')
        # check on number of photolysis reactions
        if len(self.photol_specs) > 0:
            stream_reacs.indent()
            stream_reacs.write('photol_used = (/ &\n    ' +
                         '  &\n  , '.join(self.photol_specs) +
                         '  &\n/)')
            stream_reacs.outdent()
        else:
            stream_reacs.indent()
            stream_reacs.write('! no photolysis reactions in this scheme')
            stream_reacs.outdent()
        stream_reacs.write('\nend subroutine setPhotolUsed\n')

        stream_reacs.write(dedent("""
        subroutine setChemRates()
          !integer, intent(in) :: debug_level

        """))
        stream_reacs.indent()

        for i, rate in enumerate(self.coefficients, 1):
            # LHS of coefficient
            rct = self.COEFF_ALL.format(i=i) + ' = '
            # RHS of coefficient, wrapped to prevent long lines
            if ',' in rate:
                # Always split function-like coefficients at commas
                rate_parts = expression_wrap(rate, None, ',')
            else:
                # Anything else, try to keep to less than 60 chars wide by
                # breaking at closing parentheses
                rate_parts = expression_wrap(rate, 60, ')')
            # Re-join and write the coefficient
            sep = '  &\n' + ' '*(len(rct)-2) + '& '
            stream_reacs.write(rct + sep.join(rate_parts) + '\n')

        stream_reacs.outdent()
        stream_reacs.write('\nend subroutine setChemRates\n')

        self.write_module_footer(stream_reacs, 'ChemRates_mod')

        LOG.outdent()

    # separate declaration of dimension and of array itself
    EMIS_DECLS_DIMS = dedent("""
    ! No. emission {name}s
    integer, parameter, public :: NEMIS_{name} = {count}
    """)
    EMIS = dedent("""\
    character(len={maxlen}), save, dimension(NEMIS_{name}), public :: EMIS_{name} = (/ &
      {contents} &
    /)
    """)
    # if there are no emissions
    EMIS_zero = dedent("""\
    character(len={maxlen}), save, dimension(NEMIS_{name}), public :: EMIS_{name}
    """)

    def _write_emis(self, stream_emis, name, items, stream_dims):
        #  maxlen = 10 to avoid error messages when there are no
        #  emitting species. Also, EMIS_File will be written differently when
        #  there are no emissions.
        # And separation into ChemDims and EMIS_File

        stream_dims = IndentingStreamWriter(stream_dims)
        stream_dims.indent()
        stream_dims.indent()
        if len(items) > 0:
            maxlen = max(len(_) for _ in items)
            fmt = '"{{:{maxlen}}}"'.format(maxlen=maxlen)
            stream_dims.write(self.EMIS_DECLS_DIMS.format(
                name=name, count=len(items)))
            stream_emis.write(self.EMIS.format(
                name=name, maxlen=maxlen,
                contents=' &\n, '.join(fmt.format(_) for _ in items)))
        else:
            maxlen = 10 # just took this without a special reason
            stream_dims.write(self.EMIS_DECLS_DIMS.format(
                name=name, count=len(items)))
            stream_emis.write(self.EMIS_zero.format(
                name=name, maxlen=maxlen))

    def write_emis_files(self, stream_emis, stream_dims):
        """Write emission file array definition to *stream_emis* and
        *stream_dims*."""
        LOG.info('Writing emission file array (%s)...',
                 ', '.join(self.scheme.emis_files))
        self._write_emis(stream_emis, 'File', self.scheme.emis_files, stream_dims)

    def write_emis_specs(self, stream, stream_dims):
        """Write emission species array definition to *stream*."""
        LOG.info('Writing emission species array...')
        LOG.indent()
        LOG.info(', '.join(self.emis_specs))
        self._write_emis(stream, 'Specs', self.emis_specs, stream_dims)
        LOG.outdent()

    # TMP - failed attempt at recode
    #def write_djphot_specs(self, stream):
    #    """Write djphotol species array definition to *stream*."""
    #    LOG.info('Writing djphotol species array...')
    #    LOG.info('DSDJ' + ', '.join(self.djphotol_specs))
    #    stream.write('XXXX')
    #    for f in self.djphotol_specs:
    #        stream.write('{:>10}'.format(f))
    #    #self._write_emis(stream, 'djphotol', self.djphotol_specs, stream_dims)
    #    #LOG.outdent()
    #    stream.write('\n')

    def write_femis(self, stream):
        """Write default emission control options to *stream*."""
        LOG.info('Writing femis.defaults...')
        stream.write('Name  {}  '.format(len(self.scheme.emis_files)))
        for f in self.scheme.emis_files:
            stream.write('{:>10}'.format(f))
        stream.write('\n 28   0  ')
        for f in self.scheme.emis_files:
            stream.write('{:10.1f}'.format(1.0))
        stream.write('\n')

    def write_emis_list(self, stream):
        """Write CSV list of emission files to *stream*."""
        LOG.info('Writing CM_emislist.csv...')
        stream.write(','.join(self.scheme.emis_files) + '\n')


class PrettyStreamHandler(logging.StreamHandler):
    """A :class:`logging.StreamHandler` that wraps log messages with
    severity-dependent ANSI colours."""
    #: Mapping from logging levels to ANSI colours.
    COLOURS = {
        logging.DEBUG: '\033[36m',      # Cyan
        logging.WARNING: '\033[33m',    # Yellow foreground
        logging.ERROR: '\033[31m',      # Red foreground
        logging.CRITICAL: '\033[31;7m'  # Red foreground, inverted
    }
    #: ANSI code for resetting the terminal to default colour.
    COLOUR_END = '\033[0m'

    def format(self, record):
        """Call :meth:`logging.StreamHandler.format`, and apply a colour to the
        message if output stream is a TTY."""
        msg = super(PrettyStreamHandler, self).format(record)
        if self.stream.isatty():
            colour = self.COLOURS.get(record.levelno, '')
            return colour + msg + self.COLOUR_END
        else:
            return msg

def run_genchem(name_string):
    global schemename
    schemename = name_string
    
    formatter = logging.Formatter('[%(levelname).1s%(levelname).1s] %(message)s')
    # Send logging output to stderr
    stream_handler = PrettyStreamHandler()
    #stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    # Send logging output to Log.GenOut
    file_handler = logging.FileHandler('Log.GenOut', 'w')
    file_handler.setFormatter(formatter)
    # Attach log handlers
    rootlogger = logging.getLogger('')
    rootlogger.setLevel(logging.DEBUG)
    rootlogger.addHandler(stream_handler)
    rootlogger.addHandler(file_handler)

    # Hook unhandled exception handling into logging
    def excepthook(type, value, tb):
        import traceback
        for part in traceback.format_exception(type, value, tb):
            for line in part.splitlines():
                LOG.critical(line)
        sys.exit('DS TMP2')
    sys.excepthook = excepthook

    scheme = ChemicalScheme()

    species_reader = SpeciesReader(scheme)

    with open('GenIn_Species.csv', 'r') as f:
        species_reader.read(f)

    shorthand = ShorthandMap(open('GenIn_Shorthands.txt', 'r'))
    reactions_reader = ReactionsReader(scheme, shorthand)
    with open('GenIn_Reactions.txt', 'r') as f:

        with open('CM_Reactions.log', 'w') as f2:
            reactions_reader.read(f, f2)

    # CM_ChemDims, will be written to from all Writer classes
    dimensions_writer = DimensionsWriter()
    f_dims = open('CM_ChemDims_mod.f90', 'a')
    dimensions_writer.write_beginning(f_dims)

    species_writer = SpeciesWriter(scheme)
    with open('CM_ChemSpecs_mod.f90', 'w') as f:
        species_writer.write(f, f_dims)
    with open('CM_DryDep.inc', 'w') as dry, open('CM_WetDep.inc', 'w') as wet:
        species_writer.write_depmaps(dry, wet, f_dims)

    groups_writer = GroupsWriter(scheme)
    groups_writer.write(open('CM_ChemGroups_mod.f90', 'w'))

    reactions_writer = ReactionsWriter(scheme)
    with open('CM_Reactions1.inc', 'w') as f1, open('CM_Reactions2.inc', 'w') as f2:
        reactions_writer.write_prod_loss(f1, f2)
    with open('CM_ChemRates_mod.f90', 'w') as f:
        reactions_writer.write_rates(f, f_dims)
    with open('CM_EmisFile.inc', 'w') as f:
        reactions_writer.write_emis_files(f, f_dims)
    with open('CM_EmisSpecs.inc', 'w') as f:
        reactions_writer.write_emis_specs(f, f_dims)
#    with open('CM_XXXXXXXXX.inc', 'w') as f:
#        reactions_writer.write_djphot_specs(f)
    with open('femis.defaults', 'w') as f:
        reactions_writer.write_femis(f)
    with open('CM_emislist.csv', 'w') as f:
        reactions_writer.write_emis_list(f)

    dimensions_writer.write_end(f_dims)
    f_dims.close()

if __name__ == '__main__':
    # schemename passed to GenChem
    if len(sys.argv) <= 1:
        print("""
        GenChem should be be called with the mechanism(s) as argument, or some helpful text
        E.g.:
        GenChem.py CRI_R5 MCM_v3.3
        """)
        raise ValueError('Try again!')
    #name_string = sys.argv[1]
    #for i in range(2,len(sys.argv)):
    #    name_string += ' ' + sys.argv[i]

    run_genchem(' '.join(sys.argv[1:]))

