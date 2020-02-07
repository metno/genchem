GenChem chem directory
======================

The _chem_ directory contains:

 * README_chem.md - this file

 * base_mechanisms

 * extra_mechanisms

 * scripts

 * inputs


The heart of the system is the script _GenChem.py_, which reads chemical
equations and converts them into fortran code for use in
the EMEP system of box, 1-D and 3-D codes. This is usually run via the
helper script _do.GenChem_, or even better from the box-model system and
_box/scripts/do.testChems_. (See box directory for more info, or read
the UsersGuide!)

The _base_ and _extra_ mechanism directories contain a number of  mechanisms
pre-prepared for GenChem usage. More information about the mechanism and 
role for GenChem can be found in the README files in each directory. For
each base\_mechanism we require three files, XXX_Species.csv, XXX_Reactions.txt
and XXX_Shorthands.txt. For extra\_mechanisms we require just XXX_Species.csv and XXX_Reactions.txt.


Formatting
-----------

Reactions files
...............


*   END-OF-LINE is ";" (not required in comment section)
*   Separator between rate coefficient and reaction is ":" * 
*   lines beginning with "*" are comments  
*   lines beginning with "rcemis" are emission terms  
*   lines beginning with "emisfiles" give name of  emission files, e.g. nox

        emisfile: nox, sox, voc


*   Some coefficients are defined in GenIn.shorthand, e.g. KHO2RO2 
*   Anything else is simply used as the rate coefficient.  

GenChem checks whether all shorthands used are also in
GenIn.species. For now, I have put some unused "products"
after the ";" EOL, e.g.HCOOH. 


Tracers, catalyists
...................

* Three  types of tracers allowed:

 1)   e.g. [OH] + VOC -> SOA   will put xnew(OH) into the loss rate
  of VOC, but will not change the loss rate of OH.

 2)   e.g. {O2} + OD -> OP   will ignore the O2 term. Make sure
 it is in the reaction rate though if needed!

 3)   e.g. OP + <O2> + <M> -> O3  will ignore the O2 and M term AND
  add their concentrations to the reaction rate (multiply it) -> it
  should not be there in the first place! And these species must be
  defined as own variable, e.g. O2(k) in the model itself



Base and Extra mechanisms
-------------------------

Base mechanisms
+++++++++++++++

A collection of chemical schemes adapted for GenChem usage. The schemes in
base_mechanisms usually consist of inorganic and hydrocarbon chemistry,
and for each mechanism we need the _Species.csv, _Reactions.txt and
_Shorthands.txt files.  Various reactions which can be added to these
schemes (e.g. gas-aerosol reactions) can be found in the extra_mechanisms
directory

Extra mechanisms
++++++++++++++++

A collection of chemical schemes adapted for GenChem usage. The schemes in
extra_mechanisms usually consist of reactions which can be appended to
the base_mechanisms. These extra mechanisms provide only 
 _Species.csv and _Reactions.txt files, not  the _Shorthands.txt files which
are required for base-mechanisms. 

### ShipNOx

Adds the special SHIPNOx species which is needed for global EMEP runs to
prevent over-efficient O3 generation.

### BVOC_MTERP1_EmChem19a

BVOC reactions in addition to EmChem19 basics

### BVOC_EmChem16mtx

Deprecated 

### PM_WoodFFuelInert

Provides emissions of Species such as POM_f_ffuel, EC_f_wood_age, remPPM25
and pSO4f (for those schemes that can use this). These PM compounds
are used in schemes which assume inert POM emissions (ie standard
EMEP policy runs).

### PM_VBS_EmChem19

This mechanism provides organic aerosol reactions as used in the
standard EMEP model of Simpson et al., Atmos. Chem. Physics, 2012 for
the EmChem09soa case - the latter is essentially produced by
doGenChem.py -b EmChem09 -e VBS_acp2012.

VBS denotes volatility basis set, from the work of Donahue, Robinson etc.,
and following the EMEP implementations documented in Bergström et al.,
Atmos. Chem. Physics, 2012. (The EmChem09soa case uses inert emissions
of promary organic aerosol.)


### BoxAero

Some simple aerosol reactions for box-model studies only.

### Aero2017nx

Aerosol reactions for emep and esx

### Dust

As in EMEP

### SeaSalt

As in EMEP

### FFireInert

As in EMEP


This file?
----------

This file is written in markdown, and can be converted to other formats 
with e.g. pandoc::

   pandoc README_chem.md -o tmp.pdf
