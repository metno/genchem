
.. index:: Formatting of GenChem files


..
  COMMENTED
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
  

Formatting of GenChem files
---------------------------

Three types of simple text file are needed to 
describe chemical mechanisms in the GenChem system: Reactions,
Species and Shorthands. The formatting of these is described
briefly here, and further information and background can be
found in
`Simpson et al 2020 <https://gmd.copernicus.org/preprints/gmd-2020-147/>`_.

Reactions files
+++++++++++++++

Example::

  * Some simple lines
  1.4e-12*EXP(-1310.*TINV)          : O3 + NO    = NO2 + <O2>  ; acp2004
  5.681e-34*EXP(-2.6*LogTdiv300)    : OP + <O2> + <M>  = O3 ; acp2004
  2.15e-11*EXP(110.*TINV)           : OD + <N2>  = OP       ; Updated (IUPAC 2009)

  emisfiles:sox,nox,co,voc,nh3
  rcemis(NO,KDIM)                        : = NO  ;


*   END-OF-LINE is ";". Text after this (e.g. references, or unused "products") will be ignored. 
*   Separator between rate coefficient and reaction is ":".
*   lines beginning with "*" are comments (no ";" needed here)
*   lines beginning with "rcemis" are emission terms  
*   lines beginning with "emisfiles" give name of  emission files, e.g. nox
*   Some coefficients are defined in GenIn.shorthand, e.g. TINV, LogTdiv300 
*   Anything else is simply used as the rate coefficient. (Do not add spaces!)  



Four  types of tracers/catalysts/yields are allowed, denoted by different types of parentheses:

 1) e.g. [OH] + VOC -> SOA   will put xnew(OH) into the loss rate of VOC, but will not change the loss rate of OH.

 2) e.g. {O2} + OD -> OP   will ignore the O2 term. Make sure it is in the reaction rate though if needed!

 3) e.g. OP + <O2> + <M> -> O3  will ignore the O2 and M term AND add their concentrations to the reaction rate (multiply it). This system is only used for these "special" species (O2, N2, M) as they must be pre-defined, e.g. O2(k), in boxChem and/or EMEP codes.

 4) e.g. 1.36e-11 :   [OXYL] + [OH] = \|YCOXY(0)\|  ASOC_ug1  + ...  will replace the contents of the || term with yield coefficients which will be updated each time-step in the EMEP model.  These variables (here YCOXY(0)) must be predefined in order for emep\_setup.py and the emep model to compile.



Species files
+++++++++++++++

The input to the GenChem.py script is GenIn\_Species.csv, but this
is assembled by do.GenChem from all needed  \_Species.csv files from
the base_mechanisms and extra_mechanisms sub-directories. For
example, for a typical emep run with base EmChem19a, do.GenChem
appends EmChem19a\_Species.csv, SeaSalt\_Species.csv, and many more into
one GenIn\_Species.csv. The file contains columns with species
name, type, formula, and various settings related to dry and wet deposition


The GenIn\_Species.csv file is a spreadsheet-friendly comma-separated file
where the characteristics of the chemical compounds are given::

  Spec,adv,formula,MW,DRY,WET,Groups,!Comments
  *
  RO2POOL,1,RO2POOL,xx,xx,xx,xx,!
  OD,0,O,xx,xx,xx,xx,!
  NO2,1,NO2,xx,NO2,xx,NOx;OX;OXN;daObs,!
  MACR,1,CH2=CCH3CHO,xx,MEK,xx,RCHO;carbonyl;Hstar_5p0e0;f0_0p05;DRx_2p6,!
  BSOC_ng1e2,2,C,12.,ALD,ROOH,Cstar:0.1;DeltaH:30.0;OM25;PCM;BSOA,"! semi-volatile OC from BVOC "

The meaning of the columns is:


  **Spec** -  Species name as used in model.

  **adv** -   Type of compound. CTMs usually distinguish between advected and
  non-advected (or short-lived) species, in order to minimise CPU needs
  (concentrations of short-lived compounds only need chemical reaction
  terms, not advection). In addition, the EMEP model handles semivolatile
  SOA species  through special handling (see below), and some
  species are so long-lived (e.g. CH4) that they can be accurately
  calculated without multiple iterations.  Allowed values of type are:

    0 - for short lived compounds (e.g. OH), which are not advected in the EMEP model.

    1 - for advected compounds (e.g. O3, HCHO)

    2 - for semivolatile SOA compounds (e.g. BSOC\_ng100). The EMEP model (and boxChem)
    tracks such species by compound rather than phase, and calculates
    the partitioning between the phases dynamically, based upon the
    compound's volatilty. Species labelled with
    type 2 are accounted within the list of advected species, but the
    start and end of the  semivolatile list is calculated by GenChem.py,
    to produce integer variables which demarcate these semivolatile
    compounds, e.g. FIRST_SEMIVOL=136  and LAST_SEMIVOL=176.

    3 - for compounds which react very slowly (e.g. CH4).

 
  **formula** -  If a true chemical formula is provided (e.g. CH3CHO, or 
  O=CHC(O2)(CH3)CH2OH) then GenChem will calculate the number of atoms
  (C, H, O, S or N) and the molecular weight. Such formula must use
  capital letters; lower case letters are ignored as far as processing is
  concerned, but may be used to help document the intention, e.g. nC4H10 
  is identical to C4H10, or pm25 is particulate matter but whose formula
  we do not know. For example, an entry for an organic nitrate might have
  formula 'someNO3' which mixes lower and upper case.  In this case
  the molecular weight must be given if this is needed for the chemical
  modelling. (Typically we do need the mass of emitted species, but not
  always the mass of other species since we usually use mixing ratios
  for advection and output in ppb units.  Occassionally examples occur
  where mass is not strictly required, but where one wants to know
  the nitrogen content, typically where outputs are given in terms of
  e.g. ug(N)/m3. In this case, the 'someNO3' formula would be
  enough to allow GenChem to figure out that this compound contains one
  nitrogen atom.)


  **MW** - can be dummy (xx) or a real number giving the molecular
  weight of the compound. When given, this value is used in place of
  any MW calculated from the formula. As noted above, the MW value is
  sometimes but not always needed. For some emitted compounds, usually
  connected with particulate matter where we do not know the composition,
  we have to give a dummy molecular weight.  This information is used
  internally in the model to get associated mixing ratios, but outputs
  for such compounds should always be in mass-units so that consistency
  is preserved.

  **DRY** -  dry-deposition surrogate. The EMEP and ESX models calculate
  dry-deposition explicitly for a limited number of compounds, and here
  we can choose which of these compounds can be used as a surrogate
  for the desired species.  For example, for O3 we simply use O3; for
  C2H5OOH we use the ROOH surrogate. If not dry-deposited, simply use xx.
  For the semivolatile SOA species EMEP/ESX CTMs will use this rate for
  the gas-phase fraction of the SOA.

  **WET** - wet-deposition surrogate - similar to the dry deposition
  system.  For example, for HCHO we simply use HCHO; for the semivolatile
  SOA species such as BSOC\_ng100  we specify the same wet-deposition
  as for fine-particulate matter (denoted PMf), and the EMEP/ESX CTMs
  will use this rate  for the condensed fraction of the SOA.

  **Groups** -  specifies groups which species belong to (e.g. OXN
  for oxidised nitrogen, RO2 for peroxy radicals) and allows
  surrogate species or factors to be assigned to these groups,
  e.g. Cstar:10.0;Extinc:0.4 assigns a vapour pressure Cstar (used
  in SOA modelling) to be 10 (ug/m3) and an Extinc coefficient to
  be 0.4. It is important that these groups are
  separated by semi-colons, not commas.  This rather powerful feature
  is discussed further in Simpson et al. (Submitted, 2020).


Shorthands file
+++++++++++++++

Shorthands are text-strings used in the Reactions.txt file, usually to represent commonly used rate-coefficients. The meaning of the text-string is given in \_Shorthand.txt file, e.g.  ::

  XT           temp
  FH2O         (1.0+1.4e-21*h2o*exp(2200.0*TINV))
  KHO2RO2      2.91e-13*exp(1300.*TINV) ! MCM2001 ...
  KMT12        IUPAC_troe(2.8e-31*exp(2.6*Log300divT),2.0e-12,exp(-TEMP/472.),M,0.75-1.27*(-TEMP/472.)/LOG(10.))

In these examples, XT is just a character-saving replacement for temp, FH2O gives a more complex expression, which also uses the pre-defined variable TINV = 1/temp. KHO2RO2 is a common rate-coefficient, but here we see that comments are allowed - anything afer the 2nd term. FInally, the KMT12 term shows that complex fuction calls are also allowed. IMPORTANT - avoid white space in any terms!



.. comment::

  **  NOTE !!
  This user-guide is a work-in-progress manual on the GenChem system,
  with this interim version produced for interested users, Feb. 2020.
  **
