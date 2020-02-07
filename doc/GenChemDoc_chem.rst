


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
the base_mechanisms and extra_mechanisms sub-directories. For For
example, for a typical emep run with base EmChem19a, do.GenChem
appends EmChem19a\_Species.csv, SeaSalt\_Species.csv, and many more into
one GenIn\_Species.csv. The file contains columns with species
name, type, formula, and various settings related to dry and wet deposition

Example lines::

  Spec,adv,formula,MW,DRY,WET,Groups,!Comments
  *
  RO2POOL,1,RO2POOL,xx,xx,xx,xx,!
  OD,0,O,xx,xx,xx,xx,!
  NO2,1,NO2,xx,NO2,xx,NOx;OX;OXN;daObs,!
  MACR,1,CH2=CCH3CHO,xx,MEK,xx,RCHO;carbonyl;Hstar_5p0e0;f0_0p05;DRx_2p6,!
  BSOC_ng1e2,2,C,12.,ALD,ROOH,Cstar:0.1;DeltaH:30.0;OM25;PCM;BSOA,"! semi-volatile OC from BVOC "



Further details can be found in the documentation article ...


Shorthands file
+++++++++++++++

Shorthands are text-strings used in the Reactions.txt file, usually to represent commonly used rate-coefficients. The meaning of the text-string is given in \_Shorthand.txt file, e.g.  ::

  XT           temp
  FH2O         (1.0+1.4e-21*h2o*exp(2200.0*TINV))
  KHO2RO2      2.91e-13*exp(1300.*TINV) ! MCM2001 ...
  KMT12        IUPAC_troe(2.8e-31*exp(2.6*Log300divT),2.0e-12,exp(-TEMP/472.),M,0.75-1.27*(-TEMP/472.)/LOG(10.))

In these examples, XT is just a character-saving replacement for temp, FH2O gives a more complex expression, which also uses the pre-defined variable TINV = 1/temp. KHO2RO2 is a common rate-coefficient, but here we see that comments are allowed - anything afer the 2nd term. FInally, the KMT12 term shows that complex fuction calls are also allowed. IMPORTANT - avoid white space in any terms!


.. warning

  check log300divT, H2O, ec-aging



.. warning::

  **  NOTE !!
  This user-guide is a work-in-progress manual on the GenChem system,
  with this interim version produced for interested users, Feb. 2020.
  **
