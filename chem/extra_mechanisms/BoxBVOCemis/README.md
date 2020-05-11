Biogenic/natural emissions
--------------------------

For boxChem modelling we use the simple variable SUN as driver, which
provides a reasonably realistic diurnal variation (epeak of 1.0 at noon,
zero at night), without the complicated light and temperature dependencies
used in full CTM modelling.

As example emission rates, we use the default broadleaf and coniferous
from Simpson et al, ACP, 2012:

  - broadleaf Iso with e=5ug/g/h, D=320 g/m2, Hmix=1000m has emission in #/cm3/s of 3.94e6

  - coniferous MT with emtp=1ug/g/h emtl=1ug/h, D=500 g/m2, Hmix=1000m has emission in #/cm3/s of 6.4e5

In addition, we allow factors fIsop, fMTL, fMTP and fSQT (set in config_emep.nml) which can modify these rates (eg to zero).

Results in code::
```

  4.0e6*SUN*fIsop : =  C5H8  ;
  3.0e5*SUN*fMTL  : =  APINENE  ; 
  3.0e5 *fMTP     : =  APINENE  ; 
```
