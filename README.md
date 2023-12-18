GenChem 
=======

This package contains the GenChem system, which 
includes a chemical pre-processor (**GenChem.py**) for converting chemical
equations into differential form for use in atmospheric chemical transport
models (CTMs).  Although primarily intended for users of the **EMEP MSC-W**
chemical transport model (https://github.com/metno/emep-ctm) and 
related systems, GenChem also features a simple box-model testing system
(**boxChem**), which can run as a stand-alone chemical solver, enabling for
example easy testing of chemical mechanisms against each other.

To get started, read the documentation article and guide:

   https://gmd.copernicus.org/articles/13/6447/2020/

   https://genchem.readthedocs.io/en/latest/index.html


Updates:

  2023-12: Converted system to use CloudJ photolysis, c.f. van Caspel et al. (accepted, in typesetting):

    https://gmd.copernicus.org/preprints/gmd-2023-147/

  2022-12-21: Added EmChem19X, and updated READMEs to refer to NEW report:

    Update and comparison of atmospheric chemistry mechanisms for the EMEP MSC-W model system - EmChem19a, EmChem19X, CRIv2R5Em, CB6r2Em, and MCMv3.3Em,
    Robert Bergström, Garry D. Hayman, Michael E. Jenkin  and David~Simpson,
    Technical Report MSC-W 1/2022,
    EMEP MSC-W, Norwegian Meteorological Institute, Norway
    https://emep.int/publ/reports/2022/MSCW_technical_1_2022.pdf

  2022-06-14: Revised tag: 1.1
  2022-06-13: Added -g --gnfr flags to enable GNFR emissplits, and renamed files, e.g. emissplit_snap_defaults_voc.csv

  2020-12-21: GMD article published! Link above updated 

  2020-10-22: GMD article accepted for publication. 

  2020-09-03:
    Bug fix: changed log to log10 in Shorthands file for CB6r2Em
    changed group name of PANs to PANS in Species file
    Moved chem/scripts/GenIn_Shorthands.txt to chem/generic_Shorthands.txt
     to make this file more visible.
    Added -i, -v args to getboxconcs.py in box/scripts
