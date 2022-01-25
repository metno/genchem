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

  2020-12-21: GMD article published! Link above updated 

  2020-10-22: GMD article accepted for publication. 

  2020-09-03:
    Bug fix: changed log to log10 in Shorthands file for CB6r2Em
    changed group name of PANs to PANS in Species file
    Moved chem/scripts/GenIn_Shorthands.txt to chem/generic_Shorthands.txt
     to make this file more visible.
    Added -i, -v args to getboxconcs.py in box/scripts
