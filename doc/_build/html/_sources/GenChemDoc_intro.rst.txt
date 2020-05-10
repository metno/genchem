.. index:: INTRO
  

Introduction
============

GenChem is a system to generate and test chemical mechanisms for the
EMEP MSC-W model [Simpson2012] and 1-D canopy model, ESX [SimpsonTuovinen2014].
GenChem consists of two main directories, **chem** and **box**.

The **chem** directory contains several chemical mechanisms written
in chemist-friendly format (e.g. *k1*  NO2 + OH = HNO3).
A python script *GenChem.py* can be used to convert these files
to fortran friendly input files for the EMEP model, usually with the help
of some wrapper script, either *do.GenChem*, *do.testChems*, or *emep_setup.py*.
The fortran files produced by these scripts
have the prefix "CM_" or "CMX_", where CM denotes Chemical Mechanism.

Although GenChem can be run directly from within the **chem** directory,
the strongly recommended  approach is to use the scripts available
in the **box** directory. In this approach GenChem is first applied, and
then the resulting CM files are compiled and run
as box-model simulations. Once all looks okay, a final script
can be run to add additional code, and provide an EMEP-ready
set of fortran files. This approach ensures that the CM files
compile as they should, and allows rapid testing of several chemical
mechanisms alongside each other.


Code structure
--------------

The directory structure for GenChem can be summarised as::

  XXX/chem                 # GenChem's mechanism tree
  XXX/chem/scripts         # scripts, including do.GenChem and GenChem.py
  XXX/chem/base_mechanisms # collection of main chemical schemes
  XXX/chem/extra_mechanisms # collection of extra reactions for chemical schemes
  XXX/chem/inputs           # emissplit files, see ...

  XXX/box                 # Top of ESX directory tree
  XXX/box/src             # source files
  XXX/box/scripts         # scripts 

  XXX/doc              # documentation, as .rst files plus sphinx conf system
  XXX/doc/_build       # processed documentation, as .pdf and html 
  XXX/doc/_build/html  #  .. as .html  (aim your browser at index.html here)
  XXX/doc/_build/latex # .. as .pdf  (aim your pdfreader at GenChemDoc.pdf here)

(where XXX could any suitable user-directory into which GenChem was unpacked, e.g. /home/fred/chemwork/GenChem.)


Conventions in documenentation naming
-------------------------------------

The input files to GenChem (GenIn files) as used in box or emep model
are usually built up by appending files from one *base* directory (from
base_mechanisms) and one or more (usually many!) *extra* mechanisms
from the extra_mechanisms directory. For example, GenIn_Species.csv
used for  the EMEP CTM's default EmChem19p scheme consists of  Species
files from base_mechanisms/EmChem19a, and from twelve extra_mechanisms
directories (e.g. extra_mechanisms/SeaSalt/SeaSalt_Species.csv,
extra_mechanisms/PM_VBS_EmChem19/PM_VBS_EmChem19_Species, etc.). To
avoid having to write out these names explicitly each time, we adopt
generic names, as illustrated below for the EmChem19p case::


  SCHEME               name for complete chemical mechanisms package. 
                       (currently EmChem19a, EmChem19p, CB6r2, CRIv2emep, MCM_v3.3)

  BASE_Species.csv     base_mechanisms/EmChem19a_Species.csv

  EXTRAS_Species.csv   extra_mechanisms/SeaSalt/SeaSalt_Species.csv, 
                       extra_mechanisms/Aqueous_EmChem16x/Aqueous_EmChem16x_Species.csv,
                       ....

  CMDIR_Species.csv    Either base or extras file, e.g.
                       base_mechanisms/EmChem19a_Species.csv **or**
                       extra_mechanisms/SeaSalt/SeaSalt_Species.csv, 



Requirements
------------

GenChem has been developed on Ubuntu linux systems, and
should work on any modern linux/unix computer. The code has also been
run on Windows via a virtual ubuntu environment.
The minumum requirements are a modern fortran compiler and python3 
(probably 3.5 or higher).

We have used for example

        * gfortran (gcc 4.6.1) on Linux Xubuntu PC system

        * gfortran (gcc 4.4.7) on HP supercomputer

        * ifort 13.0.1



.. comment::

  **  NOTE !!
  This user-guide is a work-in-progress manual on the GenChem system,
  with this interim version produced for interested users, Feb. 2020.
  **
