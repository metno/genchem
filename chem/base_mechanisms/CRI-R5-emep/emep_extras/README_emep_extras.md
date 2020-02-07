EMEP extras
===========

These files will be used by the EMEP MSC-W 3-D model. They provide links between GenChem's chemical
species and various externally-provided boundary conditions or biomass burning compounds. These
files have to created by hand for now. Much of this can hopefully be replaced by a namelist
config system one day, though some hard-coding will always be needed to map external names
to emep species names.

Biomass burning
---------------

So far two setups are needed, for either the FINN or GFAS sources of wildfires:

CRI-R5-emep_BiomassBurning_FINNv1.5.txt
CRI-R5-emep_BiomassBurning_GFASv1.txt


Boumdary condition mapping
--------------------------

CRI-R5-emep_BoundaryConditions.txt
