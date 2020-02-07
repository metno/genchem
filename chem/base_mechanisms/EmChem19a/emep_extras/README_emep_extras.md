EMEP extras
===========

These files will be used by the EMEP MSC-W 3-D model. They provide
links between GenChem's chemical species and various externally-provided
boundary conditions or biomass burning compounds. These files have to
created by hand for now, but some hard-coding will always be needed to
map external names to emep species names.

Biomass burning
---------------

So far two setups are needed, for either the FINN or GFAS sources of
wildfires:

EmChem19a_BiomassBurning_FINNv1.5a.txt
EmChem19a_BiomassBurning_GFASv1.txt

These files address emissions of the EmChem19a gas-phase species. BB
emissions of particles can be found in the relevant extra_mechanisms
directory, e.g.  FFireInert.


Boumdary condition mapping
--------------------------

EmChem19a_BoundaryConditions.txt
