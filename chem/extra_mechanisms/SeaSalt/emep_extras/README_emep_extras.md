EMEP extras
===========

These files will be used by the EMEP MSC-W 3-D model. They provide
links between GenChem's chemical species and various externally-provided
boundary conditions or biomass burning compounds. These files have to
created by hand for now, but some hard-coding will always be needed to
map external names to emep species names.

Biomass burning
---------------

Not relevant for sea-salt ;-)


Boumdary condition mapping
--------------------------

emep_BoundaryConditions.txt


Note that in 99% of case EMEP's BCs for sea-salt will be taken from a 3-D netcdf file.
This file just provides backup for the other 1%.
