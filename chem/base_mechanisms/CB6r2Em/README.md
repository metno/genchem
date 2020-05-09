CB6r2Em is a version of the CB6r2h code provided by Greg Yarwood, modified to
work within the EMEP/ESX/GenChem modelling systems. working.


Modifications:

   The halogen species have been moved to the extra_mechanisms directory

   Photolysis rates have bene modified to use values from the EMEP systems.

   Now just uses CB6 photolysis rate numbering. No link to MCM. (Achieved with flag
   USES%BaseChem = .... in config_box.nml. If CB6, uses new rates)

   PARLOSS added as a species, to avoid negative stiochiometry on RHS of equations.
   Fast PARLOSS + PAR = ; reaction added to achieve the desired effect in a mass
   consistent way.

   Added H2, NH3, SO4, CH3OOH for consistency with EMEP/ESX

   MW changed in CB6x for those 'artificial' species (PAR,OLE,ALD2) which are emitted in 
   emissplit. Used CB05 settings.

   Added PAN in group column for 3 species

Queries:

   SULF is said to be gaseous. Why no particulate SO4? How to treat this, since EMEP
   usually has SO4 = particulate

   Should we change species names to match MCM/EMEP. That would simplify and enable some
   of the CRI species used for emissions to have the same name in CB6 also.

