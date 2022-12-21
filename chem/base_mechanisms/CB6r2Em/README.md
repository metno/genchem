CB6r2Em is a version of the CB6r2h code provided by Greg Yarwood (Dec. 2016),
modified to work within the EMEP/ESX/GenChem modelling systems.

See:

  Update and comparison of atmospheric chemistry mechanisms for the EMEP MSC-W model system - EmChem19a, EmChem19X, CRIv2R5Em, CB6r2Em, and MCMv3.3Em,
  Robert Bergström, Garry D. Hayman, Michael E. Jenkin  and David~Simpson,
  Technical Report MSC-W 1/2022,
  EMEP MSC-W, Norwegian Meteorological Institute, Norway
  https://emep.int/publ/reports/2022/MSCW_technical_1_2022.pdf


Modifications:

   The halogen species have been moved to the extra_mechanisms directory

   Photolysis rates have been modified to use values from the EMEP systems.
   Photolysis rates for CB6 are sometimes specified using the original CB6
   numbering system, e.g. cbphot(1), and sometimes with an EMEP/boxChem
   notation, e.g. rcphot(IDNO3). The CB6r2Em_Shorthands.txt file will
   convert the cbphot variables to the equivalent EMEP/boxChem ecphot
   variables, so that a uniform system can be used across all the
   chemical-schemes in GenChem.

   PARLOSS added as a species, to avoid negative stiochiometry on RHS of
     equations.
   Fast PARLOSS + PAR = ; reaction added to achieve the desired effect in a
     mass consistent way.

   Added H2, NH3, SO4, CH3OOH for consistency with EMEP/ESX

   MW changed in CB6x for those 'artificial' species (PAR,OLE,ALD2) which are
   emitted in emissplit. 

   Added PANS in group column for 3 peroxy acyl nitrate species

   Some species names changed to match MCM/EMEP/CRI. 
      - This has been done for TERP -- name changed to APINENE (still
             representing all monoterpenes!)
      - Also done for ISOP -- name changed to C5H8

