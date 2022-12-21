CRIv2R5Em
=========

Major update of the EMEP-CRI_v2.1_R5 scheme -- many reaction rates updated 
(to new IUPAC recommendations, and largely in line with CRIv2.2 and MCM);
isoprene chemistry completely revised (Jenkin et al., 2019; The CRI v2.2 
reduced degradation scheme for isoprene, Atmospheric Environment 212 (2019) 
172–182).
 
Similar to the gas-phase scheme used for global modelling in:

   McFiggans et al, Secondary organic aerosol reduced by mixture of
   atmospheric vapours, Nature, 565, 587-593, doi:10.1029/2018JD029133,
   2019

From Robert Bergström,  28 Nov 2016, adapted by Hannah & Dave for BoxChem
[includes later revisions to some reactions after 2016]

See:

  Update and comparison of atmospheric chemistry mechanisms for the EMEP MSC-W model system - EmChem19a, EmChem19X, CRIv2R5Em, CB6r2Em, and MCMv3.3Em,
  Robert Bergström, Garry D. Hayman, Michael E. Jenkin  and David~Simpson,
  Technical Report MSC-W 1/2022,
  EMEP MSC-W, Norwegian Meteorological Institute, Norway
  https://emep.int/publ/reports/2022/MSCW_technical_1_2022.pdf


Notes from RB:
--------------

base reactions for the merged CRI_v2_R5+isoprene_C15_R1A chemical mechanism

Files taken from rv3_8_5uk, MV version, Aug 9 2011

RB: Merging scheme with C15_R1A scheme for isoprene
    Assume CH3COO2 [CH3CO3] chemistry can be treated as in the CRI_v2_R5
      scheme (somewhat simplified compared to MCM)
     also keep CRI scheme for HOCH2CO3, CH3CO2H, CH3CHO (note different rate
      coefficient in MCM!?), 
     CH3O2 (some rates different in MCM?), CH3OOH, HCHO, 

*NOTABLE* changes include using the MCM species HCC7CO instead of the CRI
      species CARB17!
         for (non-photolysis) reactions involving HCC7CO the rates have been
      taken from MCM rather than CRI
         MCM has much faster reactions for this species than CRI.
         The RN18O2 and RN18AO2 radicals (and their products RN18OOH and
      RN18NO3) are removed based on MCM-chemistry.

For a number of species the names from the isoprene scheme (C15_R1A) were
      used instead of the CRI names
    e.g. CARB11A -> MEK; CARB7 -> ACETOL; CARB6 -> MGLYOX; CARB3 -> GLYOX; 

RB -- Update August 2017:
   The reaction rates were not updated to the latest recommendations here.
   Fixing this now.
   Also changing some photolysis rates to use the same photolysis as MCM
   Changing the name of CH3COO2 to CH3CO3 (to be the same as in MCM)


Updates Jan-Feb 2019 [RB]: 
   Some changes to groups: Changing name of PAN group to PANS
                           Adding some species to OXN group
                           Removing OXN groups for some radicals (to
                           avoid problems with EMEP code that at the
                           moment does not expect radicals in the groups)

   Switching to latest IUPAC recommended rates for a number of reactions
