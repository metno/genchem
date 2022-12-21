MCMv3.3Em
=========

MCM Chemistry from web-side (http://mcm.leeds.ac.uk/MCMv3.3.1), converted by HI

 - adapted for EMEP/GenChem and with some modifications (updates)

See:

  Update and comparison of atmospheric chemistry mechanisms for the EMEP MSC-W model system - EmChem19a, EmChem19X, CRIv2R5Em, CB6r2Em, and MCMv3.3Em,
  Robert Bergström, Garry D. Hayman, Michael E. Jenkin  and David~Simpson,
  Technical Report MSC-W 1/2022,
  EMEP MSC-W, Norwegian Meteorological Institute, Norway
  https://emep.int/publ/reports/2022/MSCW_technical_1_2022.pdf


Added rcemis terms added from utils/mkrcemis.py

Added some extra species to allow testing from boxChem

Added some deposition params in Species to make ESX work

Chemistry revisions made compared to the MCM v3.3.1 
* 
*       * Fixed one bug in MCM:
*         - M3BU3ECO3 + HO2 -> C45O2 + OH + NO2 in MCM -- the NO2 is not correct (probably should have been CO2?)

*       * Updated a number of reaction rates to be in agreement with IUPAC recommended rates:
*         - the KAPHO2 reaction rate used for RCO3 + HO2 reactions was updated (set in MCM_v3.3_Shorthands.txt)
*         - SO2 + OH -> HSO3
*         - O(1D) + N2 -> O(3P)
*         - 2 NO + O2 -> 2 NO2
*         - C5H8 + NO3 -> NISOPO2
*         - C5H8 + O3 -> 0.5 CH2OOE + 0.5 HCHO + 0.3 MACR + 0.3 MACROOA + 0.2 MVK + 0.2 MVKOOA
*         - APINENE + O3 -> 0.6 APINOOA + 0.4 APINOOB
*         - APINENE + OH -> 0.572 APINAO2 + 0.353 APINBO2 + 0.075 APINCO2
*         - BPINENE + O3 -> 0.4 NOPINONE + 0.4 CH2OOF + 0.6 NOPINOOA + 0.6 HCHO
*         - BPINENE + OH -> 0.849 BPINAO2 + 0.076 BPINBO2 + 0.075 BPINCO2
*         - LIMONENE + O3 -> 0.73 LIMOOA + 0.27 LIMOOB
*         - LIMONENE + OH -> 0.408 LIMAO2 + 0.222 LIMBO2 + 0.37 LIMCO2
*         - GLYOX + NO3 -> HCOCO + HNO3
*         - MGLYOX + NO3 -> CH3CO3 + CO + HNO3
*         - CH3COCO2H + OH -> CH3CO3
*
*         - HCHO + NO3 -> HNO3 + CO + HO2 [introduced temperature dependence]
*         - CH3COCO2H photolysis -- changed products based on IUPAC
*
*       * Updated both reaction rate and products to IUPAC recommendations for one reaction:
*         - CH3CO3H + OH -> 0.5 CH3CO3 + 0.5 CH2O2CO3H [MCM only has CH3CO3 product]
*       * The CH2O2CO3H (from the above reaction) is a new peroxy radical species (not in MCM); two reactions were added for this:
*         - CH2O2CO3H + NO = NO2 + HOCH2CO3
*         - CH2O2CO3H + NO3 = NO2 + HOCH2CO3
*
*
* NOTE! The updates of the reaction rates introduced here are
*       NOT intended to be a complete revision of the MCMv3.3.1 rates!
*       Only a limited set of reactions were checked against IUPAC
*       [mostly reactions also included in the EmChem19 mechanism].
*
* Note2: some other (unintentional) differences may also exist between 
*        the EMEP MCM scheme and the MCM v3.3.1 if the latter was updated 
*        after the original downloading of the scheme.

===
