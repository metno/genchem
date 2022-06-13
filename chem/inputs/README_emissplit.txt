Emissplit files

The emissplit.defaults files are needed when running in e.g. Lagrangian
mode, or generally with anthropogenic emissions that need to to split
into species, for example NOx into NO and NO2, or (most variable) NMVOC
into C2H4, nC4H10, etc.

For sox, nox, co and nh3 the same files can be used for all
mechanisms. For VOC the file is copied (via do.GenChem) from the relevant
chem/ files.

This directory has two sub-directories:

   emissplits_gnfr_defaults  (NEW Jun 2022)

     Here we store the (usually) unchanging  emisslit files for NOx, SOx, NH3, ..
     formatted for the 19-sector GNFR_CAMS system used in EMEP

   emissplit_snap_defaults (DEPRECATED Jun 2022, was emissplit_defaults)

     Here we store the (usually) unchanging  emisslit files for NOx, SOx, NH3, ..
     formatted for the 11-sector SNAP system previously used in EMEP
