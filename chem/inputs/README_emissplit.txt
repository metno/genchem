Emissplit files

The emissplit.defaults files are needed when running in e.g. Lagrangian
mode, or generally with anthropogenic emissions that need to to split
into species, for example NOx into NO and NO2, or (most varoable) NMVOC
into C2H4, nC4H10, etc.

For sox, nox, co and nh3 the same files can be used for all
mechanisms. For VOC the file is copied (via do.GenChem) from the relevant
chem/ files.

This directory has two main sub-directories:

   emissplit_defaults

     Here we store the (usually) unchanging  emisslit files for NOx, SOx, NH3, ..

   emissplit_run

     Here we copy the files need to run the ESX model. As well as the
     files from the defaults directory, the script do.GenChem will
     copy the appopriate emissplit.defaults.voc file for the chemical
     mechanism in use.

