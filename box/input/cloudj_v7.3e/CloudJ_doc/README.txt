-----------------

22-12-2023 Willem:

This directory contains a CloudJ userguide (CloudJ_userguide_for_EMEP.pdf) for the EMEP and GenChem implementation of the 
CloudJ v7.3e photolysis code, as described in van Caspel et al. (2023) (https://doi.org/10.5194/gmd-16-7433-2023).

The userguide describes the usage and a few details of the CloudJ code. In addition, it describes the use of
the scripts (contained in this directory) which can be used to generate new cross-section input data. Such input
data can in turn be used to add new photolysis reaction rates in the EMEP/GenChem modeling systems.

- The cross-sections subdirectory contains pratmo-code fortran scripts to generate cross-section data. In addtion,
  the directory includes a few subdirectories containing data for photolysis reactions that were newly added to the
  EMEP and boxChem models.
- The stand-alone directory contains the stand-alone CloudJ v7.3e code that is shipped along with the CloudJ code 
  distribution, which is also used in practically unmodified form in the boxChem implementation of CloudJ. 

-----------------

