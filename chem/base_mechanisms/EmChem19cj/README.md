EmChem19cj chemistry scheme -- the default gas-phase mechanism of the EMEP MSC-W chemical transport model.

The EmChem19cj mechanism is a minor update of EmChem19, reflect the addition of the three explicit
glyox photolysis channels, as well as now-fixed (i.e. time-invariant) background CH4 and H2 concentrations. 
In earlier versions, these species were 'free running' after initialization, causing concentrations to go down due to loss against
oxidation over the course of a simulation (mostly impacting the global EMEP configuration).

EmChem19a and earlier versions are described below:

  Update and comparison of atmospheric chemistry mechanisms for the EMEP MSC-W model system - EmChem19a, EmChem19X, CRIv2R5Em, CB6r2Em, and MCMv3.3Em,
  Robert Bergström, Garry D. Hayman, Michael E. Jenkin  and David~Simpson,
  Technical Report MSC-W 1/2022,
  EMEP MSC-W, Norwegian Meteorological Institute, Norway
  https://emep.int/publ/reports/2022/MSCW_technical_1_2022.pdf


EmChem19a is an updated (bug fixed) version of EmChem19.

So far, the following things have been updated:

Deposition parameters added for C54CO and C5134CO2OH -- both are likely
highly soluble and are expected to deposit efficiently

For EmChem19a the BVOC_EmChem19 should be used for the monoterpene
chemistry.

Replaced PRRO2H by ACETOL (to save one advected species).


EmChem19 chemistry scheme was committed June 2019 for EMEP Report runs.

EmChem19 is an updated version of EmChem16x

  Simplified NO3 reactions for C5H8 and C3H6.

  Isoprene chemistry is largely based on the CheT2 scheme 

  Aromatic chemistry has been changed a little (but is still extremely limited)
    -- BENZENE and TOLUENE are added (only first reaction)

An extended version of EmChem19a, EmChem19X, which gives better night
time RO2 agreement with MCM is under development, but this needs more
reactions.


