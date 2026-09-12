module Aqueous_mod
  ! Aqueous reaction rates for usage in gas-phase chemistry:
  ! FAKE !!! Just to allow compilation of chem schemes if emep_setup.sh
  ! used 

integer, private, parameter :: &
  MAXK = 200,   & ! Fake, EMEP usually has 20
  NAQUEOUS = 5, & ! No. aqueous rates
  NAQRC    = 3    ! No. constant rates

!real, public, save,allocatable, dimension(:,:) :: aqrck
real, public, save, dimension(NAQUEOUS,MAXK) :: aqrck = 0.0

integer, public, parameter :: &
  ICLOHSO2  = 1, & ! for [oh] + [so2]
  ICLRC1    = 2, & ! for [h2o2] + [so2]
  ICLRC2    = 3, & ! for [o3] + [so2]
  ICLRC3    = 4, & ! for [o3] + [o2] (Fe catalytic)
  ICLHO2H2O2 = 5    ! for HO2g --> 0.5 * [H2O2]


end module Aqueous_mod
