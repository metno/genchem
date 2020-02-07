module Aqueous_mod
  ! Aqueous reaction rates for usage in gas-phase chemistry:
  ! FAKE !!! Just to allow compilation of chem schemes if emep_setup.sh
  ! used 

integer, private, parameter :: &
  MAXK = 200,   & ! Fake, EMEP usually has 20
  NAQUEOUS = 4, & ! No. aqueous rates
  NAQRC    = 3    ! No. constant rates

!real, public, save,allocatable, dimension(:,:) :: aqrck
real, public, save, dimension(NAQUEOUS,MAXK) :: aqrck = 0.0

integer, public, parameter :: &
  ICLOHSO2  = 1, & ! for [oh] + [so2]
  ICLRC1    = 2, & ! for [h2o2] + [so2]
  ICLRC2    = 3, & ! for [o3] + [so2]
  ICLRC3    = 4    ! for [o3] + [o2] (Fe catalytic)

end module Aqueous_mod
