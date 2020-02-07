! <BioNatConstants_mod.f90 - specify defined biogenic/natural constants         >!
!*****************************************************************************! 

module BioNatConstants_mod
  implicit none
  private

  ! Constants as used in chemical scheme, eg rcbio(NATBIO%TERP)
  type, private :: natbio_t
    integer :: C5H8 = 1
    integer :: TERP = 2
    integer :: Nrcbio = 2  ! No. of rcbio defined in ChemFields/Biogenics_mod
    integer :: NO   = 3    ! used for EmisNat etc
    integer :: NH3  = 4
    integer :: Rn222  = 5
  end type natbio_t
  type(natbio_t), public, parameter :: NATBIO = natbio_t()

  ! List of biogenic/natural emissions to be considered
  integer, parameter, public ::  NEMIS_BioNat  = 5
  type, private :: bio_t
    character(len=11) :: name
    real :: MW   ! molwt assumed in emission routines
  end type bio_t

  type(bio_t), save, dimension(NEMIS_BioNat), public:: &
      EMIS_BioNat =  (/ &
         bio_t(  "C5H8       ",  68.0 )  &
        ,bio_t(  "MTERP      ", 136.0 )  &
        ,bio_t(  "NO         ",  14.0 )  &
        ,bio_t(  "NH3        ",  17.0 )  &
        ,bio_t(  "RN222      ", 222.0 )  &
   /)
 ! And to match the above array:
 ! integer, public, parameter ::  BIO_C5H8=1,   BIO_MTERP=2, BIO_NO=3 ,  BIO_RN222=4 !NH3

end module BioNatConstants_mod
