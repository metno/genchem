!> Collection of 1-d (z) data used in chemical calculations, eg meteorology,
!! emission and chemical rates.

module ZchemData_mod
  use AeroConstants_mod, only: NSAREA_DEF
  use ChemDims_mod,      only: NSPEC_TOT
  use DefPhotolysis_mod, only: MAXRCPHOT
  use NumberConstants,   only: UNDEF_R
  implicit none
  private
  integer, private, parameter :: NZ = 1 ! No. vertical levels. 1 for box
  real, public, save, allocatable :: rct(:,:)
  real, public, save, dimension(NSPEC_TOT,NZ) ::&
     rcbio  = 0.0 &
    ,rcemis = 0.0
  real, public, save, dimension(MAXRCPHOT,NZ) ::&
     rcphot = 0.0 

  real, public, save, dimension(NSPEC_TOT,NZ) ::&
       xChem = 0.0  &
      ,Dchem = 0.0
  ! Some defaults   
  ! Use UNDEF_R to enforce proper initialisation
  real, public, save, dimension(NZ) :: &
            temp = UNDEF_R, rh = UNDEF_R, tinv, M = UNDEF_R, N2 = UNDEF_R &
            , O2 = UNDEF_R, H2O = UNDEF_R &
            ,log300divt, logtdiv300 &
            ,cN2O5, cHNO3, cHO2, cNO3, cNO2, cO3 & !for gas-aerosols
           !For consistency with EMEP and SOA:
            ,gamN2O5=UNDEF_R & ! for n2o5Hydrol
            ,aero_fom=UNDEF_R, aero_fss=UNDEF_R, aero_fdust=UNDEF_R &
            ,aero_fbc=UNDEF_R &! fractions
            ,xSO4, xNO3, xNH4 ! for Riemer
  real, public, save, dimension(NSAREA_DEF,NZ) :: S_m2m3    ! surface area
  integer, public, save, dimension(NZ) ::  itemp !e.g. for SOA tables
  ! variables for SOA_mod
  real, public, save, dimension(NSPEC_TOT,NZ) :: Fgas, Fpart
  !For SOA: add 1-cell standard atmosphere here:
  !CHECK FIX - needs to be in EMEP too...
  !CHECK here we assume a cell at ca. 950 hPa
  real, public, save, dimension(NZ) :: StandardP_Pa = 0.95*1.01325e5

!    real, public, allocatable, dimension(:), save :: &
!       deltaZcm             & ! layer thickness, cm
!      ,pp                     !pressure
!         ,ugdryPM             & ! for wet radius from Gerber, etc.
!   integer, public, allocatable, dimension(:), save :: &
!       itemp                  ! int of temperature
end module ZchemData_mod
!TSTEMX program tester
!TSTEMX use ZchemData_mod
!TSTEMX implicit none
!TSTEMX 
!TSTEMX end program tester

