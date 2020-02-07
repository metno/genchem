!*****************************************************************************! 
      module PhysicalConstants_mod
!----------------------------------------------------------------------------
!  Defines Physical constants 
!----------------------------------------------------------------------------
implicit none
!F private

!-- contains no subroutine:

!

  real , public, parameter ::         &
    AVOG   = 6.023e23                 & ! Avogadros number
  , ATWAIR = 28.964                   & ! mol wt of air, g/mol
  , RGAS_ATML = 0.08205               & ! Molar Gas constant (atm M-1 K-1)
  , RGAS_KG   = 287.0                 & ! Molar Gas constant (J K-1 kg-1)
  , RGAS_J    = 8.3144                  ! Molar Gas constant (J mol-1 K-1)

                                        ! NB. ( J = N m2 = kg m2 s-2 )
                                        !       M = mol l-1

  real, public, parameter  ::    &
       GRAV    = 9.807           &   ! Gravity, m s-2
    ,  EARTH_RADIUS = 6.37e6     &   ! 
    ,  CP      = 1004.0          &   ! Specific heat at const. pressure
    ,  KAPPA   = RGAS_KG/CP      &   
    ,  KARMAN  = 0.41            &   ! Von Karman  (=0.35 elsehwere in code!)
    ,  PI      = 4.0*atan(1.0)   &   ! most exact
    ,  DEG2RAD = PI/180.0        &   ! COnverts degrees to radians
    ,  RAD2DEG = 180.0/PI        &   ! COnverts radians to degrees
    ,  ROWATER = 1000.0          &   ! pw density of water kg m-3
    ,  BOLTZMANN = 1.380e-23     &   ! Boltzmann'c constant[J/deg/molec]
    ,  FREEPATH  = 6.5e-8        &   ! Mean Free Path of air [m]
    ,  VISCO     = 1.46e-5           ! Air viscosity [m2/s]   (was NU)


! Converts from mol/cm3 to nmole/m3
  real, public, parameter :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  

! Some definitions for daylight, in terms of zenith angle and cos(zen):
! (calculated from criteria that DAY_COSZEN > 1.0e-10 as daytime)

!ESX   real, public, parameter  ::  &
!ESX        DAY_ZEN   = 89.9999999942704 & !
!ESX      ,DAY_COSZEN = 1.0e-10

!=================== DEP CODE ================================Y:0

  ! CHARNOCK is used to calculate the roughness length for the 
  ! landuse category water

   real, public, parameter  :: &
       PRANDTL = 0.71,            &   ! Prandtl number (see Garratt, 1992)
       Sc_H2O  = 0.6,             &   ! Schmidt number for water
    CHARNOCK = 0.0144  !  From Garratt for k=0.41
    !CHARNOCK = 0.032   ! Charnock's alpha:
                       ! see Nordeng (1986), p.31, 
                       ! Nordeng(1991), JGR, 96, no. C4, pp. 7167-7174.
                       ! In the second of these publications, Nordeng uses
                       ! "m" to denote Charnock's alpha whilst in the first
                       ! he specifies the value 0.032.

  ! Standard temperature :

  real, public, parameter :: T0 = 273.15   ! zero degrees Celsius in Kelvin 

!=============================================================Y:0


end  module PhysicalConstants_mod
