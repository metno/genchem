! <MetFunctions.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
module MetFunctions
!____________________________________________________________________
! Miscellaneous collection of "standard" micromet functions
! Code based upon Garrett by David Simpson & Juha-Pekka Tuovinen, also
! using Stull - Jalle Hiltunen, Chalmers, 2014
! Some code from MET Norway NWP codes
!
! IMPORTANT!! - these functions are often simplified, and intended for near-
! surface application. They should not be applied for e.g. low pressure
!____________________________________________________________________
!   Language: F
!____________________________________________________________________
  use PhysicalConstants_mod, only : CP, PI, RGAS_J, RGAS_KG, AVOG, KAPPA, KARMAN, GRAV
  implicit none
  private

 !/-- Met routines - troposphere

  public :: StandardAtmos_kPa_2_km   ! US Standard Atmosphere conversion
  public :: StandardAtmos_km_2_kPa   ! US Standard Atmosphere conversion

  public :: Tpot_2_T        ! (p/P0)**KAPPA
  public :: T_2_Tpot        ! Inverse
  public :: dZfromTpot      ! Gets thickness between 2 pressure levels

 !/-- Micromet (Aerodynamic) routines

  public :: rh_to_h2ocm3
  public :: rh2vpd2  ! Conversion from fRH (frac) to VPD (Pa).  Orig.
  public :: rh2vpd  ! Conversion from fRH (frac) to VPD (Pa).  Default, uses Teten
  public :: saturated_vapour_pressure  ! for consistency with some DO3SE routines
  public :: rh2q    ! Simple conversion to specific humidity, from Jalle H.
  public :: rh2num  ! from fRH to  molec/cm3 
  public :: num2rh  ! from molec/cm3 to fRH
  public :: q2rh    ! Simple conversion from specific humidity, from Jalle H.
  public :: rh2ePa  ! Simple conversion to vapour pressure
  public :: Teten   ! Simple sat. vapour pressure calculation

  public :: AerRes
  public :: AerResM

  public :: PsiH
  public :: PsiM

  public :: phi_w   !! Added for ESX Leuning work, from J.-P.
  public :: phi_h   !! Added for ESX Leuning work, from J.-P.

  public :: Launiainen1995

  public :: Wind_at_h   !wind for given height

  public :: invMOL  ! inverse Monin-Obukhov length

 !========================================
  contains
 !=======================================================================

 elemental function StandardAtmos_km_2_kPa(h_km) result (p_kPa)
 !=======================================================================
   implicit none

  !+ Converts height (km)  to pressure (kPa) for a US standard Atmosphere
  !  Valid up to 20 km
  !
  ! pw 07/4/2010

   real :: p_kPa
   real   , intent(in)          :: h_km
!   real :: t    ! Temperature (K)

   if( h_km < 11.0 ) then   ! = p_kPa > 22.632
      ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
      !- use the power function replacament, m**n == exp(n*log m)
      p_kPa = 101.325*exp(-5.255876*log(288.15/(288.15-6.5*h_km)))
   else
      p_kPa =  22.632*exp(-0.1576884*(h_km - 11.0)  )

   end if

 end function StandardAtmos_km_2_kPa

 !=======================================================================
 elemental function StandardAtmos_kPa_2_km(p_kPa) result (h_km)
 !=======================================================================
   implicit none

  !+ Converts pressure (kPa)  to height (km) for a US standard Atmosphere
  !  Valid up to 20 km
  !
  ! ds 27/7/2003

   real, intent(in) :: p_kPa
   real             :: h_km
   real :: t    ! Temperature (K)

   if( p_kPa > 22.632 ) then   ! = 11 km height
         ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
         !- use the power function replacament, m**n == exp(n*log m)

         t = 288.15/exp(-1.0/5.255876*log(p_kPa/101.325))
         h_km = (288.15-t)/6.5
   else
         h_km = 11.0 + log( p_kPa/22.632)/(-0.1576884)
   end if

 end function StandardAtmos_kPa_2_km

 !=======================================================================

 !+
 !  Exner functions
 !  Defined here as (p/P0)**KAPPA

 ! Where KAPPA = R/CP = 0.286
 ! P0 = 1.0e5 Pa

 ! CAREFUL:  The term Exner function  can also be used for CP * (p/P0)**KAPPA
 ! Hence notation Exner_nd  - non dimensional version
 !
 ! Exner_nd returns the non-dimensional excner function (p/p0)**R/CP
 !
 ! Added 7/4/2005, Dave, based upon tpi code from tiphys
 ! Test prog at end along with results.
 !----------------------------------------------------------------------------

  elemental function Tpot_2_T(p) result(fTpot)
    ! Identical to Exner_nd
    ! Usage:   T = Tpot * Tpot_2_T(p)
     real, intent(in) :: p    ! Pressure, Pa
     real :: fTpot

      fTpot = (p/1.0e+5)**KAPPA  ! Without CP

  end function Tpot_2_T
  !-------------------------------------------------------------------
  elemental function T_2_Tpot(p) result(fT)
    ! Usage:   Tpot = T * T_2_Tpot(p)
     real, intent(in) :: p    ! Pressure, Pa
     real :: fT

      fT    = (1.0e+5/p)**KAPPA  ! Without CP

  end function T_2_Tpot

  !--------------------------------------------------------------------
  elemental function dZfromTpot(p1,p2,theta) result(dz)
     real, intent(in) :: p1, p2  ! Pressure levels, Pa
     real, intent(in) :: theta   ! potential temp, K
     real :: dz  ! thickness of layer between p1 and p2, m

     dz = theta* CP * ( Tpot_2_T(p1)-Tpot_2_T(p2) )/GRAV
  end function dZfromTpot
!____________________________________________________________________

  !=======================================================================
  !> Function Q_to_RH, returns relative humidity (RH) for Q, T and P.
  !! Code:  Jalle Hiltunen, Q_to_RH 2014

  function q2rh(Q,T,P) result(RH)

  ! Relative humidity (RH, 0-1) is related to saturated specific humidity (Qs)
  ! and temperature as follows (Stull,1995) page 86-87:
  !                   RH = Q/Qs *( 100) ????
  !                   Qs = epsilon * es/P 
  real, intent(in) :: Q,T,P
  real :: Qs, es, RH, r_s             ! es is the saturated vapor pressure
  real, PARAMETER ::  xEPSILON = 0.622  !< Ratio gas constants Rd/Rv [g vapor/g dry air]
!   ,enull = 0.61078E3  ! Saturated vapor pressure at 273.16 K

      es = Teten(T)   !! enull * exp(17.2694*(T-273.16)/(T-35.86))
      r_s = xEPSILON * es/(P-es)
      Qs = r_s/(1+r_s)
      RH = min( 1.0, Q/Qs )
  end function q2rh
  !--------------------------------------------------------------------

  function rh2q(RH,T,P) result(Q)

  ! Specific humidity (Q) is related to saturated specific humidity (Qs) and
  ! relative humidity (RH) (0-1) as follows (Stull,1995) page 86-87:
  !                   Q = Qs*RH (/100)
  !                   Qs = epsilon * es/P 
  ! Code:  Jalle Hiltunen, Q_to_RH 2014
  real, intent(in) :: RH,T,P
  real :: Qs, es, Q, r_s                 ! es is the saturated vapor pressure
  real, PARAMETER ::  xEPSILON = 0.622  !< Ratio gas constants Rd/Rv [g vapor/g dry air]
!     ,enull = 0.61078E3  ! Saturated vapor pressure at 273.16 K

      !es = enull * exp(17.2694*(T-273.16)/(T-35.86))
      es = Teten(T)   !! enull * exp(17.2694*(T-273.16)/(T-35.86))
      r_s = xEPSILON * es/(P-es)
      Qs = r_s/(1+r_s)
      Q = Qs * RH

  end function rh2q

  !--------------------------------------------------------------------

  elemental function rh2ePa(fRH,T) result(ePa)

  real, intent(in) :: fRH,T !< Fractional RH, temperature (K)
  real :: ePa               !< vapour pressure in Pascal
  !real :: es                !< es is the saturated vapor pressure (Pa)
  !real, PARAMETER :: ENULL = 0.61078E3  ! Saturated vapor pressure at 273.16 K

      !es   = ENULL * exp(17.2694*(T-273.16)/(T-35.86))
      ePa  = fRH * Teten(T)

  end function rh2ePa

  !--------------------------------------------------------------------
  !> FUNCTION Teten
  !! Calculate saturation vapour pressure. Approximate, but good from
  !! eg -40 to 40, ie for PBL

  elemental function Teten(T) result(esat)

  real, intent(in) :: T     !< Temperature (K)
  real :: esat              !< Saturation vapour pressure in Pascal
  real, PARAMETER :: ESAT0 = 611.0  ! Sat. vapour pressure (Pa) at 273K
  !real, PARAMETER :: invT0  = 1.0/273.15

    !esat = ESAT0 * exp( EFAC*(invT0 - 1/T)/ RGAS_KG )
    esat = ESAT0 * exp( 17.2694*(T-273.16)/ (T-35.86)  )

  end function Teten

  !--------------------------------------------------------------------
  !> FUNCTION rh2num
  !! Near-surface only !!

  elemental function rh2num(fRH,T) result(h2o)

  real, intent(in) :: fRH,T !< Fractional RH, temperature (K)
  real :: h2o               !< water concentration in  molec/cm3
  real :: ePa               !< vapour pressure in Pascal

!    esat = ESAT0 * exp( EFAC*(invT0 - 1/T)/ RGAS_KG )
!    h2o  = fRH* esat/Pa * M*ATWAIR/18.0
! Now, rho_air= Pa/(RGAS_KG*T) , in kg/m3, therefore:
! Now, nair = P/(RGAS_J*T)   ! moles/m3
!      Xh2o = e/P            ! vol. mixing ratio
!      [h2o] = Xh2o.nair.AVOG.1.0e-6         !  molec/cm3
!            = e/P * P/(RT) * AVOG * 1.0e-6  !  molec/cm3

    ePa = rh2ePa( fRH, T)
    h2o  = ePa/(RGAS_J * T) * AVOG * 1.0e-6

  end function rh2num
  !--------------------------------------------------------------------
  elemental function num2rh(N,T) result(fRH)

  real, intent(in) :: N,T   !< molec/cm3, temperature (K)
  real :: fRH               !< vapour pressure in Pascal

  real :: Nsat              ! molec/cm3 H2O at saturation

    Nsat = rh2num(1.0,T)
    fRH = N/Nsat

  end function num2rh

  !--------------------------------------------------------------------
  !> FUNCTION rh_to_h2ocm3

  elemental function rh_to_h2ocm3(fRH,T) result(h2o)

  real, intent(in) :: fRH,T !< Fractional RH, temperature (K)
  real :: h2o               !< water concentration in  molec/cm3
  real :: ePa               !< vapour pressure in Pascal
  real, PARAMETER :: ENULL = 0.61078E3  ! Saturated vapor pressure at 273.16 K
  real, PARAMETER :: EFAC  = 0.622*2.5e6  !...

!    esat = ESAT0 * exp( EFAC*(invT0 - 1/T)/ RGAS_KG )
!    h2o  = fRH* esat/Pa * M*ATWAIR/18.0

    ePa = rh2ePa( fRH, T)

    !Now, rho = Pa/(RGAS_KG*T) , in kg/m3, therefore:

    h2o  = ePa/(RGAS_KG * T) * 0.001/18.0 * AVOG * 1.0e-6

  end function rh_to_h2ocm3

  !--------------------------------------------------------------------

  elemental function saturated_vapour_pressure(T) result (esat)
    real, intent(in) ::  T    ! Temperature (K)
    real :: esat
    eSat = rh2ePa(1.0, T)
  end function saturated_vapour_pressure
  !--------------------------------------------------------------------
  elemental function rh2vpd2(T,fRH) result (vpd_res)
  !This function is not currently in use.

    real, intent(in) ::  T    ! Temperature (K)
    real, intent(in) ::  fRH   ! relative humidity (%)
    real :: vpd_res     ! vpd   = water vapour pressure deficit (Pa)

    !   Local:
    real :: vpSat       ! vpSat = saturated water vapour pressure (Pa)
    real :: arg

    arg   = 17.67 * (T-273.15)/(T-29.65)
    vpSat = 611.2 * exp(arg)
    vpd_res   = (1.0 - fRH) * vpSat

  end function rh2vpd2

  !--------------------------------------------------------------------
  ! rh2vpd:  Conversion from fRH (frac) to VPD (Pa).  Default, uses Teten

  elemental function rh2vpd(T,fRH) result (vpd_Pa)

    real, intent(in) ::  T    ! Temperature (K)
    real, intent(in) ::  fRH   ! relative humidity (fraction)
    real :: vpd_Pa      ! vpd   = water vapour pressure deficit (Pa)

    !   Local:
    real :: vpSat       ! vpSat = saturated water vapour pressure (Pa)

    vpSat = Teten(T)
    !arg   = 17.67 * (T-273.15)/(T-29.65)
    !vpSat = 611.2 * exp(arg)
    vpd_Pa    = (1 - fRH) * vpSat

  end function rh2vpd

  !--------------------------------------------------------------------

  function AerRes(z1,z2,uStar,Linv,Karman) result (Ra)
!...
!   Ref: Garratt, 1994, pp.55-58
!   In:
    real, intent(in) ::   z1     ! lower height (m), equivalent to h-d+1 or h-d+3
    real, intent(in) ::   z2     ! upper height (m), equivalent to z-d
    real, intent(in) ::   uStar  ! friction velocity (m/s)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)
    
    real, intent(in) ::   Karman ! von Karman's constant 
!   For AerRes, the above dummy argument is replaced by the actual argument 
!   KARMAN in the module GetMet_mod.

!   Out:
    real :: Ra      ! =  aerodynamic resistance to transfer of sensible heat
                    !from z2 to z1 (s/m)

!   uses functions:
!   PsiH   = integral flux-gradient stability function for heat 
!...

    if ( z1 > z2 ) then
      Ra = -999.0
    else
      Ra = log(z2/z1) - PsiH(z2*Linv) + PsiH(z1*Linv)
      Ra = Ra/(Karman*uStar)
    end if

  end function AerRes
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  function AerResM(z1,z2,uStar,Linv,Karman) result (Ra)
!...
!   Ref: Garratt, 1994, pp.55-58
!   In:
    real, intent(in) ::   z1     ! lower height (m), equivalent to h-d+1 or h-d+3
    real, intent(in) ::   z2     ! upper height (m), equivalent to z-d
    real, intent(in) ::   uStar  ! friction velocity (m/s)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)
    
    real, intent(in) ::   Karman ! von Karman's constant 
!   For AerRes, the above dummy argument is replaced by the actual argument 
!   KARMAN in the module GetMet_mod.

!   Out:
    real :: Ra     ! =  aerodynamic resistance to transfer of momentum
                    !from z2 to z1 (s/m)

!   uses functions:
!   PsiM   = integral flux-gradient stability function for momentum
!...

    Ra = log(z2/z1) - PsiM(z2*Linv) + PsiM(z1*Linv)
    Ra = Ra/(Karman*uStar)

  end function AerResM

  !--------------------------------------------------------------------
  elemental function PsiH(zL) result (stab_h)
    !  PsiH = integral flux-gradient stability function for heat 
    !  Ref: Garratt, 1994, pp52-54
    !  VDHH modified - use van der Hurk + Holtslag?

    ! In:
    real, intent(in) :: zL   ! surface layer stability parameter, (z-d)/L 
    
    ! Out:
    real :: stab_h         !   PsiH(zL) 
    
   ! Local
   real :: x
   real, parameter :: a=1, b=0.667, c=5.0, d=0.35
 
    if (zL <  0) then !unstable
        x    = sqrt(1.0 - 16.0 * zL)
        stab_h = 2.0 * log( (1.0 + x)/2.0 )
    else             !stable
        !ESX if ( FluxPROFILE == "Ln95" ) then
        !ESX    stab_h = -( (1+2*a/3.0*zL)**1.5 + b*(zL-c/d)* exp(-d*zL) + (b*c/d-1) )
        !ESX else 
           stab_h = -5.0 * zL
        !ESX end if
    end if

  end function PsiH

  !--------------------------------------------------------------------
  elemental function PsiM(zL) result (stab_m)
   !   Out:
   !   PsiM = integral flux-gradient stability function for momentum 
   !   Ref: Garratt, 1994, pp52-54

    real, intent(in) ::  zL    ! = surface layer stability parameter, (z-d)/L 
                               ! notation must be preserved         
    real :: stab_m
    real  :: x
   real, parameter :: a=1, b=0.667, c=5.0, d=0.35
 
    if( zL < 0) then !unstable
       x    = sqrt(sqrt(1.0 - 16.0*zL))
       stab_m = log( 0.125*(1.0+x)*(1.0+x)*(1.0+x*x) ) +  PI/2.0 - 2.0*atan(x)
    else             !stable
        !ESX if ( FluxPROFILE == "Ln95" ) then
        !ESX ? Better up to zL~10
        !ESX    stab_m = -( a*zL + b*(zl-c/d)*exp(-d*zL) + b*c/d)
        !ESX else
           stab_m = -5.0 * zL
        !ESX end if
    end if

  end function PsiM
  !--------------------------------------------------------------------

  elemental function phi_h(zL) result (phiH)
    !  PhiH = flux-gradient stability function for heat 
    real, intent(in) :: zL   ! surface layer stability parameter, (z-d)/L 
    real ::  phiH         !
 
    if (zL <  0) then !unstable
         phiH    = 1.0/sqrt(1.0 - 16.0 * zL)
    else             !stable
         phiH = 1.0 + 5 * zL
    end if

  end function phi_h

!--------------------------------------------------------------------
  elemental function phi_w(zL) result ( phiW)
    !  PhiW = flux-gradient stability function for W (water?? Check!)
    real, intent(in) :: zL   ! surface layer stability parameter, (z-d)/L 
    real ::  phiW         !
 
    if (zL <  0) then !unstable
         phiW    = 1.25*(1-3*zL)**0.3333
    else             !stable
         phiW    = 1.25*(1 + 0.2*zL)
    end if

  end function phi_w

!--------------------------------------------------------------------
  subroutine Launiainen1995 (u, z, z0m, z0mh, theta0, theta, invL)
  real, intent(in) :: u  ! winds
  real, intent(in) :: z ! mid-cell height
  real, intent(in) :: z0m ! roughness ht., momentum
  real, intent(in) :: z0mh ! ration roughness ht., momentum
  real, intent(in) :: theta0  !pot. temp at surface
  real, intent(in) :: theta   !pot. temp at ref ht.
  real, intent(out) :: invL
  real :: zeta  ! z/L
  real :: z0h, logzz0m, Rib
  z0h = z0m/ z0mh

   ! Ignoring virtual temp (and Lau has no z0):

       !Rib =      GRAV * z * (theta-theta0 ) / &
       Rib =      9.81 * z * (theta-theta0 ) / &
                    ( theta0 *  u**2 + 0.001 ) !!! EPS )

       logzz0m = log(z/z0m)
      if ( Rib <0.0 ) then

          Rib = max ( -3.0, Rib)  ! Limit used by van der Hurk + Holtslag, 1996
          zeta = ( logzz0m**2/log(z/z0h) - 0.55 ) * Rib


       else

         Rib = min ( 1.0, Rib)  ! Limit used by van der Hurk + Holtslag, 1996
         zeta =(  1.89  * logzz0m + 44.2 )*Rib**2 + ( 1.18*logzz0m -1.37) * Rib
         if( Rib > 0.08 ) zeta = zeta - 1.5*log(z0m/z0h)*Rib
       end if
       invL = zeta/z

  end subroutine Launiainen1995

!--------------------------------------------------------------------
  elemental function Wind_at_h(u_ref, z_ref, zh, d, z0, Linv) result (u_zh)
!...
!   Ref: Garratt, 1994, 
!   In:
    real, intent(in) ::   u_ref  ! windspeed at z_ref
    real, intent(in) ::   z_ref  ! centre of call, ca. 45m (m)
    real, intent(in) ::   zh     ! height required (m)
    real, intent(in) ::   d      ! displacement height (m)
    real, intent(in) ::   z0     ! roughness height (m)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)

!   Out:
    real :: u_zh    ! =   wind-speed at height h (m/s)

     u_zh = u_ref *  &
          ( log((zh-d)/z0)  -PsiM((zh-d)*Linv) + PsiM(z0*Linv)  )/ &
          ( log((z_ref-d)/z0) -PsiM((z_ref-d)*Linv) + PsiM(z0*Linv))

    !NB - COULD USE INSTEAD: Ra = log(z2/z1) - PsiM(z2*Linv) + PsiM(z1*Linv)
    ! Or could optimise with explicit PsiM, etc.

  end function Wind_at_h

!----------------------------------------------------------------------------

  elemental function invMOL(ustar_in, Hd, rho_s, t2) result(invL)
    ! taken from EMEP (as in Monteith & Unsworth, 1990, p. 237)
    real, intent(in) :: ustar_in ! friction velocity (m/s)
    real, intent(in) :: Hd       ! sensible heat flux (W/m^2)
    real, intent(in) :: rho_s    ! surface density of air (kg/m^3)
    real, intent(in) :: t2       ! 2 m temperature (K)
    
    ! out:
    real :: invL  ! inverse Monin-Obukhov length (m)
    real :: ustar

    real, parameter :: Cp = 1005 ! from esx_RadEB (J/(kg K))
    ! we limit u* to a physically plausible value
    ! to prevent numerical problems
     ustar = max( ustar_in, 0.1 )

     invL  = -1* KARMAN * GRAV * Hd & ! -Grid%Hd disliked by gfortran
            / (Cp * rho_s * ustar**3 * t2 )

    !.. we limit the range of 1/L to prevent numerical and printout problems
    !.. and because we don't trust HIRLAM or other NWPs enough.
    !   This range is very wide anyway.

    ! Grid%invL  = max( -1.0, Grid%invL ) !! limit very unstable
    ! Grid%invL  = min(  1.0, Grid%invL ) !! limit very stable
  end function invMOL

!!TSTEMX T = 273.15; Pa=50.0e3
!!TSTEMX print "(a,3f8.2,4es12.3)", "ePa ", T,fRH,kg2g*rh2q(fRH,T,Pa), Pa, ePa, xH2O, Teten(T)
!!TSTEMX T = 298.15; fRH=0.99
!!TSTEMX print *, "VPD comps", fRH, rh2vpd(T,fRH), rh2vpd(T, 1.0), saturated_vapour_pressure(T)

end module MetFunctions
!TSTEMX program testr
!TSTEMX use MetFunctions, only : Tpot_2_T, T_2_Tpot,StandardAtmos_kPa_2_km, &
!TSTEMX  StandardAtmos_km_2_kPa, rh2ePa, rh_to_h2ocm3, Teten, rh2q, rh2vpd, &
!TSTEMX  saturated_vapour_pressure, invMOL, dZfromTpot
!TSTEMX implicit none
!TSTEMX real :: p, hkm, T, fRH, kg2g=1.0e3, ePa, xH2O, Pa
!TSTEMX integer :: i, iRH
!TSTEMX fRH = 1.0; T = 273.15 + 30; Pa = 1.01325e5 ! test values
!TSTEMX ePa = rh2ePa( fRH, T )
!TSTEMX xH2O = rh_to_h2ocm3(fRH,T)
!TSTEMX print "(a,3f8.2,4es12.3)", "ePa ", T,fRH,kg2g*rh2q(fRH,T,Pa), Pa, ePa, xH2O, Teten(T)
!TSTEMX print *, 'testing invL = ', invMOL(1.5, 30., 1.1, 290.)
!TSTEMX print "(a8,4a12)", "p (hPa)", "Tpot2T-fac", "T2Tpot-fac", "h(km)"
!TSTEMX do i = 1, 20
!TSTEMX     p = 0.05 * i * 1.0e5
!TSTEMX     hkm = StandardAtmos_kPa_2_km( 0.001*p)
!TSTEMX     print "(f8.3,4f12.5)", 1.0e-2*p, Tpot_2_T(p), T_2_Tpot(p), hkm
!TSTEMX  end do
!TSTEMX do iRH = 80, 100, 5
!TSTEMX   fRH=0.01 * iRH ! 0.29
!TSTEMX   print *, "VPD comps", fRH, rh2vpd(T,fRH), rh2vpd(T, 0.9999999), saturated_vapour_pressure(T)
!TSTEMX end do
!TSTEMX print *, "TEST p to dz ", dZfromTpot(1.0e5,0.850e5,285.0)
!TSTEMX print *, StandardAtmos_km_2_kPa( 0.0)
!TSTEMX end program testr
