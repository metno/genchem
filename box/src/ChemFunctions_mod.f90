! <ChemFunctions_mod.f90 - the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
module ChemFunctions_mod
!____________________________________________________________________
! Miscellaneous collection of "standard" functions for chemical
! calculations.  Includes Troe, sine and cosine curves,  and some
! from KPP system.
!
! Where possible, reference to the EMEP documentation paper, Simpson
! et al., ACP, 2012,  are given, indicated by ACP:
! 
!____________________________________________________________________
!
!** includes
!   troe - standard chemical function
!____________________________________________________________________
!ESX use LocalVariables,     only : Grid   ! => izen, is_NWPsea
 use AeroConstants_mod,     only : AERO ! !for n2o5Hydrolysis AJ2018
 use AeroFunctions_mod,     only: UptakeRate, GammaN2O5_EJSS, GammaN2O5
 use CheckStop_mod,         only : StopAll, CheckStop
 use ChemSpecs_mod !,          only : SO4, NO3_f, NH4_f, NO3_c
 use Config_module,         only : USES
 use PhysicalConstants_mod, only : AVOG, RGAS_J !ESX, DAY_ZEN
 use NumberConstants, only : UNDEF_R !ESX, DAY_ZEN
 use ZchemData_mod,          only : x=> xChem
 use ZchemData_mod,          only : xSO4, xNO3, xNH4 & ! for RiemerN2O5
   ,gamN2O5, aero_fom,aero_fss,aero_fdust, aero_fbc  & !for n2o5Hydrolysis AJ2018
                               ,S_m2m3, cN2O5 & !for n2o5Hydrolysis AJ2018
                               ,temp, tinv, rh, M   !, amk
  implicit none
  private

!ESX  public :: troe
!ESX  public :: troeInLog  ! When log(Fc) provided
  public :: IUPAC_troe ! Using the approximate expression for F from Atkinson et al., 2006 (ACP6, 3625)
!ESX .... kept just to keep things working
  public ::  kaero !aerosol production rate
!ESX  public ::  kaero2    ! for testing

  public ::  HydrolysisN2O5
  public ::  RiemerN2O5
  public ::  ec_ageing_rate
  public ::  kmt3      ! For 3-body reactions, from Robert OCt 2009
! HI next things for KPP
  public :: Update_SUN
  public :: ARR
  public :: EP2
  public :: EP3
  public :: FALL


! weighting factor for N2O5 hydrolysis
! Mass of sulfate relative to sulfate+nitrate
! according to  Riemer N, Vogel H, Vogel B, 
! Schell B, Ackermann I, Kessler C, Hass H
! JGR 108 (D4): FEB 27 2003 

  real, parameter, public :: VOLFACSO4 = 96.0/(AVOG) * 1.2648  *0.02/0.068e-6 
  real, parameter, public :: VOLFACNO3 = 62.0/(AVOG) * 1.2648  *0.02/0.068e-6 
  real, parameter, public :: VOLFACNH4 = 18.0/(AVOG) * 1.2648  *0.02/0.068e-6 

  real, public, dimension(1:4), save :: YN2O5HYD=UNDEF_R

!HI instead of reading this from the Initialize module (and being more
!    flexible in terms of units)
!DS QUERY .... could be variable?
  real, parameter, public :: CFACTOR = 2.55e+10 ! -> everything in ppb

  real, public, save :: SUN            ! for KPP subroutines/testing
  logical, public, save :: is_day=.true.  !DS Sep2016 - used for EC_AGEING...


  !========================================
  contains
  !========================================

  ! -------------------------------------------------------------------------
  ! HI now the functions needed for KPP; START
  ! -------------------------------------------------------------------------
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_SUN - update SUN light using TIME
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Update_SUN(TIME)
      !USE my_saprc99_Parameters
      !USE my_saprc99_Global

    IMPLICIT NONE

    REAL :: TIME, SunRise, SunSet
    REAL :: Thour, Tlocal, Ttmp 
    ! PI - Value of pi
    REAL, PARAMETER :: PI = 3.14159265358979
    
    SunRise = 4.5 
    SunSet  = 19.5 
    Thour = TIME/3600.0 
    Tlocal = Thour - (INT(Thour)/24)*24

    IF ((Tlocal>=SunRise).AND.(Tlocal<=SunSet)) THEN
       Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise)
       IF (Ttmp.GT.0) THEN
          Ttmp =  Ttmp*Ttmp
       ELSE
          Ttmp = -Ttmp*Ttmp
       END IF
       SUN = ( 1.0 + COS(PI*Ttmp) )/2.0 
    ELSE
       SUN = 0.0 
    END IF

 END SUBROUTINE Update_SUN

! End of Update_SUN function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Arrhenius
! HI made functions elemental, removed DBLE() and added temperature as input
   ELEMENTAL REAL FUNCTION ARR( A0,B0,C0,TEMP)
      REAL, intent(in) :: A0,B0,C0,TEMP
      ARR =  A0 * EXP(-B0/TEMP) * (TEMP/300.0)**C0
   END FUNCTION ARR        

   ELEMENTAL REAL FUNCTION EP2(A0,C0,A2,C2,A3,C3,TEMP)
      REAL, intent(in) ::  A0,C0,A2,C2,A3,C3,TEMP
      REAL K0,K2,K3            
      K0 = A0 * EXP(-C0/TEMP)
      K2 = A2 * EXP(-C2/TEMP)
      K3 = A3 * EXP(-C3/TEMP)
      K3 = K3*CFACTOR*1.0E6
      EP2 = K0 + K3/(1.0+K3/K2 )
   END FUNCTION EP2

   ELEMENTAL REAL FUNCTION EP3(A1,C1,A2,C2,TEMP) 
      REAL, intent(in) ::  A1, C1, A2, C2,TEMP
      REAL K1, K2      
      K1 = A1 * EXP(-C1/TEMP)
      K2 = A2 * EXP(-C2/TEMP)
      EP3 = K1 + K2*(1.0E6*CFACTOR)
   END FUNCTION EP3 

   ELEMENTAL REAL FUNCTION FALL ( A0,B0,C0,A1,B1,C1,CF,TEMP)
      REAL, intent(in) ::  A0,B0,C0,A1,B1,C1,CF,TEMP
      REAL K0, K1     
      K0 = A0 * EXP(-B0/TEMP)* (TEMP/300.0)**C0
      K1 = A1 * EXP(-B1/TEMP)* (TEMP/300.0)**C1
      K0 = K0*CFACTOR*1.0E6
      K1 = K0/K1
      FALL = (K0/(1.0+K1))*   &
           CF**(1.0/(1.0+(LOG10(K1))**2))
   END FUNCTION FALL

  ! -------------------------------------------------------------------------
  ! HI now the functions needed for KPP; END
  ! -------------------------------------------------------------------------


  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! KMT3 uses air concenrtation (M) and inverse Temp (tinv) from Zmet
  !
  function kmt3(a1,c1,a3,c3,a4,c4) result (rckmt3)
     real, intent(in)  :: a1,c1,a3,c3,a4,c4
     real, dimension(size(M)) :: rckmt3
     real, dimension(size(M)) ::  k1, k3, k4
       k1 = a1 * EXP(C1*tinv)
       k3 = a3 * EXP(C3*tinv)
       k4 = a4 * EXP(C4*tinv)
       rckmt3 = k1 + (k3*M)/(1.0+(k3*M)/k4)
  end function kmt3

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental function troe(k0,kinf,Fc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! ds note - this isn't checked or optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,Fc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf    ! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!    could have Fc already as log(Fc) to save CPU, but for now
!    keep as proper Fc. Slower but less confusing

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*log(Fc))

  end function troe
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  elemental function troeInLog(k0,kinf,LogFc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! note - this isn't optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,LogFc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf    ! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!    give Fc already as log(Fc)

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*LogFc)

  end function troeInLog

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental function IUPAC_troe(k0,kinf,Fc,M,N) result (rctroe)

  !+ Calculates Troe expression 
  ! -----------------------------------------------------------
  ! note - this isn't optimised yet. Taken from
  ! Atkinson et al. ACP 2006, 6, 3625-4055. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.
  ! NOTE that in the IUPAC nomenclature k0 already contains [M] so 
  !  the k0(IUPAC)=k0*M here
  !   N=[0.75-1.27*log10(Fc)]

     real, intent(in)  :: k0,kinf,Fc,M,N
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacement, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf    ! used also below
     x    = log10(y)/N
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

     rctroe = k0M / ( 1.0 + y) * exp(x*log(Fc))

  end function IUPAC_troe
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


! N2O5 -> nitrate calculation
! EmChem16s use:
  function HydrolysisN2O5(ormethod) result(rate)
   character(len=*), intent(in) , optional:: ormethod ! overrides default method if wanted
   character(len=30), save :: method
   real, dimension(size(temp)) :: rate
   real    :: rc
   real    :: f   ! was f_riemer
   real    :: gam, gamss, gamDu, s,  s_ss, s_du, Rwet  ! for newer methods
   real, save :: g1 = 0.02, g2=0.002 ! gammas for 100% so4, 100% no3, default
  ! fixed-value gammas can be specified with e.g. gamma:0.02. we derive
  ! the numerical value, gfix, from this string
   real, save :: gfix= -999.         ! fixed-value, from gamma:xxxx values
   character(len=20) :: gtxt         ! for gamma:xxxx values
   real, parameter :: epsil = 1.0  ! one mol/cm3 to stop div by zero
   integer :: k
   integer, save :: K1, K2
   real :: xno3, xso4, fine_rate, coarse_rate, yld
   logical, save :: first_call = .true.
   character(len=*), parameter :: dtxt = 'hydroln2o5:'

   call checkstop( so4<1, "so4 not defined" )
   call checkstop( nh4_f<1, "nh4_f not defined" )
   call checkstop( no3_f<1, "no3_f not defined" )
   call checkstop( no3_c<1, "no3_c not defined" )
 
   if( first_call ) then
     method = USES%n2o5HydrolysisMethod
     if ( present(ormethod) ) method = ormethod  ! WHEN is this used?
     if( method(1:6)=="Gamma:"  ) then
       gtxt=method(7:)
       read(gtxt,*) gFix
       method='gFixed'
      end if
      K1 = 1 ! NOT: lbound(temp)
      K2 = size(rh)  !ubound(temp)
   end if
  select case ( method )
!OLD    case ( "ORIGRIEMER","OrigRiemer")
!OLD
!OLD      do k = lbound(temp,1), ubound(temp,1) ! K1, K2
!OLD       if ( rh(k)  > 0.4) then
!OLD          xNO3 = x(NO3_f,k) + x(NO3_c,k)
!OLD
!OLD         !mean molec speed of N2O5 (MW 108), m/s
!OLD         ! with density corrected for rh (moderate approx.)
!OLD          !J18 rc = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mol.speed (m/s)
!OLD          rc = sqrt(3.0 * RGAS_J * temp(k) / 0.108) & ! mol.speed (m/s)
!OLD             /(4*(2.5 - rh(k)*1.25))                   ! density
!OLD
!OLD          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )
!OLD
!OLD
!OLD          rate(k) =  (0.9*f + 0.1) * rc *  &
!OLD                !TEST   0.5 * & ! v. loosely based on Reimer 2009 
!OLD             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * xNO3  &
!OLD              + VOLFACNH4 * x(NH4_f,k) )    !SIA aerosol surface
!OLD        else
!OLD          rate(k) = 0.0
!OLD        end if
!OLD      end do ! k
  !---------------------------------------
   case ( "Smix", "SmixTen" )

!if ( DEBUG%RUNCHEM .and. DebugCell ) then
!  write(*,*) dtxt//trim(method), K1, K2, rh(K2), S_m2m3(AERO%PM_F,K2) , S_m2m3(AERO%DU_C,K2), cN2O5(1), aero_fom(1)
!stop 'NNN'
!end if
     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???
            xSO4 = x(SO4,k)
            xNO3 = x(NO3_f,k) !just fine PM
            f = 96*xSO4/( 96*xSO4 + 62* xNO3  + EPSIL )

            S = S_m2m3(AERO%PM_F,k)
            gam = GammaN2O5(temp(k),rh(k),&
                   f,aero_fom(k),aero_fss(k),aero_fdust(k),aero_fbc(k))


            rate(k) = UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM
            fine_rate = rate(k)

            !Add coarse model ! was SmixC
            S_ss = S_m2m3(AERO%SS_C,k)
            gamSS=GammaN2O5_EJSS(rh(k))
            S_du = S_m2m3(AERO%DU_C,k)
            ! gamDU=0.01 ! for dust
            ! same as UptakeRate(cN2O5,gam,S), but easier to code here:
            coarse_rate = cN2O5(k)*(gamSS*S_ss+0.01*S_du)/4
            rate(k) = rate(k) + coarse_rate
            
            ! ToDo update gam for export. Currently at fine-mod only
            !Coarse end
            if( method == "SmixTen") then
              gam = 0.1 * gam ! cf Brown et al, 2009!
              rate(k) = 0.1 * rate(k)
              fine_rate = 0.1 * fine_rate
              coarse_rate = 0.1 * coarse_rate
            end if
       else
            gam = 0.0 ! just for export
            rate(k) = 0.0
            fine_rate = 0.0
            coarse_rate = 0.0
       end if
       gamN2O5(k) = gam ! just for export
       
       ! Calculate HNO3 and ClNO2 Yields
       call HydrolysisN2O5_Yields( k, fine_rate, coarse_rate, yld )

    end do

    case ( "gFixed")  !  Fixed gammas
     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

          gam = gFix ! Found above

          S = S_m2m3(AERO%PM_F,k) !fine SIA +OM + ...
          rate(k) = UptakeRate(cN2O5(k),gam,S)
      else
         rate(k) = 0.0
         gam     = 0.0 ! just for export
      end if
      gamN2O5(k) = gam ! just for export

    end do ! k

     case default
       call StopAll("Unknown N2O5 hydrolysis"//method )

   end select
   first_call = .false.
  end function HydrolysisN2O5
  !---------------------------------------------------------------------

!===========================================================================
! N2O5 -> nitrate calculation. Some constants for
! calculation of volume fraction of sulphate aerosol, and rate of uptake
! 
!
! The first order reaction coefficient K (corrected for gas phase diffusion, 
! Schwartz, 1986) is given by
!
! K= S* alpha* v/4                               ACP:44
!    alpha=sticking coeff. for N2O5 =0.02
!    v=mean molecular speed for N2O5
!    S=aerosol surfac
!
! The surface area of the aerosols can be calculated as
! 
! S = V * surface/volume of aerosols
!     V=volume fraction of sulphate (cm3 aerosol/cm3 air)
!     (similar for nitrate and ammonium):
!
!     e.g. simplest form (not used) would be:
!     V = (so4 in moleculescm-3) x atw sulphate
!         ---------------------------------------------------------
!        AVOG X specific density of aerosols (assumed 2g/cm3*rh correction)
!
!    Or, shorter, V = C x M0/(AVOG*rho)
!
!    where C is conc. e.g. sulphate (molecule/cm3), M0 is molwt. 
!    We do not want to include  concentrations  or rho yet, so:
!
!     Let VOL =  M0/AVOG
!   
! E12:47
! The surface/volume ratio is calculated using Whitby particle distribution
! with number mean radius 0.034  and standars deviation (Sigma)=2. 
! Then surface/volume=3/r *  exp( -5/2 *(lnSigma)^2)=26.54 
! 3* exp( -5/2 *(lnSigma)^2)=0.90236
! (monodisperse aerosols; 4*pi*r^2/(4/3 pi*r^3)= 3/r =88.2)
!
! Then 
!      A = VOL * S * 0.90236 /(0.034e-6*rho) 
! and
!      K = VOL * S * 0.90236 /(0.034e-6*rho)    * alpha* v/4
! Set
!      VOLFAC= VOL*0.90236/0.034e-6 *alpha    
! Then
!      K = VOLFAC *C *v/(4*rho)
!
! rcmisc k=v/(4*rho) 
!
!      K = VOLFAC *rcmisc() *C
!
! According to Riemer et al, 2003, we weight the reaction probability
! according to the composition of the aerosol
!
! alpha(N2O5)=f*alpha1 +(1-f)alpha2                           ACP:45
!   alpha1=0.02
!   alpha2=0.002
!   f= Mso4/(Mso4+Mno3), M=aerosol mass concentration         ACP:46
 
! N2O5 -> aerosol based upon  based on Riemer 2003 and
! In testing, we had also tried a simple acounting for 
! results shown in Riemer et al., 2009.
! We did not attempt to model OC, but simply reduce the rate by
! a factor of two to loosely account for this effect. 
! June08 - changed from use of more accurate xnew to xn_2d, since
! surface area won't change so much, and anyway the uncertainties
! are large. (and xn_2d leads to fewer dependencies)

  elemental function RiemerN2O5(rh,tK,xSO4,xNO3,xNH4) result(rate) 
     !ESX real, dimension(K1:K2) :: rate
     real, intent(in) :: rh             !< rel. hum., 0-1
     real, intent(in) :: tK             !< Temp, K
     real, intent(in) :: xSO4,xNO3,xNH4 !< concentrations, vol. units, eg cm-3
     real :: rate
     real :: rc, f
     real, parameter :: EPSIL = 1.0  ! One mol/cm3 to stop div by zero

     !NB: xNO3=0.0
     ! As the partitioning between fine and coarse is so difficult
     ! we include both in the nitrate used here.

     if ( rh > 0.4 ) then 

          rc = sqrt(3 * RGAS_J * tK / 0.108) & ! mean mol. speed,m/s
             /(4*(2.5 - rh*1.25)) !density, corrected for rh (moderate approx.)

          f = 96.0*xSO4/( 96*xSO4 + 62* xNO3  + EPSIL )

           !RESTORE TO ORIG  0.5 * very loosely based on OC effects from Reimer 2009 
           !SIA aerosol surface

          rate =  (0.9*f + 0.1) * rc *  &
             ( VOLFACSO4 * xSO4 + VOLFACNO3 * xNO3 + VOLFACNH4 * xNH4 )
     else 
          rate = 0.0
     end if

  end function RiemerN2O5
  !---------------------------------------------------------------------

  elemental function kaero(rh) result(rate) 
     real, intent(in) :: rh  ! fractional RH
    ! Former rate for HNO3 -> NO3_c, not now used
     !DS real, dimension(size(rh)) :: rate
     real :: rate
     
      if ( rh  > 0.9)  then
         rate = 1.0e-4
      else
         rate = 5.0e-6
      end if

  end function kaero
  !---------------------------------------------------------------------
!ESX  function kaero2() result(rate) 
!ESX    ! New rate for HNO3 -> NO3_c, used only over sea squares
!ESX    ! as very crude simulation of sea-salt HNO3 interactions
!ESX    ! near surface (layer 16 ca. 600m).
!ESX     real, dimension(K1:K2) :: rate
!ESX     integer :: k
!ESX     
!ESX    if ( Grid%is_NWPsea) then
!ESX      rate(K1:15) = 0.0
!ESX      do k = 16, K2
!ESX        if ( rh(k)  > 0.9) then
!ESX           rate(k) = 1.0e-4
!ESX        else
!ESX           rate(k) = 5.0e-6
!ESX        end if
!ESX      end do !k
!ESX    else ! over land
!ESX      rate(K1:K2) = 0.0
!ESX    end if
!ESX  end function kaero2
 !---------------------------------------------------------------------
  function ec_ageing_rate() result(rate) 
 
   !.. Sets ageing rates for fresh EC [1/s] loosely based on Riemer etal. ACP(2004)
   !   See also Tsyro et al, JGR, 112, D23S19, 2007
   !   ---------------------------------  

    real, dimension(size(temp)) :: rate
    !D16 integer, save :: kmax
    logical, save :: first_call=.true.
    if ( first_call ) then
       !kmax = ubound(temp) !DS QUERY FAILS WHY?
     !D16  kmax=1   ! BOX MODEL FIX
       write(*,*) "ec age test BOUNDS ", lbound(temp), ubound(temp)
       first_call = .false.
    end if
    
 
    !ESXDS if ( Grid%izen <= DAY_ZEN ) then  ! daytime
    if ( is_day ) then  ! daytime

       !rate (kmax-2 : kmax)  = 3.5e-5  !  half-lifetime ~ 8h
       !rate (1 : kmax-3)     = 1.4e-4  !                ~ 2h
       rate (: )     = 1.4e-4  !                ~ 2h
      else
     !D16  rate (1 : kmax )      = 9.2e-6  !                ~ 30h
       rate (: )      = 9.2e-6  !                ~ 30h
    endif

  end function ec_ageing_rate

!-------------------------------------------------------------------  
  subroutine HydrolysisN2O5_Yields(k, fine_rate, coarse_rate, yld)
!
! Calculate the yield of HNO3, ClNO2, and particulate NO3 from 
! the hydrolysis of N2O5. The well-known N2O5 hydrolysis on wet 
! aerosols may be written as:
!
!        N2O5 --(PM_H2O)-->  2 HNO3
!
! EMEP predicts the rate of this reaction with an uptake coefficient
! (gamma) that may be calculated with several optional approaches.
! However, particulate chloride can react with the N2O5 and yield a
! different product distribution:
!
!        N2O5 + Cl-  --->  ClNO2 + NO3-
!
! Particulate nitrate ion is formed in this reaction but the ClNO2
! product efficiently photolyzes and reduces the net nitrate formed 
! from N2O5 hydrolysis.
!
! The branching ratio for these to products may be calculated as:
!
!        F_ClNO2 = [Cl-] / ([Cl-] + [H2O]/450)
!
! Thornton et al. (2010), Nature, "A large atomic chlorine source 
! inferred from mid-continental reactive nitrogen chemistry."
!
! The array YN2O5HYD is of length 4 and specifies the yield of 
! each product from the N2O5 hydrolysis:
!    YN2O5HYD(1) =: Y_HNO3 (include the factor of 2 here)
!    YN2O5HYD(2) =: Y_ClNO2 
!    YN2O5HYD(3) =: Y_NO3_f (Fine-mode nitrate)
!    YN2O5HYD(4) =: Y_NO3_c (coarse-mode nitrate) 
!-------------------------------------------------------------------
  
  use ZchemData_mod, only : x => xChem, xh2o_f, xh2o_c

  implicit none

  real, intent(in)    :: fine_rate, coarse_rate ! Hydrolysis Rates [s-1]
  integer, intent(in) :: k ! model layer
  real                :: yld
  real, parameter :: Ncm3_to_molesm3 = 1.0e6 / AVOG    ! #/cm3 to moles/m3
  real, parameter :: MWCL   = 35.453
  real :: xCl, xH2O, F_ClNO2_f, F_ClNO2_c, f_fine

  call CheckStop( SeaSalt_f<1, "CL_f not defined" )
  call CheckStop( SeaSalt_c<1, "CL_c not defined" )

  ! Fine-Mode ClNO2 Yield
  xCl = x( SeaSalt_f,k ) * Ncm3_to_molesm3 * species(SeaSalt_f)%molwt * 0.55 / MWCL
  xH2O = xh2o_f(k) / 1.0e6 / 18.0153 ! ug/m3 -> mol/m3
  if ( xCl + xH2O > 1e-15 ) then
    F_ClNO2_f = xCl / ( xCl + xH2O/450. )
  else
    F_ClNO2_f = 0.0
  end if

  ! Coarse-Mode ClNO2 Yield
  xCl = x( SeaSalt_c,k ) * Ncm3_to_molesm3 * species(SeaSalt_c)%molwt * 0.55 / MWCL
  xH2O = xh2o_c(k) / 1.0e6 / 18.0153 ! ug/m3 -> mol/m3
  if ( xCl + xH2O > 1e-15 ) then
    F_ClNO2_c = xCl / ( xCl + xH2O/450. )
  else
    F_ClNO2_c = 0.0
  end if
  
  if ( fine_rate + coarse_rate > 0.0 ) then
    ! Fraction of Rate in Fine Mode
    f_fine = fine_rate / (fine_rate + coarse_rate)

    ! Calculate total yield of ClNO2
    YN2O5HYD(2) = f_fine * F_ClNO2_f + (1.-f_fine) * F_ClNO2_c

    ! Calculate total yield of Fine and Coarse Nitrate
    YN2O5HYD(3) = f_fine * F_ClNO2_f 
    YN2O5HYD(4) = (1.-f_fine) * F_ClNO2_c

    ! Calculate total yield of HNO3
    YN2O5HYD(1) = 2.0 * ( f_fine * (1.-F_ClNO2_f) + (1.-f_fine) * (1.-F_ClNO2_c) )
  else
    YN2O5HYD(:) = 0.0
  end if

  yld = YN2O5HYD(1)

  end subroutine HydrolysisN2O5_Yields    
 
end module ChemFunctions_mod
!TSTEMX program tester
!TSTEMX use ChemFunctions_mod
!TSTEMX print *, kaero(0.5) 
!TSTEMX end program tester
