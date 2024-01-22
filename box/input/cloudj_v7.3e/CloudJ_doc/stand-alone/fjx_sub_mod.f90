!------------------------------------------------------------------------------
!     'fjx_sub_mod.F90'  for fast-JX code v7.3+                                     !
!------------------------------------------------------------------------------
!
! !MODULE: FJX
!
! !DESCRIPTION: JX version 7.2  (12/13) consistent with 7.1 data and results
!          variables in call to PHOTO_JX are same as in 7.1,
!          but a logical(out) LDARK is added to to count the number of J calcs
!          Also works with new ver 7.3+ data sets if _init_ is updated.
!
! !INTERFACE:
!
      MODULE FJX_SUB_MOD
!
! !USES:
!
      USE FJX_CMN_MOD

      IMPLICIT NONE
!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: SOLAR_JX, ACLIM_FJX, JP_ATM0, PHOTO_JX, EXITC


      CONTAINS

!-----------------------------------------------------------------------
      subroutine SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!-----------------------------------------------------------------------
!     GMTIME = UT for when J-values are wanted
!           (for implicit solver this is at the end of the time step)
!     NDAY   = integer day of the year (used for solar lat and declin)
!     YGRDJ  = laitude (radians) for grid (I,J)
!     XGDRI  = longitude (radians) for grid (I,J)
!
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in) ::   GMTIME,YGRDJ,XGRDI
      integer, intent(in) ::  NDAY
      real*8, intent(out) ::  SZA,COSSZA,SOLFX
!
      real*8  LOCT
      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
!
      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*CPI180)
      SOLDEK = asin(SINDEC)
      COSDEC = cos(SOLDEK)
      SINLAT = sin(YGRDJ)
      SOLLAT = asin(SINLAT)
      COSLAT = cos(SOLLAT)
!
      LOCT   = (((GMTIME)*15.d0)-180.d0)*CPI180 + XGRDI
      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/CPI180
!
      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*C2PI/365.d0))
!
      END SUBROUTINE SOLAR_JX


!-----------------------------------------------------------------------
      subroutine ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!-----------------------------------------------------------------------
!  Load fast-JX climatology for latitude & month given pressure grid
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)  :: YLATD
      integer, intent(in) :: MONTH, L1U
      real*8, intent(in),  dimension(L1U+1) :: PPP
      real*8, intent(out), dimension(L1U+1) :: ZZZ
      real*8, intent(out), dimension(L1U) :: TTT,DDD,OOO

      real*8, dimension(LREF) :: OREF2,TREF2
      real*8, dimension(LREF+1) :: PSTD
      integer  K, L, M, N
      real*8   DLOGP,F0,T0,PB,PC,XC,SCALEH

!  Select appropriate month
      M = max(1,min(12,MONTH))
!  Select appropriate latitudinal profiles
      N = max(1, min(18, (int(YLATD+99)/10 )))
      do K = 1,LREF
        OREF2(K) = OREF(K,N,M)
        TREF2(K) = TREF(K,N,M)
      enddo

!  Apportion O3 and T on supplied climatology z levels onto CTM levels +1
!  with mass (pressure) weighting, assuming constant mixing ratio and
!  temperature half a layer on either side of the point supplied.
!   PPP(L=1:L1_)=edge-pressure of CTM layer, PPP(L1_+1)=0 (top-of-atmos)
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
!     MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

!  Set up pressure levels for O3/T climatology - assume that value
!  given for each 2 km z* level applies from 1 km below to 1 km above,
!  so select pressures at these boundaries. Surface level values at
!  1000 mb are assumed to extend down to the actual PSURF (if > 1000)
           PSTD(1) = max(PPP(1),1000.d0)
           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
           DLOGP   = 10.d0**(-2.d0/16.d0)
      do K = 3,LREF
        PSTD(K) = PSTD(K-1)*DLOGP
      enddo
        PSTD(52)  = 0.d0
      do L = 1,L1U
        F0 = 0.d0
        T0 = 0.d0
        do K = 1,LREF
          PC   = min(PPP(L),PSTD(K))
          PB   = max(PPP(L+1),PSTD(K+1))
          if (PC .gt. PB) then
            XC = (PC-PB)/(PPP(L)-PPP(L+1))
            F0 = F0 + OREF2(K)*XC
            T0 = T0 + TREF2(K)*XC
          endif
        enddo
        TTT(L)  = T0
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
        OOO(L) = F0*1.d-6*DDD(L)
      enddo

!  Calculate effective altitudes using scale height at each level
        ZZZ(1) = 0.d0
      do L = 1,L1U-1
        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
      enddo
        ZZZ(L1U+1) = ZZZ(L1U) + ZZHT

      END SUBROUTINE ACLIM_FJX

!<<<<<<<<<<<<<<<<<<<end CTM-fastJX special call subroutines<<<<<<<<<<<<<




!<<<<<<<<<<<<<<<<<<<<<<<<begin fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<all outside calls go through PHOTO_JX<<<<<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
!  fastJX version 7.2 (f90) - Prather notes (Jan 2013)

!---calculates cloud optical depth in FJX-72 given water path & R-eff
!---assumes that clouds are 100% if in layer
!      IWP = ice water path (in layer, in cloud) in g/m**2
!      LWP = liquid water path (in layer, in cloud) in g/m**2
!      REFFI = effective radius of ice cloud (microns)
!      REFFL = effective radius of liquid cloud (microns)
!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)

!>>>R-effective determined by main code, not FJX
!   REFF determined by user - some recommendations below (from Neu & Prather)
!          REFFI is function of ice water content IWC (g/m3, 0.0001 to 0.1)
!             IWC = IWP / delta-Z (of layer in m, approx OK)
! prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
!>>>indices for cloud scattering determined by FJX core, not main code.
!   ice clouds pick heagonal or irregular:
!          NDXI = 6  ! ice hexag (cold)
!             if (TCLD .ge. 233.15) then
!          NDXI = 7  ! ice irreg
!             endif
!   liquid clouds scaled to R-eff 3 - 6 - 12 - 20 microns:
!          NDXC = 1
!            do I=2,4
!             if (REFFL .gt. 0.5*(RCC(I-1)+RCC(I))) then
!          NDXC = I
!             endif
!            enddo

!-----------------------------------------------------------------------
      SUBROUTINE PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ, PPP,ZZZ,TTT,DDD,  &
                           RRR,OOO, LWP,IWP,REFFL,REFFI,               &
                           AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)

!-----------------------------------------------------------------------
!
!  PHOTO_JX is the gateway to fast-JX calculations:
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-fo-year for sun distance, SZA (not lat or long)
!-----------------------------------------------------------------------
      implicit none

!---calling sequence variables
      integer, intent(in)                    :: L1U,ANU,NJXU
      real*8,  intent(in)                    :: U0,SZA,REFLB,SOLF
      logical, intent(in)                    :: LPRTJ
      real*8,  intent(in), dimension(L1U+1)  :: PPP,ZZZ
      real*8,  intent(in), dimension(L1U  )  :: TTT,DDD,RRR,OOO,  &
                                                LWP,IWP,REFFL,REFFI
      real*8,  intent(in), dimension(L1U,ANU):: AERSP
      integer, intent(in), dimension(L1U,ANU):: NDXAER

! reports out the JX J-values, upper level program converts to CTM chemistry J's
      real*8, intent(out), dimension(L1U-1,NJXU)::  VALJXX
      logical, intent(out)                   :: LDARK

!-----------------------------------------------------------------------
!--------key LOCAL atmospheric data needed to solve plane-parallel J----
!-----these are dimensioned JXL_, and must have JXL_ .ge. L_
      real*8, dimension(JXL1_+1) :: TTJ,DDJ,OOJ,ZZJ
      real*8, dimension(JXL1_+1) :: PPJ,RELH
      integer,dimension(JXL2_+1) :: JXTRA
      real*8, dimension(W_)      :: FJTOP,FJBOT,FSBOT,FLXD0,RFL
      real*8, dimension(JXL_, W_)   ::  AVGF, FJFLX
      real*8, dimension(JXL1_,W_)   ::  DTAUX, FLXD
      real*8, dimension(8,JXL1_,W_) ::  POMEGAX
      real*8, dimension(JXL1_)      ::  DTAU600
      real*8, dimension(8,JXL1_)    ::  POMG600
      real*8, dimension(W_,JXL1_)   ::  FFX
      real*8, dimension(W_,8)     ::  FFXNET
!---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8  FLXJ(JXL1_),FFX0,FXBOT,FABOT
      real*8             :: ODABS,ODRAY,ODI,ODL
      real*8             :: RFLECT,FREFS,FREFL,FREFI
      real*8             :: AMF2(2*JXL1_+1,2*JXL1_+1)
!------------key SCATTERING arrays for clouds+aerosols------------------
      real*8             :: OPTX(5),SSAX(5),SLEGX(8,5)
      real*8             :: OD(5,JXL1_),SSA(5,JXL1_),SLEG(8,5,JXL1_)
      real*8             :: OD600(JXL1_)
      real*8             :: PATH,RH,XTINCT
!------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,JXL_),VALJ(X_)
      real*8  FLXUP(W_),FLXDN(W_),DIRUP(W_),DIRDN(W_)
      real*8  VALJL(JXL_,X_) !2-D array of J_s returned by JRATET

      integer  LU,L2U, I,J,K,L,M,KMIE,KW,NAER,NDXI,NDXL, RATIO(W_)
      real*8   XQO3,XQO2,DTAUC,WAVE, TTTX
!-----------------------------------------------------------------------
      real*8   SWHEAT, SKPERD

      if (L1U .gt. JXL1_) then
        call EXITC(' PHOTO_JX: not enough levels in JX')
      endif

      LU = L1U - 1
      L2U = LU + LU + 2

      VALJXX(:,:) = 0.d0
      FFF(:,:) = 0.d0
      FREFI = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0

!---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
!                        or         99.0                 80 km
      if (SZA .gt. 98.d0) then
        LDARK = .true.
        return
      else
        LDARK = .false.
      endif

!---load the amtospheric column data
      do L = 1,L1U
        PPJ(L) = PPP(L)
        TTJ(L) = TTT(L)
        DDJ(L) = DDD(L)
        OOJ(L) = OOO(L)
      enddo
        PPJ(L1U+1) = 0.d0

!---calculate spherical weighting functions (AMF: Air Mass Factor)
      do L = 1,L1U+1
        ZZJ(L) = ZZZ(L)
      enddo

        RFLECT = REFLB

!-----------------------------------------------------------------------
      call SPHERE2 (U0,RAD,ZZJ,ZZHT,AMF2, L1U,JXL1_)
!-----------------------------------------------------------------------

!---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8))
!---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
      do L = 1,L1U
        do K=1,5
          OD(K,L)  = 0.d0
          SSA(K,L) = 0.d0
          do I=1,8
            SLEG(I,K,L) = 0.d0
          enddo
        enddo
      enddo

      if (LPRTJ) then
      write(6,*)'fast-JX-(7.2)---PHOTO_JX internal print: Clouds--'
      write(6,*) '  L, P1,P2, WPath, Reff, OD, Ndx'
      endif

      do L = 1,L1U
!---Liquid Water Cloud
        if (LWP(L) .gt. 1.d-5 .and. REFFL(L) .gt. 0.1d0) then
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
          ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L)
          NDXL = 1
          do I = 2,4
            if (REFFL(L) .gt. 0.5*(RCC(I-1)+RCC(I))) then
              NDXL = I
            endif
          enddo
          call OPTICL (OPTX,SSAX,SLEGX,  ODL,NDXL)
          do K = 1,5
            OD(K,L)  = OD(K,L)  + OPTX(K)
            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
            do I = 1,8
              SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
            enddo
          enddo
!>>>diagnostic print of cloud data:
          if (LPRTJ) then
            write(6,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
                  'Liq Cld',L,PPP(L),PPP(L+1),LWP(L),REFFL(L),ODL,NDXL
          endif
        endif

!---Ice Water Cloud
        if (IWP(L) .gt. 1.d-5 .and. REFFI(L) .gt. 0.1d0) then
          ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0)
          if (TTT(L) .ge. 233.15d0) then
            NDXI = 7  ! ice irreg
          else
            NDXI = 6  ! ice hexag (cold)
          endif
          call OPTICL (OPTX,SSAX,SLEGX,  ODI,NDXI)
          do K=1,5
            OD(K,L)  = OD(K,L)  + OPTX(K)
            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
            do I=1,8
              SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
            enddo
          enddo
!>>>diagnostic print of cloud data:
          if (LPRTJ) then
             write(6,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
                  'Ice Cld',L,PPP(L),PPP(L+1),IWP(L),REFFI(L),ODI,NDXI
          endif
        endif

!---aerosols in layer: check aerosol index
!---this uses data from climatology OR from current CTM (STT of aerosols)

!---FIND useful way to sum over different aerosol types!
          RH = RRR(L)
        do M = 1,ANU
          NAER = NDXAER(L,M)
          PATH = AERSP(L,M)

!---subroutines OPTICA & OPTICM return the same information:
!---  optical depth (OPTX), single-scat albedo (SSAX) and phase fn (SLEGX(8))
!---subs have slightly different inputs:
!---  PATH is the g/m2 in the layer, NAER in the cloud/aerosol index
!---  UMich aerosols use relative humidity (RH)
          if (PATH .gt. 0.d0) then
            if (NAER .gt.0) then
              call OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER)
            else
              call OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,-NAER)
            endif
            do K = 1,5
              OD(K,L)  = OD(K,L)  + OPTX(K)
              SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
              do I = 1,8
                SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
              enddo
            enddo
          endif
        enddo

        do K = 1,5
          if (OD(K,L) .gt. 0.d0) then
            SSA(K,L) = SSA(K,L)/OD(K,L)
            do I = 1,8
              SLEG(I,K,L) = SLEG(I,K,L)/OD(K,L)
            enddo
          endif
        enddo

!---Include aerosol with cloud OD at 600 nm to determine added layers:
        OD600(L) = OD(4,L)

      enddo   ! end of 'do L = 1,L1U'

!---when combining with Rayleigh and O2-O3 abs, remember the SSA and
!---  phase fn SLEG are weighted by OD and OD*SSA, respectively.
!---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
!       additonal levels at top of clouds (now uses log spacing)

!-----------------------------------------------------------------------
      call EXTRAL(OD600,L1U,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------

!---set surface reflectance
      RFL(:) = max(0.d0,min(1.d0,RFLECT))

!-----------------------------------------------------------------------
!---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
!-----------------------------------------------------------------------

      do K = 1,W_
      if (FL(K) .gt. 1.d0) then
        WAVE = WL(K)
!---Pick nearest Mie wavelength to get scattering properites------------
                               KMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255.d0 ) KMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355.d0 ) KMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) KMIE=4
        if( WAVE .gt. 800.d0 ) KMIE=5

!---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
!---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
        do L = 1,L1U
          TTTX     = TTJ(L)
          call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1) &
                        ,QO2(K,2),TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2) &
                        ,QO3(K,2),TQQ(3,2),QO3(K,3), LQQ(2))
          ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948d0
          ODRAY = DDJ(L)*QRAYL(K)

          DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

          do I = 1,8
            POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
          enddo
          POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0d0*ODRAY
          POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5d0*ODRAY
          do I = 1,8
            POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
          enddo
        enddo

      endif
      enddo

!-----------------------------------------------------------------------

      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
              AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)

!-----------------------------------------------------------------------

         FLXUP(:) = 0.d0
         DIRUP(:) = 0.d0
         FLXDN(:) = 0.d0
         DIRDN(:) = 0.d0
         FLXJ(:) = 0.d0
         FFX(:,:) = 0.d0
         FFXNET(:,:) = 0.d0

      do K = 1,W_
      if (FL(K) .gt. 1.d0) then

!direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!----     also at bottom (DN), does not include diffuse reflected flux.
        FLXUP(K) =  FJTOP(K)
        DIRUP(K) = -FLXD0(K)
        FLXDN(K) = -FJBOT(K)
        DIRDN(K) = -FSBOT(K)

        do L = 1,LU
          FFF(K,L) = SOLF*FL(K)*AVGF(L,K)
        enddo
        FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
        FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
        FREFS = FREFS + SOLF*FL(K)/WL(K)

!---for each wavelength calculate the flux budget/heating rates:
!  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L))]
!            but for spherical atmosphere!
!  FJFLX(L) = diffuse flux across top of layer L

!---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!---     need special fix at top and bottom:
!---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
        FABOT = (1.d0-RFL(K))*(FJBOT(K)+FSBOT(K))
        FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
        FLXJ(1) = FJFLX(1,K) - FXBOT
        do L = 2,LU
          FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
        enddo
        FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
!---calculate net flux deposited in each CTM layer (direct & diffuse):
        FFX0 = 0.d0
        do L=1,L1U
          FFX(K,L) = FLXD(L,K) - FLXJ(L)
          FFX0 = FFX0 + FFX(K,L)
        enddo

!  NB: the radiation level ABOVE the top CTM level is included in these budgets
!      these are the flux budget/heating terms for the column:
!  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!       these are surface fluxes to compare direct vs. diffuse:
!  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)

        FFXNET(K,1) = FLXD0(K)
        FFXNET(K,2) = FSBOT(K)
        FFXNET(K,3) = FLXD0(K) + FSBOT(K)
        FFXNET(K,4) = FJTOP(K)
        FFXNET(K,5) = FFX0
        FFXNET(K,6) = FABOT
        FFXNET(K,7) = FSBOT(K)
        FFXNET(K,8) = FJBOT(K)

!-----------------------------------------------------------------------
      endif
      enddo       ! end loop over wavelength K
!-----------------------------------------------------------------------
      FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
      FREFI = FREFI/FREFS

!---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

!-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
!-----------------------------------------------------------------------

!---mapping J-values from fast-JX onto CTM chemistry is done in main

!-----------------------------------------------------------------------
      if (LPRTJ) then
!---diagnostics below are NOT returned to the CTM code
       write(6,*)'fast-JX-(7.2)---PHOTO_JX internal print: Atmosphere--'
!---used last called values of DTAUX and POMEGAX, should be 600 nm
       do L=1,L1U
         DTAU600(L) = DTAUX(L,W_)
         do I=1,8
           POMG600(I,L) = POMEGAX(I,L,W_)
         enddo
       enddo

       call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)

!---PRINT SUMMARY of mean intensity, flux, heating rates:
       write(6,*)
       write(6,*)'fast-JX(7.2)---PHOTO_JX internal print: Mean Intens--'
       write(6,'(a,5f10.4)') &
       ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/', &
        RFLECT,SZA,U0,FREFI,FREFL

       write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
       write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
       write(6,'(a)') ' ---- 100000=Fsolar   MEAN INTENSITY per wvl bin'
         RATIO(:) = 0.d0
       do L = LU,1,-1
        do K=NW1,NW2
         if (FL(K) .gt. 1.d0) then
         RATIO(K) = (1.d5*FFF(K,L)/FL(K))
         endif
        enddo
         write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
       enddo

       write(6,*)
       write(6,*)'fast-JX(7.2)---PHOTO_JX internal print: Net Fluxes---'
       write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
       write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)

!  NB: the radiation level ABOVE the top CTM level is included in these budgets
!      these are the flux budget/heating terms for the column:
!  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!       these are surface fluxes to compare direct vs. diffuse:
!  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)


!      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
!      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
       write(6,*) ' ---NET FLUXES--- solar only to 850 nm'
       write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1)
       write(6,*) ' ---SRF FLUXES--- '
       write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1)

           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,3)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' sol TOTAL ',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,4)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' dif outtop',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,5)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' abs in atm',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,6)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' abs at srf',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,7)*FL(K)*FP(K)
         enddo
        write(6,'(a11,1p,e12.4)') ' PAR direct',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,8)*FL(K)*FP(K)
         enddo
        write(6,'(a11,1p,e12.4)') ' PAR diffus',SWHEAT

       write(6,'(2a)') ' absorb solar: W/m2/kyr K/day by wavel 10000=Fsolar', &
       '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
       do L = LU,1,-1
           SWHEAT = 0.d0
         do K=NW1,NW2
           RATIO(K) = 1.d5*FFX(K,L)
           SWHEAT = SWHEAT + FFX(K,L)*FW(K)
         enddo
           SKPERD = SWHEAT*86400.d0/(1005.d0*(PPP(L)-PPP(L+1))*100.d0/9.8d0)
         write(6,'(i9,2x,2f8.4,18i8)') L,SWHEAT,SKPERD,(RATIO(K),K=NW2,NW1,-1)
       enddo

       write(6,'(a)')
       write(6,'(a)') ' fast-JX (7.2)----J-values----'
       write(6,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
       do L = LU,1,-1
         write(6,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
       enddo

      endif

   99 continue

      END SUBROUTINE PHOTO_JX


!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
      subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
              FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in) ::   DTAUX(JXL1_,W_),POMEGAX(8,JXL1_,W_)
      real*8, intent(in) ::   AMF2(2*JXL1_+1,2*JXL1_+1)
      real*8, intent(in) ::   U0,RFL(W_)
      integer, intent(in) ::  JXTRA(JXL2_+1), LU
      real*8, intent(out) ::FJACT(JXL_,W_),FJTOP(W_),FJBOT(W_),FSBOT(W_)
      real*8, intent(out) ::  FJFLX(JXL_,W_),FLXD(JXL1_,W_),FLXD0(W_)

      integer JNDLEV(JXL_),JNELEV(JXL1_)
      integer JADDLV(JXL2_+1),JADDTO(JXL2_+1),L2LEV(JXL2_+1)
      integer JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
      integer L1U,L2U,   LZ0,LZ1,LZMID
      real*8   SUMT,SUMJ

      real*8  DTAU(JXL1_+1,W_),POMEGAJ(M2_,JXL2_+1,W_),TTAU(JXL2_+1,W_)
      real*8  FTAU2(JXL2_+1,W_),POMEGAB(M2_,W_)
      real*8  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
      real*8, dimension(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX
!--- variables used in mie code-----------------------------------------
      real*8, dimension(W_)         :: FJT,FJB
      real*8, dimension(N_,W_)      :: FJ,FZ,ZTAU
      real*8, dimension(M2_,N_,W_)  :: POMEGA
      real*8, dimension(2*JXL1_,W_)  :: FLXD2

!---there is a parallel correspondence:
!  dimension of JX arrays JXL_ .ge. dimension that CTM is using = L_
!  but calculation is done for L_=LU, L1_=L1U, L2_=L2U lengths of CTM
!
!  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
!
! in:
!     DTAUX(1:L1_,1:W_) = optical depth of each layer
!     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s abledo)
!     U0  = cos (SZA)
!     RFL(1:W_) = Lambertian albedo of surface
!     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
!        AMF2 now does both edges and middle of CTM layers
!     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
! out:
!     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
!  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
!     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
!     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
!     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
!     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
!        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
!     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
!        this should take into account sphericity, and is not just = mu0
!     FLXD0(1:W_) = sum of solar flux deposited in atmos
!        does NOT include flux on lower surface, does NOT mean absorbed!
!-----------------------------------------------------------------------
!
!     DTAU     Local optical depth of each CTM level
!     TTAU     Optical depth of air vertically above each point (to top of atm)
!     FTAU2     Attenuation of solar beam
!     POMEGAJ  Scattering phase function
!
!---new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
!   factor increase from sub-layer to sub-layer
!
!---------------------SET UP FOR MIE CODE-------------------------------
!
!-----------------wavelength independent--------------------------------
!
!  Transpose the ascending TTAU grid to a descending ZTAU grid.
!  Double the resolution - TTAU points become the odd points on the
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!  Odd point added at top of grid for unattenuated beam   (Z='inf')
!
!  The following mapping holds for JADDLV=0
!        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
!        Top:       TTAU(L2_)  ==> ZTAU(3)
!        Infinity:     0.0     ==> ZTAU(1)
!        index: 2*(L2_+1-L2)+1 ==> LZ
!
!  Mie scattering code only used from surface to level L2_
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere
!  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
!  to be incremented linearly, and the flux fz to be attenuated top-down
!    (avoiding problems where lower level fluxes are zero).
!------------------------------------------------------------------------
!
!  Ascend through atmosphere transposing grid and adding extra points
!  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
!  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
!    because we need to insert the intermediate layers (even LZ) for the
!    asymmetric scattering code.
!
!  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
!    order, expanded, doubled-level scatter grid.
!    Note that we need to deal with the expansion by JADD levels (L2L).
!      These JADDLV levels are skipped and need to be interpolated later.
!    Note that only odd LZ levels are filled,
!
!----------------------re-grid data---------------------------------------------
!  Calculate cumulative total and define levels we want J-values at.
!  Sum upwards for levels, and then downwards for Mie code readjustments.
!
!     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
!           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!     JADDLV(L2)  Number of new levels actually added at each wavelength
!            where JADDLV = 0 when there is effectively no FTAU2
!     JADDTO(L2)   Total number of new levels to add to and above level (L2)
!     JNDLEV(L) = L2 index that maps on CTM mid-layer L
!
!---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
!---    JADDLV is taken from JXTRA, which is based on visible OD.
!---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
!---these should be fixed for all wavelengths to lock-in the array sizes

      if (LU .gt. JXL_) then
        call EXITC (' OPMIE:  JXL_ .lt. L_')
      endif

      L1U = LU + 1
      L2U = 2*LU + 2


      do L2 = 1,L2U,1
        JADDLV(L2) = JXTRA(L2)
      enddo
        JADDTO(L2U+1) = 0
      do L2 = L2U,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      enddo

!---expanded grid now included CTM edge and mid layers plus expanded
!---    grid to allow for finer delta-tau at tops of clouds.
!---    DIM of new grid = L2U + JADDTO(1) + 1

!---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
!     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2U+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

!---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
!---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,LU
        JNDLEV(L) = L2LEV(2*L)
        JNELEV(L) = L2LEV(2*L+1)
      enddo
        JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

      ND = 2*L2U + 2*JADDTO(1) + 1

      if(ND .gt. N_) then
        call EXITC (' overflow of scatter arrays: ND > N_')
      endif

!----------------begin wavelength dependent set up------------------------------

!---Reinitialize arrays
      ZTAU(:,:)     = 0.d0
      FZ(:,:)       = 0.d0
      POMEGA(:,:,:) = 0.d0

      FJACT(:,:) = 0.d0
      FJTOP(:) = 0.d0
      FJBOT(:) = 0.d0
      FSBOT(:) = 0.d0
      FJFLX(:,:) = 0.d0
      FLXD(:,:) = 0.d0
      FLXD0(:) = 0.d0
      FJT(:) = 0.d0
      FJB(:) = 0.d0
      FJ(:,:) = 0.d0
      FZ(:,:) = 0.d0
      ZTAU(:,:) = 0.d0

      do K=1,W_
      if (FL(K) .gt. 1.d0) then

!---Set up optical depth DTAU(L)
       do L = 1,L1U
        DTAU(L,K) = DTAUX(L,K)
       enddo
        DTAU(L1U+1,K) = 0.d0

!---Define the total scattering phase fn for each CTM layer L=1:L_+1
!---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
!---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
        do I = 1,M2_
          POMEGAJ(I,L,K) = POMEGAX(I,L,K)
        enddo
       enddo

!---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
!---      at the middle & edges of the CTM layers L=1:2*L1_+1
!---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
!---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
        FTAU2(:,:) = 0.d0
        FTAU2(L2U+1,:) = 1.0d0
       do LL = 1,2*L1U+1
         L = (LL+1)/2
        if (AMF2(LL,LL) .gt. 0.0d0) then
           XLTAU = 0.0d0
         do II = 1,2*L1U+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5d0*DTAU(I,K)*AMF2(II,LL)
         enddo
         if (XLTAU .lt. 76.d0) then   ! zero out flux at 1e-33
          FTAU2(LL,K) = exp(-XLTAU)
         endif
        endif
       enddo

!---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
!---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:,:) = 0.d0
       do LL = 1,2*L1U
        if (AMF2(LL,LL) .gt. 0.d0) then
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
        endif
       enddo
        if (AMF2(1,1) .gt. 0.d0) then
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
        else
          FSBOT(K) = 0.d0
        endif

       do LL = 2,2*L1U,2
         L=LL/2
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       enddo

!---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
!---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0.d0
       if (AMF2(2*L1U,2*L1U) .gt. 0.d0) then
        do L=1,L1U
         FLXD0(K) = FLXD0(K) + FLXD(L,K)
        enddo
       endif

!------------------------------------------------------------------------
!  Take optical properties on CTM layers and convert to a photolysis
!  level grid corresponding to layer centres and boundaries. This is
!  required so that J-values can be calculated for the centre of CTM
!  layers; the index of these layers is kept in the JNDLEV array.
!------------------------------------------------------------------------
!---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
!---    points (1:L_) plus 1 for the mid point of added top layer.
!---combine these edge- and mid-layer points into grid of size:
!---              L2_+1 = 2*L1_+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2_+1)
!---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD

        TTAU(L2U+1,K) = 0.0d0
       do L2 = L2U,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAU(L,K)
        TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

!----solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.d0) then
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1.d0+RFL(K))
       else
        ZFLUX(K) = 0.d0
       endif

!  Calculate scattering properties, level centres then level boundaries
!>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ' upward
!     in index
       do L2 = L2U,2,-2
        L   = L2/2
        do I = 1,M2_
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
        enddo
       enddo
!---lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
         POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       enddo
!---now have POMEGAJ filled at even points from L2=3:L2_-1
!---use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
        do I = 1,M2_
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + &
                 POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
        enddo
       enddo

!---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
!---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
        do I=1,M2_
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
        enddo
       enddo

!   Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2U             ! L2 = index of CTM edge- and mid-layers
         L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
         LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
         L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

        if (L22 .gt. 0) then
          TAUBTM(K) = TTAU(L2,K)
          TAUTOP(K) = TTAU(L2+1,K)
          FBTM(K)   = FTAU2(L2,K)
          FTOP(K)   = FTAU2(L2+1,K)
         do I = 1,M2_
          POMEGAB(I,K) = POMEGAJ(I,L2,K)
         enddo

!---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
!---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
         ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
         do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

!---fraction from TAUBTM=>TAUTOP
          ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
!---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg),
!---else scale by TAU
          if (U0 .gt. 0.02d0) then
            FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
          else
            if (FTOP(K).lt.1.d-30 .or. FBTM(K).lt.1.d-30) then
              FZ(LZZ,K) = 0.d0
            else
              FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
            endif
          endif
          do I = 1,M2_
            POMEGA(I,LZZ,K) = POMEGAB(I,K) + &
                     ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
          enddo
            TAUBTM(K)    = ZTAU(LZZ,K)
            FBTM(K)      = FZ(LZZ,K)
          do I = 1,M2_
            POMEGAB(I,K) = POMEGA(I,LZZ,K)
          enddo
         enddo
        endif
       enddo

!   Now fill in the even points with simple interpolation in scatter arrays:
       do LZ = 2,ND-1,2
         ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
        do I=1,M2_
         POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
        enddo
       enddo

      endif
      enddo  ! wavelength loop!

!-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------

!---Move mean intensity from scatter array FJ(LZ=1:ND)
!---              to CTM mid-level array FJACT(L=1:L_)

      do K=1,W_
      if (FL(K) .gt. 1.d0) then

!---mean intensity at mid-layer:  4*<I> + solar
!       do L = 1,LU
!        L2L = JNDLEV(L)
!        LZ  = ND+2 - 2*L2L
!        FJACT(L,K) = 4.d0*FJ(LZ,K) + FZ(LZ,K)
!       enddo

!---mean intensity averaged throughout layer:
       do L = 1,LU
         LZ0 = ND+2 - 2*JNELEV(L)
        if (L .gt. 1) then
         LZ1 = ND+2 - 2*JNELEV(L-1)
        else
         LZ1 = ND
        endif
         SUMJ = (4.d0*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K)) &
               +(4.d0*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
         SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ+(4.d0*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
         SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
        enddo
        FJACT(L,K) = SUMJ/SUMT
       enddo

!---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
!---      average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
        L2L = JNELEV(L)
        LZ  = ND+2 - 2*L2L
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
        FJFLX(L,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0))
       enddo

!---diffuse fluxes reflected at top, incident at bottom
         FJTOP(K) = FJT(K)
         FJBOT(K) = FJB(K)

      endif
      enddo  ! wavelength loop!

      END SUBROUTINE OPMIE


!-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_) &
                             ,RFL(W_),U0,ZFLUX(W_)
      real*8, intent(out) ::  FJ(N_,W_),FJT(W_),FJB(W_)

      real*8  PM(M_,M2_),PM0(M2_)
      integer I, IM  ,K

!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!-----------------------------------------------------------------------
      do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
         PM(I,IM) = PM0(IM)
       enddo
      enddo

       call LEGND0 (-U0,PM0,M2_)
       do IM=1,M2_
         PM0(IM) = 0.25d0*PM0(IM)
       enddo

!---BLKSLV now called with all the wavelength arrays (K=1:W_)

      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

      END SUBROUTINE MIESCT


!-----------------------------------------------------------------------
      subroutine LEGND0 (X,PL,N)
!-----------------------------------------------------------------------
!---Calculates ORDINARY Legendre fns of X (real)
!---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: X
      real*8, intent(out) :: PL(N)
      integer I
      real*8  DEN
!---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo

      END SUBROUTINE LEGND0


!-----------------------------------------------------------------------
      subroutine BLKSLV &
           (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!-----------------------------------------------------------------------
!  Sets up and solves the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!  This goes back to the old, dumb, fast version 5.3
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_) &
                             ,PM(M_,M2_),PM0(M2_) &
                             ,RFL(W_),ZFLUX(W_)
      real*8, intent(out) ::  FJ(N_,W_),FJTOP(W_),FJBOT(W_)

      real*8, dimension(M_,N_,W_)    ::  A,C,H,   RR

      real*8, dimension(M_,M_,N_,W_) ::  B,AA,CC,  DD
      real*8, dimension(M_,M_) ::  E
      real*8  SUMB,SUMBX,SUMT
      integer I, J, K, L

      do K = 1,W_
      if (FL(K) .gt. 1.0d0) then
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K), &
           PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K), &
                   A(1,1,K),H(1,1,K),C(1,1,K), ND)
      endif
      enddo

      do K = 1,W_
      if (FL(K) .gt. 1.0d0) then
!-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K) &
                        -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
         enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        enddo

!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,L,K) = - E(I,J)*C(J,L,K)
         enddo
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

       enddo

!---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) &
           + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
           + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) &
           - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
           - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    +E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

!-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K) &
          + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
          + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        enddo
       enddo

!----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

!---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
!---FJBOT = scaled diffuse flux onto surface:
!---ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
!---SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2) &
            + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2) &
            + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4.d0*SUMB*RFL(K)/(1.0d0 + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4.d0*SUMT
       FJBOT(K) = 4.d0*SUMB - SUMBX

      endif
      enddo

      END SUBROUTINE BLKSLV


!-----------------------------------------------------------------------
      subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0 &
                    ,B,CC,AA,A,H,C,  ND)
!-----------------------------------------------------------------------
!  Generates coefficient matrices for the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_)
      real*8, intent(in)  ::  ZFLUX,RFL
      real*8, intent(in),dimension(N_) :: FZ,ZTAU

      real*8, intent(out),dimension(M_,M_,N_) ::  B,AA,CC
      real*8, intent(out),dimension(M_,N_) ::  A,C,H

      integer I, J, K, L1,L2,LL
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC
!
      real*8, dimension(M_,M_) :: S,T,U,V,W
!---------------------------------------------

!---------upper boundary:  2nd-order terms
       L1 = 1
       L2 = 2
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25d0*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        enddo
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1)
         A(I,L1) = 0.0d0
       enddo
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       enddo

!------------intermediate points:  can be even or odd, A & C diagonal
!---mid-layer h-points, Legendre terms 2,4,6,8
       do LL=2,ND-1,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4) &
         + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
          +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
         + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
          +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---------lower boundary:  2nd-order terms
       L1 = ND
       L2 = ND-1
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
          POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
          POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
          POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
          POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo

!------------lower boundary, 2nd-order, A-matrix is full (AA)
         DELTAU = ZTAU(L1) - ZTAU(L2)
         D2 = 0.25d0*DELTAU
         SURFAC = 4.0d0*RFL/(1.0d0 + RFL)
       do I = 1,M_
          D1 = EMU(I)/DELTAU
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
          SUM1 = SURFAC*SUM0
        do J = 1,M_
         AA(I,J,L1) = - D2*U(I,J)
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
        enddo
         H(I,L1) = H(I,L1) - 2.0d0*D2*C(I,L1) + SUM0*ZFLUX
       enddo

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0d0
       enddo

      END SUBROUTINE GEN_ID

!<<<<<<<<<<<<<<<<<<<<<<<<<<end fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<




!<<<<<begin fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!------------------------------------------------------------------------------
      subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!------------------------------------------------------------------------------
!---set CLOUD fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!      ----FJX ver 7.0+ clouds separated
! 01 W_C02 (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.55_rho=1.000
! 02 W_C04 (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.19_rho=1.000
! 03 W_C08 (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.1_rho=1.000
! 04 W_C13 (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.9_rho=1.000
! 05 W_L06 (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.9_rho=1.000
! 06 Ice-Hexagonal (Mishchencko)                  reff=50.00___G=999.9_rho=0.917
! 07 Ice-Irregular (Mishchencko)                  reff=50.00___G=999.9_rho=0.917
! 08 S-Bkg LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 09 S-Vol LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
!
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     ODCLD      ! optical depth of cloud layer @600 nm
      integer,intent(inout)::  NDCLD      ! index of cloud layer:  4:13

      integer I,J
      real*8  XTINCT, REFF,RHO

!---later versions should allow for interpolation, averaging of R-eff

!---default cloud type C1, Reff = 12 microns
      if (NDCLD .lt. 1 .or. NDCLD .gt. 9) then
         NDCLD = 3
      endif

!--rescale OD by Qext at 600 nm (J=4)
      do J=1,5
         OPTD(J) = ODCLD * QCC(J,NDCLD)/QCC(4,NDCLD)
         SSALB(J) = SCC(J,NDCLD)
        do I=1,8
         SLEG(I,J) =  PCC(I,J,NDCLD)
        enddo
      enddo

      END SUBROUTINE OPTICL


!------------------------------------------------------------------------------
      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,K)
!------------------------------------------------------------------------------
!---for the UCI aerosol data sets, calculates optical properties at fast-JX's
!              std 5 wavelengths:200-300-400-600-999nm
!
!---UCI aersols optical data  v-7.0+
!01 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435  reff=0.221___G=.0523_rho=1.630
!02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435  reff=0.386___G=.0721_rho=1.630
!03 UT-sulfate LOGN:r=0.05 s=.693 n=1.44          reff=0.166___G=.0205_rho=1.769
!04 UT-sulfate LOGN:r=0.05 s=.693 n=1.46          reff=0.166___G=.0205_rho=1.769
!05 UT-sulfatM LOGN:r=.050 s=.642 n=1.53          reff=0.140___G=.0179_rho=1.769
!06 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i    reff=0.140___G=.0179_rho=1.500
!07 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i    reff=0.150___G=.0332_rho=1.500
!08 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i reff=0.149___G=.0331_rho=1.230
!09 UM-FF04 (%BC)LOGN:r=.050 s=.642 n=1.541+0.02i reff=0.140___G=.0179_rho=1.212
!10 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i reff=0.140___G=.0179_rho=1.230
!11 MDust.15  (R.V. Martin generated phase fns)   reff=0.150___G=1.000_rho=2.600
!12 MDust.25  (R.V. Martin generated phase fns)   reff=0.250___G=1.000_rho=2.600
!13 MDust.40  (R.V. Martin generated phase fns)   reff=0.400___G=1.000_rho=2.600
!14 MDust.80  (R.V. Martin generated phase fns)  reff=0.800___G=1.000_rho=2.6001
!15 MDust1.5  (R.V. Martin generated phase fns)   reff=1.500___G=1.000_rho=2.600
!16 MDust2.5  (R.V. Martin generated phase fns)   reff=2.500___G=1.000_rho=2.600
!17 MDust4.0  (R.V. Martin generated phase fns)   reff=4.000___G=1.000_rho=2.600
!
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00+)
      integer,intent(inout)::     K       ! index of cloud/aerosols

      integer I,J
      real*8  XTINCT, REFF,RHO

      if (K.gt.NAA .or. K.lt.1)  &
        call EXITC ('OPTICA: aerosol index out-of-range')

         REFF = RAA(K)
         RHO = DAA(K)
      do J=1,5
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75d0*QAA(J,K)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SAA(J,K)
       do I=1,8
         SLEG(I,J) =  PAA(I,J,K)
       enddo
      enddo

      END SUBROUTINE OPTICA


!------------------------------------------------------------------------------
      subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,LL)
!------------------------------------------------------------------------------
!---for the U Michigan aerosol data sets, this generate fast-JX data formats.
!---NB Approximates the Legendre expansion(L) of the scattering phase fn
!---           as (2*L+1)*g**L
!---UMAER(I,J,K,L):
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
!   J=1:5 = [200, 300, 400, (550,) 600 , 1000 nm]
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC),
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00)
      integer,intent(in)::     LL         ! index of cloud/aerosols

      integer KR,J,L
      real*8  R,FRH, GCOS, XTINCT

!---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!---extrapolate phase fn from first term (g)
      L = LL

      if (L.lt.1 .or. L.gt.33)  &
          call EXITC ('OPTICM: aerosol index out-of-range')

!---pick nearest Relative Humidity
      KR =  20.d0*RELH  + 1.5d0
      KR = max(1, min(21, KR))

      do J=1,5
       SSALB(J) = UMAER(1,J,KR,L)
         XTINCT = UMAER(3,J,KR,L)
       OPTD(J) = PATH*XTINCT
         GCOS   = UMAER(2,J,KR,L)
       SLEG(1,J) =  1.d0
       SLEG(2,J) =  3.d0*GCOS
       SLEG(3,J) =  5.d0*GCOS**2
       SLEG(4,J) =  7.d0*GCOS**3
       SLEG(5,J) =  9.d0*GCOS**4
       SLEG(6,J) = 11.d0*GCOS**5
       SLEG(7,J) = 13.d0*GCOS**6
       SLEG(8,J) = 15.d0*GCOS**7
      enddo

      END SUBROUTINE OPTICM


!-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,NJXU)
!-----------------------------------------------------------------------
! in:
!        PPJ(L_+1) = pressure profile at edges
!        TTJ(L_+1) = = temperatures at mid-level
!        FFF(K=1:NW, L=1:L_) = mean actinic flux
! out:
!        VALJL(L_,JX_)  JX_ = no of dimensioned J-values in CTM code
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: LU,NJXU
      real*8, intent(in)  ::  PPJ(LU+1),TTJ(LU+1)
      real*8, intent(inout)  ::  FFF(W_,LU)
      real*8, intent(out), dimension(LU,NJXU) ::  VALJL

      real*8  VALJ(X_)
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV

      if (NJXU .lt. NJX) then
        write(6,'(A,2I5)')  'NJXU<NJX',NJXU,NJX
        call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
      endif
      do L = 1,LU
!need temperature, pressure, and density at mid-layer (for some quantum yields):
          TT   = TTJ(L)
         if (L .eq. 1) then
          PP = PPJ(1)
         else
          PP  = (PPJ(L)+PPJ(L+1))*0.5d0
         endif
          DD = 7.24e18*PP/TT

!---if W_=18/12, must zero bin-11/5 below 100 hPa:  216-222 & 287-291 nm mix
!   since O2 e-fold is too weak (by averaging the 287-291)
!   this combination of wavelengths works fine in stratosphere with O3 attenuation
        if (PP .gt. 100.d0) then
          if (W_ .eq. 18) then
            FFF(11,L) = 0.d0
          elseif (W_ .eq. 12) then
            FFF(5,L) = 0.d0
          endif
        endif

        do J = 1,NJX
          VALJ(J) = 0.d0
        enddo

        do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1), &
                         TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1), &
                         TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1), &
                         TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
            QO31D  = QO31DY*QO3TOT
           VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
           VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
           VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
        enddo

        do J = 4,NJX
         do K = 1,W_
!---also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
           call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
                     TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
           call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
                     TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          endif
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
         enddo
        enddo

        do J=1,NJX
          VALJL(L,J) = VALJ(J)
        enddo

      enddo

      END SUBROUTINE JRATET


!-----------------------------------------------------------------------
      subroutine X_interp (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!-----------------------------------------------------------------------
!  up-to-three-point linear interpolation function for X-sections
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)::  TINT,T1,T2,T3, X1,X2,X3
      integer,intent(in)::  L123
      real*8,intent(out)::  XINT

      real*8  TFACT

      if (L123 .le. 1) then
           XINT = X1
      elseif (L123 .eq. 2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
      else
        if (TINT.le. T2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
        else
             TFACT = max(0.d0,min(1.d0,(TINT-T2)/(T3-T2) ))
           XINT = X2 + TFACT*(X3 - X2)
        endif
      endif

      END SUBROUTINE X_interp


!-----------------------------------------------------------------------
      subroutine JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)

      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ,DTAU6
      real*8, intent(in), dimension(8,LU+1) :: POMEG6
      integer,intent(in), dimension(LU+LU+3) :: JXTRA
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP,DAIR,DOZO

      write(6,'(4a)') '   L z(km)     p      T   ', &
       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
       '  g(cos) CTM lyr=>'

      L = LU+2
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZZJ(L)*1.d-5,PPJ(L)

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)

        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
          DAIR = DDJ(L)/DELZ
          DOZO = OOJ(L)/DELZ

      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZKM,PPJ(L),TTJ(L),DAIR,DOZO,COLO2,COLO3,DTAU6(L), &
            POMEG6(1,L),POMEG6(2,L)/3.d0, JXTRA(L+L),JXTRA(L+L-1)

        enddo

      END SUBROUTINE JP_ATM


!-----------------------------------------------------------------------
      subroutine JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ, LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)

      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP
      write(6,'(4a)') '   L z(km)     p      T   ', &
       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
       '  g(cos) CTM lyr=>'
      L = LU+2
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZZJ(L)*1.d-5,PPJ(L)
          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)
        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
            COLO2,COLO3
        enddo

      END SUBROUTINE JP_ATM0


!-----------------------------------------------------------------------
      subroutine SPHERE2(U0,RAD,ZHL,ZZHT,AMF2, L1U,LJX1U)
!-----------------------------------------------------------------------
!----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!-----------------------------------------------------------------------
! in:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
! out:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!          these are calcualted for both layer middle and layer edge
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) ::   L1U, LJX1U
      real*8, intent(in)  ::   U0,RAD,ZHL(L1U+1),ZZHT
      real*8, intent(out) ::   AMF2(2*LJX1U+1,2*LJX1U+1)

      integer, parameter  ::  LSPH_ = 100

!     RZ      Distance from centre of Earth to each point (cm)
!     RQ      Square of radius ratios
!     SHADHT  Shadow height for the current SZA
!     XL      Slant path between points

      integer  I, J, K, II, L2
      real*8   XMU1,XMU2,XL,DIFF,SHADHT,RZ(LSPH_+1)
      real*8   RZ2(2*LSPH_+1),RQ2(2*LSPH_+1)

!--- must have top-of-atmos (NOT top-of-CTM) defined
!      ZHL(L1U+1) = ZHL(L1U) + ZZHT

      if (L1U .gt. LSPH_) then
        call EXITC(' SPHERE2: temp arrays not large enough')
      endif

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1U+1
        RZ(II)   = RAD + ZHL(II)
      enddo

!---calculate heights for edges of split CTM-layers
      L2 = 2*L1U
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1))
      enddo
        RZ2(L2+1) = RZ(L1U+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2
      enddo

!---shadow height for SZA > 90
      if (U0 .lt. 0.0d0)  then
        SHADHT = RZ2(1)/dsqrt(1.0d0 - U0**2)
      else
        SHADHT = 0.d0
      endif

!---up from the surface calculating the slant paths between each level
!---  and the level above, and deriving the appropriate Air Mass Factor
         AMF2(:,:) = 0.d0

      do 16 J = 1,2*L1U+1

!  Air Mass Factors all zero if below the tangent height
        if (RZ2(J) .lt. SHADHT) goto 16
!  Ascend from layer J calculating AMF2s
        XMU1 = abs(U0)
        do I = J,2*L1U
          XMU2     = dsqrt(1.0d0 - RQ2(I)*(1.0d0-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
        enddo
!--fix above top-of-atmos (L=L1U+1), must set DTAU(L1U+1)=0
          AMF2(2*L1U+1,J) = 1.d0
!
!  Twilight case - Emergent Beam, calc air mass factors below layer
        if (U0 .ge. 0.0d0) goto 16

!  Descend from layer J
          XMU1       = abs(U0)
         do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0d0-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
!  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ2(II))
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
            AMF2(II,J) = 2.d0*XL/(RZ2(II+1)-RZ2(II))
            XMU1      = XMU2
!  Lowest level intersected by emergent beam
          else
            XL        = RZ2(II+1)*XMU1*2.0d0
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
            goto 16
          endif
         enddo

   16 continue

      END SUBROUTINE SPHERE2


!-----------------------------------------------------------------------
      subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------
!
!    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol layers
!    this version sets up log-spaced sub-layers of increasing thickness ATAU
!
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!
!     N.B. the TTAU, etc calculated here are NOT used elsewhere

!---The log-spacing parameters have been tested for convergence and chosen
!---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!---  ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!-----------------------------------------------------------------------
!
      implicit none
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
!
      integer JTOTL,I,L,L2
      real*8  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
!
!---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
!
!---combine these edge- and mid-layer points into grid of size:
!---              L2X+1 = 2*L1X+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2X+1)
!---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
!
!---Divide thick layers to achieve better accuracy in the scattering code
!---In the original fast-J, equal sub-layers were chosen, this is wasteful
!---and this new code (ver 5.3) uses log-scale:
!---        Each succesive layer (down) increase thickness by ATAU > 1
!---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
!---        4 sub-layers with ODs = 1 - 2 - 4 - 8
!---The key parameters are:
!---        ATAU = factor increase from one layer to the next
!---        ATAUMN = the smallest OD layer desired
!---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
!---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
!---Now compute the number of log-spaced sub-layers to be added in
!---   the interval TTAU(L2) > TTAU(L2+1)
!---The objective is to have successive TAU-layers increasing by factor ATAU >1
!---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'EXTRAL: N_/L2_/L2-cutoff JXTRA:',  &
                                   NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      END SUBROUTINE EXTRAL

!<<<<<<<end fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!-----------------------------------------------------------------------
      subroutine EXITC(T_EXIT)
!-----------------------------------------------------------------------
      character(len=*), intent(in) ::  T_EXIT

      write(6,'(a)') T_EXIT
      stop

      END SUBROUTINE EXITC


      END MODULE FJX_SUB_MOD
