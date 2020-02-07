!> MODULE DefPhotolysis_mod.f90 - A  crude version of the EMEP model's.
!!***************************************************************************! 
!> Calculates Photolysis coefficients
!! Currently using "phux" method as documented in M. Kuhn et al., Atmos. 
!! Environ., 32, No.4, 693-709, 1998.
!! A second "MCM" method will be introduced soon. 
!! Compared to EMEP code, this routine calculates indices such as IDNO2, and
!! contains more descriptive structures.
!! 
!!----------------------------------------------------------------------------

  module DefPhotolysis_mod
!-----------------------------------------------------------------------------

   use CB6_Photolysis
   implicit none
   private

 !/ subroutines: 

  public  :: setphotorates
  public  :: GetPhotol      !! 
  public  :: GetPhotol2     !!  newer version to deal with array
  public  :: ZenithAngle    !! For testing this module as stand-alone
  public  :: self_test      !! For testing this module

!  private :: InitMCM        !! 
!  public  :: InitPhotol     !! calls InitPhux or InitMCM (not yet)

 !/ parameters and variables: 

  character(len=10), public, save :: PhotolMethod = "MCM"
  real, private, parameter :: PI = 4.0*atan(1.0)

  integer, public, parameter :: &
       NRCPHOT      = 35  &! Number of pre-defined photolytic reactions
      ,MAXRCPHOT    = 320   ! Biggest index. Now for CB6
     !CB6& ,MAXRCPHOT    = 57   ! Biggest index. TMP !HI 55 -> 57
     !EMEP   NRCPHOT      = 17   ! Number of pre-defined photolytic reactions

  integer, public, save :: nPhotol = 0
   
! In EMEP code, we use IDBO3 for O1D production and IDAO3 for OP production
! MCM indices:

! If not existsing in MCM (eg HO2NO2) we use one of the empty MCM indices
  integer, private, parameter :: IDNOTSET = 10

  integer, public, parameter :: & 
    IDO3_O1D   = 1,IDO3_O3P  = 2,IDH2O2    = 3,IDNO2     = 4 ,IDNO3_NO  = 5 &
   ,IDNO3_NO2  = 6,IDHONO    = 7,IDHNO3    = 8,IDHCHO_H  = 11,IDHCHO_H2 = 12&
   ,IDCH3CHO  = 13 &
   ,MCM_J17   = 17, MCM_J18   = 18, MCM_J20   = 20 &
   ,MCM_J22   = 22 , IDMEK     = 22 &
   ,MCM_J23   = 23   & ! DS for EmChem18???
   ,IDACETON  = 21   & !
   ,IDCH3COY  = 35   & ! Based on CH3COY = BIACET (in OSLO CTM3)
! The different HCOHCO rates cannot be used with EMEP yet since PHODIS didn't
! distinguish
! Will use IDHCOHCO = 31 and scale overall rate to get 1.9 CO + ...
! See NotesDJ.py
   ,IDCHOCHO_2CO  = 31 &! 
   ,IDCHOCHO_HCHO = 32 &! 
   ,IDCHOCHO_2CHO = 33 &!
!78---------------------------------------------------------------------------
   ,IDCHOCHO      = 40 &! !! FAKE: FROM scaling of 31, later in code!!
                        !  Tot/31 = 1.5 !!
   ,IDNO3         = 42 &! !! "fake" used to set a total jNO3 (since EMEP only
                          ! has one jNO3) -- assume 1.145643*IDNO3_NO2
!
   ,IDCH3COCH3  = 21, IDRCOCHO  = 34 &!  MGLOX, etc
   ,IDCH3O2H  = 41 &
   ,IDiC3H7ONO2  = 54 &! Used for ISON
   ,IDCH3COO2H  = -1 &! QUERY ???
   ,IDHO2NO2    = 25 &! IDHO2NO2 not in MCM -- use scaled MEK photolysis rate (assuming jHO2NO2 can be approximated by 1.8*jMEK)
   ,IDN2O5    = 9 ! IDN2O5 not in MCM -- use scaled H2O2 photolysis rate (assuming jN2O5 can be approximated by 7*jH2O2)

!PHUX  !MCM coeffs: (not same as EMEP)
!PHUX    integer, public, parameter ::  &
!PHUX      IDAO3    =  1 , IDBO3    =  2 , IDNO2    =  4 , &
!PHUX!TMPTEST      IDAO3    =  2 , IDBO3    =  1 , IDNO2    =  4 , &
!PHUX      IDH2O2   =  3 , IDHNO3   =  8 , IDACH2O  =  11, &
!PHUX      IDBCH2O  = 12 , IDCH3CHO = 13 , IDCH3COX = 22 , &
!PHUX      !TMP IDCH3COY = -1 , IDHCOHCO = -1 , IDRCOHCO = 34 , &
!PHUX      IDCH3COY = -1 , IDRCOHCO = 34 , &
!PHUX      IDANO3   =  6 , IDN2O5   = -1 , IDCH3O2H = 41 , &
!PHUX      IDHO2NO2 = 16 , IDACETON = 17, &
!PHUX      IDAGLYOX = 31, IDBGLYOX = 32, IDCGLYOX = 33,&!RBMCM
!PHUX      IDBNO3   =  5                         !RBMCM, ca. 10% of IDNO3
!PHUX
!PHUX  !ESXTMP- FIX LATER
!PHUX    integer, public, parameter ::  &
!PHUX      IDNO3    =  IDANO3 !!!??????   IDANO3 -> NO2+ O3P
!PHUX
!PHUX    integer, public, parameter ::  IDRCOCHO  = IDRCOHCO ! Just tmp
!PHUX    integer, public, parameter ::  IDHCOHCO  = IDRCOHCO ! Just tmp

 !!--------------------------------------------------------------------------
 !! Phux system, from chemical comparison exercise of  Kuhn et al., AE, 1998
 !! which in turn was based upon Roeths program and RADM species.

  real, private, parameter :: MINYZ=-30.
  real, private, parameter :: EMINYZ=exp(MINYZ)
!  integer, public, save :: IDO3_O3P, IDO3_O1D, &
!     IDNO3_NO, IDNO3_NO2, &
!     IDHONO, IDHCHO_H2, IDHCHO_H, IDCH3CHO, IDCH3O2H, IDH2O2, IDNO2, IDHNO3
!  integer, public, save :: IDCH3COO2H,IDCH3COCH3,IDCHOCHO,IDRCOCHO,&
!    IDN2O5, IDMEK &
!   ,IDiC3H7ONO2       &! Used for ISON
!   ,IDCHOCHO_2CO      &! .. 31
!   ,IDCHOCHO_HCHO     &! .. 32
!   ,IDCHOCHO_2CHO      ! .. 33
!  integer, public, save :: MCM_J18,MCM_J20, MCM_J22   ! new labels

 !! From MCM system:

  type, private :: j_t
    integer :: mcmJ  ! Number in MCM, see Saunders et al, 2003
    real    :: L, M, N
    real    :: exj  ! example rate, mid latitudes, summer -- based on JJA
                    ! average sza for 45N (56 degrees, from Li 2017) 
    character(len=90) :: reaction ! or comment
  end type j_t


!78chars----------------------------------------------------------------------
!MJ = corrected to Mike Jenkin's rates of Feb 2017.
  type(j_t), dimension(61) :: dj = (/ &
   j_t(  1, 6.073E-05, 1.743, 0.474, 9.447E-06, 'O3 = O1D + O2' ) &
  ,j_t(  2, 4.775E-04, 0.298, 0.080, 3.480E-04, 'O3 = O3P + O2' ) &
  ,j_t(  3, 1.041E-05, 0.723, 0.279, 4.152E-06, 'H2O2 = 2 OH' ) &
  ,j_t(  4, 1.165E-02, 0.244, 0.267, 6.271E-03, 'NO2 = NO + O3P')  &
  ,j_t(  5, 2.485E-02, 0.168, 0.108, 1.858E-02, 'NO3 = NO + O2' ) &
  ,j_t(  6, 1.747E-01, 0.155, 0.125, 1.277E-01, 'NO3 = NO2 + O3P')  &
  ,j_t(  7, 2.644E-03, 0.261, 0.288, 1.357E-03, 'HONO = OH + NO' ) &
  ,j_t(  8, 9.312E-07, 1.230, 0.307, 2.631E-07, 'HNO3 = OH + NO2')  &
  ,j_t(  9, 7.287E-05, 0.723, 0.279, 4.152E-06, 'BoxJN2O5 replaced by 7*JH2O2' ) &
  ,j_t( 10, 0.0,  0.0, 0.0, 0.0, 'NOT SET')  & !!!USED for IDNOTSET!
!78---------------------------------------------------------------------------
!carbonyls
  ,j_t( 11, 4.642E-05, 0.762, 0.353, 1.586E-05, 'HCHO = HCO + H' ) &
  ,j_t( 12, 6.853E-05, 0.477, 0.323, 2.915E-05, 'HCHO = CO + H2' ) &
  ,j_t( 13, 7.344E-06, 1.202, 0.417, 1.732E-06, 'CH3CHO = HCO + CH3')  &
  ,j_t( 14, 2.879E-05, 1.067, 0.358, 8.163E-06, 'C2H5CHO = HCO + C2H5' ) &
  ,j_t( 15, 2.792E-05, 0.805, 0.338, 9.554E-06, 'nC3H7CHO = HCO + nC3H7')  &
  ,j_t( 16, 1.675E-05, 0.805, 0.338, 5.732E-06, 'nC3H7CHO = CH3CHO + C2H4') &
  ,j_t( 17, 7.914E-05, 0.764, 0.364, 2.647E-06, 'iC3H7CHO = HCO + nC3H7')  &
!MJ,j_t(18,1.140E-05,0.396,0.298,2.44E-06,'CH2=C(CH3)CHO = CH3C=CH2 + HCO')&
!MJ,j_t(19,1.140E-05,0.396,0.298,2.44E-06,'CH2=C(CH3)CHO = CHC(CH3)CO + H')&
  ,j_t(18, 1.482e-6,  0.396, 0.298, 6.909E-7,'CH2=C(CH3)CHO = CH3C=CH2+HCO')&
!QUERy 18=19?
  ,j_t(19, 1.482e-6,  0.396, 0.298, 6.909E-7,'CH2=C(CH3)CHO = CHC(CH3)CO+H')&
  ,j_t(20, 7.600E-4,  0.396, 0.298, 3.543E-4,'C5HPALD1, MJ Oct 2015')&
  ,j_t(21, 7.992E-07, 1.578, 0.271, 1.967E-7,'CH3C(O)CH3 = CH3CO + CH3')&
  ,j_t(22, 5.804E-06, 1.092, 0.377, 1.568E-6,'CH3C(O)C2H5 = CH3CO + C2H5')&
!78---------------------------------------------------------------------------
!MJ,j_t(23,1.836E-5,0.395,0.296,3.963E-6,'CH3C(O)CH=CH2 = CH3CH=CH2 + CO')&
!MJ,j_t(24,1.836E-5,0.395,0.296,3.963E-6,'CH3C(O)CH=CH2 = CH3CO + CH=CH2')&
!QUERy 23=24?
  ,j_t(23, 2.4246e-6, 0.395, 0.296, 1.135E-6, 'CH3C(O)CH=CH2 = CH3CH=CH2+CO')&
  ,j_t(24, 2.424e-6,  0.395, 0.296, 1.135E-6, 'CH3C(O)CH=CH2 = CH3CO+CH=CH2')&
  ,j_t(25, 1.045E-05, 1.092, 0.377, 1.568E-6,'BoxJHO2NO2 replaced by 1.8*JMEK')&
  ,j_t(26, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(27, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(28, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(29, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(30, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
!alpha-dicarbonyls
! The emep code has just one compoiste rate for HCOHCO. We fake this as mcm40 below
  ,j_t(31, 6.845E-05, 0.130, 0.201, 4.430E-05, '(CHO)2 = 2 CO + H2')  &
  ,j_t(32, 1.032E-05, 0.130, 0.201, 6.680E-06, '(CHO)2 = CO + HCHO')  &
  ,j_t(33, 3.802E-05, 0.644, 0.312, 1.497E-05, '(CHO)2 = 2 HCO')  &
!78---------------------------------------------------------------------------
   ! MGLYOX, RCHOCHO.. eqn wrong!
  ,j_t(34, 1.537E-04, 0.170, 0.208, 9.599E-05, 'CH3C(O)CH3 = CH3CHO + HCO') &
  ,j_t(35, 3.326E-04, 0.148, 0.215, 2.078E-04, 'CH3C(O)C(O)CH3 = 2 CH3CHO') &
  ,j_t(36, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(37, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(38, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(39, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
!  ,j_t(40, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
!78---------------------------------------------------------------------------
  ,j_t(40, 10.2675E-05, 0.130, 0.201, 6.646E-05,'EMEP MERGE HCOHCO=1.5*J31') &
!hydroperoxides
  ,j_t(41, 7.649E-06, 0.682, 0.279, 3.124E-06, 'CH3OOH = CH3O + OH')  &
!  ,j_t( 42, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ! EMEP jNO3 assumed to be 1.145643*J6:
  ,j_t(42, 0.200144,  0.155, 0.125, 1.463E-01, &
                                     'EMEP jNO3 assumed to be 1.145643*J6') &
  ,j_t(43, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(44, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(45, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(46, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(47, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(48, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(49, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t(50, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
!organic nitrates
  ,j_t(51, 1.588E-06, 1.154, 0.318, 4.598E-07, 'CH3ONO2 = CH3O + NO2')  &
  ,j_t(52, 1.907E-06, 1.244, 0.335, 5.083E-07, 'C2H5ONO2 = C2H5O + NO2')  &
  ,j_t(53, 2.485E-06, 1.196, 0.328, 6.897E-07, 'nC3H7ONO2 = nC3H7O + NO2')  &
  ,j_t(54, 4.095E-06, 1.111, 0.316, 1.220E-06, 'iC3H7ONO2 = iC3H7O + NO2')  &
  ,j_t(55, 1.135E-05, 0.974, 0.309, 3.708E-06, 'iC4H9ONO2 = iC4H9O + NO2')  &
!78---------------------------------------------------------------------------
!MJ,j_t(56,7.549E-06,1.015,0.324,6.787E-7,'CH3C(O)CH2ONO2=CH3C(O)CH2O+NO2')&
!RBBUG!,j_t(56,4.36e-6,1.089,0.323,6.787E-7,'CH3C(O)CH2ONO2=CH3C(O)CH2O+NO2')&
  ,j_t( 56, 4.365e-5,1.089,0.323,1.301E-5,'CH3C(O)CH2ONO2 = CH3C(O)CH2O+NO2')&
   !?? Not in MCM??:
  ,j_t( 57, 3.363E-6,1.296,0.322,8.902E-7,'CH3C(O)CH2ONO2 = CH3CO + HCHO+NO2')& 
  ,j_t( 58, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t( 59, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t( 60, 0.0,   0.0,   0.0,   0.0,     'NOT SET')  &
  ,j_t( 61, 7.537E-04, 0.499, 0.266, 3.505E-04, '??')  &
  /)

 contains
 !>--------------------------------------------------------------------------

 !> Simple solar values for testing this module

    function ZenithAngle( doy, hr, lat, lon ) result (theta)
      integer, intent(in) :: doy !! Day of year (1-366)
      real, intent(in) :: hr, lat, lon
      real :: phi, lha, decl, loc_hr
      real :: theta

      loc_hr = hr + lon/15.   !! local time
      lha = (1.0+loc_hr/12)*PI    !! local hour angle
      phi = lat *PI/180.0     !! latitude in radians

     !! declination angle (varies daily)

      decl=-23.45*PI/180.0*cos(2.0*PI/365.0*(doy+10.0))

      theta = acos(cos(lha)*cos(decl)*cos(phi)+sin(decl)*sin(phi))

      !print *, "ZA ", theta*180./PI, decl*180./PI
    end function ZenithAngle

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function GetPhotol(idj,CosZ,debug_level) result (J)
    integer, intent(in) :: idj
    real,    intent(in) :: CosZ  ! Zenith angle  ! , radians
    integer, intent(in) :: debug_level
    real :: J, tmpJ= -999., Zen
    logical :: first_call = .true.
    character(len=*),parameter :: dtxt = 'GetPhotol:'
    

    if(idj==0 ) write(*,"(a,i4,2x,9g12.2)") "ZEROPHOTOL ", &
     idj, CosZ, acos(CosZ)*180/PI

    if( CosZ < 1.0e-3 ) then
       J =  0.0
       return
    end if

     if ( PhotolMethod == 'CB6' ) then 
        if ( first_call ) then
           write(*,*) dtxt//"CB6 Photolysis1 Triggered"
           call init_cb6photol()
        end if

        Zen = acos(CosZ)* 180/PI
        J   = cb6photol(1,Zen)

     else
      associate ( L=> dj(idj)%L, M=> dj(idj)%M, N => dj(idj)%N ) 

        J = L * (cosZ**M) * exp(-n/CosZ)

      end associate
      tmpJ = dj(idj)%exj
     end if

    if(debug_level>0.and.first_call ) write(*,"(a,i4,2x,9g12.2)") "GETPHOTOL ", &
     idj, dj(idj)%exj, CosZ, acos(CosZ)*180/PI, J
    first_call = .false.

  end function  GetPhotol
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function GetPhotol2(indices,CosZ,debug_level) result (Jvals)
    integer, intent(in), dimension(:) :: indices
    real,    intent(in) :: CosZ  ! Zenith angle  ! , radians
    integer, intent(in) :: debug_level
    real, dimension(size(indices)) :: Jvals
    character(len=*),parameter :: dtxt = 'GetPhotol2:'

    logical, save :: first_call = .true.
    real    ::  Zen
    integer :: idj, iphot

    if ( first_call ) then
         write(*,*) dtxt//" Photolysis1 Method" // trim(PhotolMethod)
         call init_cb6photol()
    end if
!    if(idj==0 ) write(*,"(a,i4,2x,9g12.2)") "ZEROPHOTOL ", &
!     idj, CosZ, acos(CosZ)*180/PI

     if( CosZ < 1.0e-3 ) then
       Jvals = 0.0
       return
     end if 

     Zen = acos(CosZ)* 180/PI

     if ( PhotolMethod == 'CB6' ) then 
        if ( first_call ) then
           write(*,*) "CB6 Photolysis2 Triggered, "//&
               trim(PhotolMethod)
           call init_cb6photol()
        end if

        Jvals = cb6photol(indices,Zen)
        
     else
        do iphot = 1, size(indices) 
           idj = indices(iphot)
           associate ( L=> dj(idj)%L, M=> dj(idj)%M, N => dj(idj)%N ) 

             Jvals(iphot) = L * (cosZ**M) * exp(-n/CosZ)
           end associate
      end do
     end if

    if( ( debug_level>0.and.first_call ) .or.  ( debug_level>1 ) ) then
       write(*,"(a,2x,2f7.3,9g10.2)") "GETPHOTOL2 "//trim(PhotolMethod),&
         CosZ, Zen, Jvals(1:5)
    end if
    first_call = .false.

  end function  GetPhotol2
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    function setphotorates( ind, cosZ ) result(J) ! crude for ESX
      integer, intent(in) :: ind
      real, intent(in) :: cosZ !! cosine Solar zenith angle
      real :: J    !! Photorate
      ! integer, intent(in) :: debug_level
      integer, dimension(MAXRCPHOT), save :: mapindex = -1
      integer :: imap

      if( mapindex(ind) < 1 ) then
      MAPPING:   do imap = 1, size(dj%mcmJ)
          if ( dj(imap)%mcmJ == ind ) then
             mapindex(ind) = imap
             exit MAPPING
          end if
        end do MAPPING
        print *, "DJ MAP FAILED", ind
        stop "DJ MAP FAILED"
      end if

      imap = mapindex(ind)

      if( imap < 1 ) then
        print *, "DJ MAP NEG", ind
        stop "DJ MAP NEG"
      end if
        
      associate ( L=> dj(imap)%L, M=> dj(imap)%M, N => dj(imap)%N ) 

        if( cosZ > 1.0e-30 ) then
!Bug!          J = L * cosZ**(M*exp(-n/cosZ))
              J = L * (cosZ**(M)) * exp(-n/cosZ)
        else
          J =  0.0
        end if
      end associate

    end function setphotorates
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !> Subroutine to run abovce calculations, and also demonstrate call-order

  subroutine self_test()
    real :: Za,cosZa, J  !! Solar zenith angle in radians, +cos(chi)
    integer :: idj
    integer, dimension(4), save  :: photol_mcm, photol_cb6
    real,    dimension(4), save  :: JvalsCB6, JvalsMCM

    Za =  ZenithAngle( doy=181, hr=12.0, lat=52.0, lon = 15.0 )
    cosZa = cos(Za)

   ! call InitPhotol()

    photol_mcm  = [ IDO3_O1D , IDO3_O3P , IDH2O2 , IDNO2 ] !MCM
    photol_cb6  = [ 9        , 8        , 21     , 1     ] !CB6

    print *, "Zenith Angle ",Za*180./pi

    ! 1) All MCM rates
    do idj = 1, size(dj) 
      if( dj(idj)%reaction == 'NOT SET') cycle
      J = GetPhotol(idj,CosZa,debug_level=0)
      print "(a,i3,es12.3,3x,a)", "MCM J ",idj, J, trim(dj(idj)%reaction)
    end do

    ! 2) Compare
    
    call init_cb6photol

    PhotolMethod = 'MCM'
    JvalsMCM = GetPhotol2(photol_mcm,CosZa,debug_level=0)
    PhotolMethod = 'CB6'
    JvalsCB6 = GetPhotol2(photol_cb6,CosZa,debug_level=0)

    do idj = 1, size(photol_cb6)
      print '(a,3i4,4x,a14,2x,2es10.2)', "MCM vs CB6 Jval ", idj, &
        photol_mcm(idj), photol_cb6(idj), &
        trim(dj(idj)%reaction), JvalsMCM(idj), JvalsCB6(idj)
    end do

  end subroutine self_test

end module DefPhotolysis_mod
!-----------------------------------------------------------------------------

!TSTEMX ! This code can be tested with scripts/mk.testESXmods DefPhotolysis_mod.f90
!TSTEMX program test_DefPhotolysis_mod
!TSTEMX  use DefPhotolysis_mod
!TSTEMX  implicit none
!TSTEMX  call self_test()
!TSTEMX!  print *, "NO2->NO+O3P"
!TSTEMX!  print *, "MCM ", setphotorates( IDNO2, cosZa, 10 )
!TSTEMX!!
!TSTEMX!  print *, "NO3->NO+O2"
!TSTEMX!  print *, "MCM ", setphotorates( IDBNO3, cosZa, 10 )
!TSTEMX!
!TSTEMX!  print *, "NO3->NO2+O3P"
!TSTEMX!!  print *, "PHX ", Phux( 2.73e-1, 0.29327, 0.92401, Za)
!TSTEMX!  print *, "MCM ", setphotorates( IDNO3, cosZa, 10 )
!TSTEMX!
!TSTEMX!  print *, "O3->O1D"
!TSTEMX!  !print *, "MCMA", setphotorates( IDAO3, cosZa, 10 
!TSTEMX!  print *, "MCMB", setphotorates( IDBO3, cosZa, 10 )
!TSTEMX!
!TSTEMX!  print *, "O3->O3P"
!TSTEMX end program test_DefPhotolysis_mod
