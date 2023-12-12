module Config_module
   
  use AeroConstants_mod,     only: AERO
  use AeroFunctions_mod,     only: cMolSpeed
  use ChemDims_mod
  use ChemSpecs_mod,         only: species, CM_schemes_ChemSpecs
  use Debug_module,          only: debug
  use DefPhotolysis_mod,     only: PhotolMethod
  use EmisSplit_mod,         only: NSECTORS
  use MetFunctions,          only: rh2num, num2rh
  use OutputConcs_mod,       only: OutSpecsList_t
  use OwnDataTypes_mod,      only: TXTLEN_SHORT
  use PhysicalConstants_mod, only: AVOG, RGAS_J
  use SmallUtils_mod,        only: find_index
  use ZchemData_mod                ! fIsop, etc.
   
  implicit none

  private 

  public :: Init_Box

  real, public, save :: ppb = 2.55e+10
  ! also need new types to read in initial concentrations of specific
  !  species, and emissions
  type, public :: initBox_t
    character(len=20) :: species = '-' ! name of the species
    real :: init_conc = -999. ! initial concentration of this species
  end type initBox_t

  type, private ::   box_emis_t
     character(len=4) :: poll = '-'  ! eg nox, sox, voc ... as in EMIS_File
     real             :: emis = 0.0
  end type  box_emis_t

  !! ================================= eg EMIS_File = (/ 'sox', ... /) !!
    include 'CM_EmisFile.inc'                                          !! 
  !! ==================================================================!!

  real, public, dimension(NEMIS_File), save :: box_emis = 0.0 ! emissions
  real, public, dimension(NSPEC_TOT,NSECTORS), save :: &
       emis_spec ! emissions by species and sector, molec/cm2/s
  real, public, dimension(NEMIS_File,NSECTORS), save :: &
       emis_snap ! emiss
  real, public, save ::  Hmix = 1.0e5 ! cm, for dispersion of emissions

  ! ------------------------------------------------
  ! these variables will be set through the namelist
  ! cf box_config/ tstart, tend, dt, doy, lat, lon, use_emis, &
  !    emissplit_dir, m, h2o, emis_kgm2day, wanted_spec, all_species, &
  !    initBox, dbgVOC,Temp, rh

  integer, public :: nspec_out = 0
  real, public, save :: TSTART = 0.
  real, public, save :: TEND = 24*3600.
  real, public, save :: DT = 30.
  integer, public, save :: doy = 182

  real, public, save :: lat = 45.
  real, public, save :: lon = 12.

  type, private :: box_useconfig  ! short version of emep_useconfig
 ! N2O5 hydrolysis
 ! During 2015 the aersol surface area calculation was much improved, and this
 ! leads to the need for new n2o5 hydrolysis  methods. DO NOT USE EmepReimer,
 ! but one of :; 'SmixTen' , 'Smix', 'Gamma:0.002'

    character(len=20) :: n2o5HydrolysisMethod = 'Smix'
  end type box_useconfig
  type(box_useconfig), public, save :: USES


  ! Emissions variables, also from namelist
  logical, public, save :: use_emis = .False.
  character(len=100), public, save :: emissplit_dir

  ! cloudj photolysis option
  logical, public, save :: use_cloudj = .False.
  logical, public, save :: use_hrlycj = .False.
  character(len=100), public, save :: cloudj_indir
  integer, public, save :: printJ_day = 1000

  ! From config we get e.g. 'nox', 11.4, which then needs to be
  ! transfered to box_emis with the right index 
  type(box_emis_t), dimension(NEMIS_File) :: emis_kgm2day=box_emis_t()

  ! BVOC emissions?
  ! In the box model, a default set of emissions is included in
  ! chem/extra_mechanisms/BoxBVOCemis. These can be modified 
  ! with factors (ZchemData_mod) from config_box.nml:
  ! real, public, save :: fIsop=1.0, fMTL=1.0, fMTP=1.0, fSQT=1.0
  ! The fIsop and fMTL factors will be further modulated by 
  ! light through the SUN variable.
   
  ! ------------------------------------------------

  ! number of lines with initial concentrations. This might have to become
  !  larger
  integer, public, parameter :: NSPEC_INIT = 200 
  type(initBox_t), public, dimension(NSPEC_INIT) :: initBox=initBox_t()

  ! ------------------------------------------------
  ! To allow outout of groups (eg ASOA, all PAN) we copy some ideas from ESX
  ! esx_Variables:
  !> For output species
 ! In namelist we set name and unit: (e.g. 'O3', 'ppb')
  type(OutSpecsList_t), public, dimension(NSPEC_TOT) :: & ! WAS ESX_MAXNOUT
        OutSpecs_list = OutSpecsList_t() ! NOT USED YET in BOX
  type(OutSpecsList_t), public, dimension(29) :: & ! Enough?
        OutGroups_list = OutSpecsList_t()   ! TESTING Nov 2016


   integer, public, save :: &
       nOutGroups = 0  &!> Number of output groups, eg NOx, ASOA
      ,nOutSpecs  = 0   !> Number of output specs, eg NO2, OM25_p
  ! ------------------------------------------------
!A18 for consistency with EMEP code
!=============================================================================
! We have one variable, to say if we are on master-processor
! or not: (kept here to avoid too many dependencies for box-model
! codes which don't need Par_mod

logical, public, save ::  MasterProc = .true.
 ! We allow a flexible string which can switch between different
 ! experiments called by e.g. Solver. A but crude, but
 ! it makes sure the experiments are recorded in the config
 ! system

  character(len=100), save, public :: YieldModifications = 'VBS-T10' ! Default for EmChem16mt



!=============================================================================


  contains

  subroutine Init_Box()

   integer :: i, io_config, iem, ispec
   real :: all_species = 0.    ! inital conc. for all species QUERY WHY NEEDED?
   real :: Ps = 1.01325        ! surface pressure, Pa. Can be 
   real :: tmpr

   !character(*), intent(out) :: dbgVOC
   character(80) :: msg
   character(len=*), parameter:: dtxt = 'ConfigInit:'

   namelist /box_config/ tstart, tend, dt, doy, lat, lon, use_emis, &
      USES, use_cloudj, use_hrlycj, cloudj_indir, printJ_day, & !added for consistency with EMEP
      emissplit_dir, m, h2o, Hmix, emis_kgm2day, &
      fIsop,fMTL,fMTP,fSQT, & ! BVOC factors
      all_species, initBox, &
      debug, OutSpecs_list, OutGroups_list, Temp, rh, debug

   m(:)   = -999.! Will check M and H2O later
   h2o(:) = -999.

   open(newunit=io_config, file='config_box.nml')
   read(io_config, nml=box_config)

  ! We want 3600/dt to be an integer. Also need dt >= 1.0s. Check:
  ! Naive use of modulo on dt fails due to floating-point issues.
  ! We use the integer version, and restrict dt to be > 0.1s
   tmpr = modulo(3600*1000,nint(1000*dt))
   if ( dt < 0.1 .or. tmpr /= 0 ) then
     print *, 'dt too small (<1s) or 3600/dt not integer ', 3600.0/dt,0.001*tmpr
     stop
   end if

   if ( use_emis ) then
     do i = 1, NEMIS_File
        if ( emis_kgm2day(i)%emis  > 0.0 ) then
           iem = find_index( emis_kgm2day(i)%poll, EMIS_File(:) )
           if (iem > 0 ) then
             print *, dtxt//"FOUND EMIS ", emis_kgm2day(i)%poll, iem
             box_emis(iem) = emis_kgm2day(i)%emis
           end if
        end if
     end do
   end if

   ! Decide photolysis method

!MAR2019    if( index( CM_schemes_ChemSpecs, 'CB6') > 0 ) then
!MAR2019       PhotolMethod = 'CB6'
!MAR2019       print *, dtxt//"CB6 Photolysis Triggered"
!MAR2019    end if

   ! Chemical rate coefficients:
    allocate(rct(NCHEMRATES, size(o2)))
    rct(:,:) = 0.

   ! check for M & H2O
   msg = 'ok'  !  checking M & H2O '
   if (M(1) < 0.0 ) then
       M(:) = Ps/(RGAS_J*temp(:) ) * AVOG * 1.0e-6  ! Molec/cm3
       print *, dtxt//"Sets M from Ps: ", Ps, m(1)
   else
       Ps = M(1) * RGAS_J*temp(1) / AVOG * 1.0e6  ! Pa
       print *, dtxt//"Sets Ps from M: ", Ps, m(1)
   end if
   if (h2o(1) < 0.0 ) then
      if ( rh(1) > 0.0 ) then
         h2o(:) = rh2num( rh(:), temp(:) )
         print *, dtxt//"Sets H2O from RH: ", rh(1), h2o(1)
      else
         msg = dtxt//'Error Config - no RH or H2O'
      end if
    else 
       rh(1) = num2rh( h2o(1),temp(1) )
       print *, dtxt//"Sets RH from H2O: ", rh(1), h2o(1)
    end if
    if ( msg /= 'ok' ) stop dtxt//'ERROR in namelist for M, H2O ' ! // msg)

   ! read the initial concentrations
   ! first for all species
   ppb=M(1)*1.e-9

   xChem(:,:) = all_species * ppb

      ! then specific ones
   do i = 1, NSPEC_INIT
     if (trim(initBox(i)%species) == '-') EXIT

     ispec = find_index( initBox(i)%species, species(:)%name )
     if ( ispec > 0 ) then
        xChem(ispec,:) = initBox(i)%init_conc * ppb
     else
        print *, dtxt//'init species "', trim(initBox(i)%species), &
          '" does not exist in this mechanism'
        stop dtxt//'ERROR?' ! or maybe allow
     end if
   end do

   ! set methane and hydrogen fixed background, used by EmChem19cj chem scheme
   do i = 1, NSPEC_INIT
      if (trim(initBox(i)%species) == 'CH4') & 
         methane(:) = initBox(i)%init_conc * ppb
         write(*,*) 'durpdurp', initBox(i)%init_conc
      if (trim(initBox(i)%species) == 'H2') & 
         hydrogen(:) = initBox(i)%init_conc * ppb
         write(*,*) 'durpdurp', initBox(i)%init_conc
   end do

   print "(a,2f8.2,2es12.3)", dtxt//"RH,T,xH2O(x2) ", rh(1), temp(1), rh2num(rh(1),temp(1)), h2o(1)
   o2(:) = 0.2095 * m(:)
   n2(:) = 0.79   * m(:)
   tinv(:) = 1.0 / temp(:)
   itemp(:) = nint(temp(:))

  ! from emep:
   cN2O5(:) = cMolSpeed(temp(:),108.0)
   cHNO3(:) = cMolSpeed(temp(:), 63.0)
   cHO2(:)  = cMolSpeed(temp(:), 33.0)
   cO3(:)   = cMolSpeed(temp(:), 48.0)
   cNO3(:)  = cMolSpeed(temp(:), 62.0)
   cNO2(:)  = cMolSpeed(temp(:), 46.0)

   S_m2m3(:,:) =  1.0e-4  ! crude, 100 cm2/cm3 each AERO%
   S_m2m3(AERO%PM,:) =  6.0e-4  ! crude sum
   aero_fom(:) = 0.5
   aero_fss(:) = 0.1
   aero_fdust(:) = 0.1
   aero_fbc(:) = 0.1

   ! Index of debug species
    debug%ispec = find_index( debug%Spec, species(:)%name )

  end subroutine Init_Box

 end module Config_module
