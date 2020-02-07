!> EmisSplit - module to 'split' emissions, e.g. NOx into NO and NO2
!! or VOC into many compounds.
!! (This module is hacked from emep code **** hence complexity!)
!>------------------------------------------------------------------------   

module EmisSplit_mod

!-------------------------------------------------------------------------   
!  DESCRIPTION:
!  Sets speciation for emissions to be split into inidividual species,
!  e.g. VOC, PM, NOx, per source category for each country from input files. 
!
! Input files:
!    xxxsplit.defaults   (e.g. vocsplit.defaults, pm25split.defaults)
!        where xxx can be e.g. voc, pm25
!
! Output: (module public variable)
!    EmisFrac(ie) type:
!     provides label, nsplit, MW, ...  frac(ie,isec,iland)
!    where ie is e.g. 1 for sox, 2 for nox, ... for voc...
!
!------------------------------------------------------------------------- 


use CheckStop_mod,         only: CheckStop, StopAll
use ChemDims_mod,             only: NSPEC_ADV, NSPEC_TOT, NEMIS_File
use ChemSpecs_mod,            only: species 
use ChemSpecs_mod,            only: define_chemicals ! for self_test 
use NumberConstants,      only: UNDEF_I, UNDEF_R
use PhysicalConstants_mod, only: AVOG
use SmallUtils_mod,        only: wordsplit, find_index, find_indices, trims

implicit none
private

! subroutines:
 public :: RdEmisSplit   ! Reads emisfrac, speciation of voc, pm25, etc.
 public :: doEmisSplit   ! splits e.g. XX kg/m2/h voc to molec/cm2/s C2H4, etc.
 public :: getTestEmisSplit ! Makes use of test SNAP data
 public :: self_test     ! Test of Doxygen

 ! specify EMIS_File = (/ 'sox', etc-.., and NEMIS_File
 ! e.g.  integer, parameter, public :: NEMIS_File = 5
 ! ...     EMIS_File = (/  "sox", "nox", "co ", "voc", "nh3"/)

 ! ========================
 ! ========================
  include 'CM_EmisFile.inc'
 ! ========================
 ! ========================

integer, parameter, public ::  NSECTORS=11 ! SNAP1=power stations, erc
integer, parameter, private :: NLAND = 1   ! simplified from EMEP
integer, private, parameter :: SMAX = 150 ! max number chemicals in split. MCM 3.1 had 136!

type, private :: Esplit_t
  character(len=10)        :: label = '-' ! eg nox, nmvoc
  integer                  :: nsplit = 0
  integer                  :: MolWt  = UNDEF_I
  character(len=len(species(1)%name)), &
         dimension(SMAX) :: spec = '-' ! eg C2H4, NO2
  integer, dimension(SMAX) :: index  = UNDEF_I    ! e.g. index for C2H6, ...
  real,    dimension(SMAX) :: masscorr = UNDEF_R  ! 1/mol.wt
  real,    dimension(SMAX,NSECTORS,NLAND) :: frac  = UNDEF_R
endtype Esplit_t
type(Esplit_t), public, save, dimension(NEMIS_File) :: EmisFrac = Esplit_t()

contains


!-----------------------------------------------------------------------------
!> RdEmisSplit - reads files such as emissplit.defaults.voc to work out
!! the required split of eg NMVOC into individual compounds (C2H4, etc.).
!! The splits are provided for 11 emission source sectors (SNAP) system.
!! Results are stored in EmisFrac derived type
!-----------------------------------------------------------------------------

  subroutine RdEmisSplit(idir,dbg,dbgVOC)

  character(len=*), intent(in) :: idir
  integer, intent(in) :: dbg
  ! For debug we allow a VOC species to receive 100% emissions
  character(len=*), intent(in), optional :: dbgVOC

  !-- local
  integer ::  ie             ! emission index in EMIS_FILE (1..NEMIS_FILE)
  integer ::  itot           ! Index in IX_ arrays
  integer :: io_emis, ios    ! ESX

  !-- for read-ins, dimension for max possible number of columns: 
  !-- for CRI we have 100s of VOC, hence
  character(len=10000) :: txtinput
  character(len=30), dimension(NSPEC_ADV) :: Headers
  character(len=120):: fname, errmsg
  character(len=30):: txt
  character(len=1) :: colon
  real, dimension(SMAX) :: tmp
  real     :: sumtmp
  integer  :: iland,isec,i,n, nn, nHeaders, iland_icode, skipcol
  integer  :: idbgVOC           ! For debug VOC species
  logical  :: defaults          ! Set to true for defaults, false for specials
  character(len=10):: dtxt='EmisSplit:' ! for debug
!-----------------------------------------------

  if ( dbg>0 ) write(*,*) dtxt//'idir ',trim(idir)

  do ie = 1,  NEMIS_file

    EmisFrac(ie)%label=EMIS_File(ie)    ! nox, voc, etc.
    idbgVOC = -1                        ! only used for voc case

    defaults =  .true.  !! ESX  (idef == 0)  EMEP would have idef=2

    !fname = trims( idir//"/emissplit.defaults." // EMIS_FILE(ie) )
    fname = trims( idir//"/emissplit_defaults_" // EMIS_FILE(ie) // '.csv' )
    open(newunit=io_emis,file=fname,iostat=ios)

      call CheckStop( ios/=0 , "EmisGet: ioserror: " // fname )

      if (dbg>0) write(*,*) dtxt//"defaults=", defaults,trim(fname)
 
       !/ Read text line and speciation:
       !  the following lines expect one line of a header text with the
       !  species names, followed by lines of the following format:
       !  iland, isec, tmp1, tmp2.... tmpN+1, where the N+1'th column
       !  is  optional, and for non-reactive species. These non-reactives are
       !  not used in  the rest of the program, but are sometimes needed 
       !  (e.g. VOC) !  to check mass-balance.

    n = 0

    ! emissplit file has format, e.g.:
    !# SOx splits for EmChemXX
    !# Key-word, MASS_ASSUMED= 0 by default, 
    !# change when emissions have artificial mass, e.g. 46 for NOx as NO2
    !: MASS_ASSUMED 64
    !99 99     SO2   SO4     #HEADERS
    !#DATA
    !0   1    95.00   5.00


   READ_DATA: do 

       read(unit=io_emis,fmt="(a)",iostat=ios) txtinput
       if ( ios /=  0 ) exit  READ_DATA     ! End of file
       if( txtinput(1:1) == "#" ) cycle READ_DATA
       if( txtinput(1:1) == ":" ) then  ! Assign MASS_ASSUMED to MolWt:
            read(txtinput,*) colon, txt, EmisFrac(ie)%MolWt
            cycle READ_DATA
       end if

       if ( index(txtinput,'HEADERS') > 0 ) then
            if (dbg>0) write(*,*) dtxt//"headerline "// trim(txtinput)
            !write(*,*) dtxt//"headerline "// trim(txtinput)
            call wordsplit(txtinput,NSPEC_ADV+3,Headers,nHeaders,ios)

           ! Header contains 99 99 #HEADERS and possibly UNREAC
            skipcol = 0
            if( Headers(nHeaders-1) == 'UNREAC' ) skipcol = 1
            EmisFrac(ie)%nsplit = nHeaders - 3 - skipcol

            if(dbg>0) write(*,"(a,4i4)") dtxt//"MWset, NH, nsplit, skipcol ",&
                EmisFrac(ie)%MolWt, nHeaders, EmisFrac(ie)%nsplit, skipcol

            do i = 3, nHeaders-1-skipcol
              nn = i -2
              itot= find_index( Headers(i), species(:)%name, any_case=.true. )
              if(itot > NSPEC_TOT ) then
                 call StopAll( dtxt//"itot>NSPEC_TOT" )
              else if( itot > 0 )  then  ! need to ignore UNREAC -
                 EmisFrac(ie)%index(nn) = itot 
                 EmisFrac(ie)%spec(nn) = species(itot)%name !just helper

                ! Now, get factor needed for scaling mass emissions to molec
                ! EmisFrac(ie)%MolWt non.zero if MASS_ASSUMED set
                 if ( EmisFrac(ie)%MolWt >  0 ) then
                     EmisFrac(ie)%masscorr(nn) = 1.0/EmisFrac(ie)%MolWt
                 else
                     EmisFrac(ie)%masscorr(nn) = 1.0/species(itot)%molwt
                 end if

                 if ( present(dbgVOC) ) then
                    if ( species(itot)%name == dbgVOC ) idbgVOC = nn
                 end if
                 if(dbg>2) write(*,'(a,3i4,2x,a,f8.4)') dtxt//"MassCorr ",&
                    nn, itot, idbgVOC, Headers(i), EmisFrac(ie)%masscorr(nn)

              else 
                 call StopAll ( dtxt//"SPEC not found "//Headers(i) )
              end if
            end do
            cycle READ_DATA
       end if

       n = n + 1

       read(unit=txtinput,fmt=*,iostat=ios)  iland_icode, isec,&
          (tmp(i),i=1, EmisFrac(ie)%nsplit + skipcol )

       if ( idbgVOC > 0 ) then  !! Assign all NMVOC to the debug VOC
         tmp(:) = 0.0
         tmp(idbgVOC) = 100.0  
       end if

       if(iland_icode==0)then
          iland=1!special meaning, always for ESX
       endif

       if (dbg>0 ) write(6,"(a,2i4,999f7.2)") dtxt//"Sec:"//  &
            trim(EmisFrac(ie)%label), isec, EmisFrac(ie)%nsplit, &
            tmp(1:EmisFrac(ie)%nsplit + skipcol)

       !/... some checks:
       sumtmp = sum( tmp(1:EmisFrac(ie)%nsplit+skipcol) )
       if ( ( abs( sumtmp - 100)   >  0.01   )  .or.  &
            ( defaults .and. isec  /= n      )           )   then
           write(unit=errmsg,fmt=*) "ERROR: emisfrac:"//trim(fname)//" ",&
               n, isec, sumtmp
           call CheckStop( errmsg )
       end if

       do nn = 1, EmisFrac(ie)%nsplit 

         !*** assign and convert from percent to fractions: ***

         EmisFrac(ie)%frac(nn,isec,iland) = 0.01 * tmp(nn)

         ! just a check
         !if ( dbg>0 .and. iland == 0 ) then 
         if ( dbg>0 ) then 
            itot = EmisFrac(ie)%index(nn)
            write(*,"(3x,a,4i3,i4,1x,a,f10.4)") dtxt, iland,isec, ie, nn,  &
             itot, trim(species(itot)%name), EmisFrac(ie)%frac(nn,isec,iland)
         endif
      enddo ! i

   enddo READ_DATA 
   close(io_emis)

   if(dbg>0) print *, dtxt//'CHECKS ',  n,  NSECTORS, trim(fname)
   call CheckStop(  defaults .and. n  /=  NSECTORS, &
        "ERROR: EmisGet: defaults .and. n  /=  NSECTORS: " //trim(fname))

  end do ! ie
  ios = 0
   
 end subroutine RdEmisSplit

!-----------------------------------------------------------------------------
!> doEmisSplit - converts sector-level input emissions of sox, nox etc into
!! emissions of individual compounds for an individual SNAP sector (using the
!! EmisFrac split fractions found in RdEmisSplit). 
!! Returns emissions for each species in molecules/cm2/s.
!-----------------------------------------------------------------------------

 subroutine doEmisSplit(emis_snap,units,cc,emis_spec, sec_to_0)

    real, dimension(NEMIS_File,NSECTORS), intent(in) :: &
       emis_snap ! emissions by group and sector, e.g. nox, voc, in kg/m2/h

    character(len=*),  intent(in) :: units ! kgkm2day or ...
    integer,           intent(in) :: cc    ! Country code
    integer, intent(in), optional :: sec_to_0 ! sector that is set to 0
    real, dimension(NSPEC_TOT,NSECTORS), intent(out) :: &
       emis_spec ! emissions by species and sector, molec/cm2/s

    integer :: ie, is, itot
    integer :: snap  ! SNAP emission sector (1-11, e.g. SNAP7=road traffic)

 ! from kg/km2/h to molec/cm2/s, NEEDS MW later!!!!
    real, parameter :: kgkm2h_fac= 1.0e3*AVOG/1.0e10/3600.0
    real :: units_fac
    integer  :: set_sec_to_0 = -999 ! can set all of sector to zero

    if ( units == 'kgkm2day') then
       units_fac=kgkm2h_fac/24.0
    else
       call StopAll('doEmisSplit units not defined: '//units) 
    end if
    if ( present(sec_to_0) ) then
       set_sec_to_0 = sec_to_0
    end if

    emis_spec(:,:) = 0.0
    do ie = 1, NEMIS_File   !  sox, nox, etc
      do is = 1, EmisFrac(ie)%nsplit
         itot = EmisFrac(ie)%index(is)   ! Index in total species array
         do snap = 1, NSECTORS
            if ( snap == set_sec_to_0 ) then
              emis_spec(itot,snap) = 0.
            else
              emis_spec(itot,snap) = emis_snap(ie,snap) &
                !BUG was here  * units_fac/species(itot)%molwt &
                   * units_fac * EmisFrac(ie)%masscorr(is) &
                   * EmisFrac(ie)%frac(is,snap,cc)
            end if
         end do
        !print *, "TMP", itot, species(itot)%name, emis_spec(itot)
      end do ! is
    end do ! ie
 end subroutine doEmisSplit

!-----------------------------------------------------------------------------
!> Provides test set of SNAP data - giving the mass fraction for each
!! pollutant and snap.
!! Output units same as input units, mass-based.
!! Data from emep emissions, mk.UKemisread with year=2010, 
!! EECCA/Modrun15/2015_EMEP_trends/model_input_$year

 subroutine getTestEmisSplit(emis_tot,emis_snap) ! Makes use of test SNAP data
   real, dimension(NEMIS_File), intent(in)           :: emis_tot 
   real, dimension(NEMIS_File,NSECTORS), intent(out) :: emis_snap

   real, dimension(NSECTORS)    :: f  ! Fractions in each snap 
   integer :: ie

   do ie = 1, NEMIS_File
   
     select case (EMIS_File(ie))  ! 'sox', 'nox', etc.
       case('co')    !    2181.72 (=uk total, ktonne)
  f=[ 4.55, 15.63, 14.75,  5.66,  0.51,  0.00, 41.30, 16.53, 1.06,  0.00, 0.0]
       case('nh3')   !     278.98 
  f=[ 0.30,  0.75,  0.11,  1.63,  7.02,  0.43,  3.67,  0.03, 4.59, 81.46, 0.0]
       case('voc')   !     855.20 
  f=[ 0.62,  3.23,  0.50, 15.85, 12.47, 40.40,  6.84,  5.04, 4.00, 11.03, 0.0]
  ! test for 2001
  !f=[ 0.44,  1.84,  0.40, 12.63, 18.89, 28.09, 24.21,  4.14, 2.59,  6.77, 0.0]
       case('nox')   !    1123.00 
  f=[29.99,  6.51, 10.71,  0.18,  0.03,  0.00, 33.24, 19.04, 0.31,  0.01, 0.0]
       case('pm25') !       86.88 
  f=[ 6.33, 29.43,  6.70,  8.29,  3.30,  1.96, 21.63, 14.40, 2.88,  5.09, 0.0]
       case('pmco') !       42.58 
  f=[ 6.23,  0.00,  0.91, 11.56, 12.97,  5.40, 15.93,  1.76, 1.79, 43.46, 0.0]
       case('sox')  !      427.56 
  f=[55.45, 10.86, 22.78,  5.47,  0.04,  0.00,  0.23,  4.89, 0.28,  0.00, 0.0]
       case('ivoc')   !    just test setting ivoc=0.2*voc
  f=[ 0.124,  0.646,  0.1, 3.17, 2.494, 8.08,  1.368,  1.008, 0.80, 2.206, 0.0]
       case('c5h8')   !    
  f=[ 0.00, 0.00, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00, 0.00,  0.00, 100.0]
     case default
       call StopAll('getTestEmisSplit EMIS  NOT found'//EMIS_File(ie)) 
     end select 

     emis_snap(ie,:) = emis_tot(ie) * 0.01 * f ! convert percentages to fractions and 

   end do ! ie

 end subroutine getTestEmisSplit ! Makes use of test SNAP data

!-----------------------------------------------------------------------------
!> self_test subroutine can be called by uncommenting code at bottom
!! demonstrates how this system works.
!-----------------------------------------------------------------------------
 subroutine self_test(idir)
   character(len=*), intent(in) :: idir ! input directory
   real, dimension(NSPEC_TOT,NSECTORS) :: emis_spec = 0.0  ! molec/cm2/s
   real, dimension(NSPEC_TOT,1) :: rcemis = 0.0 ! molec/cm3/s
   real, dimension(NEMIS_File)  :: emis_tot   = 15.0 ! kg/km2/day
   real, dimension(NEMIS_File,NSECTORS) :: emis_snap ! kg/km2/day
   integer :: ie, ispec, itot, ns, snap, k=1, iland=1
   real :: dz = 1.0e5  ! Mixing depth, 1km, in cm

    call define_chemicals()
    call RdEmisSplit(idir,dbg=1,dbgVOC='C3H6')

    print *, "EXAMPLE FOR TRAFFIC (SNAP 7)"
    do ie = 1, NEMIS_File
       print *, "----------"
       ns= EmisFrac(ie)%nsplit
       print "(i3,1x,a,i3,99a8)",ie, EmisFrac(ie)%label, ns, &
         adjustr(EmisFrac(ie)%spec(1:ns)),  "Total"
       print "(i3,1x,a,i3,99f8.3)",ie, EmisFrac(ie)%label, ns, &
         EmisFrac(ie)%frac(1:ns,7,iland), &
         sum( EmisFrac(ie)%frac(1:ns,7,iland) )

    end do ! ie

   ! get some test data, emis_snap:
    call getTestEmisSplit(emis_tot,emis_snap)

    ! Get emissions of each species, now in molec/cm2/s:
    call doEmisSplit(emis_snap,'kgkm2day',iland,emis_spec)

    ! And convert to molec/cm3/s (assuming here k=1):
    rcemis(:,:) = 0.0
    do snap = 1, NSECTORS
      rcemis(:,1) = rcemis(:,1) + emis_spec(:,snap)/dz
    end do
       
    do ie = 1, NEMIS_File
       ns= EmisFrac(ie)%nsplit
       do ispec = 1, ns
          itot=EmisFrac(ie)%index(ispec)
          print "(a20,es12.3)", "RCEMIS "//adjustl(species(itot)%name),&
             rcemis(itot,k)
       end do
    end do ! ie

 end subroutine self_test
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module EmisSplit_mod
!TSTEMX program tstr
!TSTEMX use EmisSplit_mod, only : self_test
!TSTEMX   call self_test('../inputs/emissplit_run')
!TSTEMX end program tstr
