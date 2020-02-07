!> MODULE OutputConcs_mod.f90 - A component of the EMEP/ESX Model System>
!!-----------------------------------------------------------------------------
!! Handles different species, groups, units where needed for output.
!! Tested so far as part of BoxChem system, but planned for ESX too.
!!-----------------------------------------------------------------------------
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!-----------------------------------------------------------------------------
module OutputConcs_mod
  use ChemDims_mod,        only: NSPEC_TOT
  use ChemSpecs_mod,       only: species, S1 => FIRST_SEMIVOL, S2=> LAST_SEMIVOL
  use ChemGroups_mod,      only: chemgroups
  use Debug_module,        only: debug  ! for %ispec
  use NumberConstants,     only : UNDEF_I, UNDEF_R
  use PhysicalConstants_mod, only: AVOG
  use SmallUtils_mod,      only: LenArray, find_index

  implicit none
  private

  public :: InitOutputConcDefs  ! finds indices, groups, units
  public :: getOutputConcs  ! gets concentrations

  integer, private, parameter :: TXTLEN=len(species(1)%name)

  !> To list wanted output species/groups:
  type, public :: OutSpecsList_t
      character(len=TXTLEN) :: name = '-' ! eg SO2
      character(len=TXTLEN) :: unit = '-' ! eg ppb
  end type OutSpecsList_t

  !> For output species
  integer, public, save :: nOutputConcs = UNDEF_I
  !integer, private, parameter :: NSPECSMAX=size(esx%OutSpecs_list(:))
  integer, private, parameter :: NMAX=1500  ! specs within groups, MCM has eg 1200 RO2!
  type, public :: OutConcs_t
      character(len=TXTLEN) :: name = '-' ! eg SO2
      character(len=TXTLEN) :: unit = '-' ! eg ppb
      integer :: nspecs   = UNDEF_I       ! number > 1 for groups
      integer :: ind      = UNDEF_I       ! index in xChem, 0 if group
      integer :: ionum    = UNDEF_I       ! io number for file
      integer,dimension(NMAX) :: specs = UNDEF_I ! indices in xChem
      logical :: group    = .false.    ! if air conc scaling needed
      logical :: timesAir = .true.     ! if air conc scaling needed, true for ppb
      logical :: semivol  = .false.    !  eg for SOA specs
  end type OutConcs_t
  ! Dimension for all species, and for now 29 groups
  type(OutConcs_t), public, dimension(NSPEC_TOT+29) :: OutConcDefs=OutConcs_t()

  logical, private, save :: dbg = .false.

 contains
 !----------------------------------------------------------------------------
  subroutine InitOutputConcDefs(wanted_specs,wanted_groups)
    type(OutSpecsList_t), dimension(:), intent(inout) :: wanted_specs
    type(OutSpecsList_t), dimension(:), intent(in)    :: wanted_groups
    integer :: ic, iS, ispec, n, igrp
    integer, save :: nspecs, ngroups
    character(len=*), parameter :: dtxt = 'getOutConcDefs:'
    integer :: gPM10

      dbg = debug%Outputs > 0
      gPM10=find_index("PM10",chemgroups%name)

 !     print *, dtxt//"NSPECS, GROUPS, gPM10 wanted ", nspecs, ngroups, gPM10

   !/ Species

     !/ special case - all
     if ( trim(wanted_specs(1)%name) == 'all' ) then

        print *, dtxt//'all species are output'
        do iS = 1, NSPEC_TOT
           wanted_specs(iS)%name = species(iS)%name
           wanted_specs(iS)%unit = 'ppb'
           if ( gPM10 < 1 ) CYCLE
           if ( find_index(iS, chemgroups(gPM10)%specs ) > 0 ) then
              wanted_specs(iS)%unit = 'ugm3'
           end if 
           if(dbg) print *, dtxt//'All->',iS,wanted_specs(iS)
        end do
     end if

     !/ now... process
     nspecs  = LenArray(wanted_specs(:)%name, "-") 
     n  = 0   ! counts valid species/groups
     do ic = 1, nspecs
       ispec = find_index( wanted_specs(ic)%name, species(:)%name )
       if ( ispec > 0 ) then
         n = n + 1
         OutConcDefs(n)%name = wanted_specs(ic)%name
         OutConcDefs(n)%unit = wanted_specs(ic)%unit
         OutConcDefs(n)%nspecs    = 1
         OutConcDefs(n)%specs(1)  = ispec
         if ( ispec >= S1 .and. ispec <= S2 ) then
            OutConcDefs(n)%semivol = .true.
         end if
         if(dbg) print *, dtxt//"PROCESS: Spec:"//&
                     wanted_specs(ic)%name, n, ispec, OutConcDefs(n)%semivol
       else
         print *, dtxt//"WARNING: OutSpec not in CM:"//&
           wanted_specs(ic)%name
       end if
     end do !ic

     print *, dtxt//"PROCESSED SPECS: ", n

      ngroups = LenArray(wanted_groups(:)%name, "-") 
      do ic = 1, ngroups
        igrp = find_index( wanted_groups(ic)%name, chemgroups(:)%name )

        if( igrp < 1  ) then
          print *, dtxt//"WARNING: OutGroup not in CM:"//&
                     wanted_groups(ic)%name
          cycle
        end if

        n = n + 1
        OutConcDefs(n)%group = .true.
        OutConcDefs(n)%name =trim( wanted_groups(ic)%name)//'s_SUM'
        OutConcDefs(n)%unit = wanted_groups(ic)%unit

        print *, dtxt//"PROCESS: Group:"//&
                     wanted_groups(ic)%name
        do iS = 1, size(chemgroups(igrp)%specs)
            ispec = chemgroups(igrp)%specs(iS)
            OutConcDefs(n)%specs(iS) = ispec
            OutConcDefs(n)%nspecs    = iS ! get last 
             if ( ispec >= S1 .and. ispec <= S2 ) then
                OutConcDefs(n)%semivol = .true.
                print *, dtxt//"SEMI:", wanted_groups(ic)%name, species(ispec)%name
             end if
        end do
     end do ! ic groups
     nOutputConcs =  n   ! count of valid species/groups
     print *, dtxt//"PROCESSED GROUPS: ", n

  end subroutine InitOutputConcDefs
 !----------------------------------------------------------------------------
  subroutine getOutputConcs(xn,Fpart,airM,xnout)
   real, dimension(:,:), intent(in) :: xn ! xn(specs,k)
   real, dimension(:,:), intent(in) :: Fpart ! fraction as particle
   real, dimension(:), intent(in)   :: airM !  air conc, m(k)
   real, dimension(:,:), intent(out) :: xnout ! xn(O3,k)
   real, dimension(size(airM)) :: outCz
   integer :: n, iS, ispec, nz, iz
   real :: tmpfac
   character(len=*), parameter :: dtxt = 'getOutputConc:'

    nz = size(airM) ! 1 for box, more for ESX

    do n = 1, nOutputConcs
      outCz(:) = 0.0

      do iS = 1, OutConcDefs(n)%nspecs
           ispec = OutConcDefs(n)%specs(iS)
           do iz = 1, nz
             tmpfac = 1.0
             if ( ispec >= S1 .and. ispec <= S2 ) then
                 tmpfac = Fpart(ispec,1)
             end if
            !-- units scaling
             select case ( OutConcDefs(n)%unit )
                 case ( "ppb" )
                   tmpfac = tmpfac * 1.0e9/airM(iz)  ! scale with air density
                 case ( "ppt" )
                   tmpfac = tmpfac * 1.0e12/airM(iz)  ! scale with air density
                 case ( "ugm3" ) ! from #/cm3 to ug/m3
                   tmpfac = tmpfac * species(ispec)%molwt/AVOG * 1.0e12
                 !case ( "molec_per_cm3", "pn_per_cm3" ) = 1.0
             end select
              !--
             outCz(iz) = outCz(iz) + tmpfac*xn(ispec,iz)
             if(ispec==debug%ispec) print "(a,3i4,2a12,2g12.3)", dtxt//"OUTCZ ", &
               n, iS, ispec, trim(OutConcDefs(n)%name), &
               trim(species(ispec)%name) , tmpfac, outCz(iz)
           end do ! iz
       end do ! iS
       xnout(n,:) = outCz(:)
       !print *, dtxt//"TMPOUT ", n, OutConcDefs(n)%name, &
       ! OutConcDefs(n)%nspecs, outCz(1), OutConcDefs(n)%unit
   end do !n
  end subroutine getOutputConcs
!=============================================================================
end module OutputConcs_mod
!TSTEMX program testr
!TSTEMX use ChemDims_mod,      only : n => NSPEC_TOT
!TSTEMX use ChemSpecs_mod,     only : define_chemicals
!TSTEMX use ChemGroups_mod, only : Init_ChemGroups
!TSTEMX use OutputConcs_mod
!TSTEMX implicit none
!TSTEMX type(OutSpecsList_t), dimension(3) :: &
!TSTEMX    wanted_specs=OutSpecsList_t(), wanted_groups=OutSpecsList_t()
!TSTEMX integer :: i
!TSTEMX real, dimension(n,2) :: xn, Fpart = 0.0
!TSTEMX real, dimension(2) :: m = 2.55e19
!TSTEMX do i = 1,n ! i is species no, so very fake here!
!TSTEMX   xn(i,:) = i*1.0
!TSTEMX end do
!TSTEMX wanted_groups(1)%name = 'OXN'; wanted_groups(1)%unit = 'ugm3'
!TSTEMX wanted_groups(2)%name = 'PAN'; wanted_groups(2)%unit = 'ppb'
!TSTEMX wanted_groups(3)%name = 'OXS'; wanted_groups(3)%unit = 'ppb'
!TSTEMX wanted_groups(1)%name = 'SO2'; wanted_groups(1)%unit = 'ppb'
!TSTEMX call define_chemicals()
!TSTEMX call Init_ChemGroups()
!TSTEMX call InitOutputConcDefs(wanted_specs,wanted_groups)
!TSTEMX call getOutputConcs(xn,Fpart,m)
!TSTEMX   print *, len(OutConcDefs(1)%name), nOutputConcs
!TSTEMX end program testr
