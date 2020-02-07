! <OwnDataTypes_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 
module OwnDataTypes_mod
use NumberConstants, only : UNDEF_I, UNDEF_R
implicit none
private

  ! depmap
  ! gtype for species groups, used in CM_ChemSpecs and Derived

  public :: print_Deriv_type
  public :: print_Asc2D
  integer, public, parameter :: &
    TXTLEN_DERIV = 34, &
    TXTLEN_SHORT = 24, &
    TXTLEN_IND   =  6

  ! Contains some user-defined data-types, and routine associated
  ! with these. Collecting them here will
  ! avoid some dependencies, and shorten some My type modules.
  !
  ! depmap used in DryDep_mod and My_WetDep_ml
  ! Deriv used in My_Derived and  Derived_mod
  ! VBST from SOA_mod

  !/-- we define a type to map indices of species to be deposited
  !   to the lesser number of species where Vg is calculated

  type, public :: depmap
    integer :: ind   ! Index of species in IXADV_ or IX_ arrays
    character(len=TXTLEN_SHORT) :: surrogate  ! Surrogate species in  calculated dep arrays
    real    :: vg    ! if CDDEP_SET, give vg in m/s
  endtype depmap

  !==================
  !/ generic groups for integers
  type, public :: typ_i2
    integer :: int1
    integer :: int2
  endtype typ_i2

  type, public :: typ_i3
    integer :: int1
    integer :: int2
    integer :: int3
  endtype typ_i3
  
  !/ generic group for two (short) strings
  type, public :: typ_ss
    character(len=TXTLEN_SHORT) :: txt1='-' ! e.g. POD1_IAM_DF
    character(len=TXTLEN_SHORT) :: txt2='-' ! e.g. POD1_IAM_DF
  endtype typ_ss

  !/ generic group for name and pointer to arrays
  type, public :: typ_sp
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: specs  ! was ptr
  endtype typ_sp

  !/ HI: generic group for name and two pointers to one integer and one
  !/  real array
  type, public :: typ_factors
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: species ! like ptr in typ_sp
    real, dimension(:), pointer :: factors 
  endtype typ_factors

  !/ HI: generic group for name and two pointers to one integer and one
  !/  character array
  type, public :: typ_maps
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: species ! like ptr in typ_sp
    character(len=TXTLEN_SHORT), dimension(:), pointer :: maps ! other 
      ! species' or variables' names to map this one to
  endtype typ_maps

  !/ generic group one (short) string & one integer
  type, public :: typ_si
    character(len=TXTLEN_SHORT) :: name = '-'
    integer :: ind = UNDEF_I
  endtype typ_si

  !/ generic group one (short) string & one real
  type, public :: typ_sr
    character(len=TXTLEN_SHORT) :: name
    real :: rval
  end type typ_sr

  !/ generic group for three (short) strings
  type, public :: typ_s3
    character(len=TXTLEN_SHORT) :: txt1,txt2,txt3
  endtype typ_s3

  !/ generic group for four (short) strings
  type, public :: typ_s4
    character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4 ! e.g. POD1_IAM_DF
  endtype typ_s4

  !/ generic group for five (short) strings
  type, public :: typ_s5i
    character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4, &
                                   txt5 ! e.g. SO2,ugS,2d,"AIR_CONC",M
    integer                     :: ind  ! e.g. IOU_DAY
  endtype typ_s5i

  !==================
  !+ Derived output type
  type, public:: Deriv
    character(len=TXTLEN_DERIV) :: name     ! e.g. DDEP_SO2_m2Conif
    character(len=TXTLEN_SHORT) :: class    ! Type of data, e.g. ADV or Mosaic
    character(len=TXTLEN_SHORT) :: subclass !  e.g. "VG", "Rns"
    character(len=TXTLEN_SHORT) :: txt      ! text where needed, e.g. "Conif"
    character(len=TXTLEN_SHORT) :: unit     ! writen in netCDF output
    integer :: index          ! index in concentation array, or other
    integer :: f2d            ! index in f_2d arrays
    logical :: dt_scale       ! used only if we need a factor on dt_advec,
    real    :: scale          !  e.g. use 100.0 to get cm/s
    logical :: avg            ! True => average data (divide by nav at end),
                              !  else accumulate over run period
    integer :: iotype         ! sets output timing
  endtype

 ! Sentinel values (moved to NumberConstants)
 ! real,    private, parameter :: UNDEF_R = -huge(0.0)
 ! integer, private, parameter :: UNDEF_I = -huge(0)

  !==================
  !+ Hourly ASCII/NetCDF output type
  type, public:: Asc2D
    character(len=TXTLEN_DERIV):: name = "NOTSET"   ! Name (no spaces!)
    character(len=TXTLEN_SHORT):: type = "NOTSET"  ! "ADVppbv" or "ADVugm3" or "SHLmcm3"
!   character(len=9) :: ofmt      ! Output format (e.g. es12.4)
    integer          :: spec = UNDEF_I   ! Species number in xn_adv or xn_shl array
                                  ! or other arrays
    integer          :: ix1 = UNDEF_I,ix2 = UNDEF_I   ! bottom-left,upper-right x
    integer          :: iy1 = UNDEF_I,iy2 = UNDEF_I   ! bottom-left,upper-right y
    integer          :: nk = UNDEF_I     ! number of vertical levels
    character(len=TXTLEN_SHORT) :: unit   ! Unit used
    real             :: unitconv = UNDEF_R  !  conv. factor
    real             :: max      = UNDEF_R        ! Max allowed value for output
  endtype

  !==================
  !+ Defines SOA, NONVOL and VBS params
  type, public :: VBST
    integer     :: index    ! just for clarity
    real        :: CiStar   ! ug/m3
   !real        :: Tref     ! Assumed 300
    real        :: DeltaH   ! kJ/mole
  endtype VBST

contains
!=========================================================================
subroutine print_Asc2D(w)
  type(Asc2D), intent(in) :: w  ! wanted
  write(*,*) "Prints Asc2D type ========================="
  write(*,"(a,a)")      "Name   :", trim(w%name)
  write(*,"(a,a)")      "type  :", trim(w%type)
  write(*,"(a,i4)")     "spec  :", w%spec
  write(*,"(a,4i5)")    "ix1,ix2,iy1,iy2  :", w%ix1, w%ix2, w%iy1,w%iy2
  write(*,"(a,a)")      "unit  :", trim(w%unit)
  write(*,"(a,i4)")     "nk  :", w%nk
  write(*,"(a,es10.3)") "unitconv  :", w%unitconv
  write(*,"(a,es10.3)") "max  :", w%max
end subroutine print_Asc2D
!=========================================================================
subroutine print_Deriv_type(w)
  type(Deriv), intent(in) :: w  ! wanted
  write(*,*) "Prints Deriv type ========================="
  write(*,"(a,a)")      "Name   :", trim(w%name)
  write(*,"(a,a)")      "class  :", trim(w%class)
  write(*,"(a,a)")    "subclass :", trim(w%subclass)
  write(*,"(a,a)")      "txt    :", trim(w%txt)
  write(*,"(a,a)")      "units  :", trim(w%unit)
  write(*,"(a,i3)")     "index  :", w%index
  write(*,"(a,i3)")     "f2d    :", w%f2d
  write(*,"(a,a10)")    "txt    :", w%txt
  write(*,"(a,es10.3)") "scale  :", w%scale
  write(*,*)          "dt_scale :", w%dt_scale
  write(*,*)               "avg :", w%avg
endsubroutine print_Deriv_type
!=========================================================================
endmodule OwnDataTypes_mod
