!----------------------------------------------------------------------------!
!> MODULE Debug_mod
!! Stores various debug settings, allowing e.. debug%chem, debug%ddep etc
!! Values cab be boolean, integer (for different levels of output) or strings
!----------------------------------------------------------------------------!
!
module Debug_module
  implicit none
  private


 !> We use DebugCell for particular cell in 1-D or 3-D systems

  logical, public, save :: DebugCell=.true.

  type, private :: debug_t
  !> Debugging flags. Use 0 for none, 1 or 2 for 
    integer :: bics  = 0   !> boundary conditions (internal or extern)
    logical :: Cell  = .false.  !> mainly for EMEP, specific grid cell
    logical :: Chem  = .false.
    integer :: config= 0   !> Extra output
    integer :: cpy   = 0   !> esx_CanopyExchange
    integer :: ddep  = 0   !> Extra output for dry-dep calcs
    integer :: do3se = 0   !> esx_do3se, eg f-factors
    integer :: driver= 0   !> Extra output
    logical :: DryRun = .false.  !> mainly for EMEP, skips chemistry
    integer :: Emis  = 0
    integer :: ispec = 0   !> will be matched against debug_spec
    integer :: Kz    = 0   !> Extra output
    integer :: mafor = 0   !> MAFOR dynamic aerosol code
    integer :: mass  = 0   !> mass budget
    integer :: megan = 0   !> MEGAN canopy model (energy balance, not BVOC)
    integer :: Outputs = 0 !> Used in OutoutConcs_mod
    integer :: Photol = 0  !> Photolysis rates
    logical :: PM    = .false.
    integer :: SOA   = 0   !> Extra output from SOA_mod
    logical :: Solver  = .false.  !> for numerical soln
    character(len=20) :: Spec = 'NOTSET'  !> Extra output, e.g for O3
    character(len=30) :: VOC   = '-'  !> for tsting emisspluit
    integer :: wdep  = 0   !> Extra output for wet-dep calcs
    integer :: Zchem = 0   !> Extra output from ChemSolver
    integer :: Zdata = 0   !> Extra output
    integer :: Zdiff = 0   !> Extra output from KdiffSolver
    integer :: Zemis = 0   !> Extra output from KdiffSolver
    integer :: Zmet  = 0   !> Extra output (ESX and EMEP only)
    integer :: Zveg  = 0   !> Extra output from 1-D , e.g. PAR
    character(len=20)     :: datetxt = '-'       ! default.
  end type debug_t
  type(debug_t), public, save :: debug = debug_t()

end module Debug_module
