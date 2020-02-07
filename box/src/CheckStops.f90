!_____________________________________________________________________________
!> In a CTM this module could include a wide range of checks, and would
!! stop all processors. (Just calling fortran "stop" is not recommended
!! for MPI systems.)

! Can call with a logical expression, e.g.
!
!   call CheckStop( t>t1, "Exceeded time")
!   call CheckStop( errmsg /= "ok", "Error: "// errmsg)
!


module CheckStop_mod
  use NumberConstants, only : UNDEF_R
  implicit none
  private
 
  public :: StopAll
  public :: CheckStop
  public :: CheckStop_ok, CheckStop_TF
  public :: checkValid

  interface CheckStop
     module procedure CheckStop_ok
     module procedure CheckStop_TF
  end interface CheckStop

contains

 !----------------------------------------------------------------------------!
  subroutine StopAll(errmsg)
      character(len=*), intent(in) :: errmsg

      if ( errmsg /= "ok" ) then
        print *,   "ERROR: ", trim(errmsg)
        !call MPI_ABORT(MPI_COMM_WORLD,9,INFO)
        stop !"ERROR: "// trim(errmsg)
      end if
  end subroutine StopAll

 !----------------------------------------------------------------------------!
 !---- Two  variations on CheckStop:

 subroutine CheckStop_ok(errmsg)    ! Test if errmsg /= "ok"
      character(len=*), intent(in) :: errmsg

      if ( errmsg /= "ok" ) call StopAll(errmsg)

 end subroutine CheckStop_ok
 !-------------------------------------------------
 subroutine CheckStop_TF(is_error, txt)
   logical, intent(in) :: is_error
   character(len=*), intent(in) :: txt

   if( is_error ) call StopAll(txt)

 end subroutine CheckStop_TF

 !----------------------------------------------------------------------------!
 !> SUBROUTINE checkValid compares for UNDEF and also for 
 !! NaN (using +0 trick). Stops code if there is a problem.

  subroutine checkValid( x, txt )
    real, intent(in) :: x
    character(len=*), intent(in) :: txt

    call CheckStop( x == UNDEF_R,  "checkValid UNDEF: "//txt )
    call CheckStop( x /= x+0, "checkValid NaN: "//txt )

  end subroutine checkValid

end module CheckStop_mod
!_____________________________________________________________________________
