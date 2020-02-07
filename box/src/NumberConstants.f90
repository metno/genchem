module NumberConstants
  implicit none
  private

  public :: is_num  ! simple function to check if string is number
  public :: is_int  ! simple function to check if string is integer

! from KPP system
! Definition of different levels of accuracy
! for REAL variables using KIND parameterization
!
! KPP SP - Single precision kind
  integer, public, parameter :: sp = selected_real_kind(6,30)
! KPP DP - Double precision kind
  integer, public, parameter :: dp = selected_real_kind(14,300)
! CYGWIN can't handle quad precision, so we re-define
! KPP QP - Quadruple precision kind
  integer, public, parameter :: qp = dp ! selected_real_kind(18,400)

!DEWS working precision can be changed here
! Typically should be dp, but qp for testing

  integer, public, parameter :: wp = qp

!integer, public, parameter :: dp = kind(0.0d0)  ! Double precision real(qp) kind

! Sentinel values
!real(qp), public, parameter :: UNDEF_D = -huge(0.0_dp)
    real,     public, parameter :: UNDEF_R = -huge(0.0)
    real(wp), public, parameter :: UNDEF_WP = -huge(0.0_wp)
    integer,  public, parameter :: UNDEF_I = -huge(0)

contains
  function is_num(txt) result (is_n)
    character(len=*), intent(in) :: txt
    integer :: i, ichar0, ichar9  ! range of 0..9
    logical :: is_n
    character(len=1) :: c1
    ichar0 = ichar('0')
    ichar9 = ichar('9')
    c1 = adjustl(txt)  ! first character

   !  do i = 40, 59 !    print *, "HCRA ", i, char(i) !  end do

     i = ichar( c1 )

    ! test for start:allowed +, -, ., number
     is_n = ( scan(c1,'.-+') > 0 )  .or. ( i >= ichar0 .and. i <= ichar9 ) 

  end function is_num

  ! is_int allows e.g 100, 1e3, etc. Only tests for decimal
  function is_int(txt) result (is_i)
    character(len=*), intent(in) :: txt
    logical :: is_i
    is_i = is_num(txt) .and. index(txt,'.') == 0
  end function is_int

end module NumberConstants
!TSTEMX program testr
!TSTEMX   use NumberConstants, only : is_num, is_int
!TSTEMX   print *, "Checks "
!TSTEMX   print *, is_num('no'), is_int('123'), is_num('+12.3'), is_num(' 12.3 ')
!TSTEMX end program testr
