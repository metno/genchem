! <boxChem_ChemSolver.f90 - derived from EMEP MSC-W  Chemical transport Model>
!-----------------------------------------------------------------------------!

  module ChemSolver

!-----------------------------------------------------------------------------!
  !=======================================================================!
  ! The following chemical solver uses variable chemical timesteps and
  ! is based on the scheme suggested in J.G. Verwer and D. Simpson (1995)
  ! "Explicit methods for stiff ODEs from atmospheric chemistry",
  ! Aplied Numerical Mathematics 18 (1995) 413.
  !
  ! Note that the exact formula used have been re-arranged for greater
  ! efficiency (Steffen Unger).
  ! Variable Dchem is used to keep track of changes from call to call.
  ! Note: decoupling of (NO3,N2O5), (PAN,CH3COO2), (MPAN,MACRO2)
  ! variable timestep (Peter Wind)
  !=======================================================================!

    use AeroConstants_mod,   only: AERO
   !Aqueous is just to allow compilation when using emep_setup
    use Aqueous_mod,         only: ICLOHSO2, ICLRC1, ICLRC2, ICLRC3, aqrck
    use BioNatConstants_mod, only: NATBIO ! for consistency with EMEP/ESX
    use CheckStop_mod,       only: CheckStop
    use ChemFields_mod, only : x, xold, xnew, cell_tinv 
    use ChemFunctions_mod         ! ESX New
    use ChemGroups_mod            ! => RO2POOL
    use ChemDims_mod              ! => NSPEC_TOT, O3, NO2, etc.
    use ChemSpecs_mod             ! => NSPEC_TOT, O3, NO2, etc.
    use Config_module,       only: YieldModifications
    use Debug_module,        only: debug
    use DefPhotolysis_mod         ! => IDHNO3, etc.
    use ModelConstants,      only: MasterProc
    use SmallUtils_mod,      only: num2str
    use YieldModifications_mod    ! eg YA_ for SOA aerosol.
    use ZchemData_mod             !rct, rcbio, rcemis for CM_Reactions

  implicit none

  private
  public  :: chemsolve        ! Runs chemical solver
  ! Just some helper routines added for boxChem
!  private :: PrintLog
  private :: checkNans

  integer, parameter:: &
      nchemMAX=1200    & ! ESX 15
     ,NUM_INITCHEM=5   & ! Number of initial time-steps with shorter dt
     ,EXTRA_ITER = 1     ! Set > 1 for even more iteration
  real, parameter:: ZERO = 1.0e-30   !  Had _dp option in some versions
  real, save     :: dt_initchem=1.0 ! shorter dt for initial time-steps
  real, public :: dbgLoss             ! for output in debugging
 

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! THIS routine can be used to print to both standard output and a log
! file. This is done in the emep, but 
! subroutine PrintLog(txt,OutputProc,ioOption)
!   character(len=*), intent(in) :: txt
!   logical, intent(in), optional :: OutputProc  !typically MasterProc, me==0
!   integer, intent(in), optional :: ioOption    !use for other files
!   logical :: ok2print
!   integer :: io
!   ok2print = .true.
!   if ( present(OutputProc) ) ok2print = OutputProc
!   if ( ok2print) then
!     io = 6                  ! default
!     if ( present(ioOption) ) io = ioOption
!     write(*,*)  trim(txt)
!     write(io,*)  trim(txt)
!   end if
! end subroutine PrintLog
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine checkNans(x,specs,txt)
    real,             dimension(:), intent(in) :: x
    character(len=*), dimension(:), intent(in) :: specs
    character(len=*) :: txt
    integer :: n
    do  n = 1, size(x)
      if ( isnan(x(n)) ) then
!      if ( x(n) /= x(n) )  then
        print *, "OOOPS!! ChemSolver:Nan found!"//txt, n, specs(n)
        stop 'OOOPS'//txt
      end if  
    end do
 end subroutine checkNans
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  subroutine chemsolve( k, dtstep, xn,Dchem,debug_level,itern)

    !.. In
    integer, intent(in) :: k   ! Needed for rates. ESX change?
    real, intent(in) :: dtstep  ! external time-step, e.g. 20 mins for dt_advec
    real, dimension(:), intent(inout) ::  xn    ! Nz, Concs., molec./cm3
    real, dimension(:), intent(inout) ::  Dchem ! Nz, Tendencies, ..'
                                                ! = d xn/dt due to chem
    integer, intent(in) :: debug_level          ! 0 for none, ..
    integer, intent(in), optional :: itern
    integer :: toiter

    logical, save ::  first_call = .true.

    real, parameter ::  CPINIT = ZERO ! 1.0e-30  ! small value for init

    !  Local
    integer :: ichem, iter,n        ! Loop indices
    integer, save ::  nchem, nzlev    ! No chem time-steps, and z-levels
    real    ::  dt2
    real    ::  P, L                ! Production, loss terms
    real    :: xextrapol   !help variable

    ! Concentrations : xold=old, x=current, xnew=predicted
    ! - dimensioned to have same size as "x"

    !real, dimension(size(xn)) :: x, xold ,xnew   ! Working array [molecules/cm3]
    real, dimension(nchemMAX), save :: dti       &! variable timestep*(c+1)/(c+2)
                           ,coeff1,coeff2,cc ! coefficients for variable timestep

!======================================================

    if ( first_call ) then
       nzlev = size( rct, dim=2 )
       print *, "ChemSolver sizes spec x nz ", size(xn), nzlev 
       !print *, "ESX RATES ", rct(1)
       Dchem=ZERO
       call makedt(dtstep,dti,nchem,coeff1,coeff2,cc)
       if ( MasterProc ) then
!         call PrintLog ('Chem dts: nchemMAX: ' //num2str(nchemMAX))
!         call PrintLog ('Chem dts: nchem: ' //num2str(nchem))
!         call PrintLog ('Chem dts: NUM_INITCHEM: ' //num2str(NUM_INITCHEM))
!         call PrintLog ('Chem dts: dt_initchem: ' //num2str( dt_initchem))
!         call PrintLog ('Chem dts: EXTRA_ITER: ' //num2str( EXTRA_ITER))
!         call PrintLog ('Chem Yields: ' //trim( YieldModifications ))
         if(debug%DryRun) write(*,*) "debug%DryRun Solver"
       end if
    endif

!======================================================


    !**  toiter gives the number of iterations used in TWOSTEP.
    !**  Use more iterations near ground:

     toiter = 3
     if ( present(itern) ) then
       if( itern > 1) then
           toiter = toiter * itern
       end if
     end if

   ! to get better accuracy if wanted (at CPU cost)
    
   !ESX toiter = toiter * EXTRA_ITER


    !** Establishment of initial conditions:
    !   Previous concentrations are estimated by the current
    !   minus Dchem because the current may be changed by
    !   processes outside the chemistry:

     xnew(:) = xn(:)

     call checkNans(xnew,species(:)%name,'posA') 

     x(:)    = xn(:) - Dchem(:)*dti(1)*1.5  ! QUERY?

     call checkNans(xnew,species(:)%name,'posB') 

     x(:)    = max (x(:), ZERO )


       !*************************************
       !     Start of integration loop      *
       !*************************************
      if ( first_call .or. YieldModificationsInUse ) then
         cell_tinv = tinv(k)  ! 1/temp for this cell
           call doYieldModifications('init')
      end if


      do ichem = 1, nchem

          do n=1,NSPEC_TOT

             xextrapol = xnew(n) + (xnew(n)-x(n)) *cc(ichem)
             xold(n) = coeff1(ichem)*xnew(n) - coeff2(ichem)*x(n)

             xold(n) = max( xold(n), ZERO )
             x(n) = xnew(n)
             xnew(n) = xextrapol

             call checkNans(xnew,species(:)%name,'posC') 


          enddo

          dt2  =  dti(ichem) !*(1.0+cc(ichem))/(1.0+2.0*cc(ichem))

          where ( xnew(:) < CPINIT  )
             xnew(:) = CPINIT
          end where

!== Here comes all chemical reactions
!=============================================================================
          if ( debug%DryRun ) then
            ! Skip fast chemistry
          else

            do iter = 1, toiter  !ESX (k)
!
! The chemistry is iterated several times, more close to the ground than aloft.
! For some reason, it proved faster for some compilers to include files as given below
! with the if statements, than to use loops.
!Just add some comments:
!At present the "difference" between My_FastReactions and My_SlowReactions
!is that in My_Reactions the products do not reacts chemically at all,
!and therefore do not need to be iterated.  We could have another class
!"slowreactions", which is not iterated or fewer times. This needs some
!work to draw a proper line ......


                   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                   include 'CM_Reactions1.inc'
                   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                   ! print *, 'N2(k) = ', N2(k)
                   ! print *, 'N2 = ', N2
                   ! print *, 'O2(k) = ', O2(k)
                   ! print *, 'O2 = ', O2
                   ! print *, 'M(k) = ', M(k)
                   ! print *, 'M = ', M
                   ! stop

                ! VBS VBS VBS VBS 

            end do !! End iterations

           !YIELDs  Allows change of gas/aerosol yield, which currently is
           ! only used for SOA species to be handled in CM_Reactions2

            if ( YieldModificationsInUse ) then
                !OLD   1/cell_tinv, iter, toiter(k)
                !OLD if( iter == toiter ) call doYieldModifications('run')
                call doYieldModifications('run')
           end if



          ! slower? species

          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
           include 'CM_Reactions2.inc'
          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

          end if ! debug%DryRun

          !if( first_call) print "(a,2es10.3)", "ChemSolver RCEMIS,BIO ", &
          !    rcemis(C5H8,10), rcbio(C5H8,10)

          if ( debug_level > 0 ) then !LOTS of output, careful. Useful
             do n=1,NSPEC_TOT
               ! write(*, "(a,2i3,1x,a,10es10.2)") "ESX TESTING CHECM ", ichem,&
                  ! n, species(n)%name, xnew(n)
               if(iter==1) write(*, "(a,2i3,1x,a,10es10.2)")"TESTCHEM ", &
                  ichem, n, species(n)%name, xnew(n)
               call checkNans(xnew,species(:)%name,'posD') 
             end do !n
          end if
       end do ! ichem

       !*************************************
       !     End of integration loop        *
       !*************************************


       !**  Saves tendencies Dchem and returns the new concentrations:

            Dchem(:) = (xnew(:) - xn(:))/dtstep

            xn(:) = xnew(:)

       first_call = .false.
  end subroutine chemsolve

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine  makedt(dtstep, dti,nchem,coeff1,coeff2,cc)

 !=====================================================================
 ! Makes coefficients for two-step (written by Peter Wind, Febr. 2003)
 ! The formulas for coeff1, coeff2 and dti can be found in:
 ! J.G. Verwer and D. Simpson, "Explicit methods for stiff ODEs from
 !  atmospheric chemistry", Aplied Numerical Mathematics 18 (1995) 413
 !
 ! Note: It is better to take first some small steps, and then
 !       larger steps, than to increase the timestep gradually.
 !=====================================================================

 implicit none

 real, intent(in) :: dtstep
 real, dimension(nchemMAX),intent(out) :: dti,coeff1,coeff2,cc
 integer,                  intent(out) :: nchem

 real    :: ttot, dt(nchemMAX), inp_dtstep
 real :: dt_init   ! time (seconds) with initially short time-steps
 integer :: i
 logical, save :: first_dt = .true.
 character(len=*), parameter :: dtxt = 'makedt:'
!_________________________

  nchem=nchemMax      ! number of chemical timesteps inside dt_advec
  inp_dtstep = dtstep ! saved for printout only

  if ( 2*NUM_INITCHEM*dt_initchem >= dtstep) then ! shorter dt_initchem needed
    dt_initchem = dtstep / ( 2*NUM_INITCHEM )
    if ( debug%Chem) write(*,'(a,4f12.3)') dtxt//"reduce dt_initchem",&
                                      inp_dtstep, dtstep, dt_initchem
  end if

   dt_init = NUM_INITCHEM*dt_initchem 

!/ ** For smaller scales, but not tested
  !emep: if(dtstep<620.0)

   nchem = NUM_INITCHEM +int((dtstep- dt_init) / dt_init )

   if (debug%Chem)write(*,'(a,4f12.3,2i3)')dtxt//"orig tstep initchem  init",& 
      inp_dtstep, dtstep, dt_initchem, dt_init, nchem, NUM_INITCHEM

   call CheckStop ( nchem == NUM_INITCHEM , dtxt//"NCHEM problem")

   dt=(dtstep - dt_init )/(nchem-NUM_INITCHEM)

   dt(1:NUM_INITCHEM)=dt_initchem     !.. first five timesteps

   if(dtstep<= dt_init )then
      nchem=int(dtstep/dt_initchem)+1
      dt=(dtstep)/(nchem)
   endif
!/ **

   call CheckStop(dtstep<dt_initchem, &
        dtxt//"Error in Solver/makedt: dtstep too small!")
   call CheckStop(nchem>nchemMAX,&
        dtxt//"Error in Solver/makedt: nchemMAX too small!")

   nchem=min(nchemMAX,nchem)

    if( MasterProc ) then

     if( first_dt ) then
      write(*,'(a,4f8.4,i4)')dtxt//'MAKEDT',dtstep,dt_init,dt(1),dt(nchem),nchem
      first_dt = .false.
     end if

     ttot=0.0
     do i=1,nchem
       ttot=ttot+dt(i)
     enddo

     !check that we are using consistent timesteps
     call CheckStop(abs(ttot-dtstep)>1.E-5, &
        dtxt// "Error in Solver/makedt: dtstep and dt not compatible")
    endif

!.. Help variables from Verwer & Simpson
    cc(1)=1.0
    coeff2(1)=1.0/(cc(1)**2+2*cc(1))
    coeff1(1)=(cc(1)+1)**2*coeff2(1)
    dti(1)=((cc(1)+1)/(cc(1)+2))*dt(1)

    do i=2,nchem
      cc(i)=dt(i-1)/dt(i)
      coeff2(i)=1.0/(cc(i)**2+2.0*cc(i))
      coeff1(i)=((cc(i)+1.0)**2)*coeff2(i)
      dti(i)=((cc(i)+1.0)/(cc(i)+2.0))*dt(i)
      cc(i)=1.0/cc(i)
    enddo

 end subroutine makedt
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module ChemSolver
