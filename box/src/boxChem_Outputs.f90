module Box_Outputs
  use Config_module,  only : dt, nspec_out, OutSpecs_list, OutGroups_list
  use ChemDims_mod,   only : NSPEC_TOT, NCHEMRATES, NPHOTOLRATES
  use ChemSpecs_mod,  only : species
  use ChemGroups_mod, only : chemgroups  ! for OM25, SIA
  use Debug_module,   only : debug
  use OutputConcs_mod,only : nOutputConcs, InitOutputConcDefs, getOutputConcs, OutConcDefs
  use PhysicalConstants_mod, only : AVOG
  use SmallUtils_mod, only : find_index, num2str, rnum2csvString, str2csvString
  use ZchemData_mod,  only : temp, xChem, tinv, m, n2, o2, h2o, rh,&
                           rcemis, Fgas, Fpart, itemp, rcbio, rct, rcphot

  implicit none

  private

  public :: boxOutputs

  contains
   subroutine boxOutputs(timeh)
    real, intent(in) :: timeh
    integer :: i, ispec, gRO2
    real, allocatable, dimension(:,:) :: xnout
    integer, save ::  io_out, old_hh
    logical, save :: first_call = .true.
    character(len=*), parameter :: dtxt='boxOut:'

    if ( first_call ) then

      old_hh = nint(timeh) 
      open(newunit=io_out, file='./box_outputs.csv',action='write')

     !------------- INIT ----------------------------------------------------
 
     call InitOutputConcDefs( OutSpecs_list(:), OutGroups_list(:) )

     !------------- INIT ----------------------------------------------------

       write(io_out,'(a)') trim( &
         str2csvString( [ character(len=len(species(1)%name)) :: &
         'time', (OutConcDefs(i)%name, i = 1, nOutputConcs), 'H2O' ]) )
       write(6,'(a)') trim( &
         str2csvString( [ character(len=len(species(1)%name)) :: &
         'time', (OutConcDefs(i)%name, i = 1, nOutputConcs), 'H2O' ]) )

       write(io_out,'(a)') trim( &
          str2csvString( [ character(len=len(OutConcDefs(1)%unit)) :: &
            '(h)', (OutConcDefs(i)%unit, i = 1, nOutputConcs), 'ppb' ] ) )

    end if  ! first_call

    allocate(xnout(nOutputConcs,1)) ! 1 for box
    call getOutputConcs(xChem(:,:),Fpart(:,:),m(:),xnout(:,:))
    !print *, dtxt//"XNOUT ", xChem(18,1), m(1), xnout(18,1)

    !write(io_out, "(a)")  trim( &
    !  rnum2csvString( [ timeh, xnout(:,1), h2o(1)*1.e9/m ], '(es18.6)') )
    write(io_out, '(f8.2,9999(:,",",es12.5))')  timeh, xnout(:,1), h2o(1)*1.e9/m

    ! EXTRA outputs, to help understand
    if ( debug%Chem .and. timeh > old_hh+1 ) then
      old_hh = old_hh + 1
      do i=1, NCHEMRATES
        write(*,'(a,i4,f9.4,es12.3)') dtxt//'RCT    ', i, timeh, rct(i,1)
      end do
      do i = 1, NPHOTOLRATES  ! size(rcphot, dim=2)
        write(*,'(a,i4, f8.1,es12.3,f12.4)') dtxt//'RCPHOT ', i, timeh, rcphot(i,1)
      end do
      gRO2=find_index("RO2",chemgroups%name)
      if ( gRO2 > 0 ) then
        do i = 1, size(chemgroups(gRO2)%specs)
          ispec=chemgroups(gRO2)%specs(i)
          write(*,'(a,i4,f8.1,es12.3)') dtxt//'RO2POOL '//species(ispec)%name, i, &
           timeh, xChem(ispec,1)*1.e12/m ! ppt
        end do
      end if
    end if
    first_call = .false.
 end subroutine boxOutputs
end module Box_Outputs

