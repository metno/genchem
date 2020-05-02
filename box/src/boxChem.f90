program Box
  use AeroFunctions_mod,    only: GammaN2O5_om
  use Config_module,        only: Init_Box, use_emis, emis_spec, emis_snap, &
                                  box_emis, emissplit_dir,Hmix, dt, tstart, &
                                  tend, doy, lat, lon
  use Box_Outputs,           only: boxOutputs
  use ChemDims_mod
  use ChemFunctions_mod
  use ChemGroups_mod
  use ChemRates_mod,         only : photol_used, setPhotolUsed, setChemRates
  use ChemSpecs_mod              !!only: define_chemicals, species, e.g. O3
  use Chemsolver,            only: chemsolve
  use Debug_module,          only: debug
  use DefPhotolysis_mod            ! ZenithAngle and Photol rates
  use EmisSplit_mod,         only: NSECTORS, RdEmisSplit, &
                                   getTestEmisSplit, doEmisSplit
  use OrganicAerosol_mod,    only: ORGANIC_AEROSOLS, & !true if SOA
                                    Init_OrganicAerosol, OrganicAerosol
  use PhysicalConstants_mod, only: AVOG
  use SmallUtils_mod,        only: num2str, find_index, to_upper
  use ZchemData_mod              ! xChem,rcphot, etc.

!UK T, RH
!   use MicroMet_mod, only : rh2num
!   use ZchemData_mod,          only : temp, h2o,rh

   implicit none
   integer :: k, doy_now
   integer :: start, end, rate, cmax, snap, iland
   real :: timeh, ZenRad, hr, time
   logical :: first_tstep=.true.       ! for BOXSOA - tmp
   real, parameter :: ppb = 2.55e+10
   real, dimension(1) :: zmid = [ 45.0 ]  ! fake array for box-model
   real :: start_time, stop_time, time_used = 0.0

   if ( digits(1.0)<50 ) stop "COMPILE ERROR:Needs double prec., e.g. f90 -r8"

!=====================================================

   call system_clock(start,rate,cmax)  ! wall time

!=====================================================

   call define_chemicals()

!=====================================================

!-----------------------------------------------------
   call Init_Box()

   call Init_ChemGroups()

   time=tstart

   call SetPhotolUsed()
   call Init_OrganicAerosol(zmid,first_tstep)

   first_tstep=.false. ! for BOXSOA  - in EMEP, this would be after all i,j

    k=1
    if (use_emis) then
     !===============================================
     ! Emissions, only needed in the beginning as long as there's no daily
     !  cycle or so.  box_emis values are set by Config_module, 
     ! eg box_emis(2) = 18.3 NOx, kg/m2/day

     ! get emission fractions in each SNAP, test data
     call getTestEmisSplit(box_emis,emis_snap)

     ! read emissplit files
     call RdEmisSplit(emissplit_dir,debug%Emis,debug%VOC)

     ! Get emissions of each species, now in molec/cm2/s:
     iland = 1 ! Box simplification
     call doEmisSplit(emis_snap,'kgkm2day',iland,emis_spec)

    ! And convert to molec/cm3/s (assuming here k=1):
      rcemis(:,:) = 0.0
      do snap = 1, NSECTORS
        rcemis(:,1) = rcemis(:,1) + emis_spec(:,snap)/Hmix
      end do

    end if

   !===============================================

   call boxOutputs(time/3600.0) ! Initial output

!print *, "DBGOD pre0 ", xChem(1,1), Dchem(1,1)

  do while ( time < tend)
   time=time + dt
   timeh = time/3600.0

  ! Many MCM/CRI tests use variable RH and T, based upon time (from
  ! midnight). Can uncomment below if wanted.
  ! Uses: tave=287.46, tamp = 5.4
  !  temp(:)= 5.4*sin((7.2722e-5*time)-1.9635)+ 287.46
  !  rh=23.*sin((7.2722E-5*time)+1.1781)+66.5
  !!  print *, "TESTING AERO GammaN2O5_om  ", GammaN2O5_om(rh)
  !  h2o=rh2num(rh,temp)


   !===============================================
   ! Photolysis
   !-----------------------------------------------
   if (NPHOTOLRATES > 0) then

     doy_now = doy + floor((timeh)/24.)
     hr = mod(timeh, 24.)
     ZenRad = ZenithAngle(doy_now, hr, lat, lon)

     !do i = 1, NPHOTOLRATES
     !  idj = photol_used(i)

      ! Photolysis rates, with correction (Ftotal) for external data if
      ! needed (we assume that rcphot scales as total Radn)
     !  rcphot(idj,:) = GetPhotol(idj,cos(ZenRad),debug_level=1 )
     !end do
     rcphot(photol_used,1) = GetPhotol2(photol_used, cos(ZenRad), &
                                debug_level=debug%Photol)

   end if

   !===============================================
   !Some BVOC rates may use SUN as scalar
    call Update_SUN(time)
   !===============================================
   ! Rate coefficients (rct)

    call  setChemRates()
    !print '(a,f7.1,3f8.2,3es12.2)', "FAKE RC ", hr, SUN, rct(1,1)

   !===============================================

    if( ORGANIC_AEROSOLS) call OrganicAerosol()


    call cpu_time(start_time) ! processor time, just the chem soln

    call chemsolve( k, dt, xChem(:,1),Dchem(:,1),debug_level=1,&
                     time=time,itern=1)

    call cpu_time(stop_time)
    time_used = time_used + stop_time - start_time

  if ( abs(timeh - nint(timeh) ) < 0.5*dt/3600 ) then  ! Print every hour
      !print *, 'CPU TIME ', start_time,stop_time, time_used,&
      ! stop_time-start_time
     call boxOutputs(timeh)

   !===============================================
    print "(a,f7.1,2f9.3,(a,f6.2),9(a,es10.2))", "RES ", hr, &
     timeh, (100.0*time)/tend, " %, O3=",xChem(O3,1)/ppb, ' ppb, OH=',&
       xChem(OH,1), ' /cm3, HO2=', xChem(HO2,1), ' /cm3'
  end if 
  end do ! time loop


  call system_clock(end)

  print '(3(a,f12.3))', 'execution of program took ', &
    float(end - start)/rate, ' seconds, chem cpu: ', time_used, '   dt=', dt

end program Box
