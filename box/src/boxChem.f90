program Box
  use AeroFunctions_mod,    only: GammaN2O5_om
  use Config_module,        only: Init_Box, use_emis, emis_spec, emis_snap, &
                                  box_emis, emissplit_dir,Hmix, dt, tstart, &
                                  tend, doy, lat, lon, & 
                                  use_cloudj, use_hrlycj, cloudj_indir, &
                                  printJ_day
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
  use cloudj_mod,            only: setup_phot_cloudj 

!UK T, RH
!   use MicroMet_mod, only : rh2num
!   use ZchemData_mod,          only : temp, h2o,rh

   implicit none
   integer :: k, doy_now
   integer :: start, end, rate, cmax, snap, iland
   real :: timeh, ZenRad, hr, time
   logical :: first_tstep=.true.       ! for BOXSOA - tmp
   logical :: first_print=.true.
   real, parameter :: ppb = 2.55e+10
   real, dimension(1) :: zmid = [ 45.0 ]  ! fake array for box-model
   real :: start_time, stop_time, time_used = 0.0
   integer, save :: photstep=-999
   integer :: hr_step

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
   hr_step = int(timeh)

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

     if(use_cloudj) then 
        if(use_hrlycj) then
          if(hr_step>photstep) then
            ! populates rcphotslice, and reassigns phot inds to match cloud on first call
            call setup_phot_cloudj(cloudj_indir,lat,lon,doy_now,hr) 
          endif
        else
          call setup_phot_cloudj(cloudj_indir,lat,lon,doy_now,hr)
        endif
        rcphot(:,1) = rcphotslice(:,1)       
     else ! use MCM rates
        ! do i = 1, NPHOTOLRATES
        !   idj = photol_used(i)
        !   Photolysis rates, with correction (Ftotal) for external data if
        !   needed (we assume that rcphot scales as total Radn)
        !   rcphot(idj,:) = GetPhotol(idj,cos(ZenRad),debug_level=1 )
        ! end do
        rcphot(photol_used,1) = GetPhotol2(photol_used, cos(ZenRad), &
                                    debug_level=debug%Photol) 
     endif

     ! 15 degrees East implies 15/360*24 = 1 hr local time offset
     ! printed in chronological order for Table 1 of GMD_Photolysis paper
     ! in the below, -12 is because the start time is at 12:00, and -1 is because of longitude
     if(hr_step .eq. printJ_day * 24 - 12 - 1 .and. first_print) then
      write(*,*) 'printing photolysis rates (sanity check)'
      write(*,*) '1. O3_O3P    ', rcphot(IDO3_O3P,1)    
      write(*,*) '2. O3_O1D    ', rcphot(IDO3_O1D,1) 
      write(*,*) '3. NO2       ', rcphot(IDNO2,1)   
      write(*,*) '4. H2CO_A    ', rcphot(IDHCHO_H,1)
      write(*,*) '5. H2CO_B    ', rcphot(IDHCHO_H2,1)
      write(*,*) '6. H2O2      ', rcphot(IDH2O2,1)
      write(*,*) '7. CH3OOH    ', rcphot(IDCH3O2H,1)
      write(*,*) '8. NO3       ', rcphot(IDNO3,1)
      write(*,*) '9. HNO2      ', rcphot(IDHONO,1)
      write(*,*) '10.HNO3      ', rcphot(IDHNO3,1)
      write(*,*) '11.HNO4      ', rcphot(IDHO2NO2,1)
      write(*,*) '12.CH3COCHO  ', rcphot(IDRCOCHO,1) ! = MGLYOX
      write(*,*) '13.GLYOX     ', rcphot(IDCHOCHO,1)
      write(*,*) '14.BIACET    ', rcphot(IDCH3COY,1)
      write(*,*) '15.MEK       ', rcphot(IDMEK,1)
      write(*,*) '16.CH3CHO    ', rcphot(IDCH3CHO,1) 
      write(*,*) '17.GLYOXA    ', rcphot(IDGLYOXA,1)
      write(*,*) '18.GLYOXB    ', rcphot(IDGLYOXB,1)
      write(*,*) '19.GLYOXC    ', rcphot(IDGLYOXC,1)
      write(*,*) '20.PAN       ', rcphot(IDPAN,1)  
      write(*,*) '21.CH3COCH3a ', rcphot(IDCH3COCH3,1)
      write(*,*) '22.MCM17     ', rcphot(MCM_J17,1)
      write(*,*) '23.MeAcr     ', rcphot(MCM_J18,1)
      write(*,*) '24.HPALD1    ', rcphot(MCM_J20,1)
      write(*,*) '25.CH3COC2H5 ', rcphot(MCM_J22,1)
      write(*,*) '26.MVK       ', rcphot(MCM_J23,1)
      write(*,*) '26.IPRNO3    ', rcphot(IDiC3H7ONO2,1)
      write(*,*) '28.CHOCHOa   ', rcphot(IDCHOCHO_2CHO,1)
      write(*,*) '29.CHOCHOb   ', rcphot(IDCHOCHO_2CO,1)
      write(*,*) '30.CHOCHOc   ', rcphot(IDCHOCHO_HCHO,1)
      write(*,*) '31.C2H5CHO   ', rcphot(IDC2H5CHO,1)
      write(*,*) '32.ACETOL    ', rcphot(IDACETOL,1)

      ! non-EmChem photolysis rates
      write(*,*) '33.CH4COCH3  ', rcphot(IDCH3COCH3,1)
      write(*,*) '34.MCM J15   ', rcphot(MCM_J15,1)
      write(*,*) '35.MCM J17   ', rcphot(MCM_J17,1)
      write(*,*) '36.MCM J18   ', rcphot(MCM_J18,1)
      write(*,*) '37.MCM J20   ', rcphot(MCM_J20,1)
      write(*,*) '38.MCM J22   ', rcphot(MCM_J22,1)
      write(*,*) '39.MCM J23   ', rcphot(MCM_J23,1)
      write(*,*) '40.iCH3H7ONO2', rcphot(IDiC3H7ONO2 ,1)
      write(*,*) '41.C2H5CHO   ', rcphot(IDC2H5CHO,1)
      write(*,*) '43.GLYALD    ', rcphot(IDGLYALD,1)
      write(*,*) '44.MCRENOL   ', rcphot(IDMCRENOL,1)
      write(*,*) '45.ETP       ', rcphot(IDETP,1)
      write(*,*) '46.ETHP      ', rcphot(IDETHP,1)
      write(*,*) '47.ATOOH     ', rcphot(IDATOOH,1)
      write(*,*) '48.R4P       ', rcphot(IDR4P,1)
      write(*,*) '49.RIPC      ', rcphot(IDRIPC,1)
      write(*,*) '50.PRALDP    ', rcphot(IDPRALDP,1)
      write(*,*) '51.IDHPE     ', rcphot(IDIDHPE,1)
      write(*,*) '52.PIP       ', rcphot(IDPIP,1)
      write(*,*) '53.ITCN      ', rcphot(IDITCN,1)
      write(*,*) '54.INPD      ', rcphot(IDINPD,1)
      write(*,*) '54.MAP       ', rcphot(IDMAP,1)
      write(*,*) '55.RP        ', rcphot(IDRP,1)
      write(*,*) '56.MCM J51   ', rcphot(MCM_J51,1)
      write(*,*) '57.MCM J52   ', rcphot(MCM_J52,1)
      write(*,*) '58.MCM J53   ', rcphot(MCM_J53,1)
      write(*,*) '59.MCM J54   ', rcphot(MCM_J54,1)
      write(*,*) '60.MCM J56   ', rcphot(MCM_J56,1)
      write(*,*) '61.R4N2      ', rcphot(IDR4N2,1)
      write(*,*) '62.MVKN      ', rcphot(IDMVKN,1)
      write(*,*) '63.INPB      ', rcphot(IDINPB,1)
      write(*,*) '64.IHN3      ', rcphot(IDIHN3,1)

      first_print = .false.
     endif

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
