!>>>>>>>>   adapted from the standalone cloudj program for use in box-chem    <<<<<<<<<<<<
!>>>>>>>>cloud-JX code (includes fractional cloud treatments) ver 7.3 (2/2015)<<<<<<<<<<<<
module CloudJ_mod
      
    USE FJX_CMN_MOD
    USE FJX_SUB_MOD
    USE FJX_INIT_MOD
    USE CLD_SUB_MOD, ONLY : CLOUD_JX
  
    use Config_module, only: cloudj_indir
    use ZchemData_mod, only: rcphotslice
    use DefPhotolysis_mod
    use ChemDims_mod,  only: NPHOTOLRATES
    use ChemRates_mod, only: photol_used, setPhotolUsed
  
    contains 
  
    subroutine setup_phot_cloudj(cloudj_dir,box_lat,box_lon,doy,utau)
    
          implicit none
  
          character(len=100), intent(in) :: cloudj_dir
          real, intent(in) :: box_lat
          real, intent(in) :: box_lon
          integer, intent(in) :: doy
          real, intent(in) :: utau
          logical, save :: first_call=.true.
  
  
          real, save, dimension(L1_) :: ETAA,ETAB, ZOFL,RI,TI,CLDP,AER1,AER2
          real, save, dimension(L_) :: WLC,WIC
          real, save :: GMTAU,PHOTAU,ALBEDO, XLNG,YLAT,XGRD,YGRD,PSURF, SCALEH
          real, save :: CF,PMID,PDEL,ZDEL,ICWC,F1,ZKM
          integer, save, dimension(L1_):: NAA1,NAA2
          integer, save :: MYEAR, MONTH, IDAY, JLAT, ILON, NHRMET
          integer, save :: I,J,K,L,N,JP, NN
          integer, save :: NRAN, RANSEED, LTOP, NJXX
          real, save, dimension(LWEPAR) :: CLDFRW,CLDIWCW,CLDLWCW   ! WCW=Cloud Water Content (g/g)
          character*6, save, dimension(JVN_)   ::  TITLJXX
          real, save, dimension(L_,8,4)  :: ZPJCLD,ZPJAVG
          real, save, dimension(21,4)    :: ERRJ, ERRJJ,ERRJ2
          integer, save :: JP04,JP09,JP11,JP15, ICLD
          integer, save, dimension(8)      :: JCNT
          character*11, save, dimension(4) ::  TITJX
    
          character*8, dimension(8), parameter :: TITERR =  &
          ['clearsky','avgcloud','cldf^3/2','ICA-beam', &
           'ICAs-ran','QCAs-mid','QCAs-avg','all ICAs']
    
          integer, save              :: CLDFLDS, CLDFLDX
    
    
    !--------------------key params sent to CLOUD_JX-------------------------
          real, save                     :: U0,SZA,REFLB,SOLF,  CLDCOR
          real, save                     :: FG0
          logical, save                  :: LPRTJ
          real,  save, dimension(L1_+1)  :: PPP,ZZZ
          real,  save, dimension(L1_  )  :: TTT,DDD,RRR,OOO
          real,  save, dimension(L1_)    :: LWP,IWP,REFFL,REFFI
          real,  save, dimension(L1_)    :: CLF
          integer, save, dimension(L1_)    :: CLDIW
          real,  save, dimension(L1_,AN_):: AERSP
          integer, save, dimension(L1_,AN_):: NDXAER
          real,  save, dimension(L_,JVN_):: VALJXX
          integer, save                    :: CLDFLAG,NRANDO,IRAN,LNRG,ICNT
          integer, save                    :: NICA,JCOUNT
    !---U0 = cos (SZA), SZA = solar zenith angle
    !---REFLB = Lambertian reflectivity at the Lower Boundary
    !---SOLF = solar flux factor for sun-earth distance
    !---FG0 = scale for asymmetry factor to get equivalent isotropic (CLDFLAG=3 only)
    !---LPRTJ = .true. = turn on internal print in both CLOUD_JX & PHOTO_JX
    !--- P = edge press (hPa), Z = edge alt (m), T = layer temp (K)
    !--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
    !--- R = layer rel.hum.(fraction)
    !---LWP/IWP = Liquid/Ice water path (g/m2)
    !---REFFL/REFFI = R-effective(microns) in liquid/ice cloud
    !---CLF = cloud fraction (0.0 to 1.0)
    !---CLDIW = integer denoting cloud in layer: 1=water, 2=ice, 3=both
    !---AERSP = aerosol path (g/m2) & NDXAER = aerosol index type
    !---  aerosols are dimensioned with up to AN_ different types in an ICA layer
    !---L1_ = parameter, dim of profile variables, L_+1 for top (non CTM) layer
    !---AN_ = parameter, dim of number of aerosols being passed
    !---VALJXX = J-values from CLOUD_JX & PHOTO_JX
    !---JVN_ = dim of max number of J-s reported out (in the order of fast-JX, not CTM)
    !---CLDFLAG = integer index for type of cloud overlap
    !---CLOUD_JX:   different cloud schemes
    !---CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
    !       CLDFLAG = 1  :  Clear sky J's
    !       CLDFLAG = 2  :  Averaged cloud cover
    !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
    !       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
    !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
    !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
    !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
    !       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)
    !---NRANDO = number of random ICAs to do J's for (CLDFLAG=4)
    !---IRAN = starting index for random number selection
    !---LNRG = flag for setting max-ran overlap groups:
    !---     =0   break max overlap groups at cloud fraction = 0
    !---     =3   else fixed 3 layers (1:9, 9:last LWcloud, LWclud+1:LTOP)
    !---     else(=6) fixed correlated length max-overlap layers
    !---NICA = total number of ICAs
  
    !---fast-JX:  INIT_JX is called only once to read in & store all fast-JX data:
    !              also sets up random sequence for cloud-JX
    !-----------------------------------------------------------------------
    if( first_call ) call INIT_FJX (TITLJXX,JVN_,NJXX)
    !-----------------------------------------------------------------------
    !--Set up atmosphere for a single column and time for J-values calculation
    !--Nominally taken from CTM, but for standalone here is read in
          open (77,file=trim(cloudj_dir)//'CTM_GrdCld.dat',status='old')
          read (77,*)
          read (77,'(i6)') MYEAR
          read (77,'(i6)') IDAY
          read (77,'(i6)') MONTH
          read (77,'(i6)') NHRMET
          read (77,'(2f6.1)') GMTAU,PHOTAU
          read (77,'(f6.1,i4)') YLAT, JLAT
          read (77,'(f6.1,i4)') XLNG, ILON
              YGRD = box_lat*CPI180
              XGRD = box_lon*CPI180
          read (77,'(f6.1)') PSURF
          read (77,'(f6.1)') ALBEDO
          read(77,'(2i5,f5.1)') LNRG,NRANDO,CLDCOR
          read(77,'(f5.2)') FG0
          read(77,'(2i5)') CLDFLDS,CLDFLDX
            ! write(6,'(a,2i4,a,f7.4)') 'LNRG = ',LNRG,NRANDO,'Cloud-Corel',CLDCOR
            ! write(6,'(a,f10.4)') 'FG0 (assym for direct beam)', FG0
            ! write(6,'(a,2i5)') '#CLDFLDS(out of) = ',CLDFLDS, CLDFLDX
            ! write(6,'(2i5,5x,a,i5)') LPAR,LWEPAR, 'LPAR / LWEPAR', L1_
          read (77,*)
          do L = 1,LPAR+1
            read (77,'(i3,1x,2f11.7,2x,f5.1,f5.2,f11.2,2(f7.3,i4))') &
                          J,ETAA(L),ETAB(L),TI(L),RI(L),ZOFL(L) &
                         ,AER1(L),NAA1(L),AER2(L),NAA2(L)
          enddo
          do L = 1,L1_
           PPP(L) = ETAA(L) + ETAB(L)*PSURF
          enddo
          
          MONTH = int(doy/30) ! first approx.; clim applies additional M = max(1,min(12,MONTH))
    !---ACLIM_FJX sets climatologies for O3, T, D & Z - overwrite with CTM data
          call ACLIM_FJX (box_lat, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
          do L = 1,L_
           TTT(L) = TI(L)
           RRR(L) = RI(L)
          enddo
           ZZZ(1)  = 16.d5*log10(1013.25d0/PPP(1))        ! zzz in cm
          do L = 1,L_
           DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
           SCALEH      = 1.3806d-19*MASFAC*TTT(L)
           ZZZ(L+1) = ZZZ(L) -( log(PPP(L+1)/PPP(L)) * SCALEH )
          enddo
           ZZZ(L1_+1) = ZZZ(L1_) + ZZHT
          REFLB = ALBEDO
          LPRTJ = .false.
    
    !---following is readin for cloud data, currently has 160 atmospheres
    !   from tropical T319 ECMWF atmosphere used in UCI CTM.
            ZPJAVG(:,:,:) = 0.d0
    
            ERRJJ(:,:) = 0.d0
            ERRJ2(:,:) = 0.d0
            JCNT(:) = 0
    
          do ICLD=1,1! CLDFLDS; default clear sky Js, no iterations over different cloud fields needed
    
            read (77,*)
           do L = LWEPAR,1,-1
            read (77,'(i3,1p,e14.5,28x,2e14.5)') J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
           enddo
  
           close(77)
    !---load CTM-based data on ozone and aerosols on top of climatology
    !---to convert kg (STT) in grid cell w/AREA (m2) to # molecules/cm^2
    !      D_O3(I,J,L) = 6.023d26*STT(I,J,L,1)/48.d0  *1.d-4 /AREAXY(I,J)
    !---to convert kg (STT) in grid cell w/AREA (m2) to PATH (g/m^2)
    !      P_AERSL(I,J,L) = STT(I,J,L,1)*1.d-3/AREAXY(I,J)
    !---this data should be available from the CTM somewhere.
    !---   ZZZ(1:L_+1) is geometric altitude (cm), approx is OK.
    !---AERSP = aerosol path (g/m2) and NDXAER() must come from CTM
    !    the index must match one in the std or UMich data sets (or add your own).
            AERSP(:,:)  = 0.d0
            NDXAER(:,:) = 0
            do L = 1,L_
              NDXAER(L,1) = NAA1(L)
              AERSP(L,1)  = AER1(L)
              NDXAER(L,2) = NAA2(L)
              AERSP(L,2)  = AER2(L)
            enddo
    !---convert cloud data from our EC met fields into cloud data for cloud-JX
    !---     init data = cloud fraction, and water content (g/m3) averaged over cell
    !---     needs: cloud fraction, ice- and liq-water path (in cloud)
    !---     and R-effective of ice and liquid clouds
            LTOP  = LWEPAR
          if (maxval(CLDFRW) .le. 0.005d0) then
            IWP(:) = 0.d0
            REFFI(:) = 0.d0
            LWP(:) = 0.d0
            REFFL(:) = 0.d0
          endif
          do L = 1,LTOP
              CLDIW(L) = 0
              CF  = CLDFRW(L)
            if (CF .gt. 0.005d0) then
              CLF(L) = CF
              WLC(L) = CLDLWCW(L) / CF
              WIC(L) = CLDIWCW(L) / CF
    !  CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
             if (WLC(L) .gt. 1.d-11) CLDIW(L) = 1
             if (WIC(L) .gt. 1.d-11) CLDIW(L) = CLDIW(L) + 2
            else
              CLF(L) = 0.d0
              WLC(L) = 0.d0
              WIC(L) = 0.d0
            endif
          enddo
    
    !---derive R-effective for clouds:  the current UCI algorithm - use your own
          do L = 1,LTOP
    !---ice clouds
            if (WIC(L) .gt. 1.d-12) then
                PDEL = PPP(L) - PPP(L+1)
                ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! m
              IWP(L) = 1000.d0*WIC(L)*PDEL*G100    ! g/m2
              ICWC =        IWP(L) / ZDEL          ! g/m3
              REFFI(L) = 164.d0 * (ICWC**0.23d0)
            else
              IWP(L) = 0.d0
              REFFI(L) = 0.d0
            endif
    !---water clouds
            if (WLC(L) .gt. 1.d-12) then
                PMID = 0.5d0*(PPP(L)+PPP(L+1))
                PDEL = PPP(L) - PPP(L+1)
              F1   = 0.005d0 * (PMID - 610.d0)
              F1   = min(1.d0, max(0.d0, F1))
              LWP(L) = 1000.d0*WLC(L)*PDEL*G100     ! g/m2
              REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
            else
              LWP(L) = 0.d0
              REFFL(L) = 0.d0
            endif
          enddo
    !---cloud input as interpreted by fast_JX
    !---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
    !           ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L)
    !           ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0)
    
    !---fast-JX:  SOLAR_JX is called  once per grid-cell to set U0, SZA, SOLF
    !--- your CTM code may have its own way of calculating and passing these quantities
    !-----------------------------------------------------------------------
          call SOLAR_JX(utau,doy,YGRD,XGRD, SZA,U0,SOLF)
    !-----------------------------------------------------------------------
  
          TTT(L1_)=230.0 ! needed to prevent series of nans, strangely (prob. compiler?)
          DDD(L1_)=2.00E+15
          OOO(L1_)=1.00E+09
  
          if (LPRTJ) then
            write(6,'(a,f8.3,3f8.5)')'solar zenith angle, solar-f' &
                   ,SZA,SOLF,U0,REFLB
            write(6,'(a,f8.3,f8.3)') 'cldj lat/lng',YLAT,XLNG
            write(6,'(a,f8.3,f8.3)') 'cldj jlat/jlng',REAL(JLAT),REAL(ILON)
            call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
          endif
          
          IDAY = doy
          GMTAU = utau
    !  locate the position of random number sequence based on year/day/hour
          IRAN = 13+ILON+3*JLAT+5*(MYEAR-1900)+7*IDAY + 11*nint(GMTAU)
    
              ZPJCLD(:,:,:) = 0.d0
    
    !  >>>>>test mode -- loop across all CLDFLAG's to accumulate errors
    
          CLDFLAG = 1 ! clear-sky = 1
          
    !=======================================================================
    
          call CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT,       &
                 DDD,RRR,OOO,   LWP,IWP,REFFL,REFFI,  CLF,CLDCOR,CLDIW, &
                 AERSP,NDXAER,L1_,AN_,VALJXX,JVN_,                      &
                 CLDFLAG,NRANDO,IRAN,LNRG,NICA,JCOUNT)
    
              JCNT(CLDFLAG) = JCNT(CLDFLAG) + JCOUNT
    !=======================================================================
  
          do J=1,NRATJ ! populate the photolysis array. NRATJ is number of phot reactions from j2j input file
            rcphotslice(J,1) = VALJXX(1,JIND(J))*JFACTA(J)
          end do
  
          if (first_call) then 
            do I=1,NRATJ
              ! name to search for must be exactly 10 characters long; first ten characters in FJX_j2j.dat
              if ('O3        ' .eq. JLABEL(I)(:10) ) IDO3_O3P  = I
              if ('O3(1D)    ' .eq. JLABEL(I)(:10) ) IDO3_O1D  = I
              if ('NO2       ' .eq. JLABEL(I)(:10) ) IDNO2     = I   
              if ('H2COa     ' .eq. JLABEL(I)(:10) ) IDHCHO_H  = I
              if ('H2COb     ' .eq. JLABEL(I)(:10) ) IDHCHO_H2 = I 
              if ('H2O2      ' .eq. JLABEL(I)(:10) ) IDH2O2    = I   
              if ('CH3OOH    ' .eq. JLABEL(I)(:10) ) IDCH3O2H  = I   
              if ('NO3c      ' .eq. JLABEL(I)(:10) ) IDNO3     = I ! lumped EMEP
              if ('HNO4      ' .eq. JLABEL(I)(:10) ) IDHO2NO2  = I
              if ('CH3COCHO  ' .eq. JLABEL(I)(:10) ) IDRCOCHO  = I
              if ('BIACET    ' .eq. JLABEL(I)(:10) ) IDCH3COY  = I
              if ('CH3COC2H5 ' .eq. JLABEL(I)(:10) ) IDMEK     = I
              if ('CH3CHO    ' .eq. JLABEL(I)(:10) ) IDCH3CHO  = I
              if ('CHOCHOa   ' .eq. JLABEL(I)(:10) ) IDGLYOXA  = I
              if ('CHOCHOb   ' .eq. JLABEL(I)(:10) ) IDGLYOXB  = I
              if ('CHOCHOc   ' .eq. JLABEL(I)(:10) ) IDGLYOXC  = I
              if ('PAN       ' .eq. JLABEL(I)(:10) ) IDPAN     = I ! only in CJX
              if ('HNO3      ' .eq. JLABEL(I)(:10) ) IDHNO3    = I
              if ('HNO2      ' .eq. JLABEL(I)(:10) ) IDHONO    = I
              if ('GLYOX     ' .eq. JLABEL(I)(:10) ) IDCHOCHO  = I ! lumped EMEP 
              if ('NO3a      ' .eq. JLABEL(I)(:10) ) IDNO3_NO  = I 
              if ('NO3b      ' .eq. JLABEL(I)(:10) ) IDNO3_NO2 = I 
              if ('CH3COCH3a ' .eq. JLABEL(I)(:10) ) IDACETON  = I ! only a-channel
              if ('N2O5      ' .eq. JLABEL(I)(:10) ) IDN2O5    = I
  
              ! duplicates with different names for historical reasons
              if ('CHOCHOa   ' .eq. JLABEL(I)(:10) ) IDCHOCHO_2CHO = I
              if ('CHOCHOb   ' .eq. JLABEL(I)(:10) ) IDCHOCHO_2CO  = I
              if ('CHOCHOc   ' .eq. JLABEL(I)(:10) ) IDCHOCHO_HCHO = I
  
              ! non-EmChem photolysis rates
              if ('CH3COCH3a ' .eq. JLABEL(I)(:10) ) IDCH3COCH3  = I
              if ('MCM15     ' .eq. JLABEL(I)(:10) ) MCM_J15     = I
              if ('MCM17     ' .eq. JLABEL(I)(:10) ) MCM_J17     = I
              if ('MeAcr     ' .eq. JLABEL(I)(:10) ) MCM_J18     = I 
              if ('HPALD1    ' .eq. JLABEL(I)(:10) ) MCM_J20     = I
              if ('CH3COC2H5 ' .eq. JLABEL(I)(:10) ) MCM_J22     = I
              if ('MVK       ' .eq. JLABEL(I)(:10) ) MCM_J23     = I
              if ('IPRNO3    ' .eq. JLABEL(I)(:10) ) IDiC3H7ONO2 = I 
              if ('C2H5CHO   ' .eq. JLABEL(I)(:10) ) IDC2H5CHO   = I
              if ('HAC       ' .eq. JLABEL(I)(:10) ) IDACETOL    = I
              if ('GLYC      ' .eq. JLABEL(I)(:10) ) IDGLYALD    = I
              if ('MCRENOL   ' .eq. JLABEL(I)(:10) ) IDMCRENOL   = I
              if ('ETP       ' .eq. JLABEL(I)(:10) ) IDETP       = I
              if ('ETHP      ' .eq. JLABEL(I)(:10) ) IDETHP      = I
              if ('ATOOH     ' .eq. JLABEL(I)(:10) ) IDATOOH     = I
              if ('R4P       ' .eq. JLABEL(I)(:10) ) IDR4P       = I
              if ('RIPC      ' .eq. JLABEL(I)(:10) ) IDRIPC      = I
              if ('PRALDP    ' .eq. JLABEL(I)(:10) ) IDPRALDP    = I
              if ('IDHPE     ' .eq. JLABEL(I)(:10) ) IDIDHPE     = I
              if ('PIP       ' .eq. JLABEL(I)(:10) ) IDPIP       = I
              if ('ITCNa     ' .eq. JLABEL(I)(:10) ) IDITCN      = I
              if ('INPDa     ' .eq. JLABEL(I)(:10) ) IDINPD      = I
              if ('MAP       ' .eq. JLABEL(I)(:10) ) IDMAP       = I
              if ('RP        ' .eq. JLABEL(I)(:10) ) IDRP        = I
              if ('MENO3     ' .eq. JLABEL(I)(:10) ) MCM_J51     = I
              if ('ETNO3     ' .eq. JLABEL(I)(:10) ) MCM_J52     = I
              if ('NPRNO3    ' .eq. JLABEL(I)(:10) ) MCM_J53     = I
              if ('IPRNO3    ' .eq. JLABEL(I)(:10) ) MCM_J54     = I
              if ('PROPNN    ' .eq. JLABEL(I)(:10) ) MCM_J56     = I
              if ('R4N2      ' .eq. JLABEL(I)(:10) ) IDR4N2      = I
              if ('MVKN      ' .eq. JLABEL(I)(:10) ) IDMVKN      = I
              if ('INPB      ' .eq. JLABEL(I)(:10) ) IDINPB      = I
              if ('IHN3      ' .eq. JLABEL(I)(:10) ) IDIHN3      = I
            enddo ! NRATJ
  
            ! place indices in the photol_used array
            call setPhotolUsed()
  
            ! verify that all phot rates were found (i.e. greater than 0) (only in EMEP)
            ! do I=1,NPHOTOLRATES
            !   if ( photol_used(I) < 0 ) write(*,*) 'Phot. rate not found in CJX.'
            !   if ( photol_used(I) < 0 ) stop 
            ! enddo
  
          endif ! first_call 
          
          first_call = .false.
    enddo ! ICLD
            
    end subroutine setup_phot_cloudj
  
  endmodule CloudJ_mod