!---addX_FJX_sss.f  -- with user supplied subroutine that supplies X-section x Q-yield
!---generates fast-JX 18-bin X-sections  revised and updated v73 (mprather,4/2015)
      implicit none
      integer, parameter :: NB_ = 100
      integer, parameter :: NS_ = 40000
      integer, parameter :: NJ_ = 18

      real*8 SRB(15,NJ_)
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN, XBIN
      real*8, dimension(NJ_) :: FFBIN,AABIN,XXBIN
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,K1,K2,NB77  ,I778
      integer INIT
      real*8 W(NS_),F(NS_)
      integer IBINJ(NS_)
      real*8 W1,W2, XT,XP,XM, MM(3), WW,XNEW, QNEW
      character*6 TITLNEW,TITLNX,TITLNQ
      character*1 ISX, ISP, ISXP
      character*72 TITLTBL(3), TITLTQ(3)
      Character*20 TITL77
      integer     TTT(3), PPP(3)

      open (1, file='wavel-bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(i5)') NB
        if (NB .gt. NB_) stop
        read(1,'(5x,f8.3)') (WBIN(I), I=1,NB+1)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
        read(1,*)
        read(1,'(5x,i5)') (IJX(I),I=16,NB)
      close (1)

      open (2, file='solar-p05nm-UCI.dat', status='OLD')
        read(2,*)
        read(2,*)
        read(2,'(f10.4,e10.3)') (W(J),F(J), J=1,NS_)
      close (2)


      do I778 = 0,1

        if (I778 .eq. 0) then
          WBIN(NB+1) = 850.
        else
          WBIN(NB+1) = 778.
        endif

!---now assign bin #(I=1:77) to each p05nm microbin J (1:40000)
        IBINJ(:) = 0
      do I=1,NB
         W1 = WBIN(I)
         W2 = WBIN(I+1)
        do J=1,NS_
          if (W(J) .gt. W1) goto 11
        enddo
          J = NS_ + 1
   11     J1 = J
        do J=J1,NS_
          if (W(J) .gt. W2) goto 12
        enddo
          J = NS_ + 1
   12     J2 = J-1
        do J=J1,J2
          IBINJ(J) = I
        enddo
      enddo
!!!!this binning does not interpolate and is OK for large bins
!     it has 7% error in the very short wavel S-R bins of pratmo.

!!!!!!!!!!!!!!!!!!initialization call to user subroutine!!!!!!!!!!!!!!

      INIT = 0
      call Q_O3 (W(1),XT,XP,XM, XNEW,
     &     INIT, MM,TTT,PPP,ISX,ISP,TITLNQ,TITLTBL)
!--store returns from Q_O3 separately
          TITLTQ(:) = TITLTBL(:)

      INIT = 0
      call X_O3 (W(1),XT,XP,XM, XNEW,
     &     INIT, MM,TTT,PPP,ISX,ISP,TITLNX,TITLTBL)

        ISP = ' '
        ISX = ' '
        ISXP = ' '

!!!!!!!call to user subroutine for XO3 and Q1D  high res and 3 Temperatures

        do K = 1,1
        write(6,*) 'T',K, TITLTBL(k)
        write(6,*) 'T',K, TITLTQ(k)
        enddo


      do K = 1,3

! - new, use the same T's for xO3 and q1D

       if (TTT(K) .gt. 0) then
        XT = TTT(K)
!        XP = PPP(K)
!        XM = MM(K)
!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
         XBIN(:) = 0.d0
! primary high-resolution wavelength loop - generate input for pratmo reference J's
        do J = 1,NS_
           I = IBINJ(J)
         if (I .gt. 0) then

      call Q_O3 (W(J),XT,XP,XM, QNEW,
     &     INIT, MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)

      call X_O3 (W(J),XT,XP,XM, XNEW,
     &      INIT, MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)

          WW = W(J)
          FBIN(I) = FBIN(I) + F(J)
          XBIN(I) = XBIN(I) + F(J)*XNEW
          ABIN(I) = ABIN(I) + F(J)*XNEW*QNEW
         endif
        enddo
        do I=1,NB
          if (XBIN(I) .gt. 0.d0) then
            ABIN(I) = ABIN(I)/XBIN(I)
          endif
          if (FBIN(I) .gt. 0.d0) then
            NB77 = I
            XBIN(I) = XBIN(I)/FBIN(I)
          endif
        enddo

!---UCI reference std 77-bin data    (only for W-end = 850 nm)
        if (I778.eq.0) then
          TITL77 = TITLTBL(1)
          write(6,'(a20,f5.0/(1p,8e10.3))') TITL77,XT,(XBIN(I),I=1,NB77)
          TITL77 = TITLTQ(1)
          write(6,'(a20,f5.0/(1p,8e10.3))') TITL77,XT,(ABIN(I),I=1,NB77)
        endif

!---secondary sum 77-bin pratmo ==> 18-bin fast-JX
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
         FFBIN(:) = 0.d0
         XXBIN(:) = 0.d0
         AABIN(:) = 0.d0
        do I=16,NB
          J = IJX(I)
          FFBIN(J) = FFBIN(J) + FBIN(I)
          XXBIN(J) = XXBIN(J) + FBIN(I)*XBIN(I)
          AABIN(J) = AABIN(J) + FBIN(I)*XBIN(I)*ABIN(I)
        enddo
        do I=1,15
         do J=1,NJ_
          FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
          XXBIN(J) = XXBIN(J) + FBIN(I)*XBIN(I)*SRB(I,J)
          AABIN(J) = AABIN(J) + FBIN(I)*XBIN(I)*ABIN(I)*SRB(I,J)
         enddo
        enddo
        do J=1,NJ_
         if (XXBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/XXBIN(J)
         if (FFBIN(J) .gt. 0.d0) XXBIN(J) = XXBIN(J)/FFBIN(J)
        enddo

! fast-JX v68 data bins for 'FJX_spec.dat'
        if (I778.eq.0) then
        TITL77 = TITLTBL(1)
        write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITL77,ISXP,TTT(K), XXBIN
        TITL77 = TITLTQ(1)
        write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITL77,ISXP,TTT(K), AABIN
        endif

! fast-J v73 data for 'FJX_spec.dat'
        if (I778.eq.1) then
        write(6,'(a)') TITLTBL(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',XXBIN(1:6),TITLNX
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'b',XXBIN(7:12),TITLNX
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',XXBIN(13:18),TITLNX

        write(6,'(a)') TITLTQ(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',AABIN(1:6),TITLNQ
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'b',AABIN(7:12),TITLNQ
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',AABIN(13:18),TITLNQ
        endif

       endif
      enddo

      write(6,*) ' WBIN(end) =  ',WBIN(NB+1)


      enddo    ! I778

      stop
      end



c-----------------------------------------------------------------------
      subroutine Q_O3 (WW,XT,XP,XM,XNEW,
     &     INIT, MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)
c-----------------------------------------------------------------------
      implicit none
      real*8, intent(in) :: WW,XT,XP,XM
      real*8, intent(out) :: XNEW
      integer, intent(inout) :: INIT
      real*8, intent(out) :: MM(3)
      integer, intent(out) :: TTT(3),PPP(3)
      character*1, intent(out) :: ISX,ISP
      character*6, intent(out) :: TITLNEW
      character*90, intent(out) :: TITLTBL(3)

      real*8  T,QO1D, Q1,Q2,Q1Q2,EX1,EX2,EX3
!---quantum yield for O3 + hv => O(1D),
!---parametric fit for range w = 306-328 nm, T = 200-320K
!---JPL_2010 standard tables:
!     Table 4A-5. Parameters for the Calculation of O(1D) Quantum Yields.
      real*8, parameter::  X1 = 304.225d0
      real*8, parameter::  X2 = 314.957d0
      real*8, parameter::  X3 = 310.737d0
      real*8, parameter::  W1 = 5.576d0
      real*8, parameter::  W2 = 6.601d0
      real*8, parameter::  W3 = 2.187d0
      real*8, parameter::  A0 = 0.90d0
      real*8, parameter::  A1 = 0.8036d0
      real*8, parameter::  A2 = 8.9061d0
      real*8, parameter::  A3 = 0.1192d0
      real*8, parameter::  A4 = 0.0765d0
      real*8, parameter::  A4x = 0.08d0
      real*8, parameter::  V1 = 0.0d0
      real*8, parameter::  V2 = 825.518d0
      real*8, parameter::  RG = 0.695d0
      character*80  ::  TABLE

      TABLE='Table 4A-5. Parameters for the Calculation of O(1D)'

      if (INIT .eq. 0) then
        TITLNEW = 'O3(1D)'
        TTT(1) = 200
        TTT(2) = 260
        TTT(3) = 320
        TITLTBL(1) ='O3(1D).Qyld O3=O(1D)+O2. JPL10 3/2013'
        write(6,'(2a)') ' species:',TITLNEW
        write(6,'(a80)') TABLE
!!!!        TITLNEW = TITLTBL(1)
        INIT = 1
        MM(:) = 0.d0
        PPP(:) = 0
        ISX = ' '
        ISP = ' '
      else

c--set QO1D = 0.0 for W.gt.340.,  = 0.48 at 193 nm, = 0.90 at 225 nm
       if (WW .lt. 306.d0) then
          QO1D = A0
       elseif (wW .gt. 328.d0) then
          QO1D = A4x
       else
         T = min (320.d0, max(200.d0, XT))
         Q1 = exp(-V1/(RG*T))
         Q2 = exp(-V2/(RG*T))
         Q1Q2 = Q1/(Q1+Q2)
         EX1 = exp( -((X1-WW)/W1)**4 )
         EX2 = exp( -((X2-WW)/W2)**2 ) * (T/300.d0)**2
         EX3 = exp( -((X3-WW)/W3)**2 ) * (T/300.d0)**1.5
        QO1D = Q1Q2*A1*EX1 + (1.d0-Q1Q2)*A2*EX2 + A3*EX3 + A4
       endif
       if (WW .lt. 220.d0) then
          QO1D = max(0.48d0, 0.48d0 + 0.42d0*(WW-190.d0)/30.d0)
       endif
       if (WW .gt. 340.d0) then
          QO1D = 0.d0
       endif
       XNEW = QO1D

      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine X_O3 (WW,XT,XP,XM,XNEW,
     &     INIT, MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)
c-----------------------------------------------------------------------
c   WW = wavelength (nm) to calc Xsection for
c   XT = temerature (K) for interpolation
c   XP = pressure (hPa) can be used in Stern-Volmer formula if need be
c   XM = air density (#/cm3), ditto
c   XNEW = cross section (cm2) as a function of WW and XT (and XP, XM)
c   INIT = initialization:
c     if INIT.eq.0 .then reads in any tables and sets Xsect name to TITLNEW
c-----------------------------------------------------------------------
      implicit none
      save T1,T2,XT1,XT2,WT1,WT2,NW

      real*8, intent(in) :: WW,XT,XP,XM
      real*8, intent(out) :: XNEW
      integer, intent(inout) :: INIT
      real*8, intent(out) :: MM(3)
      integer, intent(out) :: TTT(3),PPP(3)
      character*1, intent(out) :: ISX,ISP
      character*6, intent(out) :: TITLNEW
      character*90, intent(out) :: TITLTBL(3)

      character*80 FTBL,TABLE,FORMW
      real*8 WT1(999),WT2(999),XT1(999),XT2(999)
      real*8 XXT,XXQ,XXW,FW, T1,T2,TFACT
      integer NW,N,I,IW

      TITLNEW = 'O3    '
      FTBL = 'XO3_JPL11X.dat'

!---include the Temperatures (p, M) that you want to interpolate to and give tables.
      if (INIT .eq. 0) then

          write(6,'(2a)') ' species:',TITLNEW
        open (3, file=FTBL, status='OLD')
          read(3,'(a80)') TABLE
            write(6,'(a80)') TABLE
          read(3,'(5x,3i5)') TTT
          read(3,'(a80)') TITLTBL(1)
!!!!          TITLNEW = TITLTBL(1)
!---this is all table specific formats
          read(3,'(a80)') TABLE
           write(6,'(2a/a)') ' openfile=',FTBL, TABLE
          read(3,'(i4,1x,a)') NW,FORMW
         do N=1,NW
          read(3,FORMW) WT1(N),WT2(N),XT2(N),XT1(N)
!          write(6,'(i5,2f9.3,2f12.6)') N,WT1(N),WT2(N),XT2(N),XT1(N)
         enddo
         close(3)
        INIT = 1
         T2 = 295.d0
         T1 = 218.d0
         MM(:) = 0.d0
         PPP(:) = 0
         ISX = ' '
         ISP = ' '
         close(3)
        INIT = 1

      else

!---interpolate X-section vs. T, but use mean value in the wavelength bins
!    note that WT2(I) = WT1(I+1) -- bins do not miss any wavelengths
        XXT = min(T2, max(T1, XT))
        TFACT = (XXT - T1)/(T2 - T1)
        IW = 1
       do I=1,NW-1
        if (WW .gt. WT2(I)) IW = I+1
       enddo
         XNEW = 1.d-20*( XT1(IW) + TFACT*(XT2(IW)-XT1(IW)) )

      endif

      return
      end
