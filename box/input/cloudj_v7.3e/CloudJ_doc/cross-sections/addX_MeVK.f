!---addX_FJX_sss.f  -- with user supplied subroutine that supplies X-section x Q-yield
!---generates fast-JX 18-bin X-sections  revised and updated v73 (mprather,4/2015)
      implicit none
      integer, parameter :: NB_ = 100
      integer, parameter :: NS_ = 40000
      integer, parameter :: NJ_ = 18

      real*8 SRB(15,NJ_)
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN
      real*8, dimension(NJ_) :: FFBIN,AABIN
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,K1,K2,NB77
      integer INIT
      real*8 W(NS_),F(NS_)
      integer IBINJ(NS_)
      real*8 W1,W2, XT,XP,XM, MM(3), WW,XNEW
      character*6 TITLNEW
      character*1 ISX, ISP, ISXP
      character*90 TITLTBL(3)
      Character*20 TITL77
      integer     TTT(3),PPP(3)

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
!          write(6,'(i5,2f9.3,2i9,2f9.3)') I, W1,W2, J1,J2,W(J1),W(J2)
      enddo
!!!!this binning does not interpolate and is OK for large bins
!     it has 7% error in the very short wavel S-R bins of pratmo.

!!!!!!!!!!!!!!!!!!initialization call to user subroutine!!!!!!!!!!!!!!
      INIT = 0

      call X_MeVK (W(1),XT,XP,XM, XNEW,
     &     INIT, MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)

            ISXP = ' '
       if (ISX .eq. 'x') then
            ISXP = ISX
       else if (ISP .eq. 'p') then
            ISXP = ISP
       endif

!!!!!!!!!!!!!!!!!!call to user subroutine for Xsections at high res!!!!!!!
!  major temperature-density loop for X-sections
      do K = 1,3
       if (TTT(K) .gt. 0) then
        XT = TTT(K)
        XP = PPP(K)
        XM = MM(K)

!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
! primary high-resolution wavelength loop - generate input for pratmo reference J's
        do J = 1,NS_
         I = IBINJ(J)
         if (I .gt. 0) then

      call X_MeVK (W(J),XT,XP,XM, XNEW,
     &      INIT, MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)

           FBIN(I) = FBIN(I) + F(J)
           ABIN(I) = ABIN(I) + F(J)*XNEW
         endif
        enddo
        do I=1,NB
         if (FBIN(I) .gt. 0.d0) then
           NB77 = I
           ABIN(I) = ABIN(I)/FBIN(I)
         endif
        enddo

!---write out UCI std 77-bin data
!          write(6,'(a6,14x,f5.0/(1p,8e10.3))')
!     &       TITLNEW,XT,(ABIN(I),I=1,NB77)
!          write(6,*)

!---secondary sum 77-bin pratmo ==> 18-bin fast-JX
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
       do I=16,NB
        J = IJX(I)
        FFBIN(J) = FFBIN(J) + FBIN(I)
        AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
       enddo
       do I=1,15
        do J=1,NJ_
          FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
          AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
        enddo
       enddo
       do J=1,NJ_
        if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
       enddo

! save UCI fast-JX v68 data bins for 'FJX_spec.dat'
        if (ISP .eq. 'p') then
        write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLNEW,ISP,PPP(K), AABIN
        else
        write(6,'(a6,1x,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLNEW,   TTT(K), AABIN
        endif
        write(6,*)

! fast-J v73 data for 'FJX_spec.dat'
        write(6,'(a90)') TITLTBL(K)
        if (ISP .eq. 'p') then
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISP,PPP(K),'a',AABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',PPP(K),'b',AABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',PPP(K),'c',AABIN(13:18),TITLNEW
        else
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'a',AABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',AABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',AABIN(13:18),TITLNEW
        endif
        write(6,*)

       endif
      enddo

      stop
      end


c-------------sample subroutine for fast-JX Xsection generation---------
c-----------------------------------------------------------------------
      subroutine X_MeVK (WW,XT,XP,XM,XNEW,
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
      save W,XW,QW,QT3,QQW,NW,NB

      real*8, intent(in) :: WW,XT,XP,XM
      real*8, intent(out) :: XNEW
      integer, intent(inout) :: INIT
      real*8, intent(out) :: MM(3)
      integer, intent(out) :: TTT(3),PPP(3)
      character*1, intent(out) :: ISX,ISP
      character*6, intent(out) :: TITLNEW
      character*90, intent(out) :: TITLTBL(3)

      character*80 FTBL,TABLE,FORMW,FORMB
      real*8 W(999),XW(999),QW(999,3), QT3(3),QQW(3)
      real*8 XXT,XXQ,XXW,FW, QFACT
      integer NW,NB,  N, I,IW,IT

      TITLNEW = 'MeVK  '
      FTBL = 'XMeVK_JPL11X.dat'

!---include the Temperatures (p, M) that you want to interpolate to and give tables.
!---general format for generating all Xsections
      if (INIT .eq. 0) then

          write(6,'(2a)') ' species:',TITLNEW
        open (3, file=FTBL, status='OLD')
          read(3,'(a80)') TABLE
            write(6,'(a80)') TABLE
          read(3,'(a1)') ISX
          read(3,'(a1)') ISP
            write(6,'(4a4)') 'ISX=',ISX,'ISP=',ISP
         if (ISP .eq. 'p') then
          read(3,'(5x,3i5)') PPP
          read(3,'(5x,3f5.0)') MM
             MM(:) = MM(:)*1.e19
         endif
          read(3,'(5x,3i5)') TTT
         do I = 1,3
           if (TTT(I) .gt. 0) then
            read(3,'(a90)') TITLTBL(I)
           endif
         enddo
          TITLNEW = TITLTBL(1)
!---this is all table specific formats
          read(3,'(a80)') TABLE
           write(6,'(2a/a)') ' openfile=',FTBL, TABLE
          read(3,*)
          read(3,*)
          read(3,'(14x,3f8.0)') (QT3(I),I=1,3)
             write(6,'(a,3f8.0)') ' Temperature:',(QT3(I),I=1,3)
          read(3,'(i4,1x,a)') NW,FORMW
!---for MeVK we have Xsections and QuantumYlds (at 3 T-p's)
        do N = 1,NW
          read(3,FORMW) W(N),XW(N),(QW(N,I),I=1,3)
!         write(6,'(I3,5F12.5)') int(W(N)),XW(N),(QW(N,I),I=1,3)
        enddo
        close(3)
        INIT = 1

      else

!--Note that Xsects need to be interpolated for wavelength, but keyed to
!  Sturm-Volmer and atmos lapse rate.- evaluted for 3 specific values in table
!  can allow for interp if need be - Locate T-p interp numbers for all wavels
          XXT = min(QT3(3), max(QT3(1), XT))
        if (XXT .gt. QT3(2)) then
          IT = 2
        else
          IT = 1
        endif
          QFACT = (XXT - QT3(IT))/(QT3(IT+1) - QT3(IT))
c---interpolate X-section vs. Wavelength
        IW = 1
       do I=2,NW-1
        if (WW .gt. W(I)) IW = I
       enddo
         FW = (WW - W(Iw))/(W(IW+1) - W(IW))
         FW = min(1.d0, max(0.d0, FW))
        XXW = XW(IW) + FW*(XW(IW+1)-XW(IW))
        QQW(1) = QW(IW,1) + FW*(QW(IW+1,1)-QW(IW,1))
        QQW(2) = QW(IW,2) + FW*(QW(IW+1,2)-QW(IW,2))
        QQW(3) = QW(IW,3) + FW*(QW(IW+1,3)-QW(IW,3))
        XXQ = XXW * (QQW(IT) + QFACT*(QQW(IT+1)-QQW(IT)) )
        XNEW = XXQ * 1.d-20

      endif

      return
      end
