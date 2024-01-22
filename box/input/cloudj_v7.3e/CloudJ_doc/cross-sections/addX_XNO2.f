!---addFJX_NO2.f  -- with user supplied subroutine that supplies X-section x Q-yield
!---generates fast-JX 18-bin X-sections  revised and updated v73 (mprather,4/2015)
      implicit none
      integer, parameter :: NB_ = 100
      integer, parameter :: NS_ = 40000
      integer, parameter :: NJ_ = 18

      real*8 SRB(15,NJ_)
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN, XBIN
      real*8, dimension(NJ_) :: FFBIN,AABIN,XXBIN
      real*8 W(NS_),F(NS_)
      real*8 W1,W2,WW, XT
      real*8 XNO2, XNO2X, QNO2, QNO2X
      integer IBINJ(NS_)
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,NB77,I778
      integer     TTT(4)
      character*1 ISXP
      character*6 TITLNEW
      character*6 TITLTBL(4)
      Character*20 TITL77

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

! allow for use of WBIN(NB+1) = 778. (all v7.3+) or = 850. (all v7.2-)
      do I778 = 1,2

        if (I778 .eq. 1) then
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
! none needed for NO2 with simple tables set in data
!      call X_NO2 (WW, TT, XNO2,XNO2X)
!      call Q_NO2 (WW, TT, QNO2,QNO2X)

!!!!!!!!!!!!call to user subroutine for Xsections at high res!!!!!!!

!not extrapolated, uses Q-yld at 248K min
        TTT(1) = 220.
        TTT(2) = 294.
        TITLTBL(1) = 'NO2   '
        TITLTBL(2) = 'NO2   '
!std log extrapolation
        TTT(3) = 200.
        TTT(4) = 300.
        TITLTBL(3) = 'NO2ex '
        TITLTBL(4) = 'NO2ex '

       do K = 1,4
        XT = TTT(K)

!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
         XBIN(:) = 0.d0
! primary high-resolution wavelength loop - generate input for pratmo reference J's
        do J = 1,NS_
         I = IBINJ(J)
         if (I .gt. 0) then
          WW = W(J)

           call X_NO2 (WW, XT, XNO2,XNO2X)
           call Q_NO2 (WW, XT, QNO2,QNO2X)

           FBIN(I) = FBIN(I) + F(J)
           ABIN(I) = ABIN(I) + F(J)*XNO2 *QNO2
           XBIN(I) = XBIN(I) + F(J)*XNO2X*QNO2X

         endif
        enddo

        do I=1,NB
         if (FBIN(I) .gt. 0.d0) then
           NB77 = I
           ABIN(I) = ABIN(I)/FBIN(I)
           XBIN(I) = XBIN(I)/FBIN(I)
         endif
        enddo
!---have completed the UCI pratmo std 77-bin data
!          TITL77 = TITLTBL(K)
!          write(6,'(a20,f5.0/(1p,8e10.3))') TITL77,XT,(ABIN(I),I=1,NB77)
!          write(6,*)

!---get Fatst-JX 18 bins from secondary sum of 77-bin pratmo
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
         FFBIN(:) = 0.d0
         AABIN(:) = 0.d0
         XXBIN(:) = 0.d0
        do I=16,NB
         J = IJX(I)
         FFBIN(J) = FFBIN(J) + FBIN(I)
         AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
         XXBIN(J) = XXBIN(J) + FBIN(I)*XBIN(I)
        enddo
        do I=1,15
         do J=1,NJ_
           FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
           AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
           XXBIN(J) = XXBIN(J) + FBIN(I)*XBIN(I)*SRB(I,J)
         enddo
        enddo
        do J=1,NJ_
         if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
         if (FFBIN(J) .gt. 0.d0) XXBIN(J) = XXBIN(J)/FFBIN(J)
        enddo

! save UCI fast-JX v68 data bins for 'FJX_spec.dat'

        TITLNEW = TITLTBL(1)
        if (K .le. 2) then

         if (I778 .eq. 1) then
         write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLTBL(1),' ',TTT(K), AABIN

         else
! fast-J v73 data for 'FJX_spec.dat'
         write(6,*) TITLTBL(1)
         write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'a',AABIN(1:6),TITLNEW
         write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'b',AABIN(7:12),TITLNEW
         write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',AABIN(13:18),TITLNEW

         endif

        else

         if (I778 .eq. 1) then
         write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLTBL(3),' ',TTT(K), XXBIN

         else
! fast-J v73 data for 'FJX_spec.dat'
         write(6,*) TITLTBL(3)
         write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'a',XXBIN(1:6),TITLNEW
         write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'b',XXBIN(7:12),TITLNEW
         write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'c',XXBIN(13:18),TITLNEW

         endif

        endif

       enddo
      enddo

      stop
      end


      subroutine X_NO2 (WW, TT, XNO2,XNO2X)
c---JPL_2010 standard tables.
c---interpolates NO2 cross section vs. wavelength and temperature.
c---    XNO2 is linear, limited to range 220 to 294K
c---    XNO2X does log interpolation and extrapolates beyond 220K & 294K
c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: XNO2,XNO2X

      integer IW,I
      real*8 Q1,Q2,FT,FTX

C----JPL 2010 Table 4C-2.
      real*4, parameter, dimension(43)  :: WNM =
     &[240.964, 243.902, 246.914, 250.000, 253.165, 256.410, 259.740,
     & 263.158, 266.667, 270.270, 273.973, 277.778, 281.690, 285.714,
     & 289.855, 294.118, 298.507, 303.030, 307.692, 312.5, 317.5, 322.5,
     & 327.5, 332.5, 337.5, 342.5, 347.5, 352.5, 357.5, 362.5, 367.5,
     & 372.5, 377.5, 382.5, 387.5, 392.5, 397.5, 402.5, 407.5, 412.5,
     & 417.5, 422.5, 427.5]
      real*4, parameter, dimension(42)  :: X220K =
     &[4.14, 0.961, 0.859, 0.191, 0.496, 0.872, 1.26, 1.77, 2.36, 3.03,
     & 3.94, 5.16, 6.29, 7.72, 9.64, 11.6, 13.2, 16.0, 18.5, 20.8, 24.2,
     & 27.2, 29.4, 33.0, 37.0, 38.6, 43.5, 47.7, 49.2, 53.7, 55.2, 58.4,
     & 58.5, 59.2, 62.4, 58.5, 64.0, 57.0, 61.8, 58.3, 59.3, 56.0]
      real*4, parameter, dimension(42)  :: X294K =
     &[5.77,  2.79,  1.62,  0.998, 1.05,  1.28, 1.58, 2.05, 2.64, 3.24,
     & 4.07, 5.21, 6.23, 7.59, 9.51, 11.5, 13.2, 16.1, 18.8, 21.6, 25.3,
     & 28.7, 31.7, 35.8, 40.2, 41.8, 46.2, 49.7, 50.9, 54.9, 56.1, 59.0,
     & 59.3, 60.1, 63.0, 59.7, 64.4, 58.2, 62.4, 59.1, 59.9, 57.0]

        XNO2 = 0.d0
        XNO2X = 0.d0
      if (WW .lt. WNM(1) .or. WW .gt. WNM(43)) goto 2
      do I=1,42
       if (WW .gt. WNM(I)) IW = I
      enddo
        Q1 = X220K(IW)*1.d-20
        Q2 = X294K(IW)*1.d-20
        FTX = (TT - 220.d0) / (294.d0 - 220.d0)
        FT = max (0.d0, min (1.d0, FTX))
        XNO2 = Q1 + FT*(Q2-Q1)
        XNO2X = XNO2
      if (Q1 .gt. 0.d0) then
        XNO2X = Q1 * (Q2/Q1)**FTX
      endif
ccc        write(6,'(2f8.2,i3,4f8.4)') WW,TT,IW,Q1,Q2,XNO2,XNO2X
    2 continue
      return
      end


      subroutine Q_NO2 (WW, TT, QNO2,QNO2X)
c---JPL_2010 standard tables.
c---interpolates NO2 quantum yield vs. wavelength and temperature.
c---    QNO2 is linear, limited to range Q248 to Q298
c---    QNO2X does log interpolation and extrapolates beyond 248K & 298K
c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: QNO2,QNO2X

      integer IW
      real*8 FW,Q1,Q2,FT,FTX

C----------JPL 2010 Table 4C-3.
      real*4, parameter, dimension(26)  :: WNM =
     & [398., 399., 400., 401., 402., 403., 404., 405., 406., 407.,
     &  408., 409., 410., 411., 412., 413., 414., 415., 416., 417.,
     &  418., 419., 420., 421., 422., 423.]
      real*4, parameter, dimension(26)  :: Q298 =
     & [1.00, 0.95, 0.88, 0.75, 0.62, 0.53, 0.44, 0.37, 0.30, 0.26,
     &  0.22, 0.18, 0.15, 0.13, 0.11, 0.09, 0.08, 0.06, 0.05, 0.04,
     &  0.03, 0.02, 0.02, .015, 0.01, 0.00]
      real*4, parameter, dimension(26)  :: Q248 =
     & [1.00, 0.94, 0.86, 0.69, 0.56, 0.44, 0.34, 0.28, 0.22, 0.18,
     &  0.14, 0.12, 0.10, 0.08, 0.07, 0.06, 0.04, 0.03, 0.02, 0.02,
     &  0.02, 0.01, 0.01, 0.01, 0.01, 0.00]

      IW = (WW - 397.d0)
      IW = max(1, min(25, IW))
      FW = WW - WNM(IW)
      FW = max(0.d0, min(1.d0, FW))
      Q2 = Q298(IW) + FW*(Q298(IW+1)-Q298(IW))
      Q1 = Q248(IW) + FW*(Q248(IW+1)-Q248(IW))
      FTX = (TT - 248.d0) / (298.d0 - 248.d0)
      FT = max (0.d0, min (1.d0, FTX))
      QNO2 = Q1 + FT*(Q2-Q1)
      QNO2X = QNO2
      if (Q1 .gt. 0.d0) then
        QNO2X = Q1 * (Q2/Q1)**FTX
        QNO2X = min (1.d0, QNO2X)
      endif
ccc      write(6,'(2f8.2,i3,f7.3,4f8.4)') WW,TT,IW,FW,Q1,Q2,QNO2,QNO2X
      return
      end

