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
      character*20 TITLTBL(3)
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

!!!!!!!!!!!!!!!!!!call to user subroutine for Xsections at high res!!!!!!!
        K=1
        TTT(K) = 300
        XT = TTT(K)

!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
! primary high-resolution wavelength loop - generate input for pratmo reference J's
        do J = 1,NS_
         I = IBINJ(J)
         if (I .gt. 0) then

         call Y_PAR (W(J), XT, XNEW, TITLNEW)

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

        TITLNEW = 'Y-Par '

! save UCI fast-JX v68 data bins for 'FJX_spec.dat'
!        write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
!     &     TITLNEW,' ',TTT(K), AABIN
!        write(6,*)

! fast-J v73 data for 'FJX_spec.dat'
        write(6,'(a6)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'a',AABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'b',AABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',AABIN(13:18),TITLNEW

      write(6,*) ' Last wavelength: ',WBIN(NB+1)

      stop
      end

c>>>>>>>>>>>>>>>>>>>>>>>added Xsection<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine Y_PAR (WW, TT, XNEW, TITLNEW)
c---Photosynthetically Active Radiation: action spectrum (quantum): Y-PAR
c---traced from:
c      McCree, Keith J. (1972a). "The action spectrum, absorptance and
c        quantum yield of photosynthesis in crop plants"
c        Agricultural and Forest Meteorology 9:191-216.
c      McCree, Keith J. (1972b). Agric. & Forest Meteorology 10:443-453.
c---PAR in PAR is normally quantified as µmol photons/m2/s =? µE/m2/s
c        photosynthetic photon flux (area) density, or PPFD.

      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: XNEW
      character*6, intent(out):: TITLNEW

      integer IWW
      real*8 FWW,WWI,XX,BB

      character*6, parameter:: JNEW = 'Y-PAR '
      real*8, dimension(18), parameter :: W = [325.d0,350.d0,375.d0,
     &     400.d0,425.d0,450.d0,475.d0,500.d0,525.d0,550.d0,575.d0,
     &     600.d0,625.d0,650.d0,675.d0,700.d0,725.d0,750.d0]
      real*8, dimension(18), parameter :: Y = [0.d0,15.d-2,45.d-2,
     &     64.d-2,78.d-2,75.d-2,68.d-2,70.d-2,74.d-2,88.d-2,95.d-2,
     &     100.d-2,100.d-2,94.d-2,92.d-2,43.d-2,4.d-2,0.d0]

      TITLNEW = JNEW
      WWI = 0.04d0*(WW - 300.d0)
        IWW = WWI
        IWW = max( 1, min( 17, IWW))
        FWW = WWI - float(IWW)
        FWW = max( 0.d0, min( 1.d0, FWW))
      XNEW = Y(IWW) + (Y(IWW+1)-Y(IWW))*FWW
      return
      end
