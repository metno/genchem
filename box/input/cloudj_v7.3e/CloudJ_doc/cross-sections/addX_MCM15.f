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

      ! WBIN holds cloud-j bin edge
      ! SRB and IJX holds bin assignments
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

      ! W(J) holds solar spectrum wavelength
      ! F(J) holds solar photon count per wavelength 
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
          ! write(*,*) I, W1,W2, J1,J2,W(J1),W(J2)
      enddo
!!!!this binning does not interpolate and is OK for large bins
!     it has 7% error in the very short wavel S-R bins of pratmo.

! Wim's code down here

!!!!!!!!!!!!!!!!!!call to user subroutine for Xsections at high res!!!!!!!
        K=1
        TTT(K) = 298
        XT = TTT(K)

!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
! primary high-resolution wavelength loop - generate input for pratmo reference J's
        do J = 1,NS_

         I = IBINJ(J) ! index matching 77-pratmo bin

         if (I .gt. 0) then

          call BIACET (W(J), XNEW, TITLNEW)

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

       TITLNEW = 'MCM15 '

! save UCI fast-JX v68 data bins for 'FJX_spec.dat'
    !    write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
    ! &     TITLNEW,' ',TTT(K), AABIN
    !    write(6,*)

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

! c>>>>>>>>>>>>>>>>>>>>>>>added Xsection<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!       subroutine Y_PAR (WW, TT, XNEW, TITLNEW)
! c---Photosynthetically Active Radiation: action spectrum (quantum): Y-PAR
! c---traced from:
! c      McCree, Keith J. (1972a). "The action spectrum, absorptance and
! c        quantum yield of photosynthesis in crop plants"
! c        Agricultural and Forest Meteorology 9:191-216.
! c      McCree, Keith J. (1972b). Agric. & Forest Meteorology 10:443-453.
! c---PAR in PAR is normally quantified as �mol photons/m2/s =? �E/m2/s
! c        photosynthetic photon flux (area) density, or PPFD.

!       implicit none
!       real*8, intent(in) :: WW, TT
!       real*8, intent(out):: XNEW
!       character*6, intent(out):: TITLNEW

!       integer IWW
!       real*8 FWW,WWI,XX,BB

!       character*6, parameter:: JNEW = 'Y-PAR '
!       real*8, dimension(18), parameter :: W = [325.d0,350.d0,375.d0,
!      &     400.d0,425.d0,450.d0,475.d0,500.d0,525.d0,550.d0,575.d0,
!      &     600.d0,625.d0,650.d0,675.d0,700.d0,725.d0,750.d0]
!       real*8, dimension(18), parameter :: Y = [0.d0,15.d-2,45.d-2,
!      &     64.d-2,78.d-2,75.d-2,68.d-2,70.d-2,74.d-2,88.d-2,95.d-2,
!      &     100.d-2,100.d-2,94.d-2,92.d-2,43.d-2,4.d-2,0.d0]

!       TITLNEW = JNEW
!       WWI = 0.04d0*(WW - 300.d0)
!         IWW = WWI
!         IWW = max( 1, min( 17, IWW))
!         FWW = WWI - float(IWW)
!         FWW = max( 0.d0, min( 1.d0, FWW))
!       XNEW = Y(IWW) + (Y(IWW+1)-Y(IWW))*FWW
!       return
!       end


c>>>>>>>>>>>>>>>>>>>>>>>added Xsection<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine BIACET (WW, XNEW, TITLNEW)
c routine that returns cross-section x quantum-yield at a given input
c wavelength (WW).
c XNEW = cross section (cm2) as a function of WW.
c quantum-yield for BIACET is assumed 0.158 from 290 to 460 nm (Plum, 1983).

      implicit none
      real*8, intent(in) :: WW
      real*8, intent(out):: XNEW
      character*6, intent(out):: TITLNEW

      INTEGER :: I

      real*8  :: FB(10000), WB(10000), INPWW(10)
      real*8  :: IND(10)
      integer :: NBIACET
      integer :: IWW
      real*8  :: QYIELD
      real*8 FWW,WWI,XX,BB
      logical, save :: first_call =.true.

      ! min and max wavelength data available for BIACET
      real*8 :: BIACET_MIN=202, BIACET_MAX=364 ! markers for file length

      character*6, parameter:: JNEW = 'BIACET '

      ! set title for output
      TITLNEW = JNEW
      INPWW(1) = WW

      IF (WW .lt. 202. .or. WW .gt. 364.) THEN
        ! zero outside of this range
        XNEW = 0.
      ELSE

       open (2, file='MCM15/mcm15_in.dat', status='OLD')
         WB(:) = 0.d0
         FB(:) = 0.d0
         READ(2,'(i5)') NBIACET
         READ(2,*)     QYIELD

         ! wavelength
         DO I=1,NBIACET
           READ(2,*) WB(I)
         ENDDO
        
         ! cross-section (10**20 cm2 molecule-1)
         DO I=1,NBIACET
           READ(2,*) FB(I)
         ENDDO

c    Input, integer ND, the number of data points.
c    ND must be at least 1.
c
c    Input, double precision XD(ND), the data points.
c    Input, double precision YD(ND), the data values.
c    Input, integer NI, the number of interpolation points.
c    Input, double precision XI(NI), the interpolation points.
c    Output, double precision YI(NI), the interpolated values.

        ! find nearest neighbour input WW to BIACET wavelength
         call nearest_interp_1d ( NBIACET, WB, WB, 1, INPWW, IND )
        
         XNEW = FB(INT(IND(1)) - INT(BIACET_MIN) + 1) * 1e-21 * QYIELD
        !    if (INPWW(1) .gt. 225.) then
        !     if (first_call) then
        !     write(*,*) 'dip dupr', INPWW(1), IND(1), QYIELD
        !     if (INPWW(1) .gt.275) first_call = .false.
        !     endif
        ! end if

       close (2)

       ENDIF

       return
       end

      subroutine nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
c
c  Discussion:
c
c    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
c    constant function which interpolates the data (XD(I),YD(I)) for I = 1
c    to ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c    ND must be at least 1.
c
c    Input, double precision XD(ND), the data points.
c
c    Input, double precision YD(ND), the data values.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision YI(NI), the interpolated values.
c
      implicit none

      integer nd
      integer ni

      double precision d
      double precision d2
      integer i
      integer j
      integer k
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)

      do i = 1, ni

        k = 1
        d = abs ( xi(i) - xd(k) )

        do j = 2, nd

          d2 = abs ( xi(i) - xd(j) )

          if ( d2 .lt. d ) then
            k = j
            d = d2
          end if

        end do

        yi(i) = yd(k)

      end do

      return
      end