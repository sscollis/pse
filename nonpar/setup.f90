!==============================================================================
      subroutine setup
!==============================================================================
      use global
      implicit none

!.... local variables

      integer :: i, im, ip, j
      real    :: dm, dmm, dp, dpp, dy, aa, bb
      real    :: h_1, h_2, h_3, gam, del, alp, bet, yc
      real    :: wv1(nx), wv2(ny)

      character*1,parameter :: trans='n'
!==============================================================================

!.... make the streamwise derivative operators
!.... fourth-order-accurate central for a nonuniform mesh

      do i= 1,nx-1
         wv1(i) = x(i+1) - x(i)
      end do

      allocate( cppx(nx), cpx(nx), cx(nx), cmx(nx), cmmx(nx) )
      allocate( c2ppx(nx), c2px(nx), c2x(nx), c2mx(nx), c2mmx(nx) )

      cppx(1) = -wv1(1) / ((wv1(2) + wv1(1)) * wv1(2))
      cpx(1) = (wv1(1) + wv1(2)) / (wv1(1) * wv1(2))
      cx(1) = -1./wv1(1) - 1./(wv1(1) + wv1(2))
      cmx(1) = 0.
      cmmx(1) = 0.
     
      dm = 1. / (wv1(2) + wv1(1))
      c2ppx(1) = 2. * dm / wv1(2)
      c2x(1) = 2. * dm / wv1(1)
      c2px(1) = -c2ppx(1) - c2x(1)
      c2mx(1) = 0.
      c2mmx(1) = 0.

      cppx(nx) = 0.
      cpx(nx) = 0.
      cx(nx) = 1./wv1(nx-1) + 1./(wv1(nx-1) + wv1(nx-2))
      cmx(nx) = -(wv1(nx-1) + wv1(nx-2)) / (wv1(nx-1) * wv1(nx-2))
      cmmx(nx) = wv1(nx-1) / ((wv1(nx-1) + wv1(nx-2)) * wv1(nx-2))

      dm = 1. / (wv1(nx-1) + wv1(nx-2))
      c2x(nx) = 2. * dm / wv1(nx-1)
      c2mmx(nx) = 2. * dm / wv1(nx-2)
      c2mx(nx) = -c2x(nx) - c2mmx(nx)
      c2ppx(nx) = 0.
      c2px(nx) = 0.

      do i= 3,nx-2
         im = max(i-2,1)
         dmm = wv1(im) + wv1(i-1)
         dm = wv1(i-1)
         dp = wv1(i)
         ip = min(i+1,nx-1)
         dpp = wv1(i) + wv1(ip)

         cppx(i) = -dm*dmm*dp/(dpp**4+(-dp+dmm+dm)*dpp**3+((-dmm-dm)*dp &
             +dm*dmm)*dpp**2-dm*dmm*dp*dpp)
         cpx(i) = dm*dmm*dpp/((dp**3+(dmm+dm)*dp**2+dm*dmm*dp)*dpp-dp**4 &
             +(-dmm-dm)*dp**3-dm*dmm*dp**2)
         cmx(i) = -dmm*dp*dpp/(((dm*dmm-dm**2)*dp+dm**2*dmm-dm**3)*dpp &
             +(dm**2*dmm-dm**3)*dp+dm**3*dmm-dm**4)
         cmmx(i) = dm*dp*dpp/(((dmm**2-dm*dmm)*dp+dmm**3-dm*dmm**2)*dpp &
             +(dmm**3-dm*dmm**2)*dp+dmm**4-dm*dmm**3)
         cx(i) = -(cppx(i) + cpx(i) + cmx(i) + cmmx(i))
     
         c2ppx(i) = -2.*((dmm+dm)*dp-dm*dmm)/(dpp**4+(-dp+dmm+dm)*dpp**3 &
             +((-dmm-dm)*dp+dm*dmm)*dpp**2-dm*dmm*dp*dpp)
         c2px(i) = 2.*((dmm+dm)*dpp-dm*dmm)/((dp**3+(dmm+dm)*dp**2 &
             +dm*dmm*dp)*dpp-dp**4+(-dmm-dm)*dp**3-dm*dmm*dp**2)
         c2mx(i) = -2.*((dp-dmm)*dpp-dmm*dp)/(((dm*dmm-dm**2)*dp+dm**2*dmm &
             -dm**3)*dpp+(dm**2*dmm-dm**3)*dp+dm**3*dmm-dm**4)
         c2mmx(i) = 2.*((dp-dm)*dpp-dm*dp)/(((dmm**2-dm*dmm)*dp+dmm**3 &
             -dm*dmm**2)*dpp+(dmm**3-dm*dmm**2)*dp+dmm**4-dm*dmm**3)
         c2x(i) = -(c2ppx(i) + c2px(i) + c2mx(i) + c2mmx(i))
      end do

!.... Otto's ends

!.... i=2

      h_1=x(2)-x(1)
      h_2=x(3)-x(2)
      h_3=x(4)-x(2)
      gam=(h_1*h_3**3+h_1**2*h_3**2)/ &
            ((h_2**2+h_1*h_2)*(h_3**2-h_1*h_3) &
            -(h_2**3-h_1**2*h_2)*(h_3**2+h_1*h_3))
      del=(h_1-(h_2**2+h_1*h_2)*gam) /(h_3**2+h_1*h_3)
      alp=1.0/h_1*(h_2*gam+h_3*del-1.0)
      bet=-(alp+gam+del)
      cppx(2)=del
      cpx(2)=gam
      cx(2)=bet
      cmx(2)=alp
      cmmx(2) = 0.

      gam=2.0/h_2/((h_2+h_1)-(h_3+h_1)*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      del=-1.0/h_3*(gam*h_2*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      alp=1.0/h_1*(h_2*gam+h_3*del)
      bet=-(alp+gam+del)
      c2ppx(2)=del
      c2px(2)=gam
      c2x(2)=bet
      c2mx(2)=alp
      c2mmx(2)=0.0

!.... i=nx-1

      h_1=x(nx)-x(nx-1)
      h_2=x(nx-1)-x(nx-2)
      h_3=x(nx-1)-x(nx-3)
      gam=(h_1*h_3**3+h_1**2*h_3**2)/ &
            ((h_2**2+h_1*h_2)*(h_3**2-h_1*h_3) &
            -(h_2**3-h_1**2*h_2)*(h_3**2+h_1*h_3))
      del=(h_1-(h_2**2+h_1*h_2)*gam)/(h_3**2+h_1*h_3)
      alp=1.0/h_1*(h_2*gam+h_3*del-1.0)
      bet=-(alp+gam+del)
      cmmx(nx-1)=-del
      cmx(nx-1)=-gam
      cx(nx-1)=-bet
      cpx(nx-1)=-alp
      cppx(nx-1)=0.0

      gam=2.0/h_2/((h_2+h_1)-(h_3+h_1)*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      del=-1.0/h_3*(gam*h_2*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      alp=1.0/h_1*(h_2*gam+h_3*del)
      bet=-(alp+gam+del)
      c2mmx(nx-1)=del
      c2mx(nx-1)=gam
      c2x(nx-1)=bet
      c2px(nx-1)=alp
      c2ppx(nx-1)=0.0

!.... make the wall normal derivative operators

      allocate( opy(ny,ny), opyy(ny,ny) )

      if (.true.) then

      call opper (ny, wv2, opy)

      if (imean .ne. 2) then
	if (ystr .eq. 0.) then    ! uniform y mesh
	  do j= 1,ny
            y(j) = .5*ymax * (1. - wv2(j))
            do i= 1,ny
	      opy(i,j) = -2. * opy(i,j) / ymax
            end do
	  end do
	else                      ! stretched y mesh
	  do i= 1,ny
            yc = .5 * (1. - wv2(i))
            y(i) = ymax * ystr * yc / (1. + ystr - yc)
	  end do
	end if
      end if

!.... Use the derivative operator to compute the metrics

#ifdef CRAY
      call sgemv (trans, ny, ny, one, opy, ny, y, 1, zero, wv2, 1)
#else
      call dgemv (trans, ny, ny, one, opy, ny, y, 1, zero, wv2, 1)
#endif
     
      do j= 1,ny
	do i= 1,ny
	  opy(i,j) = opy(i,j) / wv2(i)
	end do
      end do

#ifdef CRAY
      call sgemm (trans, trans, ny, ny, ny, one, opy, ny, opy, ny, &
                  zero, opyy, ny)
#else
      call dgemm (trans, trans, ny, ny, ny, one, opy, ny, opy, ny, &
                  zero, opyy, ny)
#endif

      else

!.... Hacked FD method

	write(*,*) 'Hacked mesh in setup.f90'

	do j = 1, ny
	  y(j) = real(j-1) * ymax / real(ny-1)
	end do
	  
	aa = ymax * ystr / ( ymax - 2.0 * ystr )
	bb = 1.0 + aa / ymax
	do j = 1, ny
	  y(j) = aa * real(j-1)/real(ny-1) / ( bb - real(j-1)/real(ny-1) )
	  wv2(j) = 1.0/((bb - real(j-1)/real(ny-1))**2 / (aa * bb))
	end do

!	call opperfd4( ny, y, opy, opyy )
	call opperfd2( ny, y, opy, opyy )
	
      endif

!.... save the metric of the mapping

      allocate( dydeta(ny) )
      dydeta = wv2

!.... setup for Chebyschev integration and interpolation

      allocate( opi(ny), yint(nint,ny) )
      call chebyinit( ny, y, opi, nint, yint )

!.... print out opi for later use

      open(91,file='opi.dat')
      write(91,*) ny
      write(91,*) (opi(i),i=1,ny)
      close(91)

      open(92,file='gridy.dat')
      write(92,*) ny
      write(92,*) (y(i),i=1,ny)
      close(92) 



      return
      end subroutine setup

