!==============================================================================
      subroutine chebyinit( ny, y, opi, nint )
!==============================================================================
      implicit none

      integer :: ny, nint, i, j, ii
      real    :: y(ny), opi(ny), yint(nint,ny)
      real    :: con, pin, pi, w1(ny,ny), f(ny)

!.... Prepare for Chebyshev interpolation

      con = 2. / (ny-1)
      pin = acos(-1.) / (ny-1)
      do j= 1,ny
	do i= 1,ny
	  w1(i,j) = con * cos((i-1) * (j-1) * pin)
	end do
      end do
      do i= 1,ny
	w1(i,1) = .5 * w1(i,1)
	w1(i,ny) = .5 * w1(i,ny)
	w1(1,i) = .5 * w1 (1,i)
	w1(ny,i) = .5 * w1(ny,i)
      end do
      con = acos(-1.) / (nint-1)
      do i= 1,nint
	do j= 1,ny
	  f(j) = cos((j-1) * (i-1) * con)
	end do
	do j= 1,ny
	  yint(i,j) = 0.
	  do ii= 1,ny
	    yint(i,j) = yint(i,j) + f(ii) * w1(ii,j)
	  end do
	end do
      end do

!.... prepare for Chebyshev integration

      pi = acos(-1.)
      f(1) = -2.
      f(2) = 0.
      do i= 3,ny
	f(i) = .5 * ((cos(i*pi) - 1.) / i - (cos((i-2) * pi) - 1.) / (i-2))
      end do
      do i= 1,ny
	opi(i) = 0.
	do ii= 1,ny
	  opi(i) = opi(i) + f(ii) * w1(ii,i)
	end do
      end do
      call opper (ny, f, w1)
      do i= 1,ny
	f(i) = 0.
	do ii= 1,ny
	  f(i) = f(i) + w1(i,ii) * y(ii)
	end do
      end do
      do i= 1,ny
	opi(i) = opi(i) * f(i)
      end do

      return
      end subroutine chebyinit
