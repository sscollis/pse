!============================================================================================
! here I do not use derivatives
!============================================================================================
!                 ie  is  x  y  ulp uap  ub   opi  jx  jxx  1/rei  
!============================================================================================
subroutine jxnorm(nx, ny, x, y, u, uadj, ub, ,opi, jx, jxx,  R)  
!============================================================================================
  complex, parameter    ::  iota   = (0.0,1.0)
  complex               ::  u(ny,4,nx), uadj(ny,4,nx)
  complex               ::  a(ny), b(ny), c1(ny), dotp(ny),jx(nx),jxx(nx)
  real                  ::  ub(ny,nx), y(ny), x(nx), opi(ny), rex(nx), R  
  integer               ::  ny, nx, i, j, k
!_____________________________________________________________________________________________

  open(10, file='jx.out' )
!_____________________________________________________________________________________________

!start major do loop in i

  do i=1,nx     
     rex(i) = sqrt(x(i)*R);
!______________________________________________________________________________________________

     dotp=cmplx(0.0,0.0)
     do j=1,ny
        dotp(j)=u(j,1,i)*uadj(j,1,i)+u(j,2,i)*uadj(j,2,i)+u(j,3,i)*uadj(j,3,i)
        c1(j)=dotp(j)*Ubx(j,i)+u(j,1,i)*uadj(j,4,i)+uadj(j,1,i)*u(j,4,i)
     end do
!_______________________________________________________________________________________________

     jx(i)=cmplx(0.0,0.0)
     do j=1,ny
        jx(i)=jx(i)+opi(j)*c1(j)
     end do
     write(10,'(I4,7E16.7)') i,x(i),rex(i),real(jx(i)),aimag(jx(i)),abs(jx(i))
  end do
  close(10)

!===============================================================================================
end subroutine jxnorm
!===============================================================================================





!===============================================================================================
      subroutine chebyinit( ny, y, opi, nint, yint )
!===============================================================================================
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
!======================================================================================================
      end subroutine chebyinit
!======================================================================================================
