program inprod

!.... computes a L2 inner product using trapezoidal rule

      real, allocatable :: y(:)
      real :: yl, ur, ui, vr, vi, wr, wi, pr, pi
      complex, allocatable :: u(:,:), q(:,:)
      complex :: prod
      
      open(10,file='time.out')
      ny = 0
      do while (.true.)
	read(10,*,end=10) yl, ur, ui, vr, vi, wr, wi, pr, pi
	ny = ny + 1
      end do
10    continue
      write(*,*) 'ny = ',ny
      close(10)

      allocate( y(ny), u(ny,4), q(ny,4) )

      open(10,file='time.out')
      do j = 1, ny
	read(10,*) y(j), ur, ui, vr, vi, wr, wi, pr, pi
	u(j,1) = cmplx(ur,ui); u(j,2) = cmplx(vr,vi)
	u(j,3) = cmplx(wr,wi); u(j,4) = cmplx(pr,pi)
      end do
      close(10)

      open(10,file='tadjoint.out')
      do j = 1, ny
	read(10,*) yl, ur, ui, vr, vi, wr, wi, pr, pi
	q(j,1) = cmplx(ur,ui); q(j,2) = cmplx(vr,vi)
	q(j,3) = cmplx(wr,wi); q(j,4) = cmplx(pr,pi)
!	q(j,:) = conjg( q(j,:) )
      end do
      close(10)

!.... form the inner product
      
      prod = 0.0
      do j = 1, ny-1
	prod = prod + ( u(j,1)*q(j,1) + u(j,2)*q(j,2) + u(j,3)*q(j,3) + &
                        u(j,4)*q(j,4) + u(j+1,1)*q(j+1,1) + &
	                u(j+1,2)*q(j+1,2) + u(j+1,3)*q(j+1,3) + &
			u(j+1,4)*q(j+1,4) ) * (y(j+1) - y(j)) * 0.5
      end do

      write(*,*) prod

      stop
end program inprod
