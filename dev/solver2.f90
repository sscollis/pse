!==============================================================================
      subroutine solver
!
!     Solve the linear PSE using Chebyshev collocation
!
!     S. Scott Collis
!
!     8-20-98:  Revised to second order using midpoint rule
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, iter, il, idof
      integer :: l, ldof, l0, m, mdof, m0, info, job=0, icount=0
      integer, allocatable :: ipvt(:)

      real :: rei

      complex :: dalpha(nz,nt)
      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)
      complex :: ux(ny),  vx(ny),  wx(ny),  px(ny)
      complex :: kex
      real    :: ke(nx)

!.... midpoint rule

      complex :: alphal
      real    :: ubl(ny),  vbl(ny),  wbl(ny)
      real    :: ubxl(ny), vbxl(ny), wbxl(ny)
      real    :: ubyl(ny), vbyl(ny), wbyl(ny)
      complex :: uy(ny,ndof), uyy(ny,ndof)

      real, allocatable :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
	                   E(:,:,:), Ex(:,:,:)

      complex, allocatable :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:)

      real :: c1(ny), c2(ny)
      
      real, parameter :: fact=0.5
!==============================================================================
      rei = one / re

      allocate( u(nx,ny,nz,nt,ndof) )

      u = czero

      allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
	        C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
	        Ex(ny,ndof,ndof) )

      allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
      allocate( r(ny*ndof), ipvt(ny*ndof) )

!.... only consider only one mode for linear analysis

      i = is; n = 1; k = 1

!.... initialize the field

      do j = 1, ny
	u(i,j,k,n,1) = ubci(j,k,n)
	u(i,j,k,n,2) = vbci(j,k,n)
	u(i,j,k,n,3) = wbci(j,k,n)
	u(i,j,k,n,4) = pbci(j,k,n)
      end do

      ul = u(i,:,k,n,1)
      vl = u(i,:,k,n,2)
      wl = u(i,:,k,n,3)

      ke(i) = zero
      do j = 1, ny-1
	ke(i) = ke(i) + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                                conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j)-y(j-1))
      end do

!.... some useful debugging data files

      if (.false.) then
      do j = 1, ny
	write(50,"(4(1pe13.6,1x))") y(j),ub(j,i),vb(j,i),wb(j,i)
	write(51,"(4(1pe13.6,1x))") y(j),ubx(j,i),vbx(j,i),wbx(j,i)
	write(52,"(4(1pe13.6,1x))") y(j),uby(j,i),vby(j,i),wby(j,i)
	write(53,"(10(1pe13.6,1x))") y(j), &
                                    real(u(i,j,k,n,1)), aimag(u(i,j,k,n,1)), &
                                    real(u(i,j,k,n,2)), aimag(u(i,j,k,n,2)), &
                                    real(u(i,j,k,n,3)), aimag(u(i,j,k,n,3)), &
                                    real(u(i,j,k,n,4)), aimag(u(i,j,k,n,4))
	write(54,"(10(1pe13.6,1x))") y(j), &
                                    abs(u(i,j,k,n,1)), abs(u(i,j,k,n,2)), &
                                    abs(u(i,j,k,n,3)), abs(u(i,j,k,n,4))
      end do
      end if

!.... main marching loop

      open(9,file='pse.dat')
      do i = is+1, ie

	cpu2 = second()
	write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
	cpu = cpu2

	alpha(i,k,n) = alpha(i-1,k,n)

	c1(:) = one/(one + pt5*(cur(i-1)+cur(i))*y(:))
	c2(:) = pt5*(cur(i-1)+cur(i)) * c1(:)

!	c1(:) = one/(one + cur(i)*y(:))
!	c2(:) = cur(i) * c1(:)

	u(i,:,k,n,:) = u(i-1,:,k,n,:)

	ul = u(i-1,:,k,n,1)
	vl = u(i-1,:,k,n,2)
	wl = u(i-1,:,k,n,3)
	pl = u(i-1,:,k,n,4)

!.... compute the mode amplitude

	amp(i,k,n) = czero
	do il = is, i-1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
		       (x(il+1)-x(il))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )

	ke(i) = zero
	do j = 1, ny-1
	ke(i) = ke(i) + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                                conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
        end do

	write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
                           abs(amp(i,k,n)), ke(i), zero

	do iter = 1, niter ! loop on alpha
	
!.... compute the mode amplitude

	amp(i,k,n) = czero
	do il = is, i-1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
		       (x(il+1)-x(il))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )

!.... compute the streamwise derivative of alpha

	dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )

!.... Bertollotti sets this to zero (no advantage)

!	dalpha(k,n) = zero

!.... Evaluate the matrices.  For second order everything 
!.... should be at the midpoint of the interval

	alphal = pt5*(alpha(i-1,k,n) + alpha(i,k,n))
	ubl    = pt5*(ub(:,i-1) + ub(:,i))
	vbl    = pt5*(vb(:,i-1) + vb(:,i))
	wbl    = pt5*(wb(:,i-1) + wb(:,i))
	ubxl   = pt5*(ubx(:,i-1) + ubx(:,i))
	vbxl   = pt5*(vbx(:,i-1) + vbx(:,i))
	wbxl   = pt5*(wbx(:,i-1) + wbx(:,i))
	ubyl   = pt5*(uby(:,i-1) + uby(:,i))
	vbyl   = pt5*(vby(:,i-1) + vby(:,i))
	wbyl   = pt5*(wby(:,i-1) + wby(:,i))

!	alphal = alpha(i,k,n)
!	ubl    = ub(:,i)
!	vbl    = vb(:,i)
!	wbl    = wb(:,i)
!	ubxl   = ubx(:,i)
!	vbxl   = vbx(:,i)
!	wbxl   = wbx(:,i)
!	ubyl   = uby(:,i)
!	vbyl   = vby(:,i)
!	wbyl   = wby(:,i)

	G=zero; A=zero; B=zero; C=zero; D=zero; E=zero; Ex=zero
	Gb=czero; Ab=czero

	G(:,1,1) = one
	G(:,2,2) = one
	G(:,3,3) = one

	A(:,1,1) = c1(:) * ubl
	A(:,1,2) = -rei * two * c1(:) * c2(:)
	A(:,1,4) = c1(:)
	A(:,2,1) = rei * two * c1(:) * c2(:)
	A(:,2,2) = c1(:) * ubl
	A(:,3,3) = c1(:) * ubl
	A(:,4,1) = c1(:)

	B(:,1,1) = vbl - rei * c2(:)
	B(:,2,2) = vbl - rei * c2(:)
	B(:,2,4) = one
	B(:,3,3) = vbl - rei * c2(:)
	B(:,4,2) = one

	C(:,1,1) = wbl
	C(:,2,2) = wbl
	C(:,3,3) = wbl
	C(:,3,4) = one
	C(:,4,3) = one

	D(:,1,1) = c1(:) * ubxl + c2(:) * vbl - rei * c2(:)**2
	D(:,1,2) = ubyl + c2(:) * ubl
	D(:,2,1) = c1(:) * vbxl - two * c2(:) * ubl
	D(:,2,2) = vbyl + rei * c2(:)**2
	D(:,3,1) = c1(:) * wbxl
	D(:,3,2) = wbyl
	D(:,4,2) = c2(:)

	Ex(:,1,1) = c1(:)**2 * rei
	Ex(:,2,2) = c1(:)**2 * rei
	Ex(:,3,3) = c1(:)**2 * rei

	E(:,1,1) = rei
	E(:,2,2) = rei
	E(:,3,3) = rei

	Gb = -iota * kt(n) * omega * G + iota * alphal * A + &
	      iota * kz(k) * beta * C + D - &
	      Ex * ( iota * dalpha(k,n) - alphal**2 ) + &
	      kz(k)**2 * beta**2 * E

!.... fix the streamwise pressure gradient

	if (ipfix.eq.1 .or. (ipfix.eq.2 .and. i.eq.2) ) A(:,:,4) = zero

	Ab = A - two * iota * alphal * Ex

!.... Build the LHS matrix

	A0 = czero
	do ldof = 1, ndof
	  do mdof = 1, ndof
	    do l = 1, ny
	      l0 = (l-1)*ndof
	      do m = 1, ny
		m0 = (m-1)*ndof
		A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
		                      fact * ( B(l,ldof,mdof) * opy(l,m) - &
		                               E(l,ldof,mdof) * opyy(l,m) )
	      end do
	      m0 = (l-1)*ndof
	      A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
                                    fact * Gb(l,ldof,mdof) + &
	                            Ab(l,ldof,mdof) / (x(i)-x(i-1))
	    end do
	  end do
	end do

!.... Apply the boundary conditions to the LHS

	l = 1
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(l0+1,l0+1) = cone

	A0(l0+2,:) = czero
	A0(l0+2,l0+2) = cone

	A0(l0+3,:) = czero
	A0(l0+3,l0+3) = cone

	A0(l0+4,:) = czero
	if (ipwbc .eq. 0) then
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+4) = opy(l,m)
	  end do
	else if (ipwbc .eq. 1) then
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+2) = -rei * c2(l) * opy(l,m) - rei * opyy(l,m)
	    A0(l0+4,m0+4) = opy(l,m)
	  end do
	else
	  call error('solver$','Illegal value of ipwbc$')
	end if

	l = ny
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(l0+1,l0+1) = cone

	if (ivbc .eq. 0) then
	  A0(l0+2,:) = czero
	  A0(l0+2,l0+2) = cone
	else if (ivbc .eq. 1) then
	  A0(l0+2,:) = czero
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+2,m0+2) = opy(l,m)
	  end do
	else
	  call error('solver$','Illegal value of ivbc$')
	end if

	A0(l0+3,:) = czero
	A0(l0+3,l0+3) = cone

	A0(l0+4,:) = czero
	if (ipbc .eq. 0) then
	  A0(l0+4,l0+4) = cone	
	else if (ipbc .eq. 1) then
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+4) = opy(l,m)
	  end do
	else
	  call error('solver$','Illegal value of ipbc$')
	end if

!.... Build the RHS

	ul = u(i-1,:,k,n,1)
	vl = u(i-1,:,k,n,2)
	wl = u(i-1,:,k,n,3)
	pl = u(i-1,:,k,n,4)

	uy = zero
	uyy = zero
	do m = 1, ny
	  do l = 1, ny
	    uy(l,1) = uy(l,1) + opy(l,m) * ul(m)
	    uy(l,2) = uy(l,2) + opy(l,m) * vl(m)
	    uy(l,3) = uy(l,3) + opy(l,m) * wl(m)
	    uy(l,4) = uy(l,4) + opy(l,m) * pl(m)

	    uyy(l,1) = uyy(l,1) + opyy(l,m) * ul(m)
	    uyy(l,2) = uyy(l,2) + opyy(l,m) * vl(m)
	    uyy(l,3) = uyy(l,3) + opyy(l,m) * wl(m)
	    uyy(l,4) = uyy(l,4) + opyy(l,m) * pl(m)
	  end do
	end do

!	do l = 1, ny
!	  write(44,"(8(1pe13.6,1x))") y(l), ul(l), uy(l,1), uyy(l,1)
!	end do

!	r = czero
!	do ldof = 1, ndof
!	  do l = 1, ny
!	    l0 = (l-1)*ndof
!	    do mdof = 1, ndof
!	      r(l0+ldof) = r(l0+ldof) + &
! 	               (  B(l,ldof,mdof) *  uy(l,mdof)          - &
!	                  E(l,ldof,mdof) *  uyy(l,mdof)         + &
!                         Gb(l,ldof,mdof) *  u(i-1,l,k,n,mdof) )
!	    end do
!	  end do
!	end do
!
!	do l = 1, ny
!	  l0 = (l-1)*ndof
!	  write(45,"(10(1pe13.6,1x))") y(l), r(l0+1), r(l0+2), r(l0+3), r(l0+4)
!	end do

	r = czero
	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do mdof = 1, ndof
	      r(l0+ldof) = r(l0+ldof) + &
	  (fact-1.0) * (  B(l,ldof,mdof) *  uy(l,mdof)          - &
	                  E(l,ldof,mdof) *  uyy(l,mdof)         + &
			  Gb(l,ldof,mdof) *  u(i-1,l,k,n,mdof) ) + &
		Ab(l,ldof,mdof) / (x(i)-x(i-1)) * u(i-1,l,k,n,mdof)
	    end do
	  end do
	end do

!.... Enforce the boundary conditions on the RHS

	l = 1
	l0 = (l-1)*ndof
	r(l0+1) = czero
	r(l0+2) = czero
	r(l0+3) = czero
	r(l0+4) = czero
	
	l = ny
	l0 = (l-1)*ndof
	r(l0+1) = czero
	r(l0+2) = czero
	r(l0+3) = czero
	r(l0+4) = czero

!.... Solve the system

#ifdef CRAY
	call CGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
	if (info.ne.0) call error('solver$','Singular matrix$')
	call CGESL(A0,ny*ndof,ny*ndof,ipvt,r,0)
#else
	call ZGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
	if (info.ne.0) call error('solver$','Singular matrix$')
	call ZGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
#endif

!.... Update the solution

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    u(i,l,k,n,ldof) = r(l0+ldof)
	  end do
	end do

	ul = pt5 * ( u(i-1,:,k,n,1) + u(i,:,k,n,1) )
	vl = pt5 * ( u(i-1,:,k,n,2) + u(i,:,k,n,2) )
	wl = pt5 * ( u(i-1,:,k,n,3) + u(i,:,k,n,3) )
	pl = pt5 * ( u(i-1,:,k,n,4) + u(i,:,k,n,4) )

	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)

!.... compute some streamwise derivatives

	ux = (u(i,:,k,n,1) - u(i-1,:,k,n,1))/(x(i)-x(i-1))
	vx = (u(i,:,k,n,2) - u(i-1,:,k,n,2))/(x(i)-x(i-1))
	wx = (u(i,:,k,n,3) - u(i-1,:,k,n,3))/(x(i)-x(i-1))

!.... compute the correction for alpha (following Herbert AIAA-93-3053)

	ke(i) = zero
	kex   = zero
	do j = 1, ny-1
	  ke(i) = ke(i) + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                                  conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
	  kex = kex + pt5 * ( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
                              conjg(wl(j))*wx(j) + &
                          conjg(ul(j+1))*ux(j+1) + conjg(vl(j+1))*vx(j+1) + &
                          conjg(wl(j+1))*wx(j+1) ) * (y(j+1)-y(j))
	end do
	alpha(i,k,n) = alpha(i,k,n) - iota * kex / ke(i)

	write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                      abs(amp(i,k,n)), ke(i), abs(kex/ke(i))/abs(alpha(i,k,n))

	if (abs(kex/ke(i))/abs(alpha(i,k,n)) .lt. tol) exit

	end do    ! loop on iter

!	write(*,*) 'i, x(i), alpha(i,k,n)', i, x(i), alpha(i,k,n)

!=============================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=============================================================================

	write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                                    amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
	call flush(9)
      end do   ! loop on i

!.... output the final profile

      i = ie
      ul = u(i,:,k,n,1)
      vl = u(i,:,k,n,2)
      wl = u(i,:,k,n,3)
      pl = u(i,:,k,n,4)
!     if (ipwbc .eq. 1) pl(1) = pl(2) - (pl(3)-pl(2))/(y(3)-y(2))*y(2)
      open(10,file='pro.dat')
      write(10,"('# ',10(1pe20.13,1x))") one, zero,  real(alpha(i,k,n)), &
                                         aimag(alpha(i,k,n))

!=============================================================================
!     1    2     3     4     5     6     7     8     9
!     y,  u_r,  u_i,  v_r,  v_i,  w_r,  w_i,  p_r,  p_i
!=============================================================================

      do j = 1, ny
	write(10,"(10(1pe20.13,1x))") y(j), &
                                      real(ul(j)), aimag(ul(j)), &
                                      real(vl(j)), aimag(vl(j)), &
                                      real(wl(j)), aimag(wl(j)), &
                                      real(pl(j)), aimag(pl(j))
      end do
      close(10)

!.... save a non interpolated field file

      open  (10, file='field.pse', form='unformatted')
      write (10) ie-is+1, ny
      write (10) (((u(i,j,k,n,idof), j= 1,ny), i= is,ie), idof= 1,4)
      write (10) (x(i), i= is,ie)
      write (10) (y(j), j= 1,ny)
      write (10) (amp(i,k,n), i = is, ie)
      write (10) (((zero, j= 1,ny), i= is,ie), idof= 1,4)
      write (10) (alpha(i,k,n), i= is,ie)
      close (10)

!.... save a field file (possibly interpolated)

      call post(1)

!.... compute growth rate based on various other quantitities

      call growth(ie-is+1, ny, x(is:ie), y, amp(is:ie,k,n), &
                  alpha(is:ie,k,n), u(is:ie,:,k,n,:), opi, nint, yint)

      return
      end subroutine solver
