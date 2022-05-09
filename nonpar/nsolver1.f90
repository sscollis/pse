!==============================================================================
      subroutine nsolver
!
!     Solves the nonlinear PSE
!
!     S. Scott Collis
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, iter, il
      integer :: l, ldof, l0, m, mdof, m0, info, job=0, icount=0
      integer, allocatable :: ipvt(:)
      logical :: new(nz,mt)

      real :: rei

      complex :: dalpha(nz,mt)
      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)
      complex :: ux(ny),  vx(ny),  wx(ny),  px(ny)
      complex :: kex
      real    :: tke, ke, unorm, Fnorm, malp

      real, allocatable :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
	                   E(:,:,:), Ex(:,:,:)

      complex, allocatable :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:), F(:,:,:,:)
      complex, allocatable :: ur(:,:,:), Fl(:,:)

      real :: c1(ny), c2(ny), dum
      real, parameter :: eps = 1.0e-8
!==============================================================================

!.... Initialize the FFT's

      call fft( 0, nt, nz, dum )

      rei = one / re

      allocate( u(nx,ny,nz,mt,ndof), ur(2*mt,nz,12) )

      u = czero

      allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
	        C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
	        Ex(ny,ndof,ndof) )

      allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
      allocate( r(ny*ndof), ipvt(ny*ndof), F(mt,nz,ny,ndof), Fl(ny,ndof) )

!.... initialize the disturbance field

      i = is
      u(i,:,:,:,1) = ubci
      u(i,:,:,:,2) = vbci
      u(i,:,:,:,3) = wbci
      u(i,:,:,:,4) = pbci

!.... main marching loop

      do i = is+1, ie
	
	new = .false.

	cpu2 = second()
	write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
	cpu = cpu2

	alpha(i,:,:) = alpha(i-1,:,:)
	u(i,:,:,:,:) = u(i-1,:,:,:,:)

	c1(:) = one/(one + cur(i)*y(:))
	c2(:) = cur(i) * c1(:)

	do iter = 1, niter   ! iterate on alpha and nonlinear source term

	malp = zero

!..... form the nonlinear term

	call nonlin( i, ur, ur, F, F, c1, c2, unorm, tke )

!.... loop over the nonzero modes

	do n = 1, mt
	do k = 1, nz

        ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)

!.... compute the mode amplitude

	amp(i,k,n) = czero
	do il = 1, i-1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
		       (x(il+1)-x(il))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )

!.... compute the mode (k,n) kinetic energy

	ke = zero
	do j = 1, ny-1
	ke = ke + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                          conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
        end do

!.... localize the nonlinear term

	Fl = F(n,k,:,:)

!.... compute a norm of the (n,k) nonlinear term

	Fnorm = zero
	do j = 1, ny-1
!	  if (n.eq.1 .and. k.eq.1) then
!	    write (8,"(10(1pe13.6,1x))") y(j), real(Fl(j,1)), aimag(Fl(j,1)), &
!	                                 real(Fl(j,2)), aimag(Fl(j,2)), &
!	                                 real(Fl(j,3)), aimag(Fl(j,3)), &
!	                                 real(Fl(j,4)), aimag(Fl(j,4))
!	  end if
	  Fnorm = Fnorm + pt5 * ( abs(Fl(j,1))**2 + abs(Fl(j,2))**2 + &
                                  abs(Fl(j,3))**2 + abs(Fl(j+1,1))**2 + &
                                  abs(Fl(j+1,2))**2 + abs(Fl(j+1,3))**2 ) * &
                          (y(j+1)-y(j))
	end do
	
!	write(*,"('=> ',2(i2,1x),6(1pe13.6,1x))") kz(k), kt(n), Fnorm/unorm, &
!                                          ke/tke, abs(amp(i,k,n)), ke, Fnorm

	if ( (Fnorm/unorm .le. eps) .and. (ke/tke .le. eps) ) goto 800

	if (ke .eq. zero) then
!	  alpha(i,k,n) = cmplx(real(kt(n)*alpha(i,1,2)),zero)
	  alpha(i,k,n) = kt(n)*alpha(i,1,2)
	  alpha(i-1,k,n) = alpha(i,k,n)
	  new(k,n) = .true.
	end if

	dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )

	G=zero; A=zero; B=zero; C=zero; D=zero; E=zero; Ex=zero
	Gb=czero; Ab=czero

	G(:,1,1) = one
	G(:,2,2) = one
	G(:,3,3) = one

	A(:,1,1) = c1(:) * ub(:,i)
	A(:,1,2) = -rei * two * c1(:) * c2(:)
	A(:,1,4) = c1(:)
	A(:,2,1) = rei * two * c1(:) * c2(:)
	A(:,2,2) = c1(:) * ub(:,i)
	A(:,3,3) = c1(:) * ub(:,i)
	A(:,4,1) = c1(:)

	B(:,1,1) = vb(:,i) - rei * c2(:)
	B(:,2,2) = vb(:,i) - rei * c2(:)
	B(:,2,4) = one
	B(:,3,3) = vb(:,i) - rei * c2(:)
	B(:,4,2) = one

	C(:,1,1) = wb(:,i)
	C(:,2,2) = wb(:,i)
	C(:,3,3) = wb(:,i)
	C(:,3,4) = one
	C(:,4,3) = one

	D(:,1,1) = c1(:) * ubx(:,i) + c2(:) * vb(:,i) - rei * c2(:)**2
	D(:,1,2) = uby(:,i)
	D(:,2,1) = c1(:) * vbx(:,i) - two * c2(:) * ub(:,i)
	D(:,2,2) = vby(:,i) + rei * c2(:)**2
	D(:,3,1) = c1(:) * wbx(:,i)
	D(:,3,2) = wby(:,i)
	D(:,4,2) = c2(:)

	Ex(:,1,1) = c1(:)**2 * rei
	Ex(:,2,2) = c1(:)**2 * rei
	Ex(:,3,3) = c1(:)**2 * rei

	E(:,1,1) = rei
	E(:,2,2) = rei
	E(:,3,3) = rei

	Gb = -iota * kt(n) * omega * G + iota * alpha(i,k,n) * A + &
	      iota * kz(k) * beta * C + D - &
	      Ex * ( iota * dalpha(k,n) - alpha(i,k,n)**2 ) + &
	      kz(k)**2 * beta**2 * E

!.... fix the streamwise pressure gradient

	if (ipfix.eq.1 .or. (ipfix.eq.2 .and. i.eq.2) ) A(:,:,4) = zero

	Ab = A - two * iota * alpha(i,k,n) * Ex

!.... Build the LHS matrix

	A0 = czero

	if (new(k,n)) then

	do ldof = 1, ndof
	  do mdof = 1, ndof
	    do l = 1, ny
	      l0 = (l-1)*ndof
	      do m = 1, ny
		m0 = (m-1)*ndof
		A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
		                      B(l,ldof,mdof) * opy(l,m) - &
		                      E(l,ldof,mdof) * opyy(l,m)
	      end do
	      m0 = (l-1)*ndof
	      A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Gb(l,ldof,mdof)
	    end do
	  end do
	end do

	else

	do ldof = 1, ndof
	  do mdof = 1, ndof
	    do l = 1, ny
	      l0 = (l-1)*ndof
	      do m = 1, ny
		m0 = (m-1)*ndof
		A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
		                      B(l,ldof,mdof) * opy(l,m) - &
		                      E(l,ldof,mdof) * opyy(l,m)
	      end do
	      m0 = (l-1)*ndof
	      A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Gb(l,ldof,mdof) + &
	                            Ab(l,ldof,mdof) / (x(i)-x(i-1))
	    end do
	  end do
	end do

	end if

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

	if (k.eq.1 .and. n.eq.1) then  ! always use for 00 mode
	  A0(l0+1,:) = czero
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+1,m0+1) = opy(l,m)
	  end do
	else
	  A0(l0+1,:) = czero
	  A0(l0+1,l0+1) = cone
	endif

	if (k.eq.1 .and. n.eq.1) then  ! always use for 00 mode
	  A0(l0+2,:) = czero
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+2,m0+2) = opy(l,m)
	  end do
	else
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

	r = czero

	if (new(k,n)) then
	  do ldof = 1, ndof
	    do l = 1, ny
	      l0 = (l-1)*ndof
	      r(l0+ldof) = r(l0+ldof) + F(n,k,l,ldof) / amp(i,k,n)
	    end do
	  end do
	else
	  do ldof = 1, ndof
	    do l = 1, ny
	      l0 = (l-1)*ndof
	      do mdof = 1, ndof
		r(l0+ldof) = r(l0+ldof) + Ab(l,ldof,mdof) * &
                             u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
	      end do
	      r(l0+ldof) = r(l0+ldof) + F(n,k,l,ldof) / amp(i,k,n)
	    end do
	  end do
	end if

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

!.... update local values

	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)

!.... compute streamwise derivatives

	ux = (ul - u(i-1,:,k,n,1)) / (x(i)-x(i-1))
	vx = (vl - u(i-1,:,k,n,2)) / (x(i)-x(i-1))
	wx = (wl - u(i-1,:,k,n,3)) / (x(i)-x(i-1))
	px = (pl - u(i-1,:,k,n,4)) / (x(i)-x(i-1))

!.... compute the correction for alpha (following Herbert AIAA-93-3053)

	ke = zero
	kex = zero
	do j = 1, ny-1
	  ke = ke + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                            conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
	  kex = kex + pt5 * ( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
                              conjg(wl(j))*wx(j) + &
                          conjg(ul(j+1))*ux(j+1) + conjg(vl(j+1))*vx(j+1) + &
                          conjg(wl(j+1))*wx(j+1) ) * (y(j+1)-y(j))
	end do
	
	if (kex .ne. zero) then
	  if (new(k,n)) kex = czero
	  alpha(i,k,n) = alpha(i,k,n) - iota * kex / ke

	  write(*,"(3(i2,1x),10(1pe13.6,1x))") kz(k), kt(n), iter, &
                               alpha(i,k,n), abs(amp(i,k,n)), &
                               ke, abs(kex/ke)/&
                               max(abs(alpha(i,k,n)),1e-10)
	else
	  write(*,"(3(i2,1x),10(1pe13.6,1x))") kz(k), kt(n), iter, &
	                       alpha(i,k,n), abs(amp(i,k,n)), ke, zero
	end if

	malp = max( abs(kex/ke)/max(abs(alpha(i,k,n)),1e-10), malp )

800	continue

	end do    ! loop on k

	end do    ! loop on n

	if (malp .le. tol) goto 900 ! exit

	end do    ! loop on iter

!.... output the statistics for low wavenumbers

900	continue

	k = 1; n = 1
	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)
	ke = zero
	do j = 1, ny-1
	  ke = ke + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                            conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
	end do
	write(20,"(10(1pe20.13,1x))") &
          x(i), alpha(i,k,n), dalpha(k,n), amp(i,k,n), abs(amp(i,k,n)), &
          abs(sqrt(ke)*amp(i,k,n))

	k = 1; n = 2
	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)
	ke = zero
	do j = 1, ny-1
	  ke = ke + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                            conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
	end do
	write(21,"(10(1pe20.13,1x))") &
          x(i), alpha(i,k,n), dalpha(k,n), amp(i,k,n), abs(amp(i,k,n)), &
          abs(sqrt(ke)*amp(i,k,n))

	k = 1; n = 3
	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)
	ke = zero
	do j = 1, ny-1
	  ke = ke + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                            conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
	end do
	write(22,"(10(1pe20.13,1x))") &
          x(i), alpha(i,k,n), dalpha(k,n), amp(i,k,n), abs(amp(i,k,n)), &
          abs(sqrt(ke)*amp(i,k,n))

	k = 1; n = 4
	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)
	ke = zero
	do j = 1, ny-1
	  ke = ke + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                            conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
	end do
	write(23,"(10(1pe20.13,1x))") &
          x(i), alpha(i,k,n), dalpha(k,n), amp(i,k,n), abs(amp(i,k,n)), &
          abs(sqrt(ke)*amp(i,k,n))

!.... output intermediate profiles

	if (.false.) then
	  k = 1
	  do n = 1, 3
	    k = 1
	    ul = u(i,:,k,n,1)
	    vl = u(i,:,k,n,2)
	    wl = u(i,:,k,n,3)
	    pl = u(i,:,k,n,4)
	    do j = 1, ny
	      write (10*(k-1)+10+n,"(10(1pe13.6,1x))") &
		y(j), ul(j), vl(j), wl(j), pl(j)
	    end do
	  end do
	  stop
	end if

      end do      ! loop on i

!.... output the final profile(s)

      i = ie
      do n = 1, mt
	ul = u(i,:,1,n,1)
	vl = u(i,:,1,n,2)
	wl = u(i,:,1,n,3)
	pl = u(i,:,1,n,4)
	if (ipwbc .eq. 1) pl(1) = pl(2) - (pl(3)-pl(2))/(y(3)-y(2))*y(2)
	do j = 1, ny
	  write(9+n,"(10(1pe20.13,1x))") y(j), &
                                         real(ul(j)), aimag(ul(j)), &
                                         real(vl(j)), aimag(vl(j)), &
                                         real(wl(j)), aimag(wl(j)), &
                                         real(pl(j)), aimag(pl(j))
	end do
      end do

      return
      end subroutine nsolver

!==============================================================================
!     E X T R A   D I A G N O S T I C   S T U F F 
!==============================================================================

!	if (n.eq.1 .and. k.eq.1 .and. i.eq.2) then
!!	  alpha(i,k,n) = pt5 * (alpha(i,1,2) - conjg(alpha(i,1,2)))
!	  alpha(i,k,n) = cmplx(real(kt(n)*alpha(i,1,2)),zero)
!	  alpha(i-1,k,n) = alpha(i,k,n)
!	  new(k,n) = .true.
!	end if

!	if (n.eq.3 .and. k.eq.1 .and. alpha(i,k,n).eq.zero) then
!!	  alpha(i,k,n) = kt(n) * alpha(i,1,2)
!	  alpha(i,k,n) = cmplx(real(two*alpha(i,1,2)),zero)
!	  alpha(i-1,k,n) = alpha(i,k,n)
!	  new(k,n) = .true.
!	end if

!	if (iter.eq.1 .and. new(k,n)) &
!	if (new(k,n)) &
!	  write(*,*) "WARNING:  new mode ", kz(k), kt(n)

