!==============================================================================
      subroutine solver
!
!     Solve the linear PSE using Chebyshev collocation
!
!     Note: this is the first-order accurate scheme
!
!     S. Scott Collis
!
!     8-21-98:  1) Replaced trapezoidal integration with Chebyschev integration
!               2) Added normalization option (norm=0) Kinetic energy
!                                             (norm=1) Umax
!               3) Added continuity at wall for pressure (ipwbc=2)
!               4) Commented interpolation for ipwbc=1 on output
!     11-2-98:  1) Added field output
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, iter, il
      integer :: l, ldof, l0, m, mdof, m0, info, job=0, icount=0
      integer, allocatable :: ipvt(:)

      real :: rei

      complex :: dalpha(nz,nt)
      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)
      complex :: ux(ny),  vx(ny),  wx(ny),  px(ny)
      complex :: kex
      real    :: ke(nx)

!.... local variables

      complex :: alphal
      real    :: ubl(ny),  vbl(ny),  wbl(ny)
      real    :: ubxl(ny), vbxl(ny), wbxl(ny)
      real    :: ubyl(ny), vbyl(ny), wbyl(ny)
      complex :: uy(ny,ndof), uyy(ny,ndof), umax(nx), dumax

      real, allocatable :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
	                   E(:,:,:), Ex(:,:,:)

      complex, allocatable :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:)

      real :: c1(ny), c2(ny)

!.... wall BCs

      real :: bloc = 100.0, bcon = 0.0, bamp = -1.0e-1
      real :: coe, bl

      namelist /bc/ bloc, bcon, bamp

      real :: rdum, alphar, alphai
!==============================================================================

!.... get BC information (this doesn't work -- a weakness of PSE?)

!     open(10,file='bc.inp')
!     read(10,bc)
!     close(10)

!.... Read in alpha from a regular PSE run

      if (.false.) then
	open(20,file='pse.in')
	do i = is, ie
	  read(20,*) rdum, alphar, alphai
	  alpha(i,1,1) = cmplx(alphar, alphai)
	end do
	close(20)
      end if

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

      ke = zero
      kex = zero
      umax = zero
      dumax = zero

      if (norm.eq.0) then
	ke(i) = zero
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
                  abs(wl(j))**2)
	end do
      else if (norm.eq.1) then
	call findumax( ny, y, ul, umax(i), nint, yint )
	ke(i) = zero
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
                  abs(wl(j))**2)
	end do
      else if (norm.eq.2) then
	ke(i) = zero
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
	end do	
      else
	call error('Solver$','Illegal value for norm$')
      end if

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

      open(9,file='pse.dat')
      write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                                  amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
      call flush(9)

!.... main marching loop

      do i = is+1, ie

	cpu2 = second()
        write(*,"(/,'  i       x(i)          cpu        total cpu')")
	write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
	cpu = cpu2

	alpha(i,k,n) = alpha(i-1,k,n)

	c1(:) = one/(one + cur(i)*y(:))
	c2(:) = cur(i) * c1(:)

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
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
                                  abs(wl(j))**2)
	end do

	if (norm.eq.0) then
          write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
                   & '          ke','           dalpha')")
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
	        abs(amp(i,k,n)), ke(i), zero
	else if (norm.eq.1) then
	  call findumax( ny, y, ul, umax(i), nint, yint )
          write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
                   & '         |umax|','         dalpha')")
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(umax(i)), zero
	else if (norm.eq.2) then
	  ke(i) = zero
	  do j = 1, ny
	    ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
	  end do
          write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
                   & '          ke','           dalpha')")
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
	        abs(amp(i,k,n)), ke(i), zero	  
	end if

	do iter = 1, niter ! loop on alpha
	
!.... compute the streamwise derivative of alpha

	dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )

!.... Bertollotti sets this to zero (no advantage)

!	dalpha(k,n) = zero

!.... Evaluate the matrices (for Backward Euler)

	alphal = alpha(i,k,n)
	ubl    = ub(:,i)
	vbl    = vb(:,i)
	wbl    = wb(:,i)
	ubxl   = ubx(:,i)
	vbxl   = vbx(:,i)
	wbxl   = wbx(:,i)
	ubyl   = uby(:,i)
	vbyl   = vby(:,i)
	wbyl   = wby(:,i)

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
		                      B(l,ldof,mdof) * opy(l,m) - &
		                      E(l,ldof,mdof) * opyy(l,m)
	      end do
	      m0 = (l-1)*ndof
	      A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Gb(l,ldof,mdof) + &
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
	else if (ipwbc .eq. 2) then
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+2) = opy(l,m)
	  end do
	  A0(l0+4,m0+2) = A0(l0+4,m0+2) + c2(l)
	  A0(l0+4,m0+3) = -iota * kz(k) * beta
	  A0(l0+4,m0+1) = c1(l) / (x(i)-x(i-1))
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

	r = czero

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do mdof = 1, ndof
	      r(l0+ldof) = r(l0+ldof) + Ab(l,ldof,mdof) * &
                           u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
	    end do
	  end do
	end do

!.... Enforce the boundary conditions on the RHS

	l = 1
	l0 = (l-1)*ndof
	r(l0+1) = czero

!.... try modelling a suction/blowing slot

	if (bcon.eq.zero) then
	  r(l0+2) = czero
	else
	  coe = alog(.01) / bcon **2
	  bl = bamp * exp(coe * (x(i) - bloc) **2)
	  write(80,*) x(i), bl
	  r(l0+2) = cmplx( bl, zero )
	end if

	r(l0+3) = czero
	if (ipwbc.eq.2) then
	  r(l0+4) = c1(l) * u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
	else
	  r(l0+4) = czero
	end if
	
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

	ul = u(i,:,k,n,1)
	vl = u(i,:,k,n,2)
	wl = u(i,:,k,n,3)
	pl = u(i,:,k,n,4)

!.... compute some streamwise derivatives

	ux = (ul - u(i-1,:,k,n,1))/(x(i)-x(i-1))
	vx = (vl - u(i-1,:,k,n,2))/(x(i)-x(i-1))
	wx = (wl - u(i-1,:,k,n,3))/(x(i)-x(i-1))

	ke(i) = zero
	kex   = zero
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*( abs(ul(j))**2 + abs(vl(j))**2 + &
	          abs(wl(j))**2 )
	  kex = kex + opi(j)*( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
                conjg(wl(j))*wx(j) )
	end do

	if (norm.eq.0) then

!.... compute the correction for alpha (following Herbert AIAA-93-3053)

	  alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                              sor * (alpha(i,k,n) - iota * kex / ke(i))
	else if (norm.eq.1) then

!.... compute the correction for alpha (following Bertolotti)

	  call findumax( ny, y, ul, umax(i), nint, yint )
	  dumax = (umax(i)-umax(i-1))/(x(i)-x(i-1))
	  
!	  write(88,"(i4,1x,4(1pe13.6,1x))") i, umax(i), dumax
	  alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                              sor * (alpha(i,k,n) - iota * dumax / umax(i))
	else if (norm.eq.2) then
	  ke(i) = zero
	  kex   = zero
	  do j = 1, ny
	    ke(i) = ke(i) + opi(j)*( abs(ul(j))**2  )
	    kex = kex + opi(j)*( conjg(ul(j))*ux(j) )
	  end do
	  alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                              sor*(alpha(i,k,n) - iota * kex / ke(i))
	end if

!.... update the mode amplitude

	amp(i,k,n) = czero
	do il = is, i-1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
		       (x(il+1)-x(il))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )

!.... convergence check

	if (norm.eq.0) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), abs(kex/ke(i))/abs(alpha(i,k,n))
!	  if (abs(kex/ke(i))/abs(alpha(i,k,n)) .lt. tol) exit
	  if (abs(kex/ke(i)) .lt. tol) exit
	else if (norm.eq.1) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(umax(i)), &
                abs(dumax/umax(i))/abs(alpha(i,k,n))
!	  if (abs(dumax/umax(i))/abs(alpha(i,k,n)) .lt. tol) exit
	  if (abs(dumax/umax(i)) .lt. tol) exit
	else if (norm.eq.2) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), abs(kex/ke(i))/abs(alpha(i,k,n))
!	  if (abs(kex/ke(i))/abs(alpha(i,k,n)) .lt. tol) exit
	  if (abs(kex/ke(i)) .lt. tol) exit
	end if

	end do    ! loop on iter

!	write(*,*) 'i, x(i), alpha(i,k,n)', i, x(i), alpha(i,k,n)

!=============================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=============================================================================
	write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), &
	         abs(kex/ke(i)), amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
	call flush(9)

      end do   ! loop on i

!.... output the final profile (restart file)

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

!.... save a field file

      call post(1)

!.... compute growth rate based on various other quantitities

      call growth(ie-is+1, ny, x(is:ie), y, amp(is:ie,k,n), &
                  alpha(is:ie,k,n), u(is:ie,:,k,n,:), opi, nint, yint)


      return
      end subroutine solver
