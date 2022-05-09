!==============================================================================
      subroutine dasolver
!
!     Solve the discrete adjoint PSE using Chebyshev collocation
!
!     Note: alpha here is equal to alpha of the regular PSE
!
!     S. Scott Collis
!
!     8-21-98:  1) Replaced trapezoidal integration with Chebyschev integration
!               2) Added normalization option (norm=0) Kinetic energy
!                                             (norm=1) Umax
!               3) Added continuity at wall for pressure (ipwbc=2)
!               4) Commented interpolation for ipwbc=1 on output
!     11-2-98:  1) Added field output
!     8-13-99:  1) Added Newton iteration which works great!
!     8-16-99:  1) Fixed the Ab term in discrete adjoint along with IC
!     9-21-99:  1) Added norm=4, alpha on full J iteration
!    10-20-99:  1) Fixed the annoying "phase scramble"
!==============================================================================
      use global
      use int_str
      implicit none

      integer :: i, j, k, n, iter, il
      integer :: l, ldof, l0, m, mdof, m0, info, job=0, icount=0
      integer, allocatable :: ipvt(:)

      real :: rei

      complex :: dalpha(nz,nt)
      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)
      complex :: ux(ny),  vx(ny),  wx(ny),  px(ny)
      complex :: kex
      real    :: ke(nx), p(ny), yumax

!.... local variables

      complex :: alphal
      real    :: ubl(ny),  vbl(ny),  wbl(ny)
      real    :: ubxl(ny), vbxl(ny), wbxl(ny)
      real    :: ubyl(ny), vbyl(ny), wbyl(ny)
      complex :: uy(ny,ndof), uyy(ny,ndof), umax(nx), dumax

      real, allocatable :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
	                   E(:,:,:), Ex(:,:,:), Exp1(:,:,:), Ap1(:,:,:)

      complex, allocatable :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:)
      complex :: ua(nx,ny,ndof), uax(ny), vax(ny), wax(ny), ra(ny*ndof)
      complex :: f, dfda, A0a(ny*ndof,ny*ndof), A0b(ny*ndof,ny*ndof)

      integer :: nxi, nyi, idof
      complex :: uh(nx,ny,ndof), amph(nx), alphah(nx), uha(nx,ny,ndof)
      complex :: ut(nx,ny,ndof), ampt(nx), alphat(nx), uta(nx,ny,ndof)

      complex :: jprod(nx), jprodx, jprod1, jprod2, jprod3
      complex :: dfdah(nx), dfdat(nx), f1, f2, det, dalphah, dalphat

      real :: c1(ny), c2(ny)

      real :: rdum, alphar, alphai
      integer :: idum

      character*80 :: fname

      logical :: scramble = .false.

      real :: fact = one
      logical :: stable = .false.
!==============================================================================

!.... Don't use the zero weights on the end points

      do i = 1, ny
	p(i) = sqrt(1.0 - (cos((i-1)*acos(-1.0)/(ny-1)))**2)
      end do
      do i = 2, ny-1
	dydeta(i) = dydeta(i) * p(i)
      end do
      dydeta(1) = dydeta(1) * pt5
      dydeta(ny) = dydeta(ny) * pt5

      dydeta = opi   ! this is real integration operator

      rei = one / re

!.... Read in alpha from a regular PSE run (assumes same mesh)

      if (norm.ge.3) then
	open (10, file='field.pse', form='unformatted')
	read (10) nxi, nyi
	if (nxi.ne.(ie-is+1) .or. nyi.ne.ny) &
	  call error('dasolver$','Illegal field.pse file$')
	read (10) (((uh(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
	read (10) (rdum, i= is,ie)
	read (10) (rdum, j= 1,ny)
	read (10) (amph(i), i = is, ie)
	read (10) (((uha(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
	read (10) (alphah(i), i= is, ie)
	close (10)	
      end if

      if (ipfix.eq.2) then
	open (10, file='field.adj', form='unformatted')
	read (10) nxi, nyi
	if (nxi.ne.(ie-is+1) .or. nyi.ne.ny) &
	  call error('asolver$','Illegal field.pse file$')
	read (10) (((ut(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
	read (10) (rdum, i= is,ie)
	read (10) (rdum, j= 1,ny)
	read (10) (ampt(i), i = is, ie)
	read (10) (((uta(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
	read (10) (alphat(i), i= is, ie)
	close (10)
      end if

      allocate( u(nx,ny,nz,nt,ndof) )

      u = czero

      allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
	        C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
	        Ex(ny,ndof,ndof), Ap1(ny,ndof,ndof), Exp1(ny,ndof,ndof) )

      allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
      allocate( r(ny*ndof), ipvt(ny*ndof) )

!.... only consider only one mode for linear analysis

      i = ie; n = 1; k = 1

!.... initialize the field
      
      do j = 1, ny
	u(i,j,k,n,1) = ubci(j,k,n)
	u(i,j,k,n,2) = vbci(j,k,n)
	u(i,j,k,n,3) = wbci(j,k,n)
	u(i,j,k,n,4) = pbci(j,k,n)
      end do

!.... normalize the adjoint initial condition

      if (norm.eq.3) then
	alpha(:,k,n) = alphah
	Jprod(i) = czero
	do j = 1, ny
	  Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i-1,j,1)*u(i,j,k,n,1)+uh(i-1,j,2)*u(i,j,k,n,2)+ &
             uh(i-1,j,3)*u(i,j,k,n,3))*ub(j,i) + &
	    uh(i-1,j,1)*u(i,j,k,n,4) + &
	    uh(i-1,j,4)*u(i,j,k,n,1) - &
	    two*rei*iota*alphah(i) * &
	    (uh(i-1,j,1)*u(i,j,k,n,1)+uh(i-1,j,2)*u(i,j,k,n,2)+ &
             uh(i-1,j,3)*u(i,j,k,n,3)) )
	end do
	ubci(:,k,n) = ubci(:,k,n) / Jprod(i)
	vbci(:,k,n) = vbci(:,k,n) / Jprod(i)
	wbci(:,k,n) = wbci(:,k,n) / Jprod(i)
	pbci(:,k,n) = pbci(:,k,n) / Jprod(i)
      else if (norm.eq.4) then
	alpha(:,k,n) = alphah
	Jprod(i) = czero
	do j = 1, ny
	  Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i,j,1)*u(i,j,k,n,1)+uh(i,j,2)*u(i,j,k,n,2)+ &
             uh(i,j,3)*u(i,j,k,n,3))*ub(j,i) + &
	    uh(i,j,1)*u(i,j,k,n,4) + &
	    uh(i,j,4)*u(i,j,k,n,1) - &
	    rei*iota*(alphah(i)+alpha(i,k,n)) * &
	    (uh(i,j,1)*u(i,j,k,n,1)+uh(i,j,2)*u(i,j,k,n,2)+ &
             uh(i,j,3)*u(i,j,k,n,3)) )
	end do
	ubci(:,k,n) = ubci(:,k,n) / Jprod(i)
	vbci(:,k,n) = vbci(:,k,n) / Jprod(i)
	wbci(:,k,n) = wbci(:,k,n) / Jprod(i)
	pbci(:,k,n) = pbci(:,k,n) / Jprod(i)        
      end if

!.... re-initialize the field
      
      do j = 1, ny
	u(i,j,k,n,1) = ubci(j,k,n)
	u(i,j,k,n,2) = vbci(j,k,n)
	u(i,j,k,n,3) = wbci(j,k,n)
	u(i,j,k,n,4) = pbci(j,k,n)
      end do

      ke = zero
      kex = zero
      umax = zero
      dumax = zero

      ul = u(i,:,k,n,1)
      vl = u(i,:,k,n,2)
      wl = u(i,:,k,n,3)
      pl = u(i,:,k,n,4)

      ke(i) = zero
      do j = 1, ny
	ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
	                        abs(wl(j))**2)
      end do

      if (norm.eq.0) then
      else if (norm.eq.1) then
	call findumax( ny, y, ul, umax(i), nint, yint, yumax )
      else if (norm.eq.2) then
	ke(i) = zero
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
	end do	
      else if (norm.eq.3) then
	Jprod(i) = czero
	do j = 1, ny
	  Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i-1,j,1)*u(i,j,k,n,1)+uh(i-1,j,2)*u(i,j,k,n,2)+ &
             uh(i-1,j,3)*u(i,j,k,n,3) )*ub(j,i) + &
	    uh(i-1,j,1)*u(i,j,k,n,4) + &
	    uh(i-1,j,4)*u(i,j,k,n,1) - &
	    two*rei*iota*alphah(i) * &
	    (uh(i-1,j,1)*u(i,j,k,n,1)+uh(i-1,j,2)*u(i,j,k,n,2)+ &
             uh(i-1,j,3)*u(i,j,k,n,3)) )
	end do
      else if (norm.eq.4) then
	Jprod(i) = czero
	do j = 1, ny
	  Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+uh(i,j,3)*wl(j))*ub(j,i) + &
	    uh(i,j,1)*pl(j) + uh(i,j,4)*ul(j) - &
	    rei*iota*(alphah(i)+alpha(i,k,n)) * &
	    (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+uh(i,j,3)*wl(j)) )
	end do
      else
	call error('daSolver$','Illegal value for norm$')
      end if

      open(9,file='dadj.dat')
      write(9,"(9(1pe20.13,1x))") x(i), -alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                                  amp(i,k,n), ke(i)*abs(amp(i,k,n))**2, &
				  abs(Jprod(i))
      call flush(9)

!.... main marching loop

      do i = ie-1, is+1, -1

	cpu2 = second()
        write(*,"(/,'  i       x(i)          cpu        total cpu')")
	write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
	cpu = cpu2

        if (norm.eq.3) then
	  alpha(i,k,n) = alphah(i)
	else
	  alpha(i,k,n) = alpha(i+1,k,n)
	end if

	c1(:) = one/(one + cur(i)*y(:))
	c2(:) = cur(i) * c1(:)

	ul = u(i+1,:,k,n,1)
	vl = u(i+1,:,k,n,2)
	wl = u(i+1,:,k,n,3)
	pl = u(i+1,:,k,n,4)

!.... compute the mode amplitude

	amp(i,k,n) = czero
	do il = ie, i+1, -1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il-1,k,n)) * &
		       (x(il)-x(il-1))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )  ! to integrate backwards [SSC]

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
	  call findumax( ny, y, ul, umax(i), nint, yint, yumax )
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
	else if (norm.eq.3) then
          Jprod(i) = czero
          do j = 1, ny
	    Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	      (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+ &
               uh(i,j,3)*u(i+1,j,k,n,3))*ub(j,i+1) + &
	      uh(i,j,1)*u(i+1,j,k,n,4) + &
	      uh(i,j,4)*u(i+1,j,k,n,1) - &
	      two*rei*iota*alphah(i+1) * &
	      (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+ &
               uh(i,j,3)*u(i+1,j,k,n,3)) )
          end do
          write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
                   & '          |J|','          dalpha')")
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
	        abs(amp(i,k,n)), abs(Jprod(i)), zero
	else if (norm.eq.4) then
          Jprod(i) = czero
          do j = 1, ny
	    Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	      (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+uh(i,j,3)*wl(j))*ub(j,i) + &
	      uh(i,j,1)*pl(j) + uh(i,j,4)*ul(j) - &
	      rei*iota*(alphah(i)+alpha(i,k,n))* &
	      (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+uh(i,j,3)*wl(j)) )
          end do
          write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
                   & '          |J|','          dalpha')")
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
	        abs(amp(i,k,n)), abs(Jprod(i)), zero
	end if

	do iter = 1, niter
	
!.... compute the mode amplitude

	amp(i,k,n) = czero
	do il = ie, i+1, -1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il-1,k,n)) * &
		       (x(il)-x(il-1))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )  ! to integrate backwards [SSC]

!.... compute the streamwise derivative of alpha

        dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )

!.... Evaluate the matrices (for Backward Euler)

	c1(:) = one/(one + cur(i+1)*y(:))
	c2(:) = cur(i+1) * c1(:)

	alphal = alpha(i+1,k,n)
	ubl    = ub(:,i+1)
	vbl    = vb(:,i+1)
	wbl    = wb(:,i+1)
	ubxl   = ubx(:,i+1)
	vbxl   = vbx(:,i+1)
	wbxl   = wbx(:,i+1)
	ubyl   = uby(:,i+1)
	vbyl   = vby(:,i+1)
	wbyl   = wby(:,i+1)

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

        Ab = A
	if ( (ipfix.eq.1) .or. (ipfix.eq.2) ) then
          Ab(:,:,4) = zero
        end if

	Ab = Ab - two * iota * alphal * Ex

!.... Build the RHS

	r = czero

        if (.not. stable) then

          do j = 1, ny
            Ab(j,:,:) = transpose( Ab(j,:,:) )
          end do
          
          do ldof = 1, ndof
            do l = 1, ny
              l0 = (l-1)*ndof
              do mdof = 1, ndof
                r(l0+ldof) = r(l0+ldof) + Ab(l,ldof,mdof) * u(i+1,l,k,n,mdof)
              end do
            end do
          end do
          
        else   !.... Stabilization?

          fact = ( one + tau / (x(i)-x(i-1)) )
          
          A0 = czero
          
          do ldof = 1, ndof
            do mdof = 1, ndof
              do l = 1, ny
                l0 = (l-1)*ndof
                do m = 1, ny
                  m0 = (m-1)*ndof
                  A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
                          (fact-one) * ( B(l,ldof,mdof) * opy(l,m) - &
                                         E(l,ldof,mdof) * opyy(l,m) )
                end do
                m0 = (l-1)*ndof
                A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
                                     (fact-one) * Gb(l,ldof,mdof) + &
	                              Ab(l,ldof,mdof) / (x(i)-x(i-1))
              end do
            end do
          end do

          A0 = A0 * (x(i)-x(i-1))
          
          A0 = transpose(A0) 

          do ldof = 1, ndof
            do mdof = 1, ndof
              do l = 1, ny
                l0 = (l-1)*ndof
                do m = 1, ny
                  m0 = (m-1)*ndof
                  r(l0+ldof) = r(l0+ldof) + &
                               A0(l0+ldof,m0+mdof) * u(i+1,m,k,n,mdof)
                end do
              end do
            end do
          end do
          
        end if

	if (ipfix.eq.2) then
	  ldof = 4
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    r(l0+ldof) = r(l0+ldof) + c1(l) * &
                	 (ut(i+1,l,1) - ut(i,l,1))
	  end do
        end if

!.... Enforce the boundary conditions on the RHS

	l = 1
	l0 = (l-1)*ndof
	r(l0+1) = czero
	r(l0+2) = czero
	r(l0+3) = czero
	if (ipwbc.eq.2) then
	  call error('dasolver$','ipwbc = 2 not implemented$')
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

!.... Build the LHS matrix

	c1(:) = one/(one + cur(i)*y(:))
	c2(:) = cur(i) * c1(:)

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

        Ab = A
	if ( (ipfix.eq.1) .or. (ipfix.eq.2) ) then
          Ab(:,:,4) = zero
        end if

	Ab = Ab - two * iota * alphal * Ex

	A0 = czero

	do ldof = 1, ndof
	  do mdof = 1, ndof
	    do l = 1, ny
	      l0 = (l-1)*ndof
	      do m = 1, ny
		m0 = (m-1)*ndof
		A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + fact * ( &
		                      B(l,ldof,mdof) * opy(l,m) - &
		                      E(l,ldof,mdof) * opyy(l,m) )
	      end do
	      m0 = (l-1)*ndof
	      A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
                                    fact * Gb(l,ldof,mdof) + &
	                            Ab(l,ldof,mdof) / (x(i)-x(i-1))
	    end do
	  end do
	end do

        A0 = A0 * (x(i)-x(i-1))

!.... Apply the boundary conditions to the LHS

	l = 1
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(:,l0+1) = czero
	A0(l0+1,l0+1) = cone

	A0(l0+2,:) = czero
	A0(:,l0+2) = czero
	A0(l0+2,l0+2) = cone

	A0(l0+3,:) = czero
	A0(:,l0+3) = czero
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

!.... statically eliminate the Neumann boundary conditions

	if (.true.) then

        A0a = A0
	
	if (ipwbc .eq. 0 .or. ipwbc .eq. 1) then

	  do l = 2, ny
	    l0 = (l-1)*ndof
	    do m = 2, ny
	      m0 = (m-1)*ndof
	      A0(l0+1,m0+4) = A0(l0+1,m0+4) - A0(l0+1,4) * A0(4,m0+4) / A0(4,4)
	      A0(l0+2,m0+4) = A0(l0+2,m0+4) - A0(l0+2,4) * A0(4,m0+4) / A0(4,4)
	      A0(l0+3,m0+4) = A0(l0+3,m0+4) - A0(l0+3,4) * A0(4,m0+4) / A0(4,4)
	      A0(l0+4,m0+4) = A0(l0+4,m0+4) - A0(l0+4,4) * A0(4,m0+4) / A0(4,4)
	    end do
	    r(l0+1) = r(l0+1) - A0(l0+1,4) * r(4) / A0(4,4)
	    r(l0+2) = r(l0+2) - A0(l0+2,4) * r(4) / A0(4,4)
	    r(l0+3) = r(l0+3) - A0(l0+3,4) * r(4) / A0(4,4)
	    r(l0+4) = r(l0+4) - A0(l0+4,4) * r(4) / A0(4,4)
	  end do
	  
	  l = 1
	  l0 = (l-1)*ndof
	  do m = 2, ny
	    m0 = (m-1)*ndof
	    A0(l0+1,m0+4) = A0(l0+1,m0+4) - A0(l0+1,4) * A0(4,m0+4) / A0(4,4)
	    A0(l0+2,m0+4) = A0(l0+2,m0+4) - A0(l0+2,4) * A0(4,m0+4) / A0(4,4)
	    A0(l0+3,m0+4) = A0(l0+3,m0+4) - A0(l0+3,4) * A0(4,m0+4) / A0(4,4)
	  end do	    
	  r(l0+1) = r(l0+1) - A0(l0+1,4) * r(4) / A0(4,4)
	  r(l0+2) = r(l0+2) - A0(l0+2,4) * r(4) / A0(4,4)
	  r(l0+3) = r(l0+3) - A0(l0+3,4) * r(4) / A0(4,4)
	  
	  A0(4,:) = A0(4,:) / A0(4,4)
	  r(4) = r(4) / A0(4,4)
	  A0(:,4) = zero
	  A0(4,4) = one
	  
	  do l = 2, ny
	    l0 = (l-1)*ndof
	    do m = 1, ny
	      m0 = (m-1)*ndof
	      A0(l0+4,m0+4) = A0(l0+4,m0+4) + opy(1,l) * A0(4,m0+4) / dydeta(l)
	    end do
	    r(l0+4) = r(l0+4) + opy(1,l) * r(4) / dydeta(l)
	  end do
	  
	  l = 1
	  l0 = (l-1)*ndof
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+4) = A0(l0+4,m0+4) * opy(1,l) / dydeta(l)
	  end do
	  r(l0+4) = r(l0+4) * opy(1,l) / dydeta(l)
	else 
	  call error('lst$','Illegal value of ipwbc$')
	end if

	end if

	l = 1
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(:,l0+1) = czero
	A0(l0+1,l0+1) = cone

	A0(l0+2,:) = czero
	A0(:,l0+2) = czero
	A0(l0+2,l0+2) = cone

	A0(l0+3,:) = czero
	A0(:,l0+3) = czero
	A0(l0+3,l0+3) = cone

!.... far-field boundary conditions

	l = ny
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(:,l0+1) = czero
	A0(l0+1,l0+1) = cone

	if (ivbc .eq. 0) then
	  A0(l0+2,:) = czero
	  A0(:,l0+2) = czero
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
	A0(:,l0+3) = czero
	A0(l0+3,l0+3) = cone

	A0(l0+4,:) = czero
	if (ipbc .eq. 0) then
	  A0(:,l0+4) = czero
	  A0(l0+4,l0+4) = cone	
	else if (ipbc .eq. 1) then
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+4) = opy(l,m)
	  end do
	else
	  call error('solver$','Illegal value of ipbc$')
	end if

!.... form the discrete adjoint

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do mdof = 1, ndof
	      do m = 1, ny
		m0 = (m-1)*ndof
		A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) * dydeta(l)
	      end do
	    end do
	    r(l0+ldof) = r(l0+ldof) * dydeta(l)
	  end do
	end do

	A0 = transpose( A0 )

!.... Solve the system

#ifdef CRAY
	call CGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
	if (info.ne.0) call error('solver$','Singular matrix$')
	call CGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
#else
!	call ZGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
!	if (info.ne.0) call error('solver$','Singular matrix$')
!	call ZGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
        call ZGETRF( ny*ndof, ny*ndof, A0, ny*ndof, ipvt, info)
	if (info.ne.0) call error('solver$','Singular matrix$')
        call ZGETRS( 'N', ny*ndof, 1, A0, ny*ndof, ipvt, r, ny*ndof, info)
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

	ux = (u(i+1,:,k,n,1) - ul)/(x(i+1)-x(i))
	vx = (u(i+1,:,k,n,2) - vl)/(x(i+1)-x(i))
	wx = (u(i+1,:,k,n,3) - wl)/(x(i+1)-x(i))

!==============================================================================

!.... Build the RHS for du/dalpha

	ra = czero

        Ab = (x(i)-x(i-1))*(iota * A + two * alphal * Ex) - three * iota * Ex 

	do j = 1, ny
	  Ab(j,:,:) = transpose( Ab(j,:,:) )
	end do

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do mdof = 1, ndof
	      ra(l0+ldof) = ra(l0+ldof) - Ab(l,ldof,mdof) * u(i,l,k,n,mdof)
	    end do
	  end do
	end do

        Ab = two * iota * Exp1

	do j = 1, ny
	  Ab(j,:,:) = transpose( Ab(j,:,:) )
	end do

	do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            do mdof = 1, ndof
              ra(l0+ldof) = ra(l0+ldof) + Ab(l,ldof,mdof) * u(i+1,l,k,n,mdof)
            end do
          end do
	end do
        
!.... Enforce the boundary conditions on the RHS

	l = 1
	l0 = (l-1)*ndof
	ra(l0+1) = czero
        ra(l0+2) = czero
	ra(l0+3) = czero
        ra(l0+4) = czero
	
	l = ny
	l0 = (l-1)*ndof
	ra(l0+1) = czero
	ra(l0+2) = czero
	ra(l0+3) = czero
	ra(l0+4) = czero

	if (.true.) then

        A0b = A0 ; A0 = A0a
	
	if (ipwbc .eq. 0 .or. ipwbc .eq. 1) then

	  do l = 2, ny
	    l0 = (l-1)*ndof
	    do m = 2, ny
	      m0 = (m-1)*ndof
	      A0(l0+1,m0+4) = A0(l0+1,m0+4) - A0(l0+1,4) * A0(4,m0+4) / A0(4,4)
	      A0(l0+2,m0+4) = A0(l0+2,m0+4) - A0(l0+2,4) * A0(4,m0+4) / A0(4,4)
	      A0(l0+3,m0+4) = A0(l0+3,m0+4) - A0(l0+3,4) * A0(4,m0+4) / A0(4,4)
	      A0(l0+4,m0+4) = A0(l0+4,m0+4) - A0(l0+4,4) * A0(4,m0+4) / A0(4,4)
	    end do
	    ra(l0+1) = ra(l0+1) - A0(l0+1,4) * ra(4) / A0(4,4)
	    ra(l0+2) = ra(l0+2) - A0(l0+2,4) * ra(4) / A0(4,4)
	    ra(l0+3) = ra(l0+3) - A0(l0+3,4) * ra(4) / A0(4,4)
	    ra(l0+4) = ra(l0+4) - A0(l0+4,4) * ra(4) / A0(4,4)
	  end do
	  
	  l = 1
	  l0 = (l-1)*ndof
	  do m = 2, ny
	    m0 = (m-1)*ndof
	    A0(l0+1,m0+4) = A0(l0+1,m0+4) - A0(l0+1,4) * A0(4,m0+4) / A0(4,4)
	    A0(l0+2,m0+4) = A0(l0+2,m0+4) - A0(l0+2,4) * A0(4,m0+4) / A0(4,4)
	    A0(l0+3,m0+4) = A0(l0+3,m0+4) - A0(l0+3,4) * A0(4,m0+4) / A0(4,4)
	  end do	    
	  ra(l0+1) = ra(l0+1) - A0(l0+1,4) * ra(4) / A0(4,4)
	  ra(l0+2) = ra(l0+2) - A0(l0+2,4) * ra(4) / A0(4,4)
	  ra(l0+3) = ra(l0+3) - A0(l0+3,4) * ra(4) / A0(4,4)
	  
	  A0(4,:) = A0(4,:) / A0(4,4)
	  ra(4) = ra(4) / A0(4,4)
	  A0(:,4) = zero
	  A0(4,4) = one
	  
	  do l = 2, ny
	    l0 = (l-1)*ndof
	    do m = 1, ny
	      m0 = (m-1)*ndof
	      A0(l0+4,m0+4) = A0(l0+4,m0+4) + opy(1,l) * A0(4,m0+4) / dydeta(l)
	    end do
	    ra(l0+4) = ra(l0+4) + opy(1,l) * ra(4) / dydeta(l)
	  end do
	  
	  l = 1
	  l0 = (l-1)*ndof
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+4) = A0(l0+4,m0+4) * opy(1,l) / dydeta(l)
	  end do
	  ra(l0+4) = ra(l0+4) * opy(1,l) / dydeta(l)
	else 
	  call error('lst$','Illegal value of ipwbc$')
	end if

	end if

	l = 1
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(:,l0+1) = czero
	A0(l0+1,l0+1) = cone

	A0(l0+2,:) = czero
	A0(:,l0+2) = czero
	A0(l0+2,l0+2) = cone

	A0(l0+3,:) = czero
	A0(:,l0+3) = czero
	A0(l0+3,l0+3) = cone

!.... far-field boundary conditions

	l = ny
	l0 = (l-1)*ndof

	A0(l0+1,:) = czero
	A0(:,l0+1) = czero
	A0(l0+1,l0+1) = cone

	if (ivbc .eq. 0) then
	  A0(l0+2,:) = czero
	  A0(:,l0+2) = czero
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
	A0(:,l0+3) = czero
	A0(l0+3,l0+3) = cone

	A0(l0+4,:) = czero
	if (ipbc .eq. 0) then
	  A0(:,l0+4) = czero
	  A0(l0+4,l0+4) = cone	
	else if (ipbc .eq. 1) then
	  do m = 1, ny
	    m0 = (m-1)*ndof
	    A0(l0+4,m0+4) = opy(l,m)
	  end do
	else
	  call error('solver$','Illegal value of ipbc$')
	end if

!.... form the discrete adjoint

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do mdof = 1, ndof
	      do m = 1, ny
		m0 = (m-1)*ndof
		A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) * dydeta(l)
	      end do
	    end do
	    ra(l0+ldof) = ra(l0+ldof) * dydeta(l)
	  end do
	end do

	A0 = transpose( A0 )

        A0 = A0b

#ifdef CRAY
	call CGESL(A0,ny*ndof,ny*ndof,ipvt,ra,job)
#else
!	call ZGESL(A0,ny*ndof,ny*ndof,ipvt,ra,job)
        call ZGETRS( 'N', ny*ndof, 1, A0, ny*ndof, ipvt, ra, ny*ndof, info)
#endif

!.... Update the solution

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    ua(i,l,ldof) = ra(l0+ldof)
	  end do
	end do

!.... Note that this is d/dalpha of the discrete dudx

	uax = -ua(i,:,1)/(x(i+1)-x(i))
	vax = -ua(i,:,2)/(x(i+1)-x(i))
	wax = -ua(i,:,3)/(x(i+1)-x(i))

!.... diagnostic

        if (.false. .and. mod(i,1) .eq. 0) then
           fname = 'dduda.'
           fname = trim(fname) // trim(i2c(i))
           open(10,file=fname)
           write(10,"('# ',10(1pe20.13,1x))") one, zero, &
                            real(alphal), aimag(alphal), omega
           do j = 1, ny
              write(10,"(10(1pe20.13,1x))") y(j), &
                   real(ua(i,j,1)), aimag(ua(i,j,1)), &
                   real(ua(i,j,2)), aimag(ua(i,j,2)), &
                   real(ua(i,j,3)), aimag(ua(i,j,3)), &
                   real(ua(i,j,4)), aimag(ua(i,j,4))
           end do
        end if

!==============================================================================

!.... compute the correction for alpha

	ke(i) = zero
	kex   = zero
        dfda  = czero
	do j = 1, ny
	  ke(i) = ke(i) + opi(j)*( abs(ul(j))**2 + abs(vl(j))**2 + &
	                           abs(wl(j))**2 )
	  kex = kex + opi(j)*( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
                               conjg(wl(j))*wx(j) )
          dfda = dfda + opi(j)*(conjg(ua(i,j,1))*ux(j) + conjg(ul(j))*uax(j) &
                              + conjg(ua(i,j,2))*vx(j) + conjg(vl(j))*vax(j) &
                              + conjg(ua(i,j,3))*wx(j) + conjg(wl(j))*wax(j) )
        end do

	if (norm.eq.0) then       ! following Herbert AIAA-93-3053
          f = kex
          if (.not. newton) dfda = -iota * ke(i)
	else if (norm.eq.1) then  ! following Bertolotti
	  call findumax( ny, y, ul, umax(i), nint, yint, yumax )
	  dumax = (umax(i+1)-umax(i))/(x(i+1)-x(i))
          f = dumax
          dfda = -iota * umax(i)
	else if (norm.eq.2) then
	  ke(i) = zero
	  kex   = zero
          dfda  = czero
	  do j = 1, ny
	    ke(i) = ke(i) + opi(j)*( abs(ul(j))**2  )
	    kex = kex + opi(j)*( conjg(ul(j))*ux(j) )
            dfda = dfda + opi(j)*(conjg(ua(i,j,1))*ux(j)+conjg(ul(j))*uax(j))
	  end do
          f = kex
          if (.not. newton) dfda = -iota * ke(i)
	else if (norm.eq.3) then
          Jprod(i) = czero
          do j = 1, ny
	    Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
	      (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2) + &
               uh(i,j,3)*u(i+1,j,k,n,3) )*ub(j,i+1) + &
	      uh(i,j,1)*u(i+1,j,k,n,4) + &
	      uh(i,j,4)*u(i+1,j,k,n,1) - &
	      two*rei*iota*alphah(i+1) * &
	      (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2) + &
               uh(i,j,3)*u(i+1,j,k,n,3) ) )
          end do
	  f = zero; dfda = one
	else if (norm.eq.4) then
          Jprod(i) = czero
          dfda = czero

!.... This is yet another hack to fix the phase scramble

          if (scramble) then
          if (mod(ie,2).eq.0) then
            if (mod(i+1,2).eq.0) then
              ul = -ul; vl = -vl; wl = -wl; pl = -pl
              ua(i,:,:) = -ua(i,:,:)
            end if
          else
            if (mod(i,2).eq.0) then
              ul = -ul; vl = -vl; wl = -wl; pl = -pl
              ua(i,:,:) = -ua(i,:,:)
            end if
          end if
          end if

          do j = 1, ny
            Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
                       (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+ &
                        uh(i,j,3)*wl(j)) * ub(j,i) + &
		       uh(i,j,1)*pl(j) + uh(i,j,4)*ul(j) - &
                       rei*iota*(alphah(i)+alpha(i,k,n)) * &
                       (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+uh(i,j,3)*wl(j)) )
            dfda = dfda + opi(j) * amp(i,k,n) * amph(i) * ( &
                       (uh(i,j,1)*ua(i,j,1)+uh(i,j,2)*ua(i,j,2)+ &
                        uh(i,j,3)*ua(i,j,3)) * ub(j,i) + &
		       uh(i,j,1)*ua(i,j,4) + uh(i,j,4)*ua(i,j,1) - &
                       rei*iota*(uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+ &
                                 uh(i,j,3)*wl(j)) - &
                       rei*iota*(alphah(i)+alpha(i,k,n)) * &
                       (uh(i,j,1)*ua(i,j,1)+uh(i,j,2)*ua(i,j,2)+ &
                        uh(i,j,3)*ua(i,j,3)) ) + &
                       opi(j) * (iota) * pt5 * (x(i+1)-x(i)) * &
                       amp(i,k,n) * amph(i) * ( & 
                       (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+ &
                        uh(i,j,3)*wl(j)) * ub(j,i) + &
		       uh(i,j,1)*pl(j) + uh(i,j,4)*ul(j) - &
                       rei*iota*(alphah(i)+alpha(i,k,n)) * &
                       (uh(i,j,1)*ul(j)+uh(i,j,2)*vl(j)+uh(i,j,3)*wl(j)) )
          end do
	  Jprodx = (Jprod(i+1)-Jprod(i))/(x(i+1)-x(i))
          dfda   = - dfda / (x(i+1)-x(i))
          f      = Jprodx
	end if

        alpha(i,k,n) = (one-sor)* alpha(i,k,n) + &
                            sor * (alpha(i,k,n) - f / dfda )

!.... update the mode amplitude

	amp(i,k,n) = czero
	do il = ie, i+1, -1
	  amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il-1,k,n)) * &
		       (x(il)-x(il-1))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )  ! to integrate backwards [SSC]

!.... convergence check

	if (norm.eq.0) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), abs(f/dfda)/abs(alpha(i,k,n))
	else if (norm.eq.1) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(umax(i)), abs(f/dfda)/abs(alpha(i,k,n))
	else if (norm.eq.2) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), abs(f/dfda)/abs(alpha(i,k,n))
	else if (norm.ge.3) then
	  write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
	        abs(amp(i,k,n)), abs(Jprod(i)), abs(f/dfda)/abs(alpha(i,k,n))
	end if

	if (abs(f/dfda)/abs(alpha(i,k,n)).lt.tol) exit
!	if (abs(f/dfda).lt.tol) exit

	end do    ! loop on iter

!=============================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=============================================================================
	write(9,"(9(1pe20.13,1x))") x(i), -alpha(i,k,n), ke(i), &
	         abs(kex/ke(i)), amp(i,k,n), ke(i)*abs(amp(i,k,n))**2, &
		 abs(Jprod(i))
	call flush(9)

      end do   ! loop on i
      close(9)

!.... Quick fix for phase scramble (need to figure out why this is needed?)

      if (scramble) then
      write(*,"(/,'Fixing phase scramble...?')")
      if (mod(ie,2).eq.0) then
        do i = is+1, ie-1
          if (mod(i+1,2).eq.0) u(i,:,k,n,:) = -u(i,:,k,n,:)
        end do
      else
        do i = is+1, ie-1
          if (mod(i,2).eq.0) u(i,:,k,n,:) = -u(i,:,k,n,:)
        end do
      end if
      end if

!.... output the final profile (restart file)

      i = is+1
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

      open  (10, file='field.adj', form='unformatted')
      write (10) ie-is+1, ny
      write (10) (((u(i,j,k,n,idof), j= 1,ny), i= is,ie), idof= 1,4)
      write (10) (x(i), i= is,ie)
      write (10) (y(j), j= 1,ny)
      write (10) (amp(i,k,n), i = is, ie)
      write (10) (((ua(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
      write (10) (alpha(i,k,n), i= is,ie)
      close (10)

!.... save a few versions of the J inner product

      open(22,file='jprod.dadj')
      do i = is+1, ie-1
        Jprod(i) = czero
	Jprod1 = czero
	Jprod2 = czero
	Jprod3 = czero
	do j = 1, ny
          Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * amph(i) * ( &
            (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+&
             uh(i,j,3)*u(i+1,j,k,n,3))*ub(j,i+1) + &
            uh(i,j,1)*u(i+1,j,k,n,4) + &
            uh(i,j,4)*u(i+1,j,k,n,1) + &
            (-two)*rei*iota*alphah(i+1) * &
            (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+&
             uh(i,j,3)*u(i+1,j,k,n,3)) )
	  Jprod1 = Jprod1 + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+&
             uh(i,j,3)*u(i+1,j,k,n,3))*ub(j,i+1) + &
!	    uh(i,j,1)*u(i+1,j,k,n,4) + &
	    uh(i,j,4)*u(i+1,j,k,n,1) + &
	    (-two)*rei*iota*alphah(i+1) * &
	    (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+ &
            uh(i,j,3)*u(i+1,j,k,n,3) ) )
	  Jprod2 = Jprod2 + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2)+&
             uh(i,j,3)*u(i+1,j,k,n,3) )*ub(j,i+1) + &
	    uh(i,j,1)*u(i+1,j,k,n,4) + &
!	    uh(i,j,4)*u(i+1,j,k,n,1) + &
	    (-two)*rei*iota*alphah(i+1) * &
	    (uh(i,j,1)*u(i+1,j,k,n,1)+uh(i,j,2)*u(i+1,j,k,n,2) + &
             uh(i,j,3)*u(i+1,j,k,n,3)) )
	  Jprod3 = Jprod3 + opi(j) * amp(i,k,n) * amph(i) * ( &
	    (uh(i,j,1)*u(i,j,k,n,1)+uh(i,j,2)*u(i,j,k,n,2)+&
             uh(i,j,3)*u(i,j,k,n,3))*ub(j,i) + &
	    uh(i,j,1)*u(i,j,k,n,4) + uh(i,j,4)*u(i,j,k,n,1) - &
	    rei*iota*(alphah(i)+alpha(i,k,n)) * &
	    (uh(i,j,1)*u(i,j,k,n,1)+uh(i,j,2)*u(i,j,k,n,2)+&
             uh(i,j,3)*u(i,j,k,n,3)) )
	end do
	write(22,"(11(1pe13.6,1x))") x(i), &
                                    real(alpha(i,k,n)), aimag(alpha(i,k,n)), &
				    real(Jprod(i)), aimag(Jprod(i)), &
                                    abs(Jprod(i)), abs(Jprod1), abs(Jprod2), &
				    real(Jprod3), aimag(Jprod3), abs(Jprod3)
      end do
      close(22)

!.... save a field file

      alpha = -alpha
      call post(-1)

!.... compute growth rate based on various other quantitities

      call dgrowth(ie-(is+1)+1, ny, x(is+1:ie), y, amp(is+1:ie,k,n), &
                  alpha(is+1:ie,k,n), u(is+1:ie,:,k,n,:), opi, nint, yint, &
                  opy, rei, uby(:,is+1:ie), vby(:,is+1:ie), wby(:,is+1:ie))

      alpha = -alpha

      return
100   call error('asolver$','Error reading pse.in$')
      end subroutine dasolver
