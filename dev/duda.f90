!==============================================================================
      subroutine duda
!
!     Compute the derivative of the PSE shape function wrt alpha
!
!     Notes:  1) Validated for parallel flow
!             2)
!
!     S. Scott Collis
!
!==============================================================================
      use global
      use int_str
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
      integer :: nxl, nyl
      complex :: alphah(nx), alphat(nx)
      complex, allocatable :: uh(:,:,:), ut(:,:,:)
      complex :: Jprod, L2prod, dwda

      real :: ur, ui, vr, vi, wr, wi, prr, pri, ampr, ampi
      complex :: ampl
      character*1 :: cdum
      character*80 :: fname
!==============================================================================

      rei = one / re

!.... Read in alpha from a regular PSE run

      open(20,file='pse.in')
      do i = is, ie
        read(20,*) rdum, alphar, alphai
        alphah(i) = cmplx(alphar, alphai)
      end do
      close(20)
      write(*,"('Read pse.in')")

!.... Read in alpha from an adjoint PSE run

      open(20,file='adj.in')
      do i = is, ie
        read(20,*) rdum, alphar, alphai
        alphat(i) = cmplx(alphar, alphai)
      end do
      close(20)
      write(*,"('Read adj.in')")

      open(10,file='field.pse',form='unformatted')
      read(10) nxl, nyl 
      allocate(uh(nx,ny,ndof))
      read(10) (((uh(i,j,idof),j=1,ny),i=is,ie),idof=1,ndof)
      close(10)
      write(*,"('Read field.pse')")

      open(10,file='field.adj',form='unformatted')
      read(10) nx, ny
      allocate(ut(nx,ny,ndof))
      read(10) (((ut(i,j,idof),j=1,ny),i=is,ie),idof=1,ndof)
      close(10)
      write(*,"('Read field.adj')")

!==============================================================================
!     S o l v e   t h e   d d a l p h a   e q u a t i o n
!==============================================================================

      allocate( u(nx,ny,nz,nt,ndof) ); u = czero
      allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
	        C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
	        Ex(ny,ndof,ndof) )
      allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
      allocate( r(ny*ndof), ipvt(ny*ndof) )

!.... only consider only one mode for linear analysis

      i = is; n = 1; k = 1

      alpha(:,k,n) = alphah(:)

!.... initialize the field

      open(10, file='ddalpha.dat') 
      read(10,*) cdum, ampr, ampi
      ampl = cmplx(ampr,ampi)
      do j = 1, ny
        read(10,*) rdum, ur, ui, vr, vi, wr, wi, prr, pri
        u(i,j,k,n,1) = ampl * cmplx(ur, ui)
        u(i,j,k,n,2) = ampl * cmplx(vr, vi)
        u(i,j,k,n,3) = ampl * cmplx(wr, wi)
        u(i,j,k,n,4) = ampl * cmplx(prr, pri)
      end do
      close(10)

!.... compute the J inner product
        
      Jprod = czero
      do j = 1, ny
        Jprod = Jprod + opi(j) * ( &
                (uh(i,j,1)*ut(i,j,1) + uh(i,j,2)*ut(i,j,2))*ub(j,i) + &
                 uh(i,j,1)*ut(i,j,4) + ut(i,j,1)*uh(i,j,4) - &
                 rei*iota*(alphah(i)-alphat(i)) * &
                 (uh(i,j,1)*ut(i,j,1)+uh(i,j,2)*ut(i,j,2)) )
      end do

!.... compute the L2 inner product of the velocity

      L2prod = czero
      do j = 1, ny
        L2prod = L2prod + opi(j)*( uh(i,j,1)*ut(i,j,1) + &
                                   uh(i,j,2)*ut(i,j,2) + &
                                   uh(i,j,3)*ut(i,j,3) )
      end do

      dwda = Jprod/L2prod
      write(*,*) 'i, x(i) = ', i, x(i)
      write(*,*) 'Alphah = ', alphah(i)
      write(*,*) 'Alphat = ', alphat(i)
      write(*,*) 'Jprod = ', Jprod
      write(*,*) 'L2prod = ', L2prod
      write(*,*) 'Domega/Dalpha = ', dwda
      write(*,*)

!.... main marching loop

      do i = is+1, ie

!.... compute the J inner product
        
        Jprod = czero
        do j = 1, ny
          Jprod = Jprod + opi(j) * ( &
                (uh(i,j,1)*ut(i,j,1) + uh(i,j,2)*ut(i,j,2))*ub(j,i) + &
                 uh(i,j,1)*ut(i,j,4) + ut(i,j,1)*uh(i,j,4) - &
                 rei*iota*(alphah(i)-alphat(i)) * &
                 (uh(i,j,1)*ut(i,j,1)+uh(i,j,2)*ut(i,j,2)) )
        end do

!.... compute the L2 inner product of the velocity

        L2prod = czero
        do j = 1, ny
          L2prod = L2prod + opi(j)*( uh(i,j,1)*ut(i,j,1) + &
                                     uh(i,j,2)*ut(i,j,2) + &
                                     uh(i,j,3)*ut(i,j,3) )
        end do

        dwda = Jprod/L2prod

!.... HACK:  set d omega / d alpha = zero 

        dwda = czero


        write(*,*) 'i, x(i) = ', i, x(i)
        write(*,*) 'Alphah = ', alphah(i)
        write(*,*) 'Alphat = ', alphat(i)
        write(*,*) 'Jprod = ', Jprod
        write(*,*) 'L2prod = ', L2prod
        write(*,*) 'Domega/Dalpha = ', dwda
        write(*,*)

	c1(:) = one/(one + cur(i)*y(:))
	c2(:) = cur(i) * c1(:)

	ul = u(i-1,:,k,n,1)
	vl = u(i-1,:,k,n,2)
	wl = u(i-1,:,k,n,3)
	pl = u(i-1,:,k,n,4)
	
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

	if (ipfix.eq.1 .or. (ipfix.eq.2 .and. i.eq.2) ) then
          Ab = A; Ab(:,:,4) = zero
        end if

	Ab = Ab - two * iota * alphal * Ex

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

!.... Build the RHS for du/dalpha

	r = czero

        Ab = -iota * A - two * alphal * Ex + iota * dwda * G 

	do ldof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do mdof = 1, ndof
	      r(l0+ldof) = r(l0+ldof) + &
                   two * iota * Ex(l,ldof,mdof)* &
                   (uh(i,l,mdof)-uh(i-1,l,mdof))/(x(i)-x(i-1)) + &
                   Ab(l,ldof,mdof) * uh(i,l,mdof)
	    end do
	  end do
	end do

	if (ipfix.eq.1 .or. (ipfix.eq.2 .and. i.eq.2) ) then
           Ab = A; Ab(:,:,4) = zero
        end if
        
	Ab = Ab - two * iota * alphal * Ex
        
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
        r(l0+2) = czero
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
	call CGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
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

!.... output some ua profiles

        if ( mod(i,10) .eq. 0) then
           fname = 'ddua.'
           fname = trim(fname) // trim(i2c(i))
           open(10,file=fname)
           write(10,"('# ',10(1pe20.13,1x))") one, zero, &
                            real(alphal), aimag(alphal), &
                            omega, real(dwda), aimag(dwda)
           do j = 1, ny
              write(10,"(10(1pe20.13,1x))") y(j), &
                   real(u(i,j,k,n,1)), aimag(u(i,j,k,n,1)), &
                   real(u(i,j,k,n,2)), aimag(u(i,j,k,n,2)), &
                   real(u(i,j,k,n,3)), aimag(u(i,j,k,n,3)), &
                   real(u(i,j,k,n,4)), aimag(u(i,j,k,n,4))
           end do
        end if

      end do   ! loop on i

      return
      end subroutine duda
