!==============================================================================
      subroutine ddalpha
!
!     Solves the alpha variation equation for a parallel base flow
!
!     S. Scott Collis
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n
      integer :: l, ldof, l0, m, mdof, m0, ier

      real    :: rei

      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)
      real    :: rdum, ur, ui, vr, vi, wr, wi, ppr, ppi, alphar, alphai
      character*1 :: cdum
      complex :: Jprod, L2prod

      integer :: ipvt(ndof*ny), info
      complex :: evec(ny,ndof), r(ny*ndof), alp, aalp, avec(ny,ndof)
      complex :: r0(ny*ndof)

      real    :: G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
                 C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
                 Ex(ny,ndof,ndof)

      complex :: Ab(ny,ndof,ndof), A0(ndof*ny,ndof*ny), AA(ndof*ny,ndof*ny)

      real    :: c1(ny), c2(ny), temp(ndof), epsilon = 1.0e-10

      complex :: dwda, rdot

      integer :: lwork
      complex, allocatable :: work(:), UU(:,:), VV(:,:)
      real, allocatable :: rwork(:,:), SS(:)
!==============================================================================

!.... Read in an eigenfunction file

      open(10,file='space.reg')
      read(10,*) cdum, rdum, rdum, alphar, alphai
      alp = cmplx(alphar,alphai)
      do j = 1, ny
	read(10,*) rdum, ur, ui, vr, vi, wr, wi, ppr, ppi
	evec(j,1) = cmplx(ur,ui)
	evec(j,2) = cmplx(vr,vi)
	evec(j,3) = cmplx(wr,wi)
	evec(j,4) = cmplx(ppr,ppi)
      end do
      close(10)

!.... Read in an adjoint eigenfunction file

      open(10,file='space.adj')
      read(10,*) cdum, rdum, rdum, alphar, alphai
      aalp = cmplx(alphar,alphai)
      do j = 1, ny
	read(10,*) rdum, ur, ui, vr, vi, wr, wi, ppr, ppi
	avec(j,1) = cmplx(ur,ui)
	avec(j,2) = cmplx(vr,vi)
	avec(j,3) = cmplx(wr,wi)
	avec(j,4) = cmplx(ppr,ppi)
      end do
      close(10)

      rei = one / re

!.... Do for first station only (assume parallel flow)

      i = is
      write(*,"('Spatial DDalpha for station ',i4,' s = ',1pe13.6)") i, x(i)

!.... compute the J inner product

      Jprod = czero
      do j = 1, ny
	Jprod = Jprod + opi(j)*( &
              (evec(j,1)*avec(j,1)+evec(j,2)*avec(j,2))*ub(j,i) + &
	       evec(j,1)*avec(j,4) + avec(j,1)*evec(j,4) - &
	       rei*iota*(alp+aalp)*(evec(j,1)*avec(j,1)+evec(j,2)*avec(j,2)) )
      end do

!.... compute the L2 inner product of the velocity

      L2prod = czero
      do j = 1, ny
	L2prod = L2prod + opi(j)*( evec(j,1)*avec(j,1) + &
                                   evec(j,2)*avec(j,2) + evec(j,3)*avec(j,3) )
      end do

      dwda = Jprod/L2prod
      write(*,*) 'i, x(i) = ', i, x(i)
      write(*,*) 'Alphah = ', alp
      write(*,*) 'Alphat = ', aalp
      write(*,*) 'Jprod = ', Jprod
      write(*,*) 'L2prod = ', L2prod
      write(*,*) 'Domega/Dalpha = ', dwda

!.... normalize the adjoint eigenfunction using the J inner product

      avec = avec / Jprod

!.... recompute the J inner product as a sanity check

      Jprod = czero
      do j = 1, ny
	Jprod = Jprod + opi(j)*( &
              (evec(j,1)*avec(j,1)+evec(j,2)*avec(j,2))*ub(j,i) + &
	       evec(j,1)*avec(j,4) + avec(j,1)*evec(j,4) - &
	       rei*iota*(alp+aalp)*(evec(j,1)*avec(j,1)+evec(j,2)*avec(j,2)) )
      end do

!.... compute the L2 inner product of the velocity

      L2prod = czero
      do j = 1, ny
	L2prod = L2prod + opi(j)*( evec(j,1)*avec(j,1) + &
                                   evec(j,2)*avec(j,2) + evec(j,3)*avec(j,3) )
      end do

      dwda = Jprod/L2prod
      write(*,*) 'After normalization'
      write(*,*) 'Jprod = ', Jprod
      write(*,*) 'L2prod = ', L2prod
      write(*,*) 'Domega/Dalpha = ', dwda

!.... write out the mean profile

      open(50,file='mean.out')
      do j = 1, ny
	write(50,"(10(1pe13.6,1x))") y(j), ub(j,i), uby(j,i), wb(j,i), wby(j,i)
      end do
      close(50)	

      c1(:) = one/(one + cur(i)*y(:))
      c2(:) = cur(i) * c1(:)

      G=zero; A=zero; B=zero; C=zero; D=zero; E=zero; Ex=zero

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

      B(:,1,1) = -rei * c2(:)
      B(:,2,2) = -rei * c2(:)
      B(:,2,4) = one
      B(:,3,3) = -rei * c2(:)
      B(:,4,2) = one

      C(:,1,1) = wb(:,i)
      C(:,2,2) = wb(:,i)
      C(:,3,3) = wb(:,i)
      C(:,3,4) = one
      C(:,4,3) = one

      D(:,1,1) = -rei * c2(:)**2
      D(:,1,2) = uby(:,i)
      D(:,2,1) = -two * c2(:) * ub(:,i)
      D(:,2,2) = rei * c2(:)**2
      D(:,3,2) = wby(:,i)
      D(:,4,2) = c2(:)
      
      Ex(:,1,1) = c1(:)**2 * rei
      Ex(:,2,2) = c1(:)**2 * rei
      Ex(:,3,3) = c1(:)**2 * rei
      
      E(:,1,1) = rei
      E(:,2,2) = rei
      E(:,3,3) = rei

!==============================================================================

!.... initialize

      A0 = czero
      
      Ab = iota * alp * A + alp**2 * Ex + D + iota * beta * C + &
           beta**2 * E - iota * omega * G

!.... form the LHS matrix

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
	    A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Ab(l,ldof,mdof)
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
	  A0(l0+4,m0+2) = rei * c2(l) * opy(l,m) + rei * opyy(l,m)
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
      
!.... form the RHS vector

      Ab = iota * A + two * alp * Ex - iota * dwda * G

      r = czero

      temp(1) = one
      temp(2) = one
      temp(3) = one
      temp(4) = zero

      do ldof = 1, ndof
	do l = 1, ny
	  l0 = (l-1)*ndof
	  do mdof = 1, ndof
	    r(l0+ldof) = r(l0+ldof) - Ab(l,ldof,mdof) * &
	                 evec(l,mdof)
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

      AA = A0  ! keep around old copies
      r0 = r
      
!.... Solve the system

#ifdef CRAY
      call CGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
      if (info.ne.0) call error('solver$','Singular matrix$')
      call CGESL(A0,ny*ndof,ny*ndof,ipvt,r,0)
#else
      lwork = 32*4*ny*ndof
      allocate( work(lwork), rwork(3*ny*ndof,5*ny*ndof-4), SS(ny*ndof), &
                UU(ny*ndof,ny*ndof), VV(ny*ndof,ny*ndof) )

      call ZGESVD('A','A',ny*ndof,ny*ndof,A0,ny*ndof,SS,UU,ny*ndof,VV, &
                  ny*ndof,work,lwork,rwork,info)
      write(*,*) "info = ", info, real(work(1)), lwork
      open(22,file='svd.dat')
      do l = 1, ny*ndof
        write(22,"(i5,2(1x,1pe14.6E3))") l, SS(l)
      end do
      close(22)

      VV = conjg(transpose(VV))  ! ZGESVD return VV^H

!.... make a solution that is minimal in the 2-norm using the SVD

      r = czero
      do j = 1, ny*ndof
        if (SS(j) .gt. epsilon) &
             r = r + dot_product(UU(:,j),r0(:)) * VV(:,j) / SS(j)
      end do
#endif

!.... check the solution

      write(*,*) 'Solution tolerance = ', maxval(abs(matmul(AA,r) - r0))

!.... convert to local variables

      do l = 1, ny
	l0 = (l-1)*ndof
	ul(l) = r(l0+1)
	vl(l) = r(l0+2)
	wl(l) = r(l0+3)
	pl(l) = r(l0+4)
      end do

!.... compute the J inner product

      Jprod = czero
      do j = 1, ny
	Jprod = Jprod + opi(j)*( &
              (ul(j)*avec(j,1)+vl(j)*avec(j,2))*ub(j,i) + &
	       ul(j)*avec(j,4) + avec(j,1)*pl(j) - &
	       rei*iota*(alp+aalp)*(ul(j)*avec(j,1)+vl(j)*avec(j,2)) )
      end do
      write(*,*) 'Jprod (pre-scaling) = ', Jprod

!.... project out the homogeneous solution

      ul(:) = ul(:) - Jprod * evec(:,1)
      vl(:) = vl(:) - Jprod * evec(:,2)
      wl(:) = wl(:) - Jprod * evec(:,3)
      pl(:) = pl(:) - Jprod * evec(:,4)

!.... compute the new J inner product to be double sure

      Jprod = czero
      do j = 1, ny
	Jprod = Jprod + opi(j)*( &
              (ul(j)*avec(j,1)+vl(j)*avec(j,2))*ub(j,i) + &
	       ul(j)*avec(j,4) + avec(j,1)*pl(j) - &
	       rei*iota*(alp+aalp)*(ul(j)*avec(j,1)+vl(j)*avec(j,2)) )
      end do
      write(*,*) 'Jprod (post-scaling) = ', Jprod

!.... check the solution tolerance again to be double sure

      do l = 1, ny
	l0 = (l-1)*ndof
	r(l0+1) = ul(l)
        r(l0+2) = vl(l)
	r(l0+3) = wl(l)
	r(l0+4) = pl(l)
      end do
      write(*,*) 'Solution tolerance = ',maxval(abs(matmul(AA,r) - r0))

!.... output the solution, d/dalpha including domega/dalpha

      open(10,file='ddalpha.dat')
      write(10,"('# ',10(1pe20.13,1x))") one, zero, real(alp), aimag(alp), &
                                         omega, real(dwda), aimag(dwda)
      do j = 1, ny
	write(10,"(10(1pe20.13,1x))") y(j), &
                                      real(ul(j)), aimag(ul(j)), &
                                      real(vl(j)), aimag(vl(j)), &
                                      real(wl(j)), aimag(wl(j)), &
                                      real(pl(j)), aimag(pl(j))
      end do
      close(10)

!.... output the eigenfunction and normalized adjoint

      open(10,file='evec.dat')
      do j = 1, ny
	write(10,"(10(1pe20.13,1x))") y(j), &
                                      real(evec(j,1)), aimag(evec(j,1)), &
                                      real(evec(j,2)), aimag(evec(j,2)), &
                                      real(evec(j,3)), aimag(evec(j,3)), &
                                      real(evec(j,4)), aimag(evec(j,4))
      end do
      close(10)

      open(10,file='avec.dat')
      do j = 1, ny
	write(10,"(10(1pe20.13,1x))") y(j), &
                                      real(avec(j,1)), aimag(avec(j,1)), &
                                      real(avec(j,2)), aimag(avec(j,2)), &
                                      real(avec(j,3)), aimag(avec(j,3)), &
                                      real(avec(j,4)), aimag(avec(j,4))
      end do
      close(10)

      return
      end subroutine ddalpha
