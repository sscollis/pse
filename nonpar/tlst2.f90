!==============================================================================
      subroutine tlst
!
!     Solve the temporal eigenvalue problem
!
!     alpha = 2.3001668687353E-001  -6.7942137704083E-003
!
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, iter, i0
      integer :: l, ldof, l0, m, mdof, m0, ievec=1, ier

      real    :: rei
      integer :: ipvt(ndof*ny)
      complex :: scale

      real    :: G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
                 C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
                 Ex(ny,ndof,ndof)

      complex :: Dh(ny,ndof,ndof)
      complex :: A0(ndof*ny,ndof*ny), B0(ndof*ny,ndof*ny)
      complex :: C0(ndof*ny,ndof*ny), D0(ndof*ny,ndof*ny)
      complex :: evec(ndof*ny,ndof*ny), alp(ndof*ny), bet(ndof*ny)
      complex :: omg(ndof*ny), cs(ndof*ny)
      real    :: temp1(ndof*ny), temp2(ndof*ny), p(ny)
      integer :: index(ndof*ny)

      real    :: c1(ny), c2(ny)

!.... stuff for LAPACK eigensolver

      integer :: info, lwork
      complex, allocatable :: work(:)
      real, allocatable    :: rwork(:)

      complex :: alph = (2.3001668687353E-001,-6.7942137704083E-003)

      complex :: r(ny,ndof)

      logical :: adjoint = .true.
!==============================================================================

      rei = one / re

!.... Do for first station only (assume parallel flow)

      i = 1

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

!     alph = alpha(i,1,1)
      write (*,"('alpha = ',e13.6,1x,e13.6,)") alph

!.... write out the mean profile

      open(50,file='mean.out')
      do j = 1, ny
	write(50,"(10(1pe13.6,1x))") y(j), ub(j,i), uby(j,i), wb(j,i), wby(j,i)
      end do
      close(50)	

!.... initialize

      A0   = zero
      B0   = zero
      evec = zero
      alp  = zero
      bet  = zero
      omg  = zero

!.... B0

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    m0 = (l-1)*ndof
	    B0(l0+ldof,m0+mdof) = iota * G(l,ldof,mdof)
	  end do
	end do
      end do
      
!.... A0

      Dh = D + iota * alph * A + iota * beta * C + alph**2 * Ex + beta**2 * E

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
	    A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Dh(l,ldof,mdof)
	  end do
	end do
      end do

!.... these are the actual weights for Chebyschev

      do i = 1, ny
	p(i) = sqrt(1.0 - (cos((i-1)*acos(-1.0)/(ny-1)))**2)
      end do
      do i = 1, ny
	dydeta(i) = dydeta(i) * p(i)
      end do
      dydeta(1) = dydeta(1) * pt5
      dydeta(ny) = dydeta(ny) * pt5

!.... compute the discrete adjoint

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do m = 1, ny
	      m0 = (m-1)*ndof
	      A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) * dydeta(l)
	      B0(l0+ldof,m0+mdof) = B0(l0+ldof,m0+mdof) * dydeta(l)
	    end do
	  end do
	end do
      end do

      A0 = transpose( A0 )
      B0 = transpose( B0 )

!.... Apply the boundary conditions to B0

      l = 1
      l0 = (l-1)*ndof

      B0(l0+1,:) = czero

      B0(l0+2,:) = czero
      
      B0(l0+3,:) = czero

      B0(l0+4,:) = czero

      l = ny
      l0 = (l-1)*ndof
      
      B0(l0+1,:) = czero
      
      B0(l0+2,:) = czero
      
      B0(l0+3,:) = czero
      
      B0(l0+4,:) = czero

!.... Apply the boundary conditions to A0

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
	call error('lst$','Illegal value of ipwbc$')
      end if

      l = ny
      l0 = (l-1)*ndof
      
      A0(l0+1,:) = czero
      A0(l0+1,l0+1) = cone

      A0(l0+2,:) = czero
      A0(l0+2,l0+2) = cone

      A0(l0+3,:) = czero
      A0(l0+3,l0+3) = cone

      A0(l0+4,:) = czero
      A0(l0+4,l0+4) = cone	

!.... Lapack generalized complex eigensolver

      lwork = 8*ndof*ny
      allocate (work(lwork), rwork(lwork), STAT=ier)
      if (ier .ne. 0) then
	write(*,*) 'Error allocating work space'
	call exit(1)
      end if

#ifdef CRAY
      call CGEGV ( 'N', 'V', ndof*ny, A0, ndof*ny, B0, ndof*ny, 	&
                   alp, bet, evec, ndof*ny, evec, ndof*ny, 		&
	           work, lwork, rwork, info )
#else		     
      call ZGEGV ( 'N', 'V', ndof*ny, A0, ndof*ny, B0, ndof*ny, 	&
                   alp, bet, evec, ndof*ny, evec, ndof*ny, 		&
		   work, lwork, rwork, info )
#endif
      if (info.ne.0) write (*,*) 'Info = ', info
	
      deallocate (work, rwork)

!.... compute the frequency (temporal)

      where (bet .ne. 0) 
	omg = alp / bet
      elsewhere
	omg = zero
      end where

      alp = omg

!.... If discrete adjoint then scale by the metric of the mapping

      if (adjoint .and. .false.) then
	do j = 1, ndof*ny
	  do i = 1, ny
	    i0 = (i-1)*ndof
	    evec(i0+1,j) = evec(i0+1,j) / dydeta(i)
	    evec(i0+2,j) = evec(i0+2,j) / dydeta(i)
	    evec(i0+3,j) = evec(i0+3,j) / dydeta(i)
	    evec(i0+4,j) = evec(i0+4,j) / dydeta(i)
	  end do
	end do
      end if

!.... sort the eigenvalues by the imaginary part
      
      do j = 1, ndof*ny
	temp2(j) = aimag(alp(j))
	index(j) = j
      end do
      call PIKSR2(ndof*ny, temp2, index)
      do j = 1, ndof*ny
	temp1(j) = real(alp(index(j)))
	A0(:,j) = evec(:,index(j))
      end do
      
      alp(1:ndof*ny) = cmplx(temp1(1:ndof*ny),temp2(1:ndof*ny))
      evec(:,1:ndof*ny) = A0(:,1:ndof*ny)

!.... compute the phase speed

      cs = alp / sqrt(alph**2 + beta**2)

!.... Scale the eigenvectors in a reasonable way

      if (ievec.eq.1) then
	do j = 1, ndof*ny
	  scale = zero
	  do l = 1, ny*ndof
	    if ( abs(evec(l,j)) .gt. abs(scale) ) then
	      scale = evec(l,j)
	    end if
	  end do
	  if (scale .ne. zero) then
	    do l = 1, ny*ndof
	      evec(l,j) = evec(l,j) / scale
	    end do
	  end if
	end do
      end if

!.... output the eigenvalues and eigenfunctions to the terminal

      do j = 1, ndof*ny
	if ( aimag(cs(j)) .gt. zero) then
	  write (*,25) j, abs(alp(j)), real(alp(j)), aimag(alp(j)), &
	               real(cs(j)), aimag(cs(j))
	else
	  write (*,20) j, abs(alp(j)), real(alp(j)), aimag(alp(j)), &
		       real(cs(j)), aimag(cs(j))
	end if
      end do
100   continue

      write (*,"(/,'Which eigenfunction ==> ',$)")
      read (*,*) j
 
      if ( j .lt. 0 .or. j .gt. ndof*ny ) goto 100
      if (j .ne. 0) then
	open (unit=20,file='time.out',form='formatted',status='unknown')

        if (ipwbc .eq. 1) then
	  evec(4,j) = evec(8,j) - (evec(12,j)-evec(8,j))/(y(3)-y(2))*y(2)
	end if

	do i = 1, ny
	  i0 = (i-1)*ndof
	  write (20,50) y(i), &
			real(evec(i0+1,j)), &
			aimag(evec(i0+1,j)), &
			real(evec(i0+2,j)), &
			aimag(evec(i0+2,j)), &
			real(evec(i0+3,j)), &
			aimag(evec(i0+3,j)), &
			real(evec(i0+4,j)), &
			aimag(evec(i0+4,j))
	end do
	close (20)
	goto 100
      end if

      return
20    format(1p,i4,1x,e13.6,1x,2(e13.6,1x),2(e13.6,1x))
25    format(1p,i4,1x,e13.6,1x,2(e13.6,1x),2(e13.6,1x),' <==')
50    format(1p,11(e20.13,1x))
      end
