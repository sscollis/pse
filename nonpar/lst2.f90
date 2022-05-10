!==============================================================================
   subroutine lst
!
!     Solves the spatial eigenvalue problem for a parallel base flow
!
!     8.28.99  local solutions are added to reduce transient
!
!                              .----.   @   @
!                             / .-"-.`.  \v/
!                             | | '\ \ \_/ )
!                           ,-\ `-.' /.'  /
!                          '---`----'----'
!
!     S. Scott Collis
!     Alex Dobrinsky
!==============================================================================
   use global
   use fmax
   implicit none
   
   integer        :: i, j, k, n, iter
   integer        :: l, ldof, l0, m, mdof, m0, ievec=1, ier
   
   real           :: rei
   integer        :: ipvt(2*ndof*ny)
   complex        :: scale
   
   real           :: G(ny,ndof,ndof),  A(ny,ndof,ndof),  B(ny,ndof,ndof)
   real           :: C(ny,ndof,ndof),  D(ny,ndof,ndof),  E(ny,ndof,ndof)
   real           :: Ex(ny,ndof,ndof), F(ny,ndof,ndof),  H(ny,ndof,ndof)
   real           :: AB(ny,ndof,ndof), DB(ny,ndof,ndof), BB(ny,ndof,ndof)
   real           :: CB(ny,ndof,ndof), hm = 0.0
   
   complex        :: Dh(ny,ndof,ndof), D2(2*ndof*ny,2*ndof*ny)
   complex        :: A0(4*ndof*ny,4*ndof*ny), B0(4*ndof*ny,4*ndof*ny)
   complex        :: D0(2*ndof*ny,2*ndof*ny), D1(2*ndof*ny,2*ndof*ny)
   complex        :: evec(4*ndof*ny,4*ndof*ny), alp(4*ndof*ny)
   complex        :: omg(4*ndof*ny), cs(4*ndof*ny), bet(4*ndof*ny)
   real           :: temp1(4*ndof*ny), temp2(4*ndof*ny)
   integer        :: index(4*ndof*ny)
   
   real           :: c1(ny), c2(ny)

   character*1    :: trans = 'n'
   real           :: ubxy(ny,nx), vbxy(ny,nx), wbxy(ny,nx)

   complex        :: alphal, ke, kex, ul(ny), vl(ny), wl(ny), pl(ny), &
                                      ux(ny), vx(ny), wx(ny), px(ny)

!.... stuff for LAPACK eigensolver

   integer                   :: info, lwork
   complex, allocatable      :: work(:)
   real,    allocatable      :: rwork(:)

   real :: yumax
   complex :: umax, dumax
!==============================================================================

   rei = one / re

!.... Do for first station only (assume parallel flow)

   i = is
   write(*,"('LST for station ',i4,' s = ',1pe13.6)") i, x(i)

!.... compute cross derivatives of the mean flow

#ifdef CRAY
   call sgemm (trans, trans, ny, nx, ny, one, opy, ny, ubx, ny, &
               zero, ubxy, ny)
   call sgemm (trans, trans, ny, nx, ny, one, opy, ny, vbx, ny, &
               zero, vbxy, ny)
   call sgemm (trans, trans, ny, nx, ny, one, opy, ny, wbx, ny, &
               zero, wbxy, ny)
#else
   call dgemm (trans, trans, ny, nx, ny, one, opy, ny, ubx, ny, &
               zero, ubxy, ny)
   call dgemm (trans, trans, ny, nx, ny, one, opy, ny, vbx, ny, &
               zero, vbxy, ny)
   call dgemm (trans, trans, ny, nx, ny, one, opy, ny, wbx, ny, &
               zero, wbxy, ny)
#endif

!.... write out the mean profile

   open(50,file='mean.out')
   do j = 1, ny
     write(50,"(14(1pe13.6,1x))") y(j), &
                                  ub(j,i), ubx(j,i), uby(j,i), ubxy(j,i), &
                                  vb(j,i), vbx(j,i), vby(j,i), vbxy(j,i), &
                                  wb(j,i), wbx(j,i), wby(j,i), wbxy(j,i)
   end do
   close(50)	

!==============================================================================

   c1(:) = one/(one + cur(i)*y(:))
   c2(:) = cur(i) * c1(:)

   G=zero; A=zero; B=zero; C=zero; D=zero; E=zero; Ex=zero
   F=zero; H=zero; AB=zero; BB=zero; DB=zero; CB=zero
   
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
   
   B(:,1,1) = hm * vb(:,i) - rei * c2(:)
   B(:,2,2) = hm * vb(:,i) - rei * c2(:)
   B(:,2,4) = one
   B(:,3,3) = hm * vb(:,i) - rei * c2(:)
   B(:,4,2) = one
   
   C(:,1,1) = wb(:,i)
   C(:,2,2) = wb(:,i)
   C(:,3,3) = wb(:,i)
   C(:,3,4) = one
   C(:,4,3) = one
   
!.... added  c2(:) * ub(:,i) to the D(:,1,2). also I add
!.... the wall normal velocities, par. hm can be set to 0

   D(:,1,1) = -rei * c2(:)**2 + hm * c2(:) * vb(:,i) 
   D(:,1,2) = uby(:,i) + c2(:) * ub(:,i)  
   D(:,2,1) = -two * c2(:) * ub(:,i) 
   D(:,2,2) = rei * c2(:)**2 + hm * vby(:,i)
   D(:,3,2) = wby(:,i)
   D(:,4,2) = c2(:)
   
   Ex(:,1,1) = c1(:)**2 * rei
   Ex(:,2,2) = c1(:)**2 * rei
   Ex(:,3,3) = c1(:)**2 * rei
   
   E(:,1,1) = rei
   E(:,2,2) = rei
   E(:,3,3) = rei

   if (.true.) then

!.... adds the x derivative terms to the 1st equation 
 
   H(:,1,1) = c1(:) * ubx(:,i)
   H(:,2,1) = c1(:) * vbx(:,i)
   H(:,3,1) = c1(:) * wbx(:,i)

!.... Base flow expansion: 2nd equation, with no alpha, (u_0).

   DB(:,1,1) =  c2(:) * vbx(:,i) 
   DB(:,1,2) =  ubxy(:,i) + c2(:) * ubx(:,i)
   DB(:,2,1) = -two * c2(:) * ubx(:,i)
   DB(:,2,2) =  vbxy(:,i)
   DB(:,3,2) =  wbxy(:,i)

   BB(:,1,1) =  vbx(:,i)
   BB(:,2,2) =  vbx(:,i)
   BB(:,3,3) =  vbx(:,i)

   CB(:,1,1) =  wbx(:,i)
   CB(:,2,2) =  wbx(:,i)
   CB(:,3,3) =  wbx(:,i)

!.... Base flow expansion:  alpha terms. 2nd equation, (u_0) 

   AB(:,1,1) =  c1(:) * ubx(:,i)
   AB(:,2,2) =  c1(:) * ubx(:,i)
   AB(:,3,3) =  c1(:) * ubx(:,i)

   end if

!==============================================================================

!.... initialize

      D0 = czero
      D1 = czero
      D2 = czero

!.... D0:   no  alpha,  {u_0},  1st equation  

      Dh = D + H + iota * beta * C + beta**2 * E - iota * omega * G
      
      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    do m = 1, ny
	      m0 = (m-1)*ndof
	      D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + &
		                    B(l,ldof,mdof) * opy(l,m) - &
		                    E(l,ldof,mdof) * opyy(l,m)
	    end do
	    m0 = (l-1)*ndof
	    D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + Dh(l,ldof,mdof)
	  end do
	end do
      end do

      if (.true.) then

!....  no alpha, {u_1},  1st equation 

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    m0 = ndof*ny + (l-1)*ndof
	    D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + A(l,ldof,mdof)
	  end do
	end do
      end do

!....  no alpha, {u_0},  2nd equation 

      Dh = DB + iota * beta * CB
      
      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof + ndof*ny
	    do m = 1, ny
	      m0 = (m-1)*ndof
	      D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + &
		                    BB(l,ldof,mdof) * opy(l,m)
	    end do
	    m0 = (l-1)*ndof 
	    D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + Dh(l,ldof,mdof)
	  end do
	end do
      end do

      end if

!....  no alpha, {u_1}, 2nd equation

      Dh = D + iota * beta * C + beta**2 * E - iota * omega * G
      
      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof + ndof*ny
	    do m = 1, ny
	      m0 = (m-1)*ndof + ndof*ny
	      D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + &
		                    B(l,ldof,mdof) * opy(l,m) - &
		                    E(l,ldof,mdof) * opyy(l,m)
	    end do
	    m0 = (l-1)*ndof + ndof*ny
	    D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) + Dh(l,ldof,mdof)
	  end do
	end do
      end do
 
!.... Boundary conditions to D0           | 1 | 2 |
!.... numbering of quarters:         D0 = ---------  
!                                         | 3 | 4 |

!.... BC first quarter, beginning

      l = 1
      l0 = (l-1)*ndof

      D0(l0+1,:) = czero
      D0(l0+1,l0+1) = cone

      D0(l0+2,:) = czero
      D0(l0+2,l0+2) = cone

      D0(l0+3,:) = czero
      D0(l0+3,l0+3) = cone
      
      D0(l0+4,:) = czero
      if (ipwbc .eq. 0) then
        do m = 1, ny
          m0 = (m-1)*ndof
          D0(l0+4,m0+4) = opy(l,m)
        end do
      else if (ipwbc .eq. 1) then
        do m = 1, ny
          m0 = (m-1)*ndof
          D0(l0+4,m0+2) = -rei * c2(l) * opy(l,m) - rei * opyy(l,m)
          D0(l0+4,m0+4) = opy(l,m)
        end do
      else if (ipwbc .eq. 2) then   ! This makes D0 singular
        write(*,*) 'WARNING:  This makes D0 singular'
        do m = 1, ny
	  m0 = (m-1)*ndof
	  D0(l0+4,m0+2) = opy(l,m)
	end do
      else
	call error('lst$','Illegal value of ipwbc$')
      end if
      
!.... BC first quarter, end

      l = ny
      l0 = (l-1)*ndof

      D0(l0+1,:) = czero
      D0(l0+1,l0+1) = cone

      D0(l0+2,:) = czero
      if (ivbc .eq. 0) then
	D0(l0+2,l0+2) = cone
      else if (ivbc .eq. 1) then
	do m = 1, ny
	  m0 = (m-1)*ndof
	  D0(l0+2,m0+2) = opy(l,m)
	end do
      else
	call error('lst$','Illegal value of ivbc$')
      end if

      D0(l0+3,:) = czero
      D0(l0+3,l0+3) = cone

      D0(l0+4,:) = czero
      D0(l0+4,l0+4) = cone	

!.... BC forth quarter, beginning

      l = 1
      l0 = (l+ny-1)*ndof

      D0(l0+1,:) = czero
      D0(l0+1,l0+1) = cone
      
      D0(l0+2,:) = czero
      D0(l0+2,l0+2) = cone
      
      D0(l0+3,:) = czero
      D0(l0+3,l0+3) = cone
      
      D0(l0+4,:) = czero
      if (ipwbc .eq. 0) then
        do m = 1, ny
          m0 = (m-1)*ndof + ndof*ny
          D0(l0+4,m0+4) = opy(l,m)
        end do
      else if (ipwbc .eq. 1) then
        do m = 1, ny
          m0 = (m-1)*ndof + ndof*ny
          D0(l0+4,m0+2) = -rei * c2(l) * opy(l,m) - rei * opyy(l,m)
          D0(l0+4,m0+4) = opy(l,m)
        end do
      else if (ipwbc .eq. 2) then   ! This makes D0 singular
        write(*,*) 'WARNING:  This makes D0 singular'
        do m = 1, ny
	  m0 = (m-1)*ndof + ndof*ny
	  D0(l0+4,m0+2) = opy(l,m)
	end do
      else
	call error('lst$','Illegal value of ipwbc$')
      end if

!.... BC forth quarter, end

      l = ny
      l0 = (l+ny-1)*ndof

      D0(l0+1,:) = czero
      D0(l0+1,l0+1) = cone

      D0(l0+2,:) = czero
      if (ivbc .eq. 0) then
	D0(l0+2,l0+2) = cone
      else if (ivbc .eq. 1) then
	do m = 1, ny
	  m0 = (m-1)*ndof + ndof*ny
	  D0(l0+2,m0+2) = opy(l,m)
	end do
      else
	call error('lst$','Illegal value of ivbc$')
      end if

      D0(l0+3,:) = czero
      D0(l0+3,l0+3) = cone
      
      D0(l0+4,:) = czero
      D0(l0+4,l0+4) = cone	

!.... D1:  terms that are proportional to alpha, u_{0} 1st equation

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    m0 = (l-1)*ndof 
	    D1(l0+ldof,m0+mdof) = D1(l0+ldof,m0+mdof) + iota * A(l,ldof,mdof)
	  end do
	end do
      end do

      if (.true.) then

!.... D1:  terms that are proportional to alpha, u_{1} 1st equation

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    m0 = (l-1)*ndof + ndof*ny
	    D1(l0+ldof,m0+mdof) = D1(l0+ldof,m0+mdof) - &
                                  two * iota * Ex(l,ldof,mdof)
	  end do
	end do
      end do

      end if

!.... D1:  terms that are proportional to alpha, u_{0}, 2nd equation

      if (.true.) then   ! these are the terms that lead to the wiggles

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof + ndof*ny
	    m0 = (l-1)*ndof 
	    D1(l0+ldof,m0+mdof) = D1(l0+ldof,m0+mdof) + iota * AB(l,ldof,mdof)
	  end do
	end do
      end do

      end if

!.... D1:  terms that are proportional to alpha, u_{1}, 2nd equation

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof + ndof*ny
	    m0 = (l-1)*ndof + ndof*ny
	    D1(l0+ldof,m0+mdof) = D1(l0+ldof,m0+mdof) + iota * A(l,ldof,mdof)
	  end do
	end do
      end do

!.... Apply the boundary conditions to D1

      l = 1
      l0 = (l-1)*ndof
      
      D1(l0+1,:) = czero
      D1(l0+2,:) = czero
      D1(l0+3,:) = czero
      D1(l0+4,:) = czero
      
      l = ny
      l0 = (l-1)*ndof
      
      D1(l0+1,:) = czero
      D1(l0+2,:) = czero
      D1(l0+3,:) = czero
      D1(l0+4,:) = czero
      
      l = 1
      l0 = (l+ny-1)*ndof

      D1(l0+1,:) = czero
      D1(l0+2,:) = czero
      D1(l0+3,:) = czero
      D1(l0+4,:) = czero
      
      l = ny
      l0 = (l+ny-1)*ndof
      
      D1(l0+1,:) = czero
      D1(l0+2,:) = czero
      D1(l0+3,:) = czero
      D1(l0+4,:) = czero

!.... D2:  terms that are proportional to alpha^2, u_{0} 1st equation
      
      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof
	    m0 = (l-1)*ndof
	    D2(l0+ldof,m0+mdof) = D2(l0+ldof,m0+mdof) + Ex(l,ldof,mdof)
	  end do
	end do
      end do

!.... D2:  terms that are proportional to alpha^2,  u_{1} 1st equaiton
!....                       -------- NONE --------

!.... D2:  terms that are proportional to alpha^2,  u_{0} 2st equaiton
!....                       -------- NONE --------

!.... D2:  terms that are proportional to alpha^2,  u_{1} 2st equaiton

      do ldof = 1, ndof
	do mdof = 1, ndof
	  do l = 1, ny
	    l0 = (l-1)*ndof + ndof*ny
	    m0 = (l-1)*ndof + ndof*ny
	    D2(l0+ldof,m0+mdof) = D2(l0+ldof,m0+mdof) + Ex(l,ldof,mdof)
	  end do
	end do
      end do

!.... Apply the boundary conditions to D2

      l = 1
      l0 = (l-1)*ndof
      
      D2(l0+1,:) = czero
      D2(l0+2,:) = czero
      D2(l0+3,:) = czero
      D2(l0+4,:) = czero
      
      l = ny
      l0 = (l-1)*ndof
      
      D2(l0+1,:) = czero
      D2(l0+2,:) = czero
      D2(l0+3,:) = czero
      D2(l0+4,:) = czero

      l = 1
      l0 = (l+ny-1)*ndof

      D2(l0+1,:) = czero
      D2(l0+2,:) = czero
      D2(l0+3,:) = czero
      D2(l0+4,:) = czero
      
      l = ny
      l0 = (l+ny-1)*ndof
      
      D2(l0+1,:) = czero
      D2(l0+2,:) = czero
      D2(l0+3,:) = czero
      D2(l0+4,:) = czero

!.... form the extended system
      
      A0   = zero
      B0   = zero
      evec = zero
      alp  = zero
      bet  = zero
      omg  = zero

!.... Eigensolution method 1

#ifdef CRAY
      call CGETRF(2*ndof*ny, 2*ndof*ny, D0, 2*ndof*ny, ipvt, info)
      if (info.ne.0) then
	write(*,*) 'ERROR in CGETRF: ',info
	call exit(1)
      end if
#else
      call ZGETRF(2*ndof*ny, 2*ndof*ny, D0, 2*ndof*ny, ipvt, info)
      if (info.ne.0) then
	write(*,*) 'ERROR in ZGETRF: ',info
	call exit(1)
      end if
#endif

      D1 = -D1
      D2 = -D2
      
#ifdef CRAY
      call CGETRS('N', 2*ndof*ny, 2*ndof*ny, D0, 2*ndof*ny, &
	           ipvt, D1, 2*ndof*ny, info)
      if (info.ne.0) then
	write(*,*) 'CGETRS: ',info
	call exit(1)
      end if
      call CGETRS('N', 2*ndof*ny, 2*ndof*ny, D0, 2*ndof*ny, &
	          ipvt, D2, 2*ndof*ny, info)
      if (info.ne.0) write(*,*) 'CGETRS: ',info
#else
      call ZGETRS('N', 2*ndof*ny, 2*ndof*ny, D0, 2*ndof*ny, &
	          ipvt, D1, 2*ndof*ny, info)
      if (info.ne.0) then
	write(*,*) 'ZGETRS: ',info
	call exit(1)
      end if
      call ZGETRS('N', 2*ndof*ny, 2*ndof*ny, D0, 2*ndof*ny, &
	          ipvt, D2, 2*ndof*ny, info)
      if (info.ne.0) write(*,*) 'ZGETRS: ',info
#endif
      
      B0(1:2*ndof*ny,1:2*ndof*ny)           = D1
      B0(1:2*ndof*ny,2*ndof*ny+1:4*ndof*ny) = D2
      
!.... put in the identity matrices
      
      do l = 2*ndof*ny+1, 4*ndof*ny
	B0(l,l-2*ndof*ny) = one
      end do
      
!.... LAPACK regular complex eigensolver
      
      lwork = 16*(4*ndof*ny)
      allocate (work(lwork), rwork(2*4*ndof*ny), STAT=ier)
      if (ier .ne. 0) then
	write(*,*) 'Error allocating work space'
	call exit(1)
      end if
      cpu2 = second()
      write(*,"('Solving Eigensystem ',2(1x,1pe10.3))") cpu2-cpu, cpu2
      cpu = cpu2
#ifdef CRAY
      if (ievec.eq.1) then
	call CGEEV('N', 'V', 4*ndof*ny, B0, 4*ndof*ny, alp, evec, &
	            4*ndof*ny, evec, 4*ndof*ny, work, lwork, rwork, info)
      else
	call CGEEV('N', 'N', 4*ndof*ny, B0, 4*ndof*ny, alp, evec, &
	            4*ndof*ny, evec, 4*ndof*ny, work, lwork, rwork, info)
      end if
#else
      if (ievec.eq.1) then
	call ZGEEV('N', 'V', 4*ndof*ny, B0, 4*ndof*ny, alp, evec, &
	            4*ndof*ny, evec, 4*ndof*ny, work, lwork, rwork, info)
      else
	call ZGEEV('N', 'N', 4*ndof*ny, B0, 4*ndof*ny, alp, evec, &
	            4*ndof*ny, evec, 4*ndof*ny, work, lwork, rwork, info)
      end if
#endif
      cpu2 = second()
      write(*,"('Completed Eigensolution ',2(1x,1pe10.3))") cpu2-cpu, cpu2
      cpu = cpu2
      if (info.ne.0) then
	if (info.le.0) then
	  write(*,*) 'Error in aurgument: ',abs(info)
	else
	  write(*,*) 'Warning: Error computing eigenvalues: Error # ',info
	end if
      end if
      deallocate (work, rwork)

!.... Note that I have to invert the eigenvalue since I defined the system
!.... Backwards from Bridges and Morris!

      where (alp .ne. zero) 
	alp = one / alp
      elsewhere
	alp = zero
      end where

!.... sort the eigenvalues by the imaginary part
      
      do j = 1, 4*ndof*ny
	temp2(j) = aimag(alp(j))
	index(j) = j
      end do
      call PIKSR2(4*ndof*ny, temp2, index)
      do j = 1, 4*ndof*ny
	temp1(j) = real(alp(index(j)))
	A0(:,j) = evec(:,index(j))
      end do
      
      alp(1:4*ndof*ny) = cmplx(temp1(1:4*ndof*ny),temp2(1:4*ndof*ny))
      evec(:,1:4*ndof*ny) = A0(:,1:4*ndof*ny)

!.... compute the phase speed

      where (alp .ne. zero)
	cs = omega / alp
      elsewhere
	cs = zero
      end where

!.... Scale the eigenvectors in a reasonable way

      if (ievec.eq.1) then
	do j = 1, 4*ndof*ny
	  scale = zero
	  do l = 1, 2*ny*ndof
	    if ( abs(evec(l,j)) .gt. abs(scale) ) then
	      scale = evec(l,j)
	    end if
	  end do
	  if (abs(scale) .gt. 1.0e-10) then
	    do l = 1, 4*ny*ndof
	      evec(l,j) = evec(l,j) / scale
	    end do
	  end if
	end do
      end if

!.... output the eigenvalues and eigenfunctions to the terminal

      do j = 1, 4*ndof*ny
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
 
      if ( j .lt. 0 .or. j .gt. 4*ndof*ny ) goto 100
      if (j .ne. 0) then

	do l = 1, ny
	  l0 = (l-1)*ndof
          ul(l) = evec(l0+1,j)
          vl(l) = evec(l0+2,j)
          wl(l) = evec(l0+3,j)
          pl(l) = evec(l0+4,j)
	  l0 = (l+ny-1)*ndof
          ux(l) = evec(l0+1,j)
          vx(l) = evec(l0+2,j)
          wx(l) = evec(l0+3,j)
          px(l) = evec(l0+4,j)
        end do

!	ke  = zero
!	kex = zero
!	do l = 1, ny
!	  ke = ke + opi(l)*( abs(ul(l))**2 + abs(vl(l))**2 + &
!               abs(wl(l))**2 )
!	  kex = kex + opi(l)*( conjg(ul(l))*ux(l) + conjg(vl(l))*vx(l) + &
!                conjg(wl(l))*wx(l) )
!	end do
!        
!        alphal = alp(j) - iota * kex / ke

	call findumax( ny, y, ul, umax, nint, yint, yumax )
        dumax = cmplx( getval( ny, y,  real(ux), yumax ), &
                       getval( ny, y, aimag(ux), yumax ) )

	alphal = alp(j) - iota * dumax / umax

	open (unit=20,file='phi0.out',form='formatted',status='unknown')
	write(20,55) one, zero, real(alphal), aimag(alphal)
	do l = 1, ny
	  l0 = (l-1)*ndof
	  write (20,50) y(l), &
			real(evec(l0+1,j)), &
			aimag(evec(l0+1,j)), &
			real(evec(l0+2,j)), &
			aimag(evec(l0+2,j)), &
			real(evec(l0+3,j)), &
			aimag(evec(l0+3,j)), &
			real(evec(l0+4,j)), &
			aimag(evec(l0+4,j))
	end do
	close (20)

	open (unit=20,file='phi1.out',form='formatted',status='unknown')
	write(20,55) one, zero, real(alp(j)), aimag(alp(j))
	do l = 1, ny
	  l0 = (l+ny-1)*ndof
	  write (20,50) y(l), &
			real(evec(l0+1,j)), &
			aimag(evec(l0+1,j)), &
			real(evec(l0+2,j)), &
			aimag(evec(l0+2,j)), &
			real(evec(l0+3,j)), &
			aimag(evec(l0+3,j)), &
			real(evec(l0+4,j)), &
			aimag(evec(l0+4,j))
	end do
	close (20)

!.... Sanity check (hope it doesn't fail...)

	open (unit=20,file='phi0a.out',form='formatted',status='unknown')
	write(20,55) one, zero, real(alp(j)), aimag(alp(j))
	do l = 1, ny
	  l0 = (l+2*ny-1)*ndof
	  write (20,50) y(l), &
			real(evec(l0+1,j)/ alp(j))  , &
			aimag(evec(l0+1,j) / alp(j)), &
			real(evec(l0+2,j) / alp(j)) , &
			aimag(evec(l0+2,j) / alp(j)), &
			real(evec(l0+3,j) / alp(j)) , &
			aimag(evec(l0+3,j) / alp(j)), &
			real(evec(l0+4,j) / alp(j)) , &
			aimag(evec(l0+4,j) / alp(j))
	end do
	close (20)

	open (unit=20,file='phi1a.out',form='formatted',status='unknown')
	write(20,55) one, zero, real(alp(j)), aimag(alp(j))
	do l = 1, ny
	  l0 = (l+3*ny-1)*ndof
	  write (20,50) y(l), &
			real(evec(l0+1,j) / alp(j)), &
			aimag(evec(l0+1,j) / alp(j)), &
			real(evec(l0+2,j) / alp(j)), &
			aimag(evec(l0+2,j) / alp(j)), &
			real(evec(l0+3,j) / alp(j)), &
			aimag(evec(l0+3,j) / alp(j)), &
			real(evec(l0+4,j) / alp(j)), &
			aimag(evec(l0+4,j) / alp(j))
	end do
	close (20)

	goto 100
      end if

      return
20    format(1p,i4,1x,e13.6,1x,2(e13.6,1x),2(e13.6,1x))
25    format(1p,i4,1x,e13.6,1x,2(e13.6,1x),2(e13.6,1x),' <==')
50    format(1p,11(e20.13,1x))
55    format('# ',1p,11(e20.13,1x))
      end
