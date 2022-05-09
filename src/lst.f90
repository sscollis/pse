!==============================================================================
      subroutine lst(adjoint)
!
!     Solves the spatial eigenvalue problem for a parallel base flow
!
!     S. Scott Collis
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, iter, id
      integer :: l, ldof, l0, m, mdof, m0, ievec=1, ier
      logical :: adjoint

      real    :: rei
      integer :: ipvt(ndof*ny)
      complex :: scale

      real    :: G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
                 C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
                 Ex(ny,ndof,ndof)

      complex :: Dh(ny,ndof,ndof)
      complex :: A0(2*ndof*ny,2*ndof*ny), B0(2*ndof*ny,2*ndof*ny)
      complex :: D0(ndof*ny,ndof*ny), D1(ndof*ny,ndof*ny), D2(ndof*ny,ndof*ny)
      complex :: evec(2*ndof*ny,2*ndof*ny), alp(2*ndof*ny), bet(2*ndof*ny)
      complex :: omg(2*ndof*ny), cs(2*ndof*ny)
      real    :: temp1(2*ndof*ny), temp2(2*ndof*ny), p(ny)
      integer :: index(2*ndof*ny)

      real    :: c1(ny), c2(ny)

!.... stuff for LAPACK eigensolver

      integer :: info, lwork
      complex, allocatable :: work(:)
      real, allocatable    :: rwork(:)

      logical :: eliminate = .false.   ! set to true to eliminate BC's

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

      dydeta = opi   ! this is the integration operator

      rei = one / re

!.... Do for first station only (assume parallel flow)

      if(adjoint) then
        i = ie
        write(*,"('Spatial ALST for station ',i4,' s = ',1pe13.6)") i, x(i)
      else
        i = is
        write(*,"('Spatial LST for station ',i4,' s = ',1pe13.6)") i, x(i)
      end if

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

      D0 = czero
      D1 = czero
      D2 = czero

!.... D1:  terms that are proportional to alpha

      do ldof = 1, ndof
        do mdof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            m0 = (l-1)*ndof
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
      
!.... D2:  terms that are proportional to alpha^2
      
      do ldof = 1, ndof
        do mdof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            m0 = (l-1)*ndof
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

!.... D0:  all terms with no factor of alpha

      Dh = D + iota * beta * C + beta**2 * E - iota * omega * G
      
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

!.... Apply the boundary conditions to D0

      l = 1
      l0 = (l-1)*ndof

      D0(l0+1,:) = czero
      D0(:,l0+1) = czero
      D0(l0+1,l0+1) = cone

      D0(l0+2,:) = czero
      D0(:,l0+2) = czero
      D0(l0+2,l0+2) = cone

      D0(l0+3,:) = czero
      D0(:,l0+3) = czero
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

!.... statically eliminate the Neumann boundary conditions

      if (eliminate) then 
      if (ipwbc .eq. 0 .or. ipwbc .eq. 1) then

        do l = 2, ny
          l0 = (l-1)*ndof
          do m = 2, ny
            m0 = (m-1)*ndof
            D0(l0+1,m0+4) = D0(l0+1,m0+4) - D0(l0+1,4) * D0(4,m0+4) / D0(4,4)
            D0(l0+2,m0+4) = D0(l0+2,m0+4) - D0(l0+2,4) * D0(4,m0+4) / D0(4,4)
            D0(l0+3,m0+4) = D0(l0+3,m0+4) - D0(l0+3,4) * D0(4,m0+4) / D0(4,4)
            D0(l0+4,m0+4) = D0(l0+4,m0+4) - D0(l0+4,4) * D0(4,m0+4) / D0(4,4)
            
            D1(l0+1,m0+4) = D1(l0+1,m0+4) - D0(l0+1,4) * D1(4,m0+4) / D0(4,4)
            D1(l0+2,m0+4) = D1(l0+2,m0+4) - D0(l0+2,4) * D1(4,m0+4) / D0(4,4)
            D1(l0+3,m0+4) = D1(l0+3,m0+4) - D0(l0+3,4) * D1(4,m0+4) / D0(4,4)
            D1(l0+4,m0+4) = D1(l0+4,m0+4) - D0(l0+4,4) * D1(4,m0+4) / D0(4,4)

            D2(l0+1,m0+4) = D2(l0+1,m0+4) - D0(l0+1,4) * D2(4,m0+4) / D0(4,4)
            D2(l0+2,m0+4) = D2(l0+2,m0+4) - D0(l0+2,4) * D2(4,m0+4) / D0(4,4)
            D2(l0+3,m0+4) = D2(l0+3,m0+4) - D0(l0+3,4) * D2(4,m0+4) / D0(4,4)
            D2(l0+4,m0+4) = D2(l0+4,m0+4) - D0(l0+4,4) * D2(4,m0+4) / D0(4,4)
          end do
        end do

        l = 1
        l0 = (l-1)*ndof
        do m = 2, ny
          m0 = (m-1)*ndof
          D0(l0+1,m0+4) = D0(l0+1,m0+4) - D0(l0+1,4) * D0(4,m0+4) / D0(4,4)
          D0(l0+2,m0+4) = D0(l0+2,m0+4) - D0(l0+2,4) * D0(4,m0+4) / D0(4,4)
          D0(l0+3,m0+4) = D0(l0+3,m0+4) - D0(l0+3,4) * D0(4,m0+4) / D0(4,4)
          
          D1(l0+1,m0+4) = D1(l0+1,m0+4) - D0(l0+1,4) * D1(4,m0+4) / D0(4,4)
          D1(l0+2,m0+4) = D1(l0+2,m0+4) - D0(l0+2,4) * D1(4,m0+4) / D0(4,4)
          D1(l0+3,m0+4) = D1(l0+3,m0+4) - D0(l0+3,4) * D1(4,m0+4) / D0(4,4)

          D2(l0+1,m0+4) = D2(l0+1,m0+4) - D0(l0+1,4) * D2(4,m0+4) / D0(4,4)
          D2(l0+2,m0+4) = D2(l0+2,m0+4) - D0(l0+2,4) * D2(4,m0+4) / D0(4,4)
          D2(l0+3,m0+4) = D2(l0+3,m0+4) - D0(l0+3,4) * D2(4,m0+4) / D0(4,4)
        end do

        D0(4,:) = D0(4,:) / D0(4,4)
        D0(:,4) = zero
        D0(4,4) = one

        do l = 2, ny
          l0 = (l-1)*ndof
          do m = 1, ny
            m0 = (m-1)*ndof
            D0(l0+4,m0+4) = D0(l0+4,m0+4) + opy(1,l) * D0(4,m0+4) / dydeta(l)
            D1(l0+4,m0+4) = D1(l0+4,m0+4) + opy(1,l) * D1(4,m0+4) / dydeta(l)
            D2(l0+4,m0+4) = D2(l0+4,m0+4) + opy(1,l) * D2(4,m0+4) / dydeta(l)
          end do
        end do
        
        l = 1
        l0 = (l-1)*ndof
        do m = 1, ny
          m0 = (m-1)*ndof
          D0(l0+4,m0+4) = D0(l0+4,m0+4) * opy(1,l) / dydeta(l)
          D1(l0+4,m0+4) = D1(l0+4,m0+4) * opy(1,l) / dydeta(l)
          D2(l0+4,m0+4) = D2(l0+4,m0+4) * opy(1,l) / dydeta(l)
        end do

      else if (ipwbc .eq. 2) then
          
        do l = 2, ny
          l0 = (l-1)*ndof
          do m = 2, ny
            m0 = (m-1)*ndof
            D0(l0+1,m0+2) = D0(l0+1,m0+2) - D0(l0+1,2) * D0(4,m0+2) / D0(4,2)
            D0(l0+2,m0+2) = D0(l0+2,m0+2) - D0(l0+2,2) * D0(4,m0+2) / D0(4,2)
            D0(l0+3,m0+2) = D0(l0+3,m0+2) - D0(l0+3,2) * D0(4,m0+2) / D0(4,2)
            D0(l0+4,m0+2) = D0(l0+4,m0+2) - D0(l0+4,2) * D0(4,m0+2) / D0(4,2)
            
            D1(l0+1,m0+2) = D1(l0+1,m0+2) - D0(l0+1,2) * D1(4,m0+2) / D0(4,2)
            D1(l0+2,m0+2) = D1(l0+2,m0+2) - D0(l0+2,2) * D1(4,m0+2) / D0(4,2)
            D1(l0+3,m0+2) = D1(l0+3,m0+2) - D0(l0+3,2) * D1(4,m0+2) / D0(4,2)
            D1(l0+4,m0+2) = D1(l0+4,m0+2) - D0(l0+4,2) * D1(4,m0+2) / D0(4,2)

            D2(l0+1,m0+2) = D2(l0+1,m0+2) - D0(l0+1,2) * D2(4,m0+2) / D0(4,2)
            D2(l0+2,m0+2) = D2(l0+2,m0+2) - D0(l0+2,2) * D2(4,m0+2) / D0(4,2)
            D2(l0+3,m0+2) = D2(l0+3,m0+2) - D0(l0+3,2) * D2(4,m0+2) / D0(4,2)
            D2(l0+4,m0+2) = D2(l0+4,m0+2) - D0(l0+4,2) * D2(4,m0+2) / D0(4,2)
          end do
        end do
        
        l = 1
        l0 = (l-1)*ndof
        do m = 2, ny
          m0 = (m-1)*ndof
          D0(l0+1,m0+2) = D0(l0+1,m0+2) - D0(l0+1,2) * D0(4,m0+2) / D0(4,2)
          D0(l0+2,m0+2) = D0(l0+2,m0+2) - D0(l0+2,2) * D0(4,m0+2) / D0(4,2)
          D0(l0+3,m0+2) = D0(l0+3,m0+2) - D0(l0+3,2) * D0(4,m0+2) / D0(4,2)
          
          D1(l0+1,m0+2) = D1(l0+1,m0+2) - D0(l0+1,2) * D1(4,m0+2) / D0(4,2)
          D1(l0+2,m0+2) = D1(l0+2,m0+2) - D0(l0+2,2) * D1(4,m0+2) / D0(4,2)
          D1(l0+3,m0+2) = D1(l0+3,m0+2) - D0(l0+3,2) * D1(4,m0+2) / D0(4,2)

          D2(l0+1,m0+2) = D2(l0+1,m0+2) - D0(l0+1,2) * D2(4,m0+2) / D0(4,2)
          D2(l0+2,m0+2) = D2(l0+2,m0+2) - D0(l0+2,2) * D2(4,m0+2) / D0(4,2)
          D2(l0+3,m0+2) = D2(l0+3,m0+2) - D0(l0+3,2) * D2(4,m0+2) / D0(4,2)
        end do
        
        D0(4,:) = D0(4,:) / D0(4,2)
        D0(:,2) = zero
        D0(4,2) = one
        
        do l = 2, ny
          l0 = (l-1)*ndof
          do m = 1, ny
            m0 = (m-1)*ndof
            D0(l0+4,m0+2) = D0(l0+4,m0+2) + opy(1,l) * D0(4,m0+2) / dydeta(l)
            D1(l0+4,m0+2) = D1(l0+4,m0+2) + opy(1,l) * D1(4,m0+2) / dydeta(l)
            D2(l0+4,m0+2) = D2(l0+4,m0+2) + opy(1,l) * D2(4,m0+2) / dydeta(l)
          end do
        end do
        
        l = 1
        l0 = (l-1)*ndof
        do m = 1, ny
          m0 = (m-1)*ndof
          D0(l0+4,m0+2) = D0(l0+4,m0+2) * opy(1,l) / dydeta(l)
          D1(l0+4,m0+2) = D1(l0+4,m0+2) * opy(1,l) / dydeta(l)
          D2(l0+4,m0+2) = D2(l0+4,m0+2) * opy(1,l) / dydeta(l)
        end do

      else
        call error('lst$','Illegal value of ipwbc$')
      end if
      end if

      l = 1
      l0 = (l-1)*ndof

      D0(l0+1,:) = czero
      D0(:,l0+1) = czero
      D0(l0+1,l0+1) = cone

      D0(l0+2,:) = czero
      D0(:,l0+2) = czero
      D0(l0+2,l0+2) = cone

      D0(l0+3,:) = czero
      D0(:,l0+3) = czero
      D0(l0+3,l0+3) = cone

!.... far-field boundary conditions

      l = ny
      l0 = (l-1)*ndof
      
      D0(l0+1,:) = czero
      D0(:,l0+1) = czero
      D0(l0+1,l0+1) = cone
      
      D0(l0+2,:) = czero
      if (ivbc .eq. 0) then
        D0(:,l0+2) = czero
        D0(l0+2,l0+2) = cone
      else if (ivbc .eq. 1) then
        do m = 1, ny
          m0 = (m-1)*ndof
          D0(l0+2,m0+2) = opy(l,m)
        end do

!.... statically eliminate the boundary condition

        if (eliminate) then

        id = (ny-1)*ndof+2
        do l = 1, ny-1
          l0 = (l-1)*ndof
          do m = 1, ny-1
            m0 = (m-1)*ndof
            D0(l0+1,m0+2) = D0(l0+1,m0+2) - D0(l0+1,id)*D0(id,m0+2) / D0(id,id)
            D0(l0+2,m0+2) = D0(l0+2,m0+2) - D0(l0+2,id)*D0(id,m0+2) / D0(id,id)
            D0(l0+3,m0+2) = D0(l0+3,m0+2) - D0(l0+3,id)*D0(id,m0+2) / D0(id,id)
            D0(l0+4,m0+2) = D0(l0+4,m0+2) - D0(l0+4,id)*D0(id,m0+2) / D0(id,id)
            
            D1(l0+1,m0+2) = D1(l0+1,m0+2) - D0(l0+1,id)*D1(id,m0+2) / D0(id,id)
            D1(l0+2,m0+2) = D1(l0+2,m0+2) - D0(l0+2,id)*D1(id,m0+2) / D0(id,id)
            D1(l0+3,m0+2) = D1(l0+3,m0+2) - D0(l0+3,id)*D1(id,m0+2) / D0(id,id)
            D1(l0+4,m0+2) = D1(l0+4,m0+2) - D0(l0+4,id)*D1(id,m0+2) / D0(id,id)

            D2(l0+1,m0+2) = D2(l0+1,m0+2) - D0(l0+1,id)*D2(id,m0+2) / D0(id,id)
            D2(l0+2,m0+2) = D2(l0+2,m0+2) - D0(l0+2,id)*D2(id,m0+2) / D0(id,id)
            D2(l0+3,m0+2) = D2(l0+3,m0+2) - D0(l0+3,id)*D2(id,m0+2) / D0(id,id)
            D2(l0+4,m0+2) = D2(l0+4,m0+2) - D0(l0+4,id)*D2(id,m0+2) / D0(id,id)
          end do
        end do

        l = ny
        l0 = (l-1)*ndof
        do m = 1, ny-1
          m0 = (m-1)*ndof
          D0(l0+1,m0+2) = D0(l0+1,m0+2) - D0(l0+1,id) * D0(id,m0+2) / D0(id,id)
          D0(l0+2,m0+2) = D0(l0+2,m0+2) - D0(l0+2,id) * D0(id,m0+2) / D0(id,id)
          D0(l0+3,m0+2) = D0(l0+3,m0+2) - D0(l0+3,id) * D0(id,m0+2) / D0(id,id)
          
          D1(l0+1,m0+2) = D1(l0+1,m0+2) - D0(l0+1,id) * D1(id,m0+2) / D0(id,id)
          D1(l0+2,m0+2) = D1(l0+2,m0+2) - D0(l0+2,id) * D1(id,m0+2) / D0(id,id)
          D1(l0+3,m0+2) = D1(l0+3,m0+2) - D0(l0+3,id) * D1(id,m0+2) / D0(id,id)

          D2(l0+1,m0+2) = D2(l0+1,m0+2) - D0(l0+1,id) * D2(id,m0+2) / D0(id,id)
          D2(l0+2,m0+2) = D2(l0+2,m0+2) - D0(l0+2,id) * D2(id,m0+2) / D0(id,id)
          D2(l0+3,m0+2) = D2(l0+3,m0+2) - D0(l0+3,id) * D2(id,m0+2) / D0(id,id)
        end do

        D0(id,:) = D0(id,:) / D0(id,id)
        D0(:,id) = zero
        D0(id,id) = one

        do l = 1, ny-1
          l0 = (l-1)*ndof
          do m = 1, ny
            m0 = (m-1)*ndof
            D0(l0+2,m0+2) = D0(l0+2,m0+2) + opy(ny,l) * D0(id,m0+2) / dydeta(l)
            D1(l0+2,m0+2) = D1(l0+2,m0+2) + opy(ny,l) * D1(id,m0+2) / dydeta(l)
            D2(l0+2,m0+2) = D2(l0+2,m0+2) + opy(ny,l) * D2(id,m0+2) / dydeta(l)
          end do
        end do
        
        l = ny
        l0 = (l-1)*ndof
        do m = 1, ny
          m0 = (m-1)*ndof
          D0(l0+2,m0+2) = D0(l0+2,m0+2) * opy(ny,l) / dydeta(l)
          D1(l0+2,m0+2) = D1(l0+2,m0+2) * opy(ny,l) / dydeta(l)
          D2(l0+2,m0+2) = D2(l0+2,m0+2) * opy(ny,l) / dydeta(l)
        end do

        end if
      else
        call error('lst$','Illegal value of ivbc$')
      end if
      
      l = ny
      l0 = (l-1)*ndof

      D0(l0+3,:) = czero
      D0(:,l0+3) = czero
      D0(l0+3,l0+3) = cone
      
      D0(l0+4,:) = czero
      D0(:,l0+4) = czero
      D0(l0+4,l0+4) = cone      
      
!.... form the extended system
      
      A0   = zero
      B0   = zero
      evec = zero
      alp  = zero
      bet  = zero
      omg  = zero

!.... Eigensolution method 1

      if (.true.) then

      if (adjoint) then        ! include mapping in adjoint
        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            do m = 1, ny
              do mdof = 1, ndof
                m0 = (m-1)*ndof
                D0(l0+ldof,m0+mdof) = D0(l0+ldof,m0+mdof) * dydeta(l)
                D1(l0+ldof,m0+mdof) = D1(l0+ldof,m0+mdof) * dydeta(l)
                D2(l0+ldof,m0+mdof) = D2(l0+ldof,m0+mdof) * dydeta(l)
              end do
            end do
          end do
        end do
        
        D0 = transpose( D0 )

        call ZGETRF(ndof*ny, ndof*ny, D0, ndof*ny, ipvt, info)
        if (info.ne.0) then
          write(*,*) 'ERROR in ZGETRF: ',info
          call exit(1)
        end if

        D1 = -transpose( D1 )
        B0(ndof*ny+1:2*ndof*ny,1:ndof*ny) = -transpose( D2 )

        D2 = zero
        do l = 1, ndof*ny
          D2(l,l) = one
        end do

        call ZGETRS('N', ndof*ny, ndof*ny, D0, ndof*ny, &
                    ipvt, D1, ndof*ny, info)
        if (info.ne.0) then
          write(*,*) 'ZGETRS: ',info
          call exit(1)
        end if
        call ZGETRS('N', ndof*ny, ndof*ny, D0, ndof*ny, &
                    ipvt, D2, ndof*ny, info)
        if (info.ne.0) write(*,*) 'ZGETRS: ',info

        B0(1:ndof*ny,1:ndof*ny)           = D1
        B0(1:ndof*ny,ndof*ny+1:2*ndof*ny) = D2

      else

#ifdef CRAY
        call CGETRF(ndof*ny, ndof*ny, D0, ndof*ny, ipvt, info)
        if (info.ne.0) then
          write(*,*) 'ERROR in CGETRF: ',info
          call exit(1)
        end if
#else
        call ZGETRF(ndof*ny, ndof*ny, D0, ndof*ny, ipvt, info)
        if (info.ne.0) then
          write(*,*) 'ERROR in ZGETRF: ',info
          call exit(1)
        end if
#endif

        D1 = -D1
        D2 = -D2
      
#ifdef CRAY
        call CGETRS('N', ndof*ny, ndof*ny, D0, ndof*ny, &
                    ipvt, D1, ndof*ny, info)
        if (info.ne.0) then
          write(*,*) 'CGETRS: ',info
          call exit(1)
        end if
        call CGETRS('N', ndof*ny, ndof*ny, D0, ndof*ny, &
                    ipvt, D2, ndof*ny, info)
        if (info.ne.0) write(*,*) 'CGETRS: ',info
#else
        call ZGETRS('N', ndof*ny, ndof*ny, D0, ndof*ny, &
                    ipvt, D1, ndof*ny, info)
        if (info.ne.0) then
          write(*,*) 'ZGETRS: ',info
          call exit(1)
        end if
        call ZGETRS('N', ndof*ny, ndof*ny, D0, ndof*ny, &
                    ipvt, D2, ndof*ny, info)
        if (info.ne.0) write(*,*) 'ZGETRS: ',info
#endif
      
        B0(1:ndof*ny,1:ndof*ny)           = D1
        B0(1:ndof*ny,ndof*ny+1:2*ndof*ny) = D2
      
!.... put in the identity matrices
      
        do l = ndof*ny+1, 2*ndof*ny
          B0(l,l-ndof*ny) = one
        end do
      
      end if  ! adjoint

!.... LAPACK regular complex eigensolver
      
      lwork = 16*(2*ndof*ny)
      allocate (work(lwork), rwork(2*2*ndof*ny), STAT=ier)
      if (ier .ne. 0) then
        write(*,*) 'Error allocating work space'
        call exit(1)
      end if
      cpu2 = second()
      write(*,"('Solving Eigensystem ',2(1x,1pe10.3))") cpu2-cpu, cpu2
      cpu = cpu2
#ifdef CRAY
      if (ievec.eq.1) then
        call CGEEV('N', 'V', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                    2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
      else
        call CGEEV('N', 'N', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                    2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
      end if
#else
      if (ievec.eq.1) then
        call ZGEEV('N', 'V', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                    2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
      else
        call ZGEEV('N', 'N', 2*ndof*ny, B0, 2*ndof*ny, alp, evec, &
                    2*ndof*ny, evec, 2*ndof*ny, work, lwork, rwork, info)
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

      else

!.... Eigensolution method 2

!     A0(1:ndof*ny,1:ndof*ny)           =  D0
!     B0(1:ndof*ny,1:ndof*ny)           = -D1
!     B0(1:ndof*ny,ndof*ny+1:2*ndof*ny) = -D2 
!     do i = ndof*ny+1, 2*ndof*ny
!       A0(i,i) = one
!       B0(i,i-ndof*ny) = one
!     end do

      A0(1:ndof*ny,1:ndof*ny)           =  D2
      B0(1:ndof*ny,1:ndof*ny)           = -D1
      B0(1:ndof*ny,ndof*ny+1:2*ndof*ny) = -D0 
      do i = ndof*ny+1, 2*ndof*ny
        A0(i,i) = one
        B0(i,i-ndof*ny) = one
      end do

!.... Lapack generalized complex eigensolver

      lwork = 8*2*ndof*ny
      allocate (work(lwork), rwork(lwork), STAT=ier)
      if (ier .ne. 0) then
        write(*,*) 'Error allocating work space'
        call exit(1)
      end if
      
      cpu2 = second()
      write(*,"('Solving Eigensystem ',2(1x,1pe10.3))") cpu2-cpu, cpu2
      cpu = cpu2
#ifdef CRAY
      call CGEGV ( 'N', 'V', 2*ndof*ny, A0, 2*ndof*ny, B0, 2*ndof*ny, &
                    alp, bet, evec, 2*ndof*ny, evec, 2*ndof*ny,       &
                    work, lwork, rwork, info )           
#else
      call ZGEGV ( 'N', 'V', 2*ndof*ny, A0, 2*ndof*ny, B0, 2*ndof*ny, &
                    alp, bet, evec, 2*ndof*ny, evec, 2*ndof*ny,       &
                    work, lwork, rwork, info )           
#endif
      if (info.ne.0) write (*,*) 'Info = ', info
      deallocate (work, rwork)
      cpu2 = second()
      write(*,"('Completed Eigensolution ',2(1x,1pe10.3))") cpu2-cpu, cpu2
      cpu = cpu2

!.... compute the wave number (spatial)

!     where (bet .ne. 0) 
!       alp = alp / bet
!     elsewhere
!       alp = zero
!     end where

      where (alp .ne. 0) 
        bet = bet / alp
      elsewhere
        bet = zero
      end where
      alp = bet

      end if

!.... sort the eigenvalues by the imaginary part
      
      do j = 1, 2*ndof*ny
        temp2(j) = aimag(alp(j))
        index(j) = j
      end do
      call PIKSR2(2*ndof*ny, temp2, index)
      do j = 1, 2*ndof*ny
        temp1(j) = real(alp(index(j)))
        A0(:,j) = evec(:,index(j))
      end do
      
      alp(1:2*ndof*ny) = cmplx(temp1(1:2*ndof*ny),temp2(1:2*ndof*ny))
      evec(:,1:2*ndof*ny) = A0(:,1:2*ndof*ny)

!.... compute the phase speed

      where (alp .ne. zero)
        cs = omega / alp
      elsewhere
        cs = zero
      end where

!.... Scale the eigenvectors in a reasonable way

      if (ievec.eq.1) then
        do j = 1, 2*ndof*ny
          scale = zero
          do l = ndof+1, ny*ndof  ! start at first node off wall
            if ( abs(evec(l,j)) .gt. abs(scale) ) then
              scale = evec(l,j)
            end if
          end do
          if (abs(scale) .gt. 1.0e-10) then
            do l = 1, ny*ndof
              evec(l,j) = evec(l,j) / scale
            end do
          end if
        end do
      end if

!.... output the eigenvalues and eigenfunctions to the terminal

      do j = 1, 2*ndof*ny
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
 
      if ( j .lt. 0 .or. j .gt. 2*ndof*ny ) goto 100
      if (j .ne. 0) then
        open (unit=20,file='space.out',form='formatted',status='unknown')
        write(20,55) one, zero, real(alp(j)), aimag(alp(j)), real(omega), zero

!.... You shouldn't do this if you plan to use the eigenfunction as input
!.... to PSE

!       if (ipwbc .eq. 1) then
!         evec(4,j) = evec(8,j) - (evec(12,j)-evec(8,j))/(y(3)-y(2))*y(2)
!       end if

        do l = 1, ny
          l0 = (l-1)*ndof   !  (l+ny-1)*ndof
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

        goto 100
      end if

      return
20    format(1p,i4,1x,e13.6,1x,2(e13.6,1x),2(e13.6,1x))
25    format(1p,i4,1x,e13.6,1x,2(e13.6,1x),2(e13.6,1x),' <==')
50    format(1p,11(e20.13,1x))
55    format('# ',1p,11(e20.13,1x))
      end

!**************************************************************************
      SUBROUTINE PIKSR2 (N, ARR, BRR)
!**************************************************************************
!
!     Try the simple insertion sort.
!
!**************************************************************************
      REAL ARR(N), A
      INTEGER BRR(N), B
      
      DO J = 2, N
        A = ARR(J)
        B = BRR(J)
        DO I = J-1,1,-1
          IF(ARR(I).LE.A) GOTO 10
          ARR(I+1)=ARR(I)
          BRR(I+1)=BRR(I)
        END DO
        I = 0
  10    ARR(I+1)=A
        BRR(I+1)=B
      END DO
      
      RETURN
      END
