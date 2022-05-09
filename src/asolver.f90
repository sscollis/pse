!==============================================================================
      subroutine asolver
!
!     Solve the continuous linear adjoint PSE using Chebyshev collocation
!
!     NOTES:  1) Has the option to read in alpha from prior PSE run
!
!     S. Scott Collis
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
      complex :: kex, umax(nx), dumax
      real    :: ke(nx)

      real, allocatable :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
                           E(:,:,:), Ex(:,:,:)

      complex, allocatable :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:)

      real :: c1(ny), c2(ny)
      logical :: trap=.false.

      real :: rdum, alphar, alphai
      integer :: idum
      logical :: pse_alpha = .true.
!==============================================================================
      rei = one / re

!.... Read in alpha from a regular PSE run (assumes same mesh)

      if (pse_alpha) then
        open(20,err=100,file='pse.in')
        do i = is, ie
          read(20,*) rdum, alphar, alphai
          alpha(i,1,1) = cmplx(alphar, alphai)
        end do
        close(20)
      end if

      alpha = -alpha
      omega = -omega
      beta  = -beta

!.... Read in alpha from adjoint LNS run

      if (.false.) then
        open(20,file='alpha.out',status='old')
        do i = is, ie
          read(20,*) idum, rdum, alphar, alphai
          alpha(i,1,1) = cmplx(alphar, alphai)
        end do
        close(20)
      end if

      allocate( u(nx,ny,nz,nt,ndof) )

      u = czero

      allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
                C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
                Ex(ny,ndof,ndof) )

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

      ke = zero
      kex = zero
      umax = zero
      dumax = zero

      ul = u(i,:,k,n,1)
      vl = u(i,:,k,n,2)
      wl = u(i,:,k,n,3)

      if (trap) then
        ke(i) = zero
        do j = 1, ny-1
          ke(i) = ke(i) + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                                  conjg(wl(j))*wl(j) + &
                            conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                            conjg(wl(j+1))*wl(j+1) ) * (y(j)-y(j-1))
        end do
      else
        ke(i) = zero
        do j = 1, ny
          ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
            abs(wl(j))**2)
        end do
      end if

      if (norm.eq.0) then
      else if (norm.eq.1) then
        call findumax( ny, y, ul, umax(i), nint, yint )
      else if (norm.eq.2) then
        ke(i) = zero
        do j = 1, ny
          ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
        end do  
      else 
        call error('ASolver$','Illegal value for norm$')
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

      open(9,file='adj.dat')
      write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                                  amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
      call flush(9)

!.... main marching loop

      do i = ie-1, is, -1

        cpu2 = second()
        write(*,"(/,'  i       x(i)          cpu        total cpu')")
        write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
        cpu = cpu2

        if (.not. pse_alpha) alpha(i,k,n) = alpha(i+1,k,n)

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
        amp(i,k,n) = exp( -iota * amp(i,k,n) )  ! to integrate backwards [SSC]

        if (trap) then
          ke(i) = zero
          do j = 1, ny-1
            ke(i) = ke(i) + pt5 * ( conjg(ul(j))*ul(j) + conjg(vl(j))*vl(j) + &
                                    conjg(wl(j))*wl(j) + &
                          conjg(ul(j+1))*ul(j+1) + conjg(vl(j+1))*vl(j+1) + &
                          conjg(wl(j+1))*wl(j+1) ) * (y(j+1)-y(j))
          end do
        else
          ke(i) = zero
          do j = 1, ny
            ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
                                    abs(wl(j))**2)
          end do
        end if

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
        
!.... compute the mode amplitude

        amp(i,k,n) = czero
        do il = ie, i+1, -1
          amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il-1,k,n)) * &
                       (x(il)-x(il-1))
        end do
        amp(i,k,n) = exp( -iota * amp(i,k,n) )  ! to integrate backwards [SSC]

!.... compute the streamwise derivative of alpha

        dalpha(k,n) = ( alpha(i+1,k,n) - alpha(i,k,n) ) / ( x(i+1) - x(i) )

!.... Bertollotti sets this to zero (no advantage)

!       dalpha(k,n) = zero

        G=zero; A=zero; B=zero; C=zero; D=zero; E=zero; Ex=zero
        Gb=czero; Ab=czero

        G(:,1,1) = one
        G(:,2,2) = one
        G(:,3,3) = one

        A(:,1,1) = c1(:) * ub(:,i)
        A(:,1,2) = rei * two * c1(:) * c2(:)
        A(:,1,4) = c1(:)
        A(:,2,1) = -rei * two * c1(:) * c2(:)
        A(:,2,2) = c1(:) * ub(:,i)
        A(:,3,3) = c1(:) * ub(:,i)
        A(:,4,1) = c1(:)

        B(:,1,1) = vb(:,i) + rei * c2(:)
        B(:,2,2) = vb(:,i) + rei * c2(:)
        B(:,2,4) = one
        B(:,3,3) = vb(:,i) + rei * c2(:)
        B(:,4,2) = one

        C(:,1,1) = wb(:,i)
        C(:,2,2) = wb(:,i)
        C(:,3,3) = wb(:,i)
        C(:,3,4) = one
        C(:,4,3) = one

        D(:,1,1) = -(c1(:) * ubx(:,i) + c2(:) * vb(:,i)) + rei * c2(:)**2
        D(:,1,2) = -c1(:) * vbx(:,i) + two * c2(:) * ub(:,i)
        D(:,1,3) = -c1(:) * wbx(:,i)
        D(:,2,1) = -uby(:,i)
        D(:,2,2) = -vby(:,i) - rei * c2(:)**2
        D(:,2,3) = -wby(:,i)
        D(:,4,2) = c2(:)
      
        Ex(:,1,1) = -c1(:)**2 * rei
        Ex(:,2,2) = -c1(:)**2 * rei
        Ex(:,3,3) = -c1(:)**2 * rei
      
        E(:,1,1) = -rei
        E(:,2,2) = -rei
        E(:,3,3) = -rei

        Gb = -iota * kt(n) * omega * G + iota * alpha(i,k,n) * A + &
              iota * kz(k) * beta * C + D - &
              Ex * ( iota * dalpha(k,n) - alpha(i,k,n)**2 ) + &
              kz(k)**2 * beta**2 * E

!.... fix the streamwise pressure gradient

        if (ipfix.eq.1 .or. (ipfix.eq.2 .and. i.eq.2) ) A(:,:,4) = zero

        Ab = A - two * iota * alpha(i,k,n) * Ex

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

!.... the minus sign is here because x(i) is the unknown on the LHS

              A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Gb(l,ldof,mdof) - &
                                    Ab(l,ldof,mdof) / (x(i+1)-x(i))
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

!.... Build the RHS

        r = czero

        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            do mdof = 1, ndof

!.... the minus sign here is becuase i+1 is the known on the RHS

              r(l0+ldof) = r(l0+ldof) - Ab(l,ldof,mdof) * &
                           u(i+1,l,k,n,mdof) / (x(i+1)-x(i))
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

        ul = u(i,:,k,n,1)
        vl = u(i,:,k,n,2)
        wl = u(i,:,k,n,3)
        pl = u(i,:,k,n,4)

!.... compute some streamwise derivatives

        ux = (u(i+1,:,k,n,1) - ul)/(x(i+1)-x(i))
        vx = (u(i+1,:,k,n,2) - vl)/(x(i+1)-x(i))
        wx = (u(i+1,:,k,n,3) - wl)/(x(i+1)-x(i))

!.... compute the correction for alpha (following Herbert AIAA-93-3053)

        if (trap) then
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
        else
          ke(i) = zero
          kex   = zero
          do j = 1, ny
            ke(i) = ke(i) + opi(j)*( abs(ul(j))**2 + abs(vl(j))**2 + &
              abs(wl(j))**2 )
            kex = kex + opi(j)*( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
              conjg(wl(j))*wx(j) )
          end do
        endif

!.... update alpha

        if (.false.) then
        if (norm.eq.0) then       ! following Herbert AIAA-93-3053
          alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                              sor * (alpha(i,k,n) - iota * kex / ke(i))
        else if (norm.eq.1) then  ! following Bertolotti
          call findumax( ny, y, ul, umax(i), nint, yint )
          dumax = (umax(i+1)-umax(i))/(x(i+1)-x(i))
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
        end if

!.... update the mode amplitude

        amp(i,k,n) = czero
        do il = ie, i+1, -1
          amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il-1,k,n)) * &
                       (x(il)-x(il-1))
        end do
        amp(i,k,n) = exp( -iota * amp(i,k,n) )  ! to integrate backwards [SSC]

!.... convergence check

        if (norm.eq.0) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), abs(kex/ke(i))/abs(alpha(i,k,n))
!         if (abs(kex/ke(i))/abs(alpha(i,k,n)) .lt. tol) exit
          if (abs(kex/ke(i)) .lt. tol) exit
        else if (norm.eq.1) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(umax(i)), &
                abs(dumax/umax(i))/abs(alpha(i,k,n))
!         if (abs(dumax/umax(i))/abs(alpha(i,k,n)) .lt. tol) exit
          if (abs(dumax/umax(i)) .lt. tol) exit
        else if (norm.eq.2) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), abs(kex/ke(i))/abs(alpha(i,k,n))
!         if (abs(kex/ke(i))/abs(alpha(i,k,n)) .lt. tol) exit
          if (abs(kex/ke(i)) .lt. tol) exit
        end if

        end do    ! loop on iter

!       write(*,*) 'i, x(i), alpha(i,k,n)', i, x(i), alpha(i,k,n)

!=============================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=============================================================================

        write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), &
                                    abs(kex/ke(i)), amp(i,k,n), &
                                    ke(i)*abs(amp(i,k,n))**2
        call flush(9)

        ul = u(i,:,k,n,1)
        vl = u(i,:,k,n,2)
        wl = u(i,:,k,n,3)
        pl = u(i,:,k,n,4)
        if (ipwbc .eq. 1) pl(1) = pl(2) - (pl(3)-pl(2))/(y(3)-y(2))*y(2)
!       write(99,"(8(1pe20.13,1x))") x(i), abs(amp(i,k,n)*pl(1))

      end do   ! loop on i

!.... output the final profile

      i = is
      ul = u(i,:,k,n,1)
      vl = u(i,:,k,n,2)
      wl = u(i,:,k,n,3)
      pl = u(i,:,k,n,4)
      if (ipwbc .eq. 1) pl(1) = pl(2) - (pl(3)-pl(2))/(y(3)-y(2))*y(2)
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

      call post(-1)

!.... return back to the regular definition just in case needed elsewhere

      alpha = -alpha
      omega = -omega
      beta  = -beta

      return

100   call error('asolver$','Error reading pse.in$')

      end subroutine asolver
