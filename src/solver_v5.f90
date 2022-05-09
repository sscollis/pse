!==============================================================================
      subroutine solver
!
!     Solve the linear PSE using Chebyshev collocation
!
!     Note: this is the first-order accurate scheme in x
!
!     S. Scott Collis
!
!     8-21-98:  1) Replaced trapezoidal integration with Chebyschev integration
!               2) Added normalization option (norm=0) Kinetic energy
!                                             (norm=1) Umax
!               3) Added continuity at wall for pressure (ipwbc=2)
!               4) Commented interpolation for ipwbc=1 on output
!     11-2-98:  1) Added field output
!     8-12-99:  1) Added Newton update for alpha (backward Euler)
!     10-20:99: 1) Added midpoint in x
!               2) Fixed Newton for midpoint
!==============================================================================
      use global
      use int_str
      use fmax
      implicit none

      integer :: i, j, k, n, iter, il, idof
      integer :: l, ldof, l0, m, mdof, m0, info, job=0, icount=0
      integer, allocatable :: ipvt(:)

      real :: rei

      complex :: dalpha(nz,nt)
      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)
      complex :: ux(ny),  vx(ny),  wx(ny),  px(ny)
      complex :: kex
      real    :: ke(nx), yumax, err

!.... local variables

      complex :: alphal
      real    :: ubl(ny),  vbl(ny),  wbl(ny)
      real    :: ubxl(ny), vbxl(ny), wbxl(ny)
      real    :: ubyl(ny), vbyl(ny), wbyl(ny)
      complex :: uy(ny,ndof), uyy(ny,ndof), umax(nx), dumax

      real, allocatable :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
                           E(:,:,:), Ex(:,:,:)

      complex, allocatable :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:), ua(:,:,:)
      complex :: uax(ny), vax(ny), wax(ny)
      complex :: f, dfda

      integer :: nxi, nyi
      complex :: uh(nx,ny,ndof), amph(nx), alphah(nx), uha(nx,ny,ndof)
      complex :: ut(nx,ny,ndof), ampt(nx), alphat(nx), uta(nx,ny,ndof)

      complex :: jprod(nx), jprodx
      complex :: dfdah(nx), dfdat(nx)

      real :: c1(ny), c2(ny)

      real, parameter :: fact=0.5   ! time advancement factor

!.... wall BCs

      real :: bloc = 100.0, bcon = 0.0, bamp = -1.0e-1
      real :: coe, bl

      namelist /bc/ bloc, bcon, bamp

      real :: rdum, alphar, alphai
      character*1 :: cdum
      character*80 :: fname

!==============================================================================

!.... Read in alpha from a regular PSE run (assumes same mesh)

      if (norm.eq.3 .or. norm.eq.4) then

!.... read PSE field file

        open (10, file='field.pse', form='unformatted')
        read (10) nxi, nyi
        if (nxi.ne.(ie-is+1) .or. nyi.ne.ny) &
          call error('asolver$','Illegal field.pse file$')
        read (10) (((uh(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
        read (10) (rdum, i= is,ie)
        read (10) (rdum, j= 1,ny)
        read (10) (amph(i), i = is, ie)
        read (10) (((uha(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
        read (10) (alphah(i), i= is, ie)
        close (10)

!.... read APSE field file

        open (10, file='field.adj', form='unformatted')
        read (10) nxi, nyi
        if (nxi.ne.(ie-is+1) .or. nyi.ne.ny) &
          call error('solver$','Illegal field.adj file$')
        read (10) (((ut(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
        read (10) (rdum, i= is,ie)
        read (10) (rdum, j= 1,ny)
        read (10) (ampt(i), i = is, ie)
        read (10) (((uta(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
        read (10) (alphat(i), i = is, ie)
        close (10)
      end if

!.... read PSE field file

      if (ipfix.eq.2) then
        open (10, file='field.pse', form='unformatted')
        read (10) nxi, nyi
        if (nxi.ne.(ie-is+1) .or. nyi.ne.ny) &
          call error('asolver$','Illegal field.pse file$')
        read (10) (((uh(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
        read (10) (rdum, i= is,ie)
        read (10) (rdum, j= 1,ny)
        read (10) (amph(i), i = is, ie)
        read (10) (((uha(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
        read (10) (alphah(i), i= is, ie)
        close (10)
      end if

      rei = one / re

      allocate( u(nx,ny,nz,nt,ndof) )

      u = czero

      allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
                C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
                Ex(ny,ndof,ndof) )

      allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
      allocate( r(ny*ndof), ipvt(ny*ndof), ua(nx,ny,ndof) )

!.... only consider only one mode for linear analysis

      i = is; n = 1; k = 1

!.... initialize the field

      do j = 1, ny
        u(i,j,k,n,1) = ubci(j,k,n)
        u(i,j,k,n,2) = vbci(j,k,n)
        u(i,j,k,n,3) = wbci(j,k,n)
        u(i,j,k,n,4) = pbci(j,k,n)
      end do

      if (ipfix.eq.2) then
        alpha(:,k,n) = alphah(:)
        u(:,:,k,n,:) = uh
      end if

      ul = u(i,:,k,n,1)
      vl = u(i,:,k,n,2)
      wl = u(i,:,k,n,3)
      pl = u(i,:,k,n,4)

      ke = zero
      kex = zero
      umax = zero
      dumax = zero

      ke(i) = zero
      do j = 1, ny
        ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
          abs(wl(j))**2)
      end do

      if (norm.eq.0) then
      else if (norm.eq.1) then
        call findumax( ny, y, ul, umax(i), nint, yint, yumax )
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
      else if (norm.eq.3.or.norm.eq.4) then
        jprod(i) = czero
        do j = 1, ny
          jprod(i) = jprod(i) + opi(j) * amp(i,k,n) * ampt(i) * ( &
                     (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                      ut(i,j,3)*wl(j)) * ub(j,i) + &
                      ut(i,j,1)*pl(j) + ut(i,j,4)*ul(j) - &
                      rei*iota*(alpha(i,k,n)-alphat(i)) * &
                      (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ut(i,j,3)*wl(j)))
        end do
      else
        call error('Solver$','Illegal value for norm$')
      end if

      open(9,file='pse.dat')
      write(9,"(9(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                                  amp(i,k,n), ke(i)*abs(amp(i,k,n))**2, &
                                  abs(Jprod(i))
      call flush(9)

!.... main marching loop

      do i = is+1, ie

        cpu2 = second()
        write(*,"(/,'  i       x(i)          cpu        total cpu')")
        write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
        cpu = cpu2

        if (norm.eq.4 .or. ipfix.eq.2) then
          alpha(i,k,n) = alphah(i)
        else
          alpha(i,k,n) = alpha(i-1,k,n)
        end if

!       c1(:) = one/(one + cur(i)*y(:))
!       c2(:) = cur(i) * c1(:)

        c1(:) = one/(one + ((one-fact)*cur(i-1)+fact*cur(i))*y(:))
        c2(:) = ((one-fact)*cur(i-1)+fact*cur(i)) * c1(:)

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
        else if (norm.eq.3.or.norm.eq.4) then
          jprod(i) = czero
          do j = 1, ny
            jprod(i) = jprod(i) + opi(j) * amp(i,k,n) * ampt(i) *( &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                        ut(i,j,3)*wl(j)) *ub(j,i) + &
                        ut(i,j,1)*pl(j) + ut(i,j,4)*ul(j) - &
                        rei*iota*(alpha(i,k,n)-alphat(i)) * &
                        (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ut(i,j,3)*wl(j)))
          end do
          write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
                   & '          |J|','          dalpha')")
          write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(Jprod(i)), zero      
        end if

        do iter = 1, niter ! loop on alpha
        
!.... compute the streamwise derivative of alpha

        dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )

!.... Evaluate the matrices (for Midpoint and backward Euler)

        alphal = (one-fact)*alpha(i-1,k,n) + fact*alpha(i,k,n)
        ubl    = (one-fact)*ub(:,i-1)      + fact*ub(:,i)
        vbl    = (one-fact)*vb(:,i-1)      + fact*vb(:,i)
        wbl    = (one-fact)*wb(:,i-1)      + fact*wb(:,i)
        ubxl   = (one-fact)*ubx(:,i-1)     + fact*ubx(:,i)
        vbxl   = (one-fact)*vbx(:,i-1)     + fact*vbx(:,i)
        wbxl   = (one-fact)*wbx(:,i-1)     + fact*wbx(:,i)
        ubyl   = (one-fact)*uby(:,i-1)     + fact*uby(:,i)
        vbyl   = (one-fact)*vby(:,i-1)     + fact*vby(:,i)
        wbyl   = (one-fact)*wby(:,i-1)     + fact*wby(:,i)

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
        if ( ipfix.eq.1 ) then
          Ab(:,:,4) = zero
        else if ( ipfix.eq.2 ) then
          Ab(:,:,4) = zero
!         Ab(:,:,4) = -Ab(:,:,4) * (x(i)-x(i-1)) / (x(i+1)-x(i))
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

        r = czero

        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            do mdof = 1, 4
              r(l0+ldof) = r(l0+ldof) + &
                 (fact-one) * ( B(l,ldof,mdof) *  uy(l,mdof)          - &
                                E(l,ldof,mdof) *  uyy(l,mdof)         + &
                               Gb(l,ldof,mdof) *  u(i-1,l,k,n,mdof) ) + &
                 Ab(l,ldof,mdof) * u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
            end do
          end do
        end do

        if (ipfix.eq.2) then
          ldof = 1
          do l = 1, ny
            l0 = (l-1)*ndof
            r(l0+ldof) = r(l0+ldof) - c1(l) * &
                         (uh(i,l,4) - uh(i-1,l,4)) / &
                         (x(i)-x(i-1))
          end do
        end if

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
!       call ZGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
!       if (info.ne.0) call error('solver$','Singular matrix$')
!       call ZGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
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

        ul = pt5 * ( u(i-1,:,k,n,1) + u(i,:,k,n,1) )
        vl = pt5 * ( u(i-1,:,k,n,2) + u(i,:,k,n,2) )
        wl = pt5 * ( u(i-1,:,k,n,3) + u(i,:,k,n,3) )
        pl = pt5 * ( u(i-1,:,k,n,4) + u(i,:,k,n,4) )

!       ul = u(i,:,k,n,1)
!       vl = u(i,:,k,n,2)
!       wl = u(i,:,k,n,3)
!       pl = u(i,:,k,n,4)

!.... compute some streamwise derivatives

        ux = (u(i,:,k,n,1) - u(i-1,:,k,n,1))/(x(i)-x(i-1))
        vx = (u(i,:,k,n,2) - u(i-1,:,k,n,2))/(x(i)-x(i-1))
        wx = (u(i,:,k,n,3) - u(i-1,:,k,n,3))/(x(i)-x(i-1))
        px = (u(i,:,k,n,4) - u(i-1,:,k,n,4))/(x(i)-x(i-1))

!==============================================================================

!.... Build the RHS for du/dalpha

        if (newton) then

        r = czero

        Ab = fact * ( fact*iota*A + fact *two*alphal*Ex - &
                    iota*Ex/(x(i)-x(i-1)) ) - &
             two*iota*Ex/(x(i)-x(i-1))

        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            do mdof = 1, ndof
              r(l0+ldof) = r(l0+ldof) - Ab(l,ldof,mdof) * u(i,l,k,n,mdof)
            end do
          end do
        end do

        Ab = ( fact-one ) * &
             ( fact*iota*A + fact*two*alphal*Ex - iota*Ex/(x(i)-x(i-1)) ) - &
               two*iota*Ex/(x(i)-x(i-1))

        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            do mdof = 1, ndof
              r(l0+ldof) = r(l0+ldof) + Ab(l,ldof,mdof) * u(i-1,l,k,n,mdof) 
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

!.... Solve the system (Note that the LU factorization is already done above)

#ifdef CRAY
        call CGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
#else
!       call ZGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
        call ZGETRS( 'N', ny*ndof, 1, A0, ny*ndof, ipvt, r, ny*ndof, info)
#endif

!.... Update the solution

        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            ua(i,l,ldof) = r(l0+ldof)
          end do
        end do

!.... Note that this is d/dalpha of the discrete dudx

        uax = (ua(i,:,1))/(x(i)-x(i-1))
        vax = (ua(i,:,2))/(x(i)-x(i-1))
        wax = (ua(i,:,3))/(x(i)-x(i-1))

        else

        ua(i,:,:) = zero; uax = zero; vax = zero; wax = zero

        end if

!==============================================================================

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

!.... Update alpha

        if (norm.eq.0) then
          f = kex
          if (.not. newton) dfda = ke(i)/iota
        else if (norm.eq.1) then
          call findumax( ny, y, ul, umax(i), nint, yint, yumax )
          dumax = (umax(i)-umax(i-1))/(x(i)-x(i-1))
!
!.... Newton cannot work for u_max normalization
!
!         dfda = cmplx( getval( ny, y,  real(ua(i,:,1)), yumax ), &
!                       getval( ny, y, aimag(ua(i,:,1)), yumax ) ) / &
!                       (x(i)-x(i-1))
!         call findumax( ny, y, ua(i,:,1), dfda, nint, yint, yumax ) 
!         dfda = dfda / (x(i)-x(i-1))
          f = dumax
          dfda = umax(i) / iota
        else if (norm.eq.2) then
          ke(i) = zero
          kex   = zero
          dfda  = zero
          do j = 1, ny
            ke(i) = ke(i) + opi(j)*( abs(ul(j))**2  )
            kex = kex + opi(j)*( conjg(ul(j))*ux(j) )
            dfda = dfda + opi(j)*(conjg(ua(i,j,1))*ux(j) + conjg(ul(j))*uax(j))
          end do
          f = kex
          if (.not. newton) dfda = ke(i)/iota
        else if (norm.eq.3) then
          Jprod(i) = czero
          dfda = czero
          do j = 1, ny
            Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * ampt(i) * ( &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                        ut(i,j,3)*wl(j)) * ub(j,i) + &
                       ut(i,j,1)*pl(j) + ut(i,j,4)*ul(j) - &
                       rei*iota*(alpha(i,k,n)-alphat(i)) * &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ut(i,j,3)*wl(j)))
            dfda = dfda + opi(j) * amp(i,k,n) * ampt(i) * ( &
                       (ut(i,j,1)*ua(i,j,1)+ut(i,j,2)*ua(i,j,2)+ &
                        ut(i,j,3)*ua(i,j,3)) * ub(j,i) + &
                       ut(i,j,1)*ua(i,j,4) + ut(i,j,4)*ua(i,j,1) - &
                       rei*iota*(ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                                 ut(i,j,3)*wl(j)) - &
                       rei*iota*(alpha(i,k,n)-alphat(i)) * &
                       (ut(i,j,1)*ua(i,j,1)+ut(i,j,2)*ua(i,j,2)+&
                        ut(i,j,3)*ua(i,j,3)) )
          end do
          dfda = dfda ! + pt5*iota*Jprod(i)/(x(i)-x(i-1))
          Jprodx = (Jprod(i)-Jprod(i-1))/(x(i)-x(i-1))
          dfda = dfda / (x(i)-x(i-1))
          f = Jprodx
          if (.not. newton) dfda = -iota * Jprod(i)
        else if (norm.eq.4) then
          Jprod(i) = czero
          dfda = czero
          do j = 1, ny
            Jprod(i) = Jprod(i) + opi(j) * amp(i,k,n) * ampt(i) *( &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                        ut(i,j,3)*wl(j)) * ub(j,i) + &
                       ut(i,j,1)*pl(j) + ut(i,j,4)*ul(j) - &
                       rei*iota*(alpha(i,k,n)-alphat(i)) * &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ut(i,j,3)*wl(j)))
            dfda = dfda + opi(j) * amp(i,k,n) * ampt(i) * ( &
                       (ut(i,j,1)*ua(i,j,1)+ut(i,j,2)*ua(i,j,2)+ &
                        ut(i,j,3)*ua(i,j,3)) * ub(j,i) + &
                       ut(i,j,1)*ua(i,j,4) + ut(i,j,4)*ua(i,j,1) - &
                       rei*iota*(ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                                 ut(i,j,3)*wl(j)) - &
                       rei*iota*(alpha(i,k,n)-alphat(i)) * &
                       (ut(i,j,1)*ua(i,j,1)+ut(i,j,2)*ua(i,j,2)+ &
                        ut(i,j,3)*ua(i,j,3)) ) + &
                       opi(j) * iota * pt5 * (x(i)-x(i-1)) * &
                       amp(i,k,n) * ampt(i) * ( &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ &
                        ut(i,j,3)*wl(j)) * ub(j,i) + &
                       ut(i,j,1)*pl(j) + ut(i,j,4)*ul(j) - &
                       rei*iota*(alpha(i,k,n)-alphat(i)) * &
                       (ut(i,j,1)*ul(j)+ut(i,j,2)*vl(j)+ut(i,j,3)*wl(j)) )
          end do
          Jprodx = (Jprod(i)-Jprod(i-1))/(x(i)-x(i-1))
          dfda   = dfda / (x(i)-x(i-1))
          f      = Jprodx
        end if

        alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                            sor*(alpha(i,k,n) - f / dfda )

!.... update the mode amplitude

        amp(i,k,n) = czero
        do il = is, i-1
          amp(i,k,n) = amp(i,k,n) + &
                       pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
                       (x(il+1)-x(il))
        end do
        amp(i,k,n) = exp( iota * amp(i,k,n) )

!.... convergence check

        if (relative) then
          err = abs(f/(dfda*alpha(i,k,n)))
        else
          err = abs(f/dfda)
        end if

        if (norm.eq.0) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), ke(i), err
        else if (norm.eq.1) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(umax(i)), err
        else if (norm.eq.2) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
               abs(amp(i,k,n)), ke(i), err
        else if (norm.eq.3.or.norm.eq.4) then
          write(*,"('->',i2,1x,6(1pe13.6,1x))") iter, alpha(i,k,n), &
                abs(amp(i,k,n)), abs(Jprod(i)), err
        end if

        if (err .lt. tol) exit

        end do    ! loop on iter

!=============================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=============================================================================
        write(9,"(9(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), &
                 abs(kex/ke(i)), amp(i,k,n), ke(i)*abs(amp(i,k,n))**2, &
                 abs(Jprod(i))
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

!.... save a non interpolated field file

      open  (10, file='field.pse', form='unformatted')
      write (10) ie-is+1, ny
      write (10) (((u(i,j,k,n,idof), j= 1,ny), i= is,ie), idof= 1,4)
      write (10) (x(i), i= is,ie)
      write (10) (y(j), j= 1,ny)
      write (10) (amp(i,k,n), i = is, ie)
      write (10) (((ua(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
      write (10) (alpha(i,k,n), i= is,ie)
      close (10)

!.... save a field file (possibly interpolated)

      call post(1)

!.... compute growth rate based on various other quantitities

      call growth(ie-is+1, ny, x(is:ie), y, amp(is:ie,k,n), &
                  alpha(is:ie,k,n), u(is:ie,:,k,n,:), opi, nint, yint)


      return
100   call error('solver$','Error reading pse.in$')
      end subroutine solver
