!==============================================================================
      subroutine mfd(  )
!
!     Compute a starting profile for meanflow distortion
!
!     S. Scott Collis
!
!==============================================================================
      use global
      implicit none

      complex :: ul(ny),  vl(ny),  wl(ny),  pl(ny)      

!============================================================================
!     Initialize starting profile for Mean-Flow-Distortion (MFD)
!============================================================================

      c1(:) = one/(one + cur(i)*y(:))
      c2(:) = cur(i) * c1(:)

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
        do j = 1, ny
          ke = ke + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + abs(wl(j))**2)
        end do
 
!.... localize the nonlinear term

        Fl = F(n,k,:,:)

!.... compute a norm of the (n,k) nonlinear term

        Fnorm = zero
        do j = 1, ny
          Fnorm = Fnorm + opi(j) * ( abs(Fl(j,1))**2 + abs(Fl(j,2))**2 + &
                                     abs(Fl(j,3))**2 )
          if (n.eq.1 .and. k.eq.1) then
            write (7,"(10(1pe13.6,1x))") y(j), real(Fl(j,1)), aimag(Fl(j,1)), &
                                         real(Fl(j,2)), aimag(Fl(j,2)), &
                                         real(Fl(j,3)), aimag(Fl(j,3)), &
                                         real(Fl(j,4)), aimag(Fl(j,4))
          end if
        end do
        
        if ( (Fnorm/unorm .le. eps) .and. (ke/tke .le. eps) ) goto 80

        if (ke .eq. zero) then
!         alpha(i,k,n) = cmplx(real(kt(n)*alpha(i,1,2)),zero)
!         alpha(i,k,n) = kt(n)*alpha(i,1,2)
          if (n.eq.1) then
!           alpha(i,k,n) = cmplx(0.0,2.0*aimag(alpha(i,1,2)))
            write(*,*) 'enter alp_i =='
            read(*,*) alpi
            alpha(i,k,n) = cmplx(0.0,alpi)
          else
            alpha(i,k,n) = kt(n)*alpha(i,1,2)
          end if
          new(k,n) = .true.
        end if

        if ((.not. new(k,n)) .or. (.not. (k.eq.1 .and. n.eq.1)) ) goto 80

        dalpha(k,n) = czero

!       write(*,"('=> ',2(i2,1x),6(1pe13.6,1x))") kz(k), kt(n), alpha(i,k,n), &
!                                                 amp(i,k,n)

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

!       if (k.eq.1 .and. n.eq.1) then  ! always use for 00 mode
!         A0(l0+1,:) = czero
!         do m = 1, ny
!           m0 = (m-1)*ndof
!           A0(l0+1,m0+1) = opy(l,m)
!         end do
!       else
!         A0(l0+1,:) = czero
!         A0(l0+1,l0+1) = cone
!       endif

        if (k.eq.1 .and. n.eq.1) then  ! always use for 00 mode
!         A0(l0+2,:) = czero
!         do m = 1, ny
!           m0 = (m-1)*ndof
!           A0(l0+2,m0+2) = opy(l,m)
!         end do
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

        do ldof = 1, ndof
          do l = 1, ny
            l0 = (l-1)*ndof
            r(l0+ldof) = r(l0+ldof) + F(n,k,l,ldof) / amp(i,k,n)
          end do
        end do

!.... Enforce the boundary conditions on the RHS

        l = 1
        l0 = (l-1)*ndof
        r(l0+1) = czero
        r(l0+2) = czero
        r(l0+3) = czero
        if (ipwbc.eq.2) then
          r(l0+4) = czero
        else
          r(l0+4) = czero
        end if
        
        l = ny
        l0 = (l-1)*ndof
!       r(l0+1) = czero
!       r(l0+2) = czero
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

!.... update local values

        ul = u(i,:,k,n,1)
        vl = u(i,:,k,n,2)
        wl = u(i,:,k,n,3)
        pl = u(i,:,k,n,4)

        ke = zero
        do j = 1, ny
          ke = ke + opi(j)*( abs(ul(j))**2 + abs(vl(j))**2 + abs(wl(j))**2 )
        end do

        write(*,"(3(i2,1x),10(1pe13.6,1x))") kz(k), kt(n), 1, &
          alpha(i,k,n), abs(amp(i,k,n)), ke, zero

80      continue

        end do    ! loop on k

        end do    ! loop on n

!.... output the initial profile(s)

        i = is; k = 1
        do n = 1, mt
          ul = u(i,:,k,n,1)
          vl = u(i,:,k,n,2)
          wl = u(i,:,k,n,3)
          pl = u(i,:,k,n,4)
!         if (ipwbc .eq. 1) pl(1) = pl(2) - (pl(3)-pl(2))/(y(3)-y(2))*y(2)
          write(79+n,"('# ',4(1pe20.13,1x))") real(amp(i,k,n)),&
              aimag(amp(i,k,n)), real(alpha(i,k,n)),aimag(alpha(i,k,n))
          do j = 1, ny
            write(79+n,"(10(1pe20.13,1x))") y(j), &
              real(ul(j)), aimag(ul(j)), &
              real(vl(j)), aimag(vl(j)), &
              real(wl(j)), aimag(wl(j)), &
              real(pl(j)), aimag(pl(j))
          end do
        end do

        return
        end
