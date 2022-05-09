!==============================================================================
      subroutine nonlin( i, ur, ui, Fr, Fi, c1, c2, unorm, tke )
!
!     Compute the nonlinear term for NPSE
!
!     S. Scott Collis
!
!     9-7-98:    Changed to Chebyschev integration
!==============================================================================
      use global
      implicit none

      integer :: i 
      complex :: ui(mt,nz,12)
      real    :: ur(2*mt,nz,12), c1(ny), c2(ny)
      real    :: Fr(2*mt,nz,ny,ndof), unorm, tke
      complex :: Fi(mt,nz,ny,ndof)

      integer :: j, k, n, idof, l, il
      complex :: uy(ny,mt,nz), vy(ny,mt,nz), wy(ny,mt,nz)
      real    :: fact, ke
      complex :: ul(ny), vl(ny), wl(ny)
!==============================================================================

      write(*,*) "WARNING:  FFT routines are turned off in nonlin.f90"

!.... compute the amplitude for each mode

      do n = 1, mt
        do k = 1, nz
          amp(i,k,n) = czero
          do il = 1, i-1
            amp(i,k,n) = amp(i,k,n) + &
                         pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
                         (x(il+1)-x(il))
          end do
!         write (65,"(2(i3,1x),6(1pe20.13,1x))") kz(k),kt(n),x(i),amp(i,k,n)
          amp(i,k,n) = exp( iota * amp(i,k,n) )
!         write (66,"(2(i3,1x),6(1pe20.13,1x))") kz(k),kt(n),x(i),amp(i,k,n)
        end do
      end do

!.... compute the total energy norm used to measure the nonlinear amplitude

      tke = zero
      unorm = zero
      do n = 1, mt
        fact = two
        if (n.eq.1) fact = one
        do k = 1, nz
          ke = zero
          ul = u(i,:,k,n,1)
          vl = u(i,:,k,n,2)
          wl = u(i,:,k,n,3)
          do j = 1, ny
            ke = ke + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + abs(wl(j))**2)
          end do
          tke = tke + fact * ke
          unorm = unorm + fact * conjg(amp(i,k,n))*amp(i,k,n) * &
                                 conjg(alpha(i,k,n))*alpha(i,k,n) * ke
!         write(*,*) n, k, tke, unorm
        end do
      end do

!.... form the y-derivatives (should use BLAS here)

      do n = 1, mt
        do k = 1, nz
          do j = 1, ny
            uy(j,n,k) = czero
            vy(j,n,k) = czero
            wy(j,n,k) = czero
            do l = 1, ny
              uy(j,n,k) = uy(j,n,k) + opy(j,l) * u(i,l,k,n,1)
              vy(j,n,k) = vy(j,n,k) + opy(j,l) * u(i,l,k,n,2)
              wy(j,n,k) = wy(j,n,k) + opy(j,l) * u(i,l,k,n,3)
            end do
          end do
        end do
      end do
 
!.... form the nonlinear terms

      do j = 1, ny

        if (i.gt.1) then

        do n = 1, mt
          do k = 1, nz
            ui(n,k, 1) = amp(i,k,n) * u(i,j,k,n,1)
            ui(n,k, 2) = amp(i,k,n) * u(i,j,k,n,2)
            ui(n,k, 3) = amp(i,k,n) * u(i,j,k,n,3)
            ui(n,k, 4) = amp(i,k,n) * ( (u(i,j,k,n,1) - u(i-1,j,k,n,1) ) / &
                         (x(i) - x(i-1)) + iota * alpha(i,k,n) * u(i,j,k,n,1) )
            ui(n,k, 5) = amp(i,k,n) * ( (u(i,j,k,n,2) - u(i-1,j,k,n,2) ) / &
                         (x(i) - x(i-1)) + iota * alpha(i,k,n) * u(i,j,k,n,2) )
            ui(n,k, 6) = amp(i,k,n) * ( (u(i,j,k,n,3) - u(i-1,j,k,n,3) ) / &
                         (x(i) - x(i-1)) + iota * alpha(i,k,n) * u(i,j,k,n,3) )
            ui(n,k, 7) = amp(i,k,n) * uy(j,n,k)
            ui(n,k, 8) = amp(i,k,n) * vy(j,n,k)
            ui(n,k, 9) = amp(i,k,n) * wy(j,n,k)
            ui(n,k,10) = amp(i,k,n) * iota * kz(k) * u(i,j,k,n,1)
            ui(n,k,11) = amp(i,k,n) * iota * kz(k) * u(i,j,k,n,2)
            ui(n,k,12) = amp(i,k,n) * iota * kz(k) * u(i,j,k,n,3)
          end do
        end do

        else

        do n = 1, mt
          do k = 1, nz
            ui(n,k, 1) = amp(i,k,n) * u(i,j,k,n,1)
            ui(n,k, 2) = amp(i,k,n) * u(i,j,k,n,2)
            ui(n,k, 3) = amp(i,k,n) * u(i,j,k,n,3)
            ui(n,k, 4) = zero
            ui(n,k, 5) = zero
            ui(n,k, 6) = zero
            ui(n,k, 7) = amp(i,k,n) * uy(j,n,k)
            ui(n,k, 8) = amp(i,k,n) * vy(j,n,k)
            ui(n,k, 9) = amp(i,k,n) * wy(j,n,k)
            ui(n,k,10) = amp(i,k,n) * iota * kz(k) * u(i,j,k,n,1)
            ui(n,k,11) = amp(i,k,n) * iota * kz(k) * u(i,j,k,n,2)
            ui(n,k,12) = amp(i,k,n) * iota * kz(k) * u(i,j,k,n,3)
          end do
        end do

        end if

!       if (j.eq.28) then
!         do n = 1, mt
!           write(80,"(i3,4(1x,1pe13.6))") n, ui(n,1,1)
!         end do
!       end if

!.... convert to physical space

        do l = 1, 12
          call fft( 1, nt, nz, ur(1,1,l) )
        end do

!       if (j.eq.28) then
!         do n = 1, nt
!           write(81,"(i3,4(1x,1pe13.6))") n, ur(n,1,1), ur(n,1,2), ur(n,1,3)
!         end do
!       end if
!
!       do l = 1, 1
!         call fft( -1, nt, nz, ur(1,1,l) )
!       end do
!
!       if (j.eq.28) then
!         do n = 1, mt
!           write(82,"(i3,4(1x,1pe13.6))") n, ui(n,1,1)
!         end do
!         stop 
!       end if

        do n = 1, nt
          do k = 1, nz
            Fr(n,k,j,1) = -( ur(n,k,1) * (c1(j) * ur(n,k,4) + &
                                          c2(j) * ur(n,k,2)) + & 
                             ur(n,k,2) * ur(n,k,7) + ur(n,k,3) * ur(n,k,10) )
            Fr(n,k,j,2) = -( ur(n,k,1) * (c1(j) * ur(n,k,5) - &
                                          c2(j) * ur(n,k,1)) + &
                             ur(n,k,2) * ur(n,k,8) + ur(n,k,3) * ur(n,k,11) )
            Fr(n,k,j,3) = -( ur(n,k,1) * c1(j) * ur(n,k,6) + &
                             ur(n,k,2) * ur(n,k,9) + ur(n,k,3) * ur(n,k,12) )
            Fr(n,k,j,4) = zero
          end do
        end do

!       if (j.eq.28) then
!         do n = 1, nt
!           write(83,"(i3,4(1x,1pe13.6))") n, Fr(n,1,j,1), Fr(n,1,j,2), &
!                                           Fr(n,1,j,3),
!         end do
!       end if

!.... convert to wave space

        do l = 1, 3
          call fft( -1, nt, nz, Fr(1,1,j,l) )
        end do

!       if (j.eq.28) then
!         do n = 1, mt
!           write(84,"(i3,6(1x,1pe13.6))") n, Fi(n,1,j,1), Fi(n,1,j,2), &
!                                           Fi(n,1,j,3),
!         end do
!         stop 
!       end if

      end do

      return
      end subroutine nonlin

#ifdef CRAY
      module modfft
        real, allocatable :: table(:), work(:)
      end module modfft
!==============================================================================
      subroutine fft( sign, nt, nz, ur )
!==============================================================================
      use modfft
      implicit none
      integer :: sign, nt, nz
      real    :: ur(2*((nt+2)/2),nz)
      integer :: k, n, mt, n0
      integer :: isys = 0
      real    :: scale = 1.0
!==============================================================================
      mt = (nt+2)/2

!.... sign = -1 is the forward transform, sign = 1 is the inverse transform
!.... sign = 0 initializes the FFT coefficients

      if (sign .eq. 1) then
        scale = 1.0
        call csfft2d(-1,nt,nz,scale,ur,mt,ur,2*mt,table,work,isys)
      else if (sign .eq. -1) then
        scale = 1.0 / ( real(nt) * real(nz) )
        call scfft2d(1,nt,nz,scale,ur,2*mt,ur,mt,table,work,isys)
      else if (sign .eq. 0) then
        allocate( table(100+2*(nt+nz)), work(512*max(nt,nz)) )
        call scfft2d(0,nt,nz,scale,ur,2*mt,ur,mt,table,work,isys)
      else
        call error('fft$','Illegal value of sign$')
      end if
      return
      end
#else
      module modfft
        real, allocatable :: coeft(:), coefz(:)
      end module modfft
!==============================================================================
      subroutine fft( sign, nt, nz, ur )
!==============================================================================
      use modfft
      implicit none
      integer :: sign, nt, nz
      real    :: ur(2*((nt+2)/2),nz)
      integer :: k, n, mt, n0
!==============================================================================
      mt = (nt+2)/2

!.... sign = -1 is the forward transform, sign = 1 is the inverse transform
!.... sign = 0 initializes the FFT coefficients

      if (sign .eq. 1) then
        do n = 1, mt
          n0 = 1 + (n-1)*2
!         call zfft1d( 1, nz, ur(n0,1), mt, coefz )
        end do
        do k = 1, nz
!         call zdfft1du( -1, nt, ur(1,k), 1, coeft )  
        end do
      else if (sign .eq. -1) then
        do k = 1, nz
!         call dzfft1du( 1, nt, ur(1,k), 1, coeft )
!         call zscal1d( mt, 1.0/real(nt), ur(1,k), 1 )
        end do
        do n = 1, mt
          n0 = 1 + (n-1)*2
!         call zfft1d( -1, nz, ur(n0,1), mt, coefz )
!         call zscal1d( nz, 1.0/real(nz), ur(n0,1), mt )
        end do
      else if (sign .eq. 0) then
        allocate( coeft(nt+15), coefz(2*(nz+15)) )
!       call dzfft1dui( nt, coeft )
!       call zfft1di( nz, coefz )
      else
        call error('fft$','Illegal value of sign$')
      end if
      return
      end
#endif
