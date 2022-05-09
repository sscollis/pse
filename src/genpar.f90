!==============================================================================
      subroutine genpar
!
!     Use a parallel mean flow profile
!
!     Revised:  5-29-98  S.Collis
!               Now it interpolates the profile to a new mesh.
!
!==============================================================================
      use global
      use findmf
      implicit none

      integer :: i, j
      real :: dum

      real :: dx, xstr1, xstr2, xc, c2, aloc, bloc, alim, blim
      real, external :: fc, zeroin
      real :: s1, s2, str1, xmx
      common / xcoord / s1, s2, str1, xmx

!.... stuff for interpolation

      real, allocatable ::  uin(:), vin(:), win(:), yin(:), w(:,:), &
                            cuin(:), cvin(:), cwin(:), ct(:)

      real, external :: f
#ifdef CRAY 
      real, external :: sdot
#else
      real, external :: ddot
#endif

      character*1,parameter :: trans='n'
      integer :: ix, ii
      real :: pin, con, arg, ymaxin, ystrin, yc
!==============================================================================

      if (nx.le.2) call error('genpar$','Nx must be greater than 1$')

      allocate( x(nx), cur(nx), y(ny) )
      allocate( ub(ny,nx),  vb(ny,nx),  wb(ny,nx)  )
      allocate( uby(ny,nx), vby(ny,nx), wby(ny,nx) )
      allocate( ubx(ny,nx), vbx(ny,nx), wbx(ny,nx) )

!.... read the initial condition

      open(4, file='profile.dat', status='old', err=100)
      read (4,*,err=100,end=100) nin, ymaxin, ystrin
      write(*,"('Reading mean profile (Ny = ',i3,')')") nin
      allocate( uin(nin), vin(nin), win(nin), yin(nin), w(nin,nin), &
                cuin(nin), cvin(nin), cwin(nin), ct(nin), cyin(nin) )
      uin = zero; vin = zero; win = zero;
      do j = 1, nin
        read(4,*,err=100,end=100) dum, uin(j), vin(j), win(j) ! my format
!       read(4,*,err=100,end=100) uin(j), vin(j), win(j)      ! Streett's
      end do
      close(4)

!.... Check to see if interpolation is necessary

      if ( nin.eq.ny .and. ymaxin.eq.ymax .and. ystrin.eq.ystr) then
!     if (.true.) then
        do j = 1, ny
          ub(j,1) = uin(j)
          vb(j,1) = vin(j)
          wb(j,1) = win(j)
        end do
      else

      write(*,"('Interpolating mean flow to new mesh...')")

!.... make actual y-grid (will get redone in setup, but I need it now)

      pin = acos(-1.0) / (ny-1)
      if (ystr .eq. 0.) then  ! uniform y mesh
        do j= 1,ny
          y(j) = .5*ymax * (1. - cos((j-1) * pin))
        end do
      else                    ! stretched y mesh
        do j= 1,ny
          yc = .5 * (1. - cos((j-1) * pin))
          y(j) = ymax * ystr * yc / (1. + ystr - yc)
        end do
      end if

!.... WARNING: dirty hack

      if (.false.) then
        write(*,*) 'Hacked mesh in genpar.f90'
        do j = 1, ny
          y(j) = real(j-1) * ymax / real(ny-1)
        end do
        
        do j = 1, ny
          y(j) = ymax * ystr / ( ymax - 2.0 * ystr ) * real(j-1)/real(ny-1) &
            / ( 1.0 + ymax * ystr / ( ymax - 2.0 * ystr ) / ymax - &
            real(j-1)/real(ny-1) )
        end do
      end if

!.... make input y-grid (probably should just be input from unit 4!)

      pin = acos(-1.0) / (nin-1)
      if (ystrin .eq. 0.) then  ! uniform y mesh
        do j= 1,nin
          yin(j) = .5*ymaxin * (1. - cos((j-1) * pin))
        end do
      else                      ! stretched y mesh
        do j= 1,nin
          yc = .5 * (1. - cos((j-1) * pin))
          yin(j) = ymaxin * ystrin * yc / (1. + ystrin - yc)
        end do
      end if

!.... Interpolate to new y-grid

      ix = 1
      pin = acos(-1.0) / (nin-1)
      con = 2.0 / (nin-1)
      do j= 1,nin
        do i= 1,nin
          w(i,j) = con * cos((i-1) * (j-1) * pin)
        end do
      end do
      do i= 1,nin
        w(i,1) = .5 * w(i,1)
        w(1,i) = .5 * w(1,i)
        w(i,nin) = .5 * w(i,nin)
        w(nin,i) = .5 * w(nin,i)
      end do

#ifdef CRAY
      call sgemv ( trans, nin, nin, one, w, nin, yin, 1, zero, cyin, 1)
      call sgemv ( trans, nin, nin, one, w, nin, uin, 1, zero, cuin, 1)
      call sgemv ( trans, nin, nin, one, w, nin, vin, 1, zero, cvin, 1)
      call sgemv ( trans, nin, nin, one, w, nin, win, 1, zero, cwin, 1)
#else
      call dgemv ( trans, nin, nin, one, w, nin, yin, 1, zero, cyin, 1)
      call dgemv ( trans, nin, nin, one, w, nin, uin, 1, zero, cuin, 1)
      call dgemv ( trans, nin, nin, one, w, nin, vin, 1, zero, cvin, 1)
      call dgemv ( trans, nin, nin, one, w, nin, win, 1, zero, cwin, 1)
#endif

      do j= 1,ny
        if (y(j) .lt. yin(nin)) then
          yp = y(j)
          arg = zeroin (-1., 1., f, 0.)
          arg = acos(arg)
          do ii= 1,nin
            ct(ii) = cos((ii-1) * arg)
          end do
#ifdef CRAY
          ub(j,ix) = sdot (nin, cuin, 1, ct, 1)
          vb(j,ix) = sdot (nin, cvin, 1, ct, 1)
          wb(j,ix) = sdot (nin, cwin, 1, ct, 1)
#else
          ub(j,ix) = ddot (nin, cuin, 1, ct, 1)
          vb(j,ix) = ddot (nin, cvin, 1, ct, 1)
          wb(j,ix) = ddot (nin, cwin, 1, ct, 1)
#endif
        else
          ub(j,ix) = uin(nin)
          vb(j,ix) = vin(nin)
          wb(j,ix) = win(nin)
        end if

!       write(87,"(4(1pe13.6,1x))") y(j), ub(j,ix)
        
      end do

      deallocate( uin, vin, win, yin, w, cuin, cvin, cwin, ct, cyin )

      end if

!.... make the flow parallel

      do i= 1,nx
        do j= 1,ny
          ub(j,i) = ub(j,1)
          vb(j,i) = 0.0
          wb(j,i) = wb(j,1)
        end do
      end do

!.... make the streamwise mesh for parallel flow

!.... for the stretched mesh:  xs1 and xs2 are breakpoints, 
!....                          dx1 is the initial dx

      if (xs1 .gt. xmax) then      ! uniform mesh
        dx = xmax / (nx-1)
        do i= 1,nx
          x(i) = (i-1) * dx
        end do
      else                         ! stretched mesh
        s1 = xs1
        s2 = xs2
        xstr1 = dx1 * (nx-1)
        str1 = xstr1
        xmx = xmax

        aloc = xs1 / xstr1
        blim = xs2 / xstr1
        alim = aloc + 1.e-6
        bloc = zeroin (alim, blim, fc, 0.)
        xstr2 = (xmax - xs2) / (1. - bloc)
        c2 = .5 * (xstr2 - xstr1) / (bloc - aloc)
        do i= 1,nx
          xc = (i-1)/float(nx-1)
          if (xc .lt. aloc) then
            x(i) = xstr1 * xc
          elseif (xc .lt. bloc) then
            x(i) = xs1 + xstr1 * (xc-aloc) + c2 * (xc-aloc) **2
          else
            x(i) = xs2 + (xc - bloc) * xstr2
          endif
        end do
!       write(*,*) (x(i), i= 1,nx)
      endif

!.... assume that the curvature is zero

      write(*,"('Setting curvature to zero...')")
      cur = zero

      return
100   call error ('genpar$','Error opening profile.dat$')
      end

!==============================================================================
      real function fc(bloc)
!==============================================================================
      implicit none

      real :: aloc, bloc, xstr2, c2, xpb

      real :: xs1, xs2, xstr1, xmax
      common / xcoord / xs1, xs2, xstr1, xmax
!==============================================================================
      aloc = xs1 / xstr1
      xstr2 = (xmax - xs2) / (1. - bloc)
      c2 = .5 * (xstr2 - xstr1) / (bloc - aloc)
      xpb = xs1 + xstr1 * (bloc-aloc) + c2 * (bloc-aloc) **2
      fc = xpb - xs2
      return
      end function fc
