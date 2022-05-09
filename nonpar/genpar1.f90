!==============================================================================
      subroutine genpar
!
!     Use a parallel mean flow profile
!
!==============================================================================
      use global
      implicit none

      integer :: i, j
      real :: dum

      real :: dx, xstr1, xstr2, xc, c2, aloc, bloc, alim, blim
      real, external :: fc, zeroin
      real :: s1, s2, str1, xmx
      common / xcoord / s1, s2, str1, xmx
!==============================================================================

      if (nx.le.2) call error('genpar$','Nx must be greater than 1$')

!.... read the initial condition

      open(4, file='profile.dat', status='old', err=100)
      read(4,*) ny, ymax, ystr
      write(*,"('Reading mean profile (Ny = ',i3,')')") ny
      write(*,"('Using y-grid from profile.dat...')")

!.... Note, allways uses the y-mesh stipulated in profile.dat

      allocate( x(nx), cur(nx), y(ny) )
      allocate( ub(ny,nx),  vb(ny,nx),  wb(ny,nx)  )
      allocate( uby(ny,nx), vby(ny,nx), wby(ny,nx) )
      allocate( ubx(ny,nx), vbx(ny,nx), wbx(ny,nx) )

      do j = 1, ny
	read(4,*) ub(j,1), vb(j,1), wb(j,1)
      end do
      close(4)

!.... make the parallel flow

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
