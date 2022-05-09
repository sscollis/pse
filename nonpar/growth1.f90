!==============================================================================
      subroutine growth( nx, ny, x, y, amp, alpha, u, opi, nint, yint )
!==============================================================================
      use fmax
      implicit none
      integer, parameter :: ndof=4, japt=20
      real, parameter    :: zero=0.0, pt5=0.5, one=1.0

      integer :: nx, ny, nint
      real    :: x(nx), y(ny), yint(nint,ny), opi(ny)
      complex :: amp(nx), alpha(nx), u(nx,ny,ndof), alphat

      integer :: i, j, ii, idof, jmax
      real    :: f(ny), fr(ny), fi(ny), frmax, fimax, ymax
      real    :: q(9,nx), v(nint,ndof), vr(nint,ndof), vi(nint,ndof)
      complex :: umax(nx,ndof)

      real    :: dxb, dxf, ds, fc1, cc1, bc1

      real, external :: ddot
      integer, external :: isamax
      logical :: pse = .true.          ! method flag
!==============================================================================

      write(*,"('Writing grow.dat...',/)")

!.... make the full u-field (This is a little less accurate, but same as HLNS)

      if (.not. pse) then
	do i = 1, nx
	  u(i,:,:) = amp(i) * u(i,:,:)
	end do
      end if

!.... initialize

      q = zero

!.... loop over the streamwise direction

      do i = 1, nx

	do idof = 1, ndof

	  do j = 1, ny
	    f(j)  =   abs( u(i,j,idof) )
	    fr(j) =  real( u(i,j,idof) )
	    fi(j) = aimag( u(i,j,idof) )
	  end do

!.... interpolate to a fine mesh before finding the maximum

	  do j = 1, nint
	    v(j,idof)  = zero
	    vr(j,idof) = zero
	    vi(j,idof) = zero
	    v(j,idof)  = ddot( ny, yint(j,1), nint, f, 1 )
	    vr(j,idof) = ddot( ny, yint(j,1), nint, fr, 1 )
	    vi(j,idof) = ddot( ny, yint(j,1), nint, fi, 1 )
	  end do
	  
!.... find the maximum
	  
	  jmax = isamax(nint,v(1,idof),1)
	  q(idof,i) = v(jmax,idof)

!	  if (idof.eq.2) write(*,*) x(i), jmax, q(idof,i)

	  umax(i,idof) = cmplx( vr(jmax,idof), vi(jmax,idof) )

!.... Try the BSpline max method for streamwise velocity
!
!.... This has problems sometimes when the root is properly bracketed
!
!	  if (idof.eq.1 .or. idof.eq.2) then
!	    q(idof,i)    = findmax( ny, y, f, ymax)
!	    frmax        = getval( ny, y, fr, ymax)
!	    fimax        = getval( ny, y, fi, ymax)
!	    umax(i,idof) = cmplx( frmax, fimax )
!	  end if

!.... integrate in y to compute the total energy
	  
	  if (idof .le. 3) then
	    do j = 1, ny
	      q(5,i) = q(5,i) + opi(j) * f(j)**2
	    end do
	  end if

	end do

      end do

!.... now compute the growth rates on the interior

      do i = 2, nx-1

!.... form the derivative opperators

	dxb = x(i) - x(i-1)
	dxf = x(i+1) - x(i)
	ds = 1./(dxf+dxb)
	fc1 = ds * dxb / dxf
	cc1 = ds * (dxf/dxb - dxb/dxf)
	bc1 = -ds * dxf / dxb
	
	q(6,i) = (fc1*q(1,i+1) + cc1*q(1,i) + bc1*q(1,i-1)) / &
                  max(q(1,i),1.e-10)
	q(7,i) = .5*(fc1*q(5,i+1) + cc1*q(5,i) + bc1*q(5,i-1)) / &
                 max(q(5,i),1.e-10)
	alphat = (fc1*u(i+1,japt,1) + cc1*u(i,japt,1) + &
	         bc1*u(i-1,japt,1)) / u(i,japt,1)
	q(8,i) = -aimag(alphat)
	q(9,i) =   real(alphat)

      end do

!.... use low-order differences on the boundaries

      q(6,1) = (q(1,2)-q(1,1)) / (x(2)-x(1)) / max(q(1,1),1.e-10)
      q(7,1) = .5*(q(5,2)-q(5,1)) / (x(2)-x(1)) / max(q(5,1),1.e-10)
      alphat = (u(2,japt,1) - u(1,japt,1)) / (x(2)-x(1))
      q(8,1) = -aimag (alphat)
      q(9,1) = -real (alphat)

      q(6,nx) = (q(1,nx)-q(1,nx-1)) / (x(nx)-x(nx-1)) / max(q(1,nx),1.e-10)
!     q(6,nx) = (q(1,nx)-q(1,nx-1)) / (x(nx)-x(nx-1)) / &
!                max(pt5*(q(1,nx)+q(1,nx-1)),1.e-10)
      q(7,nx) = .5*(q(5,nx)-q(5,nx-1)) / (x(nx)-x(nx-1)) / &
                 max(q(5,nx),1.e-10)
!     q(7,nx) = .5*(q(5,nx)-q(5,nx-1)) / (x(nx)-x(nx-1)) / &
!                max(pt5*(q(5,nx)+q(5,nx-1)),1.e-10)
      alphat = (u(nx,japt,1) - u(nx-1,japt,1)) / (x(nx)-x(nx-1))
      q(8,nx) = -aimag (alphat)
      q(9,nx) = -real (alphat)

!.... correct for PSE wave term (This is the more accurate way)

      if (pse) then
	q(1,:) = q(1,:) * abs(amp)         ! umax
	q(2,:) = q(2,:) * abs(amp)         ! vmax
	q(3,:) = q(3,:) * abs(amp)         ! wmax
	q(4,:) = q(4,:) * abs(amp)         ! pmax
	q(5,:) = q(5,:) * abs(amp)**2      ! KE
	q(6,:) = q(6,:) - aimag(alpha(:))  ! growth rate (umax)
	q(7,:) = q(7,:) - aimag(alpha(:))  ! growth rate (KE)
      end if

!.... write out the results

      open (2, file="grow.dat")
      do i= 1,nx
	write (2,"(11(1pe20.13,1x))") x(i), (q(idof,i), idof= 1,9)
      end do
      close (2)

!.... write out the wall pressure

      open (2, file="pwall.dat")
      do i= 1,nx
        write (2,"(11(1pe20.13,1x))") x(i), u(i,1,4)*amp(i), &
                                      abs(u(i,1,4)*amp(i))
      end do
      close (2)

!.... if using the discrete adjoint, you must sample the pressure at the 
!.... staggered mesh location between nodes 2 and 3 (node 1 is bogus)

!!$    open (2, file="pwall.dat")
!!$    do i= 1,nx
!!$    write (2,"(11(1pe20.13,1x))") x(i), pt5*(u(i,1,4)+u(i,2,4))*amp(i), &
!!$                  abs(1./4.*u(i,1,4)+3./4.*u(i,3,4))*abs(amp(i))
!!$    end do
!!$    close (2)

!.... revert full u-field back to shape function

      if (.not. pse) then
	do i = 1, nx
	  u(i,:,:) = u(i,:,:) / amp(i)
	end do
      end if

      return
      end subroutine growth

!==============================================================================
      subroutine dgrowth( nx, ny, x, y, amp, alpha, u, opi, nint, yint )
!==============================================================================
!
!     Special interpolation for discrete adjoint
!
!==============================================================================
      use fmax
      implicit none
      integer, parameter :: ndof=4, japt=20
      real, parameter    :: zero=0.0, pt5=0.5, one=1.0

      integer :: nx, ny, nint
      real    :: x(nx), y(ny), yint(nint,ny), opi(ny)
      complex :: amp(nx), alpha(nx), u(nx,ny,ndof), alphat

      integer :: i, j, ii, idof, jmax
      real    :: f(ny), fr(ny), fi(ny), frmax, fimax, ymax
      real    :: q(9,nx), v(nint,ndof), vr(nint,ndof), vi(nint,ndof)
      complex :: umax(nx,ndof)

      real    :: dxb, dxf, ds, fc1, cc1, bc1

      real, external :: ddot
      integer, external :: isamax
      logical :: pse = .true.          ! method flag
!==============================================================================

      write(*,"('Writing grow.dat...',/)")

!.... make the full u-field (This is a little less accurate, but same as HLNS)

      if (.not. pse) then
	do i = 1, nx
	  u(i,:,:) = amp(i) * u(i,:,:)
	end do
      end if

!.... initialize

      q = zero

!.... loop over the streamwise direction

      do i = 1, nx

	do idof = 1, ndof

	  do j = 1, ny
	    f(j)  =   abs( u(i,j,idof) )
	    fr(j) =  real( u(i,j,idof) )
	    fi(j) = aimag( u(i,j,idof) )
	  end do

!.... interpolate to a fine mesh before finding the maximum

	  do j = 1, nint
	    v(j,idof)  = zero
	    vr(j,idof) = zero
	    vi(j,idof) = zero
	    v(j,idof)  = ddot( ny, yint(j,1), nint, f, 1 )
	    vr(j,idof) = ddot( ny, yint(j,1), nint, fr, 1 )
	    vi(j,idof) = ddot( ny, yint(j,1), nint, fi, 1 )
	  end do
	  
!.... find the maximum
	  
	  jmax = isamax(nint,v(1,idof),1)
	  q(idof,i) = v(jmax,idof)

!	  if (idof.eq.2) write(*,*) x(i), jmax, q(idof,i)

	  umax(i,idof) = cmplx( vr(jmax,idof), vi(jmax,idof) )

!.... Try the BSpline max method for streamwise velocity
!
!.... This has problems sometimes when the root is properly bracketed
!
!	  if (idof.eq.1 .or. idof.eq.2) then
!	    q(idof,i)    = findmax( ny, y, f, ymax)
!	    frmax        = getval( ny, y, fr, ymax)
!	    fimax        = getval( ny, y, fi, ymax)
!	    umax(i,idof) = cmplx( frmax, fimax )
!	  end if

!.... integrate in y to compute the total energy
	  
	  if (idof .le. 3) then
	    do j = 1, ny
	      q(5,i) = q(5,i) + opi(j) * f(j)**2
	    end do
	  end if

	end do

      end do

!.... now compute the growth rates on the interior

      do i = 2, nx-1

!.... form the derivative opperators

	dxb = x(i) - x(i-1)
	dxf = x(i+1) - x(i)
	ds = 1./(dxf+dxb)
	fc1 = ds * dxb / dxf
	cc1 = ds * (dxf/dxb - dxb/dxf)
	bc1 = -ds * dxf / dxb
	
	q(6,i) = (fc1*q(1,i+1) + cc1*q(1,i) + bc1*q(1,i-1)) / &
                  max(q(1,i),1.e-10)
	q(7,i) = .5*(fc1*q(5,i+1) + cc1*q(5,i) + bc1*q(5,i-1)) / &
                 max(q(5,i),1.e-10)
	alphat = (fc1*u(i+1,japt,1) + cc1*u(i,japt,1) + &
	         bc1*u(i-1,japt,1)) / u(i,japt,1)
	q(8,i) = -aimag(alphat)
	q(9,i) =   real(alphat)

      end do

!.... use low-order differences on the boundaries

      q(6,1) = (q(1,2)-q(1,1)) / (x(2)-x(1)) / max(q(1,1),1.e-10)
      q(7,1) = .5*(q(5,2)-q(5,1)) / (x(2)-x(1)) / max(q(5,1),1.e-10)
      alphat = (u(2,japt,1) - u(1,japt,1)) / (x(2)-x(1))
      q(8,1) = -aimag (alphat)
      q(9,1) = -real (alphat)

      q(6,nx) = (q(1,nx)-q(1,nx-1)) / (x(nx)-x(nx-1)) / max(q(1,nx),1.e-10)
!     q(6,nx) = (q(1,nx)-q(1,nx-1)) / (x(nx)-x(nx-1)) / &
!                max(pt5*(q(1,nx)+q(1,nx-1)),1.e-10)
      q(7,nx) = .5*(q(5,nx)-q(5,nx-1)) / (x(nx)-x(nx-1)) / &
                 max(q(5,nx),1.e-10)
!     q(7,nx) = .5*(q(5,nx)-q(5,nx-1)) / (x(nx)-x(nx-1)) / &
!                max(pt5*(q(5,nx)+q(5,nx-1)),1.e-10)
      alphat = (u(nx,japt,1) - u(nx-1,japt,1)) / (x(nx)-x(nx-1))
      q(8,nx) = -aimag (alphat)
      q(9,nx) = -real (alphat)

!.... correct for PSE wave term (This is the more accurate way)

      if (pse) then
	q(1,:) = q(1,:) * abs(amp)         ! umax
	q(2,:) = q(2,:) * abs(amp)         ! vmax
	q(3,:) = q(3,:) * abs(amp)         ! wmax
	q(4,:) = q(4,:) * abs(amp)         ! pmax
	q(5,:) = q(5,:) * abs(amp)**2      ! KE
	q(6,:) = q(6,:) - aimag(alpha(:))  ! growth rate (umax)
	q(7,:) = q(7,:) - aimag(alpha(:))  ! growth rate (KE)
      end if

!.... write out the results

      open (2, file="grow.dat")
      do i= 1,nx
	write (2,"(11(1pe20.13,1x))") x(i), (q(idof,i), idof= 1,9)
      end do
      close (2)

!.... write out the wall pressure

!!$      open (2, file="pwall.dat")
!!$      do i= 1,nx
!!$        write (2,"(11(1pe20.13,1x))") x(i), u(i,1,4)*amp(i), &
!!$                                      abs(u(i,1,4)*amp(i))
!!$      end do
!!$      close (2)

!.... if using the discrete adjoint, you must sample the pressure at the 
!.... staggered mesh location between nodes 2 and 3 (node 1 is bogus)

      open (2, file="pwall.dat")
      do i= 1,nx
        write (2,"(11(1pe20.13,1x))") x(i), &
             (1./4.*u(i,1,4)+3./4.*u(i,3,4))*amp(i), &
             abs(1./4.*u(i,1,4)+3./4.*u(i,3,4))*abs(amp(i))
      end do
      close (2)

!.... revert full u-field back to shape function

      if (.not. pse) then
	do i = 1, nx
	  u(i,:,:) = u(i,:,:) / amp(i)
	end do
      end if

      return
      end subroutine dgrowth

!==============================================================================
      subroutine findumax( ny, y, u, umax, nint, yint, ymax )
!==============================================================================
      use fmax
      implicit none
      real, parameter :: zero=0.0, pt5=0.5, one=1.0

      integer :: ny, nint
      real    :: y(ny), yint(nint,ny)
      complex :: u(ny)

      integer :: i, j, ii, jmax
      real    :: f(ny), fr(ny), fi(ny), v(nint), vr(nint), vi(nint)
      complex :: umax

      real, external :: ddot
      integer, external :: isamax

      real :: famax, frmax, fimax, ymax
!==============================================================================

      do j = 1, ny
	f(j)  =   abs( u(j) )
	fr(j) =  real( u(j) )
	fi(j) = aimag( u(j) )
      end do

!.... slow brute force technique

      if (.false.) then

!.... interpolate to a fine mesh before finding the maximum
	
	do j = 1, nint
	  v(j)  = zero
	  vr(j) = zero
	  vi(j) = zero
	  v(j)  = ddot( ny, yint(j,1), nint, f, 1 )
	  vr(j) = ddot( ny, yint(j,1), nint, fr, 1 )
	  vi(j) = ddot( ny, yint(j,1), nint, fi, 1 )
	end do
	
!.... find the maximum
	
	jmax = isamax(nint,v,1)
	umax = cmplx( vr(jmax), vi(jmax) )
	
      else

!.... Fast, accurate B-spline method
	
	famax = findmax( ny, y, f, ymax)
	frmax = getval( ny, y, fr, ymax)
	fimax = getval( ny, y, fi, ymax)
	umax  = cmplx( frmax, fimax )

      end if
      
      return
      end subroutine findumax
  
!==============================================================================
      integer function isamax(n,sx,incx)
!==============================================================================

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.

      real sx(1),smax
      integer i,incx,ix,n

      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20

!.... code for increment not equal to 1

      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      smax = abs(sx(ix))
      ix = ix + incx
      do i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
       end do
      return

!.... code for increment equal to 1

   20 smax = abs(sx(1))
      do i = 2,n
         if(abs(sx(i)).le.smax) exit
         isamax = i
         smax = abs(sx(i))
      end do
       
      return
      end function isamax
