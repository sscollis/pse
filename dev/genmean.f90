!==============================================================================
      subroutine genmean
!
!     Use a field file for the mean flow
!
!==============================================================================
      use global
      implicit none

      integer :: i, j
      real :: dum
!==============================================================================

      open (4, file='field.mean', form='unformatted', status='old', &
	    err=120)
      read (4) nx, ny
      allocate( x(nx), y(ny), cur(nx) )
      allocate( ub(ny,nx),  vb(ny,nx),  wb(ny,nx) )
      allocate( uby(ny,nx), vby(ny,nx), wby(ny,nx) )
      allocate( ubx(ny,nx), vbx(ny,nx), wbx(ny,nx) )
      read (4) ((ub(j,i), j= 1,ny), i= 1,nx), &
	       ((vb(j,i), j= 1,ny), i= 1,nx), & 
	       ((wb(j,i), j= 1,ny), i= 1,nx), &
	       ((dum, j= 1,ny), i= 1,nx)
      read (4) (x(i), i= 1,nx)
      read (4) (y(j), j= 1,ny)
      read (4) (cur(i), i= 1,nx)
      xmax = x(nx)
      close (4)
      
      return
120   call error ('genmean$','Error opening field.mean$')
      end
