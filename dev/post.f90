!==============================================================================
      subroutine post(isign)
!
!     Write out field files for the nonzero PSE modes
!
!     If int <= 1 then no interpolation is performed
!
!     isign =  1 for LNS
!     isign = -1 for ANS
!
!     S. Scott Collis
!
!     Written:  11-2-98
!
!     Revised:  1-29-99
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, idof, nxl
      character*80 :: fname, base='field'
      integer :: isign

!.... variables for interpolation to finer mesh

      integer, parameter :: korder = 3
      real, allocatable :: knot(:), bs_r(:), bs_i(:)
      complex, allocatable :: bs_alpha(:), bs_u(:,:,:)
      real :: dx
      real, allocatable :: xl(:)
      complex, allocatable :: alp(:), ul(:,:,:), ampl(:)
      integer :: il, ix, nx2
      real, external :: BSDER, BSITG
!==============================================================================

      nxl = ie-is+1

!.... Use PSE resolution for field files

      if (int.le.1) then

      write(*,"(/,'Writing field files...',/)")

      do n = 1, mt
	do k = 1, nz
	  if (maxval(abs(u(is:ie,:,k,n,:))) .eq. 0) exit
	  call makename(base,kz(k),kt(n),fname)
	  open  (10, file=fname, form='unformatted')
	  write (10) nxl, ny
	  write (10) (((amp(i,k,n)*u(i,j,k,n,idof), j= 1,ny), i= is,ie), &
                        idof= 1,4)
	  write (10) (x(i), i= is,ie)
	  write (10) (y(j), j= 1,ny)
	  close (10)
	end do
      end do

      else

!.... B-Spline interpolation to finer mesh for field files

      write(*,"(/,'Writing field files on finer mesh...',/)")
      
      nx2 = (nxl-1)*int+1

      do n = 1, mt; do k = 1, nz
      if (maxval(abs(u(is:ie,:,k,n,:))) .eq. 0) exit
      call makename(base,kz(k),kt(n),fname)

      allocate( knot(nxl+korder), bs_r(nxl), bs_i(nxl), bs_alpha(nxl), &
                bs_u(nxl,ny,4), alp(nx2), ul(nx2,ny,4), xl(nx2), ampl(nx2) )

      ul = zero

      call BSNAK( nxl, x(is:ie), korder, knot )
      call BSINT( nxl, x(is:ie), real(alpha(is:ie,k,n)), korder, knot, bs_r )
      call BSINT( nxl, x(is:ie), aimag(alpha(is:ie,k,n)), korder, knot, bs_i )
      bs_alpha = cmplx( bs_r, bs_i )
      do idof = 1, 4
	do j = 1, ny
	  call BSINT( nxl, x(is:ie), real(u(is:ie,j,k,n,idof)), korder, &
                      knot, bs_r )
	  call BSINT( nxl, x(is:ie), aimag(u(is:ie,j,k,n,idof)), korder, &
                      knot, bs_i )
	  bs_u(:,j,idof) = cmplx( bs_r, bs_i )
	end do
      end do

      ix = 0
      do i = is, ie-1
	dx = (x(i+1) - x(i))/real(int)
	do il = 1, int
	  ix = ix + 1
	  xl(ix) = x(i) + real(il-1)*dx
	  alp(ix) = cmplx(BSDER( 0,xl(ix),korder,knot,nxl,real(bs_alpha) ), &
                          BSDER( 0,xl(ix),korder,knot,nxl,aimag(bs_alpha) ) )
	  do idof = 1, 4
	    do j = 1, ny
	      ul(ix,j,idof) = cmplx( &
                      BSDER( 0,xl(ix),korder,knot,nxl,real(bs_u(:,j,idof))), &
                      BSDER( 0,xl(ix),korder,knot,nxl,aimag(bs_u(:,j,idof))) )
	    end do
	  end do

	end do
      end do
      ix = ix + 1
      xl(ix) = x(ie)
      alp(ix) = cmplx(BSDER( 0, xl(ix), korder, knot, nxl, real(bs_alpha)), &
                      BSDER( 0, xl(ix), korder, knot, nxl, aimag(bs_alpha))  )
      do idof = 1, 4
	do j = 1, ny
	  ul(ix,j,idof) = cmplx( &
                      BSDER( 0,xl(ix),korder,knot,nxl,real(bs_u(:,j,idof))), &
                      BSDER( 0,xl(ix),korder,knot,nxl,aimag(bs_u(:,j,idof))) )
	end do
      end do

      ampl = czero
      if (isign.eq.1) then
	do i = 1, nx2
	  ampl(i) = cmplx(BSITG(xl(1),xl(i),korder,knot,nxl,real(bs_alpha)), &
	    BSITG(xl(1),xl(i),korder,knot,nxl,aimag(bs_alpha)) )
	  ampl(i) = exp( iota * ampl(i) )
	end do
      else if (isign.eq.-1) then
	do i = nx2, 1, -1
	  ampl(i) = cmplx(BSITG(xl(nx2),xl(i),korder,knot,nxl,real(bs_alpha)),&
	    BSITG(xl(nx2),xl(i),korder,knot,nxl,aimag(bs_alpha)) )
	  ampl(i) = exp( iota * ampl(i) )
	end do
      else
	call error('post$','Illegal value of isign$')
      end if

      open  (10, file=fname, form='unformatted')
      write (10) nx2, ny
      write (10) (((ampl(i)*ul(i,j,idof), j= 1,ny), i= 1,nx2), idof= 1,4)
      write (10) (xl(i), i= 1,nx2)
      write (10) (y(j), j= 1,ny)
      close (10)

      deallocate( knot, bs_r, bs_i, bs_alpha, bs_u, alp, ul, xl, ampl )

      end do; end do

      end if

      return
      end subroutine post

!==============================================================================
      subroutine makename(base,kz,kt,fname)
!
!     Put a version number on a filename
!
!==============================================================================
      character*80 base, fname

      length = index(base,' ')
      fname = base
      if (abs(kz) .lt. 10 .and. abs(kt) .lt.10) then
	if (kz .ge. 0) then
	  write(fname(length:80),"('.+',i1,i1)") kz, kt
	else
	  write(fname(length:80),"('.-',i1,i1)") abs(kz), kt
	end if
      else if (abs(kz) .lt. 100 .and. abs(kt) .lt.10) then
	if (kz .ge. 0) then
	  write(fname(length:80),"('.+',i2,i1)") kz, kt
	else
	  write(fname(length:80),"('.-',i2,i1)") abs(kz), kt
	end if
      else if (abs(kz) .lt. 10 .and. abs(kt) .lt.100) then
	if (kz .ge. 0) then
	  write(fname(length:80),"('.+',i1,i2)") kz, kt
	else
	  write(fname(length:80),"('.-',i1,i2)") abs(kz), kt
	end if
      else if (abs(kz) .lt. 100 .and. abs(kt) .lt.100) then
	if (kz .ge. 0) then
	  write(fname(length:80),"('.+',i2,i2)") kz, kt
	else
	  write(fname(length:80),"('.-',i2,i2)") abs(kz), kt
	end if
      else
	call error('makename$','version number too large$')
      end if

      return
      end subroutine makename
