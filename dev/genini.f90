!==============================================================================
      subroutine genini
!
!     Generate the initial profile for PSE
!
!  Notes:
!
!     Modified so that it treats the comment on the first line correctly
!
!==============================================================================
      use global
      implicit none

      integer :: i, j, k, n, iloc
      real :: ai, ar, ur, ui, vr, vi, wr, wi, prr, pri, dum, ampr, ampi
      complex :: ampl

      character*80 :: base='inflow', name=''
      character*1  :: cdum
      logical :: status
!==============================================================================

      allocate( alpha(nx,nz,mt), amp(nx,nz,mt) )
      allocate( ubci(ny,nz,mt), vbci(ny,nz,mt), wbci(ny,nz,mt), &
                pbci(ny,nz,mt) )
      
      alpha = czero; amp = cone
      ubci = czero; vbci = czero; wbci = czero; pbci = czero

!.... read the initial condition(s)

      if (itype.eq.0 .or. itype.eq.4) then             ! linear PSE
	base='inflow'
	iloc = index(base,' ')-1
	n = 1; k = 1
	if (kz(k) .ge. 0) then
	  write(name,"(a,'.+',i1,i1)") base(1:iloc), kz(k), kt(n)
	else
	  write(name,"(a,'.-',i1,i1)") base(1:iloc), abs(kz(k)), kt(n)
	endif
	write(*,"('Reading:  ',a)") name(1:index(name,' ')-1)
	open(4, file=name,status='old',err=100)  
	read(4,*,err=100,end=100) cdum, ampr, ampi, ar, ai
	ampl = cmplx(ampr,ampi)
	alpha(is,k,n) = cmplx(ar,ai)
	do j = 1, ny
	  read(4,*,err=100,end=100) dum, ur, ui, vr, vi, wr, wi, prr, pri
	  ubci(j,k,n) = ampl * cmplx(ur, ui)
	  vbci(j,k,n) = ampl * cmplx(vr, vi)
	  wbci(j,k,n) = ampl * cmplx(wr, wi)
	  pbci(j,k,n) = ampl * cmplx(prr, pri)
	end do
	close(4)
      else if (itype.eq.2 .or. itype.eq.3) then        ! adjoint PSE
	base='outflow'
	iloc = index(base,' ')-1
	n = 1; k = 1
	if (kz(k) .ge. 0) then
	  write(name,"(a,'.+',i1,i1)") base(1:iloc), kz(k), kt(n)
	else
	  write(name,"(a,'.-',i1,i1)") base(1:iloc), abs(kz(k)), kt(n)
	endif
	write(*,"('Reading:  ',a)") name(1:index(name,' ')-1)
	open(4, file=name,status='old',err=100)  
	read(4,*,err=100,end=100) cdum, ampr, ampi, ar, ai
	ampl = cmplx(ampr,ampi)
	alpha(ie,k,n) = cmplx(ar,ai)
	do j = 1, ny
	  read(4,*,err=100,end=100) dum, ur, ui, vr, vi, wr, wi, prr, pri
	  ubci(j,k,n) = ampl * cmplx(ur, ui)
	  vbci(j,k,n) = ampl * cmplx(vr, vi)
	  wbci(j,k,n) = ampl * cmplx(wr, wi)
	  pbci(j,k,n) = ampl * cmplx(prr, pri)
	end do
	close(4)
      else if (itype.eq.1) then        ! Nonlinear PSE
	iloc = index(base,' ')-1
	do n = 1, mt
	  do k = 1, nz
	    if (kz(k) .ge. 0) then
	      write(name,"(a,'.+',i1,i1)") base(1:iloc), kz(k), kt(n)
	    else
	      write(name,"(a,'.-',i1,i1)") base(1:iloc), abs(kz(k)), kt(n)
	    endif
	    inquire(file=name, exist=status)
	    if (status) then
	      write(*,"('Reading:  ',a)") name(1:index(name,' ')-1)
	      open(4, file=name,status='old',err=100)
	      read(4,*,err=100,end=100) cdum, ampr, ampi, ar, ai
	      ampl = cmplx(ampr,ampi)
	      alpha(is,k,n) = cmplx(ar,ai)
	      do j = 1, ny
		read(4,*,err=100,end=100) dum, ur, ui, vr, vi, wr, wi, prr, pri
		ubci(j,k,n) = ampl * cmplx(ur, ui)
		vbci(j,k,n) = ampl * cmplx(vr, vi)
		wbci(j,k,n) = ampl * cmplx(wr, wi)
		pbci(j,k,n) = ampl * cmplx(prr, pri)
	      end do
	    end if
	  end do
	end do
	close(4)
      else
	call error('genini$','Illegal value for itype$')
      end if

!.... Diagnostic

      if (.false.) then
	call calcp
	k = 2; n = 2;
	open(10,file='ic.dat')
	do j = 1, ny
	  write(10,"(10(1pe13.6,1x))") y(j), &
                                       real(ubci(j,k,n)), aimag(ubci(j,k,n)), &
                                       real(vbci(j,k,n)), aimag(vbci(j,k,n)), &
                                       real(wbci(j,k,n)), aimag(wbci(j,k,n)), &
                                       real(pbci(j,k,n)), aimag(pbci(j,k,n))
	end do
	close(10)
      end if

      return
100   call error ('genini$','Error reading initial condition$')
      end

!==============================================================================
      subroutine calcp
!
!     Computes the pressure from the streamwise momentum equation assuming
!     parallel flow.
!
!==============================================================================
      use global
      implicit none
      
      integer :: i, j, k, n, l
      
      complex :: ubcy(ny), ubcyy(ny)
!==============================================================================

      i=1; k=2; n=2

      do j = 1, ny
	ubcy(j) = czero
	ubcyy(j) = czero
	do l = 1, ny
	  ubcy(j) = ubcy(j) + opy(j,l) * ubci(l,k,n)
	  ubcyy(j) = ubcyy(j) + opyy(j,l) * ubci(l,k,n)
	end do
      end do

      do j = 1, ny
	pbci(j,k,n) = one/(iota*alpha(i,k,n)) * (one/re *( -alpha(i,k,n)**2 * &
                      ubci(j,k,n) + ubcyy(j) - beta**2 * ubci(j,k,n) ) + &
                      (iota*omega - iota*alpha(i,k,n)*ub(j,i) - &
                       iota*beta*wb(j,i))*ubci(j,k,n) - uby(j,i)*vbci(j,k,n) &
                     - vb(j,i)*ubcy(j))
      end do

!      do j = 1, ny
!	write(11,"(10(1pe13.6,1x))") y(j), ubcy(j), ubcyy(j)
!      end do

      return
      end subroutine calcp
