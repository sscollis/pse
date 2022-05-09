!==============================================================================
      subroutine bc (nx, bloc, bcon, ubc, vbc, wbc, isbc)
!==============================================================================
      use ops
      use mf
      implicit none

      integer :: nx, isbc
      real    :: bloc, bcon
      complex :: ubc(nx), vbc(nx), wbc(nx)

      integer :: i, nxi
      real    :: b, coe, dum, ur, ui, vr, vi, wr, wi
      complex :: zero=(0.,0.)
!==============================================================================
      if (isbc .eq. -1) then                 ! no surface BC forcing
         do i= 1,nx
            ubc(i) = zero
            vbc(i) = zero
            wbc(i) = zero
         end do
      elseif (isbc .eq. 0) then              ! suction / blowing
         coe = alog(.01) / bcon **2
         do i= 1,nx
            b = exp(coe * (x(i) - bloc) **2)
            vbc(i) = cmplx(b, 0.)
            ubc(i) = zero
            wbc(i) = zero
         end do
      elseif (isbc .eq. 1) then              ! roughness
         coe = alog(.01) / bcon **2
         do i= 1,nx
            b = exp(coe * (x(i) - bloc) **2)
            ubc(i) = cmplx(-uby(1,i) * b, 0.)
            wbc(i) = cmplx(-wby(1,i) * b, 0.)
            vbc(i) = zero
         end do
      else                                   ! read from "surf.bc"
         open (8, file='surf.bc')
         read (8,*) nxi
	 if (nxi .ne. nx) call error('bc$','nxi .ne. nx in surf.bc$')
         do i= 1,nx
            read (8,*) dum, ur, ui, vr, vi, wr, wi
            ubc(i) = cmplx(ur,ui)
            vbc(i) = cmplx(vr,vi)
            wbc(i) = cmplx(wr,wi)
         end do
	 close (8)
      endif

      return
      end
