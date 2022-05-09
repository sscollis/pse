!==============================================================================
      subroutine mean
!==============================================================================
      use global
      implicit none

      character*1,parameter :: trans='n'
      integer :: i, j, ip1, ip2, im1, im2, nxn, nyn

      real :: ubyy(ny), div, c1, c2
!==============================================================================

      if (imean.eq.1) call genmf(1)  ! read and interpolate the bl solution

!.... compute mean flow gradients

#ifdef CRAY
      call sgemm (trans, trans, ny, nx, ny, one, opy, ny, ub, ny, &
                  zero, uby, ny)
      call sgemm (trans, trans, ny, nx, ny, one, opy, ny, vb, ny, &
                  zero, vby, ny)
      call sgemm (trans, trans, ny, nx, ny, one, opy, ny, wb, ny, &
                  zero, wby, ny)
#else
      call dgemm (trans, trans, ny, nx, ny, one, opy, ny, ub, ny, &
                  zero, uby, ny)
      call dgemm (trans, trans, ny, nx, ny, one, opy, ny, vb, ny, &
                  zero, vby, ny)
      call dgemm (trans, trans, ny, nx, ny, one, opy, ny, wb, ny, &
                  zero, wby, ny)
#endif

      do i= 1,nx
         ip2 = min(i+2,nx)
         ip1 = min(i+1,nx)
         im1 = max(i-1,1)
         im2 = max(i-2,1)
         do j= 1,ny
            ubx(j,i) = cppx(i)*ub(j,ip2) + cpx(i)*ub(j,ip1) + cx(i)*ub(j,i) &
                     + cmx(i)*ub(j,im1) + cmmx(i)*ub(j,im2)
            vbx(j,i) = cppx(i)*vb(j,ip2) + cpx(i)*vb(j,ip1) + cx(i)*vb(j,i) &
                     + cmx(i)*vb(j,im1) + cmmx(i)*vb(j,im2)
            wbx(j,i) = cppx(i)*wb(j,ip2) + cpx(i)*wb(j,ip1) + cx(i)*wb(j,i) &
                     + cmx(i)*wb(j,im1) + cmmx(i)*wb(j,im2)
         end do
      end do

      call dgemm (trans, trans, ny, 1, ny, one, opyy, ny, ub, ny, &
                  zero, ubyy, ny)

!.... diagnostic

!      i = 1
!      do j = 1, ny
!       write(88,"(4(1pe13.6,1x))") y(j), ub(j,i), uby(j,i), ubyy(j)
!      end do

!.... compute the maximum divergence of the mean-flow

      div = zero

      do i = 1, nx
        do j = 1, ny

          c1 = one/(one + cur(i)*y(j))
          c2 = cur(i) * c1

          div = max( div, c1 * ubx(j,i) + vby(j,i) + c2 * vb(j,i) )
        end do
      end do
      write(*,"('Maximum mean-flow divergence = ',1pe13.6)") div

      return
      end
