!=============================================================================!
      subroutine outf (nx, ny, u, icomp)
     
      complex :: u(ny,4,nx)
      integer :: ip(nx)
!=============================================================================!
      sc = 1.e-12
      do j= 1,ny
         do i= 1,nx
            sc = max(sc, abs(real(u(j,icomp,i))),abs(aimag(u(j,icomp,i))) )
         end do
      end do

      write(*,600) sc
 600  format (e20.5)
     
      scf = 100. / sc
      write(*,*) 'real part'
      ic = 1
 100  continue

      nmax = min(nx, ic*33)
      is = (ic-1) * 33 + 1
      do j= 1,ny
         do i= is,nmax
            up = real(u(j,icomp,i))
            ip(i) = ifix(scf * up + sign(.5, up) )
         end do
         write(601,*) (ip(i), i= is,nmax)
 601     format (33i4)
      end do

      if (nmax .ge. nx) goto 10
      ic = ic + 1
      write(*,*)
      goto 100

 10   continue

      write(*,*) 'imaginary part'
      ic = 1
 101  continue

      nmax = min0 (nx, ic*33)
      is = (ic-1) * 33 + 1
      do j= 1,ny
         do i= is,nmax
            up = aimag(u(j,icomp,i))
            ip(i) = ifix(scf * up + sign(.5, up) )
         end do
         write(601,*) (ip(i), i= is,nmax)
      end do

      if (nmax .ge. nx) goto 20
      ic = ic + 1
      write(*,*)
      goto 101

 20   continue

      return
      end
