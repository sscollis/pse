!=============================================================================!
      subroutine plot (nx, ny, omega, re, xmax, ymax, ystr, u)
!=============================================================================!
      complex u(ny,4,nx)
      dimension w1(nx,ny)
!=============================================================================!
      iter = 0

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = real(u(j,1,i))
         end do
      end do
      open (10, file='ureal.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
 1000 format (2i5, e15.5, i10, 4e15.5)
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
 1001    format (6e20.10)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = aimag(u(j,1,i))
         end do
      end do
      open (10, file='uimag.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = real(u(j,2,i))
         end do
      end do
      open (10, file='vreal.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = aimag(u(j,2,i))
         end do
      end do
      open (10, file='vimag.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = real(u(j,3,i))
         end do
      end do
      open (10, file='wreal.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = aimag(u(j,3,i))
         end do
      end do
      open (10, file='wimag.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = real(u(j,4,i))
         end do
      end do
      open (10, file='preal.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      do j= 1,ny
         do i= 1,nx
            w1(i,j) = aimag(u(j,4,i))
         end do
      end do
      open (10, file='pimag.plot')
      write (10,1000) nx, ny, omega, iter, re, xmax, ymax, ystr
      do j= 1,ny
         write (10,1001) (w1(i,j), i= 1,nx)
      end do

      return
      end
