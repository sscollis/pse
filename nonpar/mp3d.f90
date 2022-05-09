program mp3d

complex, allocatable :: u(:,:,:)

real, allocatable :: q(:,:,:), x(:), y(:)

open(10,file='field',form='unformatted')
read(10) nx, ny
allocate( u(ny,4,nx), q(nx,ny,4), x(nx), y(ny) )
read(10) (((u(j,k,i), j= 1,ny), i= 1,nx), k= 1,4)
read(10) (x(i), i= 1,nx)
read(10) (y(i), i= 1,ny)
close(10)

do k = 1, 4
  do i = 1, nx
    do j = 1, ny
      q(i,j,k) = real( u(j,k,i) )
    end do
  end do
end do

open(10,file='grid.dat',form='unformatted')
write(10) nx, ny, 1
write(10) (((x(i), i= 1,nx), j= 1,ny)), &
          (((y(j), i= 1,nx), j= 1,ny)), &
          (((0.0, i= 1,nx), j= 1,ny))
close(10)


open(10,file='q.dat',form='unformatted')
write(10) nx, ny, 1
write(10) 0.0, 0.0, 0.0, 0.0
write(10) ((1.0, i = 1,nx), j= 1, ny), &
          (((q(i,j,k), i = 1,nx), j= 1,ny), k=1,4)
close(10)

stop
end
