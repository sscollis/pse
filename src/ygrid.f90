program ygrid

!....  Test Street's wall normal mapping function

real, allocatable :: eta(:), y(:), yy(:), deta(:), d2eta(:)

ymax = 20.0
ystr = 0
ny = 64

write(*,"('Enter ny, ymax, ystr ==> ')")
read(*,*) ny, ymax, ystr

allocate( eta(ny), y(ny), yy(ny), deta(ny), d2eta(ny) )

nym = ny - 1

pin = acos(-1.) / nym

do i = 1, ny
 eta(i) = cos((i-1) * pin)
end do

do i = 1, ny
  yc = 0.5 * (1. - eta(i))
  y(i) = ymax * ystr * yc / (1. + ystr - yc)
  yy(i) = ymax * ystr * ( 1.0 - eta(i) ) / ( 1.0 + 2.0 * ystr + eta(i) )
  deta(i) = -(2.0 * ystr + 1.0 + eta(i))**2 / (2.0*ymax*ystr*(ystr+1.0))
  d2eta(i) = 0.5*(2.0 * ystr + 1.0 + eta(i))**3 / (ymax*ystr*(ystr+1.0))**2
  write(10,*) eta(i), yc, y(i), deta(i), d2eta(i)
end do

stop
end
