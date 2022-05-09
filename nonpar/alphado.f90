! compile:
!
! f90 -r8 -o alphado alphado.f90 chebyinit.f90 opper.f90
!
!========================================================================================
program alphado
!========================================================================================
  real,    external       :: maxim 
  real,    allocatable    :: x(:), y(:), rex(:) 
  real,    allocatable    :: ke(:), ke1(:), opi(:), umax(:)
  complex, allocatable    :: u(:,:,:), ul(:), vl(:), wl(:)
  complex, allocatable    :: ux(:), vx(:), wx(:), alpha(:)
  complex, parameter      :: iota=(0.0,1.0)
  real                    :: kex, dumax, pt5=0.5, re
  integer                 :: i,j,k, nx, ny, nin=20
  character(30)           :: fname
  
!========================================================================================


  print*, 'Enter the field name ===> '
  read *, fname
  print*, 'Enter reference Reynolds number ===> '
  read *, re 
  open (1, file=fname, form='unformatted')
  read (1) nx, ny
  allocate(u(ny,4,nx),x(nx),y(ny),ke(nx),ke1(nx),opi(ny),alpha(nx))
  allocate(ul(ny),vl(ny),wl(ny),ux(ny),vx(ny),wx(ny), rex(nx))
  read (1) (((u(j,k,i), j= 1,ny), i= 1,nx), k= 1,4)
  read (1) (x(i), i= 1,nx)
  read (1) (y(i), i= 1,ny)
  close (1)

!.... intialize

  alpha=cmplx(0.0,0.0)

!.... get weights
  call chebyinit(ny,y,opi,nin)

!  open(1,file='opi.dat')
!  read(1,*) (opi(j),j=1,ny)
!  close(1)
     
!.... compute the correction for alpha (following Herbert AIAA-93-3053)
!.... this is the energy alpha

  open(71,file='alpha.keadj')
  ke=0.0; i=1
  do j=1,ny
     ul(j)=u(j,1,i)
     vl(j)=u(j,2,i)
     wl(j)=u(j,3,i)
     ke(1)=ke(1)+opi(j)*(abs(ul(j))**2+abs(vl(j))**2+abs(wl(j))**2)
     ke1(1)=ke(1)+opi(j)*((real(u(j,1,i))**2+aimag(u(j,1,i))**2)+ &
                         (real(u(j,2,i))**2+aimag(u(j,2,i))**2)+ &
                         (real(u(j,3,i))**2+aimag(u(j,3,i))**2) )
  end do
  if(ke(1)-ke1(1) .ge. 1.0e-8) then 
     print*,'Hmmm........'
  end if
!  ke(1)=0.5*ke(1)
  rex(1)=sqrt(x(1)*(re))
  
  do i=2,nx
     rex(i)=sqrt(x(i)*(re))
     ke(i) = zero
     kex   = zero
     do j = 1, ny-1
        ul(j)=u(j,1,i)
        vl(j)=u(j,2,i)
        wl(j)=u(j,3,i)
        
        ux = (ul-u(j,1,i-1))/(x(i)-x(i-1))
	vx = (vl-u(j,2,i-1))/(x(i)-x(i-1))
	wx = (wl-u(j,3,i-1))/(x(i)-x(i-1))
        

        ke(i)=ke(i)+opi(j)*(abs(ul(j))**2+abs(vl(j))**2+abs(wl(j))**2)
!        kex = kex + opi(j) * ( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
!                               conjg(wl(j))*wx(j) )  
     end do
!     ke(i)=0.5*ke(i)
     kex=(ke(i)-ke(i-1))/(x(i)-x(i-1))
     alpha(i) = -iota * kex/ke(i)
     
     write(71, "(I4,6E16.7)") i, x(i), rex(i), ke(i), kex, real(alpha(i)), aimag(alpha(i))
  end do
  close(71)
  close(1)
!  allocate(umax(nx))
!  open(1,file='alpha.uadj')
!  umax(1)=maxim(u(
!  do i=2,nx
!     call maxim( ny, y, ul, umax(i), nint, yint) 
!     dumax = (umax(i-1)-umax(i))/(x(i+1)-x(i))
!     alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
!          sor * (alpha(i,k,n) - iota * dumax / umax(i))
!  end do
  close(71)
!========================================================================================
end program alphado
!========================================================================================

!===================================================================================
function maxim(vector,n)
!===================================================================================  
  real:: maxim, vector(n)
  real:: tmp
  
  integer::  i,n


  tmp=0.0
  do i=1,n
     if ( tmp .le. vector(i) ) then
        tmp=vector(i)
     end if
  end do
  
  maxim=tmp

!===================================================================================
end function maxim
!===================================================================================
