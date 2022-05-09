!=================================================
program e1
!=================================================

  real, allocatable      :: t(:), x1(:), x2(:), x(:) 
  real  allocatable      :: y1(:), y2(:), y(:)
  real                   :: k1, k2, temp(10)
  real                   :: dt, xf, yf
  integer                :: nt, i,j,k

!=================================================
  
  print*,'Enter the nt, time dimension ===> '
  read *, nt
  allocate(x1(nt),y1(nt),x2(nt),y2(nt),t(nt))
  allocate(x(nt),y(nt)) 
  t(1)=0.0
  t(nt)=1.0
  dt=1.0/real(nt-1)
  do i=1,nt-1
     t(i+1)=t(i)+dt
  end do
!... iniitial conditions

  x1(1)= -1.0
  y1(1)= 0.0
  x2(1)= -1.0
  y2(1)= 1.0

  do i=1,nt-1
     x1(i+1)=x1(i)+dt*y1(i)
     y1(i+1)=dt*(x1(i)+1)+y1(i)
     x2(i+1)=x2(i)+dt*y2(i)
     y2(i+1)=dt*(x2(i)+1)+y2(i) 
  end do


!=================================================
!Checking the final conditions for variables

  xf=sinh(1)-1
  yf=cosh(1)
 
  if ( (x1(nt)-xf) < 1.0e-4) then
     print*,'Gotch you '
     goto 100
  else if ( (x2(nt)-xf) < 1.0e-4) then
     print*,'Well, surviving '
     goto 100
  end if

  k2 = (xf-x1(nt))/(x2(nt)-x1(nt))
  k1 = 1 - k2

  x = k1*x1 + k2*x2
  y = k1*y1 + k2*y2

!... Final check 

  if ( (x(nt)-xf) < 1.0e-4) then
     print*,'Good 0 '
  end if

  if (x(1)+1 < 1.0e-4) then 
     print*,'Good 1'
  end if

100 print*,'Ho-ho-ho'
  open(1,file='e1.dat')
  do i=1,nt
     write(1,'(I4,2E16.7)') i, x(i), y(i)
  end do
  close(1)


end program e1

  
     
