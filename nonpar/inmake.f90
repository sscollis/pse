!===============================================================
program inmake
!===============================================================
!
! This little program will add two inflow files and 
! store them in the inflow file format, appropriate only 
! for lns/alns. 
!
!  f90 -r8 -o inmake inmake.f90 
! 
!===============================================================

  real, allocatable      :: y(:), u(:,:), v(:,:), w(:,:), p(:,:)
  real, allocatable      :: y1(:), u1(:,:), v1(:,:), w1(:,:), p1(:,:)
  real, allocatable      :: y2(:), u2(:,:), v2(:,:), w2(:,:), p2(:,:)
  real                   :: rdum
  integer                :: i, j, k, ny
  character(40)          :: file1, file2, file3
  character              :: cdum

!===============================================================
  
  print*, 'Enter the first file name  ==> '
  read *, file1
  print*, 'Enter the second file name ==> '
  read *, file2
  print*, 'Enter ny '
  read *, ny
  allocate(y(ny),u(ny,2),v(ny,2),w(ny,2),p(ny,2))
  allocate(y1(ny),u1(ny,2),v1(ny,2),w1(ny,2),p1(ny,2))
  allocate(y2(ny),u2(ny,2),v2(ny,2),w2(ny,2),p2(ny,2))

!===============================================================

  open(1,file=file1)
  read(1,*) cdum
  do j=1,ny
     read(1,*) y1(j),u1(j,1),u1(j,2),v1(j,1),v1(j,2),&
               w1(j,1),w1(j,2),p1(j,1),p1(j,2) 
  end do
  close(1)
  open(1,file=file2)
  read(1,*) cdum
  do j=1,ny
     read(1,*) y2(j),u2(j,1),u2(j,2),v2(j,1),v2(j,2),&
               w2(j,1),w2(j,2),p2(j,1),p2(j,2) 
  end do
  close(1)   
  
  if (y2(ny)-y1(ny) > 1.0e-6 ) then
     print*,'Incompatible y grids, exiting...'
     stop 666
  end if
  
  if (y2(20)-y1(20) > 1.0e-6 ) then
     print*,'Incompatible y grids, exiting...'
     stop 666
  end if

  do k=1,2
     do j=1,ny
        u(j,k)=u1(j,k)+u2(j,k)
        v(j,k)=v1(j,k)+v2(j,k)
        w(j,k)=w1(j,k)+w2(j,k)
        p(j,k)=p1(j,k)+p2(j,k)
     end do
  end do

  y=y1;
  print*, 'Enter the output file ==> ' 
  read *, file3
  
  open(1,file=file3)
  write(1,*) cdum
  do j=1,ny
     write(1,50) y(j),u(j,1),u(j,2),v(j,1),v(j,2),&
                 w(j,1),w(j,2),p(j,1),p(j,2) 
  end do
  close(1) 
  

50 format(1p,11(e20.13,1x))

!===============================================================
end program inmake
!===============================================================
