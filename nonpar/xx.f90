
! compile:
!
! f90 -r8 -o xx xx.f90 chebyinit.f90 opper.f90
!
!========================================================================================
program xx
!========================================================================================
  real,    external       :: maxim 
  real,    allocatable    :: x(:), y(:), rex(:) 
  real,    allocatable    :: ke(:), kex(:), kexx(:), ke1(:), opi(:), umax(:)
  complex, allocatable    :: u(:,:,:), ulxx(:), vlxx(:), wlxx(:), plxx(:)
  complex, allocatable    :: ux(:), vx(:), wx(:), alpha(:)
  complex, parameter      :: iota=(0.0,1.0)
  real                    :: dumax, pt5=0.5, re
  integer                 :: i,j,k, nx, ny, nin=20
  character(30)           :: fname
  
!========================================================================================
  integer :: i, im, ip, j
  real    :: dm, dmm, dp, dpp, dy, aa, bb
  real    :: h_1, h_2, h_3, gam, del, alp, bet, yc
  real,   allocatable    :: wv1(:), wv2(:)
  real,   allocatable    :: cppx(:), cpx(:), cx(:), cmx(:), cmmx(:) 
  real,   allocatable    :: c2ppx(:), c2px(:), c2x(:), c2mx(:), c2mmx(:) 
  character*1,parameter :: trans='n'
!=========================================================================================

!.... make the streamwise derivative operators
!.... fourth-order-accurate central for a nonuniform mesh
      
      print*, 'Enter the field name, should be the hat field ===> '
      read *, fname
      print*, 'Enter reference Reynolds number ===> '
      read *, re 
      open (1, file=fname, form='unformatted')
      read (1) nx, ny
      allocate(kexx(nx),u(ny,4,nx),x(nx),y(ny),ke(nx),ke1(nx),opi(ny),alpha(nx))
      allocate(ux(ny),vx(ny),wx(ny), rex(nx), kex(nx))
      read (1) (((u(j,k,i), j= 1,ny), i= 1,nx), k= 1,4)
      read (1) (x(i), i= 1,nx)
      read (1) (y(i), i= 1,ny)
      close (1)

      allocate(  wv1(nx), wv2(ny), cppx(nx), cpx(nx), cx(nx), cmx(nx), cmmx(nx) )
      allocate( c2ppx(nx), c2px(nx), c2x(nx), c2mx(nx), c2mmx(nx) )

      do i= 1,nx-1
         wv1(i) = x(i+1) - x(i)
      end do


      print*,'everything is allocated'

      cppx(1) = -wv1(1) / ((wv1(2) + wv1(1)) * wv1(2))
      cpx(1) = (wv1(1) + wv1(2)) / (wv1(1) * wv1(2))
      cx(1) = -1./wv1(1) - 1./(wv1(1) + wv1(2))
      cmx(1) = 0.
      cmmx(1) = 0.
     
      dm = 1. / (wv1(2) + wv1(1))
      c2ppx(1) = 2. * dm / wv1(2)
      c2x(1) = 2. * dm / wv1(1)
      c2px(1) = -c2ppx(1) - c2x(1)
      c2mx(1) = 0.
      c2mmx(1) = 0.

      cppx(nx) = 0.
      cpx(nx) = 0.
      cx(nx) = 1./wv1(nx-1) + 1./(wv1(nx-1) + wv1(nx-2))
      cmx(nx) = -(wv1(nx-1) + wv1(nx-2)) / (wv1(nx-1) * wv1(nx-2))
      cmmx(nx) = wv1(nx-1) / ((wv1(nx-1) + wv1(nx-2)) * wv1(nx-2))

      dm = 1. / (wv1(nx-1) + wv1(nx-2))
      c2x(nx) = 2. * dm / wv1(nx-1)
      c2mmx(nx) = 2. * dm / wv1(nx-2)
      c2mx(nx) = -c2x(nx) - c2mmx(nx)
      c2ppx(nx) = 0.
      c2px(nx) = 0.

      print*,'Entering the loop'

      do i= 3,nx-2
         im = max(i-2,1)
         dmm = wv1(im) + wv1(i-1)
         dm = wv1(i-1)
         dp = wv1(i)
         ip = min(i+1,nx-1)
         dpp = wv1(i) + wv1(ip)

         cppx(i) = -dm*dmm*dp/(dpp**4+(-dp+dmm+dm)*dpp**3+((-dmm-dm)*dp &
             +dm*dmm)*dpp**2-dm*dmm*dp*dpp)
         cpx(i) = dm*dmm*dpp/((dp**3+(dmm+dm)*dp**2+dm*dmm*dp)*dpp-dp**4 &
             +(-dmm-dm)*dp**3-dm*dmm*dp**2)
         cmx(i) = -dmm*dp*dpp/(((dm*dmm-dm**2)*dp+dm**2*dmm-dm**3)*dpp &
             +(dm**2*dmm-dm**3)*dp+dm**3*dmm-dm**4)
         cmmx(i) = dm*dp*dpp/(((dmm**2-dm*dmm)*dp+dmm**3-dm*dmm**2)*dpp &
             +(dmm**3-dm*dmm**2)*dp+dmm**4-dm*dmm**3)
         cx(i) = -(cppx(i) + cpx(i) + cmx(i) + cmmx(i))
     
         c2ppx(i) = -2.*((dmm+dm)*dp-dm*dmm)/(dpp**4+(-dp+dmm+dm)*dpp**3 &
             +((-dmm-dm)*dp+dm*dmm)*dpp**2-dm*dmm*dp*dpp)
         c2px(i) = 2.*((dmm+dm)*dpp-dm*dmm)/((dp**3+(dmm+dm)*dp**2 &
             +dm*dmm*dp)*dpp-dp**4+(-dmm-dm)*dp**3-dm*dmm*dp**2)
         c2mx(i) = -2.*((dp-dmm)*dpp-dmm*dp)/(((dm*dmm-dm**2)*dp+dm**2*dmm &
             -dm**3)*dpp+(dm**2*dmm-dm**3)*dp+dm**3*dmm-dm**4)
         c2mmx(i) = 2.*((dp-dm)*dpp-dm*dp)/(((dmm**2-dm*dmm)*dp+dmm**3 &
             -dm*dmm**2)*dpp+(dmm**3-dm*dmm**2)*dp+dmm**4-dm*dmm**3)
         c2x(i) = -(c2ppx(i) + c2px(i) + c2mx(i) + c2mmx(i))
      end do

!.... Otto's ends

!.... i=2

      h_1=x(2)-x(1)
      h_2=x(3)-x(2)
      h_3=x(4)-x(2)
      gam=(h_1*h_3**3+h_1**2*h_3**2)/ &
            ((h_2**2+h_1*h_2)*(h_3**2-h_1*h_3) &
            -(h_2**3-h_1**2*h_2)*(h_3**2+h_1*h_3))
      del=(h_1-(h_2**2+h_1*h_2)*gam) /(h_3**2+h_1*h_3)
      alp=1.0/h_1*(h_2*gam+h_3*del-1.0)
      bet=-(alp+gam+del)
      cppx(2)=del
      cpx(2)=gam
      cx(2)=bet
      cmx(2)=alp
      cmmx(2) = 0.

      gam=2.0/h_2/((h_2+h_1)-(h_3+h_1)*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      del=-1.0/h_3*(gam*h_2*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      alp=1.0/h_1*(h_2*gam+h_3*del)
      bet=-(alp+gam+del)
      c2ppx(2)=del
      c2px(2)=gam
      c2x(2)=bet
      c2mx(2)=alp
      c2mmx(2)=0.0

!.... i=nx-1

      h_1=x(nx)-x(nx-1)
      h_2=x(nx-1)-x(nx-2)
      h_3=x(nx-1)-x(nx-3)
      gam=(h_1*h_3**3+h_1**2*h_3**2)/ &
            ((h_2**2+h_1*h_2)*(h_3**2-h_1*h_3) &
            -(h_2**3-h_1**2*h_2)*(h_3**2+h_1*h_3))
      del=(h_1-(h_2**2+h_1*h_2)*gam)/(h_3**2+h_1*h_3)
      alp=1.0/h_1*(h_2*gam+h_3*del-1.0)
      bet=-(alp+gam+del)
      cmmx(nx-1)=-del
      cmx(nx-1)=-gam
      cx(nx-1)=-bet
      cpx(nx-1)=-alp
      cppx(nx-1)=0.0

      gam=2.0/h_2/((h_2+h_1)-(h_3+h_1)*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      del=-1.0/h_3*(gam*h_2*(h_2**2+h_1**2)/(h_3**2+h_1**2))
      alp=1.0/h_1*(h_2*gam+h_3*del)
      bet=-(alp+gam+del)
      c2mmx(nx-1)=del
      c2mx(nx-1)=gam
      c2x(nx-1)=bet
      c2px(nx-1)=alp
      c2ppx(nx-1)=0.0

!========================================================================================


      open(71,file='derxx.dat')
      !.... intialize
      call chebyinit(ny,y,opi,nin)
      ke=zero
      do i=1,nx
         rex(i)=sqrt(x(i)*(re))
         do j = 1, ny-1
            ke(i)=ke(i)+opi(j)*(abs(u(j,1,i))**2)
         end do
      end do
      kex=zeero
      kexx=zero
      do i=10,nx-10
         kex(i)=  cppx(i)*ke(i+2)+ cpx(i)*ke(i+1)+ cx(i)*ke(i)+ cmx(i)*ke(i-1)+ cmmx(i)*ke(i-2)
         kexx(i)=c2ppx(i)*ke(i+2)+c2px(i)*ke(i+1)+c2x(i)*ke(i)+c2mx(i)*ke(i-1)+c2mmx(i)*ke(i-2)
         write(71, "(I4,5E16.7)") i, x(i), rex(i), ke(i), kex(i), kexx(i)
      end do
      close(71)
      
!========================================================================================
end program xx
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
