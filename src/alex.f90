!==============================================================================
!
! This version solver3 is different from solver1 due to the organization
! and d(omega)/d(alpha) terms.
!
!
!==============================================================================

module ualpha

!==============================================================================
!
! in this module i will maintain all the necessary variables which are
! important for calculating du/d\alpha terms
!
  complex, allocatable ::  Aba(:,:,:), A0a(:,:), ra(:)
  complex, allocatable ::  ua(:,:,:), uad(:), vad(:), wad(:), pad(:)
  complex, allocatable ::  uxa(:), vxa(:), wxa(:)
  complex              ::  omal
!==============================================================================

end module ualpha

!==============================================================================
subroutine solver
!
!     Solve the linear PSE using Chebyshev collocation
!
!     Note: this is the first-order accurate scheme
!
!     S. Scott Collis
!
!     8-21-98:  1) Replaced trapezoidal integration with Chebyschev integration
!               2) Added normalization option (norm=0) Kinetic energy
!                                             (norm=1) Umax
!               3) Added continuity at wall for pressure (ipwbc=2)
!               4) Commented interpolation for ipwbc=1 on output
!     11-2-98:  1) Added field output
!
!    8-9-99     1) fixed the J normalization.
!
!==============================================================================

  use global
  use ualpha
  implicit none
  complex, allocatable     :: alpha1(:), alpha2(:), dotp(:,:), dotpx(:,:), h(:,:)
  complex, allocatable     :: dotp1(:,:), dotpx1(:,:), era(:)
  complex, allocatable     :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:)
  complex, allocatable     :: ulp(:,:,:), ulp1(:,:,:), relamp(:), amp1(:)
  complex                  :: dalpha(nz,nt), alphal
  complex                  :: ul(ny), vl(ny), wl(ny), pl(ny)
  complex                  :: ux(ny), vx(ny), wx(ny), px(ny)
  complex                  :: uy(ny,ndof), uyy(ny,ndof)
  complex                  :: kex, umax(nx), dumax
  complex                  :: temp(10), al, ala

  real, allocatable        :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
                              E(:,:,:), Ex(:,:,:), xl(:), yl(:), erar(:)

  real                     :: ke(nx), ar, ai, ar1, ai1, rdum, dxf, ampr, ampi
  real                     :: c1(ny), c2(ny), err, tmp(10), ampl
  real                     :: rei, df1, df2, df3, df4, df5, xc, xcm, xcmm
  real                     :: ur, ui, vr, vi, wr, wi, prr, pri

!_____ local variables

  real                     :: ubl(ny),  vbl(ny),  wbl(ny)
  real                     :: ubxl(ny), vbxl(ny), wbxl(ny)
  real                     :: ubyl(ny), vbyl(ny), wbyl(ny)

  integer                  :: i, j, k, n, iter, il, aflag, fnorm, idum, nxl ,nyl
  integer                  :: l, num, ldof, l0, m, mdof, m0, info, job=0, icount=0
  integer, allocatable     :: ipvt(:)

!____ wall BCs

  real                     :: bloc = 100.0, bcon = 0.0, bamp = -1.0e-1
  real                     :: coe, bl

  character                :: cdum

  namelist /bc/ bloc, bcon, bamp
!=======================================================================================
!____ allocation block

  rei = one / re

  allocate( u(nx,ny,nz,nt,ndof) )
  u = czero

  allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
            C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
            Ex(ny,ndof,ndof) )

  allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
  allocate( r(ny*ndof), ipvt(ny*ndof) )

!____ alpha variables

  allocate(Aba(ny,ndof,ndof), A0a(ny*ndof,ny*ndof) )
  allocate(ra(ny*ndof), ua(nx,ny,ndof))
  allocate(uad(ny), vad(ny), wad(ny), pad(ny))
  allocate(uxa(ny), vxa(ny), wxa(ny))
!____

  allocate(dotp(nx,ny),  dotpx(nx,ny), era(nx), erar(nx))
  allocate(dotp1(nx,ny), dotpx1(nx,ny), h(nx,10))
  allocate(relamp(nx),amp1(nx))

!____ initialize

  h  =   czero; amp1  = czero
  dotp = czero; dotpx = czero;
  relamp = czero; dotp1 = czero; dotpx1 = czero


!======================================================================================
!.... calculate d(omega)/d(alpha)

  i = is; n = 1; k = 1

!.... the adjoint eigenfunction calculated at the same x location as
!.... regular eigenfunction.

  open(1, file='adjoint.out')
  read(1,*) c, ampr, ampi, ar, ai
  ampl = cmplx(ampr,ampi)
  ala  = cmplx(ar,ai)
  do j = 1, ny
     read(4,*,err=1001,end=1001) rdum, ur, ui, vr, vi, wr, wi, prr, pri
     uad(j) = ampl * cmplx(ur, ui)
     vad(j) = ampl * cmplx(vr, vi)
     wad(j) = ampl * cmplx(wr, wi)
     pad(j) = ampl * cmplx(prr, pri)
  end do
  close(1)

  al=alpha(is,n,k)

!.... regular eigenfunction is located in ubci, vbci, wbci, pbci
!.... alpha is given at i=is, ubci are already premultiplied by ampl

!.... Computing J innerproduct for eigenfunctions
!.... everything here is 2-d, assume adjoint have not been normalize, by j product.
!_____________________________________________________________________________________

     i=is
     do j=1,ny
        temp(9)=( uad(j) * ubci(j,n,k) + vad(j) * vbci(j,n,k) )
        dotp(i,j)=temp(9) * ub(j,i) + uad(j) * pbci(j,n,k) + &
                  pad(j) * ubci(j,n,k) ! - rei * iota*(al + ala) * temp(9)
     end do

!.... compute J at the i=is point:


     h(i,3) = czero
     do j = 1,ny
        h(i,3) = h(i,3) + opi(j) * dotp(i,j)
     end do

!.... compute the cross-kinetic energy <u,u~> ( everything is in 2-d right now)

     do j=1,ny
        dotp(i,j) = ubci(j,n,k) * uad(j) + vbci(j,n,k) * vad(j)
     end do

     h(i,2) = czero
     do j=1,ny
        h(i,2) = h(i,2) + opi(j) * dotp (i,j)
     end do

!.... finaly compute the group velocity


     omal = h(i,3)/h(i,2)


!____ Read in alpha from the PSE run if norm = 3
!======================================================================================
!.... reading pse alpha

   allocate(alpha1(nx),alpha2(nx))
   alpha1=czero; alpha2=czero

  if ( norm .eq. 3 ) then
     open(1,file='alpha.pse')
     do i=is,ie
        read(1,*) idum, rdum, rdum, ar1, ai1, ar, ai
        alpha2(i)= cmplx(ar,ai)
     end do
     close(1)
  end if
!======================================================================================
!... reading adjoint pse alpha

  if ( norm .eq. 3 ) then
     open(1,file='alpha.adj')
     do i=is,ie
        read(1,*) idum, rdum, rdum, ar1, ai1, ar, ai
        amp1(i) =  cmplx(ar1,ai1)
        alpha1(i)= cmplx(ar,ai)
        alpha(i,1,1)= -alpha1(i)
     end do
     close(1)

     print*,'alpha(is,1,1) = ',alpha(is,1,1)
     write(12,'(I4,2E16.7)') (i, abs(amp1(i)), aimag(alpha1(i)),i=is,ie)
  end if

!======================================================================================
!... Calculate the difference between adjoint and pse alpha:

 write(17,'(I4,2E16.7)') (i, real(alpha1(i)+alpha2(i)), &
                             aimag(alpha1(i)+alpha2(i)), i=is,ie)

!.... only consider only one mode for linear analysis

  i = is; n = 1; k = 1

!.... initialize the field

  do j = 1, ny
     u(i,j,k,n,1) = ubci(j,k,n)
     u(i,j,k,n,2) = vbci(j,k,n)
     u(i,j,k,n,3) = wbci(j,k,n)
     u(i,j,k,n,4) = pbci(j,k,n)
  end do

  ul = u(i,:,k,n,1)
  vl = u(i,:,k,n,2)
  wl = u(i,:,k,n,3)

  ke =    zero
  kex =   zero
  umax =  zero
  dumax = zero

  ke(i) = zero
  do j=1,ny
     ke(i)=ke(i)+opi(j)*(abs(ul(j))**2+abs(vl(j))**2+abs(wl(j))**2)
  end do

!==================================================================================
!...  based on the hat fields and alpha
!...  better when one does derivatives,
!...  therefore this should be used !!!
!...  Here i is no longer neded to be is,
!... therefore is it as dummy variable.
!==================================================================================

  if ( norm .eq. 3) then
     if ( old .eq. 0 ) then
        open(1,file='afield.hat',form='unformatted')
     else if ( old .eq. 1) then
        open(1,file='afield.hato',form='unformatted')
     end if
     read(1) nxl, nyl
     allocate(ulp(nxl,nyl,4),xl(nxl),yl(nyl))
     read(1) (((ulp(i,j,num),j=1,nyl),i=is,ie),num=1,4)
     read(1) (xl(i),i=is,ie)
     read(1) (yl(j),j=1,nyl)
     close(1)

     open(1,file='field.apse',form='unformatted')
     read(1) nxl, nyl
     allocate (ulp1(nxl,nyl,4))
     read(1) (((ulp1(i,j,num),j=1,nyl),i=is,ie),num=1,4)
     read(1) (xl(i),i=is,ie)
     read(1) (yl(j),j=1,nyl)
     close(1)


!.... make sure that xl and yl are the same as x and y
!.... notice, I am only going to do second order derivatives
!__________________________________________________________________________________

     if (ie-is+1 .ne. nxl) then
        print*,'Check x mesh...'
        print*,ie,nxl
     end if
     if (ny .ne. nyl) then
        print*,'Check y mesh...'
        print*,ny,nyl
        stop 666
     end if
!... Comparing floats can be a bit ticky, I will use absolute error

     do i=is,ie
        if ( xl(i)-x(i) .ge. 1.0e-8 ) then
           print*,'xl-x = ', xl(i)-x(i)
           print*,'X meshes do not match, exiting...'
           stop 666
        end if
     end do
     do j=1,ny
        if (yl(j)-y(j) .ge. 1.0e-8 ) then
           print*,'yl-y = ', yl(j)-y(j)
           print*,'y meshes do not match, exiting...'
           stop 666
        end if
     end do
!_____________________________________________________________________________________
!.... initialize dotp

     i=is
     do j=1,ny
        temp(9)=(ulp(i,j,1)*u(i,j,1,1,1)+ulp(i,j,2)*u(i,j,1,1,2))
        dotp(i,j)=temp(9)*ub(j,i)+ulp(i,j,1)*u(i,j,1,1,4)+ &
             ulp(i,j,4)*u(i,j,1,1,1) !-rei*iota*(alpha(i,k,n)-alpha1(i))*temp(9)
     end do
     dotp1(i,:)=dotp(i,:) !... watch out....

!.... compute J at the i=is point:

     h(i,3)=czero
     do j=1,ny
        h(i,3)=h(i,3)+opi(j)*amp1(i)*dotp(i,j)
     end do
!_____________________________________________________________________________________
!.... checking amplitude at is:

     relamp(is)=amp(is,k,n)*amp1(is)


     print*, 'Amplitude at is is         ====> ',         amp(is,k,n)
     print*, 'Adjoint amplitude at is is ====> ',         amp1(is)
     print*, 'Relamplitude at ie is      ====> ',         relamp(is)
!____________________________________________________________________________________
  end if

  i=is

  if (norm.eq.0) then
     ke(i) = zero
     do j = 1, ny
        ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
             abs(wl(j))**2)
     end do
  else if (norm.eq.1) then
     call findumax( ny, y, ul, umax(i), nint, yint )
  else if (norm.eq.2) then
     ke(i) = zero
     do j = 1, ny
        ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
     end do
  end if
!__________________________________________________________________________________________
!... modified slightly to agree in format with the new version of lns solver

  open(9,file='pse.dat')
  write(9,*) ie-is+1
  write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                                amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
  call flush(9)


!=========================================================================================
!
! Initialize somehow the du/dalpha variable!
!
  ua=czero


!=========================================================================================
!
!
!.... main marching loop

  print *, ie
  open(1,file='dotp.dat')
  open(2,file='hs.dat')

  do i = is+1, ie

     cpu2 = second()
     write(*,"(/,'  i       x(i)          cpu        total cpu')")
     write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
     cpu = cpu2

     alpha(i,k,n) = alpha(i-1,k,n)

     c1(:) = one/(one + cur(i)*y(:))
     c2(:) = cur(i) * c1(:)

!.... the initial velocities are of cource taken from the inflow.+01

     ul = u(i-1,:,k,n,1)
     vl = u(i-1,:,k,n,2)
     wl = u(i-1,:,k,n,3)
     pl = u(i-1,:,k,n,4)

!.... compute the mode amplitude

     amp(i,k,n) = czero
     do il = is, i-1
        amp(i,k,n) = amp(i,k,n) + pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
                     (x(il+1)-x(il))
     end do
     amp(i,k,n) = exp( iota * amp(i,k,n) )


     ke(i) = zero
     do j = 1, ny
        ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + &
             abs(wl(j))**2)
     end do

!=========================================================================================
     if (norm.eq.0) then
        write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
             & '          ke','           dalpha')")
        write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), abs(amp(i,k,n)), ke(i), zero
     else if (norm.eq.1) then
        call findumax( ny, y, ul, umax(i), nint, yint )
        write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
             & '         |umax|','         dalpha')")
        write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), abs(amp(i,k,n)), abs(umax(i)), zero
     else if (norm.eq.2) then
        ke(i) = zero
        do j = 1, ny
           ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
        end do
        write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
             & '          ke','           dalpha')")
        write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), abs(amp(i,k,n)), ke(i), zero
     else if ( norm .eq. 3) then
        write(*,"(/,' iter','    alpha_r','      alpha_i','        |amp|', &
             & '          J','           J^,x/J^')")
        write(*,"('->',i2,1x,6(1pe13.6,1x))") 0, alpha(i,k,n), abs(amp(i,k,n)), h(i,3), zero
     end if

!=========================================================================================

     do iter = 1, niter ! loop on alpha

!.... compute the streamwise derivative of alpha

        dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )

        write(15,"(I4,2E16.7)") i, real(alpha(i,k,n)), aimag(alpha(i,k,n))

!.... Bertollotti sets this to zero (no advantage)

!       dalpha(k,n) = zero

!.... Evaluate the matrices (for Backward Euler)

        alphal = alpha(i,k,n)
        ubl    = ub(:,i)
        vbl    = vb(:,i)
        wbl    = wb(:,i)
        ubxl   = ubx(:,i)
        vbxl   = vbx(:,i)
        wbxl   = wbx(:,i)
        ubyl   = uby(:,i)
        vbyl   = vby(:,i)
        wbyl   = wby(:,i)

        G=zero; A=zero; B=zero; C=zero; D=zero; E=zero; Ex=zero
        Gb=czero; Ab=czero

        Aba=czero

        G(:,1,1) = one
        G(:,2,2) = one
        G(:,3,3) = one

        A(:,1,1) = c1(:) * ubl
        A(:,1,2) = -rei * two * c1(:) * c2(:)
        A(:,1,4) = c1(:)
        A(:,2,1) = rei * two * c1(:) * c2(:)
        A(:,2,2) = c1(:) * ubl
        A(:,3,3) = c1(:) * ubl
        A(:,4,1) = c1(:)

        B(:,1,1) = vbl - rei * c2(:)
        B(:,2,2) = vbl - rei * c2(:)
        B(:,2,4) = one
        B(:,3,3) = vbl - rei * c2(:)
        B(:,4,2) = one

        C(:,1,1) = wbl
        C(:,2,2) = wbl
        C(:,3,3) = wbl
        C(:,3,4) = one
        C(:,4,3) = one

!.. possible sign problem, +rei*c2(:)**2 ?
!.. I will write + sign

        D(:,1,1) = c1(:) * ubxl + c2(:) * vbl + rei * c2(:)**2
        D(:,1,2) = ubyl + c2(:) * ubl
        D(:,2,1) = c1(:) * vbxl - two * c2(:) * ubl
        D(:,2,2) = vbyl + rei * c2(:)**2
        D(:,3,1) = c1(:) * wbxl
        D(:,3,2) = wbyl
        D(:,4,2) = c2(:)

        Ex(:,1,1) = c1(:)**2 * rei
        Ex(:,2,2) = c1(:)**2 * rei
        Ex(:,3,3) = c1(:)**2 * rei

        E(:,1,1) = rei
        E(:,2,2) = rei
        E(:,3,3) = rei

        Gb = -iota * kt(n) * omega * G + iota * alphal * A + &
              iota * kz(k) * beta * C + D - &
              Ex * ( iota * dalpha(k,n) - alphal**2 ) + &
              kz(k)**2 * beta**2 * E


!.... fix the streamwise pressure gradient

        if (ipfix.eq.1 .or. (ipfix.eq.2 .and. i.eq.2) ) A(:,:,4) = zero

        Ab = A - two * iota * alphal * Ex
        Aba = Ab


!.... Build the LHS matrix

        A0 = czero

        do ldof = 1, ndof
           do mdof = 1, ndof
              do l = 1, ny
                 l0 = (l-1)*ndof
                 do m = 1, ny
                    m0 = (m-1)*ndof
                    A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + &
                         B(l,ldof,mdof) * opy(l,m) - &
                         E(l,ldof,mdof) * opyy(l,m)
                 end do
                 m0 = (l-1)*ndof
                 A0(l0+ldof,m0+mdof) = A0(l0+ldof,m0+mdof) + Gb(l,ldof,mdof) + &
                      Ab(l,ldof,mdof) / (x(i)-x(i-1))
              end do
           end do
        end do



!.... The boundary conditions for A0a could be a little bit different,
!.... depending on the boundary conditions on the A0 itself.
!.... Apply the boundary conditions to the LHS
!_________________________________________________________________________
!  next few lines
!  Set matrix A0 = 4ny x 4ny to the following form:
!
!  1 0 0 0 0 ...
!  0 1 0 0 0 ...
!  0 0 1 0 0 ...
!  0 0 0 0 0 ...
!  . . .
!  : : :
!
! This correspond to setting u, v, w, at the wall to the values in the
! right handside which should be zero
! then ipwbc ( pressure wall boundary conditions ) are encountered
!_________________________________________________________________________
!   small review of the boundary  conditions:
!   ipwbc = 0 ads line to a matrix A0
!
!  1 0 0 0 0 ...
!  0 1 0 0 0 ...
!  0 0 1 0 0 ...
!  0 0 0 opy(1,1) 0 0 0 opy(1,2) 0 0 0 opy(1,3) ...
!  . . .
!  : : :
!
!  this accounts to setting the dp/dy | (y=0, first node ) = r(4)
!  which is presumably zero.
!
!   ipwbc = 1 ads line to a matrix A0
!
!  1            0           0 0 0 ...
!  0            1           0 0 0 ...
!  0            0           1 0 0 ...
!  0 -rei * c2(1) * opy(1,1) - rei * opyy(1,1) 0 opy(1,1) 0 0 0 opy(1,2) ...
!  . . .
!  : : :
!
!  This accounts for dp/dy=1/R(c2 dv/dy+ d^2v/dy^2) which is direct
!  consequence of the y momentum equation at the wall where nonslip
!  condition is imposed.
!
!  ipwbc = 2 ads to the matrix A0

!  1 0 0 0 0 ...
!  0 1 0 0 0 ...
!  0 0 1 0 0 ...
!  c1(1)\dx  opy(1,1) + c2(1)  iB(-1)   0  ...
!  . . .
!  : : :
!
!  This ensures that continuity is satisfied
!
!
!_________________________________________________________________________

        l = 1
        l0 = (l-1)*ndof

        A0(l0+1,:) = czero
        A0(l0+1,l0+1) = cone

        A0(l0+2,:) = czero
        A0(l0+2,l0+2) = cone

        A0(l0+3,:) = czero
        A0(l0+3,l0+3) = cone

        A0(l0+4,:) = czero

        if (ipwbc .eq. 0) then
           do m = 1, ny
              m0 = (m-1)*ndof
              A0(l0+4,m0+4) = opy(l,m)
           end do
        else if (ipwbc .eq. 1) then
           do m = 1, ny
              m0 = (m-1)*ndof
              A0(l0+4,m0+2) = -rei * c2(l) * opy(l,m) - rei * opyy(l,m)
              A0(l0+4,m0+4) = opy(l,m)
           end do
        else if (ipwbc .eq. 2) then
           do m = 1, ny
              m0 = (m-1)*ndof
              A0(l0+4,m0+2) = opy(l,m)
           end do
           A0(l0+4,m0+2) = A0(l0+4,m0+2) + c2(l)
           A0(l0+4,m0+3) = -iota * kz(k) * beta
           A0(l0+4,m0+1) = c1(l) / (x(i)-x(i-1))
        else
           call error('solver$','Illegal value of ipwbc$')
        end if

!_________________________________________________________________________
!
!  Boundary conditions in the free stream (at the upper boundary
!
!                     4ny-4  4ny-3  4ny-2  4ny-1  4ny
!                       |
!   4ny-3     .... 0    0      1      0      0     0
!
!   4ny-2
!
!   if ivbc = 0 then
!
!
!
!                     4ny-4     4ny-3     4ny-2     4ny-1      4ny
!                       |
!   4ny-3     .... 0    0         1         0         0         0
!
!   4ny-2     .... 0    0         0         1         0         0
!
!   4ny-1     .... 0    0         0         0         1         0
!
!
!   This accounts for the homogeneous condition on the outflow
!
!  if ivbc = 1 then
!
!
!                            4ny-4  4ny-3   4ny-2   4ny-1   4ny
!                              |
!   4ny-3         ....         0      0       1       0      0    0
!
!   4ny-2 ....  opy(ny, ny-1)  0      0       0  opy(ny,ny)  0    0
!
!   4ny-1     ....             0      0       0       0      1    0
!
!   This acounts for the Neumann condition on the outflow. dv/dy=0
!
!  if ipbc = 0 then the pressure at the upper boundary is given using
!  dirichlet boundary condition, otherwise, ipbc = 1, Neumann bc are used.
!
!_________________________________________________________________________
!
!...  if the ipwbc = 0 and ivbc = 0, and ipbc = 0  then the same
!     boundary conditions will be imposed on the du/dalpha.
!_________________________________________________________________________

        l = ny
        l0 = (l-1)*ndof

        A0(l0+1,:) = czero
        A0(l0+1,l0+1) = cone

        if (ivbc .eq. 0) then
           A0(l0+2,:) = czero
           A0(l0+2,l0+2) = cone
        else if (ivbc .eq. 1) then
           A0(l0+2,:) = czero
           do m = 1, ny
              m0 = (m-1)*ndof
              A0(l0+2,m0+2) = opy(l,m)
           end do
        else
           call error('solver$','Illegal value of ivbc$')
        end if

        A0(l0+3,:) = czero
        A0(l0+3,l0+3) = cone

        A0(l0+4,:) = czero
        if (ipbc .eq. 0) then
           A0(l0+4,l0+4) = cone
        else if (ipbc .eq. 1) then
           do m = 1, ny
              m0 = (m-1)*ndof
              A0(l0+4,m0+4) = opy(l,m)
           end do
        else
           call error('solver$','Illegal value of ipbc$')
        end if

!.... Same boundary conditions are used for A0a
!.... The matrices are identical!
!.... do it before solving, the Lapack routines may modify the matrix A0

        A0a=A0

!.... Build the RHS

        r = czero

        do ldof = 1, ndof
           do l = 1, ny
              l0 = (l-1)*ndof
              do mdof = 1, ndof
                 r(l0+ldof) = r(l0+ldof) + Ab(l,ldof,mdof) * &
                      u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
              end do
           end do
        end do


!.... Enforce the boundary conditions on the RHS

        l = 1
        l0 = (l-1)*ndof
        r(l0+1) = czero

!.... try modelling a suction/blowing slot

        if (bcon.eq.zero) then
           r(l0+2) = czero
        else
           coe = alog(.01) / bcon **2
           bl = bamp * exp(coe * (x(i) - bloc) **2)
           write(80,*) x(i), bl
           r(l0+2) = cmplx( bl, zero )
        end if

        r(l0+3) = czero
        if (ipwbc.eq.2) then
           r(l0+4) = c1(l) * u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
        else
           r(l0+4) = czero
        end if

        l = ny
        l0 = (l-1)*ndof
        r(l0+1) = czero
        r(l0+2) = czero
        r(l0+3) = czero
        r(l0+4) = czero

!.... Solve the system

#ifdef CRAY
        call CGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
        if (info.ne.0) call error('solver$','Singular matrix$')
        call CGESL(A0,ny*ndof,ny*ndof,ipvt,r,0)
#else
!       call ZGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
!       if (info.ne.0) call error('solver$','Singular matrix$')
!       call ZGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
        call ZGETRF( ny*ndof, ny*ndof, A0, ny*ndof, ipvt, info)
        if (info.ne.0) call error('solver$','Singular matrix$')
        call ZGETRS( 'N', ny*ndof, 1, A0, ny*ndof, ipvt, r, ny*ndof, info)
#endif

!.... Update the solution

        do ldof = 1, ndof
           do l = 1, ny
              l0 = (l-1)*ndof
              u(i,l,k,n,ldof) = r(l0+ldof)
           end do
        end do

        ul = u(i,:,k,n,1)
        vl = u(i,:,k,n,2)
        wl = u(i,:,k,n,3)
        pl = u(i,:,k,n,4)

!.... compute some streamwise derivatives

        ux = (ul - u(i-1,:,k,n,1))/(x(i)-x(i-1))
        vx = (vl - u(i-1,:,k,n,2))/(x(i)-x(i-1))
        wx = (wl - u(i-1,:,k,n,3))/(x(i)-x(i-1))

        ke(i) = zero
        kex   = zero

!_______________________________________________________________________
!
!.... Build The RHS for the alpha equation

        ra=czero

!.... First construct the regular right hand side

        do ldof = 1, ndof
           do l = 1, ny
              l0 = (l-1)*ndof
              do mdof = 1, ndof
                 ra(l0+ldof) = ra(l0+ldof) + Aba(l,ldof,mdof) * &
                      ua(i-1,l,mdof) / (x(i)-x(i-1))
              end do
           end do
        end do

!.... Adding the extra forcing sources
!.... -iAu -2*alpha Ex u +2*iota*Ex ux + iota*d(omega)\d(alpha)G u
!.... all the derivatives are taken around alpha_0

!.... reusing matrix Aba
!.... i could just write Aba = -iota * Aba
!....  but thats will bring confusion

        Aba=czero

        Aba = -iota * A - two  * alphal * Ex + iota * omal * G

        do ldof =1, ndof
           do l =1, ny
              l0 = (l-1)*ndof
              do mdof =1,ndof
                 ra(l0+ldof) = ra(l0+ldof) + &
                      two * iota * Ex(l,ldof,mdof)* &
                      (u(i,l,k,n,mdof)-u(i-1,l,k,n,mdof))/(x(i)-x(i-1)) + &
                      Aba(l,ldof,mdof) *  u(i,l,k,n,mdof)
              end do
           end do
        end do


!.... Enforce the boundary conditions on the RHS

        l = 1
        l0 = (l-1)*ndof
        ra(l0+1) = czero
        ra(l0+2) = czero
        ra(l0+3) = czero

!.... The continuity condition can be imposed on the du/dalpha as well

        if (ipwbc.eq.2) then
           ra(l0+4) = c1(l) * u(i-1,l,k,n,mdof) / (x(i)-x(i-1))
        else
           ra(l0+4) = czero
        end if

        l = ny
        l0 = (l-1)*ndof
        ra(l0+1) = czero
        ra(l0+2) = czero
        ra(l0+3) = czero
        ra(l0+4) = czero

!.... What about trying the continuity at the outer boundary?
!.... Solve the system

#ifdef CRAY
        call CGEFA(A0a,ny*ndof,ny*ndof,ipvt,info)
        if (info.ne.0) call error('solver$','Singular matrix$')
        call CGESL(A0a,ny*ndof,ny*ndof,ipvt,r,0)
#else
!       call ZGEFA(A0a,ny*ndof,ny*ndof,ipvt,info)
!       if (info.ne.0) call error('solver$','Singular matrix$')
!       call ZGESL(A0a,ny*ndof,ny*ndof,ipvt,ra,job)
        call ZGETRF( ny*ndof, ny*ndof, A0, ny*ndof, ipvt, info)
        if (info.ne.0) call error('solver$','Singular matrix$')
        call ZGETRS( 'N', ny*ndof, 1, A0, ny*ndof, ipvt, ra, ny*ndof, info)
#endif

!.... Update the solution

        do ldof = 1, ndof
           do l = 1, ny
              l0 = (l-1)*ndof
              ua(i,l,ldof) = ra(l0+ldof)
           end do
        end do


!.... compute some streamwise derivatives

        uxa = (ua(i,:,1) - ua(i-1,:,1))/(x(i)-x(i-1))
        vxa = (ua(i,:,2) - ua(i-1,:,2))/(x(i)-x(i-1))
        wxa = (ua(i,:,3) - ua(i-1,:,3))/(x(i)-x(i-1))

!.... Finished computing ua for one iteration, now
!.... one can use it to compute the gradient.
!____________________________________________________________________________________


!====================================================================================

        do j = 1, ny
           ke(i) = ke(i) + opi(j)*( abs(ul(j))**2 + abs(vl(j))**2 + &
                abs(wl(j))**2 )
           kex = kex + opi(j)*( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
                conjg(wl(j))*wx(j) )
        end do

        if (norm.eq.0) then

!.... compute the correction for alpha (following Herbert AIAA-93-3053)

           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                sor * (alpha(i,k,n) - iota * kex / ke(i))
        else if (norm.eq.1) then

!.... compute the correction for alpha (following Bertolotti)

           call findumax( ny, y, ul, umax(i), nint, yint )
           dumax = (umax(i)-umax(i-1))/(x(i)-x(i-1))

           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                              sor * (alpha(i,k,n) - iota * dumax / umax(i))
        else if (norm.eq.2) then
!====================================================================================
!.... Computing the gradient
!
!                                                         -1
!                               < U,x U* >   | < U,x U* >|
!  alpha(new)-alpha(old) =   -  ---------- * |---------- |
!                                < U, U* >   | < U, U* > |,alpha
!
!
! Notice the conjugation operation and differentiation with respect
! to alpha can be interchanged.
!
! ke  =  < U, U* >
! kex =  < U,x U* >
!
!====================================================================================

           ke(i) = zero
           kex   = zero
           do j = 1, ny
              ke(i) = ke(i) + opi(j)*( abs(ul(j))**2  )
              kex = kex + opi(j)*( conjg(ul(j))*ux(j) )
           end do

           h(i,6) = czero
           dotp   = czero
           do j=1,ny
              dotp(i,j) = ua(i,j,1)*conjg(ul(j)) + ul(j)*conjg(ua(i,j,1))
              h(i,6)    = h(i,5) + opi(j) * dotp(i,j)
           end do


           h(i,7) = czero
           dotp   = czero
           do j=1,ny
              dotp(i,j) = uxa(j) * conjg(ul(j)) + ux(j) * conjg(ua(i,j,1))
              h(i,7)    = h(i,7) + opi(j) * dotp(i,j)
           end do

           h(i,8) = czero
           h(i,8) = cone/(h(i,7) - kex * h(i,6)/ke(i))

           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                sor*(alpha(i,k,n) - h(i,8) * kex )

        else if ( norm .eq. 3) then
!====================================================================================
!.... compute alpha making sure that jx hat  norm is constant, dj/dx=0
!.... ignore all the second order derivatives such as u,xx


           dotp(i,:) = czero
           do j=1,ny
              temp(9)   = (ulp(i,j,1)*u(i,j,1,1,1)+ulp(i,j,2)*u(i,j,1,1,2))
              dotp(i,j) = temp(9)*ub(j,i)+ulp(i,j,1)*u(i,j,1,1,4)+ &
                          ulp(i,j,4)*u(i,j,1,1,1) ! -rei*iota*(alpha(i,k,n)-alpha1(i))*temp(9)
              write(1,*) dotp(i,j)
           end do

!.... making streamwise derivatives, Always have to do one sided difference,
!.... I do not know the other side. i .ne ie, Sanity check

           dotpx(i,:)=czero
           dxf=x(i)-x(i-1)

           if (.true.) then
              dotpx(i,:)=(dotp(i,:)-dotp(i-1,:))/dxf
           end if

           h(i,1)=czero; h(i,2)=czero
           do j=1,ny
              h(i,1) = h(i,1)+opi(j)*dotp(i,j)
              h(i,2) = h(i,2)+opi(j)*dotpx(i,j)
           end do

           write(2,*) h(i,1), h(i,2)

           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                sor*(alpha(i,k,n)-iota*h(i,2)/h(i,1))

           write(16,"(I4,2E16.7)") i, real(h(i,2)/h(i,1)), aimag(h(i,2)/h(i,1))
        end if

!.... update the mode amplitude

        amp(i,k,n) = czero
        do il = is, i-1
           amp(i,k,n) = amp(i,k,n) + &
                pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
                (x(il+1)-x(il))
        end do
        amp(i,k,n) = exp( iota * amp(i,k,n) )

!.... Calculating J and Jx not hat functions
        if ( norm .eq. 3 ) then
           relamp(i)=amp(i,k,n)*amp1(i)

!... calculating real jx. this has to be the same as the one from jxmaker
!... it is exact up to a constant camp

           dotp1(i,:) = czero
           h(i,3)     = czero
           do j=1,ny
              dotp1(i,j) = relamp(i)*dotp(i,j)
           end do

           do j=1,ny
              h(i,3) = h(i,3)+opi(j)*dotp1(i,j)
           end do

!... Now calculate J,x
!... Here expression has a small error.


           dotpx1(i,:) = czero
           dotpx1(i,:) = dotpx(i,:)*relamp(i)+iota*(alpha(i,k,n)+alpha1(i))*&
                dotp(i,:)*relamp(i)
           h(i,5)=czero
           do j=1,ny
              h(i,5) = h(i,5)+opi(j)*dotpx1(i,j)
           end do

        end if

!====================================================================================

!.... convergence check
!.... for the Jx norm it does not matter
!.... since we do only one iteration

!====================================================================================
100 format ('->',i2,1x,6(1pe13.6,1x))

        if (norm.eq.0) then
!_____________________________________________________________________________________

          write(*,100) iter, alpha(i,k,n), abs(amp(i,k,n)), ke(i), &
                       abs(kex/ke(i))/abs(alpha(i,k,n))
          if (abs(kex/ke(i)) .lt. tol) exit
!_____________________________________________________________________________________

        else if (norm.eq.1) then
!_____________________________________________________________________________________

          write(*,100) iter, alpha(i,k,n), abs(amp(i,k,n)), abs(umax(i)), &
                             abs(dumax/umax(i))/abs(alpha(i,k,n))
          if (abs(dumax/umax(i)) .lt. tol) exit
!_____________________________________________________________________________________

        else if (norm.eq.2) then
!_____________________________________________________________________________________

          write(*,100) iter, alpha(i,k,n), abs(amp(i,k,n)), ke(i), &
                       abs(kex/ke(i))/abs(alpha(i,k,n))
          if (abs(kex/ke(i)) .lt. tol) exit
!_____________________________________________________________________________________

       else if ( norm .eq. 3)  then
!_____________________________________________________________________________________

          write(*,100) iter, alpha(i,k,n), abs(amp(i,k,n)), abs(h(i,3)), abs(h(i,2)/h(i,1))
          temp(4)=h(i,2)/h(i,1)
          if (abs(temp(4)) .lt. tol ) exit
!_____________________________________________________________________________________

       end if
!=====================================================================================

    end do! loop on iter
    write(15,"(I4,2E16.7)") i, real(alpha(i,k,n)), aimag(alpha(i,k,n))
    write(16,"(I4,2E16.7)") i, real(h(i,2)/h(i,1)), aimag(h(i,2)/h(i,1))


!       write(*,*) 'i, x(i), alpha(i,k,n)', i, x(i), alpha(i,k,n)
!=============================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=============================================================================
    write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), &
                                 abs(kex/ke(i)), amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
    call flush(9)

 end do! loop on i
 close(1)
 close(2)

!.... recompute J, in case norm other then 3 is selected
!.... make sure relamp is well defined, I probably repeat but anyway:
 if ( norm .eq. 3) then
    relamp(:)=amp1(:)*amp(:,k,n)

    dotp(:,:)=czero
    h(:,3)    = czero
    do i=is,ie
       do j=1,ny
          temp(9)=(ulp(i,j,1)*u(i,j,1,1,1)+ulp(i,j,2)*u(i,j,1,1,2))
          dotp(i,j)=temp(9)*ub(j,i)+ulp(i,j,1)*u(i,j,1,1,4)+ &
               ulp(i,j,4)*u(i,j,1,1,1) ! -rei*iota*(alpha(i,k,n)-alpha1(i))*temp(9)
          dotp(i,j)=relamp(i)*dotp(i,j)
          h(i,3)=h(i,3)+opi(j)*dotp(i,j)
       end do
    end do

    do i=is,ie
       write(*,*) i,abs(h(i,3)),abs(relamp(i))
    end do

    close(1)
    open(1,file='jpse.out')
    write(1,'(I4,4E16.7)') (i, x(i), real(h(i,3)), &
         aimag(h(i,3)), abs(h(i,3)), i=is,ie)
    close(1)
 end if
!.... output the final profile (restart file)

 i = ie
 ul = u(i,:,k,n,1)
 vl = u(i,:,k,n,2)
 wl = u(i,:,k,n,3)
 pl = u(i,:,k,n,4)
!     if (ipwbc .eq. 1) pl(1) = pl(2) - (pl(3)-pl(2))/(y(3)-y(2))*y(2)
 open(10,file='pro.dat')
 write(10,"('# ',10(1pe20.13,1x))") one, zero,  real(alpha(i,k,n)), &
      aimag(alpha(i,k,n))

 open(71,file='alpha.pse')
 write(71, "(I4,6E16.7)") (i, x(i), abs(amp(i,k,n)), &
      real(amp(i,k,n)), aimag(amp(i,k,n)), real(alpha(i,k,n)), aimag(alpha(i,k,n)), i=is,ie)
 close(71)

!=============================================================================
!     1    2     3     4     5     6     7     8     9
!     y,  u_r,  u_i,  v_r,  v_i,  w_r,  w_i,  p_r,  p_i
!=============================================================================

 do j = 1, ny
    write(10,"(10(1pe20.13,1x))") y(j), &
         real(ul(j)), aimag(ul(j)), &
         real(vl(j)), aimag(vl(j)), &
         real(wl(j)), aimag(wl(j)), &
         real(pl(j)), aimag(pl(j))
 end do
 close(10)

!.... save a field file

 call post(1)

!.... compute growth rate based on various other quantitities

 call growth(ie-is+1, ny, x(is:ie), y, amp(is:ie,k,n), &
      alpha(is:ie,k,n), u(is:ie,:,k,n,:), opi, nint, yint)

!... compute how close alphas are:

 if ( norm .eq. 3) then
    do i=is, ie
       era(i)  = alpha(i,k,n)+alpha1(i)
       erar(i) = abs(era(i))
    end do

!... find the maximum error:
    err=0.0
    do i=is, ie
       if (erar(i) > err ) then
          err=erar(i)
       end if
    end do
 end if
 print*, "Error in alpha is =========> ", err
 open(14,file='error.dat')
 write(14,'(3E16.7)') (x(i), aimag(era(i)),real(era(i)),i=is,ie)
 close(14)


1001 call error ('solver$','Error reading initial condition$')

 return


!====================================================================================
end subroutine solver
!====================================================================================
