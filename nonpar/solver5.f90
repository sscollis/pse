!======================================================================================
!
! This version solver3 is different from solver1 due to the organization
! and d(omega)/d(alpha) terms. 
!
!
!======================================================================================

module ualpha 

!======================================================================================
!
! in this module i will maintain all the necessary variables which are 
! important for calculating du/d\alpha terms 
! 
  
  complex, allocatable           ::  ra(:)
  complex, allocatable           ::  ua(:,:,:)
  complex, allocatable           ::  uax(:), vax(:), wax(:), pax(:) 
  complex, allocatable           ::  aduax(:), advax(:), adwax(:), adpax(:)
  complex, allocatable           ::  ual(:), val(:), wal(:), pal(:) 
  complex, allocatable           ::  adual(:), adval(:), adwal(:), adpal(:)
  complex                        ::  dwda, da

!======================================================================================

end module ualpha

!======================================================================================

!======================================================================================

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
!
!======================================================================================
!
!     .--,       .--,
!    ( (  \.---./  ) )
!     '.__/o   o\__.'               8-10-99    1) Calculate u,alpha
!        {=  ^  =}
!         >  -  <
!        /       \
!       //       \\
!      //|   .   |\\
!      "'\       /'"_.-~^`'-.
!         \  _  /--'         `
!   jgs ___)( )(___
!      (((__) (__)))
! 
!======================================================================================

  use int_str
  use global
  use ualpha
  
  
  implicit none
  
  complex, allocatable     :: dotp(:,:), dotpx(:,:)
  complex, allocatable     :: Gb(:,:,:), Ab(:,:,:), A0(:,:), r(:)
  complex, allocatable     :: adu(:,:,:), adua(:,:,:)
  complex                  :: dalpha(nz,nt), alphal

  complex                  :: ul(ny), vl(ny), wl(ny), pl(ny)
  complex                  :: ux(ny), vx(ny), wx(ny), px(ny)

  complex                  :: adul(ny), advl(ny), adwl(ny), adpl(ny), alphad(nx)
  complex                  :: adux(ny), advx(ny), adwx(ny), adpx(ny)
  complex                  :: uh(nx,ny,ndof), amph(nx), alphah(nx), uha(nx,ny,ndof)



  complex                  :: uy(ny,ndof), uyy(ny,ndof), h(nx,10), era(nx)
  complex                  :: kex, umax(nx), dumax, jprod(nx)
  complex                  :: temp(20), al, ala
  complex                  :: gu2, gumax, gke, gj 
  
  real, allocatable        :: G(:,:,:), A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), &
                              E(:,:,:), Ex(:,:,:), xl(:), yl(:), erar(:)
   
  real                     :: ke(nx), ar, ai, ar1, ai1, rdum, dxf, ampr, ampi 
  real                     :: c1(ny), c2(ny), err, tmp(10), ampl
  real                     :: rei, df1, df2, df3, df4, df5, xc, xcm, xcmm
  real                     :: ur, ui, vr, vi, wr, wi, prr, pri
  
!.... local variables

  real                     :: ubl(ny),  vbl(ny),  wbl(ny)
  real                     :: ubxl(ny), vbxl(ny), wbxl(ny)
  real                     :: ubyl(ny), vbyl(ny), wbyl(ny)
  
  integer                  :: i, j, k, n, iter, il, aflag, fnorm, idum, idof, nxl ,nyl
  integer                  :: l, num, ldof, l0, m, mdof, m0, info, job=0, icount=0
  integer, allocatable     :: ipvt(:)

!.... wall BCs

  real                     :: bloc = 100.0, bcon = 0.0, bamp = -1.0e-1
  real                     :: coe, bl

  character                :: cdum
  character(20)            :: fname

  namelist /bc/ bloc, bcon, bamp
!=======================================================================================
!  
!  Explanation of some variables
!  h - scratch variables, usually hold the J terms.
!
!  run is used when I want to do the j norm!
!
!=======================================================================================
!.... allocate variables
 

  rei = one / re
  
  allocate( u(nx,ny,nz,nt,ndof) )
  u = czero
  
  allocate( G(ny,ndof,ndof), A(ny,ndof,ndof), B(ny,ndof,ndof), &
            C(ny,ndof,ndof), D(ny,ndof,ndof), E(ny,ndof,ndof), &
            Ex(ny,ndof,ndof) )
  
  allocate( Gb(ny,ndof,ndof), Ab(ny,ndof,ndof), A0(ny*ndof,ny*ndof) )
  allocate( r(ny*ndof), ipvt(ny*ndof) )

!.... allocate alpha variables 
  
  allocate(ra(ny*ndof), ua(nx,ny,ndof))
  allocate(aduax(ny), advax(ny), adwax(ny), adpax(ny))
  allocate(uax(ny), vax(ny), wax(ny), pax(ny))
  allocate(ual(ny), val(ny), wal(ny), pal(ny)) 
  allocate(adual(ny), adval(ny), adwal(ny), adpal(ny))

  
!.... allocate J innerproduct terms 

  allocate(dotp(nx,ny),  dotpx(nx,ny), erar(nx))
  
!.... initialize some of the variables

  h  =   czero; amph  = czero
  dotp = czero; dotpx = czero; 

!.... open various debugging files


!======================================================================================
!.... consider only one mode for linear analysis

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
  pl = u(i,:,k,n,4)

  ke =    zero
  kex =   zero
  umax =  zero
  dumax = zero

!======================================================================================
!.... read the hat functions from the adjoint and pse fields to compute j norm 

  if ( run .eq. 1 ) then 
     open(10,file='afield.hat',form='unformatted')
     read(10) nxl, nyl 
     allocate(adu(nxl,nyl,4),xl(nxl),yl(nyl))
     read(10) (((adu(i,j,num),j=1,nyl),i=is,ie),num=1,4)
     read(10) (xl(i),i=is,ie)
     read(10) (yl(j),j=1,nyl) 
     close(10)
     
     open(10,file='afield.al',form='unformatted')
     read(10) nxl, nyl 
     allocate(adua(nxl,nyl,4))
     read(10) (((adua(i,j,num),j=1,nyl),i=is,ie),num=1,4)
     read(10) (xl(i),i=is,ie)
     read(10) (yl(j),j=1,nyl) 
     close(10)

!.... make sure that xl and yl are the same as x and y
  
     if (ie-is+1 .ne. nxl) then 
        print*,'Check x mesh...'
        print*,ie,nxl
     end if
     if (ny .ne. nyl) then 
        print*,'Check y mesh...'
        print*,ny,nyl
        stop 666
     end if
     
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

!.... Read pse alha

     open(10,file='alpha.pse')
     do i=is,ie
        read(10,*) idum, rdum, rdum, ar1, ai1, ar, ai
        alphah(i)= cmplx(ar,ai)
        alpha(i,n,k)=alphah(i)
     end do
     close(10)

  end if

!.... read PSE field file
  
  if (ipfix.eq.2) then
     open (10, file='field.pse', form='unformatted')
     read (10) nxl, nyl
     if (nxl.ne.(ie-is+1) .or. nyl.ne.ny) &
	  call error('solver$','Illegal field.pse file$')
     read (10) (((uh(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
     read (10) (rdum, i= is,ie)
     read (10) (rdum, j= 1,ny)
     read (10) (amph(i), i = is, ie)
     read (10) (((uha(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
     read (10) (alphah(i), i= is, ie)
     close (10)
     alpha(:,k,n) = alphah(:) 

     u(:,:,k,n,:) = uh

     ul = u(i,:,k,n,1)
     vl = u(i,:,k,n,2)
     wl = u(i,:,k,n,3)
     pl = u(i,:,k,n,4)
  
  end if

  
!=======================================================================================

  jprod = czero
  if ( norm .eq. 3) then
     i=is
     dotp(i,:) = czero
     do j=1,ny
        temp(9)   =   adu(i,j,1) * ul(j) + adu(i,j,2) * vl(j) 
        dotp(i,j) =   temp(9) * ub(j,i) + adu(i,j,1) * pl(j) + adu(i,j,4) * &
                      ul(j) - 2 * rei * iota * alpha(i,k,n) * temp(9)
        jprod(i)  =   jprod(i) + opi(j)*dotp(i,j) 
     end do
  end if

!=======================================================================================
!.... different versions of the kinetic energy, depending on the norm

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

!.... modified slightly to agree in format with the new version of lns solver

  open(9,file='pse.dat')
  write(9,*) ie-is+1
  write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), abs(kex/ke(i)),&
                              amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
  call flush(9)

!.... Initialize u alpha

  ua=czero

!=========================================================================================
!                                   main marching loop
!=========================================================================================

  do i = is+1, ie

!=========================================================================================
!.... Calculate some useful adjoint quantities

     df1 = one/(x(i) - x(i-1))
     if ( i .ne. ie ) then
        df2 = one/(x(i+1) - x(i))
     else
        df2 = one/(x(i) - x(i-1))
     end if

     if ( run .eq. 1 .and. norm .eq. 3 ) then

        adul  = adu(i,:,1)
        advl  = adu(i,:,2)
        adwl  = adu(i,:,3)
        adpl  = adu(i,:,4)

        adual  = adua(i,:,1)
        adval  = adua(i,:,2)
        adwal  = adua(i,:,3)
        adpal  = adua(i,:,4)

        if ( i.ne. ie) then 
           adux  = df2 * (adu(i+1,:,1) - adu(i,:,1))
           advx  = df2 * (adu(i+1,:,2) - adu(i,:,2))
           adwx  = df2 * (adu(i+1,:,3) - adu(i,:,3))
           adpx  = df2 * (adu(i+1,:,4) - adu(i,:,4))
        else
           
           adux  = df2 * (adu(i,:,1) - adu(i-1,:,1))
           advx  = df2 * (adu(i,:,2) - adu(i-1,:,2))
           adwx  = df2 * (adu(i,:,3) - adu(i-1,:,3))
           adpx  = df2 * (adu(i,:,4) - adu(i-1,:,4))
        end if

        aduax = -adua(i,:,1) * df2
        advax = -adua(i,:,2) * df2
        adwax = -adua(i,:,3) * df2
        adpax = -adua(i,:,4) * df2

     end if

     cpu2 = second()
     write(*,"(/,'  i       x(i)          cpu        total cpu')")
     write(*,"(i4,3(1x,1pe13.6))") i, x(i), cpu2-cpu, cpu2
     cpu = cpu2


!.... initialize alpha


     if (run .ne. 1 .and. ipfix .ne. 2 ) then
        alpha(i,k,n) = alpha(i-1,k,n)
     end if


     c1(:) = one/(one + cur(i)*y(:))
     c2(:) = cur(i) * c1(:)
     
     ul = u(i-1,:,k,n,1)
     vl = u(i-1,:,k,n,2)
     wl = u(i-1,:,k,n,3)
     pl = u(i-1,:,k,n,4)
     
!.... compute the mode amplitude 
     
     amp(i,k,n) = czero
     do il = is, i-1
        amp(i,k,n) = amp(i,k,n) + pt5 * (alpha(il,k,n) + &
                     alpha(il+1,k,n)) * (x(il+1)-x(il))
     end do
     amp(i,k,n) = exp( iota * amp(i,k,n) )
     
!=========================================================================================
!.... print statements

     if (norm.eq.0) then
        ke(i) = zero
        do j = 1, ny
           ke(i) = ke(i) + opi(j)*(abs(ul(j))**2 + & 
                   abs(vl(j))**2 + abs(wl(j))**2)
        end do
        write(*,21)
        write(*,20) 0, alpha(i,k,n), abs(amp(i,k,n)), ke(i), zero
 
     else if (norm.eq.1) then
        call findumax( ny, y, ul, umax(i), nint, yint )
        write(*,21)
        write(*,20) 0, alpha(i,k,n), abs(amp(i,k,n)), abs(umax(i)), zero

     else if (norm.eq.2) then
        ke(i) = zero
        do j = 1, ny
           ke(i) = ke(i) + opi(j)*(abs(ul(j))**2)
        end do
        write(*,21)
        write(*,20) 0, alpha(i,k,n), abs(amp(i,k,n)), abs(ke(i)), zero	  
     
     else if ( norm .eq. 3) then 
        write(*,22)
        write(*,20) 0, alpha(i,k,n), abs(amp(i,k,n)), abs(jprod(i)), abs(h(i,3)), zero
     
     end if
     
!=========================================================================================
!                                  loop on alpha
!=========================================================================================
     do iter = 1, niter 
	
!.... compute the streamwise derivative of alpha

	dalpha(k,n) = ( alpha(i,k,n) - alpha(i-1,k,n) ) / ( x(i) - x(i-1) )
        write(15,"(I4,2E16.7)") i, real(alpha(i,k,n)), aimag(alpha(i,k,n))

!.... Bertollotti sets this to zero (no advantage)

!	dalpha(k,n) = zero

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
        
	D(:,1,1) = c1(:) * ubxl + c2(:) * vbl + rei * c2(:)**2   !+rei ? -rei? 
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
        

!.... fix the streamwise pressure gradient. I save A and use Ab
!.... if ipfix = 2 then calculate the pressure derivative, but a bit later

        Ab = A
	if (ipfix.eq.1 )  then
           Ab(:,:,4) = zero
        else if ( ipfix.eq.2 ) then 
           Ab(:,:,4) = zero  
        end if

	Ab = Ab - two * iota * alphal * Ex

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

!.... subtracting the pressure gradient from the right hand side

	if (ipfix.eq.2) then
           ldof = 1
           do l = 1, ny
              l0 = (l-1)*ndof
              r(l0+ldof) = r(l0+ldof) - c1(l) * &
                   (uh(i,l,4) - uh(i-1,l,4)) / &
                   (x(i)-x(i-1))
           end do
        end if
                        

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
	call ZGEFA(A0,ny*ndof,ny*ndof,ipvt,info)
	if (info.ne.0) call error('solver$','Singular matrix$')
	call ZGESL(A0,ny*ndof,ny*ndof,ipvt,r,job)
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
        px = (pl - u(i-1,:,k,n,4))/(x(i)-x(i-1))

	ke(i) = zero
	kex   = zero
        
!=========================================================================================

!.... Build The RHS for the alpha equation

        ra=czero

!.... Adding the extra forcing sources

!.... -iota A u -2 * alpha Ex u +2 * iota * Ex ux + iota * Ex / (x(i)-x(i-1))

!.... reusing matrix Ab. Do not fix the pressure gradient here!

        Ab=czero
        Ab =  iota * A + two  * alphal * Ex - iota * Ex / (x(i)-x(i-1))
        
        do ldof =1, ndof
           do l =1, ny
              l0 = (l-1)*ndof
              do mdof =1,ndof
                 ra(l0+ldof) = ra(l0+ldof) + &
                      two * iota * Ex(l,ldof,mdof)* &
                      (u(i,l,k,n,mdof)-u(i-1,l,k,n,mdof))/(x(i)-x(i-1)) - &
                      Ab(l,ldof,mdof) *  u(i,l,k,n,mdof)
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

        r=ra

!.... What about trying the continuity at the outer boundary? 
!.... Solve the system

#ifdef CRAY
	call CGESL(A0,ny*ndof,ny*ndof,ipvt,ra,0)
#else
	call ZGESL(A0,ny*ndof,ny*ndof,ipvt,ra,job)
#endif

!.... Update the solution

	do ldof = 1, ndof
           do l = 1, ny
              l0 = (l-1)*ndof
              ua(i,l,ldof) = ra(l0+ldof)
           end do
	end do

!.... Local variables
        
	ual = ua(i,:,1)
	val = ua(i,:,2)
        wal = ua(i,:,3)
        pal = ua(i,:,4)

!.... compute some streamwise derivatives

	uax = ua(i,:,1)/(x(i)-x(i-1))
	vax = ua(i,:,2)/(x(i)-x(i-1))
        wax = ua(i,:,3)/(x(i)-x(i-1))
        pax = ua(i,:,4)/(x(i)-x(i-1))

!.... output some ua profiles

!!$        if ( mod(i,10) .eq. 0) then
!!$           fname = 'Prof/ddua.'
!!$           fname = trim(fname) // trim(i2c(i))
!!$           open(10,file=fname)
!!$           write(10,"('# ',10(1pe20.13,1x))") one, zero, &
!!$                            real(alphal), aimag(alphal), &
!!$                           omega, real(dwda), aimag(dwda)
!!$           do j = 1, ny
!!$              write(10,"(10(1pe20.13,1x))") y(j), &
!!$                   real(ua(i,j,1)), aimag(ua(i,j,1)), &
!!$                   real(ua(i,j,2)), aimag(ua(i,j,2)), &
!!$                   real(ua(i,j,3)), aimag(ua(i,j,3)), &
!!$                   real(ua(i,j,4)), aimag(ua(i,j,4))
!!$           end do
!!$        end if

!=========================================================================================
        
	do j = 1, ny
           ke(i) = ke(i) + opi(j)*( abs(ul(j))**2 + abs(vl(j))**2 + &
                abs(wl(j))**2 )
           kex = kex + opi(j)*( conjg(ul(j))*ux(j) + conjg(vl(j))*vx(j) + &
                conjg(wl(j))*wx(j) )
	end do

!=========================================================================================       
	if (norm.eq.0) then
!=========================================================================================           
!.... compute the correction for alpha (following Herbert AIAA-93-3053)
        if (newton .eq. 0) then  
          gke = czero
          do j=1,ny
            gke = gke + opi(j)*(abs(ul(j))**2 + abs(vl(j))**2 + abs(wl(j))**2)
          end do
        else if (newton .eq. 1) then
          gke = czero
          do j=1,ny
              gke =  gke + opi(j) * ( ux(j) * conjg(ua(i,j,1))+ uax(j) * conjg(ul(j)) + &
                                      vx(j) * conjg(ua(i,j,2))+ vax(j) * conjg(vl(j)) + &
                                      wx(j) * conjg(ua(i,j,3))+ wax(j) * conjg(wl(j)) ) 
          end do
        end if

        alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
             sor * (alpha(i,k,n) -  kex / gke )

!=========================================================================================
        else if (norm.eq.1) then
!=========================================================================================
!.... compute the correction for alpha (following Bertolotti)


           call findumax( ny, y, ul, umax(i), nint, yint )
           dumax = (umax(i)-umax(i-1))/(x(i)-x(i-1))
           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                sor * (alpha(i,k,n) - iota * dumax / umax(i))


!=========================================================================================
        else if (norm.eq.2) then
!=========================================================================================

           ke(i) = zero
           kex   = zero
           do j = 1, ny
              ke(i) = ke(i) + opi(j)*( abs(ul(j))**2  )
              kex = kex + opi(j)*( conjg(ul(j))*ux(j) )
           end do

           gu2 = czero
           dotp   = czero
           do j=1,ny
              dotp(i,j) =  ux(j) * conjg(ua(i,j,1))+ uax(j) * conjg(ul(j))
              gu2       = gu2 + opi(j) * dotp(i,j)   
           end do

           gu2= one/gu2

           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                           sor*(alpha(i,k,n) - gu2 * kex ) 


!==========================================================================================         
        else if ( norm .eq. 3) then
!========================================================================================== 
!.... compute alpha making \hat{J,x}\over \hat{J} -> 0 
           
           
           h(i,1)=czero 
           h(i,2)=czero
           h(i,3)=czero
           h(i,4)=czero
           h(i,5)=czero
           h(i,6)=czero
           h(i,7)=czero
           h(i,8)=czero

           dotp(i,:) = czero

!.... Computing J

           do j=1,ny
              temp(9)   =   adul(j) * ul(j) + advl(j) * vl(j) 
              dotp(i,j) =   temp(9) * ub(j,i) + adul(j) * pl(j) + adpl(j) * &
                            ul(j) - two * rei * iota * alphal * temp(9)
              h(i,1)    =   h(i,1) + opi(j) * dotp(i,j)
           end do

!.... making streamwise derivatives
           
           temp=czero
           
           do j=1,ny
              temp(1) = ux(j)   * adul(j) + ul(j) * adux(j)  + vx(j)  * advl(j) + vl(j) * advx(j)
              temp(2) = ux(j)   * adpl(j) + ul(j) * adpx(j)  + px(j)  * adul(j) + pl(j) * adux(j)
              temp(3) = adul(j) * ul(j) + advl(j) * vl(j)

              h(i,2)  = h(i,2) + opi(j) * ( ubx(j,i) * temp(3) + ub(j,i) * temp(1) + &
                        temp(2) - two * iota * dalpha(k,n) * rei * temp(3) - & 
                        two * iota * alphal * rei * temp(1) )
           end do


!.... Doing J,(x alpha)

           temp = czero

!.... Think about derivatives

           do j=1,ny

              temp(1) = ux(j)  * adul(j) + ul(j) * adux(j)  + vx(j)  * advl(j) + vl(j) * advx(j)
              temp(2) = uax(j) * adul(j) + ual(j) * adux(j)! + ul(j) * aduax(j) + ux(j) * adual(j) 
              temp(3) = vax(j) * advl(j) + val(j) * advx(j)! + vl(j) * advax(j) + vx(j) * adval(j)
              temp(4) = uax(j) * adpl(j) + ual(j) * adpx(j)! + ul(j) * adpax(j) + ux(j) * adpal(j)
              temp(5) = pax(j) * adul(j) + pal(j) * adux(j)! + pl(j) * aduax(j) + px(j) * adual(j) 
              temp(6) = ual(j) * adul(j) + val(j) * advl(j)! + vl(j) * adval(j) + ul(j) * adual(j)
              temp(9) = ul(j)  * adul(j) + vl(j) * advl(j) 

              temp(7) = temp(2) + temp(3) 
              temp(8) = temp(4) + temp(5) 

              dotp(i,j)  = ubx(j,i) * temp(6) + ub(j,i) * temp(7) + temp(8) - &
                           two * iota * dalpha(k,n) * rei * temp(6) - & 
                           two * iota * rei * temp(1) - two * iota * alphal * rei * temp(7) &
                           - two * iota * rei * temp(9)/(x(i-1)-x(i))
              
              h(i,3)     = h(i,3) + opi(j) * dotp(i,j) 
           end do
            
           alpha(i,k,n) = (one-sor)*alpha(i,k,n) + &
                sor*(alpha(i,k,n) - h(i,2)/h(i,3))

           da = -sor * h(i,2)/h(i,3)

           write(16,"(I4,4E16.7)") i, abs(h(i,1)), abs(h(i,2)), abs(h(i,3))
           call flush(16)

!.... update adjoint based on the change in alpha, using derivatives adua
           
           if ( abs(da) .le. 1.0e-7 .and. .false. ) then
              adu(i,:,1) = adu(i,:,1) + adual * da
              adu(i,:,2) = adu(i,:,2) + adval * da 
              adu(i,:,3) = adu(i,:,3) + adwal * da
              adu(i,:,4) = adu(i,:,4) + adpal * da
           end if
        
!========================================================================================== 
        end if
!========================================================================================== 
        
!.... update the mode amplitude

	amp(i,k,n) = czero
	do il = is, i-1
           amp(i,k,n) = amp(i,k,n) + &
                pt5 * (alpha(il,k,n) + alpha(il+1,k,n)) * &
                (x(il+1)-x(il))
	end do
	amp(i,k,n) = exp( iota * amp(i,k,n) )
        
!========================================================================================== 
!.... convergence check
!.... for the Jx norm it does not matter 
!.... since we do only one iteration 

!====================================================================================
100 format ('->',i2,1x,6(1pe13.6,1x))
        
	if (norm.eq.0) then
           write(*,100) iter, alpha(i,k,n), & 
                        abs(amp(i,k,n)), ke(i), &
                        abs(kex/ke(i))/abs(alpha(i,k,n))
           if (abs(kex/ke(i)) .lt. tol) exit

	else if (norm.eq.1) then
           write(*,100) iter, alpha(i,k,n), & 
                        abs(amp(i,k,n)), abs(umax(i)), &
                        abs(dumax/umax(i))/abs(alpha(i,k,n))
	  if (abs(dumax/umax(i)) .lt. tol) exit

	else if (norm.eq.2) then
           write(*,100) iter, alpha(i,k,n), abs(amp(i,k,n)), abs(kex), abs(gu2 * kex)  
           if ( abs(gu2 * kex) .lt. tol) exit


        else if ( norm .eq. 3)  then
           write(*,"('->',i2,1x,7(1pe13.6,1x))" ) iter, real(alpha(i,k,n)), aimag(alpha(i,k,n)), &
                            abs(amp(i,k,n)), abs(h(i,1)),  abs(h(i,3)),  abs(h(i,2)/h(i,3))
           temp(4)=h(i,2)/h(i,3) 
           if (abs(temp(4)) .lt. tol ) exit
     
        end if
        
!=====================================================================================
     end do! loop on iter
!=====================================================================================

     
     write(15,"(I4,2E16.7)") i, real(alpha(i,k,n)), aimag(alpha(i,k,n))
     

!=====================================================================================
!   1      2         3      4         5          6       7        8
!   x,  alpha_r,  alpha_i,  ke,  (dke/dx)/ke,  amp_r,  amp_i,  ke*|amp|^2
!=====================================================================================


     write(9,"(8(1pe20.13,1x))") x(i), alpha(i,k,n), ke(i), &
          abs(kex/ke(i)), amp(i,k,n), ke(i)*abs(amp(i,k,n))**2
     call flush(9)
    

!=====================================================================================
  end do! loop on i
!=====================================================================================

  open  (10, file='field.pse', form='unformatted')
  write (10) ie-is+1, ny
  write (10) (((u(i,j,k,n,idof), j= 1,ny), i= is,ie), idof= 1,4)
  write (10) (x(i), i= is,ie)
  write (10) (y(j), j= 1,ny)
  write (10) (amp(i,k,n), i = is, ie)
  write (10) (((ua(i,j,idof), j= 1,ny), i= is,ie), idof= 1,4)
  write (10) (alpha(i,k,n), i= is,ie)
  close (10)

!.... output the final profile (restart file)
     
     i = ie
     ul = u(i,:,k,n,1)
     vl = u(i,:,k,n,2)
     wl = u(i,:,k,n,3)
     pl = u(i,:,k,n,4)

     open(10,file='pro.dat')
     write(10,"('# ',10(1pe20.13,1x))") one, zero,  real(alpha(i,k,n)), &
          aimag(alpha(i,k,n))
     
     open(71,file='alpha.pse')
     write(71,"(i5,6e22.13)") (i, x(i), abs(amp(i,k,n)), real(amp(i,k,n)), &
          aimag(amp(i,k,n)), real(alpha(i,k,n)), aimag(alpha(i,k,n)), i=is,ie)
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
      alpha(is:ie,k,n), u(is:ie,:,k,n,:), opi, nint, yint, opy, itype, ystr, re)



!====================================================================================
!... compute how close alphas are: 

 if ( norm .eq. 3) then 
    do i=is, ie
       era(i)  = alpha(i,k,n)-alphah(i)
       erar(i) = abs(era(i))
    end do

!... find the maximum error:

    err=0.0
    do i=is, ie
       if (erar(i) > err ) then 
          err=  erar(i)
       end if
    end do
 end if
 print*, "Error in alpha is =========> ", err
 open(14,file='error.dat')
 write(14,'(3E16.7)') (x(i), aimag(era(i)),real(era(i)),i=is,ie)
 close(14)
 
  
!====================================================================================
!.... formats

19 format(i2,1x,10(1pe20.13,1x))
20 format('->',i2,1x,6(1pe13.6,1x))
21 format(/,' iter','    alpha_r','      alpha_i','        |amp|', &
        &'          ke','           dalpha')
22 format(/,' iter','    alpha_r','      alpha_i','        |amp|', &
        & '          J','           J^,x/J^')

 return


!====================================================================================
end subroutine solver
!====================================================================================





