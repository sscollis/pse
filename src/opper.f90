!==============================================================================
      subroutine opper (n, x, op)
!==============================================================================
      real x(n), op(n,n), c(n)
!==============================================================================
      nm = n - 1
      pin = acos(-1.) / nm

      do i= 1,n
         x(i) = cos((i-1) * pin)
         c(i) = 1.
      end do
      c(1) = 2.
      c(n) = 2.
      pin = .5*pin

      do j= 2,n
         ie = j - 1
         do i= 1,ie
            op(i,j) = -.5*(-1)**(i+j) * c(i)/c(j) / (sin(pin*(i+j-2)) &
                      * sin(pin*(i-j)) )
         end do
      end do

      do j= 1,nm
         is = j + 1
         do i= is,n
            op(i,j) = -op(n-i+1,n-j+1)
         end do
      end do

      do i= 2,nm
         op(i,i) = -.5 * x(i) / sin(2.*pin*(i-1)) **2
      end do
      op(1,1) = (2. * nm **2 + 1.) / 6.
      op(n,n) = -op(1,1)

      return
      end subroutine opper

!==============================================================================
      subroutine opperfd4 (nx, x, op, opp)
!==============================================================================
      integer :: nx
      real    :: x(nx), op(nx,nx), opp(nx,nx)

      integer :: i, im, ip
      real    :: dm, dmm, dp, dpp
      real    :: h_1, h_2, h_3, gam, del, alp, bet
      real    :: wv1(nx)

      real    :: cppx(nx), cpx(nx), cx(nx), cmx(nx), cmmx(nx)
      real    :: c2ppx(nx), c2px(nx), c2x(nx), c2mx(nx), c2mmx(nx)
!==============================================================================

      do i= 1,nx-1
         wv1(i) = x(i+1) - x(i)
      end do

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

      do i= 2,nx-1
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

         if (i .eq. 2) then
           cpx(i) = cpx(i) + cmmx(i)
           cx(i) = cx(i) - 3.*cmmx(i)
           cmx(i) = cmx(i) + 3.*cmmx(i)
           cmmx(i) = 0.
           c2px(i) = c2px(i) + c2mmx(i)
           c2x(i) = c2x(i) - 3.*c2mmx(i)
           c2mx(i) = c2mx(i) + 3.*c2mmx(i)
           c2mmx(i) = 0.
         endif
         if (i .eq. nx-1) then
           cmx(i) = cmx(i) + cppx(i)
           cx(i) = cx(i) - 3.*cppx(i)
           cpx(i) = cpx(i) + 3.*cppx(i)
           cppx(i) = 0.
           c2mx(i) = c2mx(i) + c2ppx(i)
           c2x(i) = c2x(i) - 3.*c2ppx(i)
           c2px(i) = c2px(i) + 3.*c2ppx(i)
           c2ppx(i) = 0.
         endif

      end do

!.... Otto's ends

      i=2
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

      i=nx-1
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

!.... Now build the matrix opperators
      
      op = 0.0
      opp = 0.0

      i = 1
      op(i,i)    = cx(i)
      op(i,i+1)  = cpx(i)
      op(i,i+2)  = cppx(i)
      opp(i,i)   = c2x(i)
      opp(i,i+1) = c2px(i)
      opp(i,i+2) = c2ppx(i)

      i = 2
      op(i,i-1)  = cmx(i)
      op(i,i)    = cx(i)
      op(i,i+1)  = cpx(i)
      op(i,i+2)  = cppx(i)
      opp(i,i-1) = c2mx(i)
      opp(i,i)   = c2x(i)
      opp(i,i+1) = c2px(i)
      opp(i,i+2) = c2ppx(i)

      do i = 3, nx-2
        op(i,i-2)  = cmmx(i)
        op(i,i-1)  = cmx(i)
        op(i,i)    = cx(i)
        op(i,i+1)  = cpx(i)
        op(i,i+2)  = cppx(i)
        opp(i,i-2) = c2mmx(i)
        opp(i,i-1) = c2mx(i)
        opp(i,i)   = c2x(i)
        opp(i,i+1) = c2px(i)
        opp(i,i+2) = c2ppx(i)
      end do

      i = nx-1
      op(i,i-2)  = cmmx(i)
      op(i,i-1)  = cmx(i)
      op(i,i)    = cx(i)
      op(i,i+1)  = cpx(i)
      opp(i,i-2) = c2mmx(i)
      opp(i,i-1) = c2mx(i)
      opp(i,i)   = c2x(i)
      opp(i,i+1) = c2px(i)

      i = nx
      op(i,i-2)  = cmmx(i)
      op(i,i-1)  = cmx(i)
      op(i,i)    = cx(i)
      opp(i,i-2) = c2mmx(i)
      opp(i,i-1) = c2mx(i)
      opp(i,i)   = c2x(i)

      return
      end subroutine opperfd4

!==============================================================================
      subroutine opperfd2 (nx, x, op, opp)
!==============================================================================
      integer :: nx
      real    :: x(nx), op(nx,nx), opp(nx,nx)

      integer :: i, im, ip
      real    :: dm, dmm, dp, dpp
      real    :: h_1, h_2, h_3, gam, del, alp, bet
      real    :: dx(nx)

      real    :: cpx(nx), cx(nx), cmx(nx)
      real    :: c2px(nx), c2x(nx), c2mx(nx)
!==============================================================================

      do i= 1,nx-1
         dx(i) = x(i+1) - x(i)
      end do

      cmx(1) =  0.0
      cx(1)  = -1.0 / dx(1)
      cpx(1) =  1.0 / dx(1)
     
      c2mx(1) = 0.0
      c2x(1)  = 0.0
      c2px(1) = 0.0

      cmx(nx) = -1.0 / dx(nx-1)
      cx(nx)  =  1.0 / dx(nx-1)
      cpx(nx) =  0.0

      c2mx(nx) = 0.0
      c2x(nx)  = 0.0
      c2px(nx) = 0.0

      do i = 2,nx-1
        alp = dx(i)/dx(i-1)
        dm  = 1.0/(alp*(alp+1.0)*dx(i))
        cmx(i) = -alp**2*dm
        cx(i)  = (alp**2 - 1.0)*dm
        cpx(i) = dm

        dmm = 2.0 / ( dx(i)**2 * (1.0 + alp) )
        c2mx(i) = dmm
        c2x(i)  = -( 1.0 + 1.0/alp ) * dmm
        c2px(i) = dmm / alp
      end do

!.... Now build the matrix opperators
      
      op = 0.0
      opp = 0.0

      i = 1
      op(i,i)    = cx(i)
      op(i,i+1)  = cpx(i)
      opp(i,i)   = c2x(i)
      opp(i,i+1) = c2px(i)

      do i = 2, nx-1
        op(i,i-1)  = cmx(i)
        op(i,i)    = cx(i)
        op(i,i+1)  = cpx(i)
        opp(i,i-1) = c2mx(i)
        opp(i,i)   = c2x(i)
        opp(i,i+1) = c2px(i)
      end do

      i = nx
      op(i,i-1)  = cmx(i)
      op(i,i)    = cx(i)
      opp(i,i-1) = c2mx(i)
      opp(i,i)   = c2x(i)

      return
      end subroutine opperfd2

!==============================================================================
      subroutine chebyinit( ny, y, opi, nint, yint )
!==============================================================================
      implicit none

      integer :: ny, nint, i, j, ii
      real    :: y(ny), opi(ny), yint(nint,ny)
      real    :: con, pin, pi, w1(ny,ny), f(ny)

!.... Prepare for Chebyshev interpolation

      con = 2. / (ny-1)
      pin = acos(-1.) / (ny-1)
      do j= 1,ny
        do i= 1,ny
          w1(i,j) = con * cos((i-1) * (j-1) * pin)
        end do
      end do
      do i= 1,ny
        w1(i,1) = .5 * w1(i,1)
        w1(i,ny) = .5 * w1(i,ny)
        w1(1,i) = .5 * w1 (1,i)
        w1(ny,i) = .5 * w1(ny,i)
      end do
      con = acos(-1.) / (nint-1)
      do i= 1,nint
        do j= 1,ny
          f(j) = cos((j-1) * (i-1) * con)
        end do
        do j= 1,ny
          yint(i,j) = 0.
          do ii= 1,ny
            yint(i,j) = yint(i,j) + f(ii) * w1(ii,j)
          end do
        end do
      end do

!.... prepare for Chebyshev integration

      pi = acos(-1.)
      f(1) = -2.
      f(2) = 0.
      do i= 3,ny
        f(i) = .5 * ((cos(i*pi) - 1.) / i - (cos((i-2) * pi) - 1.) / (i-2))
      end do
      do i= 1,ny
        opi(i) = 0.
        do ii= 1,ny
          opi(i) = opi(i) + f(ii) * w1(ii,i)
        end do
      end do
      call opper (ny, f, w1)
      do i= 1,ny
        f(i) = 0.
        do ii= 1,ny
          f(i) = f(i) + w1(i,ii) * y(ii)
        end do
      end do
      do i= 1,ny
        opi(i) = opi(i) * f(i)
      end do

      return
      end subroutine chebyinit
