!==============================================================================
      real function zeroin(ax,bx,f,tol)
!==============================================================================
      real :: ax,bx,tol
!
!  A zero of the function  f(x)  is computed in the interval ax,bx
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result ( .ge. 0.0)
!
!  output..
!
!  zeroin abcissa approximating a zero of  f  in the interval ax,bx
!
!      It is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!      This function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
      real :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      real, external :: f

!.... compute eps, the relative machine precision

      eps = 1.0
 10   eps = eps/2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) go to 10

!.... initialization

      a = ax
      b = bx
      fa = f(a)
      fb = f(b)

!.... begin step

 20   c = a
      fc = fa
      d = b - a
      e = d
 30   if (abs(fc) .ge. abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa

!.... convergence test
     
 40   tol1 = 2.0*eps*abs(b) + 0.5*tol
      xm = .5*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0) go to 90

!.... is bisection necessary

      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70

!.... is quadratic interpolation possible

      if (a .ne. c) go to 50

!.... linear interpolation

      s = fb/fa
      p = 2.0*xm*s
      q = 1.0 - s
      go to 60

!.... inverse quadratic interpolation
     
 50   q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)

!.... adjust signs

 60   if (p .gt. 0.0) q = -q
      p = abs(p)

!.... is interpolation acceptable

      if ((2.0*p) .ge. (3.0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5*e*q)) go to 70
      e = d
      d = p/q
      go to 80

!.... bisection

 70   d = xm
      e = d
     
!.... complete step

 80   a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/abs(fc))) .gt. 0.0) go to 20
      go to 30
     
!.... done
     
 90   zeroin = b
      return
      end function zeroin
