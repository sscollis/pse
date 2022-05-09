!==============================================================================
module global
!==============================================================================
!
!  Parameters, grid, and boundary condition data
!
!     switches: 
!
!     imean:  0= U0 from "profile.dat" (parallel flow)
!             1= U0 from "bl_sta.out"
!             2= U0 from "field.mean"
!
!     icurv:  0= curvature terms off
!             1= curvature terms on (if ipar .ne. 1)
!
!     ivbc:  0= v-inf homog Dirichlet
!            1= v-inf homog Neumann
!            2= continuity enforced @ ymax
!
!==============================================================================
      real, parameter    :: zero   = 0.0
      real, parameter    :: pt25   = 0.25
      real, parameter    :: pt33   = 0.333333333333333333
      real, parameter    :: pt5    = 0.5
      real, parameter    :: one    = 1.0
      real, parameter    :: onept5 = 1.5
      real, parameter    :: two    = 2.0
      real, parameter    :: three  = 3.0
      real, parameter    :: pi     = 3.14159265358979323846
      complex, parameter :: czero  = (0.0,0.0)
      complex, parameter :: cone   = (1.0,0.0)
      complex, parameter :: iota   = (0.0,1.0)
      integer, parameter :: ndof   = 4

      integer :: nx, ny, nz, nt, ns, nss, mkz, mkt, mt, mz, niter, is, ie, norm

      real :: xmax, xs1, xs2, dx1, ymax, ystr, dstar

      real :: re, beta, omega, tol
      logical :: newton, relative

      integer :: ipar, icurv, imean, ivbc, ipbc, itype, ipfix, ipwbc

      complex, allocatable :: alpha(:,:,:), amp(:,:,:)
      complex, allocatable :: ubci(:,:,:), vbci(:,:,:), wbci(:,:,:), &
                              pbci(:,:,:)

      complex, allocatable :: u(:,:,:,:,:)

      real, allocatable :: ub(:,:), ubx(:,:), uby(:,:), &
                           vb(:,:), vbx(:,:), vby(:,:), &
                           wb(:,:), wbx(:,:), wby(:,:)                  

      integer :: nint, int
      real, allocatable :: opy(:,:), opyy(:,:), opi(:), yint(:,:)

      real, allocatable :: cppx(:), cpx(:), cx(:), cmx(:), cmmx(:), &
                           c2ppx(:), c2px(:), c2x(:), c2mx(:), c2mmx(:)

      real, allocatable :: x(:), y(:), cur(:)
      integer, allocatable :: kz(:), kt(:)

      real, allocatable :: dydeta(:)

!.... stabilization parameter

      real :: tau

!.... nonlinear PSE parameters

      real :: eps, alpi, sor
      integer :: plock

#ifdef CRAY
      real :: cpu, cpu2
      real, external :: second
#else
      real*4 :: cpu, cpu2
      real*4, external :: second
#endif

end module global
