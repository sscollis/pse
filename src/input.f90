!==============================================================================
      subroutine input
!==============================================================================
      use global
      implicit none

      integer :: k, n

      namelist / parm / nx, ny, nz, nt, beta, omega, re, &
                        icurv, imean, ivbc, ipwbc, ipbc, &
                        itype, mkz, mkt, &
                        xmax, xs1, xs2, dx1, ymax, ystr, &
                        tol, ipfix, niter, is, ie, norm, nint, int, sor, &
                        newton, relative, tau

      namelist / nlin / plock, eps, alpi
!==============================================================================

!.... initialize variables

      nx=100; ny=64; nz=1; nt=1; beta=zero; omega=zero
      xmax=10.0; xs1=1.0e6; xs2=12.0; dx1=1.0
      ymax=2.0; ystr=0.0; re=100.0; tol=1.0e-10
      icurv=0; imean=0; ivbc=0; itype=0; mkz=0; mkt=1
      ipbc=0; ipfix=1; ipwbc=1; niter=20; is=1; ie=nx; norm=0; nint=1001
      sor=1.0; int=5; newton = .false.; relative=.true.; tau = zero;

      eps = 1.0e-8; plock = 0; alpi = -1.0e-3

!.... read PARM namelist

      open(10,file='pse.inp',status='old',err=1000)
      read(10,parm,err=1000)

!.... read NLIN namelist

      if (itype.eq.1) then
        read(10,nlin,err=1000)
      end if

      close(10)

!.... generate the mean flow

      if (imean.eq.0) then               ! parallel flow
        call genpar
      else if (imean.eq.1) then          ! boundary layer flow
        call genmf(0)
      else if (imean.eq.2) then          ! Navier-Stokes flow
        call genmean
      end if

      if (is.eq.0 .or. is.lt. 1) is= 1   ! for convenience
      if (ie.eq.0 .or. ie.gt.nx) ie=nx   ! for convenience

!.... generate the wavenumbers

      if (itype.ne.1) then
        nz = 1; nt = 1; mz = 1; mt = 1
        allocate(kz(nz), kt(nt))
        kz(1) = mkz
        kt(1) = mkt
      else
        if (nz.lt.3 .or. nt.lt.3) &
          call error('input$','nz or nt too small for nonlinear analysis$')
        mz = (nz+2)/2
        mt = (nt+2)/2
        allocate(kz(nz), kt(mt))
        kz(1) = 0
        do k = 2, mz
          kz(nz+2-k) = -(k-1)
          kz(k) = k-1
        end do
        do n = 1, mt
          kt(n) = n-1
        end do
      end if

      return
1000  call error('input$','Error reading pse.inp$')
      end subroutine input
