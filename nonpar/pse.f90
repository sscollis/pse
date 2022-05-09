!==============================================================================
      program pse
!
!     Solve the incompressible linear or nonlinear 
!     Parabolized Stability Equations using Chebyshev collocation
!
!     Also includes the option of generating an inflow profile by
!     solving the Orr-Sommerfeld equation locally with and without curvature
!
!     Author:   S. Scott Collis
!
!     Date:     4-29-98
!
!     Copyright:  S. Scott Collis
!                 Department of Mechanical Engineering and Materials Science
!                 Rice University, MS 321
!                 Houston, TX 77005-1892
!                 (713) 527-8101 x3617
!                 collis@rice.edu
!==============================================================================
      use global
      implicit none
!==============================================================================
      cpu = second()

!.... get user input

      call input

!.... setup the mesh and derivative operators

      call setup

!.... form the mean flow and gradients

      call mean

!.... set the curvature

      if (icurv.eq.0) cur = zero     ! no curvature

!.... linear stability theory (negative itype)

      if (itype.eq.-1) then
	call lst                     ! spatial theory
	stop
      else if (itype.eq.-2) then
	call tlst(.false.)           ! temporal theory
	stop
      else if (itype.eq.-3) then
	call adjoint                 ! spatial adjoint theory
	stop
      else if (itype.eq.-4) then
	call tadjoint                ! temporal adjoint theory
	stop
      else if (itype.eq.-5) then
	call tlst(.true.)            ! discrete temporal adjoint theory
	stop
     else if (itype .eq. -6 ) then
        call ddalpha
        stop
      end if

!.... generate the initial condition

      call genini

!.... call the appropriate PSE solver (positive itype)

      if (itype.eq.0) then
	call solver                  ! linear pse
      else if (itype.eq.1) then
	call nsolver                 ! nonlinear pse
      else if (itype.eq.2) then
	call asolver                 ! adjoint linear pse
      else if (itype.eq.3) then
	call dasolver
      else
	call error('pse$','Illegal value of linear$')
      end if

      stop
      end
