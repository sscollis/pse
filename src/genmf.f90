!==============================================================================
      module findmf
!==============================================================================
      integer :: nin
      real, allocatable :: cyin(:)
      real :: yp

      end module
      
!==============================================================================
      subroutine genmf(icode)
!
!     Use Chebyshev interpolation to convert from the BC code mesh to the 
!     current mesh in y.
!
!==============================================================================
      use global
      use findmf
      implicit none

      real, allocatable ::  uin(:), vin(:), win(:), yin(:), w(:,:), &
                            cuin(:), cvin(:), cwin(:), ct(:)

      real, external :: f, zeroin
#ifdef CRAY 
      real, external :: sdot
#else
      real, external :: ddot
#endif

      character*1,parameter :: trans='n'
      integer :: ix, i, j, ii, icode, idum
      real :: rdum, pin, con, arg
!==============================================================================

      if (icode.eq.0) then

!.... open the boundary layer file and read the x grid

        open (4, file='bl_sta.out', status='old', err=110)
        read (4,*) ns, nin, rdum, rdum, re

        allocate(x(ns))
        read (4,*) (x(i), i= 1,ns)

        do i = 1, ns
          read(4,*,end=10) idum
          read(4,*) rdum
          do j = 1, nin
            read(4,*) rdum
          end do
        end do
        nx = ns
        goto 20
10      write(*,"('Premature end of bl_sta.dat file == separation?')")
        nx = i-1
20      write(*,"('Setting Nx = ',i4)") nx
        close(4)

        open (4, file='bl_sta.out', status='old', err=110)
        read (4,*) ns, nin, rdum, rdum, re
        
        allocate( y(ny), cur(nx) )
        allocate( ub(ny,nx),  vb(ny,nx),  wb(ny,nx) )
        allocate( uby(ny,nx), vby(ny,nx), wby(ny,nx) )
        allocate( ubx(ny,nx), vbx(ny,nx), wbx(ny,nx) )

        read (4,*) (x(i), i= 1,ns)
        xmax = x(nx)

      else if (icode.eq.1) then

!.... Chebyshev interpolation of mean BL from spectral BL code nsblsp5
!.... Unit 4 is already open from previous call to genmf with icode=0

      do ix= 1,nx

        read (4,*) nin
        allocate( uin(nin), vin(nin), win(nin), yin(nin), w(nin,nin), &
                  cuin(nin), cvin(nin), cwin(nin), ct(nin), cyin(nin) )

        read (4,*) rdum, rdum, rdum, cur(ix)
        read (4,*) (rdum, yin(j), uin(j), vin(j), win(j), j= 1,nin)

        pin = pi / (nin-1)
        con = 2. / (nin-1)
        do j= 1,nin
          do i= 1,nin
            w(i,j) = con * cos((i-1) * (j-1) * pin)
          end do
        end do
        do i= 1,nin
          w(i,1) = .5 * w(i,1)
          w(1,i) = .5 * w(1,i)
          w(i,nin) = .5 * w(i,nin)
          w(nin,i) = .5 * w(nin,i)
        end do

#ifdef CRAY
        call sgemv ( trans, nin, nin, one, w, nin, yin, 1, zero, cyin, 1)
        call sgemv ( trans, nin, nin, one, w, nin, uin, 1, zero, cuin, 1)
        call sgemv ( trans, nin, nin, one, w, nin, vin, 1, zero, cvin, 1)
        call sgemv ( trans, nin, nin, one, w, nin, win, 1, zero, cwin, 1)
#else
        call dgemv ( trans, nin, nin, one, w, nin, yin, 1, zero, cyin, 1)
        call dgemv ( trans, nin, nin, one, w, nin, uin, 1, zero, cuin, 1)
        call dgemv ( trans, nin, nin, one, w, nin, vin, 1, zero, cvin, 1)
        call dgemv ( trans, nin, nin, one, w, nin, win, 1, zero, cwin, 1)
#endif

        do j= 1,ny
          if (y(j) .lt. yin(nin)) then
            yp = y(j)
            arg = zeroin (-1., 1., f, 0.)
            arg = acos(arg)
            do ii= 1,nin
              ct(ii) = cos((ii-1) * arg)
            end do
#ifdef CRAY
            ub(j,ix) = sdot (nin, cuin, 1, ct, 1)
            vb(j,ix) = sdot (nin, cvin, 1, ct, 1)
            wb(j,ix) = sdot (nin, cwin, 1, ct, 1)
#else
            ub(j,ix) = ddot (nin, cuin, 1, ct, 1)
            vb(j,ix) = ddot (nin, cvin, 1, ct, 1)
            wb(j,ix) = ddot (nin, cwin, 1, ct, 1)
#endif
          else
!           ub(j,ix) = uin(nin)
!           vb(j,ix) = vin(nin)
!           wb(j,ix) = win(nin)
            ub(j,ix) = uin(nin) + (uin(nin)-uin(nin-1))/(yin(nin)-yin(nin-1)) &
                     * (y(j) - yin(nin))
            vb(j,ix) = vin(nin) + (vin(nin)-vin(nin-1))/(yin(nin)-yin(nin-1)) &
                     * (y(j) - yin(nin))
            wb(j,ix) = win(nin) + (win(nin)-win(nin-1))/(yin(nin)-yin(nin-1)) &
                     * (y(j) - yin(nin))
          endif
        end do

        deallocate( uin, vin, win, yin, w, cuin, cvin, cwin, ct, cyin )

      end do
      close(4)

      else
        call error ('genmf$','Illegal icode$')
      end if

      return
110   call error ('genmf$','Error opening bl_sta.out$')
      end

!==============================================================================
      real function f(x)

      use findmf
!==============================================================================
      arg = acos(x)
      f = -yp
      do i= 1,nin
         f = f + cyin(i) * cos((i-1) * arg)
      end do

      return
      end
