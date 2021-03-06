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
        op(i,j) = -.5*(-1) **(i+j) * c(i)/c(j) / (sin(pin*(i+j-2)) &
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
!==============================================================================
end subroutine opper
!==============================================================================
