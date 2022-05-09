!=============================================================================
! Module int_str        transfer integer number to a string
!
! Version 1.0           Yong Chang Aug. 1, 1998
!
!  1. function i2c(int) : return a string with length 12,
!                     usage: trim(i2c(int))
!  2. function i2c_p(int) : return a character pointer, the exact length
!
!  3. function itoc(i)  : return a string with length 11,
!                     usage: trim(itoc(i))
!  4. function dig(i)   : return the digital number of i
!
!  5. function int2str(i, len_i) : return a string with length len_i
!                     usage: int2str(i,dig(i))
!
! Recommend: use function 1. or 3. or 5.
!
!==============================================================================

module int_str
  implicit none
contains

  !note: use--------:                     trim(i2c(i))
   function i2c(int) result(str)
     integer, intent(in) :: int
     character(len=12)   :: str
     
     write(str,"(i12)") int
     str = adjustl(str)
   end function i2c

   function i2c_p(int) result(str_point)
     integer, intent(in) :: int
     character, dimension(:), pointer :: str_point
     character(len=12) :: tmp
     integer:: lt,i
     
     write(tmp,"(i12)") int
     tmp = adjustl(tmp)
     lt = len_trim(tmp)
     allocate (str_point(lt))
     str_point = (/ (tmp(i:i), i = 1, lt) /)
   end function i2c_p

   function itoc(ival) result(string)               ! use :     trim(itoc(i))
   implicit none
   integer, intent(in) :: ival
   character(len=11)   :: string
   integer :: absval, k              ! local storage
   character :: temp*11
   absval = abs(ival)
   temp   = " "
   do k = 11,1,-1
      temp(k:k) = achar(48+mod(absval,10))
      absval = absval / 10
      if(absval == 0) exit
   end do
   if(ival < 0) then
      k = k - 1
      temp(k:k) = "-"
   end if
   string = adjustl(temp(k:))
   end function itoc

  function dig(i) result (d)
    integer , intent(in) :: i
    integer :: d
    integer :: absi
    absi = abs(i)
    if (absi == i) then
       d = 1
    else
       d = 2
    end if
    do
       if (absi < 10) then
          return
       else
          d = d + 1
          absi = absi / 10

       end if
    end do
  end function dig


  function int2str(i,len_i) result(s)         !  use :    int2str(i,dig(i))
    integer , intent(in) :: i,len_i
    character (len=len_i) :: s
    integer :: i_copy,j,k,is
    i_copy = abs(i)
    if ( i_copy == i ) then
       is = 1
    else
       is = 2
       s(1:1) = "-"
    end if
    do j=len_i,is,-1
       k = modulo(i_copy,10)+1
       s(j:j) = "0123456789"(k:k)
       i_copy = i_copy / 10
    end do
    return
  end function int2str
end module int_str


!!$program t
!!$  use int_str
!!$  integer :: i
!!$  write(*,*)'Input a int:'
!!$  read(*,*) i
!!$  write(*,*)'i2c:         data.'//i2c(i)//'.data'
!!$  write(*,*)'trim(i2c):   data.'//trim(i2c(i))//'.data'
!!$  write(*,*)'i2c_p:       data.',i2c_p(i),'.data'
!!$  write(*,*)'itoc:        data.'//itoc(i)//'.data'
!!$  write(*,*)'trim(itoc):  data.'//trim(itoc(i))//'.data'
!!$  write(*,*)'int2str:     data.'//int2str(i,dig(i))//'.data'
!!$
!!$end program t

!!$program test
!!$  character (len=40) :: format
!!$  integer, dimension(6) :: a = (/ 1, 3, 10, 7314, 6222, 50 /)
!!$
!!$  write( format, "(A,6(A,I1),A)" ) "(", ("I", int(log10(real(a(i))))+1, i=1,6),")"
!!$  write(*, *) format
!!$  write(*, format) (a(i), i=1,6)
!!$
!!$end program test
