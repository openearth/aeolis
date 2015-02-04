module utils_module

  use constants_module

  implicit none

contains

  
  subroutine sort(x, x2)
    
    integer :: n, i, idx
    real*8 :: tmp
    real*8, dimension(:), intent(inout) :: x
    real*8, dimension(:), intent(inout), optional :: x2

    n = size(x)

    do i=1,n
       idx = i + find_minimum(x(i:n)) - 1

       call swap(x(i), x(idx))

       if (present(x2)) then
          call swap(x2(i), x2(idx))
       end if
    end do

  end subroutine sort

  
  subroutine swap(x, y)
    
    real*8, intent(inout) :: x, y
    real*8 :: tmp

    tmp = x
    x = y
    y = tmp

  end subroutine swap

  
  subroutine split_path(str, fpath, fname)

    character(*),    intent(in) :: str
    character(slen), intent(out) :: fpath, fname
    integer :: pos, pos0

    pos0 = 1
    
    do
       pos = index(str(pos0:), '/')
       if (pos == 0) then
          fpath = str(:(pos0-1))
          fname = str(pos0:)
          exit
       else
          pos0 = pos0 + pos
       end if
    end do
    
  end subroutine split_path

  
  function find_minimum(x) result (loc)

    integer :: n, i, loc
    real*8 :: mn
    real*8, dimension(:), intent(in) :: x

    n = size(x)

    loc = 1
    mn = x(loc)
    do i=2,n
       if (x(i) < mn) then
          loc = i
          mn = x(loc)
       end if
    end do

  end function find_minimum

  
  function linear_interp(x, y, xx) result (yy)

    integer :: n
    real*8, dimension(:) :: x, y
    real*8 :: xx, yy
    real*8 :: a,  b, dyy
    integer :: j

    yy = 0.0d0

    n = size(x)

    if (n.le.0) return
    if (n.eq.1) then
       yy = y(1)
       return
    endif

    j = binary_search(x, xx)

    if (j .le. 0) then
       yy = y(1)
    elseif (j .ge. n) then
       yy = y(n)
    else
       a = x (j+1)
       b = x (j)
       if (a .eq. b) then
          dyy = 0.0d0
       else
          dyy = (y(j+1) - y(j)) / (a - b)
       endif
       yy = y(j) + (xx - x(j)) * dyy
    endif

  end function linear_interp

  
  function binary_search(xx, x) result (j)

    integer :: n
    real*8, dimension(:) :: xx
    real*8 :: x
    integer :: j
    integer :: jl, ju, jm
    logical :: l1, l2

    n = size(xx)

    jl = 0
    ju = n+1
    do while(ju-jl .gt. 1)
       jm = (ju+jl)/2
       l1 = xx(n) .gt. xx(1)
       l2 = x .gt. xx(jm)
       if ((l1 .and. l2) .or. (.not. (l1 .or. l2))) then
          jl = jm
       else
          ju = jm
       end if
    end do

    j = jl

    return

  end function binary_search

  
  function rand_normal(mean, stdev) result(c)
    
    real*8 :: mean, stdev, r, theta, c, x(2)
    
    if(stdev < 0.0d0) then
       write(*,*) "WARNING: standard deviation must be positive"
    else
       call random_number(x)
       r = (-2.0d0*log(x(1)))**0.5
       theta = 2.0d0*pi*x(2)
       c = mean + stdev*r*sin(theta)
    end if
  end function rand_normal

  
end module utils_module
