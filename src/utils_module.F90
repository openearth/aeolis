module utils_module

  use constants_module

  implicit none

contains

  subroutine assert(condition)

    logical :: condition
    
    if (.not. condition) then
       write(0,*) 'Assertion failed.'
       stop 1
    end if
    
  end subroutine assert
  
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

  subroutine split(str, arr, sep)

    character(*), intent(in) :: str
    character(*), intent(in), optional :: sep
    character(1) :: sep0
    character(10), intent(inout), dimension(:), allocatable :: arr
    integer*4 :: i, n, pos, pos0

    if (present(sep)) then
       sep0 = sep
    else
       sep0 = ' '
    end if
    
    n = 0
    i = 1

    pos0 = 1
    do
       pos = index(trim(str(pos0:)), sep0)
       if (pos == 0) then
          exit
       else if (pos .gt. 1) then
          n = n + 1
       end if
       pos0 = pos0 + pos
    end do

    allocate(arr(n+1))

    pos0 = 1
    do i = 1, n+1
       pos = index(str(pos0:), sep0)
       if (pos == 0) then
          arr(i) = trim(str(pos0:))
          exit
       end if
       arr(i) = trim(str(pos0:(pos0+pos-2)))
       pos0 = pos0 + pos
    end do

  end subroutine split

  function isin(arr, x) result (found)

    character(*), dimension(:) :: arr
    character(*) :: x
    integer*4 :: i
    logical :: found

    found = .false.
    do i = 1,size(arr)
       if (arr(i) == x) then
          found = .true.
          exit
       end if
    end do

  end function isin

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

  subroutine split_var(var_name, name, type)

    character(*), intent(in) :: var_name
    character(10), dimension(:), allocatable :: var_name_parts
    character(slen), intent(out) :: name, type

    call split(var_name, var_name_parts, '.')
    if (size(var_name_parts) == 1) then
       name = trim(var_name_parts(1))
       type = ''
    else
       name = trim(var_name_parts(1))
       type = trim(var_name_parts(2))
    end if

  end subroutine split_var
  
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
    l1 = xx(n) .gt. xx(1)
    do while(ju-jl .gt. 1)
       jm = (ju+jl)/2
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

  function dist_normal(x, mean, stdev) result (p)

    real*8 :: mean, stdev
    real*8, dimension(:) :: x
    real*8, dimension(:), allocatable :: p

    allocate(p(size(x)))

    if (stdev < 0.d0) then
       write(*,*) "WARNING: standard deviation must be positive"
    else
       p = exp((x - mean) / stdev)
    end if

  end function dist_normal
    
  function rand_normal(mean, stdev) result(c)
    
    real*8 :: mean, stdev, r, theta, c, x(2)
    
    if (stdev < 0.d0) then
       write(*,*) "WARNING: standard deviation must be positive"
    else
       call random_number(x)
       r = (-2.0d0*log(x(1)))**0.5
       theta = 2.0d0*pi*x(2)
       c = mean + stdev*r*sin(theta)
    end if
  end function rand_normal

  function first_exceedance(x, threshold) result (ith)

    real*8, dimension(:,:) :: x
    real*8 :: threshold
    integer*4 :: i, j, ith

    ith = 1
    do i = 1,size(x,1)
       do j = 1,size(x,2)
          if (x(i,j) >= threshold) then
             ith = i
             return
          end if
       end do
    end do

  end function first_exceedance

  function cumsum(x) result (xc)

    real*8, dimension(:) :: x
    real*8, dimension(:), allocatable :: xc
    integer*4 :: i, n
    
    n = size(x)
    allocate(xc(n))
    xc = 0.d0
    
    do i=1,n
       xc(i) = sum(x(1:i))
    end do

  end function cumsum
  
end module utils_module
