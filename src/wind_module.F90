module wind_module

  use constants_module
  use utils_module
  use input_module

  implicit none

contains

  subroutine generate_wind(par, gusty_wind)

    type(parameters), intent(inout) :: par
    type(windstat), dimension(:), allocatable :: wind
    type(windspeed), dimension(:), allocatable, intent(out) :: gusty_wind
    real*8, dimension(7) :: tmp
    integer*4 :: fid, ierr, n, i, it, nt, l

    if (trim(par%wind_file) /= '') then

       fid = 88

       ! count lines
       n = -1
       ierr = 0
       open(fid, file=trim(par%wind_file))
       do while (ierr == 0)
          read(fid, *, iostat=ierr)
          n = n + 1
       end do
       rewind(fid)

       ! allocate arrays
       allocate(wind(n))

       ! read data
       i = 1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr) tmp(:)
          wind(i)%t = min(par%tstop+1.d0, tmp(1))
          wind(i)%u_mean = tmp(2)
          wind(i)%u_std = tmp(3)
          wind(i)%gust_mean = tmp(4)
          wind(i)%gust_std = tmp(5)
          wind(i)%dir_mean = tmp(6) / 180.d0 * pi
          wind(i)%dir_std = tmp(7) / 180.d0 * pi

          if (wind(i)%t > par%tstop) exit
       
          i = i + 1
       end do
       close(fid)

       if (wind(i)%t < par%tstop) then
          write(*,*) "ERROR: wind definition file too short"
          stop 1
       end if
       
       ! determine time axis
       wind(1:n-1)%duration = wind(2:n)%t - wind(1:n-1)%t
       wind(n)%duration = par%tstop - wind(n)%t

       ! simulate wind gusts
       if (par%gusts) then
          call simulate_gusts(wind, gusty_wind)
       else
          allocate(gusty_wind(n))
          gusty_wind%t = wind(1:n)%t
          gusty_wind%duration = wind(1:n)%duration
          gusty_wind%dir = wind(1:n)%dir_mean
          gusty_wind%u = wind(1:n)%u_mean
       end if

    else

       allocate(gusty_wind(2))
       gusty_wind(1)%t = 0.d0
       gusty_wind(2)%t = par%tstop
       gusty_wind%duration = par%tstop
       gusty_wind%dir = 0.d0
       gusty_wind%u = 0.d0

    end if

  end subroutine generate_wind

  subroutine interpolate_wind(par, wind, t, uw, udir)

    type(parameters), intent(in) :: par
    type(windspeed), dimension(:), intent(in) :: wind
    real*8, intent(in) :: t
    real*8, dimension(:), intent(out) :: uw, udir
    integer*4 :: i
    real*8 :: tm

    ! repeat time series if necessary
!    tm = mod(t, sum(wind%duration))
    
    if (par%gusts) then
       uw = linear_interp(wind%t, wind%u, t)
       udir = linear_interp(wind%t, wind%dir, t)
    else
       i = max(1, binary_search(wind%t, t))
       uw = wind(i)%u
       udir = wind(i)%dir
    end if
    
  end subroutine interpolate_wind

  subroutine simulate_gusts(wind, gusty_wind)

    type(windstat), dimension(:), intent(in) :: wind
    type(windspeed), dimension(:), allocatable, intent(out) :: gusty_wind
    real*8, dimension(:), allocatable :: l, u, d
    integer*4 :: i, n, m
    real*8 :: p

    m = nint(2 * sum(wind%duration / wind%gust_mean))

    allocate(l(m))
    allocate(u(m))
    allocate(d(m))
    l = 0.d0
    u = 0.d0
    d = 0.d0

    n = 1
    do i = 1,size(wind)
       do while (sum(l) < sum(wind(1:i)%duration))
          l(n) = max(0.d01, rand_normal(wind(i)%gust_mean, wind(i)%gust_std))
          u(n) = max(0.d0, rand_normal(wind(i)%u_mean, wind(i)%u_std))
          d(n) = rand_normal(wind(i)%dir_mean, wind(i)%dir_std)
          n = n+1
       end do
    end do

    allocate(gusty_wind(n-1))

    do i = 1,n-1
       gusty_wind(i)%duration = l(i)
       gusty_wind(i)%u = u(i)
       gusty_wind(i)%dir = d(i)
    end do

    ! determine time axis
    gusty_wind(2:n)%t = cumsum(gusty_wind(1:n-1)%duration)

  end subroutine simulate_gusts
  
end module wind_module
