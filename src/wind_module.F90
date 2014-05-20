module wind_module

  use constants_module
  use utils_module
  use input_module

contains

  subroutine generate_wind(par, u)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, i, nt
    real*8, dimension(:), allocatable, intent(out) :: u
    real*8, dimension(:), allocatable :: duration, u_m, u_std, g_m, g_std

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
    allocate(duration(n))
    allocate(u_m(n))
    allocate(u_std(n))
    allocate(g_m(n))
    allocate(g_std(n))

    ! read data
    i = 1
    ierr = 0
    do while (ierr == 0)
       read(fid, *, iostat=ierr) duration(i), u_m(i), u_std(i), g_m(i), g_std(i)
       i = i + 1
    end do
    close(fid)
        
    ! allocate arrays
    nt = nint(sum(duration) / par%dt)
    allocate(u(nt))

    ! compute random time series
    i = 1
    it = 1
    do while (it < nt)
       l = nint(rand_normal(g_m(i), g_std(i)) / par%dt)
       u(it:min(nt,it+l)) = max(0.d0, rand_normal(u_m(i), u_std(i)))
       it = it + l
       if (it > sum(duration(1:i)) / par%dt) then
          i = i + 1
       end if
    end do

  end subroutine generate_wind

end module wind_module
