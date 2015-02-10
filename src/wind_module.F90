module wind_module

  use constants_module
  use utils_module
  use input_module

  implicit none

contains

  subroutine generate_wind(par, uout)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, i, it, nt, l
    real*8, dimension(:), pointer, intent(out) :: uout
    real*8, dimension(:), allocatable :: u
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
        
    ! checks
    if (sum(duration) < par%tstop) then
       write(*,*) "ERROR: wind definition file too short"
       stop 1
    end if

    ! maximize optimization tries
    do n=1,20

       ! allocate arrays
       nt = nint(par%tstop / par%dt)
       if (allocated(u)) deallocate(u)
       allocate(u(nt))
   
       ! compute random time series
       i = 1
       it = 1
       do while (it < nt)
          l = nint(max(par%dt, rand_normal(g_m(i), g_std(i))) / par%dt)
          u(it:min(nt,it+l)) = max(0.d0, rand_normal(u_m(i), u_std(i)))
          it = it + l
          if (it > sum(duration(1:i)) / par%dt) then
             i = i + 1
          end if
       end do
   
       ! courant check
       if (par%CFL > 0.d0) then
          if (abs(maxval(u) / par%dx * par%dt - par%CFL) < .005) then
             write(0, '(a, f6.4)') " Adapted timestep based on CFL condition: ", par%dt
             exit
          end if
          par%dt = par%CFL * par%dx / maxval(u)
       else
          exit
       end if

    end do

    ! make time step fit with output time step
    par%dt = par%tout / ceiling(par%tout / par%dt)

    allocate(uout(size(u)))
    uout = u

  end subroutine generate_wind

end module wind_module
