module supply_module

  use constants_module
  use input_module

contains

  subroutine generate_supply(par, time, supply)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, i
    real*8, dimension(:), allocatable, intent(out) :: time
    real*8, dimension(:,:), allocatable, intent(out) :: supply
    real*8 :: dt, t_prev

    fid = 99

    n = -1
    ierr = 0
    open(fid, file=trim(par%supply_file))
    do while (ierr == 0)
       read(fid, *, iostat=ierr)
       n = n + 1
    end do
    rewind(fid)

    ! allocate arrays
    allocate(time(n+1))
    allocate(supply(n, par%nx+1))

    ! read data
    i = 1
    ierr = 0
    t_prev = 0
    do while (ierr == 0)
       read(fid, *, iostat=ierr) dt, supply(i,:)
       time(i+1) = t_prev + dt
       t_prev = time(i+1)
       i = i + 1
    end do
    close(fid)

  end subroutine generate_supply

end module supply_module
