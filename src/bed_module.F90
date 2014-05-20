module bed_module

  use constants_module
  use utils_module
  use input_module

contains

  subroutine generate_bed(par, x, z)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, i
    real*8, dimension(:), allocatable, intent(out) :: x, z
    real*8, dimension(:), allocatable :: x_c, z_c

    fid = 99

    ! count lines
    n = -1
    ierr = 0
    open(fid, file=trim(par%bed_file))
    do while (ierr == 0)
       read(fid, *, iostat=ierr)
       n = n + 1
    end do
    rewind(fid)

    ! allocate arrays
    allocate(x_c(n))
    allocate(z_c(n))

    ! read data
    i = 1
    ierr = 0
    do while (ierr == 0)
       read(fid, *, iostat=ierr) x_c(i), z_c(i)
       i = i + 1
    end do
    close(fid)
        
    ! allocate arrays
    par%nx = nint(maxval(x_c) / par%dx)
    allocate(x(par%nx+1))
    allocate(z(par%nx+1))

    ! interpolate bed
    x = (/(i*par%dx, i=0, par%nx)/)
    do i=1, par%nx+1
       z(i) = linear_interp(x_c, z_c, x(i))
    end do
    
  end subroutine generate_bed

  subroutine compute_threshold_bedslope(par, x, z, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: x, z
    real*8, dimension(:), intent(inout) :: u_th
    integer :: i
    real*8 :: phi, theta

    phi = par%phi / 180.d0 * pi

    do i=1,par%nx
       theta = -atan((z(i+1) - z(i)) / (x(i+1) - x(i)))
       u_th(i) = sqrt((tan(phi) - tan(theta) / tan(phi) + cos(theta)) / \
                 (tan(phi)+1)) * u_th(i)
    end do

  end subroutine compute_threshold_bedslope

end module bed_module
