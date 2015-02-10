module moist_module

  use constants_module
  use utils_module
  use input_module

  implicit none

contains

  subroutine generate_moist(par, z, m)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, nz, i, j, nt, file_pos
    real*8, dimension(:,:), allocatable, intent(out) :: m
    real*8, dimension(:), allocatable, intent(out) :: z
    real*8, dimension(:,:), allocatable :: m0
    real*8, dimension(:), allocatable :: t
    character(2) :: fmt

    fid = 77

    ! read format
    open(fid, file=trim(par%moist_file))
    read(fid, *, iostat=ierr) fmt
    if (fmt .ne. 'TZ') then
       write(*,*) "ERROR: Unknown moist format: ", fmt
       stop 1
    end if

    ! read z values
    nz = 0
    read(fid, *, iostat=ierr) nz
    allocate(z(nz))
    do i = 1,nz
       read(fid, *, iostat=ierr) z(i)
    end do

    ! count lines
    n = -1
    ierr = 0
    do while (ierr == 0)
       read(fid, *, iostat=ierr)
       n = n + 1
    end do
    rewind(fid)

    ! allocate arrays
    allocate(t(n))
    allocate(m0(n,nz+2))

    ! read data
    ierr = 0
    do i = 1,nz+2
       read(fid, *, iostat=ierr)
    end do

    i = 1
    ierr = 0
    do while (ierr == 0)
       read(fid, *, iostat=ierr) t(i), m0(i,:)
       i = i + 1
    end do
    close(fid)

    ! checks
    if (maxval(t) < par%tstop) then
       write(*,*) "ERROR: moist definition file too short"
       stop 1
    end if

    ! allocate arrays
    nt = nint(par%tstop / par%dt)
    allocate(m(nt,nz+2))

    ! interpolate time series
    do j = 1,nz+2
       do i = 1,nt
          m(i,j) = linear_interp(t, m0(:,j), i*par%dt)
       end do
    end do

  end subroutine generate_moist

  subroutine compute_threshold_moisture(par, zm, m, z, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: zm, m, z
    real*8, dimension(par%nx+1) :: mg
    real*8, dimension(:,:), intent(inout) :: u_th
    real*8, dimension(par%nfractions, par%nx+1) :: u_th_m
    integer :: i, n

    mg = map_moisture(par, zm, m, z)

    ! convert from volumetric content (percentage of volume) to
    ! geotechnical mass content (percentage of dry mass)
    mg = mg * par%rhow / (par%rhop * (1 - par%porosity))

    u_th_m = 0.d0
    do i = 1,par%nfractions
       if (par%method_moist .eq. 'belly_johnson') then
          u_th_m(i,:) = u_th(i,:) * max(1.d0, 1.8+0.6*log10(mg))
       else if (par%method_moist .eq. 'hotta') then
          u_th_m(i,:) = u_th(i,:) + 7.5 * mg
       else
          write(*, '(a,a)') "ERROR: Unknown moisture formulation: ", par%method_moist
          stop 1
       end if

       where (mg > 0.005)
          u_th(i,:) = u_th_m(i,:)
       end where

       where (mg > 0.064)
          ! should be .04 according to Pye and Tsoar
          ! should be .64 according to Delgado-Fernandez (10% vol.)
          u_th(i,:) = 999.d0
       end where
    end do

  end subroutine compute_threshold_moisture

  function map_moisture(par, zm, m, z) result (mg)
    
    type(parameters) :: par
    real*8, dimension(:) :: zm, m, z
    real*8, dimension(par%nx+1) :: mg
    integer :: i, n

    n = size(m)
    do i = 1,par%nx+1
       if (z(i) < minval(zm)) then
          mg(i) = m(1)
       elseif (z(i) > maxval(zm)) then
          mg(i) = m(n)
       else
          mg(i) = linear_interp(zm(2:n-1), m(2:n-1), z(i))
       end if
    end do

  end function map_moisture

end module moist_module
