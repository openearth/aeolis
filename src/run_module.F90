module run_module
  
  use output_module
  use input_module
  use wind_module
  use bed_module
  use moist_module

  implicit none

contains

  subroutine run_model(par)
    
    type(parameters), intent(inout) :: par
    type(variables), dimension(:), allocatable :: var
    integer*4 :: i, j, ti, nx, nt
    real*8 :: t, dt, dx, Ta, n
    real*8 :: tstart, tlog
    real*8, pointer :: wind
    real*8, dimension(:), pointer :: x, z, u, moist_map
    real*8, dimension(:,:), pointer :: uth
    real*8, dimension(:), pointer :: rho, dist
    real*8, dimension(:,:), pointer :: Cu, Ct, supply
    real*8, dimension(:,:), pointer :: d10, d50, d90
    real*8, dimension(:,:,:), pointer :: mass
    real*8, dimension(:), allocatable :: x_tmp, z_tmp, zmoist
    real*8, dimension(:,:), allocatable :: Ct2, moist
    integer*4, parameter :: fid=20

    write(*,*) 'Initialization started...'
    
    ! bed
    write(*,*) 'Generating bed profile and composition...'
    call generate_bed(par, x_tmp, z_tmp)
    call get_pointer(var, 'x', (/par%nx+1/), x)
    call get_pointer(var, 'z', (/par%nx+1/), z)
    
    x = x_tmp
    z = z_tmp
    
    call generate_bedcomposition(par, x, z)
    open(unit=fid, file="bed.in", action="write", status="replace", form="unformatted")
    write(fid) x
    write(fid) z
    close(fid)

    ! wind
    write(*,*) 'Generating bed wind time series...'
    call get_pointer(var, 'u', (/par%nx+1/), u)
    call generate_wind(par, u)
    open(unit=fid, file="wind.in", action="write", status="replace", form="unformatted")
    write(fid) u
    close(fid)

    ! courant check
    write(0, '(a, f4.2)') " Courant condition: ", maxval(u) / par%dx * par%dt   
    if (par%dx / par%dt < maxval(u)) then
       write(0, '(a)') " Courant condition violated. Please adapt numerical parameters."
       stop 1
    end if

    ! time
    t = 0
    par%nt = int(par%tstop / par%dt)
    par%ntout = int(par%tstop / par%tout) + 1

    ! moist
    write(*,*) 'Generating moisture time series...'
    call generate_moist(par, zmoist, moist)
    open(unit=fid, file="moist.in", &
         action="write", status="replace", form="unformatted")
    do i = 1,par%nt
       write(fid) moist(i,:)
    end do
    close(fid)

    ! fractions
    call get_pointer(var, 'rho',    (/par%nfractions/), rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), dist)
    rho = par%rhop
    dist = par%grain_dist

    ! variables
    call get_pointer(var, 'Cu',     (/par%nfractions, par%nx+1/), Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%nx+1/), Ct)
    call get_pointer(var, 'uth',    (/par%nfractions, par%nx+1/), uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers+2, par%nx+1/), mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%nx+1/), supply)

    ! extra output
    call get_pointer(var, 'wind', (/0/), wind)
    call get_pointer(var, 'moist_map', (/par%nx+1/), moist_map)
    call get_pointer(var, 'd10', (/par%nlayers+2, par%nx+1/), d10)
    call get_pointer(var, 'd50', (/par%nlayers+2, par%nx+1/), d50)
    call get_pointer(var, 'd90', (/par%nlayers+2, par%nx+1/), d90)

    allocate(Ct2(par%nfractions, par%nx+1))
    Ct2 = 0.d0

    call write_dimensions(par)

    write(*,*) 'Model run started...'

    ! output
    call output_init(var, par%outputvars)
    
    ! dimensions
    nx = par%nx
    nt = par%nt
    dt = par%dt
    dx = par%dx
    Ta = par%Tp / par%dt

    tstart = get_time()
    tlog = tstart
    do ti=1,par%nt

       ! log progress
       if ( mod(dble(ti), par%nt/10.d0) < 1.d0 .or. get_time()-tlog > 60) then
          call write_progress(ti, par%nt, tstart)
          tlog = get_time()
       end if

       ! update threshold
       uth = par%u_th
       call compute_threshold_grainsize(par, uth)
       call compute_threshold_bedslope(par, x, z, uth)
       call compute_threshold_moisture(par, zmoist, moist(ti,:), z, uth)

       ! mix top layer of wet cells
       call mix_toplayer(par, uth)

       ! get available mass
       mass = get_layer_mass()

       do j=1,par%nx+1

          do i=1,par%nfractions

             ! compute transport capacity by wind, including thresholds
             Cu(i,j) = max(0.d0, 1.5e-4 * (u(ti) - uth(i,j))**3 / (u(ti) * par%VS))

             ! limit advection by available mass
             supply(i,j) = min(mass(i,1,j), par%accfac * Cu(i,j) - Ct(i,j)) / Ta

          end do

          ! scale supply to availability of fractions
          where (supply(:,j) > 0.d0)
             dist = mass(:,1,j)
          elsewhere
             dist = 0.d0
          end where
          dist = dist / max(1e-10, sum(dist))
          where (supply(:,j) > 0.d0)
             supply(:,j) = supply(:,j) * dist
          end where
             
       end do

       ! set supply in first cell to zero, since it cannot be picked up
       supply(:,1) = 0.d0

       do j=1,par%nx

          do i=1,par%nfractions

             ! compute sediment advection by wind
             Ct2(i,j+1) = max(0.d0, -par%VS * u(ti) * (Ct(i,j+1) - Ct(i,j)) / dx * dt + &
                  Ct(i,j+1) + supply(i,j+1))

          end do

       end do

       Ct = Ct2

       ! update bed elevation
       z = update_bed(z, -supply, rho, par%dt)

       ! incremental output
       call output_update(var)

       ! write output
       if ( ti == 1 .or. par%tout < dt .or. mod(ti, nint(par%tout / dt)) < 1.d0 ) then

          ! update derived variables
          if (is_output(var, 'mass')) mass = get_layer_mass()
          if (is_output(var, 'wind')) wind = u(ti)
          if (is_output(var, 'moist_map')) &
               moist_map = map_moisture(par, zmoist, moist(ti,:), z)
          if (is_output(var, 'd10')) d10 = get_layer_percentile(par, 0.1d0)
          if (is_output(var, 'd50')) d50 = get_layer_percentile(par, 0.5d0)
          if (is_output(var, 'd90')) d90 = get_layer_percentile(par, 0.9d0)

          call output_write(var)
          call output_clear(var)

       end if
       
       t = t + par%dt
    end do

    call output_close(var)
    
    write(*,*) 'Done.'

  end subroutine run_model
  
  subroutine write_progress(ti, nt, tstart)

    integer*4, intent(in) :: ti, nt
    real*8, intent(in) :: tstart
    real*8 :: p, dt1, dt2, dt3
    
    p = dble(ti) / nt
    dt1 = get_time() - tstart
    dt2 = dt1 / p
    dt3 = dt2 * (1 - p)

    write(*,'(f5.1,a2,a8,a3,a8,a3,a8)') &
         100.d0*p, '% ', &
         format_time(dt1), ' / ', &
         format_time(dt2), ' / ', &
         format_time(dt3)

  end subroutine write_progress

  function get_time() result (tm)

    integer*4 :: count,count_rate,count_max
    real*8 :: tm

    call system_clock (count,count_rate,count_max)
    tm = dble(count)/count_rate

  end function get_time

  function format_time(tm) result (str)

    real*8, intent(in) :: tm
    integer*4 :: h, m, s
    character(2) :: h_s, m_s, s_s
    character(8) :: str
    
    h = int(tm / 3600)
    m = int((tm - h * 3600) / 60)
    s = int(tm - h * 3600 - m * 60)

    if (h < 10) then
       write(h_s, '("0",i0)') h
    else
       write(h_s, '(i0)') h
    end if

    if (m < 10) then
       write(m_s, '("0",i0)') m
    else
       write(m_s, '(i0)') m
    end if

    if (s < 10) then
       write(s_s, '("0",i0)') s
    else
       write(s_s, '(i0)') s
    end if

    str = h_s // ':' // m_s // ':' // s_s

  end function format_time

end module run_module
