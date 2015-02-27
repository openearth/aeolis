module run_module
! This is a docstring of run_module
  
  use output_module
  use input_module
  use wind_module
  use moist_module
  use bed_module
  use utils_module

  implicit none

contains

  subroutine run_model(par)
    
    type(parameters), intent(inout) :: par
    type(variables), dimension(:), allocatable :: var
    type(meteorology), dimension(:), allocatable :: meteo
    integer*4 :: i, j, n, ti, nx, nt
    real*8 :: t, dt, dx, f, err, alpha
    real*8 :: tstart, tlog
    real*8, pointer :: wind
    real*8, dimension(:), pointer :: x, z
    real*8, dimension(:,:), pointer :: uth
    real*8, dimension(:), pointer :: rho, dist
    real*8, dimension(:,:), pointer :: Cu, Ct, supply, moist
    real*8, dimension(:,:), pointer :: d10, d50, d90
    real*8, dimension(:,:,:), pointer :: mass
    real*8, dimension(:), allocatable :: x_tmp, z_tmp, u, tide, frac
    real*8, dimension(:,:), allocatable :: Ct2, Ct2_prev
    integer*4, parameter :: fid=20
    
    write(*,*) 'Initialization started...'
    
    ! bed
    write(*,*) 'Generating bed profile and composition...'
    call generate_bed(par, x_tmp, z_tmp)
    call get_pointer(var, 'x', (/par%nx+1/), x)
    call get_pointer(var, 'z', (/par%nx+1/), z)
    
    x = x_tmp
    z = z_tmp
    deallocate(x_tmp)
    deallocate(z_tmp)

    call generate_bedcomposition(par, x, z)
    open(unit=fid, file=trim(par%output_dir) // "bed.in", &
         action="write", status="replace", form="unformatted")
    write(fid) x
    write(fid) z
    close(fid)

    ! wind
    write(*,*) 'Generating bed wind time series...'
    call generate_wind(par, u)
    open(unit=fid, file=trim(par%output_dir) // "wind.in", &
         action="write", status="replace", form="unformatted")
    write(fid) u
    close(fid)

    ! courant check
    if (trim(par%scheme) .eq. 'explicit') then
       write(0, '(a, f4.2)') " Courant condition: ", maxval(u) / par%dx * par%dt   
       if (par%dx / par%dt < maxval(u)) then
          write(0, '(a)') " Courant condition violated. Please adapt numerical parameters."
          stop 1
       end if
    end if
    
    ! time
    t = 0
    par%nt = int(par%tstop / par%dt)
    par%ntout = int(par%tstop / par%tout) + 1

    ! moist
!    write(*,*) 'Generating moisture time series...'
!    call generate_moist(par, zmoist, moist)
!    open(unit=fid, file=trim(par%output_dir) // "moist.in", &
!         action="write", status="replace", form="unformatted")
!    do i = 1,par%nt
!       write(fid) moist(i,:)
!    end do
    !    close(fid)

    ! tide and meteo
    write(*,*) 'Generating tide and meteo time series...'
    call generate_tide(par, tide)
    call generate_meteo(par, meteo)

    ! fractions
    call get_pointer(var, 'rho',    (/par%nfractions/), rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), dist)

    rho = par%rhom
    dist = par%grain_dist

    ! variables
    call get_pointer(var, 'Cu',     (/par%nfractions, par%nx+1/), Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%nx+1/), Ct)
    call get_pointer(var, 'uth',    (/par%nfractions, par%nx+1/), uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers+2, par%nx+1/), mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%nx+1/), supply)
    call get_pointer(var, 'moist',  (/par%nlayers+2, par%nx+1/), moist)

    ! extra output
    call get_pointer(var, 'u', (/0/), wind)
    call get_pointer(var, 'd10', (/par%nlayers+2, par%nx+1/), d10)
    call get_pointer(var, 'd50', (/par%nlayers+2, par%nx+1/), d50)
    call get_pointer(var, 'd90', (/par%nlayers+2, par%nx+1/), d90)

    allocate(frac(par%nfractions))
    allocate(Ct2(par%nfractions, par%nx+1))
    allocate(Ct2_prev(par%nfractions, par%nx+1))
    Ct2 = 0.d0
    Ct2_prev = 0.d0
    frac = 0.d0

    call write_dimensions(par)

    write(*,*) 'Model run started...'

    ! output
    call output_init(var, par%outputvars, par%output_dir)
    
    ! dimensions
    nx = par%nx
    nt = par%nt
    dt = par%dt
    dx = par%dx

    tstart = get_time()
    tlog = tstart
    do ti=1,par%nt

       ! log progress
       if ( mod(dble(ti), par%nt/10.d0) < 1.d0 .or. get_time()-tlog > 60) then
          call write_progress(ti, par%nt, tstart)
          tlog = get_time()
       end if

       ! update moisture contents
       call update_moisture(par, z, tide(ti), meteo(1), u(ti), moist)
       call mix_toplayer(par, z, tide(ti))
       
       ! update threshold
       uth = par%u_th
       call compute_threshold_grainsize(par, uth)
!       call compute_threshold_bedslope(par, x, z, uth)
       call compute_threshold_moisture(par, moist(1,:), uth)

       ! get available mass
       mass = get_layer_mass()

       ! compute transport capacity by wind, including thresholds
       alpha = (0.174 / log10(par%z0/par%k))**3
       Cu = max(0.d0, alpha * par%Cb * par%rhoa / par%g * (u(ti) - uth)**3 / (u(ti) * par%VS))

       if (trim(par%scheme) .eq. 'explicit') then
          do j=2,par%nx+1

             ! compute supply based on sediment availability
             supply(:,j) = compute_supply(par, mass(:,1,j), par%accfac * Cu(:,j), Ct(:,j))

             do i=1,par%nfractions
             
                ! compute sediment advection by wind
                Ct2(i,j) = max(0.d0, -par%VS * u(ti) * (Ct(i,j) - Ct(i,j-1)) * dt / dx + &
                     Ct(i,j) + supply(i,j))

             end do
          end do
       else
          do n=1,par%max_iter

             Ct2_prev = Ct2
             
             do j=2,par%nx+1

                ! compute supply based on sediment availability
                supply(:,j) = compute_supply(par, mass(:,1,j), par%accfac * Cu(:,j), Ct2(:,j))

                do i=1,par%nfractions
                
                   ! compute sediment advection by wind
                   Ct2(i,j) = max(0.d0, (par%VS * u(ti) * Ct2(i,j-1) * dt / dx + &
                        Ct(i,j) + supply(i,j)) / (1 + u(ti) * dt / dx))

                end do
             end do

             ! exit iteration if change is negligible
             err = sum(abs(Ct2 - Ct2_prev) / max(1e-10, Ct2_prev)) / &
                  (par%nx+1) / par%nfractions
             if (err .le. par%max_error) exit

          end do

          if (err .gt. par%max_error) &
               write(0, '(a,i4,a,f0.4,a)') &
               "WARNING: iteration not converged (i: ", ti, "; error: ", err, ")"

       end if

       Ct = Ct2

       ! add sediment deposit
       do i = 1,par%nfractions
          where (z < tide(ti))
             supply(i,:) = supply(i,:) - par%Cw * &
                  min(par%w * par%dt, tide(ti) - z) * &
                  par%grain_dist(i) / max(1e-10, sum(par%grain_dist))
          end where
       end do
       supply(:,1) = 0.d0

       ! update bed elevation
       z = update_bed(z, -supply, rho, par%dt)
       call sweep_toplayer(par)

       ! incremental output
       call output_update(var)

       ! write output
       if (ti == 1 .or. par%tout < dt .or. mod(ti, nint(par%tout / dt)) < 1.d0) then

          ! update derived variables
          if (is_output(var, 'mass')) mass = get_layer_mass()
          if (is_output(var, 'u')) wind = u(ti)
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

  function compute_supply(par, mass, Cu, Ct) result(supply)

    type(parameters), intent(in) :: par ! parameters structure
    real*8, dimension(:), intent(in) :: mass, Cu, Ct
    real*8, dimension(size(mass)) :: dist, dist2, supply

    dist = Ct / max(1e-10, Cu)

    if (sum(dist) < 1.d0) then ! deposition
    
       ! compute sediment distribution in bed
       dist2 = mass / max(1e-10, sum(mass))

       ! compute new sediment distributuion in the air
       dist = dist + dist2 * (1.d0 - sum(dist))

    end if

    ! compute distribution in air
    if (sum(dist) == 0.d0) dist = 1.d0
    dist = dist / sum(dist)

    call assert(abs(sum(dist) - 1.d0) < 1e-10)

    ! determine weighed supply
    supply = (Cu * dist - Ct) / par%Tp * par%dt

    ! limit advection by available mass
    supply = min(mass, supply)

  end function compute_supply

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
