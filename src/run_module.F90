module run_module
  
  use input_module
  use wind_module
  use bed_module
  use moist_module

  implicit none

contains

  subroutine run_model(par)
    
    type(parameters), intent(inout) :: par
    integer*4 :: i, ti, nx, nt
    real*8 :: t, dt, dx, Ta
    real*8 :: tstart
    real*8, dimension(:), allocatable :: x, z, zm, u
    real*8, dimension(:,:), allocatable :: u_th, m
    real*8, dimension(:), allocatable :: rho, dist
    real*8, dimension(:,:), allocatable :: Cu, Ct
    real*8, dimension(:,:,:), allocatable :: mass
    integer*4, parameter :: fid=20

    write(*,*) 'Initialization started...'

    ! bed
    call generate_bed(par, x, z)
    call generate_bedcomposition(par, x, z)
    open(unit=fid, file="bed.out", action="write", status="replace", form="unformatted")
    write(fid) x
    write(fid) z
    close(fid)

    ! wind
    call generate_wind(par, u)
    open(unit=fid, file="wind.out", action="write", status="replace", form="unformatted")
    write(fid) u
    close(fid)

    ! courant check
    write(0, '(a, f4.2)') " Courant condition: ", maxval(u) / par%dx * par%dt   
    if (par%dx / par%dt < maxval(u)) then
       write(0, '(a)') " Courant condition violated. Please adapt numerical parameters."
       stop 1
    end if

    ! moist
    call generate_moist(par, zm, m)
    open(unit=fid, file="moist.out", action="write", status="replace", form="unformatted")
    write(fid) m
    close(fid)

    open(unit=fid, file="moist_z.out", action="write", status="replace", form="unformatted")
    write(fid) zm
    close(fid)

    ! time
    t = 0
    par%nt = nint(par%tstop / par%dt)

    ! space
    allocate(u_th(par%nx+1, par%nfractions))
    allocate(Cu(par%nx+1, par%nfractions))
    allocate(Ct(par%nx+1, par%nfractions))
    allocate(mass(par%nx+1, par%nlayers, par%nfractions))
    u_th = 0.0
    Cu = 0.0
    Ct = 0.0
    mass = 0.0

    ! fractions
    allocate(rho(par%nfractions))
    allocate(dist(par%nfractions))

    rho = par%grain_size
    dist = par%grain_dist

    open(unit=fid, file="dims.out", action="write", status="replace", form="unformatted")
    write(fid) par%nx
    write(fid) par%dx
    write(fid) par%nt
    write(fid) par%dt
    write(fid) par%nfractions
    write(fid) par%nlayers+2
    close(fid)

    write(*,*) 'Model run started...'

    ! model
    open(unit=fid+1, file="Ct.out", action="write", status="replace", form="unformatted")
    open(unit=fid+2, file="Cu.out", action="write", status="replace", form="unformatted")
    open(unit=fid+3, file="z.out", action="write", status="replace", form="unformatted")
    open(unit=fid+4, file="u_th.out", action="write", status="replace", form="unformatted")
    open(unit=fid+5, file="mass.out", action="write", status="replace", form="unformatted")
!    open(unit=fid+6, file="mfrac.out", action="write", status="replace", form="unformatted")
!    open(unit=fid+7, file="vfrac.out", action="write", status="replace", form="unformatted")
    open(unit=fid+8, file="d10.out", action="write", status="replace", form="unformatted")
    open(unit=fid+9, file="d50.out", action="write", status="replace", form="unformatted")
    open(unit=fid+10, file="d90.out", action="write", status="replace", form="unformatted")

    ! dimensions
    nx = par%nx
    nt = par%nt
    dt = par%dt
    dx = par%dx
    Ta = par%Tp / par%dt

    tstart = get_time()
    do ti=1,par%nt

       ! log progress
       if ( mod(dble(ti), par%nt/10.d0) < 1.d0 ) then
          call write_progress(ti, par%nt, tstart)
       end if

       ! update threshold
       u_th = par%u_th
       call compute_threshold_grainsize(par, u_th)
       call compute_threshold_bedslope(par, x, z, u_th)
       call compute_threshold_moisture(par, zm, m(ti,:), z, u_th)

       ! get available mass
       mass = get_layer_mass()

       do i=1,par%nfractions

          ! compute transport capacity by wind, including thresholds
          Cu(:,i) = max(0.d0, 1.5e-4 * ((u(ti) - u_th(:,i))**3) / (u(ti) * par%VS))

          ! compute sediment advection by wind
          Ct(2:nx+1,i) = -par%VS * u(ti) * (Ct(2:nx+1,i) - Ct(1:nx,i)) / dx * dt + &
                          Cu(2:nx+1,i) / Ta + &
                          Ct(2:nx+1,i) * (1-1/Ta)

          ! limit advection by available mass
          Ct(2:nx+1,i) = min(mass(2:nx+1,1,i), Ct(2:nx+1,i))

       end do

       ! update bed elevation
       z = update_bed(z, -(Cu - Ct) / Ta, rho, par%dt, par%morfac)

       ! write output
       if ( mod(ti, nint(par%tout / dt)) < 1.d0 ) then
          write(fid+1) Ct
          write(fid+2) Cu
          write(fid+3) z
          write(fid+4) u_th
          write(fid+5) get_layer_mass()
 !         write(fid+6) get_layer_massfraction(par)
 !         write(fid+7) get_layer_volumefraction(par)
          write(fid+8) get_layer_percentile(par, 0.1d0)
          write(fid+9) get_layer_percentile(par, 0.5d0)
          write(fid+10) get_layer_percentile(par, 0.9d0)
       end if
       
       t = t + par%dt
    end do

    close(fid+1)
    close(fid+2)
    close(fid+3)
    close(fid+4)
    close(fid+5)
 !   close(fid+6)
 !   close(fid+7)
    close(fid+8)
    close(fid+9)
    close(fid+10)

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

    write(*,'(f5.1,a2,a8,a3,a8,a3,a8)') 100.d0*p, '% ', format_time(dt1), ' / ', format_time(dt2), ' / ', format_time(dt3)

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
    s = nint(tm - h * 3600 - m * 60)

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
