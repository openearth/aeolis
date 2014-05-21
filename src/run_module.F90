module run_module
  
  use input_module
  use wind_module
  use bed_module
  use supply_module

  implicit none

contains

  subroutine run_model(par)
    
    type(parameters), intent(inout) :: par
    integer*4 :: i, ti, si, n, l
    real*8 :: t, dt, dx, Ta
    real*8 :: tstart
    real*8, dimension(:), allocatable :: x, z, u_th
    real*8, dimension(:), allocatable :: u, t_supply
    real*8, dimension(:,:), allocatable :: supply
    real*8, dimension(:), allocatable :: Ct, Cu, Ca
    integer*4, parameter :: fid=20

    write(*,*) 'Initialization started...'

    ! parameters
    dt = par%dt
    dx = par%dx
    Ta = par%Tp / par%dt

    ! bed
    call generate_bed(par, x, z)
    call generate_bedcomposition(par, x, z)
    open(unit=fid, file="bed", action="write", status="replace", form="unformatted")
    write(fid) x
    write(fid) z
    close(fid)

    ! wind
    call generate_wind(par, u)
    open(unit=fid, file="wind", action="write", status="replace", form="unformatted")
    write(fid) u
    close(fid)

    ! supply
    call generate_supply(par, t_supply, supply)
    open(unit=fid, file="supply", action="write", status="replace", form="unformatted")
    write(fid) supply
    close(fid)

    ! time
    t = 0
    par%nt = nint(par%tstop / par%dt)

    ! space
    allocate(u_th(par%nx+1))
    allocate(Ct(par%nx+1))
    allocate(Cu(par%nx+1))
    allocate(Ca(par%nx+1))
    u_th = 0.d0
    Ct = 0.0
    Cu = 0.0
    Ca = par%S

    ! checks
    if (size(u) < par%nt) then
       write(*,*) "ERROR: wind definition file too short"
       stop 1
    end if

    write(*,*) 'Model run started...'

    ! update threshold
    u_th = par%u_th
    call compute_threshold_bedslope(par, x, z, u_th)

    ! model
    open(unit=fid+1, file="Ct", action="write", status="replace", form="unformatted")
    open(unit=fid+2, file="Cu", action="write", status="replace", form="unformatted")
    open(unit=fid+3, file="Ca", action="write", status="replace", form="unformatted")

    si = 1
    tstart = get_time()
    do ti=1,par%nt

       ! update supply
       if (t .gt. t_supply(si+1)) then
          si = si + 1
       end if

       ! log progress
       if ( mod(dble(ti), par%nt/10.d0) == 0 ) then
          call write_progress(ti, par%nt, tstart)
       end if

       do i=2,par%nx+1

          Cu(i) = max(0.d0, 1.5e-4 * ((u(ti) - u_th(i))**3) / (u(ti) * par%VS))
          
          Ct(i) = ((-par%VS * u(ti) * (Ct(i) - Ct(i-1)) / dx) * dt + \
                     Ct(i) + Cu(i) / Ta) / (1+1/Ta)

          if ( (Cu(i) - Ct(i)) / Ta .gt. Ca(i) ) then
             Ct(i) = ((-par%VS * u(ti) * (Ct(i) - Ct(i-1)) / dx) * dt + \
                        Ct(i) + Ca(i) / Ta) / (1+1/Ta)

             Ca(i) = Ca(i) - Ca(i) / Ta
          else
             Ca(i) = Ca(i) - (Cu(i) - Ct(i)) / Ta
          end if

          Ca(i) = Ca(i) + supply(si,i) / dx
 
       end do

       if ( mod(ti, nint(par%tout / dt)) == 0 ) then
          write(fid+1) Ct
          write(fid+2) Cu
          write(fid+3) Ca
       end if
       
       t = t + par%dt
    end do

    close(fid+1)
    close(fid+2)
    close(fid+3)

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
