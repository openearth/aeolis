module init_module

  use output_module
  use bed_module
  use wind_module
  use moist_module
  
  implicit none

  type spaceparams
     real*8, pointer :: uw, zs
     real*8, dimension(:), pointer :: x, zb
     real*8, dimension(:), pointer :: rho, dist
     real*8, dimension(:,:), pointer :: uth
     real*8, dimension(:,:), pointer :: Cu, Ct, supply, moist, thlyr
     real*8, dimension(:,:), pointer :: d10, d50, d90
     real*8, dimension(:,:,:), pointer :: mass
     type(meteorology) :: meteo
  end type spaceparams
  
contains

  subroutine init(par, var, s)

    type(parameters), intent(inout) :: par
    type(spaceparams), intent(out) :: s
    type(variables), dimension(:), allocatable, intent(out) :: var
    real*8, dimension(:), allocatable :: x_tmp, zb_tmp
    integer*4, parameter :: fid=20
    
    write(0,*) 'Initialization started...'
    
    ! bed
    write(0,*) 'Generating bed profile and composition...'
    call generate_bed(par, x_tmp, zb_tmp)
    call get_pointer(var, 'x', (/par%nx+1/), s%x)
    call get_pointer(var, 'zb', (/par%nx+1/), s%zb)
    
    s%x = x_tmp
    s%zb = zb_tmp
    deallocate(x_tmp)
    deallocate(zb_tmp)

    call generate_bedcomposition(par, s%x, s%zb)
    open(unit=fid, file=trim(par%output_dir) // "bed.in", &
         action="write", status="replace", form="unformatted")
    write(fid) s%x
    write(fid) s%zb
    close(fid)

    ! wind
    write(*,*) 'Generating bed wind time series...'
    call generate_wind(par, par%uw)
    open(unit=fid, file=trim(par%output_dir) // "wind.in", &
         action="write", status="replace", form="unformatted")
    write(fid) par%uw
    close(fid)

    ! courant check
    if (trim(par%scheme) .eq. 'explicit') then
       write(0, '(a, f4.2)') " Courant condition: ", maxval(par%uw) / par%dx * par%dt   
       if (par%dx / par%dt < maxval(par%uw)) then
          write(0, '(a)') " Courant condition violated. Please adapt numerical parameters."
          stop 1
       end if
    end if
    
    ! time
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
    call generate_tide(par, par%zs)
    call generate_meteo(par, par%meteo)

    ! fractions
    call get_pointer(var, 'rho',    (/par%nfractions/), s%rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), s%dist)

    s%rho = par%rhom
    s%dist = par%grain_dist

    ! variables
    call get_pointer(var, 'Cu',     (/par%nfractions, par%nx+1/), s%Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%nx+1/), s%Ct)
    call get_pointer(var, 'uth',    (/par%nfractions, par%nx+1/), s%uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers, par%nx+1/), s%mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%nx+1/), s%supply)
    call get_pointer(var, 'moist',  (/par%nlayers, par%nx+1/), s%moist)
    call get_pointer(var, 'thlyr',  (/par%nlayers, par%nx+1/), s%thlyr)

    ! extra output
    call get_pointer(var, 'uw', (/0/), s%uw)
    call get_pointer(var, 'zs', (/0/), s%zs)
    call get_pointer(var, 'd10', (/par%nlayers, par%nx+1/), s%d10)
    call get_pointer(var, 'd50', (/par%nlayers, par%nx+1/), s%d50)
    call get_pointer(var, 'd90', (/par%nlayers, par%nx+1/), s%d90)

  end subroutine init
  
end module init_module
