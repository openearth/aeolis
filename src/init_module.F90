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
     real*8, dimension(:,:), pointer :: uth, moist
     real*8, dimension(:,:), pointer :: Cu, Ct, p, supply, thlyr
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
    call alloc_pointer(var, 'x', (/par%nx+1/), s%x)
    call alloc_pointer(var, 'zb', (/par%nx+1/), s%zb)
    
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
    write(0,*) 'Generating wind time series...'
    call generate_wind(par, par%uw)
    open(unit=fid, file=trim(par%output_dir) // "wind.in", &
         action="write", status="replace", form="unformatted")
    write(fid) par%uw%t
    write(fid) par%uw%u
    close(fid)

    ! time
    par%ntout = int(par%tstop / par%tout) + 1

    ! moist
    write(0,*) 'Generating moisture time series...'
    call generate_moist(par, par%moist)

    ! tide and meteo
    write(0,*) 'Generating tide and meteo time series...'
    call generate_tide(par, par%zs)
    call generate_meteo(par, par%meteo)

    ! fractions
    call alloc_pointer(var, 'rho',    (/par%nfractions/), s%rho)
    call alloc_pointer(var, 'dist',   (/par%nfractions/), s%dist)

    s%rho = par%rhom
    s%dist = par%grain_dist

    ! variables
    call alloc_pointer(var, 'Cu',     (/par%nfractions, par%nx+1/), s%Cu)
    call alloc_pointer(var, 'Ct',     (/par%nfractions, par%nx+1/), s%Ct)
    call alloc_pointer(var, 'p',      (/par%nfractions, par%nx+1/), s%p)
    call alloc_pointer(var, 'uth',    (/par%nfractions, par%nx+1/), s%uth)
    call alloc_pointer(var, 'mass',   (/par%nfractions, par%nlayers, par%nx+1/), s%mass)
    call alloc_pointer(var, 'supply', (/par%nfractions, par%nx+1/), s%supply)
    call alloc_pointer(var, 'moist',  (/par%nlayers, par%nx+1/), s%moist)
    call alloc_pointer(var, 'thlyr',  (/par%nlayers, par%nx+1/), s%thlyr)

    ! extra output
    call alloc_pointer(var, 'uw', (/0/), s%uw)
    call alloc_pointer(var, 'zs', (/0/), s%zs)
    call alloc_pointer(var, 'd10', (/par%nlayers, par%nx+1/), s%d10)
    call alloc_pointer(var, 'd50', (/par%nlayers, par%nx+1/), s%d50)
    call alloc_pointer(var, 'd90', (/par%nlayers, par%nx+1/), s%d90)

  end subroutine init
  
end module init_module
