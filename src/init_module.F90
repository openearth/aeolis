module init_module

  use output_module
  use bed_module
  use wind_module
  use moist_module
  
  implicit none

  type spaceparams
     real*8, pointer :: uw, zs
     real*8, dimension(:), pointer :: rho, dist
     real*8, dimension(:,:), pointer :: x, y, zb
     real*8, dimension(:,:,:), pointer :: uth, moist
     real*8, dimension(:,:,:), pointer :: Cu, Ct, supply, thlyr
     real*8, dimension(:,:,:), pointer :: d10, d50, d90
     real*8, dimension(:,:,:,:), pointer :: mass
     type(meteorology) :: meteo
  end type spaceparams

  type spaceparams_linear
     real*8, pointer :: uw, zs
     real*8, dimension(:), pointer :: rho, dist
     real*8, dimension(:), pointer :: x, y, zb
     real*8, dimension(:,:), pointer :: uth, moist
     real*8, dimension(:,:), pointer :: Cu, Ct, p, supply, thlyr
     real*8, dimension(:,:), pointer :: d10, d50, d90
     real*8, dimension(:,:,:), pointer :: mass
  end type spaceparams_linear

contains

  subroutine init(par, var, s, sl)

    type(parameters), intent(inout) :: par
    type(spaceparams), intent(out) :: s
    type(spaceparams_linear), intent(out) :: sl
    type(variables), dimension(:), allocatable, intent(out) :: var
    real*8, dimension(:,:), allocatable :: x_tmp, y_tmp, zb_tmp
    integer*4, parameter :: fid=20

    write(0,*) 'Initialization started...'
    
    ! bed
    write(0,*) 'Generating bed profile and composition...'
    call generate_bed(par, x_tmp, y_tmp, zb_tmp)

    call alloc_variable(var, 'x',  (/par%ny+1, par%nx+1/))
    call alloc_variable(var, 'y',  (/par%ny+1, par%nx+1/))
    call alloc_variable(var, 'zb', (/par%ny+1, par%nx+1/))

    call get_pointer(var, 'x',  (/par%ny+1, par%nx+1/), s%x)
    call get_pointer(var, 'y',  (/par%ny+1, par%nx+1/), s%y)
    call get_pointer(var, 'zb', (/par%ny+1, par%nx+1/), s%zb)
    
    call get_pointer(var, 'x',  (/par%nc/), sl%x)
    call get_pointer(var, 'y',  (/par%nc/), sl%y)
    call get_pointer(var, 'zb', (/par%nc/), sl%zb)

    s%x = x_tmp
    s%y = y_tmp
    s%zb = zb_tmp
    deallocate(x_tmp)
    deallocate(y_tmp)
    deallocate(zb_tmp)

    call generate_bedcomposition(par)
    open(unit=fid, file=trim(par%output_dir) // "bed.in", &
         action="write", status="replace", form="unformatted")
    write(fid) s%x
    write(fid) s%y
    write(fid) s%zb
    close(fid)

    ! wind
    write(*,*) 'Generating wind time series...'
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
    call alloc_variable(var, 'rho',    (/par%nfractions/))
    call alloc_variable(var, 'dist',   (/par%nfractions/))

    ! variables
    call alloc_pointer(var, 'Cu',     (/par%nfractions, par%ny+1, par%nx+1/), s%Cu)
    call alloc_pointer(var, 'Ct',     (/par%nfractions, par%ny+1, par%nx+1/), s%Ct)
    call alloc_pointer(var, 'p',      (/par%nfractions, par%ny+1, par%nx+1/), s%p)
    call alloc_pointer(var, 'uth',    (/par%nfractions, par%ny+1, par%nx+1/), s%uth)
    call alloc_pointer(var, 'mass',   (/par%nfractions, par%ny+1, par%nlayers, par%nx+1/), s%mass)
    call alloc_pointer(var, 'supply', (/par%nfractions, par%ny+1, par%nx+1/), s%supply)
    call alloc_pointer(var, 'moist',  (/par%nlayers, par%ny+1, par%nx+1/), s%moist)
    call alloc_pointer(var, 'thlyr',  (/par%nlayers, par%ny+1, par%nx+1/), s%thlyr)

    ! extra output
    call alloc_variable(var, 'uw',     (/0/))
    call alloc_variable(var, 'zs',     (/0/))
    call alloc_variable(var, 'd10',    (/par%nlayers, par%ny+1, par%nx+1/))
    call alloc_variable(var, 'd50',    (/par%nlayers, par%ny+1, par%nx+1/))
    call alloc_variable(var, 'd90',    (/par%nlayers, par%ny+1, par%nx+1/))

    ! create remapping pointers
    call get_pointer(var, 'rho',    (/par%nfractions/), s%rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), s%dist)
    call get_pointer(var, 'Cu',     (/par%nfractions, par%ny+1, par%nx+1/), s%Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%ny+1, par%nx+1/), s%Ct)
    call get_pointer(var, 'uth',    (/par%nfractions, par%ny+1, par%nx+1/), s%uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers, par%ny+1, par%nx+1/), s%mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%ny+1, par%nx+1/), s%supply)
    call get_pointer(var, 'moist',  (/par%nlayers, par%ny+1, par%nx+1/), s%moist)
    call get_pointer(var, 'thlyr',  (/par%nlayers, par%ny+1, par%nx+1/), s%thlyr)
    call get_pointer(var, 'uw',     (/0/), s%uw)
    call get_pointer(var, 'zs',     (/0/), s%zs)
    call get_pointer(var, 'd10',    (/par%nlayers, par%ny+1, par%nx+1/), s%d10)
    call get_pointer(var, 'd50',    (/par%nlayers, par%ny+1, par%nx+1/), s%d50)
    call get_pointer(var, 'd90',    (/par%nlayers, par%ny+1, par%nx+1/), s%d90)

    call get_pointer(var, 'rho',    (/par%nfractions/), sl%rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), sl%dist)
    call get_pointer(var, 'Cu',     (/par%nfractions, par%nc/), sl%Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%nc/), sl%Ct)
    call get_pointer(var, 'uth',    (/par%nfractions, par%nc/), sl%uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers, par%nc/), sl%mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%nc/), sl%supply)
    call get_pointer(var, 'moist',  (/par%nlayers, par%nc/), sl%moist)
    call get_pointer(var, 'thlyr',  (/par%nlayers, par%nc/), sl%thlyr)
    call get_pointer(var, 'uw',     (/0/), sl%uw)
    call get_pointer(var, 'zs',     (/0/), sl%zs)
    call get_pointer(var, 'd10',    (/par%nlayers, par%nc/), sl%d10)
    call get_pointer(var, 'd50',    (/par%nlayers, par%nc/), sl%d50)
    call get_pointer(var, 'd90',    (/par%nlayers, par%nc/), sl%d90)

    s%rho = par%rhom
    s%dist = par%grain_dist

  end subroutine init
  
end module init_module
