module input_module

  use constants_module
  use utils_module

  implicit none

  include 'sedparams.inc'

  type windstat
     real*8 :: t = 0.d0
     real*8 :: duration = 0.d0
     real*8 :: u_mean = 0.d0
     real*8 :: u_std = 0.d0
     real*8 :: gust_mean = 0.d0
     real*8 :: gust_std = 0.d0
     real*8 :: dir_mean = 0.d0
     real*8 :: dir_std = 0.d0
  end type windstat

  type windspeed
     real*8 :: t = 0.d0
     real*8 :: duration = 0.d0
     real*8 :: u = 0.d0
     real*8 :: dir = 0.d0
  end type windspeed

  type tide
     real*8 :: t = 0.d0
     real*8 :: level = 0.d0
  end type tide

  type moisture
     real*8 :: t = 0.d0
     real*8 :: z = 0.d0
     real*8 :: moisture = 0.d0
  end type moisture
  
  type meteorology
     real*8 :: t = 0.d0 ! [s]
     real*8 :: duration = 0.d0 ! [s]
     real*8 :: solar_radiation = 1e4 ! [J/m2]
     real*8 :: air_temperature = 10.d0 ! [oC]
     real*8 :: relative_humidity = 0.4d0 ! [-]
     real*8 :: air_specific_heat = 1.0035e-3 ! [MJ/kg/K]
     real*8 :: atmospheric_pressure = 101.325 ! [kPa]
     real*8 :: latent_heat = 2.45 ! [MJ/kg]
  end type meteorology
  
  type spaceparams
     real*8, dimension(:), pointer :: rho, dist
     real*8, dimension(:,:), pointer :: uw, udir, uws, uwn
     real*8, dimension(:,:), pointer :: zb, zs
     real*8, dimension(:,:), pointer :: xz, yz, xu, yu, xv, yv
     real*8, dimension(:,:,:), pointer :: uth, moist
     real*8, dimension(:,:,:), pointer :: Cu, Ct, supply, thlyr, p
     real*8, dimension(:,:,:), pointer :: d10, d50, d90
     real*8, dimension(:,:,:,:), pointer :: mass
     
     real*8, dimension(:,:), pointer :: dsz, dnz, dsdnzi, alfaz
     real*8, dimension(:,:), pointer :: dsu, dnu, dsdnui, alfau
     real*8, dimension(:,:), pointer :: dsv, dnv, dsdnvi, alfav
     real*8, dimension(:,:), pointer :: dsc, dnc
     type(meteorology) :: meteo
  end type spaceparams

  type spaceparams_linear
     real*8, dimension(:), pointer :: rho, dist
     real*8, dimension(:), pointer :: uw, udir, uws, uwn
     real*8, dimension(:), pointer :: zb, zs
     real*8, dimension(:), pointer :: xz, yz, xu, yu, xv, yv
     real*8, dimension(:,:), pointer :: uth, moist
     real*8, dimension(:,:), pointer :: Cu, Ct, supply, thlyr, p
     real*8, dimension(:,:), pointer :: d10, d50, d90
     real*8, dimension(:,:,:), pointer :: mass
  end type spaceparams_linear

  type parameters
     logical   :: mixtoplayer   = .true.
     logical   :: sweeptoplayer = .true.
     logical   :: th_grainsize  = .true.
     logical   :: th_bedslope   = .true.
     logical   :: th_moisture   = .true.
     logical   :: th_humidity   = .true.
     logical   :: bedupdate     = .true.
     logical   :: evaporation   = .true.
     logical   :: gusts         = .true.
     
     real*8    :: Tp    = 0.d0            ! [s] adaptation time scale in transport formulation
     real*8    :: u_th  = 0.d0            ! [m/s] constant velocity threshold in transport formulation
     real*8    :: z0    = 0.d0            ! [m] height of wind measurement
     real*8    :: k     = 0.d0            ! [m] bottom friction
     real*8    :: Cb    = 0.d0            ! [-] empirical constant in transport formulation
     real*8    :: phi   = 0.d0            ! [-] angle of repose of sediment
     real*8    :: g     = 0.d0            ! [m/s^2] gravitational acceleration
     integer*4 :: nx    = 0               ! [-] number of grid cells in x direction
     integer*4 :: ny    = 0               ! [-] number of grid cells in y direction
     integer*4 :: nc    = 0               ! [-] number of grid cells
     integer*4 :: nt    = 0               ! [-] number of time steps
     real*8    :: dt    = 0.d0            ! [s] duration of time step
     real*8    :: t     = 0.d0            ! [m] current time in simulation
     real*8    :: tstop = 0.d0            ! [s] duration of simulation
     real*8    :: tout  = 0.d0            ! [s] time interval for writing model output to disk
     integer*4 :: ntout = 0               ! [-] number of time steps written to disk
     real*8    :: CFL   = 0.d0            ! [-] Courant condition
     real*8    :: accfac = 0.d0           ! [-] acceleration factor applied to transport capacity

     character(slen) :: wind_file = ''    ! wind velocity time series definition file
     character(slen) :: xgrid_file = ''   ! x grid definition file
     character(slen) :: ygrid_file = ''   ! y grid definition file
     character(slen) :: bed_file = ''     ! bed profile definition file
     character(slen) :: tide_file = ''    ! tidal time series definition file
     character(slen) :: meteo_file = ''   ! meteo time series definition file
     character(slen) :: moist_file = ''   ! moisture time series definition file
     character(slen) :: method_moist = '' ! method for conversion of moisture content to velocity threshold
     character(slen) :: output_dir = ''   ! output directory relative to working directory

     character(20) :: scheme              ! numerical scheme (euler_forward/euler_backward/maccormack)
     integer*4 :: max_iter = 0            ! [-] maximum number of point iterations per time step in euler backward scheme
     real*8    :: max_error = 0d0         ! [-] maximum relative error allowed in time step in euler backward scheme

     integer*4 :: nfractions = 0          ! [-] number of sediment fractions
     integer*4 :: nlayers = 0             ! [-] number of bed layers additional to transport and base layer
     real*8    :: layer_thickness = 0.d0  ! [m] thickness of bed layers
     real*8    :: minfrac = 0.d0          ! [-] minimum fraction allowed in top layer of bed
     real*8    :: rhop = 0.d0             ! [kg/m^3] sediment density (excluding pores)
     real*8    :: rhom = 0.d0             ! [kg/m^3] sediment density (including pores)
     real*8    :: rhow = 0.d0             ! [kg/m^3] water density
     real*8    :: rhoa = 0.d0             ! [kg/m^3] air density
     real*8    :: porosity = 0.d0         ! [-] porosity of bed
     real*8    :: A = 0.d0                ! [-] empirical constant in grain size threshold formulation
     real*8    :: Hs = 0.d0               ! [m] offshore wave height
     real*8    :: gamma = 0.d0            ! [-] maximum wave height to depth ratio
     real*8    :: facDOD = 0.d0           ! [-] wave height to depth of disturbance ratio
     real*8    :: F = 0.d0                ! [m/s] hydraulic infiltration rate
     real*8    :: Cw = 0.d0               ! [-] sediment concentration in water column
     real*8    :: w = 0.d0                ! [m/s] fall velocity of sediment in water column
     real*8    :: bi = 0.d0               ! [-] bed interaction factor
     real*8    :: L = 0.d0                ! [-] length scale in sigmoid mapping

     real*8, dimension(:), allocatable :: grain_size ! [m] median grain size for each fraction
     real*8, dimension(:), allocatable :: grain_dist ! [-] occurence of each fraction in bed

     character(10), dimension(:), allocatable  :: outputvars  ! space separated list of output variables
     character(10), dimension(:), allocatable  :: outputtypes ! space separated list of output types (sum, avg, min, max, var)

     type(windspeed), dimension(:), allocatable :: uw ! [m/s] wind speed time series
     type(tide), dimension(:), allocatable :: zs ! [m] water level elevation
     type(meteorology), dimension(:), allocatable :: meteo ! meteorological conditions
     type(moisture), dimension(:,:), allocatable :: moist ! [-] moisture content

  end type parameters

contains

  function read_cmd() result (fname)
    
    character(len=100) :: arg
    character(len=500) :: version
    character(slen)    :: fpath, fname
    integer            :: iarg, narguments
    logical            :: done

    narguments = command_argument_count()
    if (narguments > 0) then
       do iarg=1,narguments
          
          call get_command_argument(iarg, arg)
          
          if (arg == '-v') then
             write(*,*)'**********************************************************'
             write(*,*)'You are using AeoLiS version 0.1'
             write(*,*)'**********************************************************'
             stop
          elseif (arg == '-h' .or. arg == '--help') then
             write(*,*)' '
             write(*,*)'**********************************************************'
             write(*,*)'                   Welcome to AeoLiS                      '
             write(*,*)' '
             write(*,*)'Usage:'
             write(*,*)'    aeolis.exe <input>'
             write(*,*)'    aeolis.exe <input> [options]'
             write(*,*)' '
             write(*,*)'Options:'
             write(*,*)'    -v Shows the version of this AeoLiS executable'
             write(*,*)'**********************************************************'
             write(*,*)' '
             stop
          else
             call split_path(arg, fpath, fname)
             call chdir(fpath)
          end if

       end do
    end if

    write(*,*) ' '
    write(*,*) '         d8888                   888      d8b  .d8888b.   ' 
    write(*,*) '        d88888                   888      Y8P d88P  Y88b  ' 
    write(*,*) '       d88P888                   888          Y88b.       ' 
    write(*,*) '      d88P 888  .d88b.   .d88b.  888      888  "Y888b.    ' 
    write(*,*) '     d88P  888 d8P  Y8b d88""88b 888      888     "Y88b.  ' 
    write(*,*) '    d88P   888 88888888 888  888 888      888       "888  ' 
    write(*,*) '   d8888888888 Y8b.     Y88..88P 888      888 Y88b  d88P  ' 
    write(*,*) '  d88P     888  "Y8888   "Y88P"  88888888 888  "Y8888P"   '
    write(*,*) ' '

    if (len(trim(fpath)) .gt. 0) then
       write(0, '(a, a)') ' Changed working directory to: ', trim(fpath)
       write(*,*) ' '
    end if

  end function read_cmd
  
  function read_params(fname) result (par)
    
    character(len=*) :: fname
    type(parameters) :: par

    write(*,*) '**********************************************************'
    write(*,*) 'PARAMETER SETTINGS'
    write(*,*) '**********************************************************'

    par%mixtoplayer   = read_key_logical(fname, 'mixtoplayer', .true.)
    par%sweeptoplayer = read_key_logical(fname, 'sweeptoplayer', .true.)
    par%th_grainsize  = read_key_logical(fname, 'th_grainsize', .true.)
    par%th_bedslope   = read_key_logical(fname, 'th_bedslope', .false.)
    par%th_moisture   = read_key_logical(fname, 'th_moisture', .true.)
    par%th_humidity   = read_key_logical(fname, 'th_humidity', .true.)
    par%bedupdate     = read_key_logical(fname, 'bedupdate', .true.)
    par%evaporation   = read_key_logical(fname, 'evaporation', .true.)
    par%gusts         = read_key_logical(fname, 'gusts', .true.)
    
    par%Tp     = read_key_dbl(fname, 'Tp',    1.0d0)
    par%u_th   = read_key_dbl(fname, 'u_th',  0.4d0)
    par%z0     = read_key_dbl(fname, 'z0',    1.d0)
    par%k      = read_key_dbl(fname, 'k',     0.01d0)
    par%Cb     = read_key_dbl(fname, 'Cb',    1.8d0)
    par%g      = read_key_dbl(fname, 'g',     9.81d0)
    par%phi    = read_key_dbl(fname, 'phi',   40.d0)
    par%dt     = read_key_dbl(fname, 'dt',    0.05d0)
    par%nx     = read_key_dbl(fname, 'nx',    1.d0)
    par%ny     = read_key_dbl(fname, 'ny',    1.d0)
    par%tstop  = read_key_dbl(fname, 'tstop', 3600.d0)
    par%tout   = read_key_dbl(fname, 'tout',  1.d0)
    par%accfac = read_key_dbl(fname, 'accfac',  1.d0)
    par%wind_file    = read_key_str(fname, 'wind_file', '')
    par%xgrid_file   = read_key_str(fname, 'xgrid_file', '')
    par%ygrid_file   = read_key_str(fname, 'ygrid_file', '')
    par%bed_file     = read_key_str(fname, 'bed_file', '')
    par%tide_file    = read_key_str(fname, 'tide_file', '')
    par%meteo_file   = read_key_str(fname, 'meteo_file', '')
    par%moist_file   = read_key_str(fname, 'moist_file', '')
    par%method_moist = read_key_str(fname, 'method_moist', 'belly_johnson')
    par%output_dir   = read_key_str(fname, 'output_dir', '')
    par%scheme       = read_key_str(fname, 'scheme', 'euler_backward')

    if (trim(par%scheme) .ne. 'euler_forward') then
       par%max_iter  = read_key_int(fname, 'max_iter',  100)
       par%max_error = read_key_dbl(fname, 'max_error', 1d-6)
    else
       par%CFL       = read_key_dbl(fname, 'CFL',      -1.d0)
    end if

    ! bed composition
    par%nfractions      = read_key_int(fname, 'nfractions',      1)
    par%nlayers         = read_key_int(fname, 'nlayers',         3)
    par%layer_thickness = read_key_dbl(fname, 'layer_thickness', 0.0005d0)
    par%minfrac         = read_key_dbl(fname, 'minfrac',         0.0001d0)
    par%rhoa            = read_key_dbl(fname, 'rhoa',            1.25d0)
    par%rhow            = read_key_dbl(fname, 'rhow',            1025.d0)
    par%rhom            = read_key_dbl(fname, 'rhom',            1650.d0)
    par%rhop            = read_key_dbl(fname, 'rhop',            2650.d0)
    par%porosity        = read_key_dbl(fname, 'porosity',        0.4d0)
    par%A               = read_key_dbl(fname, 'A',               100.d0)
    par%Hs              = read_key_dbl(fname, 'Hs',              1.d0)
    par%gamma           = read_key_dbl(fname, 'gamma',           0.5d0)
    par%facDOD          = read_key_dbl(fname, 'facDOD',          0.1d0)
    par%F               = read_key_dbl(fname, 'F',               1.0d-4)
    par%Cw              = read_key_dbl(fname, 'Cw',              0.d0)
    par%w               = read_key_dbl(fname, 'w',               3.0d-2)
    par%bi              = read_key_dbl(fname, 'bi',              1.d0)
    par%L               = read_key_dbl(fname, 'L',               1.d0)

    if (allocated(par%grain_size)) then
       deallocate(par%grain_size)
    end if
    allocate(par%grain_size(par%nfractions))
    
    if (allocated(par%grain_dist)) then
       deallocate(par%grain_dist)
    end if
    allocate(par%grain_dist(par%nfractions))

    par%grain_size = read_key_dblvec(fname, 'grain_size', par%nfractions, 0.d0003)
    par%grain_dist = read_key_dblvec(fname, 'grain_dist', par%nfractions, 1.d0)

    par%outputvars = read_key_strvec(fname, 'outputvars', 'z')
    par%outputtypes = read_key_strvec(fname, 'outputtypes', '')

    write(*,*) '**********************************************************'
    write(*,*) ' '

    call check_params(par)

  end function read_params

  subroutine check_params(par)

    type(parameters), intent(inout) :: par
    logical :: ex

    ! check required fields
    if (trim(par%xgrid_file) == '') then
       write(*, '(a)') " No grid defined"
       stop 1
    end if

    if (trim(par%bed_file) == '') then
       write(*, '(a)') " No bathymetry defined"
       stop 1
    end if

    ! check if valid scheme is selected
    select case (trim(par%scheme))
    case ('euler_backward', 'euler_forward', 'crank_nicolson')
       ! all ok
    case default
       write(*, '(a,a)') " Unsupported scheme: ",trim(par%scheme)
       stop 1
    end select

    ! check dimensionality
    if (par%ny == 1) then
       write(*, '(a)') "Warning: using ny=1, are you sure not to use ny=0?"
    end if
    
    ! sort grain size distribution
    call sort(par%grain_size, par%grain_dist)

    ! make time step fit with output time step
    par%dt = par%tout / ceiling(par%tout / par%dt)

    ! store total number of grid cells
    par%nc = (par%nx+1) * (par%ny+1)

    ! create output directory
    if (trim(par%output_dir) .ne. '') then
       par%output_dir = trim(par%output_dir) // "/"
       inquire(file=par%output_dir, exist=ex)
       if (.not. ex) then
          write(*, '(a,a)') " Created output directory ", trim(par%output_dir)
          call system("mkdir " // trim(par%output_dir))
       end if
    end if

    ! disable gusts for implicit schemes
    if (par%scheme .ne. 'euler_forward') then
       par%gusts = .false.
    end if

  end subroutine check_params
  
  function read_key_str(fname, key, default) result (value)
    
    integer*4 :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    character(*), optional :: default

    value = read_key(fname, key)
    if (value == ' ') then
       value = default
    end if

    write(*, '(a12,a,a)') key, ' = ', trim(value)

  end function read_key_str

  function read_key_strvec(fname, key, default) result (value_arr)
    
    integer*4 :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    integer*4 :: n, i
    character(*), optional :: default
    character(10), dimension(:), allocatable :: value_arr

    value = read_key(fname, key)

    if (value == ' ') then
       value = default
    end if

    value_arr = split(value)

    write(*, '(a12,a,a)') key, ' = ', trim(value_arr(1))
    
    do i = 2,size(value_arr)
       write(*, '(a15,a)') ' ', trim(value_arr(i))
    end do
        
  end function read_key_strvec
  
  function read_key_dbl(fname, key, default) result (value_dbl)

    integer*4 :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    real*8, optional :: default
    real*8 :: value_dbl

    value = read_key(fname, key)
    if (value /= ' ') then
       read(value, '(f10.0)', iostat=ierr) value_dbl
    else
       value_dbl = default
    end if

    write(*, '(a12,a,f15.4)') key, ' = ', value_dbl
    
  end function read_key_dbl
  
  function read_key_dblvec(fname, key, n, default) result (value_dbl)
    
    integer*4 :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    integer*4 :: n, i
    real*8, optional :: default
    real*8, dimension(n):: value_dbl

    value = read_key(fname, key)

    if (value /= ' ') then
       read(value, *, iostat=ierr) value_dbl(1:n)
    else
       value_dbl = default
    end if

    write(*, '(a12,a,f15.4)') key, ' = ', value_dbl(1)

    do i = 2,n
       write(*, '(a15,f15.4)') ' ', value_dbl(i)
    end do

  end function read_key_dblvec
  
  function read_key_int(fname, key, default) result (value_int)
    
    integer*4 :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    integer*4, optional :: default
    integer*4 :: value_int

    value = read_key(fname, key)
    if (value /= ' ') then
       read(value, '(i256)', iostat=ierr) value_int
    else
       value_int = default
    end if

    write(*, '(a12,a,i15)') key, ' = ', value_int
    
  end function read_key_int

  function read_key_logical(fname, key, default) result (value_logical)
    
    integer*4 :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    logical, optional :: default
    logical :: value_logical

    value = read_key(fname, key)
    if (value /= ' ') then
       read(value, '(l256)', iostat=ierr) value_logical
    else
       value_logical = default
    end if

    write(*, '(a12,a,l1)') key, ' = ', value_logical
    
  end function read_key_logical

  function read_key(fname, key) result (value)
    
    integer*4 :: fid, ierr, idx
    character(len=*) :: fname, key
    character(slen) :: value
    character(slen) :: line
    
    value = ''

    fid = 99
    ierr = 0
    open(fid, file=fname)
    do while (ierr == 0)
       read(fid, '(a)', iostat=ierr) line
       idx = scan(line, '=')
       if (idx > 0) then
          if (key == adjustl(line(1:idx-1))) then
             value = adjustl(line(idx+1:slen))
             exit
          end if
       end if
    end do
    close(fid)
    
  end function read_key
  
end module input_module
