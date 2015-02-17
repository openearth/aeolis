module input_module

  use constants_module
  use utils_module

  implicit none

  include 'sedparams.inc'

  type parameters
     real*8    :: VS    = 0.d0
     real*8    :: Tp    = 0.d0
     real*8    :: u_th  = 0.d0
     real*8    :: z     = 0.d0
     real*8    :: S     = 0.d0
     real*8    :: phi   = 0.d0
     real*8    :: g     = 9.81d0
     integer*4 :: nx    = 0
     integer*4 :: nt    = 0
     real*8    :: dt    = 0.d0
     real*8    :: dx    = 0.d0
     real*8    :: tstop = 0.d0
     real*8    :: tout  = 0.d0
     integer*4 :: ntout = 0
     real*8    :: CFL   = 0.d0
     real*8    :: accfac = 1.d0

     character(slen) :: wind_file = ''
     character(slen) :: bed_file = ''
     character(slen) :: moist_file = ''
     character(slen) :: method_moist = ''

     integer*4 :: nfractions = 1
     integer*4 :: nlayers = 3
     integer*4 :: nmix = 3
     real*8    :: layer_thickness = 5e-4
     real*8    :: rhop = 2650.d0
     real*8    :: rhom = 1650.d0
     real*8    :: rhow = 1025.d0
     real*8    :: rhoa = 1.25d0
     real*8    :: porosity = 0.4d0
     real*8    :: A = 100.d0
     real*8, dimension(:), allocatable :: grain_size
     real*8, dimension(:), allocatable :: grain_dist

     character(10) :: scheme

     character(10), dimension(:), allocatable  :: outputvars

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
       write(0, '(a, a)') '  Changed working directory to: ', trim(fpath)
       write(*,*) ' '
    end if

  end function read_cmd
  
  function read_params(fname) result (par)
    
    character(len=*) :: fname
    type(parameters) :: par

    write(*,*) '**********************************************************'
    write(*,*) 'PARAMETER SETTINGS'
    write(*,*) '**********************************************************'
    
    par%VS    = read_key_dbl(fname, 'VS',    1.d0)
    par%Tp    = read_key_dbl(fname, 'Tp',    0.d5)
    par%u_th  = read_key_dbl(fname, 'u_th',  0.d4)
    par%z     = read_key_dbl(fname, 'z',     0.d1)
    par%S     = read_key_dbl(fname, 'S',     0.d00015)
    par%g     = read_key_dbl(fname, 'g',     9.81d0)
    par%phi   = read_key_dbl(fname, 'phi',   40.d0)
    par%dt    = read_key_dbl(fname, 'dt',    0.d05)
    par%dx    = read_key_dbl(fname, 'dx',    1.d0)
    par%tstop = read_key_dbl(fname, 'tstop', 3600.d0)
    par%tout  = read_key_dbl(fname, 'tout',  1.d0)
    par%CFL   = read_key_dbl(fname, 'CFL',   -1.d0)
    par%accfac = read_key_dbl(fname, 'accfac',  1.d0)
    par%wind_file   = read_key_str(fname, 'wind_file',  '')
    par%bed_file    = read_key_str(fname, 'bed_file',   '')
    par%moist_file  = read_key_str(fname, 'moist_file', '')
    par%method_moist = read_key_str(fname, 'method_moist', 'belly_johnson')
    par%scheme = read_key_str(fname, 'scheme', 'implicit')

    ! bed composition
    par%nfractions      = read_key_int(fname, 'nfractions',      1)
    par%nlayers         = read_key_int(fname, 'nlayers',         3)
    par%nmix            = read_key_int(fname, 'nmix',            2)
    par%layer_thickness = read_key_dbl(fname, 'layer_thickness', 0.d0005)
    par%rhoa            = read_key_dbl(fname, 'rhoa',            1.25d0)
    par%rhow            = read_key_dbl(fname, 'rhow',            1025.d0)
    par%rhom            = read_key_dbl(fname, 'rhom',            1650.d0)
    par%rhop            = read_key_dbl(fname, 'rhop',            2650.d0)
    par%porosity        = read_key_dbl(fname, 'porosity',        0.4d0)
    par%A               = read_key_dbl(fname, 'A',               100.d0)

    allocate(par%grain_size(par%nfractions))
    allocate(par%grain_dist(par%nfractions))

    par%grain_size = read_key_dblvec(fname, 'grain_size', par%nfractions, 0.d0003)
    par%grain_dist = read_key_dblvec(fname, 'grain_dist', par%nfractions, 1.d0)

    par%outputvars = read_key_strvec(fname, 'outputvars', 'z')

    write(*,*) '**********************************************************'
    write(*,*) ' '

    call check_params(par)

  end function read_params

  subroutine check_params(par)

    type(parameters), intent(inout) :: par

    ! check if valid scheme is selected
    select case (trim(par%scheme))
    case ('implicit', 'explicit')
       ! all ok
    case default
       write(0, '(a,a)') " Unsupported scheme: ",trim(par%scheme)
       stop 1
    end select

    ! sort grain size distribution
    call sort(par%grain_size, par%grain_dist)

    ! number of mix layers cannot exceed available number of layers
    par%nmix = min(par%nmix, par%nlayers+2)

    ! make time step fit with output time step
    par%dt = par%tout / ceiling(par%tout / par%dt)

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

    write(0, '(a12,a,a)') key, ' = ', trim(value)

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

    write(0, '(a12,a,a)') key, ' = ', trim(value_arr(1))
    
    do i = 2,size(value_arr)
       write(0, '(a15,a)') ' ', trim(value_arr(i))
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

    write(0, '(a12,a,f15.4)') key, ' = ', value_dbl
    
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

    write(0, '(a12,a,f15.4)') key, ' = ', value_dbl(1)

    do i = 2,n
       write(0, '(a15,f15.4)') ' ', value_dbl(i)
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

    write(0, '(a12,a,i15)') key, ' = ', value_int
    
  end function read_key_int
  
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
