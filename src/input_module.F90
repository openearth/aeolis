module input_module

  use constants_module

  implicit none

  include 'sedparams.inc'

  type parameters
     real*8    :: VS    = 0.d0
     real*8    :: Tp    = 0.d0
     real*8    :: u_th  = 0.d0
     real*8    :: z     = 0.d0
     real*8    :: S     = 0.d0
     real*8    :: phi   = 0.d0
     integer*4 :: nx    = 0.d0
     integer*4 :: nt    = 0.d0
     real*8    :: dt    = 0.d0
     real*8    :: dx    = 0.d0
     real*8    :: tstop = 0.d0
     real*8    :: tout  = 0.d0
     character(slen) :: wind_file = ''
     character(slen) :: bed_file = ''
     character(slen) :: supply_file = ''
  end type parameters

contains

  function read_params() result (fname)
    
    character(len=100) :: arg
    character(len=500) :: version
    character(slen)    :: fname
    integer            :: iarg, narguments
    logical            :: done
    
    fname = 'examples/input.txt'
    
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
             write(*,*)'    aeolis.exe'
             write(*,*)'    aeolis.exe [options]'
             write(*,*)' '
             write(*,*)'Options:'
             write(*,*)'    -v Shows the version of this AeoLiS executable'
             write(*,*)'**********************************************************'
             write(*,*)' '
             stop
          else
             fname = arg
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
    
  end function read_params
  
  function read_input(fname) result (par)
    
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
    par%phi   = read_key_dbl(fname, 'phi',   40.d0)
    par%dt    = read_key_dbl(fname, 'dt',    0.d05)
    par%dx    = read_key_dbl(fname, 'dx',    1.d0)
    par%tstop = read_key_dbl(fname, 'tstop', 3600.d0)
    par%tout  = read_key_dbl(fname, 'tout',  1.d0)
    par%wind_file   = read_key_str(fname, 'wind_file',  '')
    par%bed_file    = read_key_str(fname, 'bed_file',  '')
    par%supply_file = read_key_str(fname, 'supply_file',  '')
    
    write(*,*) '**********************************************************'
    write(*,*) ' '
    
  end function read_input

  function read_key_str(fname, key, default) result (value)
    
    integer :: ierr
    character(len=*) :: fname, key
    character(slen) :: value
    character(slen), optional :: default

    value = read_key(fname, key)
    if (value == ' ') then
       value = default
    end if

    write(0, '(a12,a,a)') key, ' = ', trim(value)
    
  end function read_key_str
  
  function read_key_dbl(fname, key, default) result (value_dbl)
    
    integer :: ierr
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
  
  function read_key_int(fname, key, default) result (value_int)
    
    integer :: ierr
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
    
    integer :: fid, ierr, idx
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
             return
          end if
       end if
    end do
    close(fid)
    
  end function read_key
  
end module input_module
