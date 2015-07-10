module logging

  use iso_c_utils
  use input_module
  
  implicit none

  abstract interface
     subroutine ilogger(level, msg)
       use iso_c_binding
       use iso_c_utils
       integer(c_int), value, intent(in) :: level !< severity
       character(c_char), intent(in) :: msg(MAXSTRINGLEN) !< c message null terminated
     end subroutine ilogger
  end interface

  procedure(ilogger), pointer :: logging_callback => null()

  ! Levels correspond to log4net and log4j
  integer, parameter, public :: LEVEL_ALL = 0
  integer, parameter, public :: LEVEL_DEBUG = 1
  integer, parameter, public :: LEVEL_INFO  = 2
  integer, parameter, public :: LEVEL_WARN  = 3
  integer, parameter, public :: LEVEL_ERROR = 4
  integer, parameter, public :: LEVEL_FATAL = 5
  integer, parameter, public :: LEVEL_OFF = 6

  ! message buffer that can be used to fill log message
  character(len=MAXSTRINGLEN) :: msgbuf

contains
  
  subroutine set_logger(c_callback) bind(C, name="set_logger")
    !DEC$ ATTRIBUTES DLLEXPORT::set_logger

    type(c_funptr), value :: c_callback

    ! set a callback that will be cauled with new messages
    call c_f_procpointer(c_callback, logging_callback)
    
  end subroutine set_logger

  subroutine log(level, msg)
    integer(c_int), intent(in) :: level
    character(len=*), intent(in) :: msg
    
    character(c_char) :: c_string(MAXSTRINGLEN)

    if (associated(logging_callback)) then
       c_string = string_to_char_array(msg)
       call logging_callback(level, c_string)
    end if
  end subroutine log

  subroutine write_progress(par, tstart)

    type(parameters), intent(in) :: par
    real*8, intent(in) :: tstart
    real*8 :: p, dt1, dt2, dt3

    p = par%t / par%tstop
    dt1 = get_time() - tstart
    dt2 = dt1 / p
    dt3 = dt2 * (1 - p)

    write(*,'(f5.1,a2,a8,a3,a8,a3,a8,a10,f4.2,a1)') &
         100.d0*p, '% ', &
         format_time(dt1), ' / ', &
         format_time(dt2), ' / ', &
         format_time(dt3), ' (avg. dt=', (par%t / par%nt), ')'

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

end module logging
