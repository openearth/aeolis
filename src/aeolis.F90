program aeolis
  ! AeoLiS main program.
  ! Reads command line options, model parameters and start model run.
  !
  ! Usage:
  !     | aeolis <input>
  !     | aeolis <input> [options]
  !
  ! Options:
  !     -v Shows the version of this AeoLiS executable

use logging
use constants_module
use input_module
use bmi_module

implicit none

character(slen) :: configfile
character(kind=c_char) :: c_configfile(slen)
real*8 :: t, tstart, tend, tlog
integer*4 :: ti

configfile = read_cmd()
c_configfile = string_to_char_array(configfile)

if (initialize(c_configfile) /= 0) &
     write(msgbuf,*) 'Initialization failed'

t = 0
ti = 1
tstart = get_time()
tlog = tstart
call get_end_time(tend)
do while (t < tend)

   ! log progress
   if ( mod(dble(ti), par%nt/10.d0) < 1.d0 .or. get_time()-tlog > 60) then
      call write_progress(ti, par%nt, tstart)
      tlog = get_time()
   end if
   
   if (update(-1.d0) /= 0) &
        write(msgbuf,*) 'Updating to timestep ', t, ' failed'

   call write_output(par, s, var)
   call get_current_time(t)

   ti = ti + 1
   
end do

if (finalize() /= 0) &
     write(msgbuf,*) 'Finalization failed'

end program
