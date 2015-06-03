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
use output_module
use bmi_module

implicit none

character(slen) :: configfile
character(kind=c_char) :: c_configfile(slen)
real*8 :: t, tstart, tend, tlog

configfile = read_cmd()
c_configfile = string_to_char_array(configfile)

if (initialize(c_configfile) /= 0) &
     write(msgbuf,*) 'Initialization failed'

! initialize output
call output_init(par, var)
call write_output(par, sl, var)

t = 0
tstart = get_time()
tlog = tstart
call get_end_time(tend)
do while (t < tend)

   ! step in time
   if (update(-1.d0) /= 0) &
        write(msgbuf,*) 'Updating to timestep ', t, ' failed'

   ! write output
   if (par%t .le. par%dt  .or. par%tout < par%dt .or. &
        mod(par%t, par%tout) < par%dt .or. par%tout - par%t < par%dt) then
      call write_output(par, sl, var)
   end if

   ! update time
   call get_current_time(t)

   ! log progress
   if ( mod(par%t, par%tstop/10.0) < par%dt .or. get_time()-tlog > 60) then
      call write_progress(par, tstart)
      tlog = get_time()
   end if

end do

!call write_dimensions(par)
call output_close(var)
if (finalize() /= 0) &
     write(msgbuf,*) 'Finalization failed'

end program
