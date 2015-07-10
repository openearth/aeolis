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

use constants_module
use input_module
use output_module
use bmi_module

implicit none

character(slen) :: configfile
character(kind=c_char) :: c_configfile(slen)
real*8 :: t, tstart, tend, tlog
integer*4 :: ierr

configfile = read_cmd()
c_configfile = string_to_char_array(configfile)

ierr = initialize(c_configfile)

par%t = 0
tstart = get_time()
tlog = tstart
call get_end_time(tend)
do while (par%t < tend)

   ! step in time
   ierr = update(-1.d0)

   ! write output
   if (par%t .le. par%dt  .or. par%tout < par%dt .or. &
        mod(par%t, par%tout) < par%dt) then

      call output_write(var)
      call output_clear(var)
   end if

   ! update time
   call get_current_time(par%t)

   ! log progress
   if ( mod(par%t, par%tstop/10.0) < par%dt .or. get_time()-tlog > 60) then
      call write_progress(par, tstart)
      tlog = get_time()
   end if

end do

call write_dimensions(par)

ierr = finalize()

end program
