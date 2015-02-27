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
use run_module

implicit none

character(slen) :: fname
type(parameters) :: par

fname = read_cmd()
par   = read_params(fname)

call run_model(par)

end program
