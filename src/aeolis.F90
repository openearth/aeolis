program aeolis

use constants_module
use input_module
use run_module

implicit none

character(slen) :: fname
type(parameters) :: par

fname = read_params()
par   = read_input(fname)

call run_model(par)

end program
