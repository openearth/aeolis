program aeolis

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
