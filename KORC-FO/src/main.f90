program main

use constants
use korc_types
use main_mpi
use initialize
use finalize
use omp_lib

implicit none
TYPE (KORC_PARAMS) :: params

call initialize_communications(params)

call initialize_korc_parameters(params)

call finalize_communications(params)

end program main
