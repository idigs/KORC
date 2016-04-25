program main

use constants
use korc_types
use main_mpi
use initialize
use omp_lib

implicit none

	CHARACTER(15) :: path_to_input = 'input_file.korc'
	TYPE (KORC_PARAMS) :: params

	write(6,*) path_to_input
	call load_korc_params(path_to_input,params)
	write(6,*) params

end program main
