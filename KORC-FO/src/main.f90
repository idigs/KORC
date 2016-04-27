program main

use korc_types
use main_mpi
use initialize
use finalize
use units
use fields

implicit none

	TYPE(KORC_PARAMS) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: ptcls
	TYPE(CHARCS_PARAMS) :: cp

	call initialize_communications(params)

	! INITIALIZATION STAGE
	call initialize_korc_parameters(params) ! Initialize korc parameters
	call initialize_particles(params,ptcls) ! Initialize particles

	call compute_charcs_plasma_params(ptcls,cp)
	! END OF INITIALIZATION STAGE


	! *** *** *** *** *** ***   *** *** *** *** *** *** ***
	! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
	! *** *** *** *** *** ***   *** *** *** *** *** *** ***



	! DEALLOCATION OF VARIABLES
	call deallocate_variables(params,ptcls)
	! DEALLOCATION OF VARIABLES

	call finalize_communications(params)

end program main
