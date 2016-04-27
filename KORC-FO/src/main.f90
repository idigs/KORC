program main

use korc_types
use main_mpi
use initialize
use finalize
use units
use emf

implicit none

	TYPE(KORC_PARAMS) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: ptcls
	TYPE(CHARCS_PARAMS) :: cp
	TYPE(FIELDS) :: EB

	call initialize_communications(params)

	! INITIALIZATION STAGE
	call initialize_korc_parameters(params) ! Initialize korc parameters

	call initialize_particles(params,ptcls) ! Initialize particles

	call initialize_fields(params,EB)

	call compute_charcs_plasma_params(ptcls,EB,cp)
	! END OF INITIALIZATION STAGE

!	write(6,*) params%num_snapshots


	! *** *** *** *** *** ***   *** *** *** *** *** *** ***
	! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
	! *** *** *** *** *** ***   *** *** *** *** *** *** ***



	! DEALLOCATION OF VARIABLES
	call deallocate_variables(params,ptcls)
	! DEALLOCATION OF VARIABLES

	call finalize_communications(params)

end program main
