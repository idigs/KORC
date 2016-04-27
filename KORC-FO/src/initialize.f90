module initialize
use korc_types
use main_mpi
use omp_lib
implicit none

INTEGER, PRIVATE :: str_length
CHARACTER(MAX_STRING_LENGTH), PRIVATE :: aux_str
contains


subroutine set_paths(params)
	implicit none
	INTEGER :: argn
	TYPE(KORC_PARAMS), INTENT(OUT) :: params

	argn = command_argument_count()
	call get_command_argument(1,params%path_to_inputs)
	call get_command_argument(2,params%path_to_outputs)
	! write(6,*) argn
	! write(6,*) TRIM(params%path_to_inputs), LEN(params%path_to_inputs)
	! write(6,*) TRIM(params%path_to_outputs), LEN(params%path_to_outputs)
end subroutine set_paths


subroutine load_korc_params(params)
	implicit none
	TYPE (KORC_PARAMS), INTENT(INOUT) :: params
	LOGICAL :: restart ! Not used, yet.
	INTEGER :: t_steps
	REAL(rp) :: DT
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	INTEGER :: output_cadence
	INTEGER :: num_species

	NAMELIST /input_parameters/ magnetic_field_model,t_steps,DT,&
				output_cadence,num_species
	
	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=input_parameters)
	close(default_unit_open)

	! params%restart = restart
	params%t_steps = t_steps
	params%DT = DT
	params%output_cadence = output_cadence
	params%num_species = num_species
	params%magnetic_field_model = TRIM(magnetic_field_model)	
end subroutine load_korc_params


subroutine initialize_korc_parameters(params)
	use korc_types
	implicit none
	TYPE(KORC_PARAMS), INTENT(OUT) :: params

	call set_paths(params)
	call load_korc_params(params)
end subroutine initialize_korc_parameters


subroutine initialize_particles(params,ptcls) 
	use korc_types
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ptcls
	REAL(rp), DIMENSION(:), ALLOCATABLE :: ppp, q, m
	INTEGER :: ii ! Iterator

	NAMELIST /plasma_species/ ppp, q, m

	! Allocate array containing variables of particles for each species
	ALLOCATE(ptcls(params%num_species))

	ALLOCATE(ppp(params%num_species))
	ALLOCATE(q(params%num_species))
	ALLOCATE(m(params%num_species))

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=plasma_species)
	close(default_unit_open)

	do ii=1,params%num_species
		ptcls(ii)%q = q(ii)
		ptcls(ii)%m = m(ii)
		ptcls(ii)%ppp = ppp(ii)
		ALLOCATE(ptcls(ii)%vars%X(3,ptcls(ii)%ppp))
		ALLOCATE(ptcls(ii)%vars%V(3,ptcls(ii)%ppp))
		ALLOCATE(ptcls(ii)%vars%gamma(ptcls(ii)%ppp))
		ALLOCATE(ptcls(ii)%vars%eta(ptcls(ii)%ppp))
	end do

	DEALLOCATE(q)
	DEALLOCATE(m)
	DEALLOCATE(ppp)
end subroutine initialize_particles


subroutine initialize_communications(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

!$OMP PARALLEL
	params%num_omp_threads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

	call initialize_mpi(params)

	call initialization_sanity_check(params) 
end subroutine initialize_communications


subroutine initialization_sanity_check(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	INTEGER :: ierr

	if (params%mpi_params%rank_topo .EQ. 0) then
		write(6,'("* * * SANITY CHECK * * *")')
	end if

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!$OMP PARALLEL
	write(6,'("MPI: ",I2," OMP thread: ", I2)') params%mpi_params%rank_topo, OMP_GET_THREAD_NUM()
!$OMP END PARALLEL

end subroutine initialization_sanity_check


end module initialize
