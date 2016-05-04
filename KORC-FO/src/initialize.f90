module initialize

use korc_types
use korc_hpc
use constants
use external_subroutines

implicit none

	PRIVATE :: set_paths, load_korc_params, initialization_sanity_check
	PUBLIC :: initialize_korc_parameters, initialize_particles, initialize_fields

contains

subroutine set_paths(params)
	implicit none
	INTEGER :: argn
	TYPE(KORC_PARAMS), INTENT(OUT) :: params

	argn = command_argument_count()
	call get_command_argument(1,params%path_to_inputs)
	call get_command_argument(2,params%path_to_outputs)
end subroutine set_paths


subroutine load_korc_params(params)
	implicit none
	TYPE (KORC_PARAMS), INTENT(INOUT) :: params
	LOGICAL :: restart ! Not used, yet.
	INTEGER(ip) :: t_steps
	REAL(rp) :: dt
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	INTEGER(ip) :: output_cadence
	INTEGER :: num_species
	INTEGER :: pic_algorithm

	NAMELIST /input_parameters/ magnetic_field_model,t_steps,dt,&
				output_cadence,num_species,pic_algorithm
	
	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=input_parameters)
	close(default_unit_open)

	! params%restart = restart
	params%t_steps = t_steps
	params%output_cadence = output_cadence
	params%num_snapshots = t_steps/output_cadence
	params%dt = dt
	params%num_species = num_species
	params%magnetic_field_model = TRIM(magnetic_field_model)
	params%pic_algorithm = pic_algorithm
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
	REAL(rp), DIMENSION(:), ALLOCATABLE :: ppp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: q
	REAL(rp), DIMENSION(:), ALLOCATABLE :: m
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Eo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Ro
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo
	LOGICAL, DIMENSION(:), ALLOCATABLE :: runaway
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Vo
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Xo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: angle, radius ! temporary vars
	REAL(rp), DIMENSION(:), ALLOCATABLE :: r ! temporary variable
	INTEGER :: ii ! Iterator

	NAMELIST /plasma_species/ ppp, q, m, Eo, runaway, Ro, Zo, r

	! Allocate array containing variables of particles for each species
	ALLOCATE(ptcls(params%num_species))

	ALLOCATE(ppp(params%num_species))
	ALLOCATE(q(params%num_species))
	ALLOCATE(m(params%num_species))
	ALLOCATE(Eo(params%num_species))
	ALLOCATE(runaway(params%num_species))
	ALLOCATE(Ro(params%num_species))
	ALLOCATE(Zo(params%num_species))

	ALLOCATE(r(params%num_species))

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=plasma_species)
	close(default_unit_open)

	do ii=1,params%num_species
		ptcls(ii)%Eo = Eo(ii)
		ptcls(ii)%runaway = runaway(ii)
		ptcls(ii)%q = q(ii)*C_E
		ptcls(ii)%m = m(ii)*C_ME
		ptcls(ii)%ppp = ppp(ii)
		ALLOCATE( ptcls(ii)%vars%X(3,ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%V(3,ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%Rgc(3,ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%Y(3,ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%E(3,ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%B(3,ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%gamma(ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%eta(ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%mu(ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%kappa(ptcls(ii)%ppp) )
		ALLOCATE( ptcls(ii)%vars%tau(ptcls(ii)%ppp) )

		ALLOCATE( Xo(3,ptcls(ii)%ppp) )
		ALLOCATE( Vo(ptcls(ii)%ppp) )
		
		ALLOCATE( angle(ptcls(ii)%ppp) )
		ALLOCATE( radius(ptcls(ii)%ppp) )

		! Initialize to zero
		ptcls(ii)%vars%X = 0.0_rp
		ptcls(ii)%vars%V = 0.0_rp
		ptcls(ii)%vars%Rgc = 0.0_rp
		ptcls(ii)%vars%Y = 0.0_rp
		ptcls(ii)%vars%E = 0.0_rp
		ptcls(ii)%vars%B = 0.0_rp
		ptcls(ii)%vars%gamma = 0.0_rp
		ptcls(ii)%vars%eta = 0.0_rp
		ptcls(ii)%vars%mu = 0.0_rp
		ptcls(ii)%vars%kappa = 0.0_rp
		ptcls(ii)%vars%tau = 0.0_rp

		! Initial condition of uniformly distributed particles on a disk in the xz-plane
		! A unique velocity direction
		call init_random_seed()
		call RANDOM_NUMBER(angle)
		angle = 2*C_PI*angle

		call init_random_seed()
		call RANDOM_NUMBER(radius)
		radius = r(ii)*radius
		
		Xo(1,:) = Ro(ii) + sqrt(radius)*cos(angle)
		Xo(2,:) = 0.0_rp
		Xo(3,:) = Zo(ii) + sqrt(radius)*sin(angle)

		ptcls(ii)%vars%X(1,:) = Xo(1,:)
		ptcls(ii)%vars%X(2,:) = Xo(2,:)
		ptcls(ii)%vars%X(3,:) = Xo(3,:)

		ptcls(ii)%vars%gamma(:) = ptcls(ii)%Eo*C_E/(ptcls(ii)%m*C_C**2)

		Vo = C_C*sqrt( 1 - 1/(ptcls(ii)%vars%gamma(:)**2) )

		ptcls(ii)%vars%V(1,:) = 0.0_rp
		ptcls(ii)%vars%V(2,:) = -Vo
		ptcls(ii)%vars%V(3,:) = 0.0_rp

		DEALLOCATE(angle)
		DEALLOCATE(radius)
		DEALLOCATE(Xo)
		DEALLOCATE(Vo)
	end do

	DEALLOCATE(ppp)
	DEALLOCATE(q)
	DEALLOCATE(m)
	DEALLOCATE(Eo)
	DEALLOCATE(runaway)
	DEALLOCATE(Ro)
	DEALLOCATE(Zo)

	DEALLOCATE(r)
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
	LOGICAL :: flag = .FALSE.

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'("* * * SANITY CHECK * * *")')
	end if

	call MPI_INITIALIZED(flag, ierr)

!$OMP PARALLEL
	write(6,'("MPI: ",I2," OMP thread: ", I2, " Initialized: ",l)') params%mpi_params%rank_topo, OMP_GET_THREAD_NUM(), flag
!$OMP END PARALLEL
end subroutine initialization_sanity_check


subroutine initialize_fields(params,EB)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(OUT) :: EB
	REAL(rp) :: Bo
	REAL(rp) :: minor_radius
	REAL(rp) :: major_radius
	REAL(rp) :: q_factor_at_separatrix
	REAL(rp) :: free_param

	NAMELIST /analytic_mag_field_params/ Bo,minor_radius,major_radius,&
			q_factor_at_separatrix,free_param

	if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
		read(default_unit_open,nml=analytic_mag_field_params)
		close(default_unit_open)

		EB%AB%Bo = Bo
		EB%AB%a = minor_radius
		EB%AB%Ro = major_radius
		EB%AB%qa = q_factor_at_separatrix
		EB%AB%co = free_param
		EB%AB%lambda = EB%AB%a / EB%AB%co
		EB%AB%Bpo = (EB%AB%a/EB%AB%Ro)*(EB%AB%Bo/EB%AB%qa)*(1+EB%AB%co**2)/EB%AB%co;

		EB%Bo = EB%AB%Bo
	else
		! Load data file containing magnetic field
	end if
end subroutine initialize_fields


end module initialize
