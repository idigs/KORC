module initialize

    use korc_types
    use constants
    use korc_hpc
    use korc_HDF5
    use korc_interp
    use rnd_numbers

    implicit none

	

	PRIVATE :: set_paths, load_korc_params, initialization_sanity_check, unitVectors
	PUBLIC :: initialize_korc_parameters, initialize_particles, initialize_fields

    contains

subroutine set_paths(params)
	implicit none
	INTEGER :: argn
	TYPE(KORC_PARAMS), INTENT(OUT) :: params

	argn = command_argument_count()
	call get_command_argument(1,params%path_to_inputs)
	call get_command_argument(2,params%path_to_outputs)

	write(6,'("* * * * * PATHS * * * * *")')
	write(6,'("The input file is:",A50)') TRIM(params%path_to_inputs)
	write(6,'("The output folder is:",A50)') TRIM(params%path_to_outputs)
end subroutine set_paths


subroutine load_korc_params(params)
	implicit none
	TYPE (KORC_PARAMS), INTENT(INOUT) :: params
	LOGICAL :: restart ! Not used, yet.
	INTEGER(ip) :: t_steps
	REAL(rp) :: dt
	LOGICAL :: radiation
	LOGICAL :: collisions
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	LOGICAL :: poloidal_flux
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename
	INTEGER(ip) :: output_cadence
	INTEGER :: num_species
	INTEGER :: num_impurity_species
	INTEGER :: pic_algorithm

	NAMELIST /input_parameters/ magnetic_field_model,poloidal_flux,&
			magnetic_field_filename,t_steps,dt,output_cadence,num_species,&
			pic_algorithm,radiation,collisions,num_impurity_species
	
	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=input_parameters)
	close(default_unit_open)

	! params%restart = restart
	params%t_steps = t_steps
	params%output_cadence = output_cadence
	params%num_snapshots = t_steps/output_cadence
	params%dt = dt
	params%num_species = num_species
	params%num_impurity_species = num_impurity_species
	params%magnetic_field_model = TRIM(magnetic_field_model)
	params%poloidal_flux = poloidal_flux
	params%magnetic_field_filename = TRIM(magnetic_field_filename)
	params%pic_algorithm = pic_algorithm
	params%radiation = radiation
	params%collisions = collisions

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'("* * * * * SIMULATION PARAMETERS * * * * *")')
		write(6,'("Number of time steps: ",I16)') params%t_steps
		write(6,'("Output cadence: ",I16)') params%output_cadence
		write(6,'("Number of outputs: ",I16)') params%num_snapshots
		write(6,'("Time step in fraction of relativistic gyro-period: ",F15.10)') params%dt
		write(6,'("Number of electron populations: ",I16)') params%num_species
		write(6,'("Magnetic field model: ",A50)') TRIM(params%magnetic_field_model)
		write(6,'("Using (JFIT) poloidal flux: ", L1)') params%poloidal_flux
		write(6,'("Magnetic field model: ",A100)') TRIM(params%magnetic_field_filename)
		write(6,'("Radiation losses included: ",L1)') params%radiation
		write(6,'("collisions losses included: ",L1)') params%collisions
	end if	
end subroutine load_korc_params


subroutine initialize_korc_parameters(params)
	use korc_types
	implicit none
	TYPE(KORC_PARAMS), INTENT(OUT) :: params

	call set_paths(params)
	call load_korc_params(params)
end subroutine initialize_korc_parameters


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_particles(params,F,spp) 
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: ppp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: q
	REAL(rp), DIMENSION(:), ALLOCATABLE :: m
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Eo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: etao
	LOGICAL, DIMENSION(:), ALLOCATABLE :: runaway
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Ro
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: r
	INTEGER :: ii,jj ! Iterator

	NAMELIST /plasma_species/ ppp, q, m, Eo, etao, runaway, Ro, Zo, r

	! Allocate array containing variables of particles for each species
	ALLOCATE(spp(params%num_species))

	ALLOCATE(ppp(params%num_species))
	ALLOCATE(q(params%num_species))
	ALLOCATE(m(params%num_species))
	ALLOCATE(Eo(params%num_species))
	ALLOCATE(etao(params%num_species))
	ALLOCATE(runaway(params%num_species))
	ALLOCATE(Ro(params%num_species))
	ALLOCATE(Zo(params%num_species))

	ALLOCATE(r(params%num_species))

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=plasma_species)
	close(default_unit_open)

	do ii=1,params%num_species
		spp(ii)%Eo = Eo(ii)*C_E
		spp(ii)%etao = etao(ii)
		spp(ii)%runaway = runaway(ii)
		spp(ii)%q = q(ii)*C_E
		spp(ii)%m = m(ii)*C_ME
		spp(ii)%ppp = ppp(ii)

		spp(ii)%Ro = Ro(ii)
		spp(ii)%Zo = Zo(ii)
		spp(ii)%r = r(ii)

		spp(ii)%gammao = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

		ALLOCATE( spp(ii)%vars%X(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%V(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Rgc(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Y(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%E(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%B(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%gamma(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%eta(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%mu(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Prad(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Pin(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%flag(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%AUX(spp(ii)%ppp) )

		! Initialize to zero
		spp(ii)%vars%X = 0.0_rp
		spp(ii)%vars%V = 0.0_rp
		spp(ii)%vars%Rgc = 0.0_rp
		spp(ii)%vars%Y = 0.0_rp
		spp(ii)%vars%E = 0.0_rp
		spp(ii)%vars%B = 0.0_rp
		spp(ii)%vars%gamma = 0.0_rp
		spp(ii)%vars%eta = 0.0_rp
		spp(ii)%vars%mu = 0.0_rp
		spp(ii)%vars%Prad = 0.0_rp
		spp(ii)%vars%Pin = 0.0_rp
		spp(ii)%vars%flag = 1_idef
		spp(ii)%vars%AUX = 0.0_rp
	end do

	DEALLOCATE(ppp)
	DEALLOCATE(q)
	DEALLOCATE(m)
	DEALLOCATE(Eo)
	DEALLOCATE(etao)
	DEALLOCATE(runaway)
	DEALLOCATE(Ro)
	DEALLOCATE(Zo)
	DEALLOCATE(r)
end subroutine initialize_particles


subroutine set_up_particles_ic(params,F,spp) 
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Vo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V1
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V2
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V3
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: b1, b2, b3
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Xo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, radius, angle ! temporary vars
	REAL(rp), DIMENSION(3) :: x = (/1.0_rp,0.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: y = (/0.0_rp,1.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: z = (/0.0_rp,0.0_rp,1.0_rp/)
	INTEGER :: ii,jj ! Iterator

	do ii=1,params%num_species
		ALLOCATE( Xo(3,spp(ii)%ppp) )
		ALLOCATE( Vo(spp(ii)%ppp) )
		ALLOCATE( V1(spp(ii)%ppp) )
		ALLOCATE( V2(spp(ii)%ppp) )
		ALLOCATE( V3(spp(ii)%ppp) )
		ALLOCATE( b1(3,spp(ii)%ppp) )
		ALLOCATE( b2(3,spp(ii)%ppp) )
		ALLOCATE( b3(3,spp(ii)%ppp) )

		
		ALLOCATE( theta(spp(ii)%ppp) )
		ALLOCATE( zeta(spp(ii)%ppp) )
		ALLOCATE( radius(spp(ii)%ppp) )
		ALLOCATE( angle(spp(ii)%ppp) )

		! Initial condition of uniformly distributed particles on a disk in the xz-plane
		! A unique velocity direction
		call init_random_seed()
		call RANDOM_NUMBER(theta)
		theta = 2.0_rp*C_PI*theta

		call init_random_seed()
		call RANDOM_NUMBER(zeta)
		zeta = 2.0_rp*C_PI*zeta

		! Uniform distribution on a disk at a fixed azimuthal theta		
		call init_random_seed()
		call RANDOM_NUMBER(radius)
		
		Xo(1,:) = ( spp(ii)%Ro + spp(ii)%r*sqrt(radius)*cos(theta) )*sin(zeta)
		Xo(2,:) = ( spp(ii)%Ro + spp(ii)%r*sqrt(radius)*cos(theta) )*cos(zeta)
		Xo(3,:) = spp(ii)%Zo + spp(ii)%r*sqrt(radius)*sin(theta)

!		do jj=1,spp(ii)%ppp
!			Xo(1,jj) = (spp(ii)%Ro - spp(ii)%r) + &
!					2.0_rp*spp(ii)%r*REAL(MODULO(jj,201_idef),rp)/201.0_rp
!			Xo(2,jj) = 0.0_rp
!			Xo(3,jj) = spp(ii)%r - &
!					2.0_rp*spp(ii)%r*FLOOR(REAL(jj,rp)/201.0_rp)/201.0_rp
!		end do

		spp(ii)%vars%X(1,:) = Xo(1,:)
		spp(ii)%vars%X(2,:) = Xo(2,:)
		spp(ii)%vars%X(3,:) = Xo(3,:)

		! Monoenergetic distribution
		spp(ii)%vars%gamma(:) = spp(ii)%gammao

		call init_random_seed()
		call RANDOM_NUMBER(angle)
		angle = 2.0_rp*C_PI*angle

		Vo = sqrt( 1.0_rp - 1.0_rp/(spp(ii)%vars%gamma(:)**2) )
        V1 = Vo*cos(C_PI*spp(ii)%etao/180.0_rp)
        V2 = Vo*sin(C_PI*spp(ii)%etao/180.0_rp)*cos(angle)
        V3 = Vo*sin(C_PI*spp(ii)%etao/180.0_rp)*sin(angle)

        call unitVectors(params,Xo,F,b1,b2,b3)

		do jj=1,spp(ii)%ppp
			if ( spp(ii)%vars%flag(jj) .EQ. 1_idef ) then
				spp(ii)%vars%V(1,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),x) + &
		                                V2(jj)*DOT_PRODUCT(b2(:,jj),x) + &
		                                V3(jj)*DOT_PRODUCT(b3(:,jj),x)

				spp(ii)%vars%V(2,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),y) + &
		                                V2(jj)*DOT_PRODUCT(b2(:,jj),y) + &
		                                V3(jj)*DOT_PRODUCT(b3(:,jj),y)

				spp(ii)%vars%V(3,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),z) + &
		                                V2(jj)*DOT_PRODUCT(b2(:,jj),z) + &
		                                V3(jj)*DOT_PRODUCT(b3(:,jj),z)
			end if
		end do

		DEALLOCATE(theta)
		DEALLOCATE(zeta)
		DEALLOCATE(radius)
		DEALLOCATE(angle)	
		DEALLOCATE(Xo)
		DEALLOCATE(Vo)
		DEALLOCATE(V1)
		DEALLOCATE(V2)
		DEALLOCATE(V3)
		DEALLOCATE(b1)
		DEALLOCATE(b2)
		DEALLOCATE(b3)
	end do
end subroutine set_up_particles_ic


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_communications(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

	call initialize_mpi(params)

!$OMP PARALLEL SHARED(params)
        params%num_omp_threads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

	call initialization_sanity_check(params) 
end subroutine initialize_communications


subroutine initialization_sanity_check(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	CHARACTER(MAX_STRING_LENGTH) :: env_variable
	INTEGER :: ierr
	LOGICAL :: flag = .FALSE.

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'("* * * SANITY CHECK * * *")')
	end if

	call GET_ENVIRONMENT_VARIABLE("OMP_PLACES",env_variable)
!	call GET_ENVIRONMENT_VARIABLE("GOMP_CPU_AFFINITY",env_variable)
	write(6,*) TRIM(env_variable)

!$OMP PARALLEL SHARED(params) PRIVATE(ierr, flag)
	call MPI_INITIALIZED(flag, ierr)
	write(6,'("MPI: ",I3," OMP/of: ",I3," / ",I3," Procs: ",I3," Init: ",l1)') &
	params%mpi_params%rank,OMP_GET_THREAD_NUM(),OMP_GET_NUM_THREADS(),OMP_GET_NUM_PROCS(),flag
!$OMP END PARALLEL
end subroutine initialization_sanity_check


subroutine initialize_fields(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(OUT) :: F
	TYPE(KORC_STRING) :: field
	REAL(rp) :: Bo
	REAL(rp) :: minor_radius
	REAL(rp) :: major_radius
	REAL(rp) :: qa
	REAL(rp) :: qo
    CHARACTER(MAX_STRING_LENGTH) :: electric_field_mode
	REAL(rp) :: Eo
    REAL(rp) :: pulse_maximum
    REAL(rp) :: pulse_duration

	NAMELIST /analytic_mag_field_params/ Bo,minor_radius,major_radius,&
			qa,qo,electric_field_mode,Eo,pulse_maximum,pulse_duration

	if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		! Load the parameters of the analytical magnetic field
		open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
		read(default_unit_open,nml=analytic_mag_field_params)
		close(default_unit_open)

		F%AB%Bo = Bo
		F%AB%a = minor_radius
		F%AB%Ro = major_radius
		F%Ro = major_radius
		F%AB%qa = qa
		F%AB%qo = qo
		F%AB%lambda = F%AB%a/sqrt(qa/qo - 1.0_rp)
		F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)

		F%Eo = Eo
		F%Bo = F%AB%Bo

        F%electric_field_mode = TRIM(electric_field_mode)
		F%to = pulse_maximum
		F%sig = pulse_duration
	else if (params%magnetic_field_model .EQ. 'EXTERNAL') then
		! Load the magnetic field from an external HDF5 file
        call load_dim_data_from_hdf5(params,F%dims)

       	call ALLOCATE_FIELDS_ARRAYS(F,params%poloidal_flux)

        call load_field_data_from_hdf5(params,F)

		if (.NOT. params%poloidal_flux) then
			field%str = 'B'
			call mean_F_field(F,F%Bo,field)
		end if
	else
		write(6,'("ERROR: when initializing fields!")')
		call korc_abort()
	end if
end subroutine initialize_fields


subroutine initialize_collision_params(params,cparams)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(COLLISION_PARAMS), INTENT(OUT) :: cparams
	REAL(rp) :: Te ! Background electron temperature in eV
	REAL(rp) :: ne! Background electron density in 1/m^3
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zj ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(:), ALLOCATABLE :: nz ! Impurity densities
	REAL(rp), DIMENSION(:), ALLOCATABLE :: IZj ! Ionization energy of impurity in eV

	NAMELIST /collision_parameters/ Te,ne,Zo,Zj,nz,IZj

	cparams%num_impurity_species = params%num_impurity_species

	ALLOCATE(Zj(params%num_impurity_species))
	ALLOCATE(Zo(params%num_impurity_species))
	ALLOCATE(nz(params%num_impurity_species))
	ALLOCATE(IZj(params%num_impurity_species))

	ALLOCATE(cparams%Zj(params%num_impurity_species))
	ALLOCATE(cparams%Zo(params%num_impurity_species))
	ALLOCATE(cparams%nz(params%num_impurity_species))
	ALLOCATE(cparams%neb(params%num_impurity_species))
	ALLOCATE(cparams%IZj(params%num_impurity_species))
	ALLOCATE(cparams%Ee_IZj(params%num_impurity_species))

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=collision_parameters)
	close(default_unit_open)

	cparams%Te = Te*C_E
	cparams%ne = ne
	cparams%nH = ne

	cparams%Zj = Zj
	cparams%Zo = Zo
	cparams%nz = nz
	cparams%IZj = C_E*IZj

	cparams%nef = ne + sum(cparams%Zj*cparams%nz)
	cparams%neb = (cparams%Zo-cparams%Zj)*cparams%nz

	cparams%rD = SQRT( C_E0*cparams%Te/(cparams%ne*C_E**2) )
	cparams%re = C_E**2/( 4.0_rp*C_PI*C_E0*C_ME*C_C**2 )
	cparams%Ee_IZj = C_ME*C_C**2/cparams%IZj


	DEALLOCATE(Zj)
	DEALLOCATE(Zo)
	DEALLOCATE(nz)
	DEALLOCATE(IZj)
end subroutine initialize_collision_params

end module initialize
