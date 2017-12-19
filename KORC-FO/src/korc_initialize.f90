module korc_initialize

    use korc_types
    use korc_constants
    use korc_hpc
    use korc_HDF5
    use korc_fields
    use korc_rnd_numbers
	use hammersley_generator
	use korc_spatial_distribution
	use korc_velocity_distribution

	use korc_avalanche ! external module
	use korc_experimental_pdf ! external module
	use korc_energy_pdfs ! external module


    implicit none
	

	PRIVATE :: set_paths,&
				load_korc_params
	PUBLIC :: initialize_korc_parameters,&
				initialize_particles,&
				initialize_fields,&
				define_time_step

    contains

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! ** SUBROUTINES FOR INITIALIZING KORC PARAMETERS ** !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine set_paths(params)
	INTEGER :: argn
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

	argn = command_argument_count()
	call get_command_argument(1,params%path_to_inputs)
	call get_command_argument(2,params%path_to_outputs)

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * PATHS * * * * *")')
		write(6,'("The input file is:",A70)') TRIM(params%path_to_inputs)
		write(6,'("The output folder is:",A70)') TRIM(params%path_to_outputs)
		write(6,'("* * * * * * * * * * * * *",/)')
	end if
end subroutine set_paths


subroutine load_korc_params(params)
	TYPE (KORC_PARAMS), INTENT(INOUT) :: params
	LOGICAL :: restart
	REAL(rp) :: simulation_time
	REAL(rp) :: snapshot_frequency
	REAL(rp) :: dt
	REAL(rp) :: minimum_particle_energy
	LOGICAL :: radiation
	LOGICAL :: collisions
	CHARACTER(MAX_STRING_LENGTH) :: collisions_model
	CHARACTER(MAX_STRING_LENGTH) :: plasma_model
	LOGICAL :: poloidal_flux
	LOGICAL :: axisymmetric
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename
	CHARACTER(MAX_STRING_LENGTH) :: outputs_list
	INTEGER :: num_species
	INTEGER :: num_impurity_species
	INTEGER :: imax,imin,ii,jj,num_outputs
	INTEGER, DIMENSION(2) :: indices

	NAMELIST /input_parameters/ restart,plasma_model,poloidal_flux,magnetic_field_filename,simulation_time,axisymmetric,&
			snapshot_frequency,dt,num_species,radiation,collisions,collisions_model,outputs_list,minimum_particle_energy
	
	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=input_parameters)
	close(default_unit_open)

 	params%restart = restart
	params%simulation_time = simulation_time
	params%snapshot_frequency = snapshot_frequency

	params%dt = dt

	if (params%restart) then
		call get_last_iteration(params)
	else
		params%ito = 1_ip
	end if

	params%num_species = num_species
	params%plasma_model = TRIM(plasma_model)
	params%poloidal_flux = poloidal_flux
	params%axisymmetric = axisymmetric
	params%magnetic_field_filename = TRIM(magnetic_field_filename)
	params%minimum_particle_energy = minimum_particle_energy*C_E
	params%radiation = radiation
	params%collisions = collisions
	params%collisions_model = TRIM(collisions_model)

	! Loading list of output parameters
	imin = SCAN(outputs_list,'{')
	imax = SCAN(outputs_list,'}')

	ii = 1_idef
	jj = 1_idef
	num_outputs = 1_idef
	do while (ii.NE.0)
		ii = SCAN(outputs_list(jj:),",")
		if (ii.NE.0) then
			jj = jj + ii
			num_outputs = num_outputs + 1_idef
		end if
	end do

	ALLOCATE(params%outputs_list(num_outputs))

	if (num_outputs.GT.1_idef) then
		indices = 0_idef
		indices(2) = SCAN(outputs_list,",")
		params%outputs_list(1) = TRIM(outputs_list(imin+1_idef:indices(2)-1_idef))
		indices(1) = indices(1) + indices(2) + 1_idef
		do ii=2_idef,num_outputs
			indices(2) = SCAN(outputs_list(indices(1):),",")
			if (indices(2).EQ.0_idef) then
				params%outputs_list(ii) = TRIM(outputs_list(indices(1):imax-1_idef))
			else
				params%outputs_list(ii) = TRIM(outputs_list(indices(1):indices(1)+indices(2)-2_idef))
				indices(1) = indices(1) + indices(2)
			end if
		end do
	else
		params%outputs_list(1) = TRIM(outputs_list(imin+1_idef:imax-1_idef))
	end if

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * SIMULATION PARAMETERS * * * * *")')
		write(6,'("Simulation time: ",E17.10," s")') params%simulation_time
		write(6,'("Output frequency: ",E17.10," s")') params%snapshot_frequency
		write(6,'("Time step in fraction of relativistic gyro-period: ",F15.10)') params%dt
		write(6,'("Number of electron populations: ",I16)') params%num_species
		write(6,'("Magnetic field model: ",A50)') TRIM(params%plasma_model)
		if (TRIM(params%plasma_model).EQ.'EXTERNAL') then
			write(6,'("USINg (JFIT) poloidal flux: ", L1)') params%poloidal_flux
			write(6,'("Axisymmetric external field: ", L1)') params%axisymmetric
			write(6,'("Magnetic field file: ",A100)') TRIM(params%magnetic_field_filename)
		end if

		write(6,'("Radiation losses included: ",L1)') params%radiation
		write(6,'("collisions losses included: ",L1)') params%collisions
		if (params%collisions) then
			write(6,'("collisions model: ",A50)') TRIM(params%collisions_model)
		end if
		write(6,'("* * * * * * * * * * * * * * * * * * * * *",/)')
	end if	
end subroutine load_korc_params


subroutine initialize_korc_parameters(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	INTEGER :: mpierr

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call set_paths(params)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call load_korc_params(params)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
end subroutine initialize_korc_parameters


subroutine define_time_step(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

! 	This definition will be changed as more species and electromagnetic fields
!	are included.

	params%dt = params%dt*(2.0_rp*C_PI*params%cpp%time_r)

	params%t_steps = CEILING(params%simulation_time/params%dt,ip)
	params%output_cadence = FLOOR(params%snapshot_frequency/params%dt,ip)
	if (params%output_cadence.EQ.0_ip) params%output_cadence = 1_ip
	params%num_snapshots = params%t_steps/params%output_cadence

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * TIME STEPPING PARAMETERS * * * * *")')
		write(6,'("Number of time steps: ",I16)') params%t_steps
		write(6,'("Output cadence: ",I16)') params%output_cadence
		write(6,'("Number of outputs: ",I16)') params%num_snapshots
		write(6,'("* * * * * * * * * ** * * * * * * * * * * * *",/)')
	end if
end subroutine define_time_step

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_particles(params,F,spp) 
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: ppp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: q
	REAL(rp), DIMENSION(:), ALLOCATABLE :: m
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Eo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: etao
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Eo_lims
	REAL(rp), DIMENSION(:), ALLOCATABLE :: etao_lims
	LOGICAL, DIMENSION(:), ALLOCATABLE :: runaway
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: spatial_distribution	
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: energy_distribution
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: pitch_distribution
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Ro
	REAL(rp), DIMENSION(:), ALLOCATABLE :: PHIo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: r_inner
	REAL(rp), DIMENSION(:), ALLOCATABLE :: r_outter
	REAL(rp), DIMENSION(:), ALLOCATABLE :: falloff_rate
	INTEGER :: ii,jj, mpierr ! Iterator

	NAMELIST /plasma_species/ ppp,q,m,Eo,etao,Eo_lims,etao_lims,runaway,spatial_distribution,&
								energy_distribution,pitch_distribution,Ro,PHIo,Zo,r_inner,r_outter,falloff_rate

	! Allocate array containing variables of particles for each species
	ALLOCATE(spp(params%num_species))

	ALLOCATE(ppp(params%num_species))
	ALLOCATE(q(params%num_species))
	ALLOCATE(m(params%num_species))
	ALLOCATE(Eo(params%num_species))
	ALLOCATE(etao(params%num_species))
	ALLOCATE(Eo_lims(2_idef*params%num_species))
	ALLOCATE(etao_lims(2_idef*params%num_species))
	ALLOCATE(runaway(params%num_species))
	ALLOCATE(spatial_distribution(params%num_species))
	ALLOCATE(energy_distribution(params%num_species))
	ALLOCATE(pitch_distribution(params%num_species))
	ALLOCATE(Ro(params%num_species))
	ALLOCATE(PHIo(params%num_species))
	ALLOCATE(Zo(params%num_species))
	ALLOCATE(r_inner(params%num_species))
	ALLOCATE(r_outter(params%num_species))
	ALLOCATE(falloff_rate(params%num_species))

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=plasma_species)
	close(default_unit_open)

	do ii=1_idef,params%num_species
		spp(ii)%runaway = runaway(ii)
		spp(ii)%spatial_distribution = TRIM(spatial_distribution(ii))		
		spp(ii)%energy_distribution = TRIM(energy_distribution(ii))
		spp(ii)%pitch_distribution = TRIM(pitch_distribution(ii))
		spp(ii)%q = q(ii)*C_E
		spp(ii)%m = m(ii)*C_ME
		spp(ii)%ppp = ppp(ii)

		spp(ii)%Ro = Ro(ii)
		spp(ii)%PHIo = C_PI*PHIo(ii)/180.0_rp
		spp(ii)%Zo = Zo(ii)
		spp(ii)%r_inner = r_inner(ii)
		spp(ii)%r_outter = r_outter(ii)
		spp(ii)%falloff_rate = falloff_rate(ii)
		

		ALLOCATE( spp(ii)%vars%X(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%V(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Rgc(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Y(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%E(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%B(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%ne(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Te(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Zeff(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%g(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%eta(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%mu(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Prad(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Pin(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%flag(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%AUX(spp(ii)%ppp) )

		SELECT CASE (TRIM(spp(ii)%energy_distribution))
			CASE ('MONOENERGETIC')
				spp(ii)%Eo = Eo(ii)*C_E
				spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

				spp(ii)%vars%g = spp(ii)%go ! Monoenergetic
				spp(ii)%Eo_lims = (/spp(ii)%Eo , spp(ii)%Eo /)
			CASE ('THERMAL')
				call thermal_distribution(params,spp(ii))

				spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 , &
									spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
			CASE ('AVALANCHE')
				call get_avalanche_PDF_params(params,spp(ii)%vars%g,spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)

				spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
				spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 , &
									spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
			CASE ('EXPERIMENTAL')
				call get_experimental_distribution(params,spp(ii)%vars%g,spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)

				spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
				spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 , &
									spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
			CASE ('GAMMA')
				call get_gamma_distribution(params,spp(ii)%vars%g,spp(ii)%go)

				spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
				spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 , &
									spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
			CASE ('UNIFORM')
				spp(ii)%Eo_lims = Eo_lims((ii-1_idef)*2_idef + 1_idef:2_idef*ii)*C_E
				spp(ii)%Eo = spp(ii)%Eo_lims(1)
				spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

				call generate_2D_hammersley_sequence(params%mpi_params%rank,params%mpi_params%nmpi,spp(ii)%vars%g,spp(ii)%vars%eta)

				spp(ii)%vars%g = (spp(ii)%Eo_lims(2) - spp(ii)%Eo_lims(1))*spp(ii)%vars%g/(spp(ii)%m*C_C**2) + &
									(spp(ii)%Eo_lims(1) + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
			CASE DEFAULT
				! Something to be done
		END SELECT

		call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

		SELECT CASE (TRIM(spp(ii)%pitch_distribution))
			CASE ('MONOPITCH')
				spp(ii)%etao = etao(ii)

				spp(ii)%vars%eta = spp(ii)%etao ! Mono-pitch-angle
				spp(ii)%etao_lims = (/spp(ii)%etao , spp(ii)%etao/)
			CASE ('THERMAL')
				spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
			CASE ('AVALANCHE')
				spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
			CASE ('EXPERIMENTAL')
				spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
			CASE ('UNIFORM')
				spp(ii)%etao_lims = etao_lims((ii-1_idef)*2_idef + 1_idef:2_idef*ii)
				spp(ii)%etao = spp(ii)%etao_lims(1)

				spp(ii)%vars%eta = (spp(ii)%etao_lims(2) - spp(ii)%etao_lims(1))*spp(ii)%vars%eta + spp(ii)%etao_lims(1)
			CASE DEFAULT
				! Something to be done
		END SELECT

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * SPECIES: ",I2," * * * * * * * * * * *")') ii
			write(6,'("Energy distribution is: ",A20)') TRIM(spp(ii)%energy_distribution)
			write(6,'("Pitch-angle distribution is: ",A20)') TRIM(spp(ii)%pitch_distribution)
			write(6,'("* * * * * * * * * * * * * * * * * * * * * *",/)')
		end if

		call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	
		! Initialize to zero
		spp(ii)%vars%X = 0.0_rp
		spp(ii)%vars%V = 0.0_rp
		spp(ii)%vars%Rgc = 0.0_rp
		spp(ii)%vars%Y = 0.0_rp
		spp(ii)%vars%E = 0.0_rp
		spp(ii)%vars%B = 0.0_rp
		spp(ii)%vars%mu = 0.0_rp
		spp(ii)%vars%Prad = 0.0_rp
		spp(ii)%vars%Pin = 0.0_rp
		spp(ii)%vars%flag = 1_is
		spp(ii)%vars%AUX = 0.0_rp
	end do

	DEALLOCATE(ppp)
	DEALLOCATE(q)
	DEALLOCATE(m)
	DEALLOCATE(Eo)
	DEALLOCATE(etao)
	DEALLOCATE(Eo_lims)
	DEALLOCATE(etao_lims)
	DEALLOCATE(runaway)
	DEALLOCATE(spatial_distribution)
	DEALLOCATE(energy_distribution)
	DEALLOCATE(pitch_distribution)
	DEALLOCATE(Ro)
	DEALLOCATE(PHIo)
	DEALLOCATE(Zo)
	DEALLOCATE(r_inner)
	DEALLOCATE(r_outter)
	DEALLOCATE(falloff_rate)
end subroutine initialize_particles


subroutine set_up_particles_ic(params,F,spp) 
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Vo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V1
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V2
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V3
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: b1, b2, b3
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta ! temporary vars
	REAL(rp), DIMENSION(3) :: x = (/1.0_rp,0.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: y = (/0.0_rp,1.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: z = (/0.0_rp,0.0_rp,1.0_rp/)
	INTEGER :: ii,jj ! Iterator

	if (params%restart) then
		call load_particles_ic(params,spp)
	else
		call intitial_spatial_distribution(params,spp)

		call initial_velocity_distribution(params,F,spp)
	end if
end subroutine set_up_particles_ic


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * *  SUBROUTINES FOR INITIALIZING FIELDS   * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_fields(params,F)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(OUT) :: F
	TYPE(KORC_STRING) :: field
	REAL(rp) :: Bo
	REAL(rp) :: minor_radius
	REAL(rp) :: major_radius
	REAL(rp) :: qa
	REAL(rp) :: qo
    CHARACTER(MAX_STRING_LENGTH) :: current_direction
    CHARACTER(MAX_STRING_LENGTH) :: electric_field_mode
	REAL(rp) :: Eo
    REAL(rp) :: pulse_maximum
    REAL(rp) :: pulse_duration

	NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
			qa,qo,electric_field_mode,Eo,pulse_maximum,pulse_duration,current_direction

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			! Load the parameters of the analytical magnetic field
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=analytical_fields_params)
			close(default_unit_open)

			F%AB%Bo = Bo
			F%AB%a = minor_radius
			F%AB%Ro = major_radius
			F%Ro = major_radius
			F%AB%qa = qa
			F%AB%qo = qo
			F%AB%lambda = F%AB%a/SQRT(qa/qo - 1.0_rp)
			F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)
			F%AB%current_direction = TRIM(current_direction)
			SELECT CASE (TRIM(F%AB%current_direction))
				CASE('PARALLEL')
					F%AB%Bp_sign = 1.0_rp
				CASE('ANTI-PARALLEL')
					F%AB%Bp_sign = -1.0_rp
				CASE DEFAULT
			END SELECT
			F%Eo = Eo
			F%Bo = F%AB%Bo

		    F%electric_field_mode = TRIM(electric_field_mode)
			F%to = pulse_maximum
			F%sig = pulse_duration
		CASE('EXTERNAL')
			! Load the magnetic field from an external HDF5 file
		    call load_dim_data_from_hdf5(params,F%dims)

			if (params%poloidal_flux) then
				call ALLOCATE_FLUX_ARRAYS(F)
			else if (params%axisymmetric) then
				call ALLOCATE_2D_FIELDS_ARRAYS(F)
			else
				call ALLOCATE_3D_FIELDS_ARRAYS(F)
			end if
		
		    call load_field_data_from_hdf5(params,F)
		CASE('UNIFORM')
			! Load the parameters of the analytical magnetic field
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=analytical_fields_params)
			close(default_unit_open)

			F%Eo = Eo
			F%Bo = Bo
		CASE DEFAULT
	END SELECT
end subroutine initialize_fields

end module korc_initialize
