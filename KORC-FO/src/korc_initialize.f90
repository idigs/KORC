!> @brief Module with subroutines to load simulation parameters and to define the time step in the simulation.
module korc_initialize
    use korc_types
    use korc_constants
    use korc_hpc
    use korc_HDF5
    use korc_fields
    use korc_rnd_numbers
	use korc_spatial_distribution
	use korc_velocity_distribution

    IMPLICIT NONE


	PRIVATE :: set_paths,&
				load_korc_params
	PUBLIC :: initialize_korc_parameters,&
				initialize_particles,&
				define_time_step

    CONTAINS

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! ** SUBROUTINES FOR INITIALIZING KORC PARAMETERS ** !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

!> @brief Subroutine that sets the input/output paths.
!!
!! @param[in,out] params Core KORC simulation parameters.
!! @param argn Number of command line inputs. The default value is two: the input files path and the outputs path.
subroutine set_paths(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
	INTEGER 							:: argn

	argn = command_argument_count()

	if (argn .EQ. 2_idef) then
		call get_command_argument(1,params%path_to_inputs)
		call get_command_argument(2,params%path_to_outputs)
	else
		call korc_abort()
	end if

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * PATHS * * * * *")')
		write(6,'("The input file is:",A70)') TRIM(params%path_to_inputs)
		write(6,'("The output folder is:",A70)') TRIM(params%path_to_outputs)
		write(6,'("* * * * * * * * * * * * *",/)')
	end if
end subroutine set_paths

!> @brief Subroutine that loads the simulation parameters from the file specified in params\%path_to_inputs
!!
!! @param[in,out] params Core KORC simulation parameters.
!! @param restart Flag to indicate if the simulations restarts (restart=T) or not (restart=F).
!! @param simulation_time Total simulation time in seconds.
!! @param snapshot_frequency Time between snapshots in time of the simulation.
!! @param dt Time step in the simulation as a fraction of the relativistic electron gyro-period @f$\tau_e = 2\pi\gamma m_e/eB_0@f$
!! @param minimum_particle_energy Minimum allowed relativistic factor @f$\gamma@f$ of simulated electrons.
!! @param radiation Flag to indicate if synchrotron radiation losses are included (radiation=T) or not (radiation=F).
!! @param collisions Flag to indicate if collisionsare included (collisions=T) or not (collisions=F).
!! @param collisions_model String with the name of the collisions model to be used in the simulation.
!! @param plasma_model String with the name of the model for the fields and plasma profiles.
!! @param magnetic_field_filename String with the name of the model for the fields and plasma profiles.
!! @param outputs_list List of electron variables to include in the outputs.
!! @param num_species Number of different populations of simulated relativistic electrons in KORC.
!! @param imax Auxiliary variable used to parse the output_list
!! @param imin Auxiliary variable used to parse the output_list
!! @param ii Iterator used to parse the output_list
!! @param jj Iterator used to parse the output_list
!! @param num_outputs Auxiliary variable used to parse the output_list
!! @param indices Auxiliary variable used to parse the output_list
!! @param HDF5_error_handling
subroutine load_korc_params(params)
	TYPE (KORC_PARAMS), INTENT(INOUT) :: params
	LOGICAL                           :: restart
	REAL(rp)                          :: simulation_time
	REAL(rp)                          :: snapshot_frequency
	REAL(rp)                          :: dt
	REAL(rp)                          :: minimum_particle_energy
	LOGICAL                           :: radiation
	LOGICAL                           :: collisions
	CHARACTER(MAX_STRING_LENGTH)      :: collisions_model
	CHARACTER(MAX_STRING_LENGTH)      :: plasma_model
	CHARACTER(MAX_STRING_LENGTH)      :: magnetic_field_filename
	CHARACTER(MAX_STRING_LENGTH)      :: outputs_list
	INTEGER                           :: num_species
	INTEGER                           :: imax
	INTEGER                           :: imin
	INTEGER                           :: ii
	INTEGER                           :: jj
	INTEGER                           :: num_outputs
	INTEGER, DIMENSION(2)             :: indices
	LOGICAL                           :: HDF5_error_handling
    INTEGER                           :: time_slice

	NAMELIST /input_parameters/ restart,plasma_model,magnetic_field_filename,simulation_time,&
								snapshot_frequency,dt,num_species,radiation,collisions,collisions_model,outputs_list,&
								minimum_particle_energy,HDF5_error_handling,time_slice

    time_slice = 0

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=input_parameters)
	close(default_unit_open)

 	params%restart = restart

	params%simulation_time = simulation_time
	params%snapshot_frequency = snapshot_frequency
	params%dt = dt

	params%num_species = num_species
	params%plasma_model = TRIM(plasma_model)
	params%magnetic_field_filename = TRIM(magnetic_field_filename)
    params%time_slice = time_slice
	params%minimum_particle_energy = minimum_particle_energy*C_E
	params%minimum_particle_g = 1.0_rp + params%minimum_particle_energy/(C_ME*C_C**2) ! Minimum value of relativistic gamma factor
	params%radiation = radiation
	params%collisions = collisions
	params%collisions_model = TRIM(collisions_model)
	if (HDF5_error_handling) then
		params%HDF5_error_handling = 1_idef
	else
		params%HDF5_error_handling = 0_idef
	end if

	! Loading list of output parameters (parsing)
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
		write(6,'("Restarting simulation: ",L1)') params%restart
		write(6,'("Number of electron populations: ",I16)') params%num_species
		write(6,'("Magnetic field model: ",A50)') TRIM(params%plasma_model)
		if (TRIM(params%plasma_model).EQ.'EXTERNAL') then
			write(6,'("Magnetic field file: ",A100)') TRIM(params%magnetic_field_filename)
		end if

		write(6,'("Radiation losses included: ",L1)') params%radiation
		write(6,'("Collisions losses included: ",L1)') params%collisions
		if (params%collisions) then
			write(6,'("Collisions model: ",A50)') TRIM(params%collisions_model)
		end if
		write(6,'("* * * * * * * * * * * * * * * * * * * * *",/)')
	end if
end subroutine load_korc_params


!> @brief Interface for calling initialization subroutines
!!
!! @param[in,out] params Core KORC simulation parameters.
!! @param mpierr MPI error status.
subroutine initialize_korc_parameters(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
	INTEGER 							:: mpierr

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call set_paths(params)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call load_korc_params(params)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
end subroutine initialize_korc_parameters


!> @brief Subroutine that defines or loads from restart file the time stepping parameters.
!!
!! @param[in,out] params Core KORC simulation parameters.
subroutine define_time_step(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

	if (params%restart) then
		call load_time_stepping_params(params)
	else
		params%ito = 1_ip

		params%dt = params%dt*(2.0_rp*C_PI*params%cpp%time_r)

		params%t_steps = CEILING(params%simulation_time/params%dt,ip)

		params%output_cadence = FLOOR(params%snapshot_frequency/params%dt,ip)

		if (params%output_cadence.EQ.0_ip) params%output_cadence = 1_ip

		params%num_snapshots = params%t_steps/params%output_cadence

		if (params%output_cadence.LT.1000_ip) then
			params%restart_output_cadence = 1000_ip
		else
			params%restart_output_cadence = params%output_cadence
		end if
	end if

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * TIME STEPPING PARAMETERS * * * * *")')
		write(6,'("Simulation time: ",E17.10," s")') params%simulation_time
		write(6,'("Output frequency: ",E17.10," s")') params%snapshot_frequency
		write(6,'("Time step in fraction of relativistic gyro-period: ",E17.10)') params%dt
		write(6,'("Number of time steps: ",I16)') params%t_steps
		write(6,'("Starting simulation at time step: ",I16)') params%ito
		write(6,'("Output cadence: ",I16)') params%output_cadence
		write(6,'("Number of outputs: ",I16)') params%num_snapshots
		write(6,'("* * * * * * * * * ** * * * * * * * * * * * *",/)')
	end if
end subroutine define_time_step

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

!> @brief Subroutine that loads the information of the initial condition of the different particle species. This subroutine calls
!! the subroutine that generates the initial energy and pitch angle distribution functions.
!!
!! @param[in,out] params Core KORC simulation parameters.
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation. See korc_types.f90 and korc_fields.f90.
!! @param[out] spp An instance of KORC's derived type SPECIES containing all the information of different electron species. See korc_types.f90.
!! @param ppp Number of computational particles used to simulate each electron species.
!! @param q Charge of each species.
!! @param m Mass of each species.
!! @param Eo Initial energy of each electron species in case of using an initial mono-energetic distribution.
!! @param etao Initial pitch-angle of each electron species in case of using an initial mono-pitch-angle distribution.
!! @param Eo_lims Minimum and maximum energy limits of a given initial non-mono-energetic distribution.
!! @param etao_lims Minimum and maximum pitch-angle limits of a given initial non-mono-pitch-angle distribution.
!! @param runaway Flag to decide whether a given electron is a runaway (runaway=T) or not (runaway=F).
!! @param spatial_distribution String describing the type of initial spatial distribution for each electron species.
!! @param energy_distribution String describing the type of initial energy distribution for each electron species.
!! @param pitch_distribution String describing the type of initial pitch-angle distribution for each electron species.
!! @param Ro Radial position of the center of the electrons' initial spatial distribution.
!! @param PHIo Azimuthal position of the electrons' initial spatial distribution, in case of using a disk at a certain poloidal section.
!! @param Zo Height of the center of the electrons' initial spatial distribution.
!! @param r_inner Minimum minor radius of the electrons' initial spatial distribution.
!! @param r_outter Maximum minor radius of the electrons' initial spatial distribution.
!! @param falloff_rate Exponential falloff or standard deviation of a non-uniform radial distribution of electrons.
!! @param shear_factor Shear factor used to generate an initial spatial distribution with an elliptic poloidal cross section.
!! See <em>Carbajal and del-Castillo-Negrete, Nuclear Fusion, submitted (2018)</em>.
!! @param ii Iterator of spp structure.
!! @param mpierr MPI error status.
subroutine initialize_particles(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) 							:: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) 	:: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: ppp
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: q
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: m
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: Eo
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: etao
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: Eo_lims
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: etao_lims
	LOGICAL, DIMENSION(:), ALLOCATABLE 						:: runaway
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: spatial_distribution
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: energy_distribution
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: pitch_distribution
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: Ro
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: PHIo
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: Zo
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: r_inner
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: r_outter
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: falloff_rate
	REAL(rp), DIMENSION(:), ALLOCATABLE 					:: shear_factor
	INTEGER 												:: ii
	INTEGER 												:: mpierr

	NAMELIST /plasma_species/ ppp,q,m,Eo,etao,Eo_lims,etao_lims,runaway,spatial_distribution,&
								energy_distribution,pitch_distribution,Ro,PHIo,Zo,r_inner,r_outter,falloff_rate,shear_factor

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
	ALLOCATE(shear_factor(params%num_species))

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
		spp(ii)%shear_factor = shear_factor(ii)

		! * * These values can change in initial_energy_pitch_dist * * !
		spp(ii)%Eo = Eo(ii)*C_E
		spp(ii)%Eo_lims = Eo_lims((ii-1_idef)*2_idef + 1_idef:2_idef*ii)*C_E
		spp(ii)%etao = etao(ii)
		spp(ii)%etao_lims = etao_lims((ii-1_idef)*2_idef + 1_idef:2_idef*ii)
		! * * These values can change in initial_energy_pitch_dist * * !

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
		ALLOCATE( spp(ii)%vars%wt(spp(ii)%ppp) )

		spp(ii)%vars%X = 0.0_rp
		spp(ii)%vars%V = 0.0_rp
		spp(ii)%vars%Rgc = 0.0_rp
		spp(ii)%vars%Y = 0.0_rp
		spp(ii)%vars%E = 0.0_rp
		spp(ii)%vars%B = 0.0_rp
		spp(ii)%vars%ne = 0.0_rp
		spp(ii)%vars%Te = 0.0_rp
		spp(ii)%vars%Zeff = 0.0_rp
		spp(ii)%vars%g = 0.0_rp
		spp(ii)%vars%eta = 0.0_rp
		spp(ii)%vars%mu = 0.0_rp
		spp(ii)%vars%Prad = 0.0_rp
		spp(ii)%vars%Pin = 0.0_rp
		spp(ii)%vars%flag = 1_is
		spp(ii)%vars%AUX = 0.0_rp
		spp(ii)%vars%wt = 0.0_rp
	end do

	call initial_energy_pitch_dist(params,spp)

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

!> @brief Subroutine with calls to subroutines to load particles' information if it is a restarting simulation, or to initialize the
!! spatial and velocity distribution of each species if it is a new simulation.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation. See korc_types.f90 and korc_fields.f90.
!! @param[in,out] spp An instance of KORC's derived type SPECIES containing all the information of different electron species. See korc_types.f90.
subroutine set_up_particles_ic(params,F,spp)
	TYPE(KORC_PARAMS), INTENT(IN) 							:: params
	TYPE(FIELDS), INTENT(IN) 								:: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp

	if (params%restart) then
		call load_particles_ic(params,spp)
	else
		call intitial_spatial_distribution(params,spp)

		call initial_gyro_distribution(params,F,spp)
	end if
end subroutine set_up_particles_ic

end module korc_initialize
