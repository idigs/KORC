!> @brief Main function of KORC.
!! @author Leopoldo Carbajal
!! @details The main function contains the calls to the main functions and subroutines. Also, it contains the variables that control
!! the behavior of the core of KORC and all other external/optional modules.
!! 
!! @param[in,out] params Contains the parameters that control the core of KORC: time steping, output list
!! @param[in,out] spp Contains the initial parameters of each species, which can be different electrons with different
!! distribution functions.
!! @param[in] F Contains the parameters of the analytical magnetic and electric fields, or in the case of using external fields it
!! contains the data used in the interpolations. See korc_fields.90 for details.
!! @param[in] P Contains the parameters of the analytical plasma profiles, or in the case of using external fields it contains the 
!! data used in the interpolations. See korc_profiles.90 for details.
!! @warning This documentation is a living document. If you find an issue please report it immediatly through the "issues" section
!! of the github repository.
program main
	use korc_types
	use korc_units
	use korc_hpc
	use korc_HDF5
	use korc_fields
	use korc_ppusher
	use korc_interp
	use korc_collisions
	use korc_initialize
	use korc_finalize
	use korc_profiles

	use korc_synthetic_camera
	use korc_binning_diagnostic

	implicit none

	TYPE(KORC_PARAMS) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: spp
	TYPE(FIELDS) :: F
	TYPE(PROFILES) :: P
	INTEGER(ip) :: it ! Iterator(s)
	REAL(rp) :: t1, t2 ! variables for timing the simulation

	call initialize_communications(params)

	t1 = MPI_WTIME()

	! * * * INITIALIZATION STAGE * * *
	call initialize_HDF5()

	call initialize_korc_parameters(params) ! Initialize korc parameters

	call initialize_fields(params,F)

	call initialize_profiles(params,P)
	
	call initialize_particles(params,F,spp) ! Initialize particles

	call initialize_collision_params(params)

	call initialize_synthetic_camera(params,F) ! Synthetic camera

	call initialize_binning_diagnostic(params) ! Binning diagnostic

	call compute_charcs_plasma_params(params,spp,F)

	call define_time_step(params)

	call initialize_particle_pusher(params)

	call normalize_variables(params,spp,F,P)

	call normalize_collisions_params(params)

	call define_collisions_time_step(params)

	! *** *** *** *** *** ***   *** *** *** *** *** *** ***
	! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
	! *** *** *** *** *** ***   *** *** *** *** *** *** ***

	call initialize_fields_interpolant(params,F)

	call initialize_profiles_interpolant(params,P)

	call set_up_particles_ic(params,F,spp)
	! * * * INITIALIZATION STAGE * * *

	! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !

	call save_simulation_parameters(params,spp,F,P)

	call save_collision_params(params)
		
	if (.NOT.params%restart) then
		call advance_particles_velocity(params,F,P,spp,0.0_rp,.TRUE.)

		! Save initial condition
		call save_simulation_outputs(params,spp,F)

		call synthetic_camera(params,spp) ! Synthetic camera!

		call binning_diagnostic(params,spp) ! Binning diagnostic
	end if
	! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !

	t2 = MPI_WTIME()

	call timing_KORC(params,t1,t2)

	t1 = MPI_WTIME()

	! Initial half-time particle push
	call advance_particles_position(params,F,spp,0.5_rp*params%dt)

	do it=params%ito,params%t_steps
        params%time = REAL(it,rp)*params%dt
		params%it = it
		if ( modulo(it,params%output_cadence) .EQ. 0_ip ) then
            call advance_particles_velocity(params,F,P,spp,params%dt,.TRUE.)
		    call advance_particles_position(params,F,spp,params%dt)

			call save_simulation_outputs(params,spp,F)
			call synthetic_camera(params,spp) ! Synthetic camera
			call binning_diagnostic(params,spp) ! Binning diagnostic
        else
            call advance_particles_velocity(params,F,P,spp,params%dt,.FALSE.)
		    call advance_particles_position(params,F,spp,params%dt)
        end if
		call save_restart_variables(params,spp,F)
	end do
	
	t2 = MPI_WTIME()

	call timing_KORC(params,t1,t2)

	! * * * FINALIZING SIMULATION * * * 
	call finalize_HDF5()

	call finalize_interpolants(params)

	! DEALLOCATION OF VARIABLES
	call deallocate_variables(params,F,spp)

	call deallocate_collisions_params(params)

	call finalize_communications(params)
	! * * * FINALIZING SIMULATION * * * 
	!...
end program main

