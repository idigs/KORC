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

	use korc_synthetic_camera

	implicit none

	TYPE(KORC_PARAMS) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: spp
	TYPE(FIELDS) :: EB
	INTEGER(ip) :: it ! Iterator(s)
	REAL(rp) :: t1, t2 ! variables for timing the simulation

	call initialize_communications(params)

	! * * * INITIALIZATION STAGE * * *
	call initialize_HDF5()

	call initialize_korc_parameters(params) ! Initialize korc parameters

	call initialize_fields(params,EB)

	call initialize_particles(params,EB,spp) ! Initialize particles

	call initialize_collision_params(params)

	call initialize_synthetic_camera(params) ! Synthetic camera

	call compute_charcs_plasma_params(params,spp,EB)

	call define_time_step(params)

	call define_collisions_time_step(params)

	call initialize_particle_pusher(params)

	call normalize_variables(params,spp,EB)

	call normalize_collisions_params(params)

	! *** *** *** *** *** ***   *** *** *** *** *** *** ***
	! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
	! *** *** *** *** *** ***   *** *** *** *** *** *** ***

	call initialize_interpolant(params,EB)

	call set_up_particles_ic(params,EB,spp)
	! * * * INITIALIZATION STAGE * * *

	call save_simulation_parameters(params,spp,EB)

	call save_collision_params(params)

	call advance_particles_velocity(params,EB,spp,0.0_rp,.TRUE.)

	! Save initial condition
	call save_simulation_outputs(params,spp,EB)

	t1 = MPI_WTIME()

	! Initial half-time particle push
	call advance_particles_position(params,EB,spp,0.5_rp*params%dt)

	do it=1,params%t_steps
        params%time = REAL(it,rp)*params%dt
		params%it = it

		if ( modulo(it,params%output_cadence) .EQ. 0_ip ) then
            call advance_particles_velocity(params,EB,spp,params%dt,.TRUE.)
		    call advance_particles_position(params,EB,spp,params%dt)

			write(6,'("Saving snapshot: ",I15)') it/params%output_cadence
			call save_simulation_outputs(params,spp,EB)

			call synthetic_camera(params,spp)
        else
            call advance_particles_velocity(params,EB,spp,params%dt,.FALSE.)
		    call advance_particles_position(params,EB,spp,params%dt)
        end if
	end do
	
	t2 = MPI_WTIME()

	call timing_KORC(params,t1,t2)

	! * * * FINALIZING SIMULATION * * * 
	call finalize_HDF5()

	call finalize_interpolant(params)

	! DEALLOCATION OF VARIABLES
	call deallocate_variables(params,EB,spp)

	call deallocate_collisions_params(params)

	call finalize_communications(params)
	! * * * FINALIZING SIMULATION * * * 
end program main
