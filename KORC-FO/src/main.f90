program main

	use korc_types
	use korc_units
	use korc_hpc
	use korc_HDF5
	use korc_fields
	use korc_ppusher
	use korc_interp
	use initialize
	use finalize

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

	call compute_charcs_plasma_params(params,spp,EB)

	call define_time_step(params)

	call normalize_variables(params,spp,EB)

	call initialize_interpolant(params,EB)

	call set_up_particles_ic(params,EB,spp)
	! * * * INITIALIZATION STAGE * * *

	call save_simulation_parameters(params,spp,EB)

	! *** *** *** *** *** ***   *** *** *** *** *** *** ***
	! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
	! *** *** *** *** *** ***   *** *** *** *** *** *** ***

	call advance_particles_velocity(params,EB,spp,0.0_rp,.FALSE.)

	! Save initial condition
	call save_simulation_outputs(params,spp,EB,0_ip)

	t1 = MPI_WTIME()

	! Initial half-time particle push
	call advance_particles_position(params,EB,spp,0.5_rp*params%dt)

	do it=1,params%t_steps
!		call advance_particles_velocity(params,EB,spp,params%dt,.FALSE.)
!		call advance_particles_position(params,EB,spp,params%dt)
		if ( modulo(it,params%output_cadence) .EQ. 0_ip ) then
            call advance_particles_velocity(params,EB,spp,params%dt,.TRUE.)
		    call advance_particles_position(params,EB,spp,params%dt)
			write(6,'("Saving snapshot: ",I15)') it/params%output_cadence
			call save_simulation_outputs(params,spp,EB,it)
        else
            call advance_particles_velocity(params,EB,spp,params%dt,.FALSE.)
		    call advance_particles_position(params,EB,spp,params%dt)
        end if
	end do
	
	t2 = MPI_WTIME()
	write(6,'("MPI: ",I2," Total time: ",F15.10)') params%mpi_params%rank, t2 - t1

	! * * * FINALIZING SIMULATION * * * 
	call finalize_HDF5()

	call finalize_interpolant(params)

	! DEALLOCATION OF VARIABLES
	call deallocate_variables(params,EB,spp)

	call finalize_communications(params)
	! * * * FINALIZING SIMULATION * * * 
end program main
