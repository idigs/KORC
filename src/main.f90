program main
  !! @note  Main function of KORC. @endnote
  !! The main program contains the calls to the main functions and subroutines. 
  !! Also, it contains the variables that control
  !! the behavior of the core of KORC and all other external/optional modules.

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
  !! Contains the parameters that control the core of KORC: 
  !! time steping, output list, etc.
  TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: spp
  !! Contains the initial parameters of each species, which 
  !! can be different electrons with different
  !! distribution functions.
  TYPE(FIELDS) :: F
  !! F: Contains the parameters of the analytical magnetic 
  !! and electric fields, or in the case of using 
  !! external fields it contains the data used in the interpolations. 
  !!See [[korc_fields(module)]] for details.
  TYPE(PROFILES) :: P
  !! P: Contains the parameters of the analytical plasma profiles, 
  !! or in the case of using external 
  !! fields it contains the data used in the interpolations. 
  !! See [[korc_profiles(module)]] for details.
  INTEGER(ip) :: it 
  !! Time iteration
    INTEGER 				:: mpierr

  call initialize_communications(params)
  !!<h2>Order of KORC operations</h2>
  !!
  !!<h3>Communication and Timing</h3>
  !! <h4>1\. Parallel Communications</h4>
  !!
  !! Subroutine [[initialize_communications]] in [[korc_hpc]] that 
  !! initializes MPI and OpenMP communications.
  
  call timing_KORC(params)
  !! <h4>2\. Timers</h4>
  !!
  !! Subroutine [[timing_KORC]] in [[korc_hpc]] that times the 
  !! execution of any parallel sections of KORC.
  
  ! * * * INITIALIZATION STAGE * * *!

  call initialize_HDF5()
  !!<h3>Initialization</h3>
  !!
  !! <h4>1\. HDF5</h4>
  !!
  !! Subroutine [[initialize_HDF5]] in [[korc_HDF5]] that initializes
  !! HDF5 library. 
  
  call initialize_korc_parameters(params)
  !! <h4>2\. Initialize korc parameters</h4>
  !!
  !! Subroutine [[initialize_korc_parameters]] in [[korc_initialize]] that 
  !! initializes paths and KORC parameters through [[load_korc_params]]
  !! on MPI processes.
  
  call initialize_fields(params,F)
  !! <h4>3\. Initialize fields</h4>
  !!
  !! Subroutine [[initialize_fields]] in [[korc_fields]] that initializes 
  !! parameters of the EM fields, either analytically or from an external HDF5
  !! file. Reads in &amp;analytical_fields_params and 
  !! &amp;externalPlasmaModel namelists from input file.
  

  call initialize_profiles(params,P,F)
  !! <h4>4\. Initialize Profiles</h4>
  !! 
  !! Subroutine [[initialize_profiles]] in [[korc_profiles]] that initializes 
  !! parameters of the plasma profiles, either analytically or from an
  !! external HDF5
  !! file. Reads in &amp;plasmaProfiles namelist from input file.
  !! Only initialized if collisions (params%collisions==T) are present.


  
  call initialize_particles(params,F,P,spp) ! Initialize particles
  !! <h4>5\. Initialize Particle Velocity Phase Space</h4>
  !! 
  !! Subroutine [[initialize_particles]] in [[korc_initialize]] that 
  !! initializes particle parameters from &amplasma_species namelist, 
  !! allocates arrays for individual particles, including location, velocity, 
  !! local EM fields and plasma profiles, etc., and 
  !! calls [[initial_energy_pitch_dist]] to assign particles' energy and pitch
  !! angle according to the chosen distribution.

!  write(6,'("init eta: ",E17.10)') spp(1)%vars%eta
  



  
  call initialize_synthetic_camera(params,F)
  !! <h4>7\. Initialize Synthetic Cameras</h4>

  call initialize_binning_diagnostic(params)
  !! <h4>8\. Initialize Binning Diagnostic</h4>

  call compute_charcs_plasma_params(params,spp,F)
  !! <h4>9\. Compute Characteristic Plasma Parameters</h4>
  !!
  !! Subroutine [[compute_charcs_plasma_params]] in [[korc_units]] calculates
  !! the characteristic plasma parameters params%cpp that are used for normalizations.
  !! Also finds the maximum non-relativistic and relativistic cyclotron frequencies
  !! to be used for setting the timstep for the time-evolution algorithms.

  call initialize_collision_params(params)
  !! <h4>6\. Initialize Collision Parameters</h4>
  !!
  !! Subroutine [[initialize_collision_params]] in [[korc_collisions]] that
  !! initializes collision parameters for the SS (single-species) and MS
  !! (multiple-species) data types, reading in namefiles from the KORC input file.
  !! MS reads in namelist &CollisionParamsMultipleSpecies while SS reads in
  !! namelist &CollisionParamsSingleSpecies.
  
  call define_time_step(params)
  !! <h4>10\. Define Time Step</h4>
  !!
  !! Subroutine [[define_time_step]] in [[korc_initialize]] either loads
  !! time-stepping parameters for a restart, or defines new parameters based
  !! on a maximum timestep
  !! set by the inverse of the relativistic cyclotron frequency.
  
  call initialize_particle_pusher(params)
  !! <h4>11\. Initialize Particle Pusher</h4>    

  if (params%SC_E) then
     call define_SC_time_step(params,F)
  end if
     
  call normalize_variables(params,spp,F,P)
  !! <h4>12\. Normalize Variables</h4>
  !!
  !! Subroutine [[normalize_variables]] in [[korc_units]] normalizes 
  !! variables consistent with characteristic plasma parameters 
  !! calculated in [[compute_charcs_plasma_params]].


  
  call normalize_collisions_params(params)
  !! <h4>13\. Normalize Collision Parameters </h4>
  !!
  !! Subroutine [[normalize_collisions_params]] in [[korc_collisions]] that
  !! normalizes collision parameters for the SS (single-species) and MS
  !! (multiple-species) data types.

  
  call define_collisions_time_step(params)
  !! <h4>14\. Define Collision Time Step</h4>
  !!
  !! Subroutine [[define_collisions_time_step]] in [[korc_collisions]] that
  !! sets subcycling iteration number for collisions based off of the collision
  !! frequency model used.

  ! *** *** *** *** *** ***   *** *** *** *** *** *** ***
  ! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
  ! *** *** *** *** *** ***   *** *** *** *** *** *** ***


  
  call initialize_fields_interpolant(params,F)
  !! <h4>15\. Initialize Fields Interpolant</h4>
  !!
  !! Subroutine [[initialize_fields_interpolant]] in [[korc_interp]] calls
  !! EZspline
  !! subroutines EZspline_init for memory allocation and boundary condition
  !! setup
  !! and EZspline_setup to compute the necessary cubic coefficients needed
  !! for subsequent
  !! field interpolations. The magnetic field can be defined in terms of an
  !! axisymmetric
  !! scalar flux function, axisymmetric field, or 3D field, while the
  !! electric field
  !! can be defined as an axisymmetric or 3D field.
  
  call initialize_profiles_interpolant(params,P)
  !! <h4>16\. Initialize Profiles Interpolant</h4>
  !!
  !! Subroutine [[initialize_profiles_interpolant]] in [[korc_interp]]
  !! calls EZspline
  !! subroutines EZlinear_init for axisymmetric (flux-surface quantities) or
  !! EZspline_init for 3D profiles for memory allocation and boundary
  !! condition setup
  !! and EZspline_setup to compute the necessary cubic coefficients needed
  !! for subsequent
  !! field interpolations. 
  !! Only initialized if collisions (params%collisions==T) are present for
  !! ne, Te, Zeff
  
  if (params%mpi_params%rank .EQ. 0) then
     write(6,'("* * * * INITIALIZING INITIAL CONDITIONS * * * *")')
  end if
  call set_up_particles_ic(params,F,spp,P)
  
  if (params%mpi_params%rank .EQ. 0) then
     write(6,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
  end if
  
!  write(6,'("post ic eta: ",E17.10)') spp(1)%vars%eta
  
  !! <h4>17\. Set Particle Initial Conditions</h4>  
  !!
  !! Subroutine [[set_up_particles_ic]] in [[korc_initialize]] calls
  !! subroutines to prescribe initial conditions or load them 
  !! from file for a restart. Initial spatial values are prescribed with 
  !! [[intitial_spatial_distribution]] in [[korc_spatial_distribution]] and 
  !! initial velocity values are prescribed with [[initial_gyro_distribution]]
  !! in [[korc_velocity_distribution]].

!  if (minval(spp(1)%vars%Y(:,1)).lt.1._rp/params%cpp%length) stop 'error with init'
  
  ! * * * INITIALIZATION STAGE * * *

  ! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !


  
  call save_simulation_parameters(params,spp,F,P)

  call save_collision_params(params)
  !! <h4>18\. Save Simulation and Collision Parameters</h4>  
  !!
  !! Subroutines [[save_simulation_parameters]] in [[korc_HDF5]] and
  !! [[save_collision_params]] in [[korc_collisions]] call
  !! subroutines to save simulation and collision parameters.

!  write(6,'("GC init eta: ",E17.10)') spp(1)%vars%eta

  if (.NOT.(params%restart.OR.params%proceed.or.params%reinit)) then
     if (params%orbit_model(1:2).eq.'FO') then

        call FO_init(params,F,spp,.true.,.false.)

     else if (params%orbit_model(1:2).eq.'GC') then

        call GC_init(params,F,spp)

     end if

     if (params%SC_E) then
     
        call init_SC_E1D(params,F,spp(1))
     
     end if

  else

     call get_fields(params,spp(1)%vars,F)

     if (params%SC_E) then

        if (params%reinit) then
           call reinit_SC_E1D(params,F)
        end if
        
     
     end if
     
  end if


  
  if (.NOT.params%restart) then
     
     call save_simulation_outputs(params,spp,F) ! Save initial condition

     call synthetic_camera(params,spp) 

     call binning_diagnostic(params,spp) 

  end if
  

  
  ! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !

!  write(6,'("pre ppusher loop eta: ",E17.10)') spp(1)%vars%eta
  
  call timing_KORC(params)


  if (params%orbit_model(1:2).eq.'FO'.and.params%field_eval.eq.'eqn') then

     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOeqn_top(params,F,P,spp)

        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call synthetic_camera(params,spp) ! Synthetic camera
        call binning_diagnostic(params,spp) ! Binning diagnostic
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'FO'.and.params%field_eval.eq.'interp') then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOinterp_top(params,F,P,spp)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call synthetic_camera(params,spp) ! Synthetic camera
        call binning_diagnostic(params,spp) ! Binning diagnostic
        call save_restart_variables(params,spp,F)
     end do
  end if

!  if (params%orbit_model.eq.'FO'.and.params%field_eval.eq.'interp') then
!     do it=params%ito,params%t_steps
!        params%time = params%init_time+REAL(it,rp)*params%dt
!        params%it = it
!        if ( modulo(it,params%output_cadence) .EQ. 0_ip ) then

!           call advance_FOinterp_vars(params,spp,params%dt, &
!                .TRUE.,.FALSE.)             

!           call save_simulation_outputs(params,spp)
!           call synthetic_camera(params,spp) ! Synthetic camera
!           call binning_diagnostic(params,spp) ! Binning diagnostic
!           call save_restart_variables(params,spp)
!        else
           
!           call advance_FOinterp_vars(params,spp,params%dt, &
!                .FALSE.,.FALSE.)

!        end if
!     end do
!  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'eqn') then
     do it=params%ito,params%t_steps,params%t_skip*params%t_it_SC
        call adv_GCeqn_top(params,F,P,spp)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip*params%t_it_SC,rp)*params%dt        
        params%it = it-1_ip+params%t_skip*params%t_it_SC

        call save_simulation_outputs(params,spp,F)
        call synthetic_camera(params,spp) ! Synthetic camera
        call binning_diagnostic(params,spp) ! Binning diagnostic
        call save_restart_variables(params,spp,F)
        
     end do
  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields) then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_psi_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call synthetic_camera(params,spp) ! Synthetic camera
        call binning_diagnostic(params,spp) ! Binning diagnostic
        call save_restart_variables(params,spp,F)
     end do
  end if
  
  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       .not.(F%axisymmetric_fields)) then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_B_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call synthetic_camera(params,spp) ! Synthetic camera
        call binning_diagnostic(params,spp) ! Binning diagnostic
        call save_restart_variables(params,spp,F)
     end do
  end if
  
  call timing_KORC(params)

  ! * * * FINALIZING SIMULATION * * * 
  call finalize_HDF5()

  
  call finalize_interpolants(params)

  
  ! DEALLOCATION OF VARIABLES
  call deallocate_variables(params,F,spp)

  
  call deallocate_collisions_params(params)

  
  call finalize_communications(params)
  ! * * * FINALIZING SIMULATION * * *

  if (params%mpi_params%rank .EQ. 0) then
     write(6,'("KORC ran successfully!")')
  end if
  
end program main
