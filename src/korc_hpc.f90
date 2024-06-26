module korc_hpc
  !! @note KORC module containing subroutines to initilize, control, 
  !! and to finalize parallel communications. @endnote
  use korc_types
  use omp_lib
  use mpi

  IMPLICIT NONE

  LOGICAL, PRIVATE :: timed_already = .FALSE. 
  !! Flag to determine if a first call to WMPI_TIME() was made already.
  REAL(rp), PRIVATE :: t1 
  !! Variable to be used in timing a parallel section of KORC.
  REAL(rp), PRIVATE :: t2 
  !! Variable to be used in timing a parallel section of KORC.

  PUBLIC :: korc_abort,&
       initialize_mpi,&
       finalize_mpi,&
       initialize_communications,&
       timing_KORC

CONTAINS

  subroutine korc_abort(errorcode)
    !! @note Subroutine that terminates the simulation. @endnote
    INTEGER,INTENT(IN) :: errorcode
    INTEGER :: mpierr 
    !! MPI error status

    !! 11: korc_hpc:set_paths
    !! 12: korc_experimental:get_Hollmann_distribution
    !! 13: korc_input:read_namelist
    !! 14: korc_hpc:load_particle_ic
    !! 15: korc_interp:interp_fields
    !! 16: korc_interp:interp_profiles
    !! 17: korc_fields:mean_F_field
    !! 18: korc_fields:initialize_fields
    !! 19: korc_spatial_distribution:initial_spatial_distribution
    !! 20: korc_experimental_pdf:load_data_from_hdf5
    !! 21: korc_interp:get_fio_ion_p
    !! 22: korc_fio_interface:initialize_nimrod
    !! 23: korc_interp:check_if_in_fields_domain
    !! 24: korc_collisions:large_angle_source
    !! 25: korc_ppusher:adv_GCinterp_psiwE_top
    !! 26: korc_collisions:define_collisions_time_step
    !! 27: main:initialize_fields_interpolant

    flush(output_unit_write)
    
    !call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    call MPI_ABORT(MPI_COMM_WORLD, errorcode, mpierr)
  end subroutine korc_abort

  subroutine set_paths(params)
    !! @note Subroutine that sets the input/output paths.@endnote
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    INTEGER 				:: argn,read_stat,exei
    !! Number of command line inputs. The default value is 
    !! two: the input files path and the outputs path.
    CHARACTER(MAX_STRING_LENGTH) :: ctmp

    argn = command_argument_count()

    if (argn .EQ. 2_idef) then
       call get_command_argument(1,params%path_to_inputs)
       call get_command_argument(2,params%path_to_outputs)
    else if (params%path_to_inputs.eq.'TEST') then
       call get_command_argument(1,params%path_to_outputs)
    else
       call korc_abort(11)
    end if

    !write(6,*) TRIM(params%path_to_outputs)
    !write(6,*) TRIM(params%path_to_inputs)
    
    if (params%mpi_params%rank .EQ. 0) then

       OPEN(UNIT=output_unit_write, &
            FILE=TRIM(params%path_to_outputs)//"output.korc", &
            STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')
       
       write(output_unit_write,'(/,"* * * * * PATHS * * * * *")')
       write(output_unit_write,*) 'The input file is: ',&
            TRIM(params%path_to_inputs)
       write(output_unit_write,'(/)')      
       write(output_unit_write,*) 'The output folder is: ',&
            TRIM(params%path_to_outputs)
       write(output_unit_write,'("* * * * * * * * * * * * *",/)')      
    end if
  end subroutine set_paths

  subroutine initialize_mpi(params)
    !! @note Subroutine that initialize MPI communications.@endnote
    !! Through this subroutine the default MPI communicator MPI_COMM_WORLD 
    !! is initialized. Also, a Cartesian
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    INTEGER 				:: mpierr
    !! MPI error status.
    INTEGER, PARAMETER 			:: NDIMS = 1
    !! Number of dimensions of non-standard topology. 
    !! NDIMS=1 for a 1-D MPI topology, NDIMS=2 for a 2-D MPI topology, 
    !! and NDIMS=3 for a 3-D MPI topology.
    INTEGER, DIMENSION(:), ALLOCATABLE 	:: DIMS
    !! Dimension of the non-standard MPI topology params::mpi_params::mpi_topo. 
    !! This is equal to the number of MPI processes in KORC.
    LOGICAL 				:: all_mpis_initialized = .FALSE.
    !! Flag to determine if all the MPI processes were initialized correctly.
    LOGICAL 				:: mpi_process_initialized = .FALSE.
    !! Flag to determine if a given MPI process was initialized correctly.
    LOGICAL, PARAMETER 			:: REORDER = .FALSE.
    !! Flag to determine if the new MPI topology params::mpi_params::mpi_topo 
    !! needs to be re-ordered.
    LOGICAL, DIMENSION(:), ALLOCATABLE 	:: PERIODS !< Something here
    !! Array of logicals determining what dimensions of the new MPI 
    !! topology params::mpi_params::mpi_topo are periodic (T) or not (F).
    INTEGER :: ii
    !! Variable to iterate over different MPI processes.
    LOGICAL :: mpiinit = .FALSE.

    !call MPI_INITIALIZED(mpiinit,mpierr)
    !write(6,*) 'initialized after',mpiinit
    
    
    call MPI_INIT(mpierr)

    !write(6,*) 'mpi_init error code',mpierr
    
    !call MPI_INITIALIZED(mpiinit,mpierr)
    !write(6,*) 'initialized after',mpiinit
    

    if (mpierr .NE. MPI_SUCCESS) then
       write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
       write(6,'(/," ERROR: Initializing MPI. Aborting... ")')
       write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
       call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
    end if
    
    call MPI_INITIALIZED(mpi_process_initialized,mpierr)

    !write(6,*) 'initialized after',mpi_process_initialized
    
    call MPI_REDUCE(mpi_process_initialized,all_mpis_initialized,1, &
         MPI_LOGICAL,MPI_LAND,0,MPI_COMM_WORLD,mpierr)

    !write(6,*) 'made it here 2'
    
    call MPI_BCAST(all_mpis_initialized,1, &
         MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, params%mpi_params%nmpi, mpierr)
    if (mpierr .NE. MPI_SUCCESS) then
       write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
       write(6,'(/," ERROR: Obtaining size of communicator. Aborting... ")')
       write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
       call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
    end if

    ! * * * Getting the rank of the MPI process in the WORLD COMMON communicator * * * !
    call MPI_COMM_RANK(MPI_COMM_WORLD, params%mpi_params%rank, mpierr)
    if (mpierr .NE. MPI_SUCCESS) then
       write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
       write(6,'(/," ERROR: Obtaining MPI rank. Aborting... ")')
       write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
       call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
    end if
    
    ! * * * Here a Cartesian topology for MPI is created * * * !
    ALLOCATE(DIMS(NDIMS))
    ALLOCATE(PERIODS(NDIMS))

    ! This loop isn't necessary but helps to do things more general in the future
    do ii=1_idef,NDIMS
       DIMS(ii) = params%mpi_params%nmpi
       PERIODS(ii) = .TRUE.
    end do

    ! * * * Here a periodic topology for MPI is created * * * !
    call MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, PERIODS, REORDER, &
         params%mpi_params%mpi_topo, mpierr)
    if (mpierr .NE. MPI_SUCCESS) then
       write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
       write(6,'(/," ERROR: Creating new MPI topology. Aborting... ")')
       write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
       call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
    end if

    ! * * * Getting the rank of the MPI process in the new topology * * * !
    call MPI_COMM_RANK(params%mpi_params%mpi_topo, params%mpi_params%rank_topo, mpierr)
    if (mpierr .NE. MPI_SUCCESS) then
       write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
       write(6,'(/," ERROR: Obtaining new MPI topology ranks. Aborting... ")')
       write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
       call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
    end if

    DEALLOCATE(DIMS)
    DEALLOCATE(PERIODS)


    if (all_mpis_initialized) then

       call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

       call set_paths(params)
       
       if (params%mpi_params%rank.EQ.0) then
          write(output_unit_write,'(/,"* * * * * COMMUNICATIONS  * * * * *")')
          write(output_unit_write,'(/,"  MPI communications initialized!  ")')
          write(output_unit_write,'(/,"  Number of MPI processes: ",I5)') params%mpi_params%nmpi
          write(output_unit_write,'(/,"* * * * * * * * * * * * * * * * * *")')
       end if

    else
       if (params%mpi_params%rank.EQ.0) then
          write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
          write(6,'(/," ERROR: MPI not initialized. Aborting... ")')
          write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
          call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
       end if
    end if
 
    !...
  end subroutine initialize_mpi


  subroutine timing_KORC(params)
    !! @note Subroutine for timing the execution of any parallel
    !! section of KORC. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
    REAL(rp) 				:: individual_runtime
    !! Execution time of each MPI process.
    REAL(rp), DIMENSION(:), ALLOCATABLE   :: runtime
    !! Execution time of KORC defined as the average of the 
    !! execution times of all MPI processes.
    INTEGER 				:: mpierr
    !! MPI error status.
    
    if (timed_already) then
       t2 = MPI_WTIME()

       ALLOCATE(runtime(params%mpi_params%nmpi))

       individual_runtime = t2 - t1

       call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

       call MPI_GATHER(individual_runtime,1,MPI_DOUBLE_PRECISION,runtime, &
            1,MPI_DOUBLE_PRECISION,0_idef, MPI_COMM_WORLD, mpierr)

       if (params%mpi_params%rank .EQ. 0_idef) then
          write(output_unit_write,'("Timing: ",F30.16," s")') &
               SUM(runtime)/REAL(params%mpi_params%nmpi,rp)
       end if

       call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

       DEALLOCATE(runtime)

       timed_already = .FALSE.
    end if

    t1 = MPI_WTIME()

    timed_already = .TRUE.
    
  end subroutine timing_KORC

  !> @brief Subroutine for finalizing MPI communications.
  !! @details This subroutine finalizes all the MPI communications and looks for errors durignt this procces.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param mpi_process_finalized Flag indicating whether an individual MPI process was finalized correctly.
  !! @param mpierr MPI error status.
  subroutine finalize_mpi(params)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    LOGICAL 						:: mpi_process_finalized = .FALSE.
    INTEGER 						:: mpierr

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    call MPI_FINALIZE(mpierr)

    call MPI_FINALIZED(mpi_process_finalized,mpierr)

    if (.NOT.mpi_process_finalized) then
       write(output_unit_write,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
       write(output_unit_write,'(/," ERROR: MPI not finalized well. MPI process: ",I5)') params%mpi_params%rank
       write(output_unit_write,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
    end if
  end subroutine finalize_mpi


  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! ** SUBROUTINES FOR INITIALIZING COMMUNICATIONS  ** !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  !! @note Subroutine for initializing MPI and open MP communications. @endnote
  !! This subroutine initializes MPI and open MP communications and looks for
  !! errors durignt this procces. The system environment
  !! variables, which are modified by the user at the moment of
  !! running/submitting a KORC simulation, are used to determine the
  !! open MP configuration. Some open MP parameters are displayed on the
  !! screen/output file.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param mpi_process_finalized Flag indicating whether an individual MPI process was finalized correctly.
  !! @param mpierr MPI error status.
  subroutine initialize_communications(params)
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    CHARACTER(MAX_STRING_LENGTH) 		:: string

    call initialize_mpi(params)

    if (params%mpi_params%rank.EQ.0) then
#ifdef ACC
      write(output_unit_write,'(/," Running on GPUs")')
#else
      write(output_unit_write,'(/," Running on CPUs")')
#endif
    endif

    params%num_omp_threads = get_num_threads()

    if (params%mpi_params%rank.EQ.0) then
       write(output_unit_write,'(/,"* * * * * * * OMP SET-UP * * * * * * *")')
       !$OMP PARALLEL
       !$OMP MASTER
       write(output_unit_write,'(/,"OMP threads per MPI process: ",I3)') params%num_omp_threads
       write(output_unit_write,'(/,"Cores available per MPI process: ",I3)') get_thread_number()
       !$OMP END MASTER
       !$OMP END PARALLEL
#ifdef GNU
       call GET_ENVIRONMENT_VARIABLE("OMP_PLACES",string)
       write(output_unit_write,'(/,"OMP places: ",A30)') TRIM(string)
       call GET_ENVIRONMENT_VARIABLE("GOMP_CPU_AFFINITY",string)
       write(output_unit_write,'(/,"OMP CPU affinity: ",A30)') TRIM(string)
#endif
       write(output_unit_write,'("* * * * * * * * * * * *  * * * * * * *",/)')
    end if
  end subroutine initialize_communications

!! @brief Get the number of threads used.
  function get_num_threads()

     implicit none

     integer :: get_num_threads
#if defined(_OPENMP)
!$OMP PARALLEL
     get_num_threads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#else
     get_num_threads = 1
#endif
  end function

!! @brief Get the maximum number of threads.
   function get_max_threads()

      implicit none

      integer :: get_max_threads
#if defined(_OPENMP)
      get_max_threads = OMP_GET_MAX_THREADS()
#else
      get_max_threads = 1
#endif
   end function

!! @brief Get the current thread number.
!!
!!  Must be called from inside a PARALLEL section.
  function get_thread_number()

     implicit none

     integer :: get_thread_number

#if defined(_OPENMP)
     get_thread_number = OMP_GET_THREAD_NUM()
#else
     get_thread_number = 0
#endif
  end function
end module korc_hpc
