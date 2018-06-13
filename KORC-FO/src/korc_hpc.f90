!> @brief KORC module containing subroutines to initilize, control, and to finalize parallel communications.
module korc_hpc
    use korc_types
    use omp_lib
    use mpi

    IMPLICIT NONE

	LOGICAL, PRIVATE :: timed_already = .FALSE. !< Flag to determine if a first call to WMPI_TIME() was made already.
	REAL(rp), PRIVATE :: t1 !< Variable to be used in timing a parallel section of KORC.
	REAL(rp), PRIVATE :: t2 !< Variable to be used in timing a parallel section of KORC.

	PUBLIC :: korc_abort,&
				initialize_mpi,&
				finalize_mpi,&
				initialize_communications,&
				timing_KORC

    CONTAINS

!> @brief Subroutine that terminates the simulation.
!!
!! @param mpierr MPI error status.
subroutine korc_abort()
	INTEGER :: mpierr !< MPI error status

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call MPI_ABORT(MPI_COMM_WORLD, -2000, mpierr)
end subroutine korc_abort

!> @brief Subroutine that initialize MPI communications.
!! @details Through this subroutine the default MPI communicator MPI_COMM_WORLD is initialized. Also, a Cartesian
!!
!! @param[in,out] params Core KORC simulation parameters.
!! @param mpierr MPI error status.
!! @param NDIMS Number of dimensions of non-standard topology. NDIMS=1 for a 1-D MPI topology, NDIMS=2 for a 2-D MPI topology, and NDIMS=3 for a 3-D MPI topology.
!! @param DIMS Dimension of the non-standard MPI topology params::mpi_params::mpi_topo. This is equal to the number of MPI processes in KORC.
!! @param all_mpis_initialized Flag to determine if all the MPI processes were initialized correctly.
!! @param mpi_process_initialized Flago to determine if a given MPI process was initialized correctly.
!! @param REORDER Flag to determine if the new MPI topology params::mpi_params::mpi_topo needs to be re-ordered.
!! @param PERIODS Array of logicals determining what dimensions of the new MPI topology params::mpi_params::mpi_topo are periodic (T) or not (F).
!! @param ii Variable to iterate over different MPI processes.
subroutine initialize_mpi(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
	INTEGER 							:: mpierr
	INTEGER, PARAMETER 					:: NDIMS = 1
	INTEGER, DIMENSION(:), ALLOCATABLE 	:: DIMS
	LOGICAL 							:: all_mpis_initialized = .FALSE.
	LOGICAL 							:: mpi_process_initialized = .FALSE.
	LOGICAL, PARAMETER 					:: REORDER = .FALSE.
	LOGICAL, DIMENSION(:), ALLOCATABLE 	:: PERIODS !< Something here
	INTEGER :: ii

	call MPI_INIT(mpierr)

	if (mpierr .NE. MPI_SUCCESS) then
		write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
		write(6,'(/," ERROR: Initializing MPI. Aborting... ")')
		write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
		call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
	end if

	call MPI_INITIALIZED(mpi_process_initialized,mpierr)

	call MPI_REDUCE(mpi_process_initialized,all_mpis_initialized,1,MPI_LOGICAL,MPI_LAND,0,MPI_COMM_WORLD,mpierr)

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
	call MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, PERIODS, REORDER, params%mpi_params%mpi_topo, mpierr)
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

	if (params%mpi_params%rank.EQ.0) then
		if (all_mpis_initialized) then
			write(6,'(/,"* * * * * COMMUNICATIONS  * * * * *")')
			write(6,'(/,"  MPI communications initialized!  ")')
			write(6,'(/,"  Number of MPI processes: ",I5)') params%mpi_params%nmpi
			write(6,'(/,"* * * * * * * * * * * * * * * * * *")')
		else
			write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
			write(6,'(/," ERROR: MPI not initialized. Aborting... ")')
			write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
			call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
		end if
	end if
!...
end subroutine initialize_mpi

!> @brief Subroutine for timing the execution of any parallel section of KORC.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param individual_runtime Execution time of each MPI process.
!! @param runtime Execution time of KORC defined as the average of the execution times of all MPI processes.
!! @param mpierr MPI error status.
subroutine timing_KORC(params)
	TYPE(KORC_PARAMS), INTENT(IN) 		:: params
	REAL(rp) 							:: individual_runtime
	REAL(rp), DIMENSION(:), ALLOCATABLE :: runtime
	INTEGER 							:: mpierr

	if (timed_already) then
		t2 = MPI_WTIME()

		ALLOCATE(runtime(params%mpi_params%nmpi))

		individual_runtime = t2 - t1

		call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

		call MPI_GATHER(individual_runtime,1,MPI_DOUBLE_PRECISION,runtime,&
				1,MPI_DOUBLE_PRECISION,0_idef, MPI_COMM_WORLD, mpierr)

		if (params%mpi_params%rank .EQ. 0_idef) then
			write(6,'("Timing: ",F30.16," s")') SUM(runtime)/REAL(params%mpi_params%nmpi,rp)
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
		write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
		write(6,'(/," ERROR: MPI not finalized well. MPI process: ",I5)') params%mpi_params%rank
		write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
	end if
end subroutine finalize_mpi


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! ** SUBROUTINES FOR INITIALIZING COMMUNICATIONS  ** !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

!> @brief Subroutine for initializing MPI and open MP communications.
!! @details This subroutine initializes MPI and open MP communications and looks for errors durignt this procces. The system environment
!! variables, which are modified by the user at the moment of running/submitting a KORC simulation, are used to determine the
!! open MP configuration. Some open MP parameters are displayed on the screen/output file.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param mpi_process_finalized Flag indicating whether an individual MPI process was finalized correctly.
!! @param mpierr MPI error status.
subroutine initialize_communications(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
	CHARACTER(MAX_STRING_LENGTH) 		:: string

	call initialize_mpi(params)

!$OMP PARALLEL SHARED(params)
        params%num_omp_threads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

	if (params%mpi_params%rank.EQ.0) then
		write(6,'(/,"* * * * * * * OMP SET-UP * * * * * * *")')
!$OMP PARALLEL
!$OMP MASTER
		write(6,'(/,"OMP threads per MPI process: ",I3)') OMP_GET_NUM_THREADS()
		write(6,'(/,"Cores available per MPI process: ",I3)') OMP_GET_NUM_PROCS()
!$OMP END MASTER
!$OMP END PARALLEL
#ifdef GNU
		call GET_ENVIRONMENT_VARIABLE("OMP_PLACES",string)
		write(6,'(/,"OMP places: ",A30)') TRIM(string)
		call GET_ENVIRONMENT_VARIABLE("GOMP_CPU_AFFINITY",string)
		write(6,'(/,"OMP CPU affinity: ",A30)') TRIM(string)
#endif
		write(6,'("* * * * * * * * * * * *  * * * * * * *",/)')
	end if
end subroutine initialize_communications

end module korc_hpc
