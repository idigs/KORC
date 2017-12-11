module korc_hpc

    use korc_types
    use omp_lib
    use mpi

    implicit none

    contains


subroutine korc_abort()
	implicit none
	INTEGER :: mpierr

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call MPI_ABORT(MPI_COMM_WORLD, -2000, mpierr)
end subroutine korc_abort


subroutine initialize_mpi(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	INTEGER :: mpierr
	INTEGER, PARAMETER :: NDIMS = 1
	INTEGER, DIMENSION(:), ALLOCATABLE :: DIMS
	LOGICAL :: all_mpis_initialized = .FALSE.
	LOGICAL :: mpi_process_initialized = .FALSE.
	LOGICAL, PARAMETER :: REORDER = .FALSE.
	LOGICAL, DIMENSION(:), ALLOCATABLE :: PERIODS
	INTEGER :: ii ! Iterator

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
end subroutine initialize_mpi


subroutine timing_KORC(params,t1,t2)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), INTENT(IN) :: t1, t2
	REAL(rp) :: individual_runtime
	REAL(rp), DIMENSION(:), ALLOCATABLE :: runtime
	INTEGER :: mpierr
	
	ALLOCATE(runtime(params%mpi_params%nmpi))	

	individual_runtime = t2 - t1

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!	write(6,'("MPI: ",I4," Total time: ",F30.16)') params%mpi_params%rank, t2 - t1

	call MPI_GATHER(individual_runtime,1,MPI_DOUBLE_PRECISION,runtime,&
			1,MPI_DOUBLE_PRECISION,0_idef, MPI_COMM_WORLD, mpierr)

	if (params%mpi_params%rank .EQ. 0_idef) then
		write(6,'("The execution time is: ",F30.16)') SUM(runtime)/REAL(params%mpi_params%nmpi,rp)
	end if

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	DEALLOCATE(runtime)	
end subroutine timing_KORC


subroutine finalize_mpi(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	LOGICAL :: mpi_process_finalized = .FALSE.
	INTEGER :: mpierr

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call MPI_FINALIZE(mpierr)

	call MPI_FINALIZED(mpi_process_finalized,mpierr)

	if (.NOT.mpi_process_finalized) then
		write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
		write(6,'(/," ERROR: MPI not finalized well. MPI process: ",I5)') params%mpi_params%rank
		write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
	end if
end subroutine finalize_mpi

end module korc_hpc
