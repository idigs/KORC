module main_mpi
use korc_types
use mpi
implicit none

contains

subroutine initialize_mpi(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	INTEGER :: ierr
	INTEGER :: NDIMS
	INTEGER, DIMENSION(:), ALLOCATABLE :: DIMS
	LOGICAL :: REORDER
	LOGICAL, DIMENSION(:), ALLOCATABLE :: PERIODS
	INTEGER(ip) :: ii ! Iterator

	call MPI_INIT(ierr)
	if (ierr .NE. MPI_SUCCESS) then
		print *,'Error starting MPI program. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD, -10, ierr)
	end if

	call MPI_COMM_SIZE(MPI_COMM_WORLD, params%mpi_params%nmpi, ierr)
	if (ierr .NE. MPI_SUCCESS) then
		print *,'Error getting the size of WORLD COMMON COMMUNICATOR. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD, -10, ierr)
	end if

	! * * * Here a Cartesian topology for MPI is created * * * !
	NDIMS = 1
	REORDER = .FALSE.
	ALLOCATE(DIMS(NDIMS))
	ALLOCATE(PERIODS(NDIMS))
	! This loop isn't necessary but helps to do things more general in the future
	do ii=1,NDIMS
		DIMS(ii) = params%mpi_params%nmpi
		PERIODS(ii) = .TRUE.
	end do

    call MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, PERIODS, REORDER, params%mpi_params%mpi_topo, ierr)
	if (ierr .NE. MPI_SUCCESS) then
		print *,'Error creating new MPI topology. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD, -10, ierr)
	end if

	! * * * Getting the rank of the MPI process in the WORLD COMMON communicator * * * !
	call MPI_COMM_RANK(MPI_COMM_WORLD, params%mpi_params%rank, ierr)
	if (ierr .NE. MPI_SUCCESS) then
		print *,'Error getting MPI rank in COMM_WORLD. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD, -10, ierr)
	end if	

	! * * * Getting the rank of the MPI process in the new topology * * * !
	call MPI_COMM_RANK(params%mpi_params%mpi_topo, params%mpi_params%rank_topo, ierr)
	if (ierr .NE. MPI_SUCCESS) then
		print *,'Error getting MPI rank in KORC topology. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD, -10, ierr)
	end if

	DEALLOCATE(DIMS)
	DEALLOCATE(PERIODS)
end subroutine initialize_mpi


subroutine finalize_mpi(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	INTEGER :: ierr

	call MPI_FINALIZE(ierr)
end subroutine finalize_mpi

end module main_mpi
