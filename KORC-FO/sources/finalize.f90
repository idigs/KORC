module finalize
use korc_types
#ifdef WITH_MPI
use mpi
#endif
implicit none

contains

subroutine finalize_communications(params)
implicit none
TYPE(KORC_PARAMS), INTENT(OUT) :: params
INTEGER :: ierr

#ifdef WITH_MPI
call MPI_FINALIZE(ierr)
#endif

end subroutine finalize_communications

end module finalize
