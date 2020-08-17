module test_hpc
  use mpi
  use fruit
  use korc_types
  implicit none

contains

  subroutine test_mpi_initialization
    use korc_hpc, only : initialize_mpi,finalize_mpi
    TYPE(KORC_PARAMS) :: params
    integer :: ierror
    integer :: size, rank
    integer :: size_k, rank_k

    params%path_to_inputs='TEST'
    
    call initialize_mpi(params)

    size_k=params%mpi_params%nmpi
    rank_k=params%mpi_params%rank

    call MPI_COMM_SIZE (MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierror)

    call assert_equals(size_k,size)
    call assert_equals(rank_k,rank)

    call finalize_mpi(params)
    
  end subroutine test_mpi_initialization

end module test_hpc



