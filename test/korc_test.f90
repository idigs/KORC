program korc_test
  use fruit
  use test_io
  use test_hpc
  implicit none
  logical ok

  CHARACTER(MX_STRING_LENGTH) :: path_to_outputs
  
  ! create output file for testing
  call set_paths(path_to_outputs)
  
  ! initialize fruit
  call init_fruit(test_unit_write=test_unit_write)
  
  ! run tests
  write(test_unit_write,*) 'Testing MPI initialization...'
  call test_mpi_initialization

  ! compile summary and finalize fruit
  call fruit_summary(test_unit_write)
  call fruit_summary_xml(path_to_outputs)
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif

  close(test_unit_write)

  
end program korc_test
