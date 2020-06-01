program korc_test
  use fruit
  use test_mpi
  implicit none
  logical ok
 
  ! initialize fruit
  call init_fruit
  
  ! run tests
  call test_mpi_initialization
  
  ! compile summary and finalize fruit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
    call exit(1)
  endif
  
end program korc_test
