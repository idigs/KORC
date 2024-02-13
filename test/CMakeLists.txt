#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.

add_library (fruit fruit/fruit.f90)

#add_library (fruit_mpi fruit/fruit_mpi.f90)
#target_link_libraries(fruit_mpi PUBLIC fruit )

set (korc_tests "")
list (APPEND korc_tests
     korc_test.f90
     test_io.f90
     test_hpc.f90
)

add_executable (xtest ${korc_tests})
set_property(TARGET xtest PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries (xtest fruit xkorc ${MPI_Fortran_LIBRARIES})
target_link_libraries (xtest ${HDF5_LIBRARIES})
target_link_libraries (xtest ${KORCSRC_LIBRARIES})
