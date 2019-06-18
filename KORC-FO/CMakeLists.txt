cmake_minimum_required (VERSION 3.14)

project (korc_proj C CXX Fortran)

#cmake_policy (VERSION 3.10)
#cmake_policy(SET CMP0074 NEW)

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)

#  Fortran specific settings.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules)

#  Set definitions.
add_definitions (-DDOUBLE_PRECISION)
add_definitions (-DPARALLEL_RANDOM)

#  Configure OpenMP
find_package (OpenMP)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")

#  Configure MPI
find_package (MPI REQUIRED)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
include_directories (${MPI_Fortran_INCLUDE_PATH})

#  Configure HDF5
find_package (HDF5 REQUIRED COMPONENTS Fortran)
include_directories (${HDF5_INCLUDE_DIRS})

#  Configure EZspline
include_directories (${EZSpline_INCLUDE_PATH})
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EZSpline_LINK_FLAGS}")

add_subdirectory (src)
