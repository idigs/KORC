cmake_minimum_required(VERSION 3.13)


cmake_policy(VERSION 3.13)

if(CMP0074)
   cmake_policy(SET CMP0074 NEW)
endif()

project(fusion-io VERSION 0.1.0 LANGUAGES C CXX Fortran)

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)

