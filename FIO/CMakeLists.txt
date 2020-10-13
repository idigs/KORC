cmake_minimum_required(VERSION 3.10)

cmake_policy(VERSION 3.10)
project(fusion-io VERSION 0.1.0 LANGUAGES C CXX Fortran)

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)

