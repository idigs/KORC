find_package(HDF5 REQUIRED COMPONENTS Fortran)

list(APPEND m3dc1_headers
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_field_lib.h
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_file.h
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_mesh.h
)
list(APPEND m3dc1_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_field_lib.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_fortran.cpp
)

add_library(m3dc1
    ${m3dc1_sources}
    ${m3dc1_headers}
)

target_include_directories(m3dc1
PUBLIC
    ${HDF5_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(m3dc1 PUBLIC ${HDF5_LIBRARIES})

set_target_properties(m3dc1
    PROPERTIES PUBLIC_HEADER "${m3dc1_headers}"
)

install(
    TARGETS m3dc1
    PUBLIC_HEADER DESTINATION include
    PRIVATE_HEADER DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
