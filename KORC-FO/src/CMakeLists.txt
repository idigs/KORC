#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.

set (korc_sources "")
list (APPEND korc_sources
    korc_collisions.f90
    korc_fields.f90
    korc_interp.f90
    korc_spatial_distribution.f90
    main.f90
    korc_constants.f90
    korc_finalize.f90
    korc_ppusher.f90
    korc_synthetic_camera.f90
    main_documentation_page.f90
    korc_HDF5.f90
    korc_coords.f90
    korc_hammersley_generator.f90
    korc_profiles.f90
    korc_types.f90
    special_fun_modules.f90
    korc_avalanche.f90
    korc_energy_pdfs.f90
    korc_hpc.f90
    korc_rnd_numbers.f90
    korc_units.f90
    special_fun_subroutines.f90
    korc_binning_diagnostic.f90
    korc_experimental_pdf.f90
    korc_initialize.f90
    korc_simple_equilibrium_pdf.f90
    korc_velocity_distribution.f90
)

add_executable (xkorc ${korc_sources})
target_link_libraries (xkorc ${MPI_Fortran_LIBRARIES})
target_link_libraries (xkorc ${HDF5_LIBRARIES})
target_link_libraries (xkorc ${EZSpline_LIBRARIES})
