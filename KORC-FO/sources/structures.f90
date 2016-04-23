module structures

use types

implicit none

type KORC_PARAMS
	 ! Control parameters for the simulation

    CHARACTER(:), ALLOCATABLE :: path_to_outputs ! Path to save the outputs. It must point to the directory where the folder outputFiles is.

	INTEGER :: numberOfCores ! The number of cores to be used in the openMP routines.

    LOGICAL :: restart
    INTEGER :: t_steps
    REAL DT; ! Time step
    CHARACTER(:), ALLOCATABLE :: magnetic_field_model;

    INTEGER output_cadence

    INTEGER :: num_species
end type KORC_PARAMS

end module structures
