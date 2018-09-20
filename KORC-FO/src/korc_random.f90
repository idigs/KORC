MODULE random
    USE, INTRINSIC :: iso_c_binding
    USE korc_types

    IMPLICIT NONE

    TYPE, PRIVATE :: RANDOM_STATES
        TYPE (C_PTR), DIMENSION(:), ALLOCATABLE :: state
    END TYPE

    TYPE(RANDOM_STATES), PRIVATE :: states

    INTERFACE
        TYPE (C_PTR) FUNCTION random_construct(seed) BIND(C, NAME='random_construct')
        USE, INTRINSIC :: iso_c_binding

        IMPLICIT NONE

        INTEGER (C_INT), VALUE :: seed

        END FUNCTION
    END INTERFACE

    INTERFACE
        REAL (C_DOUBLE) FUNCTION random_get_number(r) BIND(C, NAME='random_get_number')
        USE, INTRINSIC :: iso_c_binding

        IMPLICIT NONE

        TYPE (C_PTR), VALUE :: r

        END FUNCTION
    END INTERFACE

    INTERFACE
        SUBROUTINE random_destroy(r) BIND(C, NAME='random_destroy')
        USE, INTRINSIC :: iso_c_binding

        IMPLICIT NONE

        TYPE (C_PTR), VALUE :: r

        END SUBROUTINE
    END INTERFACE

    PUBLIC :: initialize_random

    CONTAINS

SUBROUTINE initialize_random(seed)
    USE omp_lib
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: seed
    INTEGER             :: num_threads
    INTEGER             :: thread_number

    num_threads = OMP_GET_MAX_THREADS()
    IF (.NOT. ALLOCATED(states%state)) THEN
        ALLOCATE(states%state(0:num_threads - 1))
    END IF

!$OMP PARALLEL
    thread_number = OMP_GET_THREAD_NUM()
    states%state(thread_number) = random_construct(seed + thread_number)
!$OMP END PARALLEL
END SUBROUTINE

FUNCTION get_random()
    USE omp_lib
    IMPLICIT NONE

    REAL(rp)            :: get_random

    get_random = random_get_number(states%state(OMP_GET_THREAD_NUM()))
END FUNCTION

END MODULE
