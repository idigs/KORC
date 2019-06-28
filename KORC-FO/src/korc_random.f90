MODULE random
    USE, INTRINSIC :: iso_c_binding
    USE korc_types

    IMPLICIT NONE

    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, PRIVATE :: states

    INTERFACE
        TYPE (C_PTR) FUNCTION random_construct(seed) BIND(C, NAME='random_construct')
        USE, INTRINSIC :: iso_c_binding

        IMPLICIT NONE

        INTEGER (C_INT), VALUE :: seed

        END FUNCTION
    END INTERFACE

    INTERFACE
        SUBROUTINE random_set_dist(r, low, high) BIND(C, NAME='random_set_dist')
        USE, INTRINSIC :: iso_c_binding

        IMPLICIT NONE

        TYPE (C_PTR), VALUE    :: r
        REAL (C_DOUBLE), VALUE :: low
        REAL (C_DOUBLE), VALUE :: high
        END SUBROUTINE
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

SUBROUTINE initialize_random(rank)
    USE omp_lib
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: rank

    INTEGER             :: num_threads
    INTEGER             :: thread_num

    num_threads = OMP_GET_MAX_THREADS()
    IF (.NOT. ALLOCATED(states)) THEN
        ALLOCATE(states(0:num_threads - 1))
    END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread_num)
    thread_num = OMP_GET_THREAD_NUM()
    states(thread_num) = random_construct(rank*num_threads + thread_num)
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE destroy_random
    USE omp_lib
    IMPLICIT NONE

!$OMP PARALLEL DEFAULT(SHARED)
    CALL random_destroy(states(OMP_GET_THREAD_NUM()))
!$OMP END PARALLEL

    DEALLOCATE(states)
END SUBROUTINE

FUNCTION get_random()
    USE omp_lib
    IMPLICIT NONE

    REAL(rp) :: get_random

    get_random = random_get_number(states(OMP_GET_THREAD_NUM()))
END FUNCTION

SUBROUTINE get_randoms(nums)
    USE omp_lib
    IMPLICIT NONE

    REAL(rp), DIMENSION(:), INTENT(OUT) :: nums

    INTEGER                             :: i

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    DO i = 0, SIZE(nums)
        nums(i) = get_random()
    END DO
!$OMP END PARALLEL DO
END SUBROUTINE

SUBROUTINE set_random_dist(low, high)
    USE omp_lib
    IMPLICIT NONE

    REAL(rp), INTENT(IN) :: low
    REAL(rp), INTENT(IN) :: high

!$OMP PARALLEL DEFAULT(SHARED)
    CALL random_set_dist(states(OMP_GET_THREAD_NUM()), low, high)
!$OMP END PARALLEL

END SUBROUTINE

!> @brief Gaussian random number generator.
!! @details This function returns a deviate of a Gaussian distribution @f$f_G(x;\mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}} \exp{\left( -(x-\mu)^2/2\sigma^2 \right)}@f$,
!! with mean @f$\mu@f$, and standard deviation @f$\sigma@f$.
!!
!! We use the Inverse Transform Sampling Method for sampling @f$x@f$. With this method we get @f$x = \sqrt{-2\log{(1-y)}}\cos(2\pi z)@f$,
!! where @f$y@f$ and @f$z@f$ are uniform random numbers in the interval @f$[0,1]@f$.
!!
!! @param[in] mu Mean value @f$\mu@f$ of the Gaussian distribution.
!! @param[in] mu Standard deviation @f$\sigma@f$ of the Gaussian distribution.
!! @param random_norm Sampled number @f$x@f$ from the Gaussian distribution @f$f_G(x;\mu,\sigma)@f$.
!! @param rand1 Uniform random number in the interval @f$[0,1]@f$.
!! @param rand2 Uniform random number in the interval @f$[0,1]@f$.
FUNCTION random_norm(mean,sigma)
    USE korc_constants, ONLY: C_PI
    REAL(rp), INTENT(IN) :: mean
    REAL(rp), INTENT(IN) :: sigma
    REAL(rp) :: random_norm
    REAL(rp) :: rand1, rand2

    rand1 = get_random()
    rand2 = get_random()

    random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
END FUNCTION random_norm

END MODULE
