#ifdef M3D_C1
!*******************************************************************************
!  @file korc_m3d_c1_interface.f90
!  @brief Interface for the m3d_c1 interpolation library.
!*******************************************************************************

      MODULE korc_m3d_c1
      USE, INTRINSIC :: iso_c_binding
      USE korc_types

      IMPLICIT NONE

      INTEGER (C_INT), PARAMETER :: FIO_SUCCESS        = 0
      INTEGER (C_INT), PARAMETER :: FIO_NO_DATA        = 10006

      INTEGER (C_INT), PARAMETER :: FIO_M3DC1_SOURCE   = 3

      INTEGER (C_INT), PARAMETER :: FIO_TIMESLICE      = 1

      INTEGER (C_INT), PARAMETER :: FIO_SPECIES        = 3

      INTEGER (C_INT), PARAMETER :: FIO_ELECTRON       = 1

      INTEGER (C_INT), PARAMETER :: FIO_DENSITY        = 102
      INTEGER (C_INT), PARAMETER :: FIO_TEMPERATURE    = 103

      INTEGER (C_INT), PARAMETER :: FIO_ELECTRIC_FIELD = 1001
      INTEGER (C_INT), PARAMETER :: FIO_MAGNETIC_FIELD = 1003

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_add_field(icfield, ifield, op, fac)      &
         BIND(C, NAME='fio_add_field')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: icfield
         INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
         INTEGER (C_INT), VALUE, INTENT(IN) :: op
         REAL (C_DOUBLE), VALUE, INTENT(IN) :: fac
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_allocate_search_hint(isrc, hint)         &
         BIND(C, NAME='fio_allocate_search_hint')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
         TYPE (C_PTR), INTENT(OUT)          :: hint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_close_field(ifield)                      &
         BIND(C, NAME='fio_close_field')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_close_series(iseries)                    &
         BIND(C, NAME='fio_close_series')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_close_source(isrc)                       &
         BIND(C, NAME='fio_close_source')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_create_compound_field(ifield)            &
         BIND(C, NAME='fio_create_compound_field')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), INTENT(IN) :: ifield
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_deallocate_search_hint(isrc, hint)       &
         BIND(C, NAME='fio_deallocate_search_hint')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
         TYPE (C_PTR), INTENT(INOUT)        :: hint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_eval_field(ifield, x, v, hint)         &
         BIND(C, NAME='fio_eval_field')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
         REAL (C_DOUBLE), INTENT(IN)        :: x
         REAL (C_DOUBLE), INTENT(OUT)       :: v
         TYPE (C_PTR), VALUE, INTENT(IN)    :: hint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_eval_field_deriv(ifield, x, v, hint)   &
         BIND(C, NAME='fio_eval_field_deriv')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
         REAL (C_DOUBLE), INTENT(IN)        :: x
         REAL (C_DOUBLE), INTENT(OUT)       :: v
         TYPE (C_PTR), VALUE, INTENT(IN)    :: hint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_eval_series(iseries, x, v)             &
         BIND(C, NAME='fio_eval_series')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
         REAL (C_DOUBLE), INTENT(IN)        :: x
         REAL (C_DOUBLE), INTENT(OUT)       :: v
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_get_options(isrc)                        &
         BIND(C, NAME='fio_get_options')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_get_available_fields(isrc, n, f)         &
         BIND(C, NAME='fio_get_available_fields')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN)         :: isrc
         INTEGER (C_INT), INTENT(OUT)               :: n
         INTEGER (C_INT), DIMENSION(:), INTENT(OUT) :: f
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_open_source(itype, filename, handle)     &
         BIND(C, NAME='fio_open_source')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(IN)        :: itype
         CHARACTER (kind=C_CHAR,len=1), INTENT(IN) :: filename
         INTEGER (C_INT), INTENT(OUT)              :: handle
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_get_field(isrc, type, handle)            &
         BIND(C, NAME='fio_get_field')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(in) :: isrc
         INTEGER (C_INT), VALUE, INTENT(IN) :: type
         INTEGER (C_INT), INTENT(OUT)       :: handle
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER (C_INT) FUNCTION fio_set_int_option(iopt, v)                  &
         BIND(C, NAME='fio_set_int_option')
         USE, INTRINSIC :: iso_c_binding

         IMPLICIT NONE

         INTEGER (C_INT), VALUE, INTENT(in) :: iopt
         INTEGER (C_INT), VALUE, INTENT(in) :: v
         END FUNCTION
      END INTERFACE

      CONTAINS

      SUBROUTINE initialize_m3d_c1(params, F, P, spp)

      IMPLICIT NONE

      TYPE(KORC_PARAMS), INTENT(IN)              :: params
      TYPE(FIELDS), INTENT(INOUT)                :: F
      TYPE(PROFILES), INTENT(INOUT)              :: P
      TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp

      INTEGER                                    :: ii
      INTEGER                                    :: pp
      INTEGER                                    :: status
      INTEGER                                    :: isrc

      status = fio_open_source(FIO_M3DC1_SOURCE,                               &
                               TRIM(params%magnetic_field_filename), isrc)

      status = fio_get_options(isrc)
      status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)

      status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%M3D_C1_B)
      status = fio_get_field(isrc, FIO_ELECTRIC_FIELD, F%M3D_C1_E)

      status = fio_set_int_option(FIO_SPECIES, FIO_ELECTRON);
      status = fio_get_field(isrc, FIO_DENSITY, P%M3D_C1_ne);
      status = fio_get_field(isrc, FIO_TEMPERATURE, P%M3D_C1_te);

      do ii=1_idef,params%num_species
        ALLOCATE(spp(ii)%vars%hint(spp(ii)%ppp))

        do pp = 1, spp(ii)%ppp
           status = fio_allocate_search_hint(isrc, spp(ii)%vars%hint(pp))
        end do
      end do

      END SUBROUTINE

      END MODULE
#endif
