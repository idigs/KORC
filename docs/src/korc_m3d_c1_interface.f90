#ifdef M3D_C1
!*******************************************************************************
!  @file korc_m3d_c1_interface.f90
!  @brief Interface for the m3d_c1 interpolation library.
!*******************************************************************************

MODULE korc_m3d_c1
  USE, INTRINSIC :: iso_c_binding
  USE korc_types
  USE korc_input

  IMPLICIT NONE

  INTEGER (C_INT), PARAMETER :: FIO_SUCCESS        = 0
  INTEGER (C_INT), PARAMETER :: FIO_OUT_OF_BOUNDS  = 10002
  INTEGER (C_INT), PARAMETER :: FIO_NO_DATA        = 10006

  INTEGER (C_INT), PARAMETER :: FIO_M3DC1_SOURCE   = 3

  INTEGER (C_INT), PARAMETER :: FIO_TIMESLICE      = 1

  INTEGER (C_INT), PARAMETER :: FIO_SPECIES        = 3

  INTEGER (C_INT), PARAMETER :: FIO_ELECTRON       = 1
  INTEGER (C_INT), PARAMETER :: FIO_MAIN_ION       = -1

  INTEGER (C_INT), PARAMETER :: FIO_DENSITY        = 102
  INTEGER (C_INT), PARAMETER :: FIO_TEMPERATURE    = 103

  INTEGER (C_INT), PARAMETER :: FIO_ELECTRIC_FIELD = 1001
  INTEGER (C_INT), PARAMETER :: FIO_MAGNETIC_FIELD = 1003
  INTEGER (C_INT), PARAMETER :: FIO_VECTOR_POTENTIAL = 1002

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_add_field(icfield, ifield, op, fac)      &
          BIND(C, NAME='fio_add_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: icfield
       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       INTEGER (C_INT), VALUE, INTENT(IN) :: op
       REAL (C_DOUBLE), VALUE, INTENT(IN) :: fac
     END FUNCTION fio_add_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_allocate_search_hint(isrc, hint)         &
          BIND(C, NAME='fio_allocate_search_hint')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
       TYPE (C_PTR), INTENT(OUT)          :: hint
     END FUNCTION fio_allocate_search_hint
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_field(ifield)                      &
          BIND(C, NAME='fio_close_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
     END FUNCTION fio_close_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_series(iseries)                    &
          BIND(C, NAME='fio_close_series')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
     END FUNCTION fio_close_series
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_source(isrc)                       &
          BIND(C, NAME='fio_close_source')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
     END FUNCTION fio_close_source
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_create_compound_field(ifield)            &
          BIND(C, NAME='fio_create_compound_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), INTENT(IN) :: ifield
     END FUNCTION fio_create_compound_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_deallocate_search_hint(isrc, hint)       &
          BIND(C, NAME='fio_deallocate_search_hint')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
       TYPE (C_PTR), INTENT(INOUT)        :: hint
     END FUNCTION fio_deallocate_search_hint
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
     END FUNCTION fio_eval_field
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
     END FUNCTION fio_eval_field_deriv
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_series(iseries, x, v)             &
          BIND(C, NAME='fio_eval_series')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
     END FUNCTION fio_eval_series
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_options(isrc)                        &
          BIND(C, NAME='fio_get_options')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
     END FUNCTION fio_get_options
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_available_fields(isrc, n, f)         &
          BIND(C, NAME='fio_get_available_fields')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN)         :: isrc
       INTEGER (C_INT), INTENT(OUT)               :: n
       INTEGER (C_INT), DIMENSION(:), INTENT(OUT) :: f
     END FUNCTION fio_get_available_fields
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_open_source(itype, filename, handle)     &
          BIND(C, NAME='fio_open_source')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN)        :: itype
       CHARACTER (kind=C_CHAR,len=1), INTENT(IN) :: filename
       INTEGER (C_INT), INTENT(OUT)              :: handle
     END FUNCTION fio_open_source
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_field(isrc, type, handle)            &
          BIND(C, NAME='fio_get_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(in) :: isrc
       INTEGER (C_INT), VALUE, INTENT(IN) :: type
       INTEGER (C_INT), INTENT(INOUT)     :: handle
     END FUNCTION fio_get_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_set_int_option(iopt, v)                  &
          BIND(C, NAME='fio_set_int_option')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(in) :: iopt
       INTEGER (C_INT), VALUE, INTENT(in) :: v
     END FUNCTION fio_set_int_option
  END INTERFACE

CONTAINS

  SUBROUTINE initialize_m3d_c1(params, F, P, spp)

    TYPE(KORC_PARAMS), INTENT(IN)           :: params
    TYPE(FIELDS), INTENT(INOUT)                :: F
    TYPE(PROFILES), INTENT(INOUT)              :: P
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp

    INTEGER                                    :: ii
    INTEGER                                    :: pp
    INTEGER                                    :: status
    INTEGER                                    :: isrc


    F%Efield = Efield
    F%PSIp_lim=PSIp_lim
    F%PSIp_0=PSIp_0
    
    status = fio_open_source(FIO_M3DC1_SOURCE,           &
         TRIM(params%magnetic_field_filename)            &
         &                         // C_NULL_CHAR, isrc)

    status = fio_get_options(isrc)
    status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)

    status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%M3D_C1_B)
    status = fio_get_field(isrc, FIO_ELECTRIC_FIELD, F%M3D_C1_E)
    status = fio_get_field(isrc, FIO_VECTOR_POTENTIAL, F%M3D_C1_A)

    if (.not.F%Efield) F%M3D_C1_E=-1

    status = fio_set_int_option(FIO_SPECIES, FIO_ELECTRON);
    status = fio_get_field(isrc, FIO_DENSITY, P%M3D_C1_ne);
    status = fio_get_field(isrc, FIO_TEMPERATURE, P%M3D_C1_te);

    status = fio_set_int_option(FIO_SPECIES, FIO_MAIN_ION);
    status = fio_get_field(isrc, FIO_DENSITY, P%M3D_C1_ni);
    
    !  Hardcode Bo to one for now until a better method of determining the a
    !  characteristic magnetic field value.
    F%Bo = 1.0
    F%Eo = 1.0
    F%Ro = 1.0
    F%Zo = 1.0        

    do ii = 1, params%num_species

       do pp = 1, spp(ii)%ppp
          status = fio_allocate_search_hint(isrc, spp(ii)%vars%hint(pp))
          !spp(ii)%vars%hint(pp)=c_null_ptr
       end do

       spp(ii)%vars%cart = .false.
    end do

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,*) 'Calculate B',F%M3D_C1_B
       write(output_unit_write,*) 'Calculate E',F%M3D_C1_E
       write(output_unit_write,*) 'Calculate A',F%M3D_C1_A
       write(output_unit_write,*) 'Calculate ne',P%M3D_C1_ne
       write(output_unit_write,*) 'Calculate Te',P%M3D_C1_te
       write(output_unit_write,*) 'Calculate ni',P%M3D_C1_ni
    end if
       
  END SUBROUTINE initialize_m3d_c1
  
  SUBROUTINE initialize_m3d_c1_imp(params, P,num_imp)

    TYPE(KORC_PARAMS), INTENT(IN)           :: params
    TYPE(PROFILES), INTENT(INOUT)              :: P
    INTEGER, INTENT(IN)				:: num_imp

    INTEGER                                    :: ii
    INTEGER                                    :: status
    INTEGER                                    :: isrc

    
    status = fio_open_source(FIO_M3DC1_SOURCE,           &
         TRIM(params%magnetic_field_filename)            &
         &                         // C_NULL_CHAR, isrc)

    status = fio_get_options(isrc)
    status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)

    ALLOCATE(P%M3D_C1_nimp(num_imp))
    
    do ii=1,num_imp
       status = fio_set_int_option(FIO_SPECIES, &
            FIO_MAKE_SPECIES(40, 18, 18+1-ii));
       status = fio_get_field(isrc, FIO_DENSITY, P%M3D_C1_nimp(ii));
    end do

    if (params%mpi_params%rank .EQ. 0) then
       do ii=1,num_imp
          write(output_unit_write,*) 'Calculate nimp_',ii,P%M3D_C1_nimp(ii)
       end do
    end if
       
  END SUBROUTINE initialize_m3d_c1_imp

  FUNCTION FIO_MAKE_SPECIES(m, p, e)
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(IN) :: e
    INTEGER(C_INT) :: FIO_MAKE_SPECIES
    
    FIO_MAKE_SPECIES = e + p*256 + (m-p)*65536
  END FUNCTION FIO_MAKE_SPECIES


END MODULE korc_m3d_c1
#endif
