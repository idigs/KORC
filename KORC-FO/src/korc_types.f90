module korc_types
	implicit none

! * * * * * * * * * * * * * * * * * * * * !
! * * * Real and integer precisions * * * !
! * * * * * * * * * * * * * * * * * * * * !

	INTEGER, PUBLIC, PARAMETER  :: ip = KIND(INT(1,8)) ! SELECTED_INT_KIND(10) !
	INTEGER, PUBLIC, PARAMETER  :: idef = KIND(1) !
#ifdef DOUBLE_PRECISION
	INTEGER, PUBLIC, PARAMETER :: rp = KIND(0.d0) ! Double precision
#elif SINGLE_PRECISION
	INTEGER, PUBLIC, PARAMETER :: rp = KIND(1.0) ! Single precision
#endif
	REAL(rp), PUBLIC, PARAMETER :: korc_zero = 1.0E-15

! * * * * * * * * * * * * * * * * * * * * !
! * * * Real and integer precisions * * * !
! * * * * * * * * * * * * * * * * * * * * !

	INTEGER, PUBLIC, PARAMETER :: MAX_STRING_LENGTH = 1000 ! This value can be changed, beware of string truncation
	INTEGER, PUBLIC, PARAMETER :: default_unit_open = 101
	INTEGER, PUBLIC, PARAMETER :: default_unit_write = 201

TYPE, PUBLIC :: KORC_STRING
	CHARACTER(MAX_STRING_LENGTH) :: str
END TYPE KORC_STRING

! * * * * * * * * * * * * * * * * * * * * !
! * * * Basic korc array structures * * * !
! * * * * * * * * * * * * * * * * * * * * !

TYPE, PRIVATE :: V_FIELD_3D
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: R
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PHI
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Z
END TYPE V_FIELD_3D

TYPE, PRIVATE :: V_FIELD_2D
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: R
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PHI
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z
END TYPE V_FIELD_2D

! * * * * * * * * * * * * * * * * * * * * !
! * * * Basic korc array structures * * * !
! * * * * * * * * * * * * * * * * * * * * !

TYPE, PRIVATE :: KORC_MPI
    INTEGER :: nmpi ! Number of MPI processes
	INTEGER :: rank ! Rank in WORLD COMMON communicator
	INTEGER :: rank_topo ! Rank in mpi_topo communicator
	INTEGER :: mpi_topo ! MPI communicator for certain topology
END TYPE KORC_MPI


TYPE, PUBLIC :: CHARCS_PARAMS
	REAL(rp) :: time
	REAL(rp) :: time_r
	REAL(rp) :: velocity
	REAL(rp) :: length
	REAL(rp) :: mass
	REAL(rp) :: charge
	REAL(rp) :: density
	REAL(rp) :: Eo ! Characteristic electric field
	REAL(rp) :: Bo ! Characteristic magnetic field
	REAL(rp) :: energy
	REAL(rp) :: pressure
	REAL(rp) :: temperature
END TYPE CHARCS_PARAMS


TYPE, PUBLIC :: KORC_PARAMS
	CHARACTER(MAX_STRING_LENGTH) :: path_to_inputs
	CHARACTER(MAX_STRING_LENGTH) :: path_to_outputs
	INTEGER :: num_omp_threads
	LOGICAL :: restart
	REAL(rp) :: dt
	REAL(rp) :: time
	INTEGER(ip) :: t_steps
	INTEGER(ip) :: output_cadence
	INTEGER(ip) :: num_snapshots
	INTEGER :: num_species
	INTEGER :: num_impurity_species
	INTEGER :: pic_algorithm
	LOGICAL :: radiation
	LOGICAL :: collisions
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	LOGICAL :: poloidal_flux
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename

	TYPE(KORC_MPI) :: mpi_params
	TYPE(CHARCS_PARAMS) :: cpp
END TYPE KORC_PARAMS


TYPE, PUBLIC :: PARTICLES
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X ! Position (Cartesian) dim(X) = (3,num_particles)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: V ! Velocity
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Rgc ! Guiding-center position (Cartesian)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Y ! Position in alternative coordinate system, i.e. cylindrical or toroidal coordinates.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: E ! Auxiliar vector for fields interpolations
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: B ! Auxiliar vector for fields interpolations
	REAL(rp), DIMENSION(:), ALLOCATABLE :: gamma ! Gamma relativistic
	REAL(rp), DIMENSION(:), ALLOCATABLE :: eta ! Pitch angle
	REAL(rp), DIMENSION(:), ALLOCATABLE :: mu ! Instantaneous magnetic moment
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Prad ! Radiated power (in Watts/electron)
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Pin ! Input power (in Watts/electron)
	INTEGER, DIMENSION(:), ALLOCATABLE :: flag
	REAL(rp), DIMENSION(:), ALLOCATABLE :: AUX
END TYPE PARTICLES


TYPE, PUBLIC :: SPECIES
	TYPE(PARTICLES) :: vars
	LOGICAL :: runaway
	REAL(rp) :: Eo
	REAL(rp) :: gammao
	REAL(rp) :: etao
	REAL(rp) :: wc
	REAL(rp) :: wc_r
	REAL(rp) :: q
	REAL(rp) :: m
	INTEGER :: ppp

	! Parameters for initializing spatial distribution
	REAL(rp) :: Ro
	REAL(rp) :: Zo
	REAL(rp) :: r
END TYPE SPECIES


TYPE, PRIVATE :: A_FIELD
	REAL(rp) :: Bo
	REAL(rp) :: a
	REAL(rp) :: Ro
	REAL(rp) :: qa
	REAL(rp) :: qo
	REAL(rp) :: lambda
	REAL(rp) :: Bpo
END TYPE A_FIELD


TYPE, PRIVATE :: MESH
	REAL(rp), DIMENSION(:), ALLOCATABLE :: R
	REAL(rp), DIMENSION(:), ALLOCATABLE :: PHI
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Z
END TYPE MESH


TYPE, PUBLIC :: FIELDS
	TYPE(A_FIELD) :: AB
	TYPE(V_FIELD_3D) :: E
	TYPE(V_FIELD_3D) :: B
	TYPE(MESH) :: X
	INTEGER, DIMENSION(3) :: dims ! dims(NR, NPHI, NZ)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PSIp ! Poloidal flux

	REAL(rp) :: Bo ! Characteristic magnetic field
	REAL(rp) :: Ro ! Radial position of magnetic axis

    CHARACTER(MAX_STRING_LENGTH) :: electric_field_mode
	REAL(rp) :: Eo ! Characteristic electric field
    REAL(rp) :: to
    REAL(rp) :: sig
END TYPE FIELDS


TYPE, PUBLIC :: COLLISION_PARAMS
	INTEGER :: num_impurity_species
	REAL(rp) :: Te ! Background electron temperature in eV
	REAL(rp) :: ne ! Background electron density in 1/m^3
	REAL(rp) :: nH ! Background proton density in 1/m^3
	REAL(rp) :: nef ! Free electron density in 1/m^3
	REAL(rp), DIMENSION(:), ALLOCATABLE :: neb ! Bound electron density in 1/m^3
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zi ! Atomic number of (majority) background ions
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zj ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(:), ALLOCATABLE :: nz ! Impurity densities
	REAL(rp), DIMENSION(:), ALLOCATABLE :: IZj ! Ionization energy of impurity in eV
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Ee_IZj ! me*c^2/IZj dimensionless parameter

	REAL(rp) :: rD ! Debye length
	REAL(rp) :: re ! Classical electron radius
END TYPE COLLISION_PARAMS


PUBLIC :: ALLOCATE_FIELDS_ARRAYS, DEALLOCATE_FIELDS_ARRAYS
PRIVATE :: ALLOCATE_V_FIELD_3D, DEALLOCATE_V_FIELD_3D

contains

subroutine ALLOCATE_FIELDS_ARRAYS(F,opt)
    implicit none
	TYPE(FIELDS), INTENT(INOUT) :: F
	LOGICAL, INTENT(IN) :: opt

	if (opt) then ! Using only magnetic poloidal flux
		ALLOCATE( F%PSIp(F%dims(1), F%dims(3)) )

		ALLOCATE(F%X%R(F%dims(1)))
		ALLOCATE(F%X%Z(F%dims(3)))
	else ! Using cylindrical components of magnetic field
		call ALLOCATE_V_FIELD_3D(F%B,F%dims)
	!	call ALLOCATE_V_FIELD_3D(F%E,F%dims)	
		
		ALLOCATE(F%X%R(F%dims(1)))
		ALLOCATE(F%X%PHI(F%dims(2)))
		ALLOCATE(F%X%Z(F%dims(3)))
	end if
end subroutine ALLOCATE_FIELDS_ARRAYS


subroutine DEALLOCATE_FIELDS_ARRAYS(F)
    implicit none
	TYPE(FIELDS), INTENT(INOUT) :: F

	if (ALLOCATED(F%PSIp)) DEALLOCATE(F%PSIp)
	if (ALLOCATED(F%B%R)) DEALLOCATE(F%B%R)
	if (ALLOCATED(F%B%PHI)) DEALLOCATE(F%B%PHI)
	if (ALLOCATED(F%B%Z)) DEALLOCATE(F%B%Z)
	if (ALLOCATED(F%E%R)) DEALLOCATE(F%E%R)
	if (ALLOCATED(F%E%PHI)) DEALLOCATE(F%E%PHI)
	if (ALLOCATED(F%E%Z)) DEALLOCATE(F%E%Z)

	if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
	if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
	if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)
end subroutine DEALLOCATE_FIELDS_ARRAYS


subroutine ALLOCATE_V_FIELD_3D(F,dims)
    implicit none
	TYPE(V_FIELD_3D), INTENT(INOUT) :: F
	INTEGER, DIMENSION(3), INTENT(IN) :: dims
    
    ALLOCATE(F%R(dims(1),dims(2),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(2),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(2),dims(3)))
end subroutine ALLOCATE_V_FIELD_3D


subroutine DEALLOCATE_V_FIELD_3D(F)
    implicit none
	TYPE(V_FIELD_3D), INTENT(INOUT) :: F
    
    DEALLOCATE(F%R)
    DEALLOCATE(F%PHI)
    DEALLOCATE(F%Z)
end subroutine DEALLOCATE_V_FIELD_3D


end module korc_types
