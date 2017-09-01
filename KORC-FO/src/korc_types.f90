module korc_types
	implicit none

! * * * * * * * * * * * * * * * * * * * * !
! * * * Real and integer precisions * * * !
! * * * * * * * * * * * * * * * * * * * * !

	INTEGER, PUBLIC, PARAMETER  :: ip = KIND(INT(1,8)) ! SELECTED_INT_KIND(10) !
	INTEGER, PUBLIC, PARAMETER  :: idef = KIND(1) !
	INTEGER, PUBLIC, PARAMETER  :: rdef = KIND(1.0) !
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
	REAL(rp) :: simulation_time ! Total aimed simulation time in seconds
	REAL(rp) :: snapshot_frequency ! Time between snapshots in seconds
	REAL(rp) :: dt
	REAL(rp) :: time = 0.0_rp
	INTEGER(ip) :: it = 0_ip
	INTEGER(ip) :: t_steps
	INTEGER(ip) :: output_cadence
	INTEGER(ip) :: num_snapshots
	INTEGER :: num_species
	REAL(rp) :: minimum_particle_energy ! Minimum energy of simulated particles in eV
	LOGICAL :: radiation
	LOGICAL :: collisions
	CHARACTER(MAX_STRING_LENGTH) :: collisions_model
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	LOGICAL :: poloidal_flux
	LOGICAL :: axisymmetric
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: outputs_list

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
	REAL(rp), DIMENSION(:), ALLOCATABLE :: g ! Gamma relativistic
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
	CHARACTER(MAX_STRING_LENGTH) :: energy_distribution
	CHARACTER(MAX_STRING_LENGTH) :: pitch_distribution
	REAL(rp) :: Eo
	REAL(rp) :: go
	REAL(rp) :: etao
	REAL(rp), DIMENSION(2) :: Eo_lims
	REAL(rp), DIMENSION(2) :: etao_lims
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
	TYPE(V_FIELD_3D) :: E_3D
	TYPE(V_FIELD_3D) :: B_3D
	TYPE(V_FIELD_2D) :: E_2D
	TYPE(V_FIELD_2D) :: B_2D
	TYPE(MESH) :: X
	INTEGER, DIMENSION(3) :: dims ! dims(NR, NPHI, NZ)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PSIp ! Poloidal flux
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: FLAG2D	
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: FLAG3D

	REAL(rp) :: Bo ! Characteristic magnetic field
	REAL(rp) :: Ro ! Radial position of magnetic axis
	REAL(rp) :: Zo

    CHARACTER(MAX_STRING_LENGTH) :: electric_field_mode
	REAL(rp) :: Eo ! Characteristic electric field
    REAL(rp) :: to
    REAL(rp) :: sig
END TYPE FIELDS


TYPE, PUBLIC :: PROFILES
	! Density
	REAL(rp) :: neo ! Electron density at the magnetic axis
	REAL(rp) :: nea ! Electron density at the plasma edge (if analytical)
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: ne_3D ! ne_3D(R,PHI,Z)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: ne_2D ! ne_2D(R,Z)
	
	!Temperature	
	REAL(rp) :: Teo ! Electron temperature at the magnetic axis
	REAL(rp) :: Tea ! Electron temperature at the plasma edge (if analytical)
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Te_3D ! Te_3D(R,PHI,Z)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Te_2D ! Te_2D(R,Z)
END TYPE PROFILES

PUBLIC :: ALLOCATE_FLUX_ARRAYS,ALLOCATE_2D_FIELDS_ARRAYS,ALLOCATE_3D_FIELDS_ARRAYS,DEALLOCATE_FIELDS_ARRAYS
PRIVATE :: ALLOCATE_V_FIELD_2D,ALLOCATE_V_FIELD_3D

contains

subroutine ALLOCATE_FLUX_ARRAYS(F)
    implicit none
	TYPE(FIELDS), INTENT(INOUT) :: F

	ALLOCATE( F%PSIp(F%dims(1),F%dims(3)) )

	ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))

	ALLOCATE(F%X%R(F%dims(1)))
	ALLOCATE(F%X%Z(F%dims(3)))
end subroutine ALLOCATE_FLUX_ARRAYS


subroutine ALLOCATE_2D_FIELDS_ARRAYS(F)
    implicit none
	TYPE(FIELDS), INTENT(INOUT) :: F

	call ALLOCATE_V_FIELD_2D(F%B_2D,F%dims)
!	call ALLOCATE_V_FIELD_2D(F%E_2D,F%dims)

	ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))
		
	ALLOCATE(F%X%R(F%dims(1)))
	ALLOCATE(F%X%Z(F%dims(3)))
end subroutine ALLOCATE_2D_FIELDS_ARRAYS


subroutine ALLOCATE_3D_FIELDS_ARRAYS(F)
    implicit none
	TYPE(FIELDS), INTENT(INOUT) :: F

	call ALLOCATE_V_FIELD_3D(F%B_3D,F%dims)
!	call ALLOCATE_V_FIELD_3D(F%E_3D,F%dims)	

	ALLOCATE(F%FLAG3D(F%dims(1),F%dims(2),F%dims(3)))
		
	ALLOCATE(F%X%R(F%dims(1)))
	ALLOCATE(F%X%PHI(F%dims(2)))
	ALLOCATE(F%X%Z(F%dims(3)))
end subroutine ALLOCATE_3D_FIELDS_ARRAYS


subroutine ALLOCATE_V_FIELD_2D(F,dims)
    implicit none
	TYPE(V_FIELD_2D), INTENT(INOUT) :: F
	INTEGER, DIMENSION(3), INTENT(IN) :: dims
    
    ALLOCATE(F%R(dims(1),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(3)))
end subroutine ALLOCATE_V_FIELD_2D


subroutine ALLOCATE_V_FIELD_3D(F,dims)
    implicit none
	TYPE(V_FIELD_3D), INTENT(INOUT) :: F
	INTEGER, DIMENSION(3), INTENT(IN) :: dims
    
    ALLOCATE(F%R(dims(1),dims(2),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(2),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(2),dims(3)))
end subroutine ALLOCATE_V_FIELD_3D


subroutine DEALLOCATE_FIELDS_ARRAYS(F)
    implicit none
	TYPE(FIELDS), INTENT(INOUT) :: F

	if (ALLOCATED(F%PSIp)) DEALLOCATE(F%PSIp)
	if (ALLOCATED(F%B_3D%R)) DEALLOCATE(F%B_3D%R)
	if (ALLOCATED(F%B_3D%PHI)) DEALLOCATE(F%B_3D%PHI)
	if (ALLOCATED(F%B_3D%Z)) DEALLOCATE(F%B_3D%Z)
	if (ALLOCATED(F%E_3D%R)) DEALLOCATE(F%E_3D%R)
	if (ALLOCATED(F%E_3D%PHI)) DEALLOCATE(F%E_3D%PHI)
	if (ALLOCATED(F%E_3D%Z)) DEALLOCATE(F%E_3D%Z)

	if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
	if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
	if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)

	if (ALLOCATED(F%FLAG2D)) DEALLOCATE(F%FLAG2D)
	if (ALLOCATED(F%FLAG3D)) DEALLOCATE(F%FLAG3D)
end subroutine DEALLOCATE_FIELDS_ARRAYS

end module korc_types
