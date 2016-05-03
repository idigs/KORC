module constants
implicit none
	INTEGER, PARAMETER :: ip = SELECTED_INT_KIND(10) !
#ifdef DOUBLE_PRECISION
	INTEGER, PARAMETER :: rp = KIND(0.d0) ! Double precision kind
#else
	INTEGER, PARAMETER :: rp = KIND(1.0) ! Double precision kind
#endif
	INTEGER, PARAMETER :: sp = kind(1.0) ! Single precision kind
	INTEGER, PARAMETER :: MAX_STRING_LENGTH = 1000_ip ! This value can be changed, beware of truncation errors
	INTEGER, PARAMETER :: default_unit_open = 101_ip
	INTEGER, PARAMETER :: default_unit_write = 201_ip

	REAL(rp), PARAMETER :: C_E = 1.602176E-19_rp !Electron charge in C (absolute value)
	REAL(rp), PARAMETER :: C_ME = 9.109382E-31_rp !Electron mass in kg
	REAL(rp), PARAMETER :: C_MP = 1.672621E-27_rp !Proton mass in kg
	REAL(rp), PARAMETER :: C_U = 1.660538E-27_rp !Atomic mass unit in kg
	REAL(rp), PARAMETER :: C_KB = 1.380650E-23_rp !Boltzmann constant in Joules/Kelvin
	REAL(rp), PARAMETER :: C_EPSILON = 8.854E-12_rp !Vacuum permittivity in C^2/(N*m^2)
	REAL(rp), PARAMETER :: C_C = 299792458.0_rp !Light speed in m/s
	REAL(rp), PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp)
	REAL(rp), PARAMETER :: C_MU = 4.0_rp*C_PI*1E-7_rp !Vacuum permeability in N/A^2
contains

subroutine show_consts()
	write(6,*) C_PI
end subroutine show_consts

end module constants



module korc_types
use constants
implicit none

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


TYPE, PUBLIC :: KORC_PARAMS
	CHARACTER(MAX_STRING_LENGTH) :: path_to_inputs
	CHARACTER(MAX_STRING_LENGTH) :: path_to_outputs
	INTEGER :: num_omp_threads
	LOGICAL :: restart
	INTEGER(ip) :: t_steps
	INTEGER(ip) :: output_cadence
	INTEGER(ip) :: num_snapshots
	REAL(rp) :: dt
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	INTEGER :: num_species
	INTEGER :: pic_algorithm

	TYPE(KORC_MPI) :: mpi_params
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
	REAL(rp), DIMENSION(:), ALLOCATABLE :: kappa ! Curvature
	REAL(rp), DIMENSION(:), ALLOCATABLE :: tau ! Torsion
END TYPE PARTICLES


TYPE, PUBLIC :: SPECIES
	TYPE(PARTICLES) :: vars
	LOGICAL :: runaway
	REAL(rp) :: Eo
	REAL(rp) :: wc
	REAL(rp) :: q
	REAL(rp) :: m
	INTEGER :: ppp
	! Here go the parameters for collisions, replenishment, weighting... 
END TYPE SPECIES


TYPE, PUBLIC :: CHARCS_PARAMS
	REAL(rp) :: time;
	REAL(rp) :: velocity;
	REAL(rp) :: length;
	REAL(rp) :: mass;
	REAL(rp) :: charge;
	REAL(rp) :: density;
	REAL(rp) :: electric_field;
	REAL(rp) :: magnetic_field;
	REAL(rp) :: pressure;
	REAL(rp) :: temperature;
END TYPE CHARCS_PARAMS


TYPE, PRIVATE :: A_FIELD
	REAL(rp) :: Bo
	REAL(rp) :: a
	REAL(rp) :: Ro
	REAL(rp) :: qa
	REAL(rp) :: co
	REAL(rp) :: lambda
	REAL(rp) :: Bpo
END TYPE A_FIELD


TYPE, PUBLIC :: FIELDS
	TYPE(A_FIELD) :: AB
	TYPE(V_FIELD_3D) :: E
	TYPE(V_FIELD_3D) :: B
	REAL(rp) :: Bo ! Characteristic magnetic field
	INTEGER, DIMENSION(3) :: dim ! dim(NR, NPHI, NZ)
END TYPE FIELDS

contains

end module korc_types
