module constants
implicit none
	INTEGER, PARAMETER :: MAX_STRING_LENGTH = 1000 ! This value can be changed, beware of truncation errors
	INTEGER, PARAMETER :: default_unit_open = 101

	REAL, PARAMETER :: C_E = 1.602176E-19 !Electron charge in C (absolute value)
	REAL, PARAMETER :: C_ME = 9.109382E-31 !Electron mass in kg
	REAL, PARAMETER :: C_MP = 1.672621E-27 !Proton mass in kg
	REAL, PARAMETER :: C_U = 1.660538E-27 !Atomic mass unit in kg
	REAL, PARAMETER :: C_KB = 1.380650E-23 !Boltzmann constant in Joules/Kelvin
	REAL, PARAMETER :: C_EPSILON = 8.854E-12 !Vacuum permittivity in C^2/(N*m^2)
	REAL, PARAMETER :: C_C = REAL(299792458) !Light speed in m/s
	REAL, PARAMETER :: C_PI = REAL(4)*ATAN(REAL(1))
	REAL, PARAMETER :: C_MU = (REAL(4)*C_PI)*REAL(1E-7) !Vacuum permeability in N/A^2
contains

subroutine show_consts()
	write(6,*) C_PI
end subroutine show_consts

end module constants


module korc_types
use constants
implicit none

TYPE KORC_MPI
    INTEGER :: nmpi ! Number of MPI processes
	INTEGER :: rank ! Rank in WORLD COMMON communicator
	INTEGER :: rank_topo ! Rank in mpi_topo communicator
	INTEGER :: mpi_topo ! MPI communicator for certain topology
END TYPE KORC_MPI

TYPE KORC_PARAMS 
	CHARACTER(MAX_STRING_LENGTH) :: path_to_inputs
	CHARACTER(MAX_STRING_LENGTH) :: path_to_outputs
	INTEGER :: num_omp_threads
	LOGICAL :: restart
	INTEGER :: t_steps
	REAL :: DT
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	INTEGER :: output_cadence
	INTEGER :: num_species

	TYPE(KORC_MPI) :: mpi_params
END TYPE KORC_PARAMS


TYPE PARTICLES
    REAL, DIMENSION(:,:), ALLOCATABLE :: X ! Position (Cartesian)
    REAL, DIMENSION(:,:), ALLOCATABLE :: V ! Velocity
	REAL, DIMENSION(:), ALLOCATABLE :: gamma ! Gamma relativistic
	REAL, DIMENSION(:), ALLOCATABLE :: eta ! Pitch angle
END TYPE PARTICLES

TYPE SPECIES
	TYPE(PARTICLES) :: vars
	REAL :: wc
	REAL :: q
	REAL :: m
	INTEGER :: ppp
	! Here go the parameters for collisions, replenishment, weighting... 
END TYPE SPECIES


TYPE CHARCS_PARAMS
	REAL :: time;
	REAL :: velocity;
	REAL :: length;
	REAL :: mass;
	REAL :: charge;
	REAL :: density;
	REAL :: electric_field;
	REAL :: magnetic_field;
	REAL :: pressure;
	REAL :: temperature;
END TYPE CHARCS_PARAMS


TYPE VEC_FIELD_3D
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: R
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: PHI
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: Z
END TYPE VEC_FIELD_3D

TYPE VEC_FIELD_2D
	REAL, DIMENSION(:,:), ALLOCATABLE :: R
	REAL, DIMENSION(:,:), ALLOCATABLE :: PHI
	REAL, DIMENSION(:,:), ALLOCATABLE :: Z
END TYPE VEC_FIELD_2D

TYPE EMF
	TYPE(VEC_FIELD_3D) :: E
	TYPE(VEC_FIELD_3D) :: B
END TYPE EMF

contains

end module korc_types
