module constants
implicit none
	INTEGER, PARAMETER :: rp = KIND(0.d0) ! Double precision kind
	INTEGER, PARAMETER :: sp = kind(1.0) ! Single precision kind
	INTEGER, PARAMETER :: MAX_STRING_LENGTH = 1000 ! This value can be changed, beware of truncation errors
	INTEGER, PARAMETER :: default_unit_open = 101

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
	REAL(rp) :: DT
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	INTEGER :: output_cadence
	INTEGER :: num_species

	TYPE(KORC_MPI) :: mpi_params
END TYPE KORC_PARAMS


TYPE PARTICLES
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X ! Position (Cartesian)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: V ! Velocity
	REAL(rp), DIMENSION(:), ALLOCATABLE :: gamma ! Gamma relativistic
	REAL(rp), DIMENSION(:), ALLOCATABLE :: eta ! Pitch angle
END TYPE PARTICLES

TYPE SPECIES
	TYPE(PARTICLES) :: vars
	REAL(rp) :: wc
	REAL(rp) :: q
	REAL(rp) :: m
	INTEGER :: ppp
	! Here go the parameters for collisions, replenishment, weighting... 
END TYPE SPECIES


TYPE CHARCS_PARAMS
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


TYPE VEC_FIELD_3D
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: R
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PHI
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Z
END TYPE VEC_FIELD_3D

TYPE VEC_FIELD_2D
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: R
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PHI
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z
END TYPE VEC_FIELD_2D

TYPE EMF
	TYPE(VEC_FIELD_3D) :: E
	TYPE(VEC_FIELD_3D) :: B
	INTEGER, DIMENSION(3) :: dim ! dim(NR, NPHI, NZ)
END TYPE EMF

contains

end module korc_types
