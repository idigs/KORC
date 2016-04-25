module constants
implicit none

	REAL, PARAMETER, PUBLIC :: C_E = 1.602176E-19 !Electron charge in C (absolute value)
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

TYPE KORC_PARAMS
	CHARACTER(LEN=500) :: path_to_outputs
	INTEGER :: numberOfCores
	LOGICAL :: restart
	INTEGER :: t_steps
	REAL :: DT
	CHARACTER(LEN=500) :: magnetic_field_model
	INTEGER :: output_cadence
	INTEGER :: num_species
END TYPE KORC_PARAMS

!TYPE, EXTENDS (KORC_PARAMS) MPI_KORC_PARAMS

!END TYPE MPI_KORC_PARAMS

contains

end module korc_types
