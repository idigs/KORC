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
	CHARACTER(LEN=100) :: path_to_outputs
	INTEGER :: numberOfCores
	LOGICAL :: restart
	INTEGER :: t_steps
	REAL :: DT
	CHARACTER(LEN=100) :: magnetic_field_model
	INTEGER :: output_cadence
	INTEGER :: num_species
END TYPE KORC_PARAMS

contains

subroutine load_korc_params(path_to_input,params)
	implicit none
	CHARACTER(15), INTENT(IN) :: path_to_input
	TYPE (KORC_PARAMS), INTENT(OUT) :: params

	CHARACTER(LEN=100) :: path_to_outputs
	INTEGER :: numberOfCores
	LOGICAL :: restart
	INTEGER :: t_steps
	REAL :: DT
	CHARACTER(LEN=100) :: magnetic_field_model
	INTEGER :: output_cadence
	INTEGER :: num_species

	NAMELIST /input_parameters/ path_to_outputs,numberOfCores,restart,&
			t_steps,DT,magnetic_field_model,output_cadence,num_species
	
	open(unit=101,file=path_to_input,status='OLD',form='formatted')
	read(101,nml=input_parameters)
	close(101)
	
	params%path_to_outputs = TRIM(path_to_outputs)
	params%numberOfCores = numberOfCores
	params%restart = restart
	params%t_steps = t_steps
	params%DT = DT
	params%magnetic_field_model = TRIM(magnetic_field_model)
	params%output_cadence = output_cadence
	params%num_species = num_species


end subroutine load_korc_params

end module korc_types
