module korc_constants

	USE korc_types

	IMPLICIT NONE

	REAL(rp), PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp)
	REAL(rp), PARAMETER :: C_E = 1.602176E-19_rp !Electron charge in C (absolute value)
	REAL(rp), PARAMETER :: C_ME = 9.109382E-31_rp !Electron mass in kg
	REAL(rp), PARAMETER :: C_MP = 1.672621E-27_rp !Proton mass in kg
	REAL(rp), PARAMETER :: C_U = 1.660538E-27_rp !Atomic mass unit in kg
	REAL(rp), PARAMETER :: C_KB = 1.380650E-23_rp !Boltzmann constant in Joules/Kelvin
	REAL(rp), PARAMETER :: C_C = 299792458.0_rp !Light speed in m/s
	REAL(rp), PARAMETER :: C_MU = 4.0_rp*C_PI*1E-7_rp !Vacuum permeability in N/A^2
	REAL(rp), PARAMETER :: C_E0 = 1.0_rp/(C_MU*C_C**2) !8.854E-12_rp !Vacuum permittivity in C^2/(N*m^2)
	REAL(rp), PARAMETER :: C_Ke = 1.0_rp/(4.0_rp*C_PI*C_E0) ! Coulomb constant in N*m^2/C^2
	REAL(rp), PARAMETER :: C_RE = 	C_E**2/( 4.0_rp*C_PI*C_E0*C_ME*C_C**2 ) ! Classical electron radius
	REAL(rp), PARAMETER :: C_h = 6.6261E-34_rp ! Planck constant in Joules*s
end module korc_constants
