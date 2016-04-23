module constants

implicit none

REAL, PARAMETER :: F_E = 1.602176E-19 !Electron charge in C (absolute value)
REAL, PARAMETER :: F_ME = 9.109382E-31 !Electron mass in kg
REAL, PARAMETER :: F_MP = 1.672621E-27 !Proton mass in kg
REAL, PARAMETER :: F_U = 1.660538E-27 !Atomic mass unit in kg
REAL, PARAMETER :: F_KB = 1.380650E-23 !Boltzmann constant in Joules/Kelvin
REAL, PARAMETER :: F_EPSILON = 8.854E-12 !Vacuum permittivity in C^2/(N*m^2)
REAL, PARAMETER :: F_C = REAL(299792458) !Light speed in m/s
REAL, PARAMETER :: F_PI = REAL(4)*ATAN(REAL(1))
REAL, PARAMETER :: F_MU = (REAL(4)*F_PI)*1E-7 !Vacuum permeability in N/A^2

end module constants



module types

use constants

implicit none

end module types
