!> @brief Module containing the definition of KORC derived types and KORC variables, the building blocks of the code.
module korc_types
#ifdef M3D_C1
    USE, INTRINSIC :: iso_c_binding
#endif
	implicit none

! * * * * * * * * * * * * * * * * * * * * !
! * * * Real and integer precisions * * * !
! * * * * * * * * * * * * * * * * * * * * !

	INTEGER, PUBLIC, PARAMETER 	:: is = KIND(INT(1,1)) !< Definition of 1 Byte (8 bits) Fortran KORC integer type.
	INTEGER, PUBLIC, PARAMETER 	:: ip = KIND(INT(1,8)) !< Definition of 8 Bytes (64 bits) Fortran KORC integer type.
	INTEGER, PUBLIC, PARAMETER 	:: idef = KIND(1) !< Definition of the default KORC integer type on the system where KORC is compiled.
	INTEGER, PUBLIC, PARAMETER 	:: rdef = KIND(1.0) !< Definition of the default KORC real type on the system where KORC is compiled.
#ifdef DOUBLE_PRECISION
	INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(0.d0) !< Definition of the KORC double precision real type.
#elif SINGLE_PRECISION
	INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(1.0) !< Definition of the KORC single precision real type.
#endif
	REAL(rp), PUBLIC, PARAMETER :: korc_zero = 1.0E-15 !< Definition of the zero in KORC.

! * * * * * * * * * * * * * * * * * * * * !
! * * * Real and integer precisions * * * !
! * * * * * * * * * * * * * * * * * * * * !

	INTEGER, PUBLIC, PARAMETER 	:: MAX_STRING_LENGTH = 1000 !< Default length of a KORC_STRING variable.
	INTEGER, PUBLIC, PARAMETER 	:: default_unit_open = 101 !< Default file unit for opening and reading from an external text file.
	INTEGER, PUBLIC, PARAMETER 	:: default_unit_write = 201 !< Default file unit for opening and writing to external an external text file.


TYPE, PUBLIC :: KORC_STRING !< @brief KORC string type of length MAX_STRING_LENGTH.
	CHARACTER(MAX_STRING_LENGTH) :: str
END TYPE KORC_STRING

! * * * * * * * * * * * * * * * * * * * * !
! * * * Basic korc array structures * * * !
! * * * * * * * * * * * * * * * * * * * * !

TYPE, PUBLIC :: V_FIELD_3D
!< @brief KORC 3-D vector field type
!! @details This KORC type represents a 3-D vector field varible in cylindrical coordinates. For example, this could be the 3-D magnetic
!! field, which can be written as @f$\mathbf{B}(R,\phi,Z) = B_R(R,\phi,Z) \hat{R} + B_\phi(R,\phi,Z) \hat{\phi} + B_Z(R,\phi,Z) \hat{Z}@f$.
!! All the members (components) of the V_FIELD_3D type follow the following index convention:
!! (@f\mathbf$ index,@f$\phi@f$ index,@f$Z@f$ index)

	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: R !< @f$R @f$ component of the vector field variable.
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PHI !< @f$\phi @f$ component of the vector field variable.
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Z !< @f$Z @f$ component of the vector field variable.
END TYPE V_FIELD_3D

TYPE, PUBLIC :: V_FIELD_2D
!< @brief KORC 2-D vector field type
!! @details This KORC type represents a 2-D vector field varible in cylindrical coordinates. For example, this could be the magnetic
!! field in an axisymmetric plasma, which can be written as @f$\mathbf{B}(R,Z) = B_R(R,Z) \hat{R} + B_\phi(R,Z) \hat{\phi} + B_Z(R,Z)
!! \hat{Z}@f$.
!! All the members (components) of the V_FIELD_2D type follow the following index convention:
!! (@f$R@f$ index,@f$Z@f$ index).
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: R !< @f$R @f$ component of the vector field variable.
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PHI !< @f$\phi @f$ component of the vector field variable.
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z !< @f$Z @f$ component of the vector field variable.
END TYPE V_FIELD_2D

! * * * * * * * * * * * * * * * * * * * * !
! * * * Basic korc array structures * * * !
! * * * * * * * * * * * * * * * * * * * * !

TYPE, PRIVATE :: KORC_MPI!< @brief KORC derived type to keep relevant MPI parameters.
    INTEGER :: nmpi !< Number of MPI processes.
	INTEGER :: rank !< Rank in WORLD COMMON communicator.
	INTEGER :: rank_topo !< Rank in mpi_topo communicator
	INTEGER :: mpi_topo !< MPI communicator for a certain topology.
END TYPE KORC_MPI


TYPE, PUBLIC :: CHARCS_PARAMS
!< @brief KORC derived type containing characteristic scales used in the normalization of KORC variables.
!! @details These characteristic scales are problem-dependent quantities. They are calculated in korc_units.f90 using the input
!! parameters of a given KORC simulation.

	REAL(rp) :: time
	!< @brief Characteristic non-relativistic time scale given by @f$1/\omega_{ce}@f$, where @f$\omega_{ce}=e B_0/m_e@f$ is the
	!! larger electron cyclotron frequency in the simulation.
	REAL(rp) :: time_r
	!< @brief Characteristic relativistic time scale given by @f$1/\omega_{ce}@f$, where @f$\omega_{ce}=e B_0/\gamma m_e@f$ is the
	!! larger relativistic electron cyclotron frequency in the simulation.
	REAL(rp) :: velocity !< Characteristic velocity. This is fixed to the speed of @f$c@f$.
	REAL(rp) :: length !< Characteristic length scale calculated as @f$c@f$ times the relativistic time scale.
	REAL(rp) :: mass !< Characteristic particle mass. This is equal to the electron mass @f$m_e@f$.
	REAL(rp) :: charge !< Characteristic particle charge. This is equal to the electron charge @f$q_e@f$.
	REAL(rp) :: density !< Characteristic particle density. This is equal to @f$1/l^3@f$, with @f$l@f$ the characteristic length.
	REAL(rp) :: Eo !< Characteristic electric field @f$E_0@f$. Usually @f$E_0@f$ at the magnetic axis.
	REAL(rp) :: Bo !< Characteristic magnetic field @f$B_0@f$ Usually @f$B_0@f$ at the magnetic axis.
	REAL(rp) :: energy !< Characteristic energy. This is equal to @f$m_e c^2@f$.
	REAL(rp) :: pressure !< Characteristic pressure. @todo This needs to be defined.
	REAL(rp) :: temperature !< Characteristic plasma temperature (Joules). This is equal to @f$m_e c^2@f$.
END TYPE CHARCS_PARAMS


TYPE, PUBLIC :: KORC_PARAMS
!< @brief Core KORC parameters.
!! @details This KORC derived type contains the variables that control KORC's core behavior.

	CHARACTER(MAX_STRING_LENGTH) 	:: path_to_inputs !< Absolute path to KORC's input file.
	CHARACTER(MAX_STRING_LENGTH) 	:: path_to_outputs !< Absolute path to the outputs' folder.
	INTEGER 						:: num_omp_threads !< Number of open MP threads per MPI process used in the simulation.
	LOGICAL 						:: restart !< Flag to indicate if the simulations restarts (restart=T) or not (restart=F).
	REAL(rp) 						:: simulation_time !< Total simulation time in seconds.
	REAL(rp) 						:: snapshot_frequency !< Time between snapshots in time of the simulation.
	REAL(rp) 						:: dt !< Time step in the simulation as a fraction of the relativistic electron gyro-period @f$\tau_e = 2\pi\gamma m_e/eB_0@f$.
	REAL(rp) 						:: time = 0.0_rp !< Current physical time in the simulation.
	INTEGER(ip) 					:: ito = 0_ip !< Initial time iteration in the simulation, this is different from zero in case is a restarting simulation.
	INTEGER(ip) 					:: it = 0_ip !< Current time iteration in the simulation, this is different from zero in case is a restarting simulation.
	INTEGER(ip) 					:: t_steps !< Number of time steps needed for evolving the electrons up to "simulation_time".
	INTEGER(ip) 					:: output_cadence !< Time iteration offset used to decide when the outputs are generated.
	INTEGER(ip) 					:: restart_output_cadence !< Time iteration offset used to decide when the restart files are generated.
	INTEGER(ip) 					:: num_snapshots !< Number of snapshots in time for generating the output files.
	INTEGER 						:: num_species !< Number of different populations of simulated relativistic electrons in KORC.
	REAL(rp) 						:: minimum_particle_energy !< Minimum allowed energy of simulated electrons. @todo To improve the criterium to decide when an electron will not be followed anymore in the simulation.
	REAL(rp) 						:: minimum_particle_g !< Minimum allowed relativistic factor @f$\gamma@f$ of simulated electrons.
	LOGICAL 						:: radiation !< Flag to indicate if synchrotron radiation losses are included (radiation=T) or not (radiation=F).
	LOGICAL 						:: collisions !< Flag to indicate if collisionsare included (collisions=T) or not (collisions=F).
	CHARACTER(MAX_STRING_LENGTH) 	:: collisions_model !< String with the name of the collisions model to be used in the simulation.
	CHARACTER(MAX_STRING_LENGTH) 	:: plasma_model !< String with the name of the model for the fields and plasma profiles.
	CHARACTER(MAX_STRING_LENGTH) 	:: magnetic_field_filename !< String with the name of the model for the fields and plasma profiles.
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: outputs_list !< List of electron variables to include in the outputs.
	INTEGER 						:: HDF5_error_handling !< Flag to indicate whether we allow HDF5 to show warnings during runtime (HDF5_error_handling=1) or not (HDF5_error_handling=0)

	TYPE(KORC_MPI) 					:: mpi_params !< An instance of the KORC_MPI derived type.
	TYPE(CHARCS_PARAMS) 			:: cpp !< An instance of the CHARCS_PARAMS derived type.

    INTEGER                         :: time_slice !< M3D-C1 time slice to use.
END TYPE KORC_PARAMS


TYPE, PUBLIC :: PARTICLES
!< @brief Derived type containing all the electrons' variables in the simulation.

    REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: X !< Cartesian coordinates of the electrons' position. dim(X) = (3,SPECIES::ppp).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: V !< Cartesian components of the electrons' velocity. dim(V) = dim(X).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: Rgc !< Cartesian coordinates of the electrons' guiding-center position.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: Y !< Coordinates of the electrons' position in cylindrical or toroidal coordinates.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: E !< Cartesian components of the electric field felt by each electron. dim(E) = dim(X).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: B !< Cartesian components of the magnetic field felt by each electron. dim(B) = dim(X).
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: ne !< Electron density seen by each electron. dim(ne) = (1,SPECIES::ppp).
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: Te !< Electron temperature seen by each electron. dim(ne) = (1,SPECIES::ppp).
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: Zeff !< Zeff seen by each electron. dim(ne) = (1,SPECIES::ppp).
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: g !< Instantaneous relativistic @f$\gamma = 1/\sqrt{1 - v^2/c^2}@f$ factor of each electron in the simulation.
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: eta !< Instantaneous pitch angle of each electron in the simulation.
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: mu !< Magnetic moment of each electron in the simulation.
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: Prad !< Instantaneous radiated power by each electron in the simulation.
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: Pin !< Instantaneous input power of each electron due to the electric field acceleration.
	INTEGER(is), DIMENSION(:), ALLOCATABLE :: flag !< Flag for each particle to decide whether it is being followed (flag=T) or not (flag=F).
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: AUX !< An auxiliary scalar variable for each electron.
	REAL(rp), DIMENSION(:), ALLOCATABLE    :: wt !< Weight of each electron. This is used when sampling weighted PDFs and in the synthetic camera diagnostic.
#ifdef M3D_C1
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: hint !< Hint for M3D_C1 interpolation.
    LOGICAL                                :: cart
#endif
END TYPE PARTICLES


TYPE, PUBLIC :: SPECIES
!< @brief Derived type containing the initial parameters of each electron ensemble followed in a KORC simulation.

	TYPE(PARTICLES) 				:: vars !< An instance of the PARTICLES derived type.
	LOGICAL 						:: runaway !< Flag to decide whether a given electron is a runaway (runaway=T) or not (runaway=F).
	CHARACTER(MAX_STRING_LENGTH) 	:: spatial_distribution !< String describing the type of initial spatial distribution for each electron species.
	CHARACTER(MAX_STRING_LENGTH) 	:: energy_distribution !< String describing the type of initial energy distribution for each electron species.
	CHARACTER(MAX_STRING_LENGTH) 	:: pitch_distribution !< String describing the type of initial pitch-angle distribution for each electron species.
	REAL(rp) 						:: Eo !< Initial energy of each electron species in case of using an initial mono-energetic distribution.
	REAL(rp) 						:: go !< Corresponding relativisitc factor of each electron species in case of using an initial mono-energetic distribution.
	REAL(rp) 						:: etao !< Initial pitch-angle of each electron species in case of using an initial mono-pitch-angle distribution.
	REAL(rp), DIMENSION(2) 			:: Eo_lims !< Minimum and maximum energy limits of a given initial non-mono-energetic distribution.
	REAL(rp), DIMENSION(2) 			:: etao_lims !< Minimum and maximum pitch-angle limits of a given initial non-mono-pitch-angle distribution.
	REAL(rp) 						:: wc !< The mean electron cyclotron frequency of each electron species.
	REAL(rp) 						:: wc_r !< The mean relativistic electron cyclotron frequency of each electron species.
	REAL(rp) 						:: q !< Charge of each species. @note This was left explicit to allow KORC to follow electrons and ions in the future.
	REAL(rp) 						:: m !< Mass of each species. @note This was left explicit to allow KORC to follow electrons and ions in the future.
	INTEGER 						:: ppp !< Number of computational particles used to simulate each electron species.


	REAL(rp) 						:: Ro !< Radial position of the center of the electrons' initial spatial distribution.
	REAL(rp) 						:: PHIo !< Azimuthal position of the electrons' initial spatial distribution, in case of using a disk at a certain poloidal section.
	REAL(rp) 						:: Zo !< Height of the center of the electrons' initial spatial distribution.
	REAL(rp) 						:: r_inner !< Minimum minor radius of the electrons' initial spatial distribution.
	REAL(rp) 						:: r_outter !< Maximum minor radius of the electrons' initial spatial distribution.
	REAL(rp) 						:: falloff_rate !< Exponential falloff or standard deviation of a non-uniform radial distribution of electrons.
	REAL(rp) 						:: shear_factor !< Shear factor used to generate an initial spatial distribution with an elliptic poloidal cross section. @note See <em>Carbajal and del-Castillo-Negrete, Nuclear Fusion, submitted (2018)</em>.
END TYPE SPECIES


TYPE, PRIVATE :: A_FIELD
!< @brief Derived type having all the parameters of the analytical magnetic field included in KORC.
!! @details The analytical magnetic field is given by:
!!
!! @f$\mathbf{B}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}} \left[ B_0 \hat{e}_\zeta  + B_\vartheta(r) \hat{e}_\vartheta \right]@f$,
!!
!! where @f$\eta = r/R_0@f$ is the aspect ratio, the constant @f$B_0@f$ denotes the magnitude of the toroidal magnetic field, and @f$B_\vartheta(r) = \eta B_0/q(r)@f$ is the poloidal magnetic field with safety factor @f$q(r) = q_0\left( 1 + \frac{r^2}{\lambda^2} \right)@f$. The constant @f$q_0@f$ is the safety factor at the magnetic axis and the constant @f$\lambda@f$ is obtained from the values of @f$q_0@f$ and @f$q(r)@f$ at the plasma edge @f$r=r_{edge}@f$.

	REAL(rp) 						:: Bo !< Magnitude of the toroidal magnetic field @f$B_0@f$.
	REAL(rp) 						:: a !< Plasma edge @f$r_{edge}@f$ as measured from the magnetic axis.
	REAL(rp) 						:: Ro !< Radial position of the magnetic axis @f$R_0@f$
	REAL(rp) 						:: qa !< Safety factor at the plasma edge.
	REAL(rp) 						:: qo !< Safety factor at the magnetic axis @f$q_0@f$.
	REAL(rp) 						:: lambda !< @f$\lambda@f$ parameter of @f$q(r)@f$.
	REAL(rp) 						:: Bpo !< @deprecated Parameter not used anymore. @todo Delete parameter.
	REAL(rp) 						:: Bp_sign !< Sign of @f$B_\vartheta(r)@f$. This depends on current_direction, Bp_sign=1 for current_direction='PARALLEL', and Bp_sign=-1 for current_direction='ANTI-PARALLEL'.
	CHARACTER(MAX_STRING_LENGTH) 	:: current_direction !< Direction of plasma current: PARALLEL or ANTI-PARALLEL to the toroidal magnetic field.
END TYPE A_FIELD

TYPE, PRIVATE :: MESH
!< Derived type with the cylindrical coordinates of the grid nodes at which the pre-computed plasma profiles and fields are known.
	REAL(rp), DIMENSION(:), ALLOCATABLE :: R !< Radial grid.
	REAL(rp), DIMENSION(:), ALLOCATABLE :: PHI !< Azimuthal grid.
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Z !< Z grid.
END TYPE MESH


TYPE, PUBLIC :: FIELDS
!< @brief Derived type with all the variables and data of analytical and pre-computed electric and magnetic fields.
#ifdef M3D_C1
    INTEGER                                 :: M3D_C1_B !< An M3D-C1 magnetic field.
    INTEGER                                 :: M3D_C1_E !< An M3D-C1 Electric field.
#endif

	TYPE(A_FIELD) 	 						:: AB !< An instance of the KORC derived data type A_FIELD.
	TYPE(V_FIELD_3D) 						:: E_3D !< KORC 3-D vector field of the pre-computed electric field.
	TYPE(V_FIELD_3D) 						:: B_3D !< KORC 3-D vector field of the pre-computed magnetic field.
	TYPE(V_FIELD_2D) 						:: E_2D !< KORC 2-D vector field of the pre-computed electric field.
	TYPE(V_FIELD_2D) 						:: B_2D !< KORC 3-D vector field of the pre-computed magnetic field.
	TYPE(MESH) 		 						:: X !< An instance of the KORC derived type MESH.
	INTEGER, DIMENSION(3) 					:: dims !< Dimensions of the KORC vector field. dims=(number of grid nodes along @f$R@f$, number of grid nodes along @f$\phi@f$, number of grid nodes along @f$Z@f$).
	REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: PSIp !< 2-D array for storing the data of the poloidal magnetic flux.
	REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: FLAG2D !< 2-D array defining the simulation domain where pre-computed data exist.
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: FLAG3D !< 3-D array defining the simulation domain where pre-computed data exist.

	REAL(rp) 								:: Eo !< Characteristic electric field.
	REAL(rp) 								:: Bo !< Characteristic magnetic field.
	REAL(rp) 								:: Ro !< Radial position of the magnetic axis.
	REAL(rp) 								:: Zo !< @f$Z@f$ position of the magnetic axis.

	LOGICAL 								:: Bfield !< Flag to indicate whether a pre-computed magnetic field will be used (Bfield=T) or not (Bfield=F).
	LOGICAL 								:: Bflux !< Flag to indicate whether a pre-computed poloidal magnetic flux will be used (Bflux=T) or not (Bflux=F).
	LOGICAL 								:: Efield !< Flag to indicate whether a pre-computed electric field will be used (Efield=T) or not (Efield=F).

	LOGICAL 								:: Bfield_in_file !< Flag to indicate if a pre-computed magnetic field is in the input file.
	LOGICAL 								:: Bflux_in_file !< Flag to indicate if a pre-computed poloidal magnetic flux is in the input file.
	LOGICAL 								:: Efield_in_file !< Flag to indicate if a pre-computed electric field is in the input file.

	LOGICAL 								:: axisymmetric_fields !< Flag to indicate if the pre-computed fields are axisymmetric.
END TYPE FIELDS


TYPE, PUBLIC :: PROFILES
!< @brief KORC derived data type having information about the plasma profiles. See korc_profiles.f90 for more information.
!< @details KORC can run using either analytical and pre-computed plasma profiles. Pre-computed plasma profiles, as in the case of pre-computed electric or magnetic fields, are interpolated to electrons' position in korc_profiles.f90.
!!
!! There are two types of analytical plasma profiles that can be used in KORC: 3rd degree polynomial radial plasma profiles,
!!
!! @f$f(r) = a_3r^3 + a_2r^2 +a_1r + a_0@f$,
!!
!! and radial plasma profiles with a @f$\tanh(r)@f$ dependency:
!!
!! @f$f(r) = f_0\left[1 - \tanh^n\left(\frac{2r}{a}\right)\right]@f$,
!!
!! where @f$f_0@f$ is a given plasma parameter at the magnetic axis, and @f$a@f$ is the plasma radius as measured from the magnetic axis to the last closed flux surface. Notice that the larger @f$n@f$ is, the more uniform the radial profiles are.

	TYPE(MESH) 								:: X !< An instance of the KORC derived data type MESH.
	REAL(rp) 								:: a !< Plasma radius as measured from the magnetic axis.

	INTEGER, DIMENSION(3) 					:: dims !< !< Dimensions of the arrays containing the pre-computed profiles data. dims=(number of grid nodes along @f$R@f$, number of grid nodes along @f$\phi@f$, number of grid nodes along @f$Z@f$).
	REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: FLAG2D !< 2-D array defining the simulation domain where pre-computed data exist.
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: FLAG3D !< 3-D array defining the simulation domain where pre-computed data exist.

	REAL(rp) 								:: n_ne !< @f$n@f$ used in @f$\tanh^n(r)@f$ of the electron density profile.
	REAL(rp) 								:: n_Te !< @f$n@f$ used in @f$\tanh^n(r)@f$ of the electron temperature profile.
	REAL(rp) 								:: n_Zeff !< @f$n@f$ used in @f$\tanh^n(r)@f$ of the @f$Z_{eff}@f$ profile.

	REAL(rp), DIMENSION(4) 					:: a_ne !< Coefficients of the polynomial electron density profile. See detailed description above, a_ne=(@f$a_{0}@f$,@f$a_{2}@f$,@f$a_{3}@f$,@f$a_{4}@f$).
	REAL(rp), DIMENSION(4) 					:: a_Te !< Coefficients of the polynomial electron temperature profile. See detailed description above, a_ne=(@f$a_{0}@f$,@f$a_{2}@f$,@f$a_{3}@f$,@f$a_{4}@f$).
	REAL(rp), DIMENSION(4) 					:: a_Zeff !< Coefficients of the @f$Z_{eff}@f$ profile. See detailed description above, a_ne=(@f$a_{0}@f$,@f$a_{2}@f$,@f$a_{3}@f$,@f$a_{4}@f$).

	! Zeff
	CHARACTER(MAX_STRING_LENGTH) 			:: Zeff_profile !< String containing the type of @f$Z_{eff}@f$ profile to be used in the simulation.
	REAL(rp) 								:: Zeffo !< @f$Z_{eff}@f$ at the magnetic axis.
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Zeff_3D !< 3-D array for keeping the pre-computed data of the @f$Z_{eff}@f$ profile.
	REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Zeff_2D !< 2-D array for keeping the pre-computed data of the @f$Z_{eff}@f$ profile.


	! Density
	CHARACTER(MAX_STRING_LENGTH) 			:: ne_profile !< String containing the type of electron density profile to be used in the simulation.
	REAL(rp) 								:: neo !< Electron density at the magnetic axis
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: ne_3D !< 3-D array for keeping the pre-computed data of the electron density profile.
	REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: ne_2D !< 2-D array for keeping the pre-computed data of the electron density profile.

	!Temperature
	CHARACTER(MAX_STRING_LENGTH) 			:: Te_profile !< String containing the type of electron temperature profile to be used in the simulation.
	REAL(rp) 								:: Teo !< Electron temperature at the magnetic axis
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Te_3D !< 3-D array for keeping the pre-computed data of the electron density profile.
	REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Te_2D !< 2-D array for keeping the pre-computed data of the electron density profile.

	CHARACTER(MAX_STRING_LENGTH) 			:: filename !< Full path to the HDF5 file containing the pre-computed plasma profiles.
	LOGICAL 					 			:: axisymmetric !< Flag to indicate if the plasma profiles are axisymmetric.

#ifdef M3D_C1
    INTEGER                                 :: M3D_C1_ne
    INTEGER                                 :: M3D_C1_te
    INTEGER                                 :: M3D_C1_zeff
#endif
END TYPE PROFILES

end module korc_types
