!> @brief Module containing functions and subroutines for performing interpolations using the PSPLINE library.
!! @details For a detailed documentation of the PSPLINE library we refer the user to "https://w3.pppl.gov/ntcc/PSPLINE/".
module korc_interp
    use korc_types
    use korc_coords
    use korc_rnd_numbers
	use korc_hpc

    use EZspline_obj	! psplines module
    use EZspline		! psplines module

    IMPLICIT NONE


#ifdef DOUBLE_PRECISION


!> @brief Derived type containing 3-D PSPLINE interpolants for cylindrical components of vector fields
!! @f$\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi + F_Z\hat{e}_Z@f$. Real precision of 8 bytes.
TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
	TYPE(EZspline3_r8)    :: R		!< Interpolant of @f$F_R(R,\phi,Z)@f$.
	TYPE(EZspline3_r8)    :: PHI	!< Interpolant of @f$F_\phi(R,\phi,Z)@f$.
	TYPE(EZspline3_r8)    :: Z		!< Interpolant of @f$F_Z(R,\phi,Z)@f$.

	INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NPHI !< Size of mesh containing the field data along the @f$\phi@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /) !< Periodic boundary condition for the interpolants at both ends of the @f$\phi@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE


!> @brief Derived type containing 2-D PSPLINE interpolants for cylindrical components of vector fields @f$\mathbf{F}(R,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z@f$.
!! Real precision of 8 bytes.
TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
	TYPE(EZspline2_r8)    :: A		!< Interpolant of a scalar field @f$A(R,Z)@f$.
	TYPE(EZspline2_r8)    :: R		!< Interpolant of @f$F_R(R,Z)@f$.
	TYPE(EZspline2_r8)    :: PHI	!< Interpolant of @f$F_\phi(R,Z)@f$.
	TYPE(EZspline2_r8)    :: Z		!< Interpolant of @f$F_Z(R,Z)@f$.

	INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE


!> @brief Derived type containing 3-D PSPLINE interpolants for cylindrical components of the density @f$n_e(R,\phi,Z)@f$,
!! temperature @f$T_e(R,\phi,Z)@f$, and effective charge number @f$Z_{eff}(R,\phi,Z)@f$ profiles.
!! Real precision of 8 bytes.
TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
	TYPE(EZspline3_r8)    :: ne    !< Interpolant of background electron density @f$n_e(R,\phi,Z)@f$.
	TYPE(EZspline3_r8)    :: Te    !< Interpolant of background electron temperature @f$T_e(R,\phi,Z)@f$.
	TYPE(EZspline3_r8)    :: Zeff  !< Interpolant of effective charge number @f$Z_{eff}(R,\phi,Z)@f$.

    INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NPHI !< Size of mesh containing the field data along the @f$\phi@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /) !< Periodic boundary condition for the interpolants at both ends of the @f$\phi@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE


!> @brief Derived type containing 2-D PSPLINE interpolants for cylindrical components of the density @f$n_e(R,Z)@f$,
!! temperature @f$T_e(R,Z)@f$, and effective charge number @f$Z_{eff}(R,Z)@f$ profiles.
!! Real precision of 8 bytes.
TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
	TYPE(EZspline2_r8)    :: ne	   !< Interpolant of background electron density @f$n_e(R,Z)@f$.
	TYPE(EZspline2_r8)    :: Te	   !< Interpolant of background electron temperature @f$T_e(R,Z)@f$.
	TYPE(EZspline2_r8)    :: Zeff  !< Interpolant of effective charge number @f$Z_{eff}(R,Z)@f$.

    INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE

#elif SINGLE_PRECISION


!> @brief Derived type containing 3-D PSPLINE interpolants for cylindrical components of vector fields @f$\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi + F_Z\hat{e}_Z@f$.
!! Real precision of 4 bytes.
TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
	TYPE(EZspline3_r4)    :: R     !< Interpolant of @f$F_R(R,\phi,Z)@f$.
	TYPE(EZspline3_r4)    :: PHI   !< Interpolant of @f$F_\phi(R,\phi,Z)@f$.
	TYPE(EZspline3_r4)    :: Z     !< Interpolant of @f$F_Z(R,\phi,Z)@f$.

    INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NPHI !< Size of mesh containing the field data along the @f$\phi@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /) !< Periodic boundary condition for the interpolants at both ends of the @f$\phi@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE


!> @brief Derived type containing 2-D PSPLINE interpolants for cylindrical components of vector fields @f$\mathbf{F}(R,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z@f$.
!! Real precision of 4 bytes.
TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
	TYPE(EZspline2_r4)    :: A		!< Interpolant of a scalar field @f$A(R,Z)@f$.
	TYPE(EZspline2_r4)    :: R		!< Interpolant of @f$F_R(R,Z)@f$.
	TYPE(EZspline2_r4)    :: PHI	!< Interpolant of @f$F_\phi(R,Z)@f$.
	TYPE(EZspline2_r4)    :: Z		!< Interpolant of @f$F_Z(R,Z)@f$.

    INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
    INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
    INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE


!> @brief Derived type containing 3-D PSPLINE interpolants for cylindrical components of the density @f$n_e(R,\phi,Z)@f$,
!! temperature @f$T_e(R,\phi,Z)@f$, and effective charge number @f$Z_{eff}(R,\phi,Z)@f$ profiles.
!! Real precision of 4 bytes.
TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
	TYPE(EZspline3_r4)    :: ne	   !< Interpolant of background electron density @f$n_e(R,\phi,Z)@f$.
	TYPE(EZspline3_r4)    :: Te	   !< Interpolant of background electron temperature @f$T_e(R,\phi,Z)@f$.
	TYPE(EZspline3_r4)    :: Zeff  !< Interpolant of effective charge number @f$Z_{eff}(R,\phi,Z)@f$.

    INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NPHI !< Size of mesh containing the field data along the @f$\phi@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /) !< Periodic boundary condition for the interpolants at both ends of the @f$\phi@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE


!> @brief Derived type containing 2-D PSPLINE interpolants for cylindrical components of the density @f$n_e(R,Z)@f$,
!! temperature @f$T_e(R,Z)@f$, and effective charge number @f$Z_{eff}(R,Z)@f$ profiles.
!! Real precision of 8 bytes.
TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
	TYPE(EZspline2_r4)    :: ne	   !< Interpolant of background electron density @f$n_e(R,Z)@f$.
	TYPE(EZspline2_r4)    :: Te	   !< Interpolant of background electron temperature @f$T_e(R,Z)@f$.
	TYPE(EZspline2_r4)    :: Zeff  !< Interpolant of effective charge number @f$Z_{eff}(R,Z)@f$.

    INTEGER               :: NR !< Size of mesh containing the field data along the @f$R@f$-axis.
    INTEGER               :: NZ !< Size of mesh containing the field data along the @f$Z@f$-axis.
	INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$R@f$ direction.
	INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) !< Not-a-knot boundary condition for the interpolants at both ends of the @f$Z@f$ direction.
END TYPE
#endif


!> @brief Derived type containing 2-D and 3-D arrays with the information of the spatial domain where the fields and profiles are known.
!! This info is used for detecting when a particle is lost, and therefore not followed anymore.
TYPE, PRIVATE :: KORC_INTERPOLANT_DOMAIN
	INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE      :: FLAG2D !< 2-D array with info of the spatial domain where the axisymmetric fields and plasma profiles are known.
	INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE    :: FLAG3D !< 3-D array with info of the spatial domain where the 3-D fields and plasma profiles are known.

	REAL(rp)                                          :: Ro !< Smaller radial position of the fields and profiles domain.
	REAL(rp)                                          :: Zo !< Smaller vertical position of the fields and profiles domain.

	REAL(rp)                                          :: DR    !< Separation between grid points along the radial direction.
	REAL(rp)                                          :: DPHI  !< Separation between grid points along the azimuthal direction.
	REAL(rp)                                          :: DZ    !< Separation between grid points along the vertical direction.
END TYPE


TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_2d !< An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating the magnetic field.
TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_3d !< An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating the magnetic field.
TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: efield_2d !< An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating the electric field.
TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: efield_3d !< An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating the electric field.
TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: fields_domain !< An instance of KORC_INTERPOLANT_DOMAIN used for interpolating fields.
TYPE(KORC_2D_PROFILES_INTERPOLANT), PRIVATE    :: profiles_2d !< An instance of KORC_2D_PROFILES_INTERPOLANT for interpolating plasma profiles.
TYPE(KORC_3D_PROFILES_INTERPOLANT), PRIVATE    :: profiles_3d !< An instance of KORC_3D_PROFILES_INTERPOLANT for interpolating plasma profiles.
TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: profiles_domain !< An instance of KORC_INTERPOLANT_DOMAIN used for interpolating plasma profiles.
INTEGER                                        :: ezerr !< Error status during PSPLINE interpolations.


PUBLIC :: interp_fields,&
		initialize_fields_interpolant,&
		initialize_profiles_interpolant,&
		finalize_interpolants
PRIVATE :: interp_3D_bfields,&
			interp_2D_bfields,&
			interp_3D_efields,&
			interp_2D_efields,&
			interp_2D_profiles,&
			interp_3D_profiles,&
			check_if_in_fields_domain,&
			check_if_in_profiles_domain

CONTAINS


!> @brief Subroutine that initializes fields interpolants.
!! @details This subroutine initializes either 2-D or 3-D PSPLINE interpolants using the data of fields in the KORC-dervied-type variable F.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation.
!! See korc_types.f90 and korc_fields.f90.
subroutine initialize_fields_interpolant(params,F)
	TYPE(KORC_PARAMS), INTENT(IN)  :: params
	TYPE(FIELDS), INTENT(IN)       :: F

	if (params%plasma_model .EQ. 'EXTERNAL') then

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')
		end if

		! * * * * * * * * MAGNETIC FIELD * * * * * * * * !
		if (F%Bflux) then
			bfield_2d%NR = F%dims(1)
			bfield_2d%NZ = F%dims(3)

			! Initializing poloidal flux interpolant
			call EZspline_init(bfield_2d%A,bfield_2d%NR,bfield_2d%NZ,bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			bfield_2d%A%x1 = F%X%R
			bfield_2d%A%x2 = F%X%Z

			call EZspline_setup(bfield_2d%A, F%PSIp, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))
			fields_domain%FLAG2D = F%FLAG2D

			fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
			fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
		else if (F%axisymmetric_fields) then
			bfield_2d%NR = F%dims(1)
			bfield_2d%NZ = F%dims(3)

			! Initializing R component
			call EZspline_init(bfield_2d%R,bfield_2d%NR,bfield_2d%NZ,bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			bfield_2d%R%x1 = F%X%R
			bfield_2d%R%x2 = F%X%Z

			call EZspline_setup(bfield_2d%R, F%B_2D%R, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			! Initializing PHI component
			call EZspline_init(bfield_2d%PHI,bfield_2d%NR,bfield_2d%NZ,bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			bfield_2d%PHI%x1 = F%X%R
			bfield_2d%PHI%x2 = F%X%Z

			call EZspline_setup(bfield_2d%PHI, F%B_2D%PHI, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			! Initializing Z component
			call EZspline_init(bfield_2d%Z,bfield_2d%NR,bfield_2d%NZ,bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			bfield_2d%Z%x1 = F%X%R
			bfield_2d%Z%x2 = F%X%Z

			call EZspline_setup(bfield_2d%Z, F%B_2D%Z, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))
			fields_domain%FLAG2D = F%FLAG2D

			fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
			fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
		else
			bfield_3d%NR = F%dims(1)
			bfield_3d%NPHI = F%dims(2)
			bfield_3d%NZ = F%dims(3)

			! Initializing R component of interpolant
			call EZspline_init(bfield_3d%R, bfield_3d%NR, bfield_3d%NPHI, bfield_3d%NZ,&
								bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			bfield_3d%R%x1 = F%X%R
			! bfield_3d%R%x2 = F%X%PHI
			bfield_3d%R%x3 = F%X%Z

			call EZspline_setup(bfield_3d%R, F%B_3D%R, ezerr)
			call EZspline_error(ezerr)

			! Initializing PHI component of interpolant
			call EZspline_init(bfield_3d%PHI, bfield_3d%NR, bfield_3d%NPHI, bfield_3d%NZ,&
								bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			bfield_3d%PHI%x1 = F%X%R
			! bfield_3d%PHI%x2 = F%X%PHI
			bfield_3d%PHI%x3 = F%X%Z

			call EZspline_setup(bfield_3d%PHI, F%B_3D%PHI, ezerr)
			call EZspline_error(ezerr)

			! Initializing Z component of interpolant
			call EZspline_init(bfield_3d%Z, bfield_3d%NR, bfield_3d%NPHI, bfield_3d%NZ,&
								bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			bfield_3d%Z%x1 = F%X%R
			! bfield_3d%Z%x2 = F%X%PHI
			bfield_3d%Z%x3 = F%X%Z

			call EZspline_setup(bfield_3d%Z, F%B_3D%Z, ezerr)
			call EZspline_error(ezerr)

			ALLOCATE(fields_domain%FLAG3D(bfield_3d%NR,bfield_3d%NPHI,bfield_3d%NZ))
			fields_domain%FLAG3D = F%FLAG3D

			fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
			fields_domain%DPHI = 2.0_rp*C_PI/bfield_3d%NPHI
			fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
		end if

		fields_domain%Ro = F%X%R(1)
		fields_domain%Zo = F%X%Z(1)

		! * * * * * * * * ELECTRIC FIELD * * * * * * * * !
		if (F%Efield.AND.F%Efield_in_file) then
			if (F%axisymmetric_fields) then
				efield_2d%NR = F%dims(1)
				efield_2d%NZ = F%dims(3)

				! Initializing R component
				call EZspline_init(efield_2d%R,efield_2d%NR,efield_2d%NZ,efield_2d%BCSR,efield_2d%BCSZ,ezerr)
			  	call EZspline_error(ezerr)

				efield_2d%R%x1 = F%X%R
				efield_2d%R%x2 = F%X%Z

				call EZspline_setup(efield_2d%R, F%E_2D%R, ezerr, .TRUE.)
				call EZspline_error(ezerr)

				! Initializing PHI component
				call EZspline_init(efield_2d%PHI,efield_2d%NR,efield_2d%NZ,efield_2d%BCSR,efield_2d%BCSZ,ezerr)
			  	call EZspline_error(ezerr)

				efield_2d%PHI%x1 = F%X%R
				efield_2d%PHI%x2 = F%X%Z

				call EZspline_setup(efield_2d%PHI, F%E_2D%PHI, ezerr, .TRUE.)
				call EZspline_error(ezerr)

				! Initializing Z component
				call EZspline_init(efield_2d%Z,efield_2d%NR,efield_2d%NZ,efield_2d%BCSR,efield_2d%BCSZ,ezerr)
			  	call EZspline_error(ezerr)

				efield_2d%Z%x1 = F%X%R
				efield_2d%Z%x2 = F%X%Z

				call EZspline_setup(efield_2d%Z, F%E_2D%Z, ezerr, .TRUE.)
				call EZspline_error(ezerr)
			else
				efield_3d%NR = F%dims(1)
				efield_3d%NPHI = F%dims(2)
				efield_3d%NZ = F%dims(3)

				! Initializing R component of interpolant
				call EZspline_init(efield_3d%R, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
									efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
			  	call EZspline_error(ezerr)

				efield_3d%R%x1 = F%X%R
				! efield_3d%R%x2 = F%X%PHI
				efield_3d%R%x3 = F%X%Z

				call EZspline_setup(efield_3d%R, F%E_3D%R, ezerr)
				call EZspline_error(ezerr)

				! Initializing PHI component of interpolant
				call EZspline_init(efield_3d%PHI, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
									efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
			  	call EZspline_error(ezerr)

				efield_3d%PHI%x1 = F%X%R
				! efield_3d%PHI%x2 = F%X%PHI
				efield_3d%PHI%x3 = F%X%Z

				call EZspline_setup(efield_3d%PHI, F%E_3D%PHI, ezerr)
				call EZspline_error(ezerr)

				! Initializing Z component of interpolant
				call EZspline_init(efield_3d%Z, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
									efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
			  	call EZspline_error(ezerr)

				efield_3d%Z%x1 = F%X%R
				! efield_3d%Z%x2 = F%X%PHI
				efield_3d%Z%x3 = F%X%Z

				call EZspline_setup(efield_3d%Z, F%E_3D%Z, ezerr)
				call EZspline_error(ezerr)
			end if
		end if

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * * * INTERPOLANT INITIALIZED * * * * * *")')
			write(6,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
		end if
	else if (params%plasma_model .EQ. 'ANALYTICAL') then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *")')
			write(6,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
		end if
	else if (params%plasma_model .EQ. 'UNIFORM') then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * *  * * * * * * * * * * *")')
			write(6,'("* * * * USING UNIFORM MAGNETIC FIELD * * * *")')
			write(6,'("* * * * * * * * * * *  * * * * * * * * * * *",/)')
		end if
	end if
end subroutine initialize_fields_interpolant


!> @brief Subrotuine that checks if particles in the simulation are within the spatial domain where interpolants and fields are known.
!! @details External fields and interpolants can have different spatial domains where they are defined. Therefore, it is necessary to
!! check if a given particle has left these spatial domains to stop following it, otherwise this will cause an error in the simulation.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,ou] flag Flag that determines whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param IR Variable used to localize the grid cell in the @f$(R,\phi,Z)@f$ or @f$(R,Z)@f$ grid containing the fields data that
!! corresponds to the radial position of the particles.
!! @param IPHI Variable used to localize the grid cell in the @f$(R,\phi,Z)@f$ or @f$(R,Z)@f$ grid containing the fields data that
!! corresponds to the azimuthal position of the particles.
!! @param IZ Variable used to localize the grid cell in the @f$(R,\phi,Z)@f$ or @f$(R,Z)@f$ grid containing the fields data that
!! corresponds to the vertical position of the particles.
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine check_if_in_fields_domain(Y,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: IR
    INTEGER                                                :: IPHI
    INTEGER                                                :: IZ
    INTEGER(ip)                                            :: pp
    INTEGER(ip)                                            :: ss

	ss = SIZE(Y,2)

	if (ALLOCATED(fields_domain%FLAG3D)) then
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) SHARED(Y,flag,fields_domain,bfield_3d)
		do pp=1_idef,ss
			IR = INT(FLOOR((Y(1,pp)  - fields_domain%Ro + 0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
			IPHI = INT(FLOOR((Y(2,pp)  + 0.5_rp*fields_domain%DPHI)/fields_domain%DPHI) + 1.0_rp,idef)
			IZ = INT(FLOOR((Y(3,pp)  + ABS(fields_domain%Zo) + 0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

			if ((fields_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR.((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ))) then
				flag(pp) = 0_is
			end if
		end do
!$OMP END PARALLEL DO
	else
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) SHARED(Y,flag,fields_domain,bfield_2d)
		do pp=1_idef,ss
			IR = INT(FLOOR((Y(1,pp)  - fields_domain%Ro + 0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
			IZ = INT(FLOOR((Y(3,pp)  + ABS(fields_domain%Zo) + 0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

			if ((fields_domain%FLAG2D(IR,IZ).NE.1_is).OR.((IR.GT.bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ))) then
				flag(pp) = 0_is
			end if
		end do
!$OMP END PARALLEL DO
	end if
end subroutine check_if_in_fields_domain


!> @brief Subroutine that initializes plasma profiles interpolants.
!! @details This subroutine initializes either 2-D or 3-D PSPLINE interpolants using the data of plasma profiles in the KORC-dervied-type variable P.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in] P An instance of KORC's derived type PROFILES containing all the information about the plasma profiles used in the simulation.
!! See korc_types.f90 and korc_profiles.f90.
subroutine initialize_profiles_interpolant(params,P)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PROFILES), INTENT(INOUT) :: P

	if (params%collisions) then
		if (params%plasma_model .EQ. 'EXTERNAL') then

			if (params%mpi_params%rank .EQ. 0) then
				write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * * *")')
				write(6,'("* * * * INITIALIZING PROFILES INTERPOLANT * * * *")')
			end if

			if (P%axisymmetric) then
				profiles_2d%NR = P%dims(1)
				profiles_2d%NZ = P%dims(3)

				! Initializing ne
	!			call EZspline_init(profiles_2d%ne,profiles_2d%NR,profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
				call EZlinear_Init(profiles_2d%ne,profiles_2d%NR,profiles_2d%NZ,ezerr)
			  	call EZspline_error(ezerr)

				profiles_2d%ne%x1 = P%X%R
				profiles_2d%ne%x2 = P%X%Z

				call EZspline_setup(profiles_2d%ne, P%ne_2D, ezerr, .TRUE.)
				call EZspline_error(ezerr)

				! Initializing Te
	!			call EZspline_init(profiles_2d%Te,profiles_2d%NR,profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
				call EZlinear_Init(profiles_2d%Te,profiles_2d%NR,profiles_2d%NZ,ezerr)
			  	call EZspline_error(ezerr)

				profiles_2d%Te%x1 = P%X%R
				profiles_2d%Te%x2 = P%X%Z

				call EZspline_setup(profiles_2d%Te, P%Te_2D, ezerr, .TRUE.)
				call EZspline_error(ezerr)

				! Initializing Zeff
	!			call EZspline_init(profiles_2d%Zeff,profiles_2d%NR,profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
				call EZlinear_Init(profiles_2d%Zeff,profiles_2d%NR,profiles_2d%NZ,ezerr)
			  	call EZspline_error(ezerr)

				profiles_2d%Zeff%x1 = P%X%R
				profiles_2d%Zeff%x2 = P%X%Z

				call EZspline_setup(profiles_2d%Zeff, P%Zeff_2D, ezerr, .TRUE.)
				call EZspline_error(ezerr)

				ALLOCATE(profiles_domain%FLAG2D(profiles_3d%NR,profiles_3d%NZ))
				profiles_domain%FLAG2D = P%FLAG2D

				profiles_domain%DR = ABS(P%X%R(2) - P%X%R(1))
				profiles_domain%DZ = ABS(P%X%Z(2) - P%X%Z(1))
			else
				profiles_3d%NR = P%dims(1)
				profiles_3d%NPHI = P%dims(2)
				profiles_3d%NZ = P%dims(3)

				! Initializing ne
				call EZspline_init(profiles_3d%ne, profiles_3d%NR, profiles_3d%NPHI, profiles_3d%NZ,&
									profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
			  	call EZspline_error(ezerr)

				profiles_3d%ne%x1 = P%X%R
				! profiles_3d%ne%x2 = P%X%PHI
				profiles_3d%ne%x3 = P%X%Z

				call EZspline_setup(profiles_3d%ne, P%ne_3D, ezerr)
				call EZspline_error(ezerr)

				! Initializing Te
				call EZspline_init(profiles_3d%Te, profiles_3d%NR, profiles_3d%NPHI, profiles_3d%NZ,&
									profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
			  	call EZspline_error(ezerr)

				profiles_3d%Te%x1 = P%X%R
				! profiles_3d%Te%x2 = P%X%PHI
				profiles_3d%Te%x3 = P%X%Z

				call EZspline_setup(profiles_3d%Te, P%Te_3D, ezerr)
				call EZspline_error(ezerr)

				! Initializing Zeff
				call EZspline_init(profiles_3d%Zeff, profiles_3d%NR, profiles_3d%NPHI, profiles_3d%NZ,&
									profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
			  	call EZspline_error(ezerr)

				profiles_3d%Zeff%x1 = P%X%R
				! profiles_3d%Zeff%x2 = P%X%PHI
				profiles_3d%Zeff%x3 = P%X%Z

				call EZspline_setup(profiles_3d%Zeff, P%Zeff_3D, ezerr)
				call EZspline_error(ezerr)

				ALLOCATE(profiles_domain%FLAG3D(profiles_3d%NR,profiles_3d%NPHI,profiles_3d%NZ))
				profiles_domain%FLAG3D = P%FLAG3D

				profiles_domain%DR = ABS(P%X%R(2) - P%X%R(1))
				profiles_domain%DPHI = 2.0_rp*C_PI/profiles_3d%NPHI
				profiles_domain%DZ = ABS(P%X%Z(2) - P%X%Z(1))
			end if

			profiles_domain%Ro = P%X%R(1)
			profiles_domain%Zo = P%X%Z(1)

			if (params%mpi_params%rank .EQ. 0) then
				write(6,'("* * * * * * INTERPOLANT   INITIALIZED * * * * * *")')
				write(6,'("* * * * * * * * * * * * * * * * * * * * * * * * *",/)')
			end if
		else if (params%plasma_model .EQ. 'ANALYTICAL') then
			if (params%mpi_params%rank .EQ. 0) then
				write(6,'(/,"* * * * * * * * * * * * * * * * * * * * *")')
				write(6,'("* * * * USING ANALYTICAL PROFILES * * * *")')
				write(6,'("* * * * * * * * * * * * * * * * * * * * *",/)')
			end if
		else if (params%plasma_model .EQ. 'UNIFORM') then
			if (params%mpi_params%rank .EQ. 0) then
				write(6,'(/,"* * * * * * * * * * * * *  * * * * * * * * * * *")')
				write(6,'("* * * * UNIFORM PLASMA: NO PROFILES USED * * * *")')
				write(6,'("* * * * * * * * * * * * *  * * * * * * * * * * *",/)')
			end if
		end if
	end if
end subroutine initialize_profiles_interpolant


!> @brief Subrotuine that checks if particles in the simulation are within the spatial domain where interpolants and plasma profiles are known.
!! @details External plasma profiles and interpolants can have different spatial domains where they are defined.
!! Therefore, it is necessary to check if a given particle has left these spatial domains to stop following it, otherwise this will
!! cause an error in the simulation.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,ou] flag Flag that determines whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param IR Variable used to localize the grid cell in the @f$(R,\phi,Z)@f$ or @f$(R,Z)@f$ grid containing the fields data that
!! corresponds to the radial position of the particles.
!! @param IPHI Variable used to localize the grid cell in the @f$(R,\phi,Z)@f$ or @f$(R,Z)@f$ grid containing the fields data that
!! corresponds to the azimuthal position of the particles.
!! @param IZ Variable used to localize the grid cell in the @f$(R,\phi,Z)@f$ or @f$(R,Z)@f$ grid containing the fields data that
!! corresponds to the vertical position of the particles.
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine check_if_in_profiles_domain(Y,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: IR
    INTEGER                                                :: IPHI
    INTEGER                                                :: IZ
	INTEGER(ip)                                            :: pp
    INTEGER(ip)                                            :: ss

	ss = SIZE(Y,2)

	if (ALLOCATED(profiles_domain%FLAG3D)) then
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) SHARED(Y,flag,profiles_domain,profiles_3d)
		do pp=1_idef,ss
			IR = INT(FLOOR((Y(1,pp)  - profiles_domain%Ro + 0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
			IPHI = INT(FLOOR((Y(2,pp)  + 0.5_rp*profiles_domain%DPHI)/profiles_domain%DPHI) + 1.0_rp,idef)
			IZ = INT(FLOOR((Y(3,pp)  + ABS(profiles_domain%Zo) + 0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

			if ((profiles_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR.((IR.GT.profiles_3d%NR).OR.(IZ.GT.profiles_3d%NZ))) then
				flag(pp) = 0_is
			end if
		end do
!$OMP END PARALLEL DO
	else
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) SHARED(Y,flag,profiles_domain,profiles_2d)
		do pp=1_idef,ss
			IR = INT(FLOOR((Y(1,pp)  - profiles_domain%Ro + 0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
			IZ = INT(FLOOR((Y(3,pp)  + ABS(profiles_domain%Zo) + 0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

			if ((profiles_domain%FLAG2D(IR,IZ).NE.1_is).OR.((IR.GT.profiles_2d%NR).OR.(IZ.GT.profiles_2d%NZ))) then
				flag(pp) = 0_is
			end if
		end do
!$OMP END PARALLEL DO
	end if
end subroutine check_if_in_profiles_domain


!> @brief Subroutine for interpolating the pre-computed, axisymmetric magnetic field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=@f$B_x@f$, B(2,:)=@f$B_y@f$, and B(3,:)=@f$B_z@f$.
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=@f$B_R@f$, F(2,:)=@f$B_\phi@f$, and F(3,:)=@f$B_Z@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_2D_bfields(Y,B,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
	REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,B,flag,bfield_2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(bfield_2d%R, Y(1,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(bfield_2d%PHI, Y(1,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(bfield_2d%Z, Y(1,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			B(1,pp) = F(1,pp)*COS(Y(2,pp)) - F(2,pp)*SIN(Y(2,pp))
			B(2,pp) = F(1,pp)*SIN(Y(2,pp)) + F(2,pp)*COS(Y(2,pp))
			B(3,pp) = F(3,pp)
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(F)
end subroutine interp_2D_bfields


!> @brief Subroutine for interpolating the pre-computed, 3-D magnetic field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=@f$B_x@f$, B(2,:)=@f$B_y@f$, and B(3,:)=@f$B_z@f$.
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=@f$B_R@f$, F(2,:)=@f$B_\phi@f$, and F(3,:)=@f$B_Z@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_bfields(Y,B,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
	REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,B,flag,bfield_3d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(bfield_3d%R, Y(1,pp), Y(2,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(bfield_3d%PHI, Y(1,pp), Y(2,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(bfield_3d%Z, Y(1,pp), Y(2,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			B(1,pp) = F(1,pp)*COS(Y(2,pp)) - F(2,pp)*SIN(Y(2,pp))
			B(2,pp) = F(1,pp)*SIN(Y(2,pp)) + F(2,pp)*COS(Y(2,pp))
			B(3,pp) = F(3,pp)
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(F)
end subroutine interp_3D_bfields


!> @brief Subroutine that calculates the axisymmetric magnetic field to the particles' position using the poloidal magnetic flux.
!! @details When the poloidal magnetic flux @f$\Psi(R,Z)@f$ is used in a KORC simulation, the magnetic field components are calculated as it follows:
!!
!!
!! @f$B_R = \frac{1}{R}\frac{\partial \Psi}{\partial Z}@f$,\n
!! @f$B_\phi = \frac{RoBo}{R}@f$,\n
!! @f$B_Z = -\frac{1}{R}\frac{\partial \Psi}{\partial R}@f$,
!!
!!
!! where @f$Ro@f$ and @f$Bo@f$ are the radial position of the magnetic axis and the magnetic field as measured at the magnetic axis, respectively.
!! First, the derivatives of the poloidal magnetic flux are calculated at the particles' position using the PSPLINE interpolant of
!! the poloidal magnetic flux. Then, we calculate the cylindrical components of the magnetic field, and finally we calculate its Cartesian
!! components that will be used in the particle pusher.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation.
!! See korc_types.f90 and korc_fields.f90.
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=@f$B_x@f$, B(2,:)=@f$B_y@f$, and B(3,:)=@f$B_x@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param A Variable containing the partial derivatives of the poloidal magnetic flux @f$\Psi(R,Z)@f$ and the cylindrical components
!! of the magnetic field (its value changes through the subroutine).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine calculate_magnetic_field(Y,F,B,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	TYPE(FIELDS), INTENT(IN)                               :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: A
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

	ALLOCATE(A(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,A,B,flag,bfield_2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			! FR = (dA/dZ)/R
			call EZspline_derivative(bfield_2d%A, 0, 1, Y(1,pp), Y(3,pp), A(1,pp), ezerr)
!			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			else
				A(1,pp) = A(1,pp)/Y(1,pp)

				! FPHI = Fo*Ro/R
				A(2,pp) = F%Bo*F%Ro/Y(1,pp)

				! FR = -(dA/dR)/R
				call EZspline_derivative(bfield_2d%A, 1, 0, Y(1,pp), Y(3,pp), A(3,pp), ezerr)
				call EZspline_error(ezerr)
				A(3,pp) = -A(3,pp)/Y(1,pp)

				B(1,pp) = A(1,pp)*COS(Y(2,pp)) - A(2,pp)*SIN(Y(2,pp))
				B(2,pp) = A(1,pp)*SIN(Y(2,pp)) + A(2,pp)*COS(Y(2,pp))
				B(3,pp) = A(3,pp)
			end if
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(A)
end subroutine calculate_magnetic_field


!> @brief Subroutine for interpolating the pre-computed, axisymmetric electric field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,out] E Cartesian components of interpolated electric field components. E(1,:)=@f$E_x@f$, E(2,:)=@f$E_y@f$, and E(3,:)=@f$E_z@f$.
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=@f$E_R@f$, F(2,:)=@f$E_\phi@f$, and F(3,:)=@f$E_Z@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_2D_efields(Y,E,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
	REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,E,flag,efield_2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(efield_2d%R, Y(1,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(efield_2d%PHI, Y(1,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(efield_2d%Z, Y(1,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			E(1,pp) = F(1,pp)*COS(Y(2,pp)) - F(2,pp)*SIN(Y(2,pp))
			E(2,pp) = F(1,pp)*SIN(Y(2,pp)) + F(2,pp)*COS(Y(2,pp))
			E(3,pp) = F(3,pp)
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(F)
end subroutine interp_2D_efields


!> @brief Subroutine for interpolating the pre-computed 3-D electric field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,out] E Cartesian components of interpolated electric field components. E(1,:)=@f$E_x@f$, E(2,:)=@f$E_y@f$, and E(3,:)=@f$E_z@f$.
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=@f$E_R@f$, F(2,:)=@f$E_\phi@f$, and F(3,:)=@f$E_Z@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_efields(Y,E,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
	REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,E,flag,efield_3d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(efield_3d%R, Y(1,pp), Y(2,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(efield_3d%PHI, Y(1,pp), Y(2,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(efield_3d%Z, Y(1,pp), Y(2,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			E(1,pp) = F(1,pp)*COS(Y(2,pp)) - F(2,pp)*SIN(Y(2,pp))
			E(2,pp) = F(1,pp)*SIN(Y(2,pp)) + F(2,pp)*COS(Y(2,pp))
			E(3,pp) = F(3,pp)
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(F)
end subroutine interp_3D_efields


!> @brief Subroutine that works as an interface for calling the appropriate subroutines for interpolating or calculating the electric and magnetic fields.
!!
!! @param[in,out] prtcls An instance of PARTICLES containing the variables of a given species.
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation.
!! See korc_types.f90 and korc_fields.f90.
subroutine interp_fields(prtcls,F)
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN)       :: F

	call cart_to_cyl(prtcls%X, prtcls%Y)

	call check_if_in_fields_domain(prtcls%Y, prtcls%flag)

	if (ALLOCATED(F%PSIp)) then
		call calculate_magnetic_field(prtcls%Y,F,prtcls%B,prtcls%flag)
	end if

	if (ALLOCATED(F%B_2D%R)) then
		call interp_2D_bfields(prtcls%Y,prtcls%B,prtcls%flag)
	end if

	if (ALLOCATED(F%B_3D%R)) then
		call interp_3D_bfields(prtcls%Y,prtcls%B,prtcls%flag)
	end if

	if (ALLOCATED(F%E_2D%R)) then
		call interp_2D_efields(prtcls%Y,prtcls%E,prtcls%flag)
	end if

	if (ALLOCATED(F%E_3D%R)) then
		call interp_3D_efields(prtcls%Y,prtcls%E,prtcls%flag)
	end if
end subroutine interp_fields


!> @brief Subroutine for interpolating the pre-computed, axisymmetric plasma profiles to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,out] ne Interpolated background electron density @f$n_e(R,Z)@f$.
!! @param[in,out] Te Interpolated background electron temperature @f$T_e(R,Z)@f$.
!! @param[in,out] Zeff Interpolated effective charge number @f$Z_{eff}(R,Z)@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_2D_profiles(Y,ne,Te,Zeff,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: ne
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Te
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Zeff
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(Y,ne,Te,Zeff,flag,profiles_2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(profiles_2d%ne, Y(1,pp), Y(3,pp), ne(pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(profiles_2d%Te, Y(1,pp), Y(3,pp), Te(pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(profiles_2d%Zeff, Y(1,pp), Y(3,pp), Zeff(pp), ezerr)
			call EZspline_error(ezerr)
		end if
	end do
!$OMP END PARALLEL DO
end subroutine interp_2D_profiles


!> @brief Subroutine for interpolating the pre-computed, 3-D plasma profiles to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, and Y(3,:) = @f$Z@f$.
!! @param[in,out] ne Interpolated background electron density @f$n_e(R,\phi,Z)@f$.
!! @param[in,out] Te Interpolated background electron temperature @f$T_e(R,\phi,Z)@f$.
!! @param[in,out] Zeff Interpolated effective charge number @f$Z_{eff}(R,\phi,Z)@f$.
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_profiles(Y,ne,Te,Zeff,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: ne
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Te
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Zeff
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
	INTEGER                                                :: pp
    INTEGER                                                :: ss

	ss = size(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(Y,ne,Te,Zeff,flag,profiles_2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(profiles_3d%ne, Y(1,pp), Y(2,pp), Y(3,pp), ne(pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(profiles_3d%Te, Y(1,pp), Y(2,pp), Y(3,pp), Te(pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(profiles_3d%Zeff, Y(1,pp), Y(2,pp), Y(3,pp), Zeff(pp), ezerr)
			call EZspline_error(ezerr)
		end if
	end do
!$OMP END PARALLEL DO
end subroutine interp_3D_profiles


!> @brief Subroutine that calls the appropriate subroutines for interpolating the 2-D or 3-D plasma profiles to the particles' position.
!!
!! @param[in,out] prtcls An instance of PARTICLES containing the variables of a given species.
!! @param[in] P An instance of KORC's derived type PROFILES containing all the information about the plasma profiles used in the simulation.
!!  See korc_types.f90 and korc_profiles.f90.
subroutine interp_profiles(prtcls,P)
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(PROFILES), INTENT(IN)     :: P

!	call cart_to_cyl(prtcls%X, prtcls%Y) ! This was done before when interpolating fields

	call check_if_in_profiles_domain(prtcls%Y, prtcls%flag)

	if (ALLOCATED(P%ne_2D)) then
		call interp_2D_profiles(prtcls%Y,prtcls%ne,prtcls%Te,prtcls%Zeff,prtcls%flag)
	else if (ALLOCATED(P%ne_3D)) then
		call interp_3D_profiles(prtcls%Y,prtcls%ne,prtcls%Te,prtcls%Zeff,prtcls%flag)
	else
		write(6,'("Error: NO PROFILES ALLOCATED")')
		call KORC_ABORT()
	end if
end subroutine interp_profiles


!> @brief Subroutine that frees memory allocated for PSPLINE interpolants.
!!
!! @param[in] params Core KORC simulation parameters.
subroutine finalize_interpolants(params)
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%plasma_model .EQ. 'EXTERNAL') then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * FINALIZING INTERPOLANT * * * *")')
		end if

		if (EZspline_allocated(bfield_3d%R)) call Ezspline_free(bfield_3d%R, ezerr)
		if (EZspline_allocated(bfield_3d%PHI)) call Ezspline_free(bfield_3d%PHI, ezerr)
		if (EZspline_allocated(bfield_3d%Z)) call Ezspline_free(bfield_3d%Z, ezerr)
		if (EZspline_allocated(bfield_2d%A)) call Ezspline_free(bfield_2d%A, ezerr)
		if (EZspline_allocated(bfield_2d%R)) call Ezspline_free(bfield_2d%R, ezerr)
		if (EZspline_allocated(bfield_2d%PHI)) call Ezspline_free(bfield_2d%PHI, ezerr)
		if (EZspline_allocated(bfield_2d%Z)) call Ezspline_free(bfield_2d%Z, ezerr)

		if (EZspline_allocated(profiles_3d%ne)) call Ezspline_free(profiles_3d%ne, ezerr)
		if (EZspline_allocated(profiles_3d%Te)) call Ezspline_free(profiles_3d%Te, ezerr)
		if (EZspline_allocated(profiles_3d%Zeff)) call Ezspline_free(profiles_3d%Zeff, ezerr)
		if (EZspline_allocated(profiles_2d%ne)) call Ezspline_free(profiles_2d%ne, ezerr)
		if (EZspline_allocated(profiles_2d%Te)) call Ezspline_free(profiles_2d%Te, ezerr)
		if (EZspline_allocated(profiles_2d%Zeff)) call Ezspline_free(profiles_2d%Zeff, ezerr)


		if (ALLOCATED(fields_domain%FLAG2D)) DEALLOCATE(fields_domain%FLAG2D)
		if (ALLOCATED(fields_domain%FLAG3D)) DEALLOCATE(fields_domain%FLAG3D)

		if (ALLOCATED(profiles_domain%FLAG2D)) DEALLOCATE(profiles_domain%FLAG2D)
		if (ALLOCATED(profiles_domain%FLAG3D)) DEALLOCATE(profiles_domain%FLAG3D)

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * INTERPOLANT  FINALIZED * * * *")')
		end if
	end if
end subroutine finalize_interpolants
end module korc_interp
