module korc_interp
    use korc_types
    use korc_coords
    use korc_rnd_numbers
	use korc_hpc
    use EZspline_obj	! psplines module
    use EZspline		! psplines module

    implicit none

#ifdef DOUBLE_PRECISION
	TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
		TYPE(EZspline3_r8) :: R		! 3D EZspline object
		TYPE(EZspline3_r8) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r8) :: Z		! 3D EZspline object

		INTEGER :: NR, NPHI, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
		TYPE(EZspline2_r8) :: A		! 2D EZspline object
		TYPE(EZspline2_r8) :: R		! 2D EZspline object
		TYPE(EZspline2_r8) :: PHI	! 2D EZspline object
		TYPE(EZspline2_r8) :: Z		! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
		TYPE(EZspline3_r8) :: ne	! 3D EZspline object
		TYPE(EZspline3_r8) :: Te	! 3D EZspline object
		TYPE(EZspline3_r8) :: Zeff	! 3D EZspline object

		INTEGER :: NR, NPHI, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
		TYPE(EZspline2_r8) :: ne	! 2D EZspline object
		TYPE(EZspline2_r8) :: Te	! 2D EZspline object
		TYPE(EZspline2_r8) :: Zeff	! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#elif SINGLE_PRECISION
	TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
		TYPE(EZspline3_r4) :: R		! 3D EZspline object
		TYPE(EZspline3_r4) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r4) :: Z		! 3D EZspline object

		INTEGER, PRIVATE :: NR, NPHI, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
		TYPE(EZspline2_r4) :: A		! 2D EZspline object
		TYPE(EZspline2_r4) :: R		! 2D EZspline object
		TYPE(EZspline2_r4) :: PHI	! 2D EZspline object
		TYPE(EZspline2_r4) :: Z		! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
		TYPE(EZspline3_r4) :: ne	! 3D EZspline object
		TYPE(EZspline3_r4) :: Te	! 3D EZspline object
		TYPE(EZspline3_r4) :: Zeff	! 3D EZspline object

		INTEGER :: NR, NPHI, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
		TYPE(EZspline2_r4) :: ne	! 2D EZspline object
		TYPE(EZspline2_r4) :: Te	! 2D EZspline object
		TYPE(EZspline2_r4) :: Zeff	! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#endif

	TYPE, PRIVATE :: KORC_INTERPOLANT_DOMAIN
		INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE :: FLAG2D
		INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE :: FLAG3D

		REAL(rp) :: Ro
		REAL(rp) :: Zo

		REAL(rp) :: DR
		REAL(rp) :: DPHI
		REAL(rp) :: DZ
	END TYPE 

	TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE :: bfield_2d
	TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE :: bfield_3d
	TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE :: efield_2d
	TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE :: efield_3d
	TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE :: fields_domain
	TYPE(KORC_2D_PROFILES_INTERPOLANT), PRIVATE :: profiles_2d
	TYPE(KORC_3D_PROFILES_INTERPOLANT), PRIVATE :: profiles_3d
	TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE :: profiles_domain
	INTEGER :: ezerr

    PUBLIC :: interp_fields,&
				initialize_fields_interpolant,&
				finalize_interpolants
	PRIVATE :: interp_3D_bfields,&
				interp_2D_bfields,&
				interp_3D_efields,&
				interp_2D_efields,&
				interp_2D_profiles,&
				interp_3D_profiles,&
				check_if_in_fields_domain,&
				check_if_in_profiles_domain

    contains


subroutine initialize_fields_interpolant(params,F)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F

	if (params%plasma_model .EQ. 'EXTERNAL') then
		
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')
		end if


		! * * * * * * * * MAGNETIC FIELD * * * * * * * * !
		if (params%poloidal_flux) then
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
		else if (params%axisymmetric) then
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
		if (params%axisymmetric) then
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


subroutine check_if_in_fields_domain(Y,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: IR,IPHI,IZ
	INTEGER(ip) :: pp,ss

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


subroutine initialize_profiles_interpolant(params,P)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PROFILES), INTENT(INOUT) :: P

	if (params%plasma_model .EQ. 'EXTERNAL') then
		
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * INITIALIZING PROFILES INTERPOLANT * * * *")')
		end if

		if (params%axisymmetric) then
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
end subroutine initialize_profiles_interpolant


subroutine check_if_in_profiles_domain(Y,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: IR,IPHI,IZ
	INTEGER(ip) :: pp,ss

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


subroutine interp_2D_bfields(Y,B,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

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


subroutine interp_3D_bfields(Y,B,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

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


subroutine calculate_magnetic_field(Y,F,B,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz	
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A ! A(1,:) = FR, A(2,:) = FPHI, A(3,:) = FZ
	INTEGER :: pp, ss

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


subroutine interp_2D_efields(Y,E,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: E ! E(1,:) = Ex, E(2,:) = Ey, E(3,:) = Ez
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

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


subroutine interp_3D_efields(Y,E,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: E ! E(1,:) = Ex, E(2,:) = Ey, E(3,:) = Ez
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

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


subroutine interp_fields(prtcls,F)
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F
	INTEGER :: ii, pp, ss ! Iterators

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


subroutine interp_2D_profiles(Y,ne,Te,Zeff,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ne
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Te
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Zeff
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

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


subroutine interp_3D_profiles(Y,ne,Te,Zeff,flag)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ne
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Te
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Zeff
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

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


subroutine interp_profiles(prtcls,P)
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(PROFILES), INTENT(IN) :: P
	INTEGER :: ii, pp, ss ! Iterators

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


subroutine finalize_interpolants(params)
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%plasma_model .EQ. 'EXTERNAL') then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * FINALIZING INTERPOLANT * * * *")')
		end if

		if (EZspline_allocated(bfield_3d%R)) call Ezspline_free(bfield_3d%R, ezerr)
		if (EZspline_allocated(bfield_3d%PHI))call Ezspline_free(bfield_3d%PHI, ezerr)
		if (EZspline_allocated(bfield_3d%Z)) call Ezspline_free(bfield_3d%Z, ezerr)
		if (EZspline_allocated(bfield_2d%A)) call Ezspline_free(bfield_2d%A, ezerr)
		if (EZspline_allocated(bfield_2d%R)) call Ezspline_free(bfield_2d%R, ezerr)
		if (EZspline_allocated(bfield_2d%PHI)) call Ezspline_free(bfield_2d%PHI, ezerr)
		if (EZspline_allocated(bfield_2d%Z)) call Ezspline_free(bfield_2d%Z, ezerr)

		if (EZspline_allocated(profiles_3d%ne)) call Ezspline_free(profiles_3d%ne, ezerr)
		if (EZspline_allocated(profiles_3d%Te))call Ezspline_free(profiles_3d%Te, ezerr)
		if (EZspline_allocated(profiles_3d%Zeff)) call Ezspline_free(profiles_3d%Zeff, ezerr)
		if (EZspline_allocated(profiles_2d%ne)) call Ezspline_free(profiles_2d%ne, ezerr)
		if (EZspline_allocated(profiles_2d%Te))call Ezspline_free(profiles_2d%Te, ezerr)
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
