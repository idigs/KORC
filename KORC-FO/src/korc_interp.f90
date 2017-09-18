module korc_interp

    use korc_types
    use korc_coords
    use rnd_numbers
    use EZspline_obj	! psplines module
    use EZspline		! psplines module

    implicit none

#ifdef DOUBLE_PRECISION
	TYPE, PRIVATE :: KORC_3DINTERPOLANT
		TYPE(EZspline3_r8) :: R		! 3D EZspline object
		TYPE(EZspline3_r8) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r8) :: Z		! 3D EZspline object

		INTEGER :: NR, NPHI, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2DINTERPOLANT
		TYPE(EZspline2_r8) :: A		! 2D EZspline object
		TYPE(EZspline2_r8) :: R		! 2D EZspline object
		TYPE(EZspline2_r8) :: PHI	! 2D EZspline object
		TYPE(EZspline2_r8) :: Z		! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#elif SINGLE_PRECISION
	TYPE, PRIVATE :: KORC_3DINTERPOLANT
		TYPE(EZspline3_r4) :: R		! 3D EZspline object
		TYPE(EZspline3_r4) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r4) :: Z		! 3D EZspline object

		INTEGER, PRIVATE :: NR, NPHI, NZ
		INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2DINTERPOLANT
		TYPE(EZspline2_r4) :: A		! 2D EZspline object
		TYPE(EZspline2_r4) :: R		! 2D EZspline object
		TYPE(EZspline2_r4) :: PHI	! 2D EZspline object
		TYPE(EZspline2_r4) :: Z		! 2D EZspline object

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

	TYPE(KORC_3DINTERPOLANT), PRIVATE :: interp3d
	TYPE(KORC_2DINTERPOLANT), PRIVATE :: interp2d
	TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE :: domain
	INTEGER :: ezerr


!	INTERFACE interp_magnetic_field
!	  module procedure interp_3D_B_field
!	END INTERFACE

    PUBLIC :: interp_field,initialize_interpolant,finalize_interpolant
	PRIVATE :: interp_3D_B_field,interp_2D_B_field,check_if_in_domain2D,check_if_in_domain3D,check_if_in_domain

    contains


subroutine initialize_interpolant(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F

	if (params%plasma_model .EQ. 'EXTERNAL') then
		
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * INITIALIZING  INTERPOLANT * * * *")')
		end if

		if ((.NOT.params%poloidal_flux).AND.(.NOT.params%axisymmetric)) then
			interp3d%NR = F%dims(1)
			interp3d%NPHI = F%dims(2)
			interp3d%NZ = F%dims(3)

			! Initializing R component of interpolant
			call EZspline_init(interp3d%R, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%R%x1 = F%X%R
			! interp3d%R%x2 = F%X%PHI
			interp3d%R%x3 = F%X%Z

			call EZspline_setup(interp3d%R, F%B_3D%R, ezerr)
			call EZspline_error(ezerr)

			! Initializing PHI component of interpolant
			call EZspline_init(interp3d%PHI, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%PHI%x1 = F%X%R
			! interp3d%PHI%x2 = F%X%PHI
			interp3d%PHI%x3 = F%X%Z

			call EZspline_setup(interp3d%PHI, F%B_3D%PHI, ezerr)
			call EZspline_error(ezerr)

			! Initializing Z component of interpolant
			call EZspline_init(interp3d%Z, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%Z%x1 = F%X%R
			! interp3d%Z%x2 = F%X%PHI
			interp3d%Z%x3 = F%X%Z

			call EZspline_setup(interp3d%Z, F%B_3D%Z, ezerr)
			call EZspline_error(ezerr)

			ALLOCATE(domain%FLAG3D(interp3d%NR,interp3d%NPHI,interp3d%NZ))
			domain%FLAG3D = F%FLAG3D

			domain%DR = ABS(F%X%R(2) - F%X%R(1))
			domain%DPHI = 2.0_rp*C_PI/interp3d%NPHI
			domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
		else if (params%poloidal_flux) then
			interp2d%NR = F%dims(1)
			interp2d%NZ = F%dims(3)

			! Initializing poloidal flux interpolant
			call EZspline_init(interp2d%A,interp2d%NR,interp2d%NZ,interp2d%BCSR,interp2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			interp2d%A%x1 = F%X%R
			interp2d%A%x2 = F%X%Z

			call EZspline_setup(interp2d%A, F%PSIp, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			ALLOCATE(domain%FLAG2D(interp3d%NR,interp3d%NZ))
			domain%FLAG2D = F%FLAG2D

			domain%DR = ABS(F%X%R(2) - F%X%R(1))
			domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
		else if (params%axisymmetric) then
			interp2d%NR = F%dims(1)
			interp2d%NZ = F%dims(3)

			! Initializing R component
			call EZspline_init(interp2d%R,interp2d%NR,interp2d%NZ,interp2d%BCSR,interp2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)
			
			interp2d%R%x1 = F%X%R
			interp2d%R%x2 = F%X%Z

			call EZspline_setup(interp2d%R, F%B_2D%R, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			! Initializing PHI component
			call EZspline_init(interp2d%PHI,interp2d%NR,interp2d%NZ,interp2d%BCSR,interp2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)
			
			interp2d%PHI%x1 = F%X%R
			interp2d%PHI%x2 = F%X%Z

			call EZspline_setup(interp2d%PHI, F%B_2D%PHI, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			! Initializing Z component
			call EZspline_init(interp2d%Z,interp2d%NR,interp2d%NZ,interp2d%BCSR,interp2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)
			
			interp2d%Z%x1 = F%X%R
			interp2d%Z%x2 = F%X%Z

			call EZspline_setup(interp2d%Z, F%B_2D%Z, ezerr, .TRUE.)
			call EZspline_error(ezerr)

			ALLOCATE(domain%FLAG2D(interp3d%NR,interp3d%NZ))
			domain%FLAG2D = F%FLAG2D

			domain%DR = ABS(F%X%R(2) - F%X%R(1))
			domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
		end if

		domain%Ro = F%X%R(1)
		domain%Zo = F%X%Z(1)

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * INTERPOLANT   INITIALIZED * * * *")')
			write(6,'("* * * * * * * * * * * * * * * * * * * * *",/)')
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
			write(6,'(/,"* * * * * * * * * * *  * * * * * * * * * * *")')
		end if
	end if
end subroutine initialize_interpolant


subroutine finalize_interpolant(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%plasma_model .EQ. 'EXTERNAL') then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * FINALIZING INTERPOLANT * * * *")')
		end if

		if (EZspline_allocated(interp3d%R)) call Ezspline_free(interp3d%R, ezerr)
		if (EZspline_allocated(interp3d%PHI))call Ezspline_free(interp3d%PHI, ezerr)
		if (EZspline_allocated(interp3d%Z)) call Ezspline_free(interp3d%Z, ezerr)
		if (EZspline_allocated(interp2d%A)) call Ezspline_free(interp2d%A, ezerr)
		if (EZspline_allocated(interp2d%R)) call Ezspline_free(interp2d%R, ezerr)
		if (EZspline_allocated(interp2d%PHI)) call Ezspline_free(interp2d%PHI, ezerr)
		if (EZspline_allocated(interp2d%Z)) call Ezspline_free(interp2d%Z, ezerr)
!		call EZspline_error(ezerr)

		if (ALLOCATED(domain%FLAG2D)) DEALLOCATE(domain%FLAG2D)
		if (ALLOCATED(domain%FLAG3D)) DEALLOCATE(domain%FLAG3D)

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * INTERPOLANT  FINALIZED * * * *")')
		end if
	end if
end subroutine finalize_interpolant


subroutine check_if_in_domain(Y,flag)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: IR,IPHI,IZ
	INTEGER(ip) :: pp,ss

	ss = SIZE(Y,2)

	if (ALLOCATED(domain%FLAG3D)) then
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) SHARED(Y,flag,domain,interp3d)
		do pp=1_idef,ss
			IR = INT(FLOOR((Y(1,pp)  - domain%Ro + 0.5_rp*domain%DR)/domain%DR) + 1.0_rp,idef)
			IPHI = INT(FLOOR((Y(2,pp)  + 0.5_rp*domain%DPHI)/domain%DPHI) + 1.0_rp,idef)
			IZ = INT(FLOOR((Y(3,pp)  + ABS(domain%Zo) + 0.5_rp*domain%DZ)/domain%DZ) + 1.0_rp,idef)
	
			if ((domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR.((IR.GT.interp3d%NR).OR.(IZ.GT.interp3d%NZ))) then
				flag(pp) = 0_is
			end if
		end do
!$OMP END PARALLEL DO
	else 
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) SHARED(Y,flag,domain,interp2d)
		do pp=1_idef,ss
			IR = INT(FLOOR((Y(1,pp)  - domain%Ro + 0.5_rp*domain%DR)/domain%DR) + 1.0_rp,idef)
			IZ = INT(FLOOR((Y(3,pp)  + ABS(domain%Ro) + 0.5_rp*domain%DZ)/domain%DZ) + 1.0_rp,idef)

			if ((domain%FLAG2D(IR,IZ).NE.1_is).OR.((IR.GT.interp2d%NR).OR.(IZ.GT.interp2d%NZ))) then
				flag(pp) = 0_is
			end if
		end do
!$OMP END PARALLEL DO
	end if
end subroutine check_if_in_domain


subroutine check_if_in_domain2D(Y,flag)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER(ip) :: pp,ss

    ss = SIZE(Y,2)
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(Y,flag,interp2d)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
			if (EZspline_allocated(interp2d%A)) call EZspline_isInDomain(interp2d%A, Y(1,pp), Y(3,pp), ezerr)
			if (EZspline_allocated(interp2d%R)) call EZspline_isInDomain(interp2d%R, Y(1,pp), Y(3,pp), ezerr)

			if (ezerr .NE. 0) then
				flag(pp) = 0_is
			end if
!			write(6,'("Error code is:",I3)') ezerr
        end if
	end do
!$OMP END PARALLEL DO
end subroutine check_if_in_domain2D


subroutine check_if_in_domain3D(Y,flag)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER(ip) :: pp,ss


    ss = SIZE(Y,2)
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(Y,flag,interp3d)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
			call EZspline_isInDomain3_r8(interp3d%R, Y(1,pp), Y(2,pp), Y(3,pp), ezerr)
			if (ezerr .NE. 0) then
				flag(pp) = 0_is
			end if
!			write(6,'("Error code is:",I3)') ezerr
        end if
	end do
!$OMP END PARALLEL DO
end subroutine check_if_in_domain3D


subroutine interp_2D_B_field(Y,B,flag)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

	ss = size(Y,2)


	ALLOCATE(F(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,B,flag,interp2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(interp2d%R, Y(1,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(interp2d%PHI, Y(1,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(interp2d%Z, Y(1,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			B(1,pp) = F(1,pp)*COS(Y(2,pp)) - F(2,pp)*SIN(Y(2,pp))
			B(2,pp) = F(1,pp)*SIN(Y(2,pp)) + F(2,pp)*COS(Y(2,pp))
			B(3,pp) = F(3,pp)
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(F)
end subroutine interp_2D_B_field


subroutine interp_3D_B_field(Y,B,flag)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,B,flag,interp3d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_interp(interp3d%R, Y(1,pp), Y(2,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			end if

			call EZspline_interp(interp3d%PHI, Y(1,pp), Y(2,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(interp3d%Z, Y(1,pp), Y(2,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			B(1,pp) = F(1,pp)*COS(Y(2,pp)) - F(2,pp)*SIN(Y(2,pp))
			B(2,pp) = F(1,pp)*SIN(Y(2,pp)) + F(2,pp)*COS(Y(2,pp))
			B(3,pp) = F(3,pp)
		end if
	end do
!$OMP END PARALLEL DO
	DEALLOCATE(F)
end subroutine interp_3D_B_field


subroutine calculate_magnetic_field(Y,F,B,flag)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz	
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A ! A(1,:) = FR, A(2,:) = FPHI, A(3,:) = FZ
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(A(3,ss))
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(F,Y,A,B,flag,interp2d)
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_is ) then
			call EZspline_derivative(interp2d%A, 0, 1, Y(1,pp), Y(3,pp), A(1,pp), ezerr)
!			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_is
			else
				A(1,pp) = A(1,pp)/Y(1,pp)

				A(2,pp) = - F%Bo*F%Ro/Y(1,pp)

				call EZspline_derivative(interp2d%A, 1, 0, Y(1,pp), Y(3,pp), A(3,pp), ezerr)
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


subroutine interp_field(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F
	INTEGER :: ii, pp, ss ! Iterators

	call cart_to_cyl(prtcls%X, prtcls%Y)

	call check_if_in_domain(prtcls%Y, prtcls%flag)

	if (ALLOCATED(F%PSIp)) then
		call calculate_magnetic_field(prtcls%Y,F,prtcls%B,prtcls%flag)
	else if (ALLOCATED(F%B_2D%R)) then
		call interp_2D_B_field(prtcls%Y,prtcls%B,prtcls%flag)
	else if (ALLOCATED(F%B_3D%R)) then
		call interp_3D_B_field(prtcls%Y,prtcls%B,prtcls%flag)
	else
		write(6,'("* * * * NO FIELD ARRAYS ALLOCATED * * * *")')
	end if
end subroutine interp_field

end module korc_interp
