module korc_interp

    use korc_types
    use korc_coords
    use korc_fields
    use rnd_numbers
    use EZspline_obj	! psplines module
    use EZspline		! psplines module

    implicit none

#ifdef DOUBLE_PRECISION
	TYPE, PRIVATE :: KORC_3DINTERPOLANT
		TYPE(EZspline3_r8) :: R		! 3D EZspline object
		TYPE(EZspline3_r8) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r8) :: Z		! 3D EZspline object

		INTEGER, PRIVATE :: NR, NPHI, NZ
		INTEGER, PRIVATE, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2DINTERPOLANT
		TYPE(EZspline2_r8) :: A		! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, PRIVATE, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#elif SINGLE_PRECISION
	TYPE, PRIVATE :: KORC_3DINTERPOLANT
		TYPE(EZspline3_r4) :: R		! 3D EZspline object
		TYPE(EZspline3_r4) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r4) :: Z		! 3D EZspline object

		INTEGER, PRIVATE :: NR, NPHI, NZ
		INTEGER, PRIVATE, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE

	TYPE, PRIVATE :: KORC_2DINTERPOLANT
		TYPE(EZspline2_r4) :: A		! 2D EZspline object

		INTEGER, PRIVATE :: NR, NZ
		INTEGER, PRIVATE, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#endif

	TYPE(KORC_3DINTERPOLANT), PRIVATE :: interp3d
	TYPE(KORC_2DINTERPOLANT), PRIVATE :: interp2d
	INTEGER :: ezerr


	INTERFACE interp_magnetic_field
	  module procedure interp_3D_magnetic_field
	END INTERFACE

    PUBLIC :: interp_field, interp_analytical_field, unitVectors,&
				initialize_interpolant, finalize_interpolant
	PRIVATE :: interp_magnetic_field, interp_3D_magnetic_field, cross,&
				check_if_in_domain2D, check_if_in_domain3D

    contains


function cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


subroutine initialize_interpolant(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F

	if (params%magnetic_field_model .EQ. 'EXTERNAL') then
		
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * INITIALIZING  INTERPOLANT * * * *")')
		end if

		if (.NOT. params%poloidal_flux) then
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

			call EZspline_setup(interp3d%R, F%B%R, ezerr)
			call EZspline_error(ezerr)

			! Initializing PHI component of interpolant
			call EZspline_init(interp3d%PHI, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%PHI%x1 = F%X%R
			! interp3d%PHI%x2 = F%X%PHI
			interp3d%PHI%x3 = F%X%Z

			call EZspline_setup(interp3d%PHI, F%B%PHI, ezerr)
			call EZspline_error(ezerr)

			! Initializing Z component of interpolant
			call EZspline_init(interp3d%Z, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%Z%x1 = F%X%R
			! interp3d%Z%x2 = F%X%PHI
			interp3d%Z%x3 = F%X%Z

			call EZspline_setup(interp3d%Z, F%B%Z, ezerr)
			call EZspline_error(ezerr)
		else
			interp2d%NR = F%dims(1)
			interp2d%NZ = F%dims(3)

			! Initializing poloidal flux interpolant
			call EZspline_init(interp2d%A,interp2d%NR,interp2d%NZ,interp2d%BCSR,&
								interp2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			! interp2d%A%hspline = (/2,2/)
			interp2d%A%x1 = F%X%R
			interp2d%A%x2 = F%X%Z

			call EZspline_setup(interp2d%A, F%PSIp, ezerr, .TRUE.)
			call EZspline_error(ezerr)
		end if

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("* * * * INTERPOLANT   INITIALIZED * * * *")')
			write(6,'("* * * * * * * * * * * * * * * * * * * * *",/)')
		end if
	else if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * *")')
			write(6,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *")')
			write(6,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
		end if
	else if (params%magnetic_field_model .EQ. 'UNIFORM') then
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

	if (params%magnetic_field_model .EQ. 'EXTERNAL') then
		write(6,'("* * * * FINALIZING INTERPOLANT * * * *")')

		if (.NOT. params%poloidal_flux) then
			call Ezspline_free(interp3d%R, ezerr)
			call EZspline_error(ezerr)

			call Ezspline_free(interp3d%PHI, ezerr)
			call EZspline_error(ezerr)

			call Ezspline_free(interp3d%Z, ezerr)
			call EZspline_error(ezerr)	
		else
			call Ezspline_free(interp2d%A, ezerr)
			call EZspline_error(ezerr)
		end if
		write(6,'("* * * * INTERPOLANT  FINALIZED * * * *")')
	end if
end subroutine finalize_interpolant


subroutine check_if_in_domain2D(F,Y,flag)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER(ip) :: pp,ss

    ss = SIZE(Y,2)
!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(F,Y,flag)
!$OMP DO
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_idef ) then
			call EZspline_isInDomain(interp2d%A, Y(1,pp), Y(3,pp), ezerr)
			if (ezerr .NE. 0) then
				flag(pp) = 0_idef
			end if
!			write(6,'("Error code is:",I3)') ezerr
        end if
	end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine check_if_in_domain2D


subroutine check_if_in_domain3D(F,Y,flag)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER(ip) :: pp,ss


    ss = SIZE(Y,2)
!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(F,Y,flag)
!$OMP DO
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_idef ) then
			call EZspline_isInDomain3_r8(interp3d%R, Y(1,pp), Y(2,pp), Y(3,pp), ezerr)
			if (ezerr .NE. 0) then
				flag(pp) = 0_idef
			end if
!			write(6,'("Error code is:",I3)') ezerr
        end if
	end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine check_if_in_domain3D


subroutine interp_3D_magnetic_field(Y,B,flag)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(interp3d,F,Y,B,flag)
!$OMP DO
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_idef ) then
			call EZspline_interp(interp3d%R, Y(1,pp), Y(2,pp), Y(3,pp), F(1,pp), ezerr)
			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_idef
			end if

			call EZspline_interp(interp3d%PHI, Y(1,pp), Y(2,pp), Y(3,pp), F(2,pp), ezerr)
			call EZspline_error(ezerr)

			call EZspline_interp(interp3d%Z, Y(1,pp), Y(2,pp), Y(3,pp), F(3,pp), ezerr)
			call EZspline_error(ezerr)

			B(1,pp) = F(1,pp)*cos(Y(2,pp)) - F(2,pp)*sin(Y(2,pp))
			B(2,pp) = F(1,pp)*sin(Y(2,pp)) + F(2,pp)*cos(Y(2,pp))
			B(3,pp) = F(3,pp)
		end if
	end do
!$OMP END DO
!$OMP END PARALLEL
	DEALLOCATE(F)
end subroutine interp_3D_magnetic_field


subroutine calculate_magnetic_field(Y,F,B,flag)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz	
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A ! A(1,:) = FR, A(2,:) = FPHI, A(3,:) = FZ
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(A(3,ss))
!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(interp2d,F,Y,A,B,flag)
!$OMP DO
	do pp=1_idef,ss
		if ( flag(pp) .EQ. 1_idef ) then
			call EZspline_derivative(interp2d%A, 0, 1, Y(1,pp), Y(3,pp), A(1,pp), ezerr)
!			call EZspline_error(ezerr)

			if (ezerr .NE. 0) then ! We flag the particle as lost
				flag(pp) = 0_idef
			else
				A(1,pp) = A(1,pp)/Y(1,pp)

				A(2,pp) = - F%Bo*F%Ro/Y(1,pp)

				call EZspline_derivative(interp2d%A, 1, 0, Y(1,pp), Y(3,pp), A(3,pp), ezerr)
				call EZspline_error(ezerr)
				A(3,pp) = -A(3,pp)/Y(1,pp)

				B(1,pp) = A(1,pp)*cos(Y(2,pp)) - A(2,pp)*sin(Y(2,pp))
				B(2,pp) = A(1,pp)*sin(Y(2,pp)) + A(2,pp)*cos(Y(2,pp))
				B(3,pp) = A(3,pp)
			end if
		end if
	end do
!$OMP END DO
!$OMP END PARALLEL

!	open(unit=default_unit_write,file='/home/l8c/Documents/KORC/KORC-FO/interpolation.dat',status='UNKNOWN',form='formatted')
!    do pp=1_idef,ss
!	        write(default_unit_write,'(F20.16,T22,F20.16,T44,F20.16)') Y(1,pp), Y(3,pp), A(1,pp)
!    end do
!    close(default_unit_write)
!    call korc_abort()
	DEALLOCATE(A)
end subroutine calculate_magnetic_field


subroutine interp_field(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F
	INTEGER :: ii, pp, ss ! Iterators

	call cart_to_cyl(prtcls%X, prtcls%Y)

	if (.NOT. ALLOCATED(F%PSIp)) then
		call check_if_in_domain3D(F, prtcls%Y, prtcls%flag)

		call interp_magnetic_field(prtcls%Y,prtcls%B,prtcls%flag)
	else
		call check_if_in_domain2D(F, prtcls%Y, prtcls%flag)

		call calculate_magnetic_field(prtcls%Y,F,prtcls%B,prtcls%flag)
	end if
end subroutine interp_field


subroutine interp_analytical_field(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F

	call cart_to_tor(prtcls%X, F%AB%Ro, prtcls%Y, prtcls%flag)

	call check_if_confined(F, prtcls%Y, prtcls%flag)

	call analytical_magnetic_field(F, prtcls%Y, prtcls%B, prtcls%flag)

	call analytical_electric_field(F, prtcls%Y, prtcls%E, prtcls%flag)
end subroutine interp_analytical_field


subroutine uniform_fields(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F

	call uniform_magnetic_field(F, prtcls%B)

	call uniform_electric_field(F, prtcls%E)
end subroutine uniform_fields


subroutine unitVectors(params,Xo,F,b1,b2,b3,flag)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Xo
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: b1
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: b2
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: b3
	INTEGER, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: flag
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: B
	INTEGER, DIMENSION(:), ALLOCATABLE :: local_flag
	REAL(rp), PARAMETER :: tol = korc_zero
	INTEGER :: ii, ppp

	ppp = SIZE(Xo,2) ! Number of particles

	ALLOCATE( X(3,ppp) )
	ALLOCATE( B(3,ppp) )
    ALLOCATE( local_flag(ppp) )

    local_flag = 1_idef
	
	call init_random_seed()

	if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		call cart_to_tor(Xo,F%AB%Ro,X,local_flag) ! To toroidal coords
		call analytical_magnetic_field(F,X,B,local_flag)
	else if (params%magnetic_field_model .EQ. 'UNIFORM') then
		call uniform_magnetic_field(F,B)
	else
		call cart_to_cyl(Xo, X)
        if (params%poloidal_flux) then
			call check_if_in_domain2D(F,X,local_flag)

            call calculate_magnetic_field(X,F,B,local_flag)
        else
			call check_if_in_domain3D(F,X,local_flag)

            call interp_magnetic_field(X,B,local_flag)
        end if
	end if
	
	do ii=1_idef,ppp
		if ( local_flag(ii) .EQ. 1_idef ) then
			b1(:,ii) = B(:,ii)/sqrt( DOT_PRODUCT(B(:,ii),B(:,ii)) )

		    b2(:,ii) = cross(b1(:,ii),(/0.0_rp,0.0_rp,1.0_rp/))
		    b2(:,ii) = b2(:,ii)/sqrt( DOT_PRODUCT(b2(:,ii),b2(:,ii)) )

		    b3(:,ii) = cross(b1(:,ii),b2(:,ii))
		    b3(:,ii) = b3(:,ii)/sqrt( DOT_PRODUCT(b3(:,ii),b3(:,ii)) )
		end if
	end do

	if (PRESENT(flag)) then
		flag = local_flag
	end if

	DEALLOCATE(X)
	DEALLOCATE(B)
    DEALLOCATE(local_flag)
end subroutine unitVectors


end module korc_interp
