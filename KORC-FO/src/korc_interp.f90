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
				initialize_interpolant, finalize_interpolant,&
				interp_magnetic_field_once
	PRIVATE :: interp_magnetic_field, interp_3D_magnetic_field

    contains


subroutine initialize_interpolant(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F

	if (params%magnetic_field_model .EQ. 'EXTERNAL') then

		write(6,'("* * * * * * * * * *  * * * * * * * * * *")')
		write(6,'("* * * * INITIALIZING INTERPOLANT * * * *")')
		write(6,*)

		if (.NOT. params%poloidal_flux) then
			interp3d%NR = F%dims(1)
			interp3d%NPHI = F%dims(2)
			interp3d%NZ = F%dims(3)

			write(6,'("Initializing R component of interpolant...")')
			call EZspline_init(interp3d%R, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%R%x1 = F%X%R
		!	interp3d%R%x2 = F%X%PHI
			interp3d%R%x3 = F%X%Z

			call EZspline_setup(interp3d%R, F%B%R, ezerr)
			call EZspline_error(ezerr)

			write(6,'("Initializing PHI component of interpolant...")')
			call EZspline_init(interp3d%PHI, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%PHI%x1 = F%X%R
		!	interp3d%PHI%x2 = F%X%PHI
			interp3d%PHI%x3 = F%X%Z

			call EZspline_setup(interp3d%PHI, F%B%PHI, ezerr)
			call EZspline_error(ezerr)

			write(6,'("Initializing Z component of interpolant...")')
			call EZspline_init(interp3d%Z, interp3d%NR, interp3d%NPHI, interp3d%NZ,&
								interp3d%BCSR, interp3d%BCSPHI, interp3d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp3d%Z%x1 = F%X%R
		!	interp3d%Z%x2 = F%X%PHI
			interp3d%Z%x3 = F%X%Z

			call EZspline_setup(interp3d%Z, F%B%Z, ezerr)
			call EZspline_error(ezerr)
		else
			interp2d%NR = F%dims(1)
			interp2d%NZ = F%dims(3)

			write(6,'("Initializing poloidal flux interpolant...")')
			call EZspline_init(interp2d%A, interp2d%NR, interp2d%NZ,&
								interp2d%BCSR, interp2d%BCSZ, ezerr)
		  	call EZspline_error(ezerr)

			interp2d%A%x1 = F%X%R
			interp2d%A%x2 = F%X%Z

			call EZspline_setup(interp2d%A, F%PSIp, ezerr)
			call EZspline_error(ezerr)
		end if

		write(6,'("* * * * INTERPOLANT  INITIALIZED * * * *")')
		write(6,'("* * * * * * * * * *  * * * * * * * * * *")')	
	else
		write(6,'("* * * * * * * * * * * * * * * * * * * * * * * *")')
		write(6,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *")')
		write(6,'("* * * * * * * * * * * * * * * * * * * * * * * *")')
	end if
end subroutine initialize_interpolant


subroutine finalize_interpolant(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%magnetic_field_model .EQ. 'EXTERNAL') then
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
	end if
end subroutine finalize_interpolant


subroutine interp_magnetic_field_once(F,Y,B)
	implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Bcyl ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	TYPE(KORC_3DINTERPOLANT) :: tmp_interp
	INTEGER :: pp, ss

	tmp_interp%NR = F%dims(1)
	tmp_interp%NPHI = F%dims(2)
	tmp_interp%NZ = F%dims(3)

	call EZspline_init(tmp_interp%R, tmp_interp%NR, tmp_interp%NPHI, tmp_interp%NZ,&
						tmp_interp%BCSR, tmp_interp%BCSPHI, tmp_interp%BCSZ, ezerr)
  	call EZspline_error(ezerr)

	tmp_interp%R%x1 = F%X%R
!	tmp_interp%R%x2 = F%X%PHI
	tmp_interp%R%x3 = F%X%Z

	call EZspline_setup(tmp_interp%R, F%B%R, ezerr)
	call EZspline_error(ezerr)

	call EZspline_init(tmp_interp%PHI, tmp_interp%NR, tmp_interp%NPHI, tmp_interp%NZ,&
						tmp_interp%BCSR, tmp_interp%BCSPHI, tmp_interp%BCSZ, ezerr)
  	call EZspline_error(ezerr)

	tmp_interp%PHI%x1 = F%X%R
!	tmp_interp%PHI%x2 = F%X%PHI
	tmp_interp%PHI%x3 = F%X%Z

	call EZspline_setup(tmp_interp%PHI, F%B%PHI, ezerr)
	call EZspline_error(ezerr)

	call EZspline_init(tmp_interp%Z, tmp_interp%NR, tmp_interp%NPHI, tmp_interp%NZ,&
						tmp_interp%BCSR, tmp_interp%BCSPHI, tmp_interp%BCSZ, ezerr)
  	call EZspline_error(ezerr)

	tmp_interp%Z%x1 = F%X%R
!	tmp_interp%Z%x2 = F%X%PHI
	tmp_interp%Z%x3 = F%X%Z

	call EZspline_setup(tmp_interp%Z, F%B%Z, ezerr)
	call EZspline_error(ezerr)

	! Interpolation

	ss = size(Y,2)

	ALLOCATE(Bcyl(3,ss))

	do pp=1,ss
		call EZspline_interp(tmp_interp%R, Y(1,pp), Y(2,pp), Y(3,pp), Bcyl(1,pp), ezerr)
		call EZspline_error(ezerr)

		call EZspline_interp(tmp_interp%PHI, Y(1,pp), Y(2,pp), Y(3,pp), Bcyl(2,pp), ezerr)
		call EZspline_error(ezerr)

		call EZspline_interp(tmp_interp%Z, Y(1,pp), Y(2,pp), Y(3,pp), Bcyl(3,pp), ezerr)
		call EZspline_error(ezerr)

		B(1,pp) = Bcyl(1,pp)*cos(Y(2,pp)) - Bcyl(2,pp)*sin(Y(2,pp))
		B(2,pp) = Bcyl(1,pp)*sin(Y(2,pp)) + Bcyl(2,pp)*cos(Y(2,pp))
		B(3,pp) = Bcyl(3,pp)
	end do
	
	DEALLOCATE(Bcyl)

	! Interpolation

	call Ezspline_free(tmp_interp%R, ezerr)
	call EZspline_error(ezerr)

	call Ezspline_free(tmp_interp%PHI, ezerr)
	call EZspline_error(ezerr)

	call Ezspline_free(tmp_interp%Z, ezerr)
	call EZspline_error(ezerr)
end subroutine interp_magnetic_field_once


subroutine interp_3D_magnetic_field(Y,B)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(interp3d,F,Y,B)
!$OMP DO
	do pp=1,ss
		call EZspline_interp(interp3d%R, Y(1,pp), Y(2,pp), Y(3,pp), F(1,pp), ezerr)
		call EZspline_error(ezerr)
!		if (ezerr .NE. 0) then
!			write(6,'("ezerr value: ",I15)') ezerr ! 97 outside interpolation domain
!		end if

		call EZspline_interp(interp3d%PHI, Y(1,pp), Y(2,pp), Y(3,pp), F(2,pp), ezerr)
		call EZspline_error(ezerr)

		call EZspline_interp(interp3d%Z, Y(1,pp), Y(2,pp), Y(3,pp), F(3,pp), ezerr)
		call EZspline_error(ezerr)

		B(1,pp) = F(1,pp)*cos(Y(2,pp)) - F(2,pp)*sin(Y(2,pp))
		B(2,pp) = F(1,pp)*sin(Y(2,pp)) + F(2,pp)*cos(Y(2,pp))
		B(3,pp) = F(3,pp)
	end do
!$OMP END DO
!$OMP END PARALLEL
	DEALLOCATE(F)
end subroutine interp_3D_magnetic_field


subroutine interp_field(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F
	INTEGER :: ii, pp, ss ! Iterators

	call cart_to_cyl(prtcls%X, prtcls%Y)

	if (.NOT. ALLOCATED(F%PSIp)) then
		call interp_magnetic_field(prtcls%Y,prtcls%B)
	else
		! something different
	end if

end subroutine interp_field


subroutine interp_analytical_field(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F

	call cart_to_tor(prtcls%X, F%AB%Ro, prtcls%Y)

	call analytical_magnetic_field(F,prtcls%Y,prtcls%B)
end subroutine interp_analytical_field


subroutine unitVectors(params,Xo,F,par,perp)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Xo
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: par
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: perp
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: B
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rnd_num
	REAL(rp), PARAMETER :: tol = korc_zero
	REAL(rp) :: c1, c2, c3
	REAL(rp) :: ax, ay, az
	INTEGER :: ppp
	INTEGER :: ii


	ppp = SIZE(Xo,2) ! Number of particles

	ALLOCATE( X(3,ppp) )
	ALLOCATE( B(3,ppp) )
	ALLOCATE( rnd_num(ppp) )
	
	call init_random_seed()

	if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		call cart_to_tor(Xo, F%AB%Ro, X) ! To toroidal coords
		call analytical_magnetic_field(F,X,B)
	else
		call cart_to_cyl(Xo, X)
		call interp_magnetic_field_once(F,X,B)
	end if

	call RANDOM_NUMBER(rnd_num)
	
	do ii=1,ppp
		par(:,ii) = B(:,ii)/sqrt( sum(B(:,ii)**2) )
		if ( ALL( ABS(par(:,ii)) .GE. korc_zero ) ) then
			c1 = par(1,ii)**2 + par(2,ii)**2
			
			az = sqrt(c1)*rnd_num(ii)
		
			c2 = 2.0_rp*par(1,ii)*par(3,ii)*az
			c3 = (par(2,ii)**2 + par(3,ii)**2)*az**2 -par(2,ii)**2
		
			ax = ( -c2 + sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
			ay = sqrt(1.0_rp - ax**2 - az**2)
		
			perp(:,ii) = (/ax, ay, az/)
		
			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				ax = ( -c2 - sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
				ay = sqrt(1.0_rp - ax**2 - az**2)
				
				perp(:,ii) = (/ax, ay, az/)
			end if
		
			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				az = -az
				
				c2 = 2.0_rp*par(1,ii)*par(3,ii)*az
				c3 = (par(2,ii)**2 + par(3,ii)**2)*az**2 -par(2,ii)**2
				
				ax = ( -c2 + sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
				ay = sqrt(1.0_rp - ax**2 - az**2)
				
				perp(:,ii) = (/ax, ay, az/)
				
				if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
					ax = ( -c2 - sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
					ay = sqrt(1.0_rp - ax**2 - az**2)
				    
				    perp(:,ii) = (/ax, ay, az/)
				end if
			end if
		else if ( par(1,ii) .LE. korc_zero ) then
			write(6,'("Exception!")')
		else if ( par(2,ii) .LE. korc_zero ) then
			write(6,'("Exception!")')
		else if ( par(3,ii) .LE. korc_zero ) then
			write(6,'("Exception!")')
		end if
	end do


	DEALLOCATE(X)
	DEALLOCATE(B)
	DEALLOCATE(rnd_num)
end subroutine unitVectors

end module korc_interp
