module korc_interp

    use korc_types
    use korc_coords
    use korc_fields
    use rnd_numbers
    use EZspline_obj	! psplines module
    use EZspline		! psplines module

    implicit none

#ifdef DOUBLE_PRECISION
	TYPE, PRIVATE :: KORC_INTERPOLANT
		TYPE(EZspline3_r8) :: R		! 3D EZspline object
		TYPE(EZspline3_r8) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r8) :: Z		! 3D EZspline object

		INTEGER, PRIVATE :: NR, NPHI, NZ
		INTEGER, PRIVATE, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#elif SINGLE_PRECISION
	TYPE, PRIVATE :: KORC_INTERPOLANT
		TYPE(EZspline3_r4) :: R		! 3D EZspline object
		TYPE(EZspline3_r4) :: PHI	! 3D EZspline object
		TYPE(EZspline3_r4) :: Z		! 3D EZspline object

		INTEGER, PRIVATE :: NR, NPHI, NZ
		INTEGER, PRIVATE, DIMENSION(2) :: BCSR = (/ 0, 0 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
		INTEGER, PRIVATE, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
	END TYPE
#endif
	TYPE(KORC_INTERPOLANT), PRIVATE :: interp
	INTEGER :: ezerr

    PUBLIC :: interp_field, interp_analytical_field, unitVectors,&
				initialize_interpolant

    contains


subroutine initialize_interpolant(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F

	interp%NR = F%dims(1)
	interp%NPHI = F%dims(2)
	interp%NZ = F%dims(3)

	write(6,'("* * * * * * * * * *  * * * * * * * * * *")')
	write(6,'("* * * * INITIALIZING INTERPOLANT * * * *")')

	write(6,'("Initializing R component of interpolant...")')
	call EZspline_init(interp%R, interp%NR, interp%NPHI, interp%NZ,&
						interp%BCSR, interp%BCSPHI, interp%BCSZ, ezerr)
  	call EZspline_error(ezerr)

	interp%R%x1 = F%X%R
!	interp%R%x2 = F%X%PHI
	interp%R%x3 = F%X%Z

	call EZspline_setup(interp%R, F%B%R, ezerr)
	call EZspline_error(ezerr)

	write(6,'("Initializing PHI component of interpolant...")')
	call EZspline_init(interp%PHI, interp%NR, interp%NPHI, interp%NZ,&
						interp%BCSR, interp%BCSPHI, interp%BCSZ, ezerr)
  	call EZspline_error(ezerr)

	interp%PHI%x1 = F%X%R
!	interp%PHI%x2 = F%X%PHI
	interp%PHI%x3 = F%X%Z

	call EZspline_setup(interp%PHI, F%B%PHI, ezerr)
	call EZspline_error(ezerr)

	write(6,'("Initializing Z component of interpolant...")')
	call EZspline_init(interp%Z, interp%NR, interp%NPHI, interp%NZ,&
						interp%BCSR, interp%BCSPHI, interp%BCSZ, ezerr)
  	call EZspline_error(ezerr)

	interp%Z%x1 = F%X%R
!	interp%Z%x2 = F%X%PHI
	interp%Z%x3 = F%X%Z

	call EZspline_setup(interp%Z, F%B%Z, ezerr)
	call EZspline_error(ezerr)


	write(6,'("* * * * INTERPOLANT  INITIALIZED * * * *")')
	write(6,'("* * * * * * * * * *  * * * * * * * * * *")')	

end subroutine initialize_interpolant


subroutine finalize_interpolant()
	implicit none

	call Ezspline_free(interp%R, ezerr)
	call EZspline_error(ezerr)

	call Ezspline_free(interp%PHI, ezerr)
	call EZspline_error(ezerr)

	call Ezspline_free(interp%Z, ezerr)
	call EZspline_error(ezerr)
end subroutine finalize_interpolant


subroutine interp_magnetic_field_once(F,Y,B)
	implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Bcyl ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	TYPE(KORC_INTERPOLANT) :: tmp_interp
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


subroutine interp_magnetic_field(Y,B)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: F ! F(1,:) = FR, F(2,:) = FPHI, F(3,:) = FZ
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(F(3,ss))
!	write(6,'(F20.12,F20.12,F20.12)') Y
!	write(6,*) shape(Y), shape(B), shape(F)
	do pp=1,ss
		call EZspline_interp(interp%R, Y(1,pp), Y(2,pp), Y(3,pp), F(1,pp), ezerr)
		call EZspline_error(ezerr)

		call EZspline_interp(interp%PHI, Y(1,pp), Y(2,pp), Y(3,pp), F(2,pp), ezerr)
		call EZspline_error(ezerr)

		call EZspline_interp(interp%Z, Y(1,pp), Y(2,pp), Y(3,pp), F(3,pp), ezerr)
		call EZspline_error(ezerr)

		B(1,pp) = F(1,pp)*cos(Y(2,pp)) - F(2,pp)*sin(Y(2,pp))
		B(2,pp) = F(1,pp)*sin(Y(2,pp)) + F(2,pp)*cos(Y(2,pp))
		B(3,pp) = F(3,pp)
	end do
!	write(6,'(F15.9,F15.9,F15.9)') B
	
	DEALLOCATE(F)
end subroutine interp_magnetic_field


subroutine interp_field(prtcls,EB)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER :: ii, pp, ss ! Iterators

	call cart_to_cyl(prtcls%X, prtcls%Y)

	call interp_magnetic_field(prtcls%Y,prtcls%B)

end subroutine interp_field


subroutine interp_analytical_field(prtcls,EB)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: EB

	call cart_to_tor(prtcls%X, EB%AB%Ro, prtcls%Y)

	call analytical_magnetic_field(EB,prtcls%Y,prtcls%B)
end subroutine interp_analytical_field


subroutine unitVectors(params,Xo,EB,par,perp)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Xo
	TYPE(FIELDS), INTENT(IN) :: EB
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
		call cart_to_tor(Xo, EB%AB%Ro, X) ! To toroidal coords
		call analytical_magnetic_field(EB,X,B)
	else
		call cart_to_cyl(Xo, X)
		call interp_magnetic_field_once(EB,X,B)
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
