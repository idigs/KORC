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
	PRIVATE :: interp_magnetic_field, interp_3D_magnetic_field

    contains


subroutine initialize_interpolant(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F

	if (params%magnetic_field_model .EQ. 'EXTERNAL') then

		write(6,'("* * * * INITIALIZING INTERPOLANT * * * *")')

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
			call EZspline_init(interp2d%A,interp2d%NR,interp2d%NZ,interp2d%BCSR,interp2d%BCSZ,ezerr)
		  	call EZspline_error(ezerr)

			! interp2d%A%hspline = (/2,2/)
			interp2d%A%x1 = F%X%R
			interp2d%A%x2 = F%X%Z

			call EZspline_setup(interp2d%A, F%PSIp, ezerr, .TRUE.)
			call EZspline_error(ezerr)
		end if

		write(6,'("* * * * INTERPOLANT  INITIALIZED * * * *")')
	else
		write(6,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *")')
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


subroutine calculate_magnetic_field(Y,F,B)
	implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = R, Y(2,:) = PHI, Y(3,:) = Z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A ! A(1,:) = FR, A(2,:) = FPHI, A(3,:) = FZ
	INTEGER :: pp, ss

	ss = size(Y,2)

	ALLOCATE(A(3,ss))
!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) SHARED(interp2d,F,Y,A,B)
!$OMP DO
	do pp=1,ss
		call EZspline_derivative(interp2d%A, 0, 1, Y(1,pp), Y(3,pp), A(1,pp), ezerr)
		call EZspline_error(ezerr)
		A(1,pp) = A(1,pp)/Y(1,pp)

		A(2,pp) = - F%Bo*F%Ro/Y(1,pp)

		call EZspline_derivative(interp2d%A, 1, 0, Y(1,pp), Y(3,pp), A(3,pp), ezerr)
		call EZspline_error(ezerr)
		A(3,pp) = -A(3,pp)/Y(1,pp)

		B(1,pp) = A(1,pp)*cos(Y(2,pp)) - A(2,pp)*sin(Y(2,pp))
		B(2,pp) = A(1,pp)*sin(Y(2,pp)) + A(2,pp)*cos(Y(2,pp))
		B(3,pp) = A(3,pp)
	end do
!$OMP END DO
!$OMP END PARALLEL

!	open(unit=default_unit_write,file='/home/l8c/Documents/KORC/KORC-FO/interpolation.dat',status='UNKNOWN',form='formatted')
!    do pp=1,ss
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
		call interp_magnetic_field(prtcls%Y,prtcls%B)
	else
		call calculate_magnetic_field(prtcls%Y,F,prtcls%B)
	end if
end subroutine interp_field


subroutine interp_analytical_field(prtcls,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: F

	call cart_to_tor(prtcls%X, F%AB%Ro, prtcls%Y)

	call analytical_magnetic_field(F,prtcls%Y,prtcls%B)

	call analytical_electric_field(F,prtcls%Y,prtcls%E)
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
        if (params%poloidal_flux) then
            call calculate_magnetic_field(X,F,B)
        else
            call interp_magnetic_field(X,B)
        end if
	end if

	call RANDOM_NUMBER(rnd_num)
	
	do ii=1,ppp
		par(:,ii) = B(:,ii)/sqrt( DOT_PRODUCT(B(:,ii),B(:,ii)) )
		if ( ALL( ABS(par(:,ii)) .GE. korc_zero ) ) then
			c1 = par(1,ii)**2 + par(2,ii)**2
			
			az = sqrt(c1)*rnd_num(ii)
		
			c2 = 2.0_rp*par(1,ii)*par(3,ii)*az
			c3 = (par(2,ii)**2 + par(3,ii)**2)*az**2 -par(2,ii)**2

			ax = ( -c2 + sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
			if ( (1.0_rp - (ax**2 + az**2)) .GE. korc_zero ) then
				ay = sqrt(1.0_rp - ax**2 - az**2)
			else
				ay = 0.0_rp
			end if
			perp(:,ii) = (/ax, ay, az/)

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				ax = ( -c2 - sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
				if ( (1.0_rp - (ax**2 + az**2)) .GE. korc_zero ) then
					ay = sqrt(1.0_rp - ax**2 - az**2)
				else
					ay = 0.0_rp
				end if
				perp(:,ii) = (/ax, ay, az/)
			end if
		
			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				az = -az
				
				c2 = 2.0_rp*par(1,ii)*par(3,ii)*az
				c3 = (par(2,ii)**2 + par(3,ii)**2)*az**2 -par(2,ii)**2
				
				ax = ( -c2 + sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
				if ( (1.0_rp - (ax**2 + az**2)) .GE. korc_zero ) then
					ay = sqrt(1.0_rp - ax**2 - az**2)
				else
					ay = 0.0_rp
				end if				
				perp(:,ii) = (/ax, ay, az/)
				
				if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
					ax = ( -c2 - sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
					if ( (1.0_rp - (ax**2 + az**2)) .GE. korc_zero ) then
						ay = sqrt(1.0_rp - ax**2 - az**2)
					else
						ay = 0.0_rp
				    end if
				    perp(:,ii) = (/ax, ay, az/)
				end if
			end if
		else if ( ABS(par(3,ii)) .LE. korc_zero ) then
!			write(6,'("Exception: bz .LE. 0",3F20.16)') par(:,ii)
			ax = 0.99_rp*rnd_num(ii)
			ay = -ax*par(1,ii)/par(2,ii)

			az = sqrt( 1 - ax**2 - ay**2 )
			perp(:,ii) = (/ax, ay, az/)

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				az = -sqrt( 1 - ax**2 - ay**2 )
				perp(:,ii) = (/ax, ay, az/)
			end if

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				write(6,'("Exception: bz .LE. 0",3F20.16)') par(:,ii)
			end if
		else if ( ABS(par(2,ii)) .LE. korc_zero ) then
!			write(6,'("Exception: by .LE. 0",3F20.16)') par(:,ii)
			az = 0.99_rp*rnd_num(ii)
			ax = -az*par(3,ii)/par(1,ii)

			ay = sqrt( 1 - ax**2 - az**2 )
			perp(:,ii) = (/ax, ay, az/)

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				ay = -sqrt( 1 - ax**2 - az**2 )
				perp(:,ii) = (/ax, ay, az/)
			end if

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				write(6,'("Exception: by .LE. 0",3F20.16)') par(:,ii)
			end if
		else if ( ABS(par(1,ii)) .LE. korc_zero ) then
!			write(6,'("Exception: bx .LE. 0",3F20.16)') par(:,ii)
			az = 0.99_rp*rnd_num(ii)
			ay = -az*par(3,ii)/par(2,ii)

			ax = sqrt( 1 - ay**2 - az**2 )
			perp(:,ii) = (/ax, ay, az/)

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				ax = -sqrt( 1 - ay**2 - az**2 )
				perp(:,ii) = (/ax, ay, az/)
			end if

			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				write(6,'("Exception: bx .LE. 0",3F20.16)') par(:,ii)
			end if
		end if
	end do

	DEALLOCATE(X)
	DEALLOCATE(B)
	DEALLOCATE(rnd_num)
end subroutine unitVectors

end module korc_interp
