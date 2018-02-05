MODULE korc_spatial_distribution
	USE korc_types
	USE korc_constants
	USE korc_HDF5
	USE korc_hpc
    use korc_fields
    use korc_rnd_numbers
	use korc_hammersley_generator

	IMPLICIT NONE

	PUBLIC :: intitial_spatial_distribution
	PRIVATE :: uniform,&
				disk,&
				torus,&
				elliptic_torus,&
				exponential_elliptic_torus,&
				gaussian_elliptic_torus,&
				exponential_torus,&
				gaussian_torus,&
				gaussian_energy,&
				calculate_gaussian_falloff,&
				exponential_energy,&
				calculate_exponential_falloff,&
				fzero

	CONTAINS


subroutine uniform(spp)
	TYPE(SPECIES), INTENT(INOUT) :: spp

	spp%vars%X = 0.0_rp
end subroutine uniform


subroutine disk(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, r ! temporary vars
	INTEGER :: jj ! Iterator

	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	! Uniform distribution on a disk at a fixed azimuthal theta		
	call init_random_seed()
	call RANDOM_NUMBER(r)

	r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*COS(spp%PHIo)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*SIN(spp%PHIo)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(theta)
	DEALLOCATE(r)
end subroutine disk


subroutine torus(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, r ! temporary vars
	INTEGER :: jj ! Iterator

	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta		
	call init_random_seed()
	call RANDOM_NUMBER(r)

	r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine torus


subroutine elliptic_torus(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rotation_angle, theta, zeta, r ! temporary vars
	REAL(rp), DIMENSION(:), ALLOCATABLE :: X,Y,X1,Y1
	INTEGER :: jj ! Iterator

	ALLOCATE(X1(spp%ppp))
	ALLOCATE(Y1(spp%ppp))
	ALLOCATE(X(spp%ppp))
	ALLOCATE(Y(spp%ppp))
	ALLOCATE( rotation_angle(spp%ppp) )
	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta		
	call init_random_seed()
	call RANDOM_NUMBER(r)

	r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)

	Y = r*SIN(theta)
	X = r*COS(theta) + spp%shear_factor*Y

	rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

	X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
 	Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

	spp%vars%X(1,:) = X1*SIN(zeta)
	spp%vars%X(2,:) = X1*COS(zeta)
	spp%vars%X(3,:) = Y1

	DEALLOCATE(X1)
	DEALLOCATE(Y1)
	DEALLOCATE(X)
	DEALLOCATE(Y)
	DEALLOCATE(rotation_angle)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine elliptic_torus


subroutine exponential_torus(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp) :: fl, fr, fm
	REAL(rp) :: rl, rr, rm
	REAL(rp) :: relerr
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, r ! temporary vars
	INTEGER :: pp,jj ! Iterator

	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta
	call init_random_seed()
	call RANDOM_NUMBER(r)

	do pp=1_idef,spp%ppp
		rl = 0.0_rp
		rr = spp%r_outter

		fl = fzero(rl,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))
		fr = fzero(rr,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))
		if (fl.GT.korc_zero) then
			relerr = 100*ABS(fl - fr)/fl
		else
			relerr = 100*ABS(fl - fr)/fr
		end if

		do while(relerr.GT.1.0_rp)
			rm = 0.5_rp*(rr - rl) + rl
			fm = fzero(rm,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))

			if (SIGN(1.0_rp,fm).EQ.SIGN(1.0_rp,fr)) then
				rr = rm
			else
				rl = rm
			end if

			fl = fzero(rl,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))
			fr = fzero(rr,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))

			if (fl.GT.korc_zero) then
				relerr = 100*ABS(fl - fr)/fl
			else
				relerr = 100*ABS(fl - fr)/fr
			end if
		end do
		r(pp) = rm
	end do

	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine exponential_torus


subroutine exponential_elliptic_torus(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp) :: fl, fr, fm
	REAL(rp) :: rl, rr, rm
	REAL(rp) :: relerr
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rotation_angle, theta, zeta, r ! temporary vars
	REAL(rp), DIMENSION(:), ALLOCATABLE :: X,Y,X1,Y1
	LOGICAL :: is_nan
	INTEGER :: pp,jj ! Iterator

	ALLOCATE(X1(spp%ppp))
	ALLOCATE(Y1(spp%ppp))
	ALLOCATE(X(spp%ppp))
	ALLOCATE(Y(spp%ppp))
	ALLOCATE( rotation_angle(spp%ppp) )
	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta
	call init_random_seed()
	call RANDOM_NUMBER(r)

	do pp=1_idef,spp%ppp
		rl = 0.0_rp
		rr = spp%r_outter

		fl = fzero(rl,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))
		fr = fzero(rr,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))
		if (fl.GT.korc_zero) then
			relerr = 100*ABS(fl - fr)/fl
		else
			relerr = 100*ABS(fl - fr)/fr
		end if

		do while(relerr.GT.1.0_rp)
			rm = 0.5_rp*(rr - rl) + rl
			fm = fzero(rm,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))

			if (SIGN(1.0_rp,fm).EQ.SIGN(1.0_rp,fr)) then
				rr = rm
			else
				rl = rm
			end if

			fl = fzero(rl,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))
			fr = fzero(rr,spp%r_outter,spp%vars%g(pp),spp%falloff_rate,r(pp))

			if (fl.GT.korc_zero) then
				relerr = 100*ABS(fl - fr)/fl
			else
				relerr = 100*ABS(fl - fr)/fr
			end if
		end do
		r(pp) = rm
	end do

	Y = r*SIN(theta)
	X = r*COS(theta) + spp%shear_factor*Y

	rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

	X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
 	Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

	spp%vars%X(1,:) = X1*SIN(zeta)
	spp%vars%X(2,:) = X1*COS(zeta)
	spp%vars%X(3,:) = Y1

	DEALLOCATE(X1)
	DEALLOCATE(Y1)
	DEALLOCATE(X)
	DEALLOCATE(Y)
	DEALLOCATE(rotation_angle)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine exponential_elliptic_torus


subroutine gaussian_elliptic_torus(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rotation_angle, theta, zeta, r ! temporary vars
	REAL(rp), DIMENSION(:), ALLOCATABLE :: X,Y,X1,Y1
	REAL(rp) :: sigma
	INTEGER :: pp,jj ! Iterator

	ALLOCATE(X1(spp%ppp))
	ALLOCATE(Y1(spp%ppp))
	ALLOCATE(X(spp%ppp))
	ALLOCATE(Y(spp%ppp))
	ALLOCATE( rotation_angle(spp%ppp) )
	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta
	call init_random_seed()
	call RANDOM_NUMBER(r)

	sigma = 1.0_rp/SQRT(2.0_rp*(spp%falloff_rate/params%cpp%length))
	sigma = sigma/params%cpp%length

	r = sigma*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/sigma**2))*r))
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	Y = r*SIN(theta)
	X = r*COS(theta) + spp%shear_factor*Y

	rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

	X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
 	Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

	spp%vars%X(1,:) = X1*SIN(zeta)
	spp%vars%X(2,:) = X1*COS(zeta)
	spp%vars%X(3,:) = Y1

	DEALLOCATE(X1)
	DEALLOCATE(Y1)
	DEALLOCATE(X)
	DEALLOCATE(Y)
	DEALLOCATE(rotation_angle)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine gaussian_elliptic_torus


subroutine gaussian_torus(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, r ! temporary vars
	REAL(rp) :: sigma
	INTEGER :: jj ! Iterator

	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta		
	call init_random_seed()
	call RANDOM_NUMBER(r)

	sigma = 1.0_rp/SQRT(2.0_rp*(spp%falloff_rate/params%cpp%length))
	sigma = sigma/params%cpp%length

	r = sigma*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/sigma**2))*r))
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine gaussian_torus


subroutine calculate_gaussian_falloff(params,falloff,g)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: falloff
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: g ! temporary vars
	REAL(rp) :: a, b
	
	a = -2.0768_rp
	b = 16.945_rp

	falloff = a + b/g
	falloff = falloff/params%cpp%length
end subroutine calculate_gaussian_falloff


subroutine gaussian_energy(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, r ! temporary vars
	REAL(rp), DIMENSION(:), ALLOCATABLE :: falloff_rate
	INTEGER :: jj ! Iterator

	ALLOCATE( falloff_rate(spp%ppp) )	
	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta		
	call init_random_seed()
	call RANDOM_NUMBER(r)

	call calculate_gaussian_falloff(params,falloff_rate,spp%vars%g)

	r = falloff_rate*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/falloff_rate**2))*r))
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(falloff_rate)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine gaussian_energy


SUBROUTINE calculate_exponential_falloff(g,ko)
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ko
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: g
	REAL(rp) :: a, b
	
	a = 2.0768_rp
	b = -16.945_rp

	ko = a + b/g
END SUBROUTINE calculate_exponential_falloff


FUNCTION fzero(r,a,g,ko,P) RESULT(f)
	REAL(rp) :: f
	REAL(rp), INTENT(IN) :: r	
	REAL(rp), INTENT(IN) :: g
	REAL(rp), INTENT(IN) :: a
	REAL(rp), INTENT(IN) :: ko
	REAL(rp), INTENT(IN) :: P

	f = EXP(-ko*r)*(1.0_rp + r*ko) + (1.0_rp - EXP(-ko*a)*(1.0_rp + a*ko))*P - 1.0_rp
END FUNCTION fzero


subroutine exponential_energy(params,spp)
! * The falloff coefficients are calculated as function of r with units. All these routines need to be modified so they can
! * used as function of an arbitrary plasma radius 'a'.
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, r ! temporary vars
	REAL(rp), DIMENSION(:), ALLOCATABLE :: falloff_rate
	REAL(rp) :: fo, fn
	REAL(rp) :: ro, rn, rt
	REAL(rp) :: relerr
	REAL(rp) :: P
	INTEGER :: pp,jj ! Iterator

	ALLOCATE( falloff_rate(spp%ppp) )	
	ALLOCATE( theta(spp%ppp) )
	ALLOCATE( zeta(spp%ppp) )
	ALLOCATE( r(spp%ppp) )

	! Initial condition of uniformly distributed particles on a disk in the xz-plane
	! A unique velocity direction
	call init_u_random(10986546_8)

	call init_random_seed()
	call RANDOM_NUMBER(theta)
	theta = 2.0_rp*C_PI*theta

	call init_random_seed()
	call RANDOM_NUMBER(zeta)
	zeta = 2.0_rp*C_PI*zeta

	! Uniform distribution on a disk at a fixed azimuthal theta		
	call init_random_seed()
	call RANDOM_NUMBER(r)
	
	call calculate_exponential_falloff(spp%vars%g,falloff_rate)

	do pp=1_idef,spp%ppp
		ro = 0.0_rp
		rn = spp%r_outter
		relerr = ABS(100.0_rp*(ro - rn)/rn)
		do while(relerr.GT.1.0_rp)
			fo = fzero(ro,spp%r_outter,spp%vars%g(pp),falloff_rate(pp),r(pp))
			fn = fzero(rn,spp%r_outter,spp%vars%g(pp),falloff_rate(pp),r(pp))

			rt = rn
			rn = rn - fn*(rn - ro)/(fn - fo)
			ro = rt

			relerr = ABS(100.0_rp*(rn - ro)/ro)
		end do
		r(pp) = rn
	end do

	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(falloff_rate)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine exponential_energy


subroutine intitial_spatial_distribution(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	INTEGER :: ss

	do ss=1_idef,params%num_species
		SELECT CASE (TRIM(spp(ss)%spatial_distribution))
			CASE ('UNIFORM')
				call uniform(spp(ss))
			CASE ('DISK')
				call disk(params,spp(ss))
			CASE ('TORUS')
				call torus(params,spp(ss))
			CASE ('EXPONENTIAL-TORUS')
				call exponential_torus(params,spp(ss))
			CASE ('GAUSSIAN-TORUS')
				call gaussian_torus(params,spp(ss))
			CASE ('ELLIPTIC-TORUS')
				call elliptic_torus(params,spp(ss))
			CASE ('EXPONENTIAL-ELLIPTIC-TORUS')
				call exponential_elliptic_torus(params,spp(ss))
			CASE ('GAUSSIAN-ELLIPTIC-TORUS')
				call gaussian_elliptic_torus(params,spp(ss))
			CASE ('GAUSSIAN-ENERGY')
				call gaussian_energy(params,spp(ss))
			CASE ('EXPONENTIAL-ENERGY')
				call exponential_energy(params,spp(ss))
			CASE DEFAULT
				call torus(params,spp(ss))
		END SELECT
	end do		
end subroutine intitial_spatial_distribution


END MODULE korc_spatial_distribution
