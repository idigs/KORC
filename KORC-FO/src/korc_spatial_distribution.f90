MODULE korc_spatial_distribution
	USE korc_types
	USE korc_constants
	USE korc_HDF5
	USE korc_hpc
    use korc_fields
    use korc_rnd_numbers
	use hammersley_generator

	IMPLICIT NONE

	PUBLIC :: intitial_spatial_distribution
	PRIVATE :: uniform_distribution,&
				disk_distribution,&
				torus_distribution,&
				gaussian_distribution,&
				gaussian_energy_distribution,&
				calculate_falloff,&
				exponential_energy_distribution,&
				ko,&
				fzero

	CONTAINS


subroutine uniform_distribution(spp)
	TYPE(SPECIES), INTENT(INOUT) :: spp

	spp%vars%X = 0.0_rp
end subroutine uniform_distribution


subroutine disk_distribution(params,spp)
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
end subroutine disk_distribution


subroutine torus_distribution(params,spp)
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
end subroutine torus_distribution


subroutine gaussian_distribution(params,spp)
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

	r = spp%falloff_rate*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/spp%falloff_rate**2))*r))
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine gaussian_distribution


subroutine calculate_falloff(params,falloff,g)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: falloff
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: g ! temporary vars
	REAL(rp) :: a, b
	
	a = -2.0768_rp
	b = 16.945_rp

	falloff = a + b/g
	falloff = falloff/params%cpp%length
end subroutine calculate_falloff


subroutine gaussian_energy_distribution(params,spp)
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

	call calculate_falloff(params,falloff_rate,spp%vars%g)

	r = falloff_rate*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/falloff_rate**2))*r))
	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(falloff_rate)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine gaussian_energy_distribution


FUNCTION ko(g)
	REAL(rp) :: ko
	REAL(rp), INTENT(IN) :: g
	REAL(rp) :: a, b
	
	a = 2.0768_rp
	b = -16.945_rp

	ko = a + b/g
END FUNCTION ko

FUNCTION fzero(params,r,a,g,P) RESULT(f)
	REAL(rp) :: f
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), INTENT(IN) :: r	
	REAL(rp), INTENT(IN) :: g
	REAL(rp), INTENT(IN) :: a
	REAL(rp), INTENT(IN) :: P
	REAL(rp) :: ru,au

	ru = r*params%cpp%length
	au = a*params%cpp%length

	f = EXP(-ko(g)*ru)*(1.0_rp + ru*ko(g)) + (1.0_rp - EXP(-ko(g)*au)*(1.0_rp + au*ko(g)))*P - 1.0_rp
END FUNCTION fzero


subroutine exponential_energy_distribution(params,spp)
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

!!$OMP PARALLEL DO SHARED(params,r,spp) PRIVATE(pp,ro,rn,rt,fo,fn,relerr)
	do pp=1_idef,spp%ppp
		ro = 0.0_rp
		rn = spp%r_outter
		relerr = ABS(100.0_rp*(ro - rn)/rn)
		do while(relerr.GT.1.0_rp)
			fo = fzero(params,ro,spp%r_outter,spp%vars%g(pp),r(pp))
			fn = fzero(params,rn,spp%r_outter,spp%vars%g(pp),r(pp))

			rt = rn
			rn = rn - fn*(rn - ro)/(fn - fo)
			ro = rt

			relerr = ABS(100.0_rp*(rn - ro)/ro)
		end do
		r(pp) = rn
	end do
!!$OMP END PARALLEL DO

	spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
	spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
	spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

	DEALLOCATE(falloff_rate)
	DEALLOCATE(theta)
	DEALLOCATE(zeta)
	DEALLOCATE(r)
end subroutine exponential_energy_distribution


subroutine intitial_spatial_distribution(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	INTEGER :: ss

	do ss=1_idef,params%num_species
		SELECT CASE (TRIM(spp(ss)%spatial_distribution))
			CASE ('UNIFORM')
				call uniform_distribution(spp(ss))
			CASE ('DISK')
				call disk_distribution(params,spp(ss))
			CASE ('TORUS')
				call torus_distribution(params,spp(ss))
			CASE ('GAUSSIAN')
				call gaussian_distribution(params,spp(ss))
			CASE ('GAUSSIAN-ENERGY')
				call gaussian_energy_distribution(params,spp(ss))
			CASE ('EXPONENTIAL-ENERGY')
				call exponential_energy_distribution(params,spp(ss))
			CASE DEFAULT
				call torus_distribution(params,spp(ss))
		END SELECT
	end do		
end subroutine intitial_spatial_distribution


END MODULE korc_spatial_distribution
