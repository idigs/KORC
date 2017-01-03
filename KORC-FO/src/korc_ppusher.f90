module korc_ppusher

    use korc_types
    use korc_constants
    use korc_fields
    use korc_interp
	use korc_collisions
    use korc_hpc

    implicit none

	REAL(rp), PRIVATE :: E0 ! Dimensionless vacuum permittivity

    PRIVATE :: cross,radiation_force,collision_force
    PUBLIC :: initialize_particle_pusher,advance_particles_position,advance_particles_velocity

    contains


subroutine initialize_particle_pusher(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	E0 = C_E0*(params%cpp%mass**2*params%cpp%velocity**3)/(params%cpp%charge**3*params%cpp%Bo)
end subroutine initialize_particle_pusher


function cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


subroutine radiation_force(spp,U,E,B,Frad)
    implicit none
	TYPE(SPECIES), INTENT(IN) :: spp
	REAL(rp), DIMENSION(3), INTENT(IN) :: U, E, B
	REAL(rp), DIMENSION(3), INTENT(OUT) :: Frad
	REAL(rp), DIMENSION(3) :: F1, F2, F3, V, vec
	REAL(rp) :: gamma, tmp
	
	gamma = sqrt(1.0_rp + DOT_PRODUCT(U,U))
	V = U/gamma

	tmp = spp%q**4/(6.0_rp*C_PI*E0*spp%m**2)

	F2 = tmp*( DOT_PRODUCT(E,V)*E + cross(E,B) + cross(B,cross(B,V)) )
	vec = E + cross(V,B)
	F3 = (tmp*gamma**2)*( DOT_PRODUCT(E,V)**2 - DOT_PRODUCT(vec,vec) )*V

	Frad = F2 + F3
end subroutine radiation_force


subroutine advance_particles_velocity(params,EB,spp,dt,bool)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN) :: bool
	REAL(rp), INTENT(IN) :: dt
	REAL(rp) :: Prad, B, vpar, vperp, tmp ! diagnostics and temporary variables
	REAL(rp) :: a, gammap, sigma, us, gamma, s ! variables of leapfrog of Vay, J.-L. PoP (2008)
	REAL(rp), DIMENSION(3) :: U_L, U_hs, tau, up, t
	REAL(rp), DIMENSION(3) :: U, U_RC, U_os
	REAL(rp), DIMENSION(3) :: Frad, Fcoll
	REAL(rp), DIMENSION(3) :: vec, b_unit ! diagnostics and temporary variables
	INTEGER :: ii, pp ! Iterators

	do ii = 1,params%num_species
		if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
			call interp_analytical_field(spp(ii)%vars, EB)
			if (TRIM(EB%electric_field_mode) .EQ. 'PULSE') then 
				spp(ii)%vars%E = EXP(-0.5_rp*((params%time - EB%to)/EB%sig)**2)*spp(ii)%vars%E
			end if
		else
			call interp_field(spp(ii)%vars, EB)
		end if

	    a = spp(ii)%q*dt/spp(ii)%m

!$OMP PARALLEL FIRSTPRIVATE(a,dt,bool)&
!$OMP& PRIVATE(pp,U,U_L,U_hs,tau,up,gammap,&
!$OMP& sigma,us,gamma,t,s,Frad,Fcoll,U_RC,U_os,&
!$OMP& tmp,b_unit,B,vpar,vperp,vec,Prad)&
!$OMP& SHARED(ii,spp)
!$OMP DO
		do pp=1,spp(ii)%ppp
			if ( spp(ii)%vars%flag(pp) .EQ. 1_idef ) then
				U = spp(ii)%vars%gamma(pp)*spp(ii)%vars%V(:,pp)

				!! Magnitude of magnetic field
				B = sqrt( DOT_PRODUCT(spp(ii)%vars%B(:,pp),spp(ii)%vars%B(:,pp)) )

                if (bool) then
				    !! Instantaneous guiding center
				    spp(ii)%vars%Rgc(:,pp) = spp(ii)%vars%X(:,pp)&
				    + gamma*spp(ii)%m*cross(spp(ii)%vars%V(:,pp), spp(ii)%vars%B(:,pp))/( spp(ii)%q*B**2 )
                end if

				U_L = U
				U_RC = U

				! ! ! LEAP-FROG SCHEME FOR LORENTZ FORCE ! ! !
				U_hs = U_L + 0.5_rp*a*( spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)) )		        
				tau = 0.5_rp*dt*spp(ii)%q*spp(ii)%vars%B(:,pp)/spp(ii)%m
				up = U_hs + 0.5_rp*a*spp(ii)%vars%E(:,pp)
				gammap = sqrt( 1.0_rp + DOT_PRODUCT(up,up) )
				sigma = gammap**2 - DOT_PRODUCT(tau,tau)
				us = DOT_PRODUCT(up,tau) ! variable 'u^*' in Vay, J.-L. PoP (2008)
				gamma = sqrt( 0.5_rp*(sigma + sqrt(sigma**2 + 4.0_rp*(DOT_PRODUCT(tau,tau) + us**2))) )
				t = tau/gamma
				s = 1.0_rp/(1.0_rp + DOT_PRODUCT(t,t)) ! variable 's' in Vay, J.-L. PoP (2008)
		        U_L = s*(up + DOT_PRODUCT(up,t)*t + cross(up,t))
				! ! ! LEAP-FROG SCHEME FOR LORENTZ FORCE ! ! !

				! ! ! Splitting operator for including radiation and collisions
				U_os = 0.5_rp*(U_L + U)

				if (params%radiation) then
					call radiation_force(spp(ii),U_os,spp(ii)%vars%E(:,pp),spp(ii)%vars%B(:,pp),Frad)
					U_RC = U_RC + a*Frad/spp(ii)%q
				end if

!				if (params%collisions .AND. &
!				(TRIM(params%collisions_model) .EQ. 'MULTIPLE_SPECIES')) then
!					call collision_force(spp(ii),U_os,Fcoll)
!					U_RC = U_RC + a*Fcoll/spp(ii)%q
!				end if

				U = U_L + U_RC - U

				if (params%collisions) then
					call include_collisions(params,spp(ii)%vars%B(:,pp),U)
				end if
				! ! ! Splitting operator for including radiation and collisions

				if (params%radiation .OR. params%collisions) then
					gamma = SQRT( 1.0_rp + DOT_PRODUCT(U,U) )
				end if
		        spp(ii)%vars%V(:,pp) = U/gamma
				spp(ii)%vars%gamma(pp) = gamma
		    
                if (bool) then
				    !! Parallel unit vector
				    b_unit = spp(ii)%vars%B(:,pp)/B

				    !! Parallel and perpendicular components of velocity
				    vpar = DOT_PRODUCT(spp(ii)%vars%V(:,pp), b_unit)
				    vperp =  DOT_PRODUCT(spp(ii)%vars%V(:,pp),spp(ii)%vars%V(:,pp)) - vpar**2
				    if ( vperp .GE. korc_zero ) then
					    vperp = sqrt( vperp )
				    else
					    vperp = 0.0_rp
				    end if

				    !! Pitch angle
		            spp(ii)%vars%eta(pp) = 180.0_rp*modulo(atan2(vperp,vpar), 2.0_rp*C_PI)/C_PI

				    !! Magnetic moment
				    spp(ii)%vars%mu(pp) = 0.5_rp*spp(ii)%m*gamma*vperp**2/B

					!! Radiated power
					tmp = spp(ii)%q**4/(6.0_rp*C_PI*E0*spp(ii)%m**2)

					vec = spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp))

					spp(ii)%vars%Prad(pp) = tmp*( DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%E(:,pp)) + &
					DOT_PRODUCT(cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)),spp(ii)%vars%E(:,pp)) +&
					spp(ii)%vars%gamma(pp)**2*(DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%V(:,pp))**2 - DOT_PRODUCT(vec,vec)) )

                    spp(ii)%vars%Pin(pp) = spp(ii)%q*DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%V(:,pp))
                end if
			end if
		end do
!$OMP END DO
!$OMP END PARALLEL

	end do
end subroutine advance_particles_velocity


subroutine advance_particles_position(params,EB,spp,dt)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), INTENT(IN) :: dt
	INTEGER :: ii, pp ! Iterators

    do ii = 1,params%num_species
!$OMP PARALLEL PRIVATE(pp) SHARED(ii,spp,dt,params)
!$OMP DO
	do pp = 1,spp(ii)%ppp
        if ( spp(ii)%vars%flag(pp) .EQ. 1_idef ) then
		    spp(ii)%vars%X(:,pp) = spp(ii)%vars%X(:,pp) + dt*spp(ii)%vars%V(:,pp)
        end if
	end do
!$OMP END DO
!$OMP END PARALLEL
	end do
end subroutine advance_particles_position

end module korc_ppusher
