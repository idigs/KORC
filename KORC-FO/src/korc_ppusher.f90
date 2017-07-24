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
	REAL(rp) :: g, tmp
	
	g = SQRT(1.0_rp + DOT_PRODUCT(U,U))
	V = U/g

	tmp = spp%q**4/(6.0_rp*C_PI*E0*spp%m**2)

	F2 = tmp*( DOT_PRODUCT(E,V)*E + cross(E,B) + cross(B,cross(B,V)) )
	vec = E + cross(V,B)
	F3 = (tmp*g**2)*( DOT_PRODUCT(E,V)**2 - DOT_PRODUCT(vec,vec) )*V

	Frad = F2 + F3
end subroutine radiation_force


subroutine advance_particles_velocity(params,EB,spp,dt,bool)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN) :: bool
	REAL(rp), INTENT(IN) :: dt
	REAL(rp) :: Prad, B, v, vpar, vperp, tmp ! diagnostics and temporary variables
	REAL(rp) :: a, gp, sigma, us, g, s ! variables of leapfrog of Vay, J.-L. PoP (2008)
	REAL(rp), DIMENSION(3) :: U_L, U_hs, tau, up, t
	REAL(rp), DIMENSION(3) :: U, U_RC, U_os
	REAL(rp), DIMENSION(3) :: Frad, Fcoll
	REAL(rp), DIMENSION(3) :: vec, b_unit ! diagnostics and temporary variables
	INTEGER :: ii, pp ! Iterators

	do ii = 1,params%num_species
		call get_fields(params,spp(ii)%vars, EB)

	    a = spp(ii)%q*dt/spp(ii)%m

!$OMP PARALLEL FIRSTPRIVATE(a,dt,bool)&
!$OMP& PRIVATE(pp,U,U_L,U_hs,tau,up,gp,&
!$OMP& sigma,us,g,t,s,Frad,Fcoll,U_RC,U_os,&
!$OMP& tmp,b_unit,B,vpar,v,vperp,vec,Prad)&
!$OMP& SHARED(ii,spp)

!$OMP SINGLE
	call check_collisions_params(spp(ii))
!$OMP END SINGLE

!$OMP DO
		do pp=1_idef,spp(ii)%ppp
			if ( spp(ii)%vars%flag(pp) .EQ. 1_idef ) then
				U = spp(ii)%vars%g(pp)*spp(ii)%vars%V(:,pp)

				!! Magnitude of magnetic field
				B = SQRT( DOT_PRODUCT(spp(ii)%vars%B(:,pp),spp(ii)%vars%B(:,pp)) )

				U_L = U
				U_RC = U

				! ! ! LEAP-FROG SCHEME FOR LORENTZ FORCE ! ! !
				U_hs = U_L + 0.5_rp*a*( spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)) )		        
				tau = 0.5_rp*dt*spp(ii)%q*spp(ii)%vars%B(:,pp)/spp(ii)%m
				up = U_hs + 0.5_rp*a*spp(ii)%vars%E(:,pp)
				gp = SQRT( 1.0_rp + DOT_PRODUCT(up,up) )
				sigma = gp**2 - DOT_PRODUCT(tau,tau)
				us = DOT_PRODUCT(up,tau) ! variable 'u^*' in Vay, J.-L. PoP (2008)
				g = SQRT( 0.5_rp*(sigma + SQRT(sigma**2 + 4.0_rp*(DOT_PRODUCT(tau,tau) + us**2))) )
				t = tau/g
				s = 1.0_rp/(1.0_rp + DOT_PRODUCT(t,t)) ! variable 's' in Vay, J.-L. PoP (2008)
		        U_L = s*(up + DOT_PRODUCT(up,t)*t + cross(up,t))
				! ! ! LEAP-FROG SCHEME FOR LORENTZ FORCE ! ! !

				! ! ! Splitting operator for including radiation
				U_os = 0.5_rp*(U_L + U)

				if (params%radiation) then
					call radiation_force(spp(ii),U_os,spp(ii)%vars%E(:,pp),spp(ii)%vars%B(:,pp),Frad)
					U_RC = U_RC + a*Frad/spp(ii)%q
				end if
				! ! ! Splitting operator for including radiation

				U = U_L + U_RC - U

				! ! ! Stochastic differential equations for including collisions
				if (params%collisions .AND. (TRIM(params%collisions_model) .EQ. 'SINGLE_SPECIES')) then		
					call include_CoulombCollisions(params,U)
				end if
				! ! ! Stochastic differential equations for including collisions

				if (params%radiation .OR. params%collisions) then
					g = SQRT( 1.0_rp + DOT_PRODUCT(U,U) )
				end if
		        spp(ii)%vars%V(:,pp) = U/g
				spp(ii)%vars%g(pp) = g
		    
                if (bool) then
				    !! Parallel unit vector
				    b_unit = spp(ii)%vars%B(:,pp)/B

					v = SQRT(DOT_PRODUCT(spp(ii)%vars%V(:,pp),spp(ii)%vars%V(:,pp)))
					if (v.GT.korc_zero) then
						!! Parallel and perpendicular components of velocity
						vpar = DOT_PRODUCT(spp(ii)%vars%V(:,pp), b_unit)
						vperp =  DOT_PRODUCT(spp(ii)%vars%V(:,pp),spp(ii)%vars%V(:,pp)) - vpar**2
						if ( vperp .GE. korc_zero ) then
							vperp = SQRT( vperp )
						else
							vperp = 0.0_rp
						end if

						!! Pitch angle
				        spp(ii)%vars%eta(pp) = 180.0_rp*MODULO(ATAN2(vperp,vpar), 2.0_rp*C_PI)/C_PI

						!! Magnetic moment
						spp(ii)%vars%mu(pp) = 0.5_rp*spp(ii)%m*g**2*vperp**2/B
						! See Northrop's book (The adiabatic motion of charged particles)

						!! Radiated power
						tmp = spp(ii)%q**4/(6.0_rp*C_PI*E0*spp(ii)%m**2)

						vec = spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp))

						spp(ii)%vars%Prad(pp) = tmp*( DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%E(:,pp)) + &
						DOT_PRODUCT(cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)),spp(ii)%vars%E(:,pp)) +&
						spp(ii)%vars%g(pp)**2*(DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%V(:,pp))**2 - DOT_PRODUCT(vec,vec)) )

		                spp(ii)%vars%Pin(pp) = spp(ii)%q*DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%V(:,pp))
					else
				        spp(ii)%vars%eta(pp) = 0.0_rp
						spp(ii)%vars%mu(pp) = 0.0_rp
						spp(ii)%vars%Prad(pp) = 0.0_rp
		                spp(ii)%vars%Pin(pp) = 0.0_rp
					end if
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

	if (params%magnetic_field_model .NE. 'UNIFORM') then
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
	end if
end subroutine advance_particles_position

end module korc_ppusher
