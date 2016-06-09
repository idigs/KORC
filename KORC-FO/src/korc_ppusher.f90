module korc_ppusher

    use korc_types
    use constants
    use korc_fields
    use korc_interp
    use korc_hpc

    implicit none

    PRIVATE :: cross
    PUBLIC :: advance_particles_position, advance_particles_velocity

    contains

function cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


subroutine advance_particles_velocity(params,EB,spp,dt,bool)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN) :: bool
	REAL(rp), INTENT(IN) :: dt
	REAL(rp) :: a, gammap, sigma, us, gamma, s ! variables of leapfrog of Vay, J.-L. PoP (2008)
	REAL(rp), DIMENSION(3) :: U, tau, up, t ! variables of leapfrog of Vay, J.-L. PoP (2008)
	REAL(rp), DIMENSION(3) :: U_hs
	REAL(rp), DIMENSION(3) :: vec, b_unit ! variables for diagnostics
	REAL(rp) :: B, vpar, vperp ! variables for diagnostics
	REAL(rp) :: V, kappa, Psyn, gamma_loss
	INTEGER :: ii, pp ! Iterators


	do ii = 1,params%num_species
		if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
			call interp_analytical_field(spp(ii)%vars, EB)
		else
			call interp_field(spp(ii)%vars, EB)
		end if

	a = spp(ii)%q*dt/spp(ii)%m

!$OMP PARALLEL FIRSTPRIVATE(a,dt,bool)&
!$OMP& PRIVATE(pp,U,U_hs,tau,up,gammap,sigma,us,gamma,t,s,&
!$OMP& b_unit,B,vpar,vperp,vec,V,kappa,Psyn,gamma_loss)&
!$OMP& SHARED(ii,spp)
!$OMP DO
		do pp=1,spp(ii)%ppp
			if ( spp(ii)%vars%flag(pp) .EQ. 1_idef ) then
				!! Magnitude of magnetic field
				B = sqrt( sum(spp(ii)%vars%B(:,pp)**2) )

                if (bool) then
				    !! Instantaneous guiding center
				    spp(ii)%vars%Rgc(:,pp) = spp(ii)%vars%X(:,pp)&
				    + gamma*spp(ii)%m*cross(spp(ii)%vars%V(:,pp), spp(ii)%vars%B(:,pp))&
				    /( spp(ii)%q*sum(spp(ii)%vars%B(:,pp)**2) )
                end if

                if (params%radiation_losses) then
				    !! Radiation losses operator     
		            V = sqrt( DOT_PRODUCT(spp(ii)%vars%V(:,pp), spp(ii)%vars%V(:,pp)) )
		            vec =  cross(spp(ii)%vars%V(:,pp), spp(ii)%vars%E(:,pp))&
					    + spp(ii)%vars%V(:,pp)*DOT_PRODUCT(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp))&
					    - spp(ii)%vars%B(:,pp)*V**2
		            kappa = &
				    ABS(spp(ii)%q)*sqrt( DOT_PRODUCT(vec,vec) )/(spp(ii)%vars%gamma(pp)*spp(ii)%m*V**3)
		            
		            !! Synchroton radiated power
				    Psyn = (2.0_rp/3.0_rp)*C_Ke
				    Psyn = Psyn*( (spp(ii)%q*kappa)**2 )
				    Psyn = Psyn*( (spp(ii)%vars%gamma(pp)*V)**4 )

				    spp(ii)%vars%Prad(pp) = Psyn

				    gamma_loss = - dt*Psyn/spp(ii)%m
				    !! Radiation losses operator
                end if

				!! Here we evolve V and gamma in time.
				U = spp(ii)%vars%gamma(pp)*spp(ii)%vars%V(:,pp)
				U_hs = U + &
						0.5_rp*a*( spp(ii)%vars%E(:,pp) + &
						cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)) )
		        
				tau = 0.5_rp*dt*spp(ii)%q*spp(ii)%vars%B(:,pp)/spp(ii)%m
				up = U_hs + 0.5_rp*a*spp(ii)%vars%E(:,pp)
				gammap = sqrt( 1.0_rp + sum(up**2) )
				sigma = gammap**2 - sum(tau**2)
				! us = sum(up*tau) ! variable 'u^*' in Vay, J.-L. PoP (2008)
				us = DOT_PRODUCT(up,tau) ! variable 'u^*' in Vay, J.-L. PoP (2008)
				gamma = sqrt( 0.5_rp*(sigma + sqrt( sigma**2 + 4.0_rp*(sum(tau**2) + us**2) )) )

				!! Radiation losses
				if (params%radiation_losses) then
					gamma = gamma + gamma_loss
				end if

				t = tau/gamma
				s = 1.0_rp/(1.0_rp + sum(t**2)) ! variable 's' in Vay, J.-L. PoP (2008)

		        U = s*( up + sum(up*t)*t + cross(up,t) )
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
				    spp(ii)%vars%mu(pp) = 0.5_rp*spp(ii)%m*(gamma*vperp)**2/B
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
		    spp(ii)%vars%X(:,pp) = spp(ii)%vars%X(:,pp)&
                                    + dt*spp(ii)%vars%V(:,pp)
        end if
	end do
!$OMP END DO
!$OMP END PARALLEL
	end do
end subroutine advance_particles_position

end module korc_ppusher
