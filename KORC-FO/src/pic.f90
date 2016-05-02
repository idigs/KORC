module pic
use korc_types
use emf
use main_mpi
use omp_lib
implicit none

!REAL(rp), DIMENSION(3), PRIVATE :: B, Y

PRIVATE :: cart_to_cyl, cart_to_tor, interp_field, cross
PUBLIC :: advance_particles_position, advance_particles_velocity

contains

function cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a,b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


subroutine cart_to_cyl(X,Xcyl)
implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Xcyl ! Xcyl(1,:) = R, Xcyl(2,:) = phi, Xcyl(3,:) = Z

	Xcyl(1,:) = sqrt(X(1,:)**2 + X(2,:)**2)
    Xcyl(2,:) = atan2(X(2,:), X(1,:))
    Xcyl(2,:) = modulo(Xcyl(2,:), 2.0_rp*C_PI)
    Xcyl(3,:) = X(3,:)

end subroutine cart_to_cyl


subroutine cart_to_tor(X,Ro,Xtor)
implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), INTENT(IN) :: Ro
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Xtor ! Xtor(1,:) = r, Xtor(2,:) = theta, Xtor(3,:) = zeta
	INTEGER(ip) :: pp ! Iterator(s)
	INTEGER(ip) :: ss

	ss = SIZE(X,2)

!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(X,Xtor)
!$OMP DO
	do pp=1,ss
		Xtor(1,pp) = sqrt( (sqrt(X(1,pp)**2 + X(2,pp)**2) - Ro)**2 + X(3,pp)**2 )
		Xtor(2,pp) = atan2(X(3,pp), sqrt(X(1,pp)**2 + X(2,pp)**2) - Ro)
		Xtor(2,pp) = modulo(Xtor(2,pp),2.0_rp*C_PI)
		Xtor(3,pp) = atan2(X(1,pp),X(2,pp))
		Xtor(3,pp) = modulo(Xtor(3,pp),2.0_rp*C_PI)
	end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine cart_to_tor


subroutine interp_field(prtcls,EB)
implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER(ip) ii,pp ! Iterator(s)
	INTEGER(ip) :: ss

end subroutine interp_field


subroutine interp_analytical_field(prtcls,EB)
implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: EB

	call cart_to_tor(prtcls%X, EB%AB%Ro, prtcls%Y)

	call analytical_magnetic_field(EB,prtcls%Y,prtcls%B)
end subroutine interp_analytical_field


subroutine advance_particles_velocity(params,EB,spp,dt)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), INTENT(IN) :: dt
	REAL(rp) :: a, gammap, sigma, us, gamma, s
	REAL(rp), DIMENSION(3) :: U, U_half_step, tau, up, t
	INTEGER(ip) :: ii,pp ! Iterator(s)


	do ii = 1,params%num_species
		if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
			call interp_analytical_field(spp(ii)%vars, EB)
!			write(6,'(F15.6,F15.6,F15.6)') spp(ii)%vars%B
		else
			call interp_field(spp(ii)%vars, EB)
		end if

	a = spp(ii)%q*dt/spp(ii)%m

!$OMP PARALLEL FIRSTPRIVATE(a,dt) PRIVATE(pp,U,U_half_step,tau,up,gammap,sigma,us,gamma,t,s) SHARED(spp)
!$OMP DO
		do pp=1,spp(ii)%ppp
			U = spp(ii)%vars%gamma(pp)*spp(ii)%vars%V(:,pp)
			U_half_step = U + &
					0.5_rp*a*( spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)) )
            
			tau = 0.5_rp*dt*spp(ii)%q*spp(ii)%vars%B(:,pp)/spp(ii)%m
			up = U_half_step + 0.5_rp*a*spp(ii)%vars%E(:,pp)
			gammap = sqrt( 1.0_rp + sum(up**2) )
			sigma = gammap**2 - sum(tau**2)
			us = sum(up*tau) ! variable 'u^*' in Vay, J.-L. PoP (2008)
			gamma = sqrt(0.5_rp)*sqrt( sigma + sqrt( sigma**2 + 4.0_rp*(sum(tau**2) + us**2) ) )
			t = tau/gamma
			s = 1.0_rp/(1.0_rp + sum(t**2)) ! variable 's' in Vay, J.-L. PoP (2008)

!			write(6,'("Debuggin list:")') 
!			write(6,*) dt
!			write(6,*) spp(ii)%q
!			write(6,*) spp(ii)%vars%B(:,pp)
!			write(6,*) spp(ii)%m
			
!			write(6,*) U
!			write(6,*) U_half_step
!			write(6,*) tau
!			write(6,*) up
!			write(6,*) gammap
!			write(6,*) sigma
!			write(6,*) us
!			write(6,*) gamma
!			write(6,*) t
!			write(6,*) s
            
            U = s*( up + sum(up*t)*t + cross(up,t) )
            spp(ii)%vars%V(:,pp) = U/sqrt(1 + sum(U**2));

			spp(ii)%vars%gamma(pp) = gamma

			! Instantaneous guiding center
!			spp(ii)%vars%Rgc(:,pp) = spp(ii)%vars%X(:,pp) + &
!			spp(ii)%vars%gamma(pp)*spp(ii)%m*cross(U_half_step,spp(ii)%vars%B(:,pp))/(spp(ii)%q*sum(spp(ii)%vars%B(:,pp)**2))

			! Curvature and torsion
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
	INTEGER(ip) :: ii,pp ! Iterator(s)

    do ii = 1,params%num_species
!$OMP PARALLEL PRIVATE(pp) SHARED(spp,dt,params,ii)
!$OMP DO
	do pp = 1,spp(ii)%ppp
		spp(ii)%vars%X(:,pp) = spp(ii)%vars%X(:,pp) + dt*spp(ii)%vars%V(:,pp)
	end do
!$OMP END DO
!$OMP END PARALLEL
	end do
end subroutine advance_particles_position


end module pic
