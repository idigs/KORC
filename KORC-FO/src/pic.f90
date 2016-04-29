module pic
use korc_types
use emf
use main_mpi
use omp_lib
implicit none

!REAL(rp), DIMENSION(3), PRIVATE :: B, Y

PRIVATE :: cart_to_cyl, cart_to_tor, interp_field
PUBLIC :: advance_particles_position, advance_particles_velocity

contains

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

    Xtor(1,:) = sqrt( (sqrt(X(1,:)**2 + X(2,:)**2) - Ro)**2 + X(3,:)**2 );
    Xtor(2,:) = atan2(X(3,:), sqrt(X(1,:)**2 + X(2,:)**2) - Ro);
	Xtor(2,:) = modulo(Xtor(2,:),2.0_rp*C_PI)
    Xtor(3,:) = atan2(X(1,:),X(2,:));
	Xtor(3,:) = modulo(Xtor(3,:),2.0_rp*C_PI)

end subroutine cart_to_tor


subroutine interp_field(X,F)
implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: F ! X(1,:) = x, X(2,:) = y, X(3,:) = z

end subroutine interp_field


subroutine interp_analytical_field(prtcls,EB)
implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
!    TYPE(SPECIES), INTENT(INOUT) :: spp
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER ii,pp ! Iterators

	call cart_to_tor(prtcls%X, EB%AB%Ro, prtcls%Y)
!	write(6,'(F15.6,F15.6,F15.6)') prtcls%Y

!$OMP PARALLEL PRIVATE(pp) SHARED(prtcls,EB)
!OMP DO
	do pp=1,SIZE(prtcls%X,2)
		ii = ii + 1
	end do
!OMP END DO
!$OMP END PARALLEL

end subroutine interp_analytical_field

subroutine advance_particles_velocity(params,EB,spp,dt)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), INTENT(IN) :: dt
	INTEGER :: ii,pp ! Iterator(s)


	do ii = 1,params%num_species
		if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
			call interp_analytical_field(spp(ii)%vars, EB)
!			write(6,*) spp(ii)
		else
!			call interp_field(spp(ii), EB)
		end if

	end do
	
end subroutine advance_particles_velocity


subroutine advance_particles_position(params,EB,spp,dt)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), INTENT(IN) :: dt
	INTEGER :: ii,pp ! Iterator(s)

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
