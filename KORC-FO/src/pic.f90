module pic
use korc_types
use emf
use main_mpi
use omp_lib
implicit none

PRIVATE :: cart_to_cyl, cart_to_tor
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
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Xtor ! Xtor(1,:) = r, Xtor(2,:) = theta, Xtor(3,:) = zeta

    Xtor(1,:) = sqrt( (sqrt(X(1,:)**2 + X(2,:)**2) - Ro)**2 + X(3,:)**2 );
    Xtor(2,:) = atan2(X(3,:), sqrt(X(1,:)**2 + X(2,:)**2) - Ro);
	Xtor(2,:) = modulo(Xtor(2,:),2.0_rp*C_PI)
    Xtor(3,:) = atan2(X(1,:),X(2,:));
	Xtor(3,:) = modulo(Xtor(3,:),2.0_rp*C_PI)

end subroutine cart_to_tor


subroutine advance_particles_velocity(params,EB,ptcls,dt)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ptcls
	REAL(rp), INTENT(IN) :: dt

	
end subroutine advance_particles_velocity


subroutine advance_particles_position(params,EB,ptcls,dt)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ptcls
	REAL(rp), INTENT(IN) :: dt
	INTEGER :: ii,pp ! Iterator(s)


do ii = 1,params%num_species
!$OMP PARALLEL PRIVATE(pp) SHARED(ptcls,dt,params,ii)
!$OMP DO
	do pp = 1,ptcls(ii)%ppp
		write(6,'("Something")')
	end do
!$OMP END DO
!$OMP END PARALLEL
end do


	
end subroutine advance_particles_position


end module pic
