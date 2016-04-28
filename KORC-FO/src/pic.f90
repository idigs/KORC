module pic
use korc_types
use emf
use main_mpi
implicit none

! PRIVATE ::
PUBLIC :: advance_particles

contains

subroutine cart_to_cyl(X)
implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: X

end subroutine cart_to_cyl

subroutine advance_particles(params,EB,ptcls)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ptcls

	
end subroutine advance_particles


end module pic
