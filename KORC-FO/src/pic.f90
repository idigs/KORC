module pic
use korc_types
use emf
use main_mpi
implicit none

contains

subroutine advance_particle_position(params,EB,ptcls)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ptcls

end subroutine advance_particle_position


subroutine advance_particle_velocity(params,EB,ptcls)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ptcls

end subroutine advance_particle_velocity


end module pic
