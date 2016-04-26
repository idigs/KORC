module units
use korc_types
implicit none

contains

subroutine compute_charcs_plasma_params(ptcls,cp)
implicit none
	TYPE(CHARCS_PARAMS) :: cp
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: ptcls
	INTEGER :: ii ! Iterator(s)

	do ii=1,SIZE(ptcls)

	end do

	
end subroutine compute_charcs_plasma_params

end module units
