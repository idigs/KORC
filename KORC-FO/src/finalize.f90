module finalize
use korc_types
use main_mpi
implicit none

contains

subroutine finalize_communications(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	call finalize_mpi(params)
end subroutine finalize_communications


subroutine deallocate_variables(params,ptcls)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ptcls
	INTEGER(ip) :: ii ! Iterator

	do ii=1,params%num_species
		DEALLOCATE(ptcls(ii)%vars%X)
		DEALLOCATE(ptcls(ii)%vars%V)
		DEALLOCATE(ptcls(ii)%vars%Rgc)
		DEALLOCATE(ptcls(ii)%vars%Y)
		DEALLOCATE(ptcls(ii)%vars%E)
		DEALLOCATE(ptcls(ii)%vars%B)
		DEALLOCATE(ptcls(ii)%vars%gamma)
		DEALLOCATE(ptcls(ii)%vars%eta)
		DEALLOCATE(ptcls(ii)%vars%mu)
		DEALLOCATE(ptcls(ii)%vars%kappa)
		DEALLOCATE(ptcls(ii)%vars%tau)
	end do

	DEALLOCATE(ptcls)
end subroutine deallocate_variables

end module finalize
