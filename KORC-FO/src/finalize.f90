module finalize

	use korc_types
	use korc_hpc

	implicit none

	contains

subroutine finalize_communications(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	call finalize_mpi(params)
end subroutine finalize_communications


subroutine deallocate_variables(params,EB,spp)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(INOUT) :: EB
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	INTEGER :: ii ! Iterator

	do ii=1,params%num_species
		DEALLOCATE(spp(ii)%vars%X)
		DEALLOCATE(spp(ii)%vars%V)
		DEALLOCATE(spp(ii)%vars%Rgc)
		DEALLOCATE(spp(ii)%vars%Y)
		DEALLOCATE(spp(ii)%vars%E)
		DEALLOCATE(spp(ii)%vars%B)
		DEALLOCATE(spp(ii)%vars%gamma)
		DEALLOCATE(spp(ii)%vars%eta)
		DEALLOCATE(spp(ii)%vars%mu)
		DEALLOCATE(spp(ii)%vars%kappa)
		DEALLOCATE(spp(ii)%vars%tau)
	end do

	DEALLOCATE(spp)

    if (params%magnetic_field_model .EQ. 'EXTERNAL') then
        call DEALLOCATE_FIELDS(EB)
    end if
end subroutine deallocate_variables

end module finalize
