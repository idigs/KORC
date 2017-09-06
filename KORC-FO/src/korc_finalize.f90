module korc_finalize
	use korc_types
	use korc_hpc

	implicit none

	contains

subroutine finalize_communications(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	
	call finalize_mpi(params)
end subroutine finalize_communications


subroutine deallocate_variables(params,F,spp)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	TYPE(FIELDS), INTENT(INOUT) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	INTEGER :: ii ! Iterator

	DEALLOCATE(params%outputs_list)

	do ii=1,params%num_species
		DEALLOCATE(spp(ii)%vars%X)
		DEALLOCATE(spp(ii)%vars%V)
		DEALLOCATE(spp(ii)%vars%Rgc)
		DEALLOCATE(spp(ii)%vars%Y)
		DEALLOCATE(spp(ii)%vars%E)
		DEALLOCATE(spp(ii)%vars%B)
		DEALLOCATE(spp(ii)%vars%ne)
		DEALLOCATE(spp(ii)%vars%Te)
		DEALLOCATE(spp(ii)%vars%g)
		DEALLOCATE(spp(ii)%vars%eta)
		DEALLOCATE(spp(ii)%vars%mu)
		DEALLOCATE(spp(ii)%vars%Prad)
		DEALLOCATE(spp(ii)%vars%flag)
		DEALLOCATE(spp(ii)%vars%AUX)
	end do

	DEALLOCATE(spp)

	call DEALLOCATE_FIELDS_ARRAYS(F)
end subroutine deallocate_variables

end module korc_finalize
