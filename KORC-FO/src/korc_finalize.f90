!> @brief Module containing subroutines to terminate parallel communications and free memory.
module korc_finalize
	use korc_types
	use korc_fields
	use korc_hpc

	IMPLICIT NONE

	PUBLIC :: finalize_communications,&
				deallocate_variables

	CONTAINS

!> @brief Interface to function that finalizes MPI communications. See korc_hpc.f90.
!!
!! @param[in] params Core KORC simulation parameters.
subroutine finalize_communications(params)
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	call finalize_mpi(params)
end subroutine finalize_communications


!> @brief Subroutine to free allocatable simulation variables.
!!
!! @param[in,out] params Core KORC simulation parameters.
!! @param[in,out] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation. See korc_types.f90 and
!! korc_fields.f90.
!! @param[in,out] spp An instance of KORC's derived type SPECIES containing all the information of different electron species. See korc_types.f90.
!! @param ii Iterator of the spp array
subroutine deallocate_variables(params,F,spp)
	TYPE(KORC_PARAMS), INTENT(INOUT) 						:: params
	TYPE(FIELDS), INTENT(INOUT) 							:: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	INTEGER 												:: ii

	DEALLOCATE(params%outputs_list)

	do ii=1_idef,params%num_species
		DEALLOCATE(spp(ii)%vars%X)
		DEALLOCATE(spp(ii)%vars%V)
		DEALLOCATE(spp(ii)%vars%Rgc)
		DEALLOCATE(spp(ii)%vars%Y)
		DEALLOCATE(spp(ii)%vars%E)
		DEALLOCATE(spp(ii)%vars%B)
		DEALLOCATE(spp(ii)%vars%ne)
		DEALLOCATE(spp(ii)%vars%Te)
		DEALLOCATE(spp(ii)%vars%Zeff)
		DEALLOCATE(spp(ii)%vars%g)
		DEALLOCATE(spp(ii)%vars%eta)
		DEALLOCATE(spp(ii)%vars%mu)
		DEALLOCATE(spp(ii)%vars%Prad)
		DEALLOCATE(spp(ii)%vars%flag)
		DEALLOCATE(spp(ii)%vars%AUX)
		DEALLOCATE(spp(ii)%vars%wt)
	end do

	DEALLOCATE(spp)

	call DEALLOCATE_FIELDS_ARRAYS(F)
end subroutine deallocate_variables

end module korc_finalize
