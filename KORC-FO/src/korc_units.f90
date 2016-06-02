module korc_units

    use korc_types
    use constants

    implicit none

	PUBLIC :: compute_charcs_plasma_params, define_time_step, normalize_variables

    contains

subroutine compute_charcs_plasma_params(params,spp,F)
    implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	TYPE(FIELDS), INTENT(IN) :: F
	INTEGER :: ind

	params%cpp%velocity = C_C
	params%cpp%magnetic_field = F%Bo

	! Non-relativistic cyclotron frequency
	spp(:)%wc = ( abs(spp(:)%q)/spp(:)%m )*params%cpp%magnetic_field

	! Relativistic cyclotron frequency
	spp(:)%wc_r =  abs(spp(:)%q)*params%cpp%magnetic_field/( spp(:)%gammao*spp(:)%m )


	ind = maxloc(spp(:)%wc,1) ! Index to maximum cyclotron frequency
	params%cpp%time = 2.0_rp*C_PI/spp(ind)%wc

	ind = maxloc(spp(:)%wc_r,1) ! Index to maximum relativistic cyclotron frequency
	params%cpp%time_r = 2.0_rp*C_PI/spp(ind)%wc_r

	params%cpp%mass = spp(ind)%m
	params%cpp%charge = abs(spp(ind)%q)
	params%cpp%length = params%cpp%velocity*params%cpp%time
	params%cpp%energy = params%cpp%mass*(params%cpp%velocity**2)

	params%cpp%density = 0.0_rp
	params%cpp%electric_field = 0.0_rp
	params%cpp%pressure = 0.0_rp
	params%cpp%temperature = 0.0_rp
end subroutine compute_charcs_plasma_params


subroutine define_time_step(params)
    implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

! 	This definition will be changed as more species and electromagnetic fields
!	are included.

	params%dt = params%dt*params%cpp%time_r

end subroutine define_time_step


subroutine normalize_variables(params,spp,F)
    implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	TYPE(FIELDS), INTENT(INOUT) :: F
	INTEGER :: ii ! Iterator(s)

!	Normalize params variables
	params%dt = params%dt/params%cpp%time

!	Normalize particle variables
	do ii=1,size(spp)
		spp(ii)%q = spp(ii)%q/params%cpp%charge
		spp(ii)%m = spp(ii)%m/params%cpp%mass
		spp(ii)%Eo = spp(ii)%Eo/params%cpp%energy
		spp(ii)%wc = spp(ii)%wc*params%cpp%time
		spp(ii)%vars%X = spp(ii)%vars%X/params%cpp%length
		spp(ii)%vars%V = spp(ii)%vars%V/params%cpp%velocity
		spp(ii)%vars%Rgc = spp(ii)%vars%Rgc/params%cpp%length

		spp(ii)%Ro = spp(ii)%Ro/params%cpp%length
		spp(ii)%Zo = spp(ii)%Zo/params%cpp%length
		spp(ii)%r = spp(ii)%r/params%cpp%length
	end do

!	Normalize electromagnetic fields
	if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		F%Bo = F%Bo/params%cpp%magnetic_field

		F%AB%Bo = F%AB%Bo/params%cpp%magnetic_field
		F%AB%a = F%AB%a/params%cpp%length
		F%AB%Ro = F%AB%Ro/params%cpp%length
		F%AB%lambda = F%AB%lambda/params%cpp%length
		F%AB%Bpo = F%AB%Bpo/params%cpp%magnetic_field
	else
		F%Bo = F%Bo/params%cpp%magnetic_field

		if (ALLOCATED(F%B%R)) F%B%R = F%B%R/params%cpp%magnetic_field
		if (ALLOCATED(F%B%PHI)) F%B%PHI = F%B%PHI/params%cpp%magnetic_field
		if (ALLOCATED(F%B%Z)) F%B%Z = F%B%Z/params%cpp%magnetic_field

		if (ALLOCATED(F%PSIp)) F%PSIp = &
					F%PSIp/(params%cpp%magnetic_field*(params%cpp%length**2))
		if (ALLOCATED(F%PSIp)) F%Ro = F%Ro/params%cpp%length

		F%X%R = F%X%R/params%cpp%length
		! Nothing to do for the PHI component
		F%X%Z = F%X%Z/params%cpp%length
	end if

!	Normalize other variables
		C_Ke = &
		C_Ke/( params%cpp%mass*params%cpp%length*(params%cpp%velocity**2)/(params%cpp%charge**2) )
!		.
!		.
end subroutine normalize_variables

end module korc_units
