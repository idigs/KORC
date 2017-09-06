module korc_units
    use korc_types
    use korc_constants

    implicit none

	PUBLIC :: compute_charcs_plasma_params,normalize_variables

    contains

subroutine compute_charcs_plasma_params(params,spp,F)
    implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	TYPE(FIELDS), INTENT(IN) :: F
	INTEGER :: ind

	params%cpp%velocity = C_C
	params%cpp%Bo = ABS(F%Bo)
	params%cpp%Eo = ABS(params%cpp%velocity*params%cpp%Bo)

	! Non-relativistic cyclotron frequency
	spp(:)%wc = ( ABS(spp(:)%q)/spp(:)%m )*params%cpp%Bo

	! Relativistic cyclotron frequency
	spp(:)%wc_r =  ABS(spp(:)%q)*params%cpp%Bo/( spp(:)%go*spp(:)%m )


	ind = MAXLOC(spp(:)%wc,1) ! Index to maximum cyclotron frequency
	params%cpp%time = 1.0_rp/spp(ind)%wc

	ind = MAXLOC(spp(:)%wc_r,1) ! Index to maximum relativistic cyclotron frequency
	params%cpp%time_r = 1.0_rp/spp(ind)%wc_r

	params%cpp%mass = spp(ind)%m
	params%cpp%charge = ABS(spp(ind)%q)
	params%cpp%length = params%cpp%velocity*params%cpp%time
	params%cpp%energy = params%cpp%mass*params%cpp%velocity**2

	params%cpp%density = 1.0_rp/params%cpp%length**3
	params%cpp%pressure = 0.0_rp
	params%cpp%temperature = params%cpp%energy
end subroutine compute_charcs_plasma_params


subroutine normalize_variables(params,spp,F,P)
    implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	TYPE(FIELDS), INTENT(INOUT) :: F
	TYPE(PROFILES), INTENT(INOUT) :: P
	INTEGER :: ii ! Iterator(s)

!	Normalize params variables
	params%dt = params%dt/params%cpp%time
	params%simulation_time = params%simulation_time/params%cpp%time
	params%snapshot_frequency = params%snapshot_frequency/params%cpp%time
	params%minimum_particle_energy = params%minimum_particle_energy/params%cpp%energy

!	Normalize particle variables
	do ii=1_idef,size(spp)
		spp(ii)%q = spp(ii)%q/params%cpp%charge
		spp(ii)%m = spp(ii)%m/params%cpp%mass
		spp(ii)%Eo = spp(ii)%Eo/params%cpp%energy
		spp(ii)%Eo_lims = spp(ii)%Eo_lims/params%cpp%energy
		spp(ii)%wc = spp(ii)%wc*params%cpp%time
		spp(ii)%wc_r = spp(ii)%wc_r*params%cpp%time
		spp(ii)%vars%X = spp(ii)%vars%X/params%cpp%length
		spp(ii)%vars%V = spp(ii)%vars%V/params%cpp%velocity
		spp(ii)%vars%Rgc = spp(ii)%vars%Rgc/params%cpp%length

		spp(ii)%Ro = spp(ii)%Ro/params%cpp%length
		spp(ii)%Zo = spp(ii)%Zo/params%cpp%length
		spp(ii)%r = spp(ii)%r/params%cpp%length
	end do

!	Normalize electromagnetic fields and plasma profiles
	if (params%plasma_model .EQ. 'ANALYTICAL') then
		F%Bo = F%Bo/params%cpp%Bo

		F%AB%Bo = F%AB%Bo/params%cpp%Bo
		F%AB%a = F%AB%a/params%cpp%length
		F%AB%Ro = F%AB%Ro/params%cpp%length
		F%Ro = F%Ro/params%cpp%length
		F%AB%lambda = F%AB%lambda/params%cpp%length
		F%AB%Bpo = F%AB%Bpo/params%cpp%Bo

        ! Electric field parameters
		F%Eo = F%Eo/params%cpp%Eo
        F%to = F%to/params%cpp%time
        F%sig = F%sig/params%cpp%time

		P%a = P%a/params%cpp%length
		P%neo = P%neo/params%cpp%density
		P%Teo = P%Teo/params%cpp%temperature
	else if (params%plasma_model .EQ. 'EXTERNAL') then
		F%Bo = F%Bo/params%cpp%Bo
		F%Ro = F%Ro/params%cpp%length
		F%Zo = F%Zo/params%cpp%length

		if (ALLOCATED(F%B_3D%R)) F%B_3D%R = F%B_3D%R/params%cpp%Bo
		if (ALLOCATED(F%B_3D%PHI)) F%B_3D%PHI = F%B_3D%PHI/params%cpp%Bo
		if (ALLOCATED(F%B_3D%Z)) F%B_3D%Z = F%B_3D%Z/params%cpp%Bo

		if (ALLOCATED(F%PSIp)) F%PSIp = F%PSIp/(params%cpp%Bo*params%cpp%length**2)

		if (ALLOCATED(F%B_2D%R)) F%B_2D%R = F%B_2D%R/params%cpp%Bo
		if (ALLOCATED(F%B_2D%PHI)) F%B_2D%PHI = F%B_2D%PHI/params%cpp%Bo
		if (ALLOCATED(F%B_2D%Z)) F%B_2D%Z = F%B_2D%Z/params%cpp%Bo

		F%X%R = F%X%R/params%cpp%length
		! Nothing to do for the PHI component
		F%X%Z = F%X%Z/params%cpp%length
	else if (params%plasma_model .EQ. 'UNIFORM') then
		F%Bo = F%Bo/params%cpp%Bo
		F%Eo = F%Eo/params%cpp%Eo
	end if	
end subroutine normalize_variables

end module korc_units
