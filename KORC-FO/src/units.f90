module units
use korc_types
implicit none

	PUBLIC :: compute_charcs_plasma_params
!	PRIVATE :: 
contains

subroutine compute_charcs_plasma_params(ptcls,EB,cpp)
implicit none
	TYPE(CHARCS_PARAMS), INTENT(OUT) :: cpp
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ptcls
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER :: ii ! Iterator(s)
	INTEGER :: ind
	REAL(rp), DIMENSION(SIZE(ptcls)) :: wc

	! To be defined:
	! REAL(rp) :: density;
	! REAL(rp) :: electric_field;
	! REAL(rp) :: pressure;
	! REAL(rp) :: temperature;

	cpp%velocity = C_C
	cpp%magnetic_field = EB%Bo

	wc = ( abs(ptcls(:)%q)/ptcls(:)%m )*cpp%magnetic_field
	ind = maxloc(wc,1) ! Index to maximum cyclotron frequency

	cpp%time = 2*C_PI/wc(ind)
	cpp%mass = ptcls(ind)%m
	cpp%charge = abs(ptcls(ind)%q)
	cpp%length = cpp%velocity*cpp%time

!	write(6,*) cpp
end subroutine compute_charcs_plasma_params

subroutine define_time_step(cpp,params)
implicit none
	TYPE(KORC_PARAMS), INTENT(OUT) :: params
	TYPE(CHARCS_PARAMS), INTENT(IN) :: cpp

! 	This definition will be changed as more species and electromagnetic fields
!	are included.

	params%dt = params%dt*cpp%time

end subroutine define_time_step


subroutine normalize_variables(params,ptcls,EB,cpp)
implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ptcls
	TYPE(FIELDS), INTENT(INOUT) :: EB
	TYPE(CHARCS_PARAMS), INTENT(IN) :: cpp
	INTEGER :: ii ! Iterator(s)

!	Normalize params variables
	params%dt = params%dt/cpp%time

!	Normalize particle variables
	do ii=1,size(ptcls)
		ptcls(ii)%q = ptcls(ii)%q/cpp%charge
		ptcls(ii)%m = ptcls(ii)%q/cpp%mass
		ptcls(ii)%vars%X = ptcls(ii)%vars%X/cpp%length
		ptcls(ii)%vars%V = ptcls(ii)%vars%V/cpp%velocity
		ptcls(ii)%vars%Rgc = ptcls(ii)%vars%Rgc/cpp%length
	end do

!	Normalize electromagnetic fields
	if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
		EB%AB%Bo = EB%AB%Bo/cpp%magnetic_field
		EB%AB%a = EB%AB%a/cpp%length
		EB%AB%Ro = EB%AB%Ro/cpp%length
		EB%AB%lambda = EB%AB%lambda/cpp%length
		EB%AB%Bpo = EB%AB%Bpo/cpp%magnetic_field

		EB%Bo = EB%Bo/cpp%magnetic_field
	else
		! Normalize data structures EB%E and EB%B
	end if

!	Normalize other variables
!		.
!		.
!		.
end subroutine normalize_variables

end module units
