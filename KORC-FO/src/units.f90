module units
use korc_types
implicit none

	PUBLIC :: compute_charcs_plasma_params
!	PRIVATE :: 
contains

subroutine compute_charcs_plasma_params(ptcls,EB,cp)
implicit none
	TYPE(CHARCS_PARAMS), INTENT(OUT) :: cp
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

	cp%velocity = C_C
	cp%magnetic_field = EB%Bo

	wc = ( abs(ptcls(:)%q)/ptcls(:)%m )*cp%magnetic_field
	ind = maxloc(wc,1) ! Index to maximum cyclotron frequency

	cp%time = 2*C_PI/wc(ind)
	cp%mass = ptcls(ind)%m
	cp%charge = abs(ptcls(ind)%q)
	cp%length = cp%velocity*cp%time

!	write(6,*) cp
end subroutine compute_charcs_plasma_params

end module units
