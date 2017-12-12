module korc_profiles
    use korc_types
	use korc_hpc
	use korc_coords
	use korc_interp
	use korc_HDF5

    implicit none

	PUBLIC :: get_profiles,load_profiles_data_from_hdf5
	PRIVATE :: get_analytical_profiles,uniform_profiles

    contains

subroutine get_analytical_profiles(P,Y,ne,Te,flag)
	TYPE(PROFILES), INTENT(IN) :: P
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ne ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Te ! E(1,:) = Ex, E(2,:) = Ey, E(3,:) = Ez
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
	REAL(rp) :: fr
	INTEGER(ip) pp ! Iterator(s)
	INTEGER(ip) :: ss

	ss = SIZE(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,fr) SHARED(P,Y,ne,Te,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
			fr = 1_ip - TANH(2.0_rp*Y(1,pp)/P%a)**P%n_ne
			ne(pp) = P%neo*fr

			fr = 1_ip - TANH(2.0_rp*Y(1,pp)/P%a)**P%n_Te
			Te(pp) = P%Teo*fr
        end if
	end do
!$OMP END PARALLEL DO
end subroutine get_analytical_profiles


subroutine uniform_profiles(vars,P)
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(PROFILES), INTENT(IN) :: P

	vars%ne = P%neo
	vars%Te = P%Teo
end subroutine uniform_profiles


subroutine get_profiles(params,vars,P)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(PROFILES), INTENT(IN) :: P

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			call get_analytical_profiles(P,vars%Y,vars%ne,vars%Te,vars%flag)
		CASE('EXTERNAL')
!			call interp_profiles(vars,P)
		CASE('UNIFORM')
			call uniform_profiles(vars,P)
		CASE DEFAULT
	END SELECT
end subroutine get_profiles


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! Subroutines for getting the profiles data from HDF5 files
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
subroutine load_profiles_data_from_hdf5(params,P)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PROFILES), INTENT(INOUT) :: P
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
    REAL(rp), DIMENSION(:), ALLOCATABLE :: A
	INTEGER :: h5error
    INTEGER ir, iphi, iz

end subroutine load_profiles_data_from_hdf5

end module korc_profiles
