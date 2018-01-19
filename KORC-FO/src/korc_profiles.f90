module korc_profiles
    use korc_types
	use korc_hpc
	use korc_coords
	use korc_interp
	use korc_HDF5

    implicit none

	PUBLIC :: get_profiles,&
				initialize_profiles
	PRIVATE :: get_analytical_profiles,&
				uniform_profiles,&
				load_profiles_data_from_hdf5,&
				ALLOCATE_2D_PROFILES_ARRAYS,&
				ALLOCATE_3D_PROFILES_ARRAYS

    contains

subroutine initialize_profiles(params,P)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PROFILES), INTENT(OUT) :: P
    CHARACTER(MAX_STRING_LENGTH) :: ne_profile
    CHARACTER(MAX_STRING_LENGTH) :: Te_profile
    CHARACTER(MAX_STRING_LENGTH) :: Zeff_profile
    CHARACTER(MAX_STRING_LENGTH) :: filename
	REAL(rp) :: radius_profile
	REAL(rp) :: neo
	REAL(rp) :: Teo
	REAL(rp) :: Zeffo
	REAL(rp) :: n_ne
	REAL(rp) :: n_Te
	REAL(rp) :: n_Zeff
	REAL(rp), DIMENSION(4) :: a_ne
	REAL(rp), DIMENSION(4) :: a_Te
	REAL(rp), DIMENSION(4) :: a_Zeff
	LOGICAL :: Efield, Bfield, Bflux, axisymmetric

	NAMELIST /plasmaProfiles/ radius_profile,ne_profile,neo,n_ne,a_ne,&
											Te_profile,Teo,n_Te,a_Te,&
											Zeff_profile,Zeffo,n_Zeff,a_Zeff,filename

	NAMELIST /externalPlasmaModel/ Efield, Bfield, Bflux, axisymmetric

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=plasmaProfiles)
			close(default_unit_open)

			P%a = radius_profile
			P%ne_profile = TRIM(ne_profile)
			P%neo = neo
			P%n_ne = n_ne
			P%a_ne = a_ne

			P%Te_profile = TRIM(Te_profile)
			P%Teo = Teo*C_E ! Converted to Joules
			P%n_Te = n_Te
			P%a_Te = a_Te

			P%Zeff_profile = TRIM(Zeff_profile)
			P%Zeffo = Zeffo
			P%n_Zeff = n_Zeff
			P%a_Zeff = a_Zeff
		CASE('EXTERNAL')
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=plasmaProfiles)
			close(default_unit_open)

			P%ne_profile = TRIM(ne_profile)
			P%neo = neo
			P%Te_profile = TRIM(Te_profile)
			P%Teo = Teo*C_E ! Converted to Joules
			P%Zeff_profile = TRIM(Zeff_profile)
			P%Zeffo = Zeffo

			P%filename = TRIM(filename)

			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=externalPlasmaModel)
			close(default_unit_open)

			P%axisymmetric = axisymmetric
			
			call load_profiles_data_from_hdf5(params,P)
		CASE('UNIFORM')
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=plasmaProfiles)
			close(default_unit_open)

			P%a = radius_profile
			P%ne_profile = TRIM(ne_profile)
			P%neo = neo
			P%n_ne = 0.0_rp
			P%a_ne = (/0.0_rp,0.0_rp,0.0_rp,0.0_rp/)

			P%Te_profile = TRIM(Te_profile)
			P%Teo = Teo*C_E ! Converted to Joules
			P%n_Te = 0.0_rp
			P%a_Te = (/0.0_rp,0.0_rp,0.0_rp,0.0_rp/)

			P%Zeff_profile = TRIM(Zeff_profile)
			P%Zeffo = Zeffo
			P%n_Zeff = 0.0_rp
			P%a_Zeff = (/0.0_rp,0.0_rp,0.0_rp,0.0_rp/)
		CASE DEFAULT
	END SELECT
end subroutine initialize_profiles


subroutine uniform_profiles(vars,P)
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(PROFILES), INTENT(IN) :: P

	vars%ne = P%neo
	vars%Te = P%Teo
	vars%Zeff = P%Zeffo
end subroutine uniform_profiles


subroutine get_analytical_profiles(P,Y,ne,Te,Zeff,flag)
	TYPE(PROFILES), INTENT(IN) :: P
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ne ! 
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Te ! 
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Zeff ! 
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
	REAL(rp) :: r_a
	REAL(rp) :: fr
	INTEGER(ip) pp ! Iterator(s)
	INTEGER(ip) :: ss

	ss = SIZE(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,fr,r_a) SHARED(P,Y,ne,Te,Zeff,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then

			r_a = Y(1,pp)/P%a

			SELECT CASE (TRIM(P%ne_profile))
				CASE('TANH')
					fr = 1_ip - TANH(2.0_rp*r_a)**P%n_ne
					ne(pp) = P%neo*fr
				CASE('FLAT')
					ne(pp) = P%neo
				CASE('PARABOLIC')
					fr = P%a_ne(1) + P%a_ne(2)*r_a + P%a_ne(3)*r_a**2 + P%a_ne(4)*r_a**3
					ne(pp) = P%neo*fr
				CASE DEFAULT
					ne(pp) = P%neo
			END SELECT

			SELECT CASE (TRIM(P%Te_profile))
				CASE('TANH')
					fr = 1_ip - TANH(2.0_rp*r_a)**P%n_Te
					Te(pp) = P%Teo*fr
				CASE('FLAT')
					Te(pp) = P%Teo
				CASE('PARABOLIC')
					fr = P%a_Te(1) + P%a_Te(2)*r_a + P%a_Te(3)*r_a**2 + P%a_Te(4)*r_a**3
					Te(pp) = P%Teo*fr
				CASE DEFAULT
					Te(pp) = P%Teo
			END SELECT

			SELECT CASE (TRIM(P%Zeff_profile))
				CASE('TANH')
					fr = 1_ip - TANH(2.0_rp*r_a)**P%n_Zeff
					Zeff(pp) = P%Zeffo*fr
				CASE('FLAT')
					Zeff(pp) = P%Zeffo
				CASE('PARABOLIC')
					fr = P%a_Zeff(1) + P%a_Zeff(2)*r_a + P%a_Zeff(3)*r_a**2 + P%a_Zeff(4)*r_a**3
					Zeff(pp) = P%Zeffo*fr
				CASE DEFAULT
					Zeff(pp) = P%Zeffo
			END SELECT
        end if
	end do
!$OMP END PARALLEL DO
end subroutine get_analytical_profiles


subroutine get_profiles(params,vars,P)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(PROFILES), INTENT(IN) :: P

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			call get_analytical_profiles(P,vars%Y,vars%ne,vars%Te,vars%Zeff,vars%flag)
		CASE('EXTERNAL')
			call interp_profiles(vars,P)
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
	REAL(rp) :: rdatum
	INTEGER :: h5error
    INTEGER ir, iphi, iz

	filename = TRIM(P%filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_profiles_data_from_hdf5 --> h5fopen_f")')
	end if


	dset = "/NR"
	call load_from_hdf5(h5file_id,dset,rdatum)
	P%dims(1) = INT(rdatum)

	if (P%axisymmetric) then
		P%dims(2) = 0
	else
		dset = "/NPHI"
		call load_from_hdf5(h5file_id,dset,rdatum)
		P%dims(2) = INT(rdatum)
	end if

	dset = "/NZ"
	call load_from_hdf5(h5file_id,dset,rdatum)
	P%dims(3) = INT(rdatum)

	if (P%axisymmetric) then
		call ALLOCATE_2D_PROFILES_ARRAYS(P)
	else
		call ALLOCATE_3D_PROFILES_ARRAYS(P)
	end if

	dset = "/R"
	call load_array_from_hdf5(h5file_id,dset,P%X%R)

	if (.NOT.P%axisymmetric) then
		dset = "/PHI"
		call load_array_from_hdf5(h5file_id,dset,P%X%PHI)
	end if

	dset = "/Z"
	call load_array_from_hdf5(h5file_id,dset,P%X%Z)

	dset = "/FLAG"
	if (P%axisymmetric) then
		call load_array_from_hdf5(h5file_id,dset,P%FLAG2D)
	else
		call load_array_from_hdf5(h5file_id,dset,P%FLAG3D)
	end if

	dset = "/ne"
	if (P%axisymmetric) then
		call load_array_from_hdf5(h5file_id,dset,P%ne_2D)
	else
		call load_array_from_hdf5(h5file_id,dset,P%ne_3D)
	end if

	dset = "/Te"
	if (P%axisymmetric) then
		call load_array_from_hdf5(h5file_id,dset,P%Te_2D)
		P%Te_2D = P%Te_2D*C_E
	else
		call load_array_from_hdf5(h5file_id,dset,P%Te_3D)
		P%Te_3D = P%Te_3D*C_E
	end if

	dset = "/Zeff"
	if (P%axisymmetric) then
		call load_array_from_hdf5(h5file_id,dset,P%Zeff_2D)
	else
		call load_array_from_hdf5(h5file_id,dset,P%Zeff_3D)
	end if

	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_profiles_data_from_hdf5 --> h5fclose_f")')
	end if
end subroutine load_profiles_data_from_hdf5


subroutine ALLOCATE_2D_PROFILES_ARRAYS(P)
	TYPE(PROFILES), INTENT(INOUT) :: P

	ALLOCATE(P%X%R(P%dims(1)))
	ALLOCATE(P%X%Z(P%dims(3)))
	ALLOCATE(P%FLAG2D(P%dims(1),P%dims(3)))
	ALLOCATE(P%ne_2D(P%dims(1),P%dims(3)))
	ALLOCATE(P%Te_2D(P%dims(1),P%dims(3)))
	ALLOCATE(P%Zeff_2D(P%dims(1),P%dims(3)))
end subroutine ALLOCATE_2D_PROFILES_ARRAYS


subroutine ALLOCATE_3D_PROFILES_ARRAYS(P)
	TYPE(PROFILES), INTENT(INOUT) :: P

	ALLOCATE(P%X%R(P%dims(1)))
	ALLOCATE(P%X%PHI(P%dims(2)))
	ALLOCATE(P%X%Z(P%dims(3)))
	ALLOCATE(P%FLAG3D(P%dims(1),P%dims(2),P%dims(3)))
	ALLOCATE(P%ne_3D(P%dims(1),P%dims(2),P%dims(3)))
	ALLOCATE(P%Te_3D(P%dims(1),P%dims(2),P%dims(3)))
	ALLOCATE(P%Zeff_3D(P%dims(1),P%dims(2),P%dims(3)))
end subroutine ALLOCATE_3D_PROFILES_ARRAYS

end module korc_profiles
