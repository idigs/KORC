module korc_fields
    use korc_types
	use korc_hpc
	use korc_coords
	use korc_interp
	use korc_HDF5

    implicit none

	PUBLIC :: mean_F_field,get_fields,load_field_data_from_hdf5,load_dim_data_from_hdf5
	PRIVATE :: get_analytical_fields,analytical_fields,analytical_magnetic_field,analytical_electric_field,uniform_magnetic_field,&
				uniform_electric_field,check_if_confined,uniform_fields,cross

    contains

subroutine analytical_fields(F,Y,E,B,flag)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: E ! E(1,:) = Ex, E(2,:) = Ey, E(3,:) = Ez
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
	REAL(rp) :: Ezeta, Bp, Bt, eta, q
	INTEGER(ip) pp ! Iterator(s)
	INTEGER(ip) :: ss

	ss = SIZE(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ezeta,Bp,Bt,eta,q) SHARED(F,Y,E,B,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_idef ) then
		    eta = Y(1,pp)/F%Ro
            q = F%AB%qo*(1.0_rp + (Y(1,pp)/F%AB%lambda)**2)
            Bp = eta*F%AB%Bo/(q*(1.0_rp + eta*COS(Y(2,pp))))
		    Bt = F%AB%Bo/( 1.0_rp + eta*COS(Y(2,pp)) )
			

		    B(1,pp) =  Bt*COS(Y(3,pp)) - Bp*SIN(Y(2,pp))*SIN(Y(3,pp))
		    B(2,pp) = -Bt*SIN(Y(3,pp)) - Bp*SIN(Y(2,pp))*COS(Y(3,pp))
		    B(3,pp) = Bp*COS(Y(2,pp))

			if (abs(F%Eo) > 0) then
				Ezeta = F%Eo/( 1.0_rp + eta*COS(Y(2,pp)) )

				E(1,pp) = Ezeta*COS(Y(3,pp))
				E(2,pp) = -Ezeta*SIN(Y(3,pp))
				E(3,pp) = 0.0_rp
			end if
        end if
	end do
!$OMP END PARALLEL DO
end subroutine analytical_fields

subroutine analytical_magnetic_field(F,Y,B,flag)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
	REAL(rp) :: Bp, Bt, eta, q
	INTEGER(ip) pp ! Iterator(s)
	INTEGER(ip) :: ss

	ss = SIZE(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Bp,Bt,eta,q) SHARED(F,Y,B,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_idef ) then
		    eta = Y(1,pp)/F%Ro
            q = F%AB%qo*(1.0_rp + (Y(1,pp)/F%AB%lambda)**2)
            Bp = eta*F%AB%Bo/(q*(1.0_rp + eta*COS(Y(2,pp))))
		    Bt = F%AB%Bo/( 1.0_rp + eta*COS(Y(2,pp)) )

		    B(1,pp) =  Bt*COS(Y(3,pp)) - Bp*SIN(Y(2,pp))*SIN(Y(3,pp))
		    B(2,pp) = -Bt*SIN(Y(3,pp)) - Bp*SIN(Y(2,pp))*COS(Y(3,pp))
		    B(3,pp) = Bp*COS(Y(2,pp))
        end if
	end do
!$OMP END PARALLEL DO
end subroutine analytical_magnetic_field


subroutine uniform_magnetic_field(F,B)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz

	B(1,:) = F%Bo
	B(2:3,:) = 0.0_rp
end subroutine uniform_magnetic_field


subroutine uniform_electric_field(F,E)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: E ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz

	E(1,:) = F%Eo
	E(2:3,:) = 0.0_rp
end subroutine uniform_electric_field


subroutine analytical_electric_field(F,Y,E,flag)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: E ! E(1,:) = Ex, E(2,:) = Ey, E(3,:) = Ez
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
	REAL(rp) :: Ezeta, eta
	INTEGER(ip) pp ! Iterator(s)
	INTEGER(ip) :: ss

	if (abs(F%Eo) > 0) then
		ss = SIZE(Y,2)
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ezeta,eta) SHARED(F,Y,E,flag)
		do pp=1_idef,ss
            if ( flag(pp) .EQ. 1_idef ) then
			    eta = Y(1,pp)/F%Ro		
			    Ezeta = F%Eo/( 1.0_rp + eta*COS(Y(2,pp)) )

			    E(1,pp) = Ezeta*COS(Y(3,pp))
			    E(2,pp) = -Ezeta*SIN(Y(3,pp))
			    E(3,pp) = 0.0_rp
            end if
		end do
!$OMP END PARALLEL DO
	end if
end subroutine analytical_electric_field


subroutine mean_F_field(F,Fo,op_field)
	implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), INTENT(OUT) :: Fo
	TYPE(KORC_STRING), INTENT(IN) :: op_field

	if (TRIM(op_field%str) .EQ. 'B') then
		if (ALLOCATED(F%B_3D%R)) then ! 3D field
			Fo = SUM( SQRT(F%B_3D%R**2 + F%B_3D%PHI**2 + F%B_3D%Z**2) )/SIZE(F%B_3D%R)
		else if (ALLOCATED(F%B_2D%R)) then ! Axisymmetric 2D field
			Fo = SUM( SQRT(F%B_2D%R**2 + F%B_2D%PHI**2 + F%B_2D%Z**2) )/SIZE(F%B_2D%R)
		end if
	else if (TRIM(op_field%str) .EQ. 'E') then
		if (ALLOCATED(F%E_3D%R)) then ! 3D field
			Fo = SUM( SQRT(F%E_3D%R**2 + F%E_3D%PHI**2 + F%E_3D%Z**2) )/SIZE(F%E_3D%R)
		else if (ALLOCATED(F%E_2D%R)) then ! Axisymmetric 2D field
			Fo = SUM( SQRT(F%E_2D%R**2 + F%E_2D%PHI**2 + F%E_2D%Z**2) )/SIZE(F%E_2D%R)
		end if	
	else
		write(6,'("KORC ERROR: Please enter a valid field: mean_F_field")')
		call korc_abort()
	end if
end subroutine mean_F_field


subroutine check_if_confined(F,Y,flag)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	INTEGER(ip) :: pp,ss

    ss = SIZE(Y,2)
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(F,Y,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_idef ) then
            if (Y(1,pp) .GT. F%AB%a) then
                flag(pp) = 0_idef
            end if
        end if
	end do
!$OMP END PARALLEL DO
end subroutine check_if_confined


subroutine get_analytical_fields(vars,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(FIELDS), INTENT(IN) :: F

!	call cart_to_tor(vars%X, F%AB%Ro, vars%Y, vars%flag)

!	call check_if_confined(F, vars%Y, vars%flag)

	call cart_to_tor_check_if_confined(vars%X,F,vars%Y,vars%flag)

	call analytical_fields(F,vars%Y, vars%E, vars%B, vars%flag)

!	call analytical_magnetic_field(F, vars%Y, vars%B, vars%flag)

!	call analytical_electric_field(F, vars%Y, vars%E, vars%flag)
end subroutine get_analytical_fields


subroutine uniform_fields(vars,F)
    implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(FIELDS), INTENT(IN) :: F

	call uniform_magnetic_field(F, vars%B)

	call uniform_electric_field(F, vars%E)
end subroutine uniform_fields


function cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


subroutine unitVectors(params,Xo,F,b1,b2,b3,flag)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Xo
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: b1
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: b2
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: b3
	INTEGER, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: flag
	TYPE(PARTICLES) :: vars
	REAL(rp), PARAMETER :: tol = korc_zero
	INTEGER :: ii, ppp

	ppp = SIZE(Xo,2) ! Number of particles

	ALLOCATE( vars%X(3,ppp) )
	ALLOCATE( vars%Y(3,ppp) )
	ALLOCATE( vars%B(3,ppp) )
	ALLOCATE( vars%E(3,ppp) )
    ALLOCATE( vars%flag(ppp) )

	vars%X = Xo
    vars%flag = 1_idef
	
	call init_random_seed()

	call get_fields(params,vars,F)
	
	do ii=1_idef,ppp
		if ( vars%flag(ii) .EQ. 1_idef ) then
			b1(:,ii) = vars%B(:,ii)/SQRT( DOT_PRODUCT(vars%B(:,ii),vars%B(:,ii)) )

		    b2(:,ii) = cross(b1(:,ii),(/0.0_rp,0.0_rp,1.0_rp/))
		    b2(:,ii) = b2(:,ii)/SQRT( DOT_PRODUCT(b2(:,ii),b2(:,ii)) )

		    b3(:,ii) = cross(b1(:,ii),b2(:,ii))
		    b3(:,ii) = b3(:,ii)/SQRT( DOT_PRODUCT(b3(:,ii),b3(:,ii)) )
		end if
	end do

	if (PRESENT(flag)) then
		flag = vars%flag
	end if

	DEALLOCATE( vars%X )
	DEALLOCATE( vars%Y )
	DEALLOCATE( vars%B )
	DEALLOCATE( vars%E )
    DEALLOCATE( vars%flag )
end subroutine unitVectors


subroutine get_fields(params,vars,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(FIELDS), INTENT(IN) :: F

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			call get_analytical_fields(vars, F)
		CASE('EXTERNAL')
			call interp_field(vars, F)
		CASE('UNIFORM')
			call uniform_fields(vars, F)
		CASE DEFAULT
	END SELECT
end subroutine get_fields


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! Subroutines for getting the fields data from HDF5 files
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

subroutine load_dim_data_from_hdf5(params,field_dims)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	INTEGER, DIMENSION(3), INTENT(OUT) :: field_dims
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	INTEGER :: h5error
	REAL(rp) :: rdatum
    INTEGER :: ii

	filename = TRIM(params%magnetic_field_filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fopen_f")')
	end if

	if (params%poloidal_flux.OR.params%axisymmetric) then
			dset = "/NR"
			call load_from_hdf5(h5file_id,dset,rdatum)
			field_dims(1) = INT(rdatum)

			field_dims(2) = 0

			dset = "/NZ"
			call load_from_hdf5(h5file_id,dset,rdatum)
			field_dims(3) = INT(rdatum)
	else
			dset = "/NR"
			call load_from_hdf5(h5file_id,dset,rdatum)
			field_dims(1) = INT(rdatum)

			dset = "/NPHI"
			call load_from_hdf5(h5file_id,dset,rdatum)
			field_dims(2) = INT(rdatum)

			dset = "/NZ"
			call load_from_hdf5(h5file_id,dset,rdatum)
			field_dims(3) = INT(rdatum)
	end if


	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fclose_f")')
	end if
end subroutine load_dim_data_from_hdf5


subroutine load_field_data_from_hdf5(params,F)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(INOUT) :: F
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

	filename = TRIM(params%magnetic_field_filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
	end if

	dset = "/R"
	call load_array_from_hdf5(h5file_id,dset,F%X%R)

	if ((.NOT.params%poloidal_flux).AND.(.NOT.params%axisymmetric)) then
		dset = "/PHI"
		call load_array_from_hdf5(h5file_id,dset,F%X%PHI)
	end if

	dset = "/Z"
	call load_array_from_hdf5(h5file_id,dset,F%X%Z)

	dset = '/Bo'
	call load_from_hdf5(h5file_id,dset,F%Bo)

	dset = '/Ro'
	call load_from_hdf5(h5file_id,dset,F%Ro)

	dset = '/Zo'
	call load_from_hdf5(h5file_id,dset,F%Zo)

	if ((.NOT.params%poloidal_flux).AND.(.NOT.params%axisymmetric)) then
		dset = "/BR"
		call load_array_from_hdf5(h5file_id,dset,F%B_3D%R)

		dset = "/BPHI"
		call load_array_from_hdf5(h5file_id,dset,F%B_3D%PHI)

		dset = "/BZ"
		call load_array_from_hdf5(h5file_id,dset,F%B_3D%Z)

		dset = "/FLAG"
		call load_array_from_hdf5(h5file_id,dset,F%FLAG3D)
	else if (params%poloidal_flux) then
		dset = "/PSIp"
		call load_array_from_hdf5(h5file_id,dset,F%PSIp)

		dset = "/FLAG"
		call load_array_from_hdf5(h5file_id,dset,F%FLAG2D)
	else if (params%axisymmetric) then

		dset = "/BR"
		call load_array_from_hdf5(h5file_id,dset,F%B_2D%R)

		dset = "/BPHI"
		call load_array_from_hdf5(h5file_id,dset,F%B_2D%PHI)

		dset = "/BZ"
		call load_array_from_hdf5(h5file_id,dset,F%B_2D%Z)

		dset = "/FLAG"
		call load_array_from_hdf5(h5file_id,dset,F%FLAG2D)
	end if

	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
	end if

end subroutine load_field_data_from_hdf5

end module korc_fields
