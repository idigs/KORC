!> @brief Module containing subroutines to initialize externally generated fields, and to calculate the electric and magnetic fields when using an analytical model.
module korc_fields
    use korc_types
	use korc_hpc
	use korc_coords
	use korc_interp
	use korc_HDF5

    IMPLICIT NONE

	PUBLIC :: mean_F_field,&
				get_fields,&
				initialize_fields,&
				load_field_data_from_hdf5,&
				load_dim_data_from_hdf5,&
				ALLOCATE_2D_FIELDS_ARRAYS,&
				ALLOCATE_3D_FIELDS_ARRAYS,&
				DEALLOCATE_FIELDS_ARRAYS
	PRIVATE :: get_analytical_fields,&
				analytical_fields,&
				uniform_magnetic_field,&
				uniform_electric_field,&
				uniform_fields,&
				cross,&
				analytical_electric_field_cyl,&
				ALLOCATE_V_FIELD_2D,&
				ALLOCATE_V_FIELD_3D

    CONTAINS

!> @brief Subroutine that calculates and returns the electric and magnetic field for each particle in the simulation.
!! @details The analytical magnetic field is given by:
!!
!! @f$\vec{B}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}} \left[ B_0 \hat{e}_\zeta  + B_\vartheta(r) \hat{e}_\vartheta \right]@f$,
!!
!! where @f$\eta = r/R_0@f$ is the aspect ratio, the constant @f$B_0@f$ denotes the magnitude of the toroidal magnetic field,
!! and @f$B_\vartheta(r) = \eta B_0/q(r)@f$ is the poloidal magnetic field with safety factor @f$q(r) = q_0\left( 1 + \frac{r^2}{\lambda^2} \right)@f$.
!! The constant @f$q_0@f$ is the safety factor at the magnetic axis and the constant @f$\lambda@f$ is obtained from the values of @f$q_0@f$
!! and @f$q(r)@f$ at the plasma edge @f$r=r_{edge}@f$.
!! On the other hand, the analytical electric fields is given by:
!!
!! @f$\vec{E}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}} E_0 \hat{e}_\zeta@f$,
!!
!! where @f$E_0@f$ is the electric field as measured at the mangetic axis.
!!
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in] Y Toroidal coordinates of each particle in the simulation; Y(1,:) = @f$r@f$, Y(2,:) = @f$\theta@f$, Y(3,:) = @f$\zeta@f$.
!! @paramp[in,out] B Magnetic field components in Cartesian coordinates; B(1,:) = @f$B_x@f$, B(2,:) = @f$B_y@f$, B(3,:) = @f$B_z@f$
!! @param[in,out] E Electric field components in Cartesian coordinates; E(1,:) = @f$E_x@f$, E(2,:) = @f$E_y@f$, E(3,:) = @f$E_z@f$
!! @param[in] flag Flag for each particle to decide whether it is being followed (flag=T) or not (flag=F).
!! @param Ezeta Toroidal electric field.
!! @param Bzeta Toroidal magnetic field.
!! @param Bp Poloidal magnetic field.
!! @param eta Aspect ratio.
!! @param q Safety profile.
!! @param pp Particle iterator.
!! @param ss Particle species iterator.
subroutine analytical_fields(F,Y,E,B,flag)
	TYPE(FIELDS), INTENT(IN)                               :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
	REAL(rp)                                               :: Ezeta
    REAL(rp)                                               :: Bzeta
    REAL(rp)                                               :: Bp
    REAL(rp)                                               :: eta
    REAL(rp)                                               :: q
	INTEGER(ip)                                            :: pp ! Iterator(s)
	INTEGER(ip)                                            :: ss

	ss = SIZE(Y,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ezeta,Bp,Bzeta,eta,q) SHARED(F,Y,E,B,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
		    eta = Y(1,pp)/F%Ro
            q = F%AB%qo*(1.0_rp + (Y(1,pp)/F%AB%lambda)**2)
			Bp = F%AB%Bp_sign*eta*F%AB%Bo/(q*(1.0_rp + eta*COS(Y(2,pp))))
		    Bzeta = F%AB%Bo/( 1.0_rp + eta*COS(Y(2,pp)) )


		    B(1,pp) =  Bzeta*COS(Y(3,pp)) - Bp*SIN(Y(2,pp))*SIN(Y(3,pp))
		    B(2,pp) = -Bzeta*SIN(Y(3,pp)) - Bp*SIN(Y(2,pp))*COS(Y(3,pp))
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


!> @brief Subroutine that returns the value of a uniform magnetic field.
!! @details This subroutie is used only when the simulation is ran for a 'UNIFORM' plasma. As a convention, in a uniform plasma we set
!! @f$\vec{B} = B_0 \hat{x}@f$.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in,out] B Magnetic field components in Cartesian coordinates; B(1,:) = @f$B_x@f$, B(2,:) = @f$B_y@f$, B(3,:) = @f$B_z@f$
subroutine uniform_magnetic_field(F,B)
	TYPE(FIELDS), INTENT(IN)                               :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B

	B(1,:) = F%Bo
	B(2:3,:) = 0.0_rp
end subroutine uniform_magnetic_field


!> @brief Subroutine that returns the value of a uniform electric field.
!! @details This subroutie is used only when the simulation is ran for a 'UNIFORM' plasma. As a convention, in a uniform plasma we set
!! @f$\vec{E} = E_0 \hat{x}@f$.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in,out] E Electric field components in Cartesian coordinates; E(1,:) = @f$E_x@f$, E(2,:) = @f$E_y@f$, E(3,:) = @f$E_z@f$
subroutine uniform_electric_field(F,E)
	TYPE(FIELDS), INTENT(IN)                               :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E

	E(1,:) = F%Eo
	E(2:3,:) = 0.0_rp
end subroutine uniform_electric_field


!> @brief Subrotuine that calculates and returns the electric field using the same analytical model of the 'analytical_fields' subroutine.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in] Y Cylindrical coordinates of each particle in the simulation; Y(1,:) = @f$R@f$, Y(2,:) = @f$\phi@f$, Y(3,:) = @f$Z@f$.
!! @param[in,out] E Electric field components in Cartesian coordinates; E(1,:) = @f$E_x@f$, E(2,:) = @f$E_y@f$, E(3,:) = @f$E_z@f$
!! @param[in] flag Flag for each particle to decide whether it is being followed (flag=T) or not (flag=F).
!! @param Ephi Azimuthal electric field.
!! @param pp Particle iterator.
!! @param ss Particle species iterator.
subroutine analytical_electric_field_cyl(F,Y,E,flag)
	TYPE(FIELDS), INTENT(IN)                               :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
	REAL(rp)                                               :: Ephi
	INTEGER(ip)                                            :: pp
	INTEGER(ip)                                            :: ss

	if (abs(F%Eo) > 0) then
		ss = SIZE(Y,2)
!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ephi) SHARED(F,Y,E,flag)
		do pp=1_idef,ss
            if ( flag(pp) .EQ. 1_is ) then
			    Ephi = F%Eo*F%Ro/Y(1,pp)

			    E(1,pp) = -Ephi*SIN(Y(2,pp))
			    E(2,pp) = Ephi*COS(Y(2,pp))
			    E(3,pp) = 0.0_rp
            end if
		end do
!$OMP END PARALLEL DO
	end if
end subroutine analytical_electric_field_cyl


!> @brief Subroutine that calculates the mean electric or magnetic field in case external fields are being used.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[out] Fo Mean electric or magnetic field.
!! @param[in] op_field String that specifies what mean field will be calculated. Its value can be 'B' or 'E'.
subroutine mean_F_field(F,Fo,op_field)
	TYPE(FIELDS), INTENT(IN)       :: F
	REAL(rp), INTENT(OUT)          :: Fo
	TYPE(KORC_STRING), INTENT(IN)  :: op_field

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


!> @brief Interface for calculating the analytical electric and magnetic fields for each particle in the simulation.
!! @param[in,out] vars An instance of the KORC derived type PARTICLES.
!! @param[in] F An instance of the KORC derived type FIELDS.
subroutine get_analytical_fields(vars,F)
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(FIELDS), INTENT(IN)       :: F

	call cart_to_tor_check_if_confined(vars%X,F,vars%Y,vars%flag)

	call analytical_fields(F,vars%Y, vars%E, vars%B, vars%flag)
end subroutine get_analytical_fields


!> @brief Interface for calculating the uniform electric and magnetic fields for each particle in the simulation.
!! @param[in,out] vars An instance of the KORC derived type PARTICLES.
!! @param[in] F An instance of the KORC derived type FIELDS.
subroutine uniform_fields(vars,F)
	TYPE(PARTICLES), INTENT(INOUT) :: vars
	TYPE(FIELDS), INTENT(IN)       :: F

	call uniform_magnetic_field(F, vars%B)

	call uniform_electric_field(F, vars%E)
end subroutine uniform_fields


!> @brief Function that calculates the cross product of the two vectors @f$\vec{a}@f$ and @f$\vec{b}@f$.
!! @param cross Cross product @f$\vec{a}\times \vec{b}@f$
!! @param[in] a Vector @f$\vec{a}@f$.
!! @param[in] b Vector @f$\vec{b}@f$.
function cross(a,b)
    REAL(rp), DIMENSION(3)             :: cross
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


!> @brief Subrotuine that calculates an orthonormal basis using information of the (local) magnetic field at position @f$\vec{X}_0@f$.
!! @param[in] params Core KORC simulation parameters.
!! @param[in] Xo Array with the position of the simulated particles.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in,out] b1 Basis vector pointing along the local magnetic field, that is, along @f$\vec{b} = \vec{B}/B@f$.
!! @param[in,out] b2 Basis vector perpendicular to b1
!! @param[in,out] b3 Basis vector perpendicular to b1 and b2.
!! @param[in,out] flag Flag for each particle to decide whether it is being followed (flag=T) or not (flag=F).
!! @param vars A temporary instance of the KORC derived type PARTICLES.
!! @param ii Iterator.
!! @param ppp Number of particles.
subroutine unitVectors(params,Xo,F,b1,b2,b3,flag)
	TYPE(KORC_PARAMS), INTENT(IN)                                      :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)                  :: Xo
	TYPE(FIELDS), INTENT(IN)                                           :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b1
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b2
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b3
	INTEGER(is), DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT)    :: flag
	TYPE(PARTICLES)                                                    :: vars
	INTEGER                                                            :: ii
    INTEGER                                                            :: ppp

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


!> @brief Inferface with calls to subroutines for calculating the electric and magnetic field for each particle in the simulation.
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] vars An instance of the KORC derived type PARTICLES.
!! @param[in] F An instance of the KORC derived type FIELDS.
subroutine get_fields(params,vars,F)
	TYPE(KORC_PARAMS), INTENT(IN)      :: params
	TYPE(PARTICLES), INTENT(INOUT)     :: vars
	TYPE(FIELDS), INTENT(IN)           :: F

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			call get_analytical_fields(vars, F)
		CASE('EXTERNAL')
			call interp_fields(vars, F)
			if (F%Efield.AND..NOT.F%Efield_in_file) then
				call analytical_electric_field_cyl(F,vars%Y,vars%E,vars%flag)
			end if
		CASE('UNIFORM')
			call uniform_fields(vars, F)
		CASE DEFAULT
	END SELECT
end subroutine get_fields


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * *  SUBROUTINES FOR INITIALIZING FIELDS   * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

!> @brief Subroutine that initializes the analytical or externally calculated electric and magnetic fields.
!! @details In this subroutine we load the parameters of the electric and magnetic fields from the namelists 'analytical_fields_params' and
!! 'externalPlasmaModel' in the input file.
!! @param[in] params Core KORC simulation parameters.
!! @param[out] F An instance of the KORC derived type FIELDS.
!! @param Bo Magnetic field at magnetic axis for an 'ANALITICAL' magnetic field, or the magnitude of the magnetic field for a 'UNFIROM' plasma.
!! @param minor_radius Plasma edge @f$r_{edge}@f$ as measured from the magnetic axis.
!! @param major_radius Radial position of the magnetic axis @f$R_0@f$
!! @param qa Safety factor at the plasma edge.
!! @param qo Safety factor at the magnetic axis @f$q_0@f$.
!! @param current_direction String with information about the direction of the plasma current, 'PARALLEL'  or 'ANTI-PARALLEL' to the toroidal magnetic field.
!! @param Eo Electric field at the magnetic axis.
!! @param Efield Logical variable that specifies if the electric field is going to be used on in a given simulation.
!! @param Bfield Logical variable that specifies if the magnetic field is going to be used on in a given simulation.
!! @param Bflux Logical variable that specifies if the poloidal magnetic flux is going to be used on in a given simulation.
!! @param axisymmetric_fields Logical variable that specifies if the plasma is axisymmetric.
subroutine initialize_fields(params,F)
	TYPE(KORC_PARAMS), INTENT(IN)  :: params
	TYPE(FIELDS), INTENT(OUT)      :: F
	REAL(rp)                       :: Bo
	REAL(rp)                       :: minor_radius
	REAL(rp)                       :: major_radius
	REAL(rp)                       :: qa
	REAL(rp)                       :: qo
    CHARACTER(MAX_STRING_LENGTH)   :: current_direction
	REAL(rp)                       :: Eo
	LOGICAL                        :: Efield
    LOGICAL                        :: Bfield
    LOGICAL                        :: Bflux
    LOGICAL                        :: axisymmetric_fields

	NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
			qa,qo,Eo,current_direction

	NAMELIST /externalPlasmaModel/ Efield, Bfield, Bflux, axisymmetric_fields

	SELECT CASE (TRIM(params%plasma_model))
		CASE('ANALYTICAL')
			! Load the parameters of the analytical magnetic field
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=analytical_fields_params)
			close(default_unit_open)

			F%AB%Bo = Bo
			F%AB%a = minor_radius
			F%AB%Ro = major_radius
			F%Ro = major_radius
			F%Zo = 0.0_rp
			F%AB%qa = qa
			F%AB%qo = qo
			F%AB%lambda = F%AB%a/SQRT(qa/qo - 1.0_rp)
			F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)
			F%AB%current_direction = TRIM(current_direction)
			SELECT CASE (TRIM(F%AB%current_direction))
				CASE('PARALLEL')
					F%AB%Bp_sign = 1.0_rp
				CASE('ANTI-PARALLEL')
					F%AB%Bp_sign = -1.0_rp
				CASE DEFAULT
			END SELECT
			F%Eo = Eo
			F%Bo = F%AB%Bo

		CASE('EXTERNAL')
			! Load the magnetic field from an external HDF5 file
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=externalPlasmaModel)
			close(default_unit_open)

			F%Bfield = Bfield
			F%Bflux = Bflux
			F%Efield = Efield
			F%axisymmetric_fields = axisymmetric_fields

		    call load_dim_data_from_hdf5(params,F)

			call which_fields_in_file(params,F%Bfield_in_file,F%Efield_in_file,F%Bflux_in_file)


			if (F%Bflux.AND..NOT.F%Bflux_in_file) then
				write(6,'("ERROR: Magnetic flux to be used but no data in file!")')
				call KORC_ABORT()
			end if

			if (F%Bfield.AND..NOT.F%Bfield_in_file) then
				write(6,'("ERROR: Magnetic field to be used but no data in file!")')
				call KORC_ABORT()
			end if

			if (F%Efield.AND..NOT.F%Efield_in_file) then
				if (params%mpi_params%rank.EQ.0_idef) then
					write(6,'(/,"* * * * * * * * * *  FIELDS  * * * * * * * * * *")')
					write(6,'("MESSAGE: Analytical electric field will be used.")')
					write(6,'("* * * * * * * * * * * * ** * * * * * * * * * * *",/)')
				end if
			end if

			if (F%axisymmetric_fields) then
				call ALLOCATE_2D_FIELDS_ARRAYS(F,F%Bfield,F%Bflux,F%Efield.AND.F%Efield_in_file)
			else
				call ALLOCATE_3D_FIELDS_ARRAYS(F,F%Bfield,F%Efield)
			end if

		    call load_field_data_from_hdf5(params,F)
		CASE DEFAULT
	END SELECT
end subroutine initialize_fields


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! Subroutines for getting the fields data from HDF5 files
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!> @brief Subroutine that loads the size of the arrays having the electric and magnetic field data.
!! @details All the information of externally calculated fields must be given in a rectangular, equally spaced mesh in the @f$(R,\phi,Z)@f$ space of cylindrical coordinates.
!! If the fields are axisymmetric, then the fields must be in a rectangular mesh on the @f$RZ@f$-plane.
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] F An instance of the KORC derived type FIELDS.
!! @param filename String containing the name of the HDF5 file.
!! @param gname String containing the group name of a parameter in the HDF5 file.
!! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
!! @param dset Name of data set to read from file.
!! @param h5file_id HDF5 file identifier.
!! @param group_id HDF5 group identifier.
!! @param subgroup_id HDF5 subgroup identifier.
!! @dims Array containing the size of the mesh with the data of the electric and magnetic fields. dims(1) = dimension along the @f$R@f$ coordinate,
!! dims(2) = dimension along the @f$\phi@f$ coordinate, and dims(3) = dimension along the @f$Z@f$ coordinate.
!! @param h5error HDF5 error status.
!! @param rdamum Temporary variable keeping the read data.
subroutine load_dim_data_from_hdf5(params,F)
	TYPE(KORC_PARAMS), INTENT(IN)                  :: params
	TYPE(FIELDS), INTENT(INOUT)                    :: F
	CHARACTER(MAX_STRING_LENGTH)                   :: filename
	CHARACTER(MAX_STRING_LENGTH)                   :: gname
	CHARACTER(MAX_STRING_LENGTH)                   :: subgname
	CHARACTER(MAX_STRING_LENGTH)                   :: dset
	INTEGER(HID_T)                                 :: h5file_id
	INTEGER(HID_T)                                 :: group_id
	INTEGER(HID_T)                                 :: subgroup_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE    :: dims
	INTEGER                                        :: h5error
	REAL(rp)                                       :: rdatum

	filename = TRIM(params%magnetic_field_filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fopen_f")')
	end if

	if (F%Bflux.OR.F%axisymmetric_fields) then
			dset = "/NR"
			call load_from_hdf5(h5file_id,dset,rdatum)
			F%dims(1) = INT(rdatum)

			F%dims(2) = 0

			dset = "/NZ"
			call load_from_hdf5(h5file_id,dset,rdatum)
			F%dims(3) = INT(rdatum)
	else
			dset = "/NR"
			call load_from_hdf5(h5file_id,dset,rdatum)
			F%dims(1) = INT(rdatum)

			dset = "/NPHI"
			call load_from_hdf5(h5file_id,dset,rdatum)
			F%dims(2) = INT(rdatum)

			dset = "/NZ"
			call load_from_hdf5(h5file_id,dset,rdatum)
			F%dims(3) = INT(rdatum)
	end if


	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fclose_f")')
	end if
end subroutine load_dim_data_from_hdf5


!> @brief Subroutine that queries the HDF5 file what data are present in the HDF5 input file (sanity check).
!! @param[in] params Core KORC simulation parameters.
!! @param Bfield Logical variable that specifies if the magnetic field is present in the HDF5 file.
!! @param Efield Logical variable that specifies if the electric field is present in the HDF5 file.
!! @param Bflux Logical variable that specifies if the poloidal magnetic flux is present in the HDF5 file.
!! @param filename String containing the name of the HDF5 file.
!! @param gname String containing the group name of a parameter in the HDF5 file.
!! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
!! @param dset Name of data set to read from file.
!! @param h5file_id HDF5 file identifier.
!! @param group_id HDF5 group identifier.
!! @param subgroup_id HDF5 subgroup identifier.
!! @param h5error HDF5 error status.
subroutine which_fields_in_file(params,Bfield,Efield,Bflux)
	TYPE(KORC_PARAMS), INTENT(IN)      :: params
	LOGICAL, INTENT(OUT)               :: Bfield
	LOGICAL, INTENT(OUT)               :: Efield
	LOGICAL, INTENT(OUT)               :: Bflux
	CHARACTER(MAX_STRING_LENGTH)       :: filename
	CHARACTER(MAX_STRING_LENGTH)       :: gname
	CHARACTER(MAX_STRING_LENGTH)       :: subgname
	CHARACTER(MAX_STRING_LENGTH)       :: dset
	INTEGER(HID_T)                     :: h5file_id
	INTEGER(HID_T)                     :: group_id
	INTEGER(HID_T)                     :: subgroup_id
	INTEGER                            :: h5error

	filename = TRIM(params%magnetic_field_filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
	end if

	gname = "BR"
	call h5lexists_f(h5file_id,TRIM(gname),Bfield,h5error)

	gname = "ER"
	call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

	gname = "PSIp"
	call h5lexists_f(h5file_id,TRIM(gname),Bflux,h5error)


	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
	end if
end subroutine which_fields_in_file


!> @brief Subroutine that loads the fields data from the HDF5 input file.
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
!! @param filename String containing the name of the HDF5 file.
!! @param gname String containing the group name of a parameter in the HDF5 file.
!! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
!! @param dset Name of data set to read from file.
!! @param h5file_id HDF5 file identifier.
!! @param group_id HDF5 group identifier.
!! @param subgroup_id HDF5 subgroup identifier.
!! @param h5error HDF5 error status.
subroutine load_field_data_from_hdf5(params,F)
	TYPE(KORC_PARAMS), INTENT(IN)          :: params
	TYPE(FIELDS), INTENT(INOUT)            :: F
	CHARACTER(MAX_STRING_LENGTH)           :: filename
	CHARACTER(MAX_STRING_LENGTH)           :: gname
	CHARACTER(MAX_STRING_LENGTH)           :: subgname
	CHARACTER(MAX_STRING_LENGTH)           :: dset
	INTEGER(HID_T)                         :: h5file_id
	INTEGER(HID_T)                         :: group_id
	INTEGER(HID_T)                         :: subgroup_id
	INTEGER                                :: h5error

	filename = TRIM(params%magnetic_field_filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
	end if

	dset = "/R"
	call load_array_from_hdf5(h5file_id,dset,F%X%R)

	if ((.NOT.F%Bflux).AND.(.NOT.F%axisymmetric_fields)) then
		dset = "/PHI"
		call load_array_from_hdf5(h5file_id,dset,F%X%PHI)
	end if

	dset = "/Z"
	call load_array_from_hdf5(h5file_id,dset,F%X%Z)

	dset = '/Bo'
	call load_from_hdf5(h5file_id,dset,F%Bo)

	if (F%Efield) then
		dset = '/Eo'
		call load_from_hdf5(h5file_id,dset,F%Eo)
	else
		F%Eo = 0.0_rp
	end if

	dset = '/Ro'
	call load_from_hdf5(h5file_id,dset,F%Ro)

	dset = '/Zo'
	call load_from_hdf5(h5file_id,dset,F%Zo)

	if (F%Bflux.OR.F%axisymmetric_fields) then
		dset = "/FLAG"
		call load_array_from_hdf5(h5file_id,dset,F%FLAG2D)
	else
		dset = "/FLAG"
		call load_array_from_hdf5(h5file_id,dset,F%FLAG3D)
	end if

	if (F%Bflux) then
		dset = "/PSIp"
		call load_array_from_hdf5(h5file_id,dset,F%PSIp)
	end if

	if (F%Bfield) then
		if (F%axisymmetric_fields) then
			dset = "/BR"
			call load_array_from_hdf5(h5file_id,dset,F%B_2D%R)

			dset = "/BPHI"
			call load_array_from_hdf5(h5file_id,dset,F%B_2D%PHI)

			dset = "/BZ"
			call load_array_from_hdf5(h5file_id,dset,F%B_2D%Z)
		else
			dset = "/BR"
			call load_array_from_hdf5(h5file_id,dset,F%B_3D%R)

			dset = "/BPHI"
			call load_array_from_hdf5(h5file_id,dset,F%B_3D%PHI)

			dset = "/BZ"
			call load_array_from_hdf5(h5file_id,dset,F%B_3D%Z)
		end if
	end if

	if (F%Efield.AND.F%Efield_in_file) then
		if (F%axisymmetric_fields) then
			dset = "/ER"
			call load_array_from_hdf5(h5file_id,dset,F%E_2D%R)

			dset = "/EPHI"
			call load_array_from_hdf5(h5file_id,dset,F%E_2D%PHI)

			dset = "/EZ"
			call load_array_from_hdf5(h5file_id,dset,F%E_2D%Z)
		else
			dset = "/ER"
			call load_array_from_hdf5(h5file_id,dset,F%E_3D%R)

			dset = "/EPHI"
			call load_array_from_hdf5(h5file_id,dset,F%E_3D%PHI)

			dset = "/EZ"
			call load_array_from_hdf5(h5file_id,dset,F%E_3D%Z)
		end if
	end if

	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
	end if
end subroutine load_field_data_from_hdf5


!> @brief Subroutine that allocates the variables keeping the axisymmetric fields data.
!! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
!! @param[in] bfield Logical variable that specifies if the variables that keep the magnetic field data is allocated (bfield=T) or not (bfield=F).
!! @param[in] bflux Logical variable that specifies if the variables that keep the poloidal magnetic flux data is allocated (bflux=T) or not (bflux=F).
!! @param[in] efield Logical variable that specifies if the variables that keep the electric field data is allocated (efield=T) or not (efield=F).
subroutine ALLOCATE_2D_FIELDS_ARRAYS(F,bfield,bflux,efield)
	TYPE(FIELDS), INTENT(INOUT)    :: F
	LOGICAL, INTENT(IN)            :: bfield
	LOGICAL, INTENT(IN)            :: bflux
	LOGICAL, INTENT(IN)            :: efield

	if (bfield) then
		call ALLOCATE_V_FIELD_2D(F%B_2D,F%dims)
	end if

	if (bflux) then
		ALLOCATE(F%PSIp(F%dims(1),F%dims(3)))
	end if

	if (efield) then
		call ALLOCATE_V_FIELD_2D(F%E_2D,F%dims)
	end if

	if (.NOT.ALLOCATED(F%FLAG2D)) ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))

	if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
	if (.NOT.ALLOCATED(F%X%z)) ALLOCATE(F%X%Z(F%dims(3)))
end subroutine ALLOCATE_2D_FIELDS_ARRAYS


!> @brief Subroutine that allocates the variables keeping the 3-D fields data.
!! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
!! @param[in] bfield Logical variable that specifies if the variables that keep the magnetic field data is allocated (bfield=T) or not (bfield=F).
!! @param[in] efield Logical variable that specifies if the variables that keep the electric field data is allocated (efield=T) or not (efield=F).
subroutine ALLOCATE_3D_FIELDS_ARRAYS(F,bfield,efield)
	TYPE(FIELDS), INTENT(INOUT)    :: F
	LOGICAL, INTENT(IN)            :: bfield
	LOGICAL, INTENT(IN)            :: efield

	if (bfield) then
		call ALLOCATE_V_FIELD_3D(F%B_3D,F%dims)
	end if

	if (efield) then
		call ALLOCATE_V_FIELD_3D(F%E_3D,F%dims)
	end if

	if (.NOT.ALLOCATED(F%FLAG3D)) ALLOCATE(F%FLAG3D(F%dims(1),F%dims(2),F%dims(3)))

	if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
	if (.NOT.ALLOCATED(F%X%PHI)) ALLOCATE(F%X%PHI(F%dims(2)))
	if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
end subroutine ALLOCATE_3D_FIELDS_ARRAYS


!> @brief Subroutine that allocates the cylindrical components of an axisymmetric field.
!! @param[in,out] F Vector field to be allocated.
!! @param[in] dims Dimension of the mesh containing the field data.
subroutine ALLOCATE_V_FIELD_2D(F,dims)
	TYPE(V_FIELD_2D), INTENT(INOUT)    :: F
	INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(3)))
end subroutine ALLOCATE_V_FIELD_2D


!> @brief Subroutine that allocates the cylindrical components of a 3-D field.
!! @param[in,out] F Vector field to be allocated.
!! @param[in] dims Dimension of the mesh containing the field data.
subroutine ALLOCATE_V_FIELD_3D(F,dims)
	TYPE(V_FIELD_3D), INTENT(INOUT)    :: F
	INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(2),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(2),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(2),dims(3)))
end subroutine ALLOCATE_V_FIELD_3D

!> @brief Subroutine that deallocates all the variables of the electric and magnetic fields.
!! @param[in,out] F An instance of the KORC derived type FIELDS.
subroutine DEALLOCATE_FIELDS_ARRAYS(F)
	TYPE(FIELDS), INTENT(INOUT) :: F

	if (ALLOCATED(F%PSIp)) DEALLOCATE(F%PSIp)

	if (ALLOCATED(F%B_2D%R)) DEALLOCATE(F%B_2D%R)
	if (ALLOCATED(F%B_2D%PHI)) DEALLOCATE(F%B_2D%PHI)
	if (ALLOCATED(F%B_2D%Z)) DEALLOCATE(F%B_2D%Z)

	if (ALLOCATED(F%B_3D%R)) DEALLOCATE(F%B_3D%R)
	if (ALLOCATED(F%B_3D%PHI)) DEALLOCATE(F%B_3D%PHI)
	if (ALLOCATED(F%B_3D%Z)) DEALLOCATE(F%B_3D%Z)

	if (ALLOCATED(F%E_2D%R)) DEALLOCATE(F%E_2D%R)
	if (ALLOCATED(F%E_2D%PHI)) DEALLOCATE(F%E_2D%PHI)
	if (ALLOCATED(F%E_2D%Z)) DEALLOCATE(F%E_2D%Z)

	if (ALLOCATED(F%E_3D%R)) DEALLOCATE(F%E_3D%R)
	if (ALLOCATED(F%E_3D%PHI)) DEALLOCATE(F%E_3D%PHI)
	if (ALLOCATED(F%E_3D%Z)) DEALLOCATE(F%E_3D%Z)

	if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
	if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
	if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)

	if (ALLOCATED(F%FLAG2D)) DEALLOCATE(F%FLAG2D)
	if (ALLOCATED(F%FLAG3D)) DEALLOCATE(F%FLAG3D)
end subroutine DEALLOCATE_FIELDS_ARRAYS
end module korc_fields
