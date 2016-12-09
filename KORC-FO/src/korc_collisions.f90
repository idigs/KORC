module korc_collisions
	use korc_types
	use korc_constants
	use korc_HDF5

	implicit none

	CHARACTER(LEN=*), PRIVATE, PARAMETER :: MODEL1 = 'SINGLE_SPECIES'
	CHARACTER(LEN=*), PRIVATE, PARAMETER :: MODEL2 = 'MULTIPLE_SPECIES'

	TYPE, PRIVATE :: PARAMS_MS
		INTEGER :: num_impurity_species
		REAL(rp) :: Te ! Background electron temperature in eV
		REAL(rp) :: ne ! Background electron density in 1/m^3
		REAL(rp) :: nH ! Background proton density in 1/m^3
		REAL(rp) :: nef ! Free electron density in 1/m^3
		REAL(rp), DIMENSION(:), ALLOCATABLE :: neb ! Bound electron density in 1/m^3
		REAL(rp), DIMENSION(:), ALLOCATABLE :: Zi ! Atomic number of (majority) background ions
		REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
		REAL(rp), DIMENSION(:), ALLOCATABLE :: Zj ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
		REAL(rp), DIMENSION(:), ALLOCATABLE :: nz ! Impurity densities
		REAL(rp), DIMENSION(:), ALLOCATABLE :: IZj ! Ionization energy of impurity in eV
		REAL(rp), DIMENSION(:), ALLOCATABLE :: Ee_IZj ! me*c^2/IZj dimensionless parameter

		REAL(rp) :: rD ! Debye length
		REAL(rp) :: re ! Classical electron radius
	END TYPE PARAMS_MS

	TYPE, PRIVATE :: PARAMS_SS
		REAL(rp) :: Te ! Electron temperature
		REAL(rp) :: Ti ! Ion temperature
		REAL(rp) :: ne ! Background electron density
		REAL(rp) :: Zeff ! Effective atomic number of ions
		REAL(rp) :: rD ! Debye radius
		REAL(rp) :: re ! Classical electron radius
		REAL(rp) :: CoulombLog ! Coulomb logarithm
		REAL(rp) :: VTe ! Thermal velocity of background electrons
		REAL(rp) :: delta ! delta parameter
		REAL(rp) :: Gammac ! Gamma factor
		REAL(rp) :: Tau ! Collisional time
		REAL(rp) :: ED ! Dreicer electric field
	END TYPE PARAMS_SS

	TYPE(PARAMS_MS), PRIVATE :: cparams_ms
	TYPE(PARAMS_SS), PRIVATE :: cparams_ss

	PUBLIC :: initialize_collision_params,normalize_collisions_params,&
				collision_force,deallocate_collisions_params,save_collision_params
	PRIVATE :: load_params_ms,normalize_params_ms,deallocate_params_ms

	contains

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * SUBROUTINES FOR INITIALIZING COLLISIONS PARAMS * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !


subroutine load_params_ms(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: Te ! Background electron temperature in eV
	REAL(rp) :: ne! Background electron density in 1/m^3
	INTEGER :: num_impurity_species
	REAL(rp), DIMENSION(10) :: Zo ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(10) :: Zj ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(10) :: nz ! Impurity densities
	REAL(rp), DIMENSION(10) :: IZj ! Ionization energy of impurity in eV

	NAMELIST /CollisionParamsMultipleSpecies/ num_impurity_species,Te,ne,Zo,Zj,nz,IZj


	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=CollisionParamsMultipleSpecies)
	close(default_unit_open)

!	write(*,nml=CollisionParamsMultipleSpecies)

	cparams_ms%num_impurity_species = num_impurity_species

	ALLOCATE(cparams_ms%Zj(cparams_ms%num_impurity_species))
	ALLOCATE(cparams_ms%Zo(cparams_ms%num_impurity_species))
	ALLOCATE(cparams_ms%nz(cparams_ms%num_impurity_species))
	ALLOCATE(cparams_ms%neb(cparams_ms%num_impurity_species))
	ALLOCATE(cparams_ms%IZj(cparams_ms%num_impurity_species))
	ALLOCATE(cparams_ms%Ee_IZj(cparams_ms%num_impurity_species))

	cparams_ms%Te = Te*C_E
	cparams_ms%ne = ne
	cparams_ms%nH = ne

	cparams_ms%Zj = Zj(1:cparams_ms%num_impurity_species)
	cparams_ms%Zo = Zo(1:cparams_ms%num_impurity_species)
	cparams_ms%nz = nz(1:cparams_ms%num_impurity_species)
	cparams_ms%IZj = C_E*IZj(1:cparams_ms%num_impurity_species)

	cparams_ms%nef = ne + sum(cparams_ms%Zj*cparams_ms%nz)
	cparams_ms%neb = (cparams_ms%Zo-cparams_ms%Zj)*cparams_ms%nz

	cparams_ms%rD = SQRT( C_E0*cparams_ms%Te/(cparams_ms%ne*C_E**2) )
	cparams_ms%re = C_E**2/( 4.0_rp*C_PI*C_E0*C_ME*C_C**2 )
	cparams_ms%Ee_IZj = C_ME*C_C**2/cparams_ms%IZj
end subroutine load_params_ms


subroutine load_params_ss(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: Te ! Electron temperature
	REAL(rp) :: Ti ! Ion temperature
	REAL(rp) :: ne ! Background electron density
	REAL(rp) :: Zeff ! Effective atomic number of ions

	NAMELIST /CollisionParamsSingleSpecies/ Te, Ti, ne, Zeff


	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=CollisionParamsSingleSpecies)
	close(default_unit_open)

	write(*,nml=CollisionParamsSingleSpecies)

	cparams_ss%Te = Te*C_E;
	cparams_ss%Ti = Ti*C_E;
	cparams_ss%ne = ne;
	cparams_ss%Zeff = Zeff;

	cparams_ss%rD = &
	SQRT(C_E0*cparams_ss%Te/(cparams_ss%ne*C_E**2*(cparams_ss%Te/cparams_ss%Ti)))

	cparams_ss%re = C_E**2/(4.0_rp*C_PI*C_E0*C_ME*C_C**2)
	cparams_ss%CoulombLog = 25.3_rp - 1.15_rp*LOG10(1E-3_rp*cparams_ss%ne) +&
	2.3_rp*LOG10(cparams_ss%Te/C_E)

	cparams_ss%VTe = SQRT(2.0_rp*cparams_ss%Te/C_ME)
	cparams_ss%delta = cparams_ss%VTe/C_C
	cparams_ss%Gammac = &
	cparams_ss%ne*C_E**4*cparams_ss%CoulombLog/(4.0_rp*C_PI*C_E0**2)

	cparams_ss%Tau = C_ME**2*C_C**3/cparams_ss%Gammac
	cparams_ss%ED = &
	cparams_ss%ne*C_E**3*cparams_ss%CoulombLog/(4.0_rp*C_PI*C_E0**2*cparams_ss%Te)
end subroutine load_params_ss


subroutine initialize_collision_params(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
				call load_params_ss(params)
			CASE (MODEL2)
				call load_params_ms(params)
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine initialize_collision_params


subroutine normalize_params_ms(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	cparams_ms%Te = cparams_ms%Te/params%cpp%temperature
	cparams_ms%ne = cparams_ms%ne/params%cpp%density
	cparams_ms%nH = cparams_ms%nH/params%cpp%density
	cparams_ms%nef = cparams_ms%nef/params%cpp%density
	cparams_ms%neb = cparams_ms%neb/params%cpp%density
	if (ALLOCATED(cparams_ms%nz)) cparams_ms%nz = cparams_ms%nz/params%cpp%density
	if (ALLOCATED(cparams_ms%IZj)) cparams_ms%IZj = cparams_ms%IZj/params%cpp%energy
	cparams_ms%rD = cparams_ms%rD/params%cpp%length
	cparams_ms%re = cparams_ms%re/params%cpp%length
end subroutine normalize_params_ms


subroutine normalize_collisions_params(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
				write(6,'("Something to be done")')
			CASE (MODEL2)
				call normalize_params_ms(params)
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine normalize_collisions_params


subroutine collision_force(spp,U,Fcoll)
    implicit none
	TYPE(SPECIES), INTENT(IN) :: spp
	REAL(rp), DIMENSION(3), INTENT(IN) :: U
	REAL(rp), DIMENSION(3), INTENT(OUT) :: Fcoll
	REAL(rp), DIMENSION(3) :: V, Fcolle, Fcolli
	REAL(rp) :: gamma, tmp
	REAL(rp) :: ae, ai, Clog_ef, Clog_eb, Clog_eH, Clog_eZj, Clog_eZo
	INTEGER :: ppi
	
	gamma = sqrt(1.0_rp + DOT_PRODUCT(U,U))
	V = U/gamma
	
	tmp = (gamma - 1.0_rp)*sqrt(gamma + 1.0_rp)
	Clog_ef = log(0.5_rp*tmp*(cparams_ms%rD/cparams_ms%re)/gamma)
	ae = cparams_ms%nef*Clog_ef
	do ppi=1,cparams_ms%num_impurity_species
		Clog_eb = log(tmp*cparams_ms%Ee_IZj(ppi))
		ae = ae + cparams_ms%neb(ppi)*Clog_eb
	end do

	tmp = (gamma**2 - 1.0_rp)/gamma
	Clog_eH = log( tmp*(cparams_ms%rD/cparams_ms%re) )
	ai = cparams_ms%nH*Clog_eH
	do ppi=1,cparams_ms%num_impurity_species
		Clog_eZj = log( cparams_ms%rD/(cparams_ms%Zj(ppi)*cparams_ms%re*cparams_ms%Ee_IZj(ppi)) )
		Clog_eZo = log(tmp*cparams_ms%Ee_IZj(ppi))
		ai = ai + &
			cparams_ms%nz(ppi)*(Clog_eZj*cparams_ms%Zj(ppi)**2 + Clog_eZo*cparams_ms%Zo(ppi)**2)
	end do

	tmp = gamma*(gamma + 1.0_rp)/(sqrt(DOT_PRODUCT(U,U))**3)
	Fcolle = -4.0_rp*C_PI*ae*spp%m*(cparams_ms%re**2)*tmp*U

	tmp = gamma/(sqrt(DOT_PRODUCT(U,U))**3)
	Fcolli = -4.0_rp*C_PI*ai*spp%m*(cparams_ms%re**2)*tmp*U

	Fcoll = Fcolle + Fcolli
end subroutine collision_force


subroutine save_params_ms(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER :: h5error
	REAL(rp) :: units

	if (params%mpi_params%rank .EQ. 0) then
		filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"
		call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

		gname = "collision_params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		ALLOCATE(attr_array(cparams_ms%num_impurity_species))

		dset = TRIM(gname) // "/num_impurity_species"
		attr = "Number of impurity species"
		call save_to_hdf5(h5file_id,dset,cparams_ms%num_impurity_species,attr)

		dset = TRIM(gname) // "/Te"
		attr = "Background electron temperature in eV"
		units = params%cpp%temperature/C_E
		call save_to_hdf5(h5file_id,dset,units*cparams_ms%Te,attr)

		dset = TRIM(gname) // "/ne"
		attr = "Background electron density in m^-3"
		units = params%cpp%density
		call save_to_hdf5(h5file_id,dset,units*cparams_ms%ne,attr)

		dset = TRIM(gname) // "/nH"
		attr = "Background proton density in m^-3"
		units = params%cpp%density
		call save_to_hdf5(h5file_id,dset,units*cparams_ms%nH,attr)

		dset = TRIM(gname) // "/nef"
		attr = "Free electron density in m^-3"
		units = params%cpp%density
		call save_to_hdf5(h5file_id,dset,units*cparams_ms%nef,attr)

		dset = TRIM(gname) // "/neb"
		attr_array(1) = "Bound electron density per impurity in m^-3"
		units = params%cpp%density
		call save_1d_array_to_hdf5(h5file_id,dset,units*cparams_ms%neb,attr_array)

		dset = TRIM(gname) // "/Zo"
		attr_array(1) = "Full nuclear charge of impurities"
		call save_1d_array_to_hdf5(h5file_id,dset,cparams_ms%Zo,attr_array)

		dset = TRIM(gname) // "/Zj"
		attr_array(1) = "Average charge state of impurities"
		call save_1d_array_to_hdf5(h5file_id,dset,cparams_ms%Zj,attr_array)

		dset = TRIM(gname) // "/nz"
		attr_array(1) = "Density of impurities in m^-3"
		units = params%cpp%density
		call save_1d_array_to_hdf5(h5file_id,dset,units*cparams_ms%nz,attr_array)

		dset = TRIM(gname) // "/IZj"
		attr_array(1) = " Ionization energy of impurities in eV"
		units = params%cpp%energy/C_E
		call save_1d_array_to_hdf5(h5file_id,dset,units*cparams_ms%IZj,attr_array)

		dset = TRIM(gname) // "/rD"
		attr = "Debye length in m"
		units = params%cpp%length
		call save_to_hdf5(h5file_id,dset,units*cparams_ms%rD,attr)

		dset = TRIM(gname) // "/re"
		attr = "Classical electron radius in m"
		units = params%cpp%length
		call save_to_hdf5(h5file_id,dset,units*cparams_ms%re,attr)

		DEALLOCATE(attr_array)	
	
		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if
	
end subroutine save_params_ms


subroutine save_collision_params(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
				write(6,'("Something to be done")')
			CASE (MODEL2)
				call save_params_ms(params)
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine save_collision_params


subroutine deallocate_params_ms()
	implicit none

	if (ALLOCATED(cparams_ms%Zj)) DEALLOCATE(cparams_ms%Zj)
	if (ALLOCATED(cparams_ms%Zo)) DEALLOCATE(cparams_ms%Zo)
	if (ALLOCATED(cparams_ms%nz)) DEALLOCATE(cparams_ms%nz)
	if (ALLOCATED(cparams_ms%neb)) DEALLOCATE(cparams_ms%neb)
	if (ALLOCATED(cparams_ms%IZj)) DEALLOCATE(cparams_ms%IZj)
	if (ALLOCATED(cparams_ms%Zj)) DEALLOCATE(cparams_ms%Ee_IZj)
end subroutine deallocate_params_ms


subroutine deallocate_collisions_params(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
				write(6,'("Something to be done")')
			CASE (MODEL2)
				call deallocate_params_ms()
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine deallocate_collisions_params

end module korc_collisions
