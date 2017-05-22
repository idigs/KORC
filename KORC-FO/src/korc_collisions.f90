module korc_collisions
	use korc_types
	use korc_constants
	use korc_HDF5

	IMPLICIT NONE

	CHARACTER(LEN=*), PRIVATE, PARAMETER :: MODEL1 = 'SINGLE_SPECIES'
	CHARACTER(LEN=*), PRIVATE, PARAMETER :: MODEL2 = 'MULTIPLE_SPECIES'
    REAL(rp), PRIVATE, PARAMETER :: infinity = HUGE(1.0_rp)

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
		REAL(rp) :: Ec ! Critical electric field
		REAL(rp) :: ED ! Dreicer electric field
		REAL(rp) :: dTau ! Subcycling time step in collisional time units (Tau)
        INTEGER(ip) :: subcycling_iterations
	END TYPE PARAMS_SS

	TYPE(PARAMS_MS), PRIVATE :: cparams_ms
	TYPE(PARAMS_SS), PRIVATE :: cparams_ss

	PUBLIC :: initialize_collision_params,normalize_collisions_params,&
				collision_force,deallocate_collisions_params,save_collision_params
	PRIVATE :: load_params_ms,load_params_ss,normalize_params_ms,&
				normalize_params_ss,save_params_ms,save_params_ss,&
				deallocate_params_ms,cross,unitVectors,CA,CB,CF,fun

	contains

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * SUBROUTINES FOR INITIALIZING COLLISIONS PARAMS * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !


subroutine load_params_ms(params)
	IMPLICIT NONE
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
	cparams_ms%re = C_RE
	cparams_ms%Ee_IZj = C_ME*C_C**2/cparams_ms%IZj
end subroutine load_params_ms


subroutine load_params_ss(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: Te ! Electron temperature
	REAL(rp) :: Ti ! Ion temperature
	REAL(rp) :: ne ! Background electron density
	REAL(rp) :: Zeff ! Effective atomic number of ions
	REAL(rp) :: dTau ! Subcycling time step in collisional time units (Tau)

	NAMELIST /CollisionParamsSingleSpecies/ Te, Ti, ne, Zeff, dTau


	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=CollisionParamsSingleSpecies)
	close(default_unit_open)

	cparams_ss%Te = Te*C_E
	cparams_ss%Ti = Ti*C_E
	cparams_ss%ne = ne
	cparams_ss%Zeff = Zeff
	cparams_ss%dTau = dTau

	cparams_ss%rD = &
	SQRT(C_E0*cparams_ss%Te/(cparams_ss%ne*C_E**2*(1.0_rp + cparams_ss%Te/cparams_ss%Ti)))

	cparams_ss%re = C_E**2/(4.0_rp*C_PI*C_E0*C_ME*C_C**2)
	cparams_ss%CoulombLog = 25.3_rp - 1.15_rp*LOG10(1E-6_rp*cparams_ss%ne) + 2.3_rp*LOG10(cparams_ss%Te/C_E)

	cparams_ss%VTe = SQRT(2.0_rp*cparams_ss%Te/C_ME)
	cparams_ss%delta = cparams_ss%VTe/C_C
	cparams_ss%Gammac = cparams_ss%ne*C_E**4*cparams_ss%CoulombLog/(4.0_rp*C_PI*C_E0**2)

	cparams_ss%Tau = C_ME**2*C_C**3/cparams_ss%Gammac
	cparams_ss%Ec = cparams_ss%ne*C_E**3*cparams_ss%CoulombLog/(4.0_rp*C_PI*C_E0**2*C_ME*C_C**2)
	cparams_ss%ED = cparams_ss%ne*C_E**3*cparams_ss%CoulombLog/(4.0_rp*C_PI*C_E0**2*cparams_ss%Te)
end subroutine load_params_ss


subroutine initialize_collision_params(params)
	IMPLICIT NONE
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
	IMPLICIT NONE
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


subroutine normalize_params_ss(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	cparams_ss%Te = cparams_ss%Te/params%cpp%temperature
	cparams_ss%Ti = cparams_ss%Ti/params%cpp%temperature
	cparams_ss%ne = cparams_ss%ne/params%cpp%density
	cparams_ss%rD = cparams_ss%rD/params%cpp%length
	cparams_ss%re = cparams_ss%re/params%cpp%length
	cparams_ss%VTe = cparams_ss%VTe/params%cpp%velocity
	cparams_ss%Gammac = &
	cparams_ss%Gammac*params%cpp%time/(params%cpp%mass**2*params%cpp%velocity**3)
	cparams_ss%Tau = cparams_ss%Tau/params%cpp%time
	cparams_ss%Ec = cparams_ss%Ec/params%cpp%Eo
	cparams_ss%ED = cparams_ss%ED/params%cpp%Eo
end subroutine normalize_params_ss


subroutine normalize_collisions_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
				call normalize_params_ss(params)
			CASE (MODEL2)
				call normalize_params_ms(params)
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine normalize_collisions_params


subroutine collision_force(spp,U,Fcoll)
!				if (params%collisions .AND. &
!				(TRIM(params%collisions_model) .EQ. 'MULTIPLE_SPECIES')) then
!					call collision_force(spp(ii),U_os,Fcoll)
!					U_RC = U_RC + a*Fcoll/spp(ii)%q
!				end if
    IMPLICIT NONE
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
	do ppi=1_idef,cparams_ms%num_impurity_species
		Clog_eb = log(tmp*cparams_ms%Ee_IZj(ppi))
		ae = ae + cparams_ms%neb(ppi)*Clog_eb
	end do

	tmp = (gamma**2 - 1.0_rp)/gamma
	Clog_eH = log( tmp*(cparams_ms%rD/cparams_ms%re) )
	ai = cparams_ms%nH*Clog_eH
	do ppi=1_idef,cparams_ms%num_impurity_species
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


subroutine define_collisions_time_step(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	cparams_ss%subcycling_iterations = FLOOR((cparams_ss%dTau*cparams_ss%Tau)/params%dt,ip)

	if (params%collisions .AND. (params%mpi_params%rank .EQ. 0)) then
		write(6,'(/,"* * * * * * * SUBCYCLING FOR COLLISIONS * * * * * * *")')
		write(6,'("Number of KORC iterations per collision: ",I16)') cparams_ss%subcycling_iterations
		write(6,'("* * * * * * * * * * * * * * * * * * * * * * * * * * *",/)')
	end if
end subroutine define_collisions_time_step


! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !
! * FUNCTIONS OF COLLISION OPERATOR FOR SINGLE-SPECIES PLASMAS * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !
function psi(v)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: v
	REAL(rp) :: psi
	REAL(rp) :: x

	x = v/cparams_ss%VTe
	psi = 0.5_rp*(ERF(x) - 2.0_rp*x*EXP(-x**2)/SQRT(C_PI))/x**2
end function psi


function CA(v)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: v
	REAL(rp) :: CA

	CA  = cparams_ss%Gammac*psi(v)/v
end function CA


function CF(v)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: v
	REAL(rp) :: CF

	CF  = cparams_ss%Gammac*psi(v)/cparams_ss%Te
end function CF


function CB(v)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: v
	REAL(rp) :: CB
	REAL(rp) :: x

	x = v/cparams_ss%VTe
	CB  = (0.5_rp*cparams_ss%Gammac/v)*( cparams_ss%Zeff + ERF(x) -&
			psi(v) + 0.5_rp*cparams_ss%delta**4*x**2 )
end function CB


function fun(v)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: v
	REAL(rp) :: fun
	REAL(rp) :: x

	x = v/cparams_ss%VTe
	fun = 2.0_rp*( 1.0_rp/x + x )*EXP(-x**2)/SQRT(C_PI) - ERF(x)/x**2 - psi(v)
end function fun


function cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


subroutine unitVectors(B,b1,b2,b3)
    IMPLICIT NONE
    REAL(rp), DIMENSION(3), INTENT(IN) :: B
	REAL(rp), DIMENSION(3), INTENT(OUT) :: b1
	REAL(rp), DIMENSION(3), INTENT(OUT) :: b2
	REAL(rp), DIMENSION(3), INTENT(OUT) :: b3

    b1 = B/SQRT(DOT_PRODUCT(B,B))

    b2 = cross(b1,(/0.0_rp,0.0_rp,1.0_rp/))
    b2 = b2/SQRT(DOT_PRODUCT(b2,b2))

    b3 = cross(b1,b2)
    b3 = b3/SQRT(DOT_PRODUCT(b3,b3))
end subroutine unitVectors
! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !
! * FUNCTIONS OF COLLISION OPERATOR FOR SINGLE-SPECIES PLASMAS * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !


subroutine include_collisions(params,B,U)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(3), INTENT(IN) :: B
	REAL(rp), DIMENSION(3), INTENT(INOUT) :: U
	REAL(rp), DIMENSION(3) :: b1, b2, b3
	REAL(rp), DIMENSION(3) :: x = (/1.0_rp,0.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: y = (/0.0_rp,1.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: z = (/0.0_rp,0.0_rp,1.0_rp/)
	REAL(rp) :: dt
    REAL(rp) :: rnd_num1, rnd_num2
	REAL(rp) :: U1, U2, U3
	REAL(rp) :: dU1, dU2, dU3
	REAL(rp) :: dp, dxi, dphi
	REAL(rp) :: dWp, dWphi, dWxi
	REAL(rp) :: pitch
	REAL(rp) :: xi ! cosine of pitch angle
	REAL(rp) :: g ! Gamma factor
	REAL(rp) :: v ! Velocity magnitude
	REAL(rp) :: p ! Magnitude of dimensionless momentum U

	if (MODULO(params%it+1_ip,cparams_ss%subcycling_iterations) .EQ. 0_ip) then
!		write(6,'("Collision at iteration: ",I16)') params%it
		dt = REAL(cparams_ss%subcycling_iterations,rp)*params%dt

		p = SQRT(DOT_PRODUCT(U,U))
		g = SQRT(1.0_rp + DOT_PRODUCT(U,U))
		v = p/g

		call unitVectors(B,b1,b2,b3)

		U1 = DOT_PRODUCT(U,b1);
		U2 = DOT_PRODUCT(U,b2);
		U3 = DOT_PRODUCT(U,b3);

		xi = U1/(g*v)

		call RANDOM_NUMBER(rnd_num1)
		call RANDOM_NUMBER(rnd_num2)
		dWp =  &
		SQRT(dt)*SQRT(-2.0_rp*LOG(1.0_rp-rnd_num1))*COS(2.0_rp*C_PI*rnd_num2);

		call RANDOM_NUMBER(rnd_num1)
		dWphi = 2.0_rp*C_PI*rnd_num1*SQRT(dt);

		call RANDOM_NUMBER(rnd_num1)
		dWxi = rnd_num1*SQRT(dt);


		dp = ( -CF(v) + 2.0_rp*CA(v)/p + cparams_ss%Gammac*fun(v)/(p**2*SQRT(1.0_rp + p**2)) )*dt&
			 + SQRT(2.0_rp*CA(v))*dWp

		dxi = -2.0_rp*xi*CB(v)*dt/p**2 - SQRT(2.0_rp*CB(v)*(1.0_rp - xi**2))*dWxi/p

		dphi = SQRT(2.0_rp*CB(v))*dWphi/(p*SQRT(1.0_rp - xi**2))
		if ((dphi .GE. infinity) .OR. ISNAN(dphi)) then
		    dphi = 0.0_rp
		end if

		pitch = ACOS(xi + dxi);

		dU1 = dp*COS(pitch);
		dU2 = dp*SIN(pitch)*COS(dphi);
		dU3 = dp*SIN(pitch)*SIN(dphi);

!		write(6,'("xi,pitch,dU1,dU2,dU3: ",6F25.16)') U(1),U(2),U(3),U1,U2,U3
		
		U(1) = (U1+dU1)*DOT_PRODUCT(b1,x) + (U2+dU2)*DOT_PRODUCT(b2,x)&
			+ (U3+dU3)*DOT_PRODUCT(b3,x)
		U(2) = (U1+dU1)*DOT_PRODUCT(b1,y) + (U2+dU2)*DOT_PRODUCT(b2,y)&
			+ (U3+dU3)*DOT_PRODUCT(b3,y)
		U(3) = (U1+dU1)*DOT_PRODUCT(b1,z) + (U2+dU2)*DOT_PRODUCT(b2,z)&
			+ (U3+dU3)*DOT_PRODUCT(b3,z)
	end if
end subroutine include_collisions


subroutine save_params_ms(params)
	IMPLICIT NONE
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


subroutine save_params_ss(params)
	IMPLICIT NONE
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

		dset = TRIM(gname) // "/Te"
		attr = "Background electron temperature in eV"
		units = params%cpp%temperature/C_E
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%Te,attr)

		dset = TRIM(gname) // "/Ti"
		attr = "Background ion temperature in eV"
		units = params%cpp%temperature/C_E
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%Ti,attr)

		dset = TRIM(gname) // "/ne"
		attr = "Background electron density in m^-3"
		units = params%cpp%density
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%ne,attr)

		dset = TRIM(gname) // "/Zeff"
		attr = "Effective nuclear charge of impurities"
		call save_to_hdf5(h5file_id,dset,cparams_ss%Zeff,attr)

		dset = TRIM(gname) // "/rD"
		attr = "Debye length in m"
		units = params%cpp%length
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%rD,attr)

		dset = TRIM(gname) // "/re"
		attr = "Classical electron radius in m"
		units = params%cpp%length
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%re,attr)

		dset = TRIM(gname) // "/Clog"
		attr = "Coulomb logarithm"
		call save_to_hdf5(h5file_id,dset,cparams_ss%CoulombLog,attr)

		dset = TRIM(gname) // "/VTe"
		attr = "Background electron temperature"
		units = params%cpp%velocity
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%VTe,attr)

		dset = TRIM(gname) // "/delta"
		attr = "Delta parameter VTe/C"
		call save_to_hdf5(h5file_id,dset,cparams_ss%delta,attr)

		dset = TRIM(gname) // "/Gamma"
		attr = "Gamma coefficient"
		units = (params%cpp%mass**2*params%cpp%velocity**3)/params%cpp%time
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%Gammac,attr)

		dset = TRIM(gname) // "/Tau"
		attr = "Collisional time in s"
		units = params%cpp%time
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%Tau,attr)

		dset = TRIM(gname) // "/dTau"
		attr = "Subcycling time step in s"
		units = params%cpp%time
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%dTau*cparams_ss%Tau,attr)

		dset = TRIM(gname) // "/subcycling_iterations"
		attr = "KORC iterations per collision"
		call save_to_hdf5(h5file_id,dset,cparams_ss%subcycling_iterations,attr)

		dset = TRIM(gname) // "/Ec"
		attr = "Critical electric field"
		units = params%cpp%Eo
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%Ec,attr)

		dset = TRIM(gname) // "/ED"
		attr = "Dreicer electric field"
		units = params%cpp%Eo
		call save_to_hdf5(h5file_id,dset,units*cparams_ss%ED,attr)

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if
end subroutine save_params_ss


subroutine save_collision_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
				call save_params_ss(params)
			CASE (MODEL2)
				call save_params_ms(params)
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine save_collision_params


subroutine deallocate_params_ms()
	IMPLICIT NONE

	if (ALLOCATED(cparams_ms%Zj)) DEALLOCATE(cparams_ms%Zj)
	if (ALLOCATED(cparams_ms%Zo)) DEALLOCATE(cparams_ms%Zo)
	if (ALLOCATED(cparams_ms%nz)) DEALLOCATE(cparams_ms%nz)
	if (ALLOCATED(cparams_ms%neb)) DEALLOCATE(cparams_ms%neb)
	if (ALLOCATED(cparams_ms%IZj)) DEALLOCATE(cparams_ms%IZj)
	if (ALLOCATED(cparams_ms%Zj)) DEALLOCATE(cparams_ms%Ee_IZj)
end subroutine deallocate_params_ms


subroutine deallocate_collisions_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	if (params%collisions) then
		SELECT CASE (TRIM(params%collisions_model))
			CASE (MODEL1)
!				write(6,'("Something to be done")')
			CASE (MODEL2)
				call deallocate_params_ms()
			CASE DEFAULT
				write(6,'("Default case")')
		END SELECT
	end if
end subroutine deallocate_collisions_params

end module korc_collisions
