module korc_collisions

	use korc_types
	use korc_constants

	implicit none

	PUBLIC :: initialize_collision_params

	TYPE, PUBLIC :: COLLISION_PARAMS
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
	END TYPE COLLISION_PARAMS

	contains

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * SUBROUTINES FOR INITIALIZING COLLISIONS PARAMS * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_collision_params(params,cparams)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(COLLISION_PARAMS), INTENT(OUT) :: cparams
	REAL(rp) :: Te ! Background electron temperature in eV
	REAL(rp) :: ne! Background electron density in 1/m^3
	INTEGER :: num_impurity_species
	REAL(rp), DIMENSION(10) :: Zo ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(10) :: Zj ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
	REAL(rp), DIMENSION(10) :: nz ! Impurity densities
	REAL(rp), DIMENSION(10) :: IZj ! Ionization energy of impurity in eV

	NAMELIST /CollisionParamsMultipleSpecies/ num_impurity_species,Te,ne,Zo,Zj,nz,IZj


	if (params%collisions) then
		open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
		read(default_unit_open,nml=CollisionParamsMultipleSpecies)
		close(default_unit_open)

	!	write(*,nml=CollisionParamsMultipleSpecies)

		cparams%num_impurity_species = num_impurity_species

		ALLOCATE(cparams%Zj(cparams%num_impurity_species))
		ALLOCATE(cparams%Zo(cparams%num_impurity_species))
		ALLOCATE(cparams%nz(cparams%num_impurity_species))
		ALLOCATE(cparams%neb(cparams%num_impurity_species))
		ALLOCATE(cparams%IZj(cparams%num_impurity_species))
		ALLOCATE(cparams%Ee_IZj(cparams%num_impurity_species))

		cparams%Te = Te*C_E
		cparams%ne = ne
		cparams%nH = ne

		cparams%Zj = Zj(1:cparams%num_impurity_species)
		cparams%Zo = Zo(1:cparams%num_impurity_species)
		cparams%nz = nz(1:cparams%num_impurity_species)
		cparams%IZj = C_E*IZj(1:cparams%num_impurity_species)

		cparams%nef = ne + sum(cparams%Zj*cparams%nz)
		cparams%neb = (cparams%Zo-cparams%Zj)*cparams%nz

		cparams%rD = SQRT( C_E0*cparams%Te/(cparams%ne*C_E**2) )
		cparams%re = C_E**2/( 4.0_rp*C_PI*C_E0*C_ME*C_C**2 )
		cparams%Ee_IZj = C_ME*C_C**2/cparams%IZj
	end if
end subroutine initialize_collision_params

subroutine normalize_collisions_params(params,cparams)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(COLLISION_PARAMS), INTENT(OUT) :: cparams	
	if (params%collisions) then
		cparams%Te = cparams%Te/params%cpp%temperature
		cparams%ne = cparams%ne/params%cpp%density
		cparams%nH = cparams%nH/params%cpp%density
		cparams%nef = cparams%nef/params%cpp%density
		cparams%neb = cparams%neb/params%cpp%density
		if (ALLOCATED(cparams%nz)) cparams%nz = cparams%nz/params%cpp%density
		if (ALLOCATED(cparams%IZj)) cparams%IZj = cparams%IZj/params%cpp%energy
		cparams%rD = cparams%rD/params%cpp%length
		cparams%re = cparams%re/params%cpp%length
	end if
end subroutine normalize_collisions_params


subroutine collision_force(spp,cparams,U,Fcoll)
    implicit none
	TYPE(SPECIES), INTENT(IN) :: spp
	TYPE(COLLISION_PARAMS), INTENT(IN) :: cparams
	REAL(rp), DIMENSION(3), INTENT(IN) :: U
	REAL(rp), DIMENSION(3), INTENT(OUT) :: Fcoll
	REAL(rp), DIMENSION(3) :: V, Fcolle, Fcolli
	REAL(rp) :: gamma, tmp
	REAL(rp) :: ae, ai, Clog_ef, Clog_eb, Clog_eH, Clog_eZj, Clog_eZo
	INTEGER :: ppi
	
	gamma = sqrt(1.0_rp + DOT_PRODUCT(U,U))
	V = U/gamma
	
	tmp = (gamma - 1.0_rp)*sqrt(gamma + 1.0_rp)
	Clog_ef = log(0.5_rp*tmp*(cparams%rD/cparams%re)/gamma)
	ae = cparams%nef*Clog_ef
	do ppi=1,cparams%num_impurity_species
		Clog_eb = log(tmp*cparams%Ee_IZj(ppi))
		ae = ae + cparams%neb(ppi)*Clog_eb
	end do

	tmp = (gamma**2 - 1.0_rp)/gamma
	Clog_eH = log( tmp*(cparams%rD/cparams%re) )
	ai = cparams%nH*Clog_eH
	do ppi=1,cparams%num_impurity_species
		Clog_eZj = log( cparams%rD/(cparams%Zj(ppi)*cparams%re*cparams%Ee_IZj(ppi)) )
		Clog_eZo = log(tmp*cparams%Ee_IZj(ppi))
		ai = ai + &
			cparams%nz(ppi)*(Clog_eZj*cparams%Zj(ppi)**2 + Clog_eZo*cparams%Zo(ppi)**2)
	end do

	tmp = gamma*(gamma + 1.0_rp)/(sqrt(DOT_PRODUCT(U,U))**3)
	Fcolle = -4.0_rp*C_PI*ae*spp%m*(cparams%re**2)*tmp*U

	tmp = gamma/(sqrt(DOT_PRODUCT(U,U))**3)
	Fcolli = -4.0_rp*C_PI*ai*spp%m*(cparams%re**2)*tmp*U

	Fcoll = Fcolle + Fcolli
end subroutine collision_force


subroutine deallocate_collisions_params(cparams)
	implicit none
	TYPE(COLLISION_PARAMS), INTENT(OUT) :: cparams	

	if (ALLOCATED(cparams%Zj)) DEALLOCATE(cparams%Zj)
	if (ALLOCATED(cparams%Zo)) DEALLOCATE(cparams%Zo)
	if (ALLOCATED(cparams%nz)) DEALLOCATE(cparams%nz)
	if (ALLOCATED(cparams%neb)) DEALLOCATE(cparams%neb)
	if (ALLOCATED(cparams%IZj)) DEALLOCATE(cparams%IZj)
	if (ALLOCATED(cparams%Zj)) DEALLOCATE(cparams%Ee_IZj)
end subroutine deallocate_collisions_params

end module korc_collisions
