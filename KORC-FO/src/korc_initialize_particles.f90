MODULE korc_initialize_particles
	USE korc_types
	IMPLICIT NONE

	TYPE, PRIVATE :: AVALANCHE_PDF_PARAMS
		REAL(rp) :: max_pitch_angle ! Maximum pitch angle of sampled PDF in degrees
		REAL(rp) :: min_energy ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_energy ! Maximum energy of sampled PDF in MeV
		REAL(rp) :: min_p ! Minimum momentum of sampled PDF
		REAL(rp) :: max_p ! Maximum momentum of sampled PDF
		REAL(rp) :: ne ! Background electron density in m^-3
		REAL(rp) :: Zeff ! Effective atomic number of ions
		REAL(rp) :: Ec ! Critical electric field in V/m
		REAL(rp) :: Epar ! Parallel electric field in V/m
		REAL(rp) :: Ebar ! Epar/Ec
		REAL(rp) :: Te ! Background electron temperature in eV
		REAL(rp) :: lD ! Debye length
		REAL(rp) :: bmin ! Maximum approach radius
		REAL(rp) :: CoulombLog ! Coulomb Logarithm
		REAL(rp) :: Tau ! Collisional time
	END TYPE AVALANCHE_PDF_PARAMS

	TYPE(AVALANCHE_PDF_PARAMS), PRIVATE :: avalanche_params

	CONTAINS

SUBROUTINE initialize_avalanche_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: max_pitch_angle
	REAL(rp) :: max_energy
	REAL(rp) :: ne
	REAL(rp) :: Zeff
	REAL(rp) :: Ec
	REAL(rp) :: Epar
	REAL(rp) :: Te
	NAMELIST /AvalancheGenerationPDF/ max_pitch_angle,max_energy,&
	ne,Zeff,Ec,Epar,Te

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=AvalancheGenerationPDF)
	close(default_unit_open)

!	write(*,nml=AvalancheGenerationPDF)

	avalanche_params%max_pitch_angle = max_pitch_angle
!	avalanche_params%min_energy 
	avalanche_params%max_energy = max_energy
!	avalanche_params%min_p 
!	avalanche_params%max_p 
	avalanche_params%ne = ne
	avalanche_params%Zeff = Zeff
	avalanche_params%Ec = Ec
	avalanche_params%Epar = Epar
	avalanche_params%Ebar = avalanche_params%Epar/avalanche_params%Ec
	avalanche_params%Te = Te
!	avalanche_params%lD
!	avalanche_params%bmin
!	avalanche_params%CoulombLog
!	avalanche_params%Tau
END SUBROUTINE

END MODULE korc_initialize_particles
