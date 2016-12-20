MODULE korc_avalanche
	USE korc_types
	USE korc_constants
	USE korc_HDF5
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

	TYPE(AVALANCHE_PDF_PARAMS), PRIVATE :: aval_params

	PUBLIC :: get_avalanche_PDF_params
	PRIVATE :: initialize_avalanche_params,save_avalanche_params

	CONTAINS

SUBROUTINE get_avalanche_PDF_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params

	call initialize_avalanche_params(params)

	call save_avalanche_params(params)
END SUBROUTINE get_avalanche_PDF_params

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

	aval_params%max_pitch_angle = max_pitch_angle
	aval_params%max_energy = max_energy*C_E ! In Joules
	aval_params%ne = ne
	aval_params%Zeff = Zeff
	aval_params%Ec = Ec
	aval_params%Epar = Epar
	aval_params%Ebar = aval_params%Epar/aval_params%Ec
	aval_params%Te = Te*C_E ! In Joules

	aval_params%max_p =&
	SQRT((aval_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mec^2
	aval_params%min_p = SQRT(aval_params%Ebar - 1.0_rp) ! In units of mec^2
	aval_params%min_energy = SQRT(1.0_rp + aval_params%min_p**2)*C_ME*C_C**2
	aval_params%lD = SQRT(C_E0*aval_params%Te/(aval_params%ne*C_E**2))
	aval_params%bmin =&
	aval_params%Zeff/(12.0_rp*C_PI*aval_params%ne*aval_params%lD**2)
	aval_params%CoulombLog = LOG(aval_params%lD/aval_params%bmin)
	aval_params%Tau =&
	1.0_rp/(4.0_rp*C_PI*C_C*C_RE**2*aval_params%ne*aval_params%CoulombLog)
END SUBROUTINE initialize_avalanche_params


SUBROUTINE save_avalanche_params(params)
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

	filename = TRIM(params%path_to_outputs) // "avalanche_parameters.h5"
	call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

	gname = "avalanche_pdf_params"
	call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

	dset = TRIM(gname) // "/max_pitch_angle"
	attr = "Maximum pitch angle in avalanche PDF (degrees)"
	call save_to_hdf5(h5file_id,dset,aval_params%max_pitch_angle,attr)

	dset = TRIM(gname) // "/min_energy"
	attr = "Minimum energy in avalanche PDF (eV)"
	units = 1.0_rp/C_E
	call save_to_hdf5(h5file_id,dset,units*aval_params%min_energy,attr)

	dset = TRIM(gname) // "/max_energy"
	attr = "Maximum energy in avalanche PDF (eV)"
	units = 1.0_rp/C_E
	call save_to_hdf5(h5file_id,dset,units*aval_params%max_energy,attr)

	dset = TRIM(gname) // "/max_p"
	attr = "Maximum momentum in avalanche PDF (me*c^2)"
	call save_to_hdf5(h5file_id,dset,aval_params%max_p,attr)

	dset = TRIM(gname) // "/min_p"
	attr = "Maximum momentum in avalanche PDF (me*c^2)"
	call save_to_hdf5(h5file_id,dset,aval_params%min_p,attr)

	dset = TRIM(gname) // "/ne"
	attr = "Background electron density (m^-3)"
	call save_to_hdf5(h5file_id,dset,aval_params%ne,attr)

	dset = TRIM(gname) // "/Zeff"
	attr = "Effective atomic number of ions."
	call save_to_hdf5(h5file_id,dset,aval_params%Zeff,attr)

	dset = TRIM(gname) // "/Ec"
	attr = "Critical electric field in (V/m)"
	call save_to_hdf5(h5file_id,dset,aval_params%Ec,attr)

	dset = TRIM(gname) // "/Epar"
	attr = "Parallel electric field in (V/m)"
	call save_to_hdf5(h5file_id,dset,aval_params%Epar,attr)

	dset = TRIM(gname) // "/Te"
	attr = "Background electron temperature (eV)"
	units = 1.0_rp/C_E
	call save_to_hdf5(h5file_id,dset,units*aval_params%Te,attr)

	dset = TRIM(gname) // "/lambda_D"
	attr = "Debye length (m)"
	call save_to_hdf5(h5file_id,dset,aval_params%lD,attr)

	dset = TRIM(gname) // "/bmin"
	attr = "Maximum approach radius (m)"
	call save_to_hdf5(h5file_id,dset,aval_params%bmin,attr)

	dset = TRIM(gname) // "/Clog"
	attr = "Coulomb logarithm"
	call save_to_hdf5(h5file_id,dset,aval_params%CoulombLog,attr)

	dset = TRIM(gname) // "/Tau"
	attr = "Collision time (s)"
	call save_to_hdf5(h5file_id,dset,aval_params%Tau,attr)

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_avalanche_params

END MODULE korc_avalanche
