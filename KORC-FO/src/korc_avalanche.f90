MODULE korc_avalanche
	USE korc_types
	USE korc_constants
	USE korc_HDF5
	USE korc_hpc
	IMPLICIT NONE

	TYPE, PRIVATE :: AVALANCHE_PDF_PARAMS
		REAL(rp) :: max_pitch_angle ! Maximum pitch angle of sampled PDF in degrees
		REAL(rp) :: min_pitch_angle ! Minimum pitch angle of sampled PDF in degrees
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

		REAL(rp) :: fo
		REAL(rp) :: alpha
		REAL(rp) :: cz
		REAL(rp) :: C1
		REAL(rp) :: C2
	END TYPE AVALANCHE_PDF_PARAMS

	TYPE(AVALANCHE_PDF_PARAMS), PRIVATE :: aval_params

	PUBLIC :: get_avalanche_PDF_params
	PRIVATE :: initialize_avalanche_params,save_avalanche_params,random_norm,fRE,sample_distribution

	CONTAINS

SUBROUTINE get_avalanche_PDF_params(params,g,eta)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta

	call initialize_avalanche_params(params)

	call save_avalanche_params(params)

	call sample_distribution(params,g,eta)
END SUBROUTINE get_avalanche_PDF_params


SUBROUTINE initialize_avalanche_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: max_pitch_angle
	REAL(rp) :: min_pitch_angle
	REAL(rp) :: max_energy
	REAL(rp) :: min_energy
	REAL(rp) :: ne
	REAL(rp) :: Zeff
	REAL(rp) :: Epar
	REAL(rp) :: Te
	NAMELIST /AvalancheGenerationPDF/ max_pitch_angle,min_pitch_angle,max_energy,min_energy,ne,Zeff,Epar,Te

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=AvalancheGenerationPDF)
	close(default_unit_open)

!	write(*,nml=AvalancheGenerationPDF)

	aval_params%max_pitch_angle = max_pitch_angle
	aval_params%min_pitch_angle = min_pitch_angle
	aval_params%max_energy = max_energy*C_E ! In Joules
	aval_params%ne = ne
	aval_params%Zeff = Zeff
	aval_params%Te = Te*C_E ! In Joules

	aval_params%lD = SQRT(C_E0*aval_params%Te/(aval_params%ne*C_E**2))
	aval_params%bmin = aval_params%Zeff/(12.0_rp*C_PI*aval_params%ne*aval_params%lD**2)
	aval_params%CoulombLog = LOG(aval_params%lD/aval_params%bmin)
	aval_params%Tau = 1.0_rp/(4.0_rp*C_PI*C_C*C_RE**2*aval_params%ne*aval_params%CoulombLog)

	aval_params%Ec = C_ME*C_C/(C_E*aval_params%Tau)
	aval_params%Epar = Epar
	aval_params%Ebar = aval_params%Epar/aval_params%Ec

	if (min_energy .EQ. 0.0_rp) then
		aval_params%max_p = SQRT((aval_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mec^2
		aval_params%min_p = SQRT(aval_params%Ebar - 1.0_rp) ! In units of mec^2

		aval_params%min_energy = SQRT(1.0_rp + aval_params%min_p**2)*C_ME*C_C**2
	else
		aval_params%min_energy = min_energy*C_E ! In Joules

		aval_params%max_p = SQRT((aval_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mec^2
		aval_params%min_p = SQRT((aval_params%min_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mec^2
	end if

	aval_params%alpha = (aval_params%Ebar - 1.0_rp)/(1.0_rp + aval_params%Zeff)
	aval_params%cz = SQRT(3.0_rp*(aval_params%Zeff + 5.0_rp)/C_PI)*aval_params%CoulombLog
	aval_params%fo = aval_params%alpha/aval_params%cz
	aval_params%C1 = 0.5_rp*aval_params%alpha
	aval_params%C2 = 1.0_rp/aval_params%cz - aval_params%C1
END SUBROUTINE initialize_avalanche_params


FUNCTION deg2rad(x)
	REAL(rp), INTENT(IN) :: x
	REAL(rp) :: deg2rad
	
	deg2rad = C_PI*x/180.0_rp
END FUNCTION


FUNCTION fRE(x,p)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: x ! x = cos(pitch)
	REAL(rp), INTENT(IN) :: p ! momentum
	REAL(rp) :: fRE
	
	fRE = aval_params%fo*p*EXP(-p*(aval_params%C2*x + aval_params%C1/x))/x
END FUNCTION fRE


FUNCTION random_norm(mean,sigma)
	REAL(rp), INTENT(IN) :: mean
	REAL(rp), INTENT(IN) :: sigma
	REAL(rp) :: random_norm
	REAL(rp) :: rand1, rand2

	call RANDOM_NUMBER(rand1)
	call RANDOM_NUMBER(rand2)

	random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
END FUNCTION random_norm


SUBROUTINE sample_distribution(params,g,eta)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
	REAL(rp), DIMENSION(:), ALLOCATABLE :: p
	REAL(rp) :: chi, chi_test
	REAL(rp) :: p_buffer, p_test
	REAL(rp) :: eta_buffer, eta_test
	REAL(rp) :: ratio, rand_unif
	REAL(rp), DIMENSION(:), ALLOCATABLE :: p_samples
	REAL(rp), DIMENSION(:), ALLOCATABLE :: eta_samples
	REAL(4), DIMENSION(2) :: tarray
	REAL(4) :: time_elapsed
	INTEGER :: ii,ppp,nsamples
	INTEGER :: mpierr

!	call DTIME(tarray,time_elapsed)
!	write(6,'("Time elapsed: ",F15.10)') time_elapsed

	ppp = SIZE(g)
	nsamples = ppp*params%mpi_params%nmpi
	ALLOCATE(p(ppp))

	if (params%mpi_params%rank.EQ.0_idef) then
		ALLOCATE(p_samples(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(eta_samples(nsamples))! Number of samples to distribute among all MPI processes

		call RANDOM_SEED()
		call RANDOM_NUMBER(rand_unif)
		eta_buffer = aval_params%max_pitch_angle*rand_unif
		call RANDOM_NUMBER(rand_unif)
		p_buffer = aval_params%min_p + (aval_params%max_p - aval_params%min_p)*rand_unif

		ii=2_idef
		do while (ii .LE. 1000000_idef)
			eta_test = eta_buffer + random_norm(0.0_rp,1.0_rp)
			do while ((ABS(eta_test) .GT. aval_params%max_pitch_angle).OR.(ABS(eta_test) .LT. aval_params%min_pitch_angle))
				eta_test = eta_buffer + random_norm(0.0_rp,1.0_rp)
			end do
			chi_test = COS(deg2rad(eta_test))
			chi = COS(deg2rad(eta_buffer))

			p_test = p_buffer + random_norm(0.0_rp,1.0_rp)
			do while ((p_test.LT.aval_params%min_p).OR.(p_test.GT.aval_params%max_p))
				p_test = p_buffer + random_norm(0.0_rp,1.0_rp)
			end do

			ratio = fRE(chi_test,p_test)/fRE(chi,p_buffer)

			if (ratio .GE. 1.0_rp) then
				p_buffer = p_test
				eta_buffer = eta_test
				ii = ii + 1_idef
			else 
				call RANDOM_NUMBER(rand_unif)
				if (rand_unif .LT. ratio) then
					p_buffer = p_test
					eta_buffer = eta_test
					ii = ii + 1_idef
				end if
			end if
		end do	

		eta_samples(1) = eta_buffer
		call RANDOM_SEED()
		call RANDOM_NUMBER(rand_unif)
		p_samples(1) = p_buffer

		ii=2_idef
		do while (ii .LE. nsamples)
			eta_test = eta_samples(ii-1) + random_norm(0.0_rp,1.0_rp)
			do while ((ABS(eta_test) .GT. aval_params%max_pitch_angle).OR.(ABS(eta_test) .LT. aval_params%min_pitch_angle))
				eta_test = eta_samples(ii-1) + random_norm(0.0_rp,1.0_rp)
			end do
			chi_test = COS(deg2rad(eta_test))
			chi = COS(deg2rad(eta_samples(ii-1)))

			p_test = p_samples(ii-1) + random_norm(0.0_rp,1.0_rp)
			do while ((p_test.LT.aval_params%min_p).OR.(p_test.GT.aval_params%max_p))
				p_test = p_samples(ii-1) + random_norm(0.0_rp,1.0_rp)
			end do

			ratio = fRE(chi_test,p_test)/fRE(chi,p_samples(ii-1))

			if (ratio .GE. 1.0_rp) then
				p_samples(ii) = p_test
				eta_samples(ii) = eta_test
				ii = ii + 1_idef
			else 
				call RANDOM_NUMBER(rand_unif)
				if (rand_unif .LT. ratio) then
					p_samples(ii) = p_test
					eta_samples(ii) = eta_test
					ii = ii + 1_idef
				end if
			end if
		end do	

		do ii=1_idef,ppp
			if (eta_samples(ii).LT.0.0_rp) then
				eta_samples(ii) = -eta_samples(ii)
			end if
		end do

	end if

	CALL MPI_SCATTER(p_samples,ppp,MPI_REAL8,p,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_SCATTER(eta_samples,ppp,MPI_REAL8,eta,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	g = SQRT(1.0_rp + p**2)

!	write(6,'("MPI:",I4," Minimum gamma: ",F30.16," "I6)') params%mpi_params%rank,MINVAL(p),indices(MINLOC(p))
!	call DTIME(tarray,time_elapsed)
!	write(6,'("Time elapsed: ",F15.10)') time_elapsed

	DEALLOCATE(p)
	if (params%mpi_params%rank.EQ.0_idef) then
		DEALLOCATE(p_samples)
		DEALLOCATE(eta_samples)
	end if
END SUBROUTINE sample_distribution


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

	if (params%mpi_params%rank .EQ. 0) then
		filename = TRIM(params%path_to_outputs) // "avalanche_parameters.h5"
		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		gname = "avalanche_pdf_params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/max_pitch_angle"
		attr = "Maximum pitch angle in avalanche PDF (degrees)"
		call save_to_hdf5(h5file_id,dset,aval_params%max_pitch_angle,attr)

		dset = TRIM(gname) // "/min_pitch_angle"
		attr = "Minimum pitch angle in avalanche PDF (degrees)"
		call save_to_hdf5(h5file_id,dset,aval_params%min_pitch_angle,attr)

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
	end if
END SUBROUTINE save_avalanche_params

END MODULE korc_avalanche
