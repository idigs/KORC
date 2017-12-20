MODULE korc_energy_pdfs
	USE korc_types
	USE korc_constants
	USE korc_HDF5
	USE korc_hpc
	IMPLICIT NONE

	TYPE, PRIVATE :: GAMMA_PARAMS
		REAL(rp) :: min_energy ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_energy ! Maximum energy of sampled PDF in MeV
		REAL(rp) :: min_p ! Minimum momentum of sampled PDF
		REAL(rp) :: max_p ! Maximum momentum of sampled PDF
		REAL(rp) :: k ! Shape factor of Gamma distribution
		REAL(rp) :: t ! Scale factor of Gamma distribution
	END TYPE GAMMA_PARAMS

	TYPE(GAMMA_PARAMS), PRIVATE :: gamma_pdf_params
	REAL(rp), PRIVATE, PARAMETER :: xo = (C_ME*C_C**2/C_E)/1.0E6
	REAL(rp), PRIVATE, PARAMETER :: co = (C_E*1.0E6)/(C_ME*C_C**2)
	REAL(rp), PRIVATE, PARAMETER :: minmax_buffer_size = 10.0_rp

	PUBLIC :: get_gamma_distribution
	PRIVATE :: initialize_gamma_params,&
				save_gamma_params,&
				sample_gamma_distribution,&
				deg2rad,&
				fRE,&
				random_norm,&
				fGamma

	CONTAINS

SUBROUTINE get_gamma_distribution(params,g,go)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
	REAL(rp), INTENT(OUT) :: go

	call initialize_gamma_params(params)

	call save_gamma_params(params)

	call sample_gamma_distribution(params,g,go)

END SUBROUTINE get_gamma_distribution


SUBROUTINE initialize_gamma_params(params)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: max_energy
	REAL(rp) :: min_energy
	REAL(rp) :: Zeff
	REAL(rp) :: E
	REAL(rp) :: k
	REAL(rp) :: t
	NAMELIST /EnergyGammaPDF/ max_energy,min_energy,k,t

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=EnergyGammaPDF)
	close(default_unit_open)

	gamma_pdf_params%min_energy = min_energy*C_E ! In Joules	
	gamma_pdf_params%max_energy = max_energy*C_E ! In Joules
	gamma_pdf_params%k = k
	gamma_pdf_params%t = t

	gamma_pdf_params%max_p = SQRT((gamma_pdf_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mc
	gamma_pdf_params%min_p = SQRT((gamma_pdf_params%min_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mc
END SUBROUTINE initialize_gamma_params


FUNCTION deg2rad(x)
	REAL(rp), INTENT(IN) :: x
	REAL(rp) :: deg2rad
	
	deg2rad = C_PI*x/180.0_rp
END FUNCTION


FUNCTION fGamma(x,k,t)
	REAL(rp), INTENT(IN) :: x ! Independent variable
	REAL(rp), INTENT(IN) :: k ! Shape factor
	REAL(rp), INTENT(IN) :: t ! Scale factor
	REAL(rp) :: fGamma

	fGamma = x**(k - 1.0_rp)*EXP(-x/t)/(GAMMA(k)*t**k)
END FUNCTION fGamma


FUNCTION fRE(p)
	REAL(rp), INTENT(IN) :: p ! momentum in units of mc
	REAL(rp) :: fRE
	REAL(rp) :: Eo ! In units of mc^2

	Eo = SQRT(p**2.0_rp + 1.0_rp)

	fRE = fGamma(Eo,gamma_pdf_params%k,gamma_pdf_params%t*co)
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


SUBROUTINE sample_gamma_distribution(params,g,go)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
	REAL(rp), INTENT(OUT) :: go
	REAL(rp) :: go_root
	REAL(rp), DIMENSION(:), ALLOCATABLE :: p
	REAL(rp) :: p_buffer, p_test
	REAL(rp) :: ratio, rand_unif
	REAL(rp), DIMENSION(:), ALLOCATABLE :: p_samples
	REAL(rp) :: deta
	REAL(rp) :: dp
	INTEGER :: ii,ppp,nsamples
	INTEGER :: mpierr

	ppp = SIZE(g)
	nsamples = ppp*params%mpi_params%nmpi
	ALLOCATE(p(ppp))

	dp = 1.0_rp

	if (params%mpi_params%rank.EQ.0_idef) then
		ALLOCATE(p_samples(nsamples))! Number of samples to distribute among all MPI processes

		call RANDOM_SEED()
		call RANDOM_NUMBER(rand_unif)
		p_buffer = gamma_pdf_params%min_p + (gamma_pdf_params%max_p - gamma_pdf_params%min_p)*rand_unif

		ii=2_idef
		do while (ii .LE. 1000000_idef)
			p_test = p_buffer + random_norm(0.0_rp,dp)
			do while ((p_test.LT.gamma_pdf_params%min_p).OR.(p_test.GT.gamma_pdf_params%max_p))
				p_test = p_buffer + random_norm(0.0_rp,dp)
			end do

			ratio = fRE(p_test)/fRE(p_buffer)

			if (ratio .GE. 1.0_rp) then
				p_buffer = p_test
				ii = ii + 1_idef
			else 
				call RANDOM_NUMBER(rand_unif)
				if (rand_unif .LT. ratio) then
					p_buffer = p_test
					ii = ii + 1_idef
				end if
			end if
		end do	

		call RANDOM_SEED()
		call RANDOM_NUMBER(rand_unif)
		p_samples(1) = p_buffer

		ii=2_idef
		do while (ii .LE. nsamples)
			p_test = p_samples(ii-1) + random_norm(0.0_rp,dp)
			do while ((p_test.LT.gamma_pdf_params%min_p).OR.(p_test.GT.gamma_pdf_params%max_p))
				p_test = p_samples(ii-1) + random_norm(0.0_rp,dp)
			end do

			ratio = fRE(p_test)/fRE(p_samples(ii-1))

			if (ratio .GE. 1.0_rp) then
				p_samples(ii) = p_test
				ii = ii + 1_idef
			else 
				call RANDOM_NUMBER(rand_unif)
				if (rand_unif .LT. ratio) then
					p_samples(ii) = p_test
					ii = ii + 1_idef
				end if
			end if
		end do

		go = SUM(SQRT(1.0_rp + p_samples**2))/nsamples
	end if

	CALL MPI_SCATTER(p_samples,ppp,MPI_REAL8,p,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_BCAST(go,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
	
	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	g = SQRT(1.0_rp + p**2)

	DEALLOCATE(p)
	if (params%mpi_params%rank.EQ.0_idef) then
		DEALLOCATE(p_samples)
	end if

END SUBROUTINE sample_gamma_distribution


SUBROUTINE save_gamma_params(params)
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
		filename = TRIM(params%path_to_outputs) // "gamma_distribution_parameters.h5"
		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		gname = "params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/min_energy"
		attr = "Minimum energy in avalanche PDF (eV)"
		units = 1.0_rp/C_E
		call save_to_hdf5(h5file_id,dset,units*gamma_pdf_params%min_energy,attr)

		dset = TRIM(gname) // "/max_energy"
		attr = "Maximum energy in avalanche PDF (eV)"
		units = 1.0_rp/C_E
		call save_to_hdf5(h5file_id,dset,units*gamma_pdf_params%max_energy,attr)

		dset = TRIM(gname) // "/max_p"
		attr = "Maximum momentum in avalanche PDF (me*c^2)"
		call save_to_hdf5(h5file_id,dset,gamma_pdf_params%max_p,attr)

		dset = TRIM(gname) // "/min_p"
		attr = "Maximum momentum in avalanche PDF (me*c^2)"
		call save_to_hdf5(h5file_id,dset,gamma_pdf_params%min_p,attr)

		dset = TRIM(gname) // "/k"
		attr = "Shape factor"
		call save_to_hdf5(h5file_id,dset,gamma_pdf_params%k,attr)

		dset = TRIM(gname) // "/t"
		attr = "Scale factor"
		call save_to_hdf5(h5file_id,dset,gamma_pdf_params%t,attr)

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if
END SUBROUTINE save_gamma_params

END MODULE korc_energy_pdfs
