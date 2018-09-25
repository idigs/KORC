MODULE korc_experimental_pdf
	USE korc_types
	USE korc_constants
	USE korc_HDF5
	USE korc_hpc
    USE special_functions

	IMPLICIT NONE

	TYPE, PRIVATE :: PARAMS
		REAL(rp) :: E ! Parallel electric field normalized using the critical electric field
		REAL(rp) :: Zeff ! Effective atomic number of impurities

		REAL(rp) :: max_pitch_angle ! Maximum pitch angle of sampled PDF in degrees
		REAL(rp) :: min_pitch_angle ! Minimum pitch angle of sampled PDF in degrees
		REAL(rp) :: min_energy ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_energy ! Maximum energy of sampled PDF in MeV
		REAL(rp) :: min_p ! Minimum momentum of sampled PDF
		REAL(rp) :: max_p ! Maximum momentum of sampled PDF
		REAL(rp) :: k ! Shape factor of Gamma distribution
		REAL(rp) :: t ! Scale factor of Gamma distribution
		REAL(rp) :: fGo ! Normalization factor of Gamma distribution

		REAL(rp) :: Bo
		REAL(rp) :: lambda

        REAL(rp) :: A_fact ! Multiplication factor for A in distributon.
	END TYPE PARAMS

	TYPE, PRIVATE :: HOLLMANN_PARAMS
	    CHARACTER(MAX_STRING_LENGTH) :: filename
		REAL(rp) :: E
		REAL(rp) :: Zeff
		REAL(rp) :: max_pitch_angle
		REAL(rp) :: min_pitch_angle
		REAL(rp) :: min_sampling_energy ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_sampling_energy ! Maximum energy of sampled PDF in MeV
		REAL(rp) :: min_sampling_g ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_sampling_g ! Maximum energy of sampled PDF in MeV

		REAL(rp) :: min_energy ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_energy ! Maximum energy of sampled PDF in MeV
		REAL(rp) :: min_g ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_g ! Maximum energy of sampled PDF in MeV
		REAL(rp) :: min_pitch ! Minimum energy of sampled PDF in MeV
		REAL(rp) :: max_pitch ! Maximum energy of sampled PDF in MeV

		INTEGER :: N

		REAL(rp), DIMENSION(:), ALLOCATABLE :: E_axis
		REAL(rp), DIMENSION(:), ALLOCATABLE :: g
		REAL(rp), DIMENSION(:), ALLOCATABLE :: fRE_E
		REAL(rp), DIMENSION(:), ALLOCATABLE :: fRE_pitch

		CHARACTER(MAX_STRING_LENGTH) :: current_direction
		REAL(rp) :: Bo
		REAL(rp) :: lambda

        REAL(rp) :: A_fact ! Multiplication factor for A in distributon.
	END TYPE HOLLMANN_PARAMS

	TYPE(PARAMS), PRIVATE :: pdf_params
	TYPE(HOLLMANN_PARAMS), PRIVATE :: h_params
	REAL(rp), PRIVATE, PARAMETER :: xo = (C_ME*C_C**2/C_E)/1.0E6
	REAL(rp), PRIVATE, PARAMETER :: Tol = 1.0E-5_rp
	REAL(rp), PRIVATE, PARAMETER :: minmax_buffer_size = 10.0_rp

	PUBLIC :: get_experimentalG_distribution,&
				get_Hollmann_distribution,&
				initialize_Hollmann_params,&
				sample_Hollmann_distribution
	PRIVATE :: initialize_params,&
				save_params,&
				sample_distribution,&
				deg2rad,&
				rad2deg,&
				fRE,&
				fRExPR,&
				random_norm,&
				fGamma,&
				PR,&
				P_integral,&
				IntK,&
				IntBesselK,&
				IntGamma,&
				fRE_H,&
				fRE_pitch

	CONTAINS

SUBROUTINE get_experimentalG_distribution(params,g,eta,go,etao)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
	REAL(rp), INTENT(OUT) :: go
	REAL(rp), INTENT(OUT) :: etao

	call initialize_params(params)

	call save_params(params)

	call sample_distribution(params,g,eta,go,etao)
END SUBROUTINE get_experimentalG_distribution


SUBROUTINE initialize_params(params)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: max_pitch_angle
	REAL(rp) :: min_pitch_angle
	REAL(rp) :: max_energy
	REAL(rp) :: min_energy
	REAL(rp) :: Zeff
	REAL(rp) :: E
	REAL(rp) :: k
	REAL(rp) :: t
	REAL(rp) :: Bo
	REAL(rp) :: lambda
    REAL(rp) :: A_fact

	NAMELIST /ExperimentalPDF/ max_pitch_angle,min_pitch_angle,max_energy,min_energy,Zeff,E,k,t,Bo,lambda, &
                               A_fact

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=ExperimentalPDF)
	close(default_unit_open)

	pdf_params%max_pitch_angle = max_pitch_angle
	pdf_params%min_pitch_angle = min_pitch_angle
	pdf_params%min_energy = min_energy*C_E ! In Joules
	pdf_params%max_energy = max_energy*C_E ! In Joules
	pdf_params%Zeff = Zeff
	pdf_params%E = E
	pdf_params%k = k
	pdf_params%t = t
	pdf_params%Bo = Bo
	pdf_params%lambda = lambda

	pdf_params%max_p = SQRT((pdf_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mc
	pdf_params%min_p = SQRT((pdf_params%min_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mc

	pdf_params%fGo = &
	IntGamma(SQRT(pdf_params%min_p**2.0_rp + 1.0_rp),SQRT(pdf_params%max_p**2.0_rp + 1.0_rp),pdf_params%k,pdf_params%t/xo)

    pdf_params%A_fact = A_fact
END SUBROUTINE initialize_params


FUNCTION deg2rad(x)
	REAL(rp), INTENT(IN) :: x
	REAL(rp) :: deg2rad

	deg2rad = C_PI*x/180.0_rp
END FUNCTION


FUNCTION rad2deg(x)
	REAL(rp), INTENT(IN) :: x
	REAL(rp) :: rad2deg

	rad2deg = 180.0_rp*x/C_PI
END FUNCTION


FUNCTION fGamma(x,k,t)
	REAL(rp), INTENT(IN) :: x ! Independent variable
	REAL(rp), INTENT(IN) :: k ! Shape factor
	REAL(rp), INTENT(IN) :: t ! Scale factor
	REAL(rp) :: fGamma

	fGamma = x**(k - 1.0_rp)*EXP(-x/t)/(GAMMA(k)*t**k)
END FUNCTION fGamma


FUNCTION fRE(eta,p)
	REAL(rp), INTENT(IN) 	:: eta ! pitch angle in degrees
	REAL(rp), INTENT(IN) 	:: p ! momentum in units of mc
	REAL(rp) 				:: fRE
	REAL(rp) 				:: fE
	REAL(rp) 				:: feta
	REAL(rp) 				:: A
	REAL(rp) 				:: Eo

	Eo = SQRT(p**2.0_rp + 1.0_rp)

	A = (2.0_rp*pdf_params%E/(pdf_params%Zeff + 1.0_rp))*(p**2/SQRT(p**2.0_rp + 1.0_rp))
    A = A*pdf_params%A_fact
	feta = 0.5_rp*A*EXP(A*COS(deg2rad(eta)))/SINH(A)

	fE = fGamma(Eo,pdf_params%k,pdf_params%t/xo)/pdf_params%fGo

	fRE = fE*feta
END FUNCTION fRE


FUNCTION fRExPR(eta,p)
	REAL(rp), INTENT(IN) :: eta ! pitch angle in degrees
	REAL(rp), INTENT(IN) :: p ! momentum in units of mc
	REAL(rp) :: fRExPR
	REAL(rp) :: A
	REAL(rp) :: Eo

	fRExPR = fRE(eta,p)*PR(eta,p,pdf_params%Bo,pdf_params%lambda)
END FUNCTION fRExPR


FUNCTION random_norm(mean,sigma)
	REAL(rp), INTENT(IN) :: mean
	REAL(rp), INTENT(IN) :: sigma
	REAL(rp) :: random_norm
	REAL(rp) :: rand1, rand2

	call RANDOM_NUMBER(rand1)
	call RANDOM_NUMBER(rand2)

	random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
END FUNCTION random_norm


FUNCTION IntK(v,x)
	REAL(rp) :: IntK
	REAL(rp), INTENT(IN) :: v
	REAL(rp), INTENT(IN) :: x

	IntK = (C_PI/SQRT(2.0_rp))*(1.0_rp - 0.25_rp*(4.0_rp*v**2 - 1.0_rp))*ERFC(SQRT(x))&
			 + 0.25_rp*(4.0_rp*v**2 - 1.0_rp)*SQRT(0.5_rp*C_PI/x)*EXP(-x)
END FUNCTION IntK


FUNCTION besselk(v,x)
	REAL(rp) :: besselk
	REAL(rp), INTENT(IN) :: x
	REAL(rp), INTENT(IN) :: v
	REAL(4) :: ri,rk,rip,rkp

	call bessik(REAL(x,4),REAL(v,4),ri,rk,rip,rkp)
	besselk = REAL(rk,rp)
END FUNCTION besselk


!> @brief Extended trapezoidal rule for integrating the Gamma PDF. See Sec. 4.2 of Numerical Recipies in Fortran 77.
FUNCTION IntGamma(a,b,k,t)
	REAL(rp), INTENT(IN) :: a
	REAL(rp), INTENT(IN) :: b
	REAL(rp), INTENT(IN) :: k ! shape factor
	REAL(rp), INTENT(IN) :: t ! scale factor
	REAL(rp) :: IntGamma
	REAL(rp) :: Iold
	REAL(rp) :: Inew
	REAL(rp) :: rerr
	REAL(rp) :: sum_f
	REAL(rp) :: h,z
	INTEGER :: ii,jj,npoints
	LOGICAL :: flag

	h = b - a
	sum_f = 0.5*(fGamma(a,k,t) + fGamma(b,k,t))

	Iold = 0.0_rp
	Inew = sum_f*h

	ii = 1_idef
	flag = .TRUE.
	do while (flag)
		Iold = Inew

		ii = ii + 1_idef
		npoints = 2_idef**(ii-2_idef)
		h = 0.5_rp*(b-a)/REAL(npoints,rp)
		sum_f = 0.0_rp
		do jj=1_idef,npoints
			z = a + h + 2.0_rp*(REAL(jj,rp) - 1.0_rp)*h
			sum_f = sum_f + fGamma(z,k,t)
		end do

		Inew = 0.5_rp*Iold + sum_f*h
		rerr = ABS((Inew - Iold)/Iold)
		flag = .NOT.(rerr.LT.Tol)
	end do
	IntGamma = Inew
END FUNCTION IntGamma

!> @brief Extended trapezoidal rule for integrating the modified Bessel function of second kind. See Sec. 4.2 of Numerical Recipies in Fortran 77.
FUNCTION IntBesselK(a,b)
	REAL(rp), INTENT(IN) :: a
	REAL(rp), INTENT(IN) :: b
	REAL(rp) :: IntBesselK
	REAL(rp) :: Iold
	REAL(rp) :: Inew
	REAL(rp) :: rerr
	REAL(rp) :: sum_f
	REAL(rp) :: v,h,z
	INTEGER :: ii,jj,npoints
	LOGICAL :: flag

	v = 5.0_rp/3.0_rp
	h = b - a
	sum_f = 0.5*(besselk(v,a) + besselk(v,b))

	Iold = 0.0_rp
	Inew = sum_f*h

	ii = 1_idef
	flag = .TRUE.
	do while (flag)
		Iold = Inew

		ii = ii + 1_idef
		npoints = 2_idef**(ii-2_idef)
		h = 0.5_rp*(b-a)/REAL(npoints,rp)
		sum_f = 0.0_rp
		do jj=1_idef,npoints
			z = a + h + 2.0_rp*(REAL(jj,rp) - 1.0_rp)*h
			sum_f = sum_f + besselk(v,z)
		end do

		Inew = 0.5_rp*Iold + sum_f*h
		rerr = ABS((Inew - Iold)/Iold)
		flag = .NOT.(rerr.LT.Tol)
	end do
	IntBesselK = Inew
END FUNCTION IntBesselK


SUBROUTINE P_integral(z,P)
	REAL(rp), INTENT(OUT) :: P
	REAL(rp), INTENT(IN) :: z
	REAL(rp) :: a

	P = 0.0_rp

	IF (z .LT. 0.5_rp) THEN
		a = (2.16_rp/2.0_rp**(2.0_rp/3.0_rp))*z**(1.0_rp/3.0_rp)
		P = IntBesselK(z,a) + IntK(5.0_rp/3.0_rp,a)
	ELSE IF ((z .GE. 0.5_rp).AND.(z .LT. 2.5_rp)) THEN
		a = 0.72_rp*(z + 1.0_rp)
		P = IntBesselK(z,a) + IntK(5.0_rp/3.0_rp,a)
	ELSE
		P = IntK(5.0_rp/3.0_rp,z)
	END IF
END SUBROUTINE P_integral


FUNCTION PR(eta,p,Bo,l)
	REAL(rp), INTENT(IN) :: eta ! in radians
	REAL(rp), INTENT(IN) :: p ! dimensionless (in units of mc)
	REAL(rp), INTENT(IN) :: Bo
	REAL(rp), INTENT(IN) :: l
	REAL(rp) :: PR
	REAL(rp) :: g
	REAL(rp) :: v
	REAL(rp) :: k
	REAL(rp) :: lc
	REAL(rp) :: z
	REAL(rp) :: Pi

	g = SQRT(p**2 + 1.0_rp)
	v = C_C*SQRT(1.0_rp - 1.0_rp/g**2)

	k = C_E*Bo*SIN(deg2rad(eta))/(g*C_ME*v)

	lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

	z = lc/l

	call P_integral(z,Pi)

	PR = (C_C*C_E**2)*Pi/(SQRT(3.0_rp)*C_E0*g**2*l**3)
END FUNCTION PR


SUBROUTINE sample_distribution(params,g,eta,go,etao)
	TYPE(KORC_PARAMS), INTENT(IN) 						:: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: g
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: eta
	REAL(rp), INTENT(OUT) 								:: go
	REAL(rp), INTENT(OUT) 								:: etao
	REAL(rp) 											:: go_root
	REAL(rp) 											:: etao_root
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: p
	REAL(rp) 											:: p_buffer
	REAL(rp) 											:: p_test
	REAL(rp) 											:: eta_buffer
	REAL(rp) 											:: eta_test
	REAL(rp) 											:: ratio
	REAL(rp) 											:: rand_unif
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: p_samples
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: eta_samples
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: p_tmp
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: eta_tmp
	REAL(rp) 											:: minmax
	REAL(rp) 											:: min_p
	REAL(rp) 											:: max_p
	REAL(rp) 											:: min_pitch_angle
	REAL(rp) 											:: max_pitch_angle
	REAL(rp) 											:: deta
	REAL(rp) 											:: dp
	LOGICAL 											:: lp
	LOGICAL 											:: leta
	INTEGER 											:: num_accepted
	INTEGER 											:: ii
	INTEGER 											:: jj
	INTEGER 											:: ppp
	INTEGER 											:: nsamples
	INTEGER 											:: mpierr

	ppp = SIZE(g)
	nsamples = ppp*params%mpi_params%nmpi
	ALLOCATE(p(ppp))

	deta = (pdf_params%max_pitch_angle - pdf_params%min_pitch_angle)/100.0_rp
	dp = (pdf_params%max_p - pdf_params%min_p)/100.0_rp

	do jj=1_idef,INT(minmax_buffer_size,idef)
		minmax = pdf_params%min_p - REAL(jj,rp)*dp
		if (minmax.GT.0.0_rp) then
			min_p = minmax
		end if
	end do

	max_p = pdf_params%max_p + minmax_buffer_size*dp

	if (pdf_params%min_pitch_angle.GE.korc_zero) then
		do jj=1_idef,INT(minmax_buffer_size,idef)
			minmax = pdf_params%min_pitch_angle -  REAL(jj,rp)*deta
			if (minmax.GT.0.0_rp) then
				min_pitch_angle = minmax
			end if
		end do
	else
		min_pitch_angle = 0.0_rp
	end if

	do jj=1_idef,INT(minmax_buffer_size,idef)
		minmax = pdf_params%max_pitch_angle + REAL(jj,rp)*deta
		if (minmax.LE.90.0_rp) then
			max_pitch_angle = minmax
		else
			max_pitch_angle = pdf_params%max_pitch_angle
			EXIT
		end if
	end do

	if (params%mpi_params%rank.EQ.0_idef) then
		ALLOCATE(p_samples(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(eta_samples(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(p_tmp(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(eta_tmp(nsamples))! Number of samples to distribute among all MPI processes


		!* * * Transient * * *!
		call RANDOM_SEED()

		call RANDOM_NUMBER(rand_unif)
		eta_buffer = pdf_params%min_pitch_angle + (pdf_params%max_pitch_angle - pdf_params%min_pitch_angle)*rand_unif

		call RANDOM_NUMBER(rand_unif)
		p_buffer = pdf_params%min_p + (pdf_params%max_p - pdf_params%min_p)*rand_unif

		ii=2_idef
		do while (ii .LE. 1000_idef)
			eta_test = eta_buffer + random_norm(0.0_rp,deta)
			do while ((ABS(eta_test) .GT. pdf_params%max_pitch_angle).OR.(ABS(eta_test) .LT. pdf_params%min_pitch_angle))
				eta_test = eta_buffer + random_norm(0.0_rp,deta)
			end do

			p_test = p_buffer + random_norm(0.0_rp,dp)
			do while ((p_test.LT.pdf_params%min_p).OR.(p_test.GT.pdf_params%max_p))
				p_test = p_buffer + random_norm(0.0_rp,dp)
			end do

			ratio = fRE(eta_test,p_test)/fRE(eta_buffer,p_buffer)

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
		!* * * Transient * * *!

		eta_tmp(1) = eta_buffer
		p_tmp(1) = p_buffer

		num_accepted = 0_idef
		do while(num_accepted.LT.nsamples)
			ii=2_idef
			do while (ii .LE. nsamples)
				eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
				do while ((ABS(eta_test) .GT. max_pitch_angle).OR.(ABS(eta_test) .LT. min_pitch_angle))
					eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
				end do

				p_test = p_tmp(ii-1) + random_norm(0.0_rp,dp)
				do while ((p_test.LT.min_p).OR.(p_test.GT.max_p))
					p_test = p_tmp(ii-1) + random_norm(0.0_rp,dp)
				end do

				ratio = fRE(eta_test,p_test)/fRE(eta_tmp(ii-1),p_tmp(ii-1))

				if (ratio .GE. 1.0_rp) then
					p_tmp(ii) = p_test
					eta_tmp(ii) = eta_test
					ii = ii + 1_idef
				else
					call RANDOM_NUMBER(rand_unif)
					if (rand_unif .LT. ratio) then
						p_tmp(ii) = p_test
						eta_tmp(ii) = eta_test
						ii = ii + 1_idef
					end if
				end if
			end do

			eta_tmp = ABS(eta_tmp)

			ii = 1_idef
			do while ( (ii.LT.nsamples).AND.(num_accepted.LT.nsamples) )
				lp = (p_tmp(ii).LE.pdf_params%max_p).AND.(p_tmp(ii).GE.pdf_params%min_p)
				leta = (eta_tmp(ii).LE.pdf_params%max_pitch_angle).AND.(eta_tmp(ii).GE.pdf_params%min_pitch_angle)
				if (lp.AND.leta) then
					num_accepted = num_accepted + 1_idef
					p_samples(num_accepted) = p_tmp(ii)
					eta_samples(num_accepted) = eta_tmp(ii)
				end if
				ii = ii + 1_idef
			end do
		end do

		go = SUM(SQRT(1.0_rp + p_samples**2))/nsamples
		etao = SUM(eta_samples)/nsamples
	end if

	CALL MPI_SCATTER(p_samples,ppp,MPI_REAL8,p,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_SCATTER(eta_samples,ppp,MPI_REAL8,eta,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_BCAST(go,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_BCAST(etao,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	g = SQRT(1.0_rp + p**2)

	DEALLOCATE(p)
	if (params%mpi_params%rank.EQ.0_idef) then
		DEALLOCATE(p_samples)
		DEALLOCATE(eta_samples)
		DEALLOCATE(p_tmp)
		DEALLOCATE(eta_tmp)
	end if

END SUBROUTINE sample_distribution


SUBROUTINE get_Hollmann_distribution(params,g,eta,go,etao)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
	REAL(rp), INTENT(OUT) :: go
	REAL(rp), INTENT(OUT) :: etao

	call initialize_Hollmann_params(params)

	call save_Hollmann_params(params)

	call sample_Hollmann_distribution(params,g,eta,go,etao)
END SUBROUTINE get_Hollmann_distribution


SUBROUTINE initialize_Hollmann_params(params)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
    CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: current_direction
	REAL(rp) :: E
	REAL(rp) :: Zeff
	REAL(rp) :: max_pitch_angle
	REAL(rp) :: min_pitch_angle
	REAL(rp) :: max_energy
	REAL(rp) :: min_energy
	REAL(rp) :: Bo
	REAL(rp) :: lambda
    REAL(rp) :: A_fact

	NAMELIST /HollmannPDF/ E,Zeff,max_pitch_angle,min_pitch_angle,max_energy,min_energy,filename,Bo,lambda,current_direction, &
                           A_fact


	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=HollmannPDF)
	close(default_unit_open)

	h_params%filename = TRIM(filename)

	h_params%E = E
	h_params%Zeff = Zeff
!	h_params%max_pitch_angle = max_pitch_angle                                  ! MRC
!	h_params%min_pitch_angle = min_pitch_angle                                  ! MRC

    if (TRIM(current_direction) .EQ. 'ANTICLOCKWISE') then                      ! MRC
        h_params%max_pitch_angle = 180.0_rp - min_pitch_angle                   ! MRC
        h_params%min_pitch_angle = 180.0_rp - max_pitch_angle                   ! MRC
    else                                                                        ! MRC
        h_params%max_pitch_angle = max_pitch_angle                              ! MRC
        h_params%min_pitch_angle = min_pitch_angle                              ! MRC
    end if                                                                      ! MRC

	h_params%min_sampling_energy = min_energy*C_E ! In Joules
	h_params%min_sampling_g = 1.0_rp + h_params%min_sampling_energy/(C_ME*C_C**2)
	h_params%max_sampling_energy = max_energy*C_E ! In Joules.
	h_params%max_sampling_g = 1.0_rp + h_params%max_sampling_energy/(C_ME*C_C**2)

	call load_data_from_hdf5()

	ALLOCATE(h_params%g(h_params%N))

	h_params%g = 1.0_rp + h_params%E_axis/(C_ME*C_C**2)
	h_params%max_g = MAXVAL(h_params%g)
	h_params%min_g = MINVAL(h_params%g)

	h_params%current_direction = TRIM(current_direction)

	h_params%Bo = Bo
	h_params%lambda = lambda

    h_params%A_fact = A_fact
END SUBROUTINE initialize_Hollmann_params


SUBROUTINE load_data_from_hdf5()
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	REAL(rp) :: rdatum
	INTEGER :: h5error

	filename = TRIM(h_params%filename)
	call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_data_from_hdf5 (korc_experimental) --> h5fopen_f")')
	end if

	dset = "/N"
	call load_from_hdf5(h5file_id,dset,rdatum)
	h_params%N = INT(rdatum)

	ALLOCATE(h_params%E_axis(h_params%N))
	ALLOCATE(h_params%fRE_E(h_params%N))
	ALLOCATE(h_params%fRE_pitch(h_params%N))

	dset = "/E"
	call load_array_from_hdf5(h5file_id,dset,h_params%E_axis)
	h_params%E_axis = h_params%E_axis*C_E

	dset = "/fRE_E"
	call load_array_from_hdf5(h5file_id,dset,h_params%fRE_E)

	dset = "/fRE_pitch"
	call load_array_from_hdf5(h5file_id,dset,h_params%fRE_pitch)

	call h5fclose_f(h5file_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: load_data_from_hdf5 (korc_experimental) --> h5fclose_f")')
	end if
END SUBROUTINE load_data_from_hdf5


FUNCTION fRE_H(eta,g)
	REAL(rp), INTENT(IN) 	:: eta ! pitch angle in degrees
	REAL(rp), INTENT(IN) 	:: g ! Relativistic gamma factor
	REAL(rp) 				:: fRE_H
	REAL(rp) 				:: D
	REAL(rp) 				:: g0
	REAL(rp) 				:: g1
	REAL(rp) 				:: f0
	REAL(rp) 				:: f1
	REAL(rp) 				:: m
	REAL(rp) 				:: feta
	REAL(rp) 				:: A
	INTEGER 				:: index

	index = MINLOC(ABS(h_params%g - g),1)
	D = h_params%g(index) - g

	if (D.GT.0) then
		f0 = h_params%fRE_E(index-1)
		g0 = h_params%g(index-1)

		f1 = h_params%fRE_E(index)
		g1 = h_params%g(index)
	else
		f0 = h_params%fRE_E(index)
		g0 = h_params%g(index)

		f1 = h_params%fRE_E(index+1)
		g1 = h_params%g(index+1)
	end if

	m = (f1-f0)/(g1-g0)

	fRE_H = f0 + m*(g - g0)

	A = (2.0_rp*h_params%E/(h_params%Zeff + 1.0_rp))*(g**2 - 1.0_rp)/g
    A = A*h_params%A_fact
    feta = A*EXP(-A*(1.0_rp - COS(deg2rad(eta))))/(1.0_rp - EXP(-2.0_rp*A))     ! MRC
!	feta = 0.5_rp*A*EXP(A*COS(deg2rad(eta)))/SINH(A)                            ! MRC

	fRE_H = fRE_H*feta
END FUNCTION fRE_H


FUNCTION fRE_HxPR(eta,g)
	REAL(rp), INTENT(IN) :: eta ! pitch angle in degrees
	REAL(rp), INTENT(IN) :: g ! gamma factor
	REAL(rp) :: fRE_HxPR

	fRE_HxPR = fRE_H(eta,g)*PR(eta,SQRT(g**2 - 1.0_rp),h_params%Bo,h_params%lambda)
END FUNCTION fRE_HxPR


FUNCTION fRE_pitch(g)
	REAL(rp), INTENT(IN) :: g ! Relativistic gamma factor
	REAL(rp) :: fRE_pitch
	REAL(rp) :: D
	REAL(rp) :: g0,g1,f0,f1,m
	INTEGER :: index

	index = MINLOC(ABS(h_params%g - g),1)
	D = h_params%g(index) - g

	if (D.GT.0) then
		f0 = h_params%fRE_pitch(index-1)
		g0 = h_params%g(index-1)

		f1 = h_params%fRE_pitch(index)
		g1 = h_params%g(index)
	else
		f0 = h_params%fRE_pitch(index)
		g0 = h_params%g(index)

		f1 = h_params%fRE_pitch(index+1)
		g1 = h_params%g(index+1)
	end if

	m = (f1-f0)/(g1-g0)

	fRE_pitch = f0 + m*(g - g0)
	fRE_pitch = 180.0_rp - fRE_pitch
END FUNCTION fRE_pitch


SUBROUTINE sample_Hollmann_distribution(params,g,eta,go,etao)
	TYPE(KORC_PARAMS), INTENT(IN) 						:: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: g
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: eta
	REAL(rp), INTENT(OUT) 								:: go
	REAL(rp), INTENT(OUT) 								:: etao
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: p
	REAL(rp) 											:: g_buffer
	REAL(rp) 											:: g_test
	REAL(rp) 											:: eta_buffer
	REAL(rp) 											:: eta_test
	REAL(rp) 											:: ratio
	REAL(rp) 											:: rand_unif
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: g_samples
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: eta_samples
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: g_tmp
	REAL(rp), DIMENSION(:), ALLOCATABLE 				:: eta_tmp
	REAL(rp) 											:: minmax
	REAL(rp) 											:: min_g
	REAL(rp) 											:: max_g
	REAL(rp) 											:: min_pitch_angle
	REAL(rp) 											:: max_pitch_angle
	REAL(rp) 											:: dg
	REAL(rp) 											:: deta
	LOGICAL 											:: lp
	INTEGER 											:: index_i
	INTEGER 											:: index_f
	INTEGER 											:: num_accepted
	INTEGER 											:: ii
	INTEGER 											:: jj
	INTEGER 											:: ppp
	INTEGER 											:: nsamples
	INTEGER 											:: mpierr

	ppp = SIZE(g)
	nsamples = ppp*params%mpi_params%nmpi

	index_i = MINLOC(ABS(h_params%g - h_params%min_sampling_g),1)
	index_f = MINLOC(ABS(h_params%g - h_params%max_sampling_g),1)

	deta = (h_params%max_pitch_angle - h_params%min_pitch_angle)/100.0_rp
	dg = (h_params%max_sampling_g - h_params%min_sampling_g)/100.0_rp

	do jj=1_idef,INT(minmax_buffer_size,idef)
		minmax = h_params%min_sampling_g - REAL(jj,rp)*dg
		if (minmax.GT.h_params%min_g) then
			min_g = minmax
		end if

		minmax = h_params%max_sampling_g + REAL(jj,rp)*dg
		if (minmax.LT.h_params%max_g) then
			max_g = minmax
		end if
	end do

	if (h_params%min_pitch_angle.GE.korc_zero) then
		do jj=1_idef,INT(minmax_buffer_size,idef)
			minmax = h_params%min_pitch_angle -  REAL(jj,rp)*deta
			if (minmax.GT.0.0_rp) then
				min_pitch_angle = minmax
			end if
		end do
	else
		min_pitch_angle = 0.0_rp
	end if

	do jj=1_idef,INT(minmax_buffer_size,idef)
		minmax = h_params%max_pitch_angle + REAL(jj,rp)*deta
		if (minmax.LE.90.0_rp) then
			max_pitch_angle = minmax
		else
			max_pitch_angle = h_params%max_pitch_angle
			EXIT
		end if
	end do

	if (params%mpi_params%rank.EQ.0_idef) then
		ALLOCATE(g_samples(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(eta_samples(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(g_tmp(nsamples))! Number of samples to distribute among all MPI processes
		ALLOCATE(eta_tmp(nsamples))! Number of samples to distribute among all MPI processes


		!* * * Transient * * *!
		call RANDOM_SEED()

		call RANDOM_NUMBER(rand_unif)
		eta_buffer = h_params%min_pitch_angle + (h_params%max_pitch_angle - h_params%min_pitch_angle)*rand_unif

		call RANDOM_NUMBER(rand_unif)
		g_buffer = h_params%min_sampling_g + (h_params%max_sampling_g - h_params%min_sampling_g)*rand_unif

		ii=2_idef
		do while (ii .LE. 1000_idef)
			eta_test = eta_buffer + random_norm(0.0_rp,deta)
			do while ((ABS(eta_test) .GT. h_params%max_pitch_angle).OR.(ABS(eta_test) .LT. h_params%min_pitch_angle))
				eta_test = eta_buffer + random_norm(0.0_rp,deta)
			end do

			g_test = g_buffer + random_norm(0.0_rp,dg)
			do while ((g_test.LT.h_params%min_sampling_g).OR.(g_test.GT.h_params%max_sampling_g))
				g_test = g_buffer + random_norm(0.0_rp,dg)
			end do

			ratio = fRE_H(eta_test,g_test)/fRE_H(eta_buffer,g_buffer)
			!ratio = fRE_HxPR(eta_test,g_test)/fRE_HxPR(eta_buffer,g_buffer)

			if (ratio .GE. 1.0_rp) then
				g_buffer = g_test
				eta_buffer = eta_test
				ii = ii + 1_idef
			else
				call RANDOM_NUMBER(rand_unif)
				if (rand_unif .LT. ratio) then
					g_buffer = g_test
					eta_buffer = eta_test
					ii = ii + 1_idef
				end if
			end if
		end do
		!* * * Transient * * *!

		eta_tmp(1) = eta_buffer
		g_tmp(1) = g_buffer

		num_accepted = 0_idef
		do while(num_accepted.LT.nsamples)
			ii=2_idef
			do while (ii .LE. nsamples)
				eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
				do while ((ABS(eta_test) .GT. max_pitch_angle).OR.(ABS(eta_test) .LT. min_pitch_angle))
					eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
				end do

				g_test = g_tmp(ii-1) + random_norm(0.0_rp,dg)
				do while ((g_test.LT.min_g).OR.(g_test.GT.max_g))
					g_test = g_tmp(ii-1) + random_norm(0.0_rp,dg)
				end do

				ratio = fRE_H(eta_test,g_test)/fRE_H(eta_tmp(ii-1),g_tmp(ii-1))
				!ratio = fRE_HxPR(eta_test,g_test)/fRE_HxPR(eta_tmp(ii-1),g_tmp(ii-1))

				if (ratio .GE. 1.0_rp) then
					g_tmp(ii) = g_test
					eta_tmp(ii) = eta_test
					ii = ii + 1_idef
				else
					call RANDOM_NUMBER(rand_unif)
					if (rand_unif .LT. ratio) then
						g_tmp(ii) = g_test
						eta_tmp(ii) = eta_test
						ii = ii + 1_idef
					end if
				end if
			end do

			eta_tmp = ABS(eta_tmp)

			ii = 1_idef
			do while ( (ii.LT.nsamples).AND.(num_accepted.LT.nsamples) )
				lp = (g_tmp(ii).LE.h_params%max_sampling_g).AND.(g_tmp(ii).GE.h_params%min_sampling_g)
				if (lp) then
					num_accepted = num_accepted + 1_idef
					g_samples(num_accepted) = g_tmp(ii)
					eta_samples(num_accepted) = eta_tmp(ii)
				end if
				ii = ii + 1_idef
			end do
		end do

		if (TRIM(h_params%current_direction) .EQ. 'ANTICLOCKWISE') then
			eta_samples = 180.0_rp - eta_samples
		end if

		go = SUM(g_samples)/nsamples
		etao = SUM(eta_samples)/nsamples
	end if

	CALL MPI_SCATTER(g_samples,ppp,MPI_REAL8,g,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_SCATTER(eta_samples,ppp,MPI_REAL8,eta,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_BCAST(go,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	CALL MPI_BCAST(etao,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	if (params%mpi_params%rank.EQ.0_idef) then
		DEALLOCATE(g_samples)
		DEALLOCATE(eta_samples)
		DEALLOCATE(g_tmp)
		DEALLOCATE(eta_tmp)
	end if
END SUBROUTINE sample_Hollmann_distribution


SUBROUTINE save_params(params)
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
		filename = TRIM(params%path_to_outputs) // "experimental_distribution_parameters.h5"
		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		gname = "pdf_params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/max_pitch_angle"
		attr = "Maximum pitch angle in PDF (degrees)"
		call save_to_hdf5(h5file_id,dset,pdf_params%max_pitch_angle,attr)

		dset = TRIM(gname) // "/min_pitch_angle"
		attr = "Minimum pitch angle in PDF (degrees)"
		call save_to_hdf5(h5file_id,dset,pdf_params%min_pitch_angle,attr)

		dset = TRIM(gname) // "/min_energy"
		attr = "Minimum energy in PDF (eV)"
		units = 1.0_rp/C_E
		call save_to_hdf5(h5file_id,dset,units*pdf_params%min_energy,attr)

		dset = TRIM(gname) // "/max_energy"
		attr = "Maximum energy in PDF (eV)"
		units = 1.0_rp/C_E
		call save_to_hdf5(h5file_id,dset,units*pdf_params%max_energy,attr)

		dset = TRIM(gname) // "/max_p"
		attr = "Maximum momentum in PDF (me*c)"
		call save_to_hdf5(h5file_id,dset,pdf_params%max_p,attr)

		dset = TRIM(gname) // "/min_p"
		attr = "Maximum momentum in PDF (me*c)"
		call save_to_hdf5(h5file_id,dset,pdf_params%min_p,attr)

		dset = TRIM(gname) // "/Zeff"
		attr = "Effective atomic number of ions."
		call save_to_hdf5(h5file_id,dset,pdf_params%Zeff,attr)

		dset = TRIM(gname) // "/E"
		attr = "Parallel electric field in (Ec)"
		call save_to_hdf5(h5file_id,dset,pdf_params%E,attr)

		dset = TRIM(gname) // "/k"
		attr = "Shape factor"
		call save_to_hdf5(h5file_id,dset,pdf_params%k,attr)

		dset = TRIM(gname) // "/t"
		attr = "Scale factor"
		call save_to_hdf5(h5file_id,dset,pdf_params%t,attr)

		dset = TRIM(gname) // "/fGo"
		attr = "Normalization of Gamma function"
		call save_to_hdf5(h5file_id,dset,pdf_params%fGo,attr)

		dset = TRIM(gname) // "/lambda"
		attr = "Wavelength used when PDF is weighted with the distribution of synchrotron radiation."
		call save_to_hdf5(h5file_id,dset,pdf_params%lambda,attr)

		dset = TRIM(gname) // "/Bo"
		attr = "Magnetic field used when PDF is weighted with the distribution of synchrotron radiation."
		call save_to_hdf5(h5file_id,dset,pdf_params%Bo,attr)

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if
END SUBROUTINE save_params


SUBROUTINE save_Hollmann_params(params)
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
		filename = TRIM(params%path_to_outputs) // "experimental_distribution_parameters.h5"
		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		gname = "pdf_params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/max_pitch_angle"
		attr = "Maximum pitch angle in PDF (degrees)"
		call save_to_hdf5(h5file_id,dset,h_params%max_pitch_angle,attr)

		dset = TRIM(gname) // "/min_pitch_angle"
		attr = "Minimum pitch angle in PDF (degrees)"
		call save_to_hdf5(h5file_id,dset,h_params%min_pitch_angle,attr)

		dset = TRIM(gname) // "/min_energy"
		attr = "Minimum energy in PDF (eV)"
		units = 1.0_rp/C_E
		call save_to_hdf5(h5file_id,dset,units*h_params%min_energy,attr)

		dset = TRIM(gname) // "/max_energy"
		attr = "Maximum energy in PDF (eV)"
		units = 1.0_rp/C_E
		call save_to_hdf5(h5file_id,dset,units*h_params%max_energy,attr)

		dset = TRIM(gname) // "/max_g"
		attr = "Maximum momentum in PDF (me*c)"
		call save_to_hdf5(h5file_id,dset,h_params%max_g,attr)

		dset = TRIM(gname) // "/min_g"
		attr = "Maximum momentum in PDF (me*c)"
		call save_to_hdf5(h5file_id,dset,h_params%min_g,attr)

		dset = TRIM(gname) // "/Zeff"
		attr = "Effective atomic number of ions."
		call save_to_hdf5(h5file_id,dset,h_params%Zeff,attr)

		dset = TRIM(gname) // "/E"
		attr = "Parallel electric field in (Ec)"
		call save_to_hdf5(h5file_id,dset,h_params%E,attr)

		dset = TRIM(gname) // "/lambda"
		attr = "Wavelength used when PDF is weighted with the distribution of synchrotron radiation."
		call save_to_hdf5(h5file_id,dset,h_params%lambda,attr)

		dset = TRIM(gname) // "/Bo"
		attr = "Magnetic field used when PDF is weighted with the distribution of synchrotron radiation."
		call save_to_hdf5(h5file_id,dset,h_params%Bo,attr)

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if
END SUBROUTINE save_Hollmann_params


END MODULE korc_experimental_pdf
