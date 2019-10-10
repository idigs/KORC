MODULE korc_experimental_pdf
  USE korc_types
  USE korc_constants
  USE korc_HDF5
  USE korc_hpc
  USE special_functions
  use korc_coords
  use korc_rnd_numbers
  use korc_random

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
     REAL(rp) :: sigma_E
     REAL(rp) :: Zeff
     REAL(rp) :: sigma_Z
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
       get_Hollmann_distribution_3D,&
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
  END FUNCTION deg2rad


  FUNCTION rad2deg(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: rad2deg

    rad2deg = 180.0_rp*x/C_PI
  END FUNCTION rad2deg


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

    random_norm = mean+sigma*SQRT(-2.0_rp*LOG(rand1))*COS(2.0_rp*C_PI*rand2)
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
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: g
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: eta
    REAL(rp), INTENT(OUT) 					:: go
    REAL(rp), INTENT(OUT) 					:: etao
    REAL(rp) 					:: go_root
    REAL(rp) 					:: etao_root
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: p
    REAL(rp) 					:: p_buffer
    REAL(rp) 					:: p_test
    REAL(rp) 					:: eta_buffer
    REAL(rp) 					:: eta_test
    REAL(rp) 					:: ratio
    REAL(rp) 					:: rand_unif
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: p_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: eta_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: p_tmp
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: eta_tmp
    REAL(rp) 					:: minmax
    REAL(rp) 					:: min_p
    REAL(rp) 					:: max_p
    REAL(rp) 					:: min_pitch_angle
    REAL(rp) 					:: max_pitch_angle
    REAL(rp) 					:: deta
    REAL(rp) 					:: dp
    LOGICAL 					:: lp
    LOGICAL 					:: leta
    INTEGER 					:: num_accepted
    INTEGER 					:: ii
    INTEGER 					:: jj
    INTEGER 					:: ppp
    INTEGER 					:: nsamples
    INTEGER 					:: mpierr

    ppp = SIZE(g)
    nsamples = ppp*params%mpi_params%nmpi
    ALLOCATE(p(ppp))

    deta = (pdf_params%max_pitch_angle - &
         pdf_params%min_pitch_angle)/100.0_rp
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
       ALLOCATE(p_samples(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_samples(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(p_tmp(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_tmp(nsamples))
       ! Number of samples to distribute among all MPI processes


       !* * * Transient * * *!
       call RANDOM_SEED()

       call RANDOM_NUMBER(rand_unif)
       eta_buffer = pdf_params%min_pitch_angle + &
            (pdf_params%max_pitch_angle - pdf_params%min_pitch_angle)*rand_unif

       call RANDOM_NUMBER(rand_unif)
       p_buffer = pdf_params%min_p + (pdf_params%max_p - &
            pdf_params%min_p)*rand_unif

       ii=2_idef
       do while (ii .LE. 1000_idef)
          eta_test = eta_buffer + random_norm(0.0_rp,deta)
          do while ((ABS(eta_test) .GT. pdf_params%max_pitch_angle).OR. &
               (ABS(eta_test) .LT. pdf_params%min_pitch_angle))
             eta_test = eta_buffer + random_norm(0.0_rp,deta)
          end do

          p_test = p_buffer + random_norm(0.0_rp,dp)
          do while ((p_test.LT.pdf_params%min_p).OR.(p_test.GT. &
               pdf_params%max_p))
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
             do while ((ABS(eta_test) .GT. max_pitch_angle).OR. &
                  (ABS(eta_test) .LT. min_pitch_angle))
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
             lp = (p_tmp(ii).LE.pdf_params%max_p).AND.(p_tmp(ii).GE. &
                  pdf_params%min_p)
             leta = (eta_tmp(ii).LE.pdf_params%max_pitch_angle).AND. &
                  (eta_tmp(ii).GE.pdf_params%min_pitch_angle)
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


  SUBROUTINE get_Hollmann_distribution(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    TYPE(SPECIES),  INTENT(INOUT) :: spp
!    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
!    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
!    REAL(rp), INTENT(OUT) :: go
!    REAL(rp), INTENT(OUT) :: etao

    call initialize_Hollmann_params(params)

    call save_Hollmann_params(params)

    call sample_Hollmann_distribution(params,spp)


  END SUBROUTINE get_Hollmann_distribution

  SUBROUTINE get_Hollmann_distribution_3D(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    TYPE(SPECIES),  INTENT(INOUT) :: spp

    call initialize_Hollmann_params(params)

    call save_Hollmann_params(params)

    call normalize_Hollmann_params(params)

    call sample_Hollmann_distribution_3D(params,spp)

  END SUBROUTINE get_Hollmann_distribution_3D

  SUBROUTINE initialize_Hollmann_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    CHARACTER(MAX_STRING_LENGTH) :: filename
    CHARACTER(MAX_STRING_LENGTH) :: current_direction
    REAL(rp) :: E,sigma_E
    REAL(rp) :: Zeff,sigma_Z
    REAL(rp) :: max_pitch_angle
    REAL(rp) :: min_pitch_angle
    REAL(rp) :: max_energy
    REAL(rp) :: min_energy
    REAL(rp) :: Bo
    REAL(rp) :: lambda
    REAL(rp) :: A_fact

    NAMELIST /HollmannPDF/ E,Zeff,max_pitch_angle,min_pitch_angle,max_energy, &
         min_energy,filename,Bo,lambda,current_direction,A_fact,sigma_E, &
         sigma_Z


    open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
         status='OLD',form='formatted')
    read(default_unit_open,nml=HollmannPDF)
    close(default_unit_open)

    h_params%filename = TRIM(filename)

    h_params%E = E
    h_params%sigma_E = sigma_E
    h_params%Zeff = Zeff
    h_params%sigma_Z = sigma_Z
    h_params%max_pitch_angle = max_pitch_angle
    h_params%min_pitch_angle = min_pitch_angle

    h_params%min_sampling_energy = min_energy*C_E ! In Joules
    h_params%min_sampling_g = 1.0_rp + h_params%min_sampling_energy/ &
         (C_ME*C_C**2)
    h_params%max_sampling_energy = max_energy*C_E ! In Joules.
    h_params%max_sampling_g = 1.0_rp + h_params%max_sampling_energy/ &
         (C_ME*C_C**2)

    call load_data_from_hdf5()
    ! loads h_params%E_axis 1D energy range, h_params%fRE_E
    ! energy distribution as a function of h_params%E_axis,
    ! and h_params%fRE_pitch pitch angle distribution as a
    ! function of h_params%E_axis

    ALLOCATE(h_params%g(h_params%N))

    h_params%g = 1.0_rp + h_params%E_axis/(C_ME*C_C**2)
    ! 1D range of gamma based on energy range from Hollmann input file
    h_params%max_g = MAXVAL(h_params%g)
    h_params%min_g = MINVAL(h_params%g)

    h_params%current_direction = TRIM(current_direction)

    h_params%Bo = Bo
    h_params%lambda = lambda

    h_params%A_fact = A_fact
  END SUBROUTINE initialize_Hollmann_params

  subroutine normalize_Hollmann_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    h_params%sigma_E = h_params%sigma_E/params%cpp%length
    h_params%sigma_Z = h_params%sigma_Z/params%cpp%length
    
  end subroutine normalize_Hollmann_params

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
    ! index of gamma supplied to function in Hollmann input gamma range
    D = h_params%g(index) - g

    ! linear interpolation of Hollmann input gamma range to gamma supplied
    ! to function
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
    ! end of linear interpolation, fRE_H is evaluation of input Hollmann energy
    ! distribution PDF at gamma supplied to function

    A = (2.0_rp*h_params%E/(h_params%Zeff + 1.0_rp))*(g**2 - 1.0_rp)/g
    A = A*h_params%A_fact

    feta = A*EXP(-A*(1.0_rp - COS(deg2rad(eta))))/(1.0_rp - EXP(-2.0_rp*A))     ! MRC
    !	feta = 0.5_rp*A*EXP(A*COS(deg2rad(eta)))/SINH(A)                            ! MRC

    fRE_H = fRE_H*feta
  END FUNCTION fRE_H

    FUNCTION fRE_H_3D(eta,g,R,Z,R0,Z0)
    REAL(rp), INTENT(IN) 	:: eta ! pitch angle in degrees
    REAL(rp), INTENT(IN) 	:: g ! Relativistic gamma factor
    REAL(rp), INTENT(IN)  :: R,Z,R0,Z0
    REAL(rp) 				:: fRE_H_3D
    REAL(rp) 				:: D
    REAL(rp) 				:: g0
    REAL(rp) 				:: g1
    REAL(rp) 				:: f0
    REAL(rp) 				:: f1
    REAL(rp) 				:: m
    REAL(rp) 				:: feta
    REAL(rp) 				:: A
    REAL(rp) 				:: rm,E_G,Z_G
    INTEGER 				:: index

    index = MINLOC(ABS(h_params%g - g),1)
    ! index of gamma supplied to function in Hollmann input gamma range
    D = h_params%g(index) - g

    ! linear interpolation of Hollmann input gamma range to gamma supplied
    ! to function
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

    fRE_H_3D = f0 + m*(g - g0)
    ! end of linear interpolation, fRE_H is evaluation of input Hollman energy
    ! distribution PDF at gamma supplied to function

    rm=sqrt((R-R0)**2+(Z-Z0)**2)
    E_G=exp(-(rm/h_params%sigma_E)**2/2)
    Z_G=exp(-(rm/h_params%sigma_Z)**2/2)

!    write(6,'("rm: ",E17.10)') rm
    
    A = (2.0_rp*h_params%E*E_G/(h_params%Zeff*Z_G + 1.0_rp))*(g**2 - 1.0_rp)/g
    A = A*h_params%A_fact

    feta = A*EXP(-A*(1.0_rp - COS(deg2rad(eta))))/(1.0_rp - EXP(-2.0_rp*A))     ! MRC
    !	feta = 0.5_rp*A*EXP(A*COS(deg2rad(eta)))/SINH(A)                            ! MRC

    fRE_H_3D = fRE_H_3D*feta
  END FUNCTION fRE_H_3D

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


  SUBROUTINE sample_Hollmann_distribution(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN) 			:: params
    TYPE(SPECIES), INTENT(INOUT) 		:: spp
!    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: g
!    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: eta
!    REAL(rp), INTENT(OUT) 				:: go
!    REAL(rp), INTENT(OUT) 				:: etao
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: p
    REAL(rp) 						:: g_buffer
    REAL(rp) 						:: g_test
    REAL(rp) 						:: eta_buffer
    REAL(rp) 						:: eta_test
    REAL(rp) 						:: ratio
    REAL(rp) 						:: rand_unif
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: g_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: eta_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: g_tmp
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: eta_tmp
    REAL(rp) 						:: minmax
    REAL(rp) 						:: min_g
    REAL(rp) 						:: max_g
    REAL(rp) 						:: min_pitch_angle
    REAL(rp) 						:: max_pitch_angle
    REAL(rp) 						:: dg
    REAL(rp) 						:: deta
    LOGICAL 						:: lp
    INTEGER 						:: index_i
    INTEGER 						:: index_f
    INTEGER 						:: num_accepted
    INTEGER 						:: ii
    INTEGER 						:: jj
    INTEGER 						:: ppp
    INTEGER 						:: nsamples
    INTEGER 						:: mpierr
    INTEGER,DIMENSION(33) :: seed=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
    
    nsamples = spp%ppp*params%mpi_params%nmpi

    index_i = MINLOC(ABS(h_params%g - h_params%min_sampling_g),1)
    !index of minimum gamma range desired
    index_f = MINLOC(ABS(h_params%g - h_params%max_sampling_g),1)
    !index of maximum gamma range desired
    
    deta = (h_params%max_pitch_angle - h_params%min_pitch_angle)/100.0_rp
    dg = (h_params%max_sampling_g - h_params%min_sampling_g)/100.0_rp

    ! buffer at minimum gamma boundary  
    do jj=1_idef,INT(minmax_buffer_size,idef)
       minmax = h_params%min_sampling_g - REAL(jj,rp)*dg
       if (minmax.GT.h_params%min_g) then
          min_g = minmax
       end if
       ! buffer at maximum gamma boundary  
       minmax = h_params%max_sampling_g + REAL(jj,rp)*dg
       if (minmax.LT.h_params%max_g) then
          max_g = minmax
       end if
    end do

    ! buffer at minimum pitch angle boundary  
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

    ! buffer at maximum pitch angle boundary  
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
       !! MCMC and MH algorithm perfomred on single MPI process 
       !! to sample distribution function fRE_H
       ALLOCATE(g_samples(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_samples(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(g_tmp(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_tmp(nsamples))
       ! Number of samples to distribute among all MPI processes


       !Transient!

       if (.not.params%SameRandSeed) then
          call init_random_seed()
       else
          call random_seed(put=seed)
       end if

       call RANDOM_NUMBER(rand_unif)
!       rand_unif=get_random_U()
       eta_buffer = h_params%min_pitch_angle + (h_params%max_pitch_angle &
            - h_params%min_pitch_angle)*rand_unif

       call RANDOM_NUMBER(rand_unif)
!       rand_unif=get_random_U()
       g_buffer = h_params%min_sampling_g + (h_params%max_sampling_g - &
            h_params%min_sampling_g)*rand_unif

       ii=2_idef
       do while (ii .LE. 1000_idef)
          eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
!          eta_test = eta_buffer + get_random_N()*spp%dth
          do while ((ABS(eta_test) .GT. h_params%max_pitch_angle).OR. &
               (ABS(eta_test) .LT. h_params%min_pitch_angle))
             eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
!             eta_test = eta_buffer + get_random_N()*spp%dth
          end do

          g_test = g_buffer + random_norm(0.0_rp,spp%dgam)
!          g_test = g_buffer + get_random_N()*spp%dgam
          do while ((g_test.LT.h_params%min_sampling_g).OR. &
               (g_test.GT.h_params%max_sampling_g))
             g_test = g_buffer + random_norm(0.0_rp,spp%dgam)
!             g_test = g_buffer + get_random_N()*spp%dgam
          end do

          ratio = fRE_H(eta_test,g_test)*sin(deg2rad(eta_test))/ &
               (fRE_H(eta_buffer,g_buffer)* &
                  sin(deg2rad(eta_buffer)))
          !ratio = fRE_H(eta_test,g_test)/fRE_H(eta_buffer,g_buffer)
          !ratio = fRE_HxPR(eta_test,g_test)/fRE_HxPR(eta_buffer,g_buffer)

          if (ratio .GE. 1.0_rp) then
             g_buffer = g_test
             eta_buffer = eta_test
             ii = ii + 1_idef
          else
             call RANDOM_NUMBER(rand_unif)
!             rand_unif=get_random_U()
             if (rand_unif .LT. ratio) then
                g_buffer = g_test
                eta_buffer = eta_test
                ii = ii + 1_idef
             end if
          end if
       end do

       !Transient!

       eta_tmp(1) = eta_buffer
       g_tmp(1) = g_buffer

       num_accepted = 0_idef
       do while(num_accepted.LT.nsamples)
          ii=2_idef
          do while (ii .LE. nsamples)

             if (modulo(ii,nsamples/10).eq.0) then
                write(6,'("Sample: ",I10)') ii
             end if
             
!             write(6,'("iisample",I16)') ii
             eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,spp%dth)
             !eta_test = eta_tmp(ii-1) + get_random_N()*spp%dth
!             write(6,'("max_pitch_angle: ",E17.10)') max_pitch_angle
!             write(6,'("min_pitch_angle: ",E17.10)') min_pitch_angle
             do while ((ABS(eta_test) .GT. max_pitch_angle).OR. &
                  (ABS(eta_test) .LT. min_pitch_angle))
                eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,spp%dth)
!                eta_test = eta_tmp(ii-1) + get_random_N()*spp%dth
!                write(6,'("eta_test: ",E17.10)') eta_test
             end do

             g_test = g_tmp(ii-1) + random_norm(0.0_rp,spp%dgam)
             !g_test = g_tmp(ii-1) + get_random_N()*spp%dgam

!             write(6,'("max_g: ",E17.10)') max_g
!             write(6,'("min_g: ",E17.10)') min_g
             do while ((g_test.LT.min_g).OR.(g_test.GT.max_g))
                g_test = g_tmp(ii-1) + random_norm(0.0_rp,spp%dgam)
!                g_test = g_tmp(ii-1) + get_random_N()*spp%dgam
!                write(6,'("g_test: ",E17.10)') g_test
             end do

             ratio = fRE_H(eta_test,g_test)*sin(deg2rad(eta_test))/ &
                  (fRE_H(eta_tmp(ii-1),g_tmp(ii-1))* &
                  sin(deg2rad(eta_tmp(ii-1))))
!             ratio = fRE_H(eta_test,g_test)/fRE_H(eta_tmp(ii-1),g_tmp(ii-1))
             !ratio = fRE_HxPR(eta_test,g_test)/fRE_HxPR(eta_tmp(ii-1),g_tmp(ii-1))
!             write(6,'("ratio: ",E17.10)') ratio

             if (ratio .GE. 1.0_rp) then
                g_tmp(ii) = g_test
                eta_tmp(ii) = eta_test
                ii = ii + 1_idef
             else
                call RANDOM_NUMBER(rand_unif)
!                rand_unif=get_random_U()
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
!             write(6,'("iiaccept",I16)') ii
             lp = (g_tmp(ii).LE.h_params%max_sampling_g).AND. &
                  (g_tmp(ii).GE.h_params%min_sampling_g).AND. &
                  (eta_tmp(ii).LE.h_params%max_pitch_angle).AND. &
                  (eta_tmp(ii).GE.h_params%min_pitch_angle)             
             if (lp) then
                num_accepted = num_accepted + 1_idef
                g_samples(num_accepted) = g_tmp(ii)
                eta_samples(num_accepted) = eta_tmp(ii)
             end if
             ii = ii + 1_idef
          end do
       end do


      !		if (TRIM(h_params%current_direction) .EQ. 'ANTICLOCKWISE') then
      !			eta_samples = 180.0_rp - eta_samples
      !		end if

!       go = SUM(g_samples)/nsamples
!       etao = SUM(eta_samples)/nsamples
    end if !MCMC computed on single MPI process

    CALL MPI_SCATTER(g_samples,spp%ppp,MPI_REAL8,spp%vars%g,spp%ppp,MPI_REAL8, &
         0,MPI_COMM_WORLD,mpierr)

    CALL MPI_SCATTER(eta_samples,spp%ppp,MPI_REAL8,spp%vars%eta,spp%ppp,MPI_REAL8, &
         0,MPI_COMM_WORLD,mpierr)

!    CALL MPI_BCAST(go,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

!    CALL MPI_BCAST(etao,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if (params%mpi_params%rank.EQ.0_idef) then
       DEALLOCATE(g_samples)
       DEALLOCATE(eta_samples)
       DEALLOCATE(g_tmp)
       DEALLOCATE(eta_tmp)
    end if

 !   write(6,'("sampled eta: ",E17.10)') eta
    
  END SUBROUTINE sample_Hollmann_distribution

FUNCTION PSI_ROT_exp(R,R0,sigR,Z,Z0,sigZ,theta)
  REAL(rp), INTENT(IN) 	:: R
  !! R-coordinate of MH sampled location
  REAL(rp), INTENT(IN) 	:: R0
  !! R-coordinate of center of 2D Gaussian
  REAL(rp), INTENT(IN) 	:: sigR
  !! Variance of first dimension of 2D Gaussian
  REAL(rp), INTENT(IN) 	:: Z
  !! Z-coordinate of MH sampled location
  REAL(rp), INTENT(IN) 	:: Z0
  !! Z-coordinate of center of 2D Gaussian
  REAL(rp), INTENT(IN) 	:: sigZ
  !! Variance of second dimension of 2D Gaussian
  REAL(rp), INTENT(IN) 	:: theta
  !! Angle of counter-clockwise rotation (in radians), of 2D Gaussian
  !! distribution relative to R,Z
  REAL(rp) 		:: PSI_ROT_exp
  !! Argument of exponential comprising 2D Gaussian distribution 

  PSI_ROT_exp=(R-R0)**2*((cos(theta))**2/(2*sigR**2)+ &
       (sin(theta))**2/(2*sigZ**2))+ &
       2*(R-R0)*(Z-Z0)*cos(theta)*sin(theta)*(1/(2*sigR**2)-1/(2*sigZ**2))+ &
       (Z-Z0)**2*((sin(theta))**2/(2*sigR**2)+(cos(theta))**2/(2*sigZ**2))
  
END FUNCTION PSI_ROT_exp


FUNCTION indicator_exp(psi,psi_max)
  REAL(rp), INTENT(IN)  :: psi
  REAL(rp), INTENT(IN)  :: psi_max
  REAL(rp)              :: indicator_exp

  IF (psi.LT.psi_max) THEN
     indicator_exp=1
  ELSE
     indicator_exp=0
  END IF
  
END FUNCTION indicator_exp
  
  subroutine sample_Hollmann_distribution_3D(params,spp)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.

  !! An instance of the KORC derived type FIELDS.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples,X_samples,Y_samples
  !! Major radial location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PHI_samples
  !! Azimuithal angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !! Vertical location of all samples

  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: G_samples
  !! Gamma of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: eta_samples
  !! Pitch angle of all samples

  REAL(rp) 				:: psi_max_buff
  !! Value of buffer above desired maximum argument of 2D Gaussian spatial
  !! profile
  REAL(rp) 				:: minmax
  !! Temporary variable used for setting buffers

  !! Minimum domain for momentum sampling including buffer
  REAL(rp) 				:: max_pitch_angle
  !! Maximum domain for pitch angle sampling including buffer
  REAL(rp) 				:: min_pitch_angle
  !! Minimum domain for pitch angle sampling including buffer
  REAL(rp) 				:: min_g,max_g
  REAL(rp) 				:: theta_rad
  !! Angle of rotation of 2D Gaussian spatial distribution in radians
  REAL(rp) 				:: R_buffer
  !! Previous sample of R location
  REAL(rp) 				:: Z_buffer
  !! Previous sample of Z location
 
  REAL(rp) 				:: eta_buffer
  !! Previous sample of pitch
  REAL(rp) 				:: G_buffer
  REAL(rp) 				:: R_test
  !! Present sample of R location
  REAL(rp) 				:: Z_test
  !! Present sample of Z location

  REAL(rp) 				:: eta_test
  !! Present sample of pitch angle
  REAL(rp) 				:: G_test
  REAL(rp) 				:: psi0
  !! Previous value of 2D Gaussian argument based on R_buffer, Z_buffer
  REAL(rp) 				:: psi1
  !! Present value of 2D Gaussian argument based on R_test, Z_test
  REAL(rp) 				:: f0
  !! Evaluation of Avalanche distribution with previous sample
  REAL(rp) 				:: f1
  !! Evaluation of Avalanche distribution with present sample
  REAL(rp) 				:: rand_unif
  !! Uniform random variable [0,1]
  REAL(rp) 				:: ratio
  !! MH selection criteria
  INTEGER				:: nsamples
  !! Total number of samples to be distributed over all mpi processes
  INTEGER 				:: ii
  !! Sample iterator.
  INTEGER 				:: mpierr
  !! mpi error indicator
  REAL(rp) 						:: dg,deta
  LOGICAL :: accepted
  INTEGER,DIMENSION(33) :: seed=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
  
  
  nsamples = spp%ppp*params%mpi_params%nmpi

  psi_max_buff = spp%psi_max*1.1_rp

  theta_rad=C_PI*spp%theta_gauss/180.0_rp


  deta = (h_params%max_pitch_angle - h_params%min_pitch_angle)/100.0_rp
  dg = (h_params%max_sampling_g - h_params%min_sampling_g)/100.0_rp
  
  ! buffer at minimum gamma boundary  
  do ii=1_idef,INT(minmax_buffer_size,idef)
     minmax = h_params%min_sampling_g - REAL(ii,rp)*dg
     if (minmax.GT.h_params%min_g) then
        min_g = minmax
     end if
     ! buffer at maximum gamma boundary  
     minmax = h_params%max_sampling_g + REAL(ii,rp)*dg
     if (minmax.LT.h_params%max_g) then
        max_g = minmax
     end if
  end do

  ! buffer at minimum pitch angle boundary  
  if (h_params%min_pitch_angle.GE.korc_zero) then
     do ii=1_idef,INT(minmax_buffer_size,idef)
        minmax = h_params%min_pitch_angle -  REAL(ii,rp)*deta
        if (minmax.GT.0.0_rp) then
           min_pitch_angle = minmax
        end if
     end do
  else
     min_pitch_angle = h_params%min_pitch_angle 
  end if

  ! buffer at maximum pitch angle boundary  
  do ii=1_idef,INT(minmax_buffer_size,idef)
     minmax = h_params%max_pitch_angle + REAL(ii,rp)*deta
     if (minmax.LE.180.0_rp) then
        max_pitch_angle = minmax
     end if
  end do

  
  if (params%mpi_params%rank.EQ.0_idef) then
     ALLOCATE(R_samples(nsamples))
     ALLOCATE(X_samples(nsamples))
     ALLOCATE(Y_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(PHI_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(Z_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(eta_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(G_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     
     ! Transient !

     R_buffer = spp%Ro
     Z_buffer = spp%Zo


     if (.not.params%SameRandSeed) then
        call init_random_seed()
     else
        call random_seed(put=seed)
     end if
     
!     call RANDOM_NUMBER(rand_unif)
!     eta_buffer = min_pitch_angle + (max_pitch_angle &
!          - min_pitch_angle)*rand_unif
!     eta_buffer = min_pitch_angle + (max_pitch_angle &
!          - min_pitch_angle)*get_random_mkl_U()
     eta_buffer = min_pitch_angle + (max_pitch_angle &
          - min_pitch_angle)*get_random_U()
     
!     call RANDOM_NUMBER(rand_unif)
!     G_buffer = min_g + (max_g - min_g)*rand_unif
!     G_buffer = min_g + (max_g - min_g)*get_random_mkl_U()
     G_buffer = min_g + (max_g - min_g)*get_random_U()
     
!     write(6,*) 'R_buffer',R_buffer
!     write(6,*) 'Z_buffer',Z_buffer
!     write(6,*) 'eta_buffer',eta_buffer
!     write(6,*) 'G_buffer',G_buffer
     
     !     write(6,'("length norm: ",E17.10)') params%cpp%length
     
     accepted=.false.
     ii=1_idef
     do while (ii .LE. 1000_idef)

!        write(6,'("burn:",I15)') ii
        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ
        
        !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
        !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
        eta_test = eta_buffer + get_random_N()*spp%dth
        
        !G_test = G_buffer + random_norm(0.0_rp,spp%dgam)
        !G_test = G_buffer + get_random_mkl_N(0.0_rp,spp%dgam)
        G_test = G_buffer + get_random_N()*spp%dgam


        ! Test that pitch angle and momentum are within chosen boundary
        do while ((ABS(eta_test) .GT. max_pitch_angle).OR. &
             (ABS(eta_test) .LT. min_pitch_angle))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           eta_test = eta_buffer + get_random_N()*spp%dth
        end do

        do while ((G_test.LT.min_g).OR.(G_test.GT.max_g))
           !G_test = G_buffer + random_norm(0.0_rp,spp%dgam)
           !G_test = G_buffer + get_random_mkl_N(0.0_rp,spp%dgam)
           G_test = G_buffer + get_random_N()*spp%dgam
        end do
        
        ! initialize 2D gaussian argument and distribution function, or
        ! copy from previous sample
        if (ii==1) then
           psi0=PSI_ROT_exp(R_buffer,spp%Ro,spp%sigmaR,Z_buffer,spp%Zo, &
                spp%sigmaZ,theta_rad)
           
           f0=fRE_H_3D(eta_buffer,G_buffer,R_buffer,Z_buffer,spp%Ro,spp%Zo)
!           f0=fRE_H(eta_buffer,G_buffer)
        end if

        if (accepted) then
           psi0=psi1
           f0=f1
        end if
        
        psi1=PSI_ROT_exp(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
             spp%sigmaZ,theta_rad)       
        
                
        f1=fRE_H_3D(eta_test,G_test,R_test,Z_test,spp%Ro,spp%Zo)    
!        f1=fRE_H(eta_test,G_test)
        
!        write(6,'("psi0: ",E17.10)') psi0
!        write(6,'("psi1: ",E17.10)') psi1

!        write(6,'("f0: ",E17.10)') f0
!        write(6,'("f1: ",E17.10)') f1

        
        ! Calculate acceptance ratio for MH algorithm. fRE function
        ! incorporates p^2 factor of spherical coordinate Jacobian
        ! for velocity phase space, factors of sin(pitch angle) for velocity
        ! phase space and cylindrical coordinate Jacobian R for spatial
        ! phase space incorporated here.
        ratio = indicator_exp(psi1,spp%psi_max)* &
             R_test*EXP(-psi1)*f1*sin(deg2rad(eta_test))/ &
             (R_buffer*EXP(-psi0)*f0*sin(deg2rad(eta_buffer)))

!        ratio = f1*sin(deg2rad(eta_test))/(f0*sin(deg2rad(eta_buffer)))

        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           eta_buffer = eta_test
           G_buffer = G_test
           ii = ii + 1_idef
        else
!           call RANDOM_NUMBER(rand_unif)
!           if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              eta_buffer = eta_test
              G_buffer = G_test
              ii = ii + 1_idef
           end if
        end if
     end do
     ! Transient !

     ii=1_idef
     do while (ii .LE. nsamples)

!        write(6,'("sample:",I15)') ii
        
       if (modulo(ii,nsamples/10).eq.0) then
           write(6,'("Sample: ",I10)') ii
        end if
        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ
        
        !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
        !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
        eta_test = eta_buffer + get_random_N()*spp%dth
        
        !G_test = G_buffer + random_norm(0.0_rp,spp%dgam)
        !G_test = G_buffer + get_random_mkl_N(0.0_rp,spp%dgam)
        G_test = G_buffer + get_random_N()*spp%dgam


        ! Test that pitch angle and momentum are within chosen boundary
        do while ((ABS(eta_test) .GT. max_pitch_angle).OR. &
             (ABS(eta_test) .LT. min_pitch_angle))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           eta_test = eta_buffer + get_random_N()*spp%dth
        end do

        do while ((G_test.LT.min_g).OR.(G_test.GT.max_g))
           !G_test = G_buffer + random_norm(0.0_rp,spp%dgam)
           !G_test = G_buffer + get_random_mkl_N(0.0_rp,spp%dgam)
           G_test = G_buffer + get_random_N()*spp%dgam
        end do

        if (accepted) then
           psi0=psi1
           f0=f1
        end if
        
        psi1=PSI_ROT_exp(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
             spp%sigmaZ,theta_rad)

!        write(6,'("R: ",E17.10)') R_test
!        write(6,'("R0: ",E17.10)') spp%Ro
!        write(6,'("sigma_R: ",E17.10)') spp%sigmaR
!        write(6,'("dR: ",E17.10)') spp%dR
!        write(6,'("N_dR: ",E17.10)') random_norm(0.0_rp,spp%dR)
!        write(6,'("Z: ",E17.10)') Z_test
!        write(6,'("Z0: ",E17.10)') spp%Zo
!        write(6,'("sigma_Z: ",E17.10)') spp%sigmaZ
!        write(6,'("dZ: ",E17.10)') spp%dZ
!        write(6,'("N_dR: ",Z17.10)') random_norm(0.0_rp,spp%dZ)
        
        f1=fRE_H_3D(eta_test,G_test,R_test,Z_test,spp%Ro,spp%Zo)            
!        f1=fRE_H(eta_test,G_test)
        
        ratio = indicator_exp(psi1,psi_max_buff)* &
             R_test*EXP(-psi1)*f1*sin(deg2rad(eta_test))/ &
             (R_buffer*EXP(-psi0)*f0*sin(deg2rad(eta_buffer)))

!        ratio = f1*sin(deg2rad(eta_test))/(f0*sin(deg2rad(eta_buffer)))
        
        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           eta_buffer = eta_test
           G_buffer = G_test
        else
           !call RANDOM_NUMBER(rand_unif)
           !if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              eta_buffer = eta_test
              G_buffer = G_test
           end if
        end if

!        write(6,'("R: ",E17.10)') R_buffer
!        write(6,'("Z: ",E17.10)') Z_buffer
        
        ! Only accept sample if it is within desired boundary, but
        ! add to MC above if within buffer. This helps make the boundary
        ! more defined.
        IF ((INT(indicator_exp(psi1,spp%psi_max)).EQ.1).AND. &
             (G_buffer.LE.h_params%max_sampling_g).AND. &
             (G_buffer.GE.h_params%min_sampling_g).AND. &
             (eta_buffer.LE.h_params%max_pitch_angle).AND. &
             (eta_buffer.GE.h_params%min_pitch_angle).AND. &
             ACCEPTED) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
           eta_samples(ii) = eta_buffer
           G_samples(ii) = G_buffer

!           write(6,*) 'RS',R_buffer
           
           ! Sample phi location uniformly
           !call RANDOM_NUMBER(rand_unif)
           !PHI_samples(ii) = 2.0_rp*C_PI*rand_unif
           !PHI_samples(ii) = 2.0_rp*C_PI*get_random_mkl_U()
           PHI_samples(ii) = 2.0_rp*C_PI*get_random_U()
           ii = ii + 1_idef 
        END IF
        
     end do

!  if (minval(R_samples(:)).lt.1._rp/params%cpp%length) stop 'error with sample'
!  write(6,'("R_sample: ",E17.10)') R_samples(:)*params%cpp%length

     X_samples=R_samples*cos(PHI_samples)
     Y_samples=R_samples*sin(PHI_samples)

!     write(6,*) 'R_samples',R_samples
!     write(6,*) 'PHI_samples',PHI_samples
!     write(6,*) 'Z_samples',Z_samples
!     write(6,*) 'G_samples',G_samples
!     write(6,*) 'eta_samples',eta_samples
     
  end if


  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)


  CALL MPI_SCATTER(X_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,1),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Y_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,2),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Z_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,3),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(eta_samples,spp%ppp,MPI_REAL8, &
       spp%vars%eta,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(G_samples,spp%ppp,MPI_REAL8, &
       spp%vars%g,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!  write(6,*) params%mpi_params%rank,'varX',spp%vars%X(:,1)
!  write(6,*) params%mpi_params%rank,'varY',spp%vars%X(:,2)
  
!  write(6,'("X_X: ",E17.10)') spp%vars%X(:,1)*params%cpp%length
  
  ! gamma is kept for each particle, not the momentum

  if (params%orbit_model(1:2).eq.'GC') call cart_to_cyl(spp%vars%X,spp%vars%Y)

!  write(6,*) params%mpi_params%rank,'varX',spp%vars%X(:,1)
!  write(6,*) params%mpi_params%rank,'varR',spp%vars%Y(:,1)
  

!  write(6,'("Y_R: ",E17.10)') spp%vars%Y(:,1)*params%cpp%length
!  write(6,'("Y_PHI: ",E17.10)') spp%vars%Y(:,1)*params%cpp%length
!  write(6,'("Y_Z: ",E17.10)') spp%vars%Y(:,3)*params%cpp%length
  
!  if (minval(spp%vars%Y(:,1)).lt.1._rp/params%cpp%length) stop 'error with avalanche'
  
  if (params%mpi_params%rank.EQ.0_idef) then
     DEALLOCATE(R_samples)
     DEALLOCATE(X_samples)
     DEALLOCATE(Y_samples)
     DEALLOCATE(Z_samples)
     DEALLOCATE(PHI_samples)
     DEALLOCATE(eta_samples)
     DEALLOCATE(G_samples)
  end if
  
  
end subroutine sample_Hollmann_distribution_3D

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

       dset = TRIM(gname) // "/max_sampling_g"
       attr = "Maximum momentum in PDF (me*c)"
       call save_to_hdf5(h5file_id,dset,h_params%max_sampling_g,attr)

       dset = TRIM(gname) // "/min_sampling_g"
       attr = "Maximum momentum in PDF (me*c)"
       call save_to_hdf5(h5file_id,dset,h_params%min_sampling_g,attr)
       
       dset = TRIM(gname) // "/Zeff"
       attr = "Effective atomic number of ions."
       call save_to_hdf5(h5file_id,dset,h_params%Zeff,attr)

       dset = TRIM(gname) // "/sigmaZeff"
       attr = "Effective atomic number of ions."
       call save_to_hdf5(h5file_id,dset,h_params%sigma_Z,attr)
       
       dset = TRIM(gname) // "/E"
       attr = "Parallel electric field in (Ec)"
       call save_to_hdf5(h5file_id,dset,h_params%E,attr)

       dset = TRIM(gname) // "/sigmaE"
       attr = "Parallel electric field in (Ec)"
       call save_to_hdf5(h5file_id,dset,h_params%sigma_E,attr)
       
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