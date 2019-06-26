!> @brief Module containing subroutines to initialize the velocity distribution of the simulated particles.
MODULE korc_velocity_distribution
    USE korc_types
    USE korc_constants
    USE korc_HDF5
    USE korc_hpc
    use korc_fields
    use korc_rnd_numbers
    use korc_hammersley_generator

    use korc_avalanche
    use korc_experimental_pdf
    use korc_energy_pdfs
    use korc_simple_equilibrium_pdf

    IMPLICIT NONE

    PUBLIC :: initial_gyro_distribution,&
                thermal_distribution,&
                initial_energy_pitch_dist
    PRIVATE :: fth_3V,&
                random_norm,&
                gyro_distribution

    CONTAINS

!> @brief Function used to sample the probability density function of a thermal plasma in the 3-dimensional velocity space.
!! @details This function returns @f$f_{T_e}(v) = \exp{\left( v^2/2v_{T_e}^2 \right)}@f$, where @f$v_{T_e} = \sqrt{T_e/m_e}@f$ is
!! the temperature of the thermal electrons, and @f$v = |\mathbf{v}|@f$ is the speed of the sampled electron.
!!
!! @param[in] V Velocity of the sampled electron @f$\mathbf{v}@f$.
!! @param[in] Vth Thermal velocity of the background electrons @f$v_{T_e}@f$.
!! @param fth_3V Value of @f$f_{T_e}(v)@f$.
FUNCTION fth_3V(Vth,V)
    REAL(rp), DIMENSION(3), INTENT(IN)     :: V
    REAL(rp), INTENT(IN)                 :: Vth
    REAL(rp)                             :: fth_3V

    fth_3V = EXP(-0.5_rp*DOT_PRODUCT(V,V)/Vth**2.0_rp)
END FUNCTION fth_3V


!> @brief Gaussian random number generator.
!! @details This function returns a deviate of a Gaussian distribution @f$f_G(x;\mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}} \exp{\left( -(x-\mu)^2/2\sigma^2 \right)}@f$,
!! with mean @f$\mu@f$, and standard deviation @f$\sigma@f$.
!!
!! We use the Inverse Transform Sampling Method for sampling @f$x@f$. With this method we get @f$x = \sqrt{-2\log{(1-y)}}\cos(2\pi z)@f$,
!! where @f$y@f$ and @f$z@f$ are uniform random numbers in the interval @f$[0,1]@f$.
!!
!! @param[in] mu Mean value @f$\mu@f$ of the Gaussian distribution.
!! @param[in] mu Standard deviation @f$\sigma@f$ of the Gaussian distribution.
!! @param random_norm Sampled number @f$x@f$ from the Gaussian distribution @f$f_G(x;\mu,\sigma)@f$.
!! @param rand1 Uniform random number in the interval @f$[0,1]@f$.
!! @param rand2 Uniform random number in the interval @f$[0,1]@f$.
FUNCTION random_norm(mu,sigma)
    REAL(rp), INTENT(IN)     :: mu
    REAL(rp), INTENT(IN)     :: sigma
    REAL(rp)                 :: random_norm
    REAL(rp)                 :: rand1
    REAL(rp)                 :: rand2

    call RANDOM_NUMBER(rand1)
    call RANDOM_NUMBER(rand2)

    random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
END FUNCTION random_norm


!> @brief Subroutine that samples a thermal distribution function of electrons for generating the initial condition of a set of simulated particles.
!! @details This subroutine uses the Inverse Transform Sampling Method along with the Metropolis-Hastings algorithm to generate an
!! initial condition of the velocity distribution that follows a 3-dimensional (in velocity space) thermal distribution.
!!
!! @todo Check that the gyro-distribution is initialized right in this function.
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param Vmax Velocity cutoff where we stop sampling the tail of the thermal distribution.
!! @param Vth Thermal velocity of the sampled distribution @f$v_{T_e} = \sqrt{T_e/m_e}@f$.
!! @param sv Step to sample the velocity space using the Metropolis-Hastings algorithm.
!! @param ratio Ratio used to accept or reject a sampling in the Metropolis-Hastings algorithm.
!! @param rand_unif Uniform random deviate in the interval  @f$[0,1]@f$.
!! @param V Sampled velocity.
!! @param U Sampled velocity.
!! @param b Temporary variable representing a unit vector along the @f$x@f$-axis.
!! @param ii Iterator.
!! @param ppp Number of particles per species.
subroutine thermal_distribution(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)     :: params
    TYPE(SPECIES), INTENT(INOUT)     :: spp
    REAL(rp)                         :: Vmax
    REAL(rp)                         :: Vth
    REAL(rp)                         :: sv
    REAL(rp)                         :: ratio
    REAL(rp)                         :: rand_unif
    REAL(rp), DIMENSION(3)             :: V
    REAL(rp), DIMENSION(3)             :: U
    REAL(rp), DIMENSION(3)             :: b = (/1.0_rp,0.0_rp,0.0_rp/)
    INTEGER                         :: ii
    INTEGER                         :: ppp

    Vmax = 0.9_rp
    Vth = SQRT(spp%Eo*ABS(spp%q)/spp%m)
    ppp = spp%ppp

    V = (/0.0_rp,0.0_rp,0.0_rp/)
    sv = Vth/10.0_rp

    ii=2_idef
    do while (ii .LE. 1000_idef)
        U(1) = V(1) + random_norm(0.0_rp,sv)
        do while (ABS(U(1)) .GT. Vmax)
            U(1) = V(1) + random_norm(0.0_rp,sv)
        end do

        U(2) = V(2) + random_norm(0.0_rp,sv)
        do while (ABS(U(2)) .GT. Vmax)
            U(2) = V(2) + random_norm(0.0_rp,sv)
        end do

        U(3) = V(3) + random_norm(0.0_rp,sv)
        do while (ABS(U(3)) .GT. Vmax)
            U(3) = V(3) + random_norm(0.0_rp,sv)
        end do

        ratio = fth_3V(Vth,U)/fth_3V(Vth,V)

        if (ratio .GE. 1.0_rp) then
            V = U
            ii = ii + 1_idef
        else
            call RANDOM_NUMBER(rand_unif)
            if (ratio .GT. rand_unif) then
                V = U
                ii = ii + 1_idef
            end if
        end if
    end do

    spp%vars%V(:,1) = V
    ii=2_idef
    do while (ii .LE. ppp)
        U(1) = spp%vars%V(1,ii-1) + random_norm(0.0_rp,sv)
        do while (ABS(U(1)) .GT. Vmax)
            U(1) = spp%vars%V(1,ii-1) + random_norm(0.0_rp,sv)
        end do
        U(2) = spp%vars%V(2,ii-1) + random_norm(0.0_rp,sv)
        do while (ABS(U(2)) .GT. Vmax)
            U(2) = spp%vars%V(2,ii-1) + random_norm(0.0_rp,sv)
        end do
        U(3) = spp%vars%V(3,ii-1) + random_norm(0.0_rp,sv)
        do while (ABS(U(3)) .GT. Vmax)
            U(3) = spp%vars%V(3,ii-1) + random_norm(0.0_rp,sv)
        end do

        ratio = fth_3V(Vth,U)/fth_3V(Vth,spp%vars%V(:,ii-1))

        if (ratio .GE. 1.0_rp) then
            spp%vars%V(:,ii) = U
            ii = ii + 1_idef
        else
            call RANDOM_NUMBER(rand_unif)
            if (ratio .GT. rand_unif) then
                spp%vars%V(:,ii) = U
                ii = ii + 1_idef
            end if
        end if
    end do

    do ii=1_idef,ppp
        spp%vars%g(ii) = 1.0_rp/SQRT(1.0_rp - SUM(spp%vars%V(:,ii)**2,1))
        spp%vars%eta(ii) = ACOS(DOT_PRODUCT(b,spp%vars%V(:,ii)/SQRT(SUM(spp%vars%V(:,ii)**2,1))))
    end do

    spp%go = spp%Eo/(spp%m*C_C**2)
    spp%etao = 90.0_rp
end subroutine thermal_distribution


!> @brief Subroutine that calls subroutines of different modules to initialize the energy and pitch-angle distribution in various ways.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param ii Species iterator.
!! @param mpierr MPI error status.
subroutine initial_energy_pitch_dist(params,spp)
TYPE(KORC_PARAMS), INTENT(IN)                                 :: params
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    INTEGER                                                 :: ii
    INTEGER                                                 :: mpierr

    do ii=1_idef,params%num_species
        SELECT CASE (TRIM(spp(ii)%energy_distribution))
            CASE ('MONOENERGETIC')
                spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

                spp(ii)%vars%g = spp(ii)%go ! Monoenergetic
                spp(ii)%Eo_lims = (/spp(ii)%Eo, spp(ii)%Eo /)
            CASE ('THERMAL')
                call thermal_distribution(params,spp(ii))

                spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2, &
                                    spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
            CASE ('AVALANCHE')
                call get_avalanche_distribution(params,spp(ii)%vars%g,spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)

                spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
                spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2, &
                                    spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
            CASE ('HOLLMANN')
                call get_Hollmann_distribution(params,spp(ii)%vars%g,spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)

                spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
                spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2, &
                                    spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
            CASE ('EXPERIMENTAL-GAMMA')
                call get_experimentalG_distribution(params,spp(ii)%vars%g,spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)

                spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
                spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2, &
                                    spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
            CASE ('GAMMA')
                call get_gamma_distribution(params,spp(ii)%vars%g,spp(ii)%go)

                spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
                spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 , &
                                    spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
            CASE ('UNIFORM')
                spp(ii)%Eo = spp(ii)%Eo_lims(1)
                spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

                call generate_2D_hammersley_sequence(params%mpi_params%rank,params%mpi_params%nmpi,spp(ii)%vars%g,spp(ii)%vars%eta)

                spp(ii)%vars%g = (spp(ii)%Eo_lims(2) - spp(ii)%Eo_lims(1))*spp(ii)%vars%g/(spp(ii)%m*C_C**2) + &
                                    (spp(ii)%Eo_lims(1) + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
            CASE DEFAULT
                ! Something to be done
        END SELECT

        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

        SELECT CASE (TRIM(spp(ii)%pitch_distribution))
            CASE ('MONOPITCH')
                spp(ii)%vars%eta = spp(ii)%etao ! Mono-pitch-angle
                spp(ii)%etao_lims = (/spp(ii)%etao , spp(ii)%etao/)
            CASE ('THERMAL')
                spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
            CASE ('AVALANCHE')
                spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
            CASE ('HOLLMANN')
                spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
            CASE ('EXPERIMENTAL-GAMMA')
                spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
            CASE ('UNIFORM')
                spp(ii)%etao = spp(ii)%etao_lims(1)

                spp(ii)%vars%eta = (spp(ii)%etao_lims(2) - spp(ii)%etao_lims(1))*spp(ii)%vars%eta + spp(ii)%etao_lims(1)
            CASE ('SIMPLE-EQUILIBRIUM')
                call get_equilibrium_distribution(params,spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)

                spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), MAXVAL(spp(ii)%vars%eta)/)
            CASE DEFAULT
                ! Something to be done
        END SELECT

        if (params%mpi_params%rank .EQ. 0) then
            write(6,'(/,"* * * * * SPECIES: ",I2," * * * * * * * * * * *")') ii
            write(6,'("Energy distribution is: ",A20)') TRIM(spp(ii)%energy_distribution)
            write(6,'("Pitch-angle distribution is: ",A20)') TRIM(spp(ii)%pitch_distribution)
            write(6,'("* * * * * * * * * * * * * * * * * * * * * *",/)')
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    end do
end subroutine initial_energy_pitch_dist


!> @brief Subroutine that initializes the gyro-angle distribution of the particles.
!! @details When evolving the particles in the 6-D phase space, in addition to the position (3 degrees of freedom), energy (one degree of freedom),
!! pitch angle (one degree of freedom), we need to define the gyro-angle of the particle (one degree of freedom), which is given by the
!! pitch angle and the direction of the local magnetic field. By default, this subroutine generates a uniform gyro-angle distribution.
!!
!! @note Notice that all the simulation variables are normalized here.
!! @param[in] params Core KORC simulation parameters.
!! @param[in] F An instance of the KORC derived type FIELDS. This structure has the information of the magnetic field.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param[in,out] b1 Basis vector pointing along the local magnetic field, that is, along @f$\mathbf{b} = \mathbf{B}/B@f$.
!! @param[in,out] b2 Basis vector perpendicular to b1
!! @param[in,out] b3 Basis vector perpendicular to b1 and b2.
!! @param Vo Initial particle speed.
!! @param V1 Velocity component along b1.
!! @param V2 Velocity component along b2.
!! @param V3 Velocity component along b3.
!! @param theta Uniform random number in the interval @f$[0,2\pi]@f$ representing the gyro-angle.
!! @param x Unitary vector along the @f$x@f$-axis.
!! @param y Unitary vector along the @f$y@f$-axis.
!! @param z Unitary vector along the @f$z@f$-axis.
!! @param jj Particle iterator.
subroutine gyro_distribution(params,F,spp)
    TYPE(KORC_PARAMS), INTENT(IN)         :: params
    TYPE(FIELDS), INTENT(IN)              :: F
    TYPE(SPECIES), INTENT(INOUT)          :: spp
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: b1
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: b2
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: b3
    REAL(rp), DIMENSION(:), ALLOCATABLE   :: Vo
    REAL(rp), DIMENSION(:), ALLOCATABLE   :: V1
    REAL(rp), DIMENSION(:), ALLOCATABLE   :: V2
    REAL(rp), DIMENSION(:), ALLOCATABLE   :: V3
    REAL(rp), DIMENSION(:), ALLOCATABLE   :: theta
    REAL(rp), DIMENSION(3)                :: x = (/1.0_rp,0.0_rp,0.0_rp/)
    REAL(rp), DIMENSION(3)                :: y = (/0.0_rp,1.0_rp,0.0_rp/)
    REAL(rp), DIMENSION(3)                :: z = (/0.0_rp,0.0_rp,1.0_rp/)
    INTEGER                               :: jj

    ALLOCATE( Vo(spp%ppp) )
    ALLOCATE( V1(spp%ppp) )
    ALLOCATE( V2(spp%ppp) )
    ALLOCATE( V3(spp%ppp) )
    ALLOCATE( b1(3,spp%ppp) )
    ALLOCATE( b2(3,spp%ppp) )
    ALLOCATE( b3(3,spp%ppp) )

    ALLOCATE( theta(spp%ppp) )

    ! * * * * INITIALIZE VELOCITY * * * *

    call init_random_seed()
    call RANDOM_NUMBER(theta)
    theta = 2.0_rp*C_PI*theta

    Vo = SQRT( 1.0_rp - 1.0_rp/(spp%vars%g(:)**2) )
    V1 = Vo*COS(C_PI*spp%vars%eta/180.0_rp)
    V2 = Vo*SIN(C_PI*spp%vars%eta/180.0_rp)*COS(theta)
    V3 = Vo*SIN(C_PI*spp%vars%eta/180.0_rp)*SIN(theta)

    call unitVectors(params, spp%vars, F, b1, b2, b3)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj)
    do jj = 1_idef, spp%ppp
        if ( spp%vars%flag(jj) .EQ. 1_idef ) then
            spp%vars%V(1,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),x)                  &
                             + V2(jj)*DOT_PRODUCT(b2(:,jj),x)                  &
                             + V3(jj)*DOT_PRODUCT(b3(:,jj),x)

            spp%vars%V(2,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),y)                  &
                             + V2(jj)*DOT_PRODUCT(b2(:,jj),y)                  &
                             + V3(jj)*DOT_PRODUCT(b3(:,jj),y)

            spp%vars%V(3,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),z)                  &
                             + V2(jj)*DOT_PRODUCT(b2(:,jj),z)                  &
                             + V3(jj)*DOT_PRODUCT(b3(:,jj),z)
        end if
    end do
!$OMP END PARALLEL DO

    DEALLOCATE(theta)
    DEALLOCATE(Vo)
    DEALLOCATE(V1)
    DEALLOCATE(V2)
    DEALLOCATE(V3)
    DEALLOCATE(b1)
    DEALLOCATE(b2)
    DEALLOCATE(b3)
end subroutine gyro_distribution


!> @brief Subroutine that works as an interface for initializing various gyro-angle distributions for the different simulated particle species.
!!
!! @note At this moment this subroutine only calls the subroutine to generate a uniform gyro-angle distribution. This will be modified later.
!! @param[in] params Core KORC simulation parameters.
!! @param[in] F An instance of the KORC derived type FIELDS. This structure has the information of the magnetic field.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param ss Species iterator.
subroutine initial_gyro_distribution(params,F,spp)
    TYPE(KORC_PARAMS), INTENT(IN)              :: params
    TYPE(FIELDS), INTENT(IN)                   :: F
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp
    INTEGER                                    :: ss

    do ss=1_idef,params%num_species
        SELECT CASE (TRIM(spp(ss)%energy_distribution))
            CASE ('THERMAL')
                !Nothing, all was done in initialize_particles through thermal_distribution
            CASE DEFAULT
                call gyro_distribution(params,F,spp(ss))
        END SELECT
    end do
end subroutine initial_gyro_distribution

END MODULE korc_velocity_distribution
