!> @brief Module with subroutines for generating the initial spatial distribution of the different partciles' species in the simulation.
MODULE korc_spatial_distribution
    USE korc_types
    USE korc_constants
    USE korc_HDF5
    USE korc_hpc
    use korc_fields
    use random
    use korc_hammersley_generator

    IMPLICIT NONE

    PUBLIC :: intitial_spatial_distribution
    PRIVATE :: uniform,                    &
               disk,                       &
               torus,                      &
               elliptic_torus,             &
               exponential_elliptic_torus, &
               gaussian_elliptic_torus,    &
               exponential_torus,          &
               gaussian_torus,             &
               fzero

    CONTAINS


!> @brief Initializing to zero the particles' position when simulating a 'UNIFORM' plasma.
!! @details Even though in a simulation of a uniform plasma the particles' position is not advanced, we initialize their position to zero.
!!
!! @todo Modify KORC for not allocating the particles' position spp%vars%X and to do not use it along the simulation.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
subroutine uniform(spp)
    TYPE(SPECIES), INTENT(INOUT) :: spp

    spp%vars%X = 0.0_rp
end subroutine uniform


!> @brief Subrotuine for generating a uniform disk/ring as the initial spatial condition of a given species of particles in the simulation.
!! @details This uniform disk/ring distribution is generated using the Inverse Transform Sampling method. In this case, the (toroidal) radial distribution function of the particles is:
!!
!!
!! @f$f(r) = \left\{ \begin{array}{ll} 0 & r<r_{min} \\ \frac{1}{2\pi^2(r_{max}^2-r_{min}^2)R_0} & r_{min}<r<r_{max} \\ 0 & r>r_{max} \end{array} \right. @f$,
!!
!!
!! where @f$r_{min}@f$ and @f$r_{max}@f$ are the inner and outer radius of the uniform ring distribution, and @f$R_0@f$ is the cylindrical radial position of the center of the disk/ring distribution.
!! This distribution is so that @f$\int_0^{2\pi}\int_{r_{min}}^{r_{max}} f(r) J(r,\theta) drd\theta = 1 @f$, where @f$\theta@f$ is the poloidal angle,
!! and @f$J(r,\theta)=r(R_0 + r\cos\theta)@f$ is the Jacobian of the transformation of Cartesian coordinates to toroidal coordinates.
!! Notice that in the case of a disk @f$r_{min}=0@f$. As a convention, this spatial distribution will be generated on the @f$xz@f$-plane.
!! Using the Inverse Transform Sampling method we sample @f$f(r)@f$, and obtain the radial position of the particles as @f$r = \sqrt{(r_{max}^2 - r_{min}^2)U + r_{min}^2}@f$,
!! where @f$U@f$ is a uniform deviate in @f$[0,1]@f$.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
subroutine disk(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta

    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)

!$OMP PARALLEL WORKSHARE
    spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*COS(spp%PHIo)
    spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*SIN(spp%PHIo)
    spp%vars%X(3,:) = spp%Zo + r*SIN(theta)
!$OMP END PARALLEL WORKSHARE

    DEALLOCATE(theta)
    DEALLOCATE(r)
end subroutine disk


!> @brief Subrotuine for generating a uniform torus/torus shell as the initial spatial condition of a given species of particles in the simulation.
!! @details This distribution is generated using the Inverse Transform Sampling method. This distribution follows the same radial distribution of a
!! uniform disk/ring distribution, see the documentation of the \ref korc_spatial_distribution.disk subroutine.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
!! @param zeta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform toroidal angle @f$\zeta@f$ distribution of the particles.
subroutine torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta

    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)
!
    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)

!$OMP PARALLEL WORKSHARE
    spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
    spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
    spp%vars%X(3,:) = spp%Zo + r*SIN(theta)
!$OMP END PARALLEL WORKSHARE

    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine torus


!> @brief Subroutine for generating a uniform elliptic torus as the initial spatial condition of a given particle species in the simulation.
!! @details An initial spatial distribution following the uniform distribution of \ref korc_spatial_distribution.torus is modified through
!! a shear transformation and a rotation to generate a uniform spatial distribution on tori with elliptic cross sections.
!! First, we obtain the uniform spatial distribution in a torus of minor radius @f$r_0@f$, see \ref korc_spatial_distribution.torus.
!! Then, we perform a shear transformation that changes the cross section of the torus from circular to a tilted ellipse.
!! In cylindrical coordinates this shear transformation is given by:
!!
!!
!! @f$R' = R + \alpha Z@f$,\n
!! @f$Z' = Z@f$,
!!
!!
!! where @f$\alpha@f$ is the shear factor of the transformation. Here, @f$R@f$ and @f$Z@f$ are the radial and vertical position of the particles
!! uniformly distributed in a circular torus, @f$R'@f$ and @f$Z'@f$ are their new positions when following a uniform distribution in a torus with
!! elliptic circular cross section. The center of the ellipse is @f$R_0' = R_0 + \alpha Z_0@f$, and @f$Z_0 = Z_0@f$, where @f$R_0@f$ and @f$Z_0@f$
!! is the center of the initial circular torus. The major and minor semi-axes of the tilted ellipse cross section is:
!!
!!
!! @f$a' = \left[ - \frac{2r_0^2}{\alpha \sqrt{\alpha^2 + 4} - (2+\alpha^2)}  \right]^{1/2}@f$\n
!! @f$b' = \left[ \frac{2r_0^2}{\alpha \sqrt{\alpha^2 + 4} + (2+\alpha^2)}  \right]^{1/2}@f$.
!!
!! Finally, we rotate the ellipse cross section anticlockwise along @f$(R_0',Z_0')@f$ by @f$\Theta = \cot^{-1}(\alpha/2)/2@f$, so the major semi-axis is
!! parallel to the @f$Z@f$-axis.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param rotation_angle This is the angle @f$\Theta@f$ in the equations above.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
!! @param zeta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform toroidal angle @f$\zeta@f$ distribution of the particles.
!! @param X Auxiliary vector used in the coordinate transformations.
!! @param Y Auxiliary vector used in the coordinate transformations.
!! @param X1 Auxiliary vector used in the coordinate transformations.
!! @param Y1 Auxiliary vector used in the coordinate transformations.
subroutine elliptic_torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: rotation_angle
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Y
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X1
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Y1

    ALLOCATE(X1(spp%ppp))
    ALLOCATE(Y1(spp%ppp))
    ALLOCATE(X(spp%ppp))
    ALLOCATE(Y(spp%ppp))
    ALLOCATE( rotation_angle(spp%ppp) )
    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)

    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)

    Y = r*SIN(theta)
    X = r*COS(theta) + spp%shear_factor*Y

    rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor); !> @todo Modify this approximation.

    X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
     Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

    spp%vars%X(1,:) = X1*SIN(zeta)
    spp%vars%X(2,:) = X1*COS(zeta)
    spp%vars%X(3,:) = Y1

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    DEALLOCATE(X)
    DEALLOCATE(Y)
    DEALLOCATE(rotation_angle)
    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine elliptic_torus

!> @brief Function used to find the zeros of @f$f(r)@f$ of \ref korc_spatial_distribution.exponential_torus.
!!
!! @param f Value of function.
!! @param r Guess value of radial position of the particles.
!! @param a Minor radius of the toroidal distribution @f$r_0@f$
!! @param ko Decay rate of radial distribution, see @f$f(r)@f$ of \ref korc_spatial_distribution.exponential_torus.
!! @param P Deviate of a random uniform distribution in the interval @f$[0,1]@f$.
FUNCTION fzero(r,a,ko,P) RESULT(f)
    REAL(rp)             :: f
    REAL(rp), INTENT(IN) :: r
    REAL(rp), INTENT(IN) :: a
    REAL(rp), INTENT(IN) :: ko
    REAL(rp), INTENT(IN) :: P

    f = EXP(-ko*r)*(1.0_rp + r*ko) + ( 1.0_rp - EXP(-ko*a)*(1.0_rp + a*ko) )*P - 1.0_rp
END FUNCTION fzero


!> @brief Subroutine that generates a exponentially decaying radial distribution of particles in a circular cross-section torus of
!! major and minor radi @f$R_0@f$ and @f$r_0@f$, respectively.
!! @details We generate this exponentially decaying radial distribution @f$f(r)@f$ following the same approach as in
!! \ref korc_spatial_distribution.disk, but this time, the radial distribution is given by:
!!
!!
!! @f$f(r) = \left\{ \begin{array}{ll} \frac{k_0^2}{4\pi^2 R_0}\frac{\exp{\left( -k_0 r \right)}}{1 - \exp{\left( -k_0r_0\right)}\left[ 1 + k_0 r_0 \right]} &\
!! r<r_0 \\ 0 & r>r_0 \end{array} \right. @f$.
!!
!!
!! The radial position of the particles @f$r@f$ is obtained using the Inverse Trasnform Sampling method, finding @f$r@f$ numerically
!! through the Newton-Raphson method. First, we calculate the particles' radial distribution in a disk centered at @f$(R,Z) = (0,0)@f$.
!! Then, we transfor to a new set of coordinates where the disk is centered at @f$(R,Z) = (R_0,Z_0)@f$. Finally, we generate the
!! toroidal distribution by givin each particle a toroidal angle @f$\zeta@f$ which follows a uniform distribution in the interval
!! @f$[0,2\pi]@f$.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param fl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param relerr Tolerance used to determine when to stop iterating the Newton-Raphson method for finding @f$r@f$.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
!! @param zeta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform toroidal angle @f$\zeta@f$ distribution of the particles.
!! @param pp Particle iterator.
subroutine exponential_torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp)                            :: fl
    REAL(rp)                            :: fr
    REAL(rp)                            :: fm
    REAL(rp)                            :: rl
    REAL(rp)                            :: rr
    REAL(rp)                            :: rm
    REAL(rp)                            :: relerr
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
    INTEGER                             :: pp

    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)

    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    do pp=1_idef,spp%ppp ! Newton-Raphson applied here for finding the radial distribution
        rl = 0.0_rp
        rr = spp%r_outter

        fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
        fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))
        if (fl.GT.korc_zero) then
            relerr = 100*ABS(fl - fr)/fl
        else
            relerr = 100*ABS(fl - fr)/fr
        end if

        do while(relerr.GT.1.0_rp)
            rm = 0.5_rp*(rr - rl) + rl
            fm = fzero(rm,spp%r_outter,spp%falloff_rate,r(pp))

            if (SIGN(1.0_rp,fm).EQ.SIGN(1.0_rp,fr)) then
                rr = rm
            else
                rl = rm
            end if

            fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
            fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))

            if (fl.GT.korc_zero) then
                relerr = 100*ABS(fl - fr)/fl
            else
                relerr = 100*ABS(fl - fr)/fr
            end if
        end do
        r(pp) = rm
    end do

    spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
    spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
    spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine exponential_torus


!> @brief Subroutine that generates an exponentially decaying radial distribution in an elliptic torus as the initial spatial
!! condition of a given particle species in the simulation.
!! @details As a first step, we generate an exponentially decaying radial distribution in a circular cross-section torus as in
!! \ref korc_spatial_distribution.exponential_torus. Then we transform this spatial distribution to a one in an torus with an
!! elliptic cross section, this following the same approach as in \ref korc_spatial_distribution.elliptic_torus.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param fl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param relerr Tolerance used to determine when to stop iterating the Newton-Raphson method for finding @f$r@f$.
!! @param rotation_angle This is the angle @f$\Theta@f$ in \ref korc_spatial_distribution.elliptic_torus.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
!! @param zeta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform toroidal angle @f$\zeta@f$ distribution of the particles.
!! @param X Auxiliary vector used in the coordinate transformations.
!! @param Y Auxiliary vector used in the coordinate transformations.
!! @param X1 Auxiliary vector used in the coordinate transformations.
!! @param Y1 Auxiliary vector used in the coordinate transformations.
!! @param pp Particle iterator.
subroutine exponential_elliptic_torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp)                            :: fl
    REAL(rp)                            :: fr
    REAL(rp)                            :: fm
    REAL(rp)                            :: rl
    REAL(rp)                            :: rr
    REAL(rp)                            :: rm
    REAL(rp)                            :: relerr
    REAL(rp), DIMENSION(:), ALLOCATABLE :: rotation_angle
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Y
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X1
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Y1
    INTEGER                             :: pp

    ALLOCATE(X1(spp%ppp))
    ALLOCATE(Y1(spp%ppp))
    ALLOCATE(X(spp%ppp))
    ALLOCATE(Y(spp%ppp))
    ALLOCATE( rotation_angle(spp%ppp) )
    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)

    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    do pp=1_idef,spp%ppp
        rl = 0.0_rp
        rr = spp%r_outter

        fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
        fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))
        if (fl.GT.korc_zero) then
            relerr = 100*ABS(fl - fr)/fl
        else
            relerr = 100*ABS(fl - fr)/fr
        end if

        do while(relerr.GT.1.0_rp)
            rm = 0.5_rp*(rr - rl) + rl
            fm = fzero(rm,spp%r_outter,spp%falloff_rate,r(pp))

            if (SIGN(1.0_rp,fm).EQ.SIGN(1.0_rp,fr)) then
                rr = rm
            else
                rl = rm
            end if

            fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
            fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))

            if (fl.GT.korc_zero) then
                relerr = 100*ABS(fl - fr)/fl
            else
                relerr = 100*ABS(fl - fr)/fr
            end if
        end do
        r(pp) = rm
    end do

    Y = r*SIN(theta)
    X = r*COS(theta) + spp%shear_factor*Y

    rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

    X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
     Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

    spp%vars%X(1,:) = X1*SIN(zeta)
    spp%vars%X(2,:) = X1*COS(zeta)
    spp%vars%X(3,:) = Y1

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    DEALLOCATE(X)
    DEALLOCATE(Y)
    DEALLOCATE(rotation_angle)
    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine exponential_elliptic_torus


!> @brief Subroutine that generates a Gaussian radial distribution in an elliptic torus as the initial spatial
!! condition of a given particle species in the simulation.
!! @details As a first step, we generate an Gaussian radial distribution in a circular cross-section torus as in
!! \ref korc_spatial_distribution.gaussian_torus. Then we transform this spatial distribution to a one in an torus with an
!! elliptic cross section, this following the same approach as in \ref korc_spatial_distribution.elliptic_torus.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param rotation_angle This is the angle @f$\Theta@f$ in \ref korc_spatial_distribution.elliptic_torus.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
!! @param zeta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform toroidal angle @f$\zeta@f$ distribution of the particles.
!! @param X Auxiliary vector used in the coordinate transformations.
!! @param Y Auxiliary vector used in the coordinate transformations.
!! @param X1 Auxiliary vector used in the coordinate transformations.
!! @param Y1 Auxiliary vector used in the coordinate transformations.
!! @param sigma Standard deviation @f$\sigma@f$ of the radial distribution function.
!! @param pp Particle iterator.
subroutine gaussian_elliptic_torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: rotation_angle
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Y
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X1
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Y1
    REAL(rp)                            :: sigma
    INTEGER                             :: pp

    ALLOCATE(X1(spp%ppp))
    ALLOCATE(Y1(spp%ppp))
    ALLOCATE(X(spp%ppp))
    ALLOCATE(Y(spp%ppp))
    ALLOCATE( rotation_angle(spp%ppp) )
    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)

    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    sigma = 1.0_rp/SQRT(2.0_rp*(spp%falloff_rate/params%cpp%length))
    sigma = sigma/params%cpp%length

    r = sigma*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/sigma**2))*r))
    spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
    spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
    spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

    Y = r*SIN(theta)
    X = r*COS(theta) + spp%shear_factor*Y

    rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

    X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
     Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

    spp%vars%X(1,:) = X1*SIN(zeta)
    spp%vars%X(2,:) = X1*COS(zeta)
    spp%vars%X(3,:) = Y1

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    DEALLOCATE(X)
    DEALLOCATE(Y)
    DEALLOCATE(rotation_angle)
    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine gaussian_elliptic_torus


!> @brief Subroutine that generates a Gaussian radial distribution of particles in a circular cross-section torus of
!! major and minor radi @f$R_0@f$ and @f$r_0@f$, respectively.
!! @details We generate this exponentially decaying radial distribution @f$f(r)@f$ following the same approach as in
!! \ref korc_spatial_distribution.disk, but this time, the radial distribution is given by:
!!
!!
!! @f$f(r) = \left\{ \begin{array}{ll} \frac{1}{4\pi^2 \sigma^2 R_0}\frac{\exp{\left( -\frac{r^2}{2\sigma^2} \right)}}{1 - \exp{\left( -\frac{r_0^2}{2\sigma^2} \right)}} &\
!! r<r_0 \\ 0 & r>r_0 \end{array} \right. @f$.
!!
!!
!! The radial position of the particles @f$r@f$ is obtained using the Inverse Trasnform Sampling method, finding @f$r@f$ numerically
!! through the Newton-Raphson method. First, we calculate the particles' radial distribution in a disk centered at @f$(R,Z) = (0,0)@f$.
!! Then, we transfor to a new set of coordinates where the disk is centered at @f$(R,Z) = (R_0,Z_0)@f$. Finally, we generate the
!! toroidal distribution by givin each particle a toroidal angle @f$\zeta@f$ which follows a uniform distribution in the interval
!! @f$[0,2\pi]@f$.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param fl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param relerr Tolerance used to determine when to stop iterating the Newton-Raphson method for finding @f$r@f$.
!! @param r Radial position of the particles @f$r@f$.
!! @param theta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform poloidal angle @f$\theta@f$ distribution of the particles.
!! @param zeta Uniform deviates in the range @f$[0,2\pi]@f$ representing the uniform toroidal angle @f$\zeta@f$ distribution of the particles.
!! @param pp Particle iterator.
subroutine gaussian_torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r ! temporary vars
    REAL(rp)                            :: sigma

    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

!    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta

    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)

    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    sigma = 1.0_rp/SQRT(2.0_rp*(spp%falloff_rate/params%cpp%length))
    sigma = sigma/params%cpp%length

    r = sigma*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - EXP(-0.5_rp*spp%r_outter**2/sigma**2))*r))
    spp%vars%X(1,:) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
    spp%vars%X(2,:) = ( spp%Ro + r*COS(theta) )*COS(zeta)
    spp%vars%X(3,:) = spp%Zo + r*SIN(theta)

    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine gaussian_torus


!> @brif Subroutine that contains calls to the different subroutines for initializing the simulated particles with various
!! spatial distribution functions.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param ss Species iterator.
subroutine intitial_spatial_distribution(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)              :: params
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp
    INTEGER                                    :: ss

    do ss=1_idef,params%num_species
        SELECT CASE (TRIM(spp(ss)%spatial_distribution))
            CASE ('UNIFORM')
                call uniform(spp(ss))
            CASE ('DISK')
                call disk(params,spp(ss))
            CASE ('TORUS')
                call torus(params,spp(ss))
            CASE ('EXPONENTIAL-TORUS')
                call exponential_torus(params,spp(ss))
            CASE ('GAUSSIAN-TORUS')
                call gaussian_torus(params,spp(ss))
            CASE ('ELLIPTIC-TORUS')
                call elliptic_torus(params,spp(ss))
            CASE ('EXPONENTIAL-ELLIPTIC-TORUS')
                call exponential_elliptic_torus(params,spp(ss))
            CASE ('GAUSSIAN-ELLIPTIC-TORUS')
                call gaussian_elliptic_torus(params,spp(ss))
            CASE DEFAULT
                call torus(params,spp(ss))
        END SELECT
    end do
end subroutine intitial_spatial_distribution


END MODULE korc_spatial_distribution
