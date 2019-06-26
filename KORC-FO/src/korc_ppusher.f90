!> @brief Module with subroutines for advancing the particles' position and velocity in the simulations.
module korc_ppusher
    use korc_types
    use korc_constants
    use korc_fields
    use korc_profiles
    use korc_interp
    use korc_collisions
    use korc_hpc

    IMPLICIT NONE

    REAL(rp), PRIVATE :: E0 !< Dimensionless vacuum permittivity @f$\epsilon_0 \times (m_{ch}^2 v_{ch}^3/q_{ch}^3 B_{ch})@f$, see korc_units.f90.

    PRIVATE :: cross,&
               radiation_force,&
               collision_force
    PUBLIC :: initialize_particle_pusher,&
                advance_particles_position,&
                advance_particles_velocity

    contains


!> @brief This subroutine initializes all the variables needed for advancing the particles' position and velocity.
!! @details This subroutine is specially useful when we need to define or initialize values of parameters used to calculate derived quantities.
!! The intent of this subroutine is to work as a constructor of the module.
!!
!! @param[in] params Core KORC simulation parameters.
subroutine initialize_particle_pusher(params)
    TYPE(KORC_PARAMS), INTENT(IN)  :: params

    E0 = C_E0*(params%cpp%mass**2*params%cpp%velocity**3)/(params%cpp%charge**3*params%cpp%Bo)
end subroutine initialize_particle_pusher


!> @brief Function that calculates and returns the cross product @f$\mathbf{a}\times \mathbf{b}@f$. These vectors are in Cartesian coordinates.
!! @note Notice that all the variables in this subroutine have been normalized using the characteristic scales in korc_units.f90.
!!
!! @param[in] a Vector @f$\mathbf{a}@f$.
!! @param[in] b Vector @f$\mathbf{b}@f$.
!! @param cross Value of @f$\mathbf{a}\times \mathbf{b}@f$
function cross(a,b)
    REAL(rp), DIMENSION(3), INTENT(IN) :: a
    REAL(rp), DIMENSION(3), INTENT(IN) :: b
    REAL(rp), DIMENSION(3)             :: cross

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross


!> @brief Subroutine that calculates the synchrotron radiation reaction force.
!! @details This subroutine calculates the synchrotron radiation reaction force [Carbajal et al. PoP <b>24</b>, 042512 (2017)] using the derivation of Landau-Lifshiftz of the
!! Lorentz-Abraham-Dirac radiation reaction force:\n
!! @f$\mathbf{F}_R(\mathbf{x},\mathbf{v}) = \frac{q^3}{6\pi\epsilon_0 m c^3}\left[ \mathbf{F}_1 + \mathbf{F}_2 + \mathbf{F}_3\right]@f$,
!!
!!
!! @f$\mathbf{F}_1 = \gamma \left( \frac{D \mathbf{E}}{Dt} + \mathbf{v}\times \frac{D \mathbf{B}}{Dt} \right)@f$,\n
!! @f$\mathbf{F}_2 = \frac{q}{m}\left( \frac{(\mathbf{E}\cdot\mathbf{v})}{c^2}\mathbf{E} + (\mathbf{E} + \mathbf{v}\times \mathbf{B})\times \mathbf{B} \right)@f$,\n
!! @f$\mathbf{F}_3 = -\frac{q\gamma^2}{mc^2} \left( (\mathbf{E} + \mathbf{v}\times \mathbf{B})^2 -  \frac{(\mathbf{E}\cdot\mathbf{v})^2}{c^2} \right)\mathbf{v}@f$,
!!
!!
!! where @f$\gamma = 1/\sqrt{1 - v^2/c^2}@f$ is the relativistic factor, @f$D/Dt = \partial/\partial t + \mathbf{v}\cdot\nabla@f$, @f$q@f$ and @f$m@f$ are the charge and mass of the particle,
!! and @f$\epsilon_0@f$ is the vacuum permittivity. For relativistic electrons we have @f$F_1 \ll F_2@f$ and @f$F_1 \ll F_3@f$, therefore @f$\mathbf{F}_1@f$ is not calculated here.
!!
!! @note Notice that all the variables in this subroutine have been normalized using the characteristic scales in korc_units.f90.
!! @param[in] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param[in] U @f$\mathbf{u} = \gamma \mathbf{v}@f$, where @f$\mathbf{v}@f$ is the particle's velocity.
!! @param[in] E Electric field @f$\mathbf{E}@f$ seen by each particle. This is given in Cartesian coordinates.
!! @param[in] B Magnetic field @f$\mathbf{B}@f$ seen by each particle. This is given in Cartesian coordinates.
!! @param[out] Frad The calculated synchrotron radiation reaction force @f$\mathbf{F}_R@f$.
!! @param F1 The component @f$\mathbf{F}_1@f$ of @f$\mathbf{F}_R@f$.
!! @param F2 The component @f$\mathbf{F}_2@f$ of @f$\mathbf{F}_R@f$.
!! @param F3 The component @f$\mathbf{F}_3@f$ of @f$\mathbf{F}_R@f$.
!! @param V The particle's velocity @f$\mathbf{v}@f$.
!! @param vec An auxiliary 3-D vector.
!! @param g The relativistic @f$\gamma@f$ factor of the particle.
subroutine radiation_force(spp,U,E,B,Frad)
    TYPE(SPECIES), INTENT(IN)           :: spp
    REAL(rp), DIMENSION(3), INTENT(IN)  :: U
    REAL(rp), DIMENSION(3), INTENT(IN)  :: E
    REAL(rp), DIMENSION(3), INTENT(IN)  :: B
    REAL(rp), DIMENSION(3), INTENT(OUT) :: Frad
    REAL(rp), DIMENSION(3)              :: F1
    REAL(rp), DIMENSION(3)              :: F2
    REAL(rp), DIMENSION(3)              :: F3
    REAL(rp), DIMENSION(3)              :: V
    REAL(rp), DIMENSION(3)              :: vec
    REAL(rp)                            :: g
    REAL(rp)                            :: tmp

    g = SQRT(1.0_rp + DOT_PRODUCT(U,U))
    V = U/g

    tmp = spp%q**4/(6.0_rp*C_PI*E0*spp%m**2)

    F2 = tmp*( DOT_PRODUCT(E,V)*E + cross(E,B) + cross(B,cross(B,V)) )
    vec = E + cross(V,B)
    F3 = (tmp*g**2)*( DOT_PRODUCT(E,V)**2 - DOT_PRODUCT(vec,vec) )*V

    Frad = F2 + F3
end subroutine radiation_force


!> @brief Subroutine for advancing the particles' velocity.
!! @details We are using the modified relativistic leapfrog method of J.-L. Vay, PoP <b>15</b>, 056701 (2008) for advancing the particles'
!! position and velocity. For including the synchrotron radiation reaction force we used the scheme in Tamburini et al., New J. Phys. <b>12</b>, 123005 (2010).
!! A comprehensive description of this can be found in Carbajal et al., PoP <b>24</b>, 042512 (2017).
!! The discretized equations of motion to advance the change in the position and momentum due to the Lorentz force are:
!!
!!
!! @f$\frac{\mathbf{x}^{i+1/2} - \mathbf{x}^{i-1/2}}{\Delta t}  = \mathbf{v}^i@f$\n
!! @f$\frac{\mathbf{p}^{i+1}_L - \mathbf{p}^{i}}{\Delta t} = q \left(  \mathbf{E}^{i+1/2} + \frac{\mathbf{v}^i + \mathbf{v}^{i+1}_L}{2} \times \mathbf{B}^{i+1/2} \right)@f$
!!
!!
!! where @f$\Delta t@f$ is the time step, @f$q@f$ denotes the charge,  @f$\mathbf{p}^j = m \gamma^j \mathbf{v}^j@f$, and @f$\gamma^j = 1/\sqrt{1 + v^{j2}/c^2}@f$.
!! Here @f$i@f$ and @f$i+1@f$ indicate integer time leves, while @f$i-1/2@f$ and @f$i+1/2@f$ indicate half-time steps.
!! The evolution of the relativistic @f$\gamma@f$ factor is given by @f$\gamma^{i+1} = \sqrt{1 + \left(p_L^{i+1}/mc \right)^2} = \sqrt{1 + \mathbf{p}_L^{i+1}\cdot \mathbf{p}'/m^2c^2}@f$, which can be combined with the above equations to produce:
!!
!!
!! @f$\mathbf{p}^{i+1}_L = s\left[ \mathbf{p}' + (\mathbf{p}'\cdot\mathbf{t})\mathbf{t} + \mathbf{p}'\times \mathbf{t} \right]@f$\n
!! @f$\gamma^{i+1} = \sqrt{\frac{\sigma + \sqrt{\sigma^2 + 4(\tau^2 + p^{*2})}}{2}}@f$
!!
!!
!! where we have defined @f$\mathbf{p}' = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} + \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)@f$,
!! @f$\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}@f$, @f$\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}@f$, @f$p^{*} = \mathbf{p}'\cdot \mathbf{\tau}/mc@f$, @f$\sigma = \gamma'^2 - \tau^2@f$, @f$\gamma' = \sqrt{1 + p'^2/m^2c^2}@f$, and @f$s = 1/(1+t^2)@f$.
!! The discretized equation of motion to advance the change in the momentum due to the radiation reaction force force is
!!
!!
!! @f$\frac{\mathbf{p}^{i+1}_R - \mathbf{p}^{i}}{\Delta t} = \mathbf{F}_R(\mathbf{x}^{i+1/2},\mathbf{p}^{i+1/2})@f$,
!!
!!
!! where @f$\mathbf{p}^{i+1/2} = (\mathbf{p}^{i+1}_L + \mathbf{p}^i)/2@f$. Finally, using  @f$\mathbf{p}^{i+1}_L@f$ and @f$\mathbf{p}^{i+1}_R@f$, the momentum at time level @f$i+1@f$ is given by\n
!! @f$\mathbf{p}^{i+1}  = \mathbf{p}^{i+1}_L + \mathbf{p}^{i+1}_R - \mathbf{p}^i.@f$mu
!! Collisions are included by solving the stochastic differential equation in a Cartesian coordinate system where @f$\mathbf{p}@f$ is parallel to @f$\hat{e}_z@f$:mu
!! @f$\mathbf{p} = \mathbf{A}dt + \hat{\sigma}\cdot d\mathbf{W}@f$,
!!
!!
!! where @f$\mathbf{A} = p \nu_s\hat{b}@f$, @f$\hat{b}=\mathbf{B}/B@f$ with @f$\mathbf{B}@f$ the magnetic field, and @f$\nu_s@f$ the collision frequency that corresponds to the drag force due to collisions.
!! @f$\hat{\sigma}@f$ is a diagonal 3x3 matrix with elements @f$\hat{\sigma}_{11} = p\sqrt{\nu_{\parallel}}@f$, and @f$\hat{\sigma}_{22} = \hat{\sigma}_{33} = p\sqrt{\nu_{D}}@f$, with @f$\nu_\parallel@f$ and @f$\nu_D@f$
!! the collisional frequencies producing diffusive transport along and across the direction of @f$\mathbf{p}@f$, respectively.
!! @note Notice that all the variables in this subroutine have been normalized using the characteristic scales in korc_units.f90.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in] P An instance of the KORC derived type PROFILES.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param[in] bool Logical variable used to indicate if we calculate or not quantities listed in the outputs list.
!! @param[in] dt Time step used in the leapfrog step (@f$\Delta t@f$).
!! @param Prad Total radiated power of each particle.
!! @param B Magnitude of the magnetic field seen by each particle .
!! @param v Speed of each particle.
!! @param vpar Parallel velocity @f$v_\parallel = \mathbf{v}\cdot \hat{b}@f$.
!! @param vperp Perpendicular velocity @f$v_\parallel = |\mathbf{v} - (\mathbf{v}\cdot \hat{b})\hat{b}|@f$.
!! @param tmp Temporary variable used for various computations.
!! @param a This variable is used to simplify notation in the code, and is given by @f$a=q\Delta t/m@f$,
!! @param gp This variable is @f$\gamma' = \sqrt{1 + p'^2/m^2c^2}@f$ in the above equations.
!! @param sigma This variable is @f$\sigma = \gamma'^2 - \tau^2@f$ in the above equations.
!! @param us This variable is @f$u^{*} = p^{*}/m@f$ where @f$ p^{*} = \mathbf{p}'\cdot \mathbf{\tau}/mc@f$.
!! @param g Relativistic factor @f$\gamma@f$.
!! @param s This variable is @f$s = 1/(1+t^2)@f$ in the equations above.
!! @param U_L This variable is @f$\mathbf{u}_L = \mathbf{p}_L/m@f$ where @f$\mathbf{p}^{i+1}_L = s\left[ \mathbf{p}' + (\mathbf{p}'\cdot\mathbf{t})\mathbf{t} + \mathbf{p}'\times \mathbf{t} \right]@f$.
!! @param U_hs Is @f$\mathbf{u}=\mathbf{p}/m@f$ at half-time step (@f$i+1/2@f$) in the absence of radiation losses or collisions. @f$\mathbf{u}^{i+1/2} = \mathbf{u}^i + \frac{q\Delta t}{2m}\left( \mathbf{E}^{i+1/2} + \mathbf{v}^i\times \mathbf{B}^{i+1/2} \right)@f$.
!! @param tau This variable is @f$\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}@f$.
!! @param up This variable is @f$\mathbf{u}'= \mathbf{p}'/m@f$, where @f$\mathbf{p}' = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} + \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)@f$.
!! @param t This variable is @f$\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}@f$.
!! @param U This variable is @f$\mathbf{u}^{i+1}= \mathbf{p}^{i+1}/m@f$.
!! @param U_RC This variable is @f$\mathbf{u}^{i+1}_R= \mathbf{p}^{i+1}_R/m@f$
!! @param U_os This variable is @f$\mathbf{u}^{i+1/2}= \mathbf{p}^{i+1/2}/m@f$ when radiation losses are included. Here, @f$\mathbf{p}^{i+1/2} = (\mathbf{p}^{i+1}_L + \mathbf{p}^i)/2@f$.
!! @param Frad Synchrotron radiation reaction force of each particle.
!! @param vec Auxiliary vector used in various computations.
!! @param b_unit Unitary vector pointing along the local magnetic field @f$\hat{b}@f$.
!! @param ii Species iterator.
!! @param pp Particles iterator.
!! @param ss_collisions Logical variable that indicates if collisions are included in the simulation.
subroutine advance_particles_velocity(params,F,P,spp,dt,bool)
    TYPE(KORC_PARAMS), INTENT(IN)              :: params
    TYPE(FIELDS), INTENT(IN)                   :: F
    TYPE(PROFILES), INTENT(IN)                 :: P
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN)                        :: bool
    REAL(rp), INTENT(IN)                       :: dt
    REAL(rp)                                   :: Prad
    REAL(rp)                                   :: B
    REAL(rp)                                   :: v
    REAL(rp)                                   :: vpar
    REAL(rp)                                   :: vperp
    REAL(rp)                                   :: tmp
    REAL(rp)                                   :: a
    REAL(rp)                                   :: gp
    REAL(rp)                                   :: sigma
    REAL(rp)                                   :: us
    REAL(rp)                                   :: g
    REAL(rp)                                   :: s
    REAL(rp), DIMENSION(3)                     :: U_L
    REAL(rp), DIMENSION(3)                     :: U_hs
    REAL(rp), DIMENSION(3)                     :: tau
    REAL(rp), DIMENSION(3)                     :: up
    REAL(rp), DIMENSION(3)                     :: t
    REAL(rp), DIMENSION(3)                     :: U
    REAL(rp), DIMENSION(3)                     :: U_RC
    REAL(rp), DIMENSION(3)                     :: U_os
    REAL(rp), DIMENSION(3)                     :: Frad
    REAL(rp), DIMENSION(3)                     :: vec
    REAL(rp), DIMENSION(3)                     :: b_unit
    INTEGER                                    :: ii
    INTEGER                                    :: pp
    LOGICAL                                    :: ss_collisions


    ! Determine whether we are using a single-species collision model
    ss_collisions = TRIM(params%collisions_model) .EQ. 'SINGLE_SPECIES'

    do ii = 1_idef,params%num_species

        call get_fields(params,spp(ii)%vars,F)

        call get_profiles(params,spp(ii)%vars,P)

        a = spp(ii)%q*dt/spp(ii)%m

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(pp,U,B,U_L,U_RC,U_hs,tau,up,gp,     &
!$OMP&                                     sigma,us,g,t,s,U_os,Frad,b_unit,v,  &
!$OMP&                                     vpar,vperp,tmp,vec)
        do pp=1_idef,spp(ii)%ppp
            if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then
                U = spp(ii)%vars%g(pp)*spp(ii)%vars%V(:,pp)

                !! Magnitude of magnetic field
                B = SQRT( DOT_PRODUCT(spp(ii)%vars%B(:,pp),spp(ii)%vars%B(:,pp)) )

                U_L = U
                U_RC = U

                ! ! ! LEAP-FROG SCHEME FOR LORENTZ FORCE ! ! !
                U_hs = U_L + 0.5_rp*a*( spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp)) )
                tau = 0.5_rp*dt*spp(ii)%q*spp(ii)%vars%B(:,pp)/spp(ii)%m
                up = U_hs + 0.5_rp*a*spp(ii)%vars%E(:,pp)
                gp = SQRT( 1.0_rp + DOT_PRODUCT(up,up) )
                sigma = gp**2 - DOT_PRODUCT(tau,tau)
                us = DOT_PRODUCT(up,tau) ! variable 'u^*' in Vay, J.-L. PoP (2008)
                g = SQRT( 0.5_rp*(sigma + SQRT(sigma**2 + 4.0_rp*(DOT_PRODUCT(tau,tau) + us**2))) )
                t = tau/g
                s = 1.0_rp/(1.0_rp + DOT_PRODUCT(t,t)) ! variable 's' in Vay, J.-L. PoP (2008)
                U_L = s*(up + DOT_PRODUCT(up,t)*t + cross(up,t))
                ! ! ! LEAP-FROG SCHEME FOR LORENTZ FORCE ! ! !

                ! ! ! Splitting operator for including radiation
                U_os = 0.5_rp*(U_L + U)

                if (params%radiation) then
                    call radiation_force(spp(ii),U_os,spp(ii)%vars%E(:,pp),spp(ii)%vars%B(:,pp),Frad)
                    U_RC = U_RC + a*Frad/spp(ii)%q
                end if
                ! ! ! Splitting operator for including radiation

                U = U_L + U_RC - U

                ! ! ! Stochastic differential equations for including collisions
                if (params%collisions .AND. ss_collisions) then
                    call include_CoulombCollisions(params,U,spp(ii)%vars%ne(pp),spp(ii)%vars%Te(pp),spp(ii)%vars%Zeff(pp))
                end if
                ! ! ! Stochastic differential equations for including collisions

                if (params%radiation .OR. params%collisions) then
                    g = SQRT( 1.0_rp + DOT_PRODUCT(U,U) )
                end if
                spp(ii)%vars%V(:,pp) = U/g
                spp(ii)%vars%g(pp) = g

                if (g.LT.params%minimum_particle_g) then
                    spp(ii)%vars%flag(pp) = 0_is
                end if

                if (bool) then
                    !! Parallel unit vector
                    b_unit = spp(ii)%vars%B(:,pp)/B

                    v = SQRT(DOT_PRODUCT(spp(ii)%vars%V(:,pp),spp(ii)%vars%V(:,pp)))
                    if (v.GT.korc_zero) then
                        !! Parallel and perpendicular components of velocity
                        vpar = DOT_PRODUCT(spp(ii)%vars%V(:,pp), b_unit)
                        vperp =  DOT_PRODUCT(spp(ii)%vars%V(:,pp),spp(ii)%vars%V(:,pp)) - vpar**2
                        if ( vperp .GE. korc_zero ) then
                            vperp = SQRT( vperp )
                        else
                            vperp = 0.0_rp
                        end if

                        !! Pitch angle
                        spp(ii)%vars%eta(pp) = 180.0_rp*MODULO(ATAN2(vperp,vpar), 2.0_rp*C_PI)/C_PI

                        !! Magnetic moment
                        spp(ii)%vars%mu(pp) = 0.5_rp*spp(ii)%m*g**2*vperp**2/B
                        ! See Northrop's book (The adiabatic motion of charged particles)

                        !! Radiated power
                        tmp = spp(ii)%q**4/(6.0_rp*C_PI*E0*spp(ii)%m**2)
                        vec = spp(ii)%vars%E(:,pp) + cross(spp(ii)%vars%V(:,pp),spp(ii)%vars%B(:,pp))

                        spp(ii)%vars%Prad(pp) = tmp*(DOT_PRODUCT(spp(ii)%vars%E(:,pp), spp(ii)%vars%E(:,pp)) +      &
                                                     DOT_PRODUCT(cross(spp(ii)%vars%V(:,pp), spp(ii)%vars%B(:,pp)), &
                                                                 spp(ii)%vars%E(:,pp)) +                            &
                                                     spp(ii)%vars%g(pp)**2*(DOT_PRODUCT(spp(ii)%vars%E(:,pp),       &
                                                                                        spp(ii)%vars%V(:,pp))**2 -  &
                                                                            DOT_PRODUCT(vec,vec)))

                        !! Input power due to electric field
                        spp(ii)%vars%Pin(pp) = spp(ii)%q*DOT_PRODUCT(spp(ii)%vars%E(:,pp),spp(ii)%vars%V(:,pp))
                    else
                        spp(ii)%vars%eta(pp) = 0.0_rp
                        spp(ii)%vars%mu(pp) = 0.0_rp
                        spp(ii)%vars%Prad(pp) = 0.0_rp
                        spp(ii)%vars%Pin(pp) = 0.0_rp
                    end if
                end if
            end if
        end do
!$OMP END PARALLEL DO

    end do
end subroutine advance_particles_velocity


!> @brief Subrotuine to advance particles' position.
!! @details This subroutine advances the particles position using the information of the updated velocity.\n
!! @f$\frac{\mathbf{x}^{i+1/2} - \mathbf{x}^{i-1/2}}{\Delta t}  = \mathbf{v}^i@f$\n
!!
!! @note Notice that all the variables in this subroutine have been normalized using the characteristic scales in korc_units.f90.
!! @param[in] params Core KORC simulation parameters.
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different species in the simulation.
!! @param[in] dt Time step used in the leapfrog step (@f$\Delta t@f$).
!! @param ii Species iterator.
!! @param pp Particles iterator.
subroutine advance_particles_position(params,F,spp,dt)
    TYPE(KORC_PARAMS), INTENT(IN)              :: params
    TYPE(FIELDS), INTENT(IN)                   :: F
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp
    REAL(rp), INTENT(IN)                       :: dt
    INTEGER                                    :: ii
    INTEGER                                    :: pp

    if (params%plasma_model .NE. 'UNIFORM') then
        do ii=1_idef,params%num_species
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp)
            do pp=1_idef,spp(ii)%ppp
                if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then
                    spp(ii)%vars%X(:,pp) = spp(ii)%vars%X(:,pp) + dt*spp(ii)%vars%V(:,pp)
                end if
            end do
!$OMP END PARALLEL DO
!        spp(ii)%vars%X = MERGE(spp(ii)%vars%X + dt*spp(ii)%vars%V,spp(ii)%vars%X,SPREAD(spp(ii)%vars%flag,1,3).EQ.1_idef)
        end do
    end if
end subroutine advance_particles_position

end module korc_ppusher
