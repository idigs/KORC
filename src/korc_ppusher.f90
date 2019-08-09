module korc_ppusher
  !! @note Module with subroutines for advancing the particles' position and
  !! velocity in the simulations. @endnote
  use korc_types
  use korc_constants
  use korc_fields
  use korc_profiles
  use korc_interp
  use korc_collisions
  use korc_hpc

  IMPLICIT NONE

  REAL(rp), PRIVATE :: E0
  !! Dimensionless vacuum permittivity \(\epsilon_0 \times (m_{ch}^2
  !! v_{ch}^3/q_{ch}^3 B_{ch})\), see [[korc_units]].

  PRIVATE :: cross,&
       radiation_force_p,&
       GCEoM_p,&
       aux_fields
  PUBLIC :: initialize_particle_pusher,&
       advance_FOeqn_vars,&
       advance_FOinterp_vars,&
       advance_GCeqn_vars,&
       advance_GCinterp_vars,&
       GC_init,&
       FO_init,&
       adv_GCeqn_top,&
       adv_GCinterp_top

contains



  subroutine initialize_particle_pusher(params)
    !! @note This subroutine initializes all the variables needed for advancing
    !! the particles' position and velocity. @endnote
    !! This subroutine is specially useful when we need to define or initialize
    !! values of parameters used to calculate derived quantities.
    !! The intent of this subroutine is to work as a constructor of the module.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.

    E0 = C_E0*(params%cpp%mass**2*params%cpp%velocity**3)/ &
         (params%cpp%charge**3*params%cpp%Bo)
  end subroutine initialize_particle_pusher


  pure function cross(a,b)
    !! @note Function that calculates and returns the cross product
    !! \(\mathbf{a}\times \mathbf{b}\). These vectors are in Cartesian
    !! coordinates. @endnote
    !! @note Notice that all the variables in this subroutine have been
    !! normalized using the characteristic scales in [[korc_units]]. @endnote
    REAL(rp), DIMENSION(3), INTENT(IN) :: a
    !! Vector \(\mathbf{a}\).
    REAL(rp), DIMENSION(3), INTENT(IN) :: b
    !! Vector \(\mathbf{b}\).
    REAL(rp), DIMENSION(3)             :: cross
    !!Value of \(\mathbf{a}\times \mathbf{b}\)

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function cross



  subroutine radiation_force_p(q_cache,m_cache,U_X,U_Y,U_Z,E_X,E_Y,E_Z, &
       B_X,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)

    REAL(rp), INTENT(IN)                       :: m_cache,q_cache
    
    REAL(rp), DIMENSION(8), INTENT(IN)     :: U_X,U_Y,U_Z
    !! \(\mathbf{u} = \gamma \mathbf{v}\), where \(\mathbf{v}\) is the
    !! particle's velocity.
    REAL(rp), DIMENSION(8), INTENT(IN)     :: E_X,E_Y,E_Z
    !! Electric field \(\mathbf{E}\) seen by each particle. This is given
    !! in Cartesian coordinates.
    REAL(rp), DIMENSION(8), INTENT(IN)     :: B_X,B_Y,B_Z
    !! Magnetic field \(\mathbf{B}\) seen by each particle. This is given
    !! in Cartesian coordinates.
    REAL(rp), DIMENSION(8), INTENT(OUT)    :: Frad_X,Frad_Y,Frad_Z
    !! The calculated synchrotron radiation reaction force \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(3)                 :: F1
    !! The component \(\mathbf{F}_1\) of \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(8)                 :: F2_X,F2_Y,F2_Z
    !! The component \(\mathbf{F}_2\) of \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(8)                 :: F3_X,F3_Y,F3_Z
    !! The component \(\mathbf{F}_3\) of \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(8)                 :: V_X,V_Y,V_Z
    !! The particle's velocity \(\mathbf{v}\).
    REAL(rp), DIMENSION(8)                 :: vec_X,vec_Y,vec_Z
    REAL(rp), DIMENSION(8)                 :: cross_EB_X,cross_EB_Y,cross_EB_Z
    REAL(rp), DIMENSION(8)                 :: cross_BV_X,cross_BV_Y,cross_BV_Z
    REAL(rp), DIMENSION(8)                 :: cross_BBV_X,cross_BBV_Y,cross_BBV_Z
    REAL(rp), DIMENSION(8)                 :: dot_EV,dot_vecvec
    !! An auxiliary 3-D vector.
    REAL(rp),DIMENSION(8)                               :: g
    !! The relativistic \(\gamma\) factor of the particle.
    REAL(rp)                               :: tmp
    INTEGER :: cc

    !$OMP SIMD
    !    !$OMP& aligned(g,U_X,U_Y,U_Z,V_X,V_Y,V_Z, &
    !    !$OMP& cross_EB_X,cross_EB_Y,cross_EB_Z,E_X,E_Y,E_Z,B_X,B_Y,B_Z, &
    !    !$OMP& dot_EV,cross_BV_X,cross_BV_Y,cross_BV_Z, &
    !    !$OMP& cross_BBV_X,cross_BBV_Y,cross_BBV_Z,F2_X,F2_Y,F2_Z, &
    !    !$OMP& vec_X,vec_Y,vec_Z,dot_vecvec,F3_X,F3_Y,F3_Z, &
    !    !$OMP& Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,8_idef
       g(cc) = SQRT(1.0_rp + U_X(cc)*U_X(cc)+ U_Y(cc)*U_Y(cc)+ U_Z(cc)*U_Z(cc))
       
       V_X(cc) = U_X(cc)/g(cc)
       V_Y(cc) = U_Y(cc)/g(cc)
       V_Z(cc) = U_Z(cc)/g(cc)

       tmp = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

       cross_EB_X(cc)=E_Y(cc)*B_Z(cc)-E_Z(cc)*B_Y(cc)
       cross_EB_Y(cc)=E_Z(cc)*B_X(cc)-E_X(cc)*B_Z(cc)
       cross_EB_Z(cc)=E_X(cc)*B_Y(cc)-E_Y(cc)*B_X(cc)

       dot_EV(cc)=E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc)

       cross_BV_X(cc)=B_Y(cc)*V_Z(cc)-B_Z(cc)*V_Y(cc)
       cross_BV_Y(cc)=B_Z(cc)*V_X(cc)-B_X(cc)*V_Z(cc)
       cross_BV_Z(cc)=B_X(cc)*V_Y(cc)-B_Y(cc)*V_X(cc)

       cross_BBV_X(cc)=B_Y(cc)*cross_BV_Z(cc)-B_Z(cc)*cross_BV_Y(cc)
       cross_BBV_Y(cc)=B_Z(cc)*cross_BV_X(cc)-B_X(cc)*cross_BV_Z(cc)
       cross_BBV_Z(cc)=B_X(cc)*cross_BV_Y(cc)-B_Y(cc)*cross_BV_X(cc)
       
       F2_X(cc) = tmp*( dot_EV(cc)*E_X(cc) + cross_EB_X(cc) + cross_BBV_X(cc) )
       F2_Y(cc) = tmp*( dot_EV(cc)*E_Y(cc) + cross_EB_Y(cc) + cross_BBV_Y(cc) )
       F2_Z(cc) = tmp*( dot_EV(cc)*E_Z(cc) + cross_EB_Z(cc) + cross_BBV_Z(cc) )
       
       vec_X(cc) = E_X(cc) - cross_BV_X(cc)
       vec_Y(cc) = E_Y(cc) - cross_BV_Y(cc)
       vec_Z(cc) = E_Z(cc) - cross_BV_Z(cc)

       dot_vecvec(cc)=vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+vec_Z(cc)*vec_Z(cc)
       
       F3_X(cc) = (tmp*g(cc)**2)*( dot_EV(cc)**2 - dot_vecvec(cc) )*V_X(cc)
       F3_Y(cc) = (tmp*g(cc)**2)*( dot_EV(cc)**2 - dot_vecvec(cc) )*V_Y(cc)
       F3_Z(cc) = (tmp*g(cc)**2)*( dot_EV(cc)**2 - dot_vecvec(cc) )*V_Z(cc)

       Frad_X(cc) = F2_X(cc) + F3_X(cc)
       Frad_Y(cc) = F2_Y(cc) + F3_Y(cc)
       Frad_Z(cc) = F2_Z(cc) + F3_Z(cc)
       
    end do
    !$OMP END SIMD
    
  end subroutine radiation_force_p




  subroutine FO_init(params,F,spp,output,step)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.

    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp)                                                   :: Prad
    !! Total radiated power of each particle.

    REAL(rp)                                  :: Bmag1
    !! Magnitude of the magnetic field seen by each particle .
    REAL(rp)                                                   :: v
    !! Speed of each particle.
    REAL(rp)                                                   :: vpar
    !! Parallel velocity \(v_\parallel = \mathbf{v}\cdot \hat{b}\).
    REAL(rp)                                                   :: vperp
    !! Perpendicular velocity \(v_\parallel = |\mathbf{v} - (\mathbf{v}\cdot
    !! \hat{b})\hat{b}|\).
    REAL(rp)                                                   :: tmp
    !! Temporary variable used for various computations.
    REAL(rp)                                                   :: a
    !! This variable is used to simplify notation in the code, and
    !! is given by \(a=q\Delta t/m\),

    REAL(rp), DIMENSION(3)                       :: Frad
    !! Synchrotron radiation reaction force of each particle.
    REAL(rp), DIMENSION(3)                       :: vec
    !! Auxiliary vector used in various computations.
    REAL(rp), DIMENSION(3)                       :: b_unit
    !! Unitary vector pointing along the local magnetic field \(\hat{b}\).
    INTEGER                                      :: ii
    !! Species iterator.
    INTEGER                                      :: pp
    !! Particles iterator.
    INTEGER                                      :: cc
    !! Chunk iterator.

    LOGICAL,intent(in) :: output
    LOGICAL,intent(in) :: step   
    
    INTEGER ,DIMENSION(8)                                     :: flag_cache


    do ii = 1_idef,params%num_species

       if(output) then
       
          call get_fields(params,spp(ii)%vars,F)
          !! Calls [[get_fields]] in [[korc_fields]].
          ! Interpolates fields at local particles' position and keeps in
          ! spp%vars. Fields in (R,\(\phi\),Z) coordinates.

!          write(6,'("korc_ppusher")')
!          write(6,'("B_X: ",E17.10)') spp(ii)%vars%B(:,1)
!          write(6,'("B_Z: ",E17.10)') spp(ii)%vars%B(:,2)
!          write(6,'("B_Y: ",E17.10)') spp(ii)%vars%B(:,3)
          
          !$OMP PARALLEL DO DEFAULT(none) SHARED(ii,spp) &
          !$OMP& FIRSTPRIVATE(E0) &
          !$OMP& PRIVATE(pp,b_unit,Bmag1,vpar,v,vperp,vec,tmp)
          do pp=1_idef,spp(ii)%ppp

             Bmag1 = SQRT(DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                  spp(ii)%vars%B(pp,:)))

             ! Parallel unit vector
             b_unit = spp(ii)%vars%B(pp,:)/Bmag1

             v = SQRT(DOT_PRODUCT(spp(ii)%vars%V(pp,:),spp(ii)%vars%V(pp,:)))
             if (v.GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar = DOT_PRODUCT(spp(ii)%vars%V(pp,:), b_unit)
                vperp =  DOT_PRODUCT(spp(ii)%vars%V(pp,:), &
                     spp(ii)%vars%V(pp,:)) &
                     - vpar**2
                if ( vperp .GE. korc_zero ) then
                   vperp = SQRT( vperp )
                else
                   vperp = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp) = 180.0_rp*MODULO(ATAN2(vperp,vpar), &
                     2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp) = 0.5_rp*spp(ii)%m* &
                     spp(ii)%vars%g(pp)**2*vperp**2/Bmag1
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)

                ! Radiated power
                tmp = spp(ii)%q**4/(6.0_rp*C_PI*E0*spp(ii)%m**2)
                vec = spp(ii)%vars%E(pp,:) + cross(spp(ii)%vars%V(pp,:), &
                     spp(ii)%vars%B(pp,:))

                spp(ii)%vars%Prad(pp) = tmp*( DOT_PRODUCT(spp(ii)% &
                     vars%E(pp,:), &
                     spp(ii)%vars%E(pp,:)) + &
                     DOT_PRODUCT(cross(spp(ii)%vars%V(pp,:), &
                     spp(ii)%vars%B(pp,:)),spp(ii)%vars%E(pp,:))+ &
                     spp(ii)%vars%g(pp)**2* &
                     (DOT_PRODUCT(spp(ii)%vars%E(pp,:), &
                     spp(ii)%vars%V(pp,:))**2 - DOT_PRODUCT(vec,vec)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp) = spp(ii)%q*DOT_PRODUCT( &
                     spp(ii)%vars%E(pp,:),spp(ii)%vars%V(pp,:))
             else
                spp(ii)%vars%eta(pp) = 0.0_rp
                spp(ii)%vars%mu(pp) = 0.0_rp
                spp(ii)%vars%Prad(pp) = 0.0_rp
                spp(ii)%vars%Pin(pp) = 0.0_rp
             end if


          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

       end if !(if output)

       if(step.and.(.not.params%FokPlan)) then
          dt=0.5_rp*params%dt
          
          !$OMP PARALLEL DO FIRSTPRIVATE(dt) PRIVATE(pp,cc) &
          !$OMP& SHARED(ii,spp,params)
          do pp=1_idef,spp(ii)%ppp,8

             !$OMP SIMD
             do cc=1_idef,8
                spp(ii)%vars%X(pp-1+cc,1) = spp(ii)%vars%X(pp-1+cc,1) + &
                     dt*spp(ii)%vars%V(pp-1+cc,1)
                spp(ii)%vars%X(pp-1+cc,2) = spp(ii)%vars%X(pp-1+cc,2) + &
                     dt*spp(ii)%vars%V(pp-1+cc,2)
                spp(ii)%vars%X(pp-1+cc,3) = spp(ii)%vars%X(pp-1+cc,3) + &
                     dt*spp(ii)%vars%V(pp-1+cc,3)
             end do
             !$OMP END SIMD
             
          end do
          !$OMP END PARALLEL DO

       end if !(if step)

    end do ! over species

  end subroutine FO_init

  subroutine adv_FOeqn_top(params,F,P,spp)
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(8)               :: Bmag
    REAL(rp), DIMENSION(8)               :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(8)               :: v,vpar,vperp
    REAL(rp), DIMENSION(8)               :: tmp
    REAL(rp), DIMENSION(8)               :: g
    REAL(rp), DIMENSION(8)               :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(8)               :: vec_X,vec_Y,vec_Z
    REAL(rp),DIMENSION(8) :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(8) :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(8) :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(8) :: E_X,E_Y,E_Z
    INTEGER(is),DIMENSION(8) :: flag_cache

    REAL(rp) :: B0,EF0,R0,q0,lam,ar
    REAL(rp) :: a,m_cache,q_cache
    REAL(rp) :: ne0,Te0,Zeff0


    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
 

    do ii = 1_idef,params%num_species      

       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a = q_cache*params%dt/m_cache
       
       B0=F%Bo
       EF0=F%Eo
       lam=F%AB%lambda
       R0=F%AB%Ro
       q0=F%AB%qo
       ar=F%AB%a


       
       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(E0,a,m_cache,q_cache,B0,EF0,lam,R0,q0,ar)&
       !$OMP& shared(params,ii,spp,P) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g,flag_cache)
       do pp=1_idef,spp(ii)%ppp,8

          !$OMP SIMD
          do cc=1_idef,8_idef
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             g(cc)=spp(ii)%vars%g(pp-1+cc)
             flag_cache(cc)=spp(ii)%vars%flag(pp-1+cc)
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then
             do tt=1_ip,params%t_skip

                call analytical_fields_p(B0,EF0,R0,q0,lam,ar,X_X,X_Y,X_Z, &
                     B_X,B_Y,B_Z,E_X,E_Y,E_Z,flag_cache)

                call advance_FOeqn_vars(tt,a,q_cache,m_cache,params, &
                     X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                     P,g,flag_cache)
             end do !timestep iterator

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
                spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
                spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)
                
                spp(ii)%vars%flag(pp-1+cc) = flag_cache(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
                spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
                spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)
             end do
             !$OMP END SIMD

          else

             !$OMP SIMD
             do cc=1_idef,8_idef
                B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
                E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
                E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)
             end do
             !$OMP END SIMD
             
             call advance_FP3Deqn_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
                  g,m_cache,B0,lam,R0,q0,EF0,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                  P,flag_cache)

             !$OMP SIMD
             do cc=1_idef,8_idef

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)

             end do
             !$OMP END SIMD
             
          end if
          
          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,8_idef
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))
                
                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)

                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)
                
                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)-vec_Y(cc)*vec_Y(cc)- &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD

             
       end do !particle chunk iterator
       !$OMP END PARALLEL DO
       
    end do !species iterator
    
  end subroutine adv_FOeqn_top
  
  subroutine advance_FOeqn_vars(tt,a,q_cache,m_cache,params,X_X,X_Y,X_Z, &
       V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,P,g,flag_cache)
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.

    INTEGER(ip), INTENT(IN)                                       :: tt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp), INTENT(IN)                       :: m_cache,q_cache
    !! Time step used in the leapfrog step (\(\Delta t\)).

    REAL(rp),DIMENSION(8)                                  :: Bmag



    REAL(rp),INTENT(in)                                       :: a
    !! This variable is used to simplify notation in the code, and
    !! is given by \(a=q\Delta t/m\),
    REAL(rp),DIMENSION(8)                                    :: sigma
    !! This variable is \(\sigma = \gamma'^2 - \tau^2\) in the above equations.
    REAL(rp),DIMENSION(8)                               :: us
    !! This variable is \(u^{*} = p^{*}/m\) where \( p^{*} =
    !! \mathbf{p}'\cdot \mathbf{\tau}/mc\).
    !! Variable 'u^*' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(8),INTENT(INOUT)                 :: g
    REAL(rp),DIMENSION(8) :: gp,g0
    !! Relativistic factor \(\gamma\).
    REAL(rp),DIMENSION(8)                                 :: s
    !! This variable is \(s = 1/(1+t^2)\) in the equations above.
    !! Variable 's' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(8)                            :: U_hs_X,U_hs_Y,U_hs_Z
    !! Is \(\mathbf{u}=\mathbf{p}/m\) at half-time step (\(i+1/2\)) in
    !! the absence of radiation losses or collisions. \(\mathbf{u}^{i+1/2} =
    !! \mathbf{u}^i + \frac{q\Delta t}{2m}\left( \mathbf{E}^{i+1/2} +
    !! \mathbf{v}^i\times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(8)                           :: tau_X,tau_Y,tau_Z
    !! This variable is \(\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}\).
    REAL(rp),DIMENSION(8)                            :: up_X,up_Y,up_Z
    !! This variable is \(\mathbf{u}'= \mathbf{p}'/m\), where \(\mathbf{p}'
    !! = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} +
    !! \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(8)                                     :: t_X,t_Y,t_Z
    !! This variable is \(\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}\).
    REAL(rp),DIMENSION(8),INTENT(INOUT)                     :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT)                      :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(8),INTENT(IN)                      :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(8),INTENT(IN)                     :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(8)                     :: U_L_X,U_L_Y,U_L_Z
    REAL(rp),DIMENSION(8)                     :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(8)                     :: U_RC_X,U_RC_Y,U_RC_Z
    REAL(rp),DIMENSION(8)                     :: U_os_X,U_os_Y,U_os_Z
    !! This variable is \(\mathbf{u}^{i+1}= \mathbf{p}^{i+1}/m\).
    REAL(rp),DIMENSION(8)                          :: cross_X,cross_Y,cross_Z

    REAL(rp), DIMENSION(8)                       :: Frad_X,Frad_Y,Frad_Z
    !! Synchrotron radiation reaction force of each particle.

    REAL(rp),DIMENSION(8) :: ne,Te,Zeff,Y_R,Y_PHI,Y_Z

    INTEGER                                      :: cc
    !! Chunk iterator.

    INTEGER(is),DIMENSION(8),intent(inout)             :: flag_cache

    dt=params%dt
    
    
    !$OMP SIMD
    !    !$OMP& aligned(g0,g,U_X,U_Y,U_Z,V_X,V_Y,V_Z,Bmag,B_X,B_Y,B_Z, &
    !    !$OMP& U_L_X,U_L_Y,U_L_Z,U_RC_X,U_RC_Y,U_RC_Z, &
    !    !$OMP& cross_X,cross_Y,cross_Z,U_hs_X,U_hs_Y,U_hs_Z,E_X,E_Y,E_Z, &
    !    !$OMP& tau_X,tau_Y,tau_Z,up_X,up_Y,up_Z,gp,sigma,us,t_X,t_Y,t_Z,s, &
    !    !$OMP& U_os_X,U_os_Y,U_os_Z,Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,8

       g0(cc)=g(cc)
       
       U_X(cc) = g(cc)*V_X(cc)
       U_Y(cc) = g(cc)*V_Y(cc)
       U_Z(cc) = g(cc)*V_Z(cc)
       

       ! Magnitude of magnetic field
       Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

       U_L_X(cc)=U_X(cc)
       U_L_Y(cc)=U_Y(cc)
       U_L_Z(cc)=U_Z(cc)

       U_RC_X(cc)=U_X(cc)
       U_RC_Y(cc)=U_Y(cc)
       U_RC_Z(cc)=U_Z(cc)
       
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       cross_X(cc)=V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
       cross_Y(cc)=V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
       cross_Z(cc)=V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)


       
       U_hs_X(cc) = U_L_X(cc) + 0.5_rp*a*(E_X(cc) +cross_X(cc))
       U_hs_Y(cc) = U_L_Y(cc) + 0.5_rp*a*(E_Y(cc) +cross_Y(cc))
       U_hs_Z(cc) = U_L_Z(cc) + 0.5_rp*a*(E_Z(cc) +cross_Z(cc))


       
       tau_X(cc) = 0.5_rp*a*B_X(cc)
       tau_Y(cc) = 0.5_rp*a*B_Y(cc)
       tau_Z(cc) = 0.5_rp*a*B_Z(cc)


       
       up_X(cc) = U_hs_X(cc) + 0.5_rp*a*E_X(cc)
       up_Y(cc) = U_hs_Y(cc) + 0.5_rp*a*E_Y(cc)
       up_Z(cc) = U_hs_Z(cc) + 0.5_rp*a*E_Z(cc)

       gp(cc) = SQRT( 1.0_rp + up_X(cc)*up_X(cc)+up_Y(cc)*up_Y(cc)+ &
            up_Z(cc)*up_Z(cc) )

       sigma(cc) = gp(cc)*gp(cc) - (tau_X(cc)*tau_X(cc)+ &
            tau_Y(cc)*tau_Y(cc)+tau_Z(cc)*tau_Z(cc))

       us(cc) = up_X(cc)*tau_X(cc)+up_Y(cc)*tau_Y(cc)+ &
            up_Z(cc)*tau_Z(cc)

       ! variable 'u^*' in Vay, J.-L. PoP (2008)
       g(cc) = SQRT( 0.5_rp*(sigma(cc) + SQRT(sigma(cc)*sigma(cc) + &
            4.0_rp*(tau_X(cc)*tau_X(cc)+tau_Y(cc)*tau_Y(cc)+ &
            tau_Z(cc)*tau_Z(cc) + us(cc)*us(cc)))) )

       t_X(cc) = tau_X(cc)/g(cc)
       t_Y(cc) = tau_Y(cc)/g(cc)
       t_Z(cc) = tau_Z(cc)/g(cc)

       
       s(cc) = 1.0_rp/(1.0_rp + t_X(cc)*t_X(cc)+t_Y(cc)*t_Y(cc)+ &
            t_Z(cc)*t_Z(cc))
       ! variable 's' in Vay, J.-L. PoP (2008)

       cross_X(cc)=up_Y(cc)*t_Z(cc)-up_Z(cc)*t_Y(cc)
       cross_Y(cc)=up_Z(cc)*t_X(cc)-up_X(cc)*t_Z(cc)
       cross_Z(cc)=up_X(cc)*t_Y(cc)-up_Y(cc)*t_X(cc)

       U_L_X(cc) = s(cc)*(up_X(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_X(cc) + cross_X(cc))
       U_L_Y(cc) = s(cc)*(up_Y(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Y(cc) + cross_Y(cc))
       U_L_Z(cc) = s(cc)*(up_Z(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Z(cc) + cross_Z(cc))
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !      

       U_os_X(cc) = 0.5_rp*(U_L_X(cc) + U_X(cc))
       U_os_Y(cc) = 0.5_rp*(U_L_Y(cc) + U_Y(cc))
       U_os_Z(cc) = 0.5_rp*(U_L_Z(cc) + U_Z(cc))
       ! Splitting operator for including radiation

       if (params%radiation) then
          !! Calls [[radiation_force]] in [[korc_ppusher]].
          call radiation_force_p(q_cache,m_cache,U_os_X,U_os_Y,U_os_Z, &
               E_X,E_Y,E_Z,B_Z,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)
          U_RC_X(cc) = U_RC_X(cc) + a*Frad_X(cc)/q_cache
          U_RC_Y(cc) = U_RC_Y(cc) + a*Frad_Y(cc)/q_cache
          U_RC_Z(cc) = U_RC_Z(cc) + a*Frad_Z(cc)/q_cache
       end if
       ! Splitting operator for including radiation

       U_X(cc) = U_L_X(cc) + U_RC_X(cc) - U_X(cc)
       U_Y(cc) = U_L_Y(cc) + U_RC_Y(cc) - U_Y(cc)
       U_Z(cc) = U_L_Z(cc) + U_RC_Z(cc) - U_Z(cc)
       
    end do
    !$OMP END SIMD
   

    if (params%collisions) then

       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,flag_cache)
          
    end if

    if (params%radiation.or.params%collisions) then

       !$OMP SIMD
!       !$OMP& aligned(g,U_X,U_Y,U_Z)
       do cc=1_idef,8_idef
          g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
       end do
       !$OMP END SIMD
       
    end if
    
    !$OMP SIMD
!    !$OMP& aligned(g,g0,V_X,V_Y,V_Z,U_X,U_Y,U_Z,X_X,X_Y,X_Z,flag_cache)
    do cc=1_idef,8_idef

       if (flag_cache(cc).eq.0_is) then
          g(cc)=g0(cc)
       else
          V_X(cc) = U_X(cc)/g(cc)
          V_Y(cc) = U_Y(cc)/g(cc)
          V_Z(cc) = U_Z(cc)/g(cc)
       end if

       X_X(cc) = X_X(cc) + dt*V_X(cc)*REAL(flag_cache(cc))
       X_Y(cc) = X_Y(cc) + dt*V_Y(cc)*REAL(flag_cache(cc))
       X_Z(cc) = X_Z(cc) + dt*V_Z(cc)*REAL(flag_cache(cc))
    end do
    !$OMP END SIMD
    
  end subroutine advance_FOeqn_vars

  subroutine advance_FP3Deqn_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z,g, &
       m_cache,B0,lam,R0,q0,EF0,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
       P,flag_cache)
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(8), INTENT(IN)  :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(8), INTENT(IN)  :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(8), INTENT(IN)  :: B_X,B_Y,B_Z
    INTEGER(is),DIMENSION(8), INTENT(INOUT)  :: flag_cache
    REAL(rp),DIMENSION(8) :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(8), INTENT(INOUT)  :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: g
    REAL(rp),intent(in) :: B0,EF0,R0,q0,lam,m_cache
    


!    call analytical_fields_p(B0,EF0,R0,q0,lam,X_X,X_Y,X_Z, &
!         B_X,B_Y,B_Z,E_X,E_Y,E_Z)

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,8_idef
       U_X(cc)=V_X(cc)*g(cc)
       U_Y(cc)=V_Y(cc)*g(cc)
       U_Z(cc)=V_Z(cc)*g(cc)
    end do
    !$OMP END SIMD
    
    do tt=1_ip,params%t_skip
          
       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,flag_cache)
       
    end do

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,8_idef

       g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
          
       V_X(cc)=U_X(cc)/g(cc)
       V_Y(cc)=U_Y(cc)/g(cc)
       V_Z(cc)=U_Z(cc)/g(cc)
    end do
    !$OMP END SIMD

  end subroutine advance_FP3Deqn_vars
  
  subroutine adv_FOinterp_top(params,F,P,spp)  
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(8)               :: Bmag
    REAL(rp), DIMENSION(8)               :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(8)               :: v,vpar,vperp
    REAL(rp), DIMENSION(8)               :: tmp
    REAL(rp), DIMENSION(8)               :: g
    REAL(rp), DIMENSION(8)               :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(8)               :: vec_X,vec_Y,vec_Z
    REAL(rp),DIMENSION(8) :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(8) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8) :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(8) :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(8) :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(8) :: PSIp
    INTEGER(is),DIMENSION(8) :: flag_cache
    REAL(rp) :: a,m_cache,q_cache    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
 

    do ii = 1_idef,params%num_species      

       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a = q_cache*params%dt/m_cache


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(a,m_cache,q_cache) &
       !$OMP& shared(params,ii,spp,P,F) &
       !$OMP& PRIVATE(E0,pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g, &
       !$OMP& Y_R,Y_PHI,Y_Z,flag_cache,PSIp)
       do pp=1_idef,spp(ii)%ppp,8

          !$OMP SIMD
          do cc=1_idef,8_idef
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             g(cc)=spp(ii)%vars%g(pp-1+cc)
             
             flag_cache(cc)=spp(ii)%vars%flag(pp-1+cc)
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then
             do tt=1_ip,params%t_skip

                call cart_to_cyl_p(X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

                if (params%orbit_model(3:5).eq.'new') then
                   call interp_FOfields_p(Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
                        E_X,E_Y,E_Z,PSIp,flag_cache)
                else if (params%orbit_model(3:5).eq.'old') then
                   call interp_FOfields1_p(F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
                        E_X,E_Y,E_Z,PSIp,flag_cache)
                end if

                
 !               write(6,'("B_X: ",E17.10)') B_X(1)
 !               write(6,'("B_Y: ",E17.10)') B_Y(1)
 !               write(6,'("B_Z: ",E17.10)') B_Z(1)
                
                call advance_FOinterp_vars(tt,a,q_cache,m_cache,params, &
                     X_X,X_Y,X_Z, &
                     V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,g,flag_cache,P)
             end do !timestep iterator

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
                spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
                spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)
                spp(ii)%vars%flag(pp-1+cc) = flag_cache(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
                spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
                spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)

                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
             end do
             !$OMP END SIMD

          else
             !$OMP SIMD
             do cc=1_idef,8_idef
                B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
                E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
                E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)
             end do
             !$OMP END SIMD
             
             call advance_FP3Dinterp_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
                  g,m_cache,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flag_cache,P)

             !$OMP SIMD
             do cc=1_idef,8_idef

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)

             end do
             !$OMP END SIMD
          end if

          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,8_idef
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))
                
                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)

                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)
                
                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+ &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD

             
       end do !particle chunk iterator
       !$OMP END PARALLEL DO
       
    end do !species iterator
    
  end subroutine adv_FOinterp_top
  
  subroutine advance_FOinterp_vars(tt,a,q_cache,m_cache,params,X_X,X_Y,X_Z, &
       V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,g,flag_cache,P)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    INTEGER(ip), INTENT(IN)                                       :: tt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp), INTENT(IN)                                       :: m_cache,q_cache
    !! Time step used in the leapfrog step (\(\Delta t\)).

    REAL(rp),DIMENSION(8)                                  :: Bmag



    REAL(rp),INTENT(in)                                       :: a
    !! This variable is used to simplify notation in the code, and
    !! is given by \(a=q\Delta t/m\),
    REAL(rp),DIMENSION(8)                                    :: sigma
    !! This variable is \(\sigma = \gamma'^2 - \tau^2\) in the above equations.
    REAL(rp),DIMENSION(8)                               :: us
    !! This variable is \(u^{*} = p^{*}/m\) where \( p^{*} =
    !! \mathbf{p}'\cdot \mathbf{\tau}/mc\).
    !! Variable 'u^*' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(8),INTENT(INOUT)                 :: g
    REAL(rp),DIMENSION(8) :: gp,g0
    !! Relativistic factor \(\gamma\).
    REAL(rp),DIMENSION(8)                                 :: s
    !! This variable is \(s = 1/(1+t^2)\) in the equations above.
    !! Variable 's' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(8)                            :: U_hs_X,U_hs_Y,U_hs_Z
    !! Is \(\mathbf{u}=\mathbf{p}/m\) at half-time step (\(i+1/2\)) in
    !! the absence of radiation losses or collisions. \(\mathbf{u}^{i+1/2} =
    !! \mathbf{u}^i + \frac{q\Delta t}{2m}\left( \mathbf{E}^{i+1/2} +
    !! \mathbf{v}^i\times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(8)                           :: tau_X,tau_Y,tau_Z
    !! This variable is \(\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}\).
    REAL(rp),DIMENSION(8)                            :: up_X,up_Y,up_Z
    !! This variable is \(\mathbf{u}'= \mathbf{p}'/m\), where \(\mathbf{p}'
    !! = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} +
    !! \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(8)                                     :: t_X,t_Y,t_Z
    !! This variable is \(\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}\).
    REAL(rp),DIMENSION(8),INTENT(INOUT)                     :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(8)                    :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT)                      :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(8),INTENT(IN)                      :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(8),INTENT(IN)                     :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(8)                     :: U_L_X,U_L_Y,U_L_Z
    REAL(rp),DIMENSION(8)                     :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(8)                     :: U_RC_X,U_RC_Y,U_RC_Z
    REAL(rp),DIMENSION(8)                     :: U_os_X,U_os_Y,U_os_Z
    !! This variable is \(\mathbf{u}^{i+1}= \mathbf{p}^{i+1}/m\).
    REAL(rp),DIMENSION(8)                          :: cross_X,cross_Y,cross_Z

    REAL(rp), DIMENSION(8)                       :: Frad_X,Frad_Y,Frad_Z
    !! Synchrotron radiation reaction force of each particle.

    REAL(rp),DIMENSION(8) :: ne,Te,Zeff

    INTEGER                                      :: cc
    !! Chunk iterator.

    INTEGER(is) ,DIMENSION(8),intent(inout)                   :: flag_cache

    dt=params%dt
    
    
    !$OMP SIMD
    !    !$OMP& aligned(g0,g,U_X,U_Y,U_Z,V_X,V_Y,V_Z,Bmag,B_X,B_Y,B_Z, &
    !    !$OMP& U_L_X,U_L_Y,U_L_Z,U_RC_X,U_RC_Y,U_RC_Z, &
    !    !$OMP& cross_X,cross_Y,cross_Z,U_hs_X,U_hs_Y,U_hs_Z,E_X,E_Y,E_Z, &
    !    !$OMP& tau_X,tau_Y,tau_Z,up_X,up_Y,up_Z,gp,sigma,us,t_X,t_Y,t_Z,s, &
    !    !$OMP& U_os_X,U_os_Y,U_os_Z,Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,8

       g0(cc)=g(cc)
       
       U_X(cc) = g(cc)*V_X(cc)
       U_Y(cc) = g(cc)*V_Y(cc)
       U_Z(cc) = g(cc)*V_Z(cc)
       

       ! Magnitude of magnetic field
       Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

       U_L_X(cc)=U_X(cc)
       U_L_Y(cc)=U_Y(cc)
       U_L_Z(cc)=U_Z(cc)

       U_RC_X(cc)=U_X(cc)
       U_RC_Y(cc)=U_Y(cc)
       U_RC_Z(cc)=U_Z(cc)

       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       cross_X(cc)=V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
       cross_Y(cc)=V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
       cross_Z(cc)=V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)


       
       U_hs_X(cc) = U_L_X(cc) + 0.5_rp*a*(E_X(cc) +cross_X(cc))
       U_hs_Y(cc) = U_L_Y(cc) + 0.5_rp*a*(E_Y(cc) +cross_Y(cc))
       U_hs_Z(cc) = U_L_Z(cc) + 0.5_rp*a*(E_Z(cc) +cross_Z(cc))


       
       tau_X(cc) = 0.5_rp*a*B_X(cc)
       tau_Y(cc) = 0.5_rp*a*B_Y(cc)
       tau_Z(cc) = 0.5_rp*a*B_Z(cc)


       
       up_X(cc) = U_hs_X(cc) + 0.5_rp*a*E_X(cc)
       up_Y(cc) = U_hs_Y(cc) + 0.5_rp*a*E_Y(cc)
       up_Z(cc) = U_hs_Z(cc) + 0.5_rp*a*E_Z(cc)

       gp(cc) = SQRT( 1.0_rp + up_X(cc)*up_X(cc)+up_Y(cc)*up_Y(cc)+ &
            up_Z(cc)*up_Z(cc) )

       sigma(cc) = gp(cc)*gp(cc) - (tau_X(cc)*tau_X(cc)+ &
            tau_Y(cc)*tau_Y(cc)+tau_Z(cc)*tau_Z(cc))

       us(cc) = up_X(cc)*tau_X(cc)+up_Y(cc)*tau_Y(cc)+ &
            up_Z(cc)*tau_Z(cc)

       ! variable 'u^*' in Vay, J.-L. PoP (2008)
       g(cc) = SQRT( 0.5_rp*(sigma(cc) + SQRT(sigma(cc)*sigma(cc) + &
            4.0_rp*(tau_X(cc)*tau_X(cc)+tau_Y(cc)*tau_Y(cc)+ &
            tau_Z(cc)*tau_Z(cc) + us(cc)*us(cc)))) )

       t_X(cc) = tau_X(cc)/g(cc)
       t_Y(cc) = tau_Y(cc)/g(cc)
       t_Z(cc) = tau_Z(cc)/g(cc)

       
       s(cc) = 1.0_rp/(1.0_rp + t_X(cc)*t_X(cc)+t_Y(cc)*t_Y(cc)+ &
            t_Z(cc)*t_Z(cc))
       ! variable 's' in Vay, J.-L. PoP (2008)

       cross_X(cc)=up_Y(cc)*t_Z(cc)-up_Z(cc)*t_Y(cc)
       cross_Y(cc)=up_Z(cc)*t_X(cc)-up_X(cc)*t_Z(cc)
       cross_Z(cc)=up_X(cc)*t_Y(cc)-up_Y(cc)*t_X(cc)

       U_L_X(cc) = s(cc)*(up_X(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_X(cc) + cross_X(cc))
       U_L_Y(cc) = s(cc)*(up_Y(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Y(cc) + cross_Y(cc))
       U_L_Z(cc) = s(cc)*(up_Z(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Z(cc) + cross_Z(cc))
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       U_os_X(cc) = 0.5_rp*(U_L_X(cc) + U_X(cc))
       U_os_Y(cc) = 0.5_rp*(U_L_Y(cc) + U_Y(cc))
       U_os_Z(cc) = 0.5_rp*(U_L_Z(cc) + U_Z(cc))
       ! Splitting operator for including radiation

       if (params%radiation) then
          !! Calls [[radiation_force]] in [[korc_ppusher]].
          call radiation_force_p(q_cache,m_cache,U_os_X,U_os_Y,U_os_Z, &
               E_X,E_Y,E_Z,B_Z,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)
          U_RC_X(cc) = U_RC_X(cc) + a*Frad_X(cc)/q_cache
          U_RC_Y(cc) = U_RC_Y(cc) + a*Frad_Y(cc)/q_cache
          U_RC_Z(cc) = U_RC_Z(cc) + a*Frad_Z(cc)/q_cache
       end if
       ! Splitting operator for including radiation

       U_X(cc) = U_L_X(cc) + U_RC_X(cc) - U_X(cc)
       U_Y(cc) = U_L_Y(cc) + U_RC_Y(cc) - U_Y(cc)
       U_Z(cc) = U_L_Z(cc) + U_RC_Z(cc) - U_Z(cc)
       
    end do
    !$OMP END SIMD
   
    if (params%collisions) then
       
       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,flag_cache)
       
    end if

    if (params%radiation.or.params%collisions) then

       !$OMP SIMD
       !       !$OMP& aligned(g,U_X,U_Y,U_Z)
       do cc=1_idef,8_idef
          g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
       end do
       !$OMP END SIMD
       
    end if
    
    !$OMP SIMD
    !    !$OMP& aligned(g,g0,V_X,V_Y,V_Z,U_X,U_Y,U_Z,X_X,X_Y,X_Z,flag_cache)
    do cc=1_idef,8_idef
       
       if (flag_cache(cc).eq.0_is) then
          g(cc)=g0(cc)
       else
          V_X(cc) = U_X(cc)/g(cc)
          V_Y(cc) = U_Y(cc)/g(cc)
          V_Z(cc) = U_Z(cc)/g(cc)
       end if

       X_X(cc) = X_X(cc) + dt*V_X(cc)*REAL(flag_cache(cc))
       X_Y(cc) = X_Y(cc) + dt*V_Y(cc)*REAL(flag_cache(cc))
       X_Z(cc) = X_Z(cc) + dt*V_Z(cc)*REAL(flag_cache(cc))
       
    end do
    !$OMP END SIMD
    
  end subroutine advance_FOinterp_vars

  subroutine advance_FP3Dinterp_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z,g, &
       m_cache,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flag_cache,P)    
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(8), INTENT(IN)  :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(8)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8), INTENT(IN)  :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(8), INTENT(IN)  :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(8) :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(8), INTENT(INOUT)  :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(8) :: ne,Te,Zeff
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: g
    INTEGER(is),DIMENSION(8),INTENT(INOUT) :: flag_cache
    REAL(rp),intent(in) :: m_cache
    
      
    
    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,8_idef
       U_X(cc)=V_X(cc)*g(cc)
       U_Y(cc)=V_Y(cc)*g(cc)
       U_Z(cc)=V_Z(cc)*g(cc)
    end do
    !$OMP END SIMD
    
    do tt=1_ip,params%t_skip
          
       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,flag_cache)
       
    end do

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,8_idef

       g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
          
       V_X(cc)=U_X(cc)/g(cc)
       V_Y(cc)=U_Y(cc)/g(cc)
       V_Z(cc)=U_Z(cc)/g(cc)
    end do
    !$OMP END SIMD

  end subroutine advance_FP3Dinterp_vars
  

  subroutine GC_init(params,F,spp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.

    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.

    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    REAL(rp)               :: Bmag1,pmag
    REAL(rp)               :: Bmagc
    REAL(rp)               :: rm
    REAL(rp),DIMENSION(:,:),ALLOCATABLE               :: RAphi
    REAL(rp), DIMENSION(3) :: bhat
    REAL(rp), DIMENSION(3) :: bhatc
    REAL(rp), DIMENSION(8) :: E_PHI
    REAL(rp),DIMENSION(:),ALLOCATABLE               :: RVphi

!    write(6,'("eta",E17.10)') spp(ii)%vars%eta(pp)
!    write(6,'("gam",E17.10)') spp(ii)%vars%g(pp)

    do ii = 1_idef,params%num_species


       if (spp(ii)%spatial_distribution.eq.'TRACER'.and. &
            params%FO_GC_compare) then
          call get_fields(params,spp(ii)%vars,F)
          !! Calls [[get_fields]] in [[korc_fields]].
          ! Interpolates fields at local particles' position and keeps in
          ! spp%vars. Fields in (R,\(\phi\),Z) coordinates. 

          ALLOCATE(RAphi(spp(ii)%ppp,2))
          ALLOCATE(RVphi(spp(ii)%ppp))
          RAphi=0.0_rp
          
          call cart_to_cyl(spp(ii)%vars%X,spp(ii)%vars%Y)
          
          !$OMP PARALLEL DO SHARED(params,ii,spp,F,RAphi,RVphi) &
          !$OMP&  PRIVATE(pp,Bmag1,bhat,rm)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then

                RVphi(pp)=(-sin(spp(ii)%vars%Y(pp,2))*spp(ii)%vars%V(pp,1)+ &
                     cos(spp(ii)%vars%Y(pp,2))*spp(ii)%vars%V(pp,2))* &
                     spp(ii)%vars%Y(pp,1)

                Bmag1 = SQRT(spp(ii)%vars%B(pp,1)*spp(ii)%vars%B(pp,1)+ &
                     spp(ii)%vars%B(pp,2)*spp(ii)%vars%B(pp,2)+ &
                     spp(ii)%vars%B(pp,3)*spp(ii)%vars%B(pp,3))

                !             write(6,'("pp: ",I16)') pp
                !             write(6,'("Bmag: ",E17.10)') Bmag


                bhat = spp(ii)%vars%B(pp,:)/Bmag1

                if (params%field_model.eq.'ANALYTICAL') then
                   rm=sqrt((spp(ii)%vars%Y(pp,1)-F%AB%Ro)**2+ &
                        (spp(ii)%vars%Y(pp,3))**2)

                   RAphi(pp,1)=-F%AB%lambda**2*F%AB%Bo/(2*F%AB%qo)* &
                        log(1+(rm/F%AB%lambda)**2)

                else if (params%field_model.eq.'EXTERNAL') then

                   RAphi(pp,1)=spp(ii)%vars%PSI_P(pp)/(2*C_PI)
                   
                end if

                !             write(6,'("bhat: ",E17.10)') bhat 
                !             write(6,'("V: ",E17.10)') spp(ii)%vars%V(pp,:)


                spp(ii)%vars%X(pp,:)=spp(ii)%vars%X(pp,:)- &
                     spp(ii)%m*spp(ii)%vars%g(pp)* &
                     cross(bhat,spp(ii)%vars%V(pp,:))/(spp(ii)%q*Bmag1)

                ! transforming from particle location to associated
                ! GC location

             end if ! if particle in domain, i.e. spp%vars%flag==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

          call cart_to_cyl(spp(ii)%vars%X,spp(ii)%vars%Y)
          call get_fields(params,spp(ii)%vars,F)

          !$OMP PARALLEL DO SHARED(params,ii,spp,F,RAphi,RVphi) &
          !$OMP&  PRIVATE(pp,rm)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then

                if (params%field_model.eq.'ANALYTICAL') then
                   rm=sqrt((spp(ii)%vars%Y(pp,1)-F%AB%Ro)**2+ &
                        (spp(ii)%vars%Y(pp,3))**2)
                   RAphi(pp,2)=-F%AB%lambda**2*F%AB%Bo/(2*F%AB%qo)* &
                        log(1+(rm/F%AB%lambda)**2)

                else if (params%field_model.eq.'EXTERNAL') then

                   RAphi(pp,2)=spp(ii)%vars%PSI_P(pp)/(2*C_PI)
                   
                end if

                write(6,'("RAphi1: ",E17.10)') RAphi(pp,1)
                write(6,'("RAphi2: ",E17.10)') RAphi(pp,2)
                
                spp(ii)%vars%V(pp,1)=(spp(ii)%m*spp(ii)%vars%g(pp)* &
                     RVphi(pp)+spp(ii)%q*(RAphi(pp,1)-RAphi(pp,2)))/ &
                     spp(ii)%vars%Y(pp,1)
                !GC ppar              

             end if ! if particle in domain, i.e. spp%vars%flag==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

          !$OMP PARALLEL DO SHARED(ii,spp) PRIVATE(pp,Bmagc,bhatc)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then

                Bmagc = SQRT( DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                     spp(ii)%vars%B(pp,:)))

                bhatc = spp(ii)%vars%B(pp,:)/Bmagc

                spp(ii)%vars%V(pp,1)=spp(ii)%vars%V(pp,1)/ &
                     bhatc(2)
                !GC ppar

                spp(ii)%vars%V(pp,2)=spp(ii)%m/(2*Bmagc)* &
                     (spp(ii)%vars%g(pp)**2- &
                     (1+(spp(ii)%vars%V(pp,1)/spp(ii)%m)**2))           
                !GC mu


             end if ! if particle in domain, i.e. spp%vars%flag==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

          params%GC_coords=.TRUE.
          DEALLOCATE(RAphi)
          DEALLOCATE(RVphi)

          !Preparing Output Data
          call get_fields(params,spp(ii)%vars,F)

          !$OMP PARALLEL DO shared(F,params,spp) PRIVATE(pp,E_PHI)
          do pp=1_idef,spp(ii)%ppp,8

             !$OMP SIMD
             do cc=1_idef,8_idef
                E_PHI(cc)=spp(ii)%vars%E(pp-1+cc,2)
             end do
             !$OMP END SIMD
             
             call add_analytical_E_p(params,0_ip,F,E_PHI)


             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD
                
          end do
          !$OMP END PARALLEL DO


          !$OMP PARALLEL DO SHARED(ii,spp) PRIVATE(pp,Bmag1)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then

                Bmag1 = SQRT( DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                     spp(ii)%vars%B(pp,:)))

                spp(ii)%vars%g(pp)=sqrt(1+(spp(ii)%vars%V(pp,1))**2+ &
                     2*spp(ii)%vars%V(pp,2)*Bmag1)

!                write(6,'("Bmag:",E17.10)') Bmag1
!                write(6,'("PPLL:",E17.10)') spp(ii)%vars%V(pp,1)
!                write(6,'("MU:",E17.10)') spp(ii)%vars%V(pp,2)
                
                spp(ii)%vars%eta(pp) = atan2(sqrt(2*spp(ii)%m*Bmag1* &
                     spp(ii)%vars%V(pp,2)),spp(ii)%vars%V(pp,1))*180.0_rp/C_PI

!                             write(6,'("BR",E17.10)') spp(ii)%vars%B(pp,1)
!                             write(6,'("BPHI",E17.10)') spp(ii)%vars%B(pp,2)
!                             write(6,'("BZ",E17.10)') spp(ii)%vars%B(pp,3)

                !             write(6,'("ppll",E17.10)') spp(ii)%vars%V(pp,1)
                !             write(6,'("pperp",E17.10)') sqrt(2*spp(ii)%m*Bmag1* &
                !                  spp(ii)%vars%V(pp,2))

!                             write(6,'("eta GCinit",E17.10)') spp(ii)%vars%eta(pp)
                !             write(6,'("gam",E17.10)') spp(ii)%vars%g(pp)


             end if ! if particle in domain, i.e. spp%vars%flag==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO                
       else

          if (spp(ii)%spatial_distribution.eq.'TRACER') &
               call cart_to_cyl(spp(ii)%vars%X,spp(ii)%vars%Y)
          
          params%GC_coords=.TRUE.
          call get_fields(params,spp(ii)%vars,F)
          !$OMP PARALLEL DO shared(F,params,spp) PRIVATE(pp,E_PHI)
          do pp=1_idef,spp(ii)%ppp,8

             !$OMP SIMD
             do cc=1_idef,8_idef
                E_PHI(cc)=spp(ii)%vars%E(pp-1+cc,2)
             end do
             !$OMP END SIMD
             
             call add_analytical_E_p(params,0_ip,F,E_PHI)

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD
             
          end do
          !$OMP END PARALLEL DO

          
          !$OMP PARALLEL DO SHARED(ii,spp) PRIVATE(pp,Bmag1)

          do pp=1_idef,spp(ii)%ppp
!             if ( spp(ii)%vars%flag(pp) .EQ. 1_is ) then

!                write(6,'("BR: ",E17.10)') spp(ii)%vars%B(pp,1)
!                write(6,'("BPHI: ",E17.10)') spp(ii)%vars%B(pp,2)
!                write(6,'("BZ: ",E17.10)') spp(ii)%vars%B(pp,3)
                
                Bmag1 = SQRT( DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                     spp(ii)%vars%B(pp,:)))

                pmag=sqrt(spp(ii)%vars%g(pp)**2-1)
                
                spp(ii)%vars%V(pp,1)=pmag*cos(deg2rad(spp(ii)%vars%eta(pp)))

                spp(ii)%vars%V(pp,2)=(pmag* &
                     sin(deg2rad(spp(ii)%vars%eta(pp))))**2/ &
                     (2*spp(ii)%m*Bmag1)
                
                !    write(6,'("BR",E17.10)') spp(ii)%vars%B(pp,1)
                !    write(6,'("BPHI",E17.10)') spp(ii)%vars%B(pp,2)
                !    write(6,'("BZ",E17.10)') spp(ii)%vars%B(pp,3)

                !write(6,'("ppll",E17.10)') spp(ii)%vars%V(pp,1)
                !write(6,'("mu",E17.10)') spp(ii)%vars%V(pp,2)

                !     write(6,'("eta",E17.10)') spp(ii)%vars%eta(pp)
                !     write(6,'("gam",E17.10)') spp(ii)%vars%g(pp)


!             end if ! if particle in domain, i.e. spp%vars%flag==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO  
          
       end if

    end do ! loop over particle species
    
  end subroutine GC_init

  FUNCTION deg2rad(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: deg2rad

    deg2rad = C_PI*x/180.0_rp
  END FUNCTION deg2rad

  FUNCTION rad2deg(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: rad2deg

    rad2deg = x*180.0_rp/C_PI
  END FUNCTION rad2deg
  
  subroutine adv_GCeqn_top(params,F,P,spp)
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(8)               :: Bmag
    REAL(rp),DIMENSION(8) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8) :: B_R,B_PHI,B_Z,E_PHI
    REAL(rp),DIMENSION(8) :: V_PLL,V_MU
    REAL(rp) :: B0,EF0,R0,q0,lam,ar,m_cache,q_cache,ne0,Te0,Zeff0
    INTEGER(is),DIMENSION(8)  :: flag_cache

    LOGICAL                                                    :: ss_collisions
    !! Logical variable that indicates if collisions are included in
    !! the simulation.
    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
 

    do ii = 1_idef,params%num_species      

       q_cache=spp(ii)%q
       m_cache=spp(ii)%m



       
       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(E0,q_cache,m_cache) &
       !$OMP& shared(F,P,params,ii,spp) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,flag_cache, &
       !$OMP& B_R,B_PHI,B_Z,E_PHI)
       do pp=1_idef,spp(ii)%ppp,8

          !$OMP SIMD
          do cc=1_idef,8_idef
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             flag_cache(cc)=spp(ii)%vars%flag(pp-1+cc)
          end do
          !$OMP END SIMD
          
          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCeqn_vars(tt,params,Y_R,Y_PHI, &
                     Y_Z,V_PLL,V_MU,flag_cache,q_cache,m_cache, &
                     B_R,B_PHI,B_Z,F,P)
             end do !timestep iterator

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)
                
                spp(ii)%vars%flag(pp-1+cc)=flag_cache(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)
             end do
             !$OMP END SIMD
             
          else
             
             call advance_FPeqn_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,flag_cache,m_cache, &
                  F,P)

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flag(pp-1+cc)=flag_cache(cc)
             end do
             !$OMP END SIMD
             
          end if
                  
          
          call analytical_fields_Bmag_p(F,Y_R,Y_PHI,Y_Z, &
               Bmag,E_PHI)

          !$OMP SIMD
          do cc=1_idef,8_idef
             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc)*m_cache)

             spp(ii)%vars%eta(pp-1+cc) = rad2deg(atan2(sqrt(2*m_cache* &
                  Bmag(cc)*spp(ii)%vars%V(pp-1+cc,2)), &
                  spp(ii)%vars%V(pp-1+cc,1)))
          end do
          !$OMP END SIMD
             
       end do !particle chunk iterator
       !$OMP END PARALLEL DO
       
    end do !species iterator
    
  end subroutine adv_GCeqn_top

  subroutine advance_GCeqn_vars(tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,flag_cache, &
       q_cache,m_cache,B_R,B_PHI,B_Z,F,P)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.

    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(8) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(8) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(8) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(8) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(8) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(8) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(8) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(8) :: curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(8) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,V0,E_PHI,E_Z,E_R
    REAL(rp),DIMENSION(8) :: Bmag,ne,Te,Zeff
    INTEGER(is),dimension(8), intent(inout) :: flag_cache
    
    REAL(rp) :: ar,R0
    REAL(rp),intent(IN) :: q_cache,m_cache

    ar=F%AB%a
    R0=F%AB%Ro
    
    dt=params%dt

!    write(6,'("Y_R 0: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 0: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 0: ",E17.10)') Y_Z(1)
    
    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL)
    do cc=1_idef,8_idef
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0(cc)=V_PLL(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_GC_p(F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z)

!    write(6,'("ER:",E17.10)') E_R
!    write(6,'("EPHI:",E17.10)') E_PHI
!    write(6,'("EZ:",E17.10)') E_Z
    

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k1_R,k1_PHI,k1_Z,k1_PLL)
    do cc=1_idef,8
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0(cc)   +a1*k1_PLL(cc)
    end do
    !$OMP END SIMD


!    write(6,'("Y_R 1: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 1: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 1: ",E17.10)') Y_Z(1)
    
!    write(6,'("k1R: ",E17.10)') k1_R(1)
!    write(6,'("k1PHI: ",E17.10)') k1_PHI(1)
!    write(6,'("k1Z: ",E17.10)') k1_Z(1)
!    write(6,'("k1PLL: ",E17.10)') k1_PLL(1) 
    
    call analytical_fields_GC_p(F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z)


    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k2_R,k2_PHI,k2_Z,k2_PLL)
    do cc=1_idef,8
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
    end do
    !$OMP END SIMD


!    write(6,'("Y_R 2: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 2: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 2: ",E17.10)') Y_Z(1)
    
    call analytical_fields_GC_p(F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z)


    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k3_R,k3_PHI,k3_Z,k3_PLL)
    do cc=1_idef,8
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
    end do
    !$OMP END SIMD

!    write(6,'("Y_R 3: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 3: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 3: ",E17.10)') Y_Z(1)
    
    call analytical_fields_GC_p(F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z)


    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k4_R,k4_PHI,k4_Z,k4_PLL)
    do cc=1_idef,8
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
    end do
    !$OMP END SIMD

!    write(6,'("Y_R 4: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 4: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 4: ",E17.10)') Y_Z(1)
    
    call analytical_fields_GC_p(F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z)


    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k5_R,k5_PHI,k5_Z,k5_PLL)
    do cc=1_idef,8
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
    end do
    !$OMP END SIMD

!    write(6,'("Y_R 5: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 5: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 5: ",E17.10)') Y_Z(1)
    
    call analytical_fields_GC_p(F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z)


    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k6_R,k6_PHI,k6_Z,k6_PLL)
    do cc=1_idef,8
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc) 
    end do
    !$OMP END SIMD

!    write(6,'("Y_R 6: ",E17.10)') Y_R(1)
!    write(6,'("Y_PHI 6: ",E17.10)') Y_PHI(1)
!    write(6,'("Y_Z 6: ",E17.10)') Y_Z(1)
    
    call cyl_check_if_confined_p(ar,R0,Y_R,Y_Z,flag_cache)

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,Y0_R,Y0_PHI,Y0_Z,V0)
    do cc=1_idef,8

       if (flag_cache(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0(cc)
       end if          
 
    end do
    !$OMP END SIMD
    
    if (params%collisions) then
       
       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flag_cache,F,P,E_PHI)

    end if


  end subroutine advance_GCeqn_vars

  subroutine advance_FPeqn_vars(params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,flag_cache, &
       m_cache,F,P)

    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(8), INTENT(INOUT)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8), INTENT(INOUT)  :: V_PLL,V_MU
    REAL(rp),DIMENSION(8)  :: E_PHI
    INTEGER(is),DIMENSION(8), INTENT(INOUT)  :: flag_cache
    REAL(rp),intent(in) :: m_cache

    do tt=1_ip,params%t_skip
       
       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flag_cache,F,P,E_PHI)

!       write(6,'("Collision Loop in FP")')
       
    end do


    

  end subroutine advance_FPeqn_vars


  subroutine adv_GCinterp_top(params,spp,P,F)
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(8)               :: Bmag
    REAL(rp),DIMENSION(8) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(8) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(8) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(8) :: V_PLL,V_MU
    INTEGER(is),DIMENSION(8) :: flag_cache
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar

    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
 

    do ii = 1_idef,params%num_species      

       q_cache=spp(ii)%q
       m_cache=spp(ii)%m
       
       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flag_cache,E_PHI)
       do pp=1_idef,spp(ii)%ppp,8

!          write(6,'("pp: ",I16)') pp
          
          !$OMP SIMD
          do cc=1_idef,8_idef
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             flag_cache(cc)=spp(ii)%vars%flag(pp-1+cc)           
          end do
          !$OMP END SIMD
          
          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCinterp_vars(tt,params,Y_R,Y_PHI, &
                     Y_Z,V_PLL,V_MU,q_cache,m_cache,flag_cache, &
                     F,P,B_R,B_PHI,B_Z,E_PHI)
             end do !timestep iterator

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flag(pp-1+cc)=flag_cache(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD
             
          else
             
             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flag_cache,F,P,E_PHI)

             !$OMP SIMD
             do cc=1_idef,8_idef
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flag(pp-1+cc)=flag_cache(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD
             
          end if                            
          

          !$OMP SIMD
          do cc=1_idef,8_idef
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)
             
             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))
             
             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD
             
       end do !particle chunk iterator
       !$OMP END PARALLEL DO
       
    end do !species iterator
    
  end subroutine adv_GCinterp_top
  
  subroutine advance_GCinterp_vars(tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
       q_cache,m_cache,flag_cache,F,P,B_R,B_PHI,B_Z,E_PHI)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    

    REAL(rp),DIMENSION(8)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(8) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(8) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(8) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(8) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(8) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(8) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(8) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(8) :: E_R,E_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(8) :: curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(8) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,V0
    REAL(rp),DIMENSION(8) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(8),intent(INOUT) :: flag_cache
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL)
    do cc=1_idef,8_idef

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0(cc)=V_PLL(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

    call add_analytical_E_p(params,tt,F,E_PHI)

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k1_R,k1_PHI,k1_Z,k1_PLL)
    do cc=1_idef,8
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)    
       
       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0(cc)   +a1*k1_PLL(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

    call add_analytical_E_p(params,tt,F,E_PHI)

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k2_R,k2_PHI,k2_Z,k2_PLL)
    do cc=1_idef,8
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

    call add_analytical_E_p(params,tt,F,E_PHI)

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k3_R,k3_PHI,k3_Z,k3_PLL)
    do cc=1_idef,8
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

    call add_analytical_E_p(params,tt,F,E_PHI)
    
    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k4_R,k4_PHI,k4_Z,k4_PLL)
    do cc=1_idef,8
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
    end do
    !$OMP END SIMD


    call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

    call add_analytical_E_p(params,tt,F,E_PHI)
    
    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k5_R,k5_PHI,k5_Z,k5_PLL)
    do cc=1_idef,8
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

    call add_analytical_E_p(params,tt,F,E_PHI)
    
    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k6_R,k6_PHI,k6_Z,k6_PLL)
    do cc=1_idef,8
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)


       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,Y0_R,Y0_PHI,Y0_Z,V0)
    do cc=1_idef,8

       if (flag_cache(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0(cc)
       end if          
 
    end do
    !$OMP END SIMD
    
    if (params%collisions) then       
       
       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flag_cache,F,P,E_PHI)

    end if


  end subroutine advance_GCinterp_vars

  subroutine advance_FPinterp_vars(params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
       m_cache,flag_cache,F,P,E_PHI)    
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(8), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8), INTENT(INOUT)  :: V_PLL,V_MU
    REAL(rp),DIMENSION(8), INTENT(OUT)  :: E_PHI
    REAL(rp),intent(in) :: m_cache
    INTEGER(is),DIMENSION(8),intent(INOUT) :: flag_cache

!    write(6,'("E_PHI_FP: ",E17.10)') E_PHI
    
    do tt=1_ip,params%t_skip
       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flag_cache,F,P,E_PHI)

!       write(6,'("Collision Loop in FP")')
       
    end do

!    write(6,'("V_PLL: ",E17.10)') V_PLL
!    write(6,'("V_MU: ",E17.10)') V_MU
    
  end subroutine advance_FPinterp_vars

  

  subroutine GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
       B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI, &
       gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    REAL(rp),DIMENSION(8)  :: Bmag,bhat_R,bhat_PHI,bhat_Z,Bst_R,Bst_PHI
    REAL(rp),DIMENSION(8)  :: BstdotE,BstdotgradB,EcrossB_R,EcrossB_PHI,bdotBst
    REAL(rp),DIMENSION(8)  :: bcrossgradB_R,bcrossgradB_PHI,bcrossgradB_Z,gamgc
    REAL(rp),DIMENSION(8)  :: EcrossB_Z,Bst_Z
    REAL(rp),DIMENSION(8)  :: pm,xi
    REAL(rp),DIMENSION(8),INTENT(in) :: gradB_R,gradB_PHI,gradB_Z,curlb_R
    REAL(rp),DIMENSION(8),INTENT(in) :: curlb_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL
    REAL(rp),DIMENSION(8),INTENT(IN) :: V_PLL,V_MU,Y_R,curlb_PHI
    REAL(rp),INTENT(in) :: q_cache,m_cache
    INTEGER(ip)  :: cc

    !$OMP SIMD
!    !$OMP& aligned(gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_Z, &
!    !$OMP& B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,RHS_R,RHS_PHI,RHS_Z,RHS_PLL, &
!    !$OMP& V_PLL,V_MU,Y_R,curlb_PHI)
    do cc=1_idef,8
       Bmag(cc) = SQRT(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       bhat_R(cc) = B_R(cc)/Bmag(cc)
       bhat_PHI(cc) = B_PHI(cc)/Bmag(cc)
       bhat_Z(cc) = B_Z(cc)/Bmag(cc)

       Bst_R(cc)=q_cache*B_R(cc)+V_PLL(cc)*curlb_R(cc)
       Bst_PHI(cc)=q_cache*B_PHI(cc)+V_PLL(cc)*curlb_PHI(cc)
       Bst_Z(cc)=q_cache*B_Z(cc)+V_PLL(cc)*curlb_Z(cc)

       bdotBst(cc)=bhat_R(cc)*Bst_R(cc)+bhat_PHI(cc)*Bst_PHI(cc)+ &
            bhat_Z(cc)*Bst_Z(cc)
       BstdotE(cc)=Bst_R(cc)*E_R(cc)+Bst_PHI(cc)*E_PHI(cc)+Bst_Z(cc)*E_Z(cc)   
       BstdotgradB(cc)=Bst_R(cc)*gradB_R(cc)+Bst_PHI(cc)*gradB_PHI(cc)+ &
            Bst_Z(cc)*gradB_Z(cc)

       Ecrossb_R(cc)=E_PHI(cc)*bhat_Z(cc)-E_Z(cc)*bhat_PHI(cc)
       Ecrossb_PHI(cc)=E_Z(cc)*bhat_R(cc)-E_R(cc)*bhat_Z(cc)
       Ecrossb_Z(cc)=E_R(cc)*bhat_PHI(cc)-E_PHI(cc)*bhat_R(cc)


       bcrossgradB_R(cc)=bhat_PHI(cc)*gradB_Z(cc)-bhat_Z(cc)*gradB_PHI(cc)
       bcrossgradB_PHI(cc)=bhat_Z(cc)*gradB_R(cc)-bhat_R(cc)*gradB_Z(cc)
       bcrossgradB_Z(cc)=bhat_R(cc)*gradB_PHI(cc)-bhat_PHI(cc)*gradB_R(cc)

       gamgc(cc)=sqrt(1+V_PLL(cc)*V_PLL(cc)+2*V_MU(cc)*Bmag(cc))

       pm(cc)=sqrt(gamgc(cc)**2-1)
       xi(cc)=V_PLL(cc)/pm(cc)
       
       RHS_R(cc)=(q_cache*Ecrossb_R(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_R(cc)+V_PLL(cc)*Bst_R(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PHI(cc)=(q_cache*Ecrossb_PHI(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_PHI(cc)+V_PLL(cc)*Bst_PHI(cc))/(m_cache*gamgc(cc)))/ &
            (Y_R(cc)*bdotBst(cc))
       RHS_Z(cc)=(q_cache*Ecrossb_Z(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_Z(cc)+V_PLL(cc)*Bst_Z(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PLL(cc)=(q_cache*BstdotE(cc)-V_MU(cc)*BstdotgradB(cc)/gamgc(cc))/ &
            bdotBst(cc)
       
    end do
    !$OMP END SIMD

  end subroutine GCEoM_p


  subroutine aux_fields(pp,spp,gradB,curlb,Bmag)
    TYPE(SPECIES), INTENT(IN)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(3),INTENT(INOUT) :: gradB
    REAL(rp),DIMENSION(3),INTENT(INOUT) :: curlb
    REAL(rp),INTENT(IN) :: Bmag
    REAL(rp) :: dRB
    REAL(rp) :: dPHIB
    REAL(rp) :: dZB
    INTEGER  :: pp

    dRB=(spp%vars%B(pp,1)*spp%vars%BR(pp,1)+ &
         spp%vars%B(pp,2)*spp%vars%BPHI(pp,1)+ &
         spp%vars%B(pp,3)*spp%vars%BZ(pp,1))/Bmag
    dPHIB=(spp%vars%B(pp,1)*spp%vars%BR(pp,2)+ &
         spp%vars%B(pp,2)*spp%vars%BPHI(pp,2)+ &
         spp%vars%B(pp,3)*spp%vars%BZ(pp,2))/Bmag
    dZB=(spp%vars%B(pp,1)*spp%vars%BR(pp,3)+ &
         spp%vars%B(pp,2)*spp%vars%BPHI(pp,3)+ &
         spp%vars%B(pp,3)*spp%vars%BZ(pp,3))/Bmag

    gradB(1)=dRB
    gradB(2)=dPHIB/spp%vars%Y(pp,1)
    gradB(3)=dZB

    curlb(1)=((Bmag*spp%vars%BZ(pp,2)-spp%vars%B(pp,3)*dPHIB)/spp%vars%Y(pp,1)- &
         (Bmag*spp%vars%BPHI(pp,3)-spp%vars%B(pp,2)*dZB))/Bmag**2
    curlb(2)=((Bmag*spp%vars%BR(pp,3)-spp%vars%B(pp,1)*dZB)- &
         (Bmag*spp%vars%BZ(pp,1)-spp%vars%B(pp,3)*dRB))/Bmag**2
    curlb(3)=((Bmag*spp%vars%BPHI(pp,1)-spp%vars%B(pp,2)*dRB) - &
         (Bmag*spp%vars%BPHI(pp,1)-spp%vars%B(pp,1)*dPHIB)/ &
         spp%vars%Y(pp,1))/Bmag**2+ &
         spp%vars%B(pp,2)/(Bmag*spp%vars%Y(pp,1))  

  end subroutine aux_fields


end module korc_ppusher
