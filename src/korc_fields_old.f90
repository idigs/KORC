module korc_fields
  !! @note Module containing subroutines to initialize externally
  !! generated fields, and to calculate the electric and magnetic
  !! fields when using an analytical model. @endnote
  use korc_types
  use korc_hpc
  use korc_coords
  use korc_interp
  use korc_HDF5

  IMPLICIT NONE

  PUBLIC :: mean_F_field,&
       get_fields,&
       initialize_fields,&
       load_field_data_from_hdf5,&
       load_dim_data_from_hdf5,&
       ALLOCATE_2D_FIELDS_ARRAYS,&
       ALLOCATE_3D_FIELDS_ARRAYS,&
       DEALLOCATE_FIELDS_ARRAYS
  PRIVATE :: get_analytical_fields,&
       analytical_fields,&
       analytical_fields_GC,&
       uniform_magnetic_field,&
       uniform_electric_field,&
       uniform_fields,&
       cross,&
       analytical_electric_field_cyl,&
       ALLOCATE_V_FIELD_2D,&
       ALLOCATE_V_FIELD_3D,&
       initialize_GC_fields

CONTAINS


  subroutine analytical_fields(F,Y,E,B,flag)
    !! @note Subroutine that calculates and returns the analytic electric and
    !! magnetic field for each particle in the simulation. @endnote
    !! The analytical magnetic field is given by:
    !!
    !! $$\mathbf{B}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}}
    !! \left[ B_0 \hat{e}_\zeta + B_\vartheta(r) \hat{e}_\vartheta \right],$$
    !!
    !! where \(\eta = r/R_0\) is the aspect ratio, the constant \(B_0\)
    !! denotes the magnitude of the toroidal magnetic field,
    !! and \(B_\vartheta(r) = \eta B_0/q(r)\) is the poloidal magnetic
    !! field with 
    !! safety factor \(q(r) = q_0\left( 1 + \frac{r^2}{\lambda^2} \right)\).
    !! The constant \(q_0\) is the safety factor at the magnetic axis and
    !! the constant \(\lambda\) is obtained from the values of \(q_0\)
    !! and \(q(r)\) at the plasma edge \(r=r_{edge}\). On the other hand,
    !! the analytical electric fields is given by:
    !!
    !! $$\mathbf{E}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}}
    !! E_0 \hat{e}_\zeta,$$
    !!
    !! where \(E_0\) is the electric field as measured at the mangetic axis.
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Toroidal coordinates of each particle in the simulation; 
    !! Y(1,:) = \(r\), Y(2,:) = \(\theta\), Y(3,:) = \(\zeta\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in Cartesian coordinates; 
    !! B(1,:) = \(B_x\), B(2,:) = \(B_y\), B(3,:) = \(B_z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in Cartesian coordinates; 
    !! E(1,:) = \(E_x\), E(2,:) = \(E_y\), E(3,:) = \(E_z\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp)                                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp)                                               :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp)                                               :: eta
    !! Aspect ratio \(\eta\).
    REAL(rp)                                               :: q
    !! Safety profile \(q(r)\).
    INTEGER(ip)                                            :: pp ! Iterator(s)
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.

    ss = SIZE(Y,2)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ezeta,Bp,Bzeta,eta,q) &
    !$OMP& SHARED(F,Y,E,B,flag)
    do pp=1_idef,ss
       if ( flag(pp) .EQ. 1_is ) then
          eta = Y(1,pp)/F%Ro
          q = F%AB%qo*(1.0_rp + (Y(1,pp)/F%AB%lambda)**2)
          Bp = F%AB%Bp_sign*eta*F%AB%Bo/(q*(1.0_rp + eta*COS(Y(2,pp))))
          Bzeta = F%AB%Bo/( 1.0_rp + eta*COS(Y(2,pp)) )


          B(1,pp) =  Bzeta*COS(Y(3,pp)) - Bp*SIN(Y(2,pp))*SIN(Y(3,pp))
          B(2,pp) = -Bzeta*SIN(Y(3,pp)) - Bp*SIN(Y(2,pp))*COS(Y(3,pp))
          B(3,pp) = Bp*COS(Y(2,pp))

          if (abs(F%Eo) > 0) then
             Ezeta = F%Eo/( 1.0_rp + eta*COS(Y(2,pp)) )

             E(1,pp) = Ezeta*COS(Y(3,pp))
             E(2,pp) = -Ezeta*SIN(Y(3,pp))
             E(3,pp) = 0.0_rp
          end if
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine analytical_fields


  subroutine analytical_fields_GC(params,F,Y,E,B,gradB,curlb,flag)
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Cylindrical coordinates of each particle in the simulation; 
    !! Y(1,:) = \(r\), Y(2,:) = \(\theta\), Y(3,:) = \(\zeta\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in cylindrical coordinates; 
    !! B(1,:) = \(B_R\), B(2,:) = \(B_\phi\), B(3,:) = \(B_Z\)
    REAL(rp), DIMENSION(3)   :: Btmp
    !! Placeholder for magnetic field components in cylindrical coordinates; 
    !! B(1,:) = \(B_R\), B(2,:) = \(B_\phi\), B(3,:) = \(B_Z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: gradB
    !! Gradient of magnitude of magnetic field in cylindrical coordinates; 
    !! gradB(1,:) = \(\nabla_R B\), B(2,:) = \(\nabla_\phi B_\),
    !! B(3,:) = \(\nabla_Z B\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: curlB
    !! Curl of magnetic field unit vector in cylindrical coordinates 
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in cylindricalcoordinates; 
    !! E(1,:) = \(E_R\), E(2,:) = \(E_\phi\), E(3,:) = \(E_Z\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp)                                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp)                                               :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp)                                               :: eta
    !! Aspect ratio \(\eta\).
    REAL(rp)                                               :: q
    !! Safety profile \(q(r)\).
    INTEGER(ip)                                            :: pp ! Iterator(s)
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.
    REAL(rp)    ::  dRBR
    REAL(rp)    ::  dRBPHI
    REAL(rp)    ::  dRBZ
    REAL(rp)    ::  dZBR
    REAL(rp)    ::  dZBPHI
    REAL(rp)    ::  dZBZ
    REAL(rp)    ::  Bmag
    REAL(rp)    ::  dRbhatPHI
    REAL(rp)    ::  dRbhatZ
    REAL(rp)    ::  dZbhatR
    REAL(rp)    ::  dZbhatPHI
    REAL(rp)    ::  qprof
    REAL(rp)    ::  rm

    
    ss = SIZE(Y,2)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag)
    do pp=1_idef,ss
       if ( flag(pp) .EQ. 1_is ) then

          rm=sqrt((Y(1,pp)-F%AB%Ro)**2+Y(3,pp)**2)
          qprof = 1.0_rp + (rm/F%AB%lambda)**2
          
          Btmp(1)=F%AB%Bo*Y(3,pp)/(F%AB%qo*qprof*Y(1,pp))
          Btmp(2)=-F%AB%Bo*F%AB%Ro/Y(1,pp)
          Btmp(3)=-F%AB%Bo*(Y(1,pp)-F%AB%Ro)/(F%AB%qo*qprof*Y(1,pp))

          if (.not.params%GC_coords) then
             B(1,pp) = Btmp(1)*COS(Y(2,pp)) - Btmp(2)*SIN(Y(2,pp))
             B(2,pp) = Btmp(1)*SIN(Y(2,pp)) + Btmp(2)*COS(Y(2,pp))
             B(3,pp) = Btmp(3)
          else
             B(1,pp) = Btmp(1)
             B(2,pp) = Btmp(2)
             B(3,pp) = Btmp(3)
          end if

          dRBR=-F%AB%Bo*Y(3,pp)/(F%AB%qo*qprof*Y(1,pp))*(1./Y(1,pp)+ &
               2*(Y(1,pp)-F%AB%Ro)/(F%AB%lambda**2*qprof))
          dRBPHI=F%AB%Bo*F%AB%Ro/Y(1,pp)**2
          dRBZ=F%AB%Bo/(F%AB%qo*qprof*Y(1,pp))*(-F%AB%Ro/Y(1,pp)+2*(Y(1,pp)- &
               F%AB%Ro)**2/(F%AB%lambda**2*qprof))
          dZBR=F%AB%Bo/(F%AB%qo*qprof*Y(1,pp))*(1-2*Y(3,pp)*Y(3,pp)/ &
               (F%AB%lambda**2*qprof))
          dZBPHI=0._rp
          dZBZ=F%AB%Bo*(Y(1,pp)-F%AB%Ro)/(F%AB%qo*Y(1,pp))*2*Y(3,pp)/ &
               ((F%AB%lambda*qprof)**2)

          Bmag=sqrt(dot_product(B(:,pp),B(:,pp)))

          gradB(1,pp)=(B(1,pp)*dRBR+B(2,pp)*dRBPHI+B(3,pp)*dRBZ)/Bmag
          gradB(2,pp)=0._rp
          gradB(3,pp)=(B(1,pp)*dZBR+B(2,pp)*dZBPHI+B(3,pp)*dZBZ)/Bmag
          
          dRbhatPHI=(Bmag*dRBPHI-B(2,pp)*gradB(1,pp))/Bmag**2
          dRbhatZ=(Bmag*dRBZ-B(3,pp)*gradB(1,pp))/Bmag**2
          dZbhatR=(Bmag*dZBR-B(1,pp)*gradB(3,pp))/Bmag**2
          dZbhatPHI=(Bmag*dZBPHI-B(2,pp)*gradB(3,pp))/Bmag**2

          curlb(1,pp)=-dZbhatPHI
          curlb(2,pp)=dZbhatR-dRbhatZ
          curlb(3,pp)=B(2,pp)/(Bmag*Y(1,pp))+dRbhatPHI          

          if (abs(F%Eo) > 0) then
             E(1,pp) = 0.0_rp
             E(2,pp) = -F%Eo*F%AB%Ro/Y(1,pp)
             E(3,pp) = 0.0_rp
          end if
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine analytical_fields_GC

  subroutine uniform_magnetic_field(F,B)
    !! @note Subroutine that returns the value of a uniform magnetic
    !! field. @endnote
    !! This subroutine is used only when the simulation is ran for a
    !! 'UNIFORM' plasma. As a convention, in a uniform plasma we
    !! set \(\mathbf{B} = B_0 \hat{x}\).
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in Cartesian coordinates; 
    !! B(1,:) = \(B_x\), B(2,:) = \(B_y\), B(3,:) = \(B_z\)
    B(1,:) = F%Bo
    B(2:3,:) = 0.0_rp
  end subroutine uniform_magnetic_field


  subroutine uniform_electric_field(F,E)
    !! @note Subroutine that returns the value of a uniform electric
    !! field. @endnote
    !! This subroutie is used only when the simulation is ran for a
    !! 'UNIFORM' plasma. As a convention, in a uniform plasma we set
    !! \(\mathbf{E} = E_0 \hat{x}\).
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in Cartesian coordinates; 
    !! E(1,:) = \(E_x\), E(2,:) = \(E_y\), E(3,:) = \(E_z\)

    E(1,:) = F%Eo
    E(2:3,:) = 0.0_rp
  end subroutine uniform_electric_field


  subroutine analytical_electric_field_cyl(F,Y,E,flag)
    !! @note Subrotuine that calculates and returns the electric field using the
    !! same analytical model of the 'analytical_fields' subroutine. @endnote
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Cylindrical coordinates of each particle in the simulation;
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), Y(3,:) = \(Z\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in Cartesian coordinates;
    !!  E(1,:) = \(E_x\), E(2,:) = \(E_y\), E(3,:) = \(E_z\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ephi
    !! Azimuthal electric field.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.

    if (abs(F%Eo) > 0) then
       ss = SIZE(Y,2)
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ephi) SHARED(F,Y,E,flag)
       do pp=1_idef,ss
          if ( flag(pp) .EQ. 1_is ) then
             Ephi = F%Eo*F%Ro/Y(1,pp)

             E(1,pp) = -Ephi*SIN(Y(2,pp))
             E(2,pp) = Ephi*COS(Y(2,pp))
             E(3,pp) = 0.0_rp
          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine analytical_electric_field_cyl


  subroutine mean_F_field(F,Fo,op_field)
    !! @note Subroutine that calculates the mean electric or magnetic field in
    !! case external fields are being used. @endnote
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), INTENT(OUT)          :: Fo
    !! Mean electric or magnetic field.
    TYPE(KORC_STRING), INTENT(IN)  :: op_field
    !!String that specifies what mean field will be calculated.
    !! Its value can be 'B' or 'E'.

    if (TRIM(op_field%str) .EQ. 'B') then
       if (ALLOCATED(F%B_3D%R)) then ! 3D field
          Fo = SUM( SQRT(F%B_3D%R**2 + F%B_3D%PHI**2 + F%B_3D%Z**2) )/SIZE(F%B_3D%R)
       else if (ALLOCATED(F%B_2D%R)) then ! Axisymmetric 2D field
          Fo = SUM( SQRT(F%B_2D%R**2 + F%B_2D%PHI**2 + F%B_2D%Z**2) )/SIZE(F%B_2D%R)
       end if
    else if (TRIM(op_field%str) .EQ. 'E') then
       if (ALLOCATED(F%E_3D%R)) then ! 3D field
          Fo = SUM( SQRT(F%E_3D%R**2 + F%E_3D%PHI**2 + F%E_3D%Z**2) )/SIZE(F%E_3D%R)
       else if (ALLOCATED(F%E_2D%R)) then ! Axisymmetric 2D field
          Fo = SUM( SQRT(F%E_2D%R**2 + F%E_2D%PHI**2 + F%E_2D%Z**2) )/SIZE(F%E_2D%R)
       end if
    else
       write(6,'("KORC ERROR: Please enter a valid field: mean_F_field")')
       call korc_abort()
    end if
  end subroutine mean_F_field
  

  subroutine get_analytical_fields(params,vars,F)
    !! @note Interface for calculating the analytical electric and magnetic
    !! fields for each particle in the simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT) :: vars
    !! An instance of the KORC derived type PARTICLES.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of the KORC derived type FIELDS.

    if (params%orbit_model(1:2).eq.'FO') then

       call cart_to_tor_check_if_confined(vars%X,F,vars%Y,vars%flag)
       
       call analytical_fields(F,vars%Y, vars%E, vars%B, vars%flag)

    elseif (params%orbit_model(1:2).eq.'GC') then

       if (.not.params%GC_coords) then

          call cart_to_cyl(vars%X,vars%Y)

       endif

       call cyl_check_if_confined(F,vars%Y,vars%flag)
       
       call analytical_fields_GC(params,F,vars%Y, vars%E, vars%B, vars%gradB, &
            vars%curlb, vars%flag)

    endif
    
  end subroutine get_analytical_fields

  subroutine uniform_fields(vars,F)
    !! @note Interface for calculating the uniform electric and magnetic
    !! fields for each particle in the simulation. @endnote
    TYPE(PARTICLES), INTENT(INOUT) :: vars
    !! An instance of the KORC derived type PARTICLES.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of the KORC derived type FIELDS.
    
    call uniform_magnetic_field(F, vars%B)

    call uniform_electric_field(F, vars%E)
  end subroutine uniform_fields


  function cross(a,b)
    !! @note Function that calculates the cross product of the two
    !! vectors \(\mathbf{a}\) and \(\mathbf{b}\). @endnote
    REAL(rp), DIMENSION(3)             :: cross
    !! Cross product \(\mathbf{a}\times \mathbf{b}\)
    REAL(rp), DIMENSION(3), INTENT(IN) :: a
    !!  Vector \(\mathbf{a}\).
    REAL(rp), DIMENSION(3), INTENT(IN) :: b
    !!  Vector \(\mathbf{b}\).

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function cross


  subroutine unitVectors(params,Xo,F,b1,b2,b3,flag)
    !! @note Subrotuine that calculates an orthonormal basis using information 
    !! of the (local) magnetic field at position \(\mathbf{X}_0\). @endnote
    TYPE(KORC_PARAMS), INTENT(IN)                                      :: params
    !! Core KORC simulation parameters.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)                  :: Xo
    !! Array with the position of the simulated particles.
    TYPE(FIELDS), INTENT(IN)                                           :: F
    !! F An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b1
    !! Basis vector pointing along the local magnetic field, 
    !! that is, along \(\mathbf{b} = \mathbf{B}/B\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b2
    !!  Basis vector perpendicular to b1
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b3
    !! Basis vector perpendicular to b1 and b2.
    INTEGER(is), DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT)    :: flag
    !! Flag for each particle to decide whether it is being 
    !! followed (flag=T) or not (flag=F).
    TYPE(PARTICLES)                                                    :: vars
    !! A temporary instance of the KORC derived type PARTICLES.
    INTEGER                                                            :: ii
    !! Iterator.
    INTEGER                                                            :: ppp
    !! Number of particles.

    ppp = SIZE(Xo,2) ! Number of particles

    ALLOCATE( vars%X(3,ppp) )
    ALLOCATE( vars%Y(3,ppp) )
    ALLOCATE( vars%B(3,ppp) )
    ALLOCATE( vars%gradB(3,ppp) )
    ALLOCATE( vars%curlb(3,ppp) )
    ALLOCATE( vars%E(3,ppp) )
    ALLOCATE( vars%flag(ppp) )

    vars%X = Xo
    vars%flag = 1_idef

    call init_random_seed()

    call get_fields(params,vars,F)

    do ii=1_idef,ppp
       if ( vars%flag(ii) .EQ. 1_idef ) then
          b1(:,ii) = vars%B(:,ii)/SQRT( DOT_PRODUCT(vars%B(:,ii),vars%B(:,ii)) )

          b2(:,ii) = cross(b1(:,ii),(/0.0_rp,0.0_rp,1.0_rp/))
          b2(:,ii) = b2(:,ii)/SQRT( DOT_PRODUCT(b2(:,ii),b2(:,ii)) )

          b3(:,ii) = cross(b1(:,ii),b2(:,ii))
          b3(:,ii) = b3(:,ii)/SQRT( DOT_PRODUCT(b3(:,ii),b3(:,ii)) )
       end if
    end do

    if (PRESENT(flag)) then
       flag = vars%flag
    end if

    DEALLOCATE( vars%X )
    DEALLOCATE( vars%Y )
    DEALLOCATE( vars%B )
    DEALLOCATE( vars%gradB )
    DEALLOCATE( vars%curlb )
    DEALLOCATE( vars%E )
    DEALLOCATE( vars%flag )
  end subroutine unitVectors


  subroutine get_fields(params,vars,F)
    !! @note Inferface with calls to subroutines for calculating the electric 
    !! and magnetic field for each particle in the simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    !!  An instance of the KORC derived type PARTICLES.
    TYPE(FIELDS), INTENT(IN)           :: F
    !! An instance of the KORC derived type FIELDS.

    SELECT CASE (TRIM(params%plasma_model))
    CASE('ANALYTICAL')
       if (params%field_eval.eq.'eqn') then
          call get_analytical_fields(params,vars, F)
       else
          call interp_fields(params,vars, F)          
       end if       
    CASE('EXTERNAL')

       call interp_fields(params,vars, F)
       if (F%Efield.AND..NOT.F%Efield_in_file) then
          call analytical_electric_field_cyl(F,vars%Y,vars%E,vars%flag)
       end if
    CASE('UNIFORM')
       
       call uniform_fields(vars, F)
    CASE DEFAULT
    END SELECT
  end subroutine get_fields


  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * * *  SUBROUTINES FOR INITIALIZING FIELDS   * * * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine initialize_fields(params,F)
    !! @note Subroutine that initializes the analytical or externally
    !! calculated electric and magnetic fields. @endnote
    !! In this subroutine we load the parameters of the electric and
    !! magnetic fields from the namelists 'analytical_fields_params' and
    !! 'externalPlasmaModel' in the input file.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(OUT)      :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp)                       :: Bo
    !! Magnetic field at magnetic axis for an 'ANALYTICAL' magnetic field,
    !! or the magnitude of the magnetic field for a 'UNFIROM' plasma.
    REAL(rp)                       :: minor_radius
    !! Plasma edge \(r_{edge}\) as measured from the magnetic axis.
    REAL(rp)                       :: major_radius
    !! Radial position of the magnetic axis \(R_0\)
    REAL(rp)                       :: qa
    !! Safety factor at the plasma edge.
    REAL(rp)                       :: qo
    !! Safety factor at the magnetic axis \(q_0\).
    CHARACTER(MAX_STRING_LENGTH)   :: current_direction
    !! String with information about the direction of the plasma current, 
    !! 'PARALLEL'  or 'ANTI-PARALLEL' to the toroidal magnetic field.
    REAL(rp)                       :: Eo
    !! Electric field at the magnetic axis.
    LOGICAL                        :: Efield
    !! Logical variable that specifies if the electric field is 
    !! going to be used on in a given simulation.
    LOGICAL                        :: Bfield
    !! Logical variable that specifies if the magnetic field is 
    !! going to be used on in a given simulation.
    LOGICAL                        :: Bflux
    !! Logical variable that specifies if the poloidal magnetic 
    !! flux is going to be used on in a given simulation.
    LOGICAL                        :: axisymmetric_fields
    !! Logical variable that specifies if the plasma is axisymmetric.
    INTEGER                        :: ii
    !! Iterators for creating mesh for GC model with analytic fields
    INTEGER                        :: kk
    !! Iterators for creating mesh for GC model with analytic fields
    INTEGER                        :: nR
    !! Number of mesh points in R for grid in GC model of analytical field
    INTEGER                        :: nZ
    !! Number of mesh points in Z for grid in GC model of analytical field
    real(rp)                       :: rm
    !! Minor radius at each position in the grid for
    !! GC model of analytical field
    real(rp)                       :: qr
    !! Safety factor at each position in the grid for
    !! GC model of analytical field
    real(rp)                       :: theta
    !! Poloidal angle at each position in the grid for
    !! GC model of analytical field
    
    NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
         qa,qo,Eo,current_direction,nR,nZ

    NAMELIST /externalPlasmaModel/ Efield, Bfield, Bflux, axisymmetric_fields
    
    SELECT CASE (TRIM(params%plasma_model))
       
    CASE('ANALYTICAL')
       ! Load the parameters of the analytical magnetic field
       open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
            status='OLD',form='formatted')
       read(default_unit_open,nml=analytical_fields_params)
       close(default_unit_open)

       F%AB%Bo = Bo
       F%AB%a = minor_radius
       F%AB%Ro = major_radius
       F%Ro = major_radius
       F%Zo = 0.0_rp
       F%AB%qa = qa
       F%AB%qo = qo
       F%AB%lambda = F%AB%a/SQRT(qa/qo - 1.0_rp)
       F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)
       F%AB%current_direction = TRIM(current_direction)
       SELECT CASE (TRIM(F%AB%current_direction))
       CASE('PARALLEL')
          F%AB%Bp_sign = 1.0_rp
       CASE('ANTI-PARALLEL')
          F%AB%Bp_sign = -1.0_rp
       CASE DEFAULT
       END SELECT
       F%Eo = Eo
       F%Bo = F%AB%Bo

       if (params%field_eval.eq.'interp') then
          F%dims(1) = nR
          F%dims(3) = nZ

          F%Bfield= .TRUE.
          F%axisymmetric_fields = .TRUE.
          F%Bflux=.FALSE.
          F%Efield=.FALSE.

          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,.TRUE.,.FALSE.,.FALSE.)

          do ii=1_idef,F%dims(1)
             F%X%R(ii)=(F%Ro-F%AB%a)+(ii-1)*2*F%AB%a/(F%dims(1)-1)
          end do
          do ii=1_idef,F%dims(3)
             F%X%Z(ii)=(F%Zo-F%AB%a)+(ii-1)*2*F%AB%a/(F%dims(3)-1)
          end do
          
          do ii=1_idef,F%dims(1)
             do kk=1_idef,F%dims(3)
                rm=sqrt((F%X%R(ii)-F%Ro)**2+(F%X%Z(kk)-F%Zo)**2)
                qr=F%AB%qo*(1+(rm/F%AB%lambda)**2)
                theta=atan2(F%X%Z(kk)-F%Zo,F%X%R(ii)-F%Ro)
                F%B_2D%R(ii,kk)=(rm/F%X%R(ii))* &
                     (F%AB%Bo/qr)*sin(theta)
                F%B_2D%PHI(ii,kk)=-(F%Ro/F%X%R(ii))*F%AB%Bo
                F%B_2D%Z(ii,kk)=-(rm/F%X%R(ii))* &
                     (F%AB%Bo/qr)*cos(theta)
                !! Sign convention in analytical fields corresponds to
                !! DIII-D fields with \(B_\phi<0\) and \(B_\theta<0\).
                F%FLAG2D=1.
             end do
          end do         

          if (params%orbit_model(3:5).eq.'pre') call initialize_GC_fields(F)
                    
       end if
       
    CASE('EXTERNAL')
       ! Load the magnetic field from an external HDF5 file
       open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
            status='OLD',form='formatted')
       read(default_unit_open,nml=externalPlasmaModel)
       close(default_unit_open)

       F%Bfield = Bfield
       F%Bflux = Bflux
       F%Efield = Efield
       F%axisymmetric_fields = axisymmetric_fields

       call load_dim_data_from_hdf5(params,F)
       !sets F%dims for 2D or 3D data

       call which_fields_in_file(params,F%Bfield_in_file,F%Efield_in_file, &
            F%Bflux_in_file)

       if (F%Bflux.AND..NOT.F%Bflux_in_file) then
          write(6,'("ERROR: Magnetic flux to be used but no data in file!")')
          call KORC_ABORT()
       end if

       if (F%Bfield.AND..NOT.F%Bfield_in_file) then
          write(6,'("ERROR: Magnetic field to be used but no data in file!")')
          call KORC_ABORT()
       end if

       if (F%Efield.AND..NOT.F%Efield_in_file) then
          if (params%mpi_params%rank.EQ.0_idef) then
             write(6,'(/,"* * * * * * * * * *  FIELDS  * * * * * * * * * *")')
             write(6,'("MESSAGE: Analytical electric field will be used.")')
             write(6,'("* * * * * * * * * * * * ** * * * * * * * * * * *",/)')
          end if
       end if

       if (F%axisymmetric_fields) then
          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield,F%Bflux,F%Efield.AND. &
               F%Efield_in_file)
       else
          call ALLOCATE_3D_FIELDS_ARRAYS(F,F%Bfield,F%Efield)
       end if
       !allocates 2D or 3D data arrays (fields and spatial)
       
       call load_field_data_from_hdf5(params,F)
      
       if (params%orbit_model(3:5).EQ.'pre') call initialize_GC_fields(F)


    CASE DEFAULT
    END SELECT
  end subroutine initialize_fields

  subroutine initialize_GC_fields(F)
    !! Computes the auxiliary fields \(\nabla|{\bf B}|\) and
    !! \(\nabla\times\hat{b}\) that are used in the RHS of the
    !! evolution equations for the GC orbit model.
    TYPE(FIELDS), INTENT(INOUT)      :: F
    !! An instance of the KORC derived type FIELDS.
    INTEGER                        :: ii
    !! Iterator across F%dim
    REAL(rp), DIMENSION(:,:),ALLOCATABLE :: Bmag
    !! Magnetic field magnitude
    REAL(rp), DIMENSION(:,:,:),ALLOCATABLE :: bhat
    !! Magnetic field unit vector
    
    Bmag=SQRT(F%B_2D%R**2+F%B_2D%PHI**2+F%B_2D%Z**2)

    ALLOCATE(bhat(F%dims(1),F%dims(3),3))

    bhat(:,:,1)=F%B_2D%R/Bmag
    bhat(:,:,2)=F%B_2D%PHI/Bmag
    bhat(:,:,3)=F%B_2D%Z/Bmag


    F%gradB_2D%PHI=0.
    ! No variation in phi direction

    ! Single-sided difference for axiliary fields at edge nodes
    ! Differential over R on first index, differential over Z
    ! on second index.

    ! gradB
    ! edge nodes at minimum R,Z
    F%gradB_2D%R(1,:)=(Bmag(2,:)-Bmag(1,:))/(F%X%R(2)-F%X%R(1))
    F%gradB_2D%Z(:,1)=(Bmag(:,2)-Bmag(:,1))/(F%X%Z(2)-F%X%Z(1))

    ! edge nodes at maximum R,Z
    F%gradB_2D%R(F%dims(1),:)=(Bmag(F%dims(1),:)-Bmag(F%dims(1)-1,:))/ &
         (F%X%R(F%dims(1))-F%X%R(F%dims(1)-1))
    F%gradB_2D%Z(:,F%dims(3))=(Bmag(:,F%dims(3))-Bmag(:,F%dims(3)-1))/ &
         (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1))

    ! curlb
    ! edge nodes at minimum R,Z
    ! R component has differential over Z
    F%curlb_2D%R(:,1)=-(bhat(:,2,2)-bhat(:,1,2))/ &
         (F%X%Z(2)-F%X%Z(1))

    ! PHI component has differentials over R and Z
    F%curlb_2D%PHI(1,:)=-(bhat(2,:,3)-bhat(1,:,3))/ &
         (F%X%R(2)-F%X%R(1))         

    F%curlb_2D%PHI(:,1)=F%curlb_2D%PHI(:,1)+ &
         ((bhat(:,2,1)-bhat(:,1,1))/(F%X%Z(2)-F%X%Z(1)))

    ! Z component has differentials over R
    F%curlb_2D%Z(1,:)=((bhat(2,:,2)*F%X%R(2)- &
         bhat(1,:,2)*F%X%R(1))/(F%X%R(2)-F%X%R(1)))/F%X%R(1)

    ! edge nodes at minimum R,Z
    ! R component has differential over Z
    F%curlb_2D%R(:,F%dims(3))=-(bhat(:,F%dims(3),2)- &
         bhat(:,F%dims(3)-1,2))/ &
         (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1))

    ! PHI component has differentials over R and Z
    F%curlb_2D%PHI(F%dims(1),:)=F%curlb_2D%PHI(F%dims(1),:)- &
         (bhat(F%dims(1),:,3)-bhat(F%dims(1)-1,:,3))/ &
         (F%X%R(F%dims(1))-F%X%R(F%dims(1)-1))         

    F%curlb_2D%PHI(:,F%dims(3))=F%curlb_2D%PHI(:,F%dims(3))+ &
         ((bhat(:,F%dims(3),1)-bhat(:,F%dims(3)-1,1))/ &
         (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1)))

    ! Z component has differentials over R
    F%curlb_2D%Z(F%dims(1),:)=((bhat(F%dims(1),:,2)*F%X%R(F%dims(1))- &
         bhat(F%dims(1)-1,:,2)*F%X%R(F%dims(1)-1))/(F%X%R(F%dims(1))- &
         F%X%R(F%dims(1)-1)))/F%X%R(F%dims(1))
    
    do ii=2_idef,F%dims(1)-1
       ! central difference over R for interior nodes
       F%gradB_2D%R(ii,:)=(Bmag(ii+1,:)-Bmag(ii-1,:))/ &
            (F%X%R(ii+1)-F%X%R(ii-1))
       F%curlb_2D%Z(ii,:)=((bhat(ii+1,:,2)*F%X%R(ii+1)- &
            bhat(ii-1,:,2)*F%X%R(ii-1))/(F%X%R(ii+1)-F%X%R(ii-1)))/ &
            F%X%R(ii)
       F%curlb_2D%PHI(ii,:)=F%curlb_2D%PHI(ii,:)- &
            (bhat(ii+1,:,3)-bhat(ii-1,:,3))/ &
            (F%X%R(ii+1)-F%X%R(ii-1))
    end do
    do ii=2_idef,F%dims(3)-1
       ! central difference over Z for interior nodes
       F%gradB_2D%Z(:,ii)=(Bmag(:,ii+1)-Bmag(:,ii-1))/ &
            (F%X%Z(ii+1)-F%X%Z(ii-1))
       F%curlb_2D%R(:,ii)=-(bhat(:,ii+1,2)-bhat(:,ii-1,2))/ &
            (F%X%Z(ii+1)-F%X%Z(ii-1))
       F%curlb_2D%PHI(:,ii)=F%curlb_2D%PHI(:,ii)+ &
            ((bhat(:,ii+1,1)-bhat(:,ii-1,1))/(F%X%Z(ii+1)-F%X%Z(ii-1)))
    end do
    
    DEALLOCATE(Bmag)
    DEALLOCATE(bhat) 
    
  end subroutine initialize_GC_fields

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Subroutines for getting the fields data from HDF5 files
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  !> @brief Subroutine that loads the size of the arrays having the electric and magnetic field data.
  !! @details All the information of externally calculated fields must be given in a rectangular, equally spaced mesh in the \((R,\phi,Z)\) space of cylindrical coordinates.
  !! If the fields are axisymmetric, then the fields must be in a rectangular mesh on the \(RZ\)-plane.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] F An instance of the KORC derived type FIELDS.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @dims Array containing the size of the mesh with the data of the electric and magnetic fields. dims(1) = dimension along the \(R\) coordinate,
  !! dims(2) = dimension along the \(\phi\) coordinate, and dims(3) = dimension along the \(Z\) coordinate.
  !! @param h5error HDF5 error status.
  !! @param rdamum Temporary variable keeping the read data.
  subroutine load_dim_data_from_hdf5(params,F)
    TYPE(KORC_PARAMS), INTENT(IN)                  :: params
    TYPE(FIELDS), INTENT(INOUT)                    :: F
    CHARACTER(MAX_STRING_LENGTH)                   :: filename
    CHARACTER(MAX_STRING_LENGTH)                   :: gname
    CHARACTER(MAX_STRING_LENGTH)                   :: subgname
    CHARACTER(MAX_STRING_LENGTH)                   :: dset
    INTEGER(HID_T)                                 :: h5file_id
    INTEGER(HID_T)                                 :: group_id
    INTEGER(HID_T)                                 :: subgroup_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE    :: dims
    INTEGER                                        :: h5error
    REAL(rp)                                       :: rdatum

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(6,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fopen_f")')
    end if

    if (F%Bflux.OR.F%axisymmetric_fields) then
       dset = "/NR"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(1) = INT(rdatum)

       F%dims(2) = 0

       dset = "/NZ"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(3) = INT(rdatum)
    else
       dset = "/NR"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(1) = INT(rdatum)

       dset = "/NPHI"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(2) = INT(rdatum)

       dset = "/NZ"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(3) = INT(rdatum)
    end if


    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(6,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine load_dim_data_from_hdf5


  !> @brief Subroutine that queries the HDF5 file what data are present in the HDF5 input file (sanity check).
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param Bfield Logical variable that specifies if the magnetic field is present in the HDF5 file.
  !! @param Efield Logical variable that specifies if the electric field is present in the HDF5 file.
  !! @param Bflux Logical variable that specifies if the poloidal magnetic flux is present in the HDF5 file.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @param h5error HDF5 error status.
  subroutine which_fields_in_file(params,Bfield,Efield,Bflux)
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    LOGICAL, INTENT(OUT)               :: Bfield
    LOGICAL, INTENT(OUT)               :: Efield
    LOGICAL, INTENT(OUT)               :: Bflux
    CHARACTER(MAX_STRING_LENGTH)       :: filename
    CHARACTER(MAX_STRING_LENGTH)       :: gname
    CHARACTER(MAX_STRING_LENGTH)       :: subgname
    CHARACTER(MAX_STRING_LENGTH)       :: dset
    INTEGER(HID_T)                     :: h5file_id
    INTEGER(HID_T)                     :: group_id
    INTEGER(HID_T)                     :: subgroup_id
    INTEGER                            :: h5error

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
    end if

    gname = "BR"
    call h5lexists_f(h5file_id,TRIM(gname),Bfield,h5error)

    gname = "ER"
    call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

    gname = "PSIp"
    call h5lexists_f(h5file_id,TRIM(gname),Bflux,h5error)


    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine which_fields_in_file


  !> @brief Subroutine that loads the fields data from the HDF5 input file.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @param h5error HDF5 error status.
  subroutine load_field_data_from_hdf5(params,F)
    TYPE(KORC_PARAMS), INTENT(IN)          :: params
    TYPE(FIELDS), INTENT(INOUT)            :: F
    CHARACTER(MAX_STRING_LENGTH)           :: filename
    CHARACTER(MAX_STRING_LENGTH)           :: gname
    CHARACTER(MAX_STRING_LENGTH)           :: subgname
    CHARACTER(MAX_STRING_LENGTH)           :: dset
    INTEGER(HID_T)                         :: h5file_id
    INTEGER(HID_T)                         :: group_id
    INTEGER(HID_T)                         :: subgroup_id
    INTEGER                                :: h5error

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
    end if

    dset = "/R"
    call load_array_from_hdf5(h5file_id,dset,F%X%R)

    if ((.NOT.F%Bflux).AND.(.NOT.F%axisymmetric_fields)) then
       dset = "/PHI"
       call load_array_from_hdf5(h5file_id,dset,F%X%PHI)
    end if

    dset = "/Z"
    call load_array_from_hdf5(h5file_id,dset,F%X%Z)

    dset = '/Bo'
    call load_from_hdf5(h5file_id,dset,F%Bo)

    if (F%Efield) then
       dset = '/Eo'
       call load_from_hdf5(h5file_id,dset,F%Eo)
    else
       F%Eo = 0.0_rp
    end if

    dset = '/Ro'
    call load_from_hdf5(h5file_id,dset,F%Ro)

    dset = '/Zo'
    call load_from_hdf5(h5file_id,dset,F%Zo)

    if (F%Bflux.OR.F%axisymmetric_fields) then
       dset = "/FLAG"
       call load_array_from_hdf5(h5file_id,dset,F%FLAG2D)
    else
       dset = "/FLAG"
       call load_array_from_hdf5(h5file_id,dset,F%FLAG3D)
    end if

    if (F%Bflux) then
       dset = "/PSIp"
       call load_array_from_hdf5(h5file_id,dset,F%PSIp)
    end if

    if (F%Bfield) then
       if (F%axisymmetric_fields) then
          dset = "/BR"
          call load_array_from_hdf5(h5file_id,dset,F%B_2D%R)

          dset = "/BPHI"
          call load_array_from_hdf5(h5file_id,dset,F%B_2D%PHI)

          dset = "/BZ"
          call load_array_from_hdf5(h5file_id,dset,F%B_2D%Z)
       else
          dset = "/BR"
          call load_array_from_hdf5(h5file_id,dset,F%B_3D%R)

          dset = "/BPHI"
          call load_array_from_hdf5(h5file_id,dset,F%B_3D%PHI)

          dset = "/BZ"
          call load_array_from_hdf5(h5file_id,dset,F%B_3D%Z)
       end if
    end if

    if (F%Efield.AND.F%Efield_in_file) then
       if (F%axisymmetric_fields) then
          dset = "/ER"
          call load_array_from_hdf5(h5file_id,dset,F%E_2D%R)

          dset = "/EPHI"
          call load_array_from_hdf5(h5file_id,dset,F%E_2D%PHI)

          dset = "/EZ"
          call load_array_from_hdf5(h5file_id,dset,F%E_2D%Z)
       else
          dset = "/ER"
          call load_array_from_hdf5(h5file_id,dset,F%E_3D%R)

          dset = "/EPHI"
          call load_array_from_hdf5(h5file_id,dset,F%E_3D%PHI)

          dset = "/EZ"
          call load_array_from_hdf5(h5file_id,dset,F%E_3D%Z)
       end if
    end if

    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(6,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine load_field_data_from_hdf5








  subroutine ALLOCATE_2D_FIELDS_ARRAYS(params,F,bfield,bflux,efield)
    !! @note Subroutine that allocates the variables keeping the axisymmetric
    !! fields data. @endnote
    TYPE (KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)    :: F
    !! An instance of the KORC derived type FIELDS. In this variable we keep
    !! the loaded data.
    LOGICAL, INTENT(IN)            :: bfield
    !! Logical variable that specifies if the variables that keep the magnetic
    !! field data is allocated (bfield=T) or not (bfield=F).
    LOGICAL, INTENT(IN)            :: bflux
    !! Logical variable that specifies if the variables that keep the poloidal
    !! magnetic flux data is allocated (bflux=T) or not (bflux=F).
    LOGICAL, INTENT(IN)            :: efield
    !! Logical variable that specifies if the variables that keep the electric
    !! field data is allocated (efield=T) or not (efield=F).

    if (bfield) then
       call ALLOCATE_V_FIELD_2D(F%B_2D,F%dims)

       if(params%orbit_model(3:5).EQ.'pre') then
          call ALLOCATE_V_FIELD_2D(F%curlb_2D,F%dims)
          call ALLOCATE_V_FIELD_2D(F%gradB_2D,F%dims)
       end if

    end if

    if (bflux) then
       ALLOCATE(F%PSIp(F%dims(1),F%dims(3)))
    end if

    if (efield) then
       call ALLOCATE_V_FIELD_2D(F%E_2D,F%dims)
    end if

    if (.NOT.ALLOCATED(F%FLAG2D)) ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_2D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the variables keeping the 3-D fields data.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param[in] bfield Logical variable that specifies if the variables that keep the magnetic field data is allocated (bfield=T) or not (bfield=F).
  !! @param[in] efield Logical variable that specifies if the variables that keep the electric field data is allocated (efield=T) or not (efield=F).
  subroutine ALLOCATE_3D_FIELDS_ARRAYS(F,bfield,efield)
    TYPE(FIELDS), INTENT(INOUT)    :: F
    LOGICAL, INTENT(IN)            :: bfield
    LOGICAL, INTENT(IN)            :: efield

    if (bfield) then
       call ALLOCATE_V_FIELD_3D(F%B_3D,F%dims)
    end if

    if (efield) then
       call ALLOCATE_V_FIELD_3D(F%E_3D,F%dims)
    end if

    if (.NOT.ALLOCATED(F%FLAG3D)) ALLOCATE(F%FLAG3D(F%dims(1),F%dims(2),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%PHI)) ALLOCATE(F%X%PHI(F%dims(2)))
    if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_3D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the cylindrical components of an axisymmetric field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_2D(F,dims)
    TYPE(V_FIELD_2D), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(3)))
  end subroutine ALLOCATE_V_FIELD_2D


  !> @brief Subroutine that allocates the cylindrical components of a 3-D field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_3D(F,dims)
    TYPE(V_FIELD_3D), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(2),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(2),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(2),dims(3)))
  end subroutine ALLOCATE_V_FIELD_3D

  !> @brief Subroutine that deallocates all the variables of the electric and magnetic fields.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS.
  subroutine DEALLOCATE_FIELDS_ARRAYS(F)
    TYPE(FIELDS), INTENT(INOUT) :: F

    if (ALLOCATED(F%PSIp)) DEALLOCATE(F%PSIp)

    if (ALLOCATED(F%B_2D%R)) DEALLOCATE(F%B_2D%R)
    if (ALLOCATED(F%B_2D%PHI)) DEALLOCATE(F%B_2D%PHI)
    if (ALLOCATED(F%B_2D%Z)) DEALLOCATE(F%B_2D%Z)

    if (ALLOCATED(F%gradB_2D%R)) DEALLOCATE(F%gradB_2D%R)
    if (ALLOCATED(F%gradB_2D%PHI)) DEALLOCATE(F%gradB_2D%PHI)
    if (ALLOCATED(F%gradB_2D%Z)) DEALLOCATE(F%gradB_2D%Z)

    if (ALLOCATED(F%curlb_2D%R)) DEALLOCATE(F%curlb_2D%R)
    if (ALLOCATED(F%curlb_2D%PHI)) DEALLOCATE(F%curlb_2D%PHI)
    if (ALLOCATED(F%curlb_2D%Z)) DEALLOCATE(F%curlb_2D%Z)
    
    if (ALLOCATED(F%B_3D%R)) DEALLOCATE(F%B_3D%R)
    if (ALLOCATED(F%B_3D%PHI)) DEALLOCATE(F%B_3D%PHI)
    if (ALLOCATED(F%B_3D%Z)) DEALLOCATE(F%B_3D%Z)

    if (ALLOCATED(F%E_2D%R)) DEALLOCATE(F%E_2D%R)
    if (ALLOCATED(F%E_2D%PHI)) DEALLOCATE(F%E_2D%PHI)
    if (ALLOCATED(F%E_2D%Z)) DEALLOCATE(F%E_2D%Z)

    if (ALLOCATED(F%E_3D%R)) DEALLOCATE(F%E_3D%R)
    if (ALLOCATED(F%E_3D%PHI)) DEALLOCATE(F%E_3D%PHI)
    if (ALLOCATED(F%E_3D%Z)) DEALLOCATE(F%E_3D%Z)

    if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
    if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
    if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)

    if (ALLOCATED(F%FLAG2D)) DEALLOCATE(F%FLAG2D)
    if (ALLOCATED(F%FLAG3D)) DEALLOCATE(F%FLAG3D)
  end subroutine DEALLOCATE_FIELDS_ARRAYS
end module korc_fields
