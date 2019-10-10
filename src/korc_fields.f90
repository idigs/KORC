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
       add_analytical_E_p,&
       analytical_fields_GC_p,&
       analytical_fields_Bmag_p,&
       analytical_fields_p,&
       initialize_fields,&
       load_field_data_from_hdf5,&
       load_dim_data_from_hdf5,&
       ALLOCATE_2D_FIELDS_ARRAYS,&
       ALLOCATE_3D_FIELDS_ARRAYS,&
       DEALLOCATE_FIELDS_ARRAYS,&
       calculate_SC_E1D,&
       calculate_SC_p,&
       init_SC_E1D,&
       reinit_SC_E1D,&
       define_SC_time_step
  PRIVATE :: get_analytical_fields,&
       analytical_fields,&
       analytical_fields_GC_init,&
       analytical_fields_GC,&
       uniform_magnetic_field,&
       uniform_electric_field,&
       uniform_fields,&
       cross,&
       analytical_electric_field_cyl,&
       ALLOCATE_V_FIELD_2D,&
       ALLOCATE_V_FIELD_3D,&
       initialize_GC_fields,&
       initialize_GC_fields_3D

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

    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ezeta,Bp,Bzeta,eta,q) &
    !$OMP& SHARED(F,Y,E,B,flag)
    do pp=1_idef,ss
       if ( flag(pp) .EQ. 1_is ) then
          eta = Y(pp,1)/F%Ro
          q = F%AB%qo*(1.0_rp + (Y(pp,1)/F%AB%lambda)**2)
          Bp = F%AB%Bp_sign*eta*F%AB%Bo/(q*(1.0_rp + eta*COS(Y(pp,2))))
          Bzeta = F%AB%Bo/( 1.0_rp + eta*COS(Y(pp,2)) )


          B(pp,1) =  Bzeta*COS(Y(pp,3)) - Bp*SIN(Y(pp,2))*SIN(Y(pp,3))
          B(pp,2) = -Bzeta*SIN(Y(pp,3)) - Bp*SIN(Y(pp,2))*COS(Y(pp,3))
          B(pp,3) = Bp*COS(Y(pp,2))

          if (abs(F%Eo) > 0) then
             Ezeta = -F%Eo/( 1.0_rp + eta*COS(Y(pp,2)) )

             E(pp,1) = Ezeta*COS(Y(pp,3))
             E(pp,2) = -Ezeta*SIN(Y(pp,3))
             E(pp,3) = 0.0_rp
          end if
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine analytical_fields

  subroutine analytical_fields_p(B0,E0,R0,q0,lam,ar,X_X,X_Y,X_Z, &
       B_X,B_Y,B_Z,E_X,E_Y,E_Z,flag_cache)
    REAL(rp),  INTENT(IN)      :: R0,B0,lam,q0,E0,ar
    REAL(rp),  INTENT(IN),DIMENSION(8)      :: X_X,X_Y,X_Z
    REAL(rp),  INTENT(OUT),DIMENSION(8)     :: B_X,B_Y,B_Z
    REAL(rp),  INTENT(OUT),DIMENSION(8)     :: E_X,E_Y,E_Z
    INTEGER(is),  INTENT(INOUT),DIMENSION(8)     :: flag_cache
    REAL(rp),DIMENSION(8)     :: T_R,T_T,T_Z
    REAL(rp),DIMENSION(8)                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp),DIMENSION(8)                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp),DIMENSION(8)                              :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp),DIMENSION(8)                               :: eta
    !! Aspect ratio \(\eta\).
    REAL(rp),DIMENSION(8)                                :: q
    !! Safety profile \(q(r)\).
    REAL(rp),DIMENSION(8)                             :: cT,sT,cZ,sZ
    INTEGER                                      :: cc
    !! Particle chunk iterator.

    call cart_to_tor_check_if_confined_p(ar,R0,X_X,X_Y,X_Z, &
         T_R,T_T,T_Z,flag_cache)

    !$OMP SIMD
    !    !$OMP& aligned(cT,sT,cZ,sZ,eta,q,Bp,Bzeta,B_X,B_Y,B_Z, &
    !    !$OMP& Ezeta,E_X,E_Y,E_Z,T_T,T_Z,T_R)
    do cc=1_idef,8
       cT(cc)=cos(T_T(cc))
       sT(cc)=sin(T_T(cc))
       cZ(cc)=cos(T_Z(cc))
       sZ(cc)=sin(T_Z(cc))

       eta(cc) = T_R(cc)/R0
       q(cc) = q0*(1.0_rp + (T_R(cc)*T_R(cc)/(lam*lam)))
       Bp(cc) = -eta(cc)*B0/(q(cc)*(1.0_rp + eta(cc)*cT(cc)))
       Bzeta(cc) = B0/( 1.0_rp + eta(cc)*cT(cc))


       B_X(cc) =  Bzeta(cc)*cZ(cc) - Bp(cc)*sT(cc)*sZ(cc)
       B_Y(cc) = -Bzeta(cc)*sZ(cc) - Bp(cc)*sT(cc)*cZ(cc)
       B_Z(cc) = Bp(cc)*cT(cc)

       Ezeta(cc) = -E0/( 1.0_rp + eta(cc)*cT(cc))

       E_X(cc) = Ezeta(cc)*cZ(cc)
       E_Y(cc) = -Ezeta(cc)*sZ(cc)
       E_Z(cc) = 0.0_rp
    end do
    !$OMP END SIMD

  end subroutine analytical_fields_p

  subroutine analytical_fields_GC_init(params,F,Y,E,B,gradB,curlb,flag,PSIp)
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
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: PSIp
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
    REAL(rp)    ::  rm,theta

    !    write(6,'("Y: ",E17.10)') Y

    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI, &
    !$OMP theta,PSIp) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag)
    do pp=1_idef,ss
       !       if ( flag(pp) .EQ. 1_is ) then

       rm=sqrt((Y(pp,1)-F%AB%Ro)**2+Y(pp,3)**2)
       theta=atan2(Y(pp,3),(Y(pp,1)-F%AB%Ro))
       qprof = 1.0_rp + (rm/F%AB%lambda)**2

       PSIp(pp)=Y(pp,1)*F%AB%lambda**2*F%Bo/ &
            (2*F%AB%qo*(F%AB%Ro+rm*cos(theta)))* &
            log(1+(rm/F%AB%lambda)**2)
       
       Btmp(1)=F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))
       Btmp(2)=-F%AB%Bo*F%AB%Ro/Y(pp,1)
       Btmp(3)=-F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*qprof*Y(pp,1))


       B(pp,1) = Btmp(1)*COS(Y(pp,2)) - Btmp(2)*SIN(Y(pp,2))
       B(pp,2) = Btmp(1)*SIN(Y(pp,2)) + Btmp(2)*COS(Y(pp,2))
       B(pp,3) = Btmp(3)

       dRBR=-F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))*(1./Y(pp,1)+ &
            2*(Y(pp,1)-F%AB%Ro)/(F%AB%lambda**2*qprof))
       dRBPHI=F%AB%Bo*F%AB%Ro/Y(pp,1)**2
       dRBZ=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(-F%AB%Ro/Y(pp,1)+2*(Y(pp,1)- &
            F%AB%Ro)**2/(F%AB%lambda**2*qprof))
       dZBR=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(1-2*Y(pp,3)*Y(pp,3)/ &
            (F%AB%lambda**2*qprof))
       dZBPHI=0._rp
       dZBZ=F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*Y(pp,1))*2*Y(pp,3)/ &
            ((F%AB%lambda*qprof)**2)

       Bmag=sqrt(B(pp,1)*B(pp,1)+B(pp,2)*B(pp,2)+B(pp,3)*B(pp,3))

       gradB(pp,1)=(B(pp,1)*dRBR+B(pp,2)*dRBPHI+B(pp,3)*dRBZ)/Bmag
       gradB(pp,2)=0._rp
       gradB(pp,3)=(B(pp,1)*dZBR+B(pp,2)*dZBPHI+B(pp,3)*dZBZ)/Bmag

       dRbhatPHI=(Bmag*dRBPHI-B(pp,2)*gradB(pp,1))/Bmag**2
       dRbhatZ=(Bmag*dRBZ-B(pp,3)*gradB(pp,1))/Bmag**2
       dZbhatR=(Bmag*dZBR-B(pp,1)*gradB(pp,3))/Bmag**2
       dZbhatPHI=(Bmag*dZBPHI-B(pp,2)*gradB(pp,3))/Bmag**2

       curlb(pp,1)=-dZbhatPHI
       curlb(pp,2)=dZbhatR-dRbhatZ
       curlb(pp,3)=B(pp,2)/(Bmag*Y(pp,1))+dRbhatPHI          


       
       !          if (abs(F%Eo) > 0) then
       E(pp,1) = 0.0_rp
       E(pp,2) = F%Eo*F%AB%Ro/Y(pp,1)
       E(pp,3) = 0.0_rp
       !         end if
       !      end if
    end do
    !$OMP END PARALLEL DO
  end subroutine analytical_fields_GC_init

  subroutine analytical_fields_GC(params,F,Y,E,B,gradB,curlb,flag,PSIp)
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
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: PSIp
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
    REAL(rp)    ::  rm,theta


    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI, &
    !$OMP theta) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag,PSIp)
    do pp=1_idef,ss
       !       if ( flag(pp) .EQ. 1_is ) then

       rm=sqrt((Y(pp,1)-F%AB%Ro)**2+Y(pp,3)**2)
       theta=atan2(Y(pp,3),(Y(pp,1)-F%AB%Ro))
       qprof = 1.0_rp + (rm/F%AB%lambda)**2

!       write(6,*) 'rm: ',rm
!       write(6,*) 'R0: ',F%AB%Ro
!       write(6,*) 'Y_R: ',Y(pp,1)
!       write(6,*) 'theta: ',theta
       
       PSIp(pp)=Y(pp,1)*F%AB%lambda**2*F%Bo/ &
            (2*F%AB%qo*(F%AB%Ro+rm*cos(theta)))* &
            log(1+(rm/F%AB%lambda)**2)

!       write(6,*) 'PSIp: ',PSIp(pp)
       
       Btmp(1)=F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))
       Btmp(2)=-F%AB%Bo*F%AB%Ro/Y(pp,1)
       Btmp(3)=-F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*qprof*Y(pp,1))

       B(pp,1) = Btmp(1)
       B(pp,2) = Btmp(2)
       B(pp,3) = Btmp(3)

       dRBR=-F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))*(1./Y(pp,1)+ &
            2*(Y(pp,1)-F%AB%Ro)/(F%AB%lambda**2*qprof))
       dRBPHI=F%AB%Bo*F%AB%Ro/Y(pp,1)**2
       dRBZ=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(-F%AB%Ro/Y(pp,1)+2*(Y(pp,1)- &
            F%AB%Ro)**2/(F%AB%lambda**2*qprof))
       dZBR=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(1-2*Y(pp,3)*Y(pp,3)/ &
            (F%AB%lambda**2*qprof))
       dZBPHI=0._rp
       dZBZ=F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*Y(pp,1))*2*Y(pp,3)/ &
            ((F%AB%lambda*qprof)**2)

       Bmag=sqrt(B(pp,1)*B(pp,1)+B(pp,2)*B(pp,2)+B(pp,3)*B(pp,3))

       gradB(pp,1)=(B(pp,1)*dRBR+B(pp,2)*dRBPHI+B(pp,3)*dRBZ)/Bmag
       gradB(pp,2)=0._rp
       gradB(pp,3)=(B(pp,1)*dZBR+B(pp,2)*dZBPHI+B(pp,3)*dZBZ)/Bmag

       dRbhatPHI=(Bmag*dRBPHI-B(pp,2)*gradB(pp,1))/Bmag**2
       dRbhatZ=(Bmag*dRBZ-B(pp,3)*gradB(pp,1))/Bmag**2
       dZbhatR=(Bmag*dZBR-B(pp,1)*gradB(pp,3))/Bmag**2
       dZbhatPHI=(Bmag*dZBPHI-B(pp,2)*gradB(pp,3))/Bmag**2

       curlb(pp,1)=-dZbhatPHI
       curlb(pp,2)=dZbhatR-dRbhatZ
       curlb(pp,3)=B(pp,2)/(Bmag*Y(pp,1))+dRbhatPHI          

       !          if (abs(F%Eo) > 0) then
       E(pp,1) = 0.0_rp
       E(pp,2) = F%Eo*F%AB%Ro/Y(pp,1)
       E(pp,3) = 0.0_rp
       !         end if
       !      end if
    end do
    !$OMP END PARALLEL DO

!    write(6,*) 'PSIp: ',PSIp(:)
!    write(6,*) 'B_PHI: ',B(:,2)
    
  end subroutine analytical_fields_GC

  subroutine analytical_fields_Bmag_p(F,Y_R,Y_PHI,Y_Z,Bmag,E_PHI)

    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)  :: R0,B0,lam,q0,EF0
    REAL(rp),DIMENSION(8),INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8) :: B_R,B_PHI,B_Z,rm,qprof
    REAL(rp),DIMENSION(8),INTENT(OUT) :: Bmag,E_PHI
    integer(ip) :: cc

    B0=F%Bo
    EF0=F%Eo
    lam=F%AB%lambda
    R0=F%AB%Ro
    q0=F%AB%qo

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,rm,qprof,Bmag)
    do cc=1_idef,8_idef
       rm(cc)=sqrt((Y_R(cc)-R0)*(Y_R(cc)-R0)+Y_Z(cc)*Y_Z(cc))
       qprof(cc) = 1.0_rp + (rm(cc)*rm(cc)/(lam*lam))

       B_R(cc)=B0*Y_Z(cc)/(q0*qprof(cc)*Y_R(cc))
       B_PHI(cc)=-B0*R0/Y_R(cc)
       B_Z(cc)=-B0*(Y_R(cc)-R0)/(q0*qprof(cc)*Y_R(cc))

       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       E_PHI(cc)=EF0*R0/Y_R(cc)       
    end do
    !$OMP END SIMD

  end subroutine analytical_fields_Bmag_p

  subroutine add_analytical_E_p(params,tt,F,E_PHI)

    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    TYPE(FIELDS), INTENT(IN)                                   :: F
    INTEGER(ip),INTENT(IN)  :: tt
    REAL(rp)  :: E_dyn,E_pulse,E_width,time,arg,arg1
    REAL(rp),DIMENSION(8),INTENT(INOUT) :: E_PHI
    integer(ip) :: cc

    time=(params%it+tt)*params%dt
    
    E_dyn=F%E_dyn
    E_pulse=F%E_pulse
    E_width=F%E_width

!    write(6,'("E_dyn: ",E17.10)') E_dyn
!    write(6,'("E_pulse: ",E17.10)') E_pulse
!    write(6,'("E_width: ",E17.10)') E_width
    
    !$OMP SIMD
    !    !$OMP& aligned(E_PHI)
    do cc=1_idef,8_idef

       arg=(time-E_pulse)**2/(2._rp*E_width**2)
       arg1=10._rp*(time-E_pulse)/(sqrt(2._rp)*E_width)
       
       E_PHI(cc)=E_PHI(cc)+E_dyn*exp(-arg)*(1._rp+erf(-arg1))/2._rp
    end do
    !$OMP END SIMD

  end subroutine add_analytical_E_p

  subroutine analytical_fields_GC_p(F,Y_R,Y_PHI, &
       Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
       gradB_PHI,gradB_Z,PSIp)

    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp),DIMENSION(8),INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: curlB_R,curlB_PHI,curlB_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(8)  :: dRBR,dRBPHI,dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI
    REAL(rp),DIMENSION(8)  :: dRbhatZ,dZbhatR,dZbhatPHI,qprof,rm,theta
    REAL(rp)  :: B0,E0,lam,R0,q0
    integer(ip) :: cc

    B0=F%Bo
    E0=F%Eo
    lam=F%AB%lambda
    R0=F%AB%Ro
    q0=F%AB%qo

    !$OMP SIMD
!    !$OMP& aligned(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlB_R,curlB_PHI,curlB_Z,E_R,E_PHI,E_Z,PSIp)
    do cc=1_idef,8_idef
       rm(cc)=sqrt((Y_R(cc)-R0)*(Y_R(cc)-R0)+Y_Z(cc)*Y_Z(cc))
       theta(cc)=atan2(Y_Z(cc),(Y_R(cc)-R0))
       qprof(cc) = 1.0_rp + (rm(cc)*rm(cc)/(lam*lam))

       PSIp(cc)=Y_R(cc)*lam**2*B0/ &
            (2*q0*(R0+rm(cc)*cos(theta(cc))))* &
            log(1+(rm(cc)/lam)**2)
       
       B_R(cc)=B0*Y_Z(cc)/(q0*qprof(cc)*Y_R(cc))
       B_PHI(cc)=-B0*R0/Y_R(cc)
       B_Z(cc)=-B0*(Y_R(cc)-R0)/(q0*qprof(cc)*Y_R(cc))

       dRBR(cc)=-B0*Y_Z(cc)/(q0*qprof(cc)*Y_R(cc))*(1./Y_R(cc)+ &
            2*(Y_R(cc)-R0)/(lam*lam*qprof(cc)))
       dRBPHI(cc)=B0*R0/(Y_R(cc)*Y_R(cc))
       dRBZ(cc)=B0/(q0*qprof(cc)*Y_R(cc))*(-R0/Y_R(cc)+2*(Y_R(cc)- &
            R0)*(Y_R(cc)-R0)/(lam*lam*qprof(cc)))
       dZBR(cc)=B0/(q0*qprof(cc)*Y_R(cc))*(1-2*Y_Z(cc)*Y_Z(cc)/ &
            (lam*lam*qprof(cc)))
       dZBPHI(cc)=0._rp
       dZBZ(cc)=B0*(Y_R(cc)-R0)/(q0*Y_R(cc))*2*Y_Z(cc)/ &
            (lam*lam*qprof(cc)*qprof(cc))

       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       gradB_R(cc)=(B_R(cc)*dRBR(cc)+B_PHI(cc)*dRBPHI(cc)+B_Z(cc)*dRBZ(cc))/ &
            Bmag(cc)
       gradB_PHI(cc)=0._rp
       gradB_Z(cc)=(B_R(cc)*dZBR(cc)+B_PHI(cc)*dZBPHI(cc)+B_Z(cc)*dZBZ(cc))/ &
            Bmag(cc)

       dRbhatPHI(cc)=(Bmag(cc)*dRBPHI(cc)-B_PHI(cc)*gradB_R(cc))/ &
            (Bmag(cc)*Bmag(cc))
       dRbhatZ(cc)=(Bmag(cc)*dRBZ(cc)-B_Z(cc)*gradB_R(cc))/(Bmag(cc)*Bmag(cc))
       dZbhatR(cc)=(Bmag(cc)*dZBR(cc)-B_R(cc)*gradB_Z(cc))/(Bmag(cc)*Bmag(cc))
       dZbhatPHI(cc)=(Bmag(cc)*dZBPHI(cc)-B_PHI(cc)*gradB_Z(cc))/ &
            (Bmag(cc)*Bmag(cc))

       curlb_R(cc)=-dZbhatPHI(cc)
       curlb_PHI(cc)=dZbhatR(cc)-dRbhatZ(cc)
       curlb_Z(cc)=B_PHI(cc)/(Bmag(cc)*Y_R(cc))+dRbhatPHI(cc)         

       
       E_R(cc) = 0.0_rp
       E_PHI(cc) = E0*R0/Y_R(cc)
       E_Z(cc) = 0.0_rp


       
    end do
    !$OMP END SIMD

  end subroutine analytical_fields_GC_p

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
    B(:,1) = F%Bo
    B(:,2:3) = 0.0_rp
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

    E(:,1) = F%Eo
    E(:,2:3) = 0.0_rp
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
       ss = SIZE(Y,1)
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ephi) SHARED(F,Y,E,flag)
       do pp=1_idef,ss
          if ( flag(pp) .EQ. 1_is ) then
             Ephi = F%Eo*F%Ro/Y(pp,1)

             E(pp,1) = -Ephi*SIN(Y(pp,2))
             E(pp,2) = Ephi*COS(Y(pp,2))
             E(pp,3) = 0.0_rp
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
          Fo = SUM( SQRT(F%B_3D%R**2 + F%B_3D%PHI**2 + F%B_3D%Z**2) )/ &
               SIZE(F%B_3D%R)
       else if (ALLOCATED(F%B_2D%R)) then ! Axisymmetric 2D field
          Fo = SUM( SQRT(F%B_2D%R**2 + F%B_2D%PHI**2 + F%B_2D%Z**2) )/ &
               SIZE(F%B_2D%R)
       end if
    else if (TRIM(op_field%str) .EQ. 'E') then
       if (ALLOCATED(F%E_3D%R)) then ! 3D field
          Fo = SUM( SQRT(F%E_3D%R**2 + F%E_3D%PHI**2 + F%E_3D%Z**2) )/ &
               SIZE(F%E_3D%R)
       else if (ALLOCATED(F%E_2D%R)) then ! Axisymmetric 2D field
          Fo = SUM( SQRT(F%E_2D%R**2 + F%E_2D%PHI**2 + F%E_2D%Z**2) )/ &
               SIZE(F%E_2D%R)
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

       !       call cart_to_cyl(vars%X,vars%Y)

    elseif (params%orbit_model(1:2).eq.'GC') then

       if (.not.params%GC_coords) then

          call cart_to_cyl(vars%X,vars%Y)

          call cyl_check_if_confined(F,vars%Y,vars%flag)

          call analytical_fields_GC_init(params,F,vars%Y, vars%E, vars%B, &
               vars%gradB,vars%curlb, vars%flag, vars%PSI_P)

       else
          
          call cyl_check_if_confined(F,vars%Y,vars%flag)

          call analytical_fields_GC(params,F,vars%Y, vars%E, vars%B, &
               vars%gradB,vars%curlb, vars%flag,vars%PSI_P)

       end if

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


  pure function cross(a,b)
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

    ppp = SIZE(Xo,1) ! Number of particles

    ALLOCATE( vars%X(ppp,3) )
    ALLOCATE( vars%Y(ppp,3) )
    ALLOCATE( vars%B(ppp,3) )
    ALLOCATE( vars%gradB(ppp,3) )
    ALLOCATE( vars%curlb(ppp,3) )
    ALLOCATE( vars%PSI_P(ppp) )
    ALLOCATE( vars%E(ppp,3) )
    ALLOCATE( vars%flag(ppp) )

    vars%X = Xo
    vars%flag = 1_idef
    vars%B=0._rp
    vars%PSI_P=0._rp


    call init_random_seed()

    call get_fields(params,vars,F)

    !    write(6,'("Bx: ",E17.10)') vars%B(:,1)
    !    write(6,'("By: ",E17.10)') vars%B(:,2)
    !    write(6,'("Bz: ",E17.10)') vars%B(:,3)

    do ii=1_idef,ppp
       if ( vars%flag(ii) .EQ. 1_idef ) then
          b1(ii,:) = vars%B(ii,:)/sqrt(vars%B(ii,1)*vars%B(ii,1)+ &
               vars%B(ii,2)*vars%B(ii,2)+vars%B(ii,3)*vars%B(ii,3))

          b2(ii,:) = cross(b1(ii,:),(/0.0_rp,0.0_rp,1.0_rp/))
          b2(ii,:) = b2(ii,:)/sqrt(b2(ii,1)*b2(ii,1)+b2(ii,2)*b2(ii,2)+ &
               b2(ii,3)*b2(ii,3))

          b3(ii,:) = cross(b1(ii,:),b2(ii,:))
          b3(ii,:) = b3(ii,:)/sqrt(b3(ii,1)*b3(ii,1)+b3(ii,2)*b3(ii,2)+ &
               b3(ii,3)*b3(ii,3))
       end if
    end do

    if (PRESENT(flag)) then
       flag = vars%flag
    end if

    DEALLOCATE( vars%X )
    DEALLOCATE( vars%Y )
    DEALLOCATE( vars%B )
    DEALLOCATE( vars%PSI_P )
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

    SELECT CASE (TRIM(params%field_model))
    CASE('ANALYTICAL')
       if (params%field_eval.eq.'eqn') then
          call get_analytical_fields(params,vars, F)
       else
          call interp_fields(params,vars, F)          
       end if
    CASE('EXTERNAL')

       !       write(6,'("2 size of PSI_P: ",I16)') size(vars%PSI_P)

       call interp_fields(params,vars, F)

!       write(6,'("get_fields")')
!       write(6,'("B_X: ",E17.10)') vars%B(:,1)
!       write(6,'("B_Z: ",E17.10)') vars%B(:,2)
!       write(6,'("B_Y: ",E17.10)') vars%B(:,3)
       
       !if (F%Efield.AND..NOT.F%Efield_in_file) then
       !   call analytical_electric_field_cyl(F,vars%Y,vars%E,vars%flag)
       !end if
    CASE('UNIFORM')

       call uniform_fields(vars, F)
    CASE DEFAULT
    END SELECT
  end subroutine get_fields

  subroutine calculate_SC_E1D(params,F,Vden)

    TYPE(FIELDS), INTENT(INOUT)                 :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params
    real(rp),dimension(F%dim_1D),intent(in) :: Vden
    real(rp),dimension(F%dim_1D) :: Jsamone,Jsamall,Jexp,dJdt
    real(rp),dimension(F%dim_1D) :: a,b,c,u,gam,r
    real(rp) :: dr,bet
    integer :: ii
    INTEGER 				:: mpierr

!    if (params%mpi_params%rank .EQ. 0) then
!       write(6,*) 'Calculating SC_E1D'
!    end if
    
    dr=F%r_1D(2)-F%r_1D(1)
    
    Jsamone=C_E*Vden

    ! Add sampled current densities from all MPI processes Jsamone,
    ! and output of total sampled current density Jsamall to each
    ! MPI process.
    
    call MPI_ALLREDUCE(Jsamone,Jsamall,F%dim_1D,MPI_REAL8,MPI_SUM, &
         MPI_COMM_WORLD,mpierr)

    !write(6,*) 'Jsam: ',Jsamall(1:5)
    
    Jexp=Jsamall*F%Ip0

    F%J3_SC_1D%PHI=F%J2_SC_1D%PHI
    F%J2_SC_1D%PHI=F%J1_SC_1D%PHI
    F%J1_SC_1D%PHI=Jexp    
    
    ! Calculating time-derivative of E_phi
    
    dJdt=(3*F%J1_SC_1D%PHI-4*F%J2_SC_1D%PHI+F%J3_SC_1D%PHI)/ &
         (2*F%dt_E_SC)
    
!    write(6,*) params%mpi_params%rank,'J(1)',F%J_SC_1D%PHI(1)
    
    ! Solving 1D Poisson equation with tridiagonal matrix solve

    a=0._rp
    b=-2._rp
    c=0._rp
    u=0._rp
    gam=0._rp
!    r=-2*dr**2*C_MU*Jexp
    r=2*dr**2*C_MU*dJdt

    do ii=2_idef,F%dim_1D
       a(ii)=(REAL(ii)-2._rp)/(REAL(ii)-1._rp)
       c(ii)=REAL(ii)/(REAL(ii)-1._rp)
    end do

    bet=b(2)
    u(2)=r(2)/bet
    do ii=3_idef,F%dim_1D-1
       gam(ii)=c(ii-1)/bet
       bet=b(ii)-a(ii)*gam(ii)
       if (bet.eq.0) then
          stop 'tridiag failed'         
       end if
       u(ii)=(r(ii)-a(ii)*u(ii-1))/bet
    end do

    do ii=F%dim_1D-2,2,-1
       u(ii)=u(ii)-gam(ii+1)*u(ii+1)
    end do

    u(1)=(4*u(2)-u(3))/3._rp

    ! Writing over F%A* data
    
!    F%A3_SC_1D%PHI=F%A2_SC_1D%PHI
!    F%A2_SC_1D%PHI=F%A1_SC_1D%PHI
!    F%A1_SC_1D%PHI=u
    
!    if (init) then
!       F%A3_SC_1D%PHI=F%A1_SC_1D%PHI
!       F%A2_SC_1D%PHI=F%A1_SC_1D%PHI 
!    end if

!    write(6,*) params%mpi_params%rank,'A1(1)',F%A1_SC_1D%PHI(1)
!    write(6,*) params%mpi_params%rank,'A2(1)',F%A2_SC_1D%PHI(1)
!    write(6,*) params%mpi_params%rank,'A3(1)',F%A3_SC_1D%PHI(1)
    
    ! Calculating inductive E_phi

!    F%E_SC_1D%PHI=-(3*F%A1_SC_1D%PHI-4*F%A2_SC_1D%PHI+F%A3_SC_1D%PHI)/ &
!         (2*F%dt_E_SC)

    F%E_SC_1D%PHI=u

    if (params%mpi_params%rank.eq.0) then
       write(6,*) 'J1(2)',F%J1_SC_1D%PHI(2)
       write(6,*) 'J2(2)',F%J2_SC_1D%PHI(2)
       write(6,*) 'J3(2)',F%J3_SC_1D%PHI(2)
       
       write(6,*) 'E(1)',F%E_SC_1D%PHI(1)
    end if
       
    ! Normalizing inductive E_phi
    
    F%E_SC_1D%PHI=F%E_SC_1D%PHI/params%cpp%Eo

    call initialize_SC1D_field_interpolant(params,F)
    
  end subroutine calculate_SC_E1D

  subroutine calculate_SC_p(params,F,B_R,B_PHI,B_Z,Y_R,Y_Z, &
       V_PLL,V_MU,m_cache,flag_cache,Vden)

    TYPE(FIELDS), INTENT(IN)                 :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params
    real(rp),dimension(8),intent(in) :: Y_R,Y_Z
    real(rp),dimension(8),intent(in) :: B_R,B_PHI,B_Z
    real(rp),dimension(8),intent(in) :: V_PLL,V_MU
    real(rp),intent(in) :: m_cache
    integer(is),dimension(8),intent(in) :: flag_cache
    real(rp),dimension(8) :: rm,Bmag,gam,vpll
    real(rp),dimension(F%dim_1D),intent(out) :: Vden
    real(rp),dimension(F%dim_1D) :: Vpart,Ai
    real(rp),dimension(F%dim_1D) :: r_1D
    real(rp) :: dr,sigr,ar,arg,arg1,arg2,arg3
    integer :: cc,ii,rind


    dr=F%r_1D(2)-F%r_1D(1)
    r_1D=F%r_1D
    sigr=dr

    Vpart=0._rp
    do cc=1_idef,8_idef

       ! 1D nearest grid point weighting in minor radius

       !    RR=spp%vars%Y(:,1)
       !    ZZ=spp%vars%Y(:,3)
       rm(cc)=sqrt((Y_R(cc)-F%Ro)**2+(Y_Z(cc)-F%Zo)**2)* &
            params%cpp%length

       !    write (6,*) params%mpi_params%rank,'RR',RR
       !    write (6,*) params%mpi_params%rank,'ZZ',spp%vars%Y(:,3)
       write (6,*) 'rm',rm(cc)

       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
            B_Z(cc)*B_Z(cc))
       gam(cc)=sqrt(1+V_PLL(cc)**2+ &
            2*V_MU(cc)*Bmag(cc)*m_cache)
       vpll(cc)=V_PLL(cc)/gam(cc)    
       

       ! Weighting parallel velocity

       !    write (6,*) params%mpi_params%rank,'vpll',vpll

       !   do pp=1_idef,spp%ppp
       ! NGP weighting
       rind=FLOOR((rm(cc)-dr/2)/dr)+2_ip
       Vpart(rind)=Vpart(rind)+real(flag_cache(cc))*vpll(cc)

       ! First-order weighting
!       rind=FLOOR(rm(cc)/dr)+1_ip
!       Vpart(rind)=Vpart(rind)+ &
!            vpll(cc)*(r_1D(rind+1)-rm(cc))/dr
!       Vpart(rind+1)=Vpart(rind+1)+ &
!            vpll(cc)*(rm(cc)-r_1D(rind))/dr

       ! Gaussian weighting
       
!       do ii=1_idef,F%dim_1D
!          arg=MIN((r_1D(ii)-rm(cc))**2._rp/(2._rp*sigr**2._rp),100._rp)
!          Vpart(ii)=Vpart(ii)+1/sqrt(2._rp*C_PI*sigr**2._rp)* &
!               exp(-arg)           
!       end do

    end do

    ar=F%AB%a
    ! Calculating density of minor radial annulus    
    do ii=1_idef,F%dim_1D
       ! NGP weighting
       if(ii.eq.1) then
          Vden(ii)=Vpart(ii)/(C_PI*dr**2/4)
       else
          Vden(ii)=Vpart(ii)/(2*C_PI*dr**2*(ii-1))
       end if
       ! First-order weighting
!       if(ii.eq.1) then
!          Vden(ii)=Vpart(ii)/(C_PI*dr**2/3)
!       else
!          Vden(ii)=Vpart(ii)/(2*C_PI*dr**2*(ii-1))
!       end if

       ! Gaussian weighting
!       arg=MIN(r_1D(ii)**2._rp/(2._rp*sigr**2._rp),100._rp)
!       arg1=MIN((ar-r_1D(ii))**2._rp/(2._rp*sigr**2._rp),100._rp)
!       arg2=MIN((ar-r_1D(ii))/(sqrt(2._rp)*sigr),10._rp)
!       arg3=MIN((r_1D(ii))/(sqrt(2._rp)*sigr),10._rp)

!       Ai(ii)=sqrt(C_PI*sigr)*(sqrt(2._rp)*sigr*(exp(-arg)- &
!            exp(-arg1))+r_1D(ii)*sqrt(C_PI)* &
!            (erf(arg2)-erf(-arg3)))
!       Vden(ii)=Vpart(ii)/Ai(ii)

    end do
    
  end subroutine calculate_SC_p

  subroutine init_SC_E1D(params,F,spp)
 
    TYPE(FIELDS), INTENT(INOUT)                 :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params
    TYPE(SPECIES), INTENT(IN)    :: spp
    real(rp),dimension(F%dim_1D) :: Vpart
    real(rp),dimension(spp%ppp) :: RR,ZZ,rm,vpll
    real(rp),dimension(F%dim_1D) :: Vden,Jsamone,Jsamall,Jexp,dJdt
    real(rp),dimension(F%dim_1D) :: a,b,c,u,gam,r,r_1D,Ai
    real(rp) :: dr,Isam,bet,sigr,ar,arg,arg1,arg2,arg3
    integer :: pp,ii,rind
    INTEGER 				:: mpierr

!    if (params%mpi_params%rank .EQ. 0) then
!       write(6,*) 'Calculating SC_E1D'
!    end if

    ! 1D nearest grid point weighting in minor radius
    
    RR=spp%vars%Y(:,1)
    ZZ=spp%vars%Y(:,3)
    rm=sqrt((RR-F%Ro)**2._rp+(ZZ-F%Zo)**2._rp)*params%cpp%length

!    write (6,*) params%mpi_params%rank,'RR',RR
!    write (6,*) params%mpi_params%rank,'ZZ',spp%vars%Y(:,3)
    write (6,*) 'rm',rm
    
    
    dr=F%r_1D(2)-F%r_1D(1)

    vpll=spp%vars%V(:,1)/spp%vars%g

    ! Weighting parallel velocity

!    write (6,*) 'vpll',vpll
    
    Vpart=0._rp
    r_1D=F%r_1D
    sigr=dr
    
    do pp=1_idef,spp%ppp
       ! NGP weighting
       rind=FLOOR((rm(pp)-dr/2)/dr)+2_ip
       Vpart(rind)=Vpart(rind)+vpll(pp)

       ! First-order weighting
!       rind=FLOOR(rm(pp)/dr)+1_ip
!       Vpart(rind)=Vpart(rind)+vpll(pp)*(F%r_1D(rind+1)-rm(pp))/dr
!       Vpart(rind+1)=Vpart(rind+1)+vpll(pp)*(rm(pp)-F%r_1D(rind))/dr

       ! Gaussian weighting
!       do ii=1_idef,F%dim_1D
!          arg=MIN((r_1D(ii)-rm(pp))**2._rp/(2._rp*sigr**2._rp),100._rp)          
!          Vpart(ii)=Vpart(ii)+1/sqrt(2._rp*C_PI*sigr**2._rp)* &
!               exp(-arg) 
!       end do
       
    end do

    ! Calculating density of minor radial annulus

    ar=F%AB%a
    
    do ii=1_idef,F%dim_1D
       ! NGP weighting
       if(ii.eq.1) then
          Vden(ii)=Vpart(ii)/(C_PI*dr**2/4)
       else
          Vden(ii)=Vpart(ii)/(2*C_PI*dr**2*(ii-1))
       end if
       ! First-order weighting
!       if(ii.eq.1) then
!          Vden(ii)=Vpart(ii)/(C_PI*dr**2/3)
!       else
!          Vden(ii)=Vpart(ii)/(2*C_PI*dr**2*(ii-1))
!       end if

       ! Gaussian weighting

!       arg=MIN(r_1D(ii)**2._rp/(2._rp*sigr**2._rp),100._rp)
!       arg1=MIN((ar-r_1D(ii))**2._rp/(2._rp*sigr**2._rp),100._rp)
!       arg2=MIN((ar-r_1D(ii))/(sqrt(2._rp)*sigr),10._rp)
!       arg3=MIN((r_1D(ii))/(sqrt(2._rp)*sigr),10._rp)
       
!       Ai(ii)=sqrt(C_PI*sigr)*(sqrt(2._rp)*sigr*(exp(-arg)- &
!            exp(-arg1))+r_1D(ii)*sqrt(C_PI)* &
!            (erf(arg2)-erf(-arg3)))
!       Vden(ii)=Vpart(ii)/Ai(ii)
       
    end do
    
    Jsamone=C_E*Vden

    ! Add sampled current densities from all MPI processes Jsamone,
    ! and output of total sampled current density Jsamall to each
    ! MPI process.
    
    call MPI_ALLREDUCE(Jsamone,Jsamall,F%dim_1D,MPI_REAL8,MPI_SUM, &
         MPI_COMM_WORLD,mpierr)
    
!    write(6,*) 'Jsam: ',Jsamall(1:10)

    ! Integrating current density to scale total current to
    ! experimentally determined total current


    Isam=0._rp
    do ii=1_idef,F%dim_1D
       if ((ii.eq.1).or.(ii.eq.F%dim_1D)) then
          Isam=Isam+Jsamall(ii)*r_1D(ii)/2._rp
       else 
          Isam=Isam+Jsamall(ii)*r_1D(ii)
       end if
    end do
    Isam=2._rp*C_PI*Isam*dr
!    write(6,*) params%mpi_params%rank,'Isam: ',Isam

    F%Ip0=F%Ip_exp/Isam

    
    Jexp=Jsamall*F%Ip0

    F%J3_SC_1D%PHI=Jexp 
    F%J2_SC_1D%PHI=Jexp
    F%J1_SC_1D%PHI=Jexp        
    
    ! Calculating time-derivative of E_phi
    
    dJdt=(3._rp*F%J1_SC_1D%PHI-4._rp*F%J2_SC_1D%PHI+F%J3_SC_1D%PHI)/ &
         (2._rp*F%dt_E_SC)
    
!    write(6,*) params%mpi_params%rank,'J(1)',F%J_SC_1D%PHI(1)
    
    ! Solving 1D Poisson equation with tridiagonal matrix solve

    a=0._rp
    b=-2._rp
    c=0._rp
    u=0._rp
    gam=0._rp
!    r=-2*dr**2*C_MU*Jexp
    r=2*dr**2*C_MU*dJdt

    do ii=2_idef,F%dim_1D
       a(ii)=(REAL(ii)-2._rp)/(REAL(ii)-1._rp)
       c(ii)=REAL(ii)/(REAL(ii)-1._rp)
    end do

    bet=b(2)
    u(2)=r(2)/bet
    do ii=3_idef,F%dim_1D-1
       gam(ii)=c(ii-1)/bet
       bet=b(ii)-a(ii)*gam(ii)
       if (bet.eq.0) then
          stop 'tridiag failed'         
       end if
       u(ii)=(r(ii)-a(ii)*u(ii-1))/bet
    end do

    do ii=F%dim_1D-2,2,-1
       u(ii)=u(ii)-gam(ii+1)*u(ii+1)
    end do

    u(1)=(4._rp*u(2)-u(3))/3._rp

    ! Writing over F%A* data
    
!    F%A3_SC_1D%PHI=F%A2_SC_1D%PHI
!    F%A2_SC_1D%PHI=F%A1_SC_1D%PHI
!    F%A1_SC_1D%PHI=u
    
!    if (init) then
!       F%A3_SC_1D%PHI=F%A1_SC_1D%PHI
!       F%A2_SC_1D%PHI=F%A1_SC_1D%PHI 
!    end if

!    write(6,*) params%mpi_params%rank,'A1(1)',F%A1_SC_1D%PHI(1)
!    write(6,*) params%mpi_params%rank,'A2(1)',F%A2_SC_1D%PHI(1)
!    write(6,*) params%mpi_params%rank,'A3(1)',F%A3_SC_1D%PHI(1)
    
    ! Calculating inductive E_phi

!    F%E_SC_1D%PHI=-(3*F%A1_SC_1D%PHI-4*F%A2_SC_1D%PHI+F%A3_SC_1D%PHI)/ &
!         (2*F%dt_E_SC)

    F%E_SC_1D%PHI=u

    if (params%mpi_params%rank.eq.0) then
       write(6,*) 'J1(2)',F%J1_SC_1D%PHI(2)
       write(6,*) 'J2(2)',F%J2_SC_1D%PHI(2)
       write(6,*) 'J3(2)',F%J3_SC_1D%PHI(2)
       
       write(6,*) 'E(1)',F%E_SC_1D%PHI(1)
    end if
       
    ! Normalizing inductive E_phi
    
    F%E_SC_1D%PHI=F%E_SC_1D%PHI/params%cpp%Eo

    call initialize_SC1D_field_interpolant(params,F)
    
  end subroutine init_SC_E1D
  
  subroutine reinit_SC_E1D(params,F)
 
    TYPE(FIELDS), INTENT(INOUT)                 :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params

    real(rp),dimension(F%dim_1D) :: Jsamall,Jexp,dJdt
    real(rp),dimension(F%dim_1D) :: a,b,c,u,gam,r,r_1D,Ai
    real(rp) :: dr,Isam,bet,sigr,ar,arg,arg1,arg2,arg3
    integer :: pp,ii,rind
    INTEGER 				:: mpierr

!    if (params%mpi_params%rank .EQ. 0) then
!       write(6,*) 'Calculating SC_E1D'
!    end if

    dr=F%r_1D(2)-F%r_1D(1)
    r_1D=F%r_1D
    Jsamall=F%J0_SC_1D%PHI

!    write(6,*) Jsamall
    
    Isam=0._rp
    do ii=1_idef,F%dim_1D
!       write(6,*) Isam
!       write(6,*) ii
!       write(6,*) Jsamall(ii)
!       write(6,*) (ii)
       if ((ii.eq.1_idef).or.(ii.eq.F%dim_1D)) then
          Isam=Isam+Jsamall(ii)*r_1D(ii)/2._rp
       else 
          Isam=Isam+Jsamall(ii)*r_1D(ii)
       end if
    end do
    Isam=2._rp*C_PI*Isam*dr
!    write(6,*) params%mpi_params%rank,'Isam: ',Isam

    F%Ip0=F%Ip_exp/Isam

    
    Jexp=Jsamall*F%Ip0

    F%J3_SC_1D%PHI=Jexp 
    F%J2_SC_1D%PHI=Jexp
    F%J1_SC_1D%PHI=Jexp        
    
    ! Calculating time-derivative of E_phi
    
    dJdt=(3._rp*F%J1_SC_1D%PHI-4._rp*F%J2_SC_1D%PHI+F%J3_SC_1D%PHI)/ &
         (2._rp*F%dt_E_SC)
    
!    write(6,*) params%mpi_params%rank,'J(1)',F%J_SC_1D%PHI(1)
    
    ! Solving 1D Poisson equation with tridiagonal matrix solve

    a=0._rp
    b=-2._rp
    c=0._rp
    u=0._rp
    gam=0._rp
!    r=-2*dr**2*C_MU*Jexp
    r=2*dr**2*C_MU*dJdt

    do ii=2_idef,F%dim_1D
       a(ii)=(REAL(ii)-2._rp)/(REAL(ii)-1._rp)
       c(ii)=REAL(ii)/(REAL(ii)-1._rp)
    end do

    bet=b(2)
    u(2)=r(2)/bet
    do ii=3_idef,F%dim_1D-1
       gam(ii)=c(ii-1)/bet
       bet=b(ii)-a(ii)*gam(ii)
       if (bet.eq.0) then
          stop 'tridiag failed'         
       end if
       u(ii)=(r(ii)-a(ii)*u(ii-1))/bet
    end do

    do ii=F%dim_1D-2,2,-1
       u(ii)=u(ii)-gam(ii+1)*u(ii+1)
    end do

    u(1)=(4._rp*u(2)-u(3))/3._rp

    ! Writing over F%A* data
    
!    F%A3_SC_1D%PHI=F%A2_SC_1D%PHI
!    F%A2_SC_1D%PHI=F%A1_SC_1D%PHI
!    F%A1_SC_1D%PHI=u
    
!    if (init) then
!       F%A3_SC_1D%PHI=F%A1_SC_1D%PHI
!       F%A2_SC_1D%PHI=F%A1_SC_1D%PHI 
!    end if

!    write(6,*) params%mpi_params%rank,'A1(1)',F%A1_SC_1D%PHI(1)
!    write(6,*) params%mpi_params%rank,'A2(1)',F%A2_SC_1D%PHI(1)
!    write(6,*) params%mpi_params%rank,'A3(1)',F%A3_SC_1D%PHI(1)
    
    ! Calculating inductive E_phi

!    F%E_SC_1D%PHI=-(3*F%A1_SC_1D%PHI-4*F%A2_SC_1D%PHI+F%A3_SC_1D%PHI)/ &
!         (2*F%dt_E_SC)

    F%E_SC_1D%PHI=u

!    if (params%mpi_params%rank.eq.0) then
!       write(6,*) 'J1(1)',F%J1_SC_1D%PHI(1)
!       write(6,*) 'J2(1)',F%J2_SC_1D%PHI(1)
!       write(6,*) 'J3(1)',F%J3_SC_1D%PHI(1)
       
!       write(6,*) 'E(1)',F%E_SC_1D%PHI(1)
!    end if
       
    ! Normalizing inductive E_phi
    
    F%E_SC_1D%PHI=F%E_SC_1D%PHI/params%cpp%Eo

    call initialize_SC1D_field_interpolant(params,F)
    
  end subroutine reinit_SC_E1D
  
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
    CHARACTER(MAX_STRING_LENGTH)   :: E_model
    REAL(rp)                       :: Eo,E_dyn,E_pulse,E_width   
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
    logical :: test
    integer :: res_double
    real(rp) :: RMAX,RMIN,ZMAX,ZMIN
    integer :: dim_1D
    real(rp) :: dt_E_SC,Ip_exp
    

    NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
         qa,qo,Eo,current_direction,nR,nZ,dim_1D,dt_E_SC,Ip_exp
    NAMELIST /externalPlasmaModel/ Efield, Bfield, Bflux, &
         axisymmetric_fields, Eo,E_model,E_dyn,E_pulse,E_width,res_double

    if (params%mpi_params%rank .EQ. 0) then
       write(6,'(/,"* * * * * * * * INITIALIZING FIELDS * * * * * * * *")')
    end if

    SELECT CASE (TRIM(params%field_model))

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

       F%E_model=E_model
       F%E_dyn = E_dyn
       F%E_pulse = E_pulse
       F%E_width = E_width

       F%res_double=res_double
       
       if (params%mpi_params%rank .EQ. 0) then
          write(6,'("ANALYTIC")')
          write(6,'("Magnetic field: ",E17.10)') F%Bo
          write(6,'("Electric field: ",E17.10)') F%Eo

       end if


       if (params%field_eval.eq.'interp') then
          F%dims(1) = nR
          F%dims(3) = nZ

          F%Bfield= .TRUE.
          F%axisymmetric_fields = .TRUE.
          F%Bflux=.TRUE.
          F%Efield=.TRUE.

          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield,F%Bflux,F%Efield)
          
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
                F%E_2D%R(ii,kk)=0.0_rp
                F%E_2D%PHI(ii,kk)=-(F%Ro/F%X%R(ii))*F%Eo
                F%E_2D%Z(ii,kk)=0.0_rp

                F%PSIp(ii,kk)=F%X%R(ii)*F%AB%lambda**2*F%Bo/ &
                     (2*F%AB%qo*(F%Ro+rm*cos(theta)))* &
                     log(1+(rm/F%AB%lambda)**2)
                
                !! Sign convention in analytical fields corresponds to
                !! DIII-D fields with \(B_\phi<0\) and \(B_\theta<0\).
                F%FLAG2D=1.
             end do
          end do
          
          F%FLAG2D(1:2,:)=0.
          F%FLAG2D(F%dims(1)-1:F%dims(1),:)=0.
          F%FLAG2D(:,1:2)=0.
          F%FLAG2D(:,F%dims(3)-1:F%dims(3))=0.

          if (params%orbit_model(3:5).eq.'pre') then
             if (params%mpi_params%rank .EQ. 0) then
                write(6,'("Initializing GC fields from analytic EM fields")')
             end if
             call initialize_GC_fields(F)
          end if

          F%Bfield= .FALSE.
          F%axisymmetric_fields = .TRUE.
          F%Bflux=.TRUE.
          F%Efield=.FALSE.
          
       end if

       if (params%SC_E) then

          F%dim_1D=dim_1D
          F%dt_E_SC=dt_E_SC
          F%Ip_exp=Ip_exp

          ALLOCATE(F%E_SC_1D%PHI(F%dim_1D))
          ALLOCATE(F%A1_SC_1D%PHI(F%dim_1D))
          ALLOCATE(F%A2_SC_1D%PHI(F%dim_1D))
          ALLOCATE(F%A3_SC_1D%PHI(F%dim_1D))
          ALLOCATE(F%J1_SC_1D%PHI(F%dim_1D))
          ALLOCATE(F%J2_SC_1D%PHI(F%dim_1D))
          ALLOCATE(F%J3_SC_1D%PHI(F%dim_1D))          
          ALLOCATE(F%r_1D(F%dim_1D))

          F%E_SC_1D%PHI=0._rp
          F%A1_SC_1D%PHI=0._rp
          F%A2_SC_1D%PHI=0._rp
          F%A3_SC_1D%PHI=0._rp
          F%J1_SC_1D%PHI=0._rp
          F%J2_SC_1D%PHI=0._rp
          F%J3_SC_1D%PHI=0._rp          
          F%r_1D=0._rp
                    
          do ii=1_idef,F%dim_1D
             F%r_1D(ii)=(ii-1)*F%AB%a/(F%dim_1D-1)
          end do
          
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

       F%E_model=E_model
       F%E_dyn = E_dyn
       F%E_pulse = E_pulse
       F%E_width = E_width

       F%res_double=res_double

!       write(6,'("E_dyn: ",E17.10)') E_dyn
!       write(6,'("E_pulse: ",E17.10)') E_pulse
!       write(6,'("E_width: ",E17.10)') E_width
       
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

          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
               F%Bflux,F%Efield.AND. &
               F%Efield_in_file)

       else
          call ALLOCATE_3D_FIELDS_ARRAYS(params,F,F%Bfield,F%Efield)
          
       end if
       !allocates 2D or 3D data arrays (fields and spatial)

       call load_field_data_from_hdf5(params,F)
       
       if (F%Bflux) F%Eo = Eo

       if (.not.axisymmetric_fields) then
          if (.not.F%Efield_in_file) then
             F%Eo = Eo
          
             F%E_3D%R=0._rp
             do ii=1_idef,F%dims(1)
                F%E_3D%PHI(ii,:,:)=F%Eo*F%Ro/F%X%R(ii)
             end do
             F%E_3D%Z=0._rp
          end if
       end if
       
       if (params%mpi_params%rank .EQ. 0) then

          write(6,'("EXTERNAL")')
          write(6,'("Magnetic field: ",E17.10)') F%Bo
          write(6,'("Electric field: ",E17.10)') F%Eo

       end if

       test=.true.
       
       if (F%Bflux.and.(.not.test)) then

          call initialize_fields_interpolant(params,F)

          F%Bfield=.TRUE.
          F%Efield=.TRUE.
          F%Efield_in_file=.TRUE.


          RMIN=F%X%R(1)
          RMAX=F%X%R(F%dims(1))

          ZMIN=F%X%Z(1)
          ZMAX=F%X%Z(F%dims(3))

          do ii=1_idef,res_double
             F%dims(1)=2*F%dims(1)-1
             F%dims(3)=2*F%dims(3)-1
          end do

          if (res_double>0) then
             DEALLOCATE(F%X%R)
             DEALLOCATE(F%X%Z)
             DEALLOCATE(F%PSIp)
          end if
          
          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
               F%Bflux,F%Efield.AND.F%Efield_in_file)

          do ii=1_idef,F%dims(1)
             F%X%R(ii)=RMIN+REAL(ii-1)/REAL(F%dims(1)-1)*(RMAX-RMIN)
          end do

          do ii=1_idef,F%dims(3)
             F%X%Z(ii)=ZMIN+REAL(ii-1)/REAL(F%dims(3)-1)*(ZMAX-ZMIN)
          end do            
          
          call calculate_initial_magnetic_field(F)
          
          F%E_2D%R=0._rp
          do ii=1_idef,F%dims(1)
             F%E_2D%PHI(ii,:)=F%Eo*F%Ro/F%X%R(ii)
          end do
          F%E_2D%Z=0._rp
          
       end if

       
       if (F%Bflux.and.test) then

          F%Bfield=.TRUE.
          F%Efield=.TRUE.
          F%Efield_in_file=.TRUE.
            
          
          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
               F%Bflux,F%Efield.AND.F%Efield_in_file)

          
          ! B
          ! edge nodes at minimum R,Z
          F%B_2D%Z(1,:)=-(F%PSIp(2,:)-F%PSIp(1,:))/(F%X%R(2)-F%X%R(1))/F%X%R(1)
          F%B_2D%R(:,1)=(F%PSIp(:,2)-F%PSIp(:,1))/(F%X%Z(2)-F%X%Z(1))/F%X%R(:)

          ! edge nodes at maximum R,Z
          F%B_2D%Z(F%dims(1),:)=-(F%PSIp(F%dims(1),:)-F%PSIp(F%dims(1)-1,:))/ &
               (F%X%R(F%dims(1))-F%X%R(F%dims(1)-1))/F%X%R(F%dims(1))
          F%B_2D%R(:,F%dims(3))=(F%PSIp(:,F%dims(3))-F%PSIp(:,F%dims(3)-1))/ &
               (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1))/F%X%R(:)

          do ii=2_idef,F%dims(1)-1
             ! central difference over R for interior nodes for BZ
             F%B_2D%Z(ii,:)=-(F%PSIp(ii+1,:)-F%PSIp(ii-1,:))/ &
                  (F%X%R(ii+1)-F%X%R(ii-1))/F%X%R(ii)

          end do
          do ii=2_idef,F%dims(3)-1
             ! central difference over Z for interior nodes for BR
             F%B_2D%R(:,ii)=(F%PSIp(:,ii+1)-F%PSIp(:,ii-1))/ &
                  (F%X%Z(ii+1)-F%X%Z(ii-1))/F%X%R(:)
          end do

          do ii=1_idef,F%dims(1)             
             F%B_2D%PHI(ii,:)=-F%Bo*F%Ro/F%X%R(ii)
          end do

          F%E_2D%R=0._rp
          do ii=1_idef,F%dims(1)
             F%E_2D%PHI(ii,:)=F%Eo*F%Ro/F%X%R(ii)
          end do
          F%E_2D%Z=0._rp

          F%Bfield=.FALSE.
          
       end if

       if (params%mpi_params%rank.EQ.0) then

          if (F%axisymmetric_fields) then
             if (F%Bflux) then
                write(6,'("PSIp(r=0)",E17.10)') F%PSIp(F%dims(1)/2,F%dims(3)/2)
                write(6,'("BPHI(r=0)",E17.10)') F%Bo
                write(6,'("EPHI(r=0)",E17.10)') F%Eo
             else
                write(6,'("BR(r=0)",E17.10)') F%B_2D%R(F%dims(1)/2,F%dims(3)/2)
                write(6,'("BPHI(r=0)",E17.10)') &
                     F%B_2D%PHI(F%dims(1)/2,F%dims(3)/2)
                write(6,'("BZ(r=0)",E17.10)') F%B_2D%Z(F%dims(1)/2,F%dims(3)/2)
                write(6,'("EPHI(r=0)",E17.10)') &
                     F%E_2D%PHI(F%dims(1)/2,F%dims(3)/2)
             end if
          end if


       end if

       if (params%orbit_model(3:5).EQ.'pre') then
          if (params%mpi_params%rank.eq.0) then
             write(6,'("Initializing GC fields from external EM fields")')
          end if

          if (F%axisymmetric_fields) then
             call initialize_GC_fields(F)             
          else             
             call initialize_GC_fields_3D(F)
          end if
       end if

       !       write(6,'("gradBR",E17.10)') F%gradB_2D%R(F%dims(1)/2,F%dims(3)/2)
       !       write(6,'("gradBPHI",E17.10)') F%gradB_2D%PHI(F%dims(1)/2,F%dims(3)/2)
       !       write(6,'("gradBZ",E17.10)') F%gradB_2D%Z(F%dims(1)/2,F%dims(3)/2)

    CASE DEFAULT
    END SELECT

    if (params%mpi_params%rank.eq.0) then
       write(6,'("* * * * * * * * * * * * * * * * * * * * * * * * *",/)')
    end if
       
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
  
  subroutine initialize_GC_fields_3D(F)
    !! Computes the auxiliary fields \(\nabla|{\bf B}|\) and
    !! \(\nabla\times\hat{b}\) that are used in the RHS of the
    !! evolution equations for the GC orbit model.
    TYPE(FIELDS), INTENT(INOUT)      :: F
    !! An instance of the KORC derived type FIELDS.
    INTEGER                        :: ii,jj
    !! Iterator across F%dim
    REAL(rp), DIMENSION(:,:,:),ALLOCATABLE :: Bmag
    !! Magnetic field magnitude
    REAL(rp), DIMENSION(:,:,:,:),ALLOCATABLE :: bhat
    !! Magnetic field unit vector

    Bmag=SQRT(F%B_3D%R**2+F%B_3D%PHI**2+F%B_3D%Z**2)

    ALLOCATE(bhat(F%dims(1),F%dims(2),F%dims(3),3))

    bhat(:,:,:,1)=F%B_3D%R/Bmag
    bhat(:,:,:,2)=F%B_3D%PHI/Bmag
    bhat(:,:,:,3)=F%B_3D%Z/Bmag


    F%gradB_3D%PHI=0.
    ! No variation in phi direction

    ! Single-sided difference for axiliary fields at edge nodes
    ! Differential over R on first index, differential over Z
    ! on second index.

    ! gradB
    ! edge nodes at minimum R,Z
    F%gradB_3D%R(1,:,:)=(Bmag(2,:,:)-Bmag(1,:,:))/(F%X%R(2)-F%X%R(1))
    do ii=1_idef,F%dims(1)
       F%gradB_3D%PHI(ii,1,:)=(Bmag(ii,2,:)-Bmag(ii,F%dims(2),:))/ &
            (F%X%R(ii)*(F%X%PHI(2)-F%X%PHI(F%dims(2))))
    end do
    F%gradB_3D%Z(:,:,1)=(Bmag(:,:,2)-Bmag(:,:,1))/(F%X%Z(2)-F%X%Z(1))

    ! edge nodes at maximum R,Z
    F%gradB_3D%R(F%dims(1),:,:)=(Bmag(F%dims(1),:,:)-Bmag(F%dims(1)-1,:,:))/ &
         (F%X%R(F%dims(1))-F%X%R(F%dims(1)-1))
    do ii=1_idef,F%dims(1)
       F%gradB_3D%PHI(ii,F%dims(2),:)=(Bmag(ii,1,:)-Bmag(ii,F%dims(2)-1,:))/ &
            (F%X%R(ii)*(F%X%PHI(1)-F%X%PHI(F%dims(2)-1)))
    end do       
    F%gradB_3D%Z(:,:,F%dims(3))=(Bmag(:,:,F%dims(3))-Bmag(:,:,F%dims(3)-1))/ &
         (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1))

    ! curlb
    ! edge nodes at minimum R,PHI,Z
    ! R component has differential over PHI and Z
    do ii=1_idef,F%dims(1)
       F%curlb_3D%R(ii,1,:)=(bhat(ii,2,:,3)-bhat(ii,F%dims(2),:,3))/ &
            (F%X%R(ii)*(F%X%PHI(2)-F%X%PHI(F%dims(2))))
    end do
    F%curlb_3D%R(:,:,1)=F%curlb_3D%R(:,:,1)-&
         (bhat(:,:,2,2)-bhat(:,:,1,2))/(F%X%Z(2)-F%X%Z(1))

    ! PHI component has differentials over R and Z
    F%curlb_3D%PHI(1,:,:)=-(bhat(2,:,:,3)-bhat(1,:,:,3))/ &
         (F%X%R(2)-F%X%R(1))         

    F%curlb_3D%PHI(:,:,1)=F%curlb_3D%PHI(:,:,1)+ &
         ((bhat(:,:,2,1)-bhat(:,:,1,1))/(F%X%Z(2)-F%X%Z(1)))

    ! Z component has differentials over R and PHI
    F%curlb_3D%Z(1,:,:)=((bhat(2,:,:,2)*F%X%R(2)- &
         bhat(1,:,:,2)*F%X%R(1))/(F%X%R(2)-F%X%R(1)))/F%X%R(1)

    do ii=1_idef,F%dims(1)
       F%curlb_3D%Z(ii,1,:)=F%curlb_3D%Z(ii,1,:)-&
            (bhat(ii,2,:,1)-bhat(ii,F%dims(2),:,1))/ &
            (F%X%R(ii)*(F%X%PHI(2)-F%X%PHI(F%dims(2))))
    end do

    ! edge nodes at maximum R,Z
    ! R component has differential over PHI and Z
    do ii=1_idef,F%dims(1)
       F%curlb_3D%R(ii,F%dims(2),:)=(bhat(ii,1,:,3)-bhat(ii,F%dims(2)-1,:,3))/ &
            (F%X%R(ii)*(F%X%PHI(1)-F%X%PHI(F%dims(2)-1)))
    end do
    F%curlb_3D%R(:,:,F%dims(3))=F%curlb_3D%R(:,:,F%dims(3)) &
         -(bhat(:,:,F%dims(3),2)-bhat(:,:,F%dims(3)-1,2))/ &
         (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1))

    ! PHI component has differentials over R and Z
    F%curlb_3D%PHI(F%dims(1),:,:)=F%curlb_3D%PHI(F%dims(1),:,:)- &
         (bhat(F%dims(1),:,:,3)-bhat(F%dims(1)-1,:,:,3))/ &
         (F%X%R(F%dims(1))-F%X%R(F%dims(1)-1))         

    F%curlb_3D%PHI(:,:,F%dims(3))=F%curlb_3D%PHI(:,:,F%dims(3))+ &
         ((bhat(:,:,F%dims(3),1)-bhat(:,:,F%dims(3)-1,1))/ &
         (F%X%Z(F%dims(3))-F%X%Z(F%dims(3)-1)))

    ! Z component has differentials over R and PHI
    F%curlb_3D%Z(F%dims(1),:,:)=((bhat(F%dims(1),:,:,2)*F%X%R(F%dims(1))- &
         bhat(F%dims(1)-1,:,:,2)*F%X%R(F%dims(1)-1))/(F%X%R(F%dims(1))- &
         F%X%R(F%dims(1)-1)))/F%X%R(F%dims(1))

    do ii=1_idef,F%dims(1)
       F%curlb_3D%Z(ii,F%dims(2),:)=F%curlb_3D%Z(ii,F%dims(2),:)-&
            (bhat(ii,1,:,1)-bhat(ii,F%dims(2)-1,:,1))/ &
            (F%X%R(ii)*(F%X%PHI(1)-F%X%PHI(F%dims(2)-1)))
    end do
    
    do ii=2_idef,F%dims(1)-1
       ! central difference over R for interior nodes
       F%gradB_3D%R(ii,:,:)=(Bmag(ii+1,:,:)-Bmag(ii-1,:,:))/ &
            (F%X%R(ii+1)-F%X%R(ii-1))
       
       F%curlb_3D%Z(ii,:,:)=((bhat(ii+1,:,:,2)*F%X%R(ii+1)- &
            bhat(ii-1,:,:,2)*F%X%R(ii-1))/(F%X%R(ii+1)-F%X%R(ii-1)))/ &
            F%X%R(ii)
       F%curlb_3D%PHI(ii,:,:)=-(bhat(ii+1,:,:,3)-bhat(ii-1,:,:,3))/ &
            (F%X%R(ii+1)-F%X%R(ii-1))
    end do
    do ii=2_idef,F%dims(2)-1
       ! central difference over PHI for interior nodes       
       do jj=1_idef,F%dims(1)
          F%gradB_3D%PHI(:,ii,:)=(Bmag(:,ii+1,:)-Bmag(:,ii-1,:))/ &
               (F%X%R(jj)*(F%X%PHI(ii+1)-F%X%PHI(ii-1)))
          
          F%curlb_3D%Z(jj,ii,:)=F%curlb_3D%Z(jj,ii,:)-&
               (bhat(jj,ii+1,:,1)-bhat(jj,ii-1,:,1))/ &
               (F%X%R(jj)*(F%X%PHI(ii+1)-F%X%PHI(ii-1)))
          F%curlb_3D%R(jj,ii,:)=(bhat(jj,ii+1,:,3)-bhat(jj,ii-1,:,3))/ &
            (F%X%R(jj)*(F%X%PHI(ii+1)-F%X%PHI(ii-1)))
       end do
    end do
    do ii=2_idef,F%dims(3)-1
       ! central difference over Z for interior nodes
       F%gradB_3D%Z(:,:,ii)=(Bmag(:,:,ii+1)-Bmag(:,:,ii-1))/ &
            (F%X%Z(ii+1)-F%X%Z(ii-1))
       
       F%curlb_3D%R(:,:,ii)=F%curlb_3D%R(:,:,ii)- &
            (bhat(:,:,ii+1,2)-bhat(:,:,ii-1,2))/ &
            (F%X%Z(ii+1)-F%X%Z(ii-1))
       F%curlb_3D%PHI(:,:,ii)=F%curlb_3D%PHI(:,:,ii)+ &
            ((bhat(:,:,ii+1,1)-bhat(:,:,ii-1,1))/(F%X%Z(ii+1)-F%X%Z(ii-1)))
    end do

    DEALLOCATE(Bmag)
    DEALLOCATE(bhat) 

  end subroutine initialize_GC_fields_3D

  subroutine define_SC_time_step(params,F)
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    TYPE(FIELDS), INTENT(INOUT)         :: F

    F%subcycle_E_SC = FLOOR(F%dt_E_SC/params%dt,ip)

    params%t_it_SC = params%t_skip/F%subcycle_E_SC
    params%t_skip=F%subcycle_E_SC

    F%dt_E_SC=params%t_skip*params%dt

    if (params%mpi_params%rank.EQ.0) then

     write(6,'(/,"* * * * * SC_E1D SUBCYCLING * * * * *")')
     write(6,*) "SC_E1D sybcycling iterations: ",F%subcycle_E_SC
     write(6,*) "Updated number of outputs: ", &
          params%t_steps/(params%t_skip*params%t_it_SC)

     write(6,'("* * * * * * * * * * * * * * * * * * *",/)')
    end if
    

    
  end subroutine define_SC_time_step
  
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

    if (bfield.and.(.not.ALLOCATED(F%B_2D%R))) then
       call ALLOCATE_V_FIELD_2D(F%B_2D,F%dims)

       if(params%orbit_model(3:5).EQ.'pre') then
          call ALLOCATE_V_FIELD_2D(F%curlb_2D,F%dims)
          call ALLOCATE_V_FIELD_2D(F%gradB_2D,F%dims)
       end if

    end if

    if (bflux.and.(.not.ALLOCATED(F%PSIp))) then
       ALLOCATE(F%PSIp(F%dims(1),F%dims(3)))
    end if

    if (efield.and.(.not.ALLOCATED(F%E_2D%R))) then
       call ALLOCATE_V_FIELD_2D(F%E_2D,F%dims)
    end if

    if (.NOT.ALLOCATED(F%FLAG2D)) ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_2D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the variables keeping the 3-D fields data.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param[in] bfield Logical variable that specifies if the variables that keep the magnetic field data is allocated (bfield=T) or not (bfield=F).
  !! @param[in] efield Logical variable that specifies if the variables that keep the electric field data is allocated (efield=T) or not (efield=F).
  subroutine ALLOCATE_3D_FIELDS_ARRAYS(params,F,bfield,efield)
    TYPE (KORC_PARAMS), INTENT(IN) 	:: params
    TYPE(FIELDS), INTENT(INOUT)    :: F
    LOGICAL, INTENT(IN)            :: bfield
    LOGICAL, INTENT(IN)            :: efield

    if (bfield) then
       call ALLOCATE_V_FIELD_3D(F%B_3D,F%dims)

       if(params%orbit_model(3:5).EQ.'pre') then
          call ALLOCATE_V_FIELD_3D(F%curlb_3D,F%dims)
          call ALLOCATE_V_FIELD_3D(F%gradB_3D,F%dims)
       end if
       
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

    if (ALLOCATED(F%E_SC_1D%PHI)) DEALLOCATE(F%E_SC_1D%PHI)
    if (ALLOCATED(F%J1_SC_1D%PHI)) DEALLOCATE(F%J1_SC_1D%PHI)
    if (ALLOCATED(F%J2_SC_1D%PHI)) DEALLOCATE(F%J2_SC_1D%PHI)
    if (ALLOCATED(F%J3_SC_1D%PHI)) DEALLOCATE(F%J3_SC_1D%PHI)
    if (ALLOCATED(F%A1_SC_1D%PHI)) DEALLOCATE(F%A1_SC_1D%PHI)
    if (ALLOCATED(F%A2_SC_1D%PHI)) DEALLOCATE(F%A2_SC_1D%PHI)
    if (ALLOCATED(F%A3_SC_1D%PHI)) DEALLOCATE(F%A3_SC_1D%PHI)

    if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
    if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
    if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)

    if (ALLOCATED(F%FLAG2D)) DEALLOCATE(F%FLAG2D)
    if (ALLOCATED(F%FLAG3D)) DEALLOCATE(F%FLAG3D)
  end subroutine DEALLOCATE_FIELDS_ARRAYS
end module korc_fields