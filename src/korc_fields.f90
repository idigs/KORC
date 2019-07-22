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
       get_fields_p,&
       analytical_fields_GC_p,&
       analytical_fields_Bmag_p,&
       analytical_fields_p,&
       initialize_fields,&
       load_field_data_from_hdf5,&
       load_dim_data_from_hdf5,&
       ALLOCATE_2D_FIELDS_ARRAYS,&
       ALLOCATE_3D_FIELDS_ARRAYS,&
       DEALLOCATE_FIELDS_ARRAYS
  PRIVATE :: get_analytical_fields,&
       get_analytical_fields_p,&
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

  subroutine analytical_fields_GC_init(params,F,Y,E,B,gradB,curlb,flag)
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

!    write(6,'("Y: ",E17.10)') Y
    
    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag)
    do pp=1_idef,ss
!       if ( flag(pp) .EQ. 1_is ) then

          rm=sqrt((Y(pp,1)-F%AB%Ro)**2+Y(pp,3)**2)
          qprof = 1.0_rp + (rm/F%AB%lambda)**2
          
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

    
    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag)
    do pp=1_idef,ss
!       if ( flag(pp) .EQ. 1_is ) then

          rm=sqrt((Y(pp,1)-F%AB%Ro)**2+Y(pp,3)**2)
          qprof = 1.0_rp + (rm/F%AB%lambda)**2
          
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
  end subroutine analytical_fields_GC

subroutine analytical_fields_Bmag_p(R0,B0,lam,q0,EF0,Y_R,Y_PHI,Y_Z,Bmag,E_PHI)

    REAL(rp),INTENT(IN)  :: R0,B0,lam,q0,EF0
    REAL(rp),DIMENSION(8),INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8) :: B_R,B_PHI,B_Z,rm,qprof
    REAL(rp),DIMENSION(8),INTENT(OUT) :: Bmag,E_PHI
    integer(ip) :: cc
    

    !$OMP SIMD 
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

subroutine analytical_fields_GC_p(R0,B0,lam,q0,E0,Y_R,Y_PHI, &
       Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
       gradB_PHI,gradB_Z)

    REAL(rp),INTENT(IN)  :: R0,B0,lam,q0,E0
    REAL(rp),DIMENSION(8),INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(8),INTENT(OUT) :: curlB_R,curlB_PHI,curlB_Z,E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(8)  :: dRBR,dRBPHI,dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI
    REAL(rp),DIMENSION(8)  :: dRbhatZ,dZbhatR,dZbhatPHI,qprof,rm
    integer(ip) :: cc
    

    !$OMP SIMD
    do cc=1_idef,8_idef
       rm(cc)=sqrt((Y_R(cc)-R0)*(Y_R(cc)-R0)+Y_Z(cc)*Y_Z(cc))
       qprof(cc) = 1.0_rp + (rm(cc)*rm(cc)/(lam*lam))

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
               vars%gradB,vars%curlb, vars%flag)

       else

          call cyl_check_if_confined(F,vars%Y,vars%flag)
       
          call analytical_fields_GC(params,F,vars%Y, vars%E, vars%B, &
               vars%gradB,vars%curlb, vars%flag)

    end if

    endif
    
  end subroutine get_analytical_fields

  elemental subroutine get_analytical_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI, &
       B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,R0,B0,lam,q0,E0)
    !! @note Interface for calculating the analytical electric and magnetic
    !! fields for each particle in the simulation. @endnote
    REAL(rp), INTENT(in)  :: R0,B0,E0,q0,lam
    REAL(rp),  INTENT(INOUT)   :: Y_R
    REAL(rp),  INTENT(INOUT)   :: Y_PHI
    REAL(rp),  INTENT(INOUT)   :: Y_Z
    REAL(rp),  INTENT(INOUT)   :: B_R
    REAL(rp),  INTENT(INOUT)   :: B_PHI
    REAL(rp),  INTENT(INOUT)   :: B_Z
    REAL(rp),  INTENT(INOUT)   :: gradB_R
    REAL(rp),  INTENT(INOUT)   :: gradB_PHI
    REAL(rp),  INTENT(INOUT)   :: gradB_Z
    REAL(rp),  INTENT(INOUT)   :: curlB_R
    REAL(rp),  INTENT(INOUT)   :: curlB_PHI
    REAL(rp),  INTENT(INOUT)   :: curlB_Z
    REAL(rp),  INTENT(INOUT)   :: E_R
    REAL(rp),  INTENT(INOUT)   :: E_PHI
    REAL(rp),  INTENT(INOUT)   :: E_Z


!    if (params%orbit_model(1:2).eq.'FO') then

!       call cart_to_tor_check_if_confined(vars%X,F,vars%Y,vars%flag)
       
!       call analytical_fields(F,vars%Y, vars%E, vars%B, vars%flag)

!    elseif (params%orbit_model(1:2).eq.'GC') then

!          call cyl_check_if_confined_p(F,Y_R,Y_Z,flag)
       
!          call analytical_fields_GC_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
!               curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,R0,B0,lam,q0,E0)

!    end if

    
  end subroutine get_analytical_fields_p

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
    ALLOCATE( vars%E(ppp,3) )
    ALLOCATE( vars%flag(ppp) )

    vars%X = Xo
    vars%flag = 1_idef

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

  subroutine get_fields_p(params,R0,B0,lam,q0,E0,Y_R,Y_PHI, &
                  Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flag_cache)
    !! @note Inferface with calls to subroutines for calculating the electric 
    !! and magnetic field for each particle in the simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    REAL(rp),INTENT(in)  :: R0,B0,E0,q0,lam
    REAL(rp),DIMENSION(8),INTENT(INOUT)   :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT)   :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT)   :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT)   :: curlB_R,curlB_PHI,curlB_Z
    REAL(rp),DIMENSION(8),INTENT(INOUT)   :: E_R,E_PHI,E_Z
    INTEGER(is),DIMENSION(8),INTENT(INOUT)   :: flag_cache

!    SELECT CASE (TRIM(params%plasma_model))
!    CASE('ANALYTICAL')
       if (params%field_eval.eq.'eqn') then
          call analytical_fields_GC_p(R0,B0,lam,q0,E0,Y_R,Y_PHI, &
                  Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z)
       else
          call interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
               curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)
       end if       
!    CASE('EXTERNAL')

!       call interp_fields(params,vars, F)
!       if (F%Efield.AND..NOT.F%Efield_in_file) then
!          call analytical_electric_field_cyl(F,vars%Y,vars%E,vars%flag)
!       end if
!    CASE('UNIFORM')
       
!       call uniform_fields(vars, F)
!    CASE DEFAULT
!    END SELECT
  end subroutine get_fields_p


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

    if (params%mpi_params%rank .EQ. 0) then
       write(6,'(/,"* * * * * * * * INITIALIZING FIELDS * * * * * * * *")')
    end if
       
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
          F%Bflux=.FALSE.
          F%Efield=.TRUE.

          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,.TRUE.,.FALSE.,.TRUE.)

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
                !! Sign convention in analytical fields corresponds to
                !! DIII-D fields with \(B_\phi<0\) and \(B_\theta<0\).
                F%FLAG2D=1.
             end do
          end do         

          if (params%orbit_model(3:5).eq.'pre') then
             write(6,'("Initializing GC fields from analytic EM fields")')
             call initialize_GC_fields(F)
          end if

                    
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
          
          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
               F%Bflux,F%Efield.AND. &
               F%Efield_in_file)
          
       else
          call ALLOCATE_3D_FIELDS_ARRAYS(F,F%Bfield,F%Efield)
       end if
       !allocates 2D or 3D data arrays (fields and spatial)
       
       call load_field_data_from_hdf5(params,F)
      
       if (F%Bflux.and.(params%orbit_model(1:2).EQ.'GC')) then

          call initialize_fields_interpolant(params,F)

          F%Bfield=.TRUE.
          F%Bflux=.FALSE.
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
             F%B_2D%PHI(ii,:)=F%Bo*F%Ro/F%X%R(ii)
          end do

          F%E_2D%R=0._rp
          F%E_2D%PHI=0._rp
          F%E_2D%Z=0._rp
          
       end if

       write(6,'("BR",E17.10)') F%B_2D%R(F%dims(1)/2,F%dims(3)/2)
       write(6,'("BPHI",E17.10)') F%B_2D%PHI(F%dims(1)/2,F%dims(3)/2)
       write(6,'("BZ",E17.10)') F%B_2D%Z(F%dims(1)/2,F%dims(3)/2)


       
       if (params%orbit_model(3:5).EQ.'pre') then
          write(6,'("Initializing GC fields from external EM fields")')
          call initialize_GC_fields(F)
       end if

       write(6,'("gradBR",E17.10)') F%gradB_2D%R(F%dims(1)/2,F%dims(3)/2)
       write(6,'("gradBPHI",E17.10)') F%gradB_2D%PHI(F%dims(1)/2,F%dims(3)/2)
       write(6,'("gradBZ",E17.10)') F%gradB_2D%Z(F%dims(1)/2,F%dims(3)/2)
       
    CASE DEFAULT
    END SELECT

    write(6,'("* * * * * * * * * * * * ** * * * * * * * * * * *",/)')

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
