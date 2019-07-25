module korc_interp
  !! @note Module containing functions and subroutines for performing 
  !! interpolations using the PSPLINE library. @endnote
  !! For a detailed documentation of the PSPLINE library we refer the 
  !! user to "https://w3.pppl.gov/ntcc/PSPLINE/".
  use korc_types
  use korc_coords
  use korc_rnd_numbers
  use korc_hpc

  use EZspline_obj	! psplines module
  use EZspline		! psplines module

  !$ use OMP_LIB

  IMPLICIT NONE


#ifdef DOUBLE_PRECISION


  TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for 
     !! cylindrical components of vector fields
     !! \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi + 
     !! F_Z\hat{e}_Z\). Real precision of 8 bytes. @endnote
     TYPE(EZspline3_r8)    :: R		
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3_r8)    :: PHI	
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3_r8)    :: Z		
     !! Interpolant of \(F_Z(R,\phi,Z)\).
     
     INTEGER               :: NR 
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI 
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ 
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /) 
     !! Not-a-knot boundary condition for the interpolants at both 
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /) 
     !! Periodic boundary condition for the interpolants at both 
     !! ends of the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /) 
     !! Not-a-knot boundary condition for the interpolants at both 
     !! ends of the \(Z\) direction.
  END TYPE KORC_3D_FIELDS_INTERPOLANT



  TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
     !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2_r8)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline2_r8)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline2_r8)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline2_r8)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).
     
     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2D_FIELDS_INTERPOLANT


  TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,\phi,Z)\),
     !! temperature \(T_e(R,\phi,Z)\), and effective charge number
     !! \(Z_{eff}(R,\phi,Z)\) profiles. Real precision of 8 bytes. @endnote
     TYPE(EZspline3_r8)    :: ne
     !! Interpolant of background electron density \(n_e(R,\phi,Z)\).
     TYPE(EZspline3_r8)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,\phi,Z)\).
     TYPE(EZspline3_r8)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_3D_PROFILES_INTERPOLANT


  TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,Z)\), temperature \(T_e(R,Z)\), and
     !! effective charge number \(Z_{eff}(R,Z)\) profiles.
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2_r8)    :: ne
     !! Interpolant of background electron density \(n_e(R,Z)\).
     TYPE(EZspline2_r8)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,Z)\).
     TYPE(EZspline2_r8)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2D_PROFILES_INTERPOLANT

  
#elif SINGLE_PRECISION

  
  TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of vector fields \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R +
     !! F_\phi\hat{e}_phi + F_Z\hat{e}_Z\).
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline3_r4)    :: R
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: PHI
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Z
     !! Interpolant of \(F_Z(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends of
     !! the \(Z\) direction.
  END TYPE KORC_3D_FIELDS_INTERPOLANT


  TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for cylindrical
     !! components of vector fields \(\mathbf{F}(R,Z) = F_R\hat{e}_R +
     !! F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline2_r4)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline2_r4)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline2_r4)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline2_r4)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_2D_FIELDS_INTERPOLANT



  TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,\phi,Z)\),
     !! temperature \(T_e(R,\phi,Z)\), and effective charge number
     !! \(Z_{eff}(R,\phi,Z)\) profiles.
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline3_r4)    :: ne
     !! Interpolant of background electron density \(n_e(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends of
     !! the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_3D_PROFILES_INTERPOLANT



  TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of the density \(n_e(R,Z)\),
     !! temperature \(T_e(R,Z)\), and effective charge number \(Z_{eff}(R,Z)\) profiles.
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2_r4)    :: ne
     !! Interpolant of background electron density \(n_e(R,Z)\).
     TYPE(EZspline2_r4)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,Z)\).
     TYPE(EZspline2_r4)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_2D_PROFILES_INTERPOLANT


#endif



  TYPE, PRIVATE :: KORC_INTERPOLANT_DOMAIN
     !! @note Derived type containing 2-D and 3-D arrays with the information of
     !! the spatial domain where the fields and profiles are known.
     !! This info is used for detecting when a particle is lost, and therefore not
     !! followed anymore. @endnote
     INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE      :: FLAG2D
     !! 2-D array with info of the spatial domain where the axisymmetric fields
     !! and plasma profiles are known.
     INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE    :: FLAG3D
     !! 3-D array with info of the spatial domain where the 3-D fields and plasma
     !! profiles are known.

     REAL(rp)                                          :: Ro
     !! Smaller radial position of the fields and profiles domain.
     REAL(rp)                                          :: Zo
     !! Smaller vertical position of the fields and profiles domain.

     REAL(rp)                                          :: DR
     !! Separation between grid points along the radial direction.
     REAL(rp)                                          :: DPHI  !
     ! Separation between grid points along the azimuthal direction.
     REAL(rp)                                          :: DZ
     !! Separation between grid points along the vertical direction.
  END TYPE KORC_INTERPOLANT_DOMAIN


  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_3d
  !! An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: efield_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the electric field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: efield_3d
  !! An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating
  !! the electric field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: gradB_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: curlb_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: fields_domain
  !! An instance of KORC_INTERPOLANT_DOMAIN used for interpolating fields.
  TYPE(KORC_2D_PROFILES_INTERPOLANT), PRIVATE    :: profiles_2d
  !! An instance of KORC_2D_PROFILES_INTERPOLANT for interpolating plasma
  !! profiles.
  TYPE(KORC_3D_PROFILES_INTERPOLANT), PRIVATE    :: profiles_3d
  !! An instance of KORC_3D_PROFILES_INTERPOLANT for interpolating plasma
  !! profiles.
  TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: profiles_domain
  !! An instance of KORC_INTERPOLANT_DOMAIN used for interpolating plasma
  !! profiles.
  INTEGER                                        :: ezerr
  !! Error status during PSPLINE interpolations.


  PUBLIC :: interp_fields,&
       interp_fields_p,&
       interp_FOfields_p,&
       interp_bmag_p,&
       interp_collision_p,&
!       interp_fields_FO_p,&
       interp_profiles,&
       initialize_fields_interpolant,&
       initialize_profiles_interpolant,&
       finalize_interpolants,&
       calculate_initial_magnetic_field
  PRIVATE :: interp_3D_bfields,&
       interp_2D_bfields,&
       interp_3D_efields,&
       interp_2D_efields,&
       interp_2D_profiles,&
       interp_3D_profiles,&
       check_if_in_fields_domain,&
       check_if_in_profiles_domain,&
       check_if_in_profiles_domain_p,&
       check_if_in_fields_domain_p,&
       interp_2D_gradBfields,&
       interp_2D_curlbfields,&
       gradient_2D_Bfields

CONTAINS


  subroutine initialize_fields_interpolant(params,F)
    !! @note Subroutine that initializes fields interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of fields in the KORC-dervied-type variable F.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation.
    !! See [[korc_types]] and [[korc_fields]].

    if ((params%field_model .EQ. 'EXTERNAL').or. &
         (params%field_eval.eq.'interp')) then

       if (params%mpi_params%rank .EQ. 0) then
          write(6,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')
       end if

       ! * * * * * * * * MAGNETIC FIELD * * * * * * * * !
       if (F%Bflux) then
          
          bfield_2d%NR = F%dims(1)
          bfield_2d%NZ = F%dims(3)

          ! Initializing poloidal flux interpolant
          call EZspline_init(bfield_2d%A,bfield_2d%NR,bfield_2d%NZ, &
               bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
          call EZspline_error(ezerr)

          bfield_2d%A%x1 = F%X%R
          bfield_2d%A%x2 = F%X%Z

!          write(6,'("R",E17.10)') F%X%R
!          write(6,'("Z",E17.10)') F%X%Z

          call EZspline_setup(bfield_2d%A, F%PSIp, ezerr, .TRUE.)
          call EZspline_error(ezerr)

!          write(6,'("PSIp",E17.10)') F%PSIp

          ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))
          fields_domain%FLAG2D = F%FLAG2D

          fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
          fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
       end if

       if (F%Bfield) then
          if (F%axisymmetric_fields) then
             bfield_2d%NR = F%dims(1)
             bfield_2d%NZ = F%dims(3)

             ! Initializing R component
             call EZspline_init(bfield_2d%R,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%R%x1 = F%X%R
             bfield_2d%R%x2 = F%X%Z

             call EZspline_setup(bfield_2d%R, F%B_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(bfield_2d%PHI,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%PHI%x1 = F%X%R
             bfield_2d%PHI%x2 = F%X%Z

             call EZspline_setup(bfield_2d%PHI, F%B_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(bfield_2d%Z,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%Z%x1 = F%X%R
             bfield_2d%Z%x2 = F%X%Z

             call EZspline_setup(bfield_2d%Z, F%B_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             if (params%orbit_model.eq.'GCpre') then
                gradB_2d%NR = F%dims(1)
                gradB_2d%NZ = F%dims(3)

                ! Initializing GRADBR component
                call EZspline_init(gradB_2d%R,gradB_2d%NR,gradB_2d%NZ, &
                     gradB_2d%BCSR,gradB_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_2d%R%x1 = F%X%R
                gradB_2d%R%x2 = F%X%Z

                call EZspline_setup(gradB_2d%R, F%gradB_2D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing GRADBPHI component
                call EZspline_init(gradB_2d%PHI,gradB_2d%NR,gradB_2d%NZ, &
                     gradB_2d%BCSR,gradB_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_2d%PHI%x1 = F%X%R
                gradB_2d%PHI%x2 = F%X%Z

                call EZspline_setup(gradB_2d%PHI, F%gradB_2D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing GRADBZ component
                call EZspline_init(gradB_2d%Z,gradB_2d%NR,gradB_2d%NZ, &
                     gradB_2d%BCSR,gradB_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_2d%Z%x1 = F%X%R
                gradB_2d%Z%x2 = F%X%Z

                call EZspline_setup(gradB_2d%Z, F%gradB_2D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                curlb_2d%NR = F%dims(1)
                curlb_2d%NZ = F%dims(3)

                ! Initializing CURLBR component
                call EZspline_init(curlb_2d%R,curlb_2d%NR,curlb_2d%NZ, &
                     curlb_2d%BCSR,curlb_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_2d%R%x1 = F%X%R
                curlb_2d%R%x2 = F%X%Z

                call EZspline_setup(curlb_2d%R, F%curlb_2D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing CURLBPHI component
                call EZspline_init(curlb_2d%PHI,curlb_2d%NR,curlb_2d%NZ, &
                     curlb_2d%BCSR,curlb_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_2d%PHI%x1 = F%X%R
                curlb_2d%PHI%x2 = F%X%Z

                call EZspline_setup(curlb_2d%PHI, F%curlb_2D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing CURLBZ component
                call EZspline_init(curlb_2d%Z,curlb_2d%NR,curlb_2d%NZ, &
                     curlb_2d%BCSR,curlb_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_2d%Z%x1 = F%X%R
                curlb_2d%Z%x2 = F%X%Z

                call EZspline_setup(curlb_2d%Z, F%curlb_2D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

             end if

             if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                  ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))
             fields_domain%FLAG2D = F%FLAG2D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          else 
             bfield_3d%NR = F%dims(1)
             bfield_3d%NPHI = F%dims(2)
             bfield_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(bfield_3d%R, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%R%x1 = F%X%R
             ! bfield_3d%R%x2 = F%X%PHI
             bfield_3d%R%x3 = F%X%Z

             call EZspline_setup(bfield_3d%R, F%B_3D%R, ezerr)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(bfield_3d%PHI, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%PHI%x1 = F%X%R
             ! bfield_3d%PHI%x2 = F%X%PHI
             bfield_3d%PHI%x3 = F%X%Z

             call EZspline_setup(bfield_3d%PHI, F%B_3D%PHI, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Z component of interpolant
             call EZspline_init(bfield_3d%Z, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%Z%x1 = F%X%R
             ! bfield_3d%Z%x2 = F%X%PHI
             bfield_3d%Z%x3 = F%X%Z

             call EZspline_setup(bfield_3d%Z, F%B_3D%Z, ezerr)
             call EZspline_error(ezerr)

             ALLOCATE(fields_domain%FLAG3D(bfield_3d%NR,bfield_3d%NPHI, &
                  bfield_3d%NZ))
             fields_domain%FLAG3D = F%FLAG3D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DPHI = 2.0_rp*C_PI/bfield_3d%NPHI
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          end if
       end if
       
       fields_domain%Ro = F%X%R(1)
       fields_domain%Zo = F%X%Z(1)

       ! * * * * * * * * ELECTRIC FIELD * * * * * * * * !
       if (F%Efield.AND.(F%Efield_in_file.OR. &
            params%field_eval.eq.'interp')) then
          if (F%axisymmetric_fields) then
             efield_2d%NR = F%dims(1)
             efield_2d%NZ = F%dims(3)

             ! Initializing R component
             call EZspline_init(efield_2d%R,efield_2d%NR,efield_2d%NZ, &
                  efield_2d%BCSR,efield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             efield_2d%R%x1 = F%X%R
             efield_2d%R%x2 = F%X%Z

             call EZspline_setup(efield_2d%R, F%E_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(efield_2d%PHI,efield_2d%NR,efield_2d%NZ, &
                  efield_2d%BCSR,efield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             efield_2d%PHI%x1 = F%X%R
             efield_2d%PHI%x2 = F%X%Z

             call EZspline_setup(efield_2d%PHI, F%E_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(efield_2d%Z,efield_2d%NR,efield_2d%NZ, &
                  efield_2d%BCSR,efield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             efield_2d%Z%x1 = F%X%R
             efield_2d%Z%x2 = F%X%Z

             call EZspline_setup(efield_2d%Z, F%E_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)
          else
             efield_3d%NR = F%dims(1)
             efield_3d%NPHI = F%dims(2)
             efield_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(efield_3d%R, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
                  efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             efield_3d%R%x1 = F%X%R
             ! efield_3d%R%x2 = F%X%PHI
             efield_3d%R%x3 = F%X%Z

             call EZspline_setup(efield_3d%R, F%E_3D%R, ezerr)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(efield_3d%PHI, efield_3d%NR, efield_3d%NPHI, &
                  efield_3d%NZ,efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             efield_3d%PHI%x1 = F%X%R
             ! efield_3d%PHI%x2 = F%X%PHI
             efield_3d%PHI%x3 = F%X%Z

             call EZspline_setup(efield_3d%PHI, F%E_3D%PHI, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Z component of interpolant
             call EZspline_init(efield_3d%Z, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
                  efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             efield_3d%Z%x1 = F%X%R
             ! efield_3d%Z%x2 = F%X%PHI
             efield_3d%Z%x3 = F%X%Z

             call EZspline_setup(efield_3d%Z, F%E_3D%Z, ezerr)
             call EZspline_error(ezerr)
          end if
       end if

       if (params%mpi_params%rank .EQ. 0) then
          write(6,'("* * * * * * INTERPOLANT INITIALIZED * * * * * *",/)')
       end if
    else if (params%field_model .EQ. 'ANALYTICAL') then
       if (params%mpi_params%rank .EQ. 0) then
          write(6,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *",/)')
       end if
    else if (params%field_model .EQ. 'UNIFORM') then
       if (params%mpi_params%rank .EQ. 0) then
          write(6,'("* * * * USING UNIFORM MAGNETIC FIELD * * * *",/)')
       end if
    end if
  end subroutine initialize_fields_interpolant


  subroutine check_if_in_fields_domain(Y,flag)
    !! @note Subrotuine that checks if particles in the simulation are within
    !! the spatial domain where interpolants and fields are known. @endnote
    !! External fields and interpolants can have different spatial domains where
    !! they are defined. Therefore, it is necessary to
    !! check if a given particle has left these spatial domains to stop following
    !! it, otherwise this will cause an error in the simulation.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Particles' position in cylindrical coordinates,
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed in the
    !! simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Species iterator.
    
    ss = SIZE(Y,1)

!    write(6,'("R: ",E15.10)') Y(1,1)
!    write(6,'("PHI: ",E15.10)') Y(2,1)
!    write(6,'("Z: ",E15.10)') Y(3,1)
    
    if (ALLOCATED(fields_domain%FLAG3D)) then
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) &
       !$OMP& SHARED(Y,flag,fields_domain,bfield_3d)
       do pp=1_idef,ss
          IR = INT(FLOOR((Y(pp,1)  - fields_domain%Ro + 0.5_rp* &
               fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y(pp,2)  + 0.5_rp*fields_domain%DPHI)/ &
               fields_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(fields_domain%Zo) + 0.5_rp* &
               fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

          if ((fields_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
               ((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) &
       !$OMP& SHARED(Y,flag,fields_domain,bfield_2d)
       do pp=1_idef,ss
          IR = INT(FLOOR((Y(pp,1)  - fields_domain%Ro + 0.5_rp* &
               fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(fields_domain%Zo) + 0.5_rp* &
               fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

!          write(6,'("IR: ",I10)') IR
!          write(6,'("IZ: ",I10)') IZ
          
          if ((fields_domain%FLAG2D(IR,IZ).NE.1_is).OR.((IR.GT. &
               bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine check_if_in_fields_domain
  
  subroutine check_if_in_fields_domain_p(Y_R,Y_PHI,Y_Z,flag)
    !! @note Subrotuine that checks if particles in the simulation are within
    !! the spatial domain where interpolants and fields are known. @endnote
    !! External fields and interpolants can have different spatial domains where
    !! they are defined. Therefore, it is necessary to
    !! check if a given particle has left these spatial domains to stop following
    !! it, otherwise this will cause an error in the simulation.
    REAL(rp), DIMENSION(8),  INTENT(IN)      :: Y_R,Y_PHI,Y_Z
    INTEGER(is), DIMENSION(8), INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed in the
    !! simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Species iterator.
    

!    write(6,'("YR:",E17.10)') Y_R
!    write(6,'("YPHI:",E17.10)') Y_PHI
!    write(6,'("YZ:",E17.10)') Y_Z

!    write(6,'("Ro:",E17.10)') fields_domain%Ro
!    write(6,'("Zo:",E17.10)') fields_domain%Zo
!    write(6,'("DR:",E17.10)') fields_domain%DR
!    write(6,'("DZ:",E17.10)') fields_domain%DZ
    
    if (ALLOCATED(fields_domain%FLAG3D)) then
       !$OMP SIMD
!       !$OMP&  aligned(IR,IPHI,IZ)
       do pp=1_idef,8_idef
          IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
               0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y_PHI(pp)  + 0.5_rp*fields_domain%DPHI)/ &
               fields_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
               0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

          if ((fields_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
               ((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END SIMD
    else
       !$OMP SIMD
!       !$OMP& aligned(IR,IZ)
       do pp=1_idef,8_idef
          IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
               0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
               0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

!          write(6,'("IR: ",I16)') IR
!          write(6,'("IZ: ",I16)') IZ
!          write(6,'("Size of fields_domain R: ",I16)') &
!               size(fields_domain%FLAG2D,1)
!          write(6,'("Size of fields_domain Z: ",I16)') &
!               size(fields_domain%FLAG2D,2)   
          
          if ((fields_domain%FLAG2D(IR,IZ).NE.1_is).OR. &
               ((IR.GT.bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ))) then
             flag(pp) = 0_is

!             write(6,'("Shit''s fucked.")')
          end if
       end do      
       !$OMP END SIMD
!       write(6,'("Shit''s not fucked.")')
    end if
  end subroutine check_if_in_fields_domain_p


  subroutine initialize_profiles_interpolant(params,P)
    !! @note Subroutine that initializes plasma profiles interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of plasma profiles in the KORC-dervied-type variable P.
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(INOUT) :: P
    !! An instance of KORC's derived type PROFILES containing
    !! all the information about the plasma profiles used in the simulation.
    !! See [[korc_types]] and [[korc_profiles]].
    
    if (params%collisions) then
       if ((params%profile_model .EQ. 'EXTERNAL').or. &
            (params%field_eval.eq.'interp')) then
          
          if (params%mpi_params%rank .EQ. 0) then
             write(6,'("* * * * INITIALIZING PROFILES INTERPOLANT * * * *")')
          end if

          if (P%axisymmetric) then
             profiles_2d%NR = P%dims(1)
             profiles_2d%NZ = P%dims(3)

!             write(6,'("NR",I15)') profiles_2d%NR
!             write(6,'("NZ",I15)') profiles_2d%NR
             
             ! Initializing ne
             !			call EZspline_init(profiles_2d%ne,profiles_2d%NR, &
             !profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_init(profiles_2d%ne,profiles_2d%NR,profiles_2d%NZ, &
                  profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             profiles_2d%ne%x1 = P%X%R
             profiles_2d%ne%x2 = P%X%Z

             call EZspline_setup(profiles_2d%ne, P%ne_2D, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             
             ! Initializing Te
             !			call EZspline_init(profiles_2d%Te,profiles_2d%NR, &
             !profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_init(profiles_2d%Te,profiles_2d%NR,profiles_2d%NZ, &
                  profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             profiles_2d%Te%x1 = P%X%R
             profiles_2d%Te%x2 = P%X%Z

!             write(6,'("Te_interp_R",E17.10)') profiles_2d%Te%x1
!             write(6,'("Te_interp_Z",E17.10)') profiles_2d%Te%x2

!             write(6,'("Te",E17.10)') P%Te_2D(10,:)
             
             call EZspline_setup(profiles_2d%Te, P%Te_2D, ezerr, .TRUE.)
             call EZspline_error(ezerr)
             
             ! Initializing Zeff
             !			call EZspline_init(profiles_2d%Zeff, &
             !profiles_2d%NR,profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_init(profiles_2d%Zeff,profiles_2d%NR, &
                  profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             profiles_2d%Zeff%x1 = P%X%R
             profiles_2d%Zeff%x2 = P%X%Z
             
             call EZspline_setup(profiles_2d%Zeff, P%Zeff_2D, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             
             ALLOCATE(profiles_domain%FLAG2D(profiles_2d%NR,profiles_2d%NZ))
             profiles_domain%FLAG2D = P%FLAG2D

             profiles_domain%DR = ABS(P%X%R(2) - P%X%R(1))
             profiles_domain%DZ = ABS(P%X%Z(2) - P%X%Z(1))
          else
             profiles_3d%NR = P%dims(1)
             profiles_3d%NPHI = P%dims(2)
             profiles_3d%NZ = P%dims(3)

             ! Initializing ne
             call EZspline_init(profiles_3d%ne, profiles_3d%NR, &
                  profiles_3d%NPHI, profiles_3d%NZ,&
                  profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             profiles_3d%ne%x1 = P%X%R
             ! profiles_3d%ne%x2 = P%X%PHI
             profiles_3d%ne%x3 = P%X%Z

             call EZspline_setup(profiles_3d%ne, P%ne_3D, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Te
             call EZspline_init(profiles_3d%Te, profiles_3d%NR, &
                  profiles_3d%NPHI, profiles_3d%NZ,&
                  profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             profiles_3d%Te%x1 = P%X%R
             ! profiles_3d%Te%x2 = P%X%PHI
             profiles_3d%Te%x3 = P%X%Z

             call EZspline_setup(profiles_3d%Te, P%Te_3D, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Zeff
             call EZspline_init(profiles_3d%Zeff, profiles_3d%NR, &
                  profiles_3d%NPHI, profiles_3d%NZ,&
                  profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             profiles_3d%Zeff%x1 = P%X%R
             ! profiles_3d%Zeff%x2 = P%X%PHI
             profiles_3d%Zeff%x3 = P%X%Z

             call EZspline_setup(profiles_3d%Zeff, P%Zeff_3D, ezerr)
             call EZspline_error(ezerr)

             ALLOCATE(profiles_domain%FLAG3D(profiles_3d%NR,profiles_3d%NPHI, &
                  profiles_3d%NZ))
             profiles_domain%FLAG3D = P%FLAG3D

             profiles_domain%DR = ABS(P%X%R(2) - P%X%R(1))
             profiles_domain%DPHI = 2.0_rp*C_PI/profiles_3d%NPHI
             profiles_domain%DZ = ABS(P%X%Z(2) - P%X%Z(1))
          end if

          profiles_domain%Ro = P%X%R(1)
          profiles_domain%Zo = P%X%Z(1)

          if (params%mpi_params%rank .EQ. 0) then
             write(6,'("* * * * * * INTERPOLANT   INITIALIZED * * * * * *")')
          end if
       else if (params%profile_model .EQ. 'ANALYTICAL') then
          if (params%mpi_params%rank .EQ. 0) then
             write(6,'("* * * * USING ANALYTICAL PROFILES * * * *")')
          end if
       else if (params%profile_model .EQ. 'UNIFORM') then
          if (params%mpi_params%rank .EQ. 0) then
             write(6,'("* * * * UNIFORM PLASMA: NO PROFILES USED * * * *")')
          end if
       end if
    end if
  end subroutine initialize_profiles_interpolant


  subroutine check_if_in_profiles_domain(Y,flag)
    !! @note Subrotuine that checks if particles in the simulation are
    !! within the spatial domain where interpolants and plasma profiles
    !! are known. @endnote
    !!External plasma profiles and interpolants can have different spatial
    !! domains where they are defined. Therefore, it is necessary to check
    !! if a given particle has left these spatial domains to stop following
    !! it, otherwise this will cause an error in the simulation.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Particles' position in cylindrical coordinates,
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed
    !! in the simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! @param IR Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! @param IPHI Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! @param IZ Variable used to localize the grid cell in the
    !! \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! @param pp Particle iterator.
    INTEGER(ip)                                            :: ss
    !! @param ss Species iterator.
    
    ss = SIZE(Y,1)

    if (ALLOCATED(profiles_domain%FLAG3D)) then
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) &
       !$OMP& SHARED(Y,flag,profiles_domain,profiles_3d)
       do pp=1_idef,ss
          IR = INT(FLOOR((Y(pp,1)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y(pp,2)  + 0.5_rp*profiles_domain%DPHI)/ &
               profiles_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_3d%NR).OR.(IZ.GT.profiles_3d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) &
       !$OMP& SHARED(Y,flag,profiles_domain,profiles_2d)
       do pp=1_idef,ss
          IR = INT(FLOOR((Y(pp,1)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG2D(IR,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_2d%NR).OR.(IZ.GT.profiles_2d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine check_if_in_profiles_domain

  subroutine check_if_in_profiles_domain_p(Y_R,Y_PHI,Y_Z,flag)
    REAL(rp), DIMENSION(8),  INTENT(IN)      :: Y_R,Y_PHI,Y_Z

    INTEGER(is), DIMENSION(8), INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed
    !! in the simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! @param IR Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! @param IPHI Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! @param IZ Variable used to localize the grid cell in the
    !! \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! @param pp Particle iterator.
    INTEGER(ip)                                            :: ss
    !! @param ss Species iterator.
    

    if (ALLOCATED(profiles_domain%FLAG3D)) then
       !$OMP SIMD
!       !$OMP& aligned(IR,IPHI,IZ)
       do pp=1_idef,8_idef
          IR = INT(FLOOR((Y_R(pp)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y_PHI(pp)  + 0.5_rp*profiles_domain%DPHI)/ &
               profiles_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_3d%NR).OR.(IZ.GT.profiles_3d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END SIMD
    else
       !$OMP SIMD
!       !$OMP& aligned(IR,IZ)
       do pp=1_idef,8_idef
          IR = INT(FLOOR((Y_R(pp)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG2D(IR,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_2d%NR).OR.(IZ.GT.profiles_2d%NZ))) then
             flag(pp) = 0_is

!             write(6,'("Shit''s fucked.")')
          end if
       end do      
       !$OMP END SIMD
!       write(6,'("Shit''s not fucked.")')
    end if
  end subroutine check_if_in_profiles_domain_p

subroutine interp_2D_bfields(params,Y,B,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric magnetic
  !! field to the particles' position. @endnote
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  !! Core KORC simulation parameters.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
  !! Cartesian components of interpolated magnetic field components.
  !! B(1,:)=\(B_x\), B(2,:)=\(B_y\), and B(3,:)=\(B_z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.
  
  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,B,flag,bfield_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(bfield_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(bfield_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(bfield_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

        if (.not.params%GC_coords) then
           B(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
           B(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
           B(pp,3) = F(pp,3)
        else
           B(pp,1) = F(pp,1)
           B(pp,2) = F(pp,2)
           B(pp,3) = F(pp,3)
        end if
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_bfields

subroutine gradient_2D_Bfields(Y,BR,BPHI,BZ,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! gradient of the magnitude of themagnetic field to the particles'
  !! position. Stored as cylindrical components of field. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: BR
  !! Cylindrical components of gradient of R-component of magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: BPHI
  !! Cylindrical components of gradient of R-component of magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: BZ
  !! Cylindrical components of gradient of R-component of magnetic field.  
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.
  
  ss = size(Y,1)

  ALLOCATE(F(2,ss))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(Y,BR,BPHI,BZ,flag,bfield_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_gradient(bfield_2d%R, Y(pp,1), Y(pp,3), F(:,pp), ezerr)
        call EZspline_error(ezerr)
        
        BR(pp,1) = F(pp,1)
        BR(pp,2) = 0._rp
        BR(pp,3) = F(pp,2)
        
        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if
        
        call EZspline_gradient(bfield_2d%PHI, Y(pp,1), Y(pp,3), F(:,pp), &
             ezerr)
        call EZspline_error(ezerr)

        BPHI(pp,1) = F(pp,1)
        BPHI(pp,2) = 0._rp
        BPHI(pp,3) = F(pp,2)
        
        call EZspline_gradient(bfield_2d%Z, Y(pp,1), Y(pp,3), F(:,pp), ezerr)
        call EZspline_error(ezerr)

        BZ(pp,1) = F(pp,1)
        BZ(pp,2) = 0._rp
        BZ(pp,3) = F(pp,2)
        
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine gradient_2D_Bfields

subroutine interp_2D_gradBfields(Y,gradB,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! gradient of the magnitude of themagnetic field to the particles'
  !! position. Stored as cylindrical components of field. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: gradB
  !! Cylindirical components of interpolated gradient of magnitude of
  !! magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.
  
  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,gradB,flag,gradB_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(gradB_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)
        
        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if
        
        call EZspline_interp(gradB_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)
        
        call EZspline_interp(gradB_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

!        write(6,'("PS R-gradB ",E17.10)') F(pp,1)
!        write(6,'("PS PHI-gradB ",E17.10)') F(pp,2)
!        write(6,'("PS Z-gradB ",E17.10)') F(pp,3)
        
        gradB(pp,1) = F(pp,1)
        gradB(pp,2) = F(pp,2)
        gradB(pp,3) = F(pp,3)

!        write(6,'("PHI-gradB ",E17.10)') gradB(2,1)
        
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_gradBfields

subroutine interp_2D_curlbfields(Y,curlb,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! curl of the magnetic field unit vector to the particles'
  !! position. Stored as cylindrical components of field. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: curlb
  !! Cylindirical components of interpolated curl of direction of
  !! magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.
  
  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,curlb,flag,curlb_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(curlb_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(curlb_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(curlb_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

        curlb(pp,1) = F(pp,1)
        curlb(pp,2) = F(pp,2)
        curlb(pp,3) = F(pp,3)
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_curlbfields

subroutine interp_FOfields_p(Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flag_cache)

  REAL(rp),DIMENSION(8),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(8)   :: B_R,B_PHI  
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: E_X,E_Y,E_Z
  REAL(rp),DIMENSION(8)   :: E_R,E_PHI
  REAL(rp),DIMENSION(8)   :: cP,sP  
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: cc
  !! Particle chunk iterator.
  INTEGER(is),DIMENSION(8),INTENT(INOUT)   :: flag_cache

  call check_if_in_fields_domain_p(Y_R,Y_PHI,Y_Z,flag_cache)
  
  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z,efield_2d%R, &
       efield_2d%PHI,efield_2d%Z,8,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,ezerr)
  call EZspline_error(ezerr)
 
  !$OMP SIMD
!  !$OMP& aligned (cP,sP,B_X,B_Y,E_X,E_Y,Y_PHI,B_R,B_PHI,E_R,E_PHI)
  do cc=1_idef,8_idef
     cP(cc)=cos(Y_PHI(cc))
     sP(cc)=sin(Y_PHI(cc))
     
     B_X(cc) = B_R(cc)*cP(cc) - B_PHI(cc)*sP(cc)
     B_Y(cc) = B_R(cc)*sP(cc) + B_PHI(cc)*cP(cc)

     E_X(cc) = E_R(cc)*cP(cc) - E_PHI(cc)*sP(cc)
     E_Y(cc) = E_R(cc)*sP(cc) + E_PHI(cc)*cP(cc)

  end do
  !$OMP END SIMD
  
end subroutine interp_FOfields_p

subroutine interp_FOcollision_p(Y_R,Y_PHI,Y_Z,ne,Te,Zeff)

  REAL(rp),DIMENSION(8),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: ne,Te,Zeff


!  write(6,'("YR: ",E17.10)') Y_R(1)
!  write(6,'("YPHI: ",E17.10)') Y_PHI(1)
!  write(6,'("YZ: ",E17.10)') Y_Z(1)
  
!  write(6,'("Te_interp_R",E17.10)') profiles_2d%Te%x1
!  write(6,'("Te_interp_Z",E17.10)') profiles_2d%Te%x2
  
  call EZspline_interp(profiles_2d%ne,profiles_2d%Te, &
       profiles_2d%Zeff,8,Y_R,Y_Z,ne,Te,Zeff,ezerr)
  ! this will call PSPLINE routine EZspline_interp2_bmag_cloud_r8 as there
  ! is the same number of entries
  call EZspline_error(ezerr)
 
  
end subroutine interp_FOcollision_p

subroutine interp_fields_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
     curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)

  REAL(rp),DIMENSION(8),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: curlB_R,curlB_PHI,curlB_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: E_R,E_PHI,E_Z
  INTEGER(is),DIMENSION(8),INTENT(INOUT)   :: flag_cache

  call check_if_in_fields_domain_p(Y_R,Y_PHI,Y_Z,flag_cache)
  
  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z,efield_2d%R, &
       efield_2d%PHI,efield_2d%Z,gradB_2d%R,gradB_2d%PHI,gradB_2d%Z, &
       curlb_2d%R,curlb_2d%PHI,curlb_2d%Z,8,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_PHI,curlb_Z,ezerr)
  call EZspline_error(ezerr)

end subroutine interp_fields_p

subroutine interp_collision_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
     ne,Te,Zeff,flag_cache)

  REAL(rp),DIMENSION(8),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: ne,Te,Zeff
  INTEGER(is),DIMENSION(8),INTENT(INOUT)   :: flag_cache
!  INTEGER(ip) :: ezerr

  call check_if_in_profiles_domain_p(Y_R,Y_PHI,Y_Z,flag_cache)
  
  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z,efield_2d%R, &
       efield_2d%PHI,efield_2d%Z,profiles_2d%ne,profiles_2d%Te, &
       profiles_2d%Zeff,8,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,ne,Te,Zeff,ezerr)
  call EZspline_error(ezerr)
 

end subroutine interp_collision_p

subroutine interp_bmag_p(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

  REAL(rp),DIMENSION(8),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(8),INTENT(OUT)   :: B_R,B_PHI,B_Z
!  INTEGER(ip) :: ezerr

  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
       8,Y_R,Y_Z,B_R,B_PHI,B_Z,ezerr)
  call EZspline_error(ezerr)
 

end subroutine interp_bmag_p

!> @brief Subroutine for interpolating the pre-computed, 3-D magnetic field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=\(B_x\), B(2,:)=\(B_y\), and B(3,:)=\(B_z\).
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_bfields(Y,B,flag)
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  INTEGER                                                :: pp
  INTEGER                                                :: ss

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,B,flag,bfield_3d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(bfield_3d%R, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(bfield_3d%PHI, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(bfield_3d%Z, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,3), ezerr)
        call EZspline_error(ezerr)

        B(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
        B(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
        B(pp,3) = F(pp,3)
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_3D_bfields


!> @brief Subroutine that calculates the axisymmetric magnetic field to the particles' position using the poloidal magnetic flux.
!! @details When the poloidal magnetic flux \(\Psi(R,Z)\) is used in a KORC simulation, the magnetic field components are calculated as it follows:
!!
!!
!! $$B_R = \frac{1}{R}\frac{\partial \Psi}{\partial Z},$$
!! $$B_\phi = \frac{RoBo}{R}\),$$
!! $$B_Z = -\frac{1}{R}\frac{\partial \Psi}{\partial R},$$
!!
!!
!! where \(Ro\) and \(Bo\) are the radial position of the magnetic axis and the magnetic field as measured at the magnetic axis, respectively.
!! First, the derivatives of the poloidal magnetic flux are calculated at the particles' position using the PSPLINE interpolant of
!! the poloidal magnetic flux. Then, we calculate the cylindrical components of the magnetic field, and finally we calculate its Cartesian
!! components that will be used in the particle pusher.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation.
!! See korc_types.f90 and korc_fields.f90.
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=\(B_x\), B(2,:)=\(B_y\), and B(3,:)=\(B_x\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param A Variable containing the partial derivatives of the poloidal magnetic flux \(\Psi(R,Z)\) and the cylindrical components
!! of the magnetic field (its value changes through the subroutine).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine calculate_magnetic_field(Y,F,B,PSI_P,flag)
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: PSI_P
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: A
  INTEGER                                                :: pp
  INTEGER                                                :: ss


  
  ss = size(Y,1)

  ALLOCATE(A(ss,3))
  A=0._rp
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,A,B,flag,bfield_2d,PSI_P)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then

!        write(6,'("pp: ",I16)') pp
        
!        write(6,'("Y_R: ",E17.10)') Y(:,1)
!        write(6,'("Y_PHI: ",E17.10)') Y(:,2)
!        write(6,'("Y_Z: ",E17.10)') Y(:,3)
!        write(6,'("PSI_P: ",E17.10)') PSI_P
        
        call EZspline_interp(bfield_2d%A, Y(pp,1), Y(pp,3), &
             PSI_P(pp), ezerr)
        call EZspline_error(ezerr)
        
        ! FR = (dA/dZ)/R
        call EZspline_derivative(bfield_2d%A, 0, 1, Y(pp,1), Y(pp,3), &
             A(pp,1), ezerr)
        !			call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        else
           A(pp,1) = A(pp,1)/Y(pp,1)

           ! FPHI = Fo*Ro/R
           A(pp,2) = F%Bo*F%Ro/Y(pp,1)

           ! FR = -(dA/dR)/R
           call EZspline_derivative(bfield_2d%A, 1, 0, Y(pp,1), Y(pp,3), &
                A(pp,3), ezerr)
           call EZspline_error(ezerr)
           A(pp,3) = -A(pp,3)/Y(pp,1)

!           write(6,'("B_R: ",E17.10)') A(:,1)
!           write(6,'("B_PHI: ",E17.10)') A(:,2)
!           write(6,'("B_Z: ",E17.10)') A(:,3)

!           write(6,'("B_X: ",E17.10)') B(:,1)
!           write(6,'("B_Y: ",E17.10)') B(:,2)
!           write(6,'("B_Z: ",E17.10)') B(:,3)

           
           B(pp,1) = A(pp,1)*COS(Y(pp,2)) - A(pp,2)*SIN(Y(pp,2))
           B(pp,2) = A(pp,1)*SIN(Y(pp,2)) + A(pp,2)*COS(Y(pp,2))
           B(pp,3) = A(pp,3)
        end if
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(A)
end subroutine calculate_magnetic_field

subroutine calculate_initial_magnetic_field(F)

  TYPE(FIELDS), INTENT(INOUT)                               :: F
  REAL(rp),dimension(F%dims(1),F%dims(3),2)                  :: gradA
  INTEGER                                                :: ii
  INTEGER                                                :: jj
        
        ! FR = (dA/dZ)/R
        call EZspline_gradient(bfield_2d%A,F%dims(1),F%dims(3),F%X%R, F%X%Z, &
             gradA, ezerr)
        call EZspline_error(ezerr)

        do ii=1,F%dims(1)
           F%B_2D%R(ii,:) = gradA(ii,:,2)/F%X%R(ii)
           F%B_2D%PHI(ii,:) = F%Bo*F%Ro/F%X%R(ii)
           F%B_2D%Z(ii,:) = -gradA(ii,:,1)/F%X%R(ii)
        end do
           
!        write(6,'("AR",E17.10)') gradA(1)
!        write(6,'("AZ",E17.10)') gradA(2)        
       
end subroutine calculate_initial_magnetic_field

!> @brief Subroutine for interpolating the pre-computed, axisymmetric electric field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in,out] E Cartesian components of interpolated electric field components. E(1,:)=\(E_x\), E(2,:)=\(E_y\), and E(3,:)=\(E_z\).
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=\(E_R\), F(2,:)=\(E_\phi\), and F(3,:)=\(E_Z\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_2D_efields(params,Y,E,flag)
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  !! Core KORC simulation parameters.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  INTEGER                                                :: pp
  INTEGER                                                :: ss

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,E,flag,efield_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(efield_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(efield_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(efield_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

        if (.not.params%GC_coords) then
           E(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
           E(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
           E(pp,3) = F(pp,3)
        else
           E(pp,1) = F(pp,1)
           E(pp,2) = F(pp,2)
           E(pp,3) = F(pp,3)
        end if
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_efields


!> @brief Subroutine for interpolating the pre-computed 3-D electric field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in,out] E Cartesian components of interpolated electric field components. E(1,:)=\(E_x\), E(2,:)=\(E_y\), and E(3,:)=\(E_z\).
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=\(E_R\), F(2,:)=\(E_\phi\), and F(3,:)=\(E_Z\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_efields(Y,E,flag)
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  INTEGER                                                :: pp
  INTEGER                                                :: ss

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,E,flag,efield_3d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(efield_3d%R, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(efield_3d%PHI, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(efield_3d%Z, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,3), ezerr)
        call EZspline_error(ezerr)

        E(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
        E(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
        E(pp,3) = F(pp,3)
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_3D_efields


subroutine interp_fields(params,prtcls,F)
  !! @note Subroutine that works as an interface for calling the 
  !! appropriate subroutines for interpolating or calculating the 
  !! electric and magnetic fields. @endnote
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
  TYPE(PARTICLES), INTENT(INOUT) :: prtcls
    !! An instance of PARTICLES containing the variables of a given species.
  TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of KORC's derived type FIELDS containing all the 
    !! information about the fields used in the simulation.
    !! See [[korc_types]] and [[korc_fields]].

  if (.not.params%GC_coords) call cart_to_cyl(prtcls%X,prtcls%Y)
  
!  write(6,'("BR: ",E17.10)') prtcls%BR(:,1)
!  write(6,'("Y: ",E17.10)') prtcls%X(2,1)
!  write(6,'("Z: ",E17.10)') prtcls%X(3,1)
  
  call check_if_in_fields_domain(prtcls%Y, prtcls%flag)

  if (ALLOCATED(F%PSIp).and.F%Bflux) then

!     write(6,'("3 size of PSI_P: ",I16)') size(prtcls%PSI_P)

!     write(6,'("B_X: ",E17.10)') prtcls%B(:,1)
!     write(6,'("B_Z: ",E17.10)') prtcls%B(:,3)
!     write(6,'("B_Y: ",E17.10)') prtcls%B(:,2)
!     write(6,'("PSI_P: ",E17.10)') prtcls%PSI_P
     
     call calculate_magnetic_field(prtcls%Y,F,prtcls%B,prtcls%PSI_P,prtcls%flag)
  end if

  if (ALLOCATED(F%B_2D%R)) then
     call interp_2D_bfields(params,prtcls%Y,prtcls%B,prtcls%flag)
  end if

  if (ALLOCATED(F%B_3D%R)) then
     call interp_3D_bfields(prtcls%Y,prtcls%B,prtcls%flag)
  end if

  if (ALLOCATED(F%E_2D%R)) then
     call interp_2D_efields(params,prtcls%Y,prtcls%E,prtcls%flag)
  end if

  if (ALLOCATED(F%E_3D%R)) then
     call interp_3D_efields(prtcls%Y,prtcls%E,prtcls%flag)
  end if

  if (params%GC_coords.and.ALLOCATED(F%gradB_2D%R)) then
     call interp_2D_gradBfields(prtcls%Y,prtcls%gradB,prtcls%flag)

  end if

  if (params%GC_coords.and.ALLOCATED(F%gradB_2D%R)) then
     call interp_2D_curlbfields(prtcls%Y,prtcls%curlb,prtcls%flag)
  end if

  if(params%GC_coords.and.params%orbit_model(3:6)=='grad') then
     call gradient_2D_bfields(prtcls%Y,prtcls%BR,prtcls%BPHI, &
          prtcls%BZ,prtcls%flag)
  end if
  
end subroutine interp_fields


subroutine interp_2D_profiles(Y,ne,Te,Zeff,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! plasma profiles to the particles' position. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates,
  !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: ne
  !! Interpolated background electron density !!\(n_e(R,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Te
  !! Interpolated background electron temperature \(T_e(R,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Zeff
  !! Interpolated effective charge number \(Z_{eff}(R,Z)\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the
  !! simulation (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  if (Y(2,1).eq.0) then
     ss=1_idef
  else
     ss = size(Y,1)
  end if

!  write(6,'("Also R_buffer: ",E17.10)') Y(1,ss)
  
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(Y,ne,Te,Zeff,flag,profiles_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(profiles_2d%ne, Y(pp,1), Y(pp,3), ne(pp), ezerr)
        call EZspline_error(ezerr)

!        write(6,'("Also R_buffer: ",E17.10)') Y(pp,1)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(profiles_2d%Te, Y(pp,1), Y(pp,3), Te(pp), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(profiles_2d%Zeff, Y(pp,1), Y(pp,3), Zeff(pp), ezerr)
        call EZspline_error(ezerr)
     end if
  end do
  !$OMP END PARALLEL DO
end subroutine interp_2D_profiles


subroutine interp_3D_profiles(Y,ne,Te,Zeff,flag)
  !! @note Subroutine for interpolating the pre-computed,
  !! 3-D plasma profiles to the particles' position. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates,
  !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: ne
  !! Interpolated background electron density \(n_e(R,\phi,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Te
  !! Interpolated background electron temperature \(T_e(R,\phi,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Zeff
  !! Interpolated effective charge number \(Z_{eff}(R,\phi,Z)\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in
  !! the simulation (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.
  
  ss = size(Y,1)

  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(Y,ne,Te,Zeff,flag,profiles_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(profiles_3d%ne, Y(pp,1), Y(pp,2), Y(pp,3), &
             ne(pp), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(profiles_3d%Te, Y(pp,1), Y(pp,2), Y(pp,3), &
             Te(pp), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(profiles_3d%Zeff, Y(pp,1), Y(pp,2), Y(pp,3), &
             Zeff(pp), ezerr)
        call EZspline_error(ezerr)
     end if
  end do
  !$OMP END PARALLEL DO
end subroutine interp_3D_profiles


subroutine interp_profiles(params,prtcls,P)
  !! @note Subroutine that calls the appropriate subroutines for
  !! interpolating the 2-D or 3-D plasma profiles to the particles'
  !! position. @endnote
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  !! Core KORC simulation parameters.
  TYPE(PARTICLES), INTENT(INOUT) :: prtcls
  !! An instance of PARTICLES containing the variables of a
  !! given species. Call to this subroutine generally passes spp%vars.
  TYPE(PROFILES), INTENT(IN)     :: P
  !! An instance of KORC's derived type PROFILES containing all the
  !! information about the plasma profiles used in the simulation.
  !! See[[ korc_types]] and [[korc_profiles]].

  if (.not.params%GC_coords) call cart_to_cyl(prtcls%X,prtcls%X)
  
  !write(6,'("Also R_buffer: ",E17.10)') prtcls%Y(1,1)  

  call check_if_in_profiles_domain(prtcls%Y, prtcls%flag)

  if (ALLOCATED(P%ne_2D)) then
!     write(6,'("Also R_buffer: ",E17.10)') prtcls%X(1,1)
     call interp_2D_profiles(prtcls%Y,prtcls%ne,prtcls%Te,prtcls%Zeff, &
          prtcls%flag)
  else if (ALLOCATED(P%ne_3D)) then
     call interp_3D_profiles(prtcls%Y,prtcls%ne,prtcls%Te,prtcls%Zeff, &
          prtcls%flag)
  else
     write(6,'("Error: NO PROFILES ALLOCATED")')
     call KORC_ABORT()
  end if
end subroutine interp_profiles


!> @brief Subroutine that frees memory allocated for PSPLINE interpolants.
!!
!! @param[in] params Core KORC simulation parameters.
subroutine finalize_interpolants(params)
  TYPE(KORC_PARAMS), INTENT(IN) :: params

  if (params%field_model .EQ. 'EXTERNAL') then
     if (params%mpi_params%rank .EQ. 0) then
        write(6,'("* * * * FINALIZING FIELD INTERPOLANT * * * *")')
     end if

     if (EZspline_allocated(bfield_3d%R)) call Ezspline_free(bfield_3d%R, ezerr)
     if (EZspline_allocated(bfield_3d%PHI)) &
          call Ezspline_free(bfield_3d%PHI,ezerr)
     
     if (EZspline_allocated(bfield_3d%Z)) call Ezspline_free(bfield_3d%Z, ezerr)
     if (EZspline_allocated(bfield_2d%A)) call Ezspline_free(bfield_2d%A, ezerr)
     if (EZspline_allocated(bfield_2d%R)) call Ezspline_free(bfield_2d%R, ezerr)
     if (EZspline_allocated(bfield_2d%PHI)) &
          call Ezspline_free(bfield_2d%PHI,ezerr)
     
     if (EZspline_allocated(bfield_2d%Z)) call Ezspline_free(bfield_2d%Z, ezerr)

     if (EZspline_allocated(gradB_2d%R)) call Ezspline_free(gradB_2d%R, ezerr)
     if (EZspline_allocated(gradB_2d%PHI)) &
          call Ezspline_free(gradB_2d%PHI, ezerr)
     
     if (EZspline_allocated(gradB_2d%Z)) call Ezspline_free(gradB_2d%Z, ezerr)

     if (EZspline_allocated(curlb_2d%R)) call Ezspline_free(curlb_2d%R, ezerr)
     if (EZspline_allocated(curlb_2d%PHI)) &
          call Ezspline_free(curlb_2d%PHI, ezerr)
     
     if (EZspline_allocated(curlb_2d%Z)) call Ezspline_free(curlb_2d%Z, ezerr)

     if (ALLOCATED(profiles_domain%FLAG2D)) DEALLOCATE(profiles_domain%FLAG2D)
     if (ALLOCATED(profiles_domain%FLAG3D)) DEALLOCATE(profiles_domain%FLAG3D)

     if (params%mpi_params%rank .EQ. 0) then
        write(6,'("* * * * FIELD INTERPOLANT FINALIZED * * * *")')
     end if
  end if
     
  if (params%profile_model .EQ. 'EXTERNAL') then
     if (params%mpi_params%rank .EQ. 0) then
        write(6,'("* * * * FINALIZING PROFILE INTERPOLANT * * * *")')
     end if
     
     if (EZspline_allocated(profiles_3d%ne)) &
          call Ezspline_free(profiles_3d%ne,ezerr)
     
     if (EZspline_allocated(profiles_3d%Te)) &
          call Ezspline_free(profiles_3d%Te,ezerr)
     
     if (EZspline_allocated(profiles_3d%Zeff)) call Ezspline_free( &
          profiles_3d%Zeff, ezerr)
     if (EZspline_allocated(profiles_2d%ne)) &
          call Ezspline_free(profiles_2d%ne,ezerr)
     
     if (EZspline_allocated(profiles_2d%Te)) &
          call Ezspline_free(profiles_2d%Te,ezerr)
     
     if (EZspline_allocated(profiles_2d%Zeff)) call Ezspline_free( &
          profiles_2d%Zeff, ezerr)

     if (ALLOCATED(profiles_domain%FLAG2D)) DEALLOCATE(profiles_domain%FLAG2D)
     if (ALLOCATED(profiles_domain%FLAG3D)) DEALLOCATE(profiles_domain%FLAG3D)

     if (params%mpi_params%rank .EQ. 0) then
        write(6,'("* * * * PROFILE INTERPOLANT FINALIZED * * * *")')
     end if
  end if
end subroutine finalize_interpolants
end module korc_interp
