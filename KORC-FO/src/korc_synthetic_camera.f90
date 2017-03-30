MODULE korc_synthetic_camera
	USE korc_types
	USE korc_constants
	USE korc_HDF5
    USE mpi
    USE special_functions

	IMPLICIT NONE

	TYPE, PRIVATE :: POLOIDAL_PLANE
		REAL(rp), DIMENSION(:), ALLOCATABLE :: nodes_R ! In meters
		REAL(rp), DIMENSION(:), ALLOCATABLE :: nodes_Z ! In meters

		REAL(rp) :: DR
		REAL(rp) :: DZ	
		REAL(rp) :: Rmax, Rmin
		REAL(rp) :: Zmax, Zmin

		INTEGER, DIMENSION(2) :: grid_dims ! Number of pixels (R,Z)
	END TYPE POLOIDAL_PLANE


	TYPE, PRIVATE :: CAMERA
		REAL(rp) :: aperture ! Aperture of the camera (diameter of lens) in meters
		REAL(rp) :: Riw ! Radial position of inner wall
		INTEGER, DIMENSION(2) :: num_pixels ! Number of pixels (X,Y)
		REAL(rp), DIMENSION(2) :: sensor_size ! In meters (horizontal,vertical)
		REAL(rp) :: focal_length ! Focal length in meters
		REAL(rp), DIMENSION(2) :: position ! Position of camera (R,Z)
		REAL(rp) :: incline ! Incline of camera in degrees
		REAL(rp) :: horizontal_angle_view ! Horizontal angle of view in radians
		REAL(rp) :: vertical_angle_view ! Vertical angle of view in radians
		REAL(rp), DIMENSION(:), ALLOCATABLE :: pixels_nodes_x ! In meters
		REAL(rp), DIMENSION(:), ALLOCATABLE :: pixels_nodes_y ! In meters
		REAL(rp), DIMENSION(:), ALLOCATABLE :: pixels_edges_x ! In meters
		REAL(rp), DIMENSION(:), ALLOCATABLE :: pixels_edges_y ! In meters

		REAL(rp) :: lambda_min ! Minimum wavelength in cm
		REAL(rp) :: lambda_max ! Maximum wavelength in cm
		INTEGER :: Nlambda
		REAL(rp) :: Dlambda ! In cm
		REAL(rp), DIMENSION(:), ALLOCATABLE :: lambda ! In cm
	END TYPE CAMERA

	TYPE, PRIVATE :: ANGLES
		REAL(rp), DIMENSION(:), ALLOCATABLE :: eta
		REAL(rp), DIMENSION(:), ALLOCATABLE :: beta
		REAL(rp), DIMENSION(:), ALLOCATABLE :: psi

		REAL(rp) :: threshold_angle
		REAL(rp) :: threshold_radius
	END TYPE ANGLES

	TYPE(CAMERA), PRIVATE :: cam
	TYPE(ANGLES), PRIVATE :: ang
	TYPE(POLOIDAL_PLANE), PRIVATE :: pplane

	REAL(rp), PRIVATE, PARAMETER :: Tol = 1.0E-5_rp
	REAL(rp), PRIVATE, PARAMETER :: CGS_h = 6.6261E-27_rp ! Planck constant in erg*s
	REAL(rp), PRIVATE, PARAMETER :: CGS_C = 1.0E2_rp*C_C
	REAL(rp), PRIVATE, PARAMETER :: CGS_E = 3.0E9_rp*C_E
	REAL(rp), PRIVATE, PARAMETER :: CGS_ME = 1.0E3_rp*C_ME

	INTERFACE save_snapshot_var
	  module procedure save_snapshot_var_3d,save_snapshot_var_4d
	END INTERFACE

	PRIVATE :: clockwise_rotation,anticlockwise_rotation,cross,check_if_visible,calculate_rotation_angles,&
				zeta,fx,arg,Po,P1,Psyn,chic,psic,save_synthetic_camera_params,save_snapshot,besselk,&
				spectral_angular_density,spectral_density,IntK,trapz,save_snapshot_var
	PUBLIC :: initialize_synthetic_camera,synthetic_camera

	CONTAINS

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !

SUBROUTINE test(params)
    IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
    REAL(rp), DIMENSION(3,2,4) :: A, C
    REAL(rp), DIMENSION(24) :: B, R
    INTEGER :: ii,jj,kk, ierr


    A = 0_rp

    do kk=1_idef,4
        do jj=1_idef,2
            do ii=1_idef,3        
                A(ii,jj,kk) = REAL(ii,rp) + 3_rp*(REAL(jj,rp) - 1.0_rp) + 6.0_rp*(REAL(kk,rp) - 1.0_rp)
                A(ii,jj,kk) = REAL(params%mpi_params%rank,rp)*A(ii,jj,kk)
                write(6,'("A(",I1,",",I1,",",I1,")=",F25.16)') ii,jj,kk,A(ii,jj,kk)
            end do
        end do
    end do

    B = RESHAPE(A,(/24/))

    CALL MPI_REDUCE(B,R,24,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (params%mpi_params%rank.EQ.0_idef) then
        do ii=1_idef,24
            write(6,'("R(",I2,")=",F25.16)') ii,R(ii)
        end do


        C = RESHAPE(B,(/3,2,4/))

        do kk=1_idef,4
            do jj=1_idef,2
                do ii=1_idef,3        
                    write(6,'("C(",I1,",",I1,",",I1,")=",F25.16)') ii,jj,kk,C(ii,jj,kk)
                end do
            end do
        end do
    end if
END SUBROUTINE test


SUBROUTINE testbesselkv()
	IMPLICIT NONE
	REAL(rp) :: v
	REAL(rp), DIMENSION(:), ALLOCATABLE :: x
	REAL(rp), DIMENSION(:), ALLOCATABLE :: R
	INTEGER :: nx
	INTEGER :: ii
    REAL(4) :: xnu,ri,rk,rip,rkp

	nx = 1000
	v = 1.0_rp/3.0_rp

	ALLOCATE(x(nx))
	ALLOCATE(R(nx))

	do ii=1_idef,nx
		x(ii) = REAL(ii,rp)*0.01_rp
	end do
	
	do ii=1_idef,nx
        call bessik(REAL(x(ii),4),REAL(v,4),ri,rk,rip,rkp)
		write(6,'(F25.16)') rk
	end do

	DEALLOCATE(x)
	DEALLOCATE(R)
END SUBROUTINE testbesselkv


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !


SUBROUTINE initialize_synthetic_camera(params,F)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp) :: aperture ! Aperture of the camera (diameter of lens) in meters
	REAL(rp) :: Riw ! Radial position of inner wall
	INTEGER, DIMENSION(2) :: num_pixels ! Number of pixels (X,Y)
	REAL(rp), DIMENSION(2) :: sensor_size ! (horizontal,vertical)
	REAL(rp) :: focal_length
	REAL(rp), DIMENSION(2) :: position ! Position of camera (R,Z)
	REAL(rp) :: incline
	REAL(rp) :: lambda_min ! Minimum wavelength in cm
	REAL(rp) :: lambda_max ! Maximum wavelength in cm
	INTEGER :: Nlambda
	REAL(rp) :: xmin, xmax, ymin, ymax, DX, DY
	INTEGER :: ii

	NAMELIST /SyntheticCamera/ aperture,Riw,num_pixels,sensor_size,focal_length,&
			position,incline,lambda_min,lambda_max,Nlambda

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * * * * * * * * * * * * * *")')
		write(6,'("*  Initializing synthetic camera  *")')
	end if

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=SyntheticCamera)
	close(default_unit_open)

!	write(*,nml=SyntheticCamera)

	cam%aperture = aperture
	cam%Riw = Riw
	cam%num_pixels = num_pixels
	cam%sensor_size = sensor_size
	cam%focal_length = focal_length
	cam%position = position
	cam%incline = C_PI*incline/180.0_rp
	cam%horizontal_angle_view = ATAN2(0.5_rp*cam%sensor_size(1),cam%focal_length)
	cam%vertical_angle_view = ATAN2(0.5_rp*cam%sensor_size(2),cam%focal_length)

	cam%lambda_min = 1.0E2_rp*lambda_min ! In cm
	cam%lambda_max = 1.0E2_rp*lambda_max ! In cm
	cam%Nlambda = Nlambda
	cam%Dlambda = (cam%lambda_max - cam%lambda_min)/REAL(cam%Nlambda,rp)
	ALLOCATE(cam%lambda(cam%Nlambda))
	
	do ii=1_idef,cam%Nlambda
		cam%lambda(ii) = cam%lambda_min + REAL(ii-1_idef,rp)*cam%Dlambda
	end do

	! NOTE: cam%lambda is in cm !!

	ALLOCATE(cam%pixels_nodes_x(cam%num_pixels(1)))
	ALLOCATE(cam%pixels_nodes_y(cam%num_pixels(2)))
	ALLOCATE(cam%pixels_edges_x(cam%num_pixels(1) + 1))
	ALLOCATE(cam%pixels_edges_y(cam%num_pixels(2) + 1))

	xmin = -0.5_rp*cam%sensor_size(1)
	xmax = 0.5_rp*cam%sensor_size(1)
	DX = cam%sensor_size(1)/REAL(cam%num_pixels(1),rp)

	do ii=1_idef,cam%num_pixels(1)
		cam%pixels_nodes_x(ii) = xmin + 0.5_rp*DX + REAL(ii-1_idef,rp)*DX
	end do

	do ii=1_idef,cam%num_pixels(1)+1_idef
		cam%pixels_edges_x(ii) = xmin + REAL(ii-1_idef,rp)*DX
	end do

	ymin = cam%position(2) - 0.5_rp*cam%sensor_size(2)
	ymax = cam%position(2) + 0.5_rp*cam%sensor_size(2)
	DY = cam%sensor_size(2)/REAL(cam%num_pixels(2),rp)

	do ii=1_idef,cam%num_pixels(2)
		cam%pixels_nodes_y(ii) = ymin + 0.5_rp*DY + REAL(ii-1_idef,rp)*DY
	end do

	do ii=1_idef,cam%num_pixels(2)+1_idef
		cam%pixels_edges_y(ii) = ymin + REAL(ii-1_idef,rp)*DY
	end do
	
	! Initialize ang variables
	ALLOCATE(ang%eta(cam%num_pixels(1)))
	ALLOCATE(ang%beta(cam%num_pixels(1)))
	ALLOCATE(ang%psi(cam%num_pixels(2)+1_idef))

	do ii=1_idef,cam%num_pixels(1)
		ang%eta(ii) = ABS(ATAN2(cam%pixels_nodes_x(ii),cam%focal_length))
		if (cam%pixels_edges_x(ii) .LT. 0.0_rp) then
			ang%beta(ii) = 0.5_rp*C_PI - cam%incline - ang%eta(ii)
		else
			ang%beta(ii) = 0.5_rp*C_PI - cam%incline + ang%eta(ii)
		end if
	end do

	do ii=1_idef,cam%num_pixels(2)+1_idef
		ang%psi(ii) = ATAN2(cam%pixels_edges_y(ii),cam%focal_length)
	end do

	ang%threshold_angle = ATAN2(cam%Riw,-cam%position(1))
	ang%threshold_radius = SQRT(cam%Riw**2 + cam%position(1)**2)

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'("*     Synthetic camera ready!     *")')
		write(6,'("* * * * * * * * * * * * * * * * * *",/)')
	end if


	! Initialize poloidal plane parameters
	
	pplane%grid_dims = num_pixels
	ALLOCATE(pplane%nodes_R(pplane%grid_dims(1)))
	ALLOCATE(pplane%nodes_Z(pplane%grid_dims(2)))

	! * * * * * * * ALL IN METERS * * * * * * * 

	IF (TRIM(params%magnetic_field_model) .EQ. 'ANALYTICAL') THEN
		pplane%Rmin = F%Ro - F%AB%a
		pplane%Rmax = F%Ro + F%AB%a
		pplane%Zmin = -F%AB%a
		pplane%Zmax = F%AB%a
	ELSE
		pplane%Rmin = MINVAL(F%X%R)
		pplane%Rmax = MAXVAL(F%X%R)
		pplane%Zmin = MINVAL(F%X%Z)
		pplane%Zmax = MAXVAL(F%X%Z)
	END IF	

	pplane%DR = (pplane%Rmax - pplane%Rmin)/REAL(pplane%grid_dims(1),rp)
	pplane%DZ = (pplane%Zmax - pplane%Zmin)/REAL(pplane%grid_dims(2),rp)

	do ii=1_idef,pplane%grid_dims(1)
		pplane%nodes_R(ii) = pplane%Rmin + 0.5_rp*pplane%DR + REAL(ii-1_idef,rp)*pplane%DR
	end do

	do ii=1_idef,pplane%grid_dims(2)
		pplane%nodes_Z(ii) = pplane%Zmin + 0.5_rp*pplane%DZ + REAL(ii-1_idef,rp)*pplane%DZ
	end do
	! * * * * * * * ALL IN METERS * * * * * * * 

	call save_synthetic_camera_params(params)
END SUBROUTINE initialize_synthetic_camera


! * * * * * * * * * * * * * * * !
! * * * * * FUNCTIONS * * * * * !
! * * * * * * * * * * * * * * * !

FUNCTION cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
END FUNCTION cross


FUNCTION clockwise_rotation(x,t)
	IMPLICIT NONE
	REAL(rp), DIMENSION(2), INTENT(IN) :: x
	REAL(rp), INTENT(IN) :: t ! Angle in radians
	REAL(rp), DIMENSION(2) :: clockwise_rotation

	clockwise_rotation(1) = x(1)*COS(t) + x(2)*SIN(t)
	clockwise_rotation(2) = -x(1)*SIN(t) + x(2)*COS(t)
END FUNCTION clockwise_rotation


FUNCTION anticlockwise_rotation(x,t)
	IMPLICIT NONE
	REAL(rp), DIMENSION(2), INTENT(IN) :: x
	REAL(rp), INTENT(IN) :: t ! Angle in radians
	REAL(rp), DIMENSION(2) :: anticlockwise_rotation

	anticlockwise_rotation(1) = x(1)*COS(t) - x(2)*SIN(t)
	anticlockwise_rotation(2) = x(1)*SIN(t) + x(2)*COS(t)
END FUNCTION anticlockwise_rotation


FUNCTION besselk(v,x)
	IMPLICIT NONE
	REAL(rp) :: besselk
	REAL(rp), INTENT(IN) :: x
	REAL(rp), INTENT(IN) :: v
	REAL(4) :: ri,rk,rip,rkp

	call bessik(REAL(x,4),REAL(v,4),ri,rk,rip,rkp)
	besselk = REAL(rk,rp)
END FUNCTION besselk


FUNCTION zeta(g,p,k,l)
	IMPLICIT NONE
	REAL(rp) :: zeta
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: p
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l

	zeta = (2.0_rp*C_PI/(3.0_rp*l*k*g**3))*(1.0_rp + (g*p)**2)**1.5_rp
END FUNCTION


FUNCTION fx(g,p,x)
	IMPLICIT NONE
	REAL(rp) :: fx
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: p
	REAL(rp), INTENT(IN) :: x

	fx = g*x/SQRT(1.0_rp + (g*p)**2)
END FUNCTION fx


FUNCTION Po(g,p,k,l)
	IMPLICIT NONE
	REAL(rp) :: Po
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: p
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l
	

	Po = -(4.0_rp*C_PI*CGS_C*CGS_E**2/(k*(l*g)**4))*(1.0_rp + (g*p)**2)**2
END FUNCTION Po


FUNCTION arg(g,p,k,l,x)
	IMPLICIT NONE
	REAL(rp) :: arg
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: p
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l
	REAL(rp), INTENT(IN) :: x
	REAL(rp) :: A

	A = fx(g,p,x)
	arg = 1.5_rp*zeta(g,p,k,l)*(A + (A**3)/3.0_rp)
END FUNCTION arg


FUNCTION P1(g,p,k,l,x)
	IMPLICIT NONE
	REAL(rp) :: P1
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: p
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l
	REAL(rp), INTENT(IN) :: x
	REAL(rp) :: BK13
	REAL(rp) :: BK23
	REAL(rp) :: v
	REAL(rp) :: A

	v = 1.0_rp/3.0_rp
	BK13 = besselk(v,zeta(g,p,k,l))

	v = 2.0_rp/3.0_rp
	BK23 = besselk(v,zeta(g,p,k,l))

	A = fx(g,p,x)

	P1 = ((g*p)**2)*BK13*COS(arg(g,p,k,l,x))/(1.0_rp + (g*p)**2) - 0.5_rp*BK13*(1.0_rp + A**2)*COS(arg(g,p,k,l,x))&
		+ A*BK23*SIN(arg(g,p,k,l,x))
END FUNCTION P1


FUNCTION Psyn(g,p,k,l,x)
	IMPLICIT NONE
	REAL(rp) :: Psyn
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: p
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l
	REAL(rp), INTENT(IN) :: x

	Psyn = Po(g,p,k,l)*P1(g,p,k,l,x)
END FUNCTION Psyn


FUNCTION chic(g,k,l)
	IMPLICIT NONE
	REAL(rp) :: chic
	REAL(rp), INTENT(IN) ::	g
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l
	REAL(rp) :: D
	REAL(rp) :: xi

	xi = 2.0_rp*C_PI/(3.0_rp*l*k*g**3)
	D = (0.5_rp*(SQRT(4.0_rp + (C_PI/xi)**2) - C_PI/xi))**(1.0_rp/3.0_rp)
	chic = (1.0_rp/D - D)/g
END FUNCTION chic


FUNCTION psic(k,l)
	IMPLICIT NONE
	REAL(rp) :: psic
	REAL(rp), INTENT(IN) :: k
	REAL(rp), INTENT(IN) :: l

	psic = (1.5_rp*k*l/C_PI)**(1.0_rp/3.0_rp)
END FUNCTION psic


FUNCTION IntK(v,x)
	IMPLICIT NONE
	REAL(rp) :: IntK
	REAL(rp), INTENT(IN) :: v
	REAL(rp), INTENT(IN) :: x

!	IntK = (C_PI/SQRT(2.0_rp))*(1.0_rp - 0.25_rp*(4.0_rp*v**2 - 1.0_rp))*(1 - ERF(SQRT(x))) +&
!			0.25_rp*(4.0_rp*v**2 - 1.0_rp)*SQRT(0.5_rp*C_PI/x)*EXP(-x)
	IntK = (C_PI/SQRT(2.0_rp))*(1.0_rp - 0.25_rp*(4.0_rp*v**2 - 1.0_rp))*ERFC(SQRT(x))&
			 + 0.25_rp*(4.0_rp*v**2 - 1.0_rp)*SQRT(0.5_rp*C_PI/x)*EXP(-x)
END FUNCTION IntK


FUNCTION trapz(a,b)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: a
	REAL(rp), INTENT(IN) :: b
	REAL(rp) :: trapz
	REAL(rp) :: Iold, Inew, rerr
	REAL(rp) :: v,h,z
	INTEGER :: ii,jj,npoints
	LOGICAL :: flag
	
	v = 5.0_rp/3.0_rp

	h = b - a
	trapz = 0.5*(besselk(v,a) + besselk(v,b))
	
	ii = 1_idef
	flag = .TRUE.
	do while (flag)
		Iold = trapz*h

		ii = ii + 1_idef
		npoints = 2_idef**(ii-2_idef)
		h = 0.5_rp*(b-a)/REAL(npoints,rp)

		do jj=1_idef,npoints
			z = a + h + 2.0_rp*(REAL(jj,rp) - 1.0_rp)*h
			trapz = trapz + besselk(v,z)
		end do

		Inew = trapz*h

		rerr = ABS((Inew - Iold)/Iold)

		flag = .NOT.(rerr.LT.Tol)
	end do
	trapz = Inew
END FUNCTION trapz


FUNCTION P_integral(z)
	IMPLICIT NONE
	REAL(rp) :: P_integral
	REAL(rp), INTENT(IN) :: z
	REAL(rp) :: a

	IF (z .LT. 0.5_rp) THEN
		a = (2.16_rp/2.0_rp**(2.0_rp/3.0_rp))*z**(1.0_rp/3.0_rp)
		P_integral = P_integral + trapz(z,a)
		P_integral = P_integral + IntK(5.0_rp/3.0_rp,a)
	ELSE IF ((z .GE. 0.5_rp).AND.(z .LT. 2.5_rp)) THEN
		a = 0.72_rp*(z + 1.0_rp)
		P_integral = P_integral + trapz(z,a)
		P_integral = P_integral + IntK(5.0_rp/3.0_rp,a)
	ELSE
		P_integral = IntK(5.0_rp/3.0_rp,z)
	END IF
END FUNCTION P_integral

! * * * * * * * * * * * * * * * !
! * * * * * FUNCTIONS * * * * * !
! * * * * * * * * * * * * * * * !


SUBROUTINE check_if_visible(X,V,threshold_angle,bool,angle,XC)
	IMPLICIT NONE
	REAL(rp), DIMENSION(3), INTENT(IN) :: X
	REAL(rp), DIMENSION(3), INTENT(IN) :: V
	REAL(rp), INTENT(IN) :: threshold_angle
	LOGICAL, INTENT(OUT) :: bool
	REAL(rp), INTENT(OUT) :: angle
	REAL(rp), DIMENSION(3), OPTIONAL, INTENT(OUT) :: XC
	REAL(rp), DIMENSION(3) :: vec
	REAL(rp) :: a, b, c, ciw, dis, disiw
	REAL(rp) :: sp, sn, s, psi


	a = V(1)**2 + V(2)**2
	b = 2.0_rp*(X(1)*V(1) + X(2)*V(2))
	c = X(1)**2 + X(2)**2 - cam%position(1)**2
	ciw = X(1)**2 + X(2)**2 - cam%Riw**2

	dis = b**2 - 4.0_rp*a*c
	disiw = b**2 - 4.0_rp*a*ciw
	
	if ((dis .LT. 0.0_rp).OR.(disiw .GE. 0.0_rp)) then
		bool = .FALSE. ! The particle is not visible
	else
		sp = 0.5_rp*(-b + SQRT(dis))/a
		sn = 0.5_rp*(-b - SQRT(dis))/a
		s = MAX(sp,sn)
		
		! Rotation angle along z-axis so that v is directed to the camera
		if (PRESENT(XC)) then
			XC(1) = X(1) + s*V(1)
			XC(2) = X(2) + s*V(2)
			XC(3) = X(3) + s*V(3)
			angle = ATAN2(XC(2),XC(1))
		else
			angle = ATAN2(X(2) + s*V(2),X(1) + s*V(1))
		end if
		if (angle.LT.0.0_rp) angle = angle + 2.0_rp*C_PI
	
		vec(1) = cam%position(1)*COS(angle) - X(1)
		vec(2) = cam%position(1)*SIN(angle) - X(2)
		vec(3) = cam%position(2) - X(3)

		vec = vec/SQRT(DOT_PRODUCT(vec,vec))
		
		psi = ACOS(DOT_PRODUCT(vec,V))

		if (psi.LE.threshold_angle) then
			bool = .TRUE. ! The particle is visible
		else
			bool = .FALSE. ! The particle is not visible
		end if
	end if

END SUBROUTINE check_if_visible


SUBROUTINE calculate_rotation_angles(X,bpa,apa)
	IMPLICIT NONE
	REAL(rp), DIMENSION(3), INTENT(IN) :: X
	LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: bpa
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: apa
	REAL(rp) :: R, D, psi
	REAL(rp) :: a, b, c, dis, xp, xn
	REAL(rp) :: xtmp, ytmp
	INTEGER :: ii,jj
	! bpa(:,:,1) -- > xp
	! bpa(:,:,2) -- > xn
		
	R = SQRT(SUM(X(1:2)**2))
	D = SQRT( (X(1) - cam%position(1))**2 + X(2) )
	psi = -ATAN2(X(3) - cam%position(2),D)

	bpa = .TRUE.
	
	do ii=1_idef,cam%num_pixels(1)
		a = 1.0_rp + TAN(ang%beta(ii))**2
		b = -2.0_rp*TAN(ang%beta(ii))**2*cam%position(1)
		c = (TAN(ang%beta(ii))*cam%position(1))**2 - R**2
		dis = b**2 - 4.0_rp*a*c
		
		if (dis.GT.0.0_rp) then
			do jj=1_idef,cam%num_pixels(2)

				if ((psi.GE.ang%psi(jj)).AND.(psi.LT.ang%psi(jj+1_idef))) then
					xp = 0.5_rp*(-b + SQRT(dis))/a
					xn = 0.5_rp*(-b - SQRT(dis))/a

					xtmp = xp - cam%position(1)
					ytmp = SQRT(R**2 - xp**2)

					! Check if particle is behind inner wall
					if ((ATAN2(ytmp,xtmp).GT.ang%threshold_angle).AND.(SQRT(xtmp**2+ytmp**2).GT.ang%threshold_radius)) then
						bpa(ii,jj,1) = .FALSE.
					else
						apa(ii,jj,1) = ATAN2(ytmp,xp)
						if (apa(ii,jj,1).LT.0.0_rp) apa(ii,jj,1) = apa(ii,jj,1)	+ 2.0_rp*C_PI
					end if

					xtmp = xn - cam%position(1)
					ytmp = SQRT(R**2 - xn**2)

					! Check if particle is behind inner wall
					if ((ATAN2(ytmp,xtmp).GT.ang%threshold_angle).AND.(SQRT(xtmp**2+ytmp**2).GT.ang%threshold_radius)) then
						bpa(ii,jj,2) = .FALSE.
					else
						apa(ii,jj,2) = ATAN2(ytmp,xn)
						if (apa(ii,jj,2).LT.0.0_rp) apa(ii,jj,2) = apa(ii,jj,2)	+ 2.0_rp*C_PI
					end if
				else ! Not in pixel (ii,jj)
					bpa(ii,jj,:) = .FALSE.
				end if ! Check if in pixel (ii,jj)

			end do ! NY
		else ! no real solutions
			bpa(ii,:,:) = .FALSE.
		end if ! Checking discriminant
	end do !! NX
END SUBROUTINE calculate_rotation_angles


SUBROUTINE spectral_angular_density(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	CHARACTER(MAX_STRING_LENGTH) :: var_name
	LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: bool_pixel_array
	LOGICAL :: bool
	REAL(rp), DIMENSION(3) :: binorm, n, nperp
	REAL(rp), DIMENSION(3) :: X, V, B, E, XC
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: angle_pixel_array
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: np_angular_pixel
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn_angular_pixel
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: np_lambda_pixel
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn_lambda_pixel
	REAL(rp) :: q, m, k, u, g, l, threshold_angle
	REAL(rp) :: psi, chi, beta, theta, Psyn_tmp
	REAL(rp) :: area, r, photon_energy
	REAL(rp) :: angle, clockwise
	REAL(rp) :: units
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Psyn_send_buffer,Psyn_receive_buffer, np_send_buffer, np_receive_buffer
	REAL(rp) :: lc, zeta
	INTEGER :: ii,jj,ll,ss,pp
    INTEGER :: numel, mpierr


	ALLOCATE(bool_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)
	ALLOCATE(angle_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)

	ALLOCATE(np_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))
	ALLOCATE(Psyn_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))

	ALLOCATE(np_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))
	ALLOCATE(Psyn_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))

	np_angular_pixel = 0.0_rp
	Psyn_angular_pixel = 0.0_rp

	np_lambda_pixel = 0.0_rp
	Psyn_lambda_pixel = 0.0_rp
	
	area = C_PI*(0.5_rp*cam%aperture)**2 ! In m^2

	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL FIRSTPRIVATE(q,m,area) PRIVATE(binorm,n,nperp,X,XC,V,B,E,&
!$OMP& bool_pixel_array,angle_pixel_array,k,u,g,l,threshold_angle,theta,&
!$OMP& psi,chi,beta,Psyn_tmp,bool,angle,clockwise,ii,jj,ll,pp,r,photon_energy,&
!$OMP& lc,zeta)&
!$OMP& SHARED(params,spp,ss,Psyn_angular_pixel,np_angular_pixel,np_lambda_pixel,Psyn_lambda_pixel)
!$OMP DO
		do pp=1_idef,spp(ss)%ppp
			if ( spp(ss)%vars%flag(pp) .EQ. 1_idef ) then
				V = spp(ss)%vars%V(:,pp)*params%cpp%velocity
				X = spp(ss)%vars%X(:,pp)*params%cpp%length
				g = spp(ss)%vars%g(pp)
				B = spp(ss)%vars%B(:,pp)*params%cpp%Bo
				E = spp(ss)%vars%E(:,pp)*params%cpp%Eo

				binorm = cross(V,E) + cross(V,cross(V,B))
		
				u = SQRT(DOT_PRODUCT(V,V))
				k = q*SQRT(DOT_PRODUCT(binorm,binorm))/(spp(ss)%vars%g(pp)*m*u**3)
				k = k/1.0E2_rp ! Now in cm^-1 (CGS)

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

				binorm = binorm/SQRT(DOT_PRODUCT(binorm,binorm))

				threshold_angle = (1.5_rp*k*cam%lambda_max/C_PI)**(1.0_rp/3.0_rp) ! In radians

				call check_if_visible(X,V/u,threshold_angle,bool,angle)	
			
				if (bool.EQV..TRUE.) then

					X(1:2) = clockwise_rotation(X(1:2),angle)
					V(1:2) = clockwise_rotation(V(1:2),angle)
					binorm(1:2) = clockwise_rotation(binorm(1:2),angle)

					call calculate_rotation_angles(X,bool_pixel_array,angle_pixel_array)

					clockwise = ATAN2(X(2),X(1))
					if (clockwise.LT.0.0_rp) clockwise = clockwise + 2.0_rp*C_PI

					do ii=1_idef,cam%num_pixels(1) ! NX
						do jj=1_idef,cam%num_pixels(2) ! NY
							
							if (bool_pixel_array(ii,jj,1)) then
								angle = angle_pixel_array(ii,jj,1) - clockwise

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								r = SQRT(DOT_PRODUCT(n,n))
								n = n/r

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								do ll=1_idef,cam%Nlambda ! Nlambda
									l = cam%lambda(ll)
									photon_energy = CGS_h*CGS_C/l

									if (theta .LE. threshold_angle) then
										zeta = lc/cam%lambda(ll)

										Psyn_tmp = &
										(4.0_rp*C_PI/SQRT(3.0_rp))*(CGS_C*CGS_E**2)*P_integral(zeta)/(g**2*cam%lambda(ii)**3)

										Psyn_lambda_pixel(ii,jj,ll,ss) = Psyn_lambda_pixel(ii,jj,ll,ss) + Psyn_tmp/photon_energy
										np_lambda_pixel(ii,jj,ll,ss) = np_lambda_pixel(ii,jj,ll,ss) + 1.0_rp
									end if

									if ((chi.LT.chic(g,k,l)).AND.(psi.LT.psic(k,l))) then
										Psyn_tmp = Psyn(g,psi,k,l,chi)
										if (Psyn_tmp.GT.0.0_rp) then
											Psyn_angular_pixel(ii,jj,ll,ss) = Psyn_angular_pixel(ii,jj,ll,ss) &
																			+ Psyn_tmp/photon_energy
											np_angular_pixel(ii,jj,ll,ss) = np_angular_pixel(ii,jj,ll,ss) + 1.0_rp
										end if
									end if
								end do ! Nlambda
							end if

							if (bool_pixel_array(ii,jj,2)) then
								angle = angle_pixel_array(ii,jj,2) - clockwise

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								n = n/SQRT(DOT_PRODUCT(n,n))

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								do ll=1_idef,cam%Nlambda ! Nlambda
									l = cam%lambda(ll)
									photon_energy = CGS_h*CGS_C/l

									if (theta .LE. threshold_angle) then
										zeta = lc/cam%lambda(ll)

										Psyn_tmp = &
										(4.0_rp*C_PI/SQRT(3.0_rp))*(CGS_C*CGS_E**2)*P_integral(zeta)/(g**2*cam%lambda(ii)**3)

										Psyn_lambda_pixel(ii,jj,ll,ss) = Psyn_lambda_pixel(ii,jj,ll,ss) + Psyn_tmp/photon_energy
										np_lambda_pixel(ii,jj,ll,ss) = np_lambda_pixel(ii,jj,ll,ss) + 1.0_rp
									end if

									if ((chi.LT.chic(g,k,l)).AND.(psi.LT.psic(k,l))) then
										Psyn_tmp = Psyn(g,psi,k,l,chi)
										if (Psyn_tmp.GT.0.0_rp) then
											Psyn_angular_pixel(ii,jj,ll,ss) = Psyn_angular_pixel(ii,jj,ll,ss) &
																			+ Psyn_tmp/photon_energy
											np_angular_pixel(ii,jj,ll,ss) = np_angular_pixel(ii,jj,ll,ss) + 1.0_rp
										end if
									end if
								end do ! Nlambda
							end if

						end do ! NY
					end do ! NX


				end if ! check if bool == TRUE
			end if ! if confined
		end do ! particles
!$OMP END DO
!$OMP END PARALLEL
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn has units of (erg/s)(cm^-1) 
!	or (photons/s)(cm^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

!	units = 1.0E-5_rp ! (Watts)(m^-1) Use these units when Psyn is in (erg/s*cm)
	units = 1.0E2_rp ! (Photons/s)(m^-1) Use these units when Psyn is in (Photons/s)(cm^-1)

	if (params%mpi_params%nmpi.GT.1_idef) then 
		numel = cam%num_pixels(1)*cam%num_pixels(2)*cam%Nlambda*params%num_species

		ALLOCATE(Psyn_send_buffer(numel))
		ALLOCATE(Psyn_receive_buffer(numel))
		ALLOCATE(np_send_buffer(numel))
		ALLOCATE(np_receive_buffer(numel))

		Psyn_send_buffer = RESHAPE(Psyn_angular_pixel,(/numel/))
		CALL MPI_REDUCE(Psyn_send_buffer,Psyn_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		np_send_buffer = RESHAPE(np_angular_pixel,(/numel/))
		CALL MPI_REDUCE(np_send_buffer,np_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_angular_pixel = RESHAPE(Psyn_receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species/))
		    np_angular_pixel = RESHAPE(np_receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species/))

			var_name = 'np_angular_pixel'
			call save_snapshot_var(params,np_angular_pixel,var_name)

			var_name = 'Psyn_angular_pixel'
			Psyn_angular_pixel = units*Psyn_angular_pixel
			call save_snapshot_var(params,Psyn_angular_pixel,var_name)
		end if


		Psyn_send_buffer = RESHAPE(Psyn_lambda_pixel,(/numel/))
		CALL MPI_REDUCE(Psyn_send_buffer,Psyn_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		np_send_buffer = RESHAPE(np_lambda_pixel,(/numel/))
		CALL MPI_REDUCE(np_send_buffer,np_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_lambda_pixel = RESHAPE(Psyn_receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species/))
		    np_lambda_pixel = RESHAPE(np_receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species/))

			var_name = 'np_lambda_pixel'
			call save_snapshot_var(params,np_lambda_pixel,var_name)

			var_name = 'Psyn_lambda_pixel'
			Psyn_lambda_pixel = units*Psyn_lambda_pixel
			call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
		end if

		DEALLOCATE(Psyn_send_buffer)
		DEALLOCATE(Psyn_receive_buffer)
		DEALLOCATE(np_send_buffer)
		DEALLOCATE(np_receive_buffer)

	    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	else
		var_name = 'np_angular_pixel'
		call save_snapshot_var(params,np_angular_pixel,var_name)

		var_name = 'np_lambda_pixel'
		call save_snapshot_var(params,np_lambda_pixel,var_name)

		var_name = 'Psyn_angular_pixel'
		Psyn_angular_pixel = units*Psyn_angular_pixel
		call save_snapshot_var(params,Psyn_angular_pixel,var_name)

		var_name = 'Psyn_lambda_pixel'
		Psyn_lambda_pixel = units*Psyn_lambda_pixel
		call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
	end if

	DEALLOCATE(bool_pixel_array)
	DEALLOCATE(angle_pixel_array)

	DEALLOCATE(np_angular_pixel)
    DEALLOCATE(Psyn_angular_pixel)

	DEALLOCATE(np_lambda_pixel)
    DEALLOCATE(Psyn_lambda_pixel)
END SUBROUTINE spectral_angular_density


SUBROUTINE spectral_density(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	CHARACTER(MAX_STRING_LENGTH) :: var_name
	REAL(rp), DIMENSION(3) :: binorm
	REAL(rp), DIMENSION(3) :: X, V, B, E
	REAL(rp), DIMENSION(:), ALLOCATABLE :: P
	REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: np
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PTot
	REAL(rp) :: R, Z
	REAL(rp) :: q, m, k, u, g, lc
	REAL(rp) :: photon_energy
	INTEGER :: ii,jj,ll,ss,pp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Psyn_send_buffer,Psyn_receive_buffer
	REAL(rp), DIMENSION(:), ALLOCATABLE :: np_send_buffer,np_receive_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: PTot_send_buffer,PTot_receive_buffer
    INTEGER :: numel, mpierr
	REAL(rp) :: units

	ALLOCATE(np(pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species))
	ALLOCATE(Psyn(pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species))
	ALLOCATE(PTot(pplane%grid_dims(1),pplane%grid_dims(2),params%num_species))
	ALLOCATE(P(cam%Nlambda))
	ALLOCATE(zeta(cam%Nlambda))

	np = 0.0_rp
	Psyn = 0.0_rp
	P = 0.0_rp
	PTot = 0.0_rp
	
	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL FIRSTPRIVATE(q,m) PRIVATE(binorm,X,V,B,E,&
!$OMP& k,u,g,lc,ii,jj,ll,pp,photon_energy,zeta,P,R,Z)&
!$OMP& SHARED(params,spp,ss,Psyn,PTot,np)
!$OMP DO
		do pp=1_idef,spp(ss)%ppp
			if ( spp(ss)%vars%flag(pp) .EQ. 1_idef ) then
				V = spp(ss)%vars%V(:,pp)*params%cpp%velocity
				X = spp(ss)%vars%X(:,pp)*params%cpp%length
				g = spp(ss)%vars%g(pp)
				B = spp(ss)%vars%B(:,pp)*params%cpp%Bo
				E = spp(ss)%vars%E(:,pp)*params%cpp%Eo

				binorm = cross(V,E) + cross(V,cross(V,B))
		
				u = SQRT(DOT_PRODUCT(V,V))
				k = q*SQRT(DOT_PRODUCT(binorm,binorm))/(spp(ss)%vars%g(pp)*m*u**3)
				k = k/1.0E2_rp ! Now in cm^-1 (CGS)

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3)
				zeta = lc/cam%lambda

!				jj = cam%Nlambda - MINLOC(zeta,1,mask=zeta.GT.2.5_rp) + 1_idef

				do ii=1_idef,cam%Nlambda
					P(ii) = (4.0_rp*C_PI/SQRT(3.0_rp))*(CGS_C*CGS_E**2)*P_integral(zeta(ii))/(g**2*cam%lambda(ii)**3)
				end do

				R = SQRT(SUM(X(1:2)**2))
				Z = X(3)
				
				ii = FLOOR((R - pplane%Rmin)/pplane%DR) + 1_idef
				jj = FLOOR((Z + ABS(pplane%Zmin))/pplane%DZ) + 1_idef

				Psyn(ii,jj,:,ss) = Psyn(ii,jj,:,ss) + P
				np(ii,jj,:,ss) = np(ii,jj,:,ss) + 1_idef
				PTot(ii,jj,ss) = PTot(ii,jj,ss) + spp(ss)%vars%Prad(pp);

			end if ! if confined
		end do ! particles
!$OMP END DO
!$OMP END PARALLEL
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn has units of (erg/s)(cm^-1) 
!	or (photons/s)(cm^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	if (params%mpi_params%nmpi.GT.1_idef) then 
		numel = pplane%grid_dims(1)*pplane%grid_dims(2)*cam%Nlambda*params%num_species

		ALLOCATE(Psyn_send_buffer(numel))
		ALLOCATE(Psyn_receive_buffer(numel))
		ALLOCATE(np_send_buffer(numel))
		ALLOCATE(np_receive_buffer(numel))

		Psyn_send_buffer = RESHAPE(Psyn,(/numel/))
		CALL MPI_REDUCE(Psyn_send_buffer,Psyn_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		np_send_buffer = RESHAPE(np,(/numel/))
		CALL MPI_REDUCE(np_send_buffer,np_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		numel = pplane%grid_dims(1)*pplane%grid_dims(2)*params%num_species

		ALLOCATE(PTot_send_buffer(numel))
		ALLOCATE(PTot_receive_buffer(numel))

		PTot_send_buffer = RESHAPE(PTot,(/numel/))
		CALL MPI_REDUCE(PTot_send_buffer,PTot_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn = RESHAPE(Psyn_receive_buffer,(/pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species/))

			units = 1.0E-5_rp ! (Watts)(m^-1) Use these units when Psyn is in (erg/s*cm)
!			units = 1.0E2_rp ! (Photons/s)(m^-1) Use these units when Psyn is in (Photons/s)(cm^-1)
			Psyn = units*Psyn

		    np = RESHAPE(np_receive_buffer,(/pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species/))

			var_name = 'np_pplane'
		    call save_snapshot_var(params,np,var_name)

			var_name = 'Psyn_pplane'
	    	call save_snapshot_var(params,Psyn,var_name)

			units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
			PTot = units*PTot ! (Watts)
			var_name = 'PTot_pplane'
			call save_snapshot_var(params,PTot,var_name)
		end if

		DEALLOCATE(Psyn_send_buffer)
		DEALLOCATE(Psyn_receive_buffer)
		DEALLOCATE(PTot_send_buffer)
		DEALLOCATE(PTot_receive_buffer)
		DEALLOCATE(np_send_buffer)
		DEALLOCATE(np_receive_buffer)

	    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	else
		units = 1.0E-5_rp ! (Watts)(m^-1) Use these units when Psyn is in (erg/s*cm)
!		units = 1.0E2_rp ! (Photons/s)(m^-1) Use these units when Psyn is in (Photons/s)(cm^-1)
		Psyn = units*Psyn

		var_name = 'np_pplane'
		call save_snapshot_var(params,np,var_name)

		var_name = 'Psyn_pplane'
	    call save_snapshot_var(params,Psyn,var_name)

		units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
		PTot = units*PTot ! (Watts)
		var_name = 'PTot_pplane'
		call save_snapshot_var(params,PTot,var_name)
	end if

	DEALLOCATE(np)
    DEALLOCATE(Psyn)
	DEALLOCATE(PTot)
	DEALLOCATE(P)
	DEALLOCATE(zeta)
END SUBROUTINE spectral_density


SUBROUTINE synthetic_camera(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp

	write(6,'("MPI:",I5," Synthetic camera diagnostic: ON!")') params%mpi_params%rank

	call spectral_angular_density(params,spp)

	call spectral_density(params,spp)

	write(6,'("MPI:",I5," Synthetic camera diagnostic: OFF!")') params%mpi_params%rank
END SUBROUTINE synthetic_camera


SUBROUTINE save_synthetic_camera_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	REAL(rp) :: units

	if (params%mpi_params%rank .EQ. 0) then
		filename = TRIM(params%path_to_outputs) // "synthetic_camera.h5"
		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		gname = "synthetic_camera_params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/aperture"
		attr = "Aperture of the camera (m)"
		call save_to_hdf5(h5file_id,dset,cam%aperture,attr)

		dset = TRIM(gname) // "/Riw"
		attr = "Radial position of inner wall (m)"
		call save_to_hdf5(h5file_id,dset,cam%Riw,attr)

		dset = TRIM(gname) // "/focal_length"
		attr = "Focal length of the camera (m)"
		call save_to_hdf5(h5file_id,dset,cam%focal_length,attr)

		dset = TRIM(gname) // "/incline"
		attr = "Incline of camera in degrees"
		units = 180.0_rp/C_PI
		call save_to_hdf5(h5file_id,dset,units*cam%incline,attr)

		dset = TRIM(gname) // "/horizontal_angle_view"
		attr = "Horizontal angle of view in degrees"
		units = 180.0_rp/C_PI
		call save_to_hdf5(h5file_id,dset,units*cam%horizontal_angle_view,attr)

		dset = TRIM(gname) // "/vertical_angle_view"
		attr = "Vertical angle of view in degrees"
		units = 180.0_rp/C_PI
		call save_to_hdf5(h5file_id,dset,units*cam%vertical_angle_view,attr)

		dset = TRIM(gname) // "/lambda_min"
		attr = "Minimum wavelength (m)"
		units = 1.0E-2_rp
		call save_to_hdf5(h5file_id,dset,units*cam%lambda_min,attr)

		dset = TRIM(gname) // "/lambda_max"
		attr = "Minimum wavelength (m)"
		units = 1.0E-2_rp
		call save_to_hdf5(h5file_id,dset,units*cam%lambda_max,attr)

		dset = TRIM(gname) // "/Dlambda"
		attr = "Step between finite wavelengths (m)"
		units = 1.0E-2_rp
		call save_to_hdf5(h5file_id,dset,units*cam%Dlambda,attr)

		dset = TRIM(gname) // "/Nlambda"
		attr = "Number of finite wavelengths (m)"
		call save_to_hdf5(h5file_id,dset,cam%Nlambda,attr)

	    dset = TRIM(gname) // "/lambda"
		units = 1.0E-2_rp
	    call save_1d_array_to_hdf5(h5file_id,dset,units*cam%lambda)

	    dset = TRIM(gname) // "/num_pixels"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%num_pixels)

	    dset = TRIM(gname) // "/sensor_size"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%sensor_size)

	    dset = TRIM(gname) // "/position"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%position)

	    dset = TRIM(gname) // "/pixels_nodes_x"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%pixels_nodes_x)

	    dset = TRIM(gname) // "/pixels_nodes_y"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%pixels_nodes_y)

	    dset = TRIM(gname) // "/pixels_edges_x"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%pixels_edges_x)

	    dset = TRIM(gname) // "/pixels_edges_y"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%pixels_edges_y)

		call h5gclose_f(group_id, h5error)


		gname = "poloidal_plane_params"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

	    dset = TRIM(gname) // "/grid_dims"
	    call save_1d_array_to_hdf5(h5file_id,dset,pplane%grid_dims)

	    dset = TRIM(gname) // "/nodes_R"
	    call save_1d_array_to_hdf5(h5file_id,dset,pplane%nodes_R)

	    dset = TRIM(gname) // "/nodes_Z"
	    call save_1d_array_to_hdf5(h5file_id,dset,pplane%nodes_Z)

		call h5gclose_f(group_id, h5error)


		call h5fclose_f(h5file_id, h5error)
	end if

    if (params%mpi_params%rank.EQ.0_idef) then
	    filename = TRIM(params%path_to_outputs) //"synthetic_camera_snapshots.h5"
	    call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)
	    call h5fclose_f(h5file_id, h5error)
    end if
END SUBROUTINE save_synthetic_camera_params


SUBROUTINE save_snapshot(params,part,Psyn)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(IN) :: part
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(IN) :: Psyn
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	REAL(rp) :: units

	filename = TRIM(params%path_to_outputs) //"synthetic_camera_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it'
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
    
	dset = TRIM(gname) // "/time"
	attr = "Simulation time in secs"
	call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)

		dset = "np_angular_pixel"
		call save_array_to_hdf5(subgroup_id, dset, part(:,:,:,ss))

		dset = "Psyn_angular_pixel"
!		units = 1.0E-5_rp ! (Watts)(m^-1) Use these units when Psyn is in (erg/s*cm)
		units = 1.0E2_rp ! (Photons/s)(m^-1) Use these units when Psyn is in (Photons/s)(cm^-1)
		call save_array_to_hdf5(subgroup_id, dset, units*Psyn(:,:,:,ss))

		call h5gclose_f(subgroup_id, h5error)
	end do

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot


SUBROUTINE save_snapshot_var_3d(params,var,var_name)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN) :: var
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: var_name
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	REAL(rp) :: units
	LOGICAL :: object_exists

	filename = TRIM(params%path_to_outputs) //"synthetic_camera_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it' if it doesn't exist
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

	if (.NOT.object_exists) then
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		dset = TRIM(gname) // "/time"
		attr = "Simulation time in secs"
		call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
	else
		call h5gopen_f(h5file_id, TRIM(gname), group_id, h5error)
	end if

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(group_id,TRIM(subgname),object_exists,h5error)

		if (.NOT.object_exists) then
			call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)
		else
			call h5gopen_f(group_id, TRIM(subgname), subgroup_id, h5error)
		end if

		dset = TRIM(var_name)
		call save_array_to_hdf5(subgroup_id, dset, var(:,:,ss))

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_3d


SUBROUTINE save_snapshot_var_4d(params,var,var_name)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(IN) :: var
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: var_name
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	REAL(rp) :: units
	LOGICAL :: object_exists

	filename = TRIM(params%path_to_outputs) //"synthetic_camera_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it' if it doesn't exist
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

	if (.NOT.object_exists) then
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		dset = TRIM(gname) // "/time"
		attr = "Simulation time in secs"
		call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
	else
		call h5gopen_f(h5file_id, TRIM(gname), group_id, h5error)
	end if

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(group_id,TRIM(subgname),object_exists,h5error)

		if (.NOT.object_exists) then
			call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)
		else
			call h5gopen_f(group_id, TRIM(subgname), subgroup_id, h5error)
		end if

		dset = TRIM(var_name)
		call save_array_to_hdf5(subgroup_id, dset, var(:,:,:,ss))

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_4d

END MODULE korc_synthetic_camera
