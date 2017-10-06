MODULE korc_synthetic_camera
	USE korc_types
	USE korc_constants
	USE korc_HDF5
    USE korc_hpc
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
		LOGICAL :: camera_on
		REAL(rp) :: start_at ! In seconds
!		REAL(rp) :: aperture ! Aperture of the camera (diameter of lens) in meters
		REAL(rp) :: Riw ! Radial position of inner wall
		INTEGER, DIMENSION(2) :: num_pixels ! Number of pixels (X,Y)
		REAL(rp), DIMENSION(2) :: sensor_size ! In meters (horizontal,vertical)
		REAL(rp) :: pixel_area ! Area of a single pixel of the camera sensor. This based on sensor_size and num_pixels.
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

		LOGICAL :: integrated_opt
		LOGICAL :: toroidal_sections
		LOGICAL :: photon_count
		INTEGER :: ntor_sections
	
		REAL(rp), DIMENSION(3) :: r
	END TYPE CAMERA

	TYPE, PRIVATE :: ANGLES
		REAL(rp), DIMENSION(:), ALLOCATABLE :: eta
		REAL(rp), DIMENSION(:), ALLOCATABLE :: beta
		REAL(rp), DIMENSION(:), ALLOCATABLE :: psi
		REAL(rp), DIMENSION(:), ALLOCATABLE :: ax
		REAL(rp), DIMENSION(:), ALLOCATABLE :: ay
		REAL(rp) :: ac

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
	  module procedure save_snapshot_var_1d,save_snapshot_var_2d,save_snapshot_var_3d,save_snapshot_var_4d
	END INTERFACE

	PRIVATE :: clockwise_rotation,anticlockwise_rotation,cross,check_if_visible,calculate_rotation_angles,&
				zeta,fx,arg,Po,P1,Psyn,chic,psic,save_synthetic_camera_params,besselk,&
				angular_density,spectral_density,integrated_SE_toroidal_sections,&
				IntK,nintegral_besselk,save_snapshot_var,trapz,integrated_SE_3D
	PUBLIC :: initialize_synthetic_camera,synthetic_camera

	CONTAINS

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !

SUBROUTINE test_analytical_formula()
	IMPLICIT NONE
	REAL(rp) :: z
	INTEGER :: ii

	do ii=1_idef,100
		z = 2.5_rp + REAL((ii-1_idef),rp)*0.5_rp
		write(6,*) IntK(5.0_rp/3.0_rp,z)
	end do
END SUBROUTINE test_analytical_formula


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
!	REAL(rp) :: aperture ! Aperture of the camera (diameter of lens) in meters
	REAL(rp) :: start_at ! in Seconds
	REAL(rp) :: Riw ! Radial position of inner wall
	INTEGER, DIMENSION(2) :: num_pixels ! Number of pixels (X,Y)
	REAL(rp), DIMENSION(2) :: sensor_size ! (horizontal,vertical)
	REAL(rp) :: focal_length
	REAL(rp), DIMENSION(2) :: position ! Position of camera (R,Z)
	REAL(rp) :: incline
	REAL(rp) :: lambda_min ! Minimum wavelength in cm
	REAL(rp) :: lambda_max ! Maximum wavelength in cm
	INTEGER :: Nlambda
	LOGICAL :: camera_on, integrated_opt, toroidal_sections, photon_count
	INTEGER :: ntor_sections
	REAL(rp) :: xmin, xmax, ymin, ymax, DX, DY
	INTEGER :: ii

	NAMELIST /SyntheticCamera/ camera_on,Riw,num_pixels,sensor_size,focal_length,&
								position,incline,lambda_min,lambda_max,Nlambda,integrated_opt,&
								toroidal_sections,ntor_sections,start_at,photon_count

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * * * * * * * * * * * * * *")')
		write(6,'("*  Initializing synthetic camera  *")')
	end if

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=SyntheticCamera)
	close(default_unit_open)

	!	write(*,nml=SyntheticCamera)
	
	cam%camera_on = camera_on
!	cam%aperture = aperture
	cam%start_at = start_at
	cam%Riw = Riw
	cam%num_pixels = num_pixels
	cam%sensor_size = sensor_size
	cam%pixel_area = PRODUCT(sensor_size/num_pixels)
	cam%focal_length = focal_length
	cam%position = position
	cam%incline = C_PI*incline/180.0_rp
	cam%horizontal_angle_view = ATAN2(0.5_rp*cam%sensor_size(1),cam%focal_length)
	cam%vertical_angle_view = ATAN2(0.5_rp*cam%sensor_size(2),cam%focal_length)
	cam%lambda_min = lambda_min ! In meters
	cam%lambda_max = lambda_max ! In meters
	cam%Nlambda = Nlambda
	cam%Dlambda = (cam%lambda_max - cam%lambda_min)/REAL(cam%Nlambda,rp)
	ALLOCATE(cam%lambda(cam%Nlambda))
	cam%photon_count = photon_count
	cam%integrated_opt = integrated_opt
	cam%toroidal_sections = toroidal_sections
	if (cam%toroidal_sections) then
		cam%ntor_sections = ntor_sections
	else
		cam%ntor_sections = 1_idef
	end if
	
	if (cam%camera_on) then
		cam%r = (/COS(0.5_rp*C_PI - cam%incline),-SIN(0.5_rp*C_PI - cam%incline),0.0_rp/)

		do ii=1_idef,cam%Nlambda
			cam%lambda(ii) = cam%lambda_min + REAL(ii-1_idef,rp)*cam%Dlambda
		end do

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

		ymin = -0.5_rp*cam%sensor_size(2)
		ymax = 0.5_rp*cam%sensor_size(2)
		DY = cam%sensor_size(2)/REAL(cam%num_pixels(2),rp)

		do ii=1_idef,cam%num_pixels(2)
			cam%pixels_nodes_y(ii) = ymin + 0.5_rp*DY + REAL(ii-1_idef,rp)*DY
		end do

		do ii=1_idef,cam%num_pixels(2)+1_idef
			cam%pixels_edges_y(ii) = ymin + REAL(ii-1_idef,rp)*DY
		end do
	
		! Initialize ang variables
		ALLOCATE(ang%eta(cam%num_pixels(1))) ! angle between main line of sight and a given pixel along the x-axis
		ALLOCATE(ang%beta(cam%num_pixels(1))) 
		ALLOCATE(ang%psi(cam%num_pixels(2)+1_idef)) ! angle between main line of sight and a given pixel along the y-axis

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

		ALLOCATE(ang%ax(cam%num_pixels(1))) ! angle between main line of sight and a given pixel along the x-axis
		ALLOCATE(ang%ay(cam%num_pixels(2))) ! angle between main line of sight and a given pixel along the y-axis

		do ii=1_idef,cam%num_pixels(1)
			ang%ax(ii) = ATAN2(cam%pixels_nodes_x(ii),cam%focal_length)
		end do

		do ii=1_idef,cam%num_pixels(2)
			ang%ay(ii) = ATAN2(cam%pixels_nodes_y(ii),cam%focal_length)
		end do

		if (cam%incline.GT.0.5_rp*C_PI) then
			ang%ac = cam%incline - 0.5_rp*C_PI
		else
			ang%ac = 0.5_rp*C_PI - cam%incline
		end if

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("*     Synthetic camera ready!     *")')
			write(6,'("* * * * * * * * * * * * * * * * * *",/)')
		end if


		! Initialize poloidal plane parameters
	
		pplane%grid_dims = num_pixels
		ALLOCATE(pplane%nodes_R(pplane%grid_dims(1)))
		ALLOCATE(pplane%nodes_Z(pplane%grid_dims(2)))

		! * * * * * * * ALL IN METERS * * * * * * * 

		IF (TRIM(params%plasma_model) .EQ. 'ANALYTICAL') THEN
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
	end if
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
	
	Po = -(C_C*C_E**2/(SQRT(3.0_rp)*C_E0*k*(l*g)**4))*(1.0_rp + (g*p)**2)**2
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

	IntK = (C_PI/SQRT(2.0_rp))*(1.0_rp - 0.25_rp*(4.0_rp*v**2 - 1.0_rp))*ERFC(SQRT(x))&
			 + 0.25_rp*(4.0_rp*v**2 - 1.0_rp)*SQRT(0.5_rp*C_PI/x)*EXP(-x)
END FUNCTION IntK


FUNCTION nintegral_besselk(a,b)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: a
	REAL(rp), INTENT(IN) :: b
	REAL(rp) :: nintegral_besselk
	REAL(rp) :: Iold, Inew, rerr
	REAL(rp) :: v,h,z
	INTEGER :: ii,jj,npoints
	LOGICAL :: flag
	
	v = 5.0_rp/3.0_rp

	h = b - a
	nintegral_besselk = 0.5*(besselk(v,a) + besselk(v,b))
	
	ii = 1_idef
	flag = .TRUE.
	do while (flag)
		Iold = nintegral_besselk*h

		ii = ii + 1_idef
		npoints = 2_idef**(ii-2_idef)
		h = 0.5_rp*(b-a)/REAL(npoints,rp)

		do jj=1_idef,npoints
			z = a + h + 2.0_rp*(REAL(jj,rp) - 1.0_rp)*h
			nintegral_besselk = nintegral_besselk + besselk(v,z)
		end do

		Inew = nintegral_besselk*h

		rerr = ABS((Inew - Iold)/Iold)

		flag = .NOT.(rerr.LT.Tol)
	end do
	nintegral_besselk = Inew
END FUNCTION nintegral_besselk

SUBROUTINE P_integral(z,P)
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: P
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: z
	REAL(rp) :: a
	INTEGER :: ll,ss

	ss = SIZE(z)

	P = 0.0_rp

	do ll=1_idef,ss
		IF (z(ll) .LT. 0.5_rp) THEN
			a = (2.16_rp/2.0_rp**(2.0_rp/3.0_rp))*z(ll)**(1.0_rp/3.0_rp)
			P(ll) = nintegral_besselk(z(ll),a) + IntK(5.0_rp/3.0_rp,a)
		ELSE IF ((z(ll) .GE. 0.5_rp).AND.(z(ll) .LT. 2.5_rp)) THEN
			a = 0.72_rp*(z(ll) + 1.0_rp)
			P(ll) = nintegral_besselk(z(ll),a) + IntK(5.0_rp/3.0_rp,a)
		ELSE
			P(ll) = IntK(5.0_rp/3.0_rp,z(ll))
		END IF
	end do
END SUBROUTINE P_integral


FUNCTION trapz(x,f)
	IMPLICIT NONE
	REAL(rp), DIMENSION(:), INTENT(IN) :: x
	REAL(rp), DIMENSION(:), INTENT(IN) :: f
	REAL(rp) :: trapz
	INTEGER :: N

	N = SIZE(x)

	trapz = 0.5_rp*SUM( (x(2:N) - x(1:N-1))*(f(1:N-1) + f(2:N)) )
END FUNCTION trapz


FUNCTION rad2deg(t)
	REAL(rp), INTENT(IN) :: t
	REAL(rp) :: rad2deg

	rad2deg = t*180.0_rp/C_PI

END FUNCTION rad2deg

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


SUBROUTINE is_visible(X,V,threshold_angle,bool,ii,jj)
	IMPLICIT NONE
	REAL(rp), DIMENSION(3), INTENT(IN) :: X
	REAL(rp), DIMENSION(3), INTENT(IN) :: V
	REAL(rp), INTENT(IN) :: threshold_angle
	LOGICAL, INTENT(OUT) :: bool
	INTEGER, INTENT(OUT) :: ii
	INTEGER, INTENT(OUT) :: jj
	REAL(rp), DIMENSION(3) :: vec
	REAL(rp), DIMENSION(3) :: n
	REAL(rp) :: r
	REAL(rp) :: psi
	REAL(rp) :: t,ax,ay

	vec(1) = cam%position(1) - X(1)
	vec(2) = -X(2)
	vec(3) = cam%position(2) - X(3)

	r = SQRT(DOT_PRODUCT(vec,vec))
	n = vec/r
		
	psi = ACOS(DOT_PRODUCT(n,V))

	n = (/vec(1),vec(2),0.0_rp/)
	n = n/SQRT(DOT_PRODUCT(n,n))

	if (psi.LE.threshold_angle) then
		bool = .TRUE. ! The particle is visible

		t =  ACOS(DOT_PRODUCT(n,(/1.0_rp,0.0_rp,0.0_rp/)))

		if (cam%incline.GT.0.5_rp*C_PI) then
			if (t.GT.ang%ac) then
				ax = -ACOS(DOT_PRODUCT(n,cam%r))
			else
				ax = ACOS(DOT_PRODUCT(n,cam%r))
			end if
		else
			if (t.GT.ang%ac) then
				ax = ACOS(DOT_PRODUCT(n,cam%r))
			else
				ax = -ACOS(DOT_PRODUCT(n,cam%r))
			end if
		end if

		ay = -ASIN(vec(3)/r)

		ii = MINLOC(ABS(ax - ang%ax),1)
		jj = MINLOC(ABS(ay - ang%ay),1)

		if ((ii.GT.cam%num_pixels(1)).OR.(jj.GT.cam%num_pixels(2))) bool = .FALSE.
	else
		bool = .FALSE. ! The particle is not visible
	end if
END SUBROUTINE is_visible


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
	D = SQRT( (cam%position(1) - X(1))**2 + X(2)**2 )
	psi = -ATAN2(cam%position(2) - X(3),D)

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


SUBROUTINE angular_density(params,spp)
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
	REAL(rp) :: psi, chi, beta, theta
	REAL(rp), DIMENSION(:), ALLOCATABLE :: P_lambda
	REAL(rp) :: P_angular
	REAL(rp), DIMENSION(:), ALLOCATABLE :: photon_energy
	REAL(rp) :: r,solid_angle
	REAL(rp) :: angle, clockwise
	REAL(rp) :: units
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Psyn_send_buffer,Psyn_receive_buffer, np_send_buffer, np_receive_buffer
	REAL(rp) :: lc
	REAL(rp) :: dSA ! Element of solid angle
	REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
	INTEGER :: ii,jj,ll,ss,pp
    INTEGER :: numel, mpierr


	ALLOCATE(bool_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)
	ALLOCATE(angle_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)

	ALLOCATE(np_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))
	ALLOCATE(Psyn_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))

	ALLOCATE(np_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))
	ALLOCATE(Psyn_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%Nlambda,params%num_species))

	ALLOCATE(P_lambda(cam%Nlambda))

	ALLOCATE(zeta(cam%Nlambda))

	ALLOCATE(photon_energy(cam%Nlambda))

	np_angular_pixel = 0.0_rp
	Psyn_angular_pixel = 0.0_rp

	np_lambda_pixel = 0.0_rp
	Psyn_lambda_pixel = 0.0_rp

	P_lambda = 0.0_rp
	zeta = 0.0_rp

	photon_energy = C_h*C_C/cam%lambda

	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL DO FIRSTPRIVATE(q,m,photon_energy) PRIVATE(binorm,n,nperp,X,XC,V,B,E,&
!$OMP& bool_pixel_array,angle_pixel_array,k,u,g,l,threshold_angle,theta,&
!$OMP& psi,chi,beta,P_lambda,bool,angle,clockwise,ii,jj,ll,pp,r,lc,zeta,solid_angle,dSA)&
!$OMP& SHARED(params,spp,ss,Psyn_angular_pixel,np_angular_pixel,np_lambda_pixel,Psyn_lambda_pixel)
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

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

				zeta = lc/cam%lambda

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

								dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								if (theta .LE. threshold_angle) then
									call P_integral(zeta,P_lambda)

									P_lambda = (C_C*C_E**2)*P_lambda/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)

									P_lambda = dSA*P_lambda/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))

									if (cam%photon_count) then
										Psyn_lambda_pixel(ii,jj,:,ss) = Psyn_lambda_pixel(ii,jj,:,ss) + P_lambda/photon_energy
									else
										Psyn_lambda_pixel(ii,jj,:,ss) = Psyn_lambda_pixel(ii,jj,:,ss) + P_lambda
									end if
									np_lambda_pixel(ii,jj,:,ss) = np_lambda_pixel(ii,jj,:,ss) + 1.0_rp
								end if

								do ll=1_idef,cam%Nlambda ! Nlambda
									if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
										P_angular = Psyn(g,psi,k,cam%lambda(ll),chi)
										if (P_angular.GT.0.0_rp) then
											P_angular = dSA*P_angular
											if (cam%photon_count) then
												Psyn_angular_pixel(ii,jj,ll,ss) = Psyn_angular_pixel(ii,jj,ll,ss) &
																				+ P_angular/photon_energy(ll)
											else
												Psyn_angular_pixel(ii,jj,ll,ss) = Psyn_angular_pixel(ii,jj,ll,ss) + P_angular
											end if
											np_angular_pixel(ii,jj,ll,ss) = np_angular_pixel(ii,jj,ll,ss) + 1.0_rp
										end if
									end if
								end do ! Nlambda
							end if

							if (bool_pixel_array(ii,jj,2)) then
								angle = angle_pixel_array(ii,jj,2) - clockwise

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								r = SQRT(DOT_PRODUCT(n,n))
								n = n/r

								dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								if (theta .LE. threshold_angle) then
									call P_integral(zeta,P_lambda)

									P_lambda = (C_C*C_E**2)*P_lambda/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)

									P_lambda = dSA*P_lambda/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))

									if (cam%photon_count) then
										Psyn_lambda_pixel(ii,jj,:,ss) = Psyn_lambda_pixel(ii,jj,:,ss) + P_lambda/photon_energy
									else
										Psyn_lambda_pixel(ii,jj,:,ss) = Psyn_lambda_pixel(ii,jj,:,ss) + P_lambda
									end if
									np_lambda_pixel(ii,jj,:,ss) = np_lambda_pixel(ii,jj,:,ss) + 1.0_rp
								end if

								do ll=1_idef,cam%Nlambda ! Nlambda
									if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
										P_angular = Psyn(g,psi,k,cam%lambda(ll),chi)
										if (P_angular.GT.0.0_rp) then
											P_angular = dSA*P_angular
											if (cam%photon_count) then
												Psyn_angular_pixel(ii,jj,ll,ss) = Psyn_angular_pixel(ii,jj,ll,ss) &
																				+ P_angular/photon_energy(ll)
											else
												Psyn_angular_pixel(ii,jj,ll,ss) = Psyn_angular_pixel(ii,jj,ll,ss) + P_angular
											end if
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
!$OMP END PARALLEL DO
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn has units of Watts/m 
!	or (photons/s)(m^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	units = 1.0_rp

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
		call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
	end if

	DEALLOCATE(bool_pixel_array)
	DEALLOCATE(angle_pixel_array)

	DEALLOCATE(np_angular_pixel)
    DEALLOCATE(Psyn_angular_pixel)

	DEALLOCATE(np_lambda_pixel)
    DEALLOCATE(Psyn_lambda_pixel)

	DEALLOCATE(P_lambda)

	DEALLOCATE(zeta)

	DEALLOCATE(photon_energy)
END SUBROUTINE angular_density


SUBROUTINE integrated_angular_density(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	CHARACTER(MAX_STRING_LENGTH) :: var_name
	LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: bool_pixel_array
	LOGICAL :: bool
	REAL(rp), DIMENSION(3) :: binorm, n, nperp
	REAL(rp), DIMENSION(3) :: X, V, B, E, XC
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: angle_pixel_array
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: np_angular_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Psyn_angular_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: np_lambda_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Psyn_lambda_pixel
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: P_lambda_pixel
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: P_angular_pixel
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: P_l_pixel
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: P_a_pixel
	REAL(rp), DIMENSION(:), ALLOCATABLE :: np_pixel
	REAL(rp), DIMENSION(:), ALLOCATABLE :: P_lambda, P_angular
	REAL(rp) :: q, m, k, u, g, l, threshold_angle, threshold_angle_simple_model
	REAL(rp) :: psi, chi, beta, theta
	REAL(rp), DIMENSION(:), ALLOCATABLE :: photon_energy
	REAL(rp) :: r,solid_angle
	REAL(rp) :: angle, clockwise
	REAL(rp) :: units
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer, receive_buffer
	REAL(rp) :: lc
	REAL(rp) :: dSA ! Element of solid angle
	REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
	INTEGER :: ii,jj,ll,ss,pp
    INTEGER :: numel, mpierr


	ALLOCATE(bool_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)
	ALLOCATE(angle_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)

	ALLOCATE(np_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),params%num_species))
	ALLOCATE(Psyn_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),params%num_species))

	ALLOCATE(np_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),params%num_species))
	ALLOCATE(Psyn_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),params%num_species))

	ALLOCATE(P_lambda(cam%Nlambda))
	ALLOCATE(P_angular(cam%Nlambda))

	ALLOCATE(P_l_pixel(cam%Nlambda,params%num_species))
	ALLOCATE(P_a_pixel(cam%Nlambda,params%num_species))
	ALLOCATE(np_pixel(params%num_species))

	ALLOCATE(zeta(cam%Nlambda))

	ALLOCATE(photon_energy(cam%Nlambda))

	np_angular_pixel = 0.0_rp
	Psyn_angular_pixel = 0.0_rp

	np_lambda_pixel = 0.0_rp
	Psyn_lambda_pixel = 0.0_rp

	P_l_pixel = 0.0_rp
	P_a_pixel = 0.0_rp
	np_pixel = 0.0_rp

	zeta = 0.0_rp

	photon_energy = C_h*C_C/cam%lambda

	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL DO FIRSTPRIVATE(q,m,photon_energy) PRIVATE(binorm,n,nperp,X,XC,V,B,E,&
!$OMP& bool_pixel_array,angle_pixel_array,k,u,g,l,threshold_angle,threshold_angle_simple_model,theta,&
!$OMP& psi,chi,beta,bool,angle,clockwise,ii,jj,ll,pp,r,lc,zeta,P_lambda,P_angular,solid_angle,dSA)&
!$OMP& SHARED(params,spp,ss,Psyn_angular_pixel,np_angular_pixel,np_lambda_pixel,Psyn_lambda_pixel,P_l_pixel,P_a_pixel,np_pixel)
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

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

				zeta = lc/cam%lambda

				binorm = binorm/SQRT(DOT_PRODUCT(binorm,binorm))

				threshold_angle = (1.5_rp*k*cam%lambda_max/C_PI)**(1.0_rp/3.0_rp) ! In radians

				threshold_angle_simple_model = 1.0_rp/g

				np_pixel(ss) = np_pixel(ss) + 1.0_rp ! We count all the confined particles.

				call check_if_visible(X,V/u,MAXVAL((/threshold_angle,threshold_angle_simple_model/)),bool,angle)
			
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
								angle = angle_pixel_array(ii,jj,1) - clockwise ! Here, angle is modified w.r.t. check_if_visible.

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								r = SQRT(DOT_PRODUCT(n,n))
								n = n/r

								dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								P_lambda = 0.0_rp
								P_angular = 0.0_rp

								if (theta .LE. threshold_angle_simple_model) then
									call P_integral(zeta,P_lambda)

									P_lambda = (C_C*C_E**2)*P_lambda/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)
									np_lambda_pixel(ii,jj,ss) = np_lambda_pixel(ii,jj,ss) + 1.0_rp
								end if

								do ll=1_idef,cam%Nlambda ! Nlambda
									if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
										P_angular(ll) = Psyn(g,psi,k,cam%lambda(ll),chi)
										if (P_angular(ll).GT.0.0_rp) then
											np_angular_pixel(ii,jj,ss) = np_angular_pixel(ii,jj,ss) + 1.0_rp
										else
											P_angular(ll) = 0.0_rp
										end if
									end if
								end do ! Nlambda	

								P_lambda = dSA*P_lambda/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))
								P_angular = dSA*P_angular

								P_l_pixel(:,ss)	= P_l_pixel(:,ss) + P_lambda
								P_a_pixel(:,ss)	= P_a_pixel(:,ss) + P_angular

								if (cam%photon_count) then
									P_lambda = P_lambda/photon_energy
									P_angular = P_angular/photon_energy

									Psyn_lambda_pixel(ii,jj,ss) = Psyn_lambda_pixel(ii,jj,ss) + SUM(P_lambda)
									Psyn_angular_pixel(ii,jj,ss) = Psyn_angular_pixel(ii,jj,ss) + SUM(P_angular)
								else
									Psyn_lambda_pixel(ii,jj,ss) = Psyn_lambda_pixel(ii,jj,ss) + trapz(cam%lambda,P_lambda)
									Psyn_angular_pixel(ii,jj,ss) = Psyn_angular_pixel(ii,jj,ss) + trapz(cam%lambda,P_angular)
								end if
							end if

							if (bool_pixel_array(ii,jj,2)) then
								angle = angle_pixel_array(ii,jj,2) - clockwise

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								r = SQRT(DOT_PRODUCT(n,n))
								n = n/r

								dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								P_lambda = 0.0_rp
								P_angular = 0.0_rp

								if (theta .LE. threshold_angle_simple_model) then
									call P_integral(zeta,P_lambda)

									P_lambda = (C_C*C_E**2)*P_lambda/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)
									np_lambda_pixel(ii,jj,ss) = np_lambda_pixel(ii,jj,ss) + 1.0_rp
								end if

								do ll=1_idef,cam%Nlambda ! Nlambda
									if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
										P_angular(ll) = Psyn(g,psi,k,cam%lambda(ll),chi)
										if (P_angular(ll).GT.0.0_rp) then
											np_angular_pixel(ii,jj,ss) = np_angular_pixel(ii,jj,ss) + 1.0_rp
										else
											P_angular(ll) = 0.0_rp
										end if
									end if
								end do ! Nlambda	

								P_lambda = dSA*P_lambda/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))
								P_angular = dSA*P_angular

								P_l_pixel(:,ss)	= P_l_pixel(:,ss) + P_lambda
								P_a_pixel(:,ss)	= P_a_pixel(:,ss) + P_angular

								if (cam%photon_count) then
									P_lambda = P_lambda/photon_energy
									P_angular = P_angular/photon_energy

									Psyn_lambda_pixel(ii,jj,ss) = Psyn_lambda_pixel(ii,jj,ss) + SUM(P_lambda)
									Psyn_angular_pixel(ii,jj,ss) = Psyn_angular_pixel(ii,jj,ss) + SUM(P_angular)
								else
									Psyn_lambda_pixel(ii,jj,ss) = Psyn_lambda_pixel(ii,jj,ss) + trapz(cam%lambda,P_lambda)
									Psyn_angular_pixel(ii,jj,ss) = Psyn_angular_pixel(ii,jj,ss) + trapz(cam%lambda,P_angular)
								end if
							end if

						end do ! NY
					end do ! NX
				end if ! check if bool == TRUE

			end if ! if confined
		end do ! particles
!$OMP END PARALLEL DO
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn has units of Watts/m 
!	or (photons/s)(m^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	units = 1.0_rp

	if (params%mpi_params%nmpi.GT.1_idef) then 
		numel = cam%num_pixels(1)*cam%num_pixels(2)*params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_angular_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_angular_pixel = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),params%num_species/))

			var_name = 'Psyn_angular_pixel'
			call save_snapshot_var(params,Psyn_angular_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_angular_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_angular_pixel = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),params%num_species/))

			var_name = 'np_angular_pixel'
			call save_snapshot_var(params,np_angular_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)


		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_lambda_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_lambda_pixel = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),params%num_species/))

			var_name = 'Psyn_lambda_pixel'
			call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_lambda_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_lambda_pixel = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),params%num_species/))

			var_name = 'np_lambda_pixel'
			call save_snapshot_var(params,np_lambda_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		numel = params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = np_pixel
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_pixel = receive_buffer

	        var_name = 'np_pixel'
	        call save_snapshot_var(params,np_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		numel = cam%Nlambda*params%num_species
		
		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_a_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_a_pixel = RESHAPE(receive_buffer,(/cam%Nlambda,params%num_species/))

	        var_name = 'P_a_pixel'
	        call save_snapshot_var(params,P_a_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_l_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_l_pixel = RESHAPE(receive_buffer,(/cam%Nlambda,params%num_species/))

	        var_name = 'P_l_pixel'
	        call save_snapshot_var(params,P_l_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

	    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	else
		var_name = 'np_angular_pixel'
		call save_snapshot_var(params,np_angular_pixel,var_name)

		var_name = 'np_lambda_pixel'
		call save_snapshot_var(params,np_lambda_pixel,var_name)

		var_name = 'Psyn_angular_pixel'
		call save_snapshot_var(params,Psyn_angular_pixel,var_name)

		var_name = 'Psyn_lambda_pixel'
		call save_snapshot_var(params,Psyn_lambda_pixel,var_name)


		var_name = 'np_pixel'
		call save_snapshot_var(params,np_pixel,var_name)

		var_name = 'P_a_pixel'
		call save_snapshot_var(params,P_a_pixel,var_name)

		var_name = 'P_l_pixel'
		call save_snapshot_var(params,P_l_pixel,var_name)
	end if

	DEALLOCATE(bool_pixel_array)
	DEALLOCATE(angle_pixel_array)

	DEALLOCATE(np_angular_pixel)
    DEALLOCATE(Psyn_angular_pixel)

	DEALLOCATE(np_lambda_pixel)
    DEALLOCATE(Psyn_lambda_pixel)

	DEALLOCATE(P_lambda)
	DEALLOCATE(P_angular)

	DEALLOCATE(P_l_pixel)
	DEALLOCATE(P_a_pixel)
	DEALLOCATE(np_pixel)

	DEALLOCATE(zeta)

	DEALLOCATE(photon_energy)
END SUBROUTINE integrated_angular_density


SUBROUTINE integrated_SE_toroidal_sections(params,spp)
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
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_lambda_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_angular_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_l_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_a_pixel
	REAL(rp), DIMENSION(:), ALLOCATABLE :: np_pixel
	REAL(rp), DIMENSION(:), ALLOCATABLE :: P
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: P_lambda, P_angular
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: array3D
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: array2D
	REAL(rp) :: q, m, k, u, g, l, threshold_angle, threshold_angle_simple_model
	REAL(rp) :: psi, chi, beta, theta
	REAL(rp), DIMENSION(:), ALLOCATABLE :: photon_energy
	REAL(rp) :: r,solid_angle
	REAL(rp) :: angle, clockwise
	REAL(rp) :: units
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer, receive_buffer
	REAL(rp) :: lc
	REAL(rp) :: dSA ! Element of solid angle
	REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
	REAL(rp) :: Dtor
	INTEGER :: ii,jj,ll,ss,pp
	INTEGER :: itor
    INTEGER :: numel, mpierr


	ALLOCATE(bool_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)
	ALLOCATE(angle_pixel_array(cam%num_pixels(1),cam%num_pixels(2),2)) ! (NX,NY,2)

	ALLOCATE(np_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))
	ALLOCATE(Psyn_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))

	ALLOCATE(np_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))
	ALLOCATE(Psyn_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))

	ALLOCATE(P(cam%Nlambda))
	
	ALLOCATE(P_lambda(cam%Nlambda,cam%ntor_sections))
	ALLOCATE(P_angular(cam%Nlambda,cam%ntor_sections))

	ALLOCATE(P_l_pixel(cam%Nlambda,cam%ntor_sections,params%num_species))
	ALLOCATE(P_a_pixel(cam%Nlambda,cam%ntor_sections,params%num_species))
	ALLOCATE(np_pixel(params%num_species))

	ALLOCATE(zeta(cam%Nlambda))

	ALLOCATE(photon_energy(cam%Nlambda))

	np_angular_pixel = 0.0_rp
	Psyn_angular_pixel = 0.0_rp

	np_lambda_pixel = 0.0_rp
	Psyn_lambda_pixel = 0.0_rp

	P_l_pixel = 0.0_rp
	P_a_pixel = 0.0_rp
	np_pixel = 0.0_rp

	Dtor = 2.0_rp*C_PI/REAL(cam%ntor_sections,rp)

	zeta = 0.0_rp

	photon_energy = C_h*C_C/cam%lambda

	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL DO FIRSTPRIVATE(q,m,Dtor,photon_energy) PRIVATE(binorm,n,nperp,X,XC,V,B,E,&
!$OMP& bool_pixel_array,angle_pixel_array,k,u,g,l,threshold_angle,threshold_angle_simple_model,theta,&
!$OMP& psi,chi,beta,bool,angle,clockwise,ii,jj,ll,pp,r,lc,zeta,P_lambda,P_angular,itor,solid_angle,dSA,P)&
!$OMP& SHARED(params,spp,ss,Psyn_angular_pixel,np_angular_pixel,np_lambda_pixel,Psyn_lambda_pixel,P_l_pixel,P_a_pixel,np_pixel)
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

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

				zeta = lc/cam%lambda

				binorm = binorm/SQRT(DOT_PRODUCT(binorm,binorm))

				threshold_angle = (1.5_rp*k*cam%lambda_max/C_PI)**(1.0_rp/3.0_rp) ! In radians

				threshold_angle_simple_model = 1.0_rp/g

				np_pixel(ss) = np_pixel(ss) + 1.0_rp ! We count all the confined particles.

				call check_if_visible(X,V/u,MAXVAL((/threshold_angle,threshold_angle_simple_model/)),bool,angle)
			
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
								angle = angle_pixel_array(ii,jj,1) - clockwise ! Here, angle is modified w.r.t. check_if_visible.
								itor = floor(angle_pixel_array(ii,jj,1)/Dtor) + 1_idef

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								r = SQRT(DOT_PRODUCT(n,n))
								n = n/r

								dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								P_lambda = 0.0_rp
								P_angular = 0.0_rp

								if (theta .LE. threshold_angle_simple_model) then
									call P_integral(zeta,P)

									P_lambda(:,itor) = (C_C*C_E**2)*P/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)
									np_lambda_pixel(ii,jj,itor,ss) = np_lambda_pixel(ii,jj,itor,ss) + 1.0_rp
								end if

								do ll=1_idef,cam%Nlambda ! Nlambda
									if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
										P_angular(ll,itor) = Psyn(g,psi,k,cam%lambda(ll),chi)
										if (P_angular(ll,itor).GT.0.0_rp) then
											np_angular_pixel(ii,jj,itor,ss) = np_angular_pixel(ii,jj,itor,ss) + 1.0_rp
										else
											P_angular(ll,itor) = 0.0_rp
										end if
									end if
								end do ! Nlambda	

								P_lambda(:,itor) = dSA*P_lambda(:,itor)/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))
								P_angular(:,itor) = dSA*P_angular(:,itor)

								P_l_pixel(:,itor,ss) = P_l_pixel(:,itor,ss) + P_lambda(:,itor)
								P_a_pixel(:,itor,ss) = P_a_pixel(:,itor,ss) + P_angular(:,itor)

								if (cam%photon_count) then
									P_lambda(:,itor) = P_lambda(:,itor)/photon_energy
									P_angular(:,itor) = P_angular(:,itor)/photon_energy

									Psyn_lambda_pixel(ii,jj,itor,ss) = Psyn_lambda_pixel(ii,jj,itor,ss) + SUM(P_lambda(:,itor))
									Psyn_angular_pixel(ii,jj,itor,ss) = Psyn_angular_pixel(ii,jj,itor,ss) + SUM(P_angular(:,itor))
								else
									Psyn_lambda_pixel(ii,jj,itor,ss) = Psyn_lambda_pixel(ii,jj,itor,ss) + &
																trapz(cam%lambda,P_lambda(:,itor))
									Psyn_angular_pixel(ii,jj,itor,ss) = Psyn_angular_pixel(ii,jj,itor,ss) + &
																trapz(cam%lambda,P_angular(:,itor))
								end if
							end if

							if (bool_pixel_array(ii,jj,2)) then
								angle = angle_pixel_array(ii,jj,2) - clockwise
								itor = floor(angle_pixel_array(ii,jj,2)/Dtor) + 1_idef

								XC = (/cam%position(1)*COS(angle),-cam%position(1)*SIN(angle),cam%position(2)/)

								n = XC - X
								r = SQRT(DOT_PRODUCT(n,n))
								n = n/r

								dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

								beta = ACOS(DOT_PRODUCT(n,binorm))
								if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
								if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

								nperp = n - DOT_PRODUCT(n,binorm)*binorm
								nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
								chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

								theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

								P_lambda = 0.0_rp
								P_angular = 0.0_rp

								if (theta .LE. threshold_angle_simple_model) then
									call P_integral(zeta,P)

									P_lambda(:,itor) = (C_C*C_E**2)*P/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)
									np_lambda_pixel(ii,jj,itor,ss) = np_lambda_pixel(ii,jj,itor,ss) + 1.0_rp
								end if

								do ll=1_idef,cam%Nlambda ! Nlambda
									if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
										P_angular(ll,itor) = Psyn(g,psi,k,cam%lambda(ll),chi)
										if (P_angular(ll,itor).GT.0.0_rp) then
											np_angular_pixel(ii,jj,itor,ss) = np_angular_pixel(ii,jj,itor,ss) + 1.0_rp
										else
											P_angular(ll,itor) = 0.0_rp
										end if
									end if
								end do ! Nlambda	

								P_lambda(:,itor) = dSA*P_lambda(:,itor)/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))
								P_angular(:,itor) = dSA*P_angular(:,itor)

								P_l_pixel(:,itor,ss) = P_l_pixel(:,itor,ss) + P_lambda(:,itor)
								P_a_pixel(:,itor,ss) = P_a_pixel(:,itor,ss) + P_angular(:,itor)

								if (cam%photon_count) then
									P_lambda(:,itor) = P_lambda(:,itor)/photon_energy
									P_angular(:,itor) = P_angular(:,itor)/photon_energy

									Psyn_lambda_pixel(ii,jj,itor,ss) = Psyn_lambda_pixel(ii,jj,itor,ss) + SUM(P_lambda(:,itor))
									Psyn_angular_pixel(ii,jj,itor,ss) = Psyn_angular_pixel(ii,jj,itor,ss) + SUM(P_angular(:,itor))
								else
									Psyn_lambda_pixel(ii,jj,itor,ss) = Psyn_lambda_pixel(ii,jj,itor,ss) + &
																trapz(cam%lambda,P_lambda(:,itor))
									Psyn_angular_pixel(ii,jj,itor,ss) = Psyn_angular_pixel(ii,jj,itor,ss) + &
																trapz(cam%lambda,P_angular(:,itor))
								end if
							end if

						end do ! NY
					end do ! NX
				end if ! check if bool == TRUE

			end if ! if confined
		end do ! particles
!$OMP END PARALLEL DO
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn has units of Watts/m 
!	or (photons/s)(m^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	if (params%mpi_params%rank.EQ.0_idef) then
		if (.NOT.cam%toroidal_sections) then
			ALLOCATE(array3D(cam%num_pixels(1),cam%num_pixels(2),params%num_species))
			ALLOCATE(array2D(cam%Nlambda,params%num_species))
		end if
	end if

	units = 1.0_rp

	if (params%mpi_params%nmpi.GT.1_idef) then 
		numel = cam%num_pixels(1)*cam%num_pixels(2)*params%num_species*cam%ntor_sections

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_angular_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_angular_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'Psyn_angular_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,Psyn_angular_pixel,var_name)
			else
				array3D = SUM(Psyn_angular_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if				
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_angular_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_angular_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'np_angular_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,np_angular_pixel,var_name)
			else
				array3D = SUM(np_angular_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if				
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)


		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_lambda_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_lambda_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'Psyn_lambda_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
			else
				array3D = SUM(Psyn_lambda_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_lambda_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_lambda_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'np_lambda_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,np_lambda_pixel,var_name)
			else
				array3D = SUM(np_lambda_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		numel = params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = np_pixel
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_pixel = receive_buffer

	        var_name = 'np_pixel'
	        call save_snapshot_var(params,np_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		
		numel = cam%Nlambda*params%num_species*cam%ntor_sections
		
		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_a_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_a_pixel = RESHAPE(receive_buffer,(/cam%Nlambda,cam%ntor_sections,params%num_species/))

	        var_name = 'P_a_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,P_a_pixel,var_name)
			else
				array2D = SUM(P_a_pixel,2)
				call save_snapshot_var(params,array2D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_l_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_l_pixel = RESHAPE(receive_buffer,(/cam%Nlambda,cam%ntor_sections,params%num_species/))

	        var_name = 'P_l_pixel'

			if (cam%toroidal_sections) then
		        call save_snapshot_var(params,P_l_pixel,var_name)
			else
				array2D = SUM(P_l_pixel,2)
				call save_snapshot_var(params,array2D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

	    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	else
		var_name = 'np_angular_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,np_angular_pixel,var_name)
		else
			array3D = SUM(np_angular_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	

		var_name = 'np_lambda_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,np_lambda_pixel,var_name)
		else
			array3D = SUM(np_lambda_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	

		var_name = 'Psyn_angular_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,Psyn_angular_pixel,var_name)
		else
			array3D = SUM(Psyn_angular_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	

		var_name = 'Psyn_lambda_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
		else
			array3D = SUM(Psyn_lambda_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	


		var_name = 'np_pixel'
		call save_snapshot_var(params,np_pixel,var_name)

		var_name = 'P_a_pixel'
		if (cam%toroidal_sections) then
	        call save_snapshot_var(params,P_a_pixel,var_name)
		else
			array2D = SUM(P_a_pixel,2)
			call save_snapshot_var(params,array2D,var_name)
		end if	

		var_name = 'P_l_pixel'
		if (cam%toroidal_sections) then
	        call save_snapshot_var(params,P_l_pixel,var_name)
		else
			array2D = SUM(P_l_pixel,2)
			call save_snapshot_var(params,array2D,var_name)
		end if	
	end if


	DEALLOCATE(bool_pixel_array)
	DEALLOCATE(angle_pixel_array)

	DEALLOCATE(np_angular_pixel)
    DEALLOCATE(Psyn_angular_pixel)

	DEALLOCATE(np_lambda_pixel)
    DEALLOCATE(Psyn_lambda_pixel)

	DEALLOCATE(P_lambda)
	DEALLOCATE(P_angular)

	DEALLOCATE(P_l_pixel)
	DEALLOCATE(P_a_pixel)
	DEALLOCATE(np_pixel)

	DEALLOCATE(P)

	DEALLOCATE(zeta)

	DEALLOCATE(photon_energy)

	if (ALLOCATED(array3D)) DEALLOCATE(array3D)
	if (ALLOCATED(array2D)) DEALLOCATE(array2D)
END SUBROUTINE integrated_SE_toroidal_sections


SUBROUTINE integrated_SE_3D(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	CHARACTER(MAX_STRING_LENGTH) :: var_name
	LOGICAL :: bool
	REAL(rp), DIMENSION(3) :: binorm, n, nperp
	REAL(rp), DIMENSION(3) :: X, V, B, E, XC
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: np_angular_pixel
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn_angular_pixel
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: np_lambda_pixel
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn_lambda_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_lambda_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_angular_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_l_pixel
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_a_pixel
	REAL(rp), DIMENSION(:), ALLOCATABLE :: np_pixel
	REAL(rp), DIMENSION(:), ALLOCATABLE :: P
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: P_lambda, P_angular
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: array3D
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: array2D
	REAL(rp) :: q, m, k, u, g, l, threshold_angle, threshold_angle_simple_model
	REAL(rp) :: psi, chi, beta, theta
	REAL(rp), DIMENSION(:), ALLOCATABLE :: photon_energy
	REAL(rp) :: r,solid_angle
	REAL(rp) :: angle, clockwise
	REAL(rp) :: units
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer, receive_buffer
	REAL(rp) :: lc
	REAL(rp) :: dSA ! Element of solid angle
	REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
	REAL(rp) :: azimuthal_angle,Dtor
	INTEGER :: ii,jj,ll,ss,pp
	INTEGER :: itor
    INTEGER :: numel, mpierr

	ALLOCATE(np_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))
	ALLOCATE(Psyn_angular_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))

	ALLOCATE(np_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))
	ALLOCATE(Psyn_lambda_pixel(cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species))

	ALLOCATE(P(cam%Nlambda))
	
	ALLOCATE(P_lambda(cam%Nlambda,cam%ntor_sections))
	ALLOCATE(P_angular(cam%Nlambda,cam%ntor_sections))

	ALLOCATE(P_l_pixel(cam%Nlambda,cam%ntor_sections,params%num_species))
	ALLOCATE(P_a_pixel(cam%Nlambda,cam%ntor_sections,params%num_species))
	ALLOCATE(np_pixel(params%num_species))

	ALLOCATE(zeta(cam%Nlambda))

	ALLOCATE(photon_energy(cam%Nlambda))

	np_angular_pixel = 0.0_rp
	Psyn_angular_pixel = 0.0_rp

	np_lambda_pixel = 0.0_rp
	Psyn_lambda_pixel = 0.0_rp

	P_l_pixel = 0.0_rp
	P_a_pixel = 0.0_rp
	np_pixel = 0.0_rp

	Dtor = 2.0_rp*C_PI/REAL(cam%ntor_sections,rp)

	zeta = 0.0_rp

	photon_energy = C_h*C_C/cam%lambda

	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL DO FIRSTPRIVATE(q,m,Dtor,photon_energy) PRIVATE(binorm,n,nperp,X,XC,V,B,E,&
!$OMP& k,u,g,l,threshold_angle,threshold_angle_simple_model,theta,&
!$OMP& psi,chi,beta,bool,angle,clockwise,ii,jj,ll,pp,r,lc,zeta,P_lambda,P_angular,itor,solid_angle,dSA,P,azimuthal_angle)&
!$OMP& SHARED(params,spp,ss,Psyn_angular_pixel,np_angular_pixel,np_lambda_pixel,Psyn_lambda_pixel,P_l_pixel,P_a_pixel,np_pixel)
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

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

				zeta = lc/cam%lambda

				binorm = binorm/SQRT(DOT_PRODUCT(binorm,binorm))

				threshold_angle = (1.5_rp*k*cam%lambda_max/C_PI)**(1.0_rp/3.0_rp) ! In radians

				threshold_angle_simple_model = 1.0_rp/g

				np_pixel(ss) = np_pixel(ss) + 1.0_rp ! We count all the confined particles.

				call is_visible(X,V/u,MAXVAL((/threshold_angle,threshold_angle_simple_model/)),bool,ii,jj)
			
				if (bool.EQV..TRUE.) then
					azimuthal_angle = ATAN2(X(2),X(1))

					if (azimuthal_angle.LT.0.0_rp) azimuthal_angle = 2.0_rp*C_PI + azimuthal_angle

					itor = floor(azimuthal_angle/Dtor) + 1_idef

					XC = (/cam%position(1),0.0_rp,cam%position(2)/)

					n = XC - X
					r = SQRT(DOT_PRODUCT(n,n))
					n = n/r

					dSA = DOT_PRODUCT(n,XC/SQRT(DOT_PRODUCT(XC,XC)))*cam%pixel_area/r**2

					beta = ACOS(DOT_PRODUCT(n,binorm))
					if (beta.GT.0.5_rp*C_PI) psi = beta - 0.5_rp*C_PI
					if (beta.LT.0.5_rp*C_PI) psi = 0.5_rp*C_PI - beta

					nperp = n - DOT_PRODUCT(n,binorm)*binorm
					nperp = nperp/SQRT(DOT_PRODUCT(nperp,nperp))
					chi = ABS(ACOS(DOT_PRODUCT(nperp,V/u)))

					theta = ABS(ACOS(DOT_PRODUCT(n,V/u)))

					P_lambda = 0.0_rp
					P_angular = 0.0_rp

					if (theta .LE. threshold_angle_simple_model) then
						call P_integral(zeta,P)

						P_lambda(:,itor) = (C_C*C_E**2)*P/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)
						np_lambda_pixel(ii,jj,itor,ss) = np_lambda_pixel(ii,jj,itor,ss) + REAL(cam%Nlambda,rp)
					end if

					do ll=1_idef,cam%Nlambda ! Nlambda
						if ((chi.LT.chic(g,k,cam%lambda(ll))).AND.(psi.LT.psic(k,cam%lambda(ll)))) then
							P_angular(ll,itor) = Psyn(g,psi,k,cam%lambda(ll),chi)
							if (P_angular(ll,itor).GT.0.0_rp) then
								np_angular_pixel(ii,jj,itor,ss) = np_angular_pixel(ii,jj,itor,ss) + 1.0_rp
							else
								P_angular(ll,itor) = 0.0_rp
							end if
						end if
					end do ! Nlambda	

					P_lambda(:,itor) = dSA*P_lambda(:,itor)/(2.0_rp*C_PI*(1.0_rp - COS(1.0_rp/g)))
					P_angular(:,itor) = dSA*P_angular(:,itor)

					P_l_pixel(:,itor,ss) = P_l_pixel(:,itor,ss) + P_lambda(:,itor)
					P_a_pixel(:,itor,ss) = P_a_pixel(:,itor,ss) + P_angular(:,itor)

					if (cam%photon_count) then
						P_lambda(:,itor) = P_lambda(:,itor)/photon_energy
						P_angular(:,itor) = P_angular(:,itor)/photon_energy

						Psyn_lambda_pixel(ii,jj,itor,ss) = Psyn_lambda_pixel(ii,jj,itor,ss) + SUM(P_lambda(:,itor))
						Psyn_angular_pixel(ii,jj,itor,ss) = Psyn_angular_pixel(ii,jj,itor,ss) + SUM(P_angular(:,itor))
					else
						Psyn_lambda_pixel(ii,jj,itor,ss) = Psyn_lambda_pixel(ii,jj,itor,ss) + &
													trapz(cam%lambda,P_lambda(:,itor))
						Psyn_angular_pixel(ii,jj,itor,ss) = Psyn_angular_pixel(ii,jj,itor,ss) + &
													trapz(cam%lambda,P_angular(:,itor))
					end if
				end if ! check if bool == TRUE

			end if ! if confined
		end do ! particles
!$OMP END PARALLEL DO
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn has units of Watts/m 
!	or (photons/s)(m^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	units = 1.0_rp

	if (params%mpi_params%rank.EQ.0_idef) then
		if (.NOT.cam%toroidal_sections) then
			ALLOCATE(array3D(cam%num_pixels(1),cam%num_pixels(2),params%num_species))
			ALLOCATE(array2D(cam%Nlambda,params%num_species))
		end if
	end if

	if (params%mpi_params%nmpi.GT.1_idef) then 

		numel = cam%num_pixels(1)*cam%num_pixels(2)*params%num_species*cam%ntor_sections

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_angular_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_angular_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'Psyn_angular_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,Psyn_angular_pixel,var_name)
			else
				array3D = SUM(Psyn_angular_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if				
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_angular_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_angular_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'np_angular_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,np_angular_pixel,var_name)
			else
				array3D = SUM(np_angular_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if				
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)


		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_lambda_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_lambda_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'Psyn_lambda_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
			else
				array3D = SUM(Psyn_lambda_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_lambda_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_lambda_pixel = &
								RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'np_lambda_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,np_lambda_pixel,var_name)
			else
				array3D = SUM(np_lambda_pixel,3)
				call save_snapshot_var(params,array3D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		numel = params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = np_pixel
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_pixel = receive_buffer

	        var_name = 'np_pixel'
	        call save_snapshot_var(params,np_pixel,var_name)
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		
		numel = cam%Nlambda*params%num_species*cam%ntor_sections
		
		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_a_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_a_pixel = RESHAPE(receive_buffer,(/cam%Nlambda,cam%ntor_sections,params%num_species/))

	        var_name = 'P_a_pixel'

			if (cam%toroidal_sections) then
				call save_snapshot_var(params,P_a_pixel,var_name)
			else
				array2D = SUM(P_a_pixel,2)
				call save_snapshot_var(params,array2D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_l_pixel,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_l_pixel = RESHAPE(receive_buffer,(/cam%Nlambda,cam%ntor_sections,params%num_species/))

	        var_name = 'P_l_pixel'

			if (cam%toroidal_sections) then
		        call save_snapshot_var(params,P_l_pixel,var_name)
			else
				array2D = SUM(P_l_pixel,2)
				call save_snapshot_var(params,array2D,var_name)
			end if	
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

	    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	else
		var_name = 'np_angular_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,np_angular_pixel,var_name)
		else
			array3D = SUM(np_angular_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	

		var_name = 'np_lambda_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,np_lambda_pixel,var_name)
		else
			array3D = SUM(np_lambda_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	

		var_name = 'Psyn_angular_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,Psyn_angular_pixel,var_name)
		else
			array3D = SUM(Psyn_angular_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	

		var_name = 'Psyn_lambda_pixel'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,Psyn_lambda_pixel,var_name)
		else
			array3D = SUM(Psyn_lambda_pixel,3)
			call save_snapshot_var(params,array3D,var_name)
		end if	


		var_name = 'np_pixel'
		call save_snapshot_var(params,np_pixel,var_name)

		var_name = 'P_a_pixel'
		if (cam%toroidal_sections) then
	        call save_snapshot_var(params,P_a_pixel,var_name)
		else
			array2D = SUM(P_a_pixel,2)
			call save_snapshot_var(params,array2D,var_name)
		end if	

		var_name = 'P_l_pixel'
		if (cam%toroidal_sections) then
	        call save_snapshot_var(params,P_l_pixel,var_name)
		else
			array2D = SUM(P_l_pixel,2)
			call save_snapshot_var(params,array2D,var_name)
		end if	
	end if

	DEALLOCATE(np_angular_pixel)
    DEALLOCATE(Psyn_angular_pixel)

	DEALLOCATE(np_lambda_pixel)
    DEALLOCATE(Psyn_lambda_pixel)

	DEALLOCATE(P_lambda)
	DEALLOCATE(P_angular)

	DEALLOCATE(P_l_pixel)
	DEALLOCATE(P_a_pixel)
	DEALLOCATE(np_pixel)

	DEALLOCATE(P)

	DEALLOCATE(zeta)

	DEALLOCATE(photon_energy)

	if (ALLOCATED(array3D)) DEALLOCATE(array3D)
	if (ALLOCATED(array2D)) DEALLOCATE(array2D)
END SUBROUTINE integrated_SE_3D


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
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn_lambda
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PTot
	REAL(rp) :: Rpol, Zpol
	REAL(rp) :: q, m, k, u, g, lc
	REAL(rp), DIMENSION(:), ALLOCATABLE :: photon_energy
	INTEGER :: ii,jj,ll,ss,pp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: Psyn_send_buffer,Psyn_receive_buffer
	REAL(rp), DIMENSION(:), ALLOCATABLE :: np_send_buffer,np_receive_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: PTot_send_buffer,PTot_receive_buffer
    INTEGER :: numel, mpierr
	REAL(rp) :: units

	ALLOCATE(np(pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species))
	ALLOCATE(Psyn_lambda(pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species))
	ALLOCATE(PTot(pplane%grid_dims(1),pplane%grid_dims(2),params%num_species))
	ALLOCATE(P(cam%Nlambda))
	ALLOCATE(zeta(cam%Nlambda))
	ALLOCATE(photon_energy(cam%Nlambda))

	np = 0.0_rp
	Psyn_lambda = 0.0_rp
	P = 0.0_rp
	PTot = 0.0_rp
	photon_energy = C_h*C_C/cam%lambda
	
	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass

!$OMP PARALLEL DO FIRSTPRIVATE(q,m,photon_energy) PRIVATE(binorm,X,V,B,E,k,u,g,lc,ii,jj,ll,pp,zeta,P,Rpol,Zpol)&
!$OMP& SHARED(params,spp,ss,Psyn_lambda,PTot,np)
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

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3)
				zeta = lc/cam%lambda

				call P_integral(zeta,P)

				P = (C_C*C_E**2)*P/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)

				Rpol = SQRT(SUM(X(1:2)**2))
				Zpol = X(3)
				
				ii = FLOOR((Rpol - pplane%Rmin)/pplane%DR) + 1_idef
				jj = FLOOR((Zpol + ABS(pplane%Zmin))/pplane%DZ) + 1_idef

				Psyn_lambda(ii,jj,:,ss) = Psyn_lambda(ii,jj,:,ss) + P
				np(ii,jj,:,ss) = np(ii,jj,:,ss) + 1_idef
				PTot(ii,jj,ss) = PTot(ii,jj,ss) + spp(ss)%vars%Prad(pp);
			end if ! if confined
		end do ! particles
!$OMP END PARALLEL DO
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn_lambda has units of Watts/m 
!	or (photons/s)(m^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	if (params%mpi_params%nmpi.GT.1_idef) then 
		numel = pplane%grid_dims(1)*pplane%grid_dims(2)*cam%Nlambda*params%num_species

		ALLOCATE(Psyn_send_buffer(numel))
		ALLOCATE(Psyn_receive_buffer(numel))
		ALLOCATE(np_send_buffer(numel))
		ALLOCATE(np_receive_buffer(numel))

		Psyn_send_buffer = RESHAPE(Psyn_lambda,(/numel/))
		CALL MPI_REDUCE(Psyn_send_buffer,Psyn_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		np_send_buffer = RESHAPE(np,(/numel/))
		CALL MPI_REDUCE(np_send_buffer,np_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		numel = pplane%grid_dims(1)*pplane%grid_dims(2)*params%num_species

		ALLOCATE(PTot_send_buffer(numel))
		ALLOCATE(PTot_receive_buffer(numel))

		PTot_send_buffer = RESHAPE(PTot,(/numel/))
		CALL MPI_REDUCE(PTot_send_buffer,PTot_receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)

		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_lambda = RESHAPE(Psyn_receive_buffer,(/pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species/))
		    np = RESHAPE(np_receive_buffer,(/pplane%grid_dims(1),pplane%grid_dims(2),cam%Nlambda,params%num_species/))

			var_name = 'np_pplane'
		    call save_snapshot_var(params,np,var_name)

			var_name = 'Psyn_pplane'
	    	call save_snapshot_var(params,Psyn_lambda,var_name)

			PTot = RESHAPE(PTot_receive_buffer,(/pplane%grid_dims(1),pplane%grid_dims(2),params%num_species/))

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
		var_name = 'np_pplane'
		call save_snapshot_var(params,np,var_name)

		var_name = 'Psyn_pplane'
	    call save_snapshot_var(params,Psyn_lambda,var_name)

		units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
		PTot = units*PTot ! (Watts)
		var_name = 'PTot_pplane'
		call save_snapshot_var(params,PTot,var_name)
	end if

	DEALLOCATE(np)
    DEALLOCATE(Psyn_lambda)
	DEALLOCATE(PTot)
	DEALLOCATE(P)
	DEALLOCATE(zeta)
	DEALLOCATE(photon_energy)
END SUBROUTINE spectral_density


SUBROUTINE integrated_spectral_density(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	CHARACTER(MAX_STRING_LENGTH) :: var_name
	REAL(rp), DIMENSION(3) :: binorm
	REAL(rp), DIMENSION(3) :: X, V, B, E
	REAL(rp), DIMENSION(:), ALLOCATABLE :: P
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: P_lambda
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: np_lambda
	REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: np
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: Psyn_lambda
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: PTot
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: array3D
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: array2D
	REAL(rp), DIMENSION(:), ALLOCATABLE :: array1D
	REAL(rp) :: Dtor
	REAL(rp) :: phi
	REAL(rp) :: Rpol, Zpol
	REAL(rp) :: q, m, k, u, g, lc
	REAL(rp), DIMENSION(:), ALLOCATABLE :: photon_energy
	INTEGER :: ii,jj,kk,ll,ss,pp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer, receive_buffer
    INTEGER :: numel, mpierr
	REAL(rp) :: units

	ALLOCATE(Psyn_lambda(pplane%grid_dims(1),pplane%grid_dims(2),cam%ntor_sections,params%num_species))
	ALLOCATE(np(pplane%grid_dims(1),pplane%grid_dims(2),cam%ntor_sections,params%num_species))
	ALLOCATE(PTot(pplane%grid_dims(1),pplane%grid_dims(2),cam%ntor_sections,params%num_species))

	ALLOCATE(P(cam%Nlambda))

	ALLOCATE(P_lambda(cam%Nlambda,cam%ntor_sections,params%num_species))
	ALLOCATE(np_lambda(cam%ntor_sections,params%num_species))

	ALLOCATE(zeta(cam%Nlambda))
	ALLOCATE(photon_energy(cam%Nlambda))

	np = 0.0_rp
	Psyn_lambda = 0.0_rp
	P = 0.0_rp
	P_lambda = 0.0_rp
	PTot = 0.0_rp

	np_lambda = 0.0_rp
	photon_energy = C_h*C_C/cam%lambda

	Dtor = 2.0_rp*C_PI/REAL(cam%ntor_sections,rp)
	
	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass
!$OMP PARALLEL DO FIRSTPRIVATE(q,m,photon_energy,Dtor) PRIVATE(binorm,X,V,B,E,k,u,g,lc,ii,jj,kk,ll,pp,zeta,P,Rpol,Zpol,phi)&
!$OMP& SHARED(params,spp,ss,Psyn_lambda,PTot,np,P_lambda,np_lambda)
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

				lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3)
				zeta = lc/cam%lambda

				Rpol = SQRT(SUM(X(1:2)**2))
				Zpol = X(3)
				
				ii = FLOOR((Rpol - pplane%Rmin)/pplane%DR) + 1_idef
				jj = FLOOR((Zpol + ABS(pplane%Zmin))/pplane%DZ) + 1_idef

				phi = ATAN2(X(2),X(1))
				if (phi.LT.0.0_rp) phi = 2.0_rp*C_PI + phi
				kk = floor(phi/Dtor) + 1_idef

				call P_integral(zeta,P)

				P = (C_C*C_E**2)*P/(SQRT(3.0_rp)*C_E0*g**2*cam%lambda**3)

				P_lambda(:,kk,ss) = P_lambda(:,kk,ss) + P
				np_lambda(kk,ss) = np_lambda(kk,ss) + 1.0_rp

				Psyn_lambda(ii,jj,kk,ss) = Psyn_lambda(ii,jj,kk,ss) + trapz(cam%lambda,P)
				np(ii,jj,kk,ss) = np(ii,jj,kk,ss) + 1_idef
				PTot(ii,jj,kk,ss) = PTot(ii,jj,kk,ss) + spp(ss)%vars%Prad(pp);
			end if ! if confined
		end do ! particles
!$OMP END PARALLEL DO
	end do ! species

!	* * * * * * * IMPORTANT * * * * * * *
!	* * * * * * * * * * * * * * * * * * *
!	Here Psyn_lambda has units of Watts/m 
!	or (photons/s)(m^-1). See above.
!	* * * * * * * * * * * * * * * * * * *
!	* * * * * * * IMPORTANT * * * * * * *

	if (params%mpi_params%rank.EQ.0_idef) then
		if (.NOT.cam%toroidal_sections) then
			ALLOCATE(array3D(pplane%grid_dims(1),pplane%grid_dims(2),params%num_species))
			ALLOCATE(array2D(cam%Nlambda,params%num_species))
			ALLOCATE(array1D(params%num_species))
		end if
	end if

	if (params%mpi_params%nmpi.GT.1_idef) then 
		numel = pplane%grid_dims(1)*pplane%grid_dims(2)*cam%ntor_sections*params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(Psyn_lambda,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    Psyn_lambda = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'Psyn_pplane'
			if (cam%toroidal_sections) then
				call save_snapshot_var(params,Psyn_lambda,var_name)
			else
				array3D = SUM(Psyn_lambda,3)
				call save_snapshot_var(params,array3D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			var_name = 'np_pplane'
			if (cam%toroidal_sections) then
				call save_snapshot_var(params,np,var_name)
			else
				array3D = SUM(np,3)
				call save_snapshot_var(params,array3D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(PTot,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    PTot = RESHAPE(receive_buffer,(/cam%num_pixels(1),cam%num_pixels(2),cam%ntor_sections,params%num_species/))

			units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
			PTot = units*PTot ! (Watts)
			var_name = 'PTot_pplane'
			if (cam%toroidal_sections) then
				call save_snapshot_var(params,PTot,var_name)
			else
				array3D = SUM(PTot,3)
				call save_snapshot_var(params,array3D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)


		numel = cam%ntor_sections*params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(np_lambda,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    np_lambda = RESHAPE(receive_buffer,(/cam%ntor_sections,params%num_species/))

			var_name = 'np_lambda'
			if (cam%toroidal_sections) then
				call save_snapshot_var(params,np_lambda,var_name)
			else
				array1D = SUM(np_lambda,1)
				call save_snapshot_var(params,array1D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)


		numel = cam%Nlambda*cam%ntor_sections*params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(P_lambda,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    P_lambda = RESHAPE(receive_buffer,(/cam%Nlambda,cam%ntor_sections,params%num_species/))

			var_name = 'P_lambda'
			if (cam%toroidal_sections) then
				call save_snapshot_var(params,P_lambda,var_name)
			else
				array2D = SUM(P_lambda,2)
				call save_snapshot_var(params,array2D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

	    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	else
		var_name = 'np_pplane'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,np,var_name)
		else
			array3D = SUM(np,3)
			call save_snapshot_var(params,array3D,var_name)
		end if

		var_name = 'Psyn_pplane'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,Psyn_lambda,var_name)
		else
			array3D = SUM(Psyn_lambda,3)
			call save_snapshot_var(params,array3D,var_name)
		end if

		var_name = 'P_lambda'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,P_lambda,var_name)
		else
			array2D = SUM(P_lambda,2)
			call save_snapshot_var(params,array2D,var_name)
		end if

		var_name = 'np_lambda'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,np_lambda,var_name)
		else
			array1D = SUM(np_lambda,1)
			call save_snapshot_var(params,array1D,var_name)
		end if

		units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
		PTot = units*PTot ! (Watts)
		var_name = 'PTot_pplane'
		if (cam%toroidal_sections) then
			call save_snapshot_var(params,PTot,var_name)
		else
			array3D = SUM(PTot,3)
			call save_snapshot_var(params,array3D,var_name)
		end if
	end if

	DEALLOCATE(np)
    DEALLOCATE(Psyn_lambda)
	DEALLOCATE(PTot)
	DEALLOCATE(P)
	DEALLOCATE(P_lambda)
	DEALLOCATE(np_lambda)
	DEALLOCATE(zeta)
	DEALLOCATE(photon_energy)
	
	if (ALLOCATED(array3D)) DEALLOCATE(array3D)
	if (ALLOCATED(array2D)) DEALLOCATE(array2D)
	if (ALLOCATED(array1D)) DEALLOCATE(array1D)
END SUBROUTINE integrated_spectral_density


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * SUBROUTINES TO GENERATE OUTPUTS OF THE SYNTHETIC CAMERA * * * * 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


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

!		dset = TRIM(gname) // "/aperture"
!		attr = "Aperture of the camera (m)"
!		call save_to_hdf5(h5file_id,dset,cam%aperture,attr)

		dset = TRIM(gname) // "/pixel_area"
		attr = "Pixel area (m^2)"
		call save_to_hdf5(h5file_id,dset,cam%pixel_area,attr)

		dset = TRIM(gname) // "/start_at"
		attr = "Time at which camera starts working (s)"
		call save_to_hdf5(h5file_id,dset,cam%start_at,attr)

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
		call save_to_hdf5(h5file_id,dset,cam%lambda_min,attr)

		dset = TRIM(gname) // "/lambda_max"
		attr = "Minimum wavelength (m)"
		call save_to_hdf5(h5file_id,dset,cam%lambda_max,attr)

		dset = TRIM(gname) // "/Dlambda"
		attr = "Step between finite wavelengths (m)"
		call save_to_hdf5(h5file_id,dset,cam%Dlambda,attr)

		dset = TRIM(gname) // "/Nlambda"
		attr = "Number of finite wavelengths (m)"
		call save_to_hdf5(h5file_id,dset,cam%Nlambda,attr)

	    dset = TRIM(gname) // "/lambda"
	    call save_1d_array_to_hdf5(h5file_id,dset,cam%lambda)

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

		dset = TRIM(gname) // "/photon_count"
		attr = "Logical variable: 1=Psyn is in number of photons, 0=Psyn is in Watts"
		if (cam%photon_count) then
			call save_to_hdf5(h5file_id,dset,1_idef,attr)
		else
			call save_to_hdf5(h5file_id,dset,0_idef,attr)
		end if

		dset = TRIM(gname) // "/integrated_opt"
		attr = "Logical variable: 1=integrated spectra, 0=detailed spectral info"
		if (cam%integrated_opt) then
			call save_to_hdf5(h5file_id,dset,1_idef,attr)
		else
			call save_to_hdf5(h5file_id,dset,0_idef,attr)
		end if

			dset = TRIM(gname) // "/toroidal_sections"
			attr = "Logical variable: 1=decomposed in toroidal sections, 0=no toroidal decomposition"
		if (cam%toroidal_sections) then
			call save_to_hdf5(h5file_id,dset,1_idef,attr)

			dset = TRIM(gname) // "/ntor_sections"
			attr = "Number of toroidal sections"
			call save_to_hdf5(h5file_id,dset,cam%ntor_sections,attr)
		else
			call save_to_hdf5(h5file_id,dset,0_idef,attr)
		end if

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


SUBROUTINE save_snapshot_var_1d(params,var,var_name)
IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: var
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
		call save_to_hdf5(subgroup_id,dset,var(ss))

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_1d


SUBROUTINE save_snapshot_var_2d(params,var,var_name)
IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: var
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
		call save_array_to_hdf5(subgroup_id, dset, var(:,ss))

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_2d


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


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * MAIN CALL TO SYNTHETIC CAMERA SUBROUTINES * * * * 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


SUBROUTINE synthetic_camera(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp

	if (cam%camera_on.AND.(params%time*params%cpp%time >= cam%start_at)) then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("Synthetic camera diagnostic: ON!")')
		end if

		if (cam%integrated_opt) then
!			call integrated_SE_toroidal_sections(params,spp)
			call integrated_SE_3D(params,spp)
			call integrated_spectral_density(params,spp)
		else
			call angular_density(params,spp)
			call spectral_density(params,spp)
		end if

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("Synthetic camera diagnostic: OFF!")')
		end if
	end if
END SUBROUTINE synthetic_camera

END MODULE korc_synthetic_camera
