MODULE korc_synthetic_camera
	USE korc_types
	USE korc_constants
	USE korc_HDF5
	IMPLICIT NONE

	TYPE, PRIVATE :: CAMERA
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

		REAL(rp) :: lambda_min ! Minimum wavelength in meters
		REAL(rp) :: lambda_max ! Maximum wavelength in meters
		INTEGER :: Nlambda
		REAL(rp) :: Dlambda
		REAL(rp), DIMENSION(:), ALLOCATABLE :: lambda
	END TYPE CAMERA

	TYPE(CAMERA), PRIVATE :: cam

	PRIVATE :: clockwise_rotation,anticlockwise_rotation,ajyik
	PUBLIC :: initialize_synthetic_camera,synthetic_camera

	CONTAINS


SUBROUTINE initialize_synthetic_camera(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp) :: Riw ! Radial position of inner wall
	INTEGER, DIMENSION(2) :: num_pixels ! Number of pixels (X,Y)
	REAL(rp), DIMENSION(2) :: sensor_size ! (horizontal,vertical)
	REAL(rp) :: focal_length
	REAL(rp), DIMENSION(2) :: position ! Position of camera (R,Z)
	REAL(rp) :: incline
	REAL(rp) :: lambda_min ! Minimum wavelength in meters
	REAL(rp) :: lambda_max ! Maximum wavelength in meters
	INTEGER :: Nlambda
	REAL(rp) :: xmin, xmax, ymin, ymax, DX, DY
	INTEGER :: ii

	NAMELIST /SyntheticCamera/ Riw,num_pixels,sensor_size,&
			focal_length,position,incline,lambda_min,lambda_max,Nlambda

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * * * * * * * * * * * * * *")')
		write(6,'("*  Initializing synthetic camera  *")')
	end if

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=SyntheticCamera)
	close(default_unit_open)

!	write(*,nml=SyntheticCamera)

	cam%Riw = Riw
	cam%num_pixels = num_pixels
	cam%sensor_size = sensor_size
	cam%focal_length = focal_length
	cam%position = position
	cam%incline = incline
	cam%horizontal_angle_view = ATAN2(0.5_rp*cam%sensor_size(1),cam%focal_length)
	cam%vertical_angle_view = ATAN2(0.5_rp*cam%sensor_size(2),cam%focal_length)

	cam%lambda_min = lambda_min
	cam%lambda_max = lambda_max
	cam%Nlambda = Nlambda
	cam%Dlambda = (lambda_max - lambda_min)/REAL(cam%Nlambda,rp)
	ALLOCATE(cam%lambda(cam%Nlambda))
	
	do ii=1_idef,cam%Nlambda
		cam%lambda(ii) = lambda_min + REAL(ii-1_idef,rp)*cam%Dlambda
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

	ymin = cam%position(2) - 0.5_rp*cam%sensor_size(2)
	ymax = cam%position(2) + 0.5_rp*cam%sensor_size(2)
	DY = cam%sensor_size(2)/REAL(cam%num_pixels(2),rp)

	do ii=1_idef,cam%num_pixels(2)
		cam%pixels_nodes_y(ii) = ymin + 0.5_rp*DY + REAL(ii-1_idef,rp)*DY
	end do

	do ii=1_idef,cam%num_pixels(2)+1_idef
		cam%pixels_edges_y(ii) = ymin + REAL(ii-1_idef,rp)*DY
	end do
	
	if (params%mpi_params%rank .EQ. 0) then
		write(6,'("*     Synthetic camera ready!     *")')
		write(6,'("* * * * * * * * * * * * * * * * * *",/)')
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

! * * * * * * * * * * * * * * * !
! * * * * * FUNCTIONS * * * * * !
! * * * * * * * * * * * * * * * !


SUBROUTINE isVisible(X,V,threshold_angle,bool,angle)
	IMPLICIT NONE
	REAL(rp), DIMENSION(3), INTENT(IN) :: X
	REAL(rp), DIMENSION(3), INTENT(IN) :: V
	REAL(rp), INTENT(IN) :: threshold_angle
	LOGICAL, INTENT(OUT) :: bool
	REAL(rp), INTENT(OUT) :: angle
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
		
		angle = ATAN2(X(2) + s*V(2),X(1) + s*V(1))
		if (angle.LT.0.0_rp) then
			angle = angle + 2.0_rp*C_PI
		end if
	
		vec(1) = cam%position(1)*COS(angle) - X(1)
		vec(2) = cam%position(1)*SIN(angle) - X(2)
		vec(3) = cam%position(2) - X(3)

		vec = vec/SQRT(DOT_PRODUCT(vec,vec))
		
		psi = ACOS(DOT_PRODUCT(vec,V))

		if (psi.LE.threshold_angle) then
			bool = .TRUE. ! The particle is visible
			write(6,'("Visible")')
		else
			bool = .FALSE. ! The particle is not visible
		end if
	end if

END SUBROUTINE isVisible


SUBROUTINE synthetic_camera(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	REAL(rp), DIMENSION(3) :: binormal
	REAL(rp), DIMENSION(3) :: X,V
	REAL(rp) :: q, m, k, u, threshold_angle
	LOGICAL :: bool
	REAL(rp) :: angle	
	INTEGER :: ii,pp

	do ii=1_idef,params%num_species
		q = ABS(spp(ii)%q)
		m = spp(ii)%m

		do pp=1_idef,spp(ii)%ppp
			if ( spp(ii)%vars%flag(pp) .EQ. 1_idef ) then
				V = spp(ii)%vars%V(:,pp)
				X = spp(ii)%vars%X(:,pp)

				binormal = cross(V,spp(ii)%vars%E(:,pp)) +&
				cross(V,cross(V,spp(ii)%vars%B(:,pp)))
		
				u = SQRT(DOT_PRODUCT(V,V))
				k = q*SQRT(DOT_PRODUCT(binormal,binormal))/(spp(ii)%vars%gamma(pp)*m*u**3)
				k = k/params%cpp%length ! Now with units (m)

				binormal = binormal/SQRT(DOT_PRODUCT(binormal,binormal))

				threshold_angle = (1.5_rp*k*cam%lambda_max/C_PI)**(1.0_rp/3.0_rp) ! In radians

				call isVisible(X*params%cpp%length,V/u,threshold_angle,bool,angle)

			end if
		end do

	end do
END SUBROUTINE synthetic_camera


SUBROUTINE ajyik( x, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )
!*****************************************************************************80
!
!! AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
!
!  Discussion: 
!
!    Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions 
!    Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.  X should not be zero.
!
!    Output, real ( kind = 8 ) VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2,
!    the values of J1/3(x), J2/3(x), Y1/3(x), Y2/3(x), I1/3(x), I2/3(x),
!    K1/3(x), K2/3(x).
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) b0
  real ( kind = 8 ) c0
  real ( kind = 8 ) ck
  real ( kind = 8 ) gn
  real ( kind = 8 ) gn1
  real ( kind = 8 ) gn2
  real ( kind = 8 ) gp1
  real ( kind = 8 ) gp2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  integer ( kind = 4 ) l
  real ( kind = 8 ) pi
  real ( kind = 8 ) pv1
  real ( kind = 8 ) pv2
  real ( kind = 8 ) px
  real ( kind = 8 ) qx
  real ( kind = 8 ) r
  real ( kind = 8 ) rp
  real ( kind = 8 ) rp2
  real ( kind = 8 ) rq
  real ( kind = 8 ) sk
  real ( kind = 8 ) sum
  real ( kind = 8 ) uj1
  real ( kind = 8 ) uj2
  real ( kind = 8 ) uu0
  real ( kind = 8 ) vi1
  real ( kind = 8 ) vi2
  real ( kind = 8 ) vil
  real ( kind = 8 ) vj1
  real ( kind = 8 ) vj2
  real ( kind = 8 ) vjl
  real ( kind = 8 ) vk1
  real ( kind = 8 ) vk2
  real ( kind = 8 ) vl
  real ( kind = 8 ) vsl
  real ( kind = 8 ) vv
  real ( kind = 8 ) vv0
  real ( kind = 8 ) vy1
  real ( kind = 8 ) vy2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) xk

  if ( x == 0.0D+00 ) then
    vj1 = 0.0D+00
    vj2 = 0.0D+00
    vy1 = -1.0D+300
    vy2 = 1.0D+300
    vi1 = 0.0D+00
    vi2 = 0.0D+00
    vk1 = -1.0D+300
    vk2 = -1.0D+300
    return
  end if

  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  gp1 = 0.892979511569249D+00
  gp2 = 0.902745292950934D+00
  gn1 = 1.3541179394264D+00
  gn2 = 2.678938534707747D+00
  vv0 = 0.444444444444444D+00
  uu0 = 1.1547005383793D+00
  x2 = x * x

  if ( x < 35.0D+00 ) then
    k0 = 12
  else if ( x < 50.0D+00 ) then
    k0 = 10
  else
    k0 = 8
  end if

  if ( x <= 12.0D+00 ) then

    do l = 1, 2
      vl = l / 3.0D+00
      vjl = 1.0D+00
      r = 1.0D+00
      do k = 1, 40
        r = -0.25D+00 * r * x2 / ( k * ( k + vl ) )
        vjl = vjl + r
        if ( abs ( r ) < 1.0D-15 ) then
          exit
        end if
      end do

      a0 = ( 0.5D+00 * x ) ** vl
      if ( l == 1 ) then
        vj1 = a0 / gp1 * vjl
      else
        vj2 = a0 / gp2 * vjl
      end if

    end do

  else

    do l = 1, 2

      vv = vv0 * l * l
      px = 1.0D+00
      rp = 1.0D+00

      do k = 1, k0
        rp = - 0.78125D-02 * rp &
          * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
          * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
          / ( k * ( 2.0D+00 * k - 1.0D+00 ) * x2 )
        px = px + rp
      end do

      qx = 1.0D+00
      rq = 1.0D+00
      do k = 1, k0
        rq = - 0.78125D-02 * rq &
          * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
          * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
          / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
        qx = qx + rq
      end do

      qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
      xk = x - ( 0.5D+00 * l / 3.0D+00 + 0.25D+00 ) * pi
      a0 = sqrt ( rp2 / x )
      ck = cos ( xk )
      sk = sin ( xk )
      if ( l == 1) then
        vj1 = a0 * ( px * ck - qx * sk )
        vy1 = a0 * ( px * sk + qx * ck )
      else
        vj2 = a0 * ( px * ck - qx * sk )
        vy2 = a0 * ( px * sk + qx * ck )
      end if

    end do

  end if

  if ( x <= 12.0D+00 ) then

    do l = 1, 2

      vl = l / 3.0D+00
      vjl = 1.0D+00
      r = 1.0D+00
      do k = 1, 40
        r = -0.25D+00 * r * x2 / ( k * ( k - vl ) )
        vjl = vjl + r
        if ( abs ( r ) < 1.0D-15 ) then
          exit
        end if
      end do

      b0 = ( 2.0D+00 / x ) ** vl
      if ( l == 1 ) then
        uj1 = b0 * vjl / gn1
      else
         uj2 = b0 * vjl / gn2
      end if

    end do

    pv1 = pi / 3.0D+00
    pv2 = pi / 1.5D+00
    vy1 = uu0 * ( vj1 * cos ( pv1 ) - uj1 )
    vy2 = uu0 * ( vj2 * cos ( pv2 ) - uj2 )

  end if

  if ( x <= 18.0D+00 ) then

    do l = 1, 2
      vl = l / 3.0D+00
      vil = 1.0D+00
      r = 1.0D+00
      do k = 1, 40
        r = 0.25D+00 * r * x2 / ( k * ( k + vl ) )
        vil = vil + r
        if ( abs ( r ) < 1.0D-15 ) then
          exit
        end if
      end do

      a0 = ( 0.5D+00 * x ) ** vl

      if ( l == 1 ) then
        vi1 = a0 / gp1 * vil
      else
        vi2 = a0 / gp2 * vil
      end if

    end do

  else

    c0 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )

    do l = 1, 2
      vv = vv0 * l * l
      vsl = 1.0D+00
      r = 1.0D+00
      do k = 1, k0
        r = - 0.125D+00 * r &
          * ( vv - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
        vsl = vsl + r
      end do
      if ( l == 1 ) then
        vi1 = c0 * vsl
      else
        vi2 = c0 * vsl
      end if
    end do

  end if

  if ( x <= 9.0D+00 ) then

    do l = 1, 2
      vl = l / 3.0D+00
      if ( l == 1 ) then
        gn = gn1
      else
        gn = gn2
      end if
      a0 = ( 2.0D+00 / x ) ** vl / gn
      sum = 1.0D+00
      r = 1.0D+00
      do k = 1, 60
        r = 0.25D+00 * r * x2 / ( k * ( k - vl ) )
        sum = sum + r
        if ( abs ( r ) < 1.0D-15 ) then
          exit
        end if
      end do

      if ( l == 1 ) then
        vk1 = 0.5D+00 * uu0 * pi * ( sum * a0 - vi1 )
      else
        vk2 = 0.5D+00 * uu0 * pi * ( sum * a0 - vi2 )
      end if

    end do

  else

    c0 = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )

    do l = 1, 2
      vv = vv0 * l * l
      sum = 1.0D+00
      r = 1.0D+00
      do k = 1, k0
        r = 0.125D+00 * r * ( vv - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
        sum = sum + r
      end do
      if ( l == 1 ) then
        vk1 = c0 * sum
      else
        vk2 = c0 * sum
      end if
    end do

  end if

  return
END SUBROUTINE ajyik



END MODULE korc_synthetic_camera
