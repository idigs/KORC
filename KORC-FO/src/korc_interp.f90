module korc_interp

use korc_types
use korc_coords
use korc_fields
use rnd_numbers

use EZspline_obj
use EZspline

implicit none

PUBLIC :: interp_field, interp_analytical_field, unitVectors

contains

subroutine unitVectors(Xo,EB,par,perp)
implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Xo
	TYPE(FIELDS), INTENT(IN) :: EB
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: par
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: perp
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: B
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rnd_num
	REAL(rp), PARAMETER :: tol = korc_zero
	REAL(rp) :: c1, c2, c3
	REAL(rp) :: ax, ay, az
	INTEGER :: ppp
	INTEGER :: ii


	ppp = SIZE(Xo,2) ! Number of particles

	ALLOCATE( X(3,ppp) )
	ALLOCATE( B(3,ppp) )
	ALLOCATE( rnd_num(ppp) )
	
	call cart_to_tor(Xo, EB%AB%Ro, X) ! To toroidal coords

	call init_random_seed()

	call analytical_magnetic_field(EB,X,B)

	call RANDOM_NUMBER(rnd_num)
	
	do ii=1,ppp
		par(:,ii) = B(:,ii)/sqrt( sum(B(:,ii)**2) )
		if ( ALL( ABS(par(:,ii)) .GE. korc_zero ) ) then
			c1 = par(1,ii)**2 + par(2,ii)**2
			
			az = sqrt(c1)*rnd_num(ii)
		
			c2 = 2.0_rp*par(1,ii)*par(3,ii)*az
			c3 = (par(2,ii)**2 + par(3,ii)**2)*az**2 -par(2,ii)**2
		
			ax = ( -c2 + sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
			ay = sqrt(1.0_rp - ax**2 - az**2)
		
			perp(:,ii) = (/ax, ay, az/)
		
			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				ax = ( -c2 - sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
				ay = sqrt(1.0_rp - ax**2 - az**2)
				
				perp(:,ii) = (/ax, ay, az/)
			end if
		
			if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
				az = -az
				
				c2 = 2.0_rp*par(1,ii)*par(3,ii)*az
				c3 = (par(2,ii)**2 + par(3,ii)**2)*az**2 -par(2,ii)**2
				
				ax = ( -c2 + sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
				ay = sqrt(1.0_rp - ax**2 - az**2)
				
				perp(:,ii) = (/ax, ay, az/)
				
				if ( ABS(DOT_PRODUCT(par(:,ii),perp(:,ii)) ) .GE. tol ) then
					ax = ( -c2 - sqrt(c2**2 - 4.0_rp*c1*c3) )/(2.0_rp*c1)
					ay = sqrt(1.0_rp - ax**2 - az**2)
				    
				    perp(:,ii) = (/ax, ay, az/)
				end if
			end if
		else if ( par(1,ii) .LE. korc_zero ) then
			write(6,'("Exception!")')
		else if ( par(2,ii) .LE. korc_zero ) then
			write(6,'("Exception!")')
		else if ( par(3,ii) .LE. korc_zero ) then
			write(6,'("Exception!")')
		end if
	end do


	DEALLOCATE(X)
	DEALLOCATE(B)
	DEALLOCATE(rnd_num)
end subroutine unitVectors


subroutine interp_field(prtcls,EB)
implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER :: ii, pp, ss ! Iterators

end subroutine interp_field


subroutine interp_analytical_field(prtcls,EB)
implicit none
	TYPE(PARTICLES), INTENT(INOUT) :: prtcls
	TYPE(FIELDS), INTENT(IN) :: EB

	call cart_to_tor(prtcls%X, EB%AB%Ro, prtcls%Y)

	call analytical_magnetic_field(EB,prtcls%Y,prtcls%B)
end subroutine interp_analytical_field

end module korc_interp
