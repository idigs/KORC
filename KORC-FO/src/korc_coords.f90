module korc_coords

    use korc_types
    use korc_constants

    implicit none

    PUBLIC :: cart_to_cyl,cart_to_tor,cart_to_tor_check_if_confined

    contains

subroutine cart_to_cyl(X,Xcyl)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Xcyl ! Xcyl(1,:) = R, Xcyl(2,:) = phi, Xcyl(3,:) = Z
	INTEGER :: pp, ss ! Iterators
	
	ss = SIZE(X,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(X,Xcyl)
	do pp=1_idef,ss
		Xcyl(1,pp) = SQRT(X(1,pp)**2 + X(2,pp)**2)
		Xcyl(2,pp) = ATAN2(X(2,pp), X(1,pp))
		Xcyl(2,pp) = MODULO(Xcyl(2,pp), 2.0_rp*C_PI)
		Xcyl(3,pp) = X(3,pp)
	end do
!$OMP END PARALLEL DO
end subroutine cart_to_cyl


subroutine cart_to_tor(X,Ro,Xtor,flag)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), INTENT(IN) :: Ro
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Xtor ! Xtor(1,:) = r, Xtor(2,:) = theta, Xtor(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
	INTEGER :: pp, ss ! Iterators

	ss = SIZE(X,2)

!$OMP PARALLEL DO FIRSTPRIVATE(ss,Ro) PRIVATE(pp) SHARED(X,Xtor,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
		    Xtor(1,pp) = SQRT( (SQRT(X(1,pp)**2 + X(2,pp)**2) - Ro)**2 + X(3,pp)**2 )
		    Xtor(2,pp) = ATAN2(X(3,pp), SQRT(X(1,pp)**2 + X(2,pp)**2) - Ro)
		    Xtor(2,pp) = MODULO(Xtor(2,pp),2.0_rp*C_PI)
		    Xtor(3,pp) = ATAN2(X(1,pp),X(2,pp))
		    Xtor(3,pp) = MODULO(Xtor(3,pp),2.0_rp*C_PI)
        end if
	end do
!$OMP END PARALLEL DO
end subroutine cart_to_tor


subroutine cart_to_tor_check_if_confined(X,F,Xtor,flag)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	TYPE(FIELDS), INTENT(IN) :: F
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Xtor ! Xtor(1,:) = r, Xtor(2,:) = theta, Xtor(3,:) = zeta
	INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: flag
	REAL(rp) :: a, Ro
	INTEGER :: pp, ss ! Iterators

	ss = SIZE(X,2)
	a = F%AB%a
	Ro = F%AB%Ro

!$OMP PARALLEL DO FIRSTPRIVATE(ss,a,Ro) PRIVATE(pp) SHARED(X,Xtor,flag)
	do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
		    Xtor(1,pp) = SQRT( (SQRT(X(1,pp)**2 + X(2,pp)**2) - Ro)**2 + X(3,pp)**2 )
		    Xtor(2,pp) = ATAN2(X(3,pp), SQRT(X(1,pp)**2 + X(2,pp)**2) - Ro)
		    Xtor(2,pp) = MODULO(Xtor(2,pp),2.0_rp*C_PI)
		    Xtor(3,pp) = ATAN2(X(1,pp),X(2,pp))
		    Xtor(3,pp) = MODULO(Xtor(3,pp),2.0_rp*C_PI)

			if (Xtor(3,pp) .GT. F%AB%a) then
                flag(pp) = 0_is
            end if
        end if
	end do
!$OMP END PARALLEL DO
end subroutine cart_to_tor_check_if_confined

end module korc_coords
