module korc_coords

    use korc_types
    use constants

    implicit none

    PUBLIC :: cart_to_cyl, cart_to_tor

    contains

subroutine cart_to_cyl(X,Xcyl)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Xcyl ! Xcyl(1,:) = R, Xcyl(2,:) = phi, Xcyl(3,:) = Z
	INTEGER :: pp, ss ! Iterators
	
	ss = SIZE(X,2)

!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(X,Xcyl)
!$OMP DO
	do pp=1,ss
		Xcyl(1,pp) = sqrt(X(1,pp)**2 + X(2,pp)**2)
		Xcyl(2,pp) = atan2(X(2,pp), X(1,pp))
		Xcyl(2,pp) = modulo(Xcyl(2,pp), 2.0_rp*C_PI)
		Xcyl(3,pp) = X(3,pp)
	end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine cart_to_cyl


subroutine cart_to_tor(X,Ro,Xtor)
    implicit none
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: X ! X(1,:) = x, X(2,:) = y, X(3,:) = z
	REAL(rp), INTENT(IN) :: Ro
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Xtor ! Xtor(1,:) = r, Xtor(2,:) = theta, Xtor(3,:) = zeta
	INTEGER :: pp, ss ! Iterators

	ss = SIZE(X,2)

!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(X,Xtor)
!$OMP DO
	do pp=1,ss
		Xtor(1,pp) = sqrt( (sqrt(X(1,pp)**2 + X(2,pp)**2) - Ro)**2 + X(3,pp)**2 )
		Xtor(2,pp) = atan2(X(3,pp), sqrt(X(1,pp)**2 + X(2,pp)**2) - Ro)
		Xtor(2,pp) = modulo(Xtor(2,pp),2.0_rp*C_PI)
		Xtor(3,pp) = atan2(X(1,pp),X(2,pp))
		Xtor(3,pp) = modulo(Xtor(3,pp),2.0_rp*C_PI)
	end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine cart_to_tor

end module korc_coords
