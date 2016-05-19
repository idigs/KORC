module korc_fields
    use korc_types
    implicit none

    contains

subroutine analytical_magnetic_field(EB,Y,B)
    implicit none
	TYPE(FIELDS), INTENT(IN) :: EB
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Y ! Y(1,:) = r, Y(2,:) = theta, Y(3,:) = zeta
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B ! B(1,:) = Bx, B(2,:) = By, B(3,:) = Bz
	REAL(rp) :: Bp, Br, eta
	INTEGER(ip) pp ! Iterator(s)
	INTEGER(ip) :: ss

	ss = SIZE(Y,2)

!$OMP PARALLEL FIRSTPRIVATE(ss) PRIVATE(pp,Bp,Br,eta) SHARED(EB,Y,B)
!$OMP DO
	do pp=1,ss
		Bp = EB%AB%Bpo*( Y(1,pp)/EB%AB%lambda )/( 1.0_rp + (Y(1,pp)/EB%AB%lambda)**2 )
		eta = Y(1,pp)/EB%AB%Ro
		Br = 1.0_rp/( 1.0_rp + eta*cos(Y(2,pp)) )

		B(1,pp) = Br*( EB%AB%Bo*cos(Y(3,pp)) - Bp*sin(Y(2,pp))*sin(Y(3,pp)) )
		B(2,pp) = -Br*( EB%AB%Bo*sin(Y(3,pp)) + Bp*sin(Y(2,pp))*cos(Y(3,pp)) )
		B(3,pp) = Br*Bp*cos(Y(2,pp))
	end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine analytical_magnetic_field

end module korc_fields
