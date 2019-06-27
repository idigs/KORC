!> @brief Module containing subroutines to calculate the position of the simulated particles in toroidal and cylindrical coordinates.
module korc_coords
    use korc_types
    use korc_constants

    IMPLICIT NONE

    PUBLIC :: cart_to_cyl,&
				cart_to_tor_check_if_confined

    CONTAINS

!> @brief Subroutine that converts the position of simulated particles from Cartesian @f$(x,y,z)@f$ to cylindrical @f$(R,\phi,Z)@f$ coordinates.
!! @details Here, the coordinate transformation is:
!!
!! @f$R = \sqrt{x^2 + y^2}@f$,\n
!! @f$\phi = \arctan{\left( \frac{y}{x} \right)}@f$,\n
!! @f$Z = z@f$.
!!
!! @param[in] X Particles' position in Cartesian coordinates. X(1,:) = @f$x@f$, X(2,:) = @f$y@f$, X(3,:) = @f$z@f$
!! @param[in,out] Xcyl Particles' position in cylindrical coordinates. Xcyl(1,:) = @f$R@f$, Xcyl(2,:) = @f$\phi@f$, Xcyl(3,:) = @f$Z@f$
!! @param pp Iterator.
!! @param ss Iterator.
subroutine cart_to_cyl(X,Xcyl)
    implicit none
    REAL(rp), DIMENSION(:,:), INTENT(IN)    :: X
    REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: Xcyl
    INTEGER                                 :: pp
    INTEGER                                 :: ss

    ss = SIZE(X,2)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp)
    do pp=1_idef,ss
        Xcyl(1,pp) = SQRT(X(1,pp)**2 + X(2,pp)**2)
        Xcyl(2,pp) = ATAN2(X(2,pp), X(1,pp))
        Xcyl(2,pp) = MODULO(Xcyl(2,pp), 2.0_rp*C_PI)
        Xcyl(3,pp) = X(3,pp)
    end do
!$OMP END PARALLEL DO
end subroutine cart_to_cyl


!> @brief Subroutine that converts the position of simulated particles from Cartesian @f$(x,y,z)@f$ to toroidal @f$(r,\theta, \zeta)@f$ coordinates.
!! @details In addition to performing the coordinate transformation, this subroutine checks whether a given particle is within the plasma or not.
!! A particle is not longer considered to be within the plasma if its minor radius @f$r > r_{edge}@f$, where @f$r_{edge}@f$ is the radial distance to the
!! plasma edge as measured from the magnetic axis. For more details see the analytical model of the magnetic field in korc_types.f90 and korc_fields.f90.
!!
!! The coordinate transformation is given by:
!!
!! @f$r = \sqrt{ \left[\sqrt{x^2 + y^2}-R_0\right]^2 + z^2 }@f$,\n
!! @f$\theta = \arctan{\left( \frac{z}{\sqrt{x^2 + y^2}-Ro} \right)}@f$,\n
!! @f$\zeta = \arctan{\left( \frac{x}{y} \right)}@f$,\n
!!
!! where @f$R_0@f$ is the radial position of the magnetic axis.
!!
!! @param[in] X Particles' position in Cartesian coordinates. X(1,:) = @f$x@f$, X(2,:) = @f$y@f$, X(3,:) = @f$z@f$
!! @param[in] F An instance of the KORC derived type FIELDS.
!! @param[in,out] Xtor Particles' position in cylindrical coordinates. Xtor(1,:) = @f$r@f$, Xtor(2,:) = @f$\theta@f$, Xtor(3,:) = @f$\zeta@f$
!! @param Ro Radial position of the magnetic axis.
!! @param a Distance to plasma edge as measured from the magnetic axis.
!! @param pp Iterator.
!! @param ss Iterator.
subroutine cart_to_tor_check_if_confined(X,F,Xtor,flag)
    REAL(rp), DIMENSION(:,:), INTENT(IN)     :: X
    TYPE(FIELDS), INTENT(IN)                 :: F
    REAL(rp), DIMENSION(:,:), INTENT(INOUT)  :: Xtor
    INTEGER(is), DIMENSION(:), INTENT(INOUT) :: flag
    REAL(rp)                                 :: a
    REAL(rp)                                 :: Ro
    INTEGER                                  :: pp
    INTEGER                                  :: ss

    ss = SIZE(X,2)
    a = F%AB%a
    Ro = F%AB%Ro

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp)
    do pp=1_idef,ss
        if ( flag(pp) .EQ. 1_is ) then
            Xtor(1,pp) = SQRT( (SQRT(X(1,pp)**2 + X(2,pp)**2) - Ro)**2 + X(3,pp)**2 )
            Xtor(2,pp) = ATAN2(X(3,pp), SQRT(X(1,pp)**2 + X(2,pp)**2) - Ro)
            Xtor(2,pp) = MODULO(Xtor(2,pp),2.0_rp*C_PI)
            Xtor(3,pp) = ATAN2(X(1,pp),X(2,pp))
            Xtor(3,pp) = MODULO(Xtor(3,pp),2.0_rp*C_PI)

            if (Xtor(1,pp) .GT. F%AB%a) then
                flag(pp) = 0_is
            end if
        end if
    end do
!$OMP END PARALLEL DO
end subroutine cart_to_tor_check_if_confined

end module korc_coords
