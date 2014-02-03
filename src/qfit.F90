!------------------------------------------------------------------------------
!> @brief Charge fitting library
!!
!! @author Casper Steinmann
module qfit

    use qfit_precision
    use qfit_variables
    implicit none

    private

    public :: fit_density

    contains

subroutine qfit_initialize(R, Z)

    real(dp), dimension(:,:), intent(in) :: R
    integer, dimension(:), intent(in) :: Z

    integer, 

end subroutine

!------------------------------------------------------------------------------
!> @brief fit the density to a number of points
!!
!! @details this subroutine also takes into account the nuclear charges
!!
!! @author Casper Steinmann
subroutine fit_density(density, points)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:,:), intent(in) :: points

    

end subroutine fit_density

end module qfit
