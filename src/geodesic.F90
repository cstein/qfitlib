!------------------------------------------------------------------------------
!> @brief Module to generate a geodesic grid
!!
!! @author Casper Steinmann
module geodesic

    use qfit_variables
    use qfit_utilities, only : dot

    implicit none

    private

    ! the maximum number of points
    integer, save, allocatable, dimension(:) :: point_include

contains

!------------------------------------------------------------------------------
!> @brief Initialize the geodesic grid generation module by specifying
!!        the tesselation level, nuclear coordinates and charges. Allocates
!!        memory. This memory needs to be released by a call to geodesic_finalize
!!
!! @author Casper Steinmann
!!
!! @param[in] R array of nuclear coordinates
!! @param[in] Z array of nuclear charges
subroutine geodesic_initialize( R, Z )

    !integer, intent(in) :: n
    real(dp), dimension(:,:), intent(in) :: R
    real(dp), dimension(:), intent(in) :: Z

    nnuclei = size(Z)

    allocate( Zm(nnuclei) )
    Zm = Z

    allocate( Rm(3,nnuclei) )
    Rm = R

    allocate( max_layer_points( qfit_nshell ) )
    max_layer_points = 0

    allocate( n_layer_points( qfit_nshell ) )
    n_layer_points = 0

end subroutine geodesic_initialize

!------------------------------------------------------------------------------
!> @brief Releases the memory acquired by the geodesic grid module.
!!
!! @author Casper Steinmann
subroutine geodesic_finalize

    deallocate( max_layer_points )
    deallocate( n_layer_points )
    if (allocated( point_include) ) deallocate( point_include )
    deallocate( Rm )
    deallocate( Zm )

end subroutine geodesic_finalize

!------------------------------------------------------------------------------
!> @brief Generates a geodesic grid stored in the coordinates storage. Only the
!!        first ntruepoints are usable
!!
!! @author Casper Steinmann
!! @param coordinates storage for coordinates of a geodesic grid. This must be
!!                    large enough to hold all points from all spheres
!! @param ntruepoints the coordinates 1 through ntruepoints are the true
!!                    coordinates to be used for the grid
subroutine geodesic_grid( coordinates, ntruepoints )

    real(dp), dimension(:,:), intent(out) :: coordinates
    integer, intent(out) :: ntruepoints

    integer :: ilayer
    real(dp) :: rscal

    do ilayer = 1, qfit_nshell
        rscal = qfit_vdwscale + (ilayer-1)*qfit_vdwincrement
        call geodesic_grid_layer( rscal )
    enddo

end subroutine geodesic_grid

!------------------------------------------------------------------------------
!> @brief Generates the grids for a single layer
!!
!! @author Casper Steinmann
!!
!! @param[in] rscal scaling factor for the grid for the current layer
!! @param[in] coordinates coordinates for @b all points on @b all spheres for
!!                        the current layer 
!! @param[inout] pinc list of indices whether or not a point in coordinates
!!                    should be included or not
!! @param[out] ntruepoints the number of non overlapping points.
subroutine geodesic_grid_layer( rscal ) !, coordinates, pinc, ntruepoints )

    real(dp), intent(in) :: rscal

end subroutine geodesic_grid_layer

end module geodesic
