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
    public :: qfit_initialize
    public :: qfit_finalize

    contains

subroutine qfit_initialize(R, Z)

    use connolly

    real(dp), dimension(:,:), intent(in) :: R
    real(dp), dimension(:), intent(in) :: Z

    call connolly_initialize( R, Z )

end subroutine

subroutine qfit_finalize

    use connolly

    call connolly_finalize

end subroutine

!------------------------------------------------------------------------------
!> @brief fit the density to a number of points
!!
!! @details this subroutine also takes into account the nuclear charges
!!
!! @author Casper Steinmann
subroutine fit_density(density)

    use connolly
    use qfit_integrals
    use qfit_utilities

    real(dp), dimension(:), intent(in) :: density

    real(dp), dimension(:,:), allocatable :: wrk
    real(dp), dimension(:), allocatable :: integrals
    integer :: ntotalpoints, ntruepoints, n, nnbasx, m

    real(dp), dimension(:), allocatable :: V
    real(dp), dimension(3) ::  dr
    real(dp) :: Rnm

    ! populates the max_layer_points array for memory management
    ! this requires that options have been set.
    call connolly_grid_count

    !write(luout,'(a)') 'max_layer_points:'
    !do n = 1, size( max_layer_points )
    !    write(luout,'(2i6)'), n, max_layer_points(n)
    !enddo
    !write(luout,'(a6,i6)') 'sum:', sum(max_layer_points)
    !write(luout,'(a,1i6,a,i6)') "layer_points:", max_layer_points, ', sum:',sum(max_layer_points)

    ntotalpoints = sum(max_layer_points)
    allocate( wrk( 3, ntotalpoints ) )

    call connolly_grid( wrk, ntruepoints )

    write(luout, '(i4)') ntruepoints
    write(luout,*)
    do n = 1, ntruepoints
        write(luout,'(a,f20.9,2f16.9)') 'XX',wrk(:,n)
    enddo

    nnbasx = size( density )
    allocate( integrals( nnbasx ) )
    allocate( V( ntruepoints ) )
    V = 0.0_dp

    do n = 1, ntruepoints

        do m = 1, nnuclei
            dr = wrk(:,n) - Rm(:,m)
            Rnm = sqrt(dot( dr, dr ))
            call one_electron_integrals(Zm, Rm, integrals )
            V(n) = dot( density, integrals ) / Rnm
        enddo
    enddo

    deallocate( integrals )
    deallocate( wrk )

end subroutine fit_density

end module qfit
