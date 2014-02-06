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
    integer :: ntotalpoints, ntruepoints, n, nnbasx, m, k

    real(dp), dimension(:), allocatable :: V
    real(dp), dimension(3) ::  dr
    real(dp) :: Rnk, Rmk
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: B

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
        !write(luout,'(a,f20.9,2f16.9)') 'XX',wrk(:,n)
    enddo

    allocate( A( nnuclei, nnuclei ) )
    A = zero
    allocate( B( nnuclei ) )
    B = zero

    nnbasx = size( density )
    allocate( integrals( nnbasx ) )
    integrals = zero
    allocate( V( ntruepoints ) )
    V = zero

    ! first determine the total potential
    ! todo: separate into own function
    do m = 1, nnuclei
       do k = 1, ntruepoints
           dr = Rm(:,m) - wrk(:,k)
           Rmk = sqrt( dot( dr, dr ) )

           ! call one_electron_integrals( )
           ! nuclear part + electronic part
           V(k) = V(k) + Zm(m) / Rmk
           V(k) = V(k) + dot( density, integrals )
       enddo
    enddo

    do m = 1, nnuclei
       do k = 1, ntruepoints
           dr = Rm(:,m) - wrk(:,k)
           Rmk = sqrt( dot( dr, dr ) )

           ! take care of the RHS of the linear equations
           B(m) = B(m) + V(k) / Rmk

           do n = 1, nnuclei
               dr = Rm(:,n) - wrk(:,k)
               Rnk = sqrt( dot( dr, dr ) )
               A(m,n) = A(m,n) + one / (Rmk * Rnk)
           enddo
       enddo
    enddo

    do m = 1, nnuclei
         write(luout, *) (A(m,n),n=1,nnuclei)
    enddo

    write(luout,*)
    do m = 1, nnuclei
         write(luout,*) B(m)
    enddo

    deallocate( integrals )
    deallocate( wrk )

end subroutine fit_density

end module qfit
