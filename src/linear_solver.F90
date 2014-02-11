!------------------------------------------------------------------------------
!> @brief Module to solve a system of linear equations Ax = b
!!
!! @author Casper Steinmann
module linear_solver

    use qfit_precision

    implicit none

    private

    public :: linear_solve

contains

subroutine linear_solve( A, b, x )

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(out) :: x

    integer :: lsize, info
    real(dp), dimension(:,:), allocatable :: wrk
    integer, dimension(:), allocatable :: ipiv

    lsize = size(A,1)

    ! check that A, b and x are the same dimensions
    ! and that A is square
    if (lsize /= size(A,2)) then
        stop 'ERROR [solve]: Matrix A is not square.'
    endif

    if (lsize /= size(b)) then
        stop 'ERROR [solve]: b is not the same dimension as A'
    endif

    if (lsize /= size(x)) then
        stop 'ERROR [solve]: x is not the same dimension as A'
    endif

    allocate( wrk( lsize, lsize ) )
    wrk = A

    allocate( ipiv( lsize ) )
    x = b

    ! we solve this by making an LU factorization first
    call dgetrf( lsize, lsize, wrk, lsize, ipiv, info )
    call dgetrs( 'N', lsize, 1, wrk, lsize, ipiv, x, lsize, info )

    deallocate( ipiv )
    deallocate( wrk )

end subroutine linear_solve

end module linear_solver
