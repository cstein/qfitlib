!------------------------------------------------------------------------------
!> @brief Module to solve a system of linear equations Ax = b
!!
!! @author Casper Steinmann
module linear_solver

    use qfit_precision

    implicit none

    private

    public :: linear_solve_svd

contains

!> @brief Solve system of linear equations Ax = b using an SVD approach
!!
!! @author Casper Steinmann
!! @see http://glowingpython.blogspot.com/2011/06/svd-decomposition-with-numpy.html
!! @param[in] A left hand side of linear equations to solve
!! @param[in] b right hand side of linear equations to solve
!! @param[out] x solution to the system of linear equations
subroutine linear_solve_svd( A, b, x )

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(out) :: x

    integer :: lsize, nsize, info, lwrk, i
    real(dp), dimension(:,:), allocatable :: U ! U
    real(dp), dimension(:,:), allocatable :: VT ! VT
    real(dp), dimension(:), allocatable :: S ! sigma
    real(dp), dimension(:), allocatable :: wrk ! temporary work array
    integer, dimension(:), allocatable :: ipiv

    ! temporary storage while solving the system of linear equations
    real(dp), dimension(:), allocatable :: c
    real(dp), dimension(:), allocatable :: W ! reduced 
    real(dp), dimension(:,:), allocatable :: Sm
    real(dp) :: EPS = 1.0d-5

    lsize = size(A,1)

    ! set all resulting values to zero
    x = 0.0_dp

    ! check that A, b and x are the same dimensions
    ! and that A is square
    if (lsize /= size(A,2)) then
        stop 'ERROR [solve_svd]: Matrix A is not square.'
    endif

    if (lsize /= size(b)) then
        stop 'ERROR [solve_svd]: b is not the same dimension as A'
    endif

    if (lsize /= size(x)) then
        stop 'ERROR [solve_svd]: x is not the same dimension as A'
    endif

    ! get ready to allocate work space
    lwrk = -1
    allocate(  U(lsize, lsize) )
    allocate( VT(lsize, lsize) )
    allocate( S(lsize) )
    allocate( wrk(1) )

    ! obtain correct size of work array
    call dgesvd( 'All', 'All', lsize, lsize, A, lsize, &
   &            S, U, lsize, VT, lsize, wrk, lwrk, info )

    if (info /= 0) then
        stop 'ERROR [solve_svd]: SVD memory allocation failed'
    endif

    ! now allocate the true work memory for the SVD
    lwrk = int(wrk(1))
    deallocate( wrk )
    allocate( wrk( lwrk ) )

    call dgesvd( 'All', 'All', lsize, lsize, A, lsize, &
   &            S, U, lsize, VT, lsize, wrk, lwrk, info )
    
    if (info /= 0) then
        stop 'ERROR [solve_svd]: SVD solutions not found. Aborting.'
    endif

    ! calculate modified (in full size) rhs as
    ! c = U.T * b
    allocate( c(lsize) )
    c = 0.0_dp
    call dgemv('T',lsize,lsize,1.0_dp,U,lsize,b,1,0.0_dp,c,1)

    ! reduce the solution space by removing singular values
    nsize = lsize
    do i=1,lsize
        if (abs(S(i)) .lt. EPS ) then
            nsize = nsize -1
        endif
    enddo

    allocate( ipiv( nsize ) )
    allocate( W( nsize ) )
    allocate( Sm( nsize, nsize ) )

    Sm = 0.0_dp
    W = c(1:nsize)
    do i=1,nsize
        Sm(i,i) = S(i)
    enddo

    ! solve the system of linear equations S(1:nsize) x = W
    call dgetrf( nsize, nsize, Sm, nsize, ipiv, info )
    call dgetrs( 'N', nsize, 1, Sm, nsize, ipiv, W, nsize, info )

    ! obtain the solution of the system of linear equations
    call dgemv('T',nsize,nsize,1.0_dp,VT(1:nsize,1:nsize),nsize,W,1,0.0_dp,x(1:nsize),1)

    deallocate( Sm )
    deallocate( W )
    deallocate( ipiv )
    deallocate( c )
    deallocate( wrk )
    deallocate( S )
    deallocate( VT )
    deallocate( U )

end subroutine

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
