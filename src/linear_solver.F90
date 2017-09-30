!------------------------------------------------------------------------------
!> @brief Module to solve a system of linear equations Ax = b
!!
!! @author Casper Steinmann
module linear_solver

    use qfit_precision
    use qfit_variables, only : luout, qfit_debug, qfit_verbose, qfit_eps

    implicit none

    private

    public :: linear_solve_svd
    public :: linear_solve_simple

contains

subroutine linear_solve_simple(A, b, x)

    external :: dgetrf, dgetrs

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(out) :: x

    integer :: na, info, i
    real(dp), dimension(:,:), allocatable :: wrk

    integer, dimension(:), allocatable :: o_piv

    x = 0.0_dp
    na = size(A,dim=1)
    if (na.ne.size(A,dim=2)) then
        stop 'ERROR [solve]: Matrix A is not square.'
    endif
    allocate(wrk(na,na))
    wrk = A
    do i = 1, na
        wrk(i,i) = wrk(i,i) + 0.0000001_dp
    enddo

    allocate(o_piv(na))
    call dgetrf(na, na, wrk, na, o_piv, info)
    if (info<0) then
        print *, "call to dgetrf failed. parameter",-info," wrong."
    elseif (info>0) then
        print *, "call to dgetrf failed: ",info
    endif
    x(:na) = b(:na)
    call dgetrs('N', na, 1, wrk, na, o_piv, x,na,info)
    if (info<0) then
        print *, "call to dgetrs failed. parameter",-info," wrong."
    elseif (info>0) then
        print *, "call to dgetrs failed: ",info
    endif
    deallocate(o_piv)
    deallocate(wrk)

end subroutine linear_solve_simple


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
    real(dp), dimension(:,:), allocatable :: Atemp
    real(dp), dimension(:), allocatable :: btemp

    lsize = size(A,1)
    allocate( Atemp(lsize, lsize ) )
    Atemp = A
    allocate( btemp( lsize ) )
    btemp = b

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

    call dgesvd( 'A', 'A', lsize, lsize, A, lsize, &
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
        if (abs(S(i)) .lt. qfit_eps ) then
             S(i) = 0.0_dp
        else
             S(i) = 1.0_dp / S(i)
        endif
        if (qfit_debug) then
            write(luout,'(A,F16.9)') 'QFIT_DEBUG [linear_solve_svd]: S(i)', S(i)
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

    btemp = 0.0_dp
    call dgemv('N',nsize,nsize,1.0_dp,Sm,nsize,c,1,0.0_dp,btemp,1)
    call dgemv('T',nsize,nsize,1.0_dp,VT,nsize,btemp,1,0.0_dp,x(1:nsize),1)

    deallocate( Sm )
    deallocate( W )
    deallocate( ipiv )
    deallocate( c )
    deallocate( wrk )
    deallocate( S )
    deallocate( VT )
    deallocate( U )
    deallocate( btemp )
    deallocate( Atemp )

end subroutine

end module linear_solver
