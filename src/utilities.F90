module qfit_utilities

    use qfit_precision

    contains

!------------------------------------------------------------------------------
!> @author Magnus Olsen
!
!> @brief changes case on input string. Shamelessly stolen from Polarizable Embedding module
!
!> @param[in,out] string
!> @param[in] uplo
subroutine change_case(string, uplo)

    character(len=*), intent(inout) :: string
    character(len=*), intent(in), optional :: uplo

    integer :: i, gap
    character(len=1) :: a, z, o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    end if

    gap = iachar('a') - iachar('A')

    if (o_uplo == 'u' .or. o_uplo == 'U') then
        a = 'a'
        z = 'z'
    else if (o_uplo == 'l' .or. o_uplo == 'L') then
        a = 'A'
        z = 'Z'
        gap = - gap
    else
        stop 'Unknown case specified'
    end if

    do i = 1, len_trim(string)
        if (lge(string(i:i), a) .and. lle(string(i:i), z)) then
            string(i:i) = achar(iachar(string(i:i)) - gap)
        end if
    end do

end subroutine change_case

function dot(x,y)

    real(dp), external :: ddot

    real(dp) :: dot
    real(dp), dimension(:), intent(in) :: x, y

    integer :: n, incx, incy

    incx = 1
    incy = 1

    n = size(x)

    dot = ddot(n, x, incx, y, incy)

end function dot

subroutine gemm(a, b, c, transa, transb, alpha, beta)

    external :: dgemm

    real(dp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: transa, transb
    real(dp), dimension(:,:), intent(in) :: a, b
    real(dp), dimension(:,:) , intent(inout) :: c

    integer :: m, n, k, lda, ldb, ldc
    character(len=1) :: o_transa, o_transb
    real(dp) :: o_alpha, o_beta

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0_dp
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0_dp
    end if

    if (present(transa)) then
        o_transa = transa
    else
        o_transa = 'N'
    end if

    if (present(transb)) then
        o_transb = transb
    else
        o_transb = 'N'
    end if

    if (o_transa == 'N') then
        k = size(a, 2)
    else
        k = size(a, 1)
    end if

    m = size(c, 1)
    n = size(c, 2)
    lda = max(1, size(a, 1))
    ldb = max(1, size(b, 1))
    ldc = max(1, size(c, 1))

    call dgemm(o_transa, o_transb, m, n, k, o_alpha, a(1,1), lda, b(1,1), ldb, o_beta, c(1,1), ldc)

end subroutine gemm

function trace( a )
    real(dp), intent(in), dimension(:,:) :: a
    real(dp) :: trace
    integer :: i, n
    n = size( a, 1 )
    trace = 0.0_dp
    do i = 1, n
        trace = trace + a(i,i)
    enddo
end function trace

end module qfit_utilities
