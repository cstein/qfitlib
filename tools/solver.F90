program test

    use qfit_precision
    use linear_solver

    implicit none

    real(dp), dimension(3,3) :: A
    real(dp), dimension(3) :: b
    real(dp), dimension(3) :: x

A(1,:) = (/-3.0d0, 2.0d0, -6.0d0/)
A(2,:) = (/ 5.0d0, 7.0d0, -5.0d0/)
A(3,:) = (/ 1.0d0, 4.0d0, -2.0d0/)
b(:) = (/6.0d0, 6.0d0, 8.0d0/)

call linear_solve_svd(A,b,x)

print '(3f9.4)', x

end program test

