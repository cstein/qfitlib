program test

    !use connolly
    use qfit
    use qfit_precision
    use qfit_variables

    implicit none

    real(dp), dimension(1) :: dum
    real(dp), dimension(3,3) :: R
    real(dp), dimension(3) :: Z

R(:,1) = (/0.0, 0.0, 0.0/)
R(:,2) = (/0.9, 0.0, 0.0/)
R(:,3) = (/-0.36, 0.8, 0.0/)
Z = (/8.0, 1.0, 1.0/)

R = R * aa2au

dum = 1.0_dp

!call qfit_initialize( R, Z )
!call fit_density( dum )
!call qfit_finalize()


end program test

