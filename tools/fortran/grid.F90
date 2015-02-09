program test

    use connolly
    use qfit_precision
    use qfit_variables

    implicit none

    real(dp), dimension(:,:), allocatable :: wrk
    real(dp), dimension(3,3) :: R
    real(dp), dimension(3) :: Z
    integer :: ntotalpoints, ntruepoints, k

R(:,1) = (/0.0, 0.0, 0.0/)
R(:,2) = (/0.9, 0.0, 0.0/)
R(:,3) = (/-0.41, 0.81, 0.0/)
Z = (/8.0, 1.0, 1.0/)

R = R * aa2au

qfit_nshell = 4

call connolly_initialize( R, Z )
call connolly_grid_count


ntotalpoints = sum(max_layer_points)
allocate( wrk( 3, ntotalpoints ) )
call connolly_grid( wrk, ntruepoints )

write(*,'(i4)') ntruepoints
write(*,*)

do k = 1, ntruepoints
    write(*,'(a,f20.9,2f16.9)') 'X', wrk(:,k) / aa2au
enddo

call connolly_finalize


end program test

