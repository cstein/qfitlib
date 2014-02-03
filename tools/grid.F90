program test

    use connolly
    use qfit_precision

    implicit none

    integer :: i, n, nelems, npoints, ntes
    real(dp), dimension(3,3) :: R
    integer, dimension(3) :: Z
    real(dp), allocatable, dimension(:,:) :: wrk
    real(dp), allocatable, dimension(:,:) :: coordinates

R(:,1) = (/0.0, 0.0, 0.0/)
R(:,2) = (/0.9, 0.0, 0.0/)
R(:,3) = (/-0.36, 0.8, 0.0/)
Z = (/8, 1, 1/)


ntes = 16


do n = 1, 1

    call connolly_initialize( ntes, R, Z, nelems )
    
    allocate(wrk(3,nelems))

    call connolly_grid( 1.4d0 + (n-1)*0.2d0, wrk, npoints )
    
    allocate(coordinates(3, npoints))
    
    call connolly_coordinates( wrk, coordinates )
    
    print '(i4)', npoints
    print *, ''
    do i = 1, npoints
        print '(a,3f9.4)','X', coordinates(1, i), &
                             & coordinates(2, i), &
                             & coordinates(3, i)
    enddo
    
    deallocate(coordinates)
    deallocate(wrk)

    call connolly_finalize()

enddo


end program test

