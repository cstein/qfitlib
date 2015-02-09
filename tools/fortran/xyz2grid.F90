program xyz2grid

    use connolly
    use qfit_precision
    use qfit_variables

    implicit none

    real(dp), dimension(:,:), allocatable :: wrk
    real(dp), dimension(:,:), allocatable :: R
    real(dp), dimension(:), allocatable :: nucz
    character(len=2) ,dimension(:), allocatable :: Z
    character(len=80) :: title
    integer :: ntotalpoints, ntruepoints, k
    integer :: natoms, ip

ip = 5
!open(ip, file='geom.xyz', status='old')
read(ip,*) natoms
read(ip,*) title

allocate( R(3, natoms ) )
allocate( Z(natoms) )
allocate( nucz(natoms) )
do k = 1, natoms
    read( ip, *) Z(k), R(:,k)
    if (trim(Z(k))=='H') nucz(k) =  1.0_dp
    if (trim(Z(k))=='C') nucz(k) =  6.0_dp
    if (trim(Z(k))=='N') nucz(k) =  7.0_dp
    if (trim(Z(k))=='O') nucz(k) =  8.0_dp
    if (trim(Z(k))=='Mg')nucz(k) = 12.0_dp
    if (trim(Z(k))=='S') nucz(k) = 16.0_dp
enddo

R = R * aa2au
qfit_nshell = 4
qfit_vdwscale = 1.4_dp
qfit_vdwincrement = 0.2_dp

call connolly_initialize( R, nucz )
call connolly_grid_count


ntotalpoints = sum(max_layer_points)
allocate( wrk( 3, ntotalpoints ) )
call connolly_grid( wrk, ntruepoints )

write(*,'(i4)') ntruepoints
write(*,*) 'AA'

do k = 1, ntruepoints
    write(*,'(a,f20.9,2f16.9)') 'X', wrk(:,k) / aa2au
enddo

call connolly_finalize

deallocate( R )
deallocate( Z )
deallocate( nucz )

end program xyz2grid

