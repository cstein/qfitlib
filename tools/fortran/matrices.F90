program matrices

    use qfit
    use qfit_precision
    use qfit_variables

    implicit none

    real(dp), dimension(:,:), allocatable :: wrk
    real(dp), dimension(:,:), allocatable :: R
    real(dp), dimension(:), allocatable :: nucz
    real(dp), dimension(:), allocatable :: mass
    real(dp), dimension(:), allocatable :: dummy_dens
    real(dp), dimension(3) :: cm
    real(dp), dimension(3) :: d
    character(len=2) ,dimension(:), allocatable :: Z
    character(len=80) :: title
    integer :: ntotalpoints, ntruepoints, k
    integer :: natoms, ip
    real(dp) :: tmass

ip = 5
!open(ip, file='geom.xyz', status='old')
read(ip,*) natoms
read(ip,*) title

cm = zero
d = (/0.50839, 0.47158, 0.54270/)
tmass = zero

allocate( R(3, natoms ) )
allocate( Z(natoms) )
allocate( nucz(natoms) )
allocate( mass(natoms) )

do k = 1, natoms
    read( ip, *) Z(k), R(:,k)
    if (trim(Z(k))=='H') nucz(k) =  1.0_dp
    if (trim(Z(k))=='H') mass(k) =  1.008_dp
    if (trim(Z(k))=='C') nucz(k) =  6.0_dp
    if (trim(Z(k))=='C') mass(k) =  12.011_dp
    if (trim(Z(k))=='N') nucz(k) =  7.0_dp
    if (trim(Z(k))=='N') mass(k) =  14.007_dp
    if (trim(Z(k))=='O') nucz(k) =  8.0_dp
    if (trim(Z(k))=='O') mass(k) =  15.999_dp
    if (trim(Z(k))=='S') nucz(k) = 16.0_dp
    if (trim(Z(k))=='S') mass(k) = 32.06_dp
    if (trim(Z(k))=='Mg')nucz(k) = 12.0_dp
    if (trim(Z(k))=='Mg') mass(k) = 24.06_dp
    cm = cm + R(:,k)*mass(k)
    tmass = tmass + mass(k)
enddo
cm = cm / tmass

R = R * aa2au

qfit_debug = .true.
qfit_multipole_rank = 0
qfit_vdwscale = 1.0
qfit_nshell = 4
qfit_constraint = 3
!qfit_mepfile = 'surface.mep'
call qfit_initialize( R, nucz, 0, d, cm )
call qfit_print_info()
allocate( dummy_dens(144) )
call qfit_fit( dummy_dens )
deallocate( dummy_dens )

call qfit_finalize()

deallocate( nucz )
deallocate( Z )
deallocate( R )

end program matrices

