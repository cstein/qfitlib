!------------------------------------------------------------------------------
!> @brief Module to generate a connolly grid
!!
!! @author Casper Steinmann
!!
!! @see M. L. Connolly, J. Appl. Crys. 16 (5): 548 - 558 (1983)
module connolly

    use qfit_variables
    use qfit_utilities, only : dot

    implicit none

    private

    ! constants
    real(dp), parameter :: EPS = 1.0d-8

    ! the tesselation level
    integer, save :: ntesselation

    ! the maximum number of points
    integer, save, allocatable, dimension(:) :: point_include

    public :: connolly_initialize
    public :: connolly_grid
    public :: connolly_finalize
    public :: connolly_grid_count

contains

!------------------------------------------------------------------------------
!> @brief Initialize the connolly grid generation module by specifying
!!        the tesselation level, nuclear coordinates and charges. Allocates
!!        memory. This memory needs to be released by a call to connolly_finalize
!!
!! @author Casper Steinmann
!!
!! @param[in] n tesselation level
!! @param[in] R array of nuclear coordinates
!! @param[in] Z array of nuclear charges
!! @param[out] nmaxpoints the maximum possible number of points required
subroutine connolly_initialize( R, Z )

    !integer, intent(in) :: n
    real(dp), dimension(:,:), intent(in) :: R
    real(dp), dimension(:), intent(in) :: Z

    nnuclei = size(Z)

    allocate( Zm(nnuclei) )
    Zm = Z

    allocate( Rm(3,nnuclei) )
    Rm = R

    allocate( max_layer_points( qfit_nshell ) )
    max_layer_points = 0

    allocate( n_layer_points( qfit_nshell ) )
    n_layer_points = 0

    ! points to include in the final coordinate list:
    ! we need to allocate [point_include] after settings have been read in
    ! see connolly grid count
    !allocate( point_include( max_layer_points ) )
    !point_include = 1

end subroutine

!------------------------------------------------------------------------------
!> @brief Releases the memory acquired by the connolly module.
!!
!! @author Casper Steinmann
subroutine connolly_finalize

    deallocate( max_layer_points )
    deallocate( n_layer_points )
    if (allocated( point_include) ) deallocate( point_include )
    deallocate( Rm )
    deallocate( Zm )

end subroutine connolly_finalize

subroutine connolly_grid( coordinates, ntruepoints )

    real(dp), dimension(:,:), intent(out) :: coordinates
    integer, intent(out) :: ntruepoints

    integer :: ilayer, ifrom, ito
    real(dp) :: rscal
    integer :: nlayerpoints, ioffs, ioffe

    coordinates = 0.0_dp
    ioffs = 1
    do ilayer = 1, qfit_nshell
        rscal = qfit_vdwscale + (ilayer-1)*qfit_vdwincrement
        ioffe = ioffs + max_layer_points(ilayer)

        !if (qfit_debug) then
        !    write(luout,*)
        !    write(luout,'(a,i6,f6.2,2i6)') 'calculating grid for layer [GRID]:', &
        !                               & ilayer, rscal, ioffs, ioffe-1
        !endif

        call connolly_grid_layer( rscal, coordinates(:,ioffs:ioffe-1), &
               & point_include(ioffs:ioffe-1), nlayerpoints )
        n_layer_points(ilayer) = nlayerpoints
        ioffs = ioffe
    enddo
    ntruepoints = sum( point_include )

    ! now, transfer the coordinates to keep (point_include(i) == 1)
    ! and shift them up in the coordinates array
    ito = 1
    do ifrom = 1,size(coordinates,2)
        if (point_include(ifrom) == 1) then
            coordinates(:,ito) = coordinates(:,ifrom)
            !if (qfit_debug) write(luout,'(a,3i4)') 'transfer [GRID]', ifrom, ito
            ito = ito +1
        endif
    enddo

    if ( ntruepoints /= ito-1 ) then
        write(luout,*) "ERROR: grid points could not be transferred correctly"
        stop
    endif

end subroutine

!------------------------------------------------------------------------------
!> @brief Generates the grids for the nuclei specified in connolly_init
!!
!! @author Casper Steinmann
!!
!! @param[in] rscal scaling factor for the grid
!! @param[in] coordinates coordinates for @b all points an @b all spheres. 
!! @param[out] n the number of non overlapping points. Obtained by a call to 
!!             connolly_coordinates
!! @todo convert to bohr
subroutine connolly_grid_layer( rscal, coordinates, pinc, ntruepoints )

    real(dp), intent(in) :: rscal
    real(dp), dimension(:,:), intent(out) :: coordinates
    integer, dimension(:), intent(inout) :: pinc
    integer, intent(out) :: ntruepoints


    integer :: m, i, n
    integer :: npoints, ipoint
    integer :: ncontact
    real(dp) :: rmscal, rnscal, dr2
    real(dp), dimension(3) :: dr
    real(dp), allocatable, dimension(:,:) :: points
    integer :: pp, ioffset

    ipoint = 0
    ntruepoints = 0
    ioffset = 0

    do m = 1, nnuclei
        rmscal = vdw_radii( int(Zm(m)) ) * aa2au * rscal
        pp = int(4.0_dp * pi * rmscal * rmscal * qfit_pointdensity)
        call connolly_grid_sphere_count(pp, npoints)
        !write(luout,'(a,i6,a,f6.2,3i6)') '   nuclei', m, ' Rm =', rmscal, pp, npoints, ioffset
        allocate( points(3, npoints) )

        ! generate a unit sphere, scale it and translate it
        call connolly_sphere( pp, points )
        points = points * rmscal
        do i = 1, npoints
            ! translate each point on the surface to where it belongs
            !write(luout,'(a,i4,3f12.6)') 'i',m,Rm(:,m)
            points(:,i) = Rm(:,m) - points(:,i)

            ! transfer to work array
            coordinates(:,i+ioffset) = points(:,i)
            !write(luout,'(a,2i6)') '      point', i, i+ioffset
        enddo

        do i = 1, npoints
        !    ! check that a point is not within the boundary of neighbouring nuclei
            ncontact = 1
            do n = 1, nnuclei
                if (m == n) cycle

                rnscal = vdw_radii( int(Zm(n)) ) * aa2au * rscal
                dr  = points(:,i) - Rm(:,n)
                dr2 = dot(dr, dr)
                if (dr2 < rnscal**2) then
        !            !write(luout,'(a,2i4,f9.6)') 'PP',m,n,dr2
                    ncontact = ncontact +1
                endif
            enddo

            if (ncontact > 1) then
                pinc(i+ioffset) = 0
            else
                ipoint = ipoint +1
            endif
        enddo
        ioffset = ioffset +npoints
        deallocate( points )
    enddo

    ntruepoints = ipoint

end subroutine connolly_grid_layer

!------------------------------------------------------------------------------
!> @brief Creates a connolly unitsphere at (0,0,0). Surface tesselation is
!!        input through the connolly_initialize method
!!
!! @author Casper Steinmann
!!
!! @param[out] coordinates coordinates of the connolly surface tesselation
subroutine connolly_sphere( ntes, coordinates )

    integer, intent(in) :: ntes
    real(dp), dimension(:,:), intent(out) :: coordinates

    integer :: neq, nvt, nbo
    integer :: i, io
    real(dp) :: angle
    real(dp) :: z, xy
    real(dp), allocatable, dimension(:,:) :: r_full

    neq = int(sqrt(pi * ntes ))
    nvt = neq / 2

    io = 1
    do i = 1,nvt+1
        angle = pi * (i-1) / nvt
        z  = cos( angle )
        xy = sin( angle )

        ! calculate number of points on the orthogonal circle
        nbo = int(xy*neq + EPS)
        if (nbo .lt. 1) nbo = 1

        ! allocate memory for the points
        allocate(r_full(3,nbo))

        ! assign z-value now. The subroutine connollu_orth_circle will populate x and y
        r_full(3,:) = z
        call full(xy, r_full)

        ! copy the calculated points on the sphere into the coordinate list
        coordinates(:,io:io+nbo-1) = r_full
        io = io + nbo

        deallocate(r_full)
    enddo

    ! double check that memory was not corrupted.
    if (io-nbo /= size(coordinates,2)) then
        write(*,'(a)') "ERROR: Connolly grid generation failed."
        stop
    endif

end subroutine connolly_sphere

!------------------------------------------------------------------------------
!> @brief Evaluates the number of gridpoints on a sphere given the current
!!        level of tesselation.
!!
!! @author Casper Steinmann
!!
!! @see the subroutine connolly_grid that generates the grid
!!
!! @param[out] number of points on the sphere
subroutine connolly_grid_sphere_count(npoints, nelems)

    integer, intent(in) :: npoints
    integer, intent(out) :: nelems

    integer :: nvt, neq, nbo
    integer :: i
    real(dp) :: angle, xy

    neq = int(sqrt(pi * npoints))
    nvt = neq / 2

    nelems = 0
    do i=1,nvt+1
        angle = pi * (i-1) / nvt
        xy = sin( angle )

        nbo = int(xy*neq + EPS)
        if (nbo .lt. 1) nbo = 1
        nelems = nelems + nbo
    enddo

end subroutine connolly_grid_sphere_count

subroutine connolly_grid_count

    integer :: ilayer, iatom, pp, np
    real(dp) :: rscal, radfin

    do ilayer = 1, qfit_nshell
        rscal = qfit_vdwscale + (ilayer-1)*qfit_vdwincrement
        do iatom = 1, size(Zm)
            radfin = vdw_radii( int(Zm(iatom)) ) * aa2au * rscal
            pp = int(4.0_dp * pi * radfin * radfin * qfit_pointdensity)
            call connolly_grid_sphere_count( pp, np )
            max_layer_points(ilayer) = max_layer_points(ilayer) + np
        enddo
    enddo

    allocate( point_include( sum(max_layer_points) ) )
    point_include = 1

end subroutine connolly_grid_count

subroutine full(p, r)

    real(dp), intent(in) :: p
    real(dp), dimension(:,:), intent(inout) :: r

    integer :: i, n
    real(dp) :: angle
    real(dp) :: x, y

    n = size( r, 2 )

    do i = 1, n
        angle = 2*pi*(i-1) / n
        x = cos( angle )
        y = sin( angle )

        r(1,i) = x*p
        r(2,i) = y*p
    enddo

end subroutine full

end module connolly
