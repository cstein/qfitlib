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
  !!        memory. This memory needs to be released by a call to connolly_finalize()
  !!
  !! @author Casper Steinmann
  !!
  !! @param R array of nuclear coordinates
  !! @param Z array of nuclear charges
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

end subroutine connolly_initialize

!------------------------------------------------------------------------------
!> @brief Releases the memory acquired by the connolly grid module.
!!
!! @author Casper Steinmann
subroutine connolly_finalize

    deallocate( max_layer_points )
    deallocate( n_layer_points )
    ! allocated in connolly_grid_count
    if (allocated( point_include) ) deallocate( point_include )
    deallocate( Rm )
    deallocate( Zm )

end subroutine connolly_finalize

!------------------------------------------------------------------------------
!> @brief Generates a connolly grid stored in the coordinates storage. Only the
!!        first ntruepoints are usable
!!
!! @author Casper Steinmann
!! @param coordinates storage for coordinates of a connolly grid. This must be
!!                    large enough to hold all points from all spheres
!! @param ntruepoints the coordinates 1 through ntruepoints are the true
!!                    coordinates to be used for the grid
subroutine connolly_grid( coordinates, ntruepoints )

    real(dp), dimension(:,:), intent(out) :: coordinates
    integer, intent(out) :: ntruepoints

    integer :: ilayer, ifrom, ito
    real(dp) :: rscal
    integer :: nlayerpoints, ioffs, ioffe

    coordinates = 0.0_dp

    ! loop over each layer, generating coordinates and storing them
    ! subsequently in the coordinates array.
    ioffs = 1
    do ilayer = 1, qfit_nshell
        rscal = qfit_vdwscale + (ilayer-1)*qfit_vdwincrement
        ioffe = ioffs + max_layer_points(ilayer)

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
            ito = ito +1
        endif
    enddo

    ! just to be sure that some stupid mistake did not happen
    if ( ntruepoints /= ito-1 ) then
        write(luout,*) "ERROR: grid points could not be transferred correctly"
        stop
    endif

end subroutine

!------------------------------------------------------------------------------
!> @brief Generates the grids for a single layer
!!
!! @author Casper Steinmann
!!
!! @param rscal scaling factor for the grid for the current layer
!! @param coordinates coordinates for @b all points on @b all spheres for
!!                        the current layer
!! @param pinc list of indices whether or not a point in coordinates
!!                    should be included or not
!! @param ntruepoints the number of non overlapping points.
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

    ! generate a unit sphere for each nuclei for a given layer
    do m = 1, nnuclei

        ! get ready to store the raw coordinates in temporary storage. Because
        ! of the density parameter this is not always constant.
        rmscal = vdw_radii( int(Zm(m)) ) * aa2au * rscal
        pp = int(4.0_dp * pi * rmscal * rmscal * qfit_pointdensity)
        call connolly_grid_sphere_count(pp, npoints)
        allocate( points(3, npoints) )

        ! generate a unit sphere, scale it and translate it
        call connolly_sphere( pp, points )
        points = points * rmscal
        do i = 1, npoints
            ! translate each point on the surface to where it belongs
            points(:,i) = Rm(:,m) - points(:,i)

            ! transfer to work array
            coordinates(:,i+ioffset) = points(:,i)
        enddo

        do i = 1, npoints

            ! check that a point is not within the boundary of neighbouring nuclei
            ncontact = 1
            do n = 1, nnuclei
                if (m == n) cycle

                rnscal = vdw_radii( int(Zm(n)) ) * aa2au * rscal
                dr  = points(:,i) - Rm(:,n)
                dr2 = dot(dr, dr)
                if (dr2 < rnscal**2) then
                    ncontact = ncontact +1
                endif
            enddo

            ! if the point is too close to a nuclei, discard it
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
!> @brief Creates a connolly unitsphere at (0,0,0).
!!
!! Surface tesselation is input through the connolly_initialize method
!!
!! @author Casper Steinmann
!!
!! @param ntes
!! @param coordinates coordinates of the connolly surface tesselation
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
        call complete_sphere(xy, r_full)

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
!! @param number of points on the sphere
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

!------------------------------------------------------------------------------
!> @brief evaluates the maximum amount of storage needed for the grid based on the
!!        current settings.
!!
!! @author Casper Steinmann
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

!------------------------------------------------------------------------------
subroutine complete_sphere(p, r)

    real(dp), intent(in) :: p
    real(dp), dimension(:,:), intent(inout) :: r

    integer :: i, n
    real(dp) :: angle
    real(dp) :: x, y

    n = size( r, 2 )

    do i = 1, n
        angle = 2.0_dp*pi*(i-1) / n
        x = cos( angle )
        y = sin( angle )

        r(1,i) = x*p
        r(2,i) = y*p
    enddo

end subroutine complete_sphere

end module connolly
