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

    ! the number of nuclei in the molecule
    integer, save :: nnuclei

    ! the tesselation level
    integer, save :: ntesselation

    ! the maximum number of points
    !integer, save :: nmaxpoints
    integer, save, allocatable, dimension(:) :: point_include

    public :: connolly_initialize
    public :: connolly_grid
    public :: connolly_coordinates
    public :: connolly_finalize

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
subroutine connolly_initialize( n, R, Z, nmaxpoints )

    integer, intent(in) :: n
    real(dp), dimension(:,:), intent(in) :: R
    integer, dimension(:), intent(in) :: Z
    integer, intent(out) :: nmaxpoints

    integer :: npoints

    ntesselation = n

    nnuclei = size(Z)
    call connolly_grid_count( npoints )

    nmaxpoints = nnuclei * npoints

    allocate( Zm(nnuclei) )
    Zm = Z

    allocate( Rm(3,nnuclei) )
    Rm = R

    ! points to include in the final coordinate list:
    ! 1 is to include, 0 is to not include
    allocate( point_include( nmaxpoints ) )
    point_include = 1

end subroutine

!------------------------------------------------------------------------------
!> @brief Releases the memory acquired by the connolly module.
!!
!! @author Casper Steinmann
subroutine connolly_finalize

    deallocate( point_include )
    deallocate( Rm )
    deallocate( Zm )

end subroutine connolly_finalize

!------------------------------------------------------------------------------
!> @brief Generates the grids for the nuclei specified in connolly_init
!!
!! @author Casper Steinmann
!!
!! @param[in] rscal scaling factor for the grid
!! @param[in] coordinates coordinates for @b all points an @b all spheres. 
!! @param[out] n the number of non overlapping points. Obtained by a call to 
!!             connolly_coordinates
subroutine connolly_grid( rscal, coordinates, ntruepoints )

    real(dp), intent(in) :: rscal
    real(dp), dimension(:,:), intent(out) :: coordinates
    integer, intent(out) :: ntruepoints


    integer :: m, i, n
    integer :: npoints, ipoint
    integer :: ncontact
    real(dp) :: rmscal, rnscal, dr2
    real(dp), dimension(3) :: dr
    real(dp), allocatable, dimension(:,:) :: points

    ! allocate workspace for coordinate generation over
    ! a single sphere
    call connolly_grid_count(npoints)
    allocate( points(3, npoints) )

    ipoint = 0

    do m = 1, nnuclei

        ! generate a unit sphere, scale it and translate it
        call connolly_sphere( points )
        rmscal = vdw_radii( Zm(m) ) * rscal
        points = points * rmscal
        do i = 1, npoints
            points(:,i) = Rm(:,m) - points(:,i)

            ! transfer to work array
            coordinates(:,i+(m-1)*npoints) = points(:,i)
        enddo

        do i = 1, npoints
            ! check that a point is not within the boundary of neighbouring nuclei
            ncontact = 1
            do n = 1, nnuclei
                if (m == n) cycle

                rnscal = vdw_radii( Zm(n) ) * rscal
                dr  = points(:,i) - Rm(:,n)
                dr2 = dot(dr, dr)
                if (dr2 < rnscal**2) then
                    ncontact = ncontact +1
                endif
            enddo

            if (ncontact > 1) then
                point_include(i+(m-1)*npoints) = 0
            else
                ipoint = ipoint +1
            endif
        enddo
    enddo

    deallocate( points )
    ntruepoints = ipoint

end subroutine connolly_grid

!------------------------------------------------------------------------------
!> @brief Copies the non-overlapping coordinates from wrk into coordinates
!!
!! @author Casper Steinmann
!!
!! @param[in] wrk all coordinates returned by connolly_grid
!! @param[out] coordinates the non-overlapping coordinates
subroutine connolly_coordinates( wrk, coordinates )

    real(dp), dimension(:,:), intent(in) :: wrk
    real(dp), dimension(:,:), intent(out) :: coordinates

    integer :: i,j

    j = 1
    do i = 1, size(wrk, 2)
        if (point_include(i) == 1) then
            coordinates(:,j) = wrk(:,i)
            j = j + 1
        endif
    enddo

end subroutine connolly_coordinates

!------------------------------------------------------------------------------
!> @brief Creates a connolly unitsphere at (0,0,0). Surface tesselation is
!!        input through the connolly_initialize method
!!
!! @author Casper Steinmann
!!
!! @param[out] coordinates coordinates of the connolly surface tesselation
subroutine connolly_sphere( coordinates )

    real(dp), dimension(:,:), intent(out) :: coordinates

    integer :: neq, nvt, nbo
    integer :: i, io
    real(dp) :: angle
    real(dp) :: z, xy
    real(dp), allocatable, dimension(:,:) :: r_full

    neq = int(sqrt(pi * ntesselation ))
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
subroutine connolly_grid_count(nelems)

    integer, intent(out) :: nelems

    integer :: nvt, neq, nbo
    integer :: i
    real(dp) :: angle, xy

    neq = int(sqrt(pi * ntesselation ))
    nvt = neq / 2

    nelems = 0
    do i=1,nvt+1
        angle = pi * (i-1) / nvt
        xy = sin( angle )

        nbo = int(xy*neq + EPS)
        if (nbo .lt. 1) nbo = 1
        nelems = nelems + nbo
    enddo

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
