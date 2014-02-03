!> @brief module that holds all QFITLIB-specific variables
!!
!! @author Casper Steinmann
module qfit_variables

    use qfit_precision

    implicit none

    !> whether or not we are running MBLIB.
    logical, save :: qfitrun = .false.
    !> whether or not to be more verbose with output
    logical, save :: qfit_verbose = .false.
    !> whether or not to print debug output.
    logical, save :: qfit_debug = .false.

    !
    ! Storage of nuclear coordinates and charges
    !
    !> nuclear charges
    integer, save, allocatable, dimension(:) :: Zm
    !> nuclear coordinates
    real(dp), save, allocatable, dimension(:,:) :: Rm

    !
    ! Run-time settings to be read in from input
    !
    !> unit to write to for output (default is stdout)
    integer, save :: luout = 6
    !> bitwise additive option for constraints
    !> 0: charges
    !> 1: dipole
    !> 2: quadrupole
    !> default: 0. To select charges and dipoles, use 0+1 = 1
    integer, save :: qfit_constraint = 0
    !> the scaling factor for the van der waal radii
    real(dp), save :: qfit_vdwscale = 1.4_dp
    !> the increment in scaling factor for each layer
    real(dp), save :: qfit_vdwincrement = 0.2_dp
    !> the number of layers to include
    integer, save :: qfit_nlayer = 4
    !> the density of points on the sphere
    real(dp), save :: qfit_pointdensity = 0.28_dp

    ! constants
    !> @f$ \pi @f$
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: aa2au = 1.0_dp / 0.5291772109217_dp
    ! van der waal radii for the elements
    real(dp), parameter, dimension(16) :: vdw_radii = (/ &
        & 1.2d0, 0.0d0, 0.0d0, 0.0d0, &
        & 0.0d0, 1.5d0, 1.5d0, 1.4d0, &
        & 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        & 0.0d0, 0.0d0, 0.0d0, 1.89d0 /)


end module qfit_variables
