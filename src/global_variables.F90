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
    real(dp), save, allocatable, dimension(:) :: Zm
    !> nuclear coordinates
    real(dp), save, allocatable, dimension(:,:) :: Rm
    !> the number of nuclei
    integer, save :: nnuclei
    !> total charge of the molecule
    integer, save :: total_charge
    !> total dipole of the molecule
    real(dp), dimension(3), save :: total_dipole
    !> center of mass of the molecule
    real(dp), dimension(3), save :: center_of_mass

    !
    ! Run-time settings to be read in from input
    !
    !> unit to write to for output (default is stdout)
    integer, save :: luout = 6

    !> bitwise additive option for constraints
    !> 0: nothing
    !> 1: charges
    !> 2: dipole
    !> 4: quadrupole
    !> default: 1. To select charges and dipoles, use 1+2 = 3
    integer, save :: qfit_constraint = 1
    !> the scaling factor for the van der waal radii
    real(dp), save :: qfit_vdwscale = 1.4_dp
    !> the increment in scaling factor for each layer
    real(dp), save :: qfit_vdwincrement = 0.2_dp
    !> the number of layers to include
    integer, save :: qfit_nshell = 4
    !> the density of points on the sphere
    real(dp), save :: qfit_pointdensity = 0.28_dp
    !> optional file on which we are to evaluate the mep
    character(len=80), save :: qfit_mepfile = ''
    !> whether or not to only evaluate the molecular electrostatic potential.
    !> this will skip any fitting
    logical, save :: qfit_only_calculate_mep = .false.

    !
    ! runtime variables
    !
    !> the maximum number of points in a layer
    integer, save, allocatable, dimension(:) :: max_layer_points
    !> the number of points that are in a layer
    integer, save, allocatable, dimension(:) :: n_layer_points
    !> the total number of points for all layers
    integer, save :: n_total_points
    !> the resulting charges
    real(dp), save, allocatable, dimension(:) :: fitted_charges

    ! constants
    !> @f$ \pi @f$
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: aa2au = one / 0.5291772109217_dp
    ! van der waal radii for the elements
    real(dp), parameter, dimension(0:16) :: vdw_radii = (/ &
        & 1.2d0, &
        & 1.2d0, 0.0d0, 0.0d0, 0.0d0, &
        & 0.0d0, 1.5d0, 1.5d0, 1.4d0, &
        & 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        & 0.0d0, 0.0d0, 0.0d0, 1.89d0 /)


end module qfit_variables
