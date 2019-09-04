!------------------------------------------------------------------------------
!> @brief Integral interface
!!
!! @author Casper Steinmann
module qfit_integrals

    use qfit_precision

    implicit none

    private

    public :: one_electron_integrals

    contains

!------------------------------------------------------------------------------
!> @brief evaluates the potential integrals
!!
!! @author Casper Steinmann
!!
!! @param charges Charges of the nuclei.
!! @param coord Coordinate of the surface point \f$\mathbf{R}_A\f$ where the potential integrals are evaluated.
!! @param integrals The resulting one-electron potential integrals.
!!
!! @details
!! @f[
!!   V_{\mu\nu} = \big \langle \mu \big| \frac{Z_A}{|\mathbf{r}-\mathbf{R}_A|} \big| \nu \big\rangle
!! @f]
subroutine one_electron_integrals(charges, coords, integrals)

    real(dp), intent(in), dimension(:) :: charges
    real(dp), intent(in), dimension(3) :: coords
    real(dp), dimension(:), intent(out) :: integrals

    integrals = 0.0_dp

#if defined (GEN1INT)
    call one_e_integrals_gen1int(charges, coords, integrals)
#endif

end subroutine one_electron_integrals

#if defined (GEN1INT)

subroutine one_e_integrals_gen1int(charges, coords, integrals)

    use gen1int_api

    real(dp), intent(in), dimension(:) :: charges
    real(dp), intent(in), dimension(3) :: coords
    real(dp), dimension(:), intent(out) :: integrals

    ! -- GEN1INT VARIABLES --
    integer :: ierr
    integer :: prop_sym
    integer :: num_prop
    integer :: num_ao
    integer :: num_geo_bra
    integer :: num_geo_ket
    integer :: num_geo_total
    integer :: io
    integer :: printlvl = 0
    logical :: symmetric
    logical :: triangular
    type(one_prop_t) :: prop_operator
    type(nary_tree_t) :: nary_tree_bra
    type(nary_tree_t) :: nary_tree_ket
    type(nary_tree_t) :: nary_tree_total
    type(matrix), dimension(:), allocatable :: intmats

    ! non-zero components for the operator, the first dimension is for bra and
    ! ket sub-shells, the last is the number of non-zero components, which should
    ! be 1 for non-relativistic calcualtions
    integer nnz_comp(2,1)

    ! -- SUBROUTINE VARIABLES --
    integer :: i, j
    integer :: nat
    integer, dimension(:), allocatable :: nuclei
    real(dp), dimension(:), allocatable :: real_charges
    real(dp), dimension(3,1) :: real_coords

    integrals = 0.0_dp

    ! treat the nuclei as non-nuclei because we want to also
    ! be able to evaluate other types of atomic charges
    real_coords(:,1) = coords
    nat = size(charges)
    allocate(nuclei(nat))
    allocate(real_charges(nat))
    nuclei = -1

    ! sets the non-zero components for the one-electron operator, here we only
    ! consider the (large,large) part
    nnz_comp(1,1) = 1
    nnz_comp(2,1) = 1

    ! do the correct sign according to the operator
    real_charges = charges

    call OnePropCreate(prop_name=INT_POT_ENERGY, &
                       one_prop=prop_operator,   &
                       info_prop=ierr,           &
                       idx_nuclei=nuclei,        &
                       coord_nuclei=real_coords,      &
                       charge_nuclei=real_charges,     &
                       order_geo_pot=0)

    if (ierr /= 0) stop 'ERROR: Failed to create property operator'

    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_total)

    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, &
                           num_prop=num_prop)
    call OnePropGetSymmetry(one_prop=prop_operator, &
                            prop_sym=prop_sym)

    ! gets the number of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)

    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    !write(luout,*) "gen1int: num_props", num_prop, num_ao

    if (num_prop /= 1) stop 'ERROR: Integral property failed.'

    triangular = .false.
    symmetric = (prop_sym == SYMM_INT_MAT)

    allocate(intmats(num_prop), stat=ierr)

    call MatAssociate(work_alpha=integrals(:), &
                      num_row=num_ao,          &
                      A=intmats(num_prop),     &
                      info_mat=ierr,           &
                      triangular=triangular,   &
                      symmetric=symmetric)

    if (ierr /= 0) stop 'ERROR: Failed to associate matrices.'

    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp, &
                                  one_prop=prop_operator, &
                                  nary_tree_bra=nary_tree_bra, &
                                  nary_tree_ket=nary_tree_ket, &
                                  nary_tree_total=nary_tree_total, &
                                  num_ints=num_prop, &
                                  val_ints=intmats, &
                                  num_dens=1,       &
                                  io_viewer=io,       &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)

    call MatNullify(A=intmats(num_prop))

    deallocate(intmats)
    deallocate(nuclei)
    deallocate(real_charges)

end subroutine one_e_integrals_gen1int

#endif

end module qfit_integrals
