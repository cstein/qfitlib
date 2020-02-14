!------------------------------------------------------------------------------
!> @brief Charge fitting library
!!
!! @author Casper Steinmann
module qfit

    use qfit_precision
    use qfit_variables
    use qfit_io, only : openfile


    implicit none

    external :: output

    private

    public :: qfit_initialize
    public :: qfit_finalize
    public :: qfit_set_dipole
    public :: qfit_set_transition_dipole
    public :: qfit_print_info
    public :: qfit_fit
    public :: qfit_get_results

    contains

!------------------------------------------------------------------------------
!> @brief Initializes the charge fitting library
!!
!! @author Casper Steinmann
!! @param R A \f$(3,N)\f$ array containing the coordinates of all nuclei
!! @param Z An \f$(N)\f$ array containing the nuclear charges
!! @param Q Total charge of the molecule. Default value, unless specified, is \f$Q = 0\f$.
!! @param mu Total dipole of the molecule. Default value, unless specified, is \f$\bar{\mu}=(0.0,0.0,0.0)\f$.
!! @param RCM center of mass of the molecule
subroutine qfit_initialize(R, Z, Q, mu, RCM)

    use connolly

    real(dp), dimension(:,:), intent(in) :: R
    real(dp), dimension(:), intent(in) :: Z
    integer, intent(in), optional :: Q
    real(dp), dimension(3), intent(in), optional :: mu
    real(dp), dimension(3), intent(in), optional :: RCM

    total_charge = 0
    total_dipole = zero
    center_of_mass = zero

    call connolly_initialize( R, Z )

    if (present(Q)) total_charge = Q
    if (present(mu)) total_dipole = mu
    if (present(RCM)) center_of_mass = RCM

    allocate( fitted_charges(size(Z)) )
    allocate( fitted_dipoles(3*size(Z)) )
    allocate( fitted_quadrupoles(6*size(Z)) )
    fitted_charges = zero
    fitted_dipoles = zero
    fitted_quadrupoles = zero


end subroutine

!------------------------------------------------------------------------------
!> @brief releases memory back to the program
!!
!! @author Casper Steinmann
subroutine qfit_finalize

    use connolly

    deallocate( fitted_quadrupoles )
    deallocate( fitted_dipoles )
    deallocate( fitted_charges )
    call connolly_finalize

end subroutine

!------------------------------------------------------------------------------
!> @brief sets the total charge of the system
!!
!! @author Casper Steinmann
!!
!! @details sets the constraint on the total charge of the system
!! \f$q_\mathrm{tot}=\sum_i q_i\f$
subroutine qfit_set_total_charge(charge)
    integer, intent(in) :: charge
    total_charge = charge
end subroutine qfit_set_total_charge

!------------------------------------------------------------------------------
!> @brief Assigns a dipole to QFIT which can be used to enforce constraints
!!
!! @author Casper Steinmann
!! @param mu The dipole to use as a constraint.
subroutine qfit_set_dipole(mu)
    real(dp), dimension(3), intent(in) :: mu
    total_dipole = mu
end subroutine qfit_set_dipole

!------------------------------------------------------------------------------
!> @brief Sets a transition dipole to the total dipole.
!!
!! Also performs the action of setting the total charge \f$Q\f$ to zero
!! and setting the nuclear charges to zero.
!!
!! @author Casper Steinmann
!! @param trdip The transition dipole to use as a constraint.
subroutine qfit_set_transition_dipole( trdip )
    real(dp), dimension(3), intent(in) :: trdip

    total_charge = 0
    total_dipole = -trdip
    Zm = zero

end subroutine

!------------------------------------------------------------------------------
!> @brief print out information about the calculation
!!
!! @author Casper Steinmann
subroutine qfit_print_info

    use qfit_variables

    if (qfit_multipole_rank .eq. 2) then
        write(luout, 16) 'charges, dipoles and quadrupoles'
    elseif (qfit_multipole_rank .eq. 1) then
        write(luout, 16) 'charges and dipoles'
    elseif (qfit_multipole_rank .eq. 0) then
        write(luout, 16) 'charges'
    else
        write(luout, '(/10x,a,i0,a)') 'WARNING: multipole order ',qfit_multipole_rank,' not understood.'
    endif

    if (trim(qfit_mepfile) /= '') then
        write(luout, 13) trim(qfit_mepfile)
    else
        write(luout, 10) qfit_nshell, qfit_vdwscale, qfit_vdwincrement
        write(luout, 11) qfit_pointdensity
    endif
    if (qfit_constraint .eq. -1) then
        write(luout,'(/10x,a)') 'WARNING: No constraints on charges imposed.'
    else
        if (qfit_constraint.eq.0) write(luout, 14)
        if (qfit_constraint.eq.1 .and. qfit_multipole_rank .eq. 0) write(luout, 15)
    endif
    if (qfit_verbose) write(luout,'(/10x,a)') 'Verbose mode enabled.'
    if (qfit_debug) write(luout,'(/10x,a)') 'Debug mode enabled.'

 10 format(/10x,'Generating surface using ', i2, ' layers. Each layer is'/, &
   & 10x,'scaled by', f4.1, ' plus ', f4.1,' for each successive layer.')
 11 format(/10x,'Point density is ', f5.2, ' au^-2.')
 13 format(/10x,'Will read "',a,'" for surface points.')
 14 format(/10x,'Constraining fitted charges to reproduce the', &
   &       /10x,'molecular charge.')
 15 format(/10x,'Constraining fitted charges to reproduce the', &
           /10x,'molecular dipole.')

 16 format(/10x,'Fitting ', A,' to the molecular ESP.')

end subroutine

!------------------------------------------------------------------------------
!> @brief Returns the resulting potential derived moments.
!!
!! The caller provides storage of size \f$N\f$ for the charges and
!! of size \f$3N\f$ for the dipoles where \f$N\f$ is the number of
!! atoms used upon input in the qfit_initialize() call.
!!
!! The dipole components are organized in the canonical order,
!! i.e. as \f$\mu_{1,x},\mu_{1,y},\mu_{1,z}, \mu_{2,x}, \ldots, \mu_{N,z}\f$.
!!
!! @author Casper Steinmann
!! @param charges Resulting potential derived charges.
!! @param dipoles Resulting potential derived dipoles.
!! @param quadrupoles Resulting potential derived dipoles.
subroutine qfit_get_results( charges, dipoles, quadrupoles, rmsd )

    real(dp), dimension(:), intent(out) :: charges
    real(dp), dimension(:), intent(out), optional :: dipoles
    real(dp), dimension(:), intent(out), optional :: quadrupoles
    real(dp), intent(out), optional :: rmsd

    if (size(charges) /= nnuclei) then
        write(luout,'(/A)') "ERROR: Memory allocation in input to qfit_get_results"
        write(luout,'(A,I4,A,I4)') " is wrong. Expected:", nnuclei, " got:", size(charges)
        stop
    endif
    charges = fitted_charges(1:nnuclei)
    if (qfit_multipole_rank >= 1 .and. present(dipoles)) then
        dipoles = fitted_dipoles(1:3*nnuclei)
        if (qfit_multipole_rank >= 2 .and. present(quadrupoles)) then
            quadrupoles = fitted_quadrupoles(1:6*nnuclei)
        endif
    endif
    rmsd = rmsd_esp_fit

end subroutine

!------------------------------------------------------------------------------
!> @brief Fit charges or multipole moments to the molecular esp on a number of points.
!!
!! @author Casper Steinmann
!! @param density The square AO density of the molecule.
subroutine qfit_fit(density, lupri)

    use connolly
    use qfit_integrals
    use qfit_utilities
    use linear_solver
    use auxmat, only : geom_mat, constr_mat
    use potential, only : rmsd, eval_potential

#if defined(VAR_MPI)
#if defined(USE_MPI_MOD_F90)
    use mpi
#else
    include 'mpif.h'
#endif
#endif
    ! input arguments
    real(dp), dimension(:), intent(in) :: density
    integer, optional, intent(in) :: lupri

    ! local arrays
    real(dp), dimension(:,:), allocatable :: wrk
    real(dp), dimension(:), allocatable :: integrals
    real(dp), dimension(:), allocatable :: V, Vwrk, vcharges, vdipoles, vquadrupoles
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: B

    ! local variables
    real(dp) :: Rmk
    real(dp), dimension(3) ::  dr, drhat, mu
    real(dp), dimension(6) ::  oo
    real(dp), dimension(1) :: q_one
    real(dp), dimension(:), allocatable :: moments
    real(dp), dimension(:), allocatable :: mom_pot
    real(dp), dimension(:), allocatable :: pot_diff
    real(dp) :: pot_sum
    integer :: ntotalpoints, ntruepoints
    integer :: n, m, k
    integer :: ioff, n2bas, nbas
    integer :: matsiz, matdim, condim
    logical :: lunit_open
    integer :: lumepinp
    integer :: icf, ict, idf, idt, iqf, iqt
    real(dp) :: factor
    character(len=2) :: option

#if !defined(VAR_MPI)
    integer, save :: master, myid, nprocs
    integer :: wrk_size
    integer, allocatable, dimension(:) :: wrk_sizes, displs

    myid = 0
    master = myid
    nprocs = 1
#else
    ! mpi specific variables
    integer(kind=MPI_INTEGER_KIND), save :: master, myid, nprocs, ierr
    integer(kind=MPI_INTEGER_KIND) :: wrk_size
    integer(kind=MPI_INTEGER_KIND), allocatable, dimension(:) :: wrk_sizes, displs

    if (present(lupri)) luout = lupri
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
    master = 0

    ! things to broadcast used below
    call mpi_bcast(nnuclei, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (myid.ne.master) then
        allocate(Zm(nnuclei), Rm(3,nnuclei))
    endif
    call mpi_bcast(Zm,  nnuclei, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    call mpi_bcast(Rm,3*nnuclei, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
#endif

    ! Either we generate a grid and dump it to a file so we can read
    ! it back later when fitting for many densities or we read in
    ! a supplied file by the user

    if (myid == master) then
        ! check if the user (or we) supplied a custom file to be
        ! used to evaluate the molecular electrostatic potential ...
        if (trim(qfit_mepfile) /= '' ) then

            call openfile(trim(qfit_mepfile), lumepinp, 'old', 'formatted')
            rewind( lumepinp )

            ! format is classic xyz. title line MUST hold AU or AA
            read(lumepinp,*) ntotalpoints
            ntruepoints = ntotalpoints
            read(lumepinp,*) option

            factor = one
            if (trim(option) == 'AA') then
                factor = aa2au
            endif
            allocate( wrk( 3, ntruepoints) )
            do m=1, ntruepoints
                read(lumepinp,*) option, wrk(:,m)
            enddo
            wrk = wrk * factor
        else
            !
            ! ... or create a grid around the molecule
            !

            ! populates the max_layer_points array for memory management
            ! this requires that options have been set.
            call connolly_grid_count

            ! calculate the grid surrounding the molecule
            ntotalpoints = sum(max_layer_points)
            allocate( wrk( 3, ntotalpoints ) )
            call connolly_grid( wrk, ntruepoints )

            qfit_mepfile = 'surface.mep'
            call openfile(qfit_mepfile, lumepinp, 'new', 'formatted')
            write(lumepinp,*) ntruepoints
            write(lumepinp,*) 'AU'
            do m=1,ntruepoints
                write(lumepinp,*) 'X', wrk(:,m)
            enddo
        endif

        ! allocate space on master to store entire potential
        allocate(V(ntruepoints))
        V = zero

        ! ready arrays for scatter(v) and gather(v)
        allocate(wrk_sizes(nprocs), displs(nprocs))
        wrk_sizes = ntruepoints / nprocs
        displs = 0
        do n = 1, mod(ntruepoints, nprocs)
            wrk_sizes(n) = wrk_sizes(n) +1
        enddo
        do n = 2, nprocs
            displs(n) = displs(n-1) + wrk_sizes(n-1)
        enddo
    else
        ! allocate dummy sizes on slaves
        allocate(wrk_sizes(1), displs(1))
    endif
#if defined(VAR_MPI)
    call mpi_scatter(wrk_sizes, 1, MPI_INTEGER, wrk_size, 1, MPI_INTEGER, &
                   & master, MPI_COMM_WORLD, ierr)
    if (myid.ne.master) then
        allocate(wrk(3,wrk_size))
    endif
    call mpi_scatterv(wrk(1,:), wrk_sizes, displs, MPI_DOUBLE_PRECISION, &
                    & wrk(1,:), wrk_size, MPI_DOUBLE_PRECISION, &
                    & master, MPI_COMM_WORLD, ierr)
    call mpi_scatterv(wrk(2,:), wrk_sizes, displs, MPI_DOUBLE_PRECISION, &
                    & wrk(2,:), wrk_size, MPI_DOUBLE_PRECISION, &
                    & master, MPI_COMM_WORLD, ierr)
    call mpi_scatterv(wrk(3,:), wrk_sizes, displs, MPI_DOUBLE_PRECISION, &
                    & wrk(3,:), wrk_size, MPI_DOUBLE_PRECISION, &
                    & master, MPI_COMM_WORLD, ierr)
#else
    wrk_size = ntruepoints
#endif

    ! work array for potentials
    allocate(Vwrk(wrk_size))
    Vwrk = zero

    ! integral memory
    n2bas = size(density)
    allocate(integrals(n2bas))
    integrals = zero

    ! we need integrals of electronic density, hence -1
    q_one = -one

    do k = 1, wrk_size
        call one_electron_integrals(q_one, wrk(:,k), integrals)
        Vwrk(k) = dot( density, integrals )
        do m = 1, nnuclei
            dr = Rm(:,m) - wrk(:,k)
            Rmk = sqrt( dot( dr, dr ) )
            Vwrk(k) = Vwrk(k) + Zm(m) / Rmk
        enddo
    enddo

#if defined(VAR_MPI)
    call mpi_gatherv(Vwrk, wrk_size, MPI_DOUBLE_PRECISION, &
                   & V, wrk_sizes, displs, MPI_DOUBLE_PRECISION, &
                   & master, MPI_COMM_WORLD, ierr)
#else
    V = Vwrk
#endif
    ! clean-up thread storage not needed anymore
    deallocate(integrals)
    deallocate(Vwrk)
    deallocate(displs)
    deallocate(wrk_sizes)
    if (myid.ne.master) then
        deallocate(Zm, Rm)
    endif

    ! At this point we now have the potential V(R) in all
    ! points R surrounding the molecule.

    ! Currently the matrix A is built and diagonalized
    ! on the master only

    ! get dimensions of matrices/vectors according to order of multipole moments
    ! we always assume charges and dimensionality is nnuclei
    if (myid.eq.master) then
        matdim = nnuclei
        if (qfit_multipole_rank >= 1) then
            matdim = matdim + 3*nnuclei
        endif
        if (qfit_multipole_rank >= 2) then
            matdim = matdim + 6*nnuclei
        endif

        ! set constraints based on options
        condim = 0
        if (qfit_constraint >= 0) condim = condim +1
        if (qfit_constraint >= 1) condim = condim +3

        allocate(A(matdim+condim,matdim+condim))
        allocate(B(matdim+condim))
        A = 0.0_dp
        B = 0.0_dp
        call geom_mat(Rm, wrk, V, qfit_multipole_rank, A(:matdim,:matdim), B(:matdim))
        call constr_mat(Rm, Zm, qfit_constraint, total_charge, total_dipole, matdim, condim, A, B)

        ! print the final matrix A and vector B
        if (qfit_debug) then
            matsiz = matdim+condim
            write(luout,*)
            write(luout,*) "Geometric Matrix (A):"
            call output(A,1,matsiz,1,matsiz,matsiz,matsiz,1,luout)
            write(luout,*)
            write(luout,*) "Potential Vector (b):"
            call output(B,1,matsiz,1,1,matsiz,1,1,luout)
        endif

        !
        ! solve the system of linear equations Ax = b using SVD
        !
        !call linear_solve_svd( A, b, charges )
        allocate(moments(matdim+condim))
        call linear_solve_simple(A, B, moments)


        allocate(mom_pot(ntruepoints))
        allocate(pot_diff(ntruepoints))
        call eval_potential(Rm, moments, qfit_multipole_rank, wrk, ntruepoints, mom_pot)
        pot_sum = 0.0d0
        pot_diff = 0.0d0
        do m=1, ntruepoints
          pot_diff(m) = V(m) - mom_pot(m)
        enddo
        rmsd_esp_fit = rmsd(pot_diff, ntruepoints)

        ! return the resulting charges ignoring any constraints
        icf = 1
        ict = nnuclei
        idf = ict+1
        idt = 4*nnuclei
        iqf = idt+1
        iqt = 10*nnuclei
        fitted_charges = moments(icf:ict)
        if (qfit_multipole_rank >= 1 ) then
            fitted_dipoles = moments(idf:idt)
            if (qfit_multipole_rank >= 2) then
                fitted_quadrupoles = moments(iqf:iqt)
            endif
        endif
        deallocate(pot_diff)
        deallocate(mom_pot)
        deallocate(moments)
        deallocate(B)
        deallocate(A)
        deallocate(V)
        deallocate(wrk)

        if (trim(qfit_mepfile) /= '' ) then
            inquire(unit=lumepinp, opened=lunit_open)
            if (lunit_open) close( lumepinp )
        endif
    endif ! if (myid.eq.master) then

end subroutine qfit_fit

end module qfit
