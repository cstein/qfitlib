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
        write(luout, 16) 'charges, dpioles and quadrupoles'
    elseif (qfit_multipole_rank .eq. 1) then
        write(luout, 16) 'charges and dipoles'
    elseif (qfit_multipole_rank .eq. 0) then
        write(luout, 16) 'charges'
    else
        write(luout, '(/10x,a)') 'WARNING: multipole order not understood.'
    endif

    if (trim(qfit_mepfile) /= '') then
        write(luout, 13) trim(qfit_mepfile)
    else
        write(luout, 10) qfit_nshell, qfit_vdwscale, qfit_vdwincrement
        write(luout, 11) qfit_pointdensity
    endif
    if (qfit_constraint .eq. 0) then
        write(luout,'(/10x,a)') 'WARNING: No constraints on charges imposed.'
    else
        if (iand(1,qfit_constraint).eq.1) write(luout, 14)
        if (iand(2,qfit_constraint).eq.2 .and. qfit_multipole_rank .eq. 0) write(luout, 15)
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
subroutine qfit_get_results( charges, dipoles, quadrupoles )

    real(dp), dimension(:), intent(out) :: charges
    real(dp), dimension(:), intent(out), optional :: dipoles
    real(dp), dimension(:), intent(out), optional :: quadrupoles

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
    real(dp), dimension(:), allocatable :: b

    ! local variables
    real(dp) :: Rmk
    real(dp), dimension(3) ::  dr, drhat, mu
    real(dp), dimension(6) ::  oo
    real(dp), dimension(1) :: q_one
    real(dp), dimension(:), allocatable :: moments
    integer :: ntotalpoints, ntruepoints
    integer :: n, m, k
    integer :: ioff, n2bas, nbas
    integer :: matsiz, matdim, condim
    logical :: constrain_charges, constrain_dipoles
    logical :: lunit_open
    integer :: lumepinp
    integer :: icf, ict, idf, idt, iqf, iqt
    real(dp) :: factor
    character(len=2) :: option

    ! mpi specific variables
    integer, save :: master, myid, nprocs, ierr
    integer :: wrk_size
    integer, allocatable, dimension(:) :: wrk_sizes, displs

#if !defined(VAR_MPI)
    myid = 0
    master = myid
    nprocs = 1
#else
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
    deallocate(integrals, Vwrk, displs, wrk_sizes)
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
        constrain_charges = iand(1,qfit_constraint) .ne. 0
        constrain_dipoles = iand(2,qfit_constraint) .ne. 0 ! .and. qfit_multipole_rank .eq. 0

        condim = 0
        if (constrain_charges) condim = condim +1
        if (constrain_dipoles) condim = condim +3

        allocate(A(matdim+condim,matdim+condim))
        allocate(b(matdim+condim))
        allocate(moments(matdim+condim))
        A = 0.0_dp
        b = 0.0_dp
        moments = 0.0_dp

        ! build sub-blocks of matrix A and rhs side b
        icf = 1
        ict = nnuclei
        idf = ict+1
        idt = 4*nnuclei
        iqf = idt+1
        iqt = 10*nnuclei
        call a_cc(A(icf:ict,icf:ict), b(icf:ict), V, wrk)
        if (qfit_multipole_rank >= 1) then
            call a_cd(A(idf:idt, icf:ict), V, wrk)
            call a_cd(A(icf:ict, idf:idt), V, wrk)
            call a_dd(A(idf:idt, idf:idt), b(idf:idt), V, wrk)

            if (qfit_multipole_rank >= 2) then
                !call a_cq(A(iqf:iqt, icf:ict), V, wrk)
                !call a_cq(A(icf:ict, iqf:iqt), V, wrk)

                !call a_dq(A(idf:idt, iqf:iqt), V, wrk)
                !call a_dq(A(iqf:iqt, idf:idt), V, wrk)

                !call a_qq(A(iqf:iqt, iqf:iqt), b(iqf:iqt), V, wrk)
            endif
        endif

        ! add constraints to charges and dipoles to A matrix and b vector
        ! we count from 0 (for charges) so that 1 to 3 are used for x, y
        ! and z for dipole constraints
        do m = 0, condim-1
            ioff = matdim+m+1

            ! constrain total charge to be constant
            if( m .eq. 0 .and. constrain_charges ) then
                b(ioff) = total_charge
                do n=1, nnuclei ! must be nnuclei here because we only constrain charges
                    A(ioff,n) = one
                    A(n,ioff) = one
                enddo
            endif

            ! also constrain dipole
            if ( m .gt. 0 .and. m .le. 3 .and. constrain_dipoles ) then
                b(ioff) = total_dipole(m)
                do n=1, nnuclei
                    A(ioff,n) = Rm(m, n) - center_of_mass(m)
                    A(n,ioff) = Rm(m, n) - center_of_mass(m)
                enddo
            endif
        enddo

        ! print the final matrix A and vector b
        if (qfit_debug) then
            matsiz = matdim+condim
            write(luout,*)
            write(luout,*) "Geometric Matrix (A):"
            call output(A,1,matsiz,1,matsiz,matsiz,matsiz,1,luout)
            write(luout,*)
            write(luout,*) "Potential Vector (b):"
            call output(b,1,matsiz,1,1,matsiz,1,1,luout)
        endif

        !
        ! solve the system of linear equations Ax = b using SVD
        !
        !call linear_solve_svd( A, b, charges )
        call linear_solve_simple(A, b, moments)

        ! return the resulting charges ignoring any constraints
        fitted_charges = moments(icf:ict)
        if (qfit_multipole_rank >= 1 ) then
            fitted_dipoles = moments(idf:idt)
            if (qfit_multipole_rank >= 2) then
                fitted_quadrupoles = moments(iqf:iqt)
            endif
        endif

        ! Perform statistics on the ESP computed from the electron density
        ! to the ESP we can compute using our fitted electrostatic moments
        allocate(vcharges(ntruepoints))
        allocate(vdipoles(ntruepoints))
        allocate(vquadrupoles(ntruepoints))
        vcharges = zero
        vdipoles = zero
        vquadrupoles = zero

        do k = 1, ntruepoints
            do m = 1, nnuclei
                dr = Rm(:,m) - wrk(:,k)
                Rmk = sqrt(dot(dr, dr))
                vcharges(k) = vcharges(k) + fitted_charges(m) / Rmk
                if (qfit_multipole_rank >= 1) then
                    drhat = dr / Rmk
                    mu(1) = fitted_dipoles(m)
                    mu(2) = fitted_dipoles(m+1*nnuclei)
                    mu(3) = fitted_dipoles(m+2*nnuclei)

                    vdipoles(k) = vdipoles(k) + dot(drhat, mu) / (Rmk*Rmk)
                    if (qfit_multipole_rank >= 2) then
                        ! XX
                        ! XY
                        ! XZ
                        ! YY
                        ! YZ
                        ! (ZZ = -(XX +YY))
                        oo(1) = fitted_quadrupoles(m)
                        oo(2) = fitted_quadrupoles(m+1*nnuclei)
                        oo(3) = fitted_quadrupoles(m+2*nnuclei)
                        oo(4) = fitted_quadrupoles(m+3*nnuclei)
                        oo(5) = fitted_quadrupoles(m+4*nnuclei)
                        oo(6) = fitted_quadrupoles(m+5*nnuclei)

                        mu(1) = oo(1)*drhat(1) + oo(2)*drhat(2) + oo(3)*drhat(3)
                        mu(2) = oo(2)*drhat(1) + oo(4)*drhat(2) + oo(5)*drhat(3)
                        mu(3) = oo(3)*drhat(1) + oo(5)*drhat(2) + oo(6)*drhat(3)

                        vquadrupoles(k) = vquadrupoles(k) + dot(drhat, mu) / (Rmk*Rmk*Rmk)
                    endif
                endif
            enddo

            if (qfit_debug) then
                write(luout,'(A,5F16.10)') "Vqm, Vq, Vd, Vtot = ", V(k), &
   &            vcharges(k),  vdipoles(k),  vquadrupoles(k), &
   &            vcharges(k) + vdipoles(k) + vquadrupoles(k)
            endif
            V(k) = V(k) - (vcharges(k) + vdipoles(k) + vquadrupoles(k))
            V(k) = V(k)*V(k)
        enddo
        if (qfit_verbose) then
            write(luout, '(/A,F9.6)') "@ RMSE of fitted ESP = ", sqrt(sum(V)/ntruepoints)
        endif
        deallocate(vquadrupoles)
        deallocate(vdipoles)
        deallocate(vcharges)

        deallocate(moments)
        deallocate(b)
        deallocate(A)
        deallocate(V)
        deallocate(wrk)

        if (trim(qfit_mepfile) /= '' ) then
            inquire(unit=lumepinp, opened=lunit_open)
            if (lunit_open) close( lumepinp )
        endif
    endif

end subroutine qfit_fit


!> @brief fills matrix A with charge-charge contributions
!!
!! @author Casper Steinmann
!!
!! @param A subblock of matrix A containing the charge-charge interaction
!! @param b subblock of the rhs vector B containing the QM potential
!! @param V qm potential
!! @param Rs surface points
subroutine a_cc(A, b, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(inout) :: b
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk
    real(dp), dimension(3) :: dr
    integer :: k, m, n
    integer :: nm, nn
    integer :: ntruepoints

    nm = size(A,1)
    nn = size(A,2)
    ntruepoints = size(V)

    A = 0.0_dp
    b = 0.0_dp

    do m = 1, nm
        do k = 1, ntruepoints
            dr = Rm(:,m) - Rs(:,k)
            Rmk = sqrt( dot( dr, dr ) )
            b(m) = b(m) + V(k) / Rmk
            !write(*,'(A,2I4, 6F16.6)') "B", m, k, Rm(:,m), wrk(:,k)

            do n = 1, nm
                dr = Rm(:,n) - Rs(:,k)
                Rnk = sqrt( dot( dr, dr ) )
                A(m,n) = A(m,n) + one / (Rmk * Rnk)
                !write(*,'(A,I4,5F16.10)') "    C", n, Rmk, Rnk, Rm(:,n)
            enddo
        enddo
    enddo
end subroutine a_cc

subroutine a_cd(A, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, R3, Rmnk, drmnk
    real(dp), dimension(3) :: dr
    integer :: i, k, m, n
    integer :: midx, nidx
    integer :: nm, nn, nmt, nnt
    integer :: ntruepoints
    logical :: istranspose

    ! figure out which dimension is the largest - it is the one with the x, y and z-components
    ! of the dipole
    nmt = size(A,1)
    nnt = size(A,2)
    ntruepoints = size(V)

    istranspose = .true.

    !write(*,*) "nm,nn:", nmt, nnt
    ! if nm > nn, then we have the regular matrix
    ! otherwise, it is the transpose
    nm = nmt / 3
    nn = nnt
    if (nnt.gt.nmt) then
        istranspose = .false.
        nm = nmt
        nn = nnt / 3
    endif

    !write(*,*) "is tranpose?", istranspose
    ! if transpose
    !   m is dipole
    !   n is charge
    ! otherise
    !   m is charge
    !   n is dipole

    A = 0.0_dp

    ! i is a cartesian component (x, y or z)
    do i = 1, 3

        ! loop over charges
        do n = 1, nn
            nidx = n

            do k = 1, ntruepoints
                dr = Rm(:,n) - Rs(:,k)
                Rnk = sqrt(dot( dr, dr ))
                drmnk = dr(i)

                ! loop over dipoles
                do m = 1, nm
                    dr = Rm(:,m) - Rs(:,k)
                    Rmk = sqrt(dot( dr, dr ))
                    midx = m
                    if (istranspose) then
                        midx = (i-1)*nm+m
                        R3 = Rmk**3
                        drmnk = dr(i)
                        Rmnk = Rnk
                    else
                        nidx = (i-1)*nn+n
                        R3 = Rnk**3
                        Rmnk = Rmk
                    end if
                    A(midx, nidx) = A(midx, nidx) + drmnk/(R3*Rmnk)
                    ! A(m,n) = A(m,n) + one / (Rmk * Rnk)
                    !write(*,'(i4,A,2i3,A,2I3, 3F16.5)') i, " A(", m, n, ")", midx, nidx, drmnk, R3, Rmnk
                enddo
            enddo
        enddo

    enddo
end subroutine a_cd


subroutine a_dd(A, b, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(inout) :: b
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, Rmk3, Rnk3
    real(dp), dimension(3) :: drmk, drnk
    integer :: i, j, k, m, n
    integer :: nm, nn
    integer :: midx, nidx
    integer :: ntruepoints

    nm = size(A,1) / 3
    nn = size(A,2) / 3
    ntruepoints = size(V)

    A = 0.0_dp
    b = 0.0_dp

    ! i, j are cartesian components (x, y and z) of the dipoles
    do i = 1, 3

        ! loop over dipoles
        do n = 1, nn
            nidx = (i-1)*nn + n
            do k = 1, ntruepoints
                drnk = Rm(:,n) - Rs(:,k)
                Rnk = sqrt( dot( drnk, drnk ) )
                Rnk3 = Rnk**3
                b( nidx ) = b( nidx ) + V(k) * drnk(i) / Rnk3
                !write(*,'(A,2I4, 6F16.6)') "B", m, k, Rm(:,m), wrk(:,k)

                do j = 1, 3
                    do m = 1, nm
                        midx = (j-1)*nm + m
                        drmk = Rm(:,m) - Rs(:,k)
                        Rmk = sqrt( dot( drmk, drmk ) )
                        Rmk3 = Rmk**3
                        A(midx, nidx) = A(midx, nidx) + drmk(i)*drnk(j) / (Rmk3 * Rnk3)
                        !write(*,'(2i4,A,2i3,A,2I3, 3F16.5)') j, i, " A(", m, n, ")", midx, nidx, Rmk3, Rnk3
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine a_dd

! Adds the charge <-> quadrupole interaction
! @note does not work
subroutine a_cq(A, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, R5, Rmnk
    real(dp), dimension(3) :: drnk, drmk, R
    integer :: i, j, k, m, n
    integer :: midx, nidx
    integer :: nm, nn, nmt, nnt
    integer :: ntruepoints
    logical :: istranspose

    A = 0.0_dp

    ! figure out which dimension is the largest
    nmt = size(A,1)
    nnt = size(A,2)
    ntruepoints = size(V)

    istranspose = .true.
    nm = nmt / 5
    nn = nnt
    if (nnt.gt.nmt) then
        istranspose = .false.
        nm = nmt
        nn = nnt / 5
    endif
    !write(*,*) "is tranpose?", istranspose
    ! if transpose
    !   m is quadrupole
    !   n is charge
    ! otherise
    !   m is charge
    !   n is quadrupole

    ! loop over first multipole moment
    do n = 1, nn
        nidx = n

        ! loop over surface
        do k = 1, ntruepoints
            drnk = Rm(:,n) - Rs(:,k)
            Rnk = sqrt(dot( drnk, drnk ))

            ! loop over second multipole moment
            do m = 1, nm
                midx = m
                drmk = Rm(:,m) - Rs(:,k)
                Rmk = sqrt(dot( drmk, drmk ))

                ! i is a cartesian component (x, y or z)
                do i = 1, 2 ! only over x and y
                    do j = i, 3 ! x, y and z

                        if (istranspose) then
                            midx = (i-1)*nm+m + (j-1)*nm + (i-1)*nm
                            R5 = Rmk**5
                            R = drmk
                            Rmnk = Rnk
                        else
                            nidx = (i-1)*nn+n + (j-1)*nn + (i-1)*nn
                            R5 = Rnk**5
                            R = drnk
                            Rmnk = Rmk
                        end if

                        A(midx, nidx) = A(midx, nidx) + f(R,i,j)/(R5*Rmnk)
                        ! A(m,n) = A(m,n) + one / (Rmk * Rnk)
                        !write(*,'(2i4,A,2i3,A,2I3, 3F16.5,F20.10)') i, j, " A(", m, n, ")", midx, nidx, f(R, i, j), & 
    !& R5, Rmnk, f(R,i,j)/(R5*Rmnk)
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine a_cq

! Adds the dipole <-> quadrupole interaction
! @note does not work
subroutine a_dq(A, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, R5, R3, drmnk
    real(dp), dimension(3) :: drmk, drnk, R
    integer :: i, j, k, m, n
    integer :: alpha
    integer :: midx, nidx
    integer :: nm, nn, nmt, nnt
    integer :: ntruepoints
    logical :: istranspose

    ! figure out which dimension is the largest
    nmt = size(A,1)
    nnt = size(A,2)
    ntruepoints = size(V)

    A = 0.0_dp

    istranspose = .true.
    nm = nmt / 5
    nn = nnt / 3
    if (nnt.gt.nmt) then
        istranspose = .false.
        nm = nmt / 3
        nn = nnt / 5
    endif
    !write(*,*) "nm,nn:", nmt, nnt
    !write(*,*) "is tranpose?", istranspose
    ! if transpose
    !   m is quadrupole
    !   n is dipole
    ! otherise
    !   m is dipole
    !   n is quadrupole


    ! alpha is a Cartesian component over the dipole
    do alpha = 1, 3

        ! loop over mu_n dipoles
        do n = 1, nn
            nidx = n + (alpha-1)*nn

            ! loop over surface
            do k = 1, ntruepoints
                drnk = Rm(:,n) - Rs(:,k)
                Rnk = sqrt(dot( drnk, drnk ))

                ! loop over O_m quadrupoles
                do m = 1, nm
                    drmk = Rm(:,m) - Rs(:,k)
                    Rmk = sqrt(dot( drmk, drmk ))
                    midx = m + (alpha-1)*nm

                    !write(*,*) "RNK", drnk, Rnk
                    !write(*,*) "RMK", drmk, Rmk

                    ! i is a cartesian component (x, y or z)
                    ! the i, j loops are over Cartesian components
                    ! of the quadrupole moments
                    do i = 1, 2 ! only over x and y
                        do j = i, 3 ! x, y and z

                            if (istranspose) then
                                midx = (i-1)*nm+m + (j-1)*nm + (i-1)*nm
                                R5 = Rmk**5
                                R3 = Rnk**3
                                R = drmk
                                drmnk = drnk(alpha)
                            else
                                nidx = (i-1)*nn+n + (j-1)*nn + (i-1)*nn
                                R5 = Rnk**5
                                R3 = Rmk**3
                                R = drnk
                                drmnk = drmk(alpha)
                            end if

                            A(midx, nidx) = A(midx, nidx) + drmnk*f(R,i,j)/(R5*R3)
                            ! A(m,n) = A(m,n) + one / (Rmk * Rnk)
                            !write(*,'(3i4,A,2i3,A,2I3,4F9.2,F20.10)') alpha, i, j, " A(", m, n, ")", midx, nidx, & 
                            !    & R5, R3, drmnk, f(R, i, j),  drmnk*f(R,i,j)/(R5*R3)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine a_dq

! Adds the quadrupole <-> quadrupole interaction
! @note does not work
subroutine a_qq(A, b, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(inout) :: b
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, Rmk5, Rnk5
    real(dp), dimension(3) :: drmk, drnk
    integer :: i, j, k, m, n
    integer :: alpha, beta
    integer :: midx, nidx
    integer :: nm, nn, nmt, nnt
    integer :: ntruepoints

    ! figure out which dimension is the largest
    nmt = size(A,1)
    nnt = size(A,2)
    ntruepoints = size(V)

    A = 0.0_dp
    b = 0.0_dp

    nm = nmt / 5
    nn = nnt / 5
    !write(*,*) "nm,nn:", nmt, nnt


    ! i is a cartesian component (x, y or z)
    ! the i, j loops are over Cartesian components
    ! of the quadrupole moments

    ! alpha and beta are Cartesian components over the quadrupole
    do alpha = 1, 2
        do beta = alpha, 3

            ! loop over O_n quadrupoles
            do n = 1, nn
                !nidx = n + (alpha-1)*3
                nidx = (alpha-1)*nn+n + (beta-1)*nn + (alpha-1)*nn

                ! loop over surface
                do k = 1, ntruepoints
                    drnk = Rm(:,n) - Rs(:,k)
                    Rnk = sqrt(dot( drnk, drnk ))
                    Rnk5 = Rnk**5
                    b( nidx ) = b( nidx ) + V(k) * drnk(alpha)*drnk(beta) / Rnk5

                    ! i and j are Cartesian components over the quadrupole
                    do i = 1, 2 ! only over x and y
                        do j = i, 3 ! x, y and z

                            ! loop over O_m quadrupoles
                            do m = 1, nm
                                midx = (i-1)*nm+m + (j-1)*nm + (i-1)*nm

                                drmk = Rm(:,m) - Rs(:,k)
                                Rmk = sqrt(dot( drmk, drmk ))
                                Rmk5 = Rmk**5
                                !midx = m + (alpha-1)*3

                                !if (istranspose) then
                                !    R5 = Rmk**5
                                !    R3 = Rnk**3
                                !    R = dr
                                !else
                                !    nidx = (i-1)*nn+n + (j-1)*3 + (i-1)*3
                                !    R5 = Rnk**5
                                !    R3 = Rmk**3
                                !    drmnk = dr(alpha)
                                !end if

                                A(midx, nidx) = A(midx, nidx) + f(drmk, i, j)*f(drnk, alpha, beta)/(Rnk5 * Rmk5)
                                ! A(m,n) = A(m,n) + one / (Rmk * Rnk)
                                !write(*,'(4i4,A,2i3,A,2I3)') i, j, alpha, beta, " A(", m, n, ")", midx, nidx
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine a_qq

!> @brief Switching function for traceless quadrupole moment calculations.
!!
!! @author Casper Steinmann
!!
!! @param[in] r distance vector
!! @param[in] alpha cartesian component index (1, 2 or 3)
!! @param[in] beta cartesian component index (1, 2 or 3)
!!
!! @see S. Jakobsen, F. Jensen, J. Chem. Theory Comput. (2014), 10, 5493-5504
!!      DOI: 10.1021/ct500803r
real(dp) function f(r, alpha, beta)

    real(dp), dimension(3), intent(in) :: r
    integer, intent(in) :: alpha
    integer, intent(in) :: beta

    ! let's do some sanity checks
    if (alpha < 1 .or. alpha > 3) then
        stop 'Invalid index specified for alpha.'
    endif

    if (beta < 1 .or. beta  > 3) then
        stop 'Invalid index specified for beta.'
    endif

    f = 2.0_dp * r(alpha) * r(beta)
    if (alpha .eq. beta) then
        f = r(alpha)*r(beta) - r(3)*r(3)
    endif

end function f

end module qfit
