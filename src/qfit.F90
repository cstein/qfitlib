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
    public :: qfit_print_info
    public :: qfit_fit
    public :: qfit_get_results
    public :: qfit_set_transition_dipole

    contains

!------------------------------------------------------------------------------
!> @brief initialize the charge fitting library
!!
!! @author Casper Steinmann
!! @param[in] R the coordinates of all nuclei
!! @param[in] Z the nuclear charges of the nuclei
!! @param[in] Q the total charge of the molecule
!! @param[in] D the total dipole of the molecule
!! @param[in] RCM the center of mass of the molecule
subroutine qfit_initialize(R, Z, Q, D, RCM)

    use connolly

    real(dp), dimension(:,:), intent(in) :: R
    real(dp), dimension(:), intent(in) :: Z
    integer, intent(in), optional :: Q
    real(dp), dimension(3), intent(in), optional :: D
    real(dp), dimension(3), intent(in), optional :: RCM

    total_charge = 0
    total_dipole = zero
    center_of_mass = zero

    call connolly_initialize( R, Z )

    if (present(Q)) total_charge = Q
    if (present(D)) total_dipole = D
    if (present(RCM)) center_of_mass = RCM

    allocate( fitted_charges(size(Z)) )
    fitted_charges = zero

end subroutine

!------------------------------------------------------------------------------
!> @brief releases memory back to the program
!!
!! @author Casper Steinmann
subroutine qfit_finalize

    use connolly

    deallocate( fitted_charges )
    call connolly_finalize

end subroutine

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
!> @brief returns the resulting potential derived charges
!!
!! @author Casper Steinmann
!! @param[out] charges resulting potential derived charges
subroutine qfit_get_results( charges, dipoles )

    real(dp), dimension(:), intent(out) :: charges
    real(dp), dimension(:), intent(out), optional :: dipoles

    write(luout,*) "CSS: THIS WILL FAIL"
    if (size(charges) /= nnuclei) then
        write(luout,'(/A)') "ERROR: Memory allocation in input to qfit_get_results"
        write(luout,'(A,I4,A,I4)') " is wrong. Expected:", nnuclei, " got:", size(charges)
        stop
    endif
    charges = fitted_charges(1:nnuclei)
    if (qfit_multipole_rank >= 1 .and. present(dipoles)) then
        dipoles = fitted_charges(nnuclei+1:4*nnuclei)
    endif

end subroutine

!------------------------------------------------------------------------------
!> @brief Fit charges to the molecular esp on a number of points
!!
!! @author Casper Steinmann
!! @param[in] density the density of the molecule
subroutine qfit_fit(density)

    use connolly
    use qfit_integrals
    use qfit_utilities
    use linear_solver

    real(dp), dimension(:), intent(in) :: density
    !real(dp), dimension(:), intent(out) :: charges

    ! local arrays
    real(dp), dimension(:,:), allocatable :: wrk
    real(dp), dimension(:), allocatable :: integrals
    real(dp), dimension(:), allocatable :: V, vfit
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: b

    ! local variables
    real(dp) :: Rnk, Rmk
    real(dp), dimension(3) ::  dr
    real(dp), dimension(1) :: q_one
    real(dp), dimension(:), allocatable :: charges
    integer :: ntotalpoints, ntruepoints
    integer :: n, m, k
    integer :: ioff, nconstraints, n2bas, nbas
    integer :: matsiz
    logical :: constrain_charges, constrain_dipoles
    logical :: lunit_open
    integer :: lumepinp
    real(dp) :: factor
    character(len=2) :: option

    ! get dimensions of matrices/vectors according to order of multipole moments
    ! we always assume charges and dimensionality is nnuclei
    matdim = nnuclei
    if (qfit_multipole_rank >= 1) then
        matdim = matdim + 3*nnuclei
    endif
    if (qfit_multipole_rank >= 2) then
        matdim = matdim + 5*nnuclei
    endif

    ! set constraints based on options
    constrain_charges = iand(1,qfit_constraint) .ne. 0
    constrain_dipoles = iand(2,qfit_constraint) .ne. 0 .and. qfit_multipole_rank .eq. 0

    ! obtain the number of constraints we need to take into account
    nconstraints = 0
    if (constrain_charges) nconstraints = nconstraints +1 ! charges
    if (constrain_dipoles) nconstraints = nconstraints +3 ! dipoles

    ! Either we generate a grid and dump it to a file so we can read
    ! it back later when fitting for many densities or we read in
    ! a supplied file by the user

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

    ! allocate space for evaluating the potential on those points
    allocate( V( ntruepoints ) )
    V = zero

    ! allocate space for setting up the system of linear equations
    allocate( A( matdim+nconstraints, matdim+nconstraints ) )
    allocate( b( matdim+nconstraints ) )
    A = zero
    b = zero

    ! allocate space for the charges
    allocate( charges( matdim + nconstraints ) )
    charges = zero

    ! integral memory
    n2bas = size( density )
    allocate( integrals( n2bas ) )
    integrals = zero

    ! we need integrals of electronic density, hence -1
    q_one = -one

    do k = 1, ntruepoints
        call one_electron_integrals( q_one, wrk(:,k), integrals)
        V(k) = dot( density, integrals )
        do m = 1, nnuclei
            dr = Rm(:,m) - wrk(:,k)
            Rmk = sqrt( dot( dr, dr ) )
            V(k) = V(k) + Zm(m) / Rmk
        enddo
    enddo

    call a_qq(A(1:nnuclei,1:nnuclei), b(1:nnuclei), V, wrk)
    if (qfit_multipole_rank >= 1) then
        call a_dd(A(nnuclei+1:4*nnuclei, nnuclei+1:4*nnuclei), b(nnuclei+1:4*nnuclei), V, wrk)
        call a_qd(A(nnuclei+1:4*nnuclei, 1:nnuclei), V, wrk)
        call a_qd(A(1:nnuclei, nnuclei+1:4*nnuclei), V, wrk)
    endif
    ! populate A matrix and b vector with charge <-> surface interaction

    ! add constraints to A matrix and b vector
    do m = 1, nconstraints
        ioff = m+matdim

        ! constrain total charge to be constant
        if( m .eq. 1 .and. constrain_charges ) then
            b(ioff) = total_charge
            do n=1,nnuclei ! must be nnuclei here because we only constrain charges
                A(ioff,n) = one
                A(n,ioff) = one
            enddo
        endif

        ! also constrain dipole
        if ( m .gt. 1 .and. m .le. 4 .and. constrain_dipoles ) then
            b(ioff) = total_dipole(m-1)
            do n=1,nnuclei
                A(ioff,n) = Rm(m-1,n) - center_of_mass(m-1)
                A(n,ioff) = Rm(m-1,n) - center_of_mass(m-1)
            enddo
        endif
    enddo

    if (qfit_debug) then
        matsiz = matdim+nconstraints
        write(luout,*)
        nbas = int(sqrt(real(n2bas)))
        write(luout,*) "Input Density:"
        call output(density,1,nbas,1,nbas,nbas,nbas,1,luout)

        write(luout,*) "Integrals:"
        call output(integrals,1,nbas,1,nbas,nbas,nbas,1,luout)

        write(luout,*)
        write(luout,*) "Geometric Matrix (A):"
        call output(A,1,matsiz,1,matsiz,matsiz,matsiz,1,luout)

        write(luout,*)
        write(luout,*) "Potential Vector (b):"
        call output(B,1,matsiz,1,1,matsiz,1,1,luout)
    endif

    ! solve the system of linear equations Ax = b using SVD
    call linear_solve_svd( A, b, charges )

    ! return the resulting charges ignoring any constraints
    fitted_charges = charges(1:nnuclei)


    ! lets do some statistics
    allocate( vfit( ntruepoints ) )
    vfit = zero
    do k = 1, ntruepoints
        do m = 1, nnuclei
            dr = Rm(:,m) - wrk(:,k)
            Rmk = sqrt( dot( dr, dr ) )
            vfit(k) = vfit(k) + charges(m) / Rmk
        enddo
        !if (qfit_debug) then
        !    write(luout,'(A,2F16.10)') "Vqm, Vfit = ", V(k), vfit(k)
        !endif
        V(k) = V(k) - vfit(k)
        V(k) = V(k)*V(k)
    enddo
    if (qfit_verbose) then
        write(luout, '(/1x,A,F9.6)') "@ RMS of fitted ESP = ", sqrt( sum( V ) / ntruepoints )
    endif
    deallocate( vfit )

    deallocate( integrals )
    deallocate( charges )
    deallocate( b )
    deallocate( A )
    deallocate( V )
    deallocate( wrk )

    if (trim(qfit_mepfile) /= '' ) then
        inquire(unit=lumepinp, opened=lunit_open)
        if (lunit_open) close( lumepinp )
    endif

end subroutine qfit_fit

subroutine a_qq(A, b, V, Rs)

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
end subroutine a_qq

subroutine a_qd(A, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, R, R3, Rmnk, drmnk
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
    if (nn.gt.nm) then
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
        do n = 1, nn
            nidx = n
            do k = 1, ntruepoints
                dr = Rm(:,n) - Rs(:,k)
                Rnk = sqrt(dot( dr, dr ))
                drmnk = dr(i)
                do m = 1, nm
                    dr = Rm(:,m) - Rs(:,k)
                    Rmk = sqrt(dot( dr, dr ))
                    midx = m
                    if (istranspose) then
                        midx = (i-1)*nm+m
                        R3 = Rmk * Rmk * Rmk
                        drmnk = dr(i)
                        Rmnk = Rnk
                    else
                        nidx = (i-1)*nn+n
                        R3 = Rnk * Rnk * Rnk
                        Rmnk = Rmk
                    end if
                    A(midx, nidx) = A(midx, nidx) + drmnk/(R3*Rmnk)
                    ! A(m,n) = A(m,n) + one / (Rmk * Rnk)
                    !write(*,'(i4,A,2i3,A,2I3, 3F9.4)') i, " A(", m, n, ")", midx, nidx, dr(i), R3, Rnk
                enddo
            enddo
        enddo

    enddo
end subroutine a_qd

subroutine a_dd(A, b, V, Rs)

    use qfit_utilities, only : dot

    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:), intent(inout) :: b
    real(dp), dimension(:), intent(in) ::  V
    real(dp), dimension(:,:), intent(in) :: Rs ! surface coordinates

    real(dp) :: Rmk, Rnk, Rmk3, Rnk3
    real(dp), dimension(3) :: drm, drn
    integer :: i, j, k, m, n
    integer :: nm, nn
    integer :: ntruepoints

    nm = size(A,1) / 3
    nn = size(A,2) / 3
    ntruepoints = size(V)

    A = 0.0_dp
    b = 0.0_dp

    ! i,j are cartesian components (x, y and z)
    do i = 1, 3
        do m = 1, nm
            do k = 1, ntruepoints
                drm = Rm(:,m) - Rs(:,k)
                Rmk = sqrt( dot( drm, drm ) )
                Rmk3 = Rmk*Rmk*Rmk
                b( (i-1)*nm + m ) = b( (i-1)*nm + m ) + V(k) * drm(i) / Rmk3
                !write(*,'(A,2I4, 6F16.6)') "B", m, k, Rm(:,m), wrk(:,k)

                do j = 1, 3
                    do n = 1, nm
                        drn = Rm(:,n) - Rs(:,k)
                        Rnk = sqrt( dot( drn, drn ) )
                        Rnk3 = Rnk*Rnk*Rnk
                        A( (i-1)*nm +m, (j-1)*nm + n) = A( (i-1)*nm +m, (j-1)*nm + n) + drm(i)*drn(j) / (Rmk3 * Rnk3)
                        !write(*,'(A,I4,5F16.10)') "    C", n, Rmk, Rnk, Rm(:,n)
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine a_dd

end module qfit
