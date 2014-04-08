!------------------------------------------------------------------------------
!> @brief Charge fitting library
!!
!! @author Casper Steinmann
module qfit

    use qfit_precision
    use qfit_variables
    use qfit_io, only : openfile
    implicit none

    private

    public :: qfit_initialize
    public :: qfit_finalize
    public :: qfit_print_info
    public :: qfit_fit
    public :: qfit_get_results

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

!------------------------------------------------------------------------------
!> @brief print out information about the calculation
!!
!! @author Casper Steinmann
subroutine qfit_print_info

    use qfit_variables

    if (trim(qfit_mepfile) /= '') then
        write(luout, 13) trim(qfit_mepfile)
    else
        write(luout, 10) qfit_nshell, qfit_vdwscale, qfit_vdwincrement
        write(luout, 11) qfit_pointdensity
    endif
    if (qfit_only_calculate_mep) then
        write(luout,'(/10x,a)') 'Will only print MEP. No charge-fitting will be done.'
    else
        if (qfit_constraint .eq. 0) then
            write(luout,'(/10x,a)') 'WARNING: No constraints on charges imposed.'
        else
            if (iand(1,qfit_constraint).eq.1) write(luout, 14)
            if (iand(2,qfit_constraint).eq.2) write(luout, 15)
        endif
    endif
    if (qfit_verbose) write(luout,'(/10x,a)') 'Verbose mode enabled.'
    if (qfit_debug) write(luout,'(/10x,a)') 'Debug mode enabled.'

 10 format(/10x,'Generating surface using ', i2, ' layers. Each layer is'/, &
   & 10x,'scaled by', f4.1, ' plus ', f4.1,' for each successive layer.')
 11 format(/10x,'Point density is ', f4.2, ' au^-2.')
 13 format(/10x,'Will read "',a,'" for surface points.')
 14 format(/10x,'Constraining partial charges to reproduce the', &
   &       /10x,'total molecular charge.')
 15 format(/10x,'Constraining partial charges to reproduce the', &
           /10x,'total permanet molecular dipole.')

end subroutine

!------------------------------------------------------------------------------
!> @brief returns the resulting potential derived charges
!!
!! @author Casper Steinmann
!! @param[out] charges resulting potential derived charges
subroutine qfit_get_results( charges )

    real(dp), dimension(:), intent(out) :: charges

    charges = fitted_charges(1:nnuclei)

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
    real(dp), dimension(:), allocatable :: V
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: b

    ! local variables
    real(dp) :: Rnk, Rmk
    real(dp), dimension(3) ::  dr
    real(dp), dimension(1) :: q_one
    real(dp), dimension(:), allocatable :: charges
    integer :: ntotalpoints, ntruepoints
    integer :: n, m, k
    integer :: ioff, nconstraints, nnbasx
    logical :: constrain_charges, constrain_dipoles
    logical :: lunit_open
    integer :: lumepinp
    real(dp) :: factor
    character(len=2) :: option


    ! set constraints based on options
    constrain_charges = iand(1,qfit_constraint) .ne. 0
    constrain_dipoles = iand(2,qfit_constraint) .ne. 0

    ! obtain the number of constraints we need to take into account
    nconstraints = 0
    if (constrain_charges) nconstraints = nconstraints +1 ! charges
    if (constrain_dipoles) nconstraints = nconstraints +3 ! dipoles

    ! check if the user supplied a custom file to be
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
        allocate( wrk( 3, ntotalpoints) )
        do m=1, ntotalpoints
            read(lumepinp,*) option, wrk(:,m)
        enddo
        wrk = wrk * factor

    else
        ! ... or create a grid around the molecule

        ! populates the max_layer_points array for memory management
        ! this requires that options have been set.
        call connolly_grid_count

        ! calculate the grid surrounding the molecule
        ntotalpoints = sum(max_layer_points)
        allocate( wrk( 3, ntotalpoints ) )
        call connolly_grid( wrk, ntruepoints )
    endif

    ! allocate space for evaluating the potential on those points
    allocate( V( ntruepoints ) )
    V = zero

    ! allocate space for setting up the system of linear equations
    allocate( A( nnuclei+nconstraints, nnuclei+nconstraints ) )
    allocate( b( nnuclei+nconstraints ) )
    A = zero
    b = zero

    ! allocate space for the charges
    allocate( charges( nnuclei + nconstraints ) )
    charges = zero

    ! this is the test charge we use to evaluate the electrostatic potential
    q_one = one

    ! integral memory
    nnbasx = size( density )
    allocate( integrals( nnbasx ) )
    integrals = zero


    ! determine the total potential in each point on the surface
    if (qfit_only_calculate_mep) then
        write(luout,'(/,a)') "** QFITLIB MEP **"
        write(luout,'(i4)') ntruepoints
        write(luout,*) "AU"
    endif

    do k = 1, ntruepoints
        call one_electron_integrals( q_one, wrk(:,k), integrals)
        V(k) = V(k) + dot( density, integrals )
        do m = 1, nnuclei
            dr = Rm(:,m) - wrk(:,k)
            Rmk = sqrt( dot( dr, dr ) )
            V(k) = V(k) + Zm(m) / Rmk
        enddo
    enddo

    if (qfit_only_calculate_mep) then
        do k = 1, ntruepoints
            write(luout,'(4F16.9)') wrk(:,k), V(k)
        enddo
        write(luout,'(a,/)') "*****************"
    endif

    if (.not. qfit_only_calculate_mep) then
        ! populate A matrix and b vector with charge <-> surface interaction
        do m = 1, nnuclei
            do k = 1, ntruepoints
                dr = Rm(:,m) - wrk(:,k)
                Rmk = sqrt( dot( dr, dr ) )
                b(m) = b(m) + V(k) / Rmk

                do n = 1, nnuclei
                    dr = Rm(:,n) - wrk(:,k)
                    Rnk = sqrt( dot( dr, dr ) )
                    A(m,n) = A(m,n) + one / (Rmk * Rnk)
                enddo
            enddo
        enddo

        ! add constraints to A matrix and b vector
        do m = 1, nconstraints
            ioff = m+nnuclei

            ! constrain total charge to be constant
            if( m .eq. 1 .and. constrain_charges ) then
                b(ioff) = total_charge
                do n=1,nnuclei
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
            write(luout,*)
            write(luout,*) "A-matrix"
            do m = 1, nnuclei+nconstraints
                write(luout, *) (A(m,n),n=1,nnuclei+nconstraints)
            enddo

            write(luout,*)
            write(luout,*) "b-vector"
            do m = 1, nnuclei+nconstraints
                write(luout,*) b(m)
            enddo
        endif

        ! solve the system of linear equations Ax = b using SVD
        call linear_solve_svd( A, b, charges )

        ! return the resulting charges ignoring any constraints
        fitted_charges = charges(1:nnuclei)
    endif

    deallocate( integrals )
    deallocate( charges )
    deallocate( b )
    deallocate( A )
    deallocate( V )
    deallocate( wrk )

    inquire(unit=lumepinp, opened=lunit_open)
    if (lunit_open) close( lumepinp )

end subroutine qfit_fit

end module qfit
