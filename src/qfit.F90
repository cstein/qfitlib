!------------------------------------------------------------------------------
!> @brief Charge fitting library
!!
!! @author Casper Steinmann
module qfit

    use qfit_precision
    use qfit_variables
    implicit none

    private

    public :: fit_density
    public :: qfit_initialize
    public :: qfit_finalize
    public :: qfit_print_info

    contains

!------------------------------------------------------------------------------
!> @brief initialize the charge fitting library
!!
!! @author Casper Steinmann
!! @param[in] R the coordinates of all nuclei
!! @param[in] Z the nuclear charges of the nuclei
!! @param[in] Q the total charge of the molecule
subroutine qfit_initialize(R, Z, Q)

    use connolly

    real(dp), dimension(:,:), intent(in) :: R
    real(dp), dimension(:), intent(in) :: Z
    integer, intent(in) :: Q

    call connolly_initialize( R, Z )
    total_charge = Q

end subroutine

!------------------------------------------------------------------------------
!> @brief releases memory back to the program
!!
!! @author Casper Steinmann
subroutine qfit_finalize

    use connolly

    call connolly_finalize

end subroutine

!------------------------------------------------------------------------------
!> @brief print out information about the calculation
!!
!! @author Casper Steinmann
subroutine qfit_print_info

    use qfit_variables

    write(luout, 10) qfit_nshell
    write(luout, 11) qfit_vdwscale, qfit_vdwincrement
    write(luout, 12) qfit_pointdensity
    if (qfit_verbose) write(luout,'(/10x,a)') 'Verbose mode enabled.'
    if (qfit_debug) write(luout,'(/10x,a)') 'Debug mode enabled.'

 10 format(/10x,'Uses', i2, ' shells of points around the molecule.')
 11 format(/10x,'Each layer is scaled by', f4.1, ' plus an additional', f4.1, &
   &       /10x,'for each successive layer.')
 12 format(/10x,'Point density is ', f4.2, ' au^-2.')

end subroutine

!------------------------------------------------------------------------------
!> @brief Fit charges to the molecular esp on a number of points
!!
!! @author Casper Steinmann
!! @param[in] density the density of the molecule
!! @param[out] charges the fitted charges
subroutine fit_density(density, charges)

    use connolly
    use qfit_integrals
    use qfit_utilities
    use linear_solver

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(out) :: charges

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
    integer :: ntotalpoints, ntruepoints
    integer :: n, m, k
    integer :: ioff, nconstraints, nnbasx

    ! populates the max_layer_points array for memory management
    ! this requires that options have been set.
    call connolly_grid_count

    ! calculate the grid surrounding the molecule
    ntotalpoints = sum(max_layer_points)
    allocate( wrk( 3, ntotalpoints ) )
    call connolly_grid( wrk, ntruepoints )

    ! allocate space for evaluating the potential on those points
    allocate( V( ntruepoints ) )
    V = zero

    ! obtain the number of constraints we need to take into account
    nconstraints = 1 ! currently hardcoded for charges

    ! allocate space for setting up the linear system
    allocate( A( nnuclei+nconstraints, nnuclei+nconstraints ) )
    allocate( b( nnuclei+nconstraints ) )
    A = zero
    b = zero

    ! allocate space for the charges
    allocate( fitted_charges( nnuclei + nconstraints ) )
    fitted_charges = zero
    q_one = one

    ! integral memory
    nnbasx = size( density )
    allocate( integrals( nnbasx ) )
    integrals = zero


    ! first determine the total potential
    ! todo: separate into own function
    do k = 1, ntruepoints
        call one_electron_integrals( q_one, wrk(:,k), integrals)
        V(k) = V(k) + dot( density, integrals )
        do m = 1, nnuclei
            dr = Rm(:,m) - wrk(:,k)
            Rmk = sqrt( dot( dr, dr ) )
            V(k) = V(k) + Zm(m) / Rmk
        enddo
    enddo

    ! populate A matrix and b vector
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

    ! add constraints.
    ! todo: this needs to be done way more cleverly
    do m = 1, nconstraints
        ioff = m+nnuclei
        ! constrain total charge to be constant
        if( m .eq. 1 ) then
            b(ioff) = total_charge
            do n=1,nnuclei
                A(ioff,n) = one
                A(n,ioff) = one
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

    ! solve Ax = b
    call linear_solve( A, b, fitted_charges )

    !write(luout,*)
    !write(luout,*) 'Density fitted charges'
    !do m = 1, nnuclei
        !write(luout,'(i4,f9.3)') m, fitted_charges(m)
    !enddo
    !write(luout,'(a,f9.3)') ' sum =', sum( fitted_charges(1:nnuclei) )

    charges = fitted_charges(1:nnuclei)

    deallocate( integrals )
    deallocate( V )
    deallocate( A )
    deallocate( b )
    deallocate( fitted_charges )
    deallocate( wrk )

end subroutine fit_density

end module qfit
