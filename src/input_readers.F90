!------------------------------------------------------------------------------
!> @brief Integral interface for the FMO method.
!!
!! @author Casper Steinmann
module qfit_input_readers

    use qfit_precision

    implicit none

    private

#if defined (PRG_DALTON)
    public :: dalton_input
#endif

    contains

#if defined (PRG_DALTON)
!------------------------------------------------------------------------------
!> @brief reads input section from DALTON
!!
!! @author Casper Steinmann
!!
!! @param[in] word
!! @param[in] luinp
!! @param[in] lupri
subroutine dalton_input(word, luinp, lupri)

    use qfit_variables
    use qfit_io, only : change_case

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    character(len=7) :: option

    luout = lupri

    qfitrun = .true.

    ! now parse the .dal file properly
    do
        read(luinp,'(a7)') option
        call change_case(option)

        ! read in fragment charges
        if (trim(option(2:)) == 'CONSTR') then
            read(luinp,*) option
            call change_case(option)
            if (option == 'NONE' .or. &
            &   option == 'CHARGE' .or. &
            &   option == 'DIPOLE') then

                if (option == 'NONE') then
                    qfit_constraint = 0
                endif

                ! default is 1 which means charge
                if (option == 'DIPOLE') then
                    qfit_constraint = 3 ! 1 + 2
                endif

            else
                write(luout,*) 'Constraint not recognized. Please use one of:'
                write(luout,*) '    NONE, CHARGE, DIPOLE'
            endif

        ! read atom indices in fragments. dummy read the number of atoms
        else if (trim(option(2:)) == 'VDWSCL' ) then
            read(luinp,*) qfit_vdwscale

        ! disable the electrostatic potential
        else if (trim(option(2:)) == 'VDWINC') then
            read(luinp,*) qfit_vdwincrement

        ! number of shells surrounding the molecule
        else if (trim(option(2:)) == 'NSHELL') then
            read(luinp,*) qfit_nshell

        ! surface point density
        else if (trim(option(2:)) == 'PTDENS') then
            read(luinp,*) qfit_pointdensity

        ! SVD threshold
        else if (trim(option(2:)) == 'SVDEPS') then
            read(luinp,*) qfit_eps

        ! read in the mepfile
        else if (trim(option(2:)) == 'MEPFIL') then
            read(luinp,*) qfit_mepfile

        else if (trim(option(2:)) == 'MPRANK') then
            read(luinp,*) qfit_multipole_rank

        ! verbose output
        else if (trim(option(2:)) == 'VERBOS') then
            qfit_verbose = .true.

        ! debug output
        else if (trim(option(2:)) == 'DEBUG') then
            qfit_debug = .true.
            qfit_verbose = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        else
            write(luout,*) 'Unknown option:', option, ' in *QFIT'
        end if
    end do

end subroutine dalton_input

#endif

end module qfit_input_readers
