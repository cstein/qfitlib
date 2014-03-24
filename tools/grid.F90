program test

    use connolly
    use qfit_precision
    use qfit_variables

    implicit none

    real(dp), dimension(:,:), allocatable :: wrk
    integer :: nat
    character(len=80) :: option
    character(len=2) :: atom
    real(dp) :: factor = 1.0_dp
    real(dp), dimension(:,:), allocatable :: R
    real(dp), dimension(:), allocatable :: Z
    integer :: ntotalpoints, ntruepoints, k
    integer :: ifile

qfit_nshell = 1
qfit_pointdensity = 0.28_dp
call openfile("input.xyz", ifile, 'old', 'formatted')
rewind(ifile)
read(ifile,*) nat
read(ifile,*) option
if(trim(option) == 'AA') then
    factor = aa2au
endif
allocate(R(3,nat))
allocate(Z(nat))
do k=1,nat
    read(ifile,*) atom, R(:,k)
    if (trim(atom) == 'H') Z(k) = 1.0_dp
    if (trim(atom) == 'C') Z(k) = 6.0_dp
    if (trim(atom) == 'N') Z(k) = 7.0_dp
    if (trim(atom) == 'O') Z(k) = 8.0_dp
    if (trim(atom) == 'F') Z(k) = 9.0_dp
    if (trim(atom) == 'P') Z(k) = 15.0_dp
    if (trim(atom) == 'S') Z(k) = 16.0_dp
enddo
close(ifile)

R = R * factor

call connolly_initialize( R, Z )
call connolly_grid_count

ntotalpoints = sum(max_layer_points)
allocate( wrk( 3, ntotalpoints ) )
call connolly_grid( wrk, ntruepoints )

write(*,'(i4)') ntruepoints
write(*,*)

do k = 1, ntruepoints
    write(*,'(a,f20.9,2f16.9)') 'X', wrk(:,k)
enddo

call connolly_finalize

deallocate( wrk )
deallocate( R )
deallocate( Z )


end program test

subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer, intent(out) :: lunit
    integer :: i
    logical :: lexist, lopen

    if (stat == 'old') then
        inquire(file=filename, exist=lexist)

        if (.not. lexist) then
            print *, filename, ' not found!'
            stop
        end if
    end if

    do i = 21, 99
        inquire(unit=i, opened=lopen)
        if (lopen) then
            cycle
        else
            lunit = i
            open(unit=lunit, file=filename, status=stat, form=frmt)
            exit
        end if
    end do

    return

end subroutine openfile
