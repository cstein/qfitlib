module potential

use tensor

implicit none

public :: eval_potential, rmsd

contains

subroutine eval_potential(r_mol, moments, fit_order, r_mep, nmep, mom_pot)
implicit none

double precision, dimension(:,:), intent(in) :: r_mol
double precision, dimension(:), intent(in) :: moments
integer, intent(in) :: fit_order
double precision, dimension(:,:), intent(in) :: r_mep
integer, intent(in) :: nmep
double precision, dimension(:), intent(out) :: mom_pot

integer :: nat, iat, imep, vec_idx, quad_idx
double precision :: q
double precision, dimension(3) :: r_a, r_b, vec
double precision, dimension(6) :: quad

mom_pot = 0.0d0
nat = size(r_mol, 2)
do iat = 1, nat
  ! extract moments
  q = moments(iat)
  r_a = r_mol(:,iat)

  vec_idx = nat + 3*(iat-1)+1
  if (fit_order >= 1) then
    vec = moments(vec_idx:vec_idx+2)
  end if

  quad_idx = nat + 3 * nat + 6*(iat-1) + 1
  if (fit_order >= 2) then
    quad = moments(quad_idx:quad_idx+5)
  end if

  do imep = 1, nmep
    r_b = r_mep(:,imep)

    mom_pot(imep) = mom_pot(imep) + q * t0(r_b, r_a)
    if (fit_order >= 1) then
      mom_pot(imep) = mom_pot(imep) + dot_product(vec, t1(r_b, r_a))
    endif

    if (fit_order >= 2) then
      mom_pot(imep) = mom_pot(imep) + 0.5 * dot_product(quad, t2(r_b, r_a))
    endif
  enddo
enddo

end subroutine eval_potential

function rmsd(array, n) result(var)

double precision, dimension(:), intent(in) :: array
integer, intent(in) :: n
double precision :: var

double precision, dimension(n) :: dev
double precision :: mean

mean = sum(array) / n

!allocate(dev(n))
dev = array - mean
var = sqrt(sum(dev*dev) / n)
!deallocate(dev)

end function

end module potential
