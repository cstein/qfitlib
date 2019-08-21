module auxmat

use tensor

implicit none

private

public :: geom_mat
public :: constr_mat

contains

subroutine auxvec(r_mol, r_mep, v_mep, order, A, B)
  double precision, intent(in), dimension(:,:) :: r_mol
  double precision, intent(in), dimension(:,:) :: r_mep
  double precision, intent(in), dimension(:) :: v_mep
  integer, intent(in):: order

  double precision, intent(out), dimension(:,:) :: A
  double precision, intent(out), dimension(:) :: B

  double precision :: t0_val
  double precision, dimension(3) :: t1_val
  double precision, dimension(6) :: t2_val

  integer :: nat, nmep, iat, imep, offset

  nat = size(r_mol, 2)
  nmep = size(v_mep)

  A = 0.0d0
  B = 0.0d0

  do imep = 1, nmep
    do iat = 1, nat
      offset = iat
      t0_val = t0(r_mol(:,iat), r_mep(:,imep))
      A(offset:offset, imep) = A(offset:offset, imep) + t0_val
      B(offset:offset) = B(offset:offset) + v_mep(imep) * t0_val

      if (order >= 1) then
          offset = nat + 3*(iat-1) +1
          t1_val = t1(r_mol(:,iat), r_mep(:,imep))
          A(offset:offset+2, imep) = A(offset:offset+2, imep) - 1.0d0 * t1_val(:)
          B(offset:offset+2) = B(offset:offset+2) - 1.0d0 * v_mep(imep) * t1_val
      endif

      if (order >= 2) then
          offset = nat + 3*nat + 6*(iat -1) +1
          t2_val = t2(r_mol(:,iat), r_mep(:,imep))
          A(offset:offset+5, imep) = A(offset:offset+5, imep) + 0.5d0 * t2_val(:)
          B(offset:offset+5) = B(offset:offset+5) + 0.5d0 * v_mep(imep) * t2_val
      endif

    enddo
  enddo

end subroutine auxvec



subroutine get_aux_dim(r_mol, v_mep, order, ndim, nmep)
  implicit none
  double precision, intent(in), dimension(:,:) :: r_mol
  double precision, intent(in), dimension(:) :: v_mep
  integer, intent(in) :: order
  integer, intent(out) :: ndim
  integer, intent(out) :: nmep

  integer :: nat

  nat = size(r_mol, 2)
  ndim = nat
  if (order >= 1) then
      ndim = ndim + 3 * nat
  endif
  if (order >= 2) then
      ndim = ndim + 6 * nat
  endif

  nmep = size(v_mep)

end subroutine



subroutine geom_mat(r_mol, r_mep, v_mep, order, A, B)
  implicit none
  double precision, intent(in), dimension(:,:) :: r_mol
  double precision, intent(in), dimension(:,:) :: r_mep
  double precision, intent(in), dimension(:) :: v_mep
  integer :: order
  double precision, intent(out), dimension(:,:) :: A
  double precision, intent(inout), dimension(:) :: B

  integer :: ndim
  integer :: nmep
  double precision, dimension(:,:), allocatable :: Aaux

  call get_aux_dim(r_mol, v_mep, order, ndim, nmep)
  allocate(Aaux(ndim, nmep))
  call auxvec(r_mol, r_mep, v_mep, order, Aaux, B)
  A = matmul(Aaux, transpose(Aaux))
  deallocate(Aaux)

end subroutine geom_mat

subroutine constr_mat(r_mol, z_mol, constr, constr_q, constr_d, ndim, cdim, A, B)
  use qfit_variables, only : center_of_mass
  implicit none
  double precision, intent(in), dimension(:,:) :: r_mol
  double precision, intent(in), dimension(:) :: z_mol
  integer, intent(in) :: constr
  integer, intent(in) :: constr_q
  double precision, dimension(3), intent(in) :: constr_d
  integer, intent(in) :: ndim
  integer, intent(in) :: cdim
  double precision, intent(inout), dimension(:,:) :: A
  double precision, intent(inout), dimension(:)   :: B

  integer :: nat
  integer :: i, j
  nat = size(r_mol, 2)
  do i = 1, cdim
    if (constr >= 0 .and. i.eq.1) then
      A(:nat, ndim + i) = 1.0d0
      A(ndim + i, :nat) = 1.0d0
      B(ndim + i) = constr_q
    endif
    if (constr >= 1 .and. i.gt.1) then
      do j = 1, nat
        A(j, ndim + i) = center_of_mass(i-1) - r_mol(i-1,j)
        A(ndim + i, j) = A(j, ndim +i)
      enddo
      B(ndim + i) = constr_d(i-1)
    endif
  enddo

end subroutine

end module auxmat
