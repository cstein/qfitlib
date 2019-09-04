module tensor

implicit none

public :: t0, t1, t2

contains

function t0(r_a, r_b) result(var)
  double precision, intent(in), dimension(3) :: r_a
  double precision, intent(in), dimension(3) :: r_b
  double precision :: var

  double precision :: R2
  double precision, dimension(3) :: dr

  dr = r_b - r_a
  R2 = dot_product(dr, dr)
  var = 1.0d0 / sqrt(R2)

end function t0

function t1(r_a, r_b) result(var)
  double precision, intent(in), dimension(3) :: r_a
  double precision, intent(in), dimension(3) :: r_b
  double precision, dimension(3) :: var

  double precision :: R, R2
  double precision, dimension(3) :: dr

  dr = r_b - r_a
  R2 = dot_product(dr, dr)
  R = sqrt(R2)
  var = dr / (R*R2)

end function t1

function t2(r_a, r_b) result(var)
  double precision, intent(in), dimension(3) :: r_a
  double precision, intent(in), dimension(3) :: r_b
  double precision, dimension(6) :: var

  double precision :: R, R2
  double precision, dimension(3) :: dr

  dr = r_b - r_a
  R2 = dot_product(dr, dr)
  R = sqrt(R2)

  var(1) = 3 * dr(1) * dr(1) - R2   ! xx
  var(2) = 3 * dr(1) * dr(2)        ! xy
  var(3) = 3 * dr(1) * dr(3)        ! xz
  var(4) = 3 * dr(2) * dr(2) - R2   ! yy
  var(5) = 3 * dr(2) * dr(3)        ! yz
  var(6) = 3 * dr(3) * dr(3) - R2   ! zz

  var = var / (R*R2*R2)

end function t2


end module tensor
