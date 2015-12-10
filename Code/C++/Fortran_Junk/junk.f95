module junk
implicit none
integer, parameter :: dp = kind(1.0d0)

contains

subroutine increment(x, N)
  implicit none
  integer, intent(in) :: N
  real(dp), intent(inout), dimension(N) :: x

  x(1) = x(1) + 1.0
  x(2) = x(2) + 6.0
  x(5) = x(5) - 3.0
end subroutine increment

end module junk

