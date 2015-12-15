module junk
implicit none
integer, parameter :: dp = kind(1.0d0)

contains

subroutine increment(x, M, N)
  integer :: i
  integer :: j
  integer, intent(in) :: M
  integer, intent(in) :: N
  real(dp), intent(inout), dimension(M, N) :: x

  do j = 1, N
    do i = 1, M
      x(i, j) = x(i, j) + i + 3*j
    end do
  end do
end subroutine increment

end module junk

