!==========================================================================================!
! meds_soil_solver -- the hand-rolled tridiagonal (Thomas) sweep for the implicit soil-water   !
! backward-Euler step. A subroutine with a fixed-size `intent(out)` result (never an array-      !
! valued function fed into a call -- issue #7: silently wrong at -O2, segfault at -O0 on          !
! nvfortran). The BE matrix is diagonally dominant (the C*dz/dt diagonal), so no pivoting.        !
!==========================================================================================!
module meds_soil_solver
   use meds_kinds,            only : wp, ik
   use meds_biophysics_types, only : n_soil_layer_max
   implicit none
   private

   public :: thomas_solve

contains

   !---------------------------------------------------------------------------------------!
   ! Solve  a_k x_{k-1} + b_k x_k + c_k x_{k+1} = d_k  for k = 1..n. a(1) and c(n) ignored.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine thomas_solve(a, b, c, d, x, n)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: a(n_soil_layer_max), b(n_soil_layer_max)
      real(wp),    intent(in)  :: c(n_soil_layer_max), d(n_soil_layer_max)
      real(wp),    intent(out) :: x(n_soil_layer_max)
      real(wp)    :: cp(n_soil_layer_max), dp(n_soil_layer_max), denom
      integer(ik) :: k
      x = 0.0_wp
      cp(1) = c(1) / b(1)
      dp(1) = d(1) / b(1)
      do k = 2_ik, n
         denom = b(k) - a(k) * cp(k-1)
         cp(k) = c(k) / denom
         dp(k) = (d(k) - a(k) * dp(k-1)) / denom
      end do
      x(n) = dp(n)
      do k = n - 1_ik, 1_ik, -1_ik
         x(k) = dp(k) - cp(k) * x(k+1)
      end do
   end subroutine thomas_solve

end module meds_soil_solver
