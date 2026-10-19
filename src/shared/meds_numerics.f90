! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_numerics -- small numerical helpers every layer shares: clamps, a smooth step, the   !
! smooth minimum of two limits, the step-size factor, weighted means, an ordering, and the   !
! test that an operator kept what it only moves. The linear solvers are in                   !
! meds_linear_algebra.                                                                      !
!==========================================================================================!
module meds_numerics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : safe_exp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private

   public :: quadratic_smaller_root, adaptive_step_update
   public :: logistic, clamp01, clamp, weighted_mean
   public :: ascending_order
   public :: conserved, ROUNDOFF_RTOL, pair_mean

   !----- The round-off that summing or rescaling a stored amount may leave, relative to it:  !
   !      the tolerance of conserved() for operators that only move an amount. --------------!
   real(wp), parameter :: ROUNDOFF_RTOL = 1.0e-12_wp

contains

   !----- ord(k) is the index of the k-th smallest of x(1:n), the larger index first among   !
   !      equals (the cohorts from the canopy bottom up). The cohorts are kept tallest first,  !
   !      so the model's input is reversed in O(n); any other input takes a selection sort,   !
   !      which gives the same order on a non-increasing x. -------------------------------!
   pure subroutine ascending_order(x, n, ord)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: x(n)
      integer(ik), intent(out) :: ord(n)
      logical     :: used(n)
      real(wp)    :: xmin
      integer(ik) :: j, k, imin
      if (.not. any(x(1:n-1) < x(2:n))) then
         ord = [(n - k + 1_ik, k = 1_ik, n)]
         return
      end if
      used = .false.
      do k = 1_ik, n
         imin = 0_ik
         xmin = huge(1.0_wp)
         do j = 1_ik, n
            if (.not. used(j) .and. x(j) <= xmin) then
               xmin = x(j)
               imin = j
            end if
         end do
         ord(k) = imin
         used(imin) = .true.
      end do
   end subroutine ascending_order

   !----- The smaller root of theta x^2 - (a + b) x + a b = 0: a smooth min(a, b), equal to it  !
   !      at theta = 1 (docs/science/leaf_gas_exchange.md §2). ----------------------------!
   elemental pure function quadratic_smaller_root(theta, a, b) result(x)
      real(wp), intent(in) :: theta, a, b
      real(wp)             :: x, s, disc
      s    = a + b
      disc = max(s * s - 4.0_wp * theta * a * b, 0.0_wp)
      x    = (s - sqrt(disc)) / (2.0_wp * max(theta, tiny(theta)))    ! theta in (0,1]; floor guards /0
   end function quadratic_smaller_root

   !----- The factor a step size is multiplied by: safety * err^(-1/(p_order + 1)), clamped to  !
   !      [fmin, fmax] (docs/science/numerical_scheme.md eq. (4)). --------------------------!
   elemental pure function adaptive_step_update(err, safety, fmin, fmax, p_order) result(factor)
      real(wp),    intent(in)           :: err, safety, fmin, fmax
      integer(ik), intent(in), optional :: p_order   !< order of the embedded error estimate (default 1)
      real(wp)             :: factor, expo
      integer(ik) :: p
      p = 1_ik
      if (present(p_order)) p = p_order
      expo = -1.0_wp / real(p + 1_ik, wp)
      factor = min(fmax, max(fmin, safety * max(err, tiny(err)) ** expo))   ! err = 0 gives fmax
   end function adaptive_step_update

   !----- The logistic 1 / (1 + exp(-z)), a smooth step from 0 to 1                            !
   !      (docs/science/plant_phenology.md eq. (2)). ---------------------------------------!
   elemental pure function logistic(z) result(f)
      real(wp), intent(in) :: z
      real(wp)             :: f
      f = 1.0_wp / (1.0_wp + safe_exp(-z))
   end function logistic

   !----- x clamped to [0, 1]. -----------------------------------------------------------!
   elemental pure function clamp01(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = min(1.0_wp, max(0.0_wp, x))
   end function clamp01

   !----- x clamped to [lo, hi], with lo <= hi. --------------------------------------------!
   elemental pure function clamp(x, lo, hi) result(y)
      real(wp), intent(in) :: x, lo, hi
      real(wp)             :: y
      y = min(hi, max(lo, x))
   end function clamp

   !----- The sum of x(k) w(k) over k = 1..n, in that order: a weighted mean when the weights   !
   !      sum to one. ---------------------------------------------------------------------!
   pure function weighted_mean(x, w, n) result(xbar)
      real(wp),    intent(in) :: x(:), w(:)
      integer(ik), intent(in) :: n
      real(wp) :: xbar
      integer(ik) :: k
      xbar = 0.0_wp
      do k = 1_ik, n
         xbar = xbar + x(k) * w(k)
      end do
   end function weighted_mean

   !----- True when an amount an operator only moves is unchanged to round-off:                !
   !      |after - before| <= rtol * magnitude, plus the smallest difference the number format  !
   !      resolves. A NaN fails (docs/science/numerical_scheme.md §9). ------------------------!
   elemental pure logical function conserved(before, after, magnitude, rtol) result(ok)
      real(wp), intent(in) :: before, after   !< the amount before and after the operator
      real(wp), intent(in) :: magnitude       !< sum of the absolute values of its terms
      real(wp), intent(in) :: rtol            !< relative tolerance
      ok = .not. (ieee_is_nan(before) .or. ieee_is_nan(after) .or. ieee_is_nan(magnitude))
      if (ok) ok = abs(after - before) <= rtol * magnitude + tiny(magnitude) / epsilon(magnitude)
   end function conserved

   !----- The weighted mean of two values, written so it stays between them for any weights    !
   !      however small; with no weight at all it keeps a. There is no size threshold: a       !
   !      leafless cohort's tiny leaf carbon still weighs its traits in a fusion. ------------!
   elemental pure real(wp) function pair_mean(a, wa, b, wb) result(m)
      real(wp), intent(in) :: a, wa   !< the first value and its weight (>= 0)
      real(wp), intent(in) :: b, wb   !< the second value and its weight (>= 0)
      if (wa + wb > 0.0_wp) then
         m = a + (b - a) * (wb / (wa + wb))
      else
         m = a
      end if
   end function pair_mean

end module meds_numerics
