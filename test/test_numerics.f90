!==========================================================================================!
! test_numerics -- unit tests for the shared numerical primitives + the conservation checker.  !
!   1. thomas_solve: recover a known solution of a diagonally-dominant tridiagonal system, with  !
!      the ACTUAL arrays longer than the active count n (the explicit-shape / self-contained path !
!      that let the solver move out of biophysics into shared).                                  !
!   2. quadratic_smaller_root: -> min(a,b) as theta -> 1; strictly below min for theta < 1.       !
!   3. adaptive_step_update: grow-capped (fmax), shrink-floored (fmin), identity-ish at err = 1.  !
!   4. meds_budget_check: closure predicate, imbalance algebra, accumulator fail-count, and the    !
!      no-op behaviour of the Debug hard-stop when debug = .false.                                 !
!==========================================================================================!
program test_numerics
   use meds_kinds,        only : wp, ik
   use meds_numerics,     only : thomas_solve, quadratic_smaller_root, adaptive_step_update,    &
                                 bisect_root
   use meds_budget_check, only : budget_t, closure_ok, budget_imbalance, budget_accumulate,    &
                                 budget_check_stop
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_thomas()
   call test_quadratic()
   call test_step_update()
   call test_bisect()
   call test_budget()

   if (nfail == 0_ik) then
      print '(a)', 'test_numerics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_numerics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) > atol) then
         print '(a,a,a,es15.7,a,es15.7)', '  FAIL ', name, ': got ', got, ' expected ', expect
         nfail = nfail + 1_ik
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (.not. cond) then
         print '(a,a,a,es15.7)', '  FAIL ', name, ': condition false, val = ', val
         nfail = nfail + 1_ik
      end if
   end subroutine check_true

   !----- 1. Tridiagonal solve of a 4-node system embedded in size-8 arrays. -----------------!
   subroutine test_thomas()
      integer(ik), parameter :: nmax = 8_ik, n = 4_ik
      real(wp) :: a(nmax), b(nmax), c(nmax), d(nmax), x(nmax), xtrue(n)
      integer(ik) :: k
      a = 0.0_wp ; b = 0.0_wp ; c = 0.0_wp ; d = 0.0_wp ; x = 0.0_wp
      xtrue = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      !----- Symmetric diagonally-dominant matrix: b=4, a=c=-1. d = A*xtrue. -----------------!
      do k = 1_ik, n
         a(k) = -1.0_wp ; b(k) = 4.0_wp ; c(k) = -1.0_wp
      end do
      d(1) = b(1)*xtrue(1) + c(1)*xtrue(2)
      do k = 2_ik, n-1_ik
         d(k) = a(k)*xtrue(k-1) + b(k)*xtrue(k) + c(k)*xtrue(k+1)
      end do
      d(n) = a(n)*xtrue(n-1) + b(n)*xtrue(n)
      call thomas_solve(a, b, c, d, x, n)          ! actual arrays are longer than n
      do k = 1_ik, n
         call check('thomas x', x(k), xtrue(k), 1.0e-12_wp)
      end do
   end subroutine test_thomas

   !----- 2. Smoothed-min quadratic. ---------------------------------------------------------!
   subroutine test_quadratic()
      real(wp) :: x1, x9
      x1 = quadratic_smaller_root(1.0_wp, 3.0_wp, 5.0_wp)     ! theta=1 -> exact min
      call check('quad theta=1 -> min', x1, 3.0_wp, 1.0e-12_wp)
      x9 = quadratic_smaller_root(0.9_wp, 3.0_wp, 5.0_wp)     ! co-limited: strictly below min, positive
      call check_true('quad theta<1 below min', x9 < 3.0_wp .and. x9 > 0.0_wp, x9)
   end subroutine test_quadratic

   !----- 3. Step-doubling factor clamps. ----------------------------------------------------!
   subroutine test_step_update()
      real(wp) :: f_ok, f_grow, f_shrink
      f_ok     = adaptive_step_update(1.0_wp,    0.9_wp, 0.25_wp, 4.0_wp)   ! err=1 -> safety
      call check('step err=1 -> safety', f_ok, 0.9_wp, 1.0e-12_wp)
      f_grow   = adaptive_step_update(1.0e-4_wp, 0.9_wp, 0.25_wp, 4.0_wp)   ! tiny err -> capped at fmax
      call check('step grow capped', f_grow, 4.0_wp, 1.0e-12_wp)
      f_shrink = adaptive_step_update(100.0_wp,  0.9_wp, 0.25_wp, 4.0_wp)   ! big err -> floored at fmin
      call check('step shrink floored', f_shrink, 0.25_wp, 1.0e-12_wp)
   end subroutine test_step_update

   !----- 5. bisect_root: recover a known root; a same-sign bracket returns converged=.false. -!
   subroutine test_bisect()
      real(wp) :: root
      logical  :: conv
      call bisect_root(quad_resid, 0.0_wp, 5.0_wp, 1.0e-10_wp, 100_ik, root, conv)
      call check_true('bisect converged', conv, 1.0_wp)
      call check('bisect root of x^2-4', root, 2.0_wp, 1.0e-8_wp)
      !----- Same-sign bracket: midpoint returned, converged = .false. -----------------------!
      call bisect_root(quad_resid, 3.0_wp, 5.0_wp, 1.0e-10_wp, 100_ik, root, conv)
      call check_true('bisect same-sign not converged', .not. conv, 0.0_wp)
      call check('bisect same-sign midpoint', root, 4.0_wp, 1.0e-14_wp)
   end subroutine test_bisect

   pure function quad_resid(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = x * x - 4.0_wp                 ! isolated root at x = 2 on [0, 5]
   end function quad_resid

   !----- 4. Conservation checker. -----------------------------------------------------------!
   subroutine test_budget()
      type(budget_t) :: b
      logical  :: ok_closed, ok_open
      real(wp) :: resid
      !----- Closure predicate (mixed rtol/atol). -------------------------------------------!
      ok_closed = closure_ok(1.0e-10_wp, 100.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      ok_open   = closure_ok(1.0_wp,     100.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      call check_true('closure_ok closed', ok_closed, 1.0_wp)
      call check_true('closure_ok open flagged', .not. ok_open, 0.0_wp)
      !----- Imbalance algebra: a perfectly closed store has zero residual. ------------------!
      resid = budget_imbalance(10.0_wp, 12.0_wp, 3.0_wp, 1.0_wp, 1.0_wp)     ! (12-10) - (3-1) = 0
      call check('imbalance closed = 0', resid, 0.0_wp, 1.0e-14_wp)
      !----- Accumulator: one closed + one open store -> n_fail = 1, worst > 0. --------------!
      call budget_accumulate(b, 10.0_wp, 12.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 12.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      call budget_accumulate(b, 10.0_wp, 15.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 12.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      call check_true('accumulator n_check=2', b%n_check == 2_ik, real(b%n_check, wp))
      call check_true('accumulator n_fail=1', b%n_fail == 1_ik, real(b%n_fail, wp))
      call check('accumulator worst', b%worst, 3.0_wp, 1.0e-14_wp)   ! second store off by 3
      !----- Debug hard-stop is a NO-OP when debug = .false. (must return, not halt). ---------!
      call budget_check_stop(1.0e6_wp, 1.0_wp, 1.0e-6_wp, 1.0e-8_wp, 'no-op path', .false.)
      call check_true('budget_check_stop no-op returns', .true., 1.0_wp)
   end subroutine test_budget

end program test_numerics
