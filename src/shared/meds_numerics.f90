!==========================================================================================!
! meds_numerics -- shared, dependency-free numerical PRIMITIVES (root of the library DAG:      !
! uses meds_kinds only). Consolidated per MEDS_COLUMN_DYNAMICS_DESIGN.md Part III so that the   !
! biophysics AND plant libraries -- which both link `shared` ONLY and cannot see each other --  !
! share ONE tridiagonal solve plus the small reusable scalar kernels, instead of re-deriving    !
! them behind the state/process wall.                                                          !
!                                                                                          !
! Every routine is `pure` and issue-#7-safe: `thomas_solve` is a subroutine with an intent(out) !
! array (never an array-valued function fed into a call -- silently wrong at -O2, segfault at    !
! -O0 on nvfortran). It is DIMENSION-FREE (explicit-shape dummies keyed on the passed `n`), so   !
! it needs no fixed `n_soil_layer_max` constant and thus no biophysics dependency -- the move    !
! out of the old `meds_soil_solver` that made this consolidation possible.                       !
!==========================================================================================!
module meds_numerics
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: thomas_solve, quadratic_smaller_root, adaptive_step_update

contains

   !---------------------------------------------------------------------------------------!
   ! Tridiagonal (Thomas) sweep: solve  a_k x_{k-1} + b_k x_k + c_k x_{k+1} = d_k, k = 1..n.  !
   ! a(1) and c(n) are ignored. No pivoting -- the backward-Euler matrices that call this      !
   ! (soil heat, soil water) carry a diagonally-dominant C*dz/dt diagonal. Explicit-shape       !
   ! dummies: the actual arrays may be LONGER than n (the first n elements associate), so the    !
   ! callers keep their fixed-size (n_soil_layer_max) storage and pass the active count `n`.      !
   !---------------------------------------------------------------------------------------!
   pure subroutine thomas_solve(a, b, c, d, x, n)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: a(n), b(n), c(n), d(n)
      real(wp),    intent(out) :: x(n)
      real(wp)    :: cp(n), dp(n), denom
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

   !---------------------------------------------------------------------------------------!
   ! Smaller root of the co-limitation quadratic  theta*x^2 - (a+b)*x + a*b = 0. The smaller   !
   ! root is the smooth analogue of min(a,b): it -> min(a,b) as theta -> 1 and the co-limited   !
   ! (harmonic-like) mean for smaller theta. Discriminant is guarded >= 0. `elemental` so it     !
   ! sweeps arrays. (Used by Farquhar co-limitation and pressure-volume inverses.)               !
   !---------------------------------------------------------------------------------------!
   elemental pure function quadratic_smaller_root(theta, a, b) result(x)
      real(wp), intent(in) :: theta, a, b
      real(wp)             :: x, s, disc
      s    = a + b
      disc = max(s * s - 4.0_wp * theta * a * b, 0.0_wp)
      x    = (s - sqrt(disc)) / (2.0_wp * theta)
   end function quadratic_smaller_root

   !---------------------------------------------------------------------------------------!
   ! Step-doubling size-update FACTOR (multiply the substep h by this). Classic error-scaled    !
   ! controller  factor = safety * err**(-1/2), clamped to [fmin, fmax]. The symmetric clamp     !
   ! reproduces BOTH shipped controllers: on ACCEPT (err <= 1) the raw factor >= safety >= fmin  !
   ! so `min(fmax, raw)` binds (hydrology form); on REJECT (err > 1) raw < 1 < fmax so           !
   ! `max(fmin, raw)` binds (shrink). Callers pass their OWN safety/fmin/fmax -- these DIFFER     !
   ! across kernels (hydrology fmin=0.25, hydraulics fmin=0.20), so they are arguments, never     !
   ! hard-coded here.                                                                            !
   !---------------------------------------------------------------------------------------!
   elemental pure function adaptive_step_update(err, safety, fmin, fmax) result(factor)
      real(wp), intent(in) :: err, safety, fmin, fmax
      real(wp)             :: factor
      !----- Floor err before the negative power so a perfectly-converged step (err = 0, e.g. a  !
      !      zero-flux substep) yields factor = fmax instead of 0**(-1/2) = Inf / a SIGFPE under   !
      !      -fpe0. Callers already floor err, but keep the shared primitive self-safe. -----------!
      factor = min(fmax, max(fmin, safety * max(err, tiny(err)) ** (-0.5_wp)))
   end function adaptive_step_update

end module meds_numerics
