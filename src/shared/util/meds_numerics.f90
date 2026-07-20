!==========================================================================================!
! meds_numerics -- shared numerical PRIMITIVES near the root of the library DAG (uses meds_kinds !
! + meds_constants for safe_exp; meds_constants itself uses meds_kinds only, so no cycle).       !
! Consolidated per MEDS_COLUMN_DYNAMICS_DESIGN.md Part III so that the                            !
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
   use meds_kinds,     only : wp, ik
   use meds_constants, only : safe_exp
   implicit none
   private

   public :: thomas_solve, quadratic_smaller_root, adaptive_step_update, bisect_root
   public :: gauss_legendre_7
   public :: logistic, clamp01, clamp

   !----- Interface of a pure scalar function f(x) passed to bisect_root / gauss_legendre_7. -!
   abstract interface
      pure function scalar_fn(x) result(y)
         import :: wp
         real(wp), intent(in) :: x
         real(wp)             :: y
      end function scalar_fn
   end interface

   !----- 7-point Gauss-Legendre nodes/weights on [-1,1] (fixed-order quadrature). ----------!
   integer(ik), parameter :: NG_GL7 = 7_ik
   real(wp), parameter :: gl7_x(NG_GL7) = [ -0.9491079123427585_wp, -0.7415311855993945_wp,     &
                                            -0.4058451513773972_wp,  0.0000000000000000_wp,     &
                                             0.4058451513773972_wp,  0.7415311855993945_wp,     &
                                             0.9491079123427585_wp ]
   real(wp), parameter :: gl7_w(NG_GL7) = [  0.1294849661688697_wp,  0.2797053914892766_wp,     &
                                             0.3818300505051189_wp,  0.4179591836734694_wp,     &
                                             0.3818300505051189_wp,  0.2797053914892766_wp,     &
                                             0.1294849661688697_wp ]

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
      x    = (s - sqrt(disc)) / (2.0_wp * max(theta, tiny(theta)))    ! theta in (0,1]; floor guards /0
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

   !---------------------------------------------------------------------------------------!
   ! Bracket-and-bisect root finder for a pure scalar residual f. Evaluates f at the two     !
   ! endpoints; on a SAME-SIGN bracket it returns converged = .false. and root = the bracket  !
   ! midpoint (so a caller can switch residual providers and retry). Otherwise it bisects,    !
   ! keeping the sign change, until the bracket width falls below tol or max_iter is hit;      !
   ! root is the final bracket midpoint either way. #7-safe: scalar out-args, no array temps.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine bisect_root(f, lo, hi, tol, max_iter, root, converged)
      procedure(scalar_fn)     :: f
      real(wp),    intent(in)  :: lo, hi, tol
      integer(ik), intent(in)  :: max_iter
      real(wp),    intent(out) :: root
      logical,     intent(out) :: converged
      real(wp)    :: a, b, flo, fhi, mid, fmid
      integer(ik) :: it
      a = lo ; b = hi
      flo = f(a) ; fhi = f(b)
      converged = .false. ; root = 0.5_wp * (a + b)
      if (flo * fhi <= 0.0_wp) then
         do it = 1_ik, max_iter
            mid = 0.5_wp * (a + b) ; fmid = f(mid)
            if (flo * fmid <= 0.0_wp) then ; b = mid ; else ; a = mid ; flo = fmid ; end if
            if (b - a < tol) exit
         end do
         root = 0.5_wp * (a + b) ; converged = (b - a < tol)
      end if
   end subroutine bisect_root

   !---------------------------------------------------------------------------------------!
   ! Fixed 7-point Gauss-Legendre quadrature of a pure scalar integrand f over [a,b]. Exact    !
   ! for polynomials up to degree 2*7-1 = 13; for the smooth integrands here (e.g. the           !
   ! Kirchhoff matric-flux integrand 1/(1+u^a)) it is far below any modeling tolerance. The       !
   ! [-1,1] nodes are affine-mapped to [a,b]. Pure + callback (mirrors bisect_root); intended for  !
   ! off-hot-path use (BUILD lookup tables / oracle), not a per-step kernel.                        !
   !---------------------------------------------------------------------------------------!
   pure function gauss_legendre_7(f, a, b) result(integral)
      procedure(scalar_fn) :: f
      real(wp), intent(in) :: a, b
      real(wp)    :: integral, mid, half, acc
      integer(ik) :: g
      mid  = 0.5_wp * (a + b)
      half = 0.5_wp * (b - a)
      acc  = 0.0_wp
      do g = 1_ik, NG_GL7
         acc = acc + gl7_w(g) * f(mid + half * gl7_x(g))
      end do
      integral = half * acc
   end function gauss_legendre_7

   !---------------------------------------------------------------------------------------!
   ! Numerically-safe logistic (smooth 0->1 step) -- the shared home of the 1/(1+exp(-z))     !
   ! primitive that phenology, evergreen turnover suppression, and any smooth gate reuse.      !
   ! `safe_exp` (meds_constants) clamps the argument so `-fpe0`/`-Ktrap=fp` never overflow.    !
   !---------------------------------------------------------------------------------------!
   elemental pure function logistic(z) result(f)
      real(wp), intent(in) :: z
      real(wp)             :: f
      f = 1.0_wp / (1.0_wp + safe_exp(-z))
   end function logistic

   !----- Clamp to [0,1]. ------------------------------------------------------------------!
   elemental pure function clamp01(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = min(1.0_wp, max(0.0_wp, x))
   end function clamp01

   !----- General clamp to [lo,hi] (lo <= hi assumed by the caller). -----------------------!
   elemental pure function clamp(x, lo, hi) result(y)
      real(wp), intent(in) :: x, lo, hi
      real(wp)             :: y
      y = min(hi, max(lo, x))
   end function clamp

end module meds_numerics
