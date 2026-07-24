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
   use meds_constants, only : safe_exp, tiny_num
   implicit none
   private

   public :: thomas_solve, quadratic_smaller_root, adaptive_step_update, bisect_root
   public :: gauss_legendre_7
   public :: logistic, clamp01, clamp
   public :: matrix_exp, matrix_exp_fixed, matmul_sq

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
   ! a(1) and c(n) are ignored. No PARTIAL pivoting -- the backward-Euler matrices that call    !
   ! this (soil heat, soil water) carry a diagonally-dominant C*dz/dt diagonal, so a row swap    !
   ! is never needed. A ZERO-PIVOT GUARD floors |b(1)| and each elimination |denom| away from 0   !
   ! (sign-preserving, mirroring meds_soil_biogeochem's gaussian_solve) so a degenerate/          !
   ! misconfigured caller divides by tiny_num instead of exactly 0 -- a NaN/Inf firewall, not a   !
   ! numerical solve (a genuinely singular system still returns garbage, just finite garbage).    !
   ! `pure` throughout (no error stop: this solver is called from the fast-loop RHS, which must   !
   ! stay side-effect-free/GPU-eligible -- same discipline as meds_plant_hydraulics's             !
   ! plant_water_tendency). The guard is a NO-OP for every diagonally-dominant matrix this solves  !
   ! today (|b(1)|, |denom| >> tiny_num always), so existing callers are bit-identical.            !
   ! Explicit-shape dummies: the actual arrays may be LONGER than n (the first n elements          !
   ! associate), so the callers keep their fixed-size (n_soil_layer_max) storage and pass the       !
   ! active count `n`.                                                                              !
   !---------------------------------------------------------------------------------------!
   pure subroutine thomas_solve(a, b, c, d, x, n)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: a(n), b(n), c(n), d(n)
      real(wp),    intent(out) :: x(n)
      real(wp)    :: cp(n), dp(n), denom
      integer(ik) :: k
      x = 0.0_wp
      denom = sign(max(abs(b(1)), tiny_num), b(1))
      cp(1) = c(1) / denom
      dp(1) = d(1) / denom
      do k = 2_ik, n
         denom = b(k) - a(k) * cp(k-1)
         denom = sign(max(abs(denom), tiny_num), denom)
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
   elemental pure function adaptive_step_update(err, safety, fmin, fmax, p_order) result(factor)
      real(wp),    intent(in)           :: err, safety, fmin, fmax
      !----- Embedded-pair LOWER order (MEDS_ED2_RK45_DESIGN.md sec 6): the classical I-controller  !
      !      exponent is -1/(p_order+1). Default 1 reproduces every EXISTING caller's behaviour       !
      !      byte-for-byte (ARK's ARS(2,2,2) embedded estimate is 1st-order, exponent -1/2 exactly     !
      !      as before); Cash-Karp's embedded 4th-order estimate needs p_order=4 (exponent -1/5) --      !
      !      passing 1 there would silently under/over-react to the SAME normalized error, not          !
      !      "wrong" in a way that breaks accept/reject, but not matching the method's true             !
      !      asymptotic step-size law either. --------------------------------------------------------!
      integer(ik), intent(in), optional :: p_order
      real(wp)             :: factor, expo
      integer(ik) :: p
      p = 1_ik ; if (present(p_order)) p = p_order
      expo = -1.0_wp / real(p + 1_ik, wp)
      !----- Floor err before the negative power so a perfectly-converged step (err = 0, e.g. a  !
      !      zero-flux substep) yields factor = fmax instead of 0**expo = Inf / a SIGFPE under      !
      !      -fpe0. Callers already floor err, but keep the shared primitive self-safe. -----------!
      factor = min(fmax, max(fmin, safety * max(err, tiny(err)) ** expo))
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

   !---------------------------------------------------------------------------------------!
   ! Matrix exponential exp(A) for a fixed-size m x m matrix, via scaling-and-squaring + a    !
   ! fixed-order Taylor series (issue-#7-safe: dimension-free explicit-shape dummies keyed on  !
   ! the passed `m`, no array-valued function result fed into a call). Promoted from            !
   ! meds_soil_biogeochem's CENTURY spin-up EXPM (MEDS_NUMERICS_SCOPING.md QW4) so every         !
   ! implicit-RC scheme -- soil carbon, and any future multi-node hydraulics topology -- shares   !
   ! ONE primitive instead of re-deriving it behind the state/process wall. Two entry points     !
   ! share one Taylor/squaring core:                                                             !
   !   * matrix_exp       -- ADAPTIVE squaring count s, chosen from ||A||_inf so the scaled        !
   !                         argument A/2^s is small enough for the Taylor series to converge      !
   !                         to machine precision (the general-purpose, data-dependent path).      !
   !   * matrix_exp_fixed -- a FIXED squaring count s (an explicit argument, not computed from      !
   !                         the matrix norm), so every call takes the SAME trip count regardless   !
   !                         of the data -- the GPU/warp-uniform sibling: a device kernel batching   !
   !                         this over many lanes needs one shared, data-INDEPENDENT loop bound.     !
   !                         Over-squaring is harmless (the Taylor argument only gets smaller, so    !
   !                         it is MORE accurate, not less); the caller picks s as a safe upper       !
   !                         bound for its problem class. matrix_exp delegates to this with its        !
   !                         own computed s, so the two paths share one Taylor/squaring core (and      !
   !                         are bit-identical to the pre-promotion meds_soil_biogeochem version).     !
   !---------------------------------------------------------------------------------------!
   pure subroutine matrix_exp(a_in, e, m)
      integer(ik), intent(in)  :: m
      real(wp),    intent(in)  :: a_in(m, m)
      real(wp),    intent(out) :: e(m, m)
      real(wp)    :: nrm, scale
      integer(ik) :: i, j, s

      !----- scaling: pick s so that ||A/2^s||_inf <= 1/2 (Taylor converges fast). -----------!
      nrm = 0.0_wp
      do i = 1_ik, m
         scale = 0.0_wp
         do j = 1_ik, m
            scale = scale + abs(a_in(i, j))
         end do
         nrm = max(nrm, scale)
      end do
      s = 0_ik
      scale = 1.0_wp
      do while (nrm * scale > 0.5_wp)
         s = s + 1_ik
         scale = scale * 0.5_wp
      end do
      call matrix_exp_fixed(a_in, e, m, s)
   end subroutine matrix_exp

   !----- exp(A) via a FIXED squaring count s (warp-uniform: same trip count for every call). --!
   pure subroutine matrix_exp_fixed(a_in, e, m, s)
      integer(ik), intent(in)  :: m, s
      real(wp),    intent(in)  :: a_in(m, m)
      real(wp),    intent(out) :: e(m, m)
      integer(ik), parameter :: n_taylor = 12_ik
      real(wp)    :: a(m, m), term(m, m), tmp(m, m)
      integer(ik) :: i, k, sq

      a = a_in * (2.0_wp ** (-s))                       ! A / 2^s (s=0 => A unscaled)

      !----- Taylor: E = I + A + A^2/2! + ... (term_{k} = term_{k-1} * A / k). ---------------!
      e = 0.0_wp
      term = 0.0_wp
      do i = 1_ik, m
         e(i, i)    = 1.0_wp
         term(i, i) = 1.0_wp
      end do
      do k = 1_ik, n_taylor
         call matmul_sq(term, a, tmp, m)                ! tmp = term * A
         term = tmp / real(k, wp)
         e = e + term
      end do
      !----- squaring: E = E^(2^s). ---------------------------------------------------------!
      do sq = 1_ik, s
         call matmul_sq(e, e, tmp, m)
         e = tmp
      end do
   end subroutine matrix_exp_fixed

   !----- Explicit square-matrix product C = A*B (issue-#7-safe: no array-valued function temp). --!
   pure subroutine matmul_sq(a, b, c, m)
      integer(ik), intent(in)  :: m
      real(wp),    intent(in)  :: a(m, m), b(m, m)
      real(wp),    intent(out) :: c(m, m)
      integer(ik) :: i, j, k
      real(wp)    :: acc
      do j = 1_ik, m
         do i = 1_ik, m
            acc = 0.0_wp
            do k = 1_ik, m
               acc = acc + a(i, k) * b(k, j)
            end do
            c(i, j) = acc
         end do
      end do
   end subroutine matmul_sq

end module meds_numerics
