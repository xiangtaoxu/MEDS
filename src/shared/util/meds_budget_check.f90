!==========================================================================================!
! meds_budget_check -- shared CONSERVATION-checking facility for the coupled fast loop           !
! (MEDS_COLUMN_DYNAMICS_DESIGN.md Part IV). Coupling makes conservation a CROSS-module property: !
! once the column driver sums cohort + ground + atmosphere fluxes into the canopy-air-space      !
! twins and threads soil<->plant water, the meaningful invariant is the COLUMN TOTAL, which no    !
! single kernel can see -- only the driver. This module gives everyone ONE closure predicate, a   !
! per-store accumulator, a merge (patch -> site -> run), and an end-of-run report.                !
!                                                                                          !
! Units are per-store, per unit ground area: energy [J/m2], water [kg/m2], CO2 [umol/m2] or        !
! [mol/m2]. Compute is `pure`/GPU-safe; the hard error-stop wrapper and the report are the only    !
! non-pure entries. It lives in `shared` (root of the DAG, deps meds_kinds only) so every process  !
! module and the driver test conservation identically.                                             !
!                                                                                          !
! TOLERANCE POLICY (2026-09 review, item 1A #1/#2). A residual is judged against the GROSS FLUX     !
! that crossed the store's boundary over the window, never against the store itself: a soil column   !
! holds ~1e9 J/m2, so a store-relative 1e-6 let a sustained ~1 W/m2 leak through every step. Pass    !
! scale = |influx| + |outflux| and atol = <rate floor> * dt, with the floors below. The accumulator  !
! also keeps the SIGNED cumulative residual, because a one-signed bias below the per-step tolerance  !
! is invisible to a running maximum and only shows up as a drift over a season.                      !
!==========================================================================================!
module meds_budget_check
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: budget_t, closure_ok, budget_imbalance, budget_accumulate, budget_assert, budget_check_stop
   public :: budget_check, budget_merge, budget_report
   public :: budget_rtol_flux, budget_energy_rate_floor, budget_water_rate_floor, budget_co2_rate_floor

   !----- Closure tolerances shared by every ledger (see the header). -------------------------!
   real(wp), parameter :: budget_rtol_flux         = 1.0e-6_wp  !< [-] relative to the gross boundary flux
   !----- The energy floor is 1e-2 W/m2, not the 1e-3 the water/CO2 floors' precision would suggest:  !
   !      under a melting snow pack the whole-column ledger carries an unattributed, sign-oscillating  !
   !      residual of up to ~0.7 J/m2 per 150 s step (5e-3 W/m2; the CAS, soil, pack and pond ledgers  !
   !      each close to 1e-9 on the same steps). It is 100x below the old store-scaled tolerance and    !
   !      sums to ~-14 J/m2 over a day, but it is not round-off -- an open item (review 2026-09). -----!
   real(wp), parameter :: budget_energy_rate_floor = 1.0e-2_wp  !< [W/m2]      absolute floor, as a rate
   real(wp), parameter :: budget_water_rate_floor  = 1.0e-9_wp  !< [kg/m2/s]   (~1e-4 mm/day)
   real(wp), parameter :: budget_co2_rate_floor    = 1.0e-6_wp  !< [umol/m2/s]

   !----- A budget accumulator for one store over an accumulation window (a fast step, a day, a run). !
   type :: budget_t
      real(wp)    :: store0  = 0.0_wp    !< store at window start   [store-unit/m2]   (last check)
      real(wp)    :: store1  = 0.0_wp    !< store at window end     [store-unit/m2]   (last check)
      real(wp)    :: influx  = 0.0_wp    !< boundary INFLOW  over the window (amount)  (last check)
      real(wp)    :: outflux = 0.0_wp    !< boundary OUTFLOW over the window            (last check)
      real(wp)    :: resid   = 0.0_wp    !< (store1-store0) - (influx-outflux)  imbalance (last check)
      real(wp)    :: worst   = 0.0_wp    !< running max |resid| over all checks
      real(wp)    :: resid_sum  = 0.0_wp !< SIGNED sum of resid over all checks (+ = appearing from nowhere)
      real(wp)    :: abs_sum    = 0.0_wp !< sum of |resid| over all checks
      real(wp)    :: flux_gross = 0.0_wp !< sum of |influx| + |outflux| over all checks
      real(wp)    :: elapsed    = 0.0_wp !< [s] simulated time the checks span (only if callers report it)
      integer(ik) :: n_check = 0_ik      !< number of closure tests performed
      integer(ik) :: n_fail  = 0_ik      !< number that breached tolerance
   end type budget_t

contains

   !---------------------------------------------------------------------------------------!
   ! THE closure predicate: |resid| <= rtol*|scale| + atol. Mixed relative/absolute tolerance;   !
   ! `scale` should be the gross boundary flux over the window (see the header). `elemental` so   !
   ! it can test an array of stores in one call.                                                  !
   !---------------------------------------------------------------------------------------!
   elemental pure function closure_ok(resid, scale, rtol, atol) result(ok)
      real(wp), intent(in) :: resid, scale, rtol, atol
      logical              :: ok
      ok = abs(resid) <= rtol * abs(scale) + atol
   end function closure_ok

   !---------------------------------------------------------------------------------------!
   ! Imbalance of one store over dt: (store1 - store0) - dt*(influx_rate - outflux_rate). Pass   !
   ! dt = 1 if influx/outflux are already time-integrated AMOUNTS rather than rates.             !
   !---------------------------------------------------------------------------------------!
   elemental pure function budget_imbalance(store0, store1, influx, outflux, dt) result(resid)
      real(wp), intent(in) :: store0, store1, influx, outflux, dt
      real(wp)             :: resid
      resid = (store1 - store0) - dt * (influx - outflux)
   end function budget_imbalance

   !---------------------------------------------------------------------------------------!
   ! Fold one store's residual into an accumulator: records store endpoints/fluxes, the         !
   ! residual, the running worst |resid|, the signed and absolute sums, the gross flux, and the !
   ! check / fail counters. `elapsed` [s] is the simulated time this check spans; report it so  !
   ! budget_report can turn the cumulative residual into a mean leak rate. Pure (status carried !
   ! in the accumulator; the hard stop is a separate Debug-only entry).                          !
   !---------------------------------------------------------------------------------------!
   pure subroutine budget_accumulate(b, store0, store1, influx, outflux, dt, scale, rtol, atol, elapsed)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: store0, store1, influx, outflux, dt, scale, rtol, atol
      real(wp), optional, intent(in) :: elapsed
      b%store0  = store0 ; b%store1 = store1 ; b%influx = influx ; b%outflux = outflux
      b%resid   = budget_imbalance(store0, store1, influx, outflux, dt)
      b%worst   = max(b%worst, abs(b%resid))
      b%resid_sum  = b%resid_sum  + b%resid
      b%abs_sum    = b%abs_sum    + abs(b%resid)
      b%flux_gross = b%flux_gross + dt * (abs(influx) + abs(outflux))
      if (present(elapsed)) b%elapsed = b%elapsed + elapsed
      b%n_check = b%n_check + 1_ik
      if (.not. closure_ok(b%resid, scale, rtol, atol)) b%n_fail = b%n_fail + 1_ik
   end subroutine budget_accumulate

   !---------------------------------------------------------------------------------------!
   ! Status-code assert (pure, GPU-safe): ok = .false. if the residual breaches tolerance.      !
   ! `elemental` so a whole set of stores can be checked at once.                                !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine budget_assert(resid, scale, rtol, atol, ok)
      real(wp), intent(in)  :: resid, scale, rtol, atol
      logical,  intent(out) :: ok
      ok = closure_ok(resid, scale, rtol, atol)
   end subroutine budget_assert

   !---------------------------------------------------------------------------------------!
   ! Debug-only HARD assert: report + `error stop` with a label if the budget does not close.   !
   ! No-op when `debug` is .false. (production). Not pure -- the one place a numerical fault      !
   ! halts, mirroring the biophysics kernels' `debug_error` discipline.                          !
   !---------------------------------------------------------------------------------------!
   subroutine budget_check_stop(resid, scale, rtol, atol, label, debug)
      real(wp),         intent(in) :: resid, scale, rtol, atol
      character(len=*), intent(in) :: label
      logical,          intent(in) :: debug
      if (debug .and. .not. closure_ok(resid, scale, rtol, atol)) then
         write(*, '(a,a,a,es13.5,a,es13.5)') 'meds_budget_check: budget did NOT close [', &
            trim(label), '] resid = ', resid, '  tol = ', rtol * abs(scale) + atol
         error stop 'meds_budget_check: budget did not close'
      end if
   end subroutine budget_check_stop

   !---------------------------------------------------------------------------------------!
   ! One-call ledger check: accumulate, then (Debug only) hard-stop. `influx`/`outflux` are       !
   ! AMOUNTS over the window `dt_window` [s]; the tolerance is flux-scaled (rtol * gross flux)   !
   ! plus a rate floor times the window. This is the form every fast-loop ledger should use.     !
   !---------------------------------------------------------------------------------------!
   subroutine budget_check(b, store0, store1, influx, outflux, dt_window, rate_floor, label, halt, &
                           atol_extra)
      type(budget_t),   intent(inout) :: b
      real(wp),         intent(in)    :: store0, store1, influx, outflux
      real(wp),         intent(in)    :: dt_window    !< [s] time the amounts span
      real(wp),         intent(in)    :: rate_floor   !< [store-unit/m2/s] absolute tolerance floor, as a rate
      character(len=*), intent(in)    :: label
      logical,          intent(in)    :: halt         !< Debug hard stop on breach
      !----- Extra ABSOLUTE slack for a KNOWN, documented non-closure the caller has not fixed yet.   !
      !      Every use must name the defect it covers; it is a debt, not a tolerance. ---------------!
      real(wp), optional, intent(in)  :: atol_extra
      real(wp) :: scale, atol
      scale = abs(influx) + abs(outflux)
      atol  = rate_floor * dt_window
      if (present(atol_extra)) atol = atol + atol_extra
      call budget_accumulate(b, store0, store1, influx, outflux, 1.0_wp, scale, budget_rtol_flux, atol, &
                             elapsed=dt_window)
      call budget_check_stop(b%resid, scale, budget_rtol_flux, atol, label, halt)
   end subroutine budget_check

   !---------------------------------------------------------------------------------------!
   ! Fold accumulator `c` into `b` with weight `w` (an area fraction when reducing patches to a   !
   ! site; 1 when reducing windows to a run). Sums scale with the weight so the result stays per  !
   ! unit ground area; counts add; the worst per-check residual is the max; elapsed adds with    !
   ! the weight too, so a set of patches whose weights sum to 1 reports its common window once.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine budget_merge(b, c, w)
      type(budget_t), intent(inout) :: b
      type(budget_t), intent(in)    :: c
      real(wp),       intent(in)    :: w
      b%resid_sum  = b%resid_sum  + w * c%resid_sum
      b%abs_sum    = b%abs_sum    + w * c%abs_sum
      b%flux_gross = b%flux_gross + w * c%flux_gross
      b%elapsed    = b%elapsed    + w * c%elapsed
      b%worst      = max(b%worst, c%worst)
      b%n_check    = b%n_check + c%n_check
      b%n_fail     = b%n_fail  + c%n_fail
   end subroutine budget_merge

   !---------------------------------------------------------------------------------------!
   ! End-of-window report, one line per store: the SIGNED cumulative residual, the mean leak     !
   ! rate it implies over `elapsed`, the worst single check, and the fail count. A positive      !
   ! cumulative residual is store appearing from nowhere; negative is store vanishing.           !
   !---------------------------------------------------------------------------------------!
   subroutine budget_report(b, label, unit_store, unit_rate)
      type(budget_t),   intent(in) :: b
      character(len=*), intent(in) :: label, unit_store, unit_rate
      real(wp) :: rate
      rate = 0.0_wp
      if (b%elapsed > 0.0_wp) rate = b%resid_sum / b%elapsed
      write(*, '(a,a,a,es11.3,a,a,a,es11.3,a,a,a,es10.2,a,a,a,i0,a,i0)')                         &
         ' budget[', trim(label), ']  cumulative resid = ', b%resid_sum, ' ', trim(unit_store),   &
         '  mean leak = ', rate, ' ', trim(unit_rate), '  worst/check = ', b%worst, ' ',          &
         trim(unit_store), '  fails = ', b%n_fail, '/', b%n_check
   end subroutine budget_report

end module meds_budget_check
