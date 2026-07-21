!==========================================================================================!
! meds_budget_check -- shared CONSERVATION-checking facility for the coupled fast loop           !
! (MEDS_COLUMN_DYNAMICS_DESIGN.md Part IV). Coupling makes conservation a CROSS-module property: !
! once the column driver sums cohort + ground + atmosphere fluxes into the canopy-air-space      !
! twins and threads soil<->plant water, the meaningful invariant is the COLUMN TOTAL, which no    !
! single kernel can see -- only the driver. This module gives everyone ONE closure predicate and  !
! a per-store accumulator, replacing the three bespoke inline `debug_error` guards currently       !
! re-implemented in soil_energy_step_implicit / cas_column_step_implicit / column_hydrology_flux.                    !
!                                                                                          !
! Units are per-store, per unit ground area: energy [J/m2], water [kg/m2], CO2 [umol/m2] or        !
! [mol/m2]. Compute is `pure`/GPU-safe; the hard error-stop wrapper is the only non-pure entry     !
! (Debug-only, mirroring the modules' `debug_error` discipline). It lives in `shared` (root of the !
! DAG, deps meds_kinds only) so every process module and the driver test conservation identically. !
! It is also the real, general-purpose replacement for the removed empty `src/utils/` placeholder. !
!==========================================================================================!
module meds_budget_check
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: budget_t, closure_ok, budget_imbalance, budget_accumulate, budget_assert, budget_check_stop

   !----- A budget accumulator for one store over an accumulation window (a fast step, or a day). -!
   type :: budget_t
      real(wp)    :: store0  = 0.0_wp    !< store at window start   [store-unit/m2]
      real(wp)    :: store1  = 0.0_wp    !< store at window end     [store-unit/m2]
      real(wp)    :: influx  = 0.0_wp    !< boundary INFLOW  over the window (amount or dt*rate)
      real(wp)    :: outflux = 0.0_wp    !< boundary OUTFLOW over the window
      real(wp)    :: resid   = 0.0_wp    !< (store1-store0) - (influx-outflux)   the imbalance
      real(wp)    :: worst   = 0.0_wp    !< running max |resid| over the window (for reporting)
      integer(ik) :: n_check = 0_ik      !< number of closure tests performed
      integer(ik) :: n_fail  = 0_ik      !< number that breached tolerance
   end type budget_t

contains

   !---------------------------------------------------------------------------------------!
   ! THE closure predicate: |resid| <= rtol*|scale| + atol. Mixed relative/absolute tolerance  !
   ! (matches the shipped column CO2 budget guard); `scale` is a representative store magnitude.   !
   ! `elemental` so it can test an array of stores in one call.                                  !
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
   ! residual, the running worst |resid|, and increments the check / fail counters. Pure         !
   ! (status carried in the accumulator; the hard stop is a separate Debug-only entry).           !
   !---------------------------------------------------------------------------------------!
   pure subroutine budget_accumulate(b, store0, store1, influx, outflux, dt, scale, rtol, atol)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: store0, store1, influx, outflux, dt, scale, rtol, atol
      b%store0  = store0 ; b%store1 = store1 ; b%influx = influx ; b%outflux = outflux
      b%resid   = budget_imbalance(store0, store1, influx, outflux, dt)
      b%worst   = max(b%worst, abs(b%resid))
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

end module meds_budget_check
