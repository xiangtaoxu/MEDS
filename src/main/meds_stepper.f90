!==========================================================================================!
! meds_stepper -- the master stepper (the ED2 ed_model analogue; seed of an all-process        !
! top-level loop). Lives in src/driver/, the home of the top-level utilities that wire the      !
! process modules together.                                                                    !
!                                                                                          !
! It owns only the CADENCE: each step it is told whether a month/year rolled over, and it       !
! drives the process modules on the appropriate timescale. The slow tier is the thin              !
! meds_slow_dynamics coordinator (advance_slow_dynamics), which sequences vegetation dynamics       !
! and slow soil-carbon biogeochemistry as peer domains (MEDS_SLOW_DYNAMICS_DESIGN.md Part II).     !
!==========================================================================================!
module meds_stepper
   use meds_kinds,                only : wp
   use meds_config,               only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_slow_dynamics,        only : advance_slow_dynamics
   use meds_fast_dynamics,        only : fast_context_t, fast_dynamics
   use meds_time,                 only : meds_time_t, day_of_year
   use meds_forcing_types,        only : met_driver_t
   use meds_output_types,         only : output_manager_t
   use meds_budget_check,         only : budget_t
   use meds_slow_ledger,          only : slow_ledger_t
   implicit none
   private

   public :: advance_one_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one step. The caller's calendar supplies the cadence flags; this routine passes  !
   ! them to the process drivers. When fast biophysics is on and a fast context is supplied,   !
   ! the sub-daily fast loop runs over the per-patch reservoirs BEFORE the slow loop (so a     !
   ! later fast->slow carbon handoff can hand daily-accumulated GPP to vegetation dynamics).   !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx, met_drv, step_start, mgr, &
                               run_energy_budget, run_water_budget, slow_ledger, worst_rh_seam_gap)
      type(site_t),         intent(inout) :: site
      type(meds_config_t),  intent(in)    :: cfg
      logical,              intent(in)    :: is_new_month, is_new_year
      type(fast_context_t), intent(in),    optional :: fast_ctx
      type(met_driver_t),   intent(inout), optional :: met_drv     !< live met reader (when forcing_on)
      type(meds_time_t),    intent(in),    optional :: step_start  !< calendar time at the start of this slow step
      type(output_manager_t), intent(inout), optional :: mgr       !< FAST-tier staging (forwarded to the fast loop)
      type(budget_t), intent(inout), optional :: run_energy_budget, run_water_budget !< run-level ledgers (forwarded)
      !----- The SLOW tier's own ledger (plan §10.2). Its peers above accumulate per-fast-step   !
      !      flux residuals; this one snapshots the site store across the slow step, which is    !
      !      the window neither of them can see. Same lifetime, same place in the plumbing.      !
      type(slow_ledger_t), intent(inout), optional :: slow_ledger
      !----- The soil-carbon SEAM check (design Part II section 9): |the daily pool debit - the fast    !
      !      loop's own accumulated Rh|, worst over patches. Both sides read the same frozen pool and   !
      !      the same per-pool xi integral, so it is ~0 BY CONSTRUCTION -- an assertion guard, not a    !
      !      correction. It was computed and DISCARDED on every step of every run until now, because    !
      !      nothing above the slow driver asked for it.  ----------------------------------------------!
      real(wp), intent(out), optional :: worst_rh_seam_gap

      !----- Fast loop: sub-daily biophysics over the state-hub reservoirs. When fast biophysics   !
      !      is ON a fast context MUST be supplied: the old `.and. present(fast_ctx)` SILENTLY       !
      !      skipped the loop, leaving gpp_accum=0 so carbon-mode growth ran on zero GPP (BUG1).     !
      !      Fail loud instead of silently wrong. met_drv/step_start are forwarded only when present !
      !      (forcing_on); absent -> fast_dynamics runs the constant-forcing MVP.                    !
      if (cfg%fast_biophysics_on) then
         if (.not. present(fast_ctx))                                                              &
            error stop 'advance_one_step: fast_biophysics_on=.true. but no fast_context supplied'
         if (present(met_drv) .and. present(step_start)) then
            call fast_dynamics(site, fast_ctx, cfg, met_drv=met_drv, step_start=step_start, mgr=mgr, &
                               run_energy_budget=run_energy_budget, run_water_budget=run_water_budget)
         else
            call fast_dynamics(site, fast_ctx, cfg, run_energy_budget=run_energy_budget,           &
                               run_water_budget=run_water_budget)
         end if
      end if

      !----- Slow loop: vegetation dynamics (rate assembly + demographic application), gated on the   !
      !      master slow_on freeze (holds cohort/patch state static while the fast loop still runs;   !
      !      docs/dev_plans/MEDS_SLOW_DYNAMICS_DESIGN.md Part I). Leaf phenology is the FIRST step      !
      !      INSIDE vegetation_dynamics (the folded phenology driver) and runs UNCONDITIONALLY          !
      !      whenever a step-start day-of-year is available -- pass doy whenever step_start is          !
      !      supplied; vegetation_dynamics itself no-ops the phenology advance when doy is absent        !
      !      (no calendar context, e.g. a bare test call). -----------------------------------------!
      if (cfg%slow_on) then
         !----- `rho_air` values the canopy-air store. It lives on the fast context (site-uniform,   !
         !      from the forcing), so the ledger never needs state of its own to carry it; without  !
         !      a fast context there is no canopy air worth valuing and the term stays at zero.     !
         if (present(step_start)) then
            if (present(fast_ctx)) then
               call advance_slow_dynamics(site, cfg, is_new_month, is_new_year, doy=day_of_year(step_start), &
                                          ledger=slow_ledger, rho_air=fast_ctx%rho_air,                  &
                                          worst_rh_seam_gap=worst_rh_seam_gap)
            else
               call advance_slow_dynamics(site, cfg, is_new_month, is_new_year, doy=day_of_year(step_start), &
                                          ledger=slow_ledger, worst_rh_seam_gap=worst_rh_seam_gap)
            end if
         else
            if (present(fast_ctx)) then
               call advance_slow_dynamics(site, cfg, is_new_month, is_new_year, ledger=slow_ledger,  &
                                          rho_air=fast_ctx%rho_air,                                  &
                                          worst_rh_seam_gap=worst_rh_seam_gap)
            else
               call advance_slow_dynamics(site, cfg, is_new_month, is_new_year, ledger=slow_ledger,  &
                                          worst_rh_seam_gap=worst_rh_seam_gap)
            end if
         end if
      end if
   end subroutine advance_one_step

end module meds_stepper
