! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_stepper -- the master stepper (the ED2 ed_model analogue). Lives in src/main/, beside    !
! the program that drives it.                                                                   !
!                                                                                               !
! It owns only the CADENCE: one slow step (advance_one_step) drives the fast loop and then the  !
! slow tier, the thin meds_slow_dynamics coordinator (advance_slow_dynamics), which sequences   !
! vegetation dynamics and slow soil-carbon biogeochemistry as peer domains                      !
! (docs/science/order_of_processes.md §2-§4). A calendar boundary between two steps is its own  !
! call (advance_boundary, §5): the caller makes it after the output has read the step that      !
! ended there.                                                                                  !
!==========================================================================================!
module meds_stepper
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_slow_dynamics,        only : advance_slow_dynamics, advance_boundary_dynamics, settle_tissue_water
   use meds_slow_ledger,          only : slow_ledger_t
   use meds_fast_dynamics,        only : fast_dynamics
   use meds_time,                 only : meds_time_t
   use meds_met_types,        only : met_source_t, met_cursor_t
   use meds_output_types,         only : output_buffers_t
   use meds_budget_check,         only : budget_t
   implicit none
   private

   public :: advance_one_step, advance_boundary

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one step. When fast biophysics is on, the sub-daily fast loop runs over the       !
   ! per-patch reservoirs BEFORE the slow loop, so vegetation dynamics allocates the GPP the  !
   ! fast loop accumulated over the step.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(poly, cfg, met_src, met_cur, step_start, out_bufs,                   &
                               run_energy_budget, run_water_budget, run_face_budget, slow_ledger, &
                               patch_threads, seed_canopy_air)
      type(polygon_t),      intent(inout) :: poly
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t),   intent(in),    optional :: met_src     !< the run's forcing source (when forcing_on)
      type(met_cursor_t),   intent(inout), optional :: met_cur     !< this polygon's forcing cursor
      type(meds_time_t),    intent(in),    optional :: step_start  !< calendar time at the start of this slow step
      type(output_buffers_t), intent(inout), optional :: out_bufs  !< this polygon's FAST-tier staging (to the fast loop)
      type(budget_t), intent(inout), optional :: run_energy_budget, run_water_budget !< run-level ledgers (forwarded)
      type(budget_t), intent(inout), optional :: run_face_budget   !< per-layer face closure, forwarded
      !----- The SLOW tier's own ledger (docs/science/order_of_processes.md §4). Its peers above !
      !      accumulate per-fast-step flux residuals; this one snapshots the site store across   !
      !      the slow step, which is the window neither of them can see. Same lifetime, same     !
      !      place in the plumbing.                                                              !
      !      It also keeps the soil-carbon seam's run-worst checks (docs/science/soil_carbon.md §7). !
      type(slow_ledger_t), intent(inout), optional :: slow_ledger
      integer(ik),        intent(in),    optional :: patch_threads !< the fast loop's patch threads (default 1)
      logical,            intent(inout), optional :: seed_canopy_air !< a start's first window (the fast loop's)
      type(slow_ledger_t) :: ledger   !< the caller's slow ledger, or an inactive one that books nothing

      if (present(slow_ledger)) ledger = slow_ledger

      !----- Fast loop: sub-daily biophysics over the state-hub reservoirs. Every optional is handed  !
      !      straight on: the fast loop takes the forcing only when forcing_on and the source, cursor  !
      !      and step_start are all present, and otherwise runs on the reference climate,           !
      !      met_forcing_t's defaults.                                                               !
      if (cfg%fast_biophysics_on) then
         !----- The fast loop reads each cohort's tissue water against the capacity its biomass     !
         !      allows: first shed what the slow step or the boundary left above it, booked by the  !
         !      slow ledger (docs/science/order_of_processes.md §2). ------------------------------!
         call settle_tissue_water(poly, cfg, ledger)
         call fast_dynamics(poly, cfg, met_src=met_src, met_cur=met_cur,                            &
                            step_start=step_start, out_bufs=out_bufs,                             &
                            run_energy_budget=run_energy_budget, run_water_budget=run_water_budget, &
                            run_face_budget=run_face_budget, patch_threads=patch_threads,        &
                            seed_canopy_air=seed_canopy_air)
      end if

      !----- Slow loop: vegetation dynamics (rate assembly + demographic application), gated on the  !
      !      master slow_on freeze, which holds cohort/patch state static while the fast loop still  !
      !      runs (docs/science/order_of_processes.md §2). Leaf phenology is the FIRST step INSIDE   !
      !      vegetation_dynamics and advances whenever step_start is supplied. --------------------!
      if (cfg%slow_on) call advance_slow_dynamics(poly, cfg, ledger, step_start)
      if (present(slow_ledger)) slow_ledger = ledger
   end subroutine advance_one_step

   !---------------------------------------------------------------------------------------!
   ! The calendar boundary after a step whose end (`now`) left the month -- and, with           !
   ! is_new_year, the year: the stand's monthly and yearly restructuring. The caller makes this !
   ! call after the output has read that step, so the ending period's records describe one      !
   ! stand, and what the restructuring does belongs to the period that begins. Gated on the     !
   ! master slow_on freeze, like the slow step.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine advance_boundary(poly, cfg, is_new_month, is_new_year, slow_ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      type(slow_ledger_t), intent(inout), optional :: slow_ledger
      type(slow_ledger_t) :: ledger   !< the caller's slow ledger, or an inactive one that books nothing
      if (.not. cfg%slow_on) return
      if (present(slow_ledger)) ledger = slow_ledger
      call advance_boundary_dynamics(poly, cfg, is_new_month, is_new_year, ledger)
      if (present(slow_ledger)) slow_ledger = ledger
   end subroutine advance_boundary

end module meds_stepper
