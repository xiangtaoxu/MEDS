! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_stepper -- one coupled slow step of one polygon: the fast loop, then the slow       !
! dynamics (docs/science/order_of_processes.md §2-§4); and the calendar boundary between   !
! two steps (§5), which the caller makes after the output has read the step that ended on  !
! it. A library of its own, so the tests can drive the coupled step without a run.         !
!==========================================================================================!
module meds_stepper
   use meds_kinds,              only : ik
   use meds_config,             only : meds_config_t
   use meds_time,               only : meds_time_t
   use meds_budget_check,       only : budget_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_slow_dynamics,      only : advance_slow_dynamics, advance_boundary_dynamics, settle_tissue_water
   use meds_slow_ledger,        only : slow_ledger_t
   use meds_fast_dynamics,      only : fast_dynamics
   use meds_met_types,          only : met_source_t, met_cursor_t
   use meds_output_types,       only : output_buffers_t
   implicit none
   private

   public :: advance_one_step, advance_boundary

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one slow step. With fast biophysics on, the sub-daily fast loop runs before   !
   ! the slow dynamics, so the vegetation dynamics allocates the GPP the fast loop         !
   ! accumulated over the step. Every argument after cfg is optional, so a test can step a !
   ! bare polygon; the run passes all.                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(poly, cfg, met_src, met_cur, step_start, out_bufs,                   &
                               run_energy_budget, run_water_budget, run_face_budget, slow_ledger, &
                               patch_threads, seed_canopy_air)
      type(polygon_t),        intent(inout) :: poly
      type(meds_config_t),    intent(in)    :: cfg
      type(met_source_t),     intent(in),    optional :: met_src      !< the run's forcing source
      type(met_cursor_t),     intent(inout), optional :: met_cur      !< this polygon's forcing cursor
      type(meds_time_t),      intent(in),    optional :: step_start   !< the step's start date
      type(output_buffers_t), intent(inout), optional :: out_bufs     !< this polygon's FAST-tier staging
      type(budget_t),         intent(inout), optional :: run_energy_budget, run_water_budget !< whole-column
      type(budget_t),         intent(inout), optional :: run_face_budget  !< each layer face's water
      !< the site's stores across the step, and the soil-carbon seam's worst values
      type(slow_ledger_t),    intent(inout), optional :: slow_ledger
      integer(ik),            intent(in),    optional :: patch_threads    !< the patch loop's threads (default 1)
      logical,                intent(inout), optional :: seed_canopy_air  !< a start's first window
      type(slow_ledger_t) :: ledger   !< the caller's slow ledger, or an inactive one that books nothing

      if (present(slow_ledger)) ledger = slow_ledger

      !----- The fast loop. The tissue water above what each cohort's biomass now holds is shed    !
      !      first, and the slow ledger books it (docs/science/order_of_processes.md §2). The loop !
      !      takes the forcing when forcing_on and the source, cursor and step_start are all       !
      !      given, and otherwise the reference climate. ------------------------------------------!
      if (cfg%fast_biophysics_on) then
         call settle_tissue_water(poly, cfg, ledger)
         call fast_dynamics(poly, cfg, met_src=met_src, met_cur=met_cur,                            &
                            step_start=step_start, out_bufs=out_bufs,                             &
                            run_energy_budget=run_energy_budget, run_water_budget=run_water_budget, &
                            run_face_budget=run_face_budget, patch_threads=patch_threads,        &
                            seed_canopy_air=seed_canopy_air)
      end if

      !----- The slow dynamics, unless the master slow_on freeze holds the stand while the fast    !
      !      loop runs: leaf phenology (when step_start is given), the vegetation dynamics and the !
      !      soil carbon. ---------------------------------------------------------------------------!
      if (cfg%slow_on) call advance_slow_dynamics(poly, cfg, ledger, step_start)
      if (present(slow_ledger)) slow_ledger = ledger
   end subroutine advance_one_step

   !---------------------------------------------------------------------------------------!
   ! The calendar boundary after a step whose end left the month (and, with is_new_year,   !
   ! the year): the stand's monthly and yearly restructuring, unless slow_on freezes the   !
   ! stand. The caller makes it after the output has read that step, so the ending         !
   ! period's records describe one stand, and what the restructuring does belongs to the   !
   ! period that begins.                                                                   !
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
