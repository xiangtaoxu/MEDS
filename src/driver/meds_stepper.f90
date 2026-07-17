!==========================================================================================!
! meds_stepper -- the master stepper (the ED2 ed_model analogue; seed of an all-process        !
! top-level loop). Lives in src/driver/, the home of the top-level utilities that wire the      !
! process modules together.                                                                    !
!                                                                                          !
! It owns only the CADENCE: each step it is told whether a month/year rolled over, and it       !
! drives the process modules on the appropriate timescale. Today the only process is the        !
! slow-loop VEGETATION DYNAMICS (meds_vegetation_dynamics), which assembles the demographic     !
! rates and applies them through the demography seam. A future master step would, in addition,  !
! call the fast-loop biophysics processes here.                                                 !
!==========================================================================================!
module meds_stepper
   use meds_config,               only : meds_config_t
   use meds_core_interface, only : site_t
   use meds_vegetation_dynamics,  only : vegetation_dynamics
   use meds_phenology_driver,     only : leaf_phenology
   use meds_fast_loop,            only : fast_context_t, run_fast_biophysics
   use meds_time,                 only : meds_time_t, day_of_year
   use meds_forcing_types,        only : met_driver_t
   use meds_output_types,         only : output_manager_t
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
   subroutine advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx, met_drv, step_start, mgr)
      type(site_t),         intent(inout) :: site
      type(meds_config_t),  intent(in)    :: cfg
      logical,              intent(in)    :: is_new_month, is_new_year
      type(fast_context_t), intent(in),    optional :: fast_ctx
      type(met_driver_t),   intent(inout), optional :: met_drv     !< live met reader (when forcing_on)
      type(meds_time_t),    intent(in),    optional :: step_start  !< calendar time at the start of this slow step
      type(output_manager_t), intent(inout), optional :: mgr       !< FAST-tier staging (forwarded to the fast loop)

      !----- Fast loop: sub-daily biophysics over the state-hub reservoirs. When fast biophysics   !
      !      is ON a fast context MUST be supplied: the old `.and. present(fast_ctx)` SILENTLY       !
      !      skipped the loop, leaving gpp_accum=0 so carbon-mode growth ran on zero GPP (BUG1).     !
      !      Fail loud instead of silently wrong. met_drv/step_start are forwarded only when present !
      !      (forcing_on); absent -> run_fast_biophysics runs the constant-forcing MVP.              !
      if (cfg%fast_biophysics_on) then
         if (.not. present(fast_ctx))                                                              &
            error stop 'advance_one_step: fast_biophysics_on=.true. but no fast_context supplied'
         if (present(met_drv) .and. present(step_start)) then
            call run_fast_biophysics(site, fast_ctx, cfg, met_drv=met_drv, step_start=step_start, mgr=mgr)
         else
            call run_fast_biophysics(site, fast_ctx, cfg)
         end if
      end if

      !----- Leaf phenology (opt-in): advance the per-cohort cue memory + directional status one    !
      !      daily step BETWEEN the fast loop (which supplied the daily-mean air temperature) and the !
      !      slow loop (whose carbon leaf-flush gate reads the status). Needs the step-start date for  !
      !      day-of-year; phenology_on requires forcing_on (validated), so step_start is present.      !
      if (cfg%phenology_on) then
         if (.not. present(step_start))                                                             &
            error stop 'advance_one_step: phenology_on=.true. but no step_start supplied (needs forcing_on)'
         call leaf_phenology(site, cfg, day_of_year(step_start))
      end if

      !----- Slow loop: vegetation dynamics (rate assembly + demographic application). -------!
      call vegetation_dynamics(site, cfg, is_new_month, is_new_year)
   end subroutine advance_one_step

end module meds_stepper
