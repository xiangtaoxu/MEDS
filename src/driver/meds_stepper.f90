!==========================================================================================!
! meds_stepper -- the master stepper (OUTSIDE the demography engine; seed of an all-process  !
! top-level loop, analogous to ED2's ed_model.F90 calling veg_dynamics_driver). Lives in      !
! src/driver/, the home of top-level utilities that wire the process modules together.        !
!                                                                                          !
! It drives demography ONLY through meds_demography_interface: each step it produces the     !
! demographic rate arrays (here from the phenomenological rate provider) and hands them, plus  !
! the timestep dt_yr and the structural triggers, to update_demography. The stepper OWNS the   !
! calendar and the demography on/off switch: it folds the cadence flags + cfg%demography_on    !
! into the fission/fusion triggers it passes. It touches no demography internals               !
! (kernels/fusion/sort/types) -- only the public seam + the rate module.                      !
!==========================================================================================!
module meds_stepper
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t, update_demography
   use meds_phenomenological_vital_rates, only : vital_rates
   implicit none
   private

   public :: advance_one_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one step: evaluate the phenomenological vital rates for the current site_t, then   !
   ! apply them via the demography interface. The caller's calendar supplies the cadence flags;!
   ! here we fold them, together with the demography on/off switch and the config              !
   ! fission/fusion switches, into the structural triggers handed to the engine. (A future    !
   ! master step would, in addition, call other process modules and select the rate model.)  !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(site, cfg, is_new_month, is_new_year)
      type(site_t),          intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable :: growth(:), mortality(:), recruitment(:,:)
      logical               :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse

      call vital_rates(site, cfg, growth, mortality, recruitment)

      !----- Fold cadence + demography on/off switch into the structural triggers. ---------!
      do_cohort_fissfuse  = is_new_month .and. cfg%demography_on .and.                       &
                            cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year .and. cfg%demography_on .and.                       &
                            cfg%do_patch_disturbance
      do_patch_fissfuse   = is_new_year  .and. cfg%demography_on .and.                       &
                            cfg%do_patch_fissfuse

      call update_demography(site, growth, mortality, recruitment, cfg, cfg%dt_years,       &
                             do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse)
   end subroutine advance_one_step

end module meds_stepper
