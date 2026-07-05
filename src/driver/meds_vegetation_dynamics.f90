!==========================================================================================!
! meds_vegetation_dynamics -- MEDS's slow-loop VEGETATION-DYNAMICS driver (the analogue of      !
! ED2's veg_dynamics_driver). It lives in src/driver/ (compiled into meds_aux) and drives the   !
! slow loop: it asks the demography RATE provider (meds_demography_rates) for the per-cohort     !
! demographic rates, folds the calendar cadence into the structural triggers, and applies        !
! everything through the demography engine's data seam (update_demography).                       !
!                                                                                          !
! A mechanistic CARBON path will be selected here later (growth_source): it will compute the      !
! per-cohort carbon fluxes via the plant seam and hand the resulting npp to demography, reusing   !
! the same apply-through-update_demography tail. That is the only place plant and demography meet. !
!==========================================================================================!
module meds_vegetation_dynamics
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t, update_demography
   use meds_demography_rates,     only : empirical_vital_rates
   implicit none
   private

   public :: vegetation_dynamics

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation dynamics for one step: assemble the demographic rates for the    !
   ! current stand, fold the calendar cadence + switches into the structural triggers, and    !
   ! apply everything through the demography data seam. (A future master step would call this  !
   ! on the slow cadence alongside the fast-loop biophysics processes.)                       !
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(site, cfg, is_new_month, is_new_year)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable :: growth(:), mortality(:), recruitment(:,:)
      logical               :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse

      !----- 1. Assemble demographic rates for the current stand (empirical path for now; a   !
      !         mechanistic carbon path will be selected here later via growth_source). -------!
      call empirical_vital_rates(site, cfg, growth, mortality, recruitment)

      !----- 2. Fold the calendar cadence + demography on/off + fiss/fuse switches into the   !
      !         structural triggers the engine consumes. -------------------------------------!
      do_cohort_fissfuse   = is_new_month .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_disturbance
      do_patch_fissfuse    = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_fissfuse

      !----- 3. Apply the rates through the demography engine's data seam. --------------------!
      call update_demography(site, growth, mortality, recruitment, cfg, cfg%dt_years,         &
                             do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse)
   end subroutine vegetation_dynamics

end module meds_vegetation_dynamics
