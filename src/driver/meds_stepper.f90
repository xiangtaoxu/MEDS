!==========================================================================================!
! meds_stepper -- the master stepper (OUTSIDE the demography engine; seed of an all-process  !
! top-level loop, analogous to ED2's ed_model.F90 calling veg_dynamics_driver). Lives in      !
! src/driver/, the home of top-level utilities that wire the process modules together.        !
!                                                                                          !
! It drives demography ONLY through meds_demography_interface: each step it produces the     !
! demographic rate arrays (here from the empirical evaluator) and hands them, plus the        !
! timestep dt_yr and the structural triggers, to update_demography. The stepper OWNS the      !
! calendar and the vegetation-dynamics switch: it folds the cadence flags + cfg%vegetation_   !
! dynamics_on into the fission/fusion triggers it passes. It touches no demography internals   !
! (kernels/fusion/sort/types) -- only the public seam + the rate module.                      !
!==========================================================================================!
module meds_stepper
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site, update_demography
   use meds_rates_empirical,      only : empirical_vital_rates
   implicit none
   private

   public :: advance_one_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one step: evaluate the empirical vital rates for the current site, then apply   !
   ! them via the demography interface. The caller's calendar supplies the cadence flags;     !
   ! here we fold them, together with the vegetation-dynamics switch and the config           !
   ! fission/fusion switches, into the structural triggers handed to the engine. (A future    !
   ! master step would, in addition, call other process modules and select the rate model.)  !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(comm, cfg, is_new_month, is_new_year)
      type(site),          intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable :: growth(:), mortality(:), recruitment(:,:)
      logical               :: do_cohort_fissfuse, do_patch_fissfuse

      call empirical_vital_rates(comm, cfg, growth, mortality, recruitment)

      !----- Fold cadence + vegetation-dynamics switch into the structural triggers. -------!
      do_cohort_fissfuse = is_new_month .and. cfg%vegetation_dynamics_on .and.              &
                           cfg%do_cohort_fissfuse
      do_patch_fissfuse  = is_new_year  .and. cfg%vegetation_dynamics_on .and.              &
                           cfg%do_patch_fissfuse

      call update_demography(comm, growth, mortality, recruitment, cfg, cfg%dt_years,       &
                             do_cohort_fissfuse, do_patch_fissfuse)
   end subroutine advance_one_step

end module meds_stepper
