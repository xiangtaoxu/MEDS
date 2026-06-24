!==========================================================================================!
! meds_step -- the master stepper (OUTSIDE the demography engine; seed of an all-process    !
! top-level loop, analogous to ED2's ed_model.F90 calling veg_dynamics_driver).             !
!                                                                                          !
! It drives demography ONLY through meds_demography_interface: each step it produces the     !
! demographic rate arrays (here from the empirical evaluator) and hands them, plus the        !
! cadence + fission/fusion flags from config, to update_demography. It touches no demography  !
! internals (kernels/fusion/sort/types) -- only the public seam + the rate module.           !
!==========================================================================================!
module meds_step
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
   ! them via the demography interface. Cadence flags come from the caller's calendar; the    !
   ! fission/fusion switches come from config. (A future master step would, in addition,      !
   ! call other process modules and select the rate model.)                                  !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(comm, cfg, is_new_month, is_new_year)
      type(site),          intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable :: growth(:), mortality(:), recruitment(:,:)

      call empirical_vital_rates(comm, cfg, growth, mortality, recruitment)
      call update_demography(comm, growth, mortality, recruitment, cfg, is_new_month,       &
                             is_new_year, cfg%do_cohort_fissfuse, cfg%do_patch_fissfuse)
   end subroutine advance_one_step

end module meds_step
