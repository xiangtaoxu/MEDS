!==========================================================================================!
! meds_step -- the demographic stepper: one canonical order of operations per step.        !
!                                                                                          !
! The engine knows nothing about HOW rates are computed -- it only calls the provider. The   !
! cadence flags (is_new_month, is_new_year) are supplied by the caller's calendar so the     !
! same routine serves daily (default) and monthly time-stepping.                            !
!                                                                                          !
! Order (matches the synthesised ED2-faithful sequence):                                    !
!   DAILY    competition -> growth -> mortality accumulation                                !
!   MONTHLY  patch age   -> apply mortality -> recruitment -> cohort fuse/terminate/split    !
!   ANNUAL   patch fuse  -> patch terminate -> cohort fuse/terminate (consolidate merges)    !
!==========================================================================================!
module meds_step
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mon_per_yr
   use meds_config,           only : meds_config_t
   use meds_demography_types,            only : community
   use meds_demography_interface,   only : rate_provider_t
   use meds_demography_kernels,          only : compute_competition, growth_step, mortality_accumulate,&
                                     mortality_apply_month
   use meds_recruitment,      only : recruitment_month
   use meds_cohort_dynamics,  only : new_fuse_cohorts, terminate_cohorts, split_cohorts
   use meds_patch_dynamics,   only : new_fuse_patches, terminate_patches
   use meds_sort,             only : sort_cohorts, sort_patches
   implicit none
   private

   public :: advance_one_step

contains

   subroutine advance_one_step(comm, provider, cfg, is_new_month, is_new_year)
      type(community),        intent(inout) :: comm
      class(rate_provider_t), intent(in)    :: provider
      type(meds_config_t),    intent(in)    :: cfg
      logical,                intent(in)    :: is_new_month, is_new_year
      integer(ik) :: ip

      !====================================================================================!
      ! DAILY (or per-step) biophysical-rate-driven demography.                            !
      !====================================================================================!
      call compute_competition(comm)
      call growth_step(comm, provider, cfg%dt_years)
      call mortality_accumulate(comm, provider, cfg%dt_years)

      !====================================================================================!
      ! MONTHLY structural update.                                                         !
      !====================================================================================!
      if (is_new_month) then
         do ip = 1_ik, comm%pat%n
            comm%pat%age(ip) = comm%pat%age(ip) + 1.0_wp / mon_per_yr
         end do
         call mortality_apply_month(comm, cfg%negligible_nplant)

         if (cfg%veget_dyn_on) then
            call recruitment_month(comm, cfg, provider)
            call new_fuse_cohorts(comm, cfg)
            call terminate_cohorts(comm, cfg)
            call split_cohorts(comm, cfg)
            call sort_cohorts(comm)
         end if
      end if

      !====================================================================================!
      ! ANNUAL patch dynamics (last, as in ED2).                                           !
      !====================================================================================!
      if (is_new_year .and. cfg%veget_dyn_on) then
         call sort_patches(comm)              ! oldest-first => older patch is the fusion receptor
         call new_fuse_patches(comm, cfg)
         call terminate_patches(comm, cfg)
         !----- Consolidate cohorts in the (now larger) merged patches. -------------------!
         call new_fuse_cohorts(comm, cfg)
         call terminate_cohorts(comm, cfg)
         call sort_cohorts(comm)
      end if
   end subroutine advance_one_step

end module meds_step
