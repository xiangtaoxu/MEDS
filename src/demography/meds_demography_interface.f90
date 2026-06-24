!==========================================================================================!
! meds_demography_interface -- the ONE public seam between the demographic engine and the   !
! rest of MEDS.                                                                            !
!                                                                                          !
! Anything outside src/demography/ drives demography ONLY through this module: it holds and  !
! passes the `site` state object and calls `update_demography`, supplying the demographic    !
! rates as plain DATA (three arrays: growth, mortality, recruitment) plus two fission/fusion  !
! flags. The engine never computes a rate; the master stepper (or a rate module) does, and    !
! passes the result in. The timestep (daily/weekly/monthly) is chosen outside, in config +   !
! the master stepper, and reaches the engine only as `dt_years` (via cfg) + the cadence flags.!
!                                                                                          !
! Units: growth [cm/yr] per cohort; mortality [1/yr] per cohort (total); recruitment          !
! [plant/m2/month] per (PFT, patch). The growth/mortality arrays are indexed in the current   !
! cohort order (1:site%coh%n) as seen at the start of the step.                              !
!==========================================================================================!
module meds_demography_interface
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : mon_per_yr
   use meds_config,            only : meds_config_t
   use meds_demography_types,  only : site
   use meds_demography_kernels,only : growth_step, mortality_accumulate, mortality_apply_month
   use meds_recruitment,       only : recruitment_month
   use meds_cohort_dynamics,   only : new_fuse_cohorts, terminate_cohorts, split_cohorts
   use meds_patch_dynamics,    only : new_fuse_patches, terminate_patches
   use meds_sort,              only : sort_cohorts, sort_patches
   implicit none
   private

   !----- The public demography API: the state type and the single entry point. -----------!
   public :: site
   public :: update_demography

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the demographic state by one step from supplied rates.                         !
   !                                                                                        !
   !   growth(:)        per-cohort DBH growth rate   [cm/yr]                                 !
   !   mortality(:)     per-cohort total mortality   [1/yr]                                  !
   !   recruitment(:,:) per-(PFT, patch) recruit density [plant/m2/month]                    !
   !   cfg              carries dt_years + tunables (set from the timestep mode)             !
   !   is_new_month / is_new_year  cadence flags from the master stepper's calendar          !
   !   do_cohort_fissfuse  run cohort fusion + split this month                              !
   !   do_patch_fissfuse   run patch fusion/termination this year                            !
   !---------------------------------------------------------------------------------------!
   subroutine update_demography(asite, growth, mortality, recruitment, cfg,                 &
                                is_new_month, is_new_year, do_cohort_fissfuse, do_patch_fissfuse)
      type(site),          intent(inout) :: asite
      real(wp),            intent(in)    :: growth(:)
      real(wp),            intent(in)    :: mortality(:)
      real(wp),            intent(in)    :: recruitment(:,:)
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      logical,             intent(in)    :: do_cohort_fissfuse, do_patch_fissfuse
      integer(ik) :: ip

      !====================================================================================!
      ! DAILY (or per-step) application of the supplied rates.                             !
      !====================================================================================!
      call growth_step(asite, growth, cfg%dt_years)
      call mortality_accumulate(asite, mortality, cfg%dt_years)

      !====================================================================================!
      ! MONTHLY structural update.                                                         !
      !====================================================================================!
      if (is_new_month) then
         do ip = 1_ik, asite%pat%n
            asite%pat%age(ip) = asite%pat%age(ip) + 1.0_wp / mon_per_yr
         end do
         call mortality_apply_month(asite, cfg%negligible_nplant)

         if (cfg%veget_dyn_on) then
            call recruitment_month(asite, cfg, recruitment)
            if (do_cohort_fissfuse) call new_fuse_cohorts(asite, cfg)
            call terminate_cohorts(asite, cfg)
            if (do_cohort_fissfuse) call split_cohorts(asite, cfg)
            call sort_cohorts(asite)
         end if
      end if

      !====================================================================================!
      ! ANNUAL patch dynamics (last, as in ED2).                                           !
      !====================================================================================!
      if (is_new_year .and. cfg%veget_dyn_on .and. do_patch_fissfuse) then
         call sort_patches(asite)
         call new_fuse_patches(asite, cfg)
         call terminate_patches(asite, cfg)
         !----- Consolidate cohorts in the (now larger) merged patches. -------------------!
         if (do_cohort_fissfuse) call new_fuse_cohorts(asite, cfg)
         call terminate_cohorts(asite, cfg)
         call sort_cohorts(asite)
      end if
   end subroutine update_demography

end module meds_demography_interface
