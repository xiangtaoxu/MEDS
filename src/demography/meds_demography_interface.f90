!==========================================================================================!
! meds_demography_interface -- the ONE public seam between the demographic engine and the   !
! rest of MEDS.                                                                            !
!                                                                                          !
! Anything outside src/demography/ drives demography ONLY through this module: it holds and  !
! passes the `site` state object and calls `update_demography`, supplying the demographic    !
! rates as plain DATA (three arrays: growth, mortality, recruitment) plus the timestep        !
! `dt_yr` and two structural triggers. The engine never computes a rate; the master stepper   !
! (or a rate module) does, and passes the result in. The engine knows NO calendar and NO      !
! vegetation-dynamics switch: the master stepper decides the cadence and whether structure    !
! changes, and folds both into the `do_cohort_fissfuse` / `do_patch_fissfuse` triggers.        !
!                                                                                          !
! Units: growth [cm/yr] per cohort; mortality [1/yr] per cohort (total); recruitment          !
! [plant/m2/month] per (PFT, patch); dt_yr [yr]. The growth/mortality arrays are indexed in   !
! the current cohort order (1:site%coh%n) as seen at the start of the step.                  !
!==========================================================================================!
module meds_demography_interface
   use meds_kinds,             only : wp, ik
   use meds_config,            only : meds_config_t
   use meds_demography_types,  only : site
   use meds_demography_kernels,only : growth_step, mortality_step
   use meds_recruitment,       only : recruitment_month
   use meds_cohort_dynamics,   only : new_fuse_cohorts, terminate_cohorts, split_cohorts
   use meds_patch_dynamics,    only : new_fuse_patches, terminate_patches
   use meds_disturbance,       only : apply_patch_disturbance
   use meds_sort,              only : sort_cohorts, sort_patches
   implicit none
   private

   !----- The public demography API: the state type and the single entry point. -----------!
   public :: site
   public :: update_demography

   !----- Patch structural dynamics (fusion/fission/disturbance) are an ANNUAL process: the !
   !       stepper fires do_patch_fissfuse once per year. Disturbance therefore integrates    !
   !       its yearly rate over this interval, NOT over the per-step dt_yr (which drives the   !
   !       every-step growth/mortality/ageing).                                               !
   real(wp), parameter :: patch_dynamics_interval = 1.0_wp   !< [yr]

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the demographic state by one step from supplied rates.                         !
   !                                                                                        !
   !   growth(:)        per-cohort DBH growth rate   [cm/yr]                                 !
   !   mortality(:)     per-cohort total mortality   [1/yr]                                  !
   !   recruitment(:,:) per-(PFT, patch) recruit density [plant/m2/month]                    !
   !   cfg              tunables (fusion/fission, floors, PFT traits)                        !
   !   dt_yr            timestep length [yr]                                                 !
   !   do_cohort_fissfuse  run cohort recruitment + fusion/split this step                   !
   !   do_patch_fissfuse   run patch fusion/termination this step                            !
   !                                                                                        !
   ! Growth, mortality and patch ageing advance EVERY step (length dt_yr). The structural    !
   ! triggers are decided by the master stepper (cadence + vegetation-dynamics switch folded  !
   ! in); this routine knows no calendar.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine update_demography(asite, growth, mortality, recruitment, cfg, dt_yr,          &
                                do_cohort_fissfuse, do_patch_fissfuse)
      type(site),          intent(inout) :: asite
      real(wp),            intent(in)    :: growth(:)
      real(wp),            intent(in)    :: mortality(:)
      real(wp),            intent(in)    :: recruitment(:,:)
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      logical,             intent(in)    :: do_cohort_fissfuse, do_patch_fissfuse
      integer(ik) :: ip

      !====================================================================================!
      ! PER-STEP application: growth + mortality + patch ageing (always). The kernels take  !
      ! bare arrays (offloaded via OpenMP target); the site is unpacked here on the host.   !
      !====================================================================================!
      associate (c => asite%coh)
         call growth_step(c%n, c%dbh, c%height, c%basal_area, c%agb, c%larea, c%p_dbh_crit,   &
                          c%p_wood_density, growth, dt_yr)
         call mortality_step(c%n, c%nplant, mortality, dt_yr, cfg%negligible_nplant)
      end associate
      do ip = 1_ik, asite%pat%n
         asite%pat%age(ip) = asite%pat%age(ip) + dt_yr
      end do

      !====================================================================================!
      ! Cohort restructuring (the stepper decides WHEN via the trigger).                   !
      !====================================================================================!
      if (do_cohort_fissfuse) then
         call recruitment_month(asite, cfg, recruitment)
         call new_fuse_cohorts(asite, cfg)
         call terminate_cohorts(asite, cfg)
         call split_cohorts(asite, cfg)
         call sort_cohorts(asite)
      end if

      !====================================================================================!
      ! Patch disturbance (open new age-0 gaps), then patch restructuring, then consolidate !
      ! cohorts in the merged patches (last, as ED2).                                       !
      !====================================================================================!
      if (do_patch_fissfuse) then
         call apply_patch_disturbance(asite, cfg, patch_dynamics_interval)
         call sort_patches(asite)
         call new_fuse_patches(asite, cfg)
         call terminate_patches(asite, cfg)
         call new_fuse_cohorts(asite, cfg)
         call terminate_cohorts(asite, cfg)
         call sort_cohorts(asite)
      end if
   end subroutine update_demography

end module meds_demography_interface
