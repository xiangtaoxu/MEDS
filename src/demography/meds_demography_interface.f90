!==========================================================================================!
! meds_demography_interface -- the ONE public seam between the demographic engine and the   !
! rest of MEDS.                                                                            !
!                                                                                          !
! Anything outside src/demography/ drives demography ONLY through this module: it holds and  !
! passes the `site_t` state object and calls `update_demography`, supplying the demographic    !
! rates as plain DATA (three arrays: growth, mortality, recruitment) plus the timestep        !
! `dt_yr` and two structural triggers. The engine never computes a rate; the master stepper   !
! (or a rate module) does, and passes the result in. The engine knows NO calendar and NO      !
! vegetation-dynamics switch: the master stepper decides the cadence and whether structure    !
! changes, and folds both into the `do_cohort_fissfuse` / `do_patch_fissfuse` triggers.        !
!                                                                                          !
! Units: growth [cm/yr] per cohort; mortality [1/yr] per cohort (total); recruitment          !
! [plant/m2/month] per (PFT, patch); dt_yr [yr]. The growth/mortality arrays are indexed in   !
! the current cohort order (1:site_t%cohort%n) as seen at the start of the step.                  !
!==========================================================================================!
module meds_demography_interface
   use meds_kinds,             only : wp, ik
   use meds_config,            only : meds_config_t
   use meds_demography_types,  only : site_t
   use meds_demography_dynamics,  only : growth_step, mortality_step, apply_patch_disturbance,  &
                                         apply_recruitment
   use meds_demography_structure, only : new_fuse_cohorts, terminate_cohorts, split_cohorts,   &
                                         new_fuse_patches, terminate_patches,                   &
                                         sort_cohorts, sort_patches
   implicit none
   private

   !----- The public demography API: the state type and the single entry point. -----------!
   public :: site_t
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
   !   do_cohort_fissfuse   run cohort recruitment + fusion/split this step                  !
   !   do_patch_disturbance open treefall gaps this step                                     !
   !   do_patch_fissfuse    run patch fusion/termination this step                           !
   !                                                                                        !
   ! Growth, mortality and patch ageing advance EVERY step (length dt_yr). The structural    !
   ! triggers are decided by the master stepper (cadence + vegetation-dynamics switch folded  !
   ! in); this routine knows no calendar.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine update_demography(site, growth, mortality, recruitment, cfg, dt_yr,          &
                                do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse)
      type(site_t),          intent(inout) :: site
      real(wp),            intent(in)    :: growth(:)
      real(wp),            intent(in)    :: mortality(:)
      real(wp),            intent(in)    :: recruitment(:,:)
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      logical,             intent(in)    :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse
      integer(ik) :: ip

      !====================================================================================!
      ! PER-STEP application: growth + mortality + patch ageing (always). The kernels take  !
      ! bare arrays (offloaded via OpenMP target); the site_t is unpacked here on the host.   !
      !====================================================================================!
      associate (cohort => site%cohort)
         call growth_step(cohort%n, cohort%dbh, cohort%height, cohort%basal_area, cohort%agb,    &
                          cohort%leaf_area, cohort%p_dbh_critical, cohort%p_wood_density,         &
                          growth, dt_yr)
         call mortality_step(cohort%n, cohort%nplant, mortality, dt_yr, cfg%negligible_nplant)
      end associate
      do ip = 1_ik, site%patch%n
         site%patch%age(ip) = site%patch%age(ip) + dt_yr
      end do

      !====================================================================================!
      ! Cohort restructuring (the stepper decides WHEN via the trigger).                   !
      !====================================================================================!
      if (do_cohort_fissfuse) then
         call apply_recruitment(site, cfg, recruitment)
         call new_fuse_cohorts(site, cfg)
         call terminate_cohorts(site, cfg)
         call split_cohorts(site, cfg)
         call sort_cohorts(site)
      end if

      !====================================================================================!
      ! Patch disturbance (open new age-0 gaps) and patch restructuring are gated by         !
      ! INDEPENDENT triggers; disturbance runs first so the new gap is consolidated by the   !
      ! ensuing patch + cohort fusion (last, as ED2).                                        !
      !====================================================================================!
      if (do_patch_disturbance) then
         call apply_patch_disturbance(site, cfg, patch_dynamics_interval)
      end if

      if (do_patch_fissfuse) then
         call sort_patches(site)
         call new_fuse_patches(site, cfg)
         call terminate_patches(site, cfg)
         call new_fuse_cohorts(site, cfg)
         call terminate_cohorts(site, cfg)
         call sort_cohorts(site)
      end if
   end subroutine update_demography

end module meds_demography_interface
