!==========================================================================================!
! meds_setup -- helpers to build initial communities for the demo and the test suite.      !
!==========================================================================================!
module meds_setup
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, DIST_PRIMARY
   use meds_demography_types,      only : site_t, site_alloc, cohort_ensure_capacity, rebuild_csr,  &
                                          set_cohort_size
   use meds_demography_structure, only : sort_cohorts
   implicit none
   private

   public :: init_bare_ground, add_cohort, finalize_init

contains

   !---------------------------------------------------------------------------------------!
   ! Near-bare-ground site_t: n_patch identical empty patches sharing the site_t area.     !
   !---------------------------------------------------------------------------------------!
   subroutine init_bare_ground(site, cfg, n_patch, avg_temp, min_temp)
      type(site_t),     intent(out) :: site
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch
      real(wp),            intent(in)  :: avg_temp, min_temp
      integer(ik) :: ip
      call site_alloc(site, cfg%pft%n, coh_cap = 64_ik, pat_cap = max(n_patch, 1_ik))
      site%patch%n = n_patch
      do ip = 1_ik, n_patch
         site%patch%area(ip)           = 1.0_wp / real(n_patch, wp)
         site%patch%age(ip)            = 0.0_wp
         site%patch%dist_type(ip)      = DIST_PRIMARY
         site%patch%avg_daily_temp(ip) = avg_temp
         site%patch%min_month_temp(ip) = min_temp
      end do
      site%patch%recruit_pool = 0.0_wp
      site%cohort%n = 0_ik
      call rebuild_csr(site)
   end subroutine init_bare_ground

   !---------------------------------------------------------------------------------------!
   ! Append one cohort to patch ip (caller invokes finalize_init afterwards).               !
   !---------------------------------------------------------------------------------------!
   subroutine add_cohort(site, cfg, ip, ipft, nplant, dbh)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: ip, ipft
      real(wp),            intent(in)    :: nplant, dbh
      integer(ik) :: m
      call cohort_ensure_capacity(site%cohort, site%cohort%n + 1_ik)
      m = site%cohort%n + 1_ik
      associate (cohort => site%cohort, pft => cfg%pft)
         cohort%pft(m)            = ipft
         cohort%owner_patch(m)    = ip
         cohort%nplant(m)         = nplant
         cohort%dbh(m)            = dbh
         cohort%p_dbh_critical(m) = pft%dbh_critical(ipft)
         cohort%p_wood_density(m) = pft%wood_density(ipft)
         call set_cohort_size(cohort, m)            ! height/basal_area/agb/leaf_area from dbh
      end associate
      site%cohort%n = m
   end subroutine add_cohort

   subroutine finalize_init(site)
      type(site_t), intent(inout) :: site
      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine finalize_init

end module meds_setup
