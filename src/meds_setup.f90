!==========================================================================================!
! meds_setup -- helpers to build initial communities for the demo and the test suite.      !
!==========================================================================================!
module meds_setup
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : pio4
   use meds_pft_params, only : dbh_to_height
   use meds_config,     only : meds_config_t
   use meds_demography_types,      only : site, site_alloc, cohort_ensure_capacity, rebuild_csr
   use meds_sort,       only : sort_cohorts
   implicit none
   private

   public :: init_bare_ground, add_cohort, finalize_init

contains

   !---------------------------------------------------------------------------------------!
   ! Near-bare-ground site: n_patch identical empty patches sharing the site area.     !
   !---------------------------------------------------------------------------------------!
   subroutine init_bare_ground(comm, cfg, n_patch, avg_temp, min_temp)
      type(site),     intent(out) :: comm
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch
      real(wp),            intent(in)  :: avg_temp, min_temp
      integer(ik) :: ip
      call site_alloc(comm, cfg%pft%n, coh_cap = 64_ik, pat_cap = max(n_patch, 1_ik))
      comm%pat%n = n_patch
      do ip = 1_ik, n_patch
         comm%pat%area(ip)           = 1.0_wp / real(n_patch, wp)
         comm%pat%age(ip)            = 0.0_wp
         comm%pat%dist_type(ip)      = 1_ik
         comm%pat%avg_daily_temp(ip) = avg_temp
         comm%pat%min_month_temp(ip) = min_temp
      end do
      comm%pat%recruit_pool = 0.0_wp
      comm%coh%n = 0_ik
      call rebuild_csr(comm)
   end subroutine init_bare_ground

   !---------------------------------------------------------------------------------------!
   ! Append one cohort to patch ip (caller invokes finalize_init afterwards).               !
   !---------------------------------------------------------------------------------------!
   subroutine add_cohort(comm, cfg, ip, pft, nplant, dbh)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: ip, pft
      real(wp),            intent(in)    :: nplant, dbh
      integer(ik) :: m
      call cohort_ensure_capacity(comm%coh, comm%coh%n + 1_ik)
      m = comm%coh%n + 1_ik
      associate (c => comm%coh, t => cfg%pft)
         c%pft(m)            = pft
         c%owner_patch(m)    = ip
         c%nplant(m)         = nplant
         c%dbh(m)            = dbh
         c%basal_area(m)        = pio4 * dbh * dbh
         c%height(m)           = dbh_to_height(t%height_min(pft), t%b1_height(pft), t%b2_height(pft), dbh)
         c%p_dbh_crit(m)     = t%dbh_crit(pft)
         c%p_height_min(m)      = t%height_min(pft)
         c%p_b1_height(m)         = t%b1_height(pft)
         c%p_b2_height(m)         = t%b2_height(pft)
      end associate
      comm%coh%n = m
   end subroutine add_cohort

   subroutine finalize_init(comm)
      type(site), intent(inout) :: comm
      call rebuild_csr(comm)
      call sort_cohorts(comm)
   end subroutine finalize_init

end module meds_setup
