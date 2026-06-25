!==========================================================================================!
! meds_recruitment -- monthly recruitment of new cohorts from a supplied recruit density.   !
!                                                                                          !
! Each month the per-(PFT,patch) recruit density [plant m-2 month-1] supplied from outside   !
! is accumulated into a carry-forward pool. When a pool reaches `min_recruit_size` a new      !
! cohort is spawned at the SHARED minimum reproduction height (identical for every PFT) and   !
! the pool is reset; otherwise it carries over (so rare recruiters still establish            !
! eventually). Recruitment is a HOST operation: it changes the number of cohorts.            !
!==========================================================================================!
module meds_recruitment
   use meds_kinds,          only : wp, ik
   use meds_allometry,      only : height_to_dbh
   use meds_config,         only : meds_config_t
   use meds_demography_types,  only : site_t, cohort_ensure_capacity, rebuild_csr, set_cohort_size, &
                                      assign_cohort_id
   use meds_demography_structure, only : sort_cohorts
   implicit none
   private

   public :: recruitment_month

contains

   subroutine recruitment_month(site, cfg, recruitment)
      type(site_t),          intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh

      np = site%patch%n
      if (np < 1_ik) return

      !----- All PFTs recruit at the same height -> the same diameter. --------------------!
      recruit_dbh = height_to_dbh(cfg%pft%min_reproduction_height)

      !----- Accumulate the supplied monthly recruit density into the carry-forward pool. -!
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            site%patch%recruit_pool(pf, ip) = site%patch%recruit_pool(pf, ip) + recruitment(pf, ip)
         end do
      end do

      !----- Count pools that have reached the spawn threshold. ---------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            if (site%patch%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      if (nspawn == 0_ik) return

      call cohort_ensure_capacity(site%cohort, site%cohort%n + nspawn)
      n_before = site%cohort%n
      m = site%cohort%n
      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, site%n_pft
               if (patch%recruit_pool(pf, ip) < cfg%min_recruit_size) cycle
               m = m + 1_ik
               cohort%pft(m)            = pf
               cohort%owner_patch(m)    = ip
               cohort%nplant(m)         = patch%recruit_pool(pf, ip)
               cohort%dbh(m)            = recruit_dbh
               cohort%p_dbh_critical(m)     = pft%dbh_critical(pf)
               cohort%p_wood_density(m) = pft%wood_density(pf)
               call set_cohort_size(cohort, m)         ! height/basal_area/agb/leaf_area from dbh
               patch%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         cohort%n = m
      end associate

      !----- Stamp each freshly spawned cohort with a persistent global id. ----------------!
      do k = n_before + 1_ik, site%cohort%n
         call assign_cohort_id(site, k)
      end do

      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine recruitment_month

end module meds_recruitment
