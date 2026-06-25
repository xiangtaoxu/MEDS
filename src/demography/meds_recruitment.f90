!==========================================================================================!
! meds_recruitment -- monthly recruitment of new cohorts from a supplied recruit density.   !
!                                                                                          !
! Each month the per-(PFT,patch) recruit density [plant m-2 month-1] supplied from outside   !
! is accumulated into a carry-forward pool. When a pool reaches `min_recruit_size` a new      !
! cohort is spawned at the PFT's minimum diameter and the pool is reset; otherwise it carries !
! over (so rare recruiters still establish eventually). Recruitment is a HOST operation: it   !
! changes the number of cohorts.                                                            !
!==========================================================================================!
module meds_recruitment
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : pio4
   use meds_pft_params,     only : dbh_to_height
   use meds_config,         only : meds_config_t
   use meds_demography_types,  only : site, cohort_ensure_capacity, rebuild_csr
   use meds_sort,           only : sort_cohorts
   implicit none
   private

   public :: recruitment_month

contains

   subroutine recruitment_month(comm, cfg, recruitment)
      type(site),          intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)
      integer(ik) :: ip, pf, np, m, nspawn

      np = comm%pat%n
      if (np < 1_ik) return

      !----- Accumulate the supplied monthly recruit density into the carry-forward pool. -!
      do ip = 1_ik, np
         do pf = 1_ik, comm%n_pft
            comm%pat%recruit_pool(pf, ip) = comm%pat%recruit_pool(pf, ip) + recruitment(pf, ip)
         end do
      end do

      !----- Count pools that have reached the spawn threshold. ---------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, comm%n_pft
            if (comm%pat%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      if (nspawn == 0_ik) return

      call cohort_ensure_capacity(comm%coh, comm%coh%n + nspawn)
      m = comm%coh%n
      associate (c => comm%coh, p => comm%pat, t => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, comm%n_pft
               if (p%recruit_pool(pf, ip) < cfg%min_recruit_size) cycle
               m = m + 1_ik
               c%pft(m)            = pf
               c%owner_patch(m)    = ip
               c%nplant(m)         = p%recruit_pool(pf, ip)
               c%dbh(m)            = t%dbh_min(pf)
               c%basal_area(m)        = pio4 * c%dbh(m) * c%dbh(m)
               c%height(m)           = dbh_to_height(t%height_min(pf), t%b1_height(pf), t%b2_height(pf), c%dbh(m))
               c%p_dbh_crit(m)     = t%dbh_crit(pf)
               c%p_height_min(m)      = t%height_min(pf)
               c%p_b1_height(m)         = t%b1_height(pf)
               c%p_b2_height(m)         = t%b2_height(pf)
               p%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         c%n = m
      end associate

      call rebuild_csr(comm)
      call sort_cohorts(comm)
   end subroutine recruitment_month

end module meds_recruitment
