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
   use meds_demography_types,  only : site, cohort_ensure_capacity, rebuild_csr, set_cohort_size
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
      real(wp)    :: recruit_dbh

      np = comm%pat%n
      if (np < 1_ik) return

      !----- All PFTs recruit at the same height -> the same diameter. --------------------!
      recruit_dbh = height_to_dbh(cfg%pft%min_reproduction_height)

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
               c%dbh(m)            = recruit_dbh
               c%p_dbh_crit(m)     = t%dbh_crit(pf)
               c%p_wood_density(m) = t%wood_density(pf)
               call set_cohort_size(c, m)         ! height/basal_area/agb/leaf_area from dbh
               p%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         c%n = m
      end associate

      call rebuild_csr(comm)
      call sort_cohorts(comm)
   end subroutine recruitment_month

end module meds_recruitment
