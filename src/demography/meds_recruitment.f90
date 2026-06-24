!==========================================================================================!
! meds_recruitment -- monthly recruitment of new cohorts from an external schedule.        !
!                                                                                          !
! Each month the provider's recruit density [plant m-2 month-1] is accumulated into a       !
! per-(patch,PFT) carry-forward pool. When a pool reaches `min_recruit_size` a new cohort    !
! is spawned at the PFT's minimum diameter and the pool is reset; otherwise it carries over  !
! (so rare recruiters still establish eventually). Recruitment is a HOST operation: it       !
! changes the number of cohorts.                                                            !
!==========================================================================================!
module meds_recruitment
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : pio4
   use meds_pft_params,     only : dbh2h
   use meds_config,         only : meds_config_t
   use meds_demography_types,          only : community, cohort_ensure_capacity, rebuild_csr
   use meds_demography_interface, only : rate_provider_t
   use meds_sort,           only : sort_cohorts
   implicit none
   private

   public :: recruitment_month

contains

   subroutine recruitment_month(comm, cfg, provider)
      type(community),        intent(inout) :: comm
      type(meds_config_t),    intent(in)    :: cfg
      class(rate_provider_t), intent(in)    :: provider
      integer(ik) :: ip, pf, np, m, nspawn
      real(wp)    :: rdens

      np = comm%pat%n
      if (np < 1_ik) return

      !----- Accumulate the monthly recruit density into the carry-forward pool. ----------!
      do ip = 1_ik, np
         do pf = 1_ik, comm%n_pft
            rdens = provider%recruit(pf, comm%pat%dist_type(ip), comm%pat%min_month_temp(ip))
            comm%pat%recruit_pool(pf, ip) = comm%pat%recruit_pool(pf, ip) + rdens
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
               c%basarea(m)        = pio4 * c%dbh(m) * c%dbh(m)
               c%hite(m)           = dbh2h(t%hgt_min(pf), t%b1ht(pf), t%b2ht(pf), c%dbh(m))
               c%monthly_dlnndt(m) = 0.0_wp
               c%ddbh_dt(m)        = 0.0_wp
               c%comp(m)           = 0.0_wp
               c%p_dbh_crit(m)     = t%dbh_crit(pf)
               c%p_hgt_min(m)      = t%hgt_min(pf)
               c%p_b1ht(m)         = t%b1ht(pf)
               c%p_b2ht(m)         = t%b2ht(pf)
               p%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         c%n = m
      end associate

      call rebuild_csr(comm)
      call sort_cohorts(comm)
   end subroutine recruitment_month

end module meds_recruitment
