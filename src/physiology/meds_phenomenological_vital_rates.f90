!==========================================================================================!
! meds_phenomenological_vital_rates -- the structure-only vital-rate provider (physiology layer).!
!                                                                                          !
! Produces the three demographic-rate arrays update_demography consumes -- per-cohort GROWTH   !
! and MORTALITY, and per-(PFT,patch) RECRUITMENT -- from the CURRENT demographic structure     !
! alone (size, the light a cohort sits under, and its recent growth). The relationships are    !
! PHENOMENOLOGICAL: they describe how rates vary with structure without resolving the           !
! underlying mechanism. This is the deliberate contrast with the planned mechanistic provider,  !
! which will compute the same arrays from individual carbon/water balance. Because the seam is  !
! plain data, a mechanistic module drops in by exposing the same `vital_rates` routine.         !
!                                                                                          !
! GROWTH = growth_intrinsic * suppression_competition * suppression_reproduction  [cm/yr]       !
!   growth_intrinsic    = growth_dbh_max * exp(growth_dbh_slope * max(0, ln(cap) - ln(dbh)))    !
!                         (a capped log-linear function of dbh; = growth_dbh_max for dbh>=cap)  !
!   suppression_comp    = exp(growth_lai_slope * overtopping_LAI)            in (0, 1]          !
!   suppression_repro   = 1 below min_reproduction_height, else 1 - reproduction_investment_frac!
!   Equal-height cohorts are co-dominant: they share the overtopping LAI and do not shade each  !
!   other, so the result is independent of their arbitrary within-layer order.                  !
!                                                                                          !
! MORTALITY (Camac et al. 2018, PNAS; additive hazard) = mort_gamma + mort_alpha * exp(         !
!   -mort_beta * growth_avg)  [1/yr], where growth_avg is the cohort's tracked running-mean      !
!   growth rate (a proxy for carbon budget); a freshly created cohort uses its instantaneous     !
!   growth until the average is seeded. The three parameters are wood-density-derived.           !
!                                                                                          !
! RECRUITMENT = seed_rain_recruits (baseline external rain) + a reproduction flux: the carbon    !
!   diverted from growth to reproduction (growth_intrinsic * suppression_comp *                  !
!   reproduction_investment_fraction, in carbon units via allometry) summed over reproductive    !
!   cohorts (height >= min_reproduction_height) of each PFT in the patch, times                  !
!   repro_carbon_efficiency, divided by the carbon of one minimum-size (min_cohort_height)       !
!   cohort.  Units [plant/m2/month] (the carbon flux is taken over one month).                   !
!==========================================================================================!
module meds_phenomenological_vital_rates
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mon_per_yr
   use meds_allometry,        only : dbh_to_agb, dbh_to_height, height_to_dbh
   use meds_config,           only : meds_config_t
   use meds_demography_types, only : site_t
   implicit none
   private

   public :: vital_rates

contains

   !---------------------------------------------------------------------------------------!
   ! Evaluate all three phenomenological vital rates for the current site_t.                 !
   !---------------------------------------------------------------------------------------!
   subroutine vital_rates(site, cfg, growth, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr]  per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]   per cohort, total
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)

      call growth_mortality(site, cfg, growth, mortality, recruitment)
   end subroutine vital_rates

   !---------------------------------------------------------------------------------------!
   ! The single height-sorted, per-patch overtopping-LAI sweep that produces growth and       !
   ! mortality and, in the same pass, accumulates the reproduction part of recruitment.       !
   !---------------------------------------------------------------------------------------!
   subroutine growth_mortality(site, cfg, growth, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:), mortality(:), recruitment(:,:)
      integer(ik) :: n, np, npft, i, j, k, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, layer_lai, gi, supp_comp, supp_repro, g_eff
      real(wp)    :: dbh_min, repro_dbh, new_dbh, dagb
      real(wp), allocatable :: carbon_min(:)     ! [kgC/plant] carbon of one min-size cohort, per PFT

      n    = site%cohort%n
      np   = site%patch%n
      npft = cfg%pft%n
      allocate(growth(n), mortality(n), recruitment(npft, max(np, 1_ik)))
      recruitment = 0.0_wp

      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         !----- Carbon of one minimum-size cohort (the recruit "unit"), per PFT. -----------!
         allocate(carbon_min(npft))
         dbh_min = height_to_dbh(pft%min_cohort_height)
         do pf = 1_ik, npft
            carbon_min(pf) = dbh_to_agb(dbh_min, pft%min_cohort_height, pft%wood_density(pf))
         end do

         !----- Baseline external seed rain (independent of structure). --------------------!
         do ip = 1_ik, np
            do pf = 1_ik, npft
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%seed_rain_recruits(pf)
            end do
         end do

         !----- Per-patch top-down sweep over the overtopping LAI. -------------------------!
         do ip = 1_ik, np
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping LAI accumulator
            i   = i0
            do while (i <= i1)
               !----- Equal-height cohorts are co-dominant (share the overtopping LAI). -----!
               k = i
               do while (k < i1)
                  if (cohort%height(k + 1_ik) /= cohort%height(i)) exit
                  k = k + 1_ik
               end do
               over_lai  = cum
               layer_lai = 0.0_wp
               do j = i, k
                  pf = cohort%pft(j)
                  !----- Growth: intrinsic (capped log-linear in dbh) x competition x repro. !
                  gi = pft%growth_dbh_max(pf) * exp(pft%growth_dbh_slope(pf)                  &
                       * max(0.0_wp, log(pft%growth_dbh_cap(pf)) - log(cohort%dbh(j))))
                  supp_comp = exp(pft%growth_lai_slope(pf) * over_lai)
                  if (cohort%height(j) >= pft%min_reproduction_height) then
                     supp_repro = 1.0_wp - pft%reproduction_investment_fraction(pf)
                  else
                     supp_repro = 1.0_wp
                  end if
                  growth(j) = gi * supp_comp * supp_repro

                  !----- Mortality (Camac additive hazard) from the tracked running-mean      !
                  !      growth; a not-yet-seeded cohort uses its instantaneous growth.        !
                  if (cohort%growth_avg(j) < 0.0_wp) then
                     g_eff = growth(j)
                  else
                     g_eff = cohort%growth_avg(j)
                  end if
                  mortality(j) = pft%mort_gamma(pf)                                          &
                                 + pft%mort_alpha(pf) * exp(-pft%mort_beta(pf) * g_eff)

                  !----- Reproduction flux: carbon of the growth diverted to reproduction      !
                  !      (over one month), converted to min-size recruits.                     !
                  if (cohort%height(j) >= pft%min_reproduction_height) then
                     repro_dbh = gi * supp_comp * pft%reproduction_investment_fraction(pf)   ! [cm/yr]
                     new_dbh   = cohort%dbh(j) + repro_dbh / mon_per_yr                       ! one month's worth
                     dagb      = dbh_to_agb(new_dbh, dbh_to_height(new_dbh), cohort%p_wood_density(j)) &
                                 - cohort%agb(j)                                              ! [kgC/plant/month]
                     recruitment(pf, ip) = recruitment(pf, ip)                               &
                          + cohort%nplant(j) * dagb * pft%repro_carbon_efficiency(pf) / carbon_min(pf)
                  end if

                  layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades those below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine growth_mortality

end module meds_phenomenological_vital_rates
