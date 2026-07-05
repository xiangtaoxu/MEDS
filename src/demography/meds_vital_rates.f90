!==========================================================================================!
! meds_phenomenological_vital_rates -- the structure-only vital-rate provider (ORCHESTRATOR). !
!                                                                                          !
! Produces the three demographic-rate arrays update_demography consumes -- per-cohort GROWTH   !
! and MORTALITY, and per-(PFT,patch) RECRUITMENT -- from the CURRENT demographic structure     !
! alone (size, the light a cohort sits under, and its recent growth). The relationships are    !
! PHENOMENOLOGICAL: they describe how rates vary with structure without resolving the           !
! underlying mechanism. This is the deliberate contrast with the planned mechanistic provider,  !
! which will compute the same arrays from individual carbon/water balance. Because the seam is  !
! plain data, a mechanistic module drops in by exposing the same `vital_rates` routine.         !
!                                                                                          !
! The per-individual FORMULAS live in sibling modules (meds_growth, meds_mortality,             !
! meds_recruitment); THIS routine owns the single, height-sorted, per-patch overtopping-LAI     !
! sweep that supplies each cohort its competition context and combines the three -- growth      !
! (intrinsic x competition x reproductive allocation), the Camac growth-dependent mortality     !
! hazard, and the reproduction recruit flux -- in one pass, plus the baseline seed rain.        !
!==========================================================================================!
module meds_phenomenological_vital_rates
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_demography_types, only : site_t
   use meds_growth,           only : growth_intrinsic, competition_suppression, reproduction_suppression
   use meds_mortality,        only : effective_growth, mortality_hazard
   use meds_recruitment,      only : min_cohort_carbon, reproduction_recruits
   implicit none
   private

   public :: vital_rates

contains

   !---------------------------------------------------------------------------------------!
   ! Evaluate all three phenomenological vital rates for the current site_t in ONE height-    !
   ! sorted, per-patch overtopping-LAI sweep: growth and mortality per cohort, and -- in the   !
   ! same pass -- the reproduction part of recruitment.                                        !
   !---------------------------------------------------------------------------------------!
   subroutine vital_rates(site, cfg, growth, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr]  per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]   per cohort, total
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: n, np, npft, i, j, k, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, layer_lai, gi, supp_comp, supp_repro, g_eff, repro_dbh
      real(wp), allocatable :: carbon_min(:)     ! [kgC/plant] carbon of one min-size cohort, per PFT

      n    = site%cohort%n
      np   = site%patch%n
      npft = cfg%pft%n
      allocate(growth(n), mortality(n), recruitment(npft, max(np, 1_ik)))
      recruitment = 0.0_wp

      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         !----- Carbon of one minimum-size cohort (the recruit "unit"), per PFT. -----------!
         allocate(carbon_min(npft))
         do pf = 1_ik, npft
            carbon_min(pf) = min_cohort_carbon(pft%min_cohort_height, pft%wood_density(pf))
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
                  gi         = growth_intrinsic(cohort%dbh(j), pft%growth_dbh_max(pf),          &
                                                pft%growth_dbh_slope(pf), pft%growth_dbh_cap(pf))
                  supp_comp  = competition_suppression(over_lai, pft%growth_lai_slope(pf))
                  supp_repro = reproduction_suppression(cohort%height(j), pft%min_reproduction_height, &
                                                        pft%reproduction_investment_fraction(pf))
                  growth(j) = gi * supp_comp * supp_repro

                  !----- Mortality (Camac additive hazard) from the tracked running-mean      !
                  !      growth; a not-yet-seeded cohort uses its instantaneous growth.        !
                  g_eff        = effective_growth(cohort%growth_avg(j), growth(j))
                  mortality(j) = mortality_hazard(g_eff, pft%mort_gamma(pf), pft%mort_alpha(pf), &
                                                  pft%mort_beta(pf))

                  !----- Reproduction flux: carbon of the growth diverted to reproduction      !
                  !      (over one year), converted to min-size recruits.                      !
                  if (cohort%height(j) >= pft%min_reproduction_height) then
                     repro_dbh = gi * supp_comp * pft%reproduction_investment_fraction(pf)   ! [cm/yr]
                     recruitment(pf, ip) = recruitment(pf, ip)                                    &
                          + reproduction_recruits(cohort%dbh(j), cohort%p_hgt_max(j), cohort%agb(j), &
                                                  cohort%nplant(j), cohort%p_wood_density(j),        &
                                                  repro_dbh, pft%repro_carbon_efficiency(pf),        &
                                                  carbon_min(pf))
                  end if

                  layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades those below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine vital_rates

end module meds_phenomenological_vital_rates
