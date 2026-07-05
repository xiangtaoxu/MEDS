!==========================================================================================!
! meds_demography_rates -- the DEMOGRAPHIC rate provider: it turns the current stand into the  !
! three rate arrays update_demography consumes (per-cohort growth [cm/yr] & mortality [1/yr],   !
! per-(PFT,patch) recruitment [plant/m2/yr]).                                                   !
!                                                                                          !
! Organized by DOMAIN: these are phenomenological POPULATION relationships (empirical growth,    !
! Camac mortality, reproduction->recruits), so they live in demography -- NOT in the plant        !
! physiology library. The per-individual LAWS are private here; the ONE height-sorted, per-patch  !
! overtopping-LAI sweep that supplies each cohort its competition context and reduces the           !
! reproduction flux into the (PFT,patch) recruitment array is the public entry point.               !
! Self-contained (links state + allometry only, NO plant dependency).                               !
!                                                                                          !
! A mechanistic CARBON-sourced provider (fed by the plant carbon fluxes) will be added here as a    !
! sibling, selected by growth_source, producing the same rate arrays with no engine change.         !
!==========================================================================================!
module meds_demography_rates
   use meds_kinds,            only : wp, ik
   use meds_allometry,        only : dbh_to_height, dbh_to_agb, height_to_dbh
   use meds_config,           only : meds_config_t
   use meds_demography_types, only : site_t
   implicit none
   private

   public :: empirical_vital_rates

   !----- A freshly created cohort carries a negative growth_avg sentinel until its first      !
   !       growth step seeds the running mean; until then mortality uses instantaneous growth. !
   real(wp), parameter :: growth_avg_unset = 0.0_wp   ! threshold: growth_avg < 0 means "unset"

contains

   !---------------------------------------------------------------------------------------!
   ! Evaluate the EMPIRICAL (phenomenological) vital rates for the current site_t in ONE       !
   ! height-sorted, per-patch overtopping-LAI sweep: growth and mortality per cohort, and --   !
   ! in the same pass -- the reproduction part of recruitment, reduced into the per-(PFT,patch) !
   ! array. Produces the three arrays update_demography consumes; a mechanistic (carbon)        !
   ! provider produces the same arrays with no engine change.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine empirical_vital_rates(site, cfg, growth, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr]  per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]   per cohort, total
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: n, np, npft, i, j, k, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, layer_lai
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
                  !----- Growth: empirical law (intrinsic x competition x repro allocation). --!
                  growth(j) = growth_rate_empirical(cohort%dbh(j), over_lai, cohort%height(j),  &
                                 pft%growth_dbh_max(pf), pft%growth_dbh_slope(pf),              &
                                 pft%growth_dbh_cap(pf), pft%growth_lai_slope(pf),              &
                                 pft%min_reproduction_height, pft%reproduction_investment_fraction(pf))

                  !----- Mortality (Camac additive hazard) from the tracked running-mean       !
                  !      growth; a not-yet-seeded cohort uses its instantaneous growth.        !
                  mortality(j) = mortality_rate(cohort%growth_avg(j), growth(j),               &
                                 pft%mort_gamma(pf), pft%mort_alpha(pf), pft%mort_beta(pf))

                  !----- Reproduction flux -> recruits (zero below the maturity height). -------!
                  recruitment(pf, ip) = recruitment(pf, ip)                                    &
                       + recruitment_rate(cohort%dbh(j), cohort%height(j), over_lai,           &
                                          cohort%agb(j), cohort%nplant(j), cohort%p_hgt_max(j),&
                                          cohort%p_wood_density(j), pft%growth_dbh_max(pf),    &
                                          pft%growth_dbh_slope(pf), pft%growth_dbh_cap(pf),    &
                                          pft%growth_lai_slope(pf), pft%min_reproduction_height,&
                                          pft%reproduction_investment_fraction(pf),            &
                                          pft%repro_carbon_efficiency(pf), carbon_min(pf))

                  layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades those below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine empirical_vital_rates

   !=======================================================================================!
   !  Private per-cohort LAWS (one cohort -> one rate) + their ingredients.                 !
   !=======================================================================================!

   !----- Empirical diameter growth [cm/yr]: intrinsic (capped log-linear in dbh) attenuated  !
   !       by the overtopping LAI and by the fraction of growth diverted to reproduction. -----!
   elemental pure function growth_rate_empirical(dbh, over_lai, height,                       &
              growth_dbh_max, growth_dbh_slope, growth_dbh_cap, growth_lai_slope,             &
              min_reproduction_height, reproduction_investment_fraction) result(growth)
      real(wp), intent(in) :: dbh, over_lai, height
      real(wp), intent(in) :: growth_dbh_max, growth_dbh_slope, growth_dbh_cap, growth_lai_slope
      real(wp), intent(in) :: min_reproduction_height, reproduction_investment_fraction
      real(wp)             :: growth
      growth =  growth_intrinsic(dbh, growth_dbh_max, growth_dbh_slope, growth_dbh_cap)        &
              * competition_suppression(over_lai, growth_lai_slope)                            &
              * reproduction_suppression(height, min_reproduction_height, reproduction_investment_fraction)
   end function growth_rate_empirical

   !----- Camac et al. (2018) additive mortality hazard [1/yr] on the cohort's effective       !
   !       growth (its tracked running mean, or its instantaneous growth before that is seeded).!
   elemental pure function mortality_rate(growth_avg, growth,                                 &
              mort_gamma, mort_alpha, mort_beta) result(mort)
      real(wp), intent(in) :: growth_avg, growth, mort_gamma, mort_alpha, mort_beta
      real(wp)             :: mort
      mort = mortality_hazard(effective_growth(growth_avg, growth), mort_gamma, mort_alpha, mort_beta)
   end function mortality_rate

   !----- Reproduction recruit flux from one cohort [plant/m2/yr]: below the maturity height it  !
   !       is zero; otherwise the diverted growth (repro_dbh, one year) becomes an AGB gain and  !
   !       -- times the establishment efficiency, over the min-size recruit carbon -- recruits.  !
   elemental pure function recruitment_rate(dbh, height, over_lai, agb, nplant, hgt_max,      &
              wood_density, growth_dbh_max, growth_dbh_slope, growth_dbh_cap, growth_lai_slope,&
              min_reproduction_height, reproduction_investment_fraction,                      &
              repro_carbon_efficiency, carbon_min) result(rec)
      real(wp), intent(in) :: dbh, height, over_lai, agb, nplant, hgt_max, wood_density
      real(wp), intent(in) :: growth_dbh_max, growth_dbh_slope, growth_dbh_cap, growth_lai_slope
      real(wp), intent(in) :: min_reproduction_height, reproduction_investment_fraction
      real(wp), intent(in) :: repro_carbon_efficiency, carbon_min
      real(wp)             :: rec, repro_dbh
      if (height >= min_reproduction_height) then
         repro_dbh =  growth_intrinsic(dbh, growth_dbh_max, growth_dbh_slope, growth_dbh_cap)  &
                    * competition_suppression(over_lai, growth_lai_slope)                      &
                    * reproduction_investment_fraction
         rec = reproduction_recruits(dbh, hgt_max, agb, nplant, wood_density,                 &
                                     repro_dbh, repro_carbon_efficiency, carbon_min)
      else
         rec = 0.0_wp
      end if
   end function recruitment_rate

   !----- Carbon of one minimum-size cohort [kgC/plant] (the recruit unit), for a given PFT.  !
   elemental pure function min_cohort_carbon(min_cohort_height, wood_density) result(carbon_min)
      real(wp), intent(in) :: min_cohort_height, wood_density
      real(wp)             :: carbon_min
      carbon_min = dbh_to_agb(height_to_dbh(min_cohort_height), min_cohort_height, wood_density)
   end function min_cohort_carbon

   !----- Intrinsic diameter growth [cm/yr]: capped log-linear in dbh (= growth_dbh_max for  !
   !       dbh >= growth_dbh_cap; larger for smaller stems).                                   !
   elemental pure function growth_intrinsic(dbh, growth_dbh_max, growth_dbh_slope, growth_dbh_cap) result(gi)
      real(wp), intent(in) :: dbh, growth_dbh_max, growth_dbh_slope, growth_dbh_cap
      real(wp)             :: gi
      gi = growth_dbh_max * exp(growth_dbh_slope * max(0.0_wp, log(growth_dbh_cap) - log(dbh)))
   end function growth_intrinsic

   !----- Competition suppression in (0, 1]: Beer-like attenuation by the overtopping LAI.   !
   elemental pure function competition_suppression(over_lai, growth_lai_slope) result(supp)
      real(wp), intent(in) :: over_lai, growth_lai_slope
      real(wp)             :: supp
      supp = exp(growth_lai_slope * over_lai)
   end function competition_suppression

   !----- Reproductive-allocation suppression of growth: 1 below the maturity height, else    !
   !       1 - reproduction_investment_fraction (the growth fraction diverted to reproduction).!
   elemental pure function reproduction_suppression(height, min_reproduction_height,          &
                                                    reproduction_investment_fraction) result(supp)
      real(wp), intent(in) :: height, min_reproduction_height, reproduction_investment_fraction
      real(wp)             :: supp
      if (height >= min_reproduction_height) then
         supp = 1.0_wp - reproduction_investment_fraction
      else
         supp = 1.0_wp
      end if
   end function reproduction_suppression

   !----- Effective growth for the mortality response: the tracked running-mean growth, or the  !
   !       instantaneous growth for a not-yet-seeded cohort (growth_avg < 0 sentinel). ---------!
   elemental pure function effective_growth(growth_avg, growth_instantaneous) result(g_eff)
      real(wp), intent(in) :: growth_avg, growth_instantaneous
      real(wp)             :: g_eff
      if (growth_avg < growth_avg_unset) then
         g_eff = growth_instantaneous
      else
         g_eff = growth_avg
      end if
   end function effective_growth

   !----- Camac et al. (2018) additive mortality hazard [1/yr]. --------------------------------!
   elemental pure function mortality_hazard(g_eff, mort_gamma, mort_alpha, mort_beta) result(mort)
      real(wp), intent(in) :: g_eff, mort_gamma, mort_alpha, mort_beta
      real(wp)             :: mort
      mort = mort_gamma + mort_alpha * exp(-mort_beta * g_eff)
   end function mortality_hazard

   !----- Reproduction recruit flux from one cohort [plant/m2/yr]: grow the diameter by the    !
   !       repro-allocated increment (repro_dbh, one year's worth), take the resulting AGB      !
   !       gain, and convert that reproductive carbon (times the establishment efficiency) into !
   !       minimum-size recruits. Caller has already gated this on the maturity height.  -------!
   elemental pure function reproduction_recruits(dbh, hgt_max, agb, nplant, wood_density,       &
                                repro_dbh, repro_carbon_efficiency, carbon_min) result(rec)
      real(wp), intent(in) :: dbh, hgt_max, agb, nplant, wood_density
      real(wp), intent(in) :: repro_dbh, repro_carbon_efficiency, carbon_min
      real(wp)             :: rec, new_dbh, dagb
      new_dbh = dbh + repro_dbh
      dagb    = dbh_to_agb(new_dbh, dbh_to_height(new_dbh, hgt_max), wood_density) - agb
      rec     = nplant * dagb * repro_carbon_efficiency / carbon_min
   end function reproduction_recruits

end module meds_demography_rates
