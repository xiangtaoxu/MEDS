!==========================================================================================!
! meds_plant_vital_rates -- the EMPIRICAL (phenomenological) per-cohort demographic rate      !
! laws, as stateless plant kernels. The deliberate contrast with the mechanistic carbon path  !
! (meds_plant_carbon_dynamics): these predict a cohort's rates directly from its size and the  !
! light it sits under, without resolving the underlying carbon/water balance.                  !
!                                                                                          !
! Each public kernel is a pure function of (one cohort's structure, its local environment, PFT !
! traits) -- no site_t, no cohort/patch state. The stand-scale context each cohort needs (its   !
! overtopping LAI) and the reduction of the per-cohort reproduction flux into a per-(PFT,patch)  !
! recruitment array are supplied by the caller (the vegetation-dynamics driver), which owns the  !
! height-sorted sweep over the cohort structure. That split keeps these laws stateless and       !
! Python-callable, exactly like the leaf / hydraulics / carbon kernels.                          !
!                                                                                          !
!   growth_rate_empirical  -- intrinsic (capped log-linear in dbh) x competition (overtopping    !
!                             LAI) x reproductive-allocation suppression.          [cm/yr]        !
!   mortality_rate         -- Camac et al. (2018) additive hazard on the effective growth.        !
!                             (We stay in the Camac framework.)                     [1/yr]        !
!   recruitment_rate       -- reproduction flux from one cohort: the carbon diverted from growth   !
!                             to reproduction, converted to minimum-size recruits.  [plant/m2/yr]  !
!   min_cohort_carbon      -- carbon of one minimum-size recruit (the recruit "unit"), per PFT.    !
!==========================================================================================!
module meds_plant_vital_rates
   use meds_kinds,     only : wp
   use meds_allometry, only : dbh_to_height, dbh_to_agb, height_to_dbh
   implicit none
   private

   public :: growth_rate_empirical, mortality_rate, recruitment_rate, min_cohort_carbon

   !----- A freshly created cohort carries a negative growth_avg sentinel until its first      !
   !       growth step seeds the running mean; until then mortality uses instantaneous growth. !
   real(wp), parameter :: growth_avg_unset = 0.0_wp   ! threshold: growth_avg < 0 means "unset"

contains

   !=======================================================================================!
   !  Public composite rate laws (one cohort -> one rate).                                  !
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

   !=======================================================================================!
   !  Private per-individual ingredients (moved verbatim from the former meds_growth /       !
   !  meds_mortality / meds_recruitment; combined above into the three public laws).         !
   !=======================================================================================!

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

end module meds_plant_vital_rates
