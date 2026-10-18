! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography_rates -- the per-INDIVIDUAL vital-rate LAWS, as pure scalar functions:      !
! a carbon-driven diameter growth rate, the Camac (2018) additive mortality hazard, and the    !
! reproduction-carbon -> recruit conversion. Growth, survival and fecundity: the three rates    !
! demography is made of, which is why they live in this folder and not with the physiology.     !
!                                                                                          !
! THE RULE THIS FOLDER KEEPS (src/README.md rule 4, "laws and operators do not touch"): the      !
! rate LAWS are here, and the OPERATORS beside them (meds_demography_update, cohort/patch        !
! fuse-fiss) take rate ARRAYS as arguments and never `use` this module. The engine does not      !
! compute a rate -- it APPLIES one.                                                              !
! The slow driver is the only place a rate meets its application. Python `apply_rates`, which    !
! feeds externally computed rates through those same operators, is the standing test of it.      !
!                                                                                          !
! `elemental pure`, scalars + PFT traits only (NO site_t, NO cohort SoA), over the shared        !
! allometry -- so this folder needs nothing from the plant kernels. The EMPIRICAL growth and     !
! recruitment laws are not here; they live in the Python example. Camac mortality stays here.    !
!==========================================================================================!
module meds_demography_rates
   use meds_kinds,     only : wp
   use meds_constants, only : tiny_num
   use meds_allometry, only : allometry_t, wood_to_dbh
   implicit none
   private

   public :: npp_to_growth, camac_mortality, decline_to_floor_mortality, npp_to_recruitment

contains

   !----- Prospective carbon diameter-growth rate [cm/yr]: the dbh the cohort WOULD reach if    !
   !       this step's wood NPP were added to its size-anchor wood_carbon (floored >= 0), minus   !
   !       its current dbh, over the step. Feeds both the applied growth and the mortality hazard.!
   elemental pure function npp_to_growth(allom, wood_carbon, npp_wood, dbh, dt_yr,              &
                                         rho, hgt_max, aboveground_frac) result(dbh_rate)
      type(allometry_t), intent(in) :: allom
      real(wp), intent(in) :: wood_carbon, npp_wood, dbh, dt_yr, rho, hgt_max, aboveground_frac
      real(wp)             :: dbh_rate, dbh_new
      dbh_new  = wood_to_dbh(allom, max(wood_carbon + npp_wood, 0.0_wp), rho, hgt_max, aboveground_frac)
      dbh_rate = (dbh_new - dbh) / dt_yr
   end function npp_to_growth

   !----- Camac et al. (2018) additive mortality hazard [1/yr] on the cohort's EFFECTIVE growth:  !
   !       its tracked running-mean growth (growth_avg) once it has a sample, else the            !
   !       instantaneous rate. `avg_is_set` is the caller's test (the cohort's growth_count > 0), !
   !       kept out of this law, which sees no cohort state.                                      !
   elemental pure function camac_mortality(growth_avg, growth_inst, avg_is_set,                 &
                                           mort_gamma, mort_alpha, mort_beta) result(mort)
      real(wp), intent(in) :: growth_avg, growth_inst, mort_gamma, mort_alpha, mort_beta
      logical,  intent(in) :: avg_is_set
      real(wp)             :: mort, g_eff
      g_eff = merge(growth_avg, growth_inst, avg_is_set)
      mort  = mort_gamma + mort_alpha * exp(-mort_beta * g_eff)
   end function camac_mortality

   !----- The hazard [1/yr] that takes a cohort's density down to the density floor over `window`  !
   !       [yr]: nplant * exp(-mort * window) ends exactly at nplant_floor. It is the death of a     !
   !       cohort with no leaves and no storage (docs/science/plant_carbon_allocation.md §2), and    !
   !       zero for a cohort already at the floor. ---------------------------------------------------!
   elemental pure function decline_to_floor_mortality(nplant, nplant_floor, window) result(mort)
      real(wp), intent(in) :: nplant         !< [plant/m2] the cohort's density
      real(wp), intent(in) :: nplant_floor   !< [plant/m2] the density floor (> 0)
      real(wp), intent(in) :: window         !< [yr] how long the decline takes (> 0)
      real(wp)             :: mort
      mort = log(max(nplant, nplant_floor) / max(nplant_floor, tiny_num)) / max(window, tiny_num)
   end function decline_to_floor_mortality

   !----- One cohort's reproduction-carbon recruit contribution [plant/m2/yr]: the per-plant       !
   !       reproduction NPP over the step, as an annual rate, times the density and the            !
   !       recruit carbon efficiency, over the min-size recruit carbon (min_cohort_carbon, from      !
   !       meds_allometry). Caller reduces these into the (PFT, patch) array and gates by include_pft.!
   elemental pure function npp_to_recruitment(nplant, npp_repro, dt_yr,                         &
                                              recruit_carbon_efficiency, carbon_min) result(rec)
      real(wp), intent(in) :: nplant, npp_repro, dt_yr, recruit_carbon_efficiency, carbon_min
      real(wp)             :: rec
      rec = nplant * (npp_repro / dt_yr) * recruit_carbon_efficiency / carbon_min
   end function npp_to_recruitment

end module meds_demography_rates
