! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography -- the demographic laws (docs/science/demography.md §1), as elemental    !
! pure functions of scalars and PFT traits: the prospective diameter growth rate           !
! (eq. (1)), the Camac et al. (2018) mortality hazard (eq. (3)), the decline to the        !
! density floor (eq. (4)) and a cohort's term of the recruitment rate (eq. (6)).           !
! meds_slow_update, which applies the rates, never uses this module (src/README.md rule    !
! 4). The empirical growth and recruitment laws live in the Python example.                !
!==========================================================================================!
module meds_demography
   use meds_kinds,         only : wp
   use meds_constants,     only : tiny_num
   use meds_allometry_lib, only : allometry_t, wood_to_dbh
   implicit none
   private

   public :: npp_to_growth, camac_mortality, decline_to_floor_mortality, recruit_efficiency, npp_to_recruitment

contains

   !----- The prospective diameter growth rate [cm/yr], eq. (1): the diameter the cohort would !
   !      reach if the step's wood growth were added to its wood carbon (at least 0), less its !
   !      current diameter, over the step. The mortality hazard reads it (eq. (3)). -----------!
   elemental pure function npp_to_growth(allom, wood_carbon, npp_wood, dbh, dt_yr, rho) result(dbh_rate)
      type(allometry_t), intent(in) :: allom   !< the cohort's PFT's allometry
      real(wp), intent(in) :: wood_carbon, npp_wood, dbh, dt_yr, rho
      real(wp)             :: dbh_rate, dbh_new
      dbh_new  = wood_to_dbh(allom, max(wood_carbon + npp_wood, 0.0_wp), rho)
      dbh_rate = (dbh_new - dbh) / dt_yr
   end function npp_to_growth

   !----- The Camac et al. (2018) hazard [1/yr], eq. (3), on the growth average once the cohort !
   !      has a sample (avg_is_set, the caller's test), else on the prospective rate. ----------!
   elemental pure function camac_mortality(growth_avg, growth_inst, avg_is_set,                 &
                                           mort_gamma, mort_alpha, mort_beta) result(mort)
      real(wp), intent(in) :: growth_avg, growth_inst, mort_gamma, mort_alpha, mort_beta
      logical,  intent(in) :: avg_is_set
      real(wp)             :: mort, g_eff
      g_eff = merge(growth_avg, growth_inst, avg_is_set)
      mort  = mort_gamma + mort_alpha * exp(-mort_beta * g_eff)
   end function camac_mortality

   !----- The hazard [1/yr] that takes a cohort's density to the floor over `window`, eq. (4): !
   !      the death of a cohort with no leaves and no storage (plant_carbon_allocation.md §2), !
   !      and zero for a cohort already at the floor. -----------------------------------------!
   elemental pure function decline_to_floor_mortality(nplant, nplant_floor, window) result(mort)
      real(wp), intent(in) :: nplant         !< [plant/m2] the cohort's density
      real(wp), intent(in) :: nplant_floor   !< [plant/m2] the density floor (> 0)
      real(wp), intent(in) :: window         !< [yr] how long the decline takes (> 0)
      real(wp)             :: mort
      mort = log(max(nplant, nplant_floor) / max(nplant_floor, tiny_num)) / max(window, tiny_num)
   end function decline_to_floor_mortality

   !----- The share of a cohort's reproduction carbon that becomes recruits, eq. (5): the     !
   !      PFT's efficiency in the open, reduced by the shade of its patch's LAI. ------------!
   elemental pure function recruit_efficiency(efficiency_open, shade_decay, patch_lai) result(eff)
      real(wp), intent(in) :: efficiency_open   !< [-] the efficiency in the open
      real(wp), intent(in) :: shade_decay       !< [m2/m2 per LAI] k of exp(-k LAI)
      real(wp), intent(in) :: patch_lai         !< [m2/m2] the patch's LAI
      real(wp)             :: eff
      eff = efficiency_open * exp(-shade_decay * patch_lai)
   end function recruit_efficiency

   !----- One cohort's term of the recruitment rate [plant/m2/yr], eq. (6): its reproduction  !
   !      carbon per plant as a yearly rate, times the density and the recruit efficiency     !
   !      (eq. (5)), over carbon_min. The caller sums the terms per PFT and patch and applies !
   !      include_pft. -----------------------------------------------------------------------!
   elemental pure function npp_to_recruitment(nplant, npp_repro, dt_yr,                         &
                                              recruit_carbon_efficiency, carbon_min) result(rec)
      real(wp), intent(in) :: nplant, npp_repro, dt_yr, recruit_carbon_efficiency, carbon_min
      real(wp)             :: rec
      rec = nplant * (npp_repro / dt_yr) * recruit_carbon_efficiency / carbon_min
   end function npp_to_recruitment

end module meds_demography
