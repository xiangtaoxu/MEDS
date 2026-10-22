! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_trait_dynamics -- within-lifetime CHANGES to a plant's leaf traits.            !
!                                                                                          !
! LIGHT-driven plasticity of the four cohort-level leaf traits SLA, Vcmax25, Rd25 and the leaf   !
! turnover rate, porting ED2's update_cohort_plastic_trait (Lloyd et al. 2010 canopy trait       !
! gradients). Each trait acclimates to the cumulative LAI ABOVE the cohort (its light-           !
! competition proxy) toward a target trait_toc*exp(kplastic*cum_lai_above), and the live cohort  !
! trait RELAXES toward that target only as fast as leaves are REPLACED (turnover-limited: a      !
! leaf's trait is fixed once it flushes). Thermal acclimation is site-level: a site's leaf table  !
! follows its growth temperature; per-cohort thermal acclimation is open (see docs/ROADMAP.md).  !
!                                                                                          !
!   * light_plastic_traits  -- the light-acclimated TARGET traits (pure geometry of the gradient). !
!   * update_plastic_trait  -- one replacement-weighted step of a live trait toward its target      !
!       (or an INSTANTANEOUS jump for a census restart, where no trait history exists).             !
!   * acclimate_leaf_photo_table -- a leaf table's entropy terms and Jmax:Vcmax ratio at a growth  !
!       temperature (Kattge & Knorr 2007).                                                          !
!                                                                                          !
! Stateless, elemental, scalar arithmetic (GPU/SoA-safe); the orchestration + per-cohort state   !
! live in the driver (meds_vegetation_dynamics) and the cohort SoA. docs/science/plant_traits.md !
! eqs. (1)-(2).                                                                                  !
!==========================================================================================!
module meds_plant_trait_dynamics
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : safe_exp
   use meds_pft_params, only : leaf_photo_table_t
   use meds_therm_lib,  only : kattge_knorr_entropy, kattge_knorr_jv_ratio
   implicit none
   private

   public :: light_plastic_traits, update_plastic_trait
   public :: acclimate_leaf_photo_table

   !----- Numerical guards on the gradient exponent (avoid exp overflow); NOT tunable science. --!
   real(wp), parameter :: LNEXP_MIN = -30.0_wp, LNEXP_MAX = 30.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! Light-acclimated TARGET traits for one cohort: each top-of-canopy trait scaled by the     !
   ! exponential canopy-light gradient exp(clamp(kplastic * cum_lai_above)). kplastic < 0        !
   ! (Vcmax, Rd, leaf turnover rate) => the trait falls in shade; kplastic > 0 (SLA) => it rises. !
   ! cum_lai is floored at 0 (a top-of-canopy cohort sees no shade => targets equal the _toc     !
   ! values).                                                                                    !
   ! docs/science/plant_traits.md eq. (1).                                                      !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine light_plastic_traits(cum_lai_above, sla_toc, vcmax25_toc, rd25_toc,  &
            turnover_toc, kplastic_sla, kplastic_vcmax, kplastic_rd, kplastic_turnover,              &
            sla_target, vcmax25_target, rd25_target, turnover_target)
      real(wp), intent(in)  :: cum_lai_above                                 !< [m2/m2] LAI above the cohort
      real(wp), intent(in)  :: sla_toc, vcmax25_toc, rd25_toc, turnover_toc  !< top-of-canopy trait values
      real(wp), intent(in)  :: kplastic_sla, kplastic_vcmax, kplastic_rd, kplastic_turnover  !< light-response slopes
      real(wp), intent(out) :: sla_target, vcmax25_target, rd25_target, turnover_target
      sla_target      = sla_toc      * light_gradient(kplastic_sla,      cum_lai_above)
      vcmax25_target  = vcmax25_toc  * light_gradient(kplastic_vcmax,      cum_lai_above)
      rd25_target     = rd25_toc     * light_gradient(kplastic_rd,       cum_lai_above)
      turnover_target = turnover_toc * light_gradient(kplastic_turnover, cum_lai_above)
   end subroutine light_plastic_traits

   !----- The exponential canopy-light trait gradient, guarded against exp overflow. -------------!
   elemental pure function light_gradient(kplastic, cum_lai_above) result(f)
      real(wp), intent(in) :: kplastic, cum_lai_above
      real(wp)             :: f
      f = safe_exp(max(LNEXP_MIN, min(LNEXP_MAX, kplastic * max(cum_lai_above, 0.0_wp))))
   end function light_gradient

   !---------------------------------------------------------------------------------------!
   ! One update of a live trait toward its light target. Normally REPLACEMENT-WEIGHTED: only      !
   ! newly-flushed leaves carry the new trait, so the cohort-mean moves toward the target by the   !
   ! fraction of leaves replaced this step, f = 1 - exp(-turnover_rate dt), using the cohort's     !
   ! CURRENT leaf turnover rate (a slow-turnover / shaded canopy acclimates slowly). When          !
   ! `instant` (a CENSUS restart, which carries no trait history) the trait jumps straight to the  !
   ! target. turnover_rate is per year and dt in years. docs/science/plant_traits.md eq. (2).      !
   !---------------------------------------------------------------------------------------!
   elemental pure function update_plastic_trait(current, target, turnover_rate, dt_yr, instant) result(updated)
      real(wp), intent(in) :: current, target, turnover_rate, dt_yr
      logical,  intent(in) :: instant
      real(wp)             :: updated, f
      if (instant) then
         updated = target
      else
         f       = 1.0_wp - safe_exp(-turnover_rate * dt_yr)
         updated = current + f * (target - current)
      end if
   end function update_plastic_trait

   !=======================================================================================!
   !  A leaf table's peaked-form entropy terms and Jmax:Vcmax ratio at a growth            !
   !  temperature (Kattge & Knorr 2007; docs/science/leaf_gas_exchange.md §2, Thermal      !
   !  acclimation, eq. (4)), dS = a - b T and Jmax:Vcmax = a - b T. Every other parameter   !
   !  stays the PFT table's, so applying it again at the same temperature changes nothing. !
   !=======================================================================================!
   pure subroutine acclimate_leaf_photo_table(t_growth, ds_vcmax_a, ds_vcmax_b, ds_jmax_a, ds_jmax_b,   &
                                              jv_a, jv_b, table)
      real(wp),                 intent(in)    :: t_growth     !< [K] running-mean growth temperature
      real(wp),                 intent(in)    :: ds_vcmax_a, ds_vcmax_b, ds_jmax_a, ds_jmax_b, jv_a, jv_b
      type(leaf_photo_table_t), intent(inout) :: table
      integer(ik) :: ipft
      if (t_growth <= 0.0_wp) return                          ! not yet seeded: leave the fixed values
      do ipft = 1_ik, table%n_pft
         table%pft(ipft)%ds_vcmax = kattge_knorr_entropy(ds_vcmax_a, ds_vcmax_b, t_growth)
         table%pft(ipft)%ds_jmax  = kattge_knorr_entropy(ds_jmax_a,  ds_jmax_b,  t_growth)
         table%jmax_vcmax_ratio(ipft) = kattge_knorr_jv_ratio(jv_a, jv_b, t_growth)
      end do
   end subroutine acclimate_leaf_photo_table

end module meds_plant_trait_dynamics
