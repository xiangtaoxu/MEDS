!==========================================================================================!
! meds_plant_trait_dynamics -- within-lifetime CHANGES to a plant's leaf traits.            !
!                                                                                          !
! First cut: LIGHT-driven plasticity of the four cohort-level leaf traits SLA, Vcmax25, Rd25,  !
! and leaf lifespan (llspan), porting ED2's update_cohort_plastic_trait (Lloyd et al. 2010     !
! canopy trait gradients). Each trait acclimates to the cumulative LAI ABOVE the cohort (its    !
! light-competition proxy) toward a target trait_toc*exp(kplastic*cum_lai_above), and the live   !
! cohort trait RELAXES toward that target only as fast as leaves are REPLACED (turnover-limited: !
! a leaf's trait is fixed once it flushes). Thermal acclimation is a deferred sibling.           !
!                                                                                          !
!   * light_plastic_traits -- the light-acclimated TARGET traits (pure geometry of the gradient).!
!   * relax_trait          -- one replacement-weighted step of a live trait toward its target.    !
!                                                                                          !
! Stateless, elemental, scalar arithmetic (GPU/SoA-safe); the orchestration + per-cohort state   !
! live in the driver (meds_vegetation_dynamics) and core SoA. See                                !
! docs/dev_plans/MEDS_PLANT_TRAIT_DYNAMICS_DESIGN.md.                                             !
!==========================================================================================!
module meds_plant_trait_dynamics
   use meds_kinds,     only : wp
   use meds_constants, only : safe_exp
   implicit none
   private

   public :: light_plastic_traits, relax_trait

   !----- Numerical guards on the gradient exponent (avoid exp overflow); NOT tunable science. --!
   real(wp), parameter :: LNEXP_MIN = -30.0_wp, LNEXP_MAX = 30.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! Light-acclimated TARGET traits for one cohort: each top-of-canopy trait scaled by the     !
   ! exponential canopy-light gradient exp(clamp(kplastic * cum_lai_above)). kplastic < 0        !
   ! (Vcmax/Rd) => the trait falls in shade; kplastic > 0 (SLA, llspan) => it rises. cum_lai is  !
   ! floored at 0 (a top-of-canopy cohort sees no shade => targets equal the _toc values).       !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine light_plastic_traits(cum_lai_above, sla_toc, vcmax25_toc, rd25_toc,  &
            llspan_toc, kplastic_sla, kplastic_vm0, kplastic_rd, kplastic_llspan,                  &
            sla_target, vcmax25_target, rd25_target, llspan_target)
      real(wp), intent(in)  :: cum_lai_above                                 !< [m2/m2] LAI above the cohort
      real(wp), intent(in)  :: sla_toc, vcmax25_toc, rd25_toc, llspan_toc    !< top-of-canopy trait values
      real(wp), intent(in)  :: kplastic_sla, kplastic_vm0, kplastic_rd, kplastic_llspan  !< light-response slopes
      real(wp), intent(out) :: sla_target, vcmax25_target, rd25_target, llspan_target
      sla_target     = sla_toc     * light_gradient(kplastic_sla,    cum_lai_above)
      vcmax25_target = vcmax25_toc * light_gradient(kplastic_vm0,    cum_lai_above)
      rd25_target    = rd25_toc    * light_gradient(kplastic_rd,     cum_lai_above)
      llspan_target  = llspan_toc  * light_gradient(kplastic_llspan, cum_lai_above)
   end subroutine light_plastic_traits

   !----- The exponential canopy-light trait gradient, guarded against exp overflow. -------------!
   elemental pure function light_gradient(kplastic, cum_lai_above) result(f)
      real(wp), intent(in) :: kplastic, cum_lai_above
      real(wp)             :: f
      f = safe_exp(max(LNEXP_MIN, min(LNEXP_MAX, kplastic * max(cum_lai_above, 0.0_wp))))
   end function light_gradient

   !---------------------------------------------------------------------------------------!
   ! One REPLACEMENT-WEIGHTED relaxation step of a live trait toward its light target. Only      !
   ! newly-flushed leaves carry the new trait, so the cohort-mean trait moves toward the target   !
   ! by the fraction of leaves replaced this step, f = 1 - exp(-dt/llspan), using the cohort's     !
   ! CURRENT leaf lifespan. A long-lived (shaded) canopy therefore acclimates slowly. dt and      !
   ! llspan are both in years.                                                                    !
   !---------------------------------------------------------------------------------------!
   elemental pure function relax_trait(current, target, llspan, dt_yr) result(updated)
      real(wp), intent(in) :: current, target, llspan, dt_yr
      real(wp)             :: updated, f
      f       = 1.0_wp - safe_exp(-dt_yr / max(llspan, tiny(1.0_wp)))
      updated = current + f * (target - current)
   end function relax_trait

end module meds_plant_trait_dynamics
