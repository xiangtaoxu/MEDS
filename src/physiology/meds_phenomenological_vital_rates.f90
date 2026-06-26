!==========================================================================================!
! meds_phenomenological_vital_rates -- the structure-only vital-rate provider (physiology layer).!
!                                                                                          !
! Produces the three demographic-rate arrays update_demography consumes -- per-cohort GROWTH   !
! and MORTALITY, and per-(PFT,patch) RECRUITMENT -- from the CURRENT demographic structure     !
! alone (size + the light it casts). The relationships are PHENOMENOLOGICAL: they describe how !
! rates vary with structure without resolving the underlying mechanism. This is the deliberate !
! contrast with the planned mechanistic provider, which will compute the same arrays from        !
! individual-level carbon/water balance (photosynthesis, allocation, hydraulics). Because the    !
! seam is plain data, a mechanistic module drops in by exposing the same `vital_rates` routine.  !
!                                                                                          !
! Growth & mortality (per height-sorted patch, top-down overtopping-LAI sweep):                 !
!   light        = exp(-light_ext * overtopping_LAI)                  in (0, 1]                  !
!   growth [cm/yr]   = gr_max  * light * dbh        (realized relative growth x diameter)        !
!   mortality [1/yr] = mort_base + mort_shade * (1 - light)   (baseline + shade-driven)          !
! Cohorts at EQUAL height are co-dominant -- they share the overtopping LAI and do not shade     !
! one another -- so the result is independent of their arbitrary within-layer order. Low-wood-   !
! density PFTs carry a high gr_max AND a high mort_shade: they out-grow others in the light but   !
! die fast once overtopped -- the classic ED growth-mortality trade-off.                         !
!                                                                                          !
! Recruitment: a constant per-PFT density gated only by include_pft (no environmental control) -- !
! the seam where future seed-production / seedling dynamics will live. The engine accumulates the !
! rate and spawns cohorts (meds_demography_dynamics::apply_recruitment).                          !
!==========================================================================================!
module meds_phenomenological_vital_rates
   use meds_kinds,            only : wp, ik
   use meds_allometry,        only : light_ext
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

      call growth_and_mortality(site, cfg, growth, mortality)
      call recruitment_density(site, cfg, recruitment)
   end subroutine vital_rates

   !---------------------------------------------------------------------------------------!
   ! Per-cohort growth & mortality from light competition (overtopping LAI per height-sorted !
   ! patch). Reads the site_t and computes the competition index.                            !
   !---------------------------------------------------------------------------------------!
   subroutine growth_and_mortality(site, cfg, growth, mortality)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:), mortality(:)
      integer(ik) :: n, np, i, j, k, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, light, layer_lai

      n  = site%cohort%n
      np = site%patch%n
      allocate(growth(n), mortality(n))

      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, np
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping LAI accumulator
            i   = i0
            do while (i <= i1)
               !----- Cohorts at the SAME height are CO-DOMINANT: they share the overtopping !
               !      LAI (the LAI strictly above their layer) and do not shade each other,   !
               !      so the light is independent of their (arbitrary) within-layer order.    !
               !      This removes the order/PFT bias for equal-height cohorts (e.g. recruits !
               !      born at the same height, or split daughters whose height is capped).    !
               k = i
               do while (k < i1)
                  if (cohort%height(k + 1_ik) /= cohort%height(i)) exit
                  k = k + 1_ik
               end do
               over_lai  = cum                             ! LAI strictly above this layer
               light     = exp(-light_ext * over_lai)
               layer_lai = 0.0_wp
               do j = i, k
                  pf           = cohort%pft(j)
                  growth(j)    = pft%gr_max(pf) * light * cohort%dbh(j)
                  mortality(j) = pft%mort_base(pf) + pft%mort_shade(pf) * (1.0_wp - light)
                  layer_lai    = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades the layers below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine growth_and_mortality

   !---------------------------------------------------------------------------------------!
   ! Per-(PFT, patch) potential recruitment density [plant m-2 month-1]: a constant per-PFT   !
   ! rate gated only by include_pft (no environmental control).                              !
   !---------------------------------------------------------------------------------------!
   subroutine recruitment_density(site, cfg, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: recruitment(:,:)
      integer(ik) :: ip, pf, np

      np = site%patch%n
      allocate(recruitment(cfg%pft%n, max(np, 1_ik)))
      recruitment = 0.0_wp
      associate (pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, pft%n
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%recruit_dens(pf)
            end do
         end do
      end associate
   end subroutine recruitment_density

end module meds_phenomenological_vital_rates
