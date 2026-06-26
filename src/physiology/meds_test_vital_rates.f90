!==========================================================================================!
! meds_test_vital_rates -- the TEST vital-rate evaluator (lives OUTSIDE the demography engine).!
!                                                                                          !
! Given a site_t, produce the per-cohort growth and mortality rate arrays that                 !
! `update_demography` consumes. The rates are driven by LIGHT competition: within each          !
! (height-sorted) patch a top-down sweep accumulates the OVERTOPPING LAI above each cohort, from !
! which a Beer-Lambert light fraction follows. Cohorts at EQUAL height are co-dominant -- they   !
! share the overtopping LAI and do not shade one another -- so the result is independent of      !
! their arbitrary within-layer order (no height/PFT bias). Because rates are DATA, this READS the !
! site_t (it is a support module, like setup/diagnostics) and computes the competition index.    !
! (Recruitment is a separate rate, supplied by meds_recruitment.)                              !
!                                                                                          !
!   light        = exp(-light_ext * overtopping_LAI)                  in (0, 1]                !
!   growth [cm/yr]   = gr_max  * light * dbh        (realized relative growth x diameter)        !
!   mortality [1/yr] = mort_base + mort_shade * (1 - light)   (baseline + shade-driven)          !
!                                                                                          !
! Low-wood-density PFTs carry a high gr_max AND a high mort_shade, so they out-grow others in   !
! the light but die fast once overtopped -- the classic ED growth-mortality trade-off.         !
!==========================================================================================!
module meds_test_vital_rates
   use meds_kinds,            only : wp, ik
   use meds_allometry,        only : light_ext
   use meds_config,           only : meds_config_t
   use meds_demography_types, only : site_t
   implicit none
   private

   public :: test_vital_rates

contains

   subroutine test_vital_rates(site, cfg, growth, mortality)
      type(site_t),            intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr] per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]  per cohort, total
      integer(ik) :: n, np, i, j, k, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, light, layer_lai

      n  = site%cohort%n
      np = site%patch%n
      allocate(growth(n), mortality(n))

      !----- Per-cohort growth & mortality (overtopping LAI computed per height-sorted patch). !
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
   end subroutine test_vital_rates

end module meds_test_vital_rates
