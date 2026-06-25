!==========================================================================================!
! meds_test_vital_rates -- the TEST vital-rate evaluator (lives OUTSIDE the demography engine).!
!                                                                                          !
! Given a site_t, produce the three demographic rate arrays that `update_demography` consumes:  !
! per-cohort growth, per-cohort total mortality, and per-(PFT,patch) recruitment. The rates    !
! are driven by LIGHT competition: within each (height-sorted) patch a top-down sweep          !
! accumulates the OVERTOPPING LAI above each cohort, from which a Beer-Lambert light fraction   !
! follows. Because rates are DATA, this routine READS the site_t (it is a support module, like    !
! setup/diagnostics) and computes the competition index itself.                               !
!                                                                                          !
!   light        = exp(-light_ext * overtopping_LAI)                  in (0, 1]                !
!   growth [cm/yr]   = gr_max  * light * dbh        (realized relative growth x diameter)        !
!   mortality [1/yr] = mort_base + mort_shade * (1 - light)   (baseline + shade-driven)          !
!                                                                                          !
! Low-wood-density PFTs carry a high gr_max AND a high mort_shade, so they out-grow others in   !
! the light but die fast once overtopped -- the classic ED growth-mortality trade-off.         !
! Recruitment: constant per-PFT monthly density, gated by include_pft and a min temperature.    !
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

   subroutine test_vital_rates(site, cfg, growth, mortality, recruitment)
      type(site_t),            intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr] per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]  per cohort, total
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)
      integer(ik) :: n, np, i, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, light

      n  = site%cohort%n
      np = site%patch%n
      allocate(growth(n), mortality(n))
      allocate(recruitment(cfg%pft%n, max(np, 1_ik)))

      !----- Per-cohort growth & mortality (overtopping LAI computed per height-sorted patch). !
      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, np
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping LAI accumulator
            do i = i0, i1
               over_lai  = cum                             ! LAI strictly above this cohort
               cum       = cum + cohort%nplant(i) * cohort%leaf_area(i)  ! [m2/m2] this cohort's contribution
               light     = exp(-light_ext * over_lai)
               pf        = cohort%pft(i)
               growth(i)    = pft%gr_max(pf) * light * cohort%dbh(i)
               mortality(i) = pft%mort_base(pf) + pft%mort_shade(pf) * (1.0_wp - light)
            end do
         end do

         !----- Per-(PFT, patch) recruitment density. -----------------------------------!
         recruitment = 0.0_wp
         do ip = 1_ik, np
            do pf = 1_ik, pft%n
               if (pft%include_pft(pf) == 1_ik .and.                                          &
                   patch%min_month_temp(ip) >= pft%recruit_min_temp(pf)) then
                  recruitment(pf, ip) = pft%recruit_dens(pf)
               end if
            end do
         end do
      end associate
   end subroutine test_vital_rates

end module meds_test_vital_rates
