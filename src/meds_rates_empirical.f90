!==========================================================================================!
! meds_rates_empirical -- the TEST vital-rate evaluator (lives OUTSIDE the demography engine).!
!                                                                                          !
! Given a site, produce the three demographic rate arrays that `update_demography` consumes:  !
! per-cohort growth, per-cohort total mortality, and per-(PFT,patch) recruitment. These are   !
! simple empirical relationships purely so the engine can be exercised end to end; a future   !
! mechanistic module would expose the same kind of routine and be swapped in by the master    !
! stepper. Because rates are now DATA, this routine READS the site (it is a support module,   !
! like setup/diagnostics) and computes the within-patch competition index itself.            !
!                                                                                          !
!   Growth (size + competition):  g = g0 * max(1 - dbh/dbh_crit, 0) * exp(-gcomp*comp)        !
!       comp = overtopping basal area [m2 m-2], from the height-sorted cohorts of the patch.  !
!   Mortality (TOTAL [1/yr]) = mort3 (background)                                             !
!                           + mort1*exp(-mort2*max(g,0)) (growth-dependent)                   !
!                           + treefall_rate (if dbh >= treefall_dbh)                          !
!                           + frost_mort * ramp((plant_min_temp - T)/5K) in [0,1]             !
!   Recruitment: constant per-PFT monthly density, gated by include_pft and a min temp.       !
!==========================================================================================!
module meds_rates_empirical
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : safe_exp
   use meds_pft_params,       only : pft_table_t
   use meds_config,           only : meds_config_t
   use meds_demography_types, only : site
   implicit none
   private

   public :: empirical_vital_rates

   real(wp), parameter :: frost_band = 5.0_wp     !< [K] width of the frost mortality ramp
   real(wp), parameter :: cm2_to_m2  = 1.0e-4_wp  !< basal_area [cm2/plant] * nplant -> m2/m2

contains

   subroutine empirical_vital_rates(asite, cfg, growth, mortality, recruitment)
      type(site),            intent(in)  :: asite
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr] per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]  per cohort, total
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)
      integer(ik) :: n, np, i, ip, pf, i0, i1
      real(wp)    :: cum, comp, dbh, g, chill

      n  = asite%coh%n
      np = asite%pat%n
      allocate(growth(n), mortality(n))
      allocate(recruitment(cfg%pft%n, max(np, 1_ik)))

      !----- Per-cohort growth & mortality (competition computed per height-sorted patch). -!
      associate (c => asite%coh, p => asite%pat, t => cfg%pft)
         do ip = 1_ik, np
            i0  = p%cohort_offset(ip)
            i1  = i0 + p%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping BA accumulator
            do i = i0, i1
               comp    = cum                               ! BA strictly above this cohort
               cum     = cum + c%nplant(i) * c%basal_area(i) * cm2_to_m2
               pf      = c%pft(i)
               dbh     = c%dbh(i)
               g       = t%g0(pf) * max(1.0_wp - dbh / t%dbh_crit(pf), 0.0_wp)               &
                       * safe_exp(-t%gcomp(pf) * comp)
               growth(i) = g
               chill     = (t%plant_min_temp(pf) - p%avg_daily_temp(ip)) / frost_band
               mortality(i) = t%mort3(pf)                                                    &
                            + t%mort1(pf) * safe_exp(-t%mort2(pf) * max(g, 0.0_wp))          &
                            + merge(t%treefall_rate(pf), 0.0_wp, dbh >= t%treefall_dbh(pf))  &
                            + t%frost_mort(pf) * min(max(chill, 0.0_wp), 1.0_wp)
            end do
         end do

         !----- Per-(PFT, patch) recruitment density. -----------------------------------!
         recruitment = 0.0_wp
         do ip = 1_ik, np
            do pf = 1_ik, t%n
               if (t%include_pft(pf) == 1_ik .and.                                          &
                   p%min_month_temp(ip) >= t%recruit_min_temp(pf)) then
                  recruitment(pf, ip) = t%recruit_dens(pf)
               end if
            end do
         end do
      end associate
   end subroutine empirical_vital_rates

end module meds_rates_empirical
