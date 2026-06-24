!==========================================================================================!
! meds_rates_empirical -- the TEST rate provider.                                          !
!                                                                                          !
! A concrete `rate_provider_t` implementing simple, well-behaved empirical relationships    !
! purely so the demographic engine can be exercised end-to-end. These are NOT mechanistic:  !
! a future meds_rates_mechanistic module extends the same abstract type and drops in with    !
! no change to the engine.                                                                  !
!                                                                                          !
!   Growth (size + competition):  ddbh = g0 * max(1 - dbh/dbh_crit, 0) * exp(-gcomp*comp)   !
!       -> growth slows as a stem approaches its asymptotic size AND as overtopping basal    !
!          area (comp) shades it. Pioneers are fast and shade-intolerant; climax slow and    !
!          tolerant.                                                                        !
!                                                                                          !
!   Mortality (growth-dependent, ED2 "scheme 2" analog), cause vector [1/yr]:               !
!       (1) background  = mort3                                                             !
!       (2) growth      = mort1 * exp(-mort2 * max(growth,0))   (slow growers die)          !
!       (3) treefall    = treefall_rate if dbh >= treefall_dbh                              !
!       (4) frost       = frost_mort * ramp((plant_min_temp - T)/5K) in [0,1]               !
!                                                                                          !
!   Recruitment (schedule): constant per-PFT monthly density, gated by include_pft and a    !
!       minimum monthly temperature.                                                        !
!==========================================================================================!
module meds_rates_empirical
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : safe_exp
   use meds_pft_params,     only : pft_table_t
   use meds_demography_interface, only : rate_provider_t, n_mort_cause
   implicit none
   private

   public :: empirical_provider_t, make_empirical_provider

   real(wp), parameter :: frost_band = 5.0_wp     !< [K] width of the frost mortality ramp

   type, extends(rate_provider_t) :: empirical_provider_t
      !----- Small parameter copies (self-contained; future mechanistic provider differs). !
      real(wp), allocatable :: g0(:), gcomp(:), dbh_crit(:)
      real(wp), allocatable :: mort1(:), mort2(:), mort3(:)
      real(wp), allocatable :: frost_mort(:), plant_min_temp(:)
      real(wp), allocatable :: treefall_dbh(:), treefall_rate(:)
      real(wp), allocatable :: recruit_dens(:), recruit_min_temp(:)
      integer(ik), allocatable :: include_pft(:)
   contains
      procedure :: growth    => emp_growth
      procedure :: mortality => emp_mortality
      procedure :: recruit   => emp_recruit
   end type empirical_provider_t

contains

   !---------------------------------------------------------------------------------------!
   ! Build a provider from a PFT trait table (copies only the traits it needs).            !
   !---------------------------------------------------------------------------------------!
   function make_empirical_provider(pft) result(prov)
      type(pft_table_t), intent(in) :: pft
      type(empirical_provider_t)    :: prov
      prov%g0             = pft%g0
      prov%gcomp          = pft%gcomp
      prov%dbh_crit       = pft%dbh_crit
      prov%mort1          = pft%mort1
      prov%mort2          = pft%mort2
      prov%mort3          = pft%mort3
      prov%frost_mort     = pft%frost_mort
      prov%plant_min_temp = pft%plant_min_temp
      prov%treefall_dbh   = pft%treefall_dbh
      prov%treefall_rate  = pft%treefall_rate
      prov%recruit_dens   = pft%recruit_dens
      prov%recruit_min_temp = pft%recruit_min_temp
      prov%include_pft    = pft%include_pft
   end function make_empirical_provider

   !---------------------------------------------------------------------------------------!
   pure function emp_growth(self, pft, dbh, hite, nplant, comp) result(ddbh)
      class(empirical_provider_t), intent(in) :: self
      integer(ik),                 intent(in) :: pft
      real(wp),                    intent(in) :: dbh, hite, nplant, comp
      real(wp)                                :: ddbh
      ddbh = self%g0(pft) * max(1.0_wp - dbh / self%dbh_crit(pft), 0.0_wp)                 &
           * safe_exp(-self%gcomp(pft) * comp)
   end function emp_growth

   !---------------------------------------------------------------------------------------!
   pure function emp_mortality(self, pft, dbh, hite, nplant, growth_pyr,                   &
                               avg_daily_temp, patch_age) result(mrate)
      class(empirical_provider_t), intent(in) :: self
      integer(ik),                 intent(in) :: pft
      real(wp),                    intent(in) :: dbh, hite, nplant, growth_pyr,            &
                                                 avg_daily_temp, patch_age
      real(wp)                                :: mrate(n_mort_cause)
      real(wp)                                :: chill

      mrate(1) = self%mort3(pft)
      mrate(2) = self%mort1(pft) * safe_exp(-self%mort2(pft) * max(growth_pyr, 0.0_wp))
      mrate(3) = merge(self%treefall_rate(pft), 0.0_wp, dbh >= self%treefall_dbh(pft))
      chill    = (self%plant_min_temp(pft) - avg_daily_temp) / frost_band
      mrate(4) = self%frost_mort(pft) * min(max(chill, 0.0_wp), 1.0_wp)
   end function emp_mortality

   !---------------------------------------------------------------------------------------!
   pure function emp_recruit(self, pft, dist_type, min_month_temp) result(rdens)
      class(empirical_provider_t), intent(in) :: self
      integer(ik),                 intent(in) :: pft, dist_type
      real(wp),                    intent(in) :: min_month_temp
      real(wp)                                :: rdens
      if (self%include_pft(pft) == 1_ik .and. min_month_temp >= self%recruit_min_temp(pft)) then
         rdens = self%recruit_dens(pft)
      else
         rdens = 0.0_wp
      end if
   end function emp_recruit

end module meds_rates_empirical
