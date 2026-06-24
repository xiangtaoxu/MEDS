!==========================================================================================!
! meds_demography_interface -- the ONLY seam between the demographic engine and the demographic   !
! rates (growth, mortality, recruitment).                                                  !
!                                                                                          !
! The engine never computes a rate; it asks a `rate_provider_t` for one. Today the only    !
! concrete provider is the empirical test provider (meds_rates_empirical); tomorrow a       !
! mechanistic provider extends the SAME abstract type and drops in with no engine change.   !
!                                                                                          !
! Design choices for performance + portability:                                            !
!   * The deferred bindings take SCALAR plain-data arguments and are `pure`. They are       !
!     invoked from a HOST rate-fill loop that writes per-cohort rate arrays; the heavy      !
!     per-cohort APPLICATION kernels then run as `do concurrent` on those plain arrays and   !
!     are offloadable to GPU. Polymorphism therefore never reaches the device.             !
!   * Mortality returns a CAUSE VECTOR so the empirical form is unit-testable and the       !
!     engine can record cause-resolved diagnostics.                                        !
!   * Growth depends on a competition index `comp` (overtopping basal area, m2/m2) supplied  !
!     by the engine, realising the growth-size-competition relationship.                   !
!                                                                                          !
! Units: growth [cm/yr]; mortality [1/yr]; recruitment [plant/m2/month].                    !
!==========================================================================================!
module meds_demography_interface
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: rate_provider_t
   integer(ik), parameter, public :: n_mort_cause = 4_ik   !< age, growth, treefall, frost

   type, abstract :: rate_provider_t
   contains
      procedure(growth_if),    deferred :: growth
      procedure(mortality_if), deferred :: mortality
      procedure(recruit_if),   deferred :: recruit
   end type rate_provider_t

   abstract interface
      !------------------------------------------------------------------------------------!
      ! Per-capita DBH growth rate of ONE cohort [cm/yr].                                  !
      !   comp -- overtopping basal area [m2 m-2], the competition signal.                 !
      !------------------------------------------------------------------------------------!
      pure function growth_if(self, pft, dbh, hite, nplant, comp) result(ddbh)
         import :: rate_provider_t, wp, ik
         class(rate_provider_t), intent(in) :: self
         integer(ik),            intent(in) :: pft
         real(wp),               intent(in) :: dbh, hite, nplant, comp
         real(wp)                           :: ddbh
      end function growth_if

      !------------------------------------------------------------------------------------!
      ! Mortality CAUSE VECTOR of ONE cohort [1/yr]: (age, growth, treefall, frost).       !
      !   growth_pyr -- realised annual DBH increment [cm/yr] (mortality-growth coupling). !
      !------------------------------------------------------------------------------------!
      pure function mortality_if(self, pft, dbh, hite, nplant, growth_pyr,                 &
                                 avg_daily_temp, patch_age) result(mrate)
         import :: rate_provider_t, wp, ik, n_mort_cause
         class(rate_provider_t), intent(in) :: self
         integer(ik),            intent(in) :: pft
         real(wp),               intent(in) :: dbh, hite, nplant, growth_pyr,              &
                                               avg_daily_temp, patch_age
         real(wp)                           :: mrate(n_mort_cause)
      end function mortality_if

      !------------------------------------------------------------------------------------!
      ! Monthly recruit density for ONE (patch, PFT) [plant/m2/month].                     !
      !------------------------------------------------------------------------------------!
      pure function recruit_if(self, pft, dist_type, min_month_temp) result(rdens)
         import :: rate_provider_t, wp, ik
         class(rate_provider_t), intent(in) :: self
         integer(ik),            intent(in) :: pft, dist_type
         real(wp),               intent(in) :: min_month_temp
         real(wp)                           :: rdens
      end function recruit_if
   end interface

end module meds_demography_interface
