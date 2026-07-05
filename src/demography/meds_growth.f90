!==========================================================================================!
! meds_growth -- per-individual GROWTH-rate kernels of the phenomenological vital rates.      !
!                                                                                          !
! Pure, stateless functions of (structure, environment, traits) with NO cohort/patch state.   !
! The intrinsic diameter-growth rate is a capped log-linear function of dbh, multiplicatively   !
! suppressed by neighbourhood competition (overtopping LAI) and by reproductive allocation.     !
! The vital-rates orchestrator combines them per cohort in the height-sorted per-patch sweep.   !
!==========================================================================================!
module meds_growth
   use meds_kinds, only : wp
   implicit none
   private

   public :: growth_intrinsic, competition_suppression, reproduction_suppression

contains

   !----- Intrinsic diameter growth [cm/yr]: capped log-linear in dbh (= growth_dbh_max for  !
   !       dbh >= growth_dbh_cap; larger for smaller stems).                                   !
   elemental pure function growth_intrinsic(dbh, growth_dbh_max, growth_dbh_slope, growth_dbh_cap) result(gi)
      real(wp), intent(in) :: dbh, growth_dbh_max, growth_dbh_slope, growth_dbh_cap
      real(wp)             :: gi
      gi = growth_dbh_max * exp(growth_dbh_slope * max(0.0_wp, log(growth_dbh_cap) - log(dbh)))
   end function growth_intrinsic

   !----- Competition suppression in (0, 1]: Beer-like attenuation by the overtopping LAI.   !
   elemental pure function competition_suppression(over_lai, growth_lai_slope) result(supp)
      real(wp), intent(in) :: over_lai, growth_lai_slope
      real(wp)             :: supp
      supp = exp(growth_lai_slope * over_lai)
   end function competition_suppression

   !----- Reproductive-allocation suppression of growth: 1 below the maturity height, else    !
   !       1 - reproduction_investment_fraction (the growth fraction diverted to reproduction).!
   elemental pure function reproduction_suppression(height, min_reproduction_height,          &
                                                    reproduction_investment_fraction) result(supp)
      real(wp), intent(in) :: height, min_reproduction_height, reproduction_investment_fraction
      real(wp)             :: supp
      if (height >= min_reproduction_height) then
         supp = 1.0_wp - reproduction_investment_fraction
      else
         supp = 1.0_wp
      end if
   end function reproduction_suppression

end module meds_growth
