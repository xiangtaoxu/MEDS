!==========================================================================================!
! meds_mortality -- per-individual MORTALITY-hazard kernel of the phenomenological vital rates.!
!                                                                                          !
! Pure, stateless functions of (growth, traits). The hazard follows Camac et al. (2018, PNAS): !
! an additive baseline + a low-growth term decaying with the cohort's effective growth rate     !
! (its tracked running-mean growth, or its instantaneous growth before the average is seeded).  !
!==========================================================================================!
module meds_mortality
   use meds_kinds, only : wp
   implicit none
   private

   public :: effective_growth, mortality_hazard

   !----- A freshly created cohort carries GROWTH_AVG_UNSET (< 0) until its first growth step   !
   !       seeds the running mean; until then mortality uses its instantaneous growth. --------!
   real(wp), parameter :: growth_avg_unset = 0.0_wp   ! threshold: growth_avg < 0 means "unset"

contains

   !----- Effective growth for the mortality response: the tracked running-mean growth, or the  !
   !       instantaneous growth for a not-yet-seeded cohort (growth_avg < 0 sentinel). ---------!
   elemental pure function effective_growth(growth_avg, growth_instantaneous) result(g_eff)
      real(wp), intent(in) :: growth_avg, growth_instantaneous
      real(wp)             :: g_eff
      if (growth_avg < growth_avg_unset) then
         g_eff = growth_instantaneous
      else
         g_eff = growth_avg
      end if
   end function effective_growth

   !----- Camac et al. (2018) additive mortality hazard [1/yr]. --------------------------------!
   elemental pure function mortality_hazard(g_eff, mort_gamma, mort_alpha, mort_beta) result(mort)
      real(wp), intent(in) :: g_eff, mort_gamma, mort_alpha, mort_beta
      real(wp)             :: mort
      mort = mort_gamma + mort_alpha * exp(-mort_beta * g_eff)
   end function mortality_hazard

end module meds_mortality
