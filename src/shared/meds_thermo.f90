!==========================================================================================!
! meds_thermo -- shared moist-air thermodynamic / psychrometric helpers (saturation vapour   !
! pressure and specific humidity). Placed in the shared foundation so every fast-loop          !
! biophysics process that needs humidity -- ground evaporation (meds_column_hydrology), the      !
! future canopy-air-space water balance, and the leaf/soil energy balance -- shares ONE          !
! saturation formula instead of re-deriving it. `pure`/`elemental`; the Bolton (1980)             !
! coefficients are universal empirical constants (like meds_temp_response's Arrhenius             !
! constants), not tunable model parameters, so they live in the formula.                          !
!==========================================================================================!
module meds_thermo
   use meds_kinds,     only : wp
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: sat_vapor_pressure, sat_specific_humidity

contains

   !----- Saturation vapour pressure [Pa] over liquid water (Bolton 1980). ----------------!
   elemental function sat_vapor_pressure(t_k) result(esat)
      real(wp), intent(in) :: t_k
      real(wp)             :: esat, tc
      tc   = t_k - 273.15_wp
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
   end function sat_vapor_pressure

   !----- Saturation specific humidity [kg/kg] at temperature t_k [K] and pressure p_pa [Pa]. !
   elemental function sat_specific_humidity(t_k, p_pa) result(qs)
      real(wp), intent(in) :: t_k, p_pa
      real(wp)             :: qs, esat
      esat = sat_vapor_pressure(t_k)
      qs   = 0.622_wp * esat / max(p_pa - 0.378_wp * esat, tiny_num)
   end function sat_specific_humidity

end module meds_thermo
