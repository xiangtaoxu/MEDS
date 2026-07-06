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
   use meds_constants, only : tiny_num, cp_air, cp_vap, cp_liq, cp_ice, latent_heat_fusion,    &
                              t_3ple, tsupercool_liq, tsupercool_vap
   implicit none
   private

   public :: sat_vapor_pressure, sat_specific_humidity, d_sat_vapor_pressure_dt
   public :: uext_to_temp, temp_to_uext
   public :: enthalpy_vapor, internal_energy_liquid, cp_moist, air_density
   public :: cas_enthalpy_of_temp, cas_temp_of_enthalpy

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

   !----- Clausius-Clapeyron slope d(e_sat)/dT [Pa/K] (derivative of the Bolton form). -------!
   elemental function d_sat_vapor_pressure_dt(t_k) result(desat)
      real(wp), intent(in) :: t_k
      real(wp)             :: desat, tc, esat
      tc    = t_k - 273.15_wp
      esat  = sat_vapor_pressure(t_k)
      desat = esat * 17.67_wp * 243.5_wp / (tc + 243.5_wp) ** 2
   end function d_sat_vapor_pressure_dt

   !----- Phase-change INVERTER: (internal energy, water mass, dry heat capacity) ->          !
   !      (temperature, liquid fraction). Consistent per-unit-VOLUME (soil: J/m3, kg/m3,       !
   !      J/m3/K) OR per-unit-AREA (leaf/wood: J/m2, kg/m2, J/m2/K). ED2 uextcm2tl/uint2tl.     !
   !      Continuous at u_freeze/u_melt (temp = t_3ple at both); dry (wmass=0) stores fall to    !
   !      temp = uext/dry_hcap via the ice/liquid branches (the plateau is empty).               !
   !---------------------------------------------------------------------------------------!
   elemental subroutine uext_to_temp(uext, wmass, dry_hcap, temp, fliq)
      real(wp), intent(in)  :: uext, wmass, dry_hcap
      real(wp), intent(out) :: temp, fliq
      real(wp) :: u_freeze, u_melt
      u_freeze = (dry_hcap + wmass * cp_ice) * t_3ple
      u_melt   = u_freeze + wmass * latent_heat_fusion
      if (uext <= u_freeze) then                            ! all ice
         fliq = 0.0_wp ; temp = uext / (dry_hcap + wmass * cp_ice)
      else if (uext >= u_melt) then                         ! all liquid
         fliq = 1.0_wp ; temp = (uext + wmass * cp_liq * tsupercool_liq) / (dry_hcap + wmass * cp_liq)
      else                                                  ! mixed-phase plateau (wmass > 0 here)
         temp = t_3ple ; fliq = (uext - u_freeze) / (wmass * latent_heat_fusion)
      end if
   end subroutine uext_to_temp

   !----- Forward map: (temperature, liquid fraction) -> internal energy. -------------------!
   elemental function temp_to_uext(dry_hcap, wmass, temp, fliq) result(uext)
      real(wp), intent(in) :: dry_hcap, wmass, temp, fliq
      real(wp)             :: uext
      uext = dry_hcap * temp + wmass * (fliq * cp_liq * (temp - tsupercool_liq)                &
                                        + (1.0_wp - fliq) * cp_ice * temp)
   end function temp_to_uext

   !----- Specific enthalpy of water vapour [J/kg] (thermal + phase baseline; any vapour       !
   !      flux automatically transports its latent heat, design 4b). ------------------------!
   elemental function enthalpy_vapor(t_k) result(h)
      real(wp), intent(in) :: t_k
      real(wp)             :: h
      h = cp_vap * (t_k - tsupercool_vap)
   end function enthalpy_vapor

   !----- Specific internal energy of liquid water [J/kg] (advected soil/xylem water). -------!
   elemental function internal_energy_liquid(t_k) result(u)
      real(wp), intent(in) :: t_k
      real(wp)             :: u
      u = cp_liq * (t_k - tsupercool_liq)
   end function internal_energy_liquid

   !----- Moist-air specific heat [J/kg/K]. -------------------------------------------------!
   elemental function cp_moist(shv) result(cp)
      real(wp), intent(in) :: shv
      real(wp)             :: cp
      cp = (1.0_wp - shv) * cp_air + shv * cp_vap
   end function cp_moist

   !----- Moist-air density [kg/m3] via the virtual temperature. ----------------------------!
   elemental function air_density(t_k, p_pa, shv) result(rho)
      real(wp), intent(in) :: t_k, p_pa, shv
      real(wp), parameter  :: r_dry = 287.04_wp             ! [J/kg/K] dry-air gas constant
      real(wp)             :: rho
      rho = p_pa / (r_dry * t_k * (1.0_wp + 0.608_wp * shv))
   end function air_density

   !----- Canopy-air specific enthalpy [J/kg] from temperature + specific humidity (4b). -----!
   elemental function cas_enthalpy_of_temp(t_k, shv) result(enth)
      real(wp), intent(in) :: t_k, shv
      real(wp)             :: enth
      enth = (1.0_wp - shv) * cp_air * t_k + shv * cp_vap * (t_k - tsupercool_vap)
   end function cas_enthalpy_of_temp

   !----- Canopy-air temperature [K] from specific enthalpy + specific humidity (inverse). ---!
   elemental function cas_temp_of_enthalpy(enth, shv) result(t_k)
      real(wp), intent(in) :: enth, shv
      real(wp)             :: t_k
      t_k = (enth + shv * cp_vap * tsupercool_vap) / ((1.0_wp - shv) * cp_air + shv * cp_vap)
   end function cas_temp_of_enthalpy

end module meds_thermo
