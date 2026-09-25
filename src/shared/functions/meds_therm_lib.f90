! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_therm_lib -- shared thermodynamic / thermal-property helpers for the fast loop: moist-air  !
! psychrometrics (saturation vapour pressure + specific humidity and their temperature slopes),  !
! the phase-change energy inverter, moist-air enthalpy/density, AND the SOIL thermal properties   !
! (Johansen conductivity + volumetric heat capacity -- the thermal twin of the soil retention      !
! curves in meds_hydr_lib). Placed in the shared foundation so every fast-loop biophysics       !
! process shares ONE formula instead of re-deriving it. `pure`/`elemental`; the Bolton (1980) and   !
! Johansen coefficients are universal empirical constants, not tunable model parameters.             !
!==========================================================================================!
module meds_therm_lib
   use meds_kinds,     only : wp
   use meds_constants, only : mmdry, tiny_num, cp_air, cp_vap, cp_liq, cp_ice, latent_heat_fusion,    &
                              t_3ple, tsupercool_liq, tsupercool_vap, r_dry, r_gas, r_wv,       &
                              rho_h2o, k_water, k_ice
   implicit none
   private

   public :: sat_vapor_pressure, sat_specific_humidity, sat_vapor_pressure_temp_deriv
   public :: air_vpd, specific_humidity_to_vpd

   !----- mmh2o/mmdry, for the specific-humidity -> vapour-pressure conversion. -----------!
   real(wp), parameter :: EPS_MOL = (r_gas / mmdry) / r_wv
   public :: sat_specific_humidity_temp_deriv
   public :: internal_energy_to_temp, temp_to_internal_energy
   public :: enthalpy_vapor, internal_energy_liquid, internal_energy_ice, cp_moist, air_density, cas_molar_density
   public :: temp_of_liquid_enthalpy
   public :: cas_enthalpy_of_temp, cas_temp_of_enthalpy
   !----- SOIL thermal properties (conductivity + volumetric heat capacity): thermal-property   !
   !      constitutive kernels, the thermal twin of the soil retention curves. -----------------!
   public :: soil_thermal_cond, soil_heat_cap_vol

   real(wp), parameter :: SR_FLOOR = 0.05_wp     !< Kersten log10 floor (Johansen validity limit)

contains

   !---------------------------------------------------------------------------------------!
   ! SATURATION VAPOUR PRESSURE, over liquid water OR over ice (#89).                          !
   !                                                                                          !
   ! Magnus/Bolton pair, sharing the triple-point constant so the two branches meet EXACTLY:    !
   !                                                                                          !
   !     e_liq(Tc) = 611.2 exp( 17.67 Tc / (Tc + 243.5) )      (Bolton 1980)                    !
   !     e_ice(Tc) = 611.2 exp( 21.87 Tc / (Tc + 265.5) )      (Magnus form over ice)           !
   !                                                                                          !
   ! WHY A BLEND AND NOT A BRANCH ON TEMPERATURE. `fliq` is the liquid fraction the store        !
   ! already carries prognostically (internal_energy_to_temp returns it), so the weighting is    !
   ! a state the model tracks rather than a threshold invented here. It also keeps e_sat         !
   ! CONTINUOUS through the melt plateau: a `T < 0` branch would put a step of up to 10 % of      !
   ! e_sat into the right-hand side at the triple point, and this model's adaptive controller     !
   ! has already been broken once by exactly that kind of jump (see veg_energy_balance's P6       !
   ! note). Because both forms give 611.2 Pa at Tc = 0, the blend is continuous there for any     !
   ! fliq, so a store that freezes or melts moves smoothly between the curves.                   !
   !                                                                                          !
   ! `fliq` ABSENT means pure liquid, which is what every caller did before #89 -- so omitting    !
   ! it is bit-identical. Supply it only where the surface's phase is genuinely known (snow,      !
   ! frozen soil). Dewpoint conversions and diagnostic VPD must NOT supply it: dewpoint is        !
   ! DEFINED over liquid, so an ice branch there would mis-convert the forcing.                   !
   !                                                                                          !
   ! Accuracy of the ice form against Murphy & Koop (2005): 0.1 % at -10 C, 0.4 % at -20 C,       !
   ! 0.9 % at -30 C. The error it removes is an order larger -- using the liquid curve over ice    !
   ! overstates e_sat by 10 % at -10 C, 22 % at -20 C and 34 % at -30 C.                           !
   !---------------------------------------------------------------------------------------!
   elemental function sat_vapor_pressure(t_k, fliq) result(esat)
      real(wp), intent(in)           :: t_k
      real(wp), intent(in), optional :: fliq    !< liquid fraction of the surface (absent = 1, all liquid)
      real(wp)                       :: esat, tc, fl
      tc   = t_k - 273.15_wp
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
      if (present(fliq)) then
         fl = min(1.0_wp, max(0.0_wp, fliq))
         if (fl < 1.0_wp) esat = fl * esat                                                       &
                               + (1.0_wp - fl) * 611.2_wp * exp(21.8745584_wp * tc / (tc + 265.5_wp))
      end if
   end function sat_vapor_pressure

   !----- Saturation specific humidity [kg/kg] at temperature t_k [K] and pressure p_pa [Pa]. !
   elemental function sat_specific_humidity(t_k, p_pa, fliq) result(qs)
      real(wp), intent(in)           :: t_k, p_pa
      real(wp), intent(in), optional :: fliq
      real(wp)                       :: qs, esat
      if (present(fliq)) then
         esat = sat_vapor_pressure(t_k, fliq)
      else
         esat = sat_vapor_pressure(t_k)
      end if
      qs   = 0.622_wp * esat / max(p_pa - 0.378_wp * esat, tiny_num)
   end function sat_specific_humidity

   !----- Clausius-Clapeyron slope d(e_sat)/dT [Pa/K]. `fliq` is held FIXED in the derivative:  !
   !      it is a frozen coefficient of the linearization, and on the melt plateau temperature   !
   !      is pinned at the triple point anyway, so its T-dependence has nothing to act on. ------!
   elemental function sat_vapor_pressure_temp_deriv(t_k, fliq) result(desat)
      real(wp), intent(in)           :: t_k
      real(wp), intent(in), optional :: fliq
      real(wp)                       :: desat, tc, fl, e_liq, e_ice
      tc    = t_k - 273.15_wp
      e_liq = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
      desat = e_liq * 17.67_wp * 243.5_wp / (tc + 243.5_wp) ** 2
      if (present(fliq)) then
         fl = min(1.0_wp, max(0.0_wp, fliq))
         if (fl < 1.0_wp) then
            e_ice = 611.2_wp * exp(21.8745584_wp * tc / (tc + 265.5_wp))
            desat = fl * desat                                                                   &
                  + (1.0_wp - fl) * e_ice * 21.8745584_wp * 265.5_wp / (tc + 265.5_wp) ** 2
         end if
      end if
   end function sat_vapor_pressure_temp_deriv

   !----- d(sat_specific_humidity)/dT [1/K] at temperature t_k [K], pressure p_pa [Pa] -- the  !
   !      Clausius-Clapeyron slope of qsat, folding the (p - 0.378*esat) denominator. Shared by  !
   !      every implicit latent-flux linearization (leaf/wood, snow, CAS, ground energy step).   !
   elemental function sat_specific_humidity_temp_deriv(t_k, p_pa, fliq) result(dqsdt)
      real(wp), intent(in)           :: t_k, p_pa
      real(wp), intent(in), optional :: fliq
      real(wp)                       :: dqsdt, esat, desat
      if (present(fliq)) then
         esat  = sat_vapor_pressure(t_k, fliq)
         desat = sat_vapor_pressure_temp_deriv(t_k, fliq)
      else
         esat  = sat_vapor_pressure(t_k)
         desat = sat_vapor_pressure_temp_deriv(t_k)
      end if
      dqsdt = 0.622_wp * p_pa / max((p_pa - 0.378_wp * esat) ** 2, tiny_num) * desat
   end function sat_specific_humidity_temp_deriv

   !----- Phase-change INVERTER: (internal energy, water mass, dry heat capacity) ->          !
   !      (temperature, liquid fraction). Consistent per-unit-VOLUME (soil: J/m3, kg/m3,       !
   !      J/m3/K) OR per-unit-AREA (leaf/wood: J/m2, kg/m2, J/m2/K). ED2 uextcm2tl/uint2tl.     !
   !      Continuous at u_freeze/u_melt (temp = t_3ple at both); dry (wmass=0) stores fall to    !
   !      temp = uext/dry_hcap via the ice/liquid branches (the plateau is empty).               !
   !---------------------------------------------------------------------------------------!
   elemental subroutine internal_energy_to_temp(uext, wmass, dry_hcap, temp, fliq)
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
   end subroutine internal_energy_to_temp

   !----- Forward map: (temperature, liquid fraction) -> internal energy. -------------------!
   elemental function temp_to_internal_energy(dry_hcap, wmass, temp, fliq) result(uext)
      real(wp), intent(in) :: dry_hcap, wmass, temp, fliq
      real(wp)             :: uext
      uext = dry_hcap * temp + wmass * (fliq * cp_liq * (temp - tsupercool_liq)                &
                                        + (1.0_wp - fliq) * cp_ice * temp)
   end function temp_to_internal_energy

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

   !----- EXACT inverse of internal_energy_liquid: the temperature a liquid store of specific        !
   !      internal energy u is at. Used to turn a (mass, enthalpy) pair back into the temperature    !
   !      that values it, so a paired transfer can be handed to a consumer that wants a temperature  !
   !      rather than an energy (e.g. the meltwater the snow pack sends to the ponding store). -----!
   elemental function temp_of_liquid_enthalpy(u) result(t_k)
      real(wp), intent(in) :: u
      real(wp)             :: t_k
      t_k = tsupercool_liq + u / cp_liq
   end function temp_of_liquid_enthalpy

   !----- Specific internal energy of ICE [J/kg] (frozen store: snow/frost). Shares the 0-K ice   !
   !      datum of internal_energy_to_temp's all-ice branch (u = wmass*cp_ice*T at dry_hcap=0), so a snow     !
   !      layer's energy seeded with temp_to_internal_energy(0,swe,T,0) inverts back to T exactly, and melt   !
   !      (ice at t_3ple -> liquid at t_3ple) costs latent_heat_fusion via internal_energy_liquid's !
   !      tsupercool_liq offset -- ONE datum across ice/liquid/vapour (snow<->soil<->CAS closure).  !
   elemental function internal_energy_ice(t_k) result(u)
      real(wp), intent(in) :: t_k
      real(wp)             :: u
      u = cp_ice * t_k
   end function internal_energy_ice

   !----- Moist-air specific heat [J/kg/K]. -------------------------------------------------!
   elemental function cp_moist(shv) result(cp)
      real(wp), intent(in) :: shv
      real(wp)             :: cp
      cp = (1.0_wp - shv) * cp_air + shv * cp_vap
   end function cp_moist

   !----- Moist-air density [kg/m3] via the virtual temperature. ----------------------------!
   elemental function air_density(t_k, p_pa, shv) result(rho)
      real(wp), intent(in) :: t_k, p_pa, shv
      real(wp)             :: rho
      rho = p_pa / (r_dry * t_k * (1.0_wp + 0.608_wp * shv))   ! r_dry from meds_constants (single authority)
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

   !=======================================================================================!
   !  SOIL thermal properties (design MEDS_ENERGY_BALANCE_DESIGN.md, 4c/8). Both `elemental`; !
   !  fliq = liquid fraction of the soil water, so freeze/thaw enters ice-aware.              !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Soil thermal conductivity [W/m/K] (Johansen). Kersten-number blend of dry & saturated. !
   !---------------------------------------------------------------------------------------!
   elemental function soil_thermal_cond(theta, fliq, theta_sat, k_solid, k_dry) result(kappa)
      real(wp), intent(in) :: theta, fliq, theta_sat, k_solid, k_dry
      real(wp)             :: kappa, s_r, k_e, kappa_sat
      s_r = min(max(theta / max(theta_sat, tiny_num), SR_FLOOR), 1.0_wp)      ! relative saturation
      !----- Kersten number: blend Johansen log10 form (unfrozen) with Farouki linear form (frozen)  !
      !      by liquid fraction, so the interpolation weight is ice-aware like kappa_sat below. ------!
      k_e = fliq * (log10(s_r) + 1.0_wp) + (1.0_wp - fliq) * s_r
      k_e = min(max(k_e, 0.0_wp), 1.0_wp)                                     ! Kersten number, clamped
      !----- Saturated geometric mean, ice-aware (fixes ED2's liquid-only kappa_sat). -----!
      kappa_sat = k_solid ** (1.0_wp - theta_sat)                                            &
                  * k_water ** (theta_sat * fliq) * k_ice ** (theta_sat * (1.0_wp - fliq))
      kappa = k_e * kappa_sat + (1.0_wp - k_e) * k_dry
   end function soil_thermal_cond

   !---------------------------------------------------------------------------------------!
   ! Effective volumetric heat capacity [J/m3/K] = dry matrix + phase-appropriate water.    !
   !---------------------------------------------------------------------------------------!
   elemental function soil_heat_cap_vol(theta, fliq, dry_cvol) result(c_eff)
      real(wp), intent(in) :: theta, fliq, dry_cvol
      real(wp)             :: c_eff
      c_eff = dry_cvol + theta * rho_h2o * (fliq * cp_liq + (1.0_wp - fliq) * cp_ice)
   end function soil_heat_cap_vol

   !----- Molar density of DRY air in moist air [mol/m3] = rho*(1-q)/mmdry: the CAS CO2 box     !
   !      capacity per metre of depth and the CAS<->atmosphere CO2 conductance scale. ----------!
   elemental function cas_molar_density(rho, shv) result(dmol)
      real(wp), intent(in) :: rho, shv
      real(wp)             :: dmol
      dmol = rho * (1.0_wp - shv) / mmdry
   end function cas_molar_density


   !=======================================================================================!
   !  Vapour-pressure deficit. Thermodynamics, not a diagnostic: the fast loop accumulates    !
   !  a dt-weighted CAS VPD (#264) and the output layer reads it back, so both need the SAME  !
   !  formula and it cannot live in src/io -- a kernel may not depend on the output layer.    !
   !=======================================================================================!

   !----- VPD [Pa] from temperature and ACTUAL vapour pressure. ---------------------------!
   pure elemental real(wp) function air_vpd(temp, e_vap) result(vpd)
      real(wp), intent(in) :: temp    !< [K]
      real(wp), intent(in) :: e_vap   !< [Pa] actual vapour pressure
      vpd = max(0.0_wp, sat_vapor_pressure(temp) - e_vap)
   end function air_vpd

   !----- VPD [Pa] from the canopy-air-space prognostic twins (temperature + SPECIFIC       !
   !      humidity) at a given pressure: e = q*p / (eps + (1-eps)*q), eps = mmh2o/mmdry.    !
   pure elemental real(wp) function specific_humidity_to_vpd(temp, shv, pressure) result(vpd)
      real(wp), intent(in) :: temp      !< [K]
      real(wp), intent(in) :: shv       !< [kg/kg] specific humidity
      real(wp), intent(in) :: pressure  !< [Pa]
      real(wp) :: e_vap
      e_vap = shv * pressure / max(EPS_MOL + (1.0_wp - EPS_MOL) * shv, tiny_num)
      vpd   = air_vpd(temp, e_vap)
   end function specific_humidity_to_vpd

end module meds_therm_lib
