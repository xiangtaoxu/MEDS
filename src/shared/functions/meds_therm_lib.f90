! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_therm_lib -- the thermodynamics the fast-loop stores share, so each uses one        !
! formula: tissue heat capacity and energy (docs/science/vegetation_energy_dynamics.md     !
! eqs. (8)-(10)); saturation vapour pressure and humidity (docs/science/snow_biophysics.md !
! eq. (1), docs/science/forcing.md eq. (11)); the energy -> temperature inverter with      !
! phase change, and the specific energies of liquid, ice, vapour and canopy air in one     !
! datum (MEDS_ENERGY_BALANCE_DESIGN.md §3.1); and the soil's thermal conductivity and      !
! heat capacity (docs/science/soil_biophysics.md eq. (11)). The Bolton (1980) and          !
! Johansen coefficients are fixed empirical constants, not model parameters.               !
!==========================================================================================!
module meds_therm_lib
   use meds_kinds,     only : wp
   use meds_constants, only : mmdry, tiny_num, cp_air, cp_vap, cp_liq, cp_ice, latent_heat_fusion,    &
                              t_3ple, tsupercool_liq, tsupercool_vap, r_dry,                    &
                              rho_h2o, k_water, k_ice, c2b_wood
   implicit none
   private

   public :: tissue_heat_capacity, tissue_dry_heat_capacity, tissue_energy, tissue_temp
   public :: sat_vapor_pressure, sat_specific_humidity, sat_vapor_pressure_temp_deriv
   public :: air_vpd, specific_humidity_to_vpd

   public :: sat_specific_humidity_temp_deriv
   public :: internal_energy_to_temp, temp_to_internal_energy
   public :: enthalpy_vapor, internal_energy_liquid, internal_energy_ice, air_density, cas_molar_density
   public :: temp_of_liquid_enthalpy, precip_enthalpy
   public :: cas_enthalpy_of_temp, cas_temp_of_enthalpy
   public :: soil_thermal_cond, soil_heat_cap_vol

   real(wp), parameter :: SR_FLOOR = 0.05_wp     !< [-] floor on relative saturation, where Johansen's log10 form holds

contains

   !----- A tissue's heat capacity [J/m2/K]: its dry matter plus the liquid water it holds !
   !      (docs/science/vegetation_energy_dynamics.md eq. (8)). ---------------------------!
   elemental pure real(wp) function tissue_heat_capacity(nplant, carbon, c_dry, water) result(cap)
      real(wp), intent(in) :: nplant   !< [plant/m2]
      real(wp), intent(in) :: carbon   !< [kgC/plant] the tissue's carbon
      real(wp), intent(in) :: c_dry    !< [J/kg/K]    specific heat of its dry matter
      real(wp), intent(in) :: water    !< [kg/plant]  the water it holds
      cap = tissue_dry_heat_capacity(nplant, carbon, c_dry) + water * nplant * cp_liq
   end function tissue_heat_capacity

   !----- The dry part of a tissue's heat capacity [J/m2/K]: its carbon as dry mass, at    !
   !      its dry specific heat (tissue_c_leaf or tissue_c_wood_dry). ---------------------!
   elemental pure real(wp) function tissue_dry_heat_capacity(nplant, carbon, c_dry) result(cap)
      real(wp), intent(in) :: nplant   !< [plant/m2]
      real(wp), intent(in) :: carbon   !< [kgC/plant] the tissue's carbon
      real(wp), intent(in) :: c_dry    !< [J/kg/K]    specific heat of its dry matter
      cap = carbon * nplant * c2b_wood * c_dry
   end function tissue_dry_heat_capacity

   !----- A tissue's internal energy [J per plant, or per m2 with per-m2 inputs], with its !
   !      water liquid at any temperature, since a live tissue is taken not to freeze      !
   !      (docs/science/vegetation_energy_dynamics.md eq. (9)). ---------------------------!
   elemental pure real(wp) function tissue_energy(dry_hcap, wmass, temp) result(e)
      real(wp), intent(in) :: dry_hcap   !< [J/K]  dry heat capacity
      real(wp), intent(in) :: wmass      !< [kg]   the water it holds
      real(wp), intent(in) :: temp       !< [K]
      e = dry_hcap * temp + wmass * cp_liq * (temp - tsupercool_liq)
   end function tissue_energy

   !----- The inverse (eq. (10)): the temperature of a tissue holding energy e. A tissue   !
   !      with neither dry matter nor water has no temperature of its own, and the caller  !
   !      gives it its parent's. ----------------------------------------------------------!
   elemental pure real(wp) function tissue_temp(e, wmass, dry_hcap) result(temp)
      real(wp), intent(in) :: e          !< [J]
      real(wp), intent(in) :: wmass      !< [kg]
      real(wp), intent(in) :: dry_hcap   !< [J/K]
      temp = (e + wmass * cp_liq * tsupercool_liq) / (dry_hcap + wmass * cp_liq)
   end function tissue_temp

   !----- Saturation vapour pressure [Pa] (docs/science/snow_biophysics.md eq. (1)):       !
   !      Bolton (1980) over liquid, blended with the Magnus form over ice by the liquid   !
   !      fraction fliq. Without fliq it is over liquid, as dewpoint and VPD need          !
   !      (docs/science/forcing.md §7). ---------------------------------------------------!
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

   !----- Saturation specific humidity [kg/kg]: q = 0.622 e / (p - 0.378 e), the           !
   !      conversion of docs/science/forcing.md eq. (11). ---------------------------------!
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

   !----- The slope d(e_sat)/dT [Pa/K], at fixed fliq: on the melt plateau, where fliq     !
   !      changes, the temperature is pinned at the triple point. -------------------------!
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

   !----- The slope d(q_sat)/dT [1/K] for the latent-flux linearizations                   !
   !      (docs/science/vegetation_energy_dynamics.md eq. (2)). ---------------------------!
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

   !----- Phase-change inverter: internal energy, water and dry heat capacity ->           !
   !      temperature and liquid fraction, per m3 of soil or per m2 of surface alike       !
   !      (MEDS_ENERGY_BALANCE_DESIGN.md §3.1). With no water the plateau is empty. -------!
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

   !----- The inverse: temperature and liquid fraction -> internal energy. ----------------!
   elemental function temp_to_internal_energy(dry_hcap, wmass, temp, fliq) result(uext)
      real(wp), intent(in) :: dry_hcap, wmass, temp, fliq
      real(wp)             :: uext
      uext = dry_hcap * temp + wmass * (fliq * cp_liq * (temp - tsupercool_liq)                &
                                        + (1.0_wp - fliq) * cp_ice * temp)
   end function temp_to_internal_energy

   !----- Specific enthalpy of water vapour [J/kg] in the common datum, so a vapour flux   !
   !      carries its latent heat (MEDS_ENERGY_BALANCE_DESIGN.md §4b). --------------------!
   elemental function enthalpy_vapor(t_k) result(h)
      real(wp), intent(in) :: t_k
      real(wp)             :: h
      h = cp_vap * (t_k - tsupercool_vap)
   end function enthalpy_vapor

   !----- Specific internal energy of liquid water [J/kg] (advected soil/xylem water). ----!
   elemental function internal_energy_liquid(t_k) result(u)
      real(wp), intent(in) :: t_k
      real(wp)             :: u
      u = cp_liq * (t_k - tsupercool_liq)
   end function internal_energy_liquid

   !----- The inverse of internal_energy_liquid, which turns a (mass, energy) transfer     !
   !      into the temperature a consumer takes (docs/science/soil_plant_water.md §5). ----!
   elemental function temp_of_liquid_enthalpy(u) result(t_k)
      real(wp), intent(in) :: u
      real(wp)             :: t_k
      t_k = tsupercool_liq + u / cp_liq
   end function temp_of_liquid_enthalpy

   !----- Specific internal energy of ice [J/kg], zero at 0 K: the datum of the inverter's !
   !      all-ice branch, so ice melting at t_3ple costs exactly latent_heat_fusion. ------!
   elemental function internal_energy_ice(t_k) result(u)
      real(wp), intent(in) :: t_k
      real(wp)             :: u
      u = cp_ice * t_k
   end function internal_energy_ice

   !----- The enthalpy precipitation brings, wherever it lands: rain as liquid at the air  !
   !      temperature above the canopy, snow as ice no warmer than freezing                !
   !      (docs/science/forcing.md §8). Rates give [W/m2], amounts [J/m2]. ----------------!
   elemental function precip_enthalpy(rain, snow, t_air) result(h)
      real(wp), intent(in) :: rain, snow   !< [kg/m2/s] or [kg/m2]
      real(wp), intent(in) :: t_air        !< [K] air temperature above the canopy
      real(wp)             :: h
      h = snow * internal_energy_ice(min(t_3ple, t_air)) + rain * internal_energy_liquid(t_air)
   end function precip_enthalpy

   !----- Moist-air density [kg/m3] via the virtual temperature. --------------------------!
   elemental function air_density(t_k, p_pa, shv) result(rho)
      real(wp), intent(in) :: t_k, p_pa, shv
      real(wp)             :: rho
      rho = p_pa / (r_dry * t_k * (1.0_wp + 0.608_wp * shv))
   end function air_density

   !----- Canopy-air specific enthalpy [J/kg] from temperature and specific humidity: dry  !
   !      air from 0 K, vapour in the common datum                                         !
   !      (docs/science/canopy_air_space_biophysics.md §The three CAS twins). -------------!
   elemental function cas_enthalpy_of_temp(t_k, shv) result(enth)
      real(wp), intent(in) :: t_k, shv
      real(wp)             :: enth
      enth = (1.0_wp - shv) * cp_air * t_k + shv * cp_vap * (t_k - tsupercool_vap)
   end function cas_enthalpy_of_temp

   !----- The inverse: canopy-air temperature [K] from specific enthalpy and humidity. ----!
   elemental function cas_temp_of_enthalpy(enth, shv) result(t_k)
      real(wp), intent(in) :: enth, shv
      real(wp)             :: t_k
      t_k = (enth + shv * cp_vap * tsupercool_vap) / ((1.0_wp - shv) * cp_air + shv * cp_vap)
   end function cas_temp_of_enthalpy

   !=======================================================================================!
   !  Soil thermal properties (docs/science/soil_biophysics.md §The soil thermal column).  !
   !=======================================================================================!

   !----- Soil thermal conductivity [W/m/K], Johansen's made ice-aware                     !
   !      (docs/science/soil_biophysics.md eq. (11)). -------------------------------------!
   elemental function soil_thermal_cond(theta, fliq, theta_sat, k_solid, k_dry) result(kappa)
      real(wp), intent(in) :: theta, fliq, theta_sat, k_solid, k_dry
      real(wp)             :: kappa, s_r, k_e, kappa_sat
      s_r = min(max(theta / max(theta_sat, tiny_num), SR_FLOOR), 1.0_wp)      ! relative saturation
      !----- Kersten number: Johansen's log form unfrozen, Farouki's linear frozen. -------!
      k_e = fliq * (log10(s_r) + 1.0_wp) + (1.0_wp - fliq) * s_r
      k_e = min(max(k_e, 0.0_wp), 1.0_wp)
      !----- Saturated conductivity: the geometric mean of solids, liquid and ice. --------!
      kappa_sat = k_solid ** (1.0_wp - theta_sat)                                            &
                  * k_water ** (theta_sat * fliq) * k_ice ** (theta_sat * (1.0_wp - fliq))
      kappa = k_e * kappa_sat + (1.0_wp - k_e) * k_dry
   end function soil_thermal_cond

   !----- Soil volumetric heat capacity [J/m3/K]: the solids plus the water and ice. ------!
   elemental function soil_heat_cap_vol(theta, fliq, dry_cvol) result(c_eff)
      real(wp), intent(in) :: theta, fliq, dry_cvol
      real(wp)             :: c_eff
      c_eff = dry_cvol + theta * rho_h2o * (fliq * cp_liq + (1.0_wp - fliq) * cp_ice)
   end function soil_heat_cap_vol

   !----- Molar density of dry air in moist air [mol/m3], rho (1 - q)/mmdry: the canopy    !
   !      air's CO2 capacity per metre of depth and the scale of its exchange with the     !
   !      atmosphere (docs/science/canopy_air_space_biophysics.md, the CO2 twin). ---------!
   elemental function cas_molar_density(rho, shv) result(dmol)
      real(wp), intent(in) :: rho, shv
      real(wp)             :: dmol
      dmol = rho * (1.0_wp - shv) / mmdry
   end function cas_molar_density


   !=======================================================================================!
   !  Vapour-pressure deficit. The fast loop and the output both use it, and a kernel may  !
   !  not depend on the output layer, so it lives here.                                    !
   !=======================================================================================!

   !----- VPD [Pa] from temperature and the actual vapour pressure. -----------------------!
   pure elemental real(wp) function air_vpd(temp, e_vap) result(vpd)
      real(wp), intent(in) :: temp    !< [K]
      real(wp), intent(in) :: e_vap   !< [Pa] actual vapour pressure
      vpd = max(0.0_wp, sat_vapor_pressure(temp) - e_vap)
   end function air_vpd

   !----- VPD [Pa] from temperature, specific humidity and pressure. The vapour pressure   !
   !      e = q p / (0.622 + 0.378 q) is the exact inverse of docs/science/forcing.md      !
   !      eq. (11), so a saturated forcing record reads back as VPD = 0. ------------------!
   pure elemental real(wp) function specific_humidity_to_vpd(temp, shv, pressure) result(vpd)
      real(wp), intent(in) :: temp      !< [K]
      real(wp), intent(in) :: shv       !< [kg/kg] specific humidity
      real(wp), intent(in) :: pressure  !< [Pa]
      real(wp) :: e_vap
      e_vap = shv * pressure / max(0.622_wp + 0.378_wp * shv, tiny_num)
      vpd   = air_vpd(temp, e_vap)
   end function specific_humidity_to_vpd

end module meds_therm_lib
