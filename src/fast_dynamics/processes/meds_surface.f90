! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_surface -- the ground's surface: the bare-ground skin and the one-layer snow        !
! store, blended by snowfac, each giving the soil-top boundary and the canopy-air fluxes   !
! (docs/science/snow_biophysics.md):                                                       !
!   * ground_surface_fluxes -- the bare ground's sensible and latent fluxes;               !
!   * snow_cover_fraction / snow_accumulate / snow_drain_meltwater -- the snow's mass;     !
!   * snow_surface_fluxes / snow_base_conductance / snow_surface_balance /                !
!     snow_take_exchange / snow_energy_step -- its energy.                                 !
!==========================================================================================!
module meds_surface
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan_boltzmann, &
                                       latent_heat_fusion, tsupercool_liq
   use meds_therm_lib,          only : internal_energy_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_liquid, precip_enthalpy
   use meds_patch_state_vars,   only : snow_state_t
   use meds_biophysics_params,  only : snow_params_t
   implicit none
   private

   public :: ground_surface_fluxes
   public :: snow_cover_fraction, snow_accumulate, snow_drain_meltwater
   public :: snow_surface_balance, snow_take_exchange
   public :: snow_energy_step      ! used only by the tests: a balance and its commit in one step
   public :: snow_surface_fluxes   ! used only by the tests
   public :: snow_env_t, snow_flux_t, snow_melt_t

   !----- Surface boundary conditions for the snow-surface energy balance (read-only, forced). --!
   type :: snow_env_t
      real(wp) :: abs_sw   = 0.0_wp, abs_lw = 0.0_wp    !< [W/m2] the snow-blended ground's absorbed SW and net LW
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp !< [K],[kg/kg] CAS state (forced sibling)
      real(wp) :: ggnet    = 0.0_wp                     !< [m/s] ground<->CAS conductance (heat = vapour)
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp !< [kg/m3],[Pa] CAS air
      real(wp) :: t_soil_top = 273.16_wp                !< [K] top soil-node temperature (snow-base conduction target)
      real(wp) :: k_soil_top = 0.0_wp                   !< [W/m/K] top soil layer's thermal conductivity (the caller's)
      real(wp) :: dz_soil_top = 0.05_wp                 !< [m] top soil node depth |z_node(1)| (interface spacing)
   end type snow_env_t

   !----- Snow-surface energy-balance outputs (all positive upward; per unit ground area). ------!
   type :: snow_flux_t
      real(wp) :: t_surf  = 273.16_wp, fliq = 0.0_wp    !< [K],[-] diagnosed snow-surface temperature + liquid fraction
      real(wp) :: h_snow  = 0.0_wp, le_snow = 0.0_wp    !< [W/m2] sensible + latent (sublimation/evap) to the CAS
      real(wp) :: w_flux  = 0.0_wp                      !< [kg/m2/s] vapour mass to CAS (>0 sublimation, <0 deposition)
      real(wp) :: g_base  = 0.0_wp                      !< [W/m2] conduction snow-base -> soil top (the new soil top BC)
      real(wp) :: rnet    = 0.0_wp                      !< [W/m2] net all-wave radiation absorbed by the snow surface
      real(wp) :: energy_in = 0.0_wp                    !< [W/m2] what the pack gains: rnet - h_snow - le_snow - g_base
      real(wp) :: snowfac = 0.0_wp                      !< [-] Niu-Yang07 snow-cover / burial fraction
      real(wp) :: energy_resid = 0.0_wp                 !< [J/m2] closed-budget residual (~0 by construction)
   end type snow_flux_t

   !----- Snow-mass outputs from accumulation + melt drainage (a store->store handoff to the pond). -!
   type :: snow_melt_t
      real(wp) :: melt_mass = 0.0_wp                    !< [kg/m2] free liquid drained this step -> the pond
      real(wp) :: melt_enth = 0.0_wp                    !< [J/m2]  matching enthalpy (internal_energy_liquid(t_3ple))
      real(wp) :: dump_mass = 0.0_wp                    !< [kg/m2] full-melt-out residual water dumped to the pond
      real(wp) :: dump_enth = 0.0_wp                    !< [J/m2]  full-melt-out residual energy dumped to the pond
      logical  :: melted_out = .false.                  !< .true. when the layer vanished (nlayer -> 0)
      real(wp) :: mass_resid = 0.0_wp                   !< [kg/m2] closed-budget residual (~0 by construction)
   end type snow_melt_t

contains

   !=======================================================================================!
   !  Ground / surface skin                                                                 !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Bare-ground fluxes to the canopy air, instantaneous: sensible through ggnet, latent   !
   ! from the water network's soil evaporation (already area-weighted) at the vapour       !
   ! enthalpy the canopy air receives. The caller assembles G_top and the snow blend.      !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine ground_surface_fluxes(t_ground, t_cas, ggnet, rho, soil_evap, h_bare, le_soil)
      real(wp), intent(in)  :: t_ground, t_cas, ggnet, rho, soil_evap
      real(wp), intent(out) :: h_bare, le_soil
      h_bare  = ggnet * rho * cp_air * (t_ground - t_cas)
      le_soil = soil_evap * enthalpy_vapor(t_ground)
   end subroutine ground_surface_fluxes

   !=======================================================================================!
   !  Snow -- mass side: cover fraction, accumulation, meltwater drainage.                  !
   !=======================================================================================!

   !----- Niu-Yang (2007) snow-cover fraction tanh(depth/scale): it ramps the ground optics and the  !
   !      soil-top blend and area-weights the pack's exchange (docs/science/snow_biophysics.md). ----!
   pure function snow_cover_fraction(swe, snow_depth, params, snow_rough) result(snowfac)
      real(wp),            intent(in) :: swe, snow_depth
      type(snow_params_t), intent(in) :: params
      real(wp),            intent(in) :: snow_rough   !< [m] the snow surface's roughness ([aerodynamics])
      real(wp) :: snowfac, scale
      if (swe <= params%tiny_snow_mass) then
         snowfac = 0.0_wp
      else
         scale   = max(params%ny07_a * snow_rough * (params%rho_snow / 100.0_wp) ** params%ny07_m, tiny_num)
         snowfac = tanh(snow_depth / scale)
      end if
   end function snow_cover_fraction

   !----- Accumulation of snowfall (as ice) and rain-on-snow (as liquid) with their precip_enthalpy.  !
   !      A layer is created only from min_new_snow_mass of new snow; with no layer and less, the      !
   !      caller routes the precipitation to the ground with the throughfall. -----------------------!
   pure subroutine snow_accumulate(snow, snowfall, rainf, air_temp, dt, params)
      type(snow_state_t), intent(inout) :: snow
      real(wp),            intent(in)    :: snowfall, rainf, air_temp, dt
      type(snow_params_t), intent(in)    :: params
      real(wp) :: dm_snow, dm_rain
      logical  :: has_layer
      dm_snow = max(0.0_wp, snowfall) * dt
      dm_rain = max(0.0_wp, rainf) * dt
      has_layer = (snow%nlayer >= 1_ik) .or. (dm_snow >= params%min_new_snow_mass)
      if (.not. has_layer) return                          ! bare ground + sub-threshold snow: caller handles it
      snow%swe(1)         = snow%swe(1) + dm_snow + dm_rain
      snow%snow_energy(1) = snow%snow_energy(1) + precip_enthalpy(dm_rain, dm_snow, air_temp)
      snow%snow_depth(1)  = snow%snow_depth(1) + dm_snow / params%rho_snow   ! rain fills pores / refreezes -> no bulk depth
      snow%nlayer         = 1_ik
   end subroutine snow_accumulate

   !----- Percolation: drain the free liquid above the holding capacity as a paired (mass, enthalpy) !
   !      hand-off to the pond (melt_*), and dump the residual at melt-out. ------------------------!
   pure subroutine snow_drain_meltwater(snow, params, melt)
      type(snow_state_t), intent(inout) :: snow
      type(snow_params_t), intent(in)    :: params
      type(snow_melt_t),   intent(out)   :: melt
      real(wp) :: wmass_free, hold
      melt = snow_melt_t()
      if (snow%nlayer == 0_ik) return

      !----- Vanished / trace layer: dump the whole residual to the soil, revert to bare ground. -----!
      if (snow%swe(1) <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%swe(1))
         melt%dump_enth  = snow%snow_energy(1)
         melt%melted_out = .true.
         call clear_layer(snow)
         return
      end if

      call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))

      !----- Drain the free liquid above the holding capacity (LEAF-3 1:9). ---------------------------!
      hold       = max(1.0_wp - params%liquid_holding_frac, tiny_num)
      wmass_free = max(0.0_wp, snow%swe(1) * (snow%snow_fliq(1) - params%liquid_holding_frac) / hold)
      wmass_free = min(wmass_free, snow%swe(1))
      if (wmass_free > 0.0_wp) then
         melt%melt_mass      = wmass_free
         melt%melt_enth      = wmass_free * internal_energy_liquid(t_3ple)   ! drains at the plateau
         snow%swe(1)         = snow%swe(1) - wmass_free
         snow%snow_energy(1) = snow%snow_energy(1) - melt%melt_enth
         snow%snow_depth(1)  = snow%swe(1) / params%rho_snow
      end if

      !----- If draining took the layer below the melt-out threshold, dump the remainder to soil. -----!
      if (snow%swe(1) <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%swe(1))
         melt%dump_enth  = snow%snow_energy(1)
         melt%melted_out = .true.
         call clear_layer(snow)
      else
         call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      end if
   end subroutine snow_drain_meltwater

   !----- Revert to the empty "no snow" state (nlayer = 0). ----------------------------------------!
   pure subroutine clear_layer(snow)
      type(snow_state_t), intent(inout) :: snow
      snow%swe(1)         = 0.0_wp
      snow%snow_energy(1) = 0.0_wp
      snow%snow_depth(1)  = 0.0_wp
      snow%snow_fliq(1)   = 0.0_wp
      snow%snow_temp(1)   = t_3ple
      snow%nlayer         = 0_ik
   end subroutine clear_layer

   !=======================================================================================!
   !  Snow -- energy side: surface fluxes, base conduction, backward-Euler energy advance.  !
   !=======================================================================================!

   !----- Sensible and latent fluxes at a snow-surface temperature (positive up), the latent at the  !
   !      vapour enthalpy the canopy air receives, and saturated over ice by the pack's own liquid  !
   !      fraction (docs/science/snow_biophysics.md eq. (1)). ---------------------------------------!
   pure subroutine snow_surface_fluxes(t_surf, fliq, env, h_flux, w_flux, le_flux)
      real(wp),         intent(in)  :: t_surf
      real(wp),         intent(in)  :: fliq      !< [-] liquid fraction of the snow surface
      type(snow_env_t), intent(in)  :: env
      real(wp),         intent(out) :: h_flux, w_flux, le_flux
      h_flux = env%ggnet * env%rho_air * cp_air * (t_surf - env%can_temp)
      w_flux = env%ggnet * env%rho_air * (sat_specific_humidity(t_surf, env%press, fliq) - env%can_shv)
      le_flux = w_flux * enthalpy_vapor(t_surf)
   end subroutine snow_surface_fluxes

   !----- Snow-base to soil-top conductance [W/m2/K]: the half snow layer and the top soil node in    !
   !      series, so a deepening pack decouples the soil (docs/science/snow_biophysics.md). ----------!
   pure function snow_base_conductance(snow_depth, env, params) result(gcond)
      real(wp),            intent(in) :: snow_depth
      type(snow_env_t),    intent(in) :: env
      type(snow_params_t), intent(in) :: params
      real(wp) :: gcond, r_series
      r_series = 0.5_wp * max(snow_depth, tiny_num) / max(params%k_snow, tiny_num)                  &
               + max(env%dz_soil_top, tiny_num) / max(env%k_soil_top, tiny_num)
      gcond    = 1.0_wp / max(r_series, tiny_num)
   end function snow_base_conductance

   !----- The pack's exchange with its surroundings over dt, from a linearized backward-Euler step on  !
   !      its enthalpy that it does not commit (docs/science/snow_biophysics.md, "The snow store"):     !
   !      with every flux linearized about T^n, swe e(T) = E^n + dt af R(T) has one root, found in      !
   !      closed form on the plateau or the ice or liquid line, and the fluxes are taken at that T.     !
   !      Sublimating mass leaves with the pack's specific energy, so the rest keeps T. The march      !
   !      calls this at every evaluation, against the canopy air and soil top it holds then; the pack   !
   !      takes the exchange's time integral afterwards (snow_take_exchange). ------------------------!
   pure subroutine snow_surface_balance(snow, env, params, dt, area_frac, flux)
      type(snow_state_t),  intent(in)  :: snow
      type(snow_env_t),    intent(in)  :: env
      type(snow_params_t), intent(in)  :: params
      real(wp),            intent(in)  :: dt
      real(wp),            intent(in)  :: area_frac    !< snow-cover fraction: scales all boundary exchange (sub-column)
      type(snow_flux_t),   intent(out) :: flux
      real(wp) :: t_n, fliq_n, e_old, swe0, gcond, af, dtaf
      real(wp) :: h_n, w_n, le_n, g_n, dh_dt, dw_dt, demit_dt, l_vap, r0, drdt
      real(wp) :: u_freeze, u_melt, e_bal, t_new, dtemp, h_s, w_s, le_s, g_s, emit

      flux = snow_flux_t()
      af = max(0.0_wp, min(1.0_wp, area_frac))
      if (snow%nlayer == 0_ik .or. snow%swe(1) <= params%tiny_snow_mass .or. af <= 0.0_wp) then
         flux%t_surf = env%t_soil_top                              ! no snow surface (or zero cover) to balance
         return
      end if

      swe0  = snow%swe(1)
      e_old = snow%snow_energy(1)
      dtaf  = dt * af
      flux%snowfac = af
      call internal_energy_to_temp(e_old, swe0, 0.0_wp, t_n, fliq_n)
      gcond = snow_base_conductance(snow%snow_depth(1), env, params)

      !----- The fluxes at T^n and their slopes. A kilogram of vapour carries enthalpy_vapor; less the    !
      !      pack's specific energy, that is the latent heat it takes from the pack (sublimation on ice,  !
      !      evaporation on water). ------------------------------------------------------------------!
      call snow_surface_fluxes(t_n, fliq_n, env, h_n, w_n, le_n)
      g_n      = gcond * (t_n - env%t_soil_top)
      dh_dt    = env%ggnet * env%rho_air * cp_air
      dw_dt    = env%ggnet * env%rho_air * sat_specific_humidity_temp_deriv(t_n, env%press, fliq_n)
      demit_dt = 4.0_wp * params%snow_emiss * stefan_boltzmann * t_n ** 3
      l_vap    = enthalpy_vapor(t_n) - e_old / swe0
      r0       = env%abs_sw + env%abs_lw - h_n - l_vap * w_n - g_n
      drdt     = -(demit_dt + dh_dt + l_vap * dw_dt + gcond)                ! < 0

      !----- swe e(T) = E^n + dt af (r0 + drdt (T - T^n)): the plateau first, else the ice or liquid line. --!
      u_freeze = swe0 * cp_ice * t_3ple
      u_melt   = u_freeze + swe0 * latent_heat_fusion
      t_new    = t_3ple
      e_bal    = e_old + dtaf * (r0 + drdt * (t_3ple - t_n))
      if (e_bal < u_freeze) then                                             ! frozen: e = cp_ice T
         t_new = (e_old + dtaf * (r0 - drdt * t_n)) / (swe0 * cp_ice - dtaf * drdt)
         e_bal = swe0 * cp_ice * t_new
      else if (e_bal > u_melt) then                                          ! liquid: e = cp_liq (T - T_sc)
         t_new = (e_old + swe0 * cp_liq * tsupercool_liq + dtaf * (r0 - drdt * t_n))                   &
                 / (swe0 * cp_liq - dtaf * drdt)
         e_bal = swe0 * internal_energy_liquid(t_new)
      end if

      !----- The fluxes at T, on the linearization. ------------------------------------------------!
      dtemp = t_new - t_n
      h_s   = h_n + dh_dt * dtemp
      w_s   = w_n + dw_dt * dtemp
      g_s   = g_n + gcond * dtemp
      emit  = demit_dt * dtemp                     ! the emission added as the surface moves from T^n to T
      if (w_s * dtaf < swe0) then
         le_s = w_s * (l_vap + e_bal / swe0)       ! the latent heat plus the leaving mass's own energy
      else                                         ! the whole pack sublimates, with all the energy left
         w_s  = swe0 / dtaf
         le_s = e_old / dtaf + env%abs_sw + env%abs_lw - emit - h_s - g_s
      end if

      !----- The area-weighted fluxes, positive away from the pack: h_snow + le_snow -> the canopy air's !
      !      enthalpy, w_flux -> its vapour, g_base -> the soil top; energy_in is what the pack gains. ---!
      flux%t_surf    = t_new
      flux%h_snow    = af * h_s
      flux%le_snow   = af * le_s
      flux%w_flux    = af * w_s
      flux%g_base    = af * g_s
      flux%rnet      = af * (env%abs_sw + env%abs_lw - emit)   ! snow's area-weighted net radiation into the store
      flux%energy_in = flux%rnet - flux%h_snow - flux%le_snow - flux%g_base
   end subroutine snow_surface_balance

   !----- The pack takes an exchange's amounts: energy_amount [J/m2] in, vapour_amount [kg/m2] out,  !
   !      keeping e(T) per kilogram, so its temperature reads off the new energy. A pack whose mass    !
   !      is gone is left at t_3ple for the drain kernel to empty. --------------------------------!
   pure subroutine snow_take_exchange(snow, energy_amount, vapour_amount, params)
      type(snow_state_t),  intent(inout) :: snow
      real(wp),            intent(in)    :: energy_amount, vapour_amount
      type(snow_params_t), intent(in)    :: params
      if (snow%nlayer == 0_ik) return
      snow%snow_energy(1) = snow%snow_energy(1) + energy_amount
      snow%swe(1)         = max(0.0_wp, snow%swe(1) - vapour_amount)
      snow%snow_depth(1)  = snow%swe(1) / params%rho_snow
      if (snow%swe(1) > params%tiny_snow_mass) then
         call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      else
         snow%snow_temp(1) = t_3ple                               ! vanished; drain kernel dumps residual
         snow%snow_fliq(1) = 0.0_wp
      end if
   end subroutine snow_take_exchange

   !----- One whole step at a fixed environment: the balance, then the pack takes it. The unit the    !
   !      snow tests check; energy_resid is zero by construction. --------------------------------!
   pure subroutine snow_energy_step(snow, env, params, dt, area_frac, flux)
      type(snow_state_t),  intent(inout) :: snow
      type(snow_env_t),    intent(in)    :: env
      type(snow_params_t), intent(in)    :: params
      real(wp),            intent(in)    :: dt
      real(wp),            intent(in)    :: area_frac
      type(snow_flux_t),   intent(out)   :: flux
      real(wp) :: e_old
      e_old = snow%snow_energy(1)
      call snow_surface_balance(snow, env, params, dt, area_frac, flux)
      if (flux%snowfac <= 0.0_wp) return
      call snow_take_exchange(snow, dt * flux%energy_in, dt * flux%w_flux, params)
      flux%t_surf = merge(snow%snow_temp(1), flux%t_surf, snow%swe(1) > params%tiny_snow_mass)
      flux%fliq   = snow%snow_fliq(1)
      flux%energy_resid = (snow%snow_energy(1) - e_old) - dt * flux%energy_in
   end subroutine snow_energy_step

end module meds_surface
