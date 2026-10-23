! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_ground_biophysics -- the stateless land-surface interface kernels: the bare-ground  !
! skin and the one-layer snow / temporary-surface-water store, blended by snowfac, each    !
! giving the soil-top boundary and the canopy-air fluxes (docs/science/snow_biophysics.md):!
!   * ground_surface_fluxes -- the bare ground's sensible and latent fluxes;               !
!   * snow_cover_fraction / snow_accumulate / snow_drain_meltwater -- the snow's mass;     !
!   * snow_surface_fluxes / snow_base_conductance / snow_energy_step -- its energy.        !
!==========================================================================================!
module meds_ground_biophysics
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan_boltzmann, &
                                       latent_heat_fusion, tsupercool_liq
   use meds_therm_lib,          only : internal_energy_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_liquid, precip_enthalpy
   use meds_patch_state_vars,   only : snow_state_t
   use meds_soil_types, only : snow_env_t, snow_flux_t, snow_melt_t
   use meds_biophysics_params,  only : snow_params_t
   implicit none
   private

   public :: ground_surface_fluxes
   public :: snow_cover_fraction, snow_accumulate, snow_drain_meltwater
   public :: snow_energy_step
   public :: snow_surface_fluxes   ! used only by the tests

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

   !----- Percolation after the energy step: drain the free liquid above the holding capacity as a   !
   !      paired (mass, enthalpy) hand-off to the pond (melt_*), and dump the residual at melt-out. --!
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

   !----- Advance the snow store one dt by a linearized backward-Euler step on its enthalpy           !
   !      (docs/science/snow_biophysics.md, "The snow store"): with every flux linearized about T^n,   !
   !      swe e(T) = E^n + dt af R(T) has one root, found in closed form on the plateau or the ice or  !
   !      liquid line. Sublimating mass leaves with the pack's specific energy, so the rest keeps T. -!
   pure subroutine snow_energy_step(snow, env, params, dt, area_frac, flux)
      type(snow_state_t), intent(inout) :: snow
      type(snow_env_t),    intent(in)    :: env
      type(snow_params_t), intent(in)    :: params
      real(wp),            intent(in)    :: dt
      real(wp),            intent(in)    :: area_frac    !< snow-cover fraction: scales all boundary exchange (sub-column)
      type(snow_flux_t),   intent(out)   :: flux
      real(wp) :: t_n, fliq_n, e_old, swe0, swe1, gcond, af, dtaf
      real(wp) :: h_n, w_n, le_n, g_n, dh_dt, dw_dt, demit_dt, l_vap, r0, drdt
      real(wp) :: u_freeze, u_melt, e_bal, t_new, dtemp, h_s, w_s, le_s, g_s, emit, net_flux

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
         swe1 = swe0 - w_s * dtaf
      else                                         ! the whole pack sublimates, with all the energy left
         w_s  = swe0 / dtaf
         le_s = e_old / dtaf + env%abs_sw + env%abs_lw - emit - h_s - g_s
         swe1 = 0.0_wp
      end if

      !----- Energy and mass update: the pack keeps e(T) per kilogram, so it reads off at T. ---------!
      net_flux = af * (env%abs_sw + env%abs_lw - emit - h_s - le_s - g_s)   ! sub-column area weight
      snow%snow_energy(1) = e_old + dt * net_flux
      snow%swe(1)         = swe1
      snow%snow_depth(1)  = max(0.0_wp, snow%swe(1)) / params%rho_snow

      if (snow%swe(1) > params%tiny_snow_mass) then
         call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      else
         snow%snow_temp(1) = t_3ple                               ! vanished; drain kernel dumps residual
         snow%snow_fliq(1) = 0.0_wp
      end if

      !----- Report the turbulent + conduction fluxes for the CAS + soil coupling. h_snow + le_snow -> !
      !      CAS enthalpy; w_flux -> CAS vapour; g_base -> soil top BC. ------------------------------!
      flux%t_surf  = merge(snow%snow_temp(1), t_new, snow%swe(1) > params%tiny_snow_mass)
      flux%fliq    = snow%snow_fliq(1)
      flux%h_snow  = af * h_s
      flux%le_snow = af * le_s
      flux%w_flux  = af * w_s
      flux%g_base  = af * g_s
      flux%rnet    = af * (env%abs_sw + env%abs_lw - emit)   ! snow's area-weighted net radiation into the store
      flux%energy_resid = (snow%snow_energy(1) - e_old) - dt * net_flux
   end subroutine snow_energy_step

end module meds_ground_biophysics
