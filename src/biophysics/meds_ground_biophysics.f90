!==========================================================================================!
! meds_ground_biophysics -- the stateless LAND-SURFACE-INTERFACE kernels: the bare-ground skin  !
! energy balance and the snow / temporary-surface-water store. These are the two mutually-        !
! exclusive MODES of one interface (bare vs snow-covered, blended by snowfac), both producing the  !
! soil-top thermal BC and the surface->CAS fluxes.                                                  !
!                                                                                          !
! GROUND skin:                                                                                     !
!   * ground_surface_fluxes -- the BARE-ground sensible + latent fluxes to the CAS (soil_evap is     !
!     the frozen hydrology-authority mass flux); the caller assembles G_top and the snow blend.      !
! SNOW store (design MEDS_SNOW_DESIGN.md P0, single bulk layer; temp + liquid fraction are a         !
! read-off of the shared inverter uext_to_temp, so MELT/refreeze is the internal-energy plateau):    !
!   * snow_cover_fraction / snow_accumulate / snow_drain_meltwater  -- the MASS side.                 !
!   * snow_surface_fluxes / snow_base_conductance / snow_energy_step -- the ENERGY side.               !
!==========================================================================================!
module meds_ground_biophysics
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan,           &
                                       latent_heat_vap
   use meds_therm_lib,          only : uext_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_ice, internal_energy_liquid
   use meds_column_state_types, only : snow_column_t
   use meds_biophysics_types,   only : leaf_energy_env_t, snow_params_t, snow_env_t, snow_flux_t,  &
                                       snow_melt_t
   implicit none
   private

   public :: ground_surface_fluxes
   public :: snow_cover_fraction, snow_accumulate, snow_drain_meltwater
   public :: snow_surface_fluxes, snow_base_conductance, snow_energy_step

contains

   !=======================================================================================!
   !  GROUND / surface skin                                                                 !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Bare-ground surface fluxes to the CAS -- INSTANTANEOUS, no state advance. Sensible from   !
   ! the ground<->CAS conductance ggnet; latent from the FROZEN hydrology-authority mass flux   !
   ! soil_evap (already area-weighted by the caller via chydro_forcing_t%snow_free_frac), on the SAME !
   ! enthalpy twin the CAS receives. t_ground = soil_temp(1) is supplied by the soil kernel.      !
   ! The caller assembles G_top (net into the soil surface) and the snow-fraction blend, since    !
   ! those diverge between the bare and snow-covered paths.                                         !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine ground_surface_fluxes(t_ground, t_cas, ggnet, rho, soil_evap, h_bare, le_soil)
      real(wp), intent(in)  :: t_ground, t_cas, ggnet, rho, soil_evap
      real(wp), intent(out) :: h_bare, le_soil
      h_bare  = ggnet * rho * cp_air * (t_ground - t_cas)
      le_soil = soil_evap * enthalpy_vapor(t_ground)
   end subroutine ground_surface_fluxes

   !=======================================================================================!
   !  SNOW -- MASS side: cover fraction, accumulation, meltwater drainage.                  !
   !=======================================================================================!

   !----- Niu-Yang (2007) snow-cover / burial fraction from geometric depth. Monotone tanh in      !
   !      [0,1]; the (rho_snow/rho_fresh) normalization folds into the tuned scale. Used to RAMP     !
   !      the ground optics, the soil-top-BC blend, and the aerodynamic roughness (design §4f/§4g/§6). !
   pure function snow_cover_fraction(swe, snow_depth, params) result(snowfac)
      real(wp),            intent(in) :: swe, snow_depth
      type(snow_params_t), intent(in) :: params
      real(wp) :: snowfac, scale
      if (swe <= params%tiny_snow_mass) then
         snowfac = 0.0_wp
      else
         scale   = max(params%ny07_a * params%z0_snow * (params%rho_snow / 100.0_wp) ** params%ny07_m, tiny_num)
         snowfac = tanh(snow_depth / scale)
      end if
   end function snow_cover_fraction

   !----- Accumulation: frozen precip (snowf) + rain-on-snow (rainf) onto the single bulk layer.    !
   !      snowf lands as ICE at min(t_3ple, tair); rain-on-snow lands as LIQUID at tair and refreezes !
   !      automatically via the inverter downstream. A layer is CREATED only when the total new snow  !
   !      mass reaches min_new_snow_mass; below that (and on bare ground, nlayer=0) NOTHING is added   !
   !      here -- the caller folds sub-threshold snowfall into the soil-top store and routes rain to    !
   !      infiltration (design §4a). Mass + enthalpy are conserved on the shared datum.                 !
   pure subroutine snow_accumulate(snow, snowf, rainf, tair, dt, params)
      type(snow_column_t), intent(inout) :: snow
      real(wp),            intent(in)    :: snowf, rainf, tair, dt
      type(snow_params_t), intent(in)    :: params
      real(wp) :: dm_snow, dm_rain
      logical  :: has_layer
      dm_snow = max(0.0_wp, snowf) * dt
      dm_rain = max(0.0_wp, rainf) * dt
      has_layer = (snow%nlayer >= 1_ik) .or. (dm_snow >= params%min_new_snow_mass)
      if (.not. has_layer) return                          ! bare ground + sub-threshold snow: caller handles it
      snow%swe(1)         = snow%swe(1) + dm_snow + dm_rain
      snow%snow_energy(1) = snow%snow_energy(1)                                                     &
                          + dm_snow * internal_energy_ice(min(t_3ple, tair))                        &
                          + dm_rain * internal_energy_liquid(tair)           ! rain-on-snow: liquid enthalpy, refreezes later
      snow%snow_depth(1)  = snow%snow_depth(1) + dm_snow / params%rho_snow   ! rain fills pores / refreezes -> no bulk depth
      snow%nlayer         = 1_ik
   end subroutine snow_accumulate

   !----- Percolation: after the energy step, drain the FREE liquid (above the holding capacity) as a !
   !      PAIRED (mass, enthalpy) transfer to the soil (melt%melt_*), and handle full melt-out (residual !
   !      dump%). The caller adds melt_mass to infiltration availability AND melt_enth to soil_energy(1)  !
   !      (design §4e). temp/fliq are re-diagnosed from the drained (swe, energy). Mass conserved.        !
   pure subroutine snow_drain_meltwater(snow, params, melt)
      type(snow_column_t), intent(inout) :: snow
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

      call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))

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
         call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      end if
   end subroutine snow_drain_meltwater

   !----- Revert to the empty "no snow" state (nlayer = 0). ----------------------------------------!
   pure subroutine clear_layer(snow)
      type(snow_column_t), intent(inout) :: snow
      snow%swe(1)         = 0.0_wp
      snow%snow_energy(1) = 0.0_wp
      snow%snow_depth(1)  = 0.0_wp
      snow%snow_fliq(1)   = 0.0_wp
      snow%snow_temp(1)   = t_3ple
      snow%nlayer         = 0_ik
   end subroutine clear_layer

   !=======================================================================================!
   !  SNOW -- ENERGY side: surface fluxes, base conduction, backward-Euler energy advance.  !
   !=======================================================================================!

   !----- Turbulent + latent surface fluxes at a given snow-surface temperature (positive UPWARD).  !
   !      Latent uses the SAME vapour-enthalpy twin the CAS receives: because the store is           !
   !      ice-referenced, removing enthalpy_vapor from an ice layer debits sublimation (= vaporization !
   !      + fusion) automatically; a wet layer's liquid reference supplies only vaporization (§4d).    !
   pure subroutine snow_surface_fluxes(t_surf, env, h_flux, w_flux, le_flux)
      real(wp),         intent(in)  :: t_surf
      type(snow_env_t), intent(in)  :: env
      real(wp),         intent(out) :: h_flux, w_flux, le_flux
      h_flux = env%ggnet * env%rho_air * cp_air * (t_surf - env%can_temp)
      w_flux = env%ggnet * env%rho_air * (sat_specific_humidity(t_surf, env%press) - env%can_shv)
      le_flux = w_flux * enthalpy_vapor(t_surf)
   end subroutine snow_surface_fluxes

   !----- Snow-base -> soil-top conduction CONDUCTANCE [W/m2/K] (series resistance of the half snow   !
   !      layer + the top soil node). k_snow << k_soil throttles this as the pack deepens -- the       !
   !      physical decoupling that caps the winter soil surface (design §4g). ------------------------!
   pure function snow_base_conductance(snow_depth, env, params) result(gcond)
      real(wp),            intent(in) :: snow_depth
      type(snow_env_t),    intent(in) :: env
      type(snow_params_t), intent(in) :: params
      real(wp) :: gcond, r_series
      r_series = 0.5_wp * max(snow_depth, tiny_num) / max(params%k_snow, tiny_num)                  &
               + max(env%dz_soil_top, tiny_num) / max(env%k_soil_top, tiny_num)
      gcond    = 1.0_wp / max(r_series, tiny_num)
   end function snow_base_conductance

   !----- Advance the snow store one dt by a linearized-implicit (backward-Euler) energy step.        !
   !      Prognostic = extensive internal energy (J/m2) + water-equivalent mass (swe); temperature +   !
   !      liquid fraction are read-offs. The store update captures WARMING and MELT together (the       !
   !      inverter re-partitions), the emission response is made consistent with the linearization slope !
   !      (so the update stays bounded -- the wood-store lesson), and sublimation debits vapour enthalpy  !
   !      while removing mass. Reports the turbulent/conduction fluxes for the CAS + soil coupling.       !
   pure subroutine snow_energy_step(snow, env, params, dt, area_frac, flux)
      type(snow_column_t), intent(inout) :: snow
      type(snow_env_t),    intent(in)    :: env
      type(snow_params_t), intent(in)    :: params
      real(wp),            intent(in)    :: dt
      real(wp),            intent(in)    :: area_frac    !< snow-cover fraction: scales ALL boundary exchange (sub-column)
      type(snow_flux_t),   intent(out)   :: flux
      real(wp) :: t_n, fliq_n, cap, e_old, swe0, gcond, af
      real(wp) :: h_n, w_n, le_n, g_n, r_n, drdt, dqsatdt, t_star
      real(wp) :: h_s, w_s, le_s, g_s, emiss_corr, net_flux

      flux = snow_flux_t()
      af = max(0.0_wp, min(1.0_wp, area_frac))
      if (snow%nlayer == 0_ik .or. snow%swe(1) <= params%tiny_snow_mass .or. af <= 0.0_wp) then
         flux%t_surf = env%t_soil_top                              ! no snow surface (or zero cover) to balance
         return
      end if

      swe0  = snow%swe(1)
      e_old = snow%snow_energy(1)
      flux%snowfac = af
      call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, t_n, fliq_n)
      gcond = snow_base_conductance(snow%snow_depth(1), env, params)

      !----- Fluxes + linearization slope at T^n (all drdt terms <= 0). ----------------------------!
      call snow_surface_fluxes(t_n, env, h_n, w_n, le_n)
      g_n     = gcond * (t_n - env%t_soil_top)
      r_n     = env%abs_sw + env%abs_lw - h_n - le_n - g_n
      dqsatdt = sat_specific_humidity_temp_deriv(t_n, env%press)
      drdt    = -4.0_wp * params%snow_emiss * stefan * t_n ** 3                                     &
              - env%ggnet * env%rho_air * cp_air                                                    &
              - latent_heat_vap * env%rho_air * env%ggnet * dqsatdt                                 &
              - gcond
      cap     = snow%swe(1) * (fliq_n * cp_liq + (1.0_wp - fliq_n) * cp_ice)   ! sensible heat capacity [J/m2/K]

      !----- Plateau-aware implicit surface temperature. On the melt plateau (0<fliq<1) the           !
      !      temperature is PINNED at t_3ple; a warming/cooling step that would cross t_3ple is pinned  !
      !      there too (the layer starts to melt/freeze -- the excess energy goes to fliq via the        !
      !      inverter, not to temperature). t_star stays within one Newton step of t_n, so bounded.      !
      !      Only the AREA-FRACTION af of the surface exchanges (sub-column), so af scales the flux AND   !
      !      its slope in the implicit estimate (cap = full pack heat capacity is NOT scaled): af -> 0    !
      !      makes the store a no-op (a thin patchy pack barely exchanges -> stable + continuous).        !
      if (fliq_n > 0.0_wp .and. fliq_n < 1.0_wp) then
         t_star = t_3ple
      else
         t_star = t_n + af * r_n / (cap / dt - af * drdt)
         if ((t_n < t_3ple .and. t_star > t_3ple) .or. (t_n > t_3ple .and. t_star < t_3ple)) t_star = t_3ple
      end if

      !----- Fluxes at t_star; cap sublimation to the available mass (af*w_s*dt <= swe; deposition free). !
      call snow_surface_fluxes(t_star, env, h_s, w_s, le_s)
      g_s = gcond * (t_star - env%t_soil_top)
      if (w_s > 0.0_wp) w_s = min(w_s, snow%swe(1) / (af * dt))
      le_s = w_s * enthalpy_vapor(t_star)

      !----- Emission response consistent with drdt's -4*eps*sigma*T^3 term: abs_lw is NET at T^n, so   !
      !      the ADDED emission as the surface moves to t_star is eps*sigma*(t_star^4 - t_n^4). Including !
      !      it makes the energy update match the linearization (bounded, no wood-style overshoot).      !
      emiss_corr = params%snow_emiss * stefan * (t_star ** 4 - t_n ** 4)
      net_flux   = af * (env%abs_sw + env%abs_lw - h_s - le_s - g_s - emiss_corr)   ! sub-column area weight

      !----- ENERGY + MASS update. The inverter downstream re-partitions warming vs melt; sublimation    !
      !      removes af-scaled mass and its vapour enthalpy (net_flux already debits af*le_s). ---------!
      snow%snow_energy(1) = e_old + dt * net_flux
      snow%swe(1)         = swe0 - af * w_s * dt
      snow%snow_depth(1)  = max(0.0_wp, snow%swe(1)) / params%rho_snow

      if (snow%swe(1) > params%tiny_snow_mass) then
         call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      else
         snow%snow_temp(1) = t_3ple                               ! vanished; drain kernel dumps residual
         snow%snow_fliq(1) = 0.0_wp
      end if

      !----- Report the turbulent + conduction fluxes for the CAS + soil coupling. h_snow + le_snow -> !
      !      CAS enthalpy; w_flux -> CAS vapour; g_base -> soil top BC. The residual closes by           !
      !      construction: (energy change) - dt*(net_flux) = 0. -------------------------------------!
      flux%t_surf  = merge(snow%snow_temp(1), t_star, snow%swe(1) > params%tiny_snow_mass)
      flux%fliq    = snow%snow_fliq(1)
      flux%h_snow  = af * h_s
      flux%le_snow = af * le_s
      flux%w_flux  = af * w_s
      flux%g_base  = af * g_s
      flux%rnet    = af * (env%abs_sw + env%abs_lw - emiss_corr)   ! snow's area-weighted net radiation into the store
      flux%energy_resid = (snow%snow_energy(1) - e_old) - dt * net_flux
   end subroutine snow_energy_step

end module meds_ground_biophysics
