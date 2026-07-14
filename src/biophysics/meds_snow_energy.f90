!==========================================================================================!
! meds_snow_energy -- the ENERGY side of the MEDS temporary-surface-water / snow store (design    !
! MEDS_SNOW_DESIGN.md P0, single bulk layer). Stateless per-store kernels: the snow-SURFACE energy  !
! balance (net SW with snow albedo + net LW − sensible − latent[sublimation/evap] − base            !
! conduction), the snow-BASE conduction that becomes the soil top BC, and a backward-Euler energy    !
! advance of the store. Snow temperature + liquid fraction are a read-off of the shared inverter     !
! uext_to_temp (dry_hcap = 0), so MELT/refreeze is the internal-energy plateau, captured for free.    !
! Links meds_shared ONLY. Advances snow_energy (extensive J/m2) directly; the plateau handles phase.  !
!==========================================================================================!
module meds_snow_energy
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan,            &
                                      latent_heat_vap
   use meds_thermo,            only : uext_to_temp, sat_specific_humidity, sat_vapor_pressure,     &
                                      d_sat_vapor_pressure_dt, enthalpy_vapor
   use meds_column_state_types, only : snow_column_t
   use meds_biophysics_types,  only : snow_params_t, snow_env_t, snow_flux_t
   use meds_snow_mass,         only : snow_cover_fraction
   implicit none
   private

   public :: snow_surface_fluxes, snow_base_conductance, snow_energy_step

contains

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
      real(wp) :: h_n, w_n, le_n, g_n, r_n, drdt, dqsatdt, esat, t_star
      real(wp) :: h_s, w_s, le_s, g_s, emiss_corr, net_flux

      flux = snow_flux_t()
      af = max(0.0_wp, min(1.0_wp, area_frac))
      if (snow%nlayer == 0_ik .or. snow%swe(1) <= params%tiny_snow_mass .or. af <= 0.0_wp) then
         flux%t_surf = env%t_soil_top ; return                     ! no snow surface (or zero cover) to balance
      end if

      swe0  = snow%swe(1) ; e_old = snow%snow_energy(1)
      flux%snowfac = af
      call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, t_n, fliq_n)
      gcond = snow_base_conductance(snow%snow_depth(1), env, params)

      !----- Fluxes + linearization slope at T^n (all drdt terms <= 0). ----------------------------!
      call snow_surface_fluxes(t_n, env, h_n, w_n, le_n)
      g_n     = gcond * (t_n - env%t_soil_top)
      r_n     = env%abs_sw + env%abs_lw - h_n - le_n - g_n
      esat    = sat_vapor_pressure(t_n)
      dqsatdt = 0.622_wp * env%press / max((env%press - 0.378_wp * esat) ** 2, tiny_num)            &
                * d_sat_vapor_pressure_dt(t_n)
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
         snow%snow_temp(1) = t_3ple ; snow%snow_fliq(1) = 0.0_wp   ! vanished; drain kernel dumps residual
      end if

      !----- Report the turbulent + conduction fluxes for the CAS + soil coupling. h_snow + le_snow -> !
      !      CAS enthalpy; w_flux -> CAS vapour; g_base -> soil top BC. The residual closes by           !
      !      construction: (energy change) - dt*(net_flux) = 0. -------------------------------------!
      flux%t_surf  = merge(snow%snow_temp(1), t_star, snow%swe(1) > params%tiny_snow_mass)
      flux%fliq    = snow%snow_fliq(1)
      flux%h_snow  = af * h_s ; flux%le_snow = af * le_s ; flux%w_flux = af * w_s ; flux%g_base = af * g_s
      flux%rnet    = af * (env%abs_sw + env%abs_lw - emiss_corr)   ! snow's area-weighted net radiation into the store
      flux%energy_resid = (snow%snow_energy(1) - e_old) - dt * net_flux
   end subroutine snow_energy_step

end module meds_snow_energy
