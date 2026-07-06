!==========================================================================================!
! meds_surface_energy -- the three per-store SURFACE energy kernels of the energy family        !
! (design MEDS_ENERGY_BALANCE_DESIGN.md, 4a/4b/4d): the per-cohort leaf/wood energy balance, the  !
! per-patch ground/surface balance, and the per-patch canopy-air-space enthalpy update. All         !
! STATELESS: each takes the sibling stores' temperatures as forced inputs (the coupled leaf<->CAS   !
! <->ground<->soil fixed point is deferred to P3). Prognostic INTERNAL ENERGY / specific enthalpy;   !
! temperature diagnosed via the meds_thermo inverter. P1 is liquid-only, one linearized step.        !
!==========================================================================================!
module meds_surface_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : pi, stefan, cp_air, cp_liq, latent_heat_vap, tiny_num
   use meds_biophysics_types, only : leaf_energy_env_t, leaf_energy_flux_t, veg_thermal_params_t
   use meds_thermo,           only : sat_specific_humidity, sat_vapor_pressure,                &
                                     d_sat_vapor_pressure_dt, enthalpy_vapor, uext_to_temp,     &
                                     cas_temp_of_enthalpy
   implicit none
   private

   public :: veg_energy_balance, ground_surface_balance, canopy_air_update

   real(wp), parameter :: LEAF_MAXWHC = 0.11_wp     !< [kg/m2 leaf] film-holding capacity (wetted fraction)

contains

   !---------------------------------------------------------------------------------------!
   ! Leaf OR wood cohort energy balance (design 4a). One L-stable linearized BE step on the  !
   ! prognostic store_energy; fluxes go OUT to the CAS as twins. is_leaf toggles transpiration !
   ! and the sensible area basis (2*LAI flat plates vs pi*WAI cylinders).                       !
   !---------------------------------------------------------------------------------------!
   pure subroutine veg_energy_balance(store_energy, env, tparams, dt, is_leaf, flux)
      real(wp),                   intent(inout) :: store_energy
      type(leaf_energy_env_t),    intent(in)    :: env
      type(veg_thermal_params_t), intent(in)    :: tparams
      real(wp),                   intent(in)    :: dt
      logical,                    intent(in)    :: is_leaf
      type(leaf_energy_flux_t),   intent(out)   :: flux

      real(wp) :: t_n, fliq_n, cap, area_h, esat, dqsatdt, drdt, t_star, r_n, r_star, e_old
      real(wp) :: h_n, w_n, qw_n, tr_n, qt_n, gv_n
      real(wp) :: h_s, w_s, qw_s, tr_s, qt_s, gv_s

      area_h = tparams%effarea_heat ; if (.not. is_leaf) area_h = pi          ! flat plates vs cylinders
      cap    = env%dry_hcap + env%wmass * cp_liq                             ! [J/m2/K] total heat capacity
      e_old  = store_energy
      call uext_to_temp(store_energy, env%wmass, env%dry_hcap, t_n, fliq_n)

      !----- Fluxes + linearization slope at T^n. ------------------------------------------!
      call veg_surface_fluxes(t_n, env, tparams, is_leaf, area_h, h_n, w_n, qw_n, tr_n, qt_n, gv_n)
      r_n  = env%abs_sw + env%abs_lw - h_n - qw_n - qt_n                     ! [W/m2] net into store
      esat = sat_vapor_pressure(t_n)
      dqsatdt = 0.622_wp * env%press / max((env%press - 0.378_wp * esat) ** 2, tiny_num)      &
                * d_sat_vapor_pressure_dt(t_n)
      drdt = -8.0_wp * tparams%leaf_emiss * stefan * t_n ** 3 * env%area_index                &
             - area_h * env%area_index * env%gbh * env%rho_air * cp_air                       &
             - latent_heat_vap * env%rho_air * gv_n * dqsatdt                 ! all terms <= 0

      !----- Implicit (linearized) temperature estimate, then conservative energy update. ---!
      t_star = t_n + r_n / (cap / dt - drdt)
      call veg_surface_fluxes(t_star, env, tparams, is_leaf, area_h, h_s, w_s, qw_s, tr_s, qt_s, gv_s)
      r_star = env%abs_sw + env%abs_lw - h_s - qw_s - qt_s
      store_energy = store_energy + dt * r_star

      call uext_to_temp(store_energy, env%wmass, env%dry_hcap, flux%temp, flux%fliq)
      flux%h_flux   = h_s ; flux%qw_flux = qw_s ; flux%q_transp = qt_s
      flux%w_flux   = w_s ; flux%transp  = tr_s
      flux%energy_resid = (store_energy - e_old) - dt * r_star               ! = 0 by construction
   end subroutine veg_energy_balance

   !----- Sensible + film-evaporation + transpiration fluxes at temperature t (helper). -----!
   pure subroutine veg_surface_fluxes(t, env, tparams, is_leaf, area_h, h_flux, w_flux, qw_flux, &
                                      transp, q_transp, g_vapor)
      real(wp),                   intent(in)  :: t, area_h
      type(leaf_energy_env_t),    intent(in)  :: env
      type(veg_thermal_params_t), intent(in)  :: tparams
      logical,                    intent(in)  :: is_leaf
      real(wp),                   intent(out) :: h_flux, w_flux, qw_flux, transp, q_transp, g_vapor
      real(wp) :: grad, w_max, sigma_w, sigma_eff, g_series, g_ev, g_tr
      grad   = sat_specific_humidity(t, env%press) - env%can_shv
      w_max  = LEAF_MAXWHC * max(env%area_index, tiny_num)
      sigma_w = 0.0_wp
      if (env%leaf_water > 0.0_wp) sigma_w = min(1.0_wp, (env%leaf_water / w_max) ** (2.0_wp/3.0_wp))
      sigma_eff = sigma_w ; if (grad < 0.0_wp) sigma_eff = 1.0_wp            ! dew wets the full surface
      !----- Sensible. --------------------------------------------------------------------!
      h_flux = area_h * env%area_index * env%gbh * env%rho_air * cp_air * (t - env%can_temp)
      !----- Interception-film evaporation (both leaf and wood). ---------------------------!
      g_ev   = tparams%effarea_evap * env%area_index * env%gbw * sigma_eff
      w_flux  = g_ev * env%rho_air * grad
      qw_flux = w_flux * enthalpy_vapor(t)
      !----- Transpiration (leaf only; boundary layer in series with stomata). -------------!
      transp = 0.0_wp ; g_tr = 0.0_wp
      if (is_leaf .and. env%gbw + env%gsw > tiny_num) then
         g_series = env%gbw * env%gsw / (env%gbw + env%gsw)
         g_tr     = tparams%effarea_transp * env%area_index * g_series * env%fs_open
         transp   = g_tr * env%rho_air * grad
      end if
      q_transp = transp * enthalpy_vapor(t)
      g_vapor  = g_ev + g_tr                                                 ! total vapour conductance (slope)
   end subroutine veg_surface_fluxes

   !---------------------------------------------------------------------------------------!
   ! Ground/surface energy balance (design 4d) -- INSTANTANEOUS, no state advance. Returns    !
   ! G_top = Rn − H − LE (net into the soil surface) + the sensible/latent diagnostics.       !
   ! t_ground = soil_temp(1) is supplied by the soil kernel. env carries the ground conductance !
   ! (gbh = ggnet), vapour conductance (gbw), CAS state, and absorbed radiation (abs_sw/lw).    !
   !---------------------------------------------------------------------------------------!
   pure subroutine ground_surface_balance(t_ground, env, g_top, h_ground, le_ground)
      real(wp),                intent(in)  :: t_ground
      type(leaf_energy_env_t), intent(in)  :: env
      real(wp),                intent(out) :: g_top, h_ground, le_ground
      real(wp) :: rn, w_flux_gc
      rn        = env%abs_sw + env%abs_lw
      h_ground  = env%gbh * env%rho_air * cp_air * (t_ground - env%can_temp)
      w_flux_gc = env%gbw * env%rho_air * (sat_specific_humidity(t_ground, env%press) - env%can_shv)
      le_ground = w_flux_gc * enthalpy_vapor(t_ground)                       ! same enthalpy twin the CAS receives
      g_top     = rn - h_ground - le_ground
   end subroutine ground_surface_balance

   !---------------------------------------------------------------------------------------!
   ! Canopy-air-space enthalpy + humidity update (design 4b). Advances BOTH prognostic twins  !
   ! one step from the summed cohort/ground fluxes and the atmospheric exchange (implicit in    !
   ! the atm term for L-stability). Returns the closed-budget residual (~0).                    !
   !---------------------------------------------------------------------------------------!
   pure subroutine canopy_air_update(cas_enthalpy, cas_shv, cas_temp, can_depth,               &
                                     coh_h_flux, coh_qw_flux, coh_w_flux, coh_transp,           &
                                     ground_h_flux, ground_qw_flux, ground_w_flux, dew,         &
                                     ustar, enthalpy_atm, w_flux_ac, rho_air, dt, resid)
      real(wp), intent(inout) :: cas_enthalpy, cas_shv, cas_temp
      real(wp), intent(in)    :: can_depth
      real(wp), intent(in)    :: coh_h_flux, coh_qw_flux, coh_w_flux, coh_transp
      real(wp), intent(in)    :: ground_h_flux, ground_qw_flux, ground_w_flux, dew
      real(wp), intent(in)    :: ustar, enthalpy_atm, w_flux_ac, rho_air, dt
      real(wp), intent(out)   :: resid
      real(wp) :: wcapcan, wci, f_sens, gatm, enth_new, shv_new
      wcapcan = rho_air * can_depth                                         ! [kg/m2] CAS air mass per ground area
      wci     = 1.0_wp / max(wcapcan, tiny_num)
      f_sens  = coh_h_flux + coh_qw_flux + ground_h_flux + ground_qw_flux   ! [W/m2] into CAS from surfaces
      gatm    = rho_air * ustar                                            ! [kg/m2/s] atm<->CAS exchange
      !----- Enthalpy: implicit in the atmospheric-exchange term. --------------------------!
      enth_new = (cas_enthalpy + dt * wci * (f_sens + gatm * enthalpy_atm)) / (1.0_wp + dt * wci * gatm)
      !----- Specific humidity twin. -------------------------------------------------------!
      shv_new  = cas_shv + dt * wci * (coh_w_flux + coh_transp + ground_w_flux - dew + w_flux_ac)
      resid    = wcapcan * (enth_new - cas_enthalpy)                                            &
                 - dt * (f_sens + gatm * (enthalpy_atm - enth_new))         ! = 0 by construction
      cas_enthalpy = enth_new
      cas_shv      = shv_new
      cas_temp     = cas_temp_of_enthalpy(enth_new, shv_new)
   end subroutine canopy_air_update

end module meds_surface_energy
