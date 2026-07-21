!==========================================================================================!
! meds_vegetation_biophysics -- the stateless per-cohort VEGETATION-surface kernels: the        !
! leaf/wood tissue energy balance and the canopy interception film they share. Coupled through   !
! the wetted-fraction film (`leaf_water` -> sigma_w -> leaf latent flux), so they live together.   !
!                                                                                          !
!   * veg_energy_step_implicit -- one L-stable linearized BE step on the prognostic tissue energy  !
!                                 (leaf OR wood); fluxes go OUT to the CAS as twins.                 !
!   * intercept_canopy_layer   -- per-cohort capacity-limited interception bucket (Beer fraction);   !
!                                 the caller sweeps height-sorted cohorts top->bottom.                !
!==========================================================================================!
module meds_vegetation_biophysics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : pi, tiny_num, stefan, cp_air, cp_liq, latent_heat_vap
   use meds_therm_lib,        only : uext_to_temp, sat_specific_humidity,                        &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor
   use meds_biophysics_types, only : leaf_energy_env_t, leaf_energy_flux_t, veg_thermal_params_t
   implicit none
   private

   public :: veg_energy_step_implicit, intercept_canopy_layer

   real(wp), parameter :: LEAF_MAXWHC = 0.11_wp     !< [kg/m2 leaf] film-holding capacity (wetted fraction)

contains

   !---------------------------------------------------------------------------------------!
   ! Leaf OR wood cohort energy balance (design 4a). One L-stable linearized BE step on the  !
   ! prognostic store_energy; fluxes go OUT to the CAS as twins. is_leaf toggles transpiration !
   ! and the sensible area basis (2*LAI flat plates vs pi*WAI cylinders).                       !
   !---------------------------------------------------------------------------------------!
   pure subroutine veg_energy_step_implicit(store_energy, env, tparams, dt, is_leaf, flux)
      real(wp),                   intent(inout) :: store_energy
      type(leaf_energy_env_t),    intent(in)    :: env
      type(veg_thermal_params_t), intent(in)    :: tparams
      real(wp),                   intent(in)    :: dt
      logical,                    intent(in)    :: is_leaf
      type(leaf_energy_flux_t),   intent(out)   :: flux

      real(wp) :: t_n, fliq_n, cap, area_h, dqsatdt, drdt, t_star, r_n, e_old
      real(wp) :: delta_e, turb_star, turb_avg, scale
      real(wp) :: h_n, w_n, qw_n, tr_n, qt_n, gv_n
      real(wp) :: h_s, w_s, qw_s, tr_s, qt_s, gv_s

      area_h = tparams%effarea_heat
      if (.not. is_leaf) area_h = pi                                        ! flat plates vs cylinders
      cap    = env%dry_hcap + env%wmass * cp_liq                             ! [J/m2/K] total heat capacity
      e_old  = store_energy
      call uext_to_temp(store_energy, env%wmass, env%dry_hcap, t_n, fliq_n)

      !----- Fluxes + linearization slope at T^n. ------------------------------------------!
      call veg_surface_fluxes(t_n, env, tparams, is_leaf, area_h, h_n, w_n, qw_n, tr_n, qt_n, gv_n)
      r_n  = env%abs_sw + env%abs_lw - h_n - qw_n - qt_n                     ! [W/m2] net into store
      dqsatdt = sat_specific_humidity_temp_deriv(t_n, env%press)
      drdt = -8.0_wp * tparams%leaf_emiss * stefan * t_n ** 3 * env%area_index                &
             - area_h * env%area_index * env%gbh * env%rho_air * cp_air                       &
             - latent_heat_vap * env%rho_air * gv_n * dqsatdt                 ! all terms <= 0

      !----- Implicit (linearized) temperature estimate. t_star stays within one Newton step of  !
      !      t_n (drdt<=0 => cap/dt-drdt>0), so the energy update below is BOUNDED even when the   !
      !      store is stiff (tau = cap/|drdt| << dt) and far from equilibrium. Applying the        !
      !      ENDPOINT flux dt*r_star over the whole step instead over-removes energy there and can  !
      !      drive T negative (a young-stand wood cohort: tau~35 s vs dt_fast=1800 s).             !
      t_star  = t_n + r_n / (cap / dt - drdt)
      delta_e = cap * (t_star - t_n)                                         ! [J/m2] L-stable energy change
      store_energy = e_old + delta_e
      call uext_to_temp(store_energy, env%wmass, env%dry_hcap, flux%temp, flux%fliq)

      !----- Conservative flux report: the radiation the store did NOT keep leaves as turbulent + !
      !      latent flux (step-AVERAGED), split by the instantaneous shares at t_star. For a linear !
      !      flux delta_e == dt*r_star, so this reduces to the endpoint report (scale = 1).         !
      call veg_surface_fluxes(t_star, env, tparams, is_leaf, area_h, h_s, w_s, qw_s, tr_s, qt_s, gv_s)
      turb_star = h_s + qw_s + qt_s
      turb_avg  = env%abs_sw + env%abs_lw - delta_e / dt
      scale     = 1.0_wp
      if (abs(turb_star) > tiny_num) scale = turb_avg / turb_star
      flux%h_flux   = h_s * scale
      flux%qw_flux  = qw_s * scale
      flux%q_transp = qt_s * scale
      flux%w_flux   = w_s * scale
      flux%transp   = tr_s * scale
      flux%energy_resid = delta_e - dt * (env%abs_sw + env%abs_lw - flux%h_flux                 &
                          - flux%qw_flux - flux%q_transp)                    ! = 0 by construction
   end subroutine veg_energy_step_implicit

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
      sigma_eff = sigma_w
      if (grad < 0.0_wp) sigma_eff = 1.0_wp                                 ! dew wets the full surface
      !----- Sensible. --------------------------------------------------------------------!
      h_flux = area_h * env%area_index * env%gbh * env%rho_air * cp_air * (t - env%can_temp)
      !----- Interception-film evaporation (both leaf and wood). ---------------------------!
      g_ev   = tparams%effarea_evap * env%area_index * env%gbw * sigma_eff
      w_flux  = g_ev * env%rho_air * grad
      qw_flux = w_flux * enthalpy_vapor(t)
      !----- Transpiration (leaf only; boundary layer in series with stomata). -------------!
      transp = 0.0_wp
      g_tr   = 0.0_wp
      if (is_leaf .and. env%gbw + env%gsw > tiny_num) then
         g_series = env%gbw * env%gsw / (env%gbw + env%gsw)
         g_tr     = tparams%effarea_transp * env%area_index * g_series * env%fs_open
         transp   = g_tr * env%rho_air * grad
      end if
      q_transp = transp * enthalpy_vapor(t)
      g_vapor  = g_ev + g_tr                                                 ! total vapour conductance (slope)
   end subroutine veg_surface_fluxes

   !---------------------------------------------------------------------------------------!
   ! Per-cohort canopy interception (design 3c). One canopy layer: capacity-limited bucket    !
   ! with a Beer interception fraction; the caller sweeps the height-sorted cohorts top->bottom !
   ! feeding each `throughfall` as the next cohort's `rain_above`. `leaf_water` is the per-cohort !
   ! prognostic film [kg/m2 ground]; `sigma_w` (wetted fraction) is exported for the CAS-space     !
   ! evaporation. Not `elemental`: the cascade is a sequential recurrence.                          !
   !---------------------------------------------------------------------------------------!
   pure subroutine intercept_canopy_layer(leaf_water, rain_above, lai, sai, e_canopy, dt,     &
                                          dewmx, k_int, alpha_pi, throughfall, drip, sigma_w)
      real(wp), intent(inout) :: leaf_water
      real(wp), intent(in)    :: rain_above, lai, sai, e_canopy, dt, dewmx, k_int, alpha_pi
      real(wp), intent(out)   :: throughfall, drip, sigma_w
      real(wp) :: pai, f_pi, w_max, q_grab, room, q_intr
      pai    = lai + sai
      f_pi   = alpha_pi * (1.0_wp - exp(-k_int * pai))
      w_max  = dewmx * pai
      q_grab = f_pi * rain_above
      room   = max(0.0_wp, w_max - leaf_water) / dt
      q_intr = min(q_grab, room + e_canopy)                      ! bounded by capacity + evap headroom
      leaf_water  = min(max(leaf_water + (q_intr - e_canopy) * dt, 0.0_wp), w_max)
      drip        = max(0.0_wp, q_grab - q_intr)
      throughfall = (rain_above - q_grab) + drip                 ! gap-throughfall + drip -> below
      if (w_max > tiny_num) then
         sigma_w = min(1.0_wp, (leaf_water / w_max) ** (2.0_wp / 3.0_wp))
      else
         sigma_w = 0.0_wp
      end if
   end subroutine intercept_canopy_layer

end module meds_vegetation_biophysics
