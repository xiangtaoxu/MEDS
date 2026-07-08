!==========================================================================================!
! meds_column_energy -- THE stateless energy-balance kernels for the WHOLE soil-veg-air column   !
! (design MEDS_ENERGY_BALANCE_DESIGN.md). "Column" here is the full vertical land-surface column,  !
! not just soil: the internal-energy balance of the four coupled stores that carry temperature --   !
!   * soil thermal column   -- `soil_energy_flux`     (per patch, per layer)                          !
!   * ground / surface skin  -- `ground_surface_balance` (per patch; the top soil layer)               !
!   * canopy air space        -- `canopy_air_update`    (per patch)                                      !
!   * leaf / wood tissue       -- `veg_energy_balance`   (per cohort, leaf OR wood)                        !
!                                                                                          !
! Prognostic INTERNAL ENERGY / specific enthalpy (phase-safe); temperature diagnosed via the        !
! meds_thermo inverter. All STATELESS: each kernel takes the sibling stores' temperatures as forced   !
! inputs -- the coupled leaf<->CAS<->ground<->soil fixed point is deferred to P3. P1 is liquid-only,   !
! one implicit/linearized step per call. The soil sweep reuses meds_soil_solver (Thomas) + the          !
! negative-z geometry; the inner `soil_heat_be_step` is bare-array + device-eligible (growth_step).      !
!==========================================================================================!
module meds_column_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, tiny_num, pi, stefan, cp_air, cp_liq, latent_heat_vap
   use meds_biophysics_types, only : soil_energy_column_t, energy_forcing_t, soil_thermal_params_t, &
                                     soil_params_t, energy_opts_t, energy_flux_t, n_soil_layer_max, &
                                     ENERGY_PHASE_OFF, leaf_energy_env_t, leaf_energy_flux_t,        &
                                     veg_thermal_params_t
   use meds_thermo,           only : uext_to_temp, internal_energy_liquid, sat_specific_humidity,   &
                                     sat_vapor_pressure, d_sat_vapor_pressure_dt, enthalpy_vapor,    &
                                     cas_temp_of_enthalpy
   use meds_soil_thermal,     only : soil_thermal_cond, soil_heat_cap_vol
   use meds_numerics,         only : thomas_solve
   implicit none
   private

   public :: soil_energy_flux, veg_energy_balance, ground_surface_balance, canopy_air_update

   real(wp), parameter :: LEAF_MAXWHC = 0.11_wp     !< [kg/m2 leaf] film-holding capacity (wetted fraction)

contains

   !=======================================================================================!
   !  Soil thermal column                                                                   !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! The soil seam: advance one patch's soil thermal column over dt.                        !
   !---------------------------------------------------------------------------------------!
   subroutine soil_energy_flux(col, forcing, therm, soil, opts, dt, flux)
      type(soil_energy_column_t),  intent(inout) :: col
      type(energy_forcing_t),      intent(in)    :: forcing
      type(soil_thermal_params_t), intent(in)    :: therm
      type(soil_params_t),         intent(in)    :: soil
      type(energy_opts_t),         intent(in)    :: opts
      real(wp),                    intent(in)    :: dt
      type(energy_flux_t),         intent(out)   :: flux

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max) :: t_n, fl_n, kappa, c_eff, q_src, t_new, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: fliq_use, wmass, e0, e1, div

      n = soil%n_active

      !----- Invert energy -> (temp, fliq) and node-wise thermal properties at state^n. ----!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o                          ! [kg/m3] water mass per volume
         call uext_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k), t_n(k), fl_n(k))
         fliq_use = fl_n(k)
         if (opts%phase_change == ENERGY_PHASE_OFF) fliq_use = 1.0_wp     ! liquid-only in P1
         kappa(k) = soil_thermal_cond(forcing%soil_water(k), fliq_use, soil%theta_sat(k),     &
                                      therm%soil_solid_conductivity(k), therm%soil_dry_conductivity(k))
         c_eff(k) = soil_heat_cap_vol(forcing%soil_water(k), fliq_use, therm%soil_dry_heat_capacity(k))
         q_src(k) = -forcing%root_heat_sink(k) / soil%dz(k)               ! [W/m3] source (sink is negative)
      end do

      !----- Implicit BE conduction solve for temperature^{n+1}. ----------------------------!
      call soil_heat_be_step(t_n, soil%dz, soil%dz_node, kappa, c_eff, q_src, forcing%g_top,   &
                             forcing%geothermal, dt, n, t_new)

      !----- Conservative energy update: conductive faces from T^{n+1}, advective upwind. ---!
      hf(0)  = -forcing%g_top                                             ! top face (positive up)
      qwf(0) = 0.0_wp
      do k = 1_ik, n - 1_ik
         kf(k)  = (soil%dz(k) + soil%dz(k+1)) / (soil%dz(k) / kappa(k) + soil%dz(k+1) / kappa(k+1))
         hf(k)  = -kf(k) * (t_new(k) - t_new(k+1)) / soil%dz_node(k)
         !----- Upwind the liquid enthalpy on the SOURCE layer. w_flux is UPWARD-positive, so     !
         !      <= 0 (downward flow) draws from layer k (above face k); > 0 (upward flow) from      !
         !      k+1 (below). The sign is carried by w_flux (no extra minus), matching hf and ED2    !
         !      rk4_derivs qw_flux_g. (This convention was verified correct -- BUG6 was refuted.)   !
         if (forcing%w_flux(k) <= 0.0_wp) then
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_new(k))
         else
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_new(k+1))
         end if
      end do
      hf(n)  = forcing%geothermal                                        ! bottom geothermal (positive up)
      qwf(n) = 0.0_wp

      e0 = 0.0_wp ; e1 = 0.0_wp
      do k = 1_ik, n
         e0  = e0 + col%soil_energy(k) * soil%dz(k)
         div = (hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))
         col%soil_energy(k) = col%soil_energy(k) + dt / soil%dz(k) * div + dt * q_src(k)
         e1  = e1 + col%soil_energy(k) * soil%dz(k)
      end do

      !----- Re-diagnose temperature + liquid fraction from the committed energy. ----------!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o
         call uext_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k),         &
                           col%soil_temp(k), col%soil_fliq(k))
      end do

      !----- Diagnostics + closed energy budget. -------------------------------------------!
      flux%ground_heat  = forcing%g_top
      flux%bottom_heat  = -forcing%geothermal
      flux%energy_resid = (e1 - e0) - dt * (forcing%g_top - flux%bottom_heat                   &
                          - sum(forcing%root_heat_sink(1:n)))
      flux%nsub      = 1_ik
      flux%converged = .true.
      if (opts%debug_error .and. abs(flux%energy_resid) > opts%atol * sum(c_eff(1:n)*soil%dz(1:n))) then
         error stop 'soil_energy_flux: energy budget did not close'
      end if
   end subroutine soil_energy_flux

   !---------------------------------------------------------------------------------------!
   ! Bare-array device-eligible inner step: implicit BE conduction (dz-weighted harmonic     !
   ! face kappa) + explicit source q_src, top Neumann g_top, bottom geothermal geo. Returns   !
   ! temperature^{n+1}. (Advection is applied in the conservative energy update, not here.)   !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_heat_be_step(t_n, dz, dz_node, kappa, c_eff, q_src, g_top, geo, dt, nzg, t_new)
      integer(ik), intent(in)  :: nzg
      real(wp),    intent(in)  :: t_n(n_soil_layer_max), dz(n_soil_layer_max)
      real(wp),    intent(in)  :: dz_node(n_soil_layer_max), kappa(n_soil_layer_max)
      real(wp),    intent(in)  :: c_eff(n_soil_layer_max), q_src(n_soil_layer_max)
      real(wp),    intent(in)  :: g_top, geo, dt
      real(wp),    intent(out) :: t_new(n_soil_layer_max)
      real(wp) :: kf(n_soil_layer_max), a(n_soil_layer_max), b(n_soil_layer_max)
      real(wp) :: c(n_soil_layer_max), r(n_soil_layer_max)
      integer(ik) :: k
      do k = 1_ik, nzg - 1_ik
         kf(k) = (dz(k) + dz(k+1)) / (dz(k) / kappa(k) + dz(k+1) / kappa(k+1))     ! series-resistor face
      end do
      do k = 1_ik, nzg
         if (k >= 2_ik) then ; a(k) = -kf(k-1) / dz_node(k-1) ; else ; a(k) = 0.0_wp ; end if
         if (k <= nzg-1_ik) then ; c(k) = -kf(k) / dz_node(k) ; else ; c(k) = 0.0_wp ; end if
         b(k) = c_eff(k) * dz(k) / dt - a(k) - c(k)
         r(k) = c_eff(k) * dz(k) / dt * t_n(k) + dz(k) * q_src(k)
      end do
      r(1)   = r(1)   + g_top                                            ! top Neumann surface flux
      r(nzg) = r(nzg) + geo                                              ! bottom geothermal source
      call thomas_solve(a, b, c, r, t_new, nzg)
   end subroutine soil_heat_be_step

   !=======================================================================================!
   !  Leaf / wood tissue                                                                    !
   !=======================================================================================!

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

   !=======================================================================================!
   !  Ground / surface skin  +  Canopy air space                                            !
   !=======================================================================================!

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
                                     ustar, temp1, enthalpy_atm, w_flux_ac, rho_air, dt, resid)
      real(wp), intent(inout) :: cas_enthalpy, cas_shv, cas_temp
      real(wp), intent(in)    :: can_depth
      real(wp), intent(in)    :: coh_h_flux, coh_qw_flux, coh_w_flux, coh_transp
      real(wp), intent(in)    :: ground_h_flux, ground_qw_flux, ground_w_flux, dew
      real(wp), intent(in)    :: ustar, temp1, enthalpy_atm, w_flux_ac, rho_air, dt
      real(wp), intent(out)   :: resid
      real(wp) :: wcapcan, wci, f_sens, gatm, enth_new, shv_new
      wcapcan = rho_air * can_depth                                         ! [kg/m2] CAS air mass per ground area
      wci     = 1.0_wp / max(wcapcan, tiny_num)
      f_sens  = coh_h_flux + coh_qw_flux + ground_h_flux + ground_qw_flux   ! [W/m2] into CAS from surfaces
      !----- atm<->CAS scalar conductance = rho*ustar*c3, where c3 (temp1) is the dimensionless   !
      !      Monin-Obukhov scalar-transfer coefficient from the aerodynamics solver (aero_out%      !
      !      temp1). Dropping it (c3=1) over-couples the CAS to the free atmosphere ~1/c3 (BUG8).   !
      !      The vapour twin's temp2 (== temp1 while z0q==z0h) rides in the caller-formed w_flux_ac !
      !      (= rho*ustar*temp2*(shv_atm - can_shv)). -----------------------------------------------!
      gatm    = rho_air * ustar * temp1                                    ! [kg/m2/s] atm<->CAS exchange
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

end module meds_column_energy
