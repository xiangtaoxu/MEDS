!==========================================================================================!
! meds_vegetation_biophysics -- the stateless per-cohort VEGETATION-surface kernels: the        !
! leaf/wood tissue energy balance and the canopy interception film they share. Coupled through   !
! the wetted-fraction film (`leaf_water` -> sigma_w -> leaf latent flux), so they live together.   !
!                                                                                          !
!   * veg_energy_diagnostic    -- the quasi-steady leaf OR wood surface temperature + its CAS/budget  !
!                                 flux contributions (the operator-split / ARK-diagnostic closure).    !
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

   public :: veg_energy_diagnostic, veg_energy_step_implicit, intercept_canopy_layer
   public :: sensible_heat_coeff, lw_emission_slope, le_conductance_flux, leaf_transp_coeff
   public :: leaf_film_coeff

   real(wp), parameter :: LEAF_MAXWHC = 0.11_wp     !< [kg/m2 leaf] film-holding capacity (wetted fraction)

contains

   !---------------------------------------------------------------------------------------!
   ! Diagnostic (quasi-steady) leaf OR wood surface energy balance -- the closure BOTH the      !
   ! operator-split and the ARK surface path share. Solves the linearized steady balance for the  !
   ! temperature offset from the CAS and returns the flux contributions the caller accumulates.    !
   ! ONE kernel serves leaf and wood: WOOD is the le_slope = le_ref = 0 case (no transpiration).    !
   !                                                                                          !
   !   t_emit   -- LW emission base: t_cas for the split sweep (reduces to the current form exactly)  !
   !               or the start-of-sub-step leaf temperature for the Picard/prognostic path.          !
   !   a_store  -- backward-Euler storage conductance cap/dt for a PROGNOSTIC leaf (0 = diagnostic);   !
   !               t_store0 is the store's start-of-sub-step temperature it relaxes from.              !
   ! With t_emit = t_cas and a_store = 0 (the ARK-diagnostic call) the two extra terms are exactly     !
   ! -0.0 / +0.0, so the result is bit-identical to the bare (abs_sw+abs_lw-le_ref) diagnostic.        !
   ! drnet groups (t_cas - t_emit) + dt_temp so the split's te = t_cas case cancels to 0.0 with NO      !
   ! rounding ulp (single authority for that ordering trick).                                           !
   !                                                                                          !
   ! WETTED-CANOPY EXTENSION (MEDS_ED2_RK45_DESIGN.md sec 3.4, P1, all four OPTIONAL): f_wet splits    !
   ! the single latent pathway into a (1-f_wet) DRY (stomatal, conductance le_slope/le_ref, unchanged)  !
   ! share that transpires and an f_wet WET (boundary-layer film, conductance le_slope_wet/le_ref_wet)   !
   ! share that evaporates/condenses as film_evap -- i.e. a wet leaf transpires less and evaporates its  !
   ! film instead, exactly the "compete for the leaf" partition sec 3.4 specifies. Any caller that omits  !
   ! f_wet (every existing call site: ARK's surface_derivs, and the split's own WOOD branch) gets f_wet=0,  !
   ! which makes the wet share vanish identically and reduces every line below to EXACTLY the pre-P1      !
   ! formula (dry share = (1-0)*le_slope = le_slope) -- so this is a BEHAVIOR-PRESERVING extension for     !
   ! every caller that doesn't opt in. -----------------------------------------------------------------!
   elemental pure subroutine veg_energy_diagnostic(abs_sw, abs_lw, h_coeff, le_slope, lw_slope, le_ref, &
                                         t_cas, t_emit, a_store, t_store0,                          &
                                         dt_temp, t_store, transp, dh, drnet,                        &
                                         f_wet, le_slope_wet, le_ref_wet, film_evap, q_extra)
      real(wp), intent(in)  :: abs_sw, abs_lw, h_coeff, le_slope, lw_slope, le_ref
      real(wp), intent(in)  :: t_cas, t_emit, a_store, t_store0
      real(wp), intent(out) :: dt_temp    !< temperature offset from the CAS [K]
      real(wp), intent(out) :: t_store    !< diagnosed store temperature (t_cas + dt_temp) [K]
      real(wp), intent(out) :: transp     !< transpiration mass flux (0 for wood) [kg/m2/s]
      real(wp), intent(out) :: dh         !< sensible-flux contribution to coh_h [W/m2]
      real(wp), intent(out) :: drnet      !< net-radiation contribution to coh_rnet [W/m2]
      real(wp), intent(in),  optional :: f_wet         !< [-] wetted canopy fraction (absent/0 = today's dry-only path)
      real(wp), intent(in),  optional :: le_slope_wet  !< [W/m2/K] film-evap latent slope (boundary layer only, no stomata)
      real(wp), intent(in),  optional :: le_ref_wet    !< [W/m2] film-evap latent reference term
      real(wp), intent(out), optional :: film_evap     !< [kg/m2/s] film evaporation (dew if negative)
      !----- NON-radiative advected-enthalpy source (P2 qloss/qwflux_wl): enters the dt_temp balance   !
      !      exactly like abs_sw (it shifts the equilibrium temperature and hence dh/transp), but is    !
      !      EXCLUDED from drnet on purpose -- it is an internal soil<->leaf transfer already debited    !
      !      from the soil's own root_heat_sink, not a whole-column boundary radiative input. Folding    !
      !      it into drnet (as an earlier version of this call did, by adding it to abs_sw itself)        !
      !      double-counts it in the coh_rnet-derived e_in ledger: dh/transp already carry qx into the    !
      !      CAS via the shifted dt_temp, so drnet must stay q_extra-free for the whole-column identity    !
      !      (state change == b-weighted boundary flux) to close. Absent/0 for every caller but the P2    !
      !      advective-enthalpy pre-pass (build_column_frozen), so this is a no-op elsewhere. -------------!
      real(wp), intent(in),  optional :: q_extra
      real(wp) :: fw, les_dry, ler_dry, les_wet, ler_wet, qx, denom, denom_true, g_slave
      !----- Coupling ("heat-capacity") FLOOR: h_coeff/les_*/lw_slope all scale with the tissue's own    !
      !      area index (LAI/WAI), so a just-recruited cohort (or a cohort whose wood area is smaller     !
      !      still) can present a denominator that is small in absolute W/m2/K terms without ever          !
      !      tripping the bare divide-by-zero guard (tiny_num). If the numerator carries ANY term that     !
      !      does not collapse in step (an upstream q_extra/abs_lw glitch, or simply roundoff at these      !
      !      magnitudes), the ratio is unbounded even though no single input is individually pathological.  !
      !      The floor bounds that ratio. It is applied to the DENOMINATOR (see below) so the kernel stays   !
      !      continuous: for any established canopy (les_dry alone typically reaches O(1-10) once stomata    !
      !      are open) the floor is inactive and the result is the exact unmodified balance. ----------------!
      real(wp), parameter :: veg_coupling_floor = 1.0_wp   !< [W/m2/K]
      fw = 0.0_wp ; les_wet = 0.0_wp ; ler_wet = 0.0_wp ; qx = 0.0_wp
      if (present(f_wet)) fw = f_wet
      if (present(le_slope_wet)) les_wet = fw * le_slope_wet
      if (present(le_ref_wet))   ler_wet = fw * le_ref_wet
      if (present(q_extra))      qx = q_extra
      les_dry = (1.0_wp - fw) * le_slope
      ler_dry = (1.0_wp - fw) * le_ref
      !----- The floor is applied to the DENOMINATOR, never as a branch on the result (P6). An earlier   !
      !      version special-cased the sub-floor case (dt_temp/transp/dh forced to 0), which put a JUMP   !
      !      DISCONTINUITY in this kernel -- and therefore in the whole-column ODE right-hand side -- at   !
      !      denom = veg_coupling_floor. That is fatal for an embedded-error adaptive integrator: les_dry  !
      !      carries the live tcas-dependent dqsat/qsat, so denom moves BETWEEN the stages of a single      !
      !      RK45 step, and a cohort whose coupling sits near the threshold has stages landing on opposite  !
      !      sides of the jump. The 5th- and 4th-order solutions then differ by an O(1) amount that does    !
      !      NOT shrink with dt, so the controller can never reach err <= 1: it rejects and halves down to   !
      !      the sub-step floor forever (observed: a growing canopy crossing the threshold turned a 6-min    !
      !      30-yr run into one that could not finish a single month). max() keeps every output CONTINUOUS   !
      !      in the state while still bounding the ratio: for denom >= floor this is EXACTLY the original    !
      !      formula (no bias for any real canopy), and below it dt_temp saturates at numer/floor instead     !
      !      of blowing up. Each flux below still carries its own conductance, which vanishes with LAI/WAI,   !
      !      so a negligible cohort contributes negligible fluxes without needing a special case. ------------!
      !----- SLAVING CONDUCTANCE (g_slave): the floor's deficit, made REAL and routed to the CAS.       !
      !      Clamping the denominator alone bounds dt_temp but leaves every FLUX on its own true         !
      !      conductance, so the tissue's balance no longer closes: it absorbs                            !
      !      N*(1 - denom_true/veg_coupling_floor) that leaves as nothing.  A diagnostic tissue holds      !
      !      no energy, so the caller books coh_rnet in and coh_h/coh_qw out with NOTHING re-checking      !
      !      the tissue between -- the gap is energy created from nowhere, and it showed up as the only     !
      !      unclosed term in the whole-column ledger of every vegetated run (~1 W/m2 on a 7-cohort          !
      !      stand; wood dominates, since WAI ~ 0.2*LAI puts it under the floor first).                     !
      !                                                                                          !
      !      The deficit is instead given to the tissue as an explicit extra conductance TO THE CANOPY      !
      !      AIR, and dh carries it.  That is exactly what "unresolvable" means physically, and it is       !
      !      ED2's own treatment (leaf_resolvable/wood_resolvable, stable_cohorts.f90) made continuous:      !
      !      a tissue too weakly coupled to resolve is thermally SLAVED to the canopy air, so whatever it    !
      !      absorbs passes straight through to the CAS instead of raising its temperature.  ED2 zeroes      !
      !      both the fluxes and the offset (trivially conservative); this keeps both and stays exact.       !
      !                                                                                          !
      !      Properties, all required and all held: EXACTLY conserving by construction (every outgoing       !
      !      term now sums to dt_temp*denom == the numerator, whichever branch max() took); CONTINUOUS in    !
      !      the state (g_slave -> 0 as denom_true -> the floor from below, so the P6 jump-discontinuity     !
      !      that broke the RK45 controller does not return); and IDENTICALLY ZERO for any established       !
      !      canopy (denom_true >= floor => g_slave = 0 => bit-identical to the pre-fix kernel).             !
      denom_true = h_coeff + les_dry + les_wet + lw_slope + a_store
      denom      = max(denom_true, veg_coupling_floor)
      g_slave    = denom - denom_true
      dt_temp = (abs_sw + abs_lw + qx - (ler_dry + ler_wet) - lw_slope * (t_cas - t_emit)             &
                 + a_store * (t_store0 - t_cas))                                                    &
                / denom
      t_store = t_cas + dt_temp
      transp  = (ler_dry + les_dry * dt_temp) / latent_heat_vap
      dh      = (h_coeff + g_slave) * dt_temp
      drnet   = abs_sw + abs_lw - lw_slope * ((t_cas - t_emit) + dt_temp)
      if (present(film_evap)) film_evap = (ler_wet + les_wet * dt_temp) / latent_heat_vap
   end subroutine veg_energy_diagnostic

   !---------------------------------------------------------------------------------------!
   ! The four `veg_energy_diagnostic` linearization coefficients (leaf OR wood, one authority     !
   ! each -- were duplicated inline between the split and ARK-frozen pre-pass call sites).        !
   !---------------------------------------------------------------------------------------!
   pure function sensible_heat_coeff(area_basis, gbh, rho, cp) result(h_coeff)
      real(wp), intent(in) :: area_basis   !< effarea_heat*LAI (leaf) or pi*WAI (wood) [m2/m2]
      real(wp), intent(in) :: gbh, rho, cp
      real(wp) :: h_coeff
      h_coeff = area_basis * gbh * rho * cp
   end function sensible_heat_coeff

   pure function leaf_transp_coeff(effarea_transp, lai, gbw, gsw) result(g_tr)
      real(wp), intent(in) :: effarea_transp, lai, gbw, gsw
      real(wp) :: g_tr
      g_tr = 0.0_wp
      if (gbw + gsw > tiny_num) g_tr = effarea_transp * lai * gbw * gsw / (gbw + gsw)
   end function leaf_transp_coeff

   !----- Film-evaporation conductance (sec 3.4, P1): the WET-fraction latent pathway, boundary   !
   !      layer ONLY (no stomatal resistance -- water sits directly on the surface), vs             !
   !      leaf_transp_coeff's dry-fraction boundary+stomata SERIES conductance just above. Same       !
   !      area_index/effarea sidedness convention as leaf_transp_coeff, so leaf and wood share one     !
   !      call (wood has no stomata anyway; its transp branch is gated off by the caller). -----------!
   pure function leaf_film_coeff(effarea_evap, area_index, gbw) result(g_ev)
      real(wp), intent(in) :: effarea_evap, area_index, gbw
      real(wp) :: g_ev
      g_ev = effarea_evap * area_index * gbw
   end function leaf_film_coeff

   pure function lw_emission_slope(emiss, te, area_index) result(lw_slope)
      real(wp), intent(in) :: emiss, te, area_index
      real(wp) :: lw_slope
      lw_slope = 4.0_wp * emiss * stefan * te**3 * area_index
   end function lw_emission_slope

   pure function le_conductance_flux(rho, g_tr, dq) result(le)
      real(wp), intent(in) :: rho, g_tr, dq
      real(wp) :: le
      le = latent_heat_vap * rho * g_tr * dq
   end function le_conductance_flux

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
