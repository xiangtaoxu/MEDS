! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_biophysics -- the stateless per-cohort VEGETATION-surface kernels: the        !
! leaf/wood tissue energy balance and the canopy interception film they share. Coupled     !
! through the wetted-fraction film (`leaf_water` -> sigma_w -> leaf latent flux), so they  !
! live together.                                                                           !
!                                                                                          !
!   * veg_energy_balance     -- a leaf OR wood tissue's energy balance, relaxed exactly    !
!                               over the step (docs/science/vegetation_energy_dynamics.md  !
!                               eqs. (1)-(7)).                                             !
!   * intercept_canopy_layer -- per-cohort capacity-limited interception bucket (Beer      !
!                               fraction); the caller sweeps height-sorted cohorts         !
!                               top->bottom (docs/science/soil_biophysics.md, "Canopy      !
!                               interception").                                            !
!==========================================================================================!
module meds_plant_biophysics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : pi, tiny_num, stefan, cp_air, cp_liq, latent_heat_vap
   use meds_therm_lib,        only : internal_energy_to_temp, sat_specific_humidity,                        &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor
   use meds_plant_types, only : veg_thermal_params_t
   implicit none
   private

   public :: veg_energy_balance, intercept_canopy_layer
   public :: sensible_heat_coeff, lw_emission_slope, leaf_transp_coeff
   public :: leaf_film_coeff


contains

   !---------------------------------------------------------------------------------------!
   ! Leaf OR wood surface energy balance (docs/science/vegetation_energy_dynamics.md eqs.  !
   ! (1)-(7)): every flux linearized about the canopy-air temperature, and the tissue      !
   ! relaxed EXACTLY over the step toward its zero-inertia offset with its own heat        !
   ! capacity. Returns the endpoint temperature and the step-average flux contributions    !
   ! the caller accumulates. ONE kernel serves leaf and wood: WOOD is the le_slope =       !
   ! le_ref = 0 case (no transpiration).                                                   !
   !                                                                                       !
   !   t_emit            -- LW emission base, the temperature lw_slope linearizes emission !
   !                        about. Every caller passes t_cas, so the (t_cas - t_emit) terms !
   !                        are exactly 0.                                                 !
   !   store_hcap_per_dt -- the tissue's heat capacity over the step, cap/dt (eq. (8) sizes !
   !                        cap); 0 gives the zero-inertia limit. t_store0 is the          !
   !                        start-of-step tissue temperature the store relaxes from.       !
   ! drnet groups (t_cas - t_emit) + dt_avg so the t_emit = t_cas case cancels to 0.0 with !
   ! NO rounding ulp (single authority for that ordering trick).                           !
   !                                                                                       !
   ! WETTED CANOPY (all four OPTIONAL): f_wet splits the latent pathway into a (1-f_wet)   !
   ! DRY share that transpires through boundary layer and stomata in series                !
   ! (le_slope/le_ref) and an f_wet WET share that evaporates or condenses the intercepted !
   ! film through the boundary layer alone (le_slope_wet/le_ref_wet), returned as          !
   ! film_evap: a wet leaf transpires less and evaporates its film instead. A caller that  !
   ! omits f_wet gets f_wet = 0, which makes the wet share vanish identically. ------------!
   elemental pure subroutine veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope, lw_slope, le_ref, &
                                         t_cas, t_emit, store_hcap_per_dt, t_store0,                          &
                                         dt_temp, t_store, transp, dh, drnet,                        &
                                         f_wet, le_slope_wet, le_ref_wet, film_evap, q_extra,      &
                                         h_evap, h_evap_wet)
      real(wp), intent(in)  :: abs_sw, abs_lw, h_coeff, le_slope, lw_slope, le_ref
      real(wp), intent(in)  :: t_cas, t_emit, store_hcap_per_dt, t_store0
      real(wp), intent(out) :: dt_temp    !< temperature offset from the CAS [K]
      real(wp), intent(out) :: t_store    !< diagnosed store temperature (t_cas + dt_temp) [K]
      real(wp), intent(out) :: transp     !< transpiration mass flux (0 for wood) [kg/m2/s]
      real(wp), intent(out) :: dh         !< sensible-flux contribution to coh_h [W/m2]
      real(wp), intent(out) :: drnet      !< net-radiation contribution to coh_rnet [W/m2]
      real(wp), intent(in),  optional :: f_wet         !< [-] wetted canopy fraction (absent/0 = today's dry-only path)
      real(wp), intent(in),  optional :: le_slope_wet  !< [W/m2/K] film-evap latent slope (boundary layer only, no stomata)
      real(wp), intent(in),  optional :: le_ref_wet    !< [W/m2] film-evap latent reference term
      real(wp), intent(out), optional :: film_evap     !< [kg/m2/s] film evaporation (dew if negative)
      !----- NON-radiative advected-enthalpy source (Q_adv of eq. (1): the enthalpy the sapflow        !
      !      carries in, built by build_column_frozen): enters the dt_temp balance exactly like abs_sw !
      !      (it shifts the equilibrium temperature and hence dh/transp), but is EXCLUDED from drnet   !
      !      on purpose -- it is an internal soil<->leaf transfer already debited from the soil's own  !
      !      root_heat_sink, not a whole-column boundary radiative input. Folding it into drnet would  !
      !      double-count it in the coh_rnet-derived e_in ledger: dh/transp already carry qx into the  !
      !      CAS via the shifted dt_temp, so drnet must stay q_extra-free for the whole-column         !
      !      identity (state change == b-weighted boundary flux) to close. Absent = 0. ----------------!
      real(wp), intent(in),  optional :: q_extra
      !----- Energy the store pays per kg of water it evaporates -- the factor le_slope/le_ref were  !
      !      built with, needed here only to turn the energy flux back into a MASS flux (h_v of eq.  !
      !      (2)). Default latent_heat_vap; the fast loop passes the full vapour enthalpy at the     !
      !      canopy-air temperature the balance is linearized about, so the store pays what the CAS  !
      !      receives. ------------------------------------------------------------------------------!
      real(wp), intent(in),  optional :: h_evap      !< [J/kg] for transpiration
      real(wp), intent(in),  optional :: h_evap_wet  !< [J/kg] for film evaporation
      real(wp) :: fw, les_dry, ler_dry, les_wet, ler_wet, qx, denom, denom_true, g_slave, he, hw
      real(wp) :: numer, dt_diag, dt_prev, dt_avg, x, w_end, w_avg
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
      he = latent_heat_vap ; if (present(h_evap))     he = h_evap
      hw = latent_heat_vap ; if (present(h_evap_wet)) hw = h_evap_wet
      if (present(f_wet)) fw = f_wet
      if (present(le_slope_wet)) les_wet = fw * le_slope_wet
      if (present(le_ref_wet))   ler_wet = fw * le_ref_wet
      if (present(q_extra))      qx = q_extra
      les_dry = (1.0_wp - fw) * le_slope
      ler_dry = (1.0_wp - fw) * le_ref
      !----- The floor is applied to the DENOMINATOR, never as a branch on the result: a branch        !
      !      (dt_temp/transp/dh forced to 0 below the floor) would put a JUMP DISCONTINUITY in this    !
      !      kernel -- and so in the whole-column ODE right-hand side -- at denom =                    !
      !      veg_coupling_floor. That is fatal for an embedded-error adaptive integrator: les_dry      !
      !      carries the live tcas-dependent dqsat/qsat, so denom moves BETWEEN the stages of a single !
      !      RK45 step, and a cohort whose coupling sits near the threshold has stages landing on      !
      !      opposite sides of the jump. The 5th- and 4th-order solutions then differ by an O(1)       !
      !      amount that does NOT shrink with dt, so the controller can never reach err <= 1 and       !
      !      halves down to the sub-step floor forever: a growing canopy crossing the threshold stalls !
      !      the run. max() keeps every output CONTINUOUS in the state while still bounding the ratio: !
      !      for denom >= floor this is EXACTLY the unfloored formula (no bias for any real canopy),   !
      !      and below it dt_temp saturates at numer/floor instead of blowing up. Each flux below      !
      !      still carries its own conductance, which vanishes with LAI/WAI, so a negligible cohort    !
      !      contributes negligible fluxes without needing a special case. ----------------------------!
      !----- SLAVING CONDUCTANCE (g_slave, docs/science/vegetation_energy_dynamics.md eq. (1)): the    !
      !      floor's deficit, made REAL and routed to the CAS. Clamping the denominator alone bounds   !
      !      dt_temp but leaves every FLUX on its own true conductance, so the tissue's balance would  !
      !      not close: it would absorb N*(1 - denom_true/veg_coupling_floor) that leaves as nothing   !
      !      -- energy created from nowhere in the whole-column ledger (~1 W/m2 on a 7-cohort stand;   !
      !      wood dominates, since WAI ~ 0.2*LAI puts it under the floor first).                       !
      !                                                                                                !
      !      The deficit is instead given to the tissue as an explicit extra conductance TO THE CANOPY !
      !      AIR, and dh carries it. That is exactly what "unresolvable" means physically, and it is   !
      !      ED2's own treatment (leaf_resolvable/wood_resolvable, stable_cohorts.f90) made            !
      !      continuous: a tissue too weakly coupled to resolve is thermally SLAVED to the canopy air, !
      !      so whatever it absorbs passes straight through to the CAS instead of raising its          !
      !      temperature. ED2 zeroes both the fluxes and the offset (trivially conservative); this     !
      !      keeps both and stays exact.                                                               !
      !                                                                                                !
      !      Properties, all required and all held: EXACTLY conserving by construction (the outgoing   !
      !      terms and the storage change sum to the numerator, eq. (7), whichever branch max() took); !
      !      CONTINUOUS in the state (g_slave -> 0 as denom_true -> the floor from below, so there is  !
      !      no jump); and IDENTICALLY ZERO for any established canopy (denom_true >= floor => g_slave !
      !      = 0).                                                                                     !
      !----- COUPLING conductance D0 of eq. (1): everything that ties the tissue to its surroundings.  !
      !      It EXCLUDES store_hcap_per_dt, which is a property of the tissue's own inertia, not of    !
      !      the coupling.                                                                             !
      denom_true = h_coeff + les_dry + les_wet + lw_slope
      denom      = max(denom_true, veg_coupling_floor)
      g_slave    = denom - denom_true

      !----- The DIAGNOSTIC (zero-inertia) offset: where the tissue would sit if it had no memory. ---!
      numer   = abs_sw + abs_lw + qx - (ler_dry + ler_wet) - lw_slope * (t_cas - t_emit)
      dt_diag = numer / denom
      dt_prev = t_store0 - t_cas

      !=========================================================================================!
      ! EXACT EXPONENTIAL RELAXATION (docs/science/vegetation_energy_dynamics.md eqs. (3)-(7),  !
      ! and "Why not backward Euler"). Under the Category-0 coefficient freeze the tissue       !
      ! equation cap*dT/dt = numer - denom*(T - t_cas) is LINEAR with tau = cap/denom, so the   !
      ! step has a closed form, and using it removes the time-stepping error entirely rather    !
      ! than bounding it:                                                                       !
      !                                                                                         !
      !   x     = dt/tau = denom/store_hcap_per_dt   (store_hcap_per_dt = cap/dt)               !
      !   w_end = exp(-x)                            weight on the OLD state at the ENDPOINT    !
      !   w_avg = (1 - exp(-x))/x                    weight on the OLD state, STEP-AVERAGED     !
      !                                                                                         !
      ! The two weights are NOT interchangeable: the committed STATE is the endpoint, while the !
      ! flux the canopy air receives is the step AVERAGE, and pairing them closes               !
      ! store_hcap_per_dt*(dt_end - dt_prev) + denom*dt_avg == numer identically.               !
      ! store_hcap_per_dt -> 0 gives x -> infinity, both weights -> 0 and the pure steady       !
      ! balance: "diagnostic" is this kernel's zero-inertia LIMIT, not a separate branch.       !
      !=========================================================================================!
      if (store_hcap_per_dt > tiny_num) then
         x = denom / store_hcap_per_dt
         !----- exp(-x) and (1-exp(-x))/x, both accurate across the whole range. Fortran has no      !
         !      expm1, so the small-x branch uses the series (relative error ~ x^4/120 < 1e-18 at    !
         !      the threshold) -- without it, 1 - exp(-x) cancels catastrophically as x -> 0, which  !
         !      is exactly the large-capacity limit a big cohort lives in. The large-x branch avoids  !
         !      underflow and gives the correct 1/x tail for w_avg. -------------------------------!
         if (x < 1.0e-4_wp) then
            w_end = 1.0_wp - x*(1.0_wp - 0.5_wp*x*(1.0_wp - x/3.0_wp))
            w_avg = 1.0_wp - 0.5_wp*x*(1.0_wp - x/3.0_wp*(1.0_wp - 0.25_wp*x))
         else if (x > 500.0_wp) then
            w_end = 0.0_wp
            w_avg = 1.0_wp / x
         else
            w_end = exp(-x)
            w_avg = (1.0_wp - w_end) / x
         end if
      else
         w_end = 0.0_wp ; w_avg = 0.0_wp            ! zero heat capacity => pure diagnostic
      end if

      dt_temp = (1.0_wp - w_end) * dt_diag + w_end * dt_prev      ! ENDPOINT  -> the committed state
      dt_avg  = (1.0_wp - w_avg) * dt_diag + w_avg * dt_prev      ! AVERAGE   -> every reported flux

      t_store = t_cas + dt_temp
      transp  = (ler_dry + les_dry * dt_avg) / he
      dh      = (h_coeff + g_slave) * dt_avg
      drnet   = abs_sw + abs_lw - lw_slope * ((t_cas - t_emit) + dt_avg)
      if (present(film_evap)) film_evap = (ler_wet + les_wet * dt_avg) / hw
   end subroutine veg_energy_balance


   !---------------------------------------------------------------------------------------!
   ! The four `veg_energy_balance` linearization coefficients (leaf OR wood), one          !
   ! authority each (docs/science/vegetation_energy_dynamics.md eqs. (1)-(2)).             !
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

   !----- Film-evaporation conductance: the WET-fraction latent pathway, boundary layer ONLY (no   !
   !      stomatal resistance -- water sits directly on the surface), vs leaf_transp_coeff's       !
   !      dry-fraction boundary+stomata SERIES conductance just above. Same area_index/effarea     !
   !      sidedness convention as leaf_transp_coeff, so leaf and wood share one call (wood has no  !
   !      stomata anyway; its transp branch is gated off by the caller). --------------------------!
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


   !---------------------------------------------------------------------------------------!
   ! Per-cohort canopy interception (docs/science/soil_biophysics.md, "Canopy              !
   ! interception"). One canopy layer: capacity-limited bucket (leaf_max*LAI +             !
   ! wood_max*SAI, the PFT's film capacities per unit leaf and wood area) with a Beer      !
   ! interception fraction; the caller sweeps the height-sorted cohorts top->bottom        !
   ! feeding each `throughfall` as the next cohort's `rain_above`. `leaf_water` is the     !
   ! per-cohort prognostic film [kg/m2 ground]; `sigma_w` (wetted fraction) is exported    !
   ! for the CAS-space evaporation. Not `elemental`: the cascade is a sequential           !
   ! recurrence.                                                                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine intercept_canopy_layer(leaf_water, rain_above, lai, sai, e_canopy, dt,     &
                                          leaf_max, wood_max, k_int, alpha_pi, throughfall, drip, sigma_w)
      real(wp), intent(inout) :: leaf_water
      real(wp), intent(in)    :: rain_above, lai, sai, e_canopy, dt, k_int, alpha_pi
      real(wp), intent(in)    :: leaf_max, wood_max   !< [kg/m2 leaf], [kg/m2 wood] film capacities
      real(wp), intent(out)   :: throughfall, drip, sigma_w
      real(wp) :: pai, f_pi, w_max, q_grab, room, q_intr, w_new, over
      pai    = lai + sai
      f_pi   = alpha_pi * (1.0_wp - exp(-k_int * pai))
      w_max  = leaf_max * lai + wood_max * sai
      q_grab = f_pi * rain_above
      room   = max(0.0_wp, w_max - leaf_water) / dt
      q_intr = min(q_grab, room + e_canopy)                      ! bounded by capacity + evap headroom
      !----- A film already ABOVE capacity -- the canopy lost leaf area under a full film, e.g.    !
      !      the daily slow step shedding leaves -- DRIPS its excess to the layer below. Clipping  !
      !      it instead would send the water to neither the ground nor any flux, and the           !
      !      whole-column ledger would lose exactly that much in heavy rain. ----------------------!
      w_new = leaf_water + (q_intr - e_canopy) * dt
      over  = max(0.0_wp, w_new - w_max)
      leaf_water  = min(max(w_new, 0.0_wp), w_max)
      drip        = max(0.0_wp, q_grab - q_intr) + over / dt
      throughfall = (rain_above - q_grab) + drip                 ! gap-throughfall + drip -> below
      if (w_max > tiny_num) then
         sigma_w = min(1.0_wp, (leaf_water / w_max) ** (2.0_wp / 3.0_wp))
      else
         sigma_w = 0.0_wp
      end if
   end subroutine intercept_canopy_layer

end module meds_plant_biophysics
