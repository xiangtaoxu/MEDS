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
   use meds_constants,        only : tiny_num, stefan_boltzmann, latent_heat_vap
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
   !   store_hcap_per_dt -- the tissue's heat capacity over the step, cap/dt (eq. (8) sizes !
   !                        cap); 0 gives the zero-inertia limit. t_store0 is the          !
   !                        start-of-step tissue temperature the store relaxes from.       !
   !   t_parent          -- the temperature of what the tissue grows from: its wood for a  !
   !                        leaf, the canopy air for wood. A tissue with no area has no    !
   !                        coupling and no heat capacity, and takes this temperature      !
   !                        (docs/science/vegetation_energy_dynamics.md, "A tissue with no !
   !                        area"); nothing weights it, since every flux it carries is 0.  !
   !                                                                                       !
   ! WETTED CANOPY (all four OPTIONAL): f_wet splits the latent pathway into a (1-f_wet)   !
   ! DRY share that transpires through boundary layer and stomata in series                !
   ! (le_slope/le_ref) and an f_wet WET share that evaporates or condenses the intercepted !
   ! film through the boundary layer alone (le_slope_wet/le_ref_wet), returned as          !
   ! film_evap: a wet leaf transpires less and evaporates its film instead. A caller that  !
   ! omits f_wet gets f_wet = 0, which makes the wet share vanish identically. ------------!
   elemental pure subroutine veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope, lw_slope, le_ref, &
                                         t_cas, t_parent, store_hcap_per_dt, t_store0,              &
                                         dt_temp, t_store, transp, dh, drnet,                        &
                                         f_wet, le_slope_wet, le_ref_wet, film_evap, q_extra,      &
                                         h_evap, h_evap_wet)
      real(wp), intent(in)  :: abs_sw, abs_lw, h_coeff, le_slope, lw_slope, le_ref
      real(wp), intent(in)  :: t_cas, t_parent, store_hcap_per_dt, t_store0
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
      real(wp) :: fw, les_dry, ler_dry, les_wet, ler_wet, qx, denom, he, hw
      real(wp) :: numer, imbalance, dt_prev, dt_avg, x, w_avg, w_flux
      fw = 0.0_wp
      les_wet = 0.0_wp
      ler_wet = 0.0_wp
      qx = 0.0_wp
      he = latent_heat_vap
      if (present(h_evap))     he = h_evap
      hw = latent_heat_vap
      if (present(h_evap_wet)) hw = h_evap_wet
      if (present(f_wet)) fw = f_wet
      if (present(le_slope_wet)) les_wet = fw * le_slope_wet
      if (present(le_ref_wet))   ler_wet = fw * le_ref_wet
      if (present(q_extra))      qx = q_extra
      les_dry = (1.0_wp - fw) * le_slope
      ler_dry = (1.0_wp - fw) * le_ref

      !----- COUPLING conductance D0 and FORCING of eq. (1). D0 is everything that ties the tissue to   !
      !      its surroundings; it EXCLUDES store_hcap_per_dt, which is the tissue's own inertia. Every   !
      !      term of both scales with the tissue's own area, as its heat capacity does (eq. (8)), so a   !
      !      thin cohort is simply a small one: splitting a cohort changes no temperature. -------------!
      denom   = h_coeff + les_dry + les_wet + lw_slope
      numer   = abs_sw + abs_lw + qx - (ler_dry + ler_wet)
      dt_prev = t_store0 - t_cas

      !=========================================================================================!
      ! EXACT EXPONENTIAL RELAXATION (docs/science/vegetation_energy_dynamics.md eqs. (3)-(7),  !
      ! and "Why not backward Euler"). Under the Category-0 coefficient freeze the tissue       !
      ! equation cap*dT/dt = numer - denom*(T - t_cas) is LINEAR with tau = cap/denom, so the   !
      ! step has a closed form, and using it removes the time-stepping error entirely rather    !
      ! than bounding it:                                                                       !
      !                                                                                         !
      !   x      = dt/tau = denom/store_hcap_per_dt   (store_hcap_per_dt = cap/dt)              !
      !   w_avg  = (1 - exp(-x))/x                   the ENDPOINT's share of the start-of-step  !
      !                                              imbalance, per unit of capacity           !
      !   w_flux = (1 - w_avg)/x                     the same share for the step AVERAGE        !
      !                                                                                         !
      ! The committed STATE is the endpoint, while the flux the canopy air receives is the step  !
      ! AVERAGE, and pairing them closes store_hcap_per_dt*(dt_end - dt_prev) + denom*dt_avg ==  !
      ! numer identically (eq. (7)). Written on the start-of-step imbalance numer -              !
      ! denom*dt_prev, the step holds for every x >= 0 with no division by denom: a store with   !
      ! no coupling (x = 0) warms by numer*dt/cap, and the heat store bounds the step whatever   !
      ! the coupling. As x grows both shares tend to the zero-inertia balance numer/denom. -----!
      !=========================================================================================!
      if (store_hcap_per_dt > tiny_num) then
         x = denom / store_hcap_per_dt
         !----- (1-exp(-x))/x and (x - 1 + exp(-x))/x**2, both accurate across the whole range.        !
         !      Fortran has no expm1, so the small-x branch uses their series (relative error ~ x**4/120  !
         !      < 1e-18 at the threshold) -- without it, 1 - exp(-x) cancels catastrophically as         !
         !      x -> 0, which is exactly the large-capacity limit a big cohort lives in. The large-x     !
         !      branch avoids underflow and gives the 1/x tails. ------------------------------------!
         if (x < 1.0e-4_wp) then
            w_avg  = 1.0_wp - 0.5_wp*x*(1.0_wp - x/3.0_wp*(1.0_wp - 0.25_wp*x))
            w_flux = 0.5_wp*(1.0_wp - x/3.0_wp*(1.0_wp - 0.25_wp*x*(1.0_wp - 0.2_wp*x)))
         else if (x > 500.0_wp) then
            w_avg  = 1.0_wp / x
            w_flux = (1.0_wp - w_avg) / x
         else
            w_avg  = (1.0_wp - exp(-x)) / x
            w_flux = (1.0_wp - w_avg) / x
         end if
         imbalance = numer - denom * dt_prev
         dt_temp   = dt_prev + w_avg  * imbalance / store_hcap_per_dt   ! ENDPOINT -> the committed state
         dt_avg    = dt_prev + w_flux * imbalance / store_hcap_per_dt   ! AVERAGE  -> every reported flux
      else if (denom > tiny_num) then
         dt_temp = numer / denom          ! zero heat capacity: the steady balance, at once
         dt_avg  = dt_temp
      else
         dt_temp = t_parent - t_cas       ! no area: no coupling and no heat (eq. (8))
         dt_avg  = dt_temp
      end if

      t_store = t_cas + dt_temp
      transp  = (ler_dry + les_dry * dt_avg) / he
      dh      = h_coeff * dt_avg
      drnet   = abs_sw + abs_lw - lw_slope * dt_avg
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
      lw_slope = 4.0_wp * emiss * stefan_boltzmann * te**3 * area_index
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
