! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_biophysics -- the stateless per-cohort vegetation-surface kernels, coupled    !
! through the interception film's wetted fraction:                                         !
!   * veg_energy_balance     -- a leaf or wood tissue's energy balance, relaxed exactly    !
!                               (docs/science/vegetation_energy_dynamics.md eqs. (1)-(7)), !
!                               and its four linearization coefficients;                   !
!   * intercept_canopy_layer -- one cohort's interception bucket                           !
!                               (docs/science/soil_biophysics.md, "Canopy interception").  !
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
   ! Leaf or wood energy balance (docs/science/vegetation_energy_dynamics.md eqs. (1)-(7)):!
   ! fluxes linearized about the canopy-air temperature, the tissue relaxed exactly over   !
   ! the step from t_store0 with heat capacity store_hcap_per_dt = cap/dt (0: zero         !
   ! inertia). Wood is the le_slope = le_ref = 0 case. A tissue with no area takes         !
   ! t_parent: its wood's temperature for a leaf, the canopy air's for wood. The optional  !
   ! f_wet splits the latent pathway into a dry share (stomata) and a wet share that       !
   ! evaporates the film (film_evap); absent, f_wet = 0.                                   !
   !---------------------------------------------------------------------------------------!
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
      real(wp), intent(in),  optional :: f_wet         !< [-] wetted canopy fraction (absent = 0, dry only)
      real(wp), intent(in),  optional :: le_slope_wet  !< [W/m2/K] film-evap latent slope (boundary layer only, no stomata)
      real(wp), intent(in),  optional :: le_ref_wet    !< [W/m2] film-evap latent reference term
      real(wp), intent(out), optional :: film_evap     !< [kg/m2/s] film evaporation (dew if negative)
      !----- The sapflow's advected enthalpy, Q_adv of eq. (1): it shifts the balance like abs_sw but     !
      !      stays out of drnet, being an internal soil-to-plant transfer (eq. (6)). Absent = 0. -------!
      real(wp), intent(in),  optional :: q_extra
      !----- h_v of eq. (2), the energy per kg evaporated that le_slope/le_ref were built with, to turn  !
      !      the latent flux back into mass. Default latent_heat_vap. ----------------------------------!
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

      !----- Coupling conductance (denom, without the tissue's own inertia) and forcing (numer) of eq. (1). !
      !      Both scale with the tissue's area, so splitting a cohort changes no temperature. ------------!
      denom   = h_coeff + les_dry + les_wet + lw_slope
      numer   = abs_sw + abs_lw + qx - (ler_dry + ler_wet)
      dt_prev = t_store0 - t_cas

      !=========================================================================================!
      ! Exact exponential relaxation, eqs. (3)-(7), with x = dt/tau = denom/store_hcap_per_dt:  !
      ! the endpoint (the committed state) takes the share w_avg = (1-exp(-x))/x of the         !
      ! start-of-step imbalance and the step average (every reported flux) w_flux = (1-w_avg)/x;!
      ! pairing them closes eq. (7). On the imbalance, it holds at x = 0 with no division by    !
      ! denom.                                                                                  !
      !=========================================================================================!
      if (store_hcap_per_dt > tiny_num) then
         x = denom / store_hcap_per_dt
         !----- Both weights without cancellation: their series at small x (Fortran has no expm1; a big  !
         !      cohort lives there), the 1/x tails at large x. ------------------------------------------!
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
         dt_temp   = dt_prev + w_avg  * imbalance / store_hcap_per_dt   ! endpoint: the committed state
         dt_avg    = dt_prev + w_flux * imbalance / store_hcap_per_dt   ! average: every reported flux
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
   ! The four `veg_energy_balance` linearization coefficients (leaf or wood), one          !
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

   !----- Film-evaporation conductance: the wet share's latent pathway, through the boundary layer   !
   !      only; same sidedness convention as leaf_transp_coeff, so leaf and wood share it. ---------!
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
   ! One cohort's canopy interception (docs/science/soil_biophysics.md, "Canopy            !
   ! interception"): a capacity-limited film bucket [kg/m2 ground] with a Beer interception!
   ! fraction, returning the wetted fraction sigma_w. The caller sweeps height-sorted      !
   ! cohorts top to bottom, each throughfall the next rain_above, so it is not elemental.  !
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
      !----- A film above capacity (leaf area shed under a full film) drips its excess to the     !
      !      layer below, so the column ledger keeps the water. -------------------------------------!
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
