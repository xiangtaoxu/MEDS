!==========================================================================================!
! meds_fast_time_derivs -- the pure whole-column tendency (RHS) for the fast-loop IMEX-ARK        !
! integrator (design docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md, phase P0). The overhaul replaces the       !
! operator-split + Picard fast step with ONE additive Runge-Kutta advance; an ARK needs a         !
! side-effect-free f(y) = dy/dt for the whole column, whereas today every kernel is               !
! *advance-and-commit* (forms a flux, hides it inside a backward-Euler denominator, and mutates   !
! its `intent(inout)` store). This module gives the ARK stepper that RHS.                          !
!                                                                                          !
! TWO ENTRY POINTS:                                                                               !
!   * surface_derivs -- the CAS surface block (leaf-energy diagnostic + ground skin + the three    !
!     CAS twins). Depends only on meds_therm_lib; validated bit-for-bit against the split.            !
!   * column_derivs  -- the WHOLE column: surface_derivs + the soil-heat column, the soil-water    !
!     (Richards) column, and the per-cohort plant-hydraulics 2x2, assembled from the tendency-      !
!     exposing kernels soil_energy_time_deriv / soil_water_time_deriv / plant_water_tendency (each a    !
!     side-effect-free sibling of its BE kernel, reusing the SAME flux helpers -- no re-derivation).!
!                                                                                          !
! Faithfulness: each reservoir's tendency is the EXPLICIT RHS whose backward-Euler advance over    !
! dt reproduces the corresponding split kernel as dt -> 0 (soil heat/water) or exactly (the CAS     !
! implicit-in-atm twins; the 2x2 matrix-exponential). Verified per-reservoir in test_column_derivs.!
! The surface hydrology BCs (q_top, root uptake) and the leaf gas-exchange pre-pass are            !
! the frozen, explicit part of the additive split (held constant across the ARK macro-step).       !
!==========================================================================================!
module meds_fast_time_derivs
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : latent_heat_vap, stefan, cp_air, tiny_num, rho_h2o, mmdry
   use meds_therm_lib,           only : cas_molar_density, cas_temp_of_enthalpy, sat_specific_humidity,                    &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_to_temp,       &
                                     internal_energy_liquid
   use meds_biophysics_types, only : energy_forcing_t
   use meds_column_constants, only : n_soil_layer_max
   use meds_column_reservoirs, only : soil_energy_column_t
   use meds_column_params, only : soil_thermal_params_t, soil_params_t
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t
   use meds_soil_energy,      only : soil_energy_time_deriv
   use meds_soil_water,       only : soil_water_time_deriv
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_time_deriv
   use meds_ground_biophysics, only : ground_surface_fluxes
   use meds_canopy_aerodynamics, only : mo_surface_layer, cas_atm_conductances
   use meds_plant_biophysics, only : veg_energy_balance, lw_emission_slope
   use meds_column_state_ops, only : assemble_soil_energy_forcing
   use meds_fast_snow,        only : snow_stage_t
   use meds_fast_types,       only : surface_state_t, surface_tend_t, cas_boundary_t,             &
                                     tissue_coefficients_t, canopy_film_capacity_t, ground_boundary_t, &
                                     column_state_t, column_frozen_t, column_tend_t,               &
                                     stage_bflux_t, column_bflux_t
   implicit none
   private

   public :: surface_derivs, column_derivs, cas_conductances
   !----- exported so the split path relaxes on the SAME timescale rather than keeping a copy    !
   !      that could drift out of step with this one. -------------------------------------------!
   public :: TAU_COND

   !----- CAS supersaturation relaxation timescale [s] for the smooth condensation sink: q relaxes to    !
   !      qsat with this e-folding time (dq/dt|_cond = -(q-qsat)/TAU_COND). Short vs dt_fast=1800s so a    !
   !      supersaturated canopy air condenses within the step, but SMOOTH (no clamp) so the ARK2 embedded  !
   !      error stays small. Physical dew/fog forms fast; the exact value is not critical (any strong     !
   !      restoring rate pins q near qsat). ----------------------------------------------------------!
   real(wp), parameter :: TAU_COND = 300.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! refresh_cas_conductances -- N2a.  Re-solve ONLY the bulk Monin-Obukhov surface layer at a live   !
   ! canopy-air state and return the three CAS<->atmosphere conductances.  Deliberately NOT a full     !
   ! canopy_aerodynamics call: the per-cohort boundary layers (leaf_gbw/leaf_gbh, and therefore        !
   ! h_coeff_leaf/g_transp_leaf) stay frozen, because the sec-1g decomposition puts <= 0.7% of the lag gain on    !
   ! them while they are what makes a full refresh expensive.  Canopy geometry (displacement,           !
   ! roughness) and the reference-level state are frozen inputs -- only (T_cas, q_cas) is live, which   !
   ! is exactly the loop that was oscillating.                                                          !
   !                                                                                          !
   ! temp2 == temp1 in canopy_aerodynamics (z0q = z0h), so vapour reuses the heat transfer factor;      !
   ! the CO2 conductance rides the same factor on the MOLAR capacity, computed at the live humidity.    !
   !---------------------------------------------------------------------------------------!
   pure subroutine refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp) :: tcas, ustar, temp1, zeta, rib, obu, can_dmol
      tcas = cas_temp_of_enthalpy(cas_enthalpy, cas_shv)
      !----- column_prepass sets aenv%can_theta = tcas, so passing tcas here reproduces exactly    !
      !      what the frozen pre-pass fed this same kernel -- the ONLY difference is the state. ---!
      call mo_surface_layer(cas%aero_cfg, cas%mo_u_ref, cas%mo_zref, cas%mo_displace, cas%mo_rough, &
                            cas%mo_theta_atm, cas%mo_shv_atm, tcas, cas_shv,                      &
                            ustar, temp1, zeta, rib, obu)
      can_dmol = cas_molar_density(cas%mo_rho, cas_shv)
      call cas_atm_conductances(cas%mo_rho, can_dmol, ustar, temp1, temp1, g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine refresh_cas_conductances

   !---------------------------------------------------------------------------------------!
   ! cas_conductances -- the ONE place that decides whether a surface-layer re-solve is needed.       !
   ! Returns the CAS<->atmosphere conductances for the canopy-air state (cas_enthalpy, cas_shv): a     !
   ! fresh solve when the bundle carries live MO inputs, otherwise the bundle's own values verbatim.   !
   ! Single-sourcing the rule is what keeps a tendency and the boundary flux charged against it from   !
   ! ever disagreeing -- see stage_bnd (meds_fast_rk45) for the failure this prevents.                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      if (cas%mo_live) then
         call refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      else
         g_atm_heat = cas%g_atm_heat ; g_atm_vapour = cas%g_atm_vapour ; g_atm_co2 = cas%g_atm_co2
      end if
   end subroutine cas_conductances



   !---------------------------------------------------------------------------------------!
   ! surface_derivs -- the CAS surface-block RHS (leaf-energy diagnostic + ground skin + the three  !
   ! CAS twins). Faithful, side-effect-free transcription of column_fast_step's surface path         !
   ! (build_column_frozen). The longwave emission base is the current tcas, so the   !
   ! split's `tcas - te` term is identically zero. Integrating d_cas_* with a backward-Euler-in-the-  !
   ! atmosphere step reproduces the split's committed enth1/shv1/co21 exactly.                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine surface_derivs(y, cas, tissue, film, ground, snow, t_ground, n, f)
      type(surface_state_t),        intent(in)  :: y
      type(cas_boundary_t),         intent(in)  :: cas      !< CAS <-> atmosphere boundary (a stage may pass its own conductances)
      type(tissue_coefficients_t),  intent(in)  :: tissue   !< per-cohort leaf/wood energy coefficients
      type(canopy_film_capacity_t), intent(in)  :: film     !< interception film
      type(ground_boundary_t),      intent(in)  :: ground   !< bare-ground boundary
      type(snow_stage_t),           intent(in)  :: snow     !< the snow stage's outcome
      real(wp),               intent(in)  :: t_ground   !< [K] soil-top temperature at THIS evaluation (a live input, not frozen)
      integer(ik),            intent(in)  :: n
      type(surface_tend_t),   intent(out) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat
      real(wp)    :: lw_slope, le_slope, le_ref, dtl, tl, transp_i, dh, drnet
      real(wp)    :: lw_slope_w, dtw, tw, transp_w                    !< diagnostic WOOD balance (own store)
      real(wp)    :: coh_h, coh_qw, coh_transp, coh_rnet, coh_film_evap
      real(wp)    :: h_evap_l, h_film_l, h_evap_w, h_film_w   !< [J/kg] energy per kg evaporated (leaf/wood; transp/film)
      real(wp)    :: h_bare, le_soil   !< bare-soil half of the snowfac blend (C4)
      !----- Canopy-SURFACE water (sec 3.4, P2c): the wetted-fraction film-evap latent terms, using the  !
      !      FROZEN conductance (film%g_film_leaf/w, sec 3.4/P1's leaf_film_coeff, precomputed once in the    !
      !      Act-1 pre-pass) but state-dependent dqdt/qsat_c-qcas -- mirrors le_slope/le_ref's own          !
      !      frozen-conductance/live-state split for the dry pathway just above. Harmless when              !
      !      canopy_water_on is off: film%f_wet_c(i) stays 0.0 there, which makes veg_energy_balance's     !
      !      wet pathway vanish identically (proven in its own docstring), so these are pure no-ops on        !
      !      the byte-identical default path. --------------------------------------------------------------!
      real(wp)    :: le_slope_wet, le_ref_wet, le_slope_wet_w, le_ref_wet_w
      real(wp)    :: gah_l, gaw_l, gac_l      !< N2a: the conductances this evaluation uses
      type(cas_source_t) :: cas_src
      type(cas_column_t) :: cas_col
      integer(ik) :: i

      allocate(f%leaf_temp(n), f%wood_temp(n), f%transp_c(n), f%film_evap_leaf(n), f%film_evap_wood(n))

      tcas   = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, cas%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, cas%press)

      coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
      coh_film_evap = 0.0_wp
      do i = 1_ik, n
         lw_slope = lw_emission_slope(tissue%leaf_emiss, tcas, tissue%lai(i))
         !----- The leaf pays the FULL specific enthalpy of the vapour it sheds, h_evap = enthalpy_vapor  !
         !      at the canopy-air temperature the balance is linearized around (the same reference        !
         !      qsat_c/dqdt use; the cp_vap*(t_leaf - t_cas) difference is ~0.2% of h_evap), and the      !
         !      CAS receives the SAME number -- one flux, both sides. It used to pay latent_heat_vap      !
         !      only, with the liquid part (h_evap - L) charged to the SOIL through coh_qsoil. Once P2    !
         !      added the explicit soil -> wood -> leaf advected-enthalpy chain (qloss, qwflux_wl) that   !
         !      proxy charged the soil TWICE and heated the leaf for free (~30 W/m2 at 3 mm/day). ED2's   !
         !      rk4_derivs charges the leaf tq2enthalpy(T_leaf) and the soil uint_water once; so does    !
         !      this now (2026-09 review, item 1A #10). The FILM pays h_evap minus the liquid enthalpy    !
         !      the intercepted water arrived with (film_liquid_enthalpy), so film store + leaf + CAS close with no !
         !      slack term (item 1A #2). ------------------------------------------------------------------!
         h_evap_l = enthalpy_vapor(tcas)
         h_film_l = h_evap_l - film%film_liquid_enthalpy
         le_slope = h_evap_l * cas%rho * tissue%g_transp_leaf(i) * dqdt
         le_ref   = h_evap_l * cas%rho * tissue%g_transp_leaf(i) * (qsat_c - qcas)
         le_slope_wet = h_film_l * cas%rho * film%g_film_leaf(i) * dqdt
         le_ref_wet   = h_film_l * cas%rho * film%g_film_leaf(i) * (qsat_c - qcas)
         !----- ARK-diagnostic leaf: emission base = t_cas, no storage (t_emit = tcas, store_hcap_per_dt = 0).   !
         !      qwflux_wl (sapflow's advected enthalpy, sec 2/6, P2) folds in via q_extra -- it shifts   !
         !      the equilibrium temperature (and hence dh/transp) like any other absorbed energy, but    !
         !      is kept OUT of drnet (an internal soil<->leaf transfer, not a boundary radiative input;   !
         !      see veg_energy_balance's own doc-comment). 0.0 when unset (every existing caller), so  !
         !      this is a no-op unless build_column_frozen populates it. --------------------------------!
         call veg_energy_balance(tissue%abs_sw(i), tissue%abs_lw(i), tissue%h_coeff_leaf(i), le_slope,          &
                                    lw_slope, le_ref, tcas, tcas, tissue%leaf_hcap_per_dt(i), tissue%t_leaf0(i),  &
                                    dtl, tl, transp_i, dh, drnet, q_extra=tissue%qwflux_wl(i),        &
                                    f_wet=film%f_wet_c(i), le_slope_wet=le_slope_wet,               &
                                    le_ref_wet=le_ref_wet, film_evap=f%film_evap_leaf(i),          &
                                    h_evap=h_evap_l, h_evap_wet=h_film_l)
         f%leaf_temp(i) = tl
         f%transp_c(i)  = transp_i                                          ! per-cohort demand (pre src_frac)
         coh_h      = coh_h      + dh
         coh_qw     = coh_qw     + (transp_i + f%film_evap_leaf(i)) * h_evap_l   ! what the leaf paid
         coh_transp = coh_transp + transp_i
         coh_film_evap = coh_film_evap + f%film_evap_leaf(i)
         coh_rnet   = coh_rnet   + drnet
         !----- Diagnostic WOOD balance (own store; emission base = tcas, no transpiration). Wood        !
         !      sensible + net-LW join coh_h / coh_rnet; a diagnostic wood has no storage so the two     !
         !      wood terms are equal (h_coeff_w*dtw) and telescope in the ledger. Frozen wood inputs are !
         !      zero when wood is not diagnostic (build_column_frozen), making this a no-op then.        !
         lw_slope_w = lw_emission_slope(tissue%leaf_emiss, tcas, tissue%wai(i))
         h_evap_w = h_evap_l
         h_film_w = h_evap_w - film%film_liquid_enthalpy
         le_slope_wet_w = h_film_w * cas%rho * film%g_film_w(i) * dqdt
         le_ref_wet_w   = h_film_w * cas%rho * film%g_film_w(i) * (qsat_c - qcas)
         !----- Diagnostic WOOD = the le_slope = le_ref = 0 case of the same kernel (no transp).       !
         !      q_wood_net (qloss - qwflux_wl, sec 2/6, P2) folds in via q_extra the same way qwflux_wl   !
         !      does for leaf above (kept out of drnet) -- 0.0 when unset. ------------------------------!
         call veg_energy_balance(tissue%abs_sw_wood(i), tissue%abs_lw_wood(i), tissue%h_coeff_w(i),          &
                                    0.0_wp, lw_slope_w, 0.0_wp, tcas, tcas, tissue%wood_hcap_per_dt(i),        &
                                    tissue%t_wood0(i),                                                &
                                    dtw, tw, transp_w, dh, drnet, q_extra=tissue%q_wood_net(i),       &
                                    f_wet=film%f_wet_c(i), le_slope_wet=le_slope_wet_w,             &
                                    le_ref_wet=le_ref_wet_w, film_evap=f%film_evap_wood(i),        &
                                    h_evap=h_evap_w, h_evap_wet=h_film_w)
         f%wood_temp(i) = tw
         coh_h    = coh_h    + dh
         coh_qw   = coh_qw   + f%film_evap_wood(i) * h_evap_w   ! wood film-evap -> CAS, what the wood paid
         coh_film_evap = coh_film_evap + f%film_evap_wood(i)
         coh_rnet = coh_rnet + drnet
      end do

      !----- GROUND SURFACE = snowfac-blended snow + (1-snowfac) bare soil (C4, issue #76). The snow  !
      !      terms come from the shared pre-column stage (meds_fast_snow) and are ALREADY snowfac-     !
      !      weighted; the bare-soil sensible is scaled by (1-snowfac). soil_evap already carries its  !
      !      own (1-snowfac) AREA factor via chydro_forcing_t%snow_free_frac, so le_soil is an         !
      !      area-integrated tile flux and is added whole -- do NOT scale it again.                    !
      !                                                                                                !
      !      With no snow every added term is 0 and snowfac = 0, so this reduces to                    !
      !        h_ground = h_bare, le_ground = le_soil, g_top = rad - h_bare - le_soil                  !
      !      which is the pre-C4 expression EXACTLY -- snow-off bit-identity is structural here, not   !
      !      a property to re-verify. snow%ground_rad is seeded to abs_sw_ground + abs_lw_ground by     !
      !      build_column_frozen for the same reason. -------------------------------------------------!
      call ground_surface_fluxes(t_ground, tcas, ground%ggnet, cas%rho, ground%soil_evap, h_bare, le_soil)
      f%h_ground  = snow%h_snow  + (1.0_wp - snow%snowfac) * h_bare
      f%le_ground = snow%le_snow + le_soil
      f%g_top     = snow%g_base                                                                   &
                    + (1.0_wp - snow%snowfac) * (ground%abs_sw_ground + ground%abs_lw_ground - h_bare)       &
                    - le_soil

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + coh_film_evap + ground%soil_evap + snow%subl_rate
      f%coh_rnet   = coh_rnet
      f%coh_transp = coh_transp

      !----- SMOOTH condensation sink (dew/fog): when the CAS is supersaturated, condense the excess as a  !
      !      continuous relaxation of q toward qsat with timescale TAU_COND -- C1-smooth except a mild kink !
      !      at RH=1 (vs the removed HARD state clamp). The condensed vapour leaves as liquid at Tcas: the  !
      !      CAS loses vapour mass (src_vap) AND only the LIQUID enthalpy (src_enth -= cond*u_liq); the     !
      !      latent heat (enthalpy_vapor - u_liq) STAYS in the air (warming it, as real condensation does). !
      !      Replaces the discontinuous clamp -> keeps the ARK2 embedded error smooth through saturation.   !
      !      §8g: gated so a split-vs-ARK comparison can hold the MODEL fixed (the split never reaches
      !      this routine, so leaving it always-on made the two schemes different models, not just
      !      different integrators).
      if (cas%cas_condensation) then
         f%cond     = (cas%cas_mass_capacity / TAU_COND) * max(0.0_wp, y%cas_shv - qsat_c)
      else
         f%cond     = 0.0_wp
      end if
      f%src_vap  = f%src_vap  - f%cond
      f%cond_enth = f%cond * internal_energy_liquid(tcas)
      f%src_enth = f%src_enth - f%cond_enth

      !----- CAS box tendencies (shared kernel; the condensation adjustment above is folded into    !
      !      src_enth/src_vap, so the box math is scheme-agnostic). ---------------------------------!
      cas_src%surface_enthalpy_source = f%src_enth
      cas_src%surface_vapor_source    = f%src_vap
      cas_src%biotic_co2_source       = cas%nee_biotic
      cas_col%air_mass_capacity        = cas%cas_mass_capacity
      cas_col%air_molar_capacity       = cas%cas_molar_capacity
      !----- The CAS<->atm conductances follow the LIVE canopy-air state.  This is the RK45 path's     !
      !      whole story -- its CAS tendency is built here, so this call is what unfreezes the           !
      !      Monin-Obukhov feedback for that scheme.  ARK arrives with mo_live already cleared on its   !
      !      stage copy, so it reuses the value it solved once for this stage rather than re-solving on  !
      !      every one of the Newton's residual evaluations. ---------------------------------------------!
      call cas_conductances(cas, y%cas_enthalpy, y%cas_shv, gah_l, gaw_l, gac_l)
      cas_col%atm_conductance_enthalpy = gah_l
      cas_col%atm_conductance_vapor    = gaw_l
      cas_col%atm_conductance_co2      = gac_l
      cas_col%atm_enthalpy             = cas%enthalpy_atm
      cas_col%atm_specific_humidity    = cas%shv_atm
      cas_col%atm_co2                  = cas%co2_atm
      call cas_column_time_deriv(y%cas_enthalpy, y%cas_shv, y%cas_co2, cas_src, cas_col,        &
                                 f%d_cas_enthalpy, f%d_cas_shv, f%d_cas_co2)
   end subroutine surface_derivs

   !---------------------------------------------------------------------------------------!
   ! column_derivs -- the WHOLE-column RHS. Diagnoses the soil-top temperature from the state (so    !
   ! the ground skin couples to the current soil-top energy), runs the surface block, then assembles !
   ! the soil-heat, soil-water, and per-cohort plant-water-MASS tendencies from the frozen forcing +   !
   ! the surface couplings (qloss -> root heat sink, the FROZEN aggregate uptake -> root water     !
   ! sink, transp_c -> the per-cohort transpiration demand). Commits nothing. This mirrors             !
   ! column_fast_step's operator sequence WITHOUT the backward-Euler denominators -- it is the         !
   ! tendency an IMEX-ARK/RK45 stage evaluates.                                                        !
   !                                                                                          !
   ! PLANT WATER MASS (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2): unlike the retired plant_water_tendency  !
   ! (which re-evaluated the full nonlinear PV-curve/conductance system -- psi's own stiff ODE -- at    !
   ! every stage), the mass ODE is a TRIVIAL affine expression: frozen%plant%sapflow_frozen/uptake_frozen are     !
   ! ONE frozen pair of numbers (the Act-1 pre-pass's time-averaged solve_plant_water output, held       !
   ! constant across the whole macro-step); only the per-plant transpiration demand is REFRESHED each    !
   ! stage, from the CURRENT surface_derivs evaluation -- exactly the sec 6 stability argument (mass     !
   ! adds no stiff mode because its inflow is frozen and its outflow moves at the CAS timescale). ------!
   pure subroutine column_derivs(y, frozen, n, nsl, f, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n      !< number of cohorts
      integer(ik),           intent(in)  :: nsl    !< number of active soil layers
      type(column_tend_t),   intent(out) :: f
      !----- this stage's OWN surface tendencies (coh_rnet, src_enth/src_vap, ...), for a fully-        !
      !      explicit multi-stage integrator (RK45) to b-weight into its OWN whole-column boundary-      !
      !      flux ledger -- mirrors column_be_stage's sf_out (meds_fast_ark.f90), same rationale. --------!
      type(surface_tend_t),  intent(out), optional :: sf_out

      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_column_t) :: soil_e
      type(energy_forcing_t)     :: eforc
      real(wp)                   :: t_ground, fliq1, wmass1, root_uptake(n_soil_layer_max)
      real(wp)                   :: transp_i, qloss_total, e_drain
      real(wp)                   :: qface_own(n_soil_layer_max)
      integer(ik)                :: k, i

      allocate(f%d_leaf_water_mass(n), f%d_wood_water_mass(n), f%leaf_temp(n))
      allocate(f%d_leaf_surf_water(n), f%d_wood_surf_water(n))

      !----- Diagnose the soil-top temperature from the current state so the ground skin sees the   !
      !      prognostic soil-top energy (the coupling the surface block needs). ---------------------!
      wmass1   = y%theta(1) * rho_h2o
      call internal_energy_to_temp(y%soil_energy(1), wmass1, frozen%params%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      !----- 1. Surface block (leaf + ground + CAS twins). ------------------------------------!
      y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,    &
                          t_ground, n, surf_tend)
      f%d_cas_enthalpy = surf_tend%d_cas_enthalpy ; f%d_cas_shv = surf_tend%d_cas_shv ; f%d_cas_co2 = surf_tend%d_cas_co2
      f%g_top = surf_tend%g_top ; f%leaf_temp(1:n) = surf_tend%leaf_temp(1:n)

      !----- 3. Soil-water (Richards) column: frozen top flux + the FROZEN aggregate root-uptake     !
      !      sink (frozen%roots%uptake, sec 3 -- the Act-1 pre-pass's plant-side request, already rescaled by    !
      !      the soil's own realized supply), NOT the stage-refreshed surf_tend%coh_transp -- the soil forcing  !
      !      must be the SAME frozen number the mass ODE below debits from wood_water_mass, or the two    !
      !      sides of the wood<->soil interface no longer cancel to machine precision. -----------------!
      !      frozen%roots%uptake is advance_soil_water_column's uptake_total, which ALREADY carries the psi-wilting  !
      !      ramp (face_and_sink applied f_wilt_ramp inside the scratch solve). Passing it back through  !
      !      the ramp here limited it a second time whenever psi_soil < psi_open, so the soil lost        !
      !      frozen%roots%uptake*fwilt while wood_water_mass gained frozen%roots%uptake -- water created from nothing on    !
      !      dry soil. apply_wilt_limit=.false. takes the sink as-is. ---------------------------------!
      do k = 1_ik, nsl
         root_uptake(k) = frozen%roots%uptake * frozen%roots%root_share(k)
      end do
      call soil_water_time_deriv(y%theta, frozen%params%soil, frozen%params%hydro_opts, nsl,                 &
                                 frozen%hydrology%q_top,                                                    &
                               root_uptake, f%dtheta_dt, f%drainage_rate, f%uptake_rate, qface_own, &
                               apply_wilt_limit=.false.)

      !----- 2. Soil-heat column: g_top from the surface, root heat sink = qloss (the advected enthalpy  !
      !      of the water the roots extract, at the root-weighted soil temperature, sec 2/6 P2),          !
      !      distributed by the static root_share profile. The transpired water's liquid enthalpy is      !
      !      NOT charged to the soil any more (the old coh_qsoil proxy): the leaf pays the full vapour    !
      !      enthalpy itself in surface_derivs, as ED2 does. qloss_frozen sums to 0 when the P2 wiring    !
      !      is unset (every caller besides build_column_frozen), so this is a no-op there. --------------!
      soil_e%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      !----- Faces, drainage and the bottom-face enthalpy all ride THIS stage's OWN water tendency     !
      !      (issue #78 item 3, C2): RK45 integrates its own theta, so the scratch solve's time-mean     !
      !      faces (w_flux_frozen) and its clip/floor corrections are the WRONG numbers here -- on a      !
      !      saturated column the two trajectories' interior faces differed by ~2.8 kg/m2 per step and   !
      !      the borrowed clip cooling happened to cancel it, which is how the two defects hid each other !
      !      while the soil surface sat at 345 K. The whole-column ledger cannot see a purely VERTICAL    !
      !      misplacement; only the faces' provenance protects against it. --------------------------!
      qloss_total = sum(frozen%roots%qloss_frozen(1:n))
      e_drain = f%drainage_rate * internal_energy_liquid(frozen%hydrology%t_bot)
      call assemble_soil_energy_forcing(eforc, nsl, surf_tend%g_top, frozen%hydrology%geothermal, y%theta,          &
                                        frozen%roots%root_share, qloss_total, qface_own,                        &
                                        frozen%hydrology%infiltration, frozen%hydrology%t_infil, e_drain)
      call soil_energy_time_deriv(soil_e, eforc, frozen%params%therm, frozen%params%soil, frozen%params%energy_opts, f%dedt)

      !----- 4. Per-cohort plant WATER MASS: frozen sapflow/uptake (Act 1) in, REFRESHED per-plant   !
      !      transpiration demand out -- see the header. transp_i mirrors exactly the conversion the    !
      !      retired plant_water_tendency call used (surf_tend%transp_c is per-m2-ground; per-plant divides     !
      !      by nplant). the aggregate uptake is the soil's realized supply (the instantaneous supply           !
      !      throttle is retired, matching the split path's own P0 design -- the plant's mass STORAGE     !
      !      absorbs any soil-supply/demand mismatch instead of throttling transp itself). --------------!
      do i = 1_ik, n
         transp_i = surf_tend%transp_c(i) / max(frozen%plant%nplant(i), tiny_num)
         f%d_leaf_water_mass(i) = frozen%plant%sapflow_frozen(i) - transp_i
         f%d_wood_water_mass(i) = frozen%plant%uptake_frozen(i)  - frozen%plant%sapflow_frozen(i)
      end do

      !----- 5. Canopy-SURFACE water (sec 3.4, P2c): frozen interception (Act 1) in, REFRESHED         !
      !      per-stage film evaporation out -- the surface-film analogue of the internal mass ODE        !
      !      just above. Already per-m2-ground on both sides (unlike leaf_water_mass, no nplant          !
      !      division needed: frozen%film%intercept_leaf/wood and surf_tend%film_evap_leaf/wood share that convention). !
      !      Zero when canopy_water_on is off (frozen%film%intercept_leaf/wood and surf_tend%film_evap_leaf/wood all      !
      !      0.0 then), so this is a no-op unless build_column_frozen populates it. --------------------!
      f%d_leaf_surf_water(1:n) = frozen%film%intercept_leaf(1:n) - surf_tend%film_evap_leaf(1:n)
      f%d_wood_surf_water(1:n) = frozen%film%intercept_wood(1:n) - surf_tend%film_evap_wood(1:n)
      if (present(sf_out)) sf_out = surf_tend
   end subroutine column_derivs

end module meds_fast_time_derivs
