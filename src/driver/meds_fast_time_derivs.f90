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
! The surface hydrology BCs (q_top, psi_e, soil_psi_root) and the leaf gas-exchange pre-pass are    !
! the frozen, explicit part of the additive split (held constant across the ARK macro-step).       !
!==========================================================================================!
module meds_fast_time_derivs
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : latent_heat_vap, stefan, cp_air, tiny_num, rho_h2o
   use meds_therm_lib,           only : cas_temp_of_enthalpy, sat_specific_humidity,                    &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, uext_to_temp,       &
                                     internal_energy_liquid
   use meds_biophysics_types, only : n_soil_layer_max, soil_energy_column_t, energy_forcing_t,   &
                                     soil_thermal_params_t, soil_params_t, soil_opts_t, energy_opts_t
   use meds_soil_energy,      only : soil_energy_time_deriv
   use meds_soil_water,       only : soil_water_time_deriv
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_time_deriv
   use meds_ground_biophysics, only : ground_surface_fluxes
   use meds_vegetation_biophysics, only : veg_energy_diagnostic
   use meds_fast_types,       only : surface_state_t, surface_frozen_t, surface_tend_t,           &
                                     column_state_t, column_frozen_t, column_tend_t,               &
                                     stage_bflux_t, column_bflux_t
   implicit none
   private

   public :: surface_derivs, column_derivs, root_weighted_psi

   !----- CAS supersaturation relaxation timescale [s] for the smooth condensation sink: q relaxes to    !
   !      qsat with this e-folding time (dq/dt|_cond = -(q-qsat)/TAU_COND). Short vs dt_fast=1800s so a    !
   !      supersaturated canopy air condenses within the step, but SMOOTH (no clamp) so the ARK2 embedded  !
   !      error stays small. Physical dew/fog forms fast; the exact value is not critical (any strong     !
   !      restoring rate pins q near qsat). ----------------------------------------------------------!
   real(wp), parameter :: TAU_COND = 300.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! root-weighted mean of a per-layer quantity (e.g. psi_soil) by the static root_frac profile   !
   ! (assumed to sum to 1) -- the one authority for a formula the split and the ARK frozen        !
   ! pre-pass both compute after their own column_hydrology_flux call.                            !
   !---------------------------------------------------------------------------------------!
   pure function root_weighted_psi(psi_soil, root_frac, nsl) result(psi_root)
      real(wp),    intent(in) :: psi_soil(:), root_frac(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: psi_root
      integer(ik) :: k
      psi_root = 0.0_wp
      do k = 1_ik, nsl
         psi_root = psi_root + psi_soil(k) * root_frac(k)
      end do
   end function root_weighted_psi

   !---------------------------------------------------------------------------------------!
   ! surface_derivs -- the CAS surface-block RHS (leaf-energy diagnostic + ground skin + the three  !
   ! CAS twins). Faithful, side-effect-free transcription of column_fast_step's surface path         !
   ! (meds_fast_split.f90). The longwave emission base is the current tcas, so the   !
   ! split's `tcas - te` term is identically zero. Integrating d_cas_* with a backward-Euler-in-the-  !
   ! atmosphere step reproduces the split's committed enth1/shv1/co21 exactly.                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine surface_derivs(y, fro, n, f)
      type(surface_state_t),  intent(in)  :: y
      type(surface_frozen_t), intent(in)  :: fro
      integer(ik),            intent(in)  :: n
      type(surface_tend_t),   intent(out) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat
      real(wp)    :: lw_slope, le_slope, le_ref, dtl, tl, transp_i, dh, drnet
      real(wp)    :: lw_slope_w, dtw, tw, transp_w                    !< diagnostic WOOD balance (own store)
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet
      type(cas_source_t) :: cas_src
      type(cas_column_t) :: cas_col
      integer(ik) :: i

      allocate(f%leaf_temp(n), f%wood_temp(n), f%transp_c(n))

      tcas   = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, fro%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, fro%press)

      coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
      do i = 1_ik, n
         lw_slope = 4.0_wp * fro%leaf_emiss * stefan * tcas ** 3 * fro%lai(i)
         le_slope = latent_heat_vap * fro%rho * fro%g_tr_f(i) * dqdt
         le_ref   = latent_heat_vap * fro%rho * fro%g_tr_f(i) * (qsat_c - qcas)
         !----- ARK-diagnostic leaf: emission base = t_cas, no storage (t_emit = tcas, a_store = 0).   !
         !      qwflux_wl (sapflow's advected enthalpy, sec 2/6, P2) folds in via q_extra -- it shifts   !
         !      the equilibrium temperature (and hence dh/transp) like any other absorbed energy, but    !
         !      is kept OUT of drnet (an internal soil<->leaf transfer, not a boundary radiative input;   !
         !      see veg_energy_diagnostic's own doc-comment). 0.0 when unset (every existing caller), so  !
         !      this is a no-op unless build_column_frozen populates it. --------------------------------!
         call veg_energy_diagnostic(fro%abs_sw(i), fro%abs_lw(i), fro%h_coeff_f(i), le_slope,          &
                                    lw_slope, le_ref, tcas, tcas, 0.0_wp, tcas,                   &
                                    dtl, tl, transp_i, dh, drnet, q_extra=fro%qwflux_wl(i))
         f%leaf_temp(i) = tl
         f%transp_c(i)  = transp_i                                          ! per-cohort demand (pre src_frac)
         coh_h      = coh_h      + dh
         coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl)
         coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)
         coh_transp = coh_transp + transp_i
         coh_rnet   = coh_rnet   + drnet
         !----- Diagnostic WOOD balance (own store; emission base = tcas, no transpiration). Wood        !
         !      sensible + net-LW join coh_h / coh_rnet; a diagnostic wood has no storage so the two     !
         !      wood terms are equal (h_coeff_w*dtw) and telescope in the ledger. Frozen wood inputs are !
         !      zero when wood is not diagnostic (build_column_frozen), making this a no-op then.        !
         lw_slope_w = 4.0_wp * fro%leaf_emiss * stefan * tcas ** 3 * fro%wai(i)
         !----- Diagnostic WOOD = the le_slope = le_ref = 0 case of the same kernel (no transp).       !
         !      q_wood_net (qloss - qwflux_wl, sec 2/6, P2) folds in via q_extra the same way qwflux_wl   !
         !      does for leaf above (kept out of drnet) -- 0.0 when unset. ------------------------------!
         call veg_energy_diagnostic(fro%abs_sw_wood(i), fro%abs_lw_wood(i), fro%h_coeff_w(i),          &
                                    0.0_wp, lw_slope_w, 0.0_wp, tcas, tcas, 0.0_wp, tcas,         &
                                    dtw, tw, transp_w, dh, drnet, q_extra=fro%q_wood_net(i))
         f%wood_temp(i) = tw
         coh_h    = coh_h    + dh
         coh_rnet = coh_rnet + drnet
      end do
      coh_qw     = coh_qw     * fro%src_frac
      coh_qsoil  = coh_qsoil  * fro%src_frac
      coh_transp = coh_transp * fro%src_frac

      call ground_surface_fluxes(fro%t_ground, tcas, fro%ggnet, fro%rho, fro%soil_evap, f%h_ground, f%le_ground)
      f%g_top     = fro%abs_sw_ground + fro%abs_lw_ground - f%h_ground - f%le_ground

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + fro%soil_evap
      f%coh_rnet   = coh_rnet
      f%coh_qsoil  = coh_qsoil
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
      if (fro%cas_condensation) then
         f%cond     = (fro%wcap / TAU_COND) * max(0.0_wp, y%cas_shv - qsat_c)
      else
         f%cond     = 0.0_wp
      end if
      f%src_vap  = f%src_vap  - f%cond
      f%src_enth = f%src_enth - f%cond * internal_energy_liquid(tcas)

      !----- CAS box tendencies (shared kernel; the condensation adjustment above is folded into    !
      !      src_enth/src_vap, so the box math is scheme-agnostic). ---------------------------------!
      cas_src%surface_enthalpy_source = f%src_enth
      cas_src%surface_vapor_source    = f%src_vap
      cas_src%biotic_co2_source       = fro%nee_biotic
      cas_col%air_mass_capacity        = fro%wcap
      cas_col%air_molar_capacity       = fro%ccap
      cas_col%atm_conductance_enthalpy = fro%gah
      cas_col%atm_conductance_vapor    = fro%gaw
      cas_col%atm_conductance_co2      = fro%gac
      cas_col%atm_enthalpy             = fro%enth_atm
      cas_col%atm_specific_humidity    = fro%shv_atm
      cas_col%atm_co2                  = fro%co2_atm
      call cas_column_time_deriv(y%cas_enthalpy, y%cas_shv, y%cas_co2, cas_src, cas_col,        &
                                 f%d_cas_enthalpy, f%d_cas_shv, f%d_cas_co2)
   end subroutine surface_derivs

   !---------------------------------------------------------------------------------------!
   ! column_derivs -- the WHOLE-column RHS. Diagnoses the soil-top temperature from the state (so    !
   ! the ground skin couples to the current soil-top energy), runs the surface block, then assembles !
   ! the soil-heat, soil-water, and per-cohort plant-water-MASS tendencies from the frozen forcing +   !
   ! the surface couplings (coh_qsoil -> root heat sink, the FROZEN aggregate uptake -> root water     !
   ! sink, transp_c -> the per-cohort transpiration demand). Commits nothing. This mirrors             !
   ! column_fast_step's operator sequence WITHOUT the backward-Euler denominators -- it is the         !
   ! tendency an IMEX-ARK/RK45 stage evaluates.                                                        !
   !                                                                                          !
   ! PLANT WATER MASS (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2): unlike the retired plant_water_tendency  !
   ! (which re-evaluated the full nonlinear PV-curve/conductance system -- psi's own stiff ODE -- at    !
   ! every stage), the mass ODE is a TRIVIAL affine expression: fro%sapflow_frozen/uptake_frozen are     !
   ! ONE frozen pair of numbers (the Act-1 pre-pass's time-averaged solve_plant_water output, held       !
   ! constant across the whole macro-step); only the per-plant transpiration demand is REFRESHED each    !
   ! stage, from the CURRENT surface_derivs evaluation -- exactly the sec 6 stability argument (mass     !
   ! adds no stiff mode because its inflow is frozen and its outflow moves at the CAS timescale). ------!
   pure subroutine column_derivs(y, fro, n, nsl, f, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n      !< number of cohorts
      integer(ik),           intent(in)  :: nsl    !< number of active soil layers
      type(column_tend_t),   intent(out) :: f
      !----- this stage's OWN surface tendencies (coh_rnet, src_enth/src_vap, ...), for a fully-        !
      !      explicit multi-stage integrator (RK45) to b-weight into its OWN whole-column boundary-      !
      !      flux ledger -- mirrors column_be_stage's sf_out (meds_fast_ark.f90), same rationale. --------!
      type(surface_tend_t),  intent(out), optional :: sf_out

      type(surface_state_t)      :: ys
      type(surface_frozen_t)     :: fs
      type(surface_tend_t)       :: sf
      type(soil_energy_column_t) :: soil_e
      type(energy_forcing_t)     :: eforc
      real(wp)                   :: t_ground, fliq1, wmass1, root_uptake(n_soil_layer_max)
      real(wp)                   :: transp_i, qloss_total, e_infil, e_runof, e_drain
      integer(ik)                :: k, i

      allocate(f%d_leaf_water_mass(n), f%d_wood_water_mass(n), f%leaf_temp(n))

      !----- Diagnose the soil-top temperature from the current state so the ground skin sees the   !
      !      prognostic soil-top energy (the coupling the surface block needs). ---------------------!
      wmass1   = y%theta(1) * rho_h2o
      call uext_to_temp(y%soil_energy(1), wmass1, fro%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      !----- 1. Surface block (leaf + ground + CAS twins). ------------------------------------!
      ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
      fs = fro%surf ; fs%t_ground = t_ground
      call surface_derivs(ys, fs, n, sf)
      f%d_cas_enthalpy = sf%d_cas_enthalpy ; f%d_cas_shv = sf%d_cas_shv ; f%d_cas_co2 = sf%d_cas_co2
      f%g_top = sf%g_top ; f%leaf_temp(1:n) = sf%leaf_temp(1:n)

      !----- 2. Soil-heat column: g_top from the surface, root heat sink from the shed enthalpy         !
      !      (transpiration's coh_qsoil, pre-existing) PLUS qloss (uptake's advected enthalpy, sec        !
      !      2/6, P2) -- both distributed by the SAME static root_frac profile (matching how              !
      !      coh_qsoil already was); qloss_frozen sums to 0 when the P2 advective-enthalpy wiring is        !
      !      unset (every caller besides build_column_frozen), so this is a no-op there. ------------------!
      soil_e%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top      = sf%g_top
      eforc%geothermal = fro%geothermal
      qloss_total      = sum(fro%qloss_frozen(1:n))
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         eforc%root_heat_sink(k) = (sf%coh_qsoil + qloss_total) * fro%soil%root_frac(k)
         eforc%w_flux(k)         = 0.0_wp                       ! interior advection lumped (baseline)
      end do
      !----- boundary water-enthalpy advection (mirrors meds_fast_ark.f90's column_be_stage exactly):  !
      !      infiltrating throughfall carries its liquid enthalpy INTO the soil top, runoff/drainage    !
      !      carry it OUT (state^n temps, frozen). Without this, e_in/e_out (rk45_column_step) count     !
      !      infiltration/runoff/drainage as whole-column boundary flux, but the soil STATE never pays    !
      !      or gets paid for it -- a residual exactly equal to the omitted term (caught via RK45's        !
      !      whole-energy budget test with a nonzero background drainage). root_heat_sink is a SINK,        !
      !      so q_src = -sink/dz: subtract an inflow, add an outflow. -----------------------------------!
      e_infil = fro%infiltration * internal_energy_liquid(fro%rain_temp)
      e_runof = fro%runoff_surf  * internal_energy_liquid(fro%surf%t_ground)
      e_drain = fro%drainage     * internal_energy_liquid(fro%t_bot)
      eforc%root_heat_sink(1)   = eforc%root_heat_sink(1)   - e_infil + e_runof
      eforc%root_heat_sink(nsl) = eforc%root_heat_sink(nsl) + e_drain
      call soil_energy_time_deriv(soil_e, eforc, fro%therm, fro%soil, fro%energy_opts, f%dedt)

      !----- 3. Soil-water (Richards) column: frozen top flux + the FROZEN aggregate root-uptake     !
      !      sink (fro%uptake, sec 3 -- the Act-1 pre-pass's plant-side request, already rescaled by    !
      !      the soil's own realized supply), NOT the stage-refreshed sf%coh_transp -- the soil forcing  !
      !      must be the SAME frozen number the mass ODE below debits from wood_water_mass, or the two    !
      !      sides of the wood<->soil interface no longer cancel to machine precision. -----------------!
      do k = 1_ik, nsl
         root_uptake(k) = fro%uptake * fro%soil%root_frac(k)
      end do
      call soil_water_time_deriv(y%theta, fro%soil, fro%hydro_opts, nsl, fro%q_top, fro%psi_e,     &
                               root_uptake, f%dtheta_dt, f%drainage_rate, f%uptake_rate)

      !----- 4. Per-cohort plant WATER MASS: frozen sapflow/uptake (Act 1) in, REFRESHED per-plant   !
      !      transpiration demand out -- see the header. transp_i mirrors exactly the conversion the    !
      !      retired plant_water_tendency call used (sf%transp_c is per-m2-ground; per-plant divides     !
      !      by nplant). fro%surf%src_frac stays at its 1.0 default here (the instantaneous supply       !
      !      throttle is retired, matching the split path's own P0 design -- the plant's mass STORAGE     !
      !      absorbs any soil-supply/demand mismatch instead of throttling transp itself). --------------!
      do i = 1_ik, n
         transp_i = sf%transp_c(i) * fro%surf%src_frac / max(fro%nplant(i), tiny_num)
         f%d_leaf_water_mass(i) = fro%sapflow_frozen(i) - transp_i
         f%d_wood_water_mass(i) = fro%uptake_frozen(i)  - fro%sapflow_frozen(i)
      end do
      if (present(sf_out)) sf_out = sf
   end subroutine column_derivs

end module meds_fast_time_derivs
