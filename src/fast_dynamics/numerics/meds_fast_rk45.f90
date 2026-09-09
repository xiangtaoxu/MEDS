!==========================================================================================!
! meds_fast_rk45 -- the ED2-faithful adaptive Cash-Karp RK45 fast-loop integrator            !
! (docs/dev_plans/MEDS_ED2_RK45_DESIGN.md, phase P2). Peer of meds_fast_ark: dispatched by     !
! meds_fast_step%column_fast_step when cfg%time_integrator == INTEG_RK45. Unlike the ARK's       !
! IMEX-ESDIRK stages (implicit CAS+soil, operator-split mass), RK45 is FULLY EXPLICIT over the   !
! SAME pure RHS meds_fast_time_derivs%column_derivs used by the test-only RK4 oracle -- CAS,      !
! soil energy, soil water, and plant water mass are ALL genuinely integrated by the Cash-Karp      !
! stages (no operator split at all), which is why mass rides the embedded-error WRMS here          !
! (with_mass=.true.) where it is deliberately excluded on the ARK path.                             !
!                                                                                          !
! Stability (design doc sec 6): with the Act-1 pre-pass's fluxes frozen, the stiffest              !
! integrated mode is the CAS (tau ~ 130 s), so an explicit method needs dt <~ 360 s -- a few         !
! adaptive substeps at dt_fast = 900 s, not the microsecond steps a naively-explicit stiff solve      !
! would otherwise need. The plant water mass ODE adds no additional stiffness (its inflow is         !
! frozen and its outflow moves at the CAS timescale).                                               !
!                                                                                          !
! Shares build_column_frozen with the ARK path (meds_fast_ark) -- the Act-1 pre-pass (leaf gas       !
! exchange, aerodynamics, plant-hydraulics frozen sapflow/uptake, advective enthalpy) does not        !
! depend on which stepper advances the macro-step, so it is reused verbatim, not re-derived.          !
!==========================================================================================!
module meds_fast_rk45
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, cp_liq, cp_air
   use meds_therm_lib,           only : cas_temp_of_enthalpy, internal_energy_liquid, internal_energy_to_temp
   use meds_soil_water,       only : pond_overflow
   use meds_fast_time_derivs, only : surface_derivs, column_derivs, cas_conductances
   use meds_fast_types,       only : column_state_t, column_frozen_t, column_tend_t, error_control_t, &
                                     surface_state_t, surface_tend_t,          &
                                     column_config_t, column_cohort_t, column_forcing_t,         &
                                     column_budget_t, mask_is_full
   use meds_fast_ark,         only : build_column_frozen
   use meds_column_state_ops, only : state_init, state_axpy, state_accum, state_sub, zero_like, apply_process_mask,     &
                                     clamp_theta, clamp_cas, clamp_soil_energy, soil_water_store,   &
                                     soil_energy_store, plant_water_store, canopy_film_store,       &
                                     deposit_condensate, clamp_canopy_film, unpack_column_state,     &
                                     diagnose_soil_temps
   use meds_fast_control,     only : state_wrms_grouped,                                         &
                                     step_control_factor
   use meds_config,           only : meds_config_t, CTRL_L2_STRICT
   use meds_biophysics_types, only : aero_env_t, aero_geom_t, aero_out_t, patch_biophys_t
   use meds_biophysics_opts, only : SOIL_BC_AQUIFER
   use meds_budget_check,     only : budget_check, budget_energy_rate_floor, budget_water_rate_floor
   implicit none
   private

   public :: rk45_column_step, adaptive_rk45_march, column_fast_step_rk45
   public :: rk45_state_railed

   !----- Cash-Karp embedded 5(4) tableau (Cash & Karp 1990, ACM TOMS 16:201; the SAME          !
   !      coefficients as Numerical Recipes' rkck). c_i (stage times) are documentation only --   !
   !      the frozen forcing frozen has no explicit time dependence, so they never enter the code.   !
   real(wp), parameter :: A21 = 0.2_wp
   real(wp), parameter :: A31 = 3.0_wp/40.0_wp,      A32 = 9.0_wp/40.0_wp
   real(wp), parameter :: A41 = 0.3_wp,               A42 = -0.9_wp,             A43 = 1.2_wp
   real(wp), parameter :: A51 = -11.0_wp/54.0_wp,     A52 = 2.5_wp,              A53 = -70.0_wp/27.0_wp, &
                          A54 = 35.0_wp/27.0_wp
   real(wp), parameter :: A61 = 1631.0_wp/55296.0_wp, A62 = 175.0_wp/512.0_wp,   A63 = 575.0_wp/13824.0_wp, &
                          A64 = 44275.0_wp/110592.0_wp, A65 = 253.0_wp/4096.0_wp
   !----- 5th-order (committed, "local extrapolation") and embedded 4th-order b-vectors; b2=b*2=0,  !
   !      b5=0 (5th order does not use k5). -----------------------------------------------------!
   real(wp), parameter :: B1 = 37.0_wp/378.0_wp,    B3 = 250.0_wp/621.0_wp, B4 = 125.0_wp/594.0_wp, &
                          B6 = 512.0_wp/1771.0_wp
   real(wp), parameter :: BS1 = 2825.0_wp/27648.0_wp, BS3 = 18575.0_wp/48384.0_wp,                &
                          BS4 = 13525.0_wp/55296.0_wp, BS5 = 277.0_wp/14336.0_wp, BS6 = 0.25_wp
   !----- Embedded-pair lower order (for step_control_factor's -1/(p+1) exponent, sec 6). ---------!
   integer(ik), parameter :: RK45_P_ORDER = 4_ik

   !----- WORK BUDGET per dt_fast (P6, MEDS_ED2_RK45_DESIGN.md). Section 6's stability estimate puts   !
   !      normal operation at ~3 accepted sub-steps; a healthy step never approaches this. When the     !
   !      explicit surface enters its stiff (dense-canopy + cold) regime it instead thrashes to         !
   !      hundreds of sub-steps and STILL rails -- work that is pure waste, since the dispatcher then    !
   !      discards the step and redoes it on the implicit-CAS split path anyway. Bailing out at the cap  !
   !      turns that runaway into an early, cheap "this step is stiff" signal: it caps the wasted        !
   !      explicit work at ~20x normal instead of ~150x, which is the difference between a 30-yr run     !
   !      finishing and a single cold month taking longer than the rest of the run combined. -----------!
   integer(ik), parameter, public :: RK45_WORK_CAP = 64_ik

contains

   !---------------------------------------------------------------------------------------!
   ! stage_bnd -- the whole-column BOUNDARY-flux rates at one stage's state (analogous to           !
   ! column_be_stage's stage_bflux_t in meds_fast_ark, but for a fully-explicit stage: no BE          !
   ! solve, so the "state" IS just the stage's own input y_stage, and surf_tend is column_derivs' own surface     !
   ! diagnostic at y_stage). cond_enth is pre-multiplied by u_liq at THIS stage's own tcas (mirrors ARK's   !
   ! t_cas1 reference, sec 3.4/9's "one flux, both sides" -- the CAS's own reference, not a frozen     !
   ! one, since condensation genuinely happens at the evolving CAS temperature). ---------------------!
   pure subroutine stage_bnd(y_stage, frozen, surf_tend, rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i, &
                             atm_heat_i)
      type(column_state_t),  intent(in)  :: y_stage
      type(column_frozen_t), intent(in)  :: frozen
      type(surface_tend_t),  intent(in)  :: surf_tend
      real(wp),               intent(out) :: rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i
      real(wp),               intent(out) :: atm_heat_i   !< [W/m2] sensible export, the reported H
      real(wp) :: gah_i, gaw_i, gac_i
      rnet_i      = surf_tend%coh_rnet
      !----- The boundary flux must be charged at the conductance the TENDENCY used. column_derivs   !
      !      built this stage's CAS tendency from a live-state surface-layer re-solve, so reading the !
      !      state^n frozen%cas%gah here would book a boundary flux the state update never took -- the  !
      !      "borrow one solve's flux while committing another's state" defect class this project has !
      !      already paid for three times.                                                             !
      !                                                                                          !
      !      RE-SOLVED rather than read back off surf_tend. Reporting it on surface_tend_t is the more        !
      !      obviously-safe design (one number, both sides) and was tried first, but adding fields to  !
      !      that type perturbs ifx's inlining/FMA choices inside surface_derivs, which cost the       !
      !      then-default frozen path its bit-for-bit identity (~1e-12) and would have broken the      !
      !      git-stash/cmp verification protocol this project relies on. cas_conductances is `pure`    !
      !      and both call sites pass the same (frozen, y_stage) pair, so the recompute returns the same       !
      !      number. -----------------------------------------------------------------------------------!
      call cas_conductances(frozen%cas, y_stage%cas_enthalpy, y_stage%cas_shv, gah_i, gaw_i, gac_i)
      atm_enth_i  = gah_i * (y_stage%cas_enthalpy - frozen%cas%enth_atm)
      atm_vap_i   = gaw_i * (y_stage%cas_shv      - frozen%cas%shv_atm)
      atm_heat_i  = gah_i * cp_air * (cas_temp_of_enthalpy(y_stage%cas_enthalpy, y_stage%cas_shv)          &
                                       - frozen%cas%mo_theta_atm)
      cond_i      = surf_tend%cond
      cond_enth_i = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
   end subroutine stage_bnd

   !---------------------------------------------------------------------------------------!
   ! rk45_column_step -- one Cash-Karp 5(4) step of size dt from y, over the pure RHS            !
   ! column_derivs. Commits the 5TH-order solution (local extrapolation, matching how              !
   ! adaptive_ark_march/adaptive_imex_march both commit their higher-order result); y_err is the    !
   ! (5th - 4th) embedded difference for the adaptive controller. w_out/e_in/e_out are the           !
   ! whole-column boundary-flux AMOUNTS over dt, b-weighted by the SAME 5th-order b-vector as the    !
   ! state commit (the consistent quadrature for a boundary integral over this step). e_in's          !
   ! infiltration term (frozen%hydrology%infiltration*u_liq(rain_temp)) IS the whole-column precip-energy          !
   ! input -- the caller must NOT also add a separate forc%precip term on top (double-counts nearly   !
   ! the full infiltrating share whenever infiltration ~= precip); mirrors ARK's own bf%whole_enth_in, !
   ! which folds e_infil in the same way with no further outer addition. -----------------------------!
   ! STAGE CLAMPING (mirrors ark2_column_step's base3 clamp_theta/clamp_cas): every stage's y_stage is       !
   ! explicit-only here (no ESDIRK stabilization), so a stiff surface<->soil coupling under a too-large  !
   ! trial dt (a sparse/near-bare patch: small ground-CAS heat capacity relative to its coupling          !
   ! conductance) can drive theta/soil_energy far outside their physical domain WITHIN a single stage      !
   ! evaluation -- ground_evaporation's fractional pow() then hits a domain error (negative base) rather   !
   ! than merely a large-but-finite value the controller could reject. Clamping y_stage before each             !
   ! column_derivs call keeps every stage evaluation finite (an in-range y_stage is untouched, so this is a      !
   ! no-op on any well-resolved step); the resulting k_i is then a legitimate derivative at a physical      !
   ! (if boundary-pinned) state, and the normal 5th/4th embedded-error comparison rejects and shrinks dt     !
   ! exactly as it would for any other oversized step -- no separate detection logic needed. ---------------!
   pure subroutine rk45_column_step(y, frozen, n, nsl, dt, y_out, y_err, w_out, e_in, e_out,          &
                                    clamp_stage_n, clamp_commit_n, clamp_mass, clamp_energy, cond_out, &
                                    tissue_leaf_int, tissue_wood_int, cond_enth_out,                   &
                                    atm_heat_out, atm_vap_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      real(wp),               intent(out) :: w_out, e_in, e_out
      !----- CLAMP telemetry (see column_budget_t%clamp_*). STAGE and COMMIT are reported separately    !
      !      because they mean opposite things: a stage clamp bounds a throwaway input and is expected   !
      !      on an oversized trial, while a commit clamp edits the kept state with no ledger entry.      !
      !      The caller decides what to do with each -- notably, it must discard the COMMIT tallies of   !
      !      a REJECTED trial, since that state is thrown away. -------------------------------------!
      integer(ik), optional,  intent(inout) :: clamp_stage_n
      integer(ik), optional,  intent(out)   :: clamp_commit_n
      real(wp),    optional,  intent(out)   :: clamp_mass, clamp_energy
      real(wp),    optional,  intent(out)   :: cond_out   !< [kg/m2] condensate to deposit (row 1b)
      real(wp),    optional,  intent(out)   :: cond_enth_out !< [J/m2] its liquid enthalpy, b-weighted at the stage CAS temps
      real(wp), optional,     intent(out) :: atm_heat_out, atm_vap_out   !< [J/m2],[kg/m2] net CAS->atm sensible/vapour over dt
      real(wp),    optional,  intent(out)   :: tissue_leaf_int(n), tissue_wood_int(n)  !< [K*s]

      type(column_tend_t)  :: k1, k2, k3, k4, k5, k6
      type(column_state_t) :: y_stage, y_4th
      type(surface_tend_t) :: surf_tend
      real(wp) :: rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6), atm_heat(6)
      real(wp) :: bw_rnet, bw_atm_enth, bw_atm_vap, bw_cond, bw_cond_enth, bw_drain, bw_atm_heat
      !----- Per-stage tissue temperatures, kept so the store's energy can be set from their TIME    !
      !      INTEGRAL rather than from the last stage -- see column_bflux_t's note. RK45 is the       !
      !      accuracy baseline, so paying n x 6 here is the right trade. ------------------------------!
      real(wp) :: tleaf_s(n,6), twood_s(n,6)

      call column_derivs(y, frozen, n, nsl, k1, sf_out=surf_tend)
      call stage_bnd(y, frozen, surf_tend, rnet(1), atm_enth(1), atm_vap(1), cond(1), cond_enth(1), atm_heat(1))
      tleaf_s(1:n,1) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,1) = surf_tend%wood_temp(1:n)

      call state_init(y, n, nsl, y_stage) ; call state_accum(y_stage, dt*A21, k1, n, nsl)
      call clamp_theta(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call clamp_soil_energy(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, n, nsl, k2, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(2), atm_enth(2), atm_vap(2), cond(2), cond_enth(2), atm_heat(2))
      tleaf_s(1:n,2) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,2) = surf_tend%wood_temp(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A31, k1, n, nsl) ; call state_accum(y_stage, dt*A32, k2, n, nsl)
      call clamp_theta(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call clamp_soil_energy(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, n, nsl, k3, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(3), atm_enth(3), atm_vap(3), cond(3), cond_enth(3), atm_heat(3))
      tleaf_s(1:n,3) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,3) = surf_tend%wood_temp(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A41, k1, n, nsl) ; call state_accum(y_stage, dt*A42, k2, n, nsl)
      call state_accum(y_stage, dt*A43, k3, n, nsl)
      call clamp_theta(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call clamp_soil_energy(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, n, nsl, k4, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(4), atm_enth(4), atm_vap(4), cond(4), cond_enth(4), atm_heat(4))
      tleaf_s(1:n,4) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,4) = surf_tend%wood_temp(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A51, k1, n, nsl) ; call state_accum(y_stage, dt*A52, k2, n, nsl)
      call state_accum(y_stage, dt*A53, k3, n, nsl) ; call state_accum(y_stage, dt*A54, k4, n, nsl)
      call clamp_theta(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call clamp_soil_energy(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, n, nsl, k5, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(5), atm_enth(5), atm_vap(5), cond(5), cond_enth(5), atm_heat(5))
      tleaf_s(1:n,5) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,5) = surf_tend%wood_temp(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A61, k1, n, nsl) ; call state_accum(y_stage, dt*A62, k2, n, nsl)
      call state_accum(y_stage, dt*A63, k3, n, nsl) ; call state_accum(y_stage, dt*A64, k4, n, nsl)
      call state_accum(y_stage, dt*A65, k5, n, nsl)
      call clamp_theta(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call clamp_soil_energy(y_stage, frozen, nsl, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, n, nsl, k6, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6), atm_heat(6))
      tleaf_s(1:n,6) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,6) = surf_tend%wood_temp(1:n)

      !----- y_out = y + dt*(B1*k1 + B3*k3 + B4*k4 + B6*k6)  [5th order; b2=b5=0]. -----------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt*B1, k1, n, nsl)
      call state_accum(y_out, dt*B3, k3, n, nsl)
      call state_accum(y_out, dt*B4, k4, n, nsl)
      call state_accum(y_out, dt*B6, k6, n, nsl)
      !----- The COMMITTED state is deliberately NOT clamped (C1, MEDS_INTEGRATOR_PARITY.md [RETIRED] row 7).    !
      !      It used to be. The argument for clamping it was that a clamp which bites shows up as a      !
      !      large 5th-vs-4th discrepancy and the controller then rejects the step -- but that does NOT  !
      !      cover the accept test, which is `err <= 1 .OR. dt <= dt_floor`: at the sub-step floor a      !
      !      clamped state was committed unconditionally, with the controller out of moves, and whatever !
      !      mass or energy the clamp moved was kept with no ledger entry. Measured on the saturated      !
      !      column: 142.7 kg/m2 of water and (on nvfortran) 1.5e5 J/m2 of energy over 96 sub-steps.     !
      !                                                                                                  !
      !      The STAGE clamps above stay -- they bound throwaway inputs so column_derivs stays evaluable  !
      !      (ground_evaporation's fractional pow() needs a non-negative base), and a discarded stage     !
      !      state breaks no books. What is gone is editing the state that is kept.                      !
      !                                                                                                  !
      !      An out-of-range commit is now handled the way every other bad step is: the embedded error    !
      !      grows, the controller rejects and shrinks; and if a floor-forced step still commits a railed !
      !      state, rk45_state_railed catches it at the dispatch and the hybrid rescue redoes the whole   !
      !      dt_fast on the split path. That path already existed -- the clamp was pre-empting it with a  !
      !      silent correction instead of letting it fire. -----------------------------------------------!
      if (present(clamp_commit_n)) clamp_commit_n = 0_ik
      if (present(clamp_mass))     clamp_mass     = 0.0_wp
      if (present(clamp_energy))   clamp_energy   = 0.0_wp

      !----- y_4th = y + dt*(BS1*k1 + BS3*k3 + BS4*k4 + BS5*k5 + BS6*k6)  [embedded 4th order]. --!
      !      y_4th is a SEPARATE named temporary, not y_err itself: state_sub's `out` dummy has     !
      !      intent(out), so aliasing it with the `b` argument (both = y_err) is undefined -- the    !
      !      compiler may treat intent(out) as non-aliased and clobber b's value (incl. deallocating !
      !      its allocatable components) before the subtraction reads it. Confirmed empirically:      !
      !      the aliased form gave a y_err that did not shrink with dt at all (a dt-INDEPENDENT        !
      !      garbage value), instead of the expected O(dt^5) embedded-error scaling. -------------------!
      call state_init(y, n, nsl, y_4th)
      call state_accum(y_4th, dt*BS1, k1, n, nsl)
      call state_accum(y_4th, dt*BS3, k3, n, nsl)
      call state_accum(y_4th, dt*BS4, k4, n, nsl)
      call state_accum(y_4th, dt*BS5, k5, n, nsl)
      call state_accum(y_4th, dt*BS6, k6, n, nsl)
      call state_sub(y_out, y_4th, n, nsl, y_err)   ! y_err := y_5th - y_4th (the embedded estimate)

      !----- whole-column boundary AMOUNTS over dt: b-weighted (5th-order vector) sums of the      !
      !      state-dependent per-stage rates, PLUS the frozen constants (added once, undiluted --   !
      !      sum(b)=1 for any consistent RK b-vector, so a CONSTANT rate integrates to rate*dt        !
      !      regardless of how it is spread across the b-weighted sum). --------------------------!
      !----- C2: DRAINAGE b-weighted from the stage tendencies, so the ledger books the water that     !
      !      RK45's OWN theta trajectory actually shed through the bottom face. It used to book          !
      !      frozen%hydrology%drainage -- the Act-1 scratch solve's frozen value -- while committing its own theta,   !
      !      so the two disagreed by exactly the amount the RK trajectory departed from the scratch      !
      !      solve. Invisible while unsaturated (the two nearly coincide) and the dominant term once     !
      !      the column saturates. Same b-vector as the state commit, for the same reason. -------------!
      bw_drain     = B1*k1%drainage_rate + B3*k3%drainage_rate + B4*k4%drainage_rate                  &
                     + B6*k6%drainage_rate
      bw_rnet      = B1*rnet(1)      + B3*rnet(3)      + B4*rnet(4)      + B6*rnet(6)
      bw_atm_enth  = B1*atm_enth(1)  + B3*atm_enth(3)  + B4*atm_enth(4)  + B6*atm_enth(6)
      bw_atm_vap   = B1*atm_vap(1)   + B3*atm_vap(3)   + B4*atm_vap(4)   + B6*atm_vap(6)
      bw_atm_heat  = B1*atm_heat(1)  + B3*atm_heat(3)  + B4*atm_heat(4)  + B6*atm_heat(6)
      bw_cond      = B1*cond(1)      + B3*cond(3)      + B4*cond(4)      + B6*cond(6)
      bw_cond_enth = B1*cond_enth(1) + B3*cond_enth(3) + B4*cond_enth(4) + B6*cond_enth(6)
      !----- Tissue-temperature TIME INTEGRAL over this sub-step, same b-vector as the state commit.  !
      !      sum(B) = 1, so dividing by the march length recovers a temperature. --------------------!
      if (present(tissue_leaf_int)) then
         tissue_leaf_int(1:n) = dt * (B1*tleaf_s(1:n,1) + B3*tleaf_s(1:n,3)                       &
                                      + B4*tleaf_s(1:n,4) + B6*tleaf_s(1:n,6))
         tissue_wood_int(1:n) = dt * (B1*twood_s(1:n,1) + B3*twood_s(1:n,3)                       &
                                      + B4*twood_s(1:n,4) + B6*twood_s(1:n,6))
      end if

      !----- NO clip / theta_res-floor term any more (#78 item 3). Those were the SCRATCH solve's        !
      !      post-solve corrections, and they entered here as whole-column BOUNDARY flux -- clipped      !
      !      water leaving for a pond that could hold no heat, floored water created from nowhere. Both  !
      !      halves are gone: the saturation excess is shed by the caller's commit guard, which routes  !
      !      it to the pond as an internal transfer (crossing no boundary), and the theta_res floor is    !
      !      likewise applied by the caller on its OWN committed state. --------------------------------!
      !----- ground_rad is the snowfac-BLENDED radiative input (= abs_sw+abs_lw when bare, C4). --!
      !----- #78 item 4: the infiltration enthalpy is a POND -> SOIL transfer between two tracked      !
      !      stores now, so it is no longer a boundary input here -- the caller's pond store absorbs    !
      !      the other half. Leaving it produced exactly the infiltration enthalpy as a spurious        !
      !      surplus (7.4e4 J/m2 on the wet fixture, which is precip*dt*u_liq to three digits). --------!
      e_in  = (bw_rnet + frozen%snow%ground_rad) * dt
      e_out = bw_atm_enth * dt                                                                    &
              + bw_drain * dt * internal_energy_liquid(frozen%hydrology%t_bot)
      !----- cond is EXCLUDED from w_out (row 1b): it is deposited into soil layer 1 by the caller, !
      !      not lost across the boundary. cond_out returns the amount for that deposit. ----------!
      !----- runoff is NOT booked here (#75): the caller rebuilds the ponding store from RK45's own !
      !      trajectory and derives its own overflow, so taking the frozen scratch runoff too would  !
      !      double-count it. Drainage IS RK45's own (b-weighted above). --------------------------!
      w_out = bw_atm_vap * dt + bw_drain * dt
      if (present(cond_out))      cond_out      = bw_cond * dt
      if (present(cond_enth_out)) cond_enth_out = bw_cond_enth * dt
      if (present(atm_heat_out))  atm_heat_out  = bw_atm_heat * dt
      if (present(atm_vap_out))   atm_vap_out   = bw_atm_vap  * dt
   end subroutine rk45_column_step

   !---------------------------------------------------------------------------------------!
   ! adaptive_rk45_march -- integrate to t_end with the Cash-Karp embedded error driving the        !
   ! step controller (reusing meds_fast_control's WRMS + I/PI controller + warm start, sec 6 --      !
   ! "no new controller"), mirroring adaptive_ark_march's accept/reject structure exactly, just       !
   ! over rk45_column_step instead of ark2_column_step. with_mass=.true. (the default): RK45 has no   !
   ! operator split at all, so mass genuinely differs between the 5th/4th solutions -- a live error    !
   ! signal, unlike ARK where it is structurally zero. -----------------------------------------------!
   subroutine adaptive_rk45_march(y0, frozen, n, nsl, t_end, ec, dt_init, y_out, nsteps, nrej,       &
                                  w_out_acc, e_in_acc, e_out_acc, dt_warm_out,                    &
                                  clamp_stage_n, clamp_commit_n, clamp_mass, clamp_energy, cond_acc,  &
                                  cond_enth_acc, atm_heat_acc, atm_vap_acc,                            &
                                  ood_max, tissue_leaf_acc, tissue_wood_acc)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec
      type(column_state_t),  intent(out) :: y_out
      integer(ik),            intent(out) :: nsteps, nrej
      real(wp),               intent(out) :: w_out_acc, e_in_acc, e_out_acc
      real(wp),    optional,  intent(out) :: dt_warm_out
      !----- CLAMP telemetry over the whole march. The two kinds accumulate on DIFFERENT populations:  !
      !      stage clamps count on every trial (a rejected trial's clamp still says the controller      !
      !      probed too far -- that is the signal), while commit clamps count only on ACCEPTED steps,   !
      !      because a rejected step's clamped state is discarded and never breaks any book. ----------!
      integer(ik), optional,  intent(inout) :: clamp_stage_n, clamp_commit_n
      real(wp),    optional,  intent(inout) :: clamp_mass, clamp_energy
      real(wp),    optional,  intent(inout) :: cond_acc   !< [kg/m2] accumulated condensate (row 1b)
      real(wp),    optional,  intent(inout) :: cond_enth_acc !< [J/m2] ...and its stage-valued liquid enthalpy
      real(wp), optional,     intent(inout) :: atm_heat_acc, atm_vap_acc   !< [J/m2],[kg/m2] net CAS->atm export, accepted sub-steps
      !----- [K*s] per-cohort tissue-temperature time integrals over the accepted march. Caller zeroes.!
      real(wp),    optional,  intent(inout) :: tissue_leaf_acc(n), tissue_wood_acc(n)
      !----- [m3/m3] running MAX theta excursion outside [theta_res, theta_sat] at a stage-1 RHS input   !
      !      (issue #78 item 2 telemetry -- see column_budget_t%theta_ood_max). Measured on EVERY trial, !
      !      accepted or not: a rejected trial's out-of-domain probe is still a place the RHS was        !
      !      evaluated, which is exactly what this is reporting. ---------------------------------------!
      real(wp),    optional,  intent(inout) :: ood_max

      type(column_state_t) :: y, y_new, y_err, y_zero
      real(wp) :: t, dt, err, err_prev, fac, dt_floor
      real(wp) :: w_out, e_in, e_out, dt_try, dt_warm
      real(wp) :: cmass_i, cenergy_i, cond_i, cond_enth_i
      real(wp)    :: atm_heat_i, atm_vap_i
      real(wp) :: tl_int_i(n), tw_int_i(n)
      integer(ik) :: ccommit_i, kood
      logical  :: clamped

      w_out_acc = 0.0_wp ; e_in_acc = 0.0_wp ; e_out_acc = 0.0_wp
      !----- substep FLOOR: bound the worst case to ~t_end/64 sub-steps (mirrors adaptive_ark_       !
      !      march's own floor and rationale -- a floor step is not guaranteed L-stable here          !
      !      (RK45 is explicit), but the sec 6 stability estimate already bounds normal operation      !
      !      to ~3 substeps, so a stiffness-driven runaway signals a genuinely pathological step,      !
      !      not routine behaviour; degrade gracefully rather than grind at a tiny floor forever). ---!
      dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)

      call state_init(y0, n, nsl, y)
      t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik
      err_prev = -1.0_wp
      dt_warm = dt
      do
         if (t >= t_end - tiny_num) exit
         dt_try = dt
         dt = min(dt, t_end - t)
         clamped = dt < dt_try - tiny_num
         !----- stage-1 reads y UNCLAMPED (clamp_theta guards stages 2-6 only), so this is the one      !
         !      place a constitutive kernel sees theta outside its domain. Record how far. -------------!
         if (present(ood_max)) then
            do kood = 1_ik, nsl
               ood_max = max(ood_max, y%theta(kood) - frozen%params%soil%theta_sat(kood),                       &
                             frozen%params%soil%theta_res(kood) - y%theta(kood))
            end do
         end if
         call rk45_column_step(y, frozen, n, nsl, dt, y_new, y_err, w_out, e_in, e_out,               &
                               clamp_stage_n=clamp_stage_n, clamp_commit_n=ccommit_i,              &
                               clamp_mass=cmass_i, clamp_energy=cenergy_i, cond_out=cond_i,       &
                               tissue_leaf_int=tl_int_i, tissue_wood_int=tw_int_i,                &
                               cond_enth_out=cond_enth_i, atm_heat_out=atm_heat_i, atm_vap_out=atm_vap_i)
         !----- named temporary: never pass a derived-type-valued function result straight into a  !
         !      call (the nvfortran whole-program-optimizer trap documented in CLAUDE.md). --------!
         y_zero = zero_like(y_err, n, nsl)
         err = state_wrms_grouped(y_err, y_zero, y_new, n, nsl, ec%tols)
         if (err /= err .or. dt /= dt) then    ! non-finite: too big a step, reject deterministically
            if (dt <= dt_floor) then
               !----- floor-forced commit of a non-finite trial: this is the WORST case for the commit  !
               !      clamp (the controller is out of moves), so it must be tallied like any accept. ---!
               if (present(clamp_commit_n)) clamp_commit_n = clamp_commit_n + ccommit_i
               if (present(clamp_mass))     clamp_mass     = clamp_mass     + cmass_i
               if (present(clamp_energy))   clamp_energy   = clamp_energy   + cenergy_i
               if (present(cond_acc))       cond_acc       = cond_acc       + cond_i
               if (present(cond_enth_acc))  cond_enth_acc  = cond_enth_acc  + cond_enth_i
            if (present(atm_heat_acc))   atm_heat_acc   = atm_heat_acc   + atm_heat_i
            if (present(atm_vap_acc))    atm_vap_acc    = atm_vap_acc    + atm_vap_i
               call state_init(y_new, n, nsl, y) ; t = t + dt_floor ; nsteps = nsteps + 1_ik ; exit
            end if
            nrej = nrej + 1_ik ; dt = max(dt * ec%fmin, dt_floor) ; cycle
         end if
         fac = step_control_factor(err, err_prev, ec)
         if (err <= 1.0_wp .or. dt <= dt_floor) then
            if (ec%level == CTRL_L2_STRICT .and. err > 1.0_wp) &
               error stop 'adaptive_rk45_march: L2 strict -- floor step cannot meet tolerance'
            call state_init(y_new, n, nsl, y)
            !----- this step's clamped state is now the committed state, so its unbookkept mass/energy  !
            !      joins the running tally (a REJECTED trial's is discarded with the trial). -----------!
            if (present(clamp_commit_n)) clamp_commit_n = clamp_commit_n + ccommit_i
            if (present(clamp_mass))     clamp_mass     = clamp_mass     + cmass_i
            if (present(clamp_energy))   clamp_energy   = clamp_energy   + cenergy_i
            if (present(cond_acc))       cond_acc       = cond_acc       + cond_i
            if (present(cond_enth_acc))  cond_enth_acc  = cond_enth_acc  + cond_enth_i
            if (present(atm_heat_acc))   atm_heat_acc   = atm_heat_acc   + atm_heat_i
            if (present(atm_vap_acc))    atm_vap_acc    = atm_vap_acc    + atm_vap_i
            w_out_acc = w_out_acc + w_out ; e_in_acc = e_in_acc + e_in ; e_out_acc = e_out_acc + e_out
            !----- Accumulate over ACCEPTED sub-steps only, exactly like every amount above. --------!
            if (present(tissue_leaf_acc)) then
               tissue_leaf_acc(1:n) = tissue_leaf_acc(1:n) + tl_int_i(1:n)
               tissue_wood_acc(1:n) = tissue_wood_acc(1:n) + tw_int_i(1:n)
            end if
            t = t + dt ; nsteps = nsteps + 1_ik
            err_prev = err
            if (.not. clamped) dt_warm = dt
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         !----- WORK-BUDGET bail-out (P6): stop burning explicit sub-steps on a step the dispatcher    !
         !      is going to discard anyway. y_out is left at the partial-time state, which is FINE --   !
         !      the caller returns immediately without committing it (stiff_bail). ---------------------!
         if (nsteps + nrej >= RK45_WORK_CAP) exit
      end do
      call state_init(y, n, nsl, y_out)
      if (present(dt_warm_out)) dt_warm_out = dt_warm
   end subroutine adaptive_rk45_march


   !=======================================================================================!
   !  INTEG_RK45 path: the ED2-faithful adaptive Cash-Karp fast step. Shares build_column_frozen  !
   !  with the ARK path (meds_fast_ark) -- the Act-1 pre-pass does not depend on which stepper      !
   !  advances the macro-step. -----------------------------------------------------------------!
   !=======================================================================================!
   subroutine column_fast_step_rk45(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,  &
                                    gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, &
                                    stiff_bail, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: col_config
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters
      !----- P6: .true. when the explicit march hit its work budget (RK45_WORK_CAP) without resolving   !
      !      the step -- the stiff-regime signal. biophys is then left UNTOUCHED at state^n and the caller    !
      !      must redo this dt_fast on the implicit-CAS split path. -------------------------------------!
      logical,     optional,   intent(out)   :: stiff_bail
      real(wp),    optional,   intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) per-cohort diagnostic capture

      type(column_frozen_t)  :: frozen
      type(column_state_t)   :: y, y_out
      type(surface_state_t)  :: y_stage
      type(surface_tend_t)   :: surf_tend
      real(wp)    :: dt0, wcap, enth0, shv0, enth1, shv1
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, w_plant0, w_plant1, w_surface0
      real(wp)    :: w_out_acc, e_in_acc, e_out_acc, w_in, w_out, e_in, e_out
      real(wp)    :: tg, fl, dt_warm_next, cond_dep, cond_dep_enth
      real(wp)    :: clip_mass_rk, clip_enth_rk, dm_clip, w_pond_rk, runoff_rk
      real(wp)    :: floor_mass_rk, floor_enth_rk, dm_floor  ! #78 item 3: the theta_res floor at commit
      real(wp)    :: e_pond0, e_pond_rk, over_enth_rk   ! #78 item 4
      !----- Canopy-SURFACE water (sec 3.4, P2c) ledger scratch. --------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, leaf_cap_i, wood_cap_i, intercept_total
      real(wp)    :: tissue_store0, tissue_store1
      real(wp)    :: tl_int_acc(col_cohort%n), tw_int_acc(col_cohort%n)
      real(wp)    :: cap_leaf_a(col_cohort%n), cap_wood_a(col_cohort%n)
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, k, i, nsteps, nrej
      logical     :: halt_budgets

      n = col_cohort%n ; nsl = col_config%soil%n_active
      !----- The bottom-BC guard is GONE (Phase 0/3). The aquifer BC no longer carries a storage      !
      !      bucket or a water-table state for this path to borrow: it is a head-driven boundary flux,  !
      !      and RK45 gets it through soil_water_time_deriv on its OWN theta like every other face. ----!
      !----- CLAMP counters accumulate down the call chain, so this sub-step's tally starts clean. ----!
      budget%clamp_stage_n = 0_ik ; budget%clamp_commit_n = 0_ik
      budget%clamp_mass    = 0.0_wp ; budget%clamp_energy = 0.0_wp
      budget%theta_ood_max = 0.0_wp
      !----- state^n ponding store, captured BEFORE the unpack below overwrites it. ------------------!
      w_surface0 = biophys%soil_w%w_surface
      e_pond0    = biophys%soil_w%w_surface_enth   ! #78 item 4
      halt_budgets = col_config%energy%debug_error .and. mask_is_full(col_config%mask)
      if (present(stiff_bail)) stiff_bail = .false.

      cond_dep = 0.0_wp ; cond_dep_enth = 0.0_wp
      budget%atm_heat_export = 0.0_wp ; budget%atm_vap_export = 0.0_wp
      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, n, nsl, &
                               frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      dt0 = dt_fast
      if (biophys%adapt_dt_last > tiny_num) dt0 = min(biophys%adapt_dt_last, dt_fast)
      if (col_config%integrator%dt_init > tiny_num) dt0 = min(col_config%integrator%dt_init, dt_fast)
      ec = col_config%integrator%error_control
      ec%p_order = RK45_P_ORDER
      tl_int_acc(1:n) = 0.0_wp ; tw_int_acc(1:n) = 0.0_wp
      call adaptive_rk45_march(y, frozen, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                              w_out_acc, e_in_acc, e_out_acc, dt_warm_out=dt_warm_next,           &
                              clamp_stage_n=budget%clamp_stage_n, clamp_commit_n=budget%clamp_commit_n, &
                              clamp_mass=budget%clamp_mass, clamp_energy=budget%clamp_energy,          &
                              cond_acc=cond_dep, cond_enth_acc=cond_dep_enth,                     &
                              atm_heat_acc=budget%atm_heat_export, atm_vap_acc=budget%atm_vap_export, &
                              ood_max=budget%theta_ood_max,                                         &
                              tissue_leaf_acc=tl_int_acc, tissue_wood_acc=tw_int_acc)
      biophys%adapt_dt_last = dt_warm_next
      budget%integ_nsteps = nsteps ; budget%integ_nrej = nrej

      !----- P6 WORK-BUDGET bail: the march gave up (stiff regime). Return NOW, before committing      !
      !      anything into biophys or running the budget ledgers -- y_out is a partial-time state, so       !
      !      committing it would corrupt the column and trip the conservation checks. biophys is untouched  !
      !      at state^n, so the caller can cleanly redo this dt_fast on the split path. -----------------!
      if (nsteps + nrej >= RK45_WORK_CAP) then
         if (present(stiff_bail)) stiff_bail = .true.
         if (present(converged))  converged  = .false.
         if (present(iters))      iters      = nsteps
         return
      end if

      !----- §5.1 PROCESS MASK: a masked-off component is restored to state^n (y). RK45 genuinely   !
      !      integrates soil water AND plant mass (unlike ARK, where soil water is fully operator-    !
      !      split and mass is separately split too), so BOTH honour the mask here directly. -----------!
      call apply_process_mask(col_config%mask, y, y_out, n, nsl)

      !----- Canopy-SURFACE water (sec 3.4, P2c): capacity clamp + overflow bookkeeping (mirrors the    !
      !      split path's own post-hoc treatment, sec 9's "clamp, don't silently over-apply"; the ODE     !
      !      itself is left unclamped mid-integration, matching the general lesson that the closure-safe   !
      !      clamp point is a store's BOUNDARY, not a flux feeding an already-self-consistent solve).       !
      !      Gated behind canopy_water_on per the P1 nvfortran lesson (gate new arithmetic behind its        !
      !      own flag from the start), not just relied on to telescope to a no-op. ------------------------!
      !      DEFICIT (the floor, symmetric to overflow): film_evap uses a state^n-frozen conductance      !
      !      (rescaled by availability in build_column_frozen, but only an approximation -- evaporative     !
      !      demand can still grow through the step as tcas/leaf_temp evolve), so the store can still        !
      !      transiently overdraw below 0. Left uncorrected, a negative store corrupts the NEXT dt_fast's      !
      !      frozen interception rate (intercept_canopy_layer's own internal floor silently "fixes" a           !
      !      negative starting bucket, fabricating mass that was never really lost). Floor at 0 and bookkeep     !
      !      the shortfall as a NEGATIVE addition to w_out/e_out (the store APPEARED to gain from clamping        !
      !      up to 0, so the ledger's outflow must shrink by the same amount to match) -- the exact mirror         !
      !      of surf_overflow's sign. ---------------------------------------------------------------------------!
      surf_overflow = 0.0_wp ; surf_deficit = 0.0_wp
      if (col_config%canopy_water_on) then
         call clamp_canopy_film(y_out, col_cohort%lai, col_cohort%wai, col_config%hydro%dewmx, n, surf_overflow, surf_deficit)
      end if

      !----- unpack into biophys + re-derive the diagnostic soil/leaf/wood temperatures. -----------!
      call unpack_column_state(y_out, n, nsl, biophys)
      !----- Ponding / aquifer / water-table stores, mirroring column_fast_step_ark. RK45 was DROPPING  !
      !      them: the scratch advance_soil_water_column computes the end-of-step pond, but nothing wrote    !
      !      it back and the whole_water ledger below carried no w_surface term either, so any water     !
      !      that ponded left the tracked stores without appearing in any flux. Invisible while the      !
      !      pond stays empty (every other test here is free-draining at theta = 0.30) and exactly the   !
      !      size of the pond once it fills -- 7.5 kg/m2 in the saturated test that caught it. Gated on  !
      !      the same §5.1 soil-water mask ARK uses, for the same reason: these ARE the soil-water       !
      !      store, so they must freeze with theta or the two schemes run different reduced systems. ----!
      if (col_config%mask%soil_water) then
         biophys%soil_w%w_surface = frozen%hydrology%w_surface1
      end if
      !----- ROW 1b: DEPOSIT THE CONDENSATE (same routing as column_fast_step_ark, for the full         !
      !      rationale). Dew/fog landed on a surface inside the column; it used to be booked into        !
      !      w_out and vanish. Paired mass + enthalpy into soil layer 1 at the CAS temperature, so the   !
      !      whole-column ledger closes with no boundary term. Same destination on all three paths --    !
      !      routing it anywhere else on one path would reopen a scheme asymmetry. --------------------!
      if (cond_dep > 0.0_wp) then
         !----- at the b-weighted STAGE enthalpy the CAS was debited, not u_liq(T_end) (item 1A #6). --!
         call deposit_condensate(y_out, col_config%soil%dz(1), cond_dep, cond_dep_enth)
         biophys%soil_e%soil_energy(1) = y_out%soil_energy(1) ; biophys%soil_w%theta(1) = y_out%theta(1)
      end if
      call diagnose_soil_temps(y_out, col_config%soil_thermal%soil_dry_heat_capacity, nsl, biophys%soil_e%soil_temp, &
            biophys%soil_e%soil_fliq)
      !----- CONSTITUTIVE-DOMAIN GUARDS on RK45's committed theta, both paired transfers.               !
      !                                                                                                 !
      !      An explicit method has no post-solve hook, so unlike the implicit sibling these guards are  !
      !      the only place theta can be brought back inside [theta_res, theta_sat] before it is          !
      !      committed -- which is what keeps the van Genuchten curves in domain on the next step.        !
      !                                                                                                 !
      !      A saturation-RELIEF TENDENCY (-max(0, theta-theta_sat)/tau inside the stages, so the mass    !
      !      never oversaturates in the first place) was built and measured as the alternative, and       !
      !      REJECTED: with the interior faces below wired correctly it changed nothing observable --     !
      !      soil-surface peak 295.10 vs 295.05 K, whole suite green, zero stage clamps either way. Its   !
      !      only real effect was a ~27x smaller in-stage excursion past theta_sat, which is issue #78    !
      !      item 2's concern (out-of-domain constitutive evaluation), not item 3's. Shipping a new       !
      !      timescale parameter and a new tendency term for no measurable gain is how the borrowed       !
      !      frozen%hydrology%clip_enth got here in the first place; if a production run ever shows a large commit     !
      !      clip, add relief then, with that measurement in hand.                                       !
      !                                                                                                 !
      !      1. SATURATION (upper): the excess moves to the pond with its enthalpy at the layer's own    !
      !         temperature. Now that the pond has a thermal state (#78 item 4) this is a transfer       !
      !         between two tracked stores, not a boundary loss.                                         !
      !      2. RESIDUAL (lower): theta_res is the OTHER edge of the constitutive domain, and RK45 had   !
      !         no guard on it at all -- it compensated the scratch solve's floor ENTHALPY               !
      !         (frozen%hydrology%floor_enth) without ever applying the floor to its own theta. Deleting that         !
      !         compensation without adding the floor would leave Se < 0 reachable. Unlike the clip      !
      !         this one CREATES water, which cannot be a transfer from anywhere, so it is booked as a   !
      !         boundary INPUT (floor_mass_rk/floor_enth_rk -> w_in/e_in) rather than hidden. The split  !
      !         path leaves the same correction in its mass RESIDUAL; making it an explicit ledger term  !
      !         here means a run that leans on it shows up in the books instead of in the noise.         !
      !                                                                                                 !
      !      Dunne runoff needs no term: f_sat is nonzero only under SOIL_BC_AQUIFER, which this path    !
      !      hard-errors on (C5), so RK45's runoff is purely pond overflow. Frozen infiltration is       !
      !      likewise correct rather than a compromise -- advance_soil_water_column computes infl from       !
      !      state^n BEFORE its own solve, so the split path freezes it identically. -------------------!
      clip_mass_rk = 0.0_wp ; clip_enth_rk = 0.0_wp
      floor_mass_rk = 0.0_wp ; floor_enth_rk = 0.0_wp
      if (col_config%mask%soil_water) then
         do k = 1_ik, nsl
            if (y_out%theta(k) > col_config%soil%theta_sat(k)) then
               dm_clip = (y_out%theta(k) - col_config%soil%theta_sat(k)) * col_config%soil%dz(k) * rho_h2o
               clip_mass_rk = clip_mass_rk + dm_clip
               clip_enth_rk = clip_enth_rk + dm_clip * internal_energy_liquid(biophys%soil_e%soil_temp(k))
               y_out%soil_energy(k) = y_out%soil_energy(k)                                            &
                                    - dm_clip * internal_energy_liquid(biophys%soil_e%soil_temp(k))       &
                                      / col_config%soil%dz(k)
               y_out%theta(k)       = col_config%soil%theta_sat(k)
               biophys%soil_w%theta(k)       = y_out%theta(k)
               biophys%soil_e%soil_energy(k) = y_out%soil_energy(k)
            else if (y_out%theta(k) < col_config%soil%theta_res(k)) then
               dm_floor = (col_config%soil%theta_res(k) - y_out%theta(k)) * col_config%soil%dz(k) * rho_h2o
               floor_mass_rk = floor_mass_rk + dm_floor
               floor_enth_rk = floor_enth_rk + dm_floor * internal_energy_liquid(biophys%soil_e%soil_temp(k))
               y_out%soil_energy(k) = y_out%soil_energy(k)                                            &
                                    + dm_floor * internal_energy_liquid(biophys%soil_e%soil_temp(k))      &
                                      / col_config%soil%dz(k)
               y_out%theta(k)       = col_config%soil%theta_res(k)
               biophys%soil_w%theta(k)       = y_out%theta(k)
               biophys%soil_e%soil_energy(k) = y_out%soil_energy(k)
            end if
         end do
         !----- REBUILD the pond from RK45's OWN trajectory (#75), rather than inheriting            !
         !      frozen%hydrology%w_surface1. That frozen value is the SCRATCH solve's end-of-step pond and already !
         !      contains the scratch's own saturation clip -- mass RK45's theta never shed. Adding    !
         !      RK45's own clip on top of it counted that water twice. The composition below is       !
         !      advance_soil_water_column's own, evaluated on this path's numbers: what could not         !
         !      infiltrate, plus what this trajectory's own theta had to shed at the saturation guard. !
         !      q_over (Dunne) is identically 0 here (the aquifer BC is head-driven, no saturated area). !
         w_pond_rk   = w_surface0 + (frozen%hydrology%precip_ground - frozen%hydrology%infiltration) * dt_fast + clip_mass_rk
         !----- ...and its ENTHALPY on the SAME trajectory, term for term (#78 item 4). Composing the    !
         !      pond's mass from one trajectory and its enthalpy from another is the defect class this    !
         !      whole issue is about, so the enthalpy mirrors the mass line above term for term: precip   !
         !      in at the temperature the kernel used, infiltration out at the pond temperature it        !
         !      reported, and the commit clip in at the layer temperatures it was valued from. -----------!
         e_pond_rk   = e_pond0                                                                        &
                     + frozen%hydrology%precip_ground * dt_fast * internal_energy_liquid(frozen%hydrology%t_precip)             &
                     - frozen%hydrology%infiltration  * dt_fast * internal_energy_liquid(frozen%hydrology%t_infil)              &
                     + clip_enth_rk
         !----- overflow + empty-pond reset through the SAME kernel advance_soil_water_column uses (step 4):  !
         !      the overflow carries the pond's mean specific enthalpy, not u_liq of the plateau-pinned  !
         !      read-off temperature (2026-09 winter residual). The kernel speaks in RATES over dt_fast;  !
         !      this ledger books AMOUNTS, hence the *dt_fast. -----------------------------------------!
         call pond_overflow(w_pond_rk, e_pond_rk, dt_fast, col_config%hydro%w_pond_max, frozen%hydrology%t_precip,   &
                            runoff_rk, over_enth_rk)
         runoff_rk    = runoff_rk    * dt_fast
         over_enth_rk = over_enth_rk * dt_fast
         biophys%soil_w%w_surface      = w_pond_rk
         biophys%soil_w%w_surface_enth = e_pond_rk
      else
         w_pond_rk = frozen%hydrology%w_surface1 ; runoff_rk = 0.0_wp
         e_pond_rk = frozen%hydrology%w_surface_enth1 ; over_enth_rk = 0.0_wp
      end if

      !----- Commit the tissue temperatures from the b-weighted stage time integrals (below); the      !
      !      frozen store relaxed them from frozen%tissue%t_leaf0/t_wood0 against a = cap/dt inside every    !
      !      stage, so nothing is corrected afterwards. ----------------------------------------------!
      tissue_store0 = 0.0_wp ; tissue_store1 = 0.0_wp
      do i = 1_ik, n
         !----- Derive the capacity from the SAME a_store the kernel relaxed against, not by         !
         !      recomputing it from dry_hcap + wmass. The two agree by construction today, but only    !
         !      this form guarantees that zeroing a_leaf/a_wood zeroes the ledger's store term too --  !
         !      i.e. that "no capacity" is a clean no-op end to end rather than a state change the     !
         !      fluxes never paid for. --------------------------------------------------------------!
         cap_leaf_a(i) = frozen%tissue%a_leaf(i) * dt_fast
         cap_wood_a(i) = frozen%tissue%a_wood(i) * dt_fast
         tissue_store0 = tissue_store0 + cap_leaf_a(i) * frozen%tissue%t_leaf0(i)                          &
                                       + cap_wood_a(i) * frozen%tissue%t_wood0(i)
      end do
      !----- Commit the TIME-AVERAGE over the march, not the last stage -- identical treatment to    !
      !      the ARK path, and for the identical reason (see column_bflux_t's note). -----------------!
      if (col_config%mask%veg_energy) then
         biophys%leaf_temp(1:n) = tl_int_acc(1:n) / dt_fast
         biophys%wood_temp(1:n) = tw_int_acc(1:n) / dt_fast
      else
         biophys%leaf_temp(1:n) = frozen%tissue%t_leaf0(1:n) ; biophys%wood_temp(1:n) = frozen%tissue%t_wood0(1:n)
      end if
      do i = 1_ik, n
         tissue_store1 = tissue_store1 + cap_leaf_a(i) * biophys%leaf_temp(i)                             &
                                       + cap_wood_a(i) * biophys%wood_temp(i)
      end do

      !----- WHOLE-COLUMN CONSERVATION LEDGER (design doc sec 8 gates 2/3 -- the headline           !
      !      deliverable): unlike ARK's per-kernel + whole ledger, RK45 checks the WHOLE-COLUMN        !
      !      water/energy budgets only (soil water and plant mass are genuinely integrated states,      !
      !      not operator-split, so there is no separate "soil_water (rk45)"/frozen-flux kernel          !
      !      check the way ARK needs one -- the whole-column ledger IS the individual-store ledger      !
      !      here, since every store advances through the SAME column_derivs RHS). -------------------!
      wcap = frozen%cas%wcap
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv   ! AFTER the prognostic-wood CAS credit above
      e_soil0 = soil_energy_store(y%soil_energy,     col_config%soil%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, col_config%soil%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     col_config%soil%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, col_config%soil%dz, nsl)
      w_plant0 = plant_water_store(col_cohort%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(col_cohort%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- Canopy-SURFACE water (sec 3.4, P2c): already ground-area-referenced (no nplant factor,     !
      !      unlike w_plant0/1 above). Valued at u_liq(rain_temp) = frozen%film%film_u_ref; the tissue pays   !
      !      enthalpy_vapor - film_u_ref per kg of film it evaporates (surface_derivs), so the store       !
      !      closes exactly against the CAS credit. All zero when canopy_water_on is off. ----------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%rain_temp)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%rain_temp)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))

      !----- e_in is e_in_acc ALONE -- NOT e_in_acc + a separate forc%precip energy term. Precip's       !
      !      energy already enters the ledger via rk45_column_step's OWN per-substep e_infil            !
      !      (frozen%hydrology%infiltration*u_liq(rain_temp), b-weighted into e_in_acc), which is the SAME frozen      !
      !      quantity feeding column_derivs' root_heat_sink(1) -- i.e. what the SOIL state actually        !
      !      receives. Adding a second, independent forc%precip*u_liq(cas_temp) term here (as an           !
      !      earlier version of this line did) double-counts nearly the full infiltrating share            !
      !      whenever infiltration ~= precip (the common, non-runoff case) -- mirrors ARK's own             !
      !      whole_energy ledger, which uses acc%whole_enth_in (e_infil baked in via bf%whole_enth_in)       !
      !      directly, with no further outer precip addition. w_in stays forc%precip*dt_fast (unlike        !
      !      e_in, w_out_acc has no infiltration-side counterpart to double against). The INTERCEPTED       !
      !      share (intercept_total) needs its own e_in term at the SAME rain_temp reference, mirroring      !
      !      the split path's own intercepted_total treatment -- 0 when canopy_water_on is off. w_out_acc/    !
      !      e_*_acc and w_in are all AMOUNTS over the whole dt_fast (budget_accumulate below uses dt=1),      !
      !      so intercept_total (a RATE) needs *dt_fast to match, while surf_overflow/surf_deficit (already    !
      !      amounts, in y_out's own units) need no such scaling. surf_deficit SUBTRACTS (the exact mirror       !
      !      of surf_overflow's sign -- flooring a negative store UP to 0 makes it appear to gain, so the         !
      !      ledger's outflow must shrink by the same amount to match). ------------------------------------------!
      !----- forc%snowf is a boundary water input that lands in the PACK (C4). It was absent here --   !
      !      split has always carried it -- so with a pack the ledger saw mass appear with no source.  !
      !----- floor_mass_rk (#78 item 3) is water CREATED by the theta_res guard above. It has no source  !
      !      inside the column, so it enters as a boundary input rather than as a silent correction --    !
      !      0 on any column that never dried past theta_res. --------------------------------------------!
      w_in  = (forc%precip + forc%snowf + biophys%shed_water_rate) * dt_fast + floor_mass_rk ! P4: shed water is a boundary
                                                               ! input too; its energy needs NO separate
                                                               ! term here, for the SAME reason precip's
                                                               ! doesn't -- it rides e_in_acc via
                                                               ! rk45_column_step's own e_infil, once mixed
                                                               ! into hforc%precip_ground (build_column_frozen,
                                                               ! shared with ARK).
      !----- C2: RK45's OWN pond overflow (from its own clip), not the frozen scratch's runoff. ----!
      w_out = w_out_acc + surf_overflow - surf_deficit + runoff_rk
      e_in  = e_in_acc + intercept_total * dt_fast * internal_energy_liquid(frozen%hydrology%rain_temp)        &
              + frozen%snow%acc_enth + floor_enth_rk                                             &
              + (frozen%hydrology%precip_ground - frozen%snow%melt_rate) * dt_fast                             &
                * internal_energy_liquid(frozen%hydrology%t_precip)
      !----- #78 items 3+4: the commit clip's enthalpy does NOT appear here -- it is a soil -> pond      !
      !      transfer between two tracked stores, so it telescopes inside the ledger rather than         !
      !      crossing its boundary, and the SCRATCH solve's clip (which used to leave as boundary flux   !
      !      for mass this trajectory never shed) is gone entirely. What does leave is the pond           !
      !      OVERFLOW, at the pond's own temperature -- runoff carries real energy now that the water it  !
      !      drains has a temperature to carry. ------------------------------------------------------!
      e_out = e_out_acc + (surf_overflow - surf_deficit) * internal_energy_liquid(frozen%hydrology%rain_temp)     &
              + over_enth_rk

      !----- FLUX-scaled tolerances (meds_budget_check header), same rule as the ARK ledgers. --------!
      call budget_check(budget%whole_water,                                                            &
                        w_soil0 + wcap*shv0 + w_surface0 + w_plant0 + surf_water0                     &
                        + frozen%snow%swe0,                                                          &
                        w_soil1 + wcap*shv1 + w_pond_rk + w_plant1 + surf_water1                        &
                        + frozen%snow%swe1,                                                          &
                        w_in, w_out, dt_fast, budget_water_rate_floor, 'whole_water (rk45)', halt_budgets)
      !----- snow store + its accumulated precip enthalpy join the ledger (C4); 0 without snow. -----!
      call budget_check(budget%whole_energy,                                                           &
                        !----- No melt rebase any more (#78 item 4): the pack sends its meltwater to  !
                        !      the POND, not to soil layer 1, so the pack/pond pair telescopes on its !
                        !      own (same as the ARK path). ----------------------------------------!
                        e_soil0                           + wcap*enth0 + surf_enth0                    &
                        + frozen%snow%enth0 + e_pond0 + tissue_store0,                              &
                        e_soil1 + wcap*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond_rk            &
                        + tissue_store1,                                                                &
                        e_in, e_out, dt_fast, budget_energy_rate_floor, 'whole_energy (rk45)', halt_budgets)
      !----- NOT YET CHECKED: a per-kernel cas_co2 closure (ARK's own budget%cas_co2 check) would need  !
      !      a b-weighted per-stage CO2 atmospheric-exchange accumulation this first pass does not      !
      !      track (only rnet/atm_enth/atm_vap/cond are tracked in rk45_column_step) -- deferred; the   !
      !      whole-column water/energy ledger above are the design doc's headline gates (sec 8          !
      !      gates 2/3), not the per-kernel CAS/soil individual checks ARK additionally provides. -------!

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_rk45

   !----- Did the explicit march commit a CLAMP-PINNED state? The 5th/4th embedded pair can rail  !
   !      TOGETHER, which makes the error controller see err~0 and accept physically impossible   !
   !      temperatures; this is the second rescue trigger at the dispatch (meds_fast_step). Lives !
   !      here, with the scheme it describes, rather than in the dispatcher. --------------------!
   pure function rk45_state_railed(biophys, nsl) result(railed)
      type(patch_biophys_t), intent(in) :: biophys
      integer(ik),           intent(in) :: nsl
      logical :: railed
      real(wp), parameter :: T_LO = 185.0_wp, T_HI = 345.0_wp
      railed =      biophys%cas%can_temp <= T_LO .or. biophys%cas%can_temp >= T_HI                          &
               .or. any(biophys%soil_e%soil_temp(1:nsl) <= T_LO)                                        &
               .or. any(biophys%soil_e%soil_temp(1:nsl) >= T_HI)
   end function rk45_state_railed

end module meds_fast_rk45
