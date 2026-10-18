! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_rk45 -- the ED2-style adaptive Cash-Karp RK45 fast-loop integrator, the        !
! accuracy baseline (docs/science/numerical_scheme.md §3 `rk45`). Peer of meds_fast_ark:   !
! dispatched by meds_fast_step%column_fast_step when cfg%time_integrator == INTEG_RK45.    !
! Unlike the ARK's ESDIRK stages (implicit CAS + soil heat, operator-split water), RK45 is !
! FULLY EXPLICIT over the pure RHS meds_fast_time_derivs%column_derivs that the test-only  !
! RK4 oracle also uses: CAS, soil energy, soil water and plant water mass all ride the     !
! Cash-Karp stages (the water at the network's constant rates), so the error norm scores   !
! all of them.                                                                             !
!                                                                                          !
! Stability: with the frozen record's fluxes held, the stiffest integrated mode is the CAS !
! (tau ~ 130 s), so an explicit method needs dt <~ 360 s -- a few adaptive sub-steps at    !
! dt_fast = 900 s, not the tiny steps an explicit solve of the unfrozen column would need. !
! The plant water mass adds no stiffness: its rates are the network's, constant over the   !
! step. A step the explicit march cannot resolve is redone on ARK by the dispatcher        !
! (rk45_state_railed, RK45_WORK_CAP).                                                      !
!                                                                                          !
! Shares build_column_frozen with the ARK path: the frozen record (leaf gas exchange,      !
! aerodynamics, the soil-plant water network's step, advective enthalpy) does not depend   !
! on which stepper advances the step, so it is reused verbatim.                            !
!==========================================================================================!
module meds_fast_rk45
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, cp_liq, cp_air, n_soil_layer_max
   use meds_therm_lib,           only : cas_temp_of_enthalpy, internal_energy_liquid, internal_energy_to_temp
   use meds_fast_time_derivs, only : surface_derivs, column_derivs, cas_conductances
   use meds_fast_types,       only : column_state_t, column_frozen_t, column_tend_t, error_control_t, &
                                     surface_state_t, surface_tend_t,          &
                                     column_config_t, column_cohort_t, column_forcing_t,         &
                                     column_budget_t, mask_is_full
   use meds_fast_frozen,      only : build_column_frozen
   use meds_column_state_ops, only : state_init, state_axpy, state_accum, state_sub, zero_like, apply_process_mask,     &
                                     clamp_cas, soil_water_store,                                   &
                                     soil_energy_store, plant_water_store, canopy_film_store,       &
                                     deposit_condensate, clamp_canopy_film, unpack_column_state,     &
                                     diagnose_soil_temps, transp_per_leaf_area
   use meds_site_diag_types,  only : CD_TRANSP
   use meds_fast_control,     only : state_wrms_grouped,                                         &
                                     step_control_factor
   use meds_config,           only : meds_config_t, CTRL_L2_STRICT
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_fast_types, only : patch_biophys_t, fast_workspace_t
   use meds_biophysics_opts, only : SOIL_BC_AQUIFER
   use meds_budget_check,     only : budget_check, last_check_closed, budget_energy_rate_floor,      &
                                     budget_water_rate_floor
   implicit none
   private

   public :: rk45_column_step, adaptive_rk45_march, column_fast_step_rk45
   public :: rk45_state_railed, rk45_ledgers_stop

   !----- Cash-Karp embedded 5(4) tableau (Cash & Karp 1990, ACM TOMS 16:201; the SAME                 !
   !      coefficients as Numerical Recipes' rkck). The stage times c_i are not needed: the            !
   !      frozen record has no explicit time dependence, so they never enter the code. ----------------!
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
   !----- Embedded-pair lower order, the p in step_control_factor's -1/(p+1) exponent                  !
   !      (docs/science/numerical_scheme.md eq. (4)). -------------------------------------------------!
   integer(ik), parameter :: RK45_P_ORDER = 4_ik

   !----- WORK BUDGET per dt_fast. The stability estimate in the module header puts normal operation     !
   !      at ~3 accepted sub-steps; a healthy step never approaches this. In its stiff (dense-canopy     !
   !      + cold) regime the explicit surface instead thrashes to hundreds of sub-steps and STILL        !
   !      rails -- wasted work, since the dispatcher then discards the step and redoes it on ARK.        !
   !      Bailing out at the cap turns that runaway into an early, cheap "this step is stiff"            !
   !      signal: it caps the wasted explicit work at ~20x normal instead of ~150x, the difference       !
   !      between a 30-yr run finishing and a single cold month taking longer than the rest of the       !
   !      run combined. ---------------------------------------------------------------------------------!
   integer(ik), parameter, public :: RK45_WORK_CAP = 64_ik

contains

   !---------------------------------------------------------------------------------------!
   ! stage_bnd -- the whole-column BOUNDARY-flux rates at one stage's state (the explicit              !
   ! analogue of column_be_stage's stage_bflux_t in meds_fast_be_stage): no backward-Euler             !
   ! solve, so the "state" IS the stage's input y_stage, and surf_tend is column_derivs' own           !
   ! surface diagnostic at y_stage. cond_enth is what surface_derivs debited from the CAS at           !
   ! THIS stage's own tcas, since condensation happens at the evolving CAS temperature. ---------------!
   pure subroutine stage_bnd(y_stage, frozen, surf_tend, rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i, &
                             atm_heat_i)
      type(column_state_t),  intent(in)  :: y_stage
      type(column_frozen_t), intent(in)  :: frozen
      type(surface_tend_t),  intent(in)  :: surf_tend
      real(wp),               intent(out) :: rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i
      real(wp),               intent(out) :: atm_heat_i   !< [W/m2] sensible export, the reported H
      real(wp) :: gah_i, gaw_i, gac_i
      rnet_i      = surf_tend%coh_rnet
      !----- The boundary flux must be charged at the conductance the TENDENCY used. column_derivs       !
      !      built this stage's CAS tendency from a live-state surface-layer re-solve, so reading        !
      !      the state^n frozen%cas%g_atm_heat here would book a boundary flux the state update          !
      !      never took -- the "borrow one solve's flux while committing another's state" defect         !
      !      class (docs/science/numerical_scheme.md §5b).                                               !
      !                                                                                                  !
      !      RE-SOLVED rather than read back off surf_tend: cas_conductances is `pure` and both call     !
      !      sites pass the same (frozen, y_stage) pair, so the recompute returns the same number,       !
      !      whereas adding fields to surface_tend_t perturbs ifx's inlining and FMA choices inside      !
      !      surface_derivs (~1e-12 in every result). ---------------------------------------------------!
      call cas_conductances(frozen%cas, y_stage%cas_enthalpy, y_stage%cas_shv, gah_i, gaw_i, gac_i)
      atm_enth_i  = gah_i * (y_stage%cas_enthalpy - frozen%cas%enthalpy_atm)
      atm_vap_i   = gaw_i * (y_stage%cas_shv      - frozen%cas%shv_atm)
      atm_heat_i  = gah_i * cp_air * (cas_temp_of_enthalpy(y_stage%cas_enthalpy, y_stage%cas_shv)          &
                                       - frozen%cas%mo_theta_atm)
      cond_i      = surf_tend%cond
      cond_enth_i = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
   end subroutine stage_bnd

   !---------------------------------------------------------------------------------------!
   ! rk45_column_step -- one Cash-Karp 5(4) step of size dt from y, over the pure RHS                   !
   ! column_derivs. Commits the 5TH-order solution (local extrapolation, as adaptive_ark_march          !
   ! commits its higher-order result); y_err is the (5th - 4th) embedded difference for the             !
   ! adaptive controller. w_out/e_in/e_out are the whole-column boundary-flux AMOUNTS over dt,          !
   ! b-weighted by the SAME 5th-order b-vector as the state commit (the consistent quadrature           !
   ! for a boundary integral over this step). They carry no infiltration term: infiltration             !
   ! moves water and enthalpy between two tracked stores, the pond and the soil, and the                !
   ! caller books the rainfall into the pond once per dt_fast.                                          !
   !                                                                                                    !
   ! STAGE CLAMPING (the explicit counterpart of ark2_column_step's base3 clamp_cas): every             !
   ! y_stage is explicit-only here (no ESDIRK stabilization), so a stiff leaf<->canopy-air coupling     !
   ! under a too-large trial dt can drive the canopy-air enthalpy or humidity far outside its           !
   ! physical domain WITHIN a single stage evaluation, where cas_temp_of_enthalpy diagnoses a wild      !
   ! temperature and qsat overflows rather than returning a large-but-finite value the controller      !
   ! could reject. clamp_cas bounds each y_stage's canopy air before its column_derivs call; the        !
   ! resulting k_i is then a legitimate derivative at a physical (if boundary-pinned) state, and the    !
   ! normal 5th/4th embedded-error comparison rejects and shrinks dt exactly as for any other           !
   ! oversized step -- no separate detection logic needed. The soil energy needs no such clamp: a       !
   ! stage reads it only through the temperature read-off, which any finite energy inverts with its     !
   ! liquid fraction in [0, 1], and through terms linear in that temperature (conduction, the ground's  !
   ! sensible heat, the enthalpy of its evaporation); the evaporation rate itself is frozen from the    !
   ! water network for the whole dt_fast. --------------------------------------------------------------!
   pure subroutine rk45_column_step(y, frozen, col_config, n, nsl, dt, y_out, y_err, w_out, e_in, e_out,          &
                                    clamp_stage_n, cond_out,                                         &
                                    tissue_leaf_int, tissue_wood_int, transp_int, cond_enth_out,       &
                                    atm_heat_out, atm_vap_out, face_resid_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      real(wp),               intent(out) :: w_out, e_in, e_out
      !----- STAGE-clamp telemetry (see column_budget_t%clamp_stage_n): a stage clamp bounds a       !
      !      throwaway input and is expected on an oversized trial. -----------------------------------!
      integer(ik), optional,  intent(inout) :: clamp_stage_n
      real(wp),    optional,  intent(out)   :: cond_out   !< [kg/m2] condensate to deposit into soil layer 1
      real(wp),    optional,  intent(out)   :: cond_enth_out !< [J/m2] its liquid enthalpy, b-weighted at the stage CAS temps
      real(wp), optional,     intent(out) :: atm_heat_out, atm_vap_out   !< [J/m2],[kg/m2] net CAS->atm sensible/vapour over dt
      real(wp),    optional,  intent(out)   :: tissue_leaf_int(n), tissue_wood_int(n)  !< [K*s]
      real(wp),    optional,  intent(out)   :: transp_int(n)   !< [kg/m2] each cohort's transpiration over dt
      !----- PER-LAYER FACE CLOSURE: sum over interior layers of |the mass this step's committed       !
      !      theta actually moved - the mass the SOIL-ENERGY equation was charged for|. Zero by        !
      !      construction while the energy advects on its own stage's water faces, which is the        !
      !      invariant; nonzero the moment a face is borrowed from elsewhere -- the defect class both  !
      !      whole-column ledgers are blind to, because a misplacement between LAYERS still sums       !
      !      correctly against the boundary (docs/science/numerical_scheme.md §5b). -------------------!
      real(wp),    optional,  intent(out)   :: face_resid_out   !< [kg/m2]

      type(column_tend_t)  :: k1, k2, k3, k4, k5, k6
      type(column_state_t) :: y_stage, y_4th
      type(surface_tend_t) :: surf_tend
      real(wp) :: rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6), atm_heat(6)
      real(wp) :: bw_rnet, bw_atm_enth, bw_atm_vap, bw_cond, bw_cond_enth, bw_drain, bw_atm_heat, bw_g_bottom
      !----- Per-stage tissue temperatures, kept so the store's energy can be set from their TIME    !
      !      INTEGRAL rather than from the last stage -- see column_bflux_t's note. RK45 is the       !
      !      accuracy baseline, so paying n x 6 here is the right trade. ------------------------------!
      real(wp) :: tleaf_s(n,6), twood_s(n,6), transp_s(n,6)
      real(wp) :: bw_face(n_soil_layer_max), f_in, f_out, uptk, face_resid
      integer(ik) :: kf

      call column_derivs(y, frozen, col_config, n, nsl, k1, sf_out=surf_tend)
      call stage_bnd(y, frozen, surf_tend, rnet(1), atm_enth(1), atm_vap(1), cond(1), cond_enth(1), atm_heat(1))
      tleaf_s(1:n,1) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,1) = surf_tend%wood_temp(1:n)
      transp_s(1:n,1) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage) ; call state_accum(y_stage, dt*A21, k1, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_config, n, nsl, k2, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(2), atm_enth(2), atm_vap(2), cond(2), cond_enth(2), atm_heat(2))
      tleaf_s(1:n,2) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,2) = surf_tend%wood_temp(1:n)
      transp_s(1:n,2) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A31, k1, n, nsl) ; call state_accum(y_stage, dt*A32, k2, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_config, n, nsl, k3, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(3), atm_enth(3), atm_vap(3), cond(3), cond_enth(3), atm_heat(3))
      tleaf_s(1:n,3) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,3) = surf_tend%wood_temp(1:n)
      transp_s(1:n,3) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A41, k1, n, nsl) ; call state_accum(y_stage, dt*A42, k2, n, nsl)
      call state_accum(y_stage, dt*A43, k3, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_config, n, nsl, k4, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(4), atm_enth(4), atm_vap(4), cond(4), cond_enth(4), atm_heat(4))
      tleaf_s(1:n,4) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,4) = surf_tend%wood_temp(1:n)
      transp_s(1:n,4) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A51, k1, n, nsl) ; call state_accum(y_stage, dt*A52, k2, n, nsl)
      call state_accum(y_stage, dt*A53, k3, n, nsl) ; call state_accum(y_stage, dt*A54, k4, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_config, n, nsl, k5, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(5), atm_enth(5), atm_vap(5), cond(5), cond_enth(5), atm_heat(5))
      tleaf_s(1:n,5) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,5) = surf_tend%wood_temp(1:n)
      transp_s(1:n,5) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A61, k1, n, nsl) ; call state_accum(y_stage, dt*A62, k2, n, nsl)
      call state_accum(y_stage, dt*A63, k3, n, nsl) ; call state_accum(y_stage, dt*A64, k4, n, nsl)
      call state_accum(y_stage, dt*A65, k5, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_config, n, nsl, k6, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6), atm_heat(6))
      tleaf_s(1:n,6) = surf_tend%leaf_temp(1:n) ; twood_s(1:n,6) = surf_tend%wood_temp(1:n)
      transp_s(1:n,6) = surf_tend%transp_c(1:n)

      !----- y_out = y + dt*(B1*k1 + B3*k3 + B4*k4 + B6*k6)  [5th order; b2=b5=0]. -----------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt*B1, k1, n, nsl)
      call state_accum(y_out, dt*B3, k3, n, nsl)
      call state_accum(y_out, dt*B4, k4, n, nsl)
      call state_accum(y_out, dt*B6, k6, n, nsl)
      !----- PER-LAYER FACE CLOSURE, measured on the b-weighted commit.                                 !
      !                                                                                                 !
      !      INTERIOR LAYERS ONLY (2 .. nsl-1). At layer 1 the water column's top face is q_top         !
      !      (infiltration MINUS soil evaporation) while the energy column's is infiltration alone      !
      !      -- deliberately different quantities, since the evaporated water's enthalpy is paid at     !
      !      the surface, not advected across the soil top face. Layer nsl is excluded for the          !
      !      mirror reason at the drainage face. Both boundaries are already covered by the             !
      !      whole-column ledgers; the interior is precisely what those ledgers cannot see. ------------!
      if (present(face_resid_out)) then
         do kf = 1_ik, nsl
            bw_face(kf) = B1*k1%soil_face(kf) + B3*k3%soil_face(kf)                                &
                        + B4*k4%soil_face(kf) + B6*k6%soil_face(kf)
         end do
         face_resid = 0.0_wp
         do kf = 2_ik, nsl - 1_ik
            f_in  = bw_face(kf-1_ik) * rho_h2o * dt
            f_out = bw_face(kf)      * rho_h2o * dt
            uptk  = frozen%roots%uptake * frozen%roots%root_share(kf) * dt
            face_resid = face_resid                                                                &
               + abs((y_out%theta(kf) - y%theta(kf)) * col_config%soil%dz(kf) * rho_h2o         &
                     - (f_in - f_out - uptk))
         end do
         face_resid_out = face_resid
      end if

      !----- The COMMITTED state is not clamped. The accept test is `err <= 1 .or. dt <= dt_floor`,       !
      !      so at the sub-step floor a step is committed with the controller out of moves, and a         !
      !      clamp on the kept state would keep whatever mass or energy it moved with no ledger           !
      !      entry. The STAGE clamps above bound throwaway inputs so column_derivs stays evaluable,       !
      !      and a discarded stage state breaks no books. An out-of-range commit is handled like          !
      !      every other bad step: the embedded error grows and the controller rejects and shrinks;       !
      !      a floor-forced step that still commits a railed state is caught by rk45_state_railed at      !
      !      the dispatch, and the rescue redoes the whole dt_fast on ARK. -------------------------------!

      !----- y_4th = y + dt*(BS1*k1 + BS3*k3 + BS4*k4 + BS5*k5 + BS6*k6)  [embedded 4th order].         !
      !      y_4th is a SEPARATE named temporary, not y_err itself: passing y_err as both               !
      !      state_sub's `b` and its `out` aliases two dummies, one of them modified, which Fortran     !
      !      forbids -- the compiler may clobber b (incl. deallocating its allocatable components)      !
      !      before the subtraction reads it. The aliased form gives a y_err that does not shrink       !
      !      with dt at all, instead of the expected O(dt^5) embedded-error scaling. -------------------!
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
      !----- DRAINAGE b-weighted from the stage tendencies, with the same b-vector as the state commit, !
      !      so the ledger books the water the committed theta shed through the bottom face. Every      !
      !      stage carries the water network's drainage, so this is that rate times dt. ---------------!
      bw_drain     = B1*k1%drainage_rate + B3*k3%drainage_rate + B4*k4%drainage_rate                  &
                     + B6*k6%drainage_rate
      !----- The heat the soil tendency moved across the column's bottom face, b-weighted like the     !
      !      drainage: zero on the adiabatic base, the conduction to deep_temp under the Dirichlet     !
      !      anchor, which the soil store carries and so the ledger must too. -------------------------!
      bw_g_bottom  = B1*k1%g_bottom + B3*k3%g_bottom + B4*k4%g_bottom + B6*k6%g_bottom
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
      !----- ...and each cohort's transpiration, the same b-weighted stage fluxes its leaf water was    !
      !      debited and the canopy air credited. -----------------------------------------------------!
      if (present(transp_int)) transp_int(1:n) = dt * (B1*transp_s(1:n,1) + B3*transp_s(1:n,3)       &
                                                     + B4*transp_s(1:n,4) + B6*transp_s(1:n,6))

      !----- No clip or floor term: the water network keeps every layer inside its curve, so             !
      !      nothing corrects theta after the stages (docs/science/numerical_scheme.md §5b,              !
      !      "Saturation is the water network's, and nothing is clipped"). ------------------------------!
      !----- ground_rad is the snowfac-BLENDED radiative input (= abs_sw+abs_lw when bare). -------------!
      !----- The infiltration enthalpy is a POND -> SOIL transfer between two tracked stores, so it      !
      !      is not a boundary input here; the caller's pond store carries the other half. Booking       !
      !      it here would add exactly the infiltration enthalpy as a spurious surplus. -----------------!
      e_in  = (bw_rnet + frozen%snow%ground_rad + bw_g_bottom) * dt
      e_out = bw_atm_enth * dt                                                                    &
              + bw_drain * dt * internal_energy_liquid(frozen%hydrology%t_bot)
      !----- cond is EXCLUDED from w_out: the caller deposits it into soil layer 1, so it is not       !
      !      lost across the boundary; cond_out returns the amount for that deposit. Runoff is not     !
      !      booked here either: the caller adds the network's runoff once per dt_fast. Drainage is    !
      !      RK45's own (b-weighted above). -----------------------------------------------------------!
      w_out = bw_atm_vap * dt + bw_drain * dt
      if (present(cond_out))      cond_out      = bw_cond * dt
      if (present(cond_enth_out)) cond_enth_out = bw_cond_enth * dt
      if (present(atm_heat_out))  atm_heat_out  = bw_atm_heat * dt
      if (present(atm_vap_out))   atm_vap_out   = bw_atm_vap  * dt
   end subroutine rk45_column_step

   !---------------------------------------------------------------------------------------!
   ! adaptive_rk45_march -- integrate to t_end with the Cash-Karp embedded error driving the           !
   ! step controller, with the same grouped WRMS, I/PI controller and warm start as the ARK            !
   ! (docs/science/numerical_scheme.md eqs. (3)-(4), p = 4) and adaptive_ark_march's                   !
   ! accept/reject structure, over rk45_column_step. Every state RK45 integrates is scored except      !
   ! the canopy films, which carry no group; soil water and plant water are scored too, though with    !
   ! their rates constant over the step their embedded difference is round-off. -----------------------!
   subroutine adaptive_rk45_march(y0, frozen, col_config, n, nsl, t_end, ec, dt_init, y_out, nsteps, nrej,       &
                                  w_out_acc, e_in_acc, e_out_acc, dt_warm_out,                    &
                                  clamp_stage_n, cond_acc,                                            &
                                  cond_enth_acc, atm_heat_acc, atm_vap_acc, face_resid_acc,           &
                                  ood_max, tissue_leaf_acc, tissue_wood_acc, transp_acc)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec
      type(column_state_t),  intent(out) :: y_out
      integer(ik),            intent(out) :: nsteps, nrej
      real(wp),               intent(out) :: w_out_acc, e_in_acc, e_out_acc
      real(wp),    optional,  intent(out) :: dt_warm_out
      !----- STAGE-clamp telemetry over the whole march, counted on every trial: a rejected trial's   !
      !      clamp still says the controller probed too far, which is the signal. ----------------------!
      integer(ik), optional,  intent(inout) :: clamp_stage_n
      !----- Per-layer face closure, accumulated over ACCEPTED sub-steps only -- a rejected             !
      !      trial's state is thrown away, so its residual is meaningless. -----------------------------!
      real(wp),    optional,  intent(inout) :: face_resid_acc
      real(wp),    optional,  intent(inout) :: cond_acc   !< [kg/m2] accumulated condensate
      real(wp),    optional,  intent(inout) :: cond_enth_acc !< [J/m2] ...and its stage-valued liquid enthalpy
      real(wp), optional,     intent(inout) :: atm_heat_acc, atm_vap_acc   !< [J/m2],[kg/m2] net CAS->atm export, accepted sub-steps
      !----- [K*s] per-cohort tissue-temperature time integrals over the accepted march. Caller zeroes.!
      real(wp),    optional,  intent(inout) :: tissue_leaf_acc(n), tissue_wood_acc(n)
      real(wp),    optional,  intent(inout) :: transp_acc(n)   !< [kg/m2] per-cohort transpiration, accepted march
      !----- [m3/m3] running MAX theta excursion outside [theta_res, theta_sat] at a stage-1 RHS         !
      !      input (telemetry: see column_budget_t%theta_ood_max). Measured on EVERY trial,              !
      !      accepted or not: a rejected trial's out-of-domain probe is still a place the RHS was        !
      !      evaluated, which is exactly what this reports. ---------------------------------------------!
      real(wp),    optional,  intent(inout) :: ood_max

      type(column_state_t) :: y, y_new, y_err, y_zero
      real(wp) :: t, dt, err, err_prev, fac, dt_floor
      real(wp) :: w_out, e_in, e_out, dt_try, dt_warm
      real(wp) :: cond_i, cond_enth_i, face_i
      real(wp)    :: atm_heat_i, atm_vap_i
      real(wp) :: tl_int_i(n), tw_int_i(n), tr_int_i(n)
      integer(ik) :: kood
      logical  :: clamped

      w_out_acc = 0.0_wp ; e_in_acc = 0.0_wp ; e_out_acc = 0.0_wp
      !----- Sub-step FLOOR: t_end/64, as in adaptive_ark_march. A floor step is not guaranteed        !
      !      stable here (RK45 is explicit), but the module header's stability estimate bounds         !
      !      normal operation to ~3 sub-steps, so a stiffness-driven runaway signals a genuinely       !
      !      pathological step; degrade gracefully rather than grind at a tiny floor. -----------------!
      dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)   ! clamp-ok: the smallest sub-step [s] (above)

      call state_init(y0, n, nsl, y)
      t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik
      err_prev = -1.0_wp
      dt_warm = dt
      do
         if (t >= t_end - tiny_num) exit
         dt_try = dt
         dt = min(dt, t_end - t)
         clamped = dt < dt_try - tiny_num
         !----- theta is never clamped: record how far the step's starting theta lies outside the       !
         !      constitutive domain (the water network allows a little storage above saturation). -----!
         if (present(ood_max)) then
            do kood = 1_ik, nsl
               ood_max = max(ood_max, y%theta(kood) - col_config%soil%theta_sat(kood),                       &
                             col_config%soil%theta_res(kood) - y%theta(kood))
            end do
         end if
         call rk45_column_step(y, frozen, col_config, n, nsl, dt, y_new, y_err, w_out, e_in, e_out,               &
                               clamp_stage_n=clamp_stage_n, cond_out=cond_i,                      &
                               tissue_leaf_int=tl_int_i, tissue_wood_int=tw_int_i, transp_int=tr_int_i, &
                               cond_enth_out=cond_enth_i, atm_heat_out=atm_heat_i, atm_vap_out=atm_vap_i, &
                               face_resid_out=face_i)
         !----- named temporary: never pass a derived-type-valued function result straight into a  !
         !      call (the nvfortran whole-program-optimizer trap documented in CLAUDE.md). --------!
         y_zero = zero_like(y_err, n, nsl)
         err = state_wrms_grouped(y_err, y_zero, y_new, n, nsl, ec%tols)
         if (err /= err .or. dt /= dt) then    ! non-finite: too big a step, reject deterministically
            if (dt <= dt_floor) then
               !----- floor-forced commit of a non-finite trial: tallied like any accept. ------------!
               if (present(face_resid_acc)) face_resid_acc = face_resid_acc + face_i
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
            if (present(face_resid_acc)) face_resid_acc = face_resid_acc + face_i
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
            if (present(transp_acc)) transp_acc(1:n) = transp_acc(1:n) + tr_int_i(1:n)
            t = t + dt ; nsteps = nsteps + 1_ik
            err_prev = err
            if (.not. clamped) dt_warm = dt
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         !----- WORK-BUDGET bail-out: stop burning explicit sub-steps on a step the dispatcher is       !
         !      going to discard anyway. y_out is left at the partial-time state, which is FINE --      !
         !      the caller returns immediately without committing it (stiff_bail). ---------------------!
         if (nsteps + nrej >= RK45_WORK_CAP) exit
      end do
      call state_init(y, n, nsl, y_out)
      if (present(dt_warm_out)) dt_warm_out = dt_warm
   end subroutine adaptive_rk45_march


   !=======================================================================================!
   !  INTEG_RK45 path: the ED2-style adaptive Cash-Karp fast step. Shares build_column_frozen       !
   !  with the ARK path (meds_fast_ark): the frozen record does not depend on which stepper         !
   !  advances the step. ---------------------------------------------------------------------------!
   !=======================================================================================!
   subroutine column_fast_step_rk45(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,  &
                                    ws, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, &
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
      type(fast_workspace_t),  intent(inout) :: ws       !< this thread's reusable step storage
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters
      !----- .true. when the explicit march hit its work budget (RK45_WORK_CAP) without resolving        !
      !      the step -- the stiff-regime signal. biophys is then left UNTOUCHED at state^n and the      !
      !      caller redoes this dt_fast on ARK. ---------------------------------------------------------!
      logical,     optional,   intent(out)   :: stiff_bail
      real(wp),    optional,   intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) per-cohort diagnostic capture

      type(column_frozen_t)  :: frozen
      type(column_state_t)   :: y, y_out
      type(surface_state_t)  :: y_stage
      type(surface_tend_t)   :: surf_tend
      real(wp)    :: dt0, cas_mass_capacity, enth0, shv0, enth1, shv1
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, w_plant0, w_plant1, w_surface0
      real(wp)    :: w_out_acc, e_in_acc, e_out_acc, w_in, w_out, e_in, e_out
      real(wp)    :: tg, fl, dt_warm_next, cond_dep, cond_dep_enth, face_acc
      real(wp)    :: w_pond_rk, runoff_rk
      real(wp)    :: e_pond0, e_pond_rk, over_enth_rk   ! [J/m2] pond enthalpy (start, end); runoff enthalpy
      !----- Canopy-SURFACE water ledger scratch. -----------------------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, leaf_cap_i, wood_cap_i, intercept_total
      real(wp)    :: tissue_store0, tissue_store1
      real(wp)    :: tl_int_acc(col_cohort%n), tw_int_acc(col_cohort%n), tr_int_acc(col_cohort%n)
      real(wp)    :: cap_leaf_a(col_cohort%n), cap_wood_a(col_cohort%n)
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, k, i, nsteps, nrej

      n = col_cohort%n ; nsl = col_config%soil%n_active
      !----- The stage-clamp counter accumulates down the call chain, so this sub-step's tally starts  !
      !      clean. ------------------------------------------------------------------------------------!
      budget%clamp_stage_n = 0_ik ; budget%theta_ood_max = 0.0_wp
      !----- state^n ponding store, captured BEFORE the unpack below overwrites it. ------------------!
      w_surface0 = biophys%soil_w%w_surface
      e_pond0    = biophys%soil_w%w_surface_enth
      if (present(stiff_bail)) stiff_bail = .false.

      cond_dep = 0.0_wp ; cond_dep_enth = 0.0_wp ; face_acc = 0.0_wp
      budget%atm_heat_export = 0.0_wp ; budget%atm_vap_export = 0.0_wp
      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws%spac, &
                               n, nsl, frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      dt0 = dt_fast
      if (biophys%adapt_dt_last > tiny_num) dt0 = min(biophys%adapt_dt_last, dt_fast)
      if (col_config%integrator%dt_init > tiny_num) dt0 = min(col_config%integrator%dt_init, dt_fast)
      ec = col_config%integrator%error_control
      ec%p_order = RK45_P_ORDER
      tl_int_acc(1:n) = 0.0_wp ; tw_int_acc(1:n) = 0.0_wp ; tr_int_acc(1:n) = 0.0_wp
      call adaptive_rk45_march(y, frozen, col_config, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                              w_out_acc, e_in_acc, e_out_acc, dt_warm_out=dt_warm_next,           &
                              clamp_stage_n=budget%clamp_stage_n,                                  &
                              cond_acc=cond_dep, cond_enth_acc=cond_dep_enth,                     &
                              atm_heat_acc=budget%atm_heat_export, atm_vap_acc=budget%atm_vap_export, &
                              ood_max=budget%theta_ood_max, face_resid_acc=face_acc,                 &
                              tissue_leaf_acc=tl_int_acc, tissue_wood_acc=tw_int_acc, transp_acc=tr_int_acc)
      biophys%adapt_dt_last = dt_warm_next
      budget%integ_nsteps = nsteps ; budget%integ_nrej = nrej
      !----- PER-LAYER FACE CLOSURE: the theta the stages committed against the face flows and        !
      !      root uptake they were given, layer by layer. --------------------------------------------!
      budget%soil_face_mass%resid   = face_acc
      budget%soil_face_mass%worst   = max(budget%soil_face_mass%worst, abs(face_acc))
      budget%soil_face_mass%abs_sum = budget%soil_face_mass%abs_sum + abs(face_acc)
      budget%soil_face_mass%n_check = budget%soil_face_mass%n_check + 1_ik

      !----- WORK-BUDGET bail: the march gave up (stiff regime). Return NOW, before committing           !
      !      anything into biophys or running the budget ledgers -- y_out is a partial-time state,       !
      !      so committing it would corrupt the column and trip the conservation checks. biophys         !
      !      is untouched at state^n, so the caller can cleanly redo this dt_fast on ARK. ---------------!
      if (nsteps + nrej >= RK45_WORK_CAP) then
         if (present(stiff_bail)) stiff_bail = .true.
         if (present(converged))  converged  = .false.
         if (present(iters))      iters      = nsteps
         return
      end if

      !----- PROCESS MASK (process_mask_t): a masked-off component is restored to state^n (y).         !
      !      Soil water and plant water ride the stages at the network's rates, so both honour the     !
      !      mask here. -------------------------------------------------------------------------------!
      call apply_process_mask(col_config%mask, y, y_out, n, nsl)

      !----- Canopy-SURFACE water: capacity clamp and overflow bookkeeping, at the store's boundary        !
      !      after the march; the ODE itself is left unclamped mid-integration, since the                  !
      !      closure-safe clamp point is a store's BOUNDARY, not a flux feeding an already                 !
      !      self-consistent solve. Gated behind canopy_water_on, so with the store off no film            !
      !      arithmetic runs. DEFICIT (the floor, symmetric to overflow): film_evap uses a                 !
      !      state^n-frozen conductance (rescaled by availability in build_column_frozen, but only         !
      !      approximately -- the evaporative demand can grow through the step as tcas and                 !
      !      leaf_temp evolve), so the store can transiently overdraw below 0. Left negative, it           !
      !      would corrupt the NEXT dt_fast's frozen interception rate (intercept_canopy_layer's           !
      !      own floor on a negative starting bucket creates water). Floor at 0 and book the               !
      !      shortfall as a NEGATIVE addition to w_out/e_out (flooring up to 0 makes the store             !
      !      appear to gain, so the outflow must shrink by the same amount) -- the mirror of               !
      !      surf_overflow's sign. ------------------------------------------------------------------------!
      surf_overflow = 0.0_wp ; surf_deficit = 0.0_wp
      if (col_config%canopy_water_on) then
         call clamp_canopy_film(y_out, col_cohort%lai, col_cohort%wai,                                  &
                                col_config%leaf_surf_water_max(col_cohort%pft(1:n)),                    &
                                col_config%wood_surf_water_max(col_cohort%pft(1:n)), n, surf_overflow, surf_deficit)
      end if

      !----- unpack into biophys + re-derive the diagnostic soil/leaf/wood temperatures. -----------!
      call unpack_column_state(y_out, n, nsl, biophys)
      !----- DEPOSIT THE CONDENSATE (the same routing as column_fast_step_ark, which has the full        !
      !      rationale): dew and fog land on a surface inside the column, so the condensate goes,        !
      !      mass and enthalpy paired, into soil layer 1 at the CAS temperature, and the                 !
      !      whole-column ledger closes with no boundary term. Both schemes use the same                 !
      !      destination. -------------------------------------------------------------------------------!
      if (cond_dep > 0.0_wp) then
         !----- at the b-weighted STAGE enthalpy the CAS was debited, not u_liq(T_end). ---------------!
         call deposit_condensate(y_out, col_config%soil%dz(1), cond_dep, cond_dep_enth)
         biophys%soil_e%soil_energy(1) = y_out%soil_energy(1) ; biophys%soil_w%theta(1) = y_out%theta(1)
      end if
      call diagnose_soil_temps(y_out, col_config%soil_thermal%soil_dry_heat_capacity, nsl, biophys%soil_e%soil_temp, &
            biophys%soil_e%soil_fliq)
      !----- The pond is the water network's, as theta is (docs/science/soil_plant_water.md): the         !
      !      stages moved theta at the network's steady rate, so the committed theta is the               !
      !      network's, and the pond, its runoff and their enthalpy are the network's too. The            !
      !      network keeps every layer inside its curve (above residual; above saturation only            !
      !      through its specific storage), so the commit needs no guards. Under the soil_water           !
      !      mask the pond keeps its state^n value, like theta. ------------------------------------------!
      if (col_config%mask%soil_water) then
         w_pond_rk    = frozen%hydrology%w_surface1
         e_pond_rk    = frozen%hydrology%w_surface_enth1
         runoff_rk    = frozen%hydrology%runoff_surf * dt_fast
         over_enth_rk = frozen%hydrology%runoff_enth * dt_fast
         biophys%soil_w%w_surface      = w_pond_rk
         biophys%soil_w%w_surface_enth = e_pond_rk
      else
         w_pond_rk = w_surface0 ; e_pond_rk = e_pond0 ; runoff_rk = 0.0_wp ; over_enth_rk = 0.0_wp
      end if

      !----- Commit the tissue temperatures from the b-weighted stage time integrals (below); the      !
      !      frozen store relaxed them from frozen%tissue%t_leaf0/t_wood0 against a = cap/dt inside every    !
      !      stage, so nothing is corrected afterwards. ----------------------------------------------!
      tissue_store0 = 0.0_wp ; tissue_store1 = 0.0_wp
      do i = 1_ik, n
         !----- Capacity from the SAME *_hcap_per_dt the kernel relaxed against, as on the ARK           !
         !      path (column_fast_step_ark says why). ---------------------------------------------------!
         cap_leaf_a(i) = frozen%tissue%leaf_hcap_per_dt(i) * dt_fast
         cap_wood_a(i) = frozen%tissue%wood_hcap_per_dt(i) * dt_fast
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
      !----- The transpiration each cohort's leaves lost over the step, as on the ARK path. ---------!
      if (present(cdiag)) cdiag(CD_TRANSP, 1:n) = transp_per_leaf_area(tr_int_acc(1:n), dt_fast,          &
                                                                       col_cohort%nplant(1:n), col_cohort%leaf_area(1:n))
      do i = 1_ik, n
         tissue_store1 = tissue_store1 + cap_leaf_a(i) * biophys%leaf_temp(i)                             &
                                       + cap_wood_a(i) * biophys%wood_temp(i)
      end do
      !----- ...plus what the tissue WATER stores gained as their water mass changed, valued at each    !
      !      tissue's start-of-step temperature (build_column_frozen). The tissue temperatures never saw  !
      !      that part, so cap*T alone would leave it unbooked. -------------------------------------------!
      tissue_store1 = tissue_store1 + frozen%tissue%water_store_enth * dt_fast

      !----- WHOLE-COLUMN CONSERVATION LEDGER: unlike ARK's per-kernel + whole ledgers, RK45 checks     !
      !      the WHOLE-COLUMN water and energy budgets only -- every store advances through the         !
      !      SAME column_derivs RHS, so the whole-column ledger is the individual-store ledger here. ---!
      cas_mass_capacity = frozen%cas%cas_mass_capacity
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv
      e_soil0 = soil_energy_store(y%soil_energy,     col_config%soil%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, col_config%soil%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     col_config%soil%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, col_config%soil%dz, nsl)
      w_plant0 = plant_water_store(col_cohort%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(col_cohort%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- Canopy-SURFACE water, already per m2 of ground (no nplant factor, unlike w_plant0/1).         !
      !      Valued at frozen%film%film_liquid_enthalpy = u_liq(t_film_valuation); the tissue pays         !
      !      enthalpy_vapor - film_liquid_enthalpy per kg of film it evaporates (surface_derivs),          !
      !      so the store closes exactly against the CAS credit. All zero when canopy_water_on is          !
      !      off. -----------------------------------------------------------------------------------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))

      !----- e_in: the b-weighted stage inputs (e_in_acc: net radiation, ground radiation, the             !
      !      bottom face), the intercepted rain at the film's valuation, the rainfall enthalpy the         !
      !      pack took (acc_enth), and the water reaching the pond other than meltwater, at the            !
      !      temperature that values it (t_pond_inflow; the meltwater's enthalpy is already inside         !
      !      the pack/pond pair). Infiltration is a pond -> soil transfer, not an input. w_out_acc,        !
      !      e_*_acc and w_in are AMOUNTS over the whole dt_fast, so the rate intercept_total takes        !
      !      *dt_fast, while surf_overflow/surf_deficit are amounts already. surf_deficit SUBTRACTS        !
      !      (the mirror of surf_overflow's sign: flooring a negative store up to 0 makes it appear        !
      !      to gain, so the outflow shrinks by the same amount). -----------------------------------------!
      !----- w_in: rainfall, snowfall (into the pack, or below the new-snow threshold into the             !
      !      throughfall) and the shed water, all boundary water inputs. The shed water's enthalpy         !
      !      rides the precip_ground term of e_in (build_column_frozen mixes it into                       !
      !      hforc%precip_ground). ------------------------------------------------------------------------!
      w_in  = (forc%rainfall + forc%snowfall + biophys%shed_water_rate) * dt_fast
      !----- The pond's overflow over the step, from the water network. ------------------------------!
      w_out = w_out_acc + surf_overflow - surf_deficit + runoff_rk
      e_in  = e_in_acc + intercept_total * dt_fast * internal_energy_liquid(frozen%hydrology%t_film_valuation)        &
              + frozen%snow%acc_enth                                                             &
              + (frozen%hydrology%precip_ground - frozen%snow%melt_rate) * dt_fast                             &
                * internal_energy_liquid(frozen%hydrology%t_pond_inflow)
      !----- What leaves is the pond OVERFLOW, at the pond's own temperature                            !
      !      (docs/science/numerical_scheme.md §5b, "The pond holds heat"). ----------------------------!
      e_out = e_out_acc + (surf_overflow - surf_deficit) * internal_energy_liquid(frozen%hydrology%t_film_valuation)     &
              + over_enth_rk

      !----- FLUX-scaled tolerances (meds_budget_check header), same rule as the ARK ledgers. The     !
      !      checks never stop the run here: column_fast_step may still discard this step for the ARK !
      !      rescue, and calls rk45_ledgers_stop only for a step it keeps. ------------------------!
      call budget_check(budget%whole_water,                                                            &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0                     &
                        + frozen%snow%swe0,                                                          &
                        w_soil1 + cas_mass_capacity*shv1 + w_pond_rk + w_plant1 + surf_water1                        &
                        + frozen%snow%swe1,                                                          &
                        w_in, w_out, dt_fast, budget_water_rate_floor, 'whole_water (rk45)', .false.)
      !----- The snow store and its accumulated rainfall enthalpy join the ledger; 0 without snow. ----!
      call budget_check(budget%whole_energy,                                                           &
                        !----- The pack hands its meltwater, mass and enthalpy, to the POND: both     !
                        !      are tracked stores, so the pair telescopes here (as on ARK). ----------!
                        e_soil0                           + cas_mass_capacity*enth0 + surf_enth0                    &
                        + frozen%snow%enth0 + e_pond0 + tissue_store0,                              &
                        e_soil1 + cas_mass_capacity*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond_rk            &
                        + tissue_store1,                                                                &
                        e_in, e_out, dt_fast, budget_energy_rate_floor, 'whole_energy (rk45)', .false.)
      !----- CAS CO2 is not checked on this path: a per-kernel cas_co2 closure (ARK's budget%cas_co2)    !
      !      needs a b-weighted per-stage CO2 atmospheric-exchange accumulation, which                   !
      !      rk45_column_step does not keep (it tracks rnet, atm_enth, atm_vap and cond only). ----------!

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_rk45

   !----- Under [fast].debug_error, stop the run if a step column_fast_step KEEPS breached a        !
   !      whole-column ledger. column_fast_step_rk45 records the checks without stopping, because   !
   !      a step handed to the ARK rescue is discarded and must not end the run. -------------------!
   subroutine rk45_ledgers_stop(budget, dt_fast, col_config)
      type(column_budget_t), intent(in) :: budget
      real(wp),              intent(in) :: dt_fast
      type(column_config_t), intent(in) :: col_config
      if (.not. (col_config%energy%debug_error .and. mask_is_full(col_config%mask))) return
      if (.not. last_check_closed(budget%whole_water, dt_fast, budget_water_rate_floor)) then
         write(*, '(a,es13.5)') 'meds_fast_rk45: whole_water (rk45) did not close, resid = ',               &
                                budget%whole_water%resid
         error stop 'meds_fast_rk45: a kept RK45 step breached its whole-column water ledger'
      end if
      if (.not. last_check_closed(budget%whole_energy, dt_fast, budget_energy_rate_floor)) then
         write(*, '(a,es13.5)') 'meds_fast_rk45: whole_energy (rk45) did not close, resid = ',              &
                                budget%whole_energy%resid
         error stop 'meds_fast_rk45: a kept RK45 step breached its whole-column energy ledger'
      end if
   end subroutine rk45_ledgers_stop

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
