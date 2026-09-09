!==========================================================================================!
! meds_fast_ark -- the production fast-loop integrator: an L-stable ESDIRK2 (ARS(2,2,2) tableau;   !
! the explicit part is empty, so despite the historical "IMEX-ARK" name this is a diagonally       !
! implicit scheme). Design: docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md.                                !
!                                                                                          !
! Hosts the TABLEAU AND THE MARCH, and nothing else: the scheme entry `column_fast_step_ark`   !
! (called from `meds_fast_step`), the two-stage step `ark2_column_step`, and its embedded-error !
! controller `adaptive_ark_march`.                                                              !
!                                                                                          !
! What used to live here and no longer does, because none of it is ARK's:                     !
!   * the frozen work record builder      -> meds_fast_frozen    (RK45 builds the same record) !
!   * the implicit stage + its Newton     -> meds_fast_be_stage  (the RK4 oracle uses it too)  !
!   * the water-mass and canopy-film advance -> meds_fast_be_stage                             !
!   * the column_state_t algebra + ledgers   -> meds_column_state_ops (moved by PR #120)       !
! Importing those from a module named for one scheme made the two schemes look related when    !
! they are only siblings. Nothing outside this module imports `meds_fast_ark` now except the    !
! dispatcher and the RHS test.                                                                  !
!                                                                                          !
! Soil water is committed ONCE per dt_fast from the scratch advance_soil_water_column solve (the ARK   !
! stages pass theta through); the pond is carried in column_state_t but committed the same way.    !
! Whole-column water and energy close to round-off on every bottom BC (free-drain, bedrock,        !
! aquifer); see meds_budget_check and the ledgers at the end of column_fast_step_ark.              !
!==========================================================================================!
module meds_fast_ark
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o
   use meds_hydr_lib, only : water_content
   use meds_config, only : meds_config_t, INTEG_ARK, CTRL_L2_STRICT
   use meds_fast_control, only : state_wrms_grouped, step_control_factor
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_soil_types, only : snow_melt_t
   use meds_fast_types, only : patch_biophys_t
   use meds_fast_time_derivs, only : surface_derivs
   use meds_column_state_ops, only : state_init, state_extrap, state_err_diff, state_sub, bflux_zero, bflux_add, &
                                     bflux_bweight, clamp_cas, clamp_theta, clamp_soil_energy, soil_water_store, &
                                     soil_energy_store, plant_water_store, canopy_film_store, deposit_condensate, &
                                     clamp_canopy_film, unpack_column_state, diagnose_soil_temps, &
                                     assemble_soil_energy_forcing, apply_process_mask
   use meds_fast_be_stage, only : column_be_stage, advance_water_mass_full, advance_surf_water_full
   use meds_fast_frozen, only : build_column_frozen
   use meds_fast_types, only : column_config_t, column_cohort_t, column_forcing_t, column_budget_t, alloc_column_cohort, &
                               column_state_t, column_frozen_t, surface_state_t, cas_boundary_t, surface_tend_t, &
                               stage_bflux_t, column_bflux_t, error_control_t, column_tend_t, mask_is_full
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction, ground_surface_fluxes
   use meds_therm_lib, only : internal_energy_liquid, internal_energy_to_temp, internal_energy_ice, temp_of_liquid_enthalpy
   use meds_budget_check, only : budget_check, budget_energy_rate_floor, budget_water_rate_floor, budget_co2_rate_floor
   implicit none
   private

   !----- `niter` on column_be_stage / ark2_column_step / adaptive_ark_march is a two-valued port:   !
   !      <= 1 selects the uncoupled single-BE pass, anything > 1 selects the coupled 2x2 Newton      !
   !      (whose own iteration cap is NEWT_MAX, not this number). Passing a named sentinel rather      !
   !      than a bare 8 stops it reading as a tunable count it never was -- see meds_config's          !
   !      ark_coupled. -------------------------------------------------------------------------------!
   integer(ik), parameter :: NEWT_COUPLED = 2_ik

   public :: column_fast_step_ark
   public :: ark2_column_step, adaptive_ark_march

contains



   !---------------------------------------------------------------------------------------!
   ! ark2_column_step -- one 2nd-order L-stable IMEX step via the ARS(2,2,2) additive Runge-Kutta      !
   ! (Ascher-Ruuth-Spiteri 1997, Appl.Numer.Math. 25:151; identical gamma in Giraldo et al. 2013     !
   ! "ARK2"). Stiffly accurate on both tableaux, so the two ESDIRK stages map onto imex_euler's       !
   ! "last BE solve = committed state" structure with dt -> gamma*dt. The biotic CO2 source is folded !
   ! IMPLICIT (stays in the CO2 BE numerator), so f_E == 0 and the scheme reduces to a clean 2-solve  !
   ! ESDIRK2. PLANT WATER MASS IS OPERATOR-SPLIT OUT of the tableau, exactly like psi was before it     !
   ! (MEDS_ED2_RK45_DESIGN.md sec 1/4/5): mass is frozen through the stages, then advanced once over    !
   ! the full dt (advance_water_mass_full -- now a trivial closed-form Euler step, no iteration, since   !
   ! the frozen sapflow/uptake pre-pass already absorbed the only stiff physics), and excluded from the  !
   ! embedded error. y_err = (Y3-base3)-(Y2-y_n) is the free embedded 1st-order estimate for the         !
   ! adaptive controller (2 solves/step vs step-doubling's 3). Hydraulics WORK counters (section 5.3)     !
   ! now come from the Act-1 pre-pass's solve_plant_water_batch call (build_column_frozen), not from       !
   ! this per-stage endpoint update -- there is no more per-stage hydraulics solve to count. -------------!
   subroutine ark2_column_step(y, frozen, n, nsl, dt, y_out, y_err, niter, bf, clamp_n)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      integer(ik), optional, intent(in)  :: niter
      type(column_bflux_t), optional, intent(out) :: bf   !< b-weighted boundary-flux AMOUNTS over dt (ledger)
      !----- STAGE-clamp activations (see column_budget_t%clamp_stage_n). The ARK clamps its ARS       !
      !      stage-3 extrapolation base and NOTHING else -- y_out is committed unclamped -- so this     !
      !      path contributes to the stage counter only, and its commit counter stays 0 by             !
      !      construction. That asymmetry against RK45 is a result, not an omission. ------------------!
      integer(ik), optional, intent(inout) :: clamp_n
      real(wp), parameter :: GAMMA = 0.2928932188134524_wp   ! 1 - 1/sqrt(2)
      real(wp), parameter :: BETA  = 2.4142135623730951_wp   ! (1-gamma)/gamma = 1 + sqrt(2)
      type(column_state_t)  :: Y2, base3, Y3
      type(stage_bflux_t)   :: bf2, bf3
      type(surface_tend_t)  :: sf2, sf3
      real(wp)              :: transp_bw(n)
      real(wp)              :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      integer(ik) :: np
      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      !----- Stage 2: gamma*dt BE stage from y_n (CAS+soil only; mass frozen -- it is split out). -----!
      call column_be_stage(y, frozen, n, nsl, GAMMA*dt, Y2, niter=np, bf=bf2, sf_out=sf2)
      !----- Stage 3: extrapolated base. The BETA=2.414 extrapolation can overshoot BOTH the vG theta   !
      !      range AND the CAS enthalpy into a wild temperature where qsat(T) overflows to NaN; clamp     !
      !      both to physical ranges so the stage stays FINITE. This only bites on a genuinely oversized  !
      !      step (which the adaptive controller then rejects and shrinks normally) -- in range it is an  !
      !      identity, so no accuracy cost. Without the CAS clamp a big transient poisons the whole march !
      !      with NaN. --------------------------------------------------------------------------------!
      call state_extrap(y, BETA, Y2, n, nsl, base3)
      call clamp_theta(base3, frozen, nsl, nfire=clamp_n)
      call clamp_cas(base3, nfire=clamp_n)
      call column_be_stage(base3, frozen, n, nsl, GAMMA*dt, Y3, niter=np, bf=bf3, sf_out=sf3)
      call state_init(Y3, n, nsl, y_out)
      !----- operator-split mass: closed-form Euler over the FULL dt from y_n, using the SAME b-weighted !
      !      (1-gamma, gamma) per-cohort transp the CAS's own vapour balance used (bf2/bf3's ledger is    !
      !      b-weighted identically, sec 1/3/4/5) -- NOT a separate endpoint evaluation, so the mass       !
      !      debit and the CAS credit agree to within the tableau's own stage algebra. --------------------!
      transp_bw(1:n) = (1.0_wp - GAMMA)*sf2%transp_c(1:n) + GAMMA*sf3%transp_c(1:n)
      call advance_water_mass_full(y, frozen, n, nsl, dt, transp_bw, y_out)
      !----- operator-split canopy-SURFACE water (sec 3.4, P2c): same b-weighting discipline, using the  !
      !      SAME sf2/sf3 (already captured above for transp_bw) -- film_evap_leaf/wood are zero when     !
      !      canopy_water_on is off, so this is a no-op then. --------------------------------------------!
      film_evap_leaf_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_leaf(1:n) + GAMMA*sf3%film_evap_leaf(1:n)
      film_evap_wood_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_wood(1:n) + GAMMA*sf3%film_evap_wood(1:n)
      call advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
      !----- embedded 1st-order error estimate (mass is split out -> zeroed). -------------------------!
      call state_err_diff(Y3, base3, Y2, y, n, nsl, y_err)
      !----- b-weighted boundary-flux amounts: b^I = (0, 1-gamma, gamma) -> exact telescoping.         !
      !      (Water closure is exact only when clamp_theta is inactive; it barely moves over gamma*dt.) !
      if (present(bf)) then
         call bflux_bweight(bf, bf2, bf3, dt, GAMMA)
         !----- b-weighted TIME INTEGRAL of the tissue temperatures, same b^I = (0, 1-gamma, gamma)   !
         !      weights and the same telescoping argument as every other amount in bf. This is what    !
         !      the store's energy is set from -- see column_bflux_t's own note. --------------------!
         allocate(bf%tissue_leaf_int(n), bf%tissue_wood_int(n))
         bf%tissue_leaf_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%leaf_temp(1:n) + GAMMA*sf3%leaf_temp(1:n))
         bf%tissue_wood_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%wood_temp(1:n) + GAMMA*sf3%wood_temp(1:n))
      end if
   end subroutine ark2_column_step













   !---------------------------------------------------------------------------------------!
   ! adaptive_ark_march -- integrate to t_end with the ARK2 embedded error estimate driving the       !
   ! step controller (2 solves/step, no step-doubling). y_err (already the 2nd-1st order difference)   !
   ! is the local error; the WRMS of it vs tolerance drives accept/reject via adaptive_step_update     !
   ! (p=1 embedded -> exponent -1/2). Reports the step + reject count.                                 !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_ark_march(y0, frozen, n, nsl, t_end, ec, dt_init, y_out, nsteps, nrej, niter, acc, &
                                 dt_warm_out, clamp_n)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec       !< tolerances + controller + strictness (meds_fast_control)
      type(column_state_t),  intent(out) :: y_out
      integer(ik),           intent(out) :: nsteps, nrej
      integer(ik), optional, intent(in)  :: niter    !< coupled leaf<->CAS Newton cap (default 8)
      !----- section 8e WARM START: the controller proposal to seed the NEXT call with. Only steps that  !
      !      were not truncated by the end of the interval update it -- the last step of a march is       !
      !      usually a short remainder, and seeding from that would bias the next call small and undo     !
      !      the saving. Absent => caller does not want a warm start. ------------------------------------!
      real(wp),    optional, intent(out) :: dt_warm_out
      type(column_bflux_t), optional, intent(out) :: acc  !< accumulated boundary-flux amounts (ledger)
      !----- STAGE-clamp activations over the whole march, REJECTED trials included: a clamp on a step  !
      !      that was then thrown away still says the controller was probing too far, which is the      !
      !      signal wanted. Accumulated (intent(inout)), so the caller zeroes it. -----------------------!
      integer(ik), optional, intent(inout) :: clamp_n

      type(column_state_t) :: y, y_new, y_err, y_lo
      type(column_bflux_t) :: bfsub
      real(wp)             :: t, dt, err, err_prev, fac, dt_floor
      integer(ik)          :: np
      real(wp)             :: dt_try, dt_warm
      logical              :: clamped

      np = 8_ik ; if (present(niter)) np = max(1_ik, niter)
      if (present(acc)) call bflux_zero(acc, n)
      !----- substep FLOOR: bound the worst case to ~t_end/DT_FLOOR sub-steps. The ARK2 BE stages are    !
      !      L-stable, so a floor step is STABLE (bounded) even when the embedded error stays above tol   !
      !      -- e.g. a stiff transient the tolerance can't resolve. A tiny absolute floor (the old 1e-2s) !
      !      let a pathological step balloon to ~1.8e5 sub-steps and stall the march; t_end/64 caps it at !
      !      64 and degrades gracefully. (Also surfaces a genuine non-finite state promptly rather than   !
      !      grinding at the floor forever.) -----------------------------------------------------------!
      dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)

      call state_init(y0, n, nsl, y)
      t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik
      err_prev = -1.0_wp                                          ! < 0 => first step uses the I-controller
      dt_warm = dt                                                ! fallback if every step is end-clamped
      do
         if (t >= t_end - tiny_num) exit
         dt_try = dt
         dt = min(dt, t_end - t)
         clamped = dt < dt_try - tiny_num
         call ark2_column_step(y, frozen, n, nsl, dt, y_new, y_err, niter=np, bf=bfsub, clamp_n=clamp_n)
         call state_sub(y_new, y_err, n, nsl, y_lo)               ! the 1st-order embedded solution
         !----- per-group WRMS over the WHOLE column state (see state_wrms_grouped's header). The ARK's  !
         !      theta and water-mass terms are structurally zero here -- both ride operator-split maps     !
         !      outside the ESDIRK tableau -- so they dilute rather than inform. That is accepted          !
         !      deliberately: one norm, no per-scheme opt-outs, and the measured cost is 0-1% in accuracy  !
         !      against a 9-15% FALL in sub-steps (MEDS_INTEGRATOR_PARITY.md [RETIRED] §3f). --------------------!
         err = state_wrms_grouped(y_new, y_lo, y, n, nsl, ec%tols)
         !----- ROBUSTNESS: a non-finite err (a stage -- typically the BETA=2.414 stage-3 extrapolation    !
         !      base3 -- overshot the CAS enthalpy into a region where qsat(T) overflows) is a step that   !
         !      is simply TOO BIG: REJECT it and shrink dt deterministically (the NaN poisons the adaptive !
         !      fac, so use fmin directly). At a smaller dt, Y2 ~ y and base3 no longer overshoots, so the !
         !      step becomes finite and the march recovers -- the correct adaptive response, not a force-  !
         !      accept. Only if even a floor-sized step is non-finite do we commit + bail so meds_main's    !
         !      has_nan check reports it cleanly instead of the march hanging.                              !
         if (err /= err .or. dt /= dt) then
            if (dt <= dt_floor) then
               call state_init(y_new, n, nsl, y) ; t = t + dt_floor ; nsteps = nsteps + 1_ik ; exit
            end if
            nrej = nrej + 1_ik ; dt = max(dt * ec%fmin, dt_floor) ; cycle
         end if
         fac = step_control_factor(err, err_prev, ec)             ! I (default) or PI (Gustafsson) controller
         if (err <= 1.0_wp .or. dt <= dt_floor) then
            !----- L2 STRICT: a floor-forced accept that still breaches tolerance is a FAILURE to meet the  !
            !      requested accuracy -- fail hard rather than silently commit an under-resolved step (L1    !
            !      degrades gracefully; L2 is the faithful/validation mode). -------------------------------!
            if (ec%level == CTRL_L2_STRICT .and. err > 1.0_wp) &
               error stop 'adaptive_ark_march: L2 strict -- floor step cannot meet tolerance'
            call state_init(y_new, n, nsl, y)
            if (present(acc)) call bflux_add(acc, bfsub)          ! accumulate ONLY accepted substeps
            t = t + dt ; nsteps = nsteps + 1_ik
            err_prev = err                                        ! remember for the PI controller
            !----- WARM-START SEED = the step that was just ACCEPTED, recorded BEFORE the controller's    !
            !      growth factor is applied. Seeding the next call with the GROWN proposal (dt*fac, fac    !
            !      up to fmax=5) re-imports the very over-estimate that gets rejected -- measured, it left !
            !      the rejection rate at 26-29% instead of collapsing it. The last accepted size is the    !
            !      one with evidence behind it. ---------------------------------------------------------!
            if (.not. clamped) dt_warm = dt
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         if (nsteps + nrej > 4096_ik) exit                        ! hard backstop (should never trigger)
      end do
      call state_init(y, n, nsl, y_out)
      if (present(dt_warm_out)) dt_warm_out = dt_warm
   end subroutine adaptive_ark_march
   !=======================================================================================!
   !  INTEG_ARK path: the coupled IMEX-ARK fast step (docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md). Shares the   !
   !  split's frozen pre-pass (build_column_frozen), packs the state into the pure column vector,     !
   !  advances one dt_fast with the ARK stepper, then unpacks. PARTIAL rainfall>0 guard-lift: the ARK   !
   !  now carries the split's soil-boundary water-enthalpy advection (rain/runoff/drainage liquid      !
   !  enthalpy, in column_be_stage) and persists the scratch hydrology's ponding/aquifer/water-table   !
   !  (column_state_t still doesn't advance them prognostically -> a lagged operator split, so the      !
   !  whole-WATER budget closes only to the split-error tolerance, not machine). STILL restricted to   !
   !  free-drain + no Zeng-Decker: those bottom BCs need prognostic aquifer/z_wt in the state vector.  !
   !=======================================================================================!
   subroutine column_fast_step_ark(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,  &
                                   gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, cdiag)
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
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) per-cohort diagnostic capture
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters

      type(column_frozen_t)  :: frozen
      type(column_state_t)   :: y, y_out, ycur, ytmp, yerr
      type(surface_state_t)  :: y_stage
      type(surface_tend_t)   :: surf_tend
      type(column_bflux_t)   :: acc, bfsub
      real(wp)    :: tg, fl, dt0, cas_mass_capacity, cas_molar_capacity, enth0, shv0, co20, enth1, shv1, co21, e_soil0, e_soil1, &
           w_soil0, w_soil1
      real(wp)    :: e_pond0, e_pond1   !< pond enthalpy store, start/end (#78 item 4)
      real(wp)    :: w_plant0, w_plant1
      real(wp)    :: w_surface0, dt_warm_next
      !----- Canopy-SURFACE water (sec 3.4, P2c) ledger scratch. --------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, leaf_cap_i, wood_cap_i, intercept_total
      !----- Row-1b condensate actually DEPOSITED into soil layer 1 this step: mass [kg/m2] and the      !
      !      paired enthalpy [J/m2]. Held in locals rather than recomputed at the budget calls so the      !
      !      inflow term below is the IDENTICAL number the state update applied -- the "one flux, both    !
      !      sides" discipline the whole design rests on. --------------------------------------------!
      real(wp)    :: cond_dep_mass, cond_dep_enth
      real(wp)    :: tissue_store0, tissue_store1
      real(wp)    :: cap_leaf_a(col_cohort%n), cap_wood_a(col_cohort%n)
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, k, i, isub, nsub, nsteps, nrej
      logical     :: halt_budgets     !< §5.1: hard-stop on a non-closing budget (full column + debug only)

      n = col_cohort%n ; nsl = col_config%soil%n_active
      !----- CLAMP counters are ACCUMULATED down the call chain (intent(inout) all the way to the      !
      !      kernels), so this sub-step's tally starts clean here. The commit counters stay 0 on this   !
      !      path: the ARK clamps its stage-3 extrapolation base only, never y_out. -------------------!
      budget%clamp_stage_n = 0_ik ; budget%clamp_commit_n = 0_ik ; budget%theta_ood_max = 0.0_wp
      budget%clamp_mass    = 0.0_wp ; budget%clamp_energy = 0.0_wp

      !----- §5.1: a reduced column freezes a store while its fluxes still act on the neighbours, so it  !
      !      cannot conserve by construction -- suppress the HARD stops (soft n_fail counters still run). !
      halt_budgets = col_config%energy%debug_error .and. mask_is_full(col_config%mask)

      !----- The bottom-BC guard is GONE (Phase 0/3). All three BCs are now pure boundary conditions   !
      !      with no prognostic state behind them: bedrock seals the face, free drainage takes the       !
      !      unit-gradient limit, and the aquifer is head-driven against a saturated zone at the column  !
      !      base. The ARK commits the scratch advance_soil_water_column theta verbatim, so it inherits all  !
      !      three unchanged. -------------------------------------------------------------------------!
      w_surface0 = biophys%soil_w%w_surface

      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, n, nsl, &
                               frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      !----- STATE^n pond enthalpy, snapshotted HERE (#78 item 4). It must be read before anything      !
      !      commits to biophys%soil_w: the pond commit below overwrites w_surface_enth with the scratch's  !
      !      end-of-step value, so a snapshot taken alongside e_soil0 further down reads the COMMITTED  !
      !      value and the store change collapses to exactly zero. That is what it did -- the whole-    !
      !      column residual came out equal to -e_pond0 to six figures, which is what identified it.    !
      !      build_column_frozen mutates biophys%snow (and nothing in biophys%soil_w), so state^n is intact.    !
      e_pond0 = biophys%soil_w%w_surface_enth

      !----- advance one dt_fast: adaptive (embedded-error) or GPU-warp-uniform fixed substeps. ----!
      if (col_config%integrator%adaptive) then
         !----- section 8e WARM START: seed from the step the controller converged to on the PREVIOUS      !
         !      dt_fast call for this patch. Column stiffness barely changes call to call, so the cold     !
         !      start was re-discovering the same step size every call and paying ~1 rejection for it.     !
         !      An explicit [fast].ark_dt_init still wins (it is an operator override); with neither, the  !
         !      first call for a patch cold-starts exactly as before. ---------------------------------!
         dt0 = dt_fast
         if (biophys%adapt_dt_last > tiny_num) dt0 = min(biophys%adapt_dt_last, dt_fast)
         if (col_config%integrator%dt_init > tiny_num) dt0 = min(col_config%integrator%dt_init, dt_fast)
         !----- The UNIFIED error-control bundle (§8c Layer 1): build_error_control seeds every tolerance !
         !      group from the setting that governs it today (and honours the [fast].rtol_all master      !
         !      dial), plus the controller + strictness. Defaults (CTRL_I, CTRL_L1, rtol_all unset)       !
         !      reproduce the legacy march byte-for-byte. ------------------------------------------------!
         ec = col_config%integrator%error_control
         call adaptive_ark_march(y, frozen, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                                 niter=merge(NEWT_COUPLED, 1_ik, col_config%integrator%coupled_newton), acc=acc, &
                                 dt_warm_out=dt_warm_next,                                          &
                                 clamp_n=budget%clamp_stage_n)
         biophys%adapt_dt_last = dt_warm_next
      else
         nsub = max(1_ik, col_config%integrator%fixed_substeps) ; nrej = 0_ik ; ycur = y ; call bflux_zero(acc, n)
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, frozen, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr,          &
                                  niter=merge(NEWT_COUPLED, 1_ik, col_config%integrator%coupled_newton), bf=bfsub, &
                                  clamp_n=budget%clamp_stage_n)
            call bflux_add(acc, bfsub)
            ycur = ytmp
         end do
         y_out = ycur ; nsteps = nsub
      end if
      !----- section 5.3 WORK counters: record what the march actually cost. hydro_nsub/hydro_nonconv  !
      !      are set ONCE by the Act-1 pre-pass's solve_plant_water_batch call (build_column_frozen),     !
      !      NOT here -- there is no more per-stage hydraulics solve to accumulate over sub-steps. -------!
      budget%integ_nsteps = nsteps ; budget%integ_nrej = nrej

      !----- SOIL WATER is operator-split out: the ESDIRK stages passed theta through unchanged (=theta^n); !
      !      commit the AUTHORITATIVE end-of-step theta from the scratch advance_soil_water_column HERE, once,  !
      !      so a single consistent theta feeds the state commit, the soil_temp read-off, and BOTH the      !
      !      soil_water and whole_water storage terms (w_soil1 below). ------------------------------------!
      y_out%theta(1:nsl) = frozen%hydrology%theta1(1:nsl)
      !----- pond: same treatment as theta -- committed from the scratch solve, but THROUGH the state !
      !      vector so there is one authority. Phase 1 replaces these two lines with a stage RHS. -----!
      y_out%w_surface      = frozen%hydrology%w_surface1
      y_out%w_surface_enth = frozen%hydrology%w_surface_enth1     ! #78 item 4: paired with the mass

      !----- §5.1 PROCESS MASK. The mask must mean the same thing under every scheme, so it is applied  !
      !      at the ARK's single state-commit point: a masked-off component is restored to state^n (y),  !
      !      leaving the ODE one dimension smaller while its couplings still acted during the march.     !
      !      This mirrors the split path's freeze exactly. mask%veg_energy needs no case here -- the ARK !
      !      error-stops on prognostic leaf/wood, so no vegetation energy store exists on this path. ----!
      call apply_process_mask(col_config%mask, y, y_out, n, nsl)

      !----- Canopy-SURFACE water (sec 3.4, P2c): capacity clamp + overflow/deficit bookkeeping (mirrors   !
      !      the split path's own post-hoc treatment, sec 9's "clamp, don't silently over-apply"). DEFICIT  !
      !      (the floor, symmetric to overflow): film_evap uses a state^n-frozen conductance (rescaled by     !
      !      availability in build_column_frozen, but only an approximation), so the store can still           !
      !      transiently overdraw below 0 -- left uncorrected, this corrupts the NEXT dt_fast's frozen           !
      !      interception rate (intercept_canopy_layer's own internal floor silently "fixes" a negative           !
      !      starting bucket, fabricating mass). Floor at 0 and bookkeep the shortfall as a NEGATIVE               !
      !      addition to the outflow ledger (the exact mirror of surf_overflow's sign). Gated behind                !
      !      canopy_water_on per the P1 nvfortran lesson. -------------------------------------------------------!
      surf_overflow = 0.0_wp ; surf_deficit = 0.0_wp
      if (col_config%canopy_water_on) then
         call clamp_canopy_film(y_out, col_cohort%lai, col_cohort%wai, col_config%soil_water_opts%dewmx, n, surf_overflow, &
                                surf_deficit)
      end if

      !----- unpack into biophys + re-derive the diagnostic soil temperatures + leaf temperatures. -----!
      call unpack_column_state(y_out, n, nsl, biophys)
      !----- persist the scratch hydrology's ponding/aquifer/water-table (lagged operator split).       !
      !      §5.1: these are part of the SOIL-WATER store, so they must obey the same freeze as theta   !
      !      -- otherwise mask%soil_water=.false. means something different on this path than on the    !
      !      split path (which restores the whole soil_column_t), and the two schemes are no longer     !
      !      running the same reduced system. They are the ONLY writes to biophys%soil_w besides theta, and !
      !      the hydrology ran on soil_w_scratch, so skipping them leaves the store at state^n. --------!
      if (col_config%mask%soil_water) then
         biophys%soil_w%w_surface      = y_out%w_surface
         biophys%soil_w%w_surface_enth = y_out%w_surface_enth
      end if
      call diagnose_soil_temps(y_out, col_config%soil_thermal%soil_dry_heat_capacity, nsl, biophys%soil_e%soil_temp, &
            biophys%soil_e%soil_fliq)
      call internal_energy_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        col_config%soil_thermal%soil_dry_heat_capacity(1), tg, fl)
      y_stage%cas_enthalpy = y_out%cas_enthalpy ; y_stage%cas_shv = y_out%cas_shv ; y_stage%cas_co2 = y_out%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, tg, n, surf_tend)
      !----- Commit the tissue temperatures the frozen store already produced. surf_tend%leaf_temp/wood_temp !
      !      ARE the dt_fast endpoints: surface_derivs solved the balance with leaf_hcap_per_dt/wood_hcap_per_dt = cap/dt   !
      !      relaxing from frozen%tissue%t_leaf0/t_wood0, so the store is already inside the CAS solve and    !
      !      there is nothing to correct afterwards. An earlier attempt applied the store as a POST-     !
      !      COMMIT lump on the committed CAS; that is unstable, because the tissue and the canopy air   !
      !      have comparable heat capacities (cap_wood/cap_cas ~ 0.5 here) and the tissue relaxes far    !
      !      inside one step, so dumping its demand on an already-solved CAS drove a growing oscillation !
      !      (measured: 49 K CAS kicks per step, wood diverging to 247 K against a 331 K balance).       !
      !      §5.1 mask: freezing veg_energy holds the tissue at its entry temperature while the fluxes   !
      !      it drives into the CAS are kept. -----------------------------------------------------------!
      tissue_store0 = 0.0_wp ; tissue_store1 = 0.0_wp
      do i = 1_ik, n
         !----- Derive the capacity from the SAME store_hcap_per_dt the kernel relaxed against, not by         !
         !      recomputing it from dry_hcap + wmass. The two agree by construction today, but only    !
         !      this form guarantees that zeroing leaf_hcap_per_dt/wood_hcap_per_dt zeroes the ledger's store term too --  !
         !      i.e. that "no capacity" is a clean no-op end to end rather than a state change the     !
         !      fluxes never paid for. --------------------------------------------------------------!
         cap_leaf_a(i) = frozen%tissue%leaf_hcap_per_dt(i) * dt_fast
         cap_wood_a(i) = frozen%tissue%wood_hcap_per_dt(i) * dt_fast
         tissue_store0 = tissue_store0 + cap_leaf_a(i) * frozen%tissue%t_leaf0(i)                          &
                                       + cap_wood_a(i) * frozen%tissue%t_wood0(i)
      end do
      !----- Commit the TIME-AVERAGE of the tissue temperature over the march, not the final stage's  !
      !      value. That is the temperature whose store energy equals the b-weighted complement of    !
      !      the flux the canopy air actually received, so state and ledger are the same number and   !
      !      the whole-column energy closes to machine precision. Taking surf_tend%leaf_temp (the last       !
      !      stage) instead leaves ~3.9e4 J unaccounted per step -- 4e-4 relative, against a 1e-6     !
      !      tolerance -- because the CAS moves across the stages. -----------------------------------!
      if (col_config%mask%veg_energy) then
         if (allocated(acc%tissue_leaf_int)) then
            biophys%leaf_temp(1:n) = acc%tissue_leaf_int(1:n) / dt_fast
            biophys%wood_temp(1:n) = acc%tissue_wood_int(1:n) / dt_fast
         else
            biophys%leaf_temp(1:n) = surf_tend%leaf_temp(1:n) ; biophys%wood_temp(1:n) = surf_tend%wood_temp(1:n)
         end if
      else
         biophys%leaf_temp(1:n) = frozen%tissue%t_leaf0(1:n) ; biophys%wood_temp(1:n) = frozen%tissue%t_wood0(1:n)
      end if
      do i = 1_ik, n
         tissue_store1 = tissue_store1 + cap_leaf_a(i) * biophys%leaf_temp(i)                             &
                                       + cap_wood_a(i) * biophys%wood_temp(i)
      end do

      !----- WHOLE-COLUMN CONSERVATION LEDGER: close the same 7 budgets the split closes, using the     !
      !      b-weighted boundary-flux AMOUNTS accumulated over the substeps (acc). The flux-form CAS    !
      !      commits + the energy_resid=0 soil-heat column make the identity exact -> machine-precision !
      !      closure for ENERGY (incl. the frozen rain/runoff/drainage advection, a fixed source). dt=1  !
      !      because acc holds AMOUNTS, not rates. Whole-WATER carries the lagged ponding split, so it   !
      !      closes only to the operator-split tolerance. ---------------------------------------------!
      !----- ROW 1b: DEPOSIT THE CONDENSATE (same routing on RK45; see the rationale below). Full      !
      !      rationale). Dew/fog landed on a surface inside the column; it used to be booked into        !
      !      whole_wat_out/whole_enth_out and vanish. Paired mass + enthalpy into soil layer 1 at the    !
      !      CAS temperature it condensed at, so the whole-column ledger closes with no boundary term.   !
      !      Applied to y_out BEFORE the store snapshot below, and mirrored into biophys (already unpacked). !
      !      Same destination on all three paths -- routing it elsewhere here would reopen a scheme      !
      !      asymmetry. acc%whole_cond is 0 whenever cas_condensation is off. ------------------------!
      cond_dep_mass = 0.0_wp ; cond_dep_enth = 0.0_wp
      if (acc%whole_cond > 0.0_wp) then
         cond_dep_mass = acc%whole_cond
         cond_dep_enth = acc%whole_cond_enth       ! b-weighted at the stage CAS temperatures (see column_bflux_t)
         call deposit_condensate(y_out, col_config%soil%dz(1), cond_dep_mass, cond_dep_enth)
         biophys%soil_w%theta(1)       = y_out%theta(1)
         biophys%soil_e%soil_energy(1) = y_out%soil_energy(1)
      end if

      !----- No wood_rnet term any more: the wood's absorbed radiation is ALWAYS carried by            !
      !      surface_derivs' own balance (its frozen diagnostic inputs are now filled unconditionally), !
      !      so it is already inside coh_rnet / acc%whole_enth_in. The old term existed only because    !
      !      the prognostic-wood branch zeroed those inputs and owned the radiation separately. -------!

      cas_mass_capacity = frozen%cas%cas_mass_capacity ; cas_molar_capacity = frozen%cas%cas_molar_capacity
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv ; co21 = y_out%cas_co2
      e_pond1 = frozen%hydrology%w_surface_enth1
      e_soil0 = soil_energy_store(y%soil_energy,     col_config%soil%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, col_config%soil%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     col_config%soil%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, col_config%soil%dz, nsl)
      !----- Plant internal water MASS is now a genuine store the march evolves (MEDS_ED2_RK45_       !
      !      DESIGN.md sec 1/3/4/5, P2): it absorbs exactly the transp<->uptake mismatch the OLD        !
      !      tolerance inflation below used to paper over (the SAME closure gain split's P0 already      !
      !      has). Omitting it here would make the store's own real change read as a leak. -------------!
      w_plant0 = plant_water_store(col_cohort%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(col_cohort%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- Canopy-SURFACE water (sec 3.4, P2c): already ground-area-referenced (no nplant factor,     !
      !      unlike w_plant0/1 above). Valued at u_liq(t_film_valuation) = frozen%film%film_liquid_enthalpy, the liquid       !
      !      enthalpy the intercepted water arrived with; the tissue pays enthalpy_vapor - film_liquid_enthalpy per  !
      !      kg it evaporates (surface_derivs), so this store closes exactly against the CAS credit. All   !
      !      zero when canopy_water_on is off. -------------------------------------------------------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))
      !----- L2/debug_error mode (col_config%energy%debug_error) promotes a non-closing budget from a       !
      !      silently-counted n_fail to a hard `error stop` -- the enforced half of the conservation   !
      !      check (plan MEDS_NUMERICS_SCOPING.md sec 4/QW2), mirroring the split path; off by         !
      !      default so production behaviour is unchanged. Each check reuses budget%*%resid, which        !
      !      budget_accumulate just set as a side effect. --------------------------------------------!
      !----- Tolerances are FLUX-scaled (meds_budget_check header): rtol * gross boundary flux over  !
      !      the step plus a rate floor * dt_fast. Store-scaled tolerances let a ~1 W/m2 leak through. !
      budget%atm_heat_export = acc%atm_heat_out ; budget%atm_vap_export = acc%atm_vap_out
      call budget_check(budget%cas_energy, cas_mass_capacity*enth0, cas_mass_capacity*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                        dt_fast, budget_energy_rate_floor, 'cas_energy (ark)', halt_budgets)
      call budget_check(budget%cas_water,  cas_mass_capacity*shv0,  cas_mass_capacity*shv1,  acc%cas_vap_in,  acc%cas_vap_out, &
                        dt_fast, budget_water_rate_floor, 'cas_water (ark)', halt_budgets)
      call budget_check(budget%cas_co2,    cas_molar_capacity*co20,  cas_molar_capacity*co21,  acc%cas_co2_in,  acc%cas_co2_out, &
                        dt_fast, budget_co2_rate_floor, 'cas_co2 (ark)', halt_budgets)
      !----- cond_dep_enth is a boundary INPUT to the SOIL store, and it has to be said here even though  !
      !      the whole-column ledger needs no term for it.  The row-1b deposit moves condensate CAS ->    !
      !      soil layer 1 AFTER the march: whole-column sees an internal transfer between two stores it   !
      !      already tracks, so it telescopes there; but this per-kernel budget sees only the soil, for    !
      !      which the same transfer is an inflow.  Omitting it made the soil budget short by exactly the  !
      !      deposit -- measured +2.95e3 J/m2 (winter) / +4.65e3 (summer) on a forced month, with the      !
      !      soil_water twin short by the paired mass and the two residuals' RATIO landing exactly on      !
      !      internal_energy_liquid(T_cas).  ARK is the only path that noticed: RK45 keeps no per-kernel   !
      !      soil budgets (every store rides the same column_derivs RHS, so its whole-column ledger IS     !
      !      the per-store one) and split fills these from the KERNEL's own residual, which for           !
      !      soil_energy_step_implicit is zero by construction and so cannot see a post-solve deposit. ----!
      call budget_check(budget%soil_energy, e_soil0, e_soil1, acc%soil_enth_in + cond_dep_enth,            &
                        acc%soil_enth_out, dt_fast, budget_energy_rate_floor, 'soil_energy (ark)', halt_budgets)
      !----- SOIL WATER (fully frozen now): storage theta^n -> theta1 (w_soil0 -> w_soil1, both from the    !
      !      scratch solve), inflow q_top*rho, outflow drainage + realized uptake -- all from the frozen    !
      !      hflux, which closed its OWN mass budget to machine precision inside advance_soil_water_column. -----!
      !----- ...and the paired MASS, for the same reason (see the soil_energy note just above). ---------!
      !----- ...plus the scratch's own two post-solve corrections: clipped water LEFT the soil for the  !
      !      pond, floored water was CREATED in it -- both are in the committed theta1. -----------------!
      call budget_check(budget%soil_water,  w_soil0, w_soil1,                                              &
                        (frozen%hydrology%q_top*rho_h2o + frozen%hydrology%floor_mass)*dt_fast + cond_dep_mass,                    &
                        (frozen%hydrology%drainage + frozen%roots%uptake                                    &
                         + frozen%hydrology%clip_mass)*dt_fast,                                             &
                        dt_fast, budget_water_rate_floor, 'soil_water (ark)', halt_budgets)
      !----- whole-WATER: rainfall IN; drainage + runoff + CAS-vapour OUT; ponding + plant internal water !
      !      MASS in the store. The soil + ponding + drainage/runoff/rainfall terms are frozen fast-step    !
      !      amounts; the CAS-vapour exchange g_atm_vapour*(shv-shv_atm) is the ARK-accumulated part (acc%           !
      !      whole_wat_out). The plant's OWN w_plant0/1 store term (above) now absorbs the transp<->        !
      !      uptake mismatch the OLD comment here described (the ARK re-evaluates transpiration per          !
      !      ESDIRK stage as the CAS VPD evolves, while the committed soil theta lost the FROZEN scratch      !
      !      uptake_total) -- exactly the closure split's P0 already has, so the tolerance inflation that     !
      !      used to paper over that gap is gone; this closes to machine precision like the other 6 now. ----!
      !----- SNOW joins the ledgers (C4): the pack is a real mass + energy STORE, its accumulated      !
      !      rainfall enthalpy is a boundary INPUT, and SNOWFALL is a boundary water input -- forc%snowfall  !
      !      was absent from w_in here (split has always carried it), so with a pack the ledger saw     !
      !      mass appear with no source and leaked exactly snowfall*dt every step. Sublimation already     !
      !      leaves via the CAS vapour term and meltwater already moved pack -> soil as a paired        !
      !      transfer, so neither needs a term. All zero without snow. --------------------------------!
      !----- frozen%hydrology%floor_mass: the theta_res floor's water is created inside the column and enters as a   !
      !      boundary INPUT, exactly as its enthalpy (e_floor -> acc%whole_enth_in) already did; the      !
      !      mass half was missing (2026-09 review, item 1A #5). -------------------------------------!
      call budget_check(budget%whole_water,                                                                &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0 + frozen%snow%swe0,  &
                        w_soil1 + cas_mass_capacity*shv1 + frozen%hydrology%w_surface1 + w_plant1 + surf_water1 &
                        + frozen%snow%swe1, &
                        acc%whole_wat_in + (forc%rainfall + forc%snowfall + biophys%shed_water_rate                &
                                            + frozen%hydrology%floor_mass)*dt_fast,                                    &
                        acc%whole_wat_out                                                                   &
                        + (frozen%hydrology%runoff_surf + frozen%hydrology%drainage)*dt_fast                &
                                          + surf_overflow - surf_deficit,                                 &
                        dt_fast, budget_water_rate_floor, 'whole_water (ark)', halt_budgets)
      call budget_check(budget%whole_energy,                                                               &
                             !----- No melt rebase any more (#78 item 4): the pack hands its meltwater to  !
                             !      the POND, not to soil layer 1, so e_soil0 no longer contains the melt   !
                             !      enthalpy and the pack/pond pair telescopes on its own. -------------!
                        e_soil0                           + cas_mass_capacity*enth0 + surf_enth0                     &
                        + frozen%snow%enth0 + e_pond0 + tissue_store0,                               &
                        e_soil1 + cas_mass_capacity*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond1                &
                        + tissue_store1,                                                                 &
                        acc%whole_enth_in + intercept_total*dt_fast*internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                          + frozen%snow%acc_enth                                       &
                                          + (frozen%hydrology%precip_ground - frozen%snow%melt_rate)*dt_fast         &
                                            * internal_energy_liquid(frozen%hydrology%t_pond_inflow),                      &
                        acc%whole_enth_out + (surf_overflow - surf_deficit)                                 &
                        * internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                           + frozen%hydrology%runoff_enth*dt_fast,                                    &
                        dt_fast, budget_energy_rate_floor, 'whole_energy (ark)', halt_budgets)

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_ark

end module meds_fast_ark
