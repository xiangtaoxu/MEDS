! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_ark -- the production fast-loop integrator: an L-stable two-stage ESDIRK2      !
! (ARS(2,2,2) tableau). The config name is "ark", but the explicit part of the tableau is  !
! empty, so this is a diagonally implicit scheme, not an IMEX one. The step, its error     !
! estimate, the error norm and the step controller: docs/science/numerical_scheme.md §3    !
! `ark`, eqs. (1)-(4).                                                                     !
!                                                                                          !
! Hosts the TABLEAU AND THE MARCH, and nothing else: the scheme entry                      !
! `column_fast_step_ark` (called from `meds_fast_step`), the two-stage step                !
! `ark2_column_step`, and its embedded-error controller `adaptive_ark_march`.              !
!                                                                                          !
! ARK-ONLY. Everything the two schemes share lives elsewhere, because importing it from a  !
! module named for one scheme makes siblings look related:                                 !
!   * the frozen work record builder         -> meds_fast_frozen    (RK45 builds it too)   !
!   * the implicit stage + its Newton        -> meds_fast_be_stage  (the RHS test uses it) !
!   * the water-mass and canopy-film advance -> meds_fast_be_stage                         !
!   * the column_state_t algebra + ledgers   -> meds_column_state_ops                      !
! Nothing outside this module imports `meds_fast_ark` except the dispatcher and the RHS    !
! test (test_column_derivs).                                                               !
!                                                                                          !
! Soil and plant water are the soil-plant water network's (meds_spac_solver, run in        !
! build_column_frozen): the stages move theta at its steady rate, plant water is advanced  !
! at its rates once per step, and the end-of-step theta and pond are committed from it     !
! once per dt_fast. Whole-column water and energy close to round-off on every bottom BC    !
! (drainage, free to sealed, or aquifer; docs/science/numerical_scheme.md §5b); see        !
! meds_budget_check and the ledgers at the end of column_fast_step_ark.                    !
!==========================================================================================!
module meds_fast_ark
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o
   use meds_plant_hydr_lib,  only : tissue_psi_to_water
   use meds_config, only : meds_config_t, INTEG_ARK, CTRL_L2_STRICT
   use meds_fast_control, only : state_wrms_grouped, step_control_factor
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_soil_types, only : snow_melt_t
   use meds_fast_types, only : patch_biophys_t
   use meds_fast_time_derivs, only : surface_derivs
   use meds_column_state_ops, only : state_init, state_extrap, state_err_diff, state_sub, bflux_zero, bflux_add, &
                                     bflux_bweight, clamp_cas, soil_water_store,                    &
                                     soil_energy_store, plant_water_store, canopy_film_store, deposit_condensate, &
                                     clamp_canopy_film, unpack_column_state, diagnose_soil_temps, &
                                     assemble_soil_energy_forcing, apply_process_mask, transp_per_leaf_area, &
                                     commit_tissue_energy
   use meds_fast_be_stage, only : column_be_stage, advance_plant_water, advance_surf_water_full
   use meds_fast_frozen, only : build_column_frozen
   use meds_diag_vars,       only : CO_TRANSP
   use meds_fast_types, only : column_config_t, column_cohort_t, column_forcing_t, column_budget_t, alloc_column_cohort, &
                               column_state_t, column_frozen_t, surface_state_t, cas_boundary_t, surface_tend_t, &
                               stage_bflux_t, column_bflux_t, error_control_t, column_tend_t, mask_is_full, &
                               ark_stages_t, fast_workspace_t
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction, ground_surface_fluxes
   use meds_therm_lib, only : internal_energy_liquid, internal_energy_to_temp, internal_energy_ice, liquid_energy_to_temp
   use meds_budget_check, only : budget_check, BUDGET_ENERGY_RATE_FLOOR, BUDGET_WATER_RATE_FLOOR, BUDGET_CO2_RATE_FLOOR
   implicit none
   private

   !----- `niter` on column_be_stage / ark2_column_step / adaptive_ark_march is a two-valued port:     !
   !      <= 1 selects the uncoupled single-BE pass, anything > 1 selects the coupled 2x2 Newton       !
   !      (whose own iteration cap is NEWT_MAX, not this number). NEWT_COUPLED names the second        !
   !      value, so it does not read as a tunable count; [fast].ark_coupled picks between the two. ----!
   integer(ik), parameter :: NEWT_COUPLED = 2_ik

   public :: column_fast_step_ark
   public :: ark2_column_step, adaptive_ark_march

contains



   !---------------------------------------------------------------------------------------!
   ! ark2_column_step -- one second-order, L-stable ESDIRK2 step, docs/science/numerical_scheme.md     !
   ! eq. (1): the ARS(2,2,2) tableau (Ascher-Ruuth-Spiteri 1997, Appl.Numer.Math. 25:151; the same     !
   ! gamma as Giraldo et al. 2013 "ARK2"). Each stage is one backward-Euler solve over gamma*dt        !
   ! (column_be_stage), and the method is stiffly accurate, so its last stage, Y3, is the step. The    !
   ! biotic CO2 source is folded into the implicit CO2 solve, so the explicit tableau is empty. Soil   !
   ! and plant WATER are the soil-plant water network's: the stages move soil water at its steady      !
   ! rate and pass plant water through, which is then advanced once over the full dt at the            !
   ! network's rates (advance_plant_water); neither adds to the embedded error. y_err =                !
   ! (Y3-base3)-(Y2-y_n) is the free first-order error estimate, eq. (2) (2 solves per step against    !
   ! step doubling's 3). ------------------------------------------------------------------------------!
   subroutine ark2_column_step(y, frozen, col_config, n, nsl, dt, y_out, y_err, stages, niter, bf, clamp_n)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(inout) :: y_out, y_err   !< filled in place
      type(ark_stages_t),    intent(inout) :: stages    !< this thread's stage storage, filled in place
      integer(ik), optional, intent(in)  :: niter
      type(column_bflux_t), optional, intent(out) :: bf   !< b-weighted boundary-flux AMOUNTS over dt (ledger)
      !----- STAGE-clamp activations (see column_budget_t%clamp_stage_n). The ARK clamps its ARS       !
      !      stage-3 extrapolation base and nothing else; y_out is committed unclamped. ---------------!
      integer(ik), optional, intent(inout) :: clamp_n
      real(wp), parameter :: GAMMA = 0.2928932188134524_wp   ! 1 - 1/sqrt(2)
      real(wp), parameter :: BETA  = 2.4142135623730951_wp   ! (1-gamma)/gamma = 1 + sqrt(2)
      type(stage_bflux_t)   :: bf2, bf3
      real(wp)              :: transp_bw(n)
      real(wp)              :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      integer(ik) :: np
      associate (Y2 => stages%y2, base3 => stages%base3, Y3 => stages%y3, sf2 => stages%sf2, sf3 => stages%sf3)
         np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

         !----- Stage 2: gamma*dt BE stage from y_n (CAS+soil only; mass frozen -- it is split out). -----!
         call column_be_stage(y, frozen, col_config, n, nsl, GAMMA*dt, Y2, sf2, niter=np, bf=bf2)
         !----- Stage 3: extrapolated base. The BETA=2.414 extrapolation can overshoot the CAS enthalpy   !
         !      into a wild temperature where qsat(T) overflows to NaN; clamp it to a physical range so the  !
         !      stage stays FINITE. This only bites on a genuinely oversized step (which the adaptive        !
         !      controller then rejects and shrinks normally) -- in range it is an identity, so no accuracy  !
         !      cost. Without the CAS clamp a big transient poisons the whole march with NaN. Soil water      !
         !      needs no clamp: it moves at the water network's steady rate, so the extrapolation lands      !
         !      between its start and end values. ---------------------------------------------------!
         call state_extrap(y, BETA, Y2, n, nsl, base3)
         call clamp_cas(base3, nfire=clamp_n)
         call column_be_stage(base3, frozen, col_config, n, nsl, GAMMA*dt, Y3, sf3, niter=np, bf=bf3)
         call state_init(Y3, n, nsl, y_out)
         !----- Plant water over the FULL dt at the network's rates, debited the SAME b-weighted          !
         !      (1-gamma, gamma) per-cohort transpiration the CAS's own vapour balance was credited        !
         !      (bf2/bf3's ledger is b-weighted identically), so debit and credit are one number. --------!
         transp_bw(1:n) = (1.0_wp - GAMMA)*sf2%transp_c(1:n) + GAMMA*sf3%transp_c(1:n)
         call advance_plant_water(y, frozen, n, dt, transp_bw, y_out)
         !----- Operator-split canopy-SURFACE water: the same b-weighting, from the same sf2/sf3. With    !
         !      canopy_water_on off, film_evap_leaf/wood are zero and this is a no-op. -------------------!
         film_evap_leaf_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_leaf(1:n) + GAMMA*sf3%film_evap_leaf(1:n)
         film_evap_wood_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_wood(1:n) + GAMMA*sf3%film_evap_wood(1:n)
         call advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
         !----- Embedded first-order error estimate, eq. (2); the plant water, films and pond are not     !
         !      counted (state_entry_rules). -------------------------------------------------------------!
         call state_err_diff(Y3, base3, Y2, y, n, nsl, y_err)
         !----- b-weighted boundary-flux amounts: b^I = (0, 1-gamma, gamma) -> exact telescoping.         !
         if (present(bf)) then
            call bflux_bweight(bf, bf2, bf3, dt, GAMMA)
            !----- b-weighted TIME INTEGRAL of the tissue temperatures, same b^I = (0, 1-gamma, gamma)   !
            !      weights and the same telescoping argument as every other amount in bf. This is what    !
            !      the store's energy is set from -- see column_bflux_t's own note. --------------------!
            allocate(bf%tissue_leaf_int(n), bf%tissue_wood_int(n))
            bf%tissue_leaf_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%leaf_temp(1:n) + GAMMA*sf3%leaf_temp(1:n))
            bf%tissue_wood_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%wood_temp(1:n) + GAMMA*sf3%wood_temp(1:n))
            !----- ...and each cohort's transpiration, the transp_bw the leaf water was debited. -------!
            allocate(bf%transp_coh(n))
            bf%transp_coh(1:n) = dt * transp_bw(1:n)
         end if
      end associate
   end subroutine ark2_column_step













   !---------------------------------------------------------------------------------------!
   ! adaptive_ark_march -- integrate to t_end with the embedded error estimate driving the step        !
   ! controller (2 solves per step, no step doubling). y_err (already the second-minus-first-order     !
   ! difference) is the local error; its grouped weighted RMS against the tolerances                   !
   ! (docs/science/numerical_scheme.md eq. (3)) drives accept/reject, and step_control_factor          !
   ! sizes the next step, eq. (4) with p = 1. Reports the step and reject counts.                      !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_ark_march(y0, frozen, col_config, n, nsl, t_end, ec, dt_init, y_out, ws, nsteps, nrej, niter, acc, &
                                 dt_warm_out, clamp_n)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec       !< tolerances + controller + strictness (meds_fast_control)
      type(column_state_t),  intent(inout) :: y_out  !< filled in place
      type(fast_workspace_t), intent(inout) :: ws    !< this thread's march and stage storage
      integer(ik),           intent(out) :: nsteps, nrej
      integer(ik), optional, intent(in)  :: niter    !< <= 1 uncoupled pass, > 1 (default) coupled Newton
      !----- WARM START (docs/science/numerical_scheme.md §3 `ark`): the step size to seed the NEXT       !
      !      call with, the last accepted step not truncated by the end of the interval -- the last       !
      !      step of a march is usually a short remainder, and seeding from it would bias the next        !
      !      call small. Absent => the caller does not want a warm start. --------------------------------!
      real(wp),    optional, intent(out) :: dt_warm_out
      type(column_bflux_t), optional, intent(out) :: acc  !< accumulated boundary-flux amounts (ledger)
      !----- STAGE-clamp activations over the whole march, REJECTED trials included: a clamp on a step  !
      !      that was then thrown away still says the controller was probing too far, which is the      !
      !      signal wanted. Accumulated (intent(inout)), so the caller zeroes it. -----------------------!
      integer(ik), optional, intent(inout) :: clamp_n

      type(column_bflux_t) :: bfsub
      real(wp)             :: t, dt, err, err_prev, fac, dt_floor
      integer(ik)          :: np
      real(wp)             :: dt_try, dt_warm
      logical              :: clamped
      associate (y => ws%y_march, y_new => ws%y_new, y_err => ws%y_err, y_lo => ws%y_lo)
         np = 8_ik ; if (present(niter)) np = max(1_ik, niter)
         if (present(acc)) call bflux_zero(acc, n)
         !----- Sub-step FLOOR, t_end/64 (never below 1e-2 s). The BE stages are L-stable, so a floor-size   !
         !      step is STABLE (bounded) even when the embedded error stays above tolerance -- e.g. a        !
         !      stiff transient the tolerance cannot resolve. The 1e-2 s bound alone would let a             !
         !      pathological step balloon to ~1.8e5 sub-steps and stall the march; t_end/64 caps it at 64    !
         !      and degrades gracefully, and surfaces a genuine non-finite state promptly rather than        !
         !      grinding at the floor. ----------------------------------------------------------------------!
         dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)   ! clamp-ok: the smallest sub-step [s] (above)

         call state_init(y0, n, nsl, y)
         t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik
         err_prev = -1.0_wp                                          ! < 0 => first step uses the I-controller
         dt_warm = dt                                                ! fallback if every step is end-clamped
         do
            if (t >= t_end - tiny_num) exit
            dt_try = dt
            dt = min(dt, t_end - t)
            clamped = dt < dt_try - tiny_num
            call ark2_column_step(y, frozen, col_config, n, nsl, dt, y_new, y_err, ws%stages, niter=np, bf=bfsub,          &
                                  clamp_n=clamp_n)
            call state_sub(y_new, y_err, n, nsl, y_lo)               ! the 1st-order embedded solution
            !----- Grouped WRMS over the WHOLE column state, eq. (3) (state_wrms_grouped). The ARK's           !
            !      theta and plant-water terms are structurally zero here -- both ride operator-split maps     !
            !      outside the ESDIRK tableau -- so they dilute rather than inform. Accepted deliberately:     !
            !      one norm, no per-scheme opt-outs; docs/science/numerical_scheme.md §7, item 3, gives        !
            !      the measured cost. -------------------------------------------------------------------------!
            err = state_wrms_grouped(y_new, y_lo, y, n, nsl, ec%tols)
            !----- ROBUSTNESS: a non-finite err (a stage -- typically the BETA=2.414 stage-3 extrapolation     !
            !      base3 -- overshot the CAS enthalpy into a region where qsat(T) overflows) is a step that    !
            !      is simply TOO BIG: REJECT it and shrink dt deterministically (the NaN poisons the           !
            !      adaptive fac, so use fmin directly). At a smaller dt, Y2 ~ y and base3 no longer            !
            !      overshoots, so the step becomes finite and the march recovers -- the correct adaptive       !
            !      response, not a force-accept. Only if even a floor-sized step is non-finite do we commit    !
            !      and bail, so the driver's has_nan check (meds_polygon) reports it instead of the march      !
            !      hanging.                                                                                    !
            if (err /= err .or. dt /= dt) then
               if (dt <= dt_floor) then
                  call state_init(y_new, n, nsl, y) ; t = t + dt_floor ; nsteps = nsteps + 1_ik
                  exit
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
               !----- WARM-START SEED = the step just ACCEPTED, recorded BEFORE the controller's growth       !
               !      factor is applied. Seeding the next call with the GROWN proposal (dt*fac, fac up to     !
               !      fmax = 5) re-imports the over-estimate that then gets rejected (measured: a 26-29%      !
               !      rejection rate, against under 20% this way). The last accepted size is the one with     !
               !      evidence behind it. --------------------------------------------------------------------!
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
      end associate
   end subroutine adaptive_ark_march
   !=======================================================================================!
   !  INTEG_ARK path: the coupled ESDIRK2 fast step (docs/science/numerical_scheme.md §3 `ark`).        !
   !  Builds the frozen record (build_column_frozen, which also runs the step's water network),         !
   !  packs the state into the column vector, advances one dt_fast with the ARK stepper, commits        !
   !  the network's water, and closes the ledgers. -----------------------------------------------------!
   !=======================================================================================!
   subroutine column_fast_step_ark(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,  &
                                   ws, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, cdiag)
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
      type(fast_workspace_t),  intent(inout) :: ws    !< this thread's reusable step storage
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CO_FAST_ROWS, ncoh) per-cohort diagnostic capture
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters

      type(column_frozen_t)  :: frozen
      type(column_state_t)   :: y, y_out, ycur, ytmp, yerr
      type(surface_state_t)  :: y_stage
      type(surface_tend_t)   :: surf_tend
      type(column_bflux_t)   :: acc, bfsub
      real(wp)    :: tg, fl, dt0, cas_mass_capacity, cas_molar_capacity, enth0, shv0, co20, enth1, shv1, co21, e_soil0, e_soil1, &
           w_soil0, w_soil1
      real(wp)    :: e_pond0, e_pond1   !< [J/m2] pond enthalpy, start and end of the step
      real(wp)    :: w_plant0, w_plant1
      real(wp)    :: w_surface0, dt_warm_next
      !----- Canopy-SURFACE water ledger scratch. -----------------------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, leaf_cap_i, wood_cap_i, intercept_total
      !----- Condensate actually DEPOSITED into soil layer 1 this step: mass [kg/m2] and the paired      !
      !      enthalpy [J/m2]. Held in locals rather than recomputed at the budget calls so the inflow    !
      !      term below is the IDENTICAL number the state update applied ("one flux, both sides"). ------!
      real(wp)    :: cond_dep_mass, cond_dep_enth
      real(wp)    :: tissue_store0, tissue_store1
      real(wp)    :: t_leaf_avg(col_cohort%n), t_wood_avg(col_cohort%n)
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, k, i, isub, nsub, nsteps, nrej
      logical     :: halt_budgets     !< hard-stop on a non-closing budget (debug_error, full column only)

      n = col_cohort%n ; nsl = col_config%soil%n_active
      !----- The clamp counter is ACCUMULATED down the call chain (intent(inout) to the kernels), so    !
      !      this sub-step's tally starts clean here. The ARK clamps only its stage-3 extrapolation     !
      !      base, never y_out; theta_ood_max stays 0 on this path, which has no explicit stage that    !
      !      could carry theta out of its domain. ------------------------------------------------------!
      budget%clamp_stage_n = 0_ik ; budget%theta_ood_max = 0.0_wp

      !----- A reduced column (process_mask_t) freezes a store while its fluxes still act on the          !
      !      neighbours, so it cannot conserve by construction: suppress the HARD stops (the soft         !
      !      n_fail counters still run). -----------------------------------------------------------------!
      halt_budgets = col_config%energy%debug_error .and. mask_is_full(col_config%mask)

      !----- The bottom BCs are pure boundary conditions with no prognostic state behind them (drainage, !
      !      free to sealed, or a head-driven aquifer); the water network takes each as a boundary       !
      !      flow, and the ARK commits its theta. -----------------------------------------------------!
      w_surface0 = biophys%soil_water%w_surface

      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws%spac, &
                               n, nsl, frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      !----- STATE^n pond enthalpy, snapshotted HERE, before anything commits to biophys%soil_water: the    !
      !      pond commit below overwrites w_surface_energy with the network's end-of-step value, so a     !
      !      snapshot taken alongside e_soil0 further down would read the COMMITTED value, the store    !
      !      change would collapse to zero and the whole-column residual would equal -e_pond0.          !
      !      build_column_frozen mutates biophys%snow and nothing in biophys%soil_water, so state^n is      !
      !      intact here. ------------------------------------------------------------------------------!
      e_pond0 = biophys%soil_water%w_surface_energy

      !----- advance one dt_fast: adaptive (embedded-error) or GPU-warp-uniform fixed substeps. ----!
      if (col_config%integrator%adaptive) then
         !----- WARM START (docs/science/numerical_scheme.md §3 `ark`): seed from the step the             !
         !      controller settled on in the PREVIOUS dt_fast for this patch. Column stiffness barely      !
         !      changes from call to call, so a cold start would rediscover the same step size every       !
         !      call and pay about one rejection for it. An explicit [fast].ark_dt_init wins (an           !
         !      operator override); with neither, a patch's first call starts from dt_fast. ---------------!
         dt0 = dt_fast
         if (biophys%adapt_dt_last > tiny_num) dt0 = min(biophys%adapt_dt_last, dt_fast)
         if (col_config%integrator%dt_init > tiny_num) dt0 = min(col_config%integrator%dt_init, dt_fast)
         !----- The error-control bundle, built once per run by build_error_control (meds_fast_config):    !
         !      every tolerance group's setting ([fast].rtol, atol_scale and water_atol), plus the         !
         !      controller and its strictness (docs/science/numerical_scheme.md eqs. (3)-(4)). ------------!
         ec = col_config%integrator%error_control
         call adaptive_ark_march(y, frozen, col_config, n, nsl, dt_fast, ec, dt0, y_out, ws, nsteps, nrej,         &
                                 niter=merge(NEWT_COUPLED, 1_ik, col_config%integrator%coupled_newton), acc=acc, &
                                 dt_warm_out=dt_warm_next,                                          &
                                 clamp_n=budget%clamp_stage_n)
         biophys%adapt_dt_last = dt_warm_next
      else
         nsub = max(1_ik, col_config%integrator%fixed_substeps) ; nrej = 0_ik ; ycur = y ; call bflux_zero(acc, n)
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, frozen, col_config, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr, ws%stages, &
                                  niter=merge(NEWT_COUPLED, 1_ik, col_config%integrator%coupled_newton), bf=bfsub, &
                                  clamp_n=budget%clamp_stage_n)
            call bflux_add(acc, bfsub)
            ycur = ytmp
         end do
         y_out = ycur ; nsteps = nsub
      end if
      !----- WORK counters: what the march cost (the water network's own work is recorded where it     !
      !      ran, in build_column_frozen). ------------------------------------------------------------!
      budget%integ_nsteps = nsteps ; budget%integ_nrej = nrej

      !----- SOIL WATER is the water network's: the ESDIRK stages moved theta at its steady rate; commit !
      !      its end-of-step theta HERE, once, so a single consistent theta feeds the state commit, the    !
      !      soil_temp read-off, and BOTH the soil_water and whole_water storage terms (w_soil1 below). ---!
      y_out%theta(1:nsl) = frozen%hydrology%theta1(1:nsl)
      !----- pond: same treatment as theta -- the network's end-of-step pond, committed THROUGH the     !
      !      state vector so there is one authority. ------------------------------------------------!
      y_out%w_surface      = frozen%hydrology%w_surface1
      y_out%w_surface_energy = frozen%hydrology%w_surface_energy1     ! paired with the mass

      !----- PROCESS MASK (process_mask_t, meds_fast_types). The mask must mean the same thing under     !
      !      every scheme, so it is applied at the ARK's single state-commit point: a masked-off         !
      !      component is restored to state^n (y), leaving the ODE one dimension smaller while its       !
      !      couplings still acted during the march. mask%veg_energy is applied below, where the         !
      !      tissue temperatures are committed: they are not part of column_state_t. --------------------!
      call apply_process_mask(col_config%mask, y, y_out, n, nsl)

      !----- Canopy-SURFACE water: capacity clamp and overflow/deficit bookkeeping, at the store's         !
      !      boundary after the march rather than inside it. DEFICIT (the floor, symmetric to              !
      !      overflow): film_evap uses a state^n-frozen conductance (rescaled by availability in           !
      !      build_column_frozen, but only approximately), so the store can transiently overdraw below     !
      !      0. Left negative, it would corrupt the NEXT dt_fast's frozen interception rate                !
      !      (intercept_canopy_layer's own floor on a negative starting bucket creates water). Floor       !
      !      at 0 and book the shortfall as a NEGATIVE addition to the outflow (the mirror of              !
      !      surf_overflow's sign). Gated behind canopy_water_on, so with the store off no film            !
      !      arithmetic runs. -----------------------------------------------------------------------------!
      surf_overflow = 0.0_wp ; surf_deficit = 0.0_wp
      if (col_config%canopy_water_on) then
         call clamp_canopy_film(y_out, col_cohort%lai, col_cohort%wai,                                  &
                                col_config%leaf_surf_water_max(col_cohort%pft(1:n)),                    &
                                col_config%wood_surf_water_max(col_cohort%pft(1:n)), n, surf_overflow, surf_deficit)
      end if

      !----- unpack into biophys + re-derive the diagnostic soil temperatures + leaf temperatures. -----!
      call unpack_column_state(y_out, n, nsl, biophys)
      !----- The pond is part of the SOIL-WATER store, so it obeys the same soil_water mask as theta.     !
      !      It is the only write to biophys%soil_water besides theta, so skipping it leaves the store at     !
      !      state^n. ------------------------------------------------------------------------------------!
      if (col_config%mask%soil_water) then
         biophys%soil_water%w_surface      = y_out%w_surface
         biophys%soil_water%w_surface_energy = y_out%w_surface_energy
      end if
      call diagnose_soil_temps(y_out, col_config%soil_thermal%soil_dry_heat_capacity, nsl, biophys%soil_energy%soil_temp, &
            biophys%soil_energy%soil_fliq)
      call internal_energy_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        col_config%soil_thermal%soil_dry_heat_capacity(1), tg, fl)
      y_stage%cas_enthalpy = y_out%cas_enthalpy ; y_stage%cas_shv = y_out%cas_shv ; y_stage%cas_co2 = y_out%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, tg, n, surf_tend)
      !----- The tissue heat store is already inside the CAS solve: surface_derivs relaxes each tissue     !
      !      from frozen%tissue%t_leaf0/t_wood0 against leaf_hcap_per_dt/wood_hcap_per_dt = cap/dt_fast    !
      !      (docs/science/vegetation_energy_dynamics.md §One formulation), so nothing is corrected        !
      !      afterwards. Applying the store as a lump on the already-committed CAS is unstable: tissue     !
      !      and canopy air have comparable heat capacities (cap_wood/cap_cas ~ 0.5) and the tissue        !
      !      relaxes far inside one step, so the lump drives a growing oscillation (measured: 49 K CAS     !
      !      kicks per step, wood diverging to 247 K against a 331 K balance). The veg_energy mask         !
      !      holds the tissue at its entry temperature while the fluxes it drives into the CAS are         !
      !      kept. ----------------------------------------------------------------------------------------!
      !----- Commit the TIME-AVERAGE of where each tissue relaxed to over the march, not the final    !
      !      stage's value. That is the temperature whose heat equals the b-weighted complement of the!
      !      flux the canopy air actually received, so the tissue's energy (commit_tissue_energy) and !
      !      the ledger are the same number and the whole-column energy closes to machine precision.  !
      !      Taking the last stage instead leaves ~3.9e4 J unaccounted per step -- 4e-4 relative,     !
      !      against a 1e-6 tolerance -- because the CAS moves across the stages. --------------------!
      if (allocated(acc%tissue_leaf_int)) then
         t_leaf_avg(1:n) = acc%tissue_leaf_int(1:n) / dt_fast
         t_wood_avg(1:n) = acc%tissue_wood_int(1:n) / dt_fast
      else
         t_leaf_avg(1:n) = surf_tend%leaf_temp(1:n)
         t_wood_avg(1:n) = surf_tend%wood_temp(1:n)
      end if
      call commit_tissue_energy(frozen%tissue, col_cohort%nplant, col_config%mask%veg_energy,                  &
                                col_config%mask%hydraulics, dt_fast, t_leaf_avg, t_wood_avg, n, biophys,      &
                                tissue_store0, tissue_store1)
      !----- The transpiration each cohort's leaves lost over the step: the b-weighted amount over the  !
      !      accepted sub-steps, which the leaf water was debited and the canopy air credited. --------!
      if (present(cdiag) .and. allocated(acc%transp_coh))                                           &
         cdiag(CO_TRANSP, 1:n) = transp_per_leaf_area(acc%transp_coh(1:n), dt_fast, col_cohort%nplant(1:n), &
                                                      col_cohort%leaf_area(1:n))

      !----- WHOLE-COLUMN CONSERVATION LEDGER: the column's seven budgets, from the b-weighted           !
      !      boundary-flux AMOUNTS accumulated over the sub-steps (acc). The flux-form CAS commits and   !
      !      the soil-heat column's energy_resid = 0 make the identity exact, so ENERGY closes to        !
      !      machine precision (including the frozen rain/runoff/drainage advection, a fixed source),    !
      !      and water closes to round-off (docs/science/numerical_scheme.md §5b). acc holds AMOUNTS,    !
      !      not rates. ---------------------------------------------------------------------------------!
      !----- DEPOSIT THE CONDENSATE (the same routing on RK45). Dew and fog land on a surface inside     !
      !      the column, so the condensate goes, mass and enthalpy paired, into soil layer 1 at the CAS  !
      !      temperature it condensed at, and the whole-column ledger closes with no boundary term       !
      !      (docs/science/numerical_scheme.md §7, item 4, says why not onto the foliage). Applied to    !
      !      y_out BEFORE the store snapshot below, and mirrored into biophys (already unpacked).        !
      !      acc%whole_cond is 0 whenever cas_condensation is off. --------------------------------------!
      cond_dep_mass = 0.0_wp ; cond_dep_enth = 0.0_wp
      if (acc%whole_cond > 0.0_wp) then
         cond_dep_mass = acc%whole_cond
         cond_dep_enth = acc%whole_cond_enth       ! b-weighted at the stage CAS temperatures (see column_bflux_t)
         call deposit_condensate(y_out, col_config%soil%dz(1), cond_dep_mass, cond_dep_enth)
         biophys%soil_water%theta(1)       = y_out%theta(1)
         biophys%soil_energy%soil_energy(1) = y_out%soil_energy(1)
      end if

      !----- No separate wood-radiation term: surface_derivs' balance carries the wood's absorbed       !
      !      radiation, so it is already inside coh_rnet / acc%whole_enth_in. --------------------------!

      cas_mass_capacity = frozen%cas%cas_mass_capacity ; cas_molar_capacity = frozen%cas%cas_molar_capacity
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv ; co21 = y_out%cas_co2
      e_pond1 = frozen%hydrology%w_surface_energy1
      e_soil0 = soil_energy_store(y%soil_energy,     col_config%soil%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, col_config%soil%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     col_config%soil%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, col_config%soil%dz, nsl)
      !----- Plant internal water MASS is a store the step evolves: it absorbs the step's mismatch      !
      !      between transpiration and root uptake (docs/science/plant_hydraulics.md eqs. (1)-(2)),     !
      !      so leaving it out would make that real change read as a leak. -----------------------------!
      w_plant0 = plant_water_store(col_cohort%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(col_cohort%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- Canopy-SURFACE water, already per m2 of ground (no nplant factor, unlike w_plant0/1). It      !
      !      is valued at frozen%film%film_liquid_enthalpy = u_liq(t_film_valuation), the liquid           !
      !      enthalpy the intercepted water arrived with; the tissue pays enthalpy_vapor -                 !
      !      film_liquid_enthalpy per kg it evaporates (surface_derivs), so this store closes exactly      !
      !      against the CAS credit. All zero when canopy_water_on is off. --------------------------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))
      !----- [fast].debug_error promotes a non-closing budget from a counted n_fail to a hard           !
      !      `error stop` (docs/science/numerical_scheme.md §5b); off by default. Tolerances are        !
      !      FLUX-scaled (meds_budget_check header): rtol * gross boundary flux over the step plus a    !
      !      rate floor * dt_fast, because a store-scaled tolerance lets a ~1 W/m2 leak through. -------!
      budget%atm_heat_export = acc%atm_heat_out ; budget%atm_vap_export = acc%atm_vap_out
      call budget_check(budget%cas_energy, cas_mass_capacity*enth0, cas_mass_capacity*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'cas_energy (ark)', halt_budgets)
      call budget_check(budget%cas_water,  cas_mass_capacity*shv0,  cas_mass_capacity*shv1,  acc%cas_vap_in,  acc%cas_vap_out, &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'cas_water (ark)', halt_budgets)
      call budget_check(budget%cas_co2,    cas_molar_capacity*co20,  cas_molar_capacity*co21,  acc%cas_co2_in,  acc%cas_co2_out, &
                        dt_fast, BUDGET_CO2_RATE_FLOOR, 'cas_co2 (ark)', halt_budgets)
      !----- cond_dep_enth is a boundary INPUT to the SOIL store, though the whole-column ledger needs     !
      !      no term for it. The deposit moves condensate CAS -> soil layer 1 AFTER the march: the         !
      !      whole-column ledger sees a transfer between two stores it tracks, which telescopes, but       !
      !      this soil budget sees only the soil, for which the same transfer is an inflow. Leaving it     !
      !      out makes the soil budget short by exactly the deposit (and the soil_water twin by the        !
      !      paired mass; the two residuals' ratio is internal_energy_liquid(T_cas)). RK45 keeps no        !
      !      per-kernel soil budgets: every store rides the same column_derivs RHS, so its                 !
      !      whole-column ledger is the per-store one. ----------------------------------------------------!
      call budget_check(budget%soil_energy, e_soil0, e_soil1, acc%soil_enth_in + cond_dep_enth,            &
                        acc%soil_enth_out, dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'soil_energy (ark)', halt_budgets)
      !----- SOIL WATER: storage theta^n -> theta1 (both the water network's), in through the top face   !
      !      net of evaporation and exfiltration (q_top), out through the bottom and to the roots -- the   !
      !      network's own flows, which moved its water exactly -- plus the condensate deposit (see the    !
      !      soil_energy note just above). ---------------------------------------------------------!
      call budget_check(budget%soil_water,  w_soil0, w_soil1,                                              &
                        frozen%hydrology%q_top*rho_h2o*dt_fast + cond_dep_mass,                            &
                        (frozen%hydrology%drainage + frozen%roots%uptake)*dt_fast,                         &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'soil_water (ark)', halt_budgets)
      !----- whole-WATER: rainfall IN; drainage + runoff + CAS-vapour OUT; soil, ponding, plant internal  !
      !      water and the films in the store. The soil, pond, plant, drainage and runoff terms are the     !
      !      water network's step; the CAS-vapour exchange g_atm_vapour*(shv-shv_atm) is the ARK-          !
      !      accumulated part (acc%whole_wat_out). -------------------------------------------------!
      !----- SNOW: the pack is a mass + energy STORE, its accumulated rainfall enthalpy a boundary       !
      !      INPUT, and SNOWFALL a boundary water input (it lands in the pack, or, below the new-snow    !
      !      threshold, in the throughfall). Sublimation leaves via the CAS vapour term and meltwater    !
      !      moves pack -> pond as a paired transfer, so neither needs a term. All zero without snow. ---!
      call budget_check(budget%whole_water,                                                                &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0 + frozen%snow%swe0,  &
                        w_soil1 + cas_mass_capacity*shv1 + frozen%hydrology%w_surface1 + w_plant1 + surf_water1 &
                        + frozen%snow%swe1, &
                        acc%whole_wat_in + (forc%rainfall + forc%snowfall + biophys%shed_water_rate)*dt_fast,       &
                        acc%whole_wat_out                                                                   &
                        + (frozen%hydrology%runoff_surf + frozen%hydrology%drainage)*dt_fast                &
                                          + surf_overflow - surf_deficit,                                 &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'whole_water (ark)', halt_budgets)
      call budget_check(budget%whole_energy,                                                               &
                             !----- The pack hands its meltwater, mass and enthalpy, to the POND:           !
                             !      both are tracked stores, so the pair telescopes here. ------------------!
                        e_soil0                           + cas_mass_capacity*enth0 + surf_enth0                     &
                        + frozen%snow%enth0 + e_pond0 + tissue_store0,                               &
                        e_soil1 + cas_mass_capacity*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond1                &
                        + tissue_store1,                                                                 &
                        acc%whole_enth_in + intercept_total*dt_fast*internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                          + frozen%snow%acc_enth                                       &
                                          + frozen%hydrology%ground_inflow_enth*dt_fast,                              &
                        acc%whole_enth_out + (surf_overflow - surf_deficit)                                 &
                        * internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                           + frozen%hydrology%runoff_enth*dt_fast,                                    &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'whole_energy (ark)', halt_budgets)

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_ark

end module meds_fast_ark
