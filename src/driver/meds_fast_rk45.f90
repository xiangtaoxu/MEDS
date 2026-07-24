!==========================================================================================!
! meds_fast_rk45 -- the ED2-faithful adaptive Cash-Karp RK45 fast-loop integrator            !
! (docs/dev_plans/MEDS_ED2_RK45_DESIGN.md, phase P2). Peer of meds_fast_ark: dispatched by     !
! meds_fast_split%column_fast_step when cfg%time_integrator == INTEG_RK4. Unlike the ARK's      !
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
   use meds_constants,        only : tiny_num, rho_h2o
   use meds_therm_lib,           only : cas_temp_of_enthalpy, internal_energy_liquid, uext_to_temp
   use meds_fast_time_derivs, only : surface_derivs, column_derivs
   use meds_fast_types,       only : column_state_t, column_frozen_t, column_tend_t,             &
                                     surface_state_t, surface_frozen_t, surface_tend_t,          &
                                     column_config_t, column_cohort_t, column_forcing_t,         &
                                     column_budget_t, mask_is_full,                              &
                                     WOODEN_DIAGNOSTIC
   use meds_fast_ark,         only : state_init, state_axpy, state_accum, state_sub,             &
                                     build_column_frozen
   use meds_fast_control,     only : error_control_t, build_error_control, state_wrms_grouped,   &
                                     step_control_factor
   use meds_config,           only : meds_config_t, CTRL_L2_STRICT
   use meds_biophysics_types, only : aero_env_t, aero_geom_t, aero_out_t, patch_biophys_t
   use meds_budget_check,     only : budget_accumulate, budget_check_stop
   implicit none
   private

   public :: rk45_column_step, adaptive_rk45_march, column_fast_step_rk45

   !----- Cash-Karp embedded 5(4) tableau (Cash & Karp 1990, ACM TOMS 16:201; the SAME          !
   !      coefficients as Numerical Recipes' rkck). c_i (stage times) are documentation only --   !
   !      the frozen forcing fro has no explicit time dependence, so they never enter the code.   !
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

contains

   !---------------------------------------------------------------------------------------!
   ! stage_bnd -- the whole-column BOUNDARY-flux rates at one stage's state (analogous to           !
   ! column_be_stage's stage_bflux_t in meds_fast_ark, but for a fully-explicit stage: no BE          !
   ! solve, so the "state" IS just the stage's own input ys, and sf is column_derivs' own surface     !
   ! diagnostic at ys). cond_enth is pre-multiplied by u_liq at THIS stage's own tcas (mirrors ARK's   !
   ! t_cas1 reference, sec 3.4/9's "one flux, both sides" -- the CAS's own reference, not a frozen     !
   ! one, since condensation genuinely happens at the evolving CAS temperature). ---------------------!
   pure subroutine stage_bnd(ys, fro, sf, rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i)
      type(column_state_t),  intent(in)  :: ys
      type(column_frozen_t), intent(in)  :: fro
      type(surface_tend_t),  intent(in)  :: sf
      real(wp),               intent(out) :: rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i
      real(wp) :: tcas_i
      tcas_i      = cas_temp_of_enthalpy(ys%cas_enthalpy, ys%cas_shv)
      rnet_i      = sf%coh_rnet
      atm_enth_i  = fro%surf%gah * (ys%cas_enthalpy - fro%surf%enth_atm)
      atm_vap_i   = fro%surf%gaw * (ys%cas_shv      - fro%surf%shv_atm)
      cond_i      = sf%cond
      cond_enth_i = sf%cond * internal_energy_liquid(tcas_i)
   end subroutine stage_bnd

   !---------------------------------------------------------------------------------------!
   ! rk45_column_step -- one Cash-Karp 5(4) step of size dt from y, over the pure RHS            !
   ! column_derivs. Commits the 5TH-order solution (local extrapolation, matching how              !
   ! adaptive_ark_march/adaptive_imex_march both commit their higher-order result); y_err is the    !
   ! (5th - 4th) embedded difference for the adaptive controller. w_out/e_in/e_out are the           !
   ! whole-column boundary-flux AMOUNTS over dt, b-weighted by the SAME 5th-order b-vector as the    !
   ! state commit (the consistent quadrature for a boundary integral over this step). e_in's          !
   ! infiltration term (fro%infiltration*u_liq(rain_temp)) IS the whole-column precip-energy          !
   ! input -- the caller must NOT also add a separate forc%precip term on top (double-counts nearly   !
   ! the full infiltrating share whenever infiltration ~= precip); mirrors ARK's own bf%whole_enth_in, !
   ! which folds e_infil in the same way with no further outer addition. -----------------------------!
   pure subroutine rk45_column_step(y, fro, n, nsl, dt, y_out, y_err, w_out, e_in, e_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      real(wp),               intent(out) :: w_out, e_in, e_out

      type(column_tend_t)  :: k1, k2, k3, k4, k5, k6
      type(column_state_t) :: ys, y_4th
      type(surface_tend_t) :: sf
      real(wp) :: rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6)
      real(wp) :: bw_rnet, bw_atm_enth, bw_atm_vap, bw_cond, bw_cond_enth

      call column_derivs(y, fro, n, nsl, k1, sf_out=sf)
      call stage_bnd(y, fro, sf, rnet(1), atm_enth(1), atm_vap(1), cond(1), cond_enth(1))

      call state_init(y, n, nsl, ys) ; call state_accum(ys, dt*A21, k1, n, nsl)
      call column_derivs(ys, fro, n, nsl, k2, sf_out=sf)
      call stage_bnd(ys, fro, sf, rnet(2), atm_enth(2), atm_vap(2), cond(2), cond_enth(2))

      call state_init(y, n, nsl, ys)
      call state_accum(ys, dt*A31, k1, n, nsl) ; call state_accum(ys, dt*A32, k2, n, nsl)
      call column_derivs(ys, fro, n, nsl, k3, sf_out=sf)
      call stage_bnd(ys, fro, sf, rnet(3), atm_enth(3), atm_vap(3), cond(3), cond_enth(3))

      call state_init(y, n, nsl, ys)
      call state_accum(ys, dt*A41, k1, n, nsl) ; call state_accum(ys, dt*A42, k2, n, nsl)
      call state_accum(ys, dt*A43, k3, n, nsl)
      call column_derivs(ys, fro, n, nsl, k4, sf_out=sf)
      call stage_bnd(ys, fro, sf, rnet(4), atm_enth(4), atm_vap(4), cond(4), cond_enth(4))

      call state_init(y, n, nsl, ys)
      call state_accum(ys, dt*A51, k1, n, nsl) ; call state_accum(ys, dt*A52, k2, n, nsl)
      call state_accum(ys, dt*A53, k3, n, nsl) ; call state_accum(ys, dt*A54, k4, n, nsl)
      call column_derivs(ys, fro, n, nsl, k5, sf_out=sf)
      call stage_bnd(ys, fro, sf, rnet(5), atm_enth(5), atm_vap(5), cond(5), cond_enth(5))

      call state_init(y, n, nsl, ys)
      call state_accum(ys, dt*A61, k1, n, nsl) ; call state_accum(ys, dt*A62, k2, n, nsl)
      call state_accum(ys, dt*A63, k3, n, nsl) ; call state_accum(ys, dt*A64, k4, n, nsl)
      call state_accum(ys, dt*A65, k5, n, nsl)
      call column_derivs(ys, fro, n, nsl, k6, sf_out=sf)
      call stage_bnd(ys, fro, sf, rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6))

      !----- y_out = y + dt*(B1*k1 + B3*k3 + B4*k4 + B6*k6)  [5th order; b2=b5=0]. -----------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt*B1, k1, n, nsl)
      call state_accum(y_out, dt*B3, k3, n, nsl)
      call state_accum(y_out, dt*B4, k4, n, nsl)
      call state_accum(y_out, dt*B6, k6, n, nsl)

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
      bw_rnet      = B1*rnet(1)      + B3*rnet(3)      + B4*rnet(4)      + B6*rnet(6)
      bw_atm_enth  = B1*atm_enth(1)  + B3*atm_enth(3)  + B4*atm_enth(4)  + B6*atm_enth(6)
      bw_atm_vap   = B1*atm_vap(1)   + B3*atm_vap(3)   + B4*atm_vap(4)   + B6*atm_vap(6)
      bw_cond      = B1*cond(1)      + B3*cond(3)      + B4*cond(4)      + B6*cond(6)
      bw_cond_enth = B1*cond_enth(1) + B3*cond_enth(3) + B4*cond_enth(4) + B6*cond_enth(6)

      e_in  = (bw_rnet + fro%surf%abs_sw_ground + fro%surf%abs_lw_ground) * dt                    &
              + fro%infiltration * dt * internal_energy_liquid(fro%rain_temp)
      e_out = bw_atm_enth * dt + bw_cond_enth * dt                                                &
              + fro%drainage * dt * internal_energy_liquid(fro%t_bot)                             &
              + fro%runoff_surf * dt * internal_energy_liquid(fro%surf%t_ground)
      w_out = bw_atm_vap * dt + bw_cond * dt + fro%drainage * dt + fro%runoff_surf * dt
   end subroutine rk45_column_step

   !---------------------------------------------------------------------------------------!
   ! adaptive_rk45_march -- integrate to t_end with the Cash-Karp embedded error driving the        !
   ! step controller (reusing meds_fast_control's WRMS + I/PI controller + warm start, sec 6 --      !
   ! "no new controller"), mirroring adaptive_ark_march's accept/reject structure exactly, just       !
   ! over rk45_column_step instead of ark2_column_step. with_mass=.true. (the default): RK45 has no   !
   ! operator split at all, so mass genuinely differs between the 5th/4th solutions -- a live error    !
   ! signal, unlike ARK where it is structurally zero. -----------------------------------------------!
   subroutine adaptive_rk45_march(y0, fro, n, nsl, t_end, ec, dt_init, y_out, nsteps, nrej,       &
                                  w_out_acc, e_in_acc, e_out_acc, dt_warm_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec
      type(column_state_t),  intent(out) :: y_out
      integer(ik),            intent(out) :: nsteps, nrej
      real(wp),               intent(out) :: w_out_acc, e_in_acc, e_out_acc
      real(wp),    optional,  intent(out) :: dt_warm_out

      type(column_state_t) :: y, y_new, y_err, y_zero
      real(wp) :: t, dt, err, err_prev, fac, dt_floor
      real(wp) :: w_out, e_in, e_out, dt_try, dt_warm
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
         call rk45_column_step(y, fro, n, nsl, dt, y_new, y_err, w_out, e_in, e_out)
         !----- named temporary: never pass a derived-type-valued function result straight into a  !
         !      call (the nvfortran whole-program-optimizer trap documented in CLAUDE.md). --------!
         y_zero = zero_like(y_err, n, nsl)
         err = state_wrms_grouped(y_err, y_zero, y_new, n, nsl, ec%tols, with_mass=.true.)
         if (err /= err .or. dt /= dt) then    ! non-finite: too big a step, reject deterministically
            if (dt <= dt_floor) then
               call state_init(y_new, n, nsl, y) ; t = t + dt_floor ; nsteps = nsteps + 1_ik ; exit
            end if
            nrej = nrej + 1_ik ; dt = max(dt * ec%fmin, dt_floor) ; cycle
         end if
         fac = step_control_factor(err, err_prev, ec)
         if (err <= 1.0_wp .or. dt <= dt_floor) then
            if (ec%level == CTRL_L2_STRICT .and. err > 1.0_wp) &
               error stop 'adaptive_rk45_march: L2 strict -- floor step cannot meet tolerance'
            call state_init(y_new, n, nsl, y)
            w_out_acc = w_out_acc + w_out ; e_in_acc = e_in_acc + e_in ; e_out_acc = e_out_acc + e_out
            t = t + dt ; nsteps = nsteps + 1_ik
            err_prev = err
            if (.not. clamped) dt_warm = dt
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         if (nsteps + nrej > 4096_ik) exit   ! hard backstop (should never trigger)
      end do
      call state_init(y, n, nsl, y_out)
      if (present(dt_warm_out)) dt_warm_out = dt_warm
   end subroutine adaptive_rk45_march

   !----- a column_state_t of the SAME shape as `ref`, every field zeroed -- lets state_wrms_grouped   !
   !      (which takes two STATES to difference) read a already-a-difference y_err directly, without    !
   !      a bespoke "WRMS of one state" variant. Trivial and allocation-only; not a hot path (once      !
   !      per accept/reject trial, not per stage). ------------------------------------------------------!
   pure function zero_like(ref, n, nsl) result(z)
      type(column_state_t), intent(in) :: ref
      integer(ik),           intent(in) :: n, nsl
      type(column_state_t) :: z
      z%cas_enthalpy = 0.0_wp ; z%cas_shv = 0.0_wp ; z%cas_co2 = 0.0_wp
      z%soil_energy = 0.0_wp ; z%theta = 0.0_wp
      allocate(z%leaf_water_mass(n), z%wood_water_mass(n))
      z%leaf_water_mass = 0.0_wp ; z%wood_water_mass = 0.0_wp
   end function zero_like

   !=======================================================================================!
   !  INTEG_RK4 path: the ED2-faithful adaptive Cash-Karp fast step. Shares build_column_frozen  !
   !  with the ARK path (meds_fast_ark) -- the Act-1 pre-pass does not depend on which stepper      !
   !  advances the macro-step. -----------------------------------------------------------------!
   !=======================================================================================!
   subroutine column_fast_step_rk45(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,  &
                                    gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budg
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters

      type(column_frozen_t)  :: fro
      type(column_state_t)   :: y, y_out
      type(surface_state_t)  :: ys
      type(surface_frozen_t) :: fs
      type(surface_tend_t)   :: sf
      real(wp)    :: dt0, wcap, enth0, shv0, enth1, shv1
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, w_plant0, w_plant1
      real(wp)    :: w_out_acc, e_in_acc, e_out_acc, w_in, w_out, e_in, e_out
      real(wp)    :: tg, fl, dt_warm_next
      !----- Canopy-SURFACE water (sec 3.4, P2c) ledger scratch. --------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, leaf_cap_i, wood_cap_i, intercept_total
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, k, i, nsteps, nrej
      logical     :: halt_budgets

      n = coh%n ; nsl = ccfg%soil%n_active
      halt_budgets = ccfg%energy%debug_error .and. mask_is_full(ccfg%mask)

      call build_column_frozen(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, n, nsl, &
                               fro, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)

      dt0 = dt_fast
      if (bio%adapt_dt_last > tiny_num) dt0 = min(bio%adapt_dt_last, dt_fast)
      if (cfg%ark_dt_init  > tiny_num)  dt0 = min(cfg%ark_dt_init,   dt_fast)
      ec = build_error_control(cfg)
      ec%p_order = RK45_P_ORDER
      call adaptive_rk45_march(y, fro, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                              w_out_acc, e_in_acc, e_out_acc, dt_warm_out=dt_warm_next)
      bio%adapt_dt_last = dt_warm_next
      budg%integ_nsteps = nsteps ; budg%integ_nrej = nrej

      !----- §5.1 PROCESS MASK: a masked-off component is restored to state^n (y). RK45 genuinely   !
      !      integrates soil water AND plant mass (unlike ARK, where soil water is fully operator-    !
      !      split and mass is separately split too), so BOTH honour the mask here directly. -----------!
      if (.not. ccfg%mask%cas_energy) y_out%cas_enthalpy       = y%cas_enthalpy
      if (.not. ccfg%mask%cas_vapour) y_out%cas_shv            = y%cas_shv
      if (.not. ccfg%mask%cas_co2)    y_out%cas_co2            = y%cas_co2
      if (.not. ccfg%mask%soil_heat)  y_out%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      if (.not. ccfg%mask%soil_water) y_out%theta(1:nsl)       = y%theta(1:nsl)
      if (.not. ccfg%mask%hydraulics) then
         y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
         y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)
         !----- Canopy-SURFACE water (sec 3.4, P2c) rides the SAME hydraulics mask entry as internal   !
         !      water mass -- both are "plant water" stores; a dedicated mask field is deferred (no      !
         !      test scenario needs surf_water reduced independently of internal mass yet). -------------!
         y_out%leaf_surf_water(1:n) = y%leaf_surf_water(1:n)
         y_out%wood_surf_water(1:n) = y%wood_surf_water(1:n)
      end if

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
      if (ccfg%canopy_water_on) then
         do i = 1_ik, n
            leaf_cap_i = ccfg%hydro%dewmx * coh%lai(i) ; wood_cap_i = ccfg%hydro%dewmx * coh%wai(i)
            surf_overflow = surf_overflow + max(0.0_wp, y_out%leaf_surf_water(i) - leaf_cap_i)          &
                                           + max(0.0_wp, y_out%wood_surf_water(i) - wood_cap_i)
            surf_deficit  = surf_deficit  + max(0.0_wp, -y_out%leaf_surf_water(i))                      &
                                           + max(0.0_wp, -y_out%wood_surf_water(i))
            y_out%leaf_surf_water(i) = min(max(y_out%leaf_surf_water(i), 0.0_wp), leaf_cap_i)
            y_out%wood_surf_water(i) = min(max(y_out%wood_surf_water(i), 0.0_wp), wood_cap_i)
         end do
      end if

      !----- unpack into bio + re-derive the diagnostic soil/leaf/wood temperatures. -----------!
      bio%cas%can_enthalpy = y_out%cas_enthalpy ; bio%cas%can_shv = y_out%cas_shv ; bio%cas%can_co2 = y_out%cas_co2
      bio%cas%can_temp = cas_temp_of_enthalpy(y_out%cas_enthalpy, y_out%cas_shv)
      bio%soil_e%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      bio%soil_w%theta(1:nsl)       = y_out%theta(1:nsl)
      bio%leaf_water_mass(1:n) = y_out%leaf_water_mass(1:n)
      bio%wood_water_mass(1:n) = y_out%wood_water_mass(1:n)
      bio%leaf_surf_water(1:n) = y_out%leaf_surf_water(1:n)
      bio%wood_surf_water(1:n) = y_out%wood_surf_water(1:n)
      do k = 1_ik, nsl
         call uext_to_temp(y_out%soil_energy(k), y_out%theta(k)*rho_h2o,                          &
                           ccfg%soil_thermal%soil_dry_heat_capacity(k), bio%soil_e%soil_temp(k), bio%soil_e%soil_fliq(k))
      end do
      call uext_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        ccfg%soil_thermal%soil_dry_heat_capacity(1), tg, fl)
      fs = fro%surf ; fs%t_ground = tg
      ys%cas_enthalpy = y_out%cas_enthalpy ; ys%cas_shv = y_out%cas_shv ; ys%cas_co2 = y_out%cas_co2
      call surface_derivs(ys, fs, n, sf)
      bio%leaf_temp(1:n) = sf%leaf_temp(1:n)
      if (ccfg%wood_energy_model == WOODEN_DIAGNOSTIC) bio%wood_temp(1:n) = sf%wood_temp(1:n)

      !----- WHOLE-COLUMN CONSERVATION LEDGER (design doc sec 8 gates 2/3 -- the headline           !
      !      deliverable): unlike ARK's per-kernel + whole ledger, RK45 checks the WHOLE-COLUMN        !
      !      water/energy budgets only (soil water and plant mass are genuinely integrated states,      !
      !      not operator-split, so there is no separate "soil_water (rk45)"/frozen-flux kernel          !
      !      check the way ARK needs one -- the whole-column ledger IS the individual-store ledger      !
      !      here, since every store advances through the SAME column_derivs RHS). -------------------!
      wcap = fro%surf%wcap
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv
      e_soil0 = 0.0_wp ; e_soil1 = 0.0_wp ; w_soil0 = 0.0_wp ; w_soil1 = 0.0_wp
      do k = 1_ik, nsl
         e_soil0 = e_soil0 + y%soil_energy(k)     * ccfg%soil%dz(k)
         e_soil1 = e_soil1 + y_out%soil_energy(k) * ccfg%soil%dz(k)
         w_soil0 = w_soil0 + y%theta(k)     * ccfg%soil%dz(k) * rho_h2o
         w_soil1 = w_soil1 + y_out%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      w_plant0 = sum(coh%nplant(1:n) * (y%leaf_water_mass(1:n)     + y%wood_water_mass(1:n)))
      w_plant1 = sum(coh%nplant(1:n) * (y_out%leaf_water_mass(1:n) + y_out%wood_water_mass(1:n)))
      !----- Canopy-SURFACE water (sec 3.4, P2c): already ground-area-referenced (no nplant factor,     !
      !      unlike w_plant0/1 above). Valued at the SAME fixed rain_temp reference the split path's       !
      !      own surf_enth0/1 uses (KNOWN DEFERRED IMPRECISION, mirrors the P1/P0 root_heat_sink notes:      !
      !      a film_evap*(enthalpy_vapor(tl)-u_liq(rain_temp)) mismatch survives, since film evaporation      !
      !      is credited to the CAS at the LEAF temperature but the store is booked at rain_temp -- see        !
      !      the looser whole_energy tolerance this same imprecision earns on the split path). All zero        !
      !      when canopy_water_on is off, so this is a no-op on the byte-identical default path. --------------!
      surf_water0 = sum(y%leaf_surf_water(1:n)     + y%wood_surf_water(1:n))
      surf_water1 = sum(y_out%leaf_surf_water(1:n) + y_out%wood_surf_water(1:n))
      surf_enth0  = surf_water0 * internal_energy_liquid(fro%rain_temp)
      surf_enth1  = surf_water1 * internal_energy_liquid(fro%rain_temp)
      intercept_total = sum(fro%intercept_leaf(1:n) + fro%intercept_wood(1:n))

      !----- e_in is e_in_acc ALONE -- NOT e_in_acc + a separate forc%precip energy term. Precip's       !
      !      energy already enters the ledger via rk45_column_step's OWN per-substep e_infil            !
      !      (fro%infiltration*u_liq(rain_temp), b-weighted into e_in_acc), which is the SAME frozen      !
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
      w_in  = (forc%precip + forc%shed_water_rate) * dt_fast   ! P4: shed water is a boundary input too;
                                                                ! its energy needs NO separate term here,
                                                                ! for the SAME reason precip's doesn't --
                                                                ! it rides e_in_acc via rk45_column_step's
                                                                ! own e_infil, once mixed into hforc%precip_ground
                                                                ! (build_column_frozen, shared with ARK).
      w_out = w_out_acc + surf_overflow - surf_deficit
      e_in  = e_in_acc + intercept_total * dt_fast * internal_energy_liquid(fro%rain_temp)
      e_out = e_out_acc + (surf_overflow - surf_deficit) * internal_energy_liquid(fro%rain_temp)

      call budget_accumulate(budg%whole_water, w_soil0 + wcap*shv0 + w_plant0 + surf_water0,      &
                             w_soil1 + wcap*shv1 + w_plant1 + surf_water1, w_in, w_out, 1.0_wp,     &
                             max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
      call budget_check_stop(budg%whole_water%resid, max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp, &
                             1.0e-4_wp, 'whole_water (rk45)', halt_budgets)
      call budget_accumulate(budg%whole_energy, e_soil0 + wcap*enth0 + surf_enth0,                &
                             e_soil1 + wcap*enth1 + surf_enth1,                                    &
                             e_in, e_out, 1.0_wp, abs(e_soil1 + wcap*enth1), 1.0e-6_wp,             &
                             merge(5.0e6_wp, 1.0e0_wp, ccfg%canopy_water_on))
      call budget_check_stop(budg%whole_energy%resid, abs(e_soil1 + wcap*enth1), 1.0e-6_wp,       &
                             merge(5.0e6_wp, 1.0e0_wp, ccfg%canopy_water_on), 'whole_energy (rk45)', halt_budgets)
      !----- NOT YET CHECKED: a per-kernel cas_co2 closure (ARK's own budg%cas_co2 check) would need  !
      !      a b-weighted per-stage CO2 atmospheric-exchange accumulation this first pass does not      !
      !      track (only rnet/atm_enth/atm_vap/cond are tracked in rk45_column_step) -- deferred; the   !
      !      whole-column water/energy ledger above are the design doc's headline gates (sec 8          !
      !      gates 2/3), not the per-kernel CAS/soil individual checks ARK additionally provides. -------!

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_rk45

end module meds_fast_rk45
