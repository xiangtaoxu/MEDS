!==========================================================================================!
! meds_fast_ark -- the IMEX-ARK fast-loop scheme backend: the peer of meds_fast_split. Hosts        !
! the ARK dispatch (column_fast_step_ark, called from meds_fast_split%column_fast_step when          !
! cfg%time_integrator=="ark") + its frozen pre-pass (build_column_frozen) + the production ARS(2,2,2) !
! time-integrator machinery (design docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md): the L-stable ESDIRK      !
! stage solve (column_be_stage + the 2x2 leaf<->CAS Newton arrowhead newton_surface_solve/jac_surface), !
! the ark2_column_step 2nd-order step + its embedded-error adaptive controller adaptive_ark_march,      !
! and the shared state-vector building blocks (state_init/state_wrms/state_extrap/state_sub/            !
! state_err_diff/clamp_cas/clamp_theta) + the operator-split plant-hydraulics advance                    !
! (advance_hydraulics_full) + the boundary-flux conservation ledger (bflux_*).                            !
!                                                                                          !
! column_be_stage/advance_hydraulics_full/state_init/state_wrms are ALSO the shared building blocks       !
! the test-only RK4/IMEX-Euler oracle (meds_fast_rk4_oracle) calls cross-module -- they are PUBLIC here   !
! (previously private, since their only caller lived in the same file before this split) purely as a      !
! consequence of the file separation, not a behaviour change.                                              !
!                                                                                          !
! INTEG_ARK path (column_fast_step_ark): shares the split's frozen pre-pass (build_column_frozen),        !
! packs the state into the pure column vector, advances one dt_fast with the ARK stepper, then unpacks.    !
! PARTIAL precip>0 guard-lift: the ARK now carries the split's soil-boundary water-enthalpy advection      !
! (rain/runoff/drainage liquid enthalpy, in column_be_stage) and persists the scratch hydrology's          !
! ponding/aquifer/water-table (column_state_t still doesn't advance them prognostically -> a lagged        !
! operator split, so the whole-WATER budget closes only to the split-error tolerance, not machine).        !
! STILL restricted to free-drain + no Zeng-Decker: those bottom BCs need prognostic aquifer/z_wt in the    !
! state vector.                                                                                             !
!==========================================================================================!
module meds_fast_ark
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, latent_heat_vap, rho_h2o, r_gas, pi, &
                                     tsupercool_liq, grav_head
   use meds_plant_hydraulics, only : rhizosphere_cond, solve_plant_water
   use meds_hydr_lib, only : soil_hydr_cond_from_theta
   use meds_config,           only : meds_config_t, hydraulics_config_t,                          &
                                     SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED,               &
                                     INTEG_SPLIT, INTEG_ARK, CTRL_L2_STRICT
   use meds_fast_control,     only : error_control_t, build_error_control, state_wrms_grouped,   &
                                     step_control_factor
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,        &
                                     energy_forcing_t, energy_opts_t, energy_flux_t,           &
                                     soil_column_t, soil_energy_column_t, chydro_forcing_t, chydro_flux_t, &
                                     leaf_energy_env_t, leaf_energy_flux_t, SOIL_BC_FREE_DRAIN, &
                                     snow_params_t, snow_env_t, snow_flux_t, snow_melt_t
   use meds_fast_time_derivs, only : surface_derivs, root_weighted_psi
   use meds_fast_types,       only : column_config_t, column_cohort_t, column_forcing_t,       &
                                     column_budget_t, alloc_column_cohort,                      &
                                     LEAFEN_DIAGNOSTIC, LEAFEN_PROGNOSTIC,                       &
                                     WOODEN_DIAGNOSTIC, WOODEN_PROGNOSTIC,                       &
                                     SOILH2O_LAGGED, SOILH2O_COUPLED,                            &
                                     column_state_t, column_frozen_t, surface_state_t,          &
                                     surface_frozen_t, surface_tend_t, stage_bflux_t, column_bflux_t, &
                                     mask_is_full
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_soil_energy,      only : soil_energy_step_implicit
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_vegetation_biophysics, only : veg_energy_diagnostic, veg_energy_step_implicit,      &
                                     sensible_heat_coeff, leaf_transp_coeff
   use meds_soil_water,       only : column_hydrology_flux
   use meds_ground_biophysics, only : snow_energy_step, snow_base_conductance,                  &
                                     snow_accumulate, snow_drain_meltwater, snow_cover_fraction, &
                                     ground_surface_fluxes
   use meds_plant_interface,  only : leaf_gas_exchange_batch,                                  &
                                     stem_maintenance_respiration,                             &
                                     fine_root_maintenance_respiration,                        &
                                     hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_soil_biogeochem,  only : heterotrophic_respiration_flux, heterotrophic_respiration_matrix, &
                                     assemble_env_scalar, assemble_transfer_matrix
   use meds_biogeochem_types, only : co2_opts_t, n_soil_pool
   use meds_therm_lib,           only : cas_temp_of_enthalpy, cas_enthalpy_of_temp, sat_specific_humidity, &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_liquid,  &
                                     sat_vapor_pressure, uext_to_temp, temp_to_uext
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok, budget_check_stop
   implicit none
   private

   public :: column_fast_step_ark, aero_bottom_to_top, column_prepass
   public :: ark2_column_step, adaptive_ark_march, bflux_zero, bflux_add
   public :: column_be_stage, advance_hydraulics_full, state_init

contains

   !---------------------------------------------------------------------------------------!
   ! column_be_stage -- ONE backward-Euler stage of the STIFF, backward-Euler-integrable block only:  !
   ! the CAS twins (BE-in-atm; np>1 -> the coupled 2x2 Newton arrowhead) + soil heat + soil water     !
   ! (BE-Thomas), driven by the frozen surface sources. Plant hydraulics is DELIBERATELY EXCLUDED --  !
   ! solve_plant_water is an EXACT matrix exponential, not a backward-Euler stage, so it cannot ride   !
   ! an ESDIRK accumulation (it would drop the order + overshoot psi); psi is PASSED THROUGH here and  !
   ! advanced separately by advance_hydraulics_full over the full step. This is the reusable ESDIRK    !
   ! stage primitive: for the CAS+soil block it solves Y = base + dt*f_I(Y), so both imex_euler_column_!
   ! step (gamma=1) and each ark2 stage (gamma*dt) are just a column_be_stage call. Reuses the         !
   ! validated production kernels -- no new numerics.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine column_be_stage(y, fro, n, nsl, dt, y_out, niter, bf)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter    !< 1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux RATES for the ARK ledger

      type(surface_state_t)      :: ys
      type(surface_frozen_t)     :: fs
      type(surface_tend_t)       :: sf
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      real(wp)    :: t_ground, fliq1, wmass1, wcap, ccap, gah, gaw, gac
      real(wp)    :: enth1, shv1, e_infil, e_runof, e_drain, t_cas1
      integer(ik) :: k, np, nfeval
      logical     :: ok

      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call uext_to_temp(y%soil_energy(1), wmass1, fro%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      wcap = fro%surf%wcap ; ccap = fro%surf%ccap
      gah  = fro%surf%gah  ; gaw  = fro%surf%gaw ; gac = fro%surf%gac
      fs = fro%surf ; fs%t_ground = t_ground

      !----- CAS enthalpy + humidity. np==1: the uncoupled single-BE-pass baseline. np>1: a DIRECT 2x2  !
      !      Newton solve of the coupled backward-Euler surface block (the arrowhead). The FINAL sf     !
      !      drives the soil sinks (single-flux-per-interface).                                         !
      if (np <= 1_ik) then
         ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
         call surface_derivs(ys, fs, n, sf)
         enth1 = (wcap*y%cas_enthalpy + dt*(sf%src_enth + gah*fro%surf%enth_atm)) / (wcap + dt*gah)
         shv1  = (wcap*y%cas_shv      + dt*(sf%src_vap  + gaw*fro%surf%shv_atm )) / (wcap + dt*gaw)
      else
         call newton_surface_solve(y, fs, n, dt, wcap, gah, gaw, enth1, shv1, sf, nfeval, ok)
      end if
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      y_out%cas_co2      = (ccap*y%cas_co2 + dt*(fro%surf%nee_biotic + gac*fro%surf%co2_atm)) / (ccap + dt*gac)

      !----- soil-heat column: implicit BE-Thomas (soil_energy_step_implicit). ---------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top = sf%g_top ; eforc%geothermal = fro%geothermal
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         eforc%root_heat_sink(k) = sf%coh_qsoil * fro%soil%root_frac(k)
         eforc%w_flux(k)         = 0.0_wp
      end do
      !----- boundary water-enthalpy advection (guard-lift): infiltrating throughfall carries its liquid !
      !      enthalpy INTO the soil top, runoff/drainage carry it OUT (state^n temps, frozen -> the fixed !
      !      source telescopes exactly under the b-weighted ledger). Matches the split :436-439. root_    !
      !      heat_sink is a SINK, so q_src = -sink/dz: subtract an inflow, add an outflow. --------------!
      e_infil = fro%infiltration * internal_energy_liquid(fro%rain_temp)
      e_runof = fro%runoff_surf  * internal_energy_liquid(fro%surf%t_ground)   ! t_ground frozen @ state^n
      e_drain = fro%drainage     * internal_energy_liquid(fro%t_bot)
      eforc%root_heat_sink(1)   = eforc%root_heat_sink(1)   - e_infil + e_runof
      eforc%root_heat_sink(nsl) = eforc%root_heat_sink(nsl) + e_drain
      call soil_energy_step_implicit(se, eforc, fro%therm, fro%soil, fro%energy_opts, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- soil water is OPERATOR-SPLIT OUT of the ESDIRK stages: theta is PASSED THROUGH (held at the   !
      !      stage input = theta^n) and the AUTHORITATIVE end-of-step theta is committed once, from the     !
      !      scratch column_hydrology_flux (fro%theta1), in column_fast_step_ark. Re-solving it here with a !
      !      relief-free single-BE Richards drifted to saturation over long wet runs (no ponding/runoff),   !
      !      then hung the next scratch solve; the robust ponding/runoff/free-drain solve is the SOLE       !
      !      soil-water authority now (the ED2 "single soil-water authority" principle). theta feeds only   !
      !      the t_ground diagnosis + the soil-energy thermal property above, both correctly at theta^n. ---!
      y_out%theta(1:nsl) = y%theta(1:nsl)

      !----- plant hydraulics PASSED THROUGH (advanced by advance_hydraulics_full, not here). ----!
      y_out%psi(:, 1:n) = y%psi(:, 1:n)

      !----- emit this stage's boundary-flux RATES for the conservation ledger (§2.3). The b-weight   !
      !      + cross-substep accumulation happens in ark2_column_step / adaptive_ark_march. Every      !
      !      quantity is the committed-state flux, so the accumulated amounts telescope to closure.    !
      if (present(bf)) then
         t_cas1 = cas_temp_of_enthalpy(enth1, shv1)         ! committed CAS temp for the dew liquid enthalpy
         associate (fs2 => fro%surf)
            bf%cas_enth_in  = sf%src_enth + gah*fs2%enth_atm    ; bf%cas_enth_out = gah*enth1
            bf%cas_vap_in   = sf%src_vap  + gaw*fs2%shv_atm     ; bf%cas_vap_out  = gaw*shv1
            bf%cas_co2_in   = fs2%nee_biotic + gac*fs2%co2_atm  ; bf%cas_co2_out  = gac*y_out%cas_co2
            bf%soil_enth_in = sf%g_top + fro%geothermal + e_infil
            bf%soil_enth_out= sf%coh_qsoil * sum(fro%soil%root_frac(1:nsl)) + e_runof + e_drain
            !----- soil water is out of the ARK: its storage delta + q_top/drainage/uptake fluxes are     !
            !      re-sourced once/step from the frozen hflux in column_fast_step_ark, so the per-stage    !
            !      bf carries ONLY the CAS-vapour exchange (drainage/runoff/precip are frozen fast-step).  !
            bf%soil_wat_in  = 0.0_wp                            ; bf%soil_wat_out = 0.0_wp
            !----- condensation (dew) leaves the CAS as liquid at Tcas -> a whole-column water + liquid-   !
            !      enthalpy OUTPUT (the CAS-side loss is already in src_vap/src_enth, so cas_water/energy   !
            !      close automatically). -----------------------------------------------------------------!
            bf%whole_enth_in= sf%coh_rnet + fs2%abs_sw_ground + fs2%abs_lw_ground + e_infil
            bf%whole_enth_out= gah*(enth1 - fs2%enth_atm) + e_runof + e_drain                             &
                             + sf%cond*internal_energy_liquid(t_cas1)
            bf%whole_wat_in = 0.0_wp                            ; bf%whole_wat_out = gaw*(shv1 - fs2%shv_atm) + sf%cond
         end associate
      end if
   end subroutine column_be_stage
   !---------------------------------------------------------------------------------------!
   ! newton_surface_solve -- the ARROWHEAD: a direct 2x2 Newton solve of the coupled backward-Euler   !
   ! CAS surface block { R_H, R_q } = 0 for (cas_enthalpy H1, cas_shv q1), where the surface sources   !
   ! src_enth/src_vap depend nonlinearly on (H1,q1) through tcas, qcas and qsat. Replaces the leaf<->  !
   ! CAS Picard iteration: quadratic convergence (~1 step) + robust near saturation. Numerical         !
   ! Jacobian (finite-difference surface_derivs) -- captures the strong VPD self-limiting d src_vap/dq  !
   ! with no derivation risk. Singular-Jacobian guard + line search + supersaturation clamp + eval cap; !
   ! never error stops (GPU-safe). Commits the CAS via the FLUX form so budgets close for ANY sf.      !
   !---------------------------------------------------------------------------------------!
   subroutine newton_surface_solve(y, fs, n, dt, wcap, gah, gaw, enth1, shv1, sf, nfeval, ok)
      type(column_state_t),   intent(in)    :: y
      type(surface_frozen_t), intent(in)    :: fs
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, wcap, gah, gaw
      real(wp),               intent(out)   :: enth1, shv1
      type(surface_tend_t),   intent(out)   :: sf
      integer(ik),            intent(out)   :: nfeval
      logical,                intent(out)   :: ok

      type(surface_state_t) :: ys
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt
      integer(ik) :: it, ls
      real(wp),    parameter :: RTOL_N = 1.0e-7_wp, ATOL_H = 5.0e1_wp, ATOL_Q = 1.0e-6_wp
      real(wp),    parameter :: DETEPS = 1.0e-30_wp
      integer(ik), parameter :: NEWT_MAX = 4_ik, LS_MAX = 6_ik, FEVAL_CAP = 24_ik

      H0 = y%cas_enthalpy ; q0 = y%cas_shv
      Hk = H0 ; qk = q0 ; nfeval = 0_ik ; ok = .false.
      ys%cas_co2 = y%cas_co2
      ys%cas_enthalpy = Hk ; ys%cas_shv = qk
      call surface_derivs(ys, fs, n, sf) ; nfeval = nfeval + 1_ik
      R_H = wcap*(Hk - H0)/dt - sf%src_enth - gah*(fs%enth_atm - Hk)
      R_q = wcap*(qk - q0)/dt - sf%src_vap  - gaw*(fs%shv_atm  - qk)

      do it = 1_ik, NEWT_MAX
         if ( abs(R_H)*dt/wcap <= ATOL_H + RTOL_N*abs(Hk) .and.                                  &
              abs(R_q)*dt/wcap <= ATOL_Q + RTOL_N*abs(qk) ) then
            ok = .true. ; exit
         end if
         call jac_surface(Hk, qk, y%cas_co2, fs, sf, n, wcap, gah, gaw, dt, J11, J12, J21, J22, nfeval)
         detJ = J11*J22 - J12*J21
         if (detJ <= DETEPS*abs(J11*J22) .or. detJ <= 0.0_wp) then       ! singular / sign-flipped guard
            delH = -R_H / max(J11, tiny_num)                             ! damped-diagonal (Picard-like) fallback
            delq = -R_q / max(J22, tiny_num)
         else
            delH = (-R_H*J22 + R_q*J12)/detJ
            delq = (-R_q*J11 + R_H*J21)/detJ
         end if
         lam = 1.0_wp ; rn0 = R_H*R_H + R_q*R_q ; Ht = Hk ; qt = qk ; RHt = R_H ; Rqt = R_q
         do ls = 1_ik, LS_MAX
            Ht = Hk + lam*delH ; qt = qk + lam*delq
            !----- NO supersaturation STATE clamp: ED2 never clamps can_shv to SAT*qsat -- doing so is a   !
            !      harsh state discontinuity that inflates the ARK2 embedded error in cas_shv AND cas_     !
            !      enthalpy near RH=1 and thrashes the adaptive controller. Like ED2 we TOLERATE transient !
            !      supersaturation; the smooth condensation SINK in surface_derivs relaxes it physically.  !
            ys%cas_enthalpy = Ht ; ys%cas_shv = qt
            call surface_derivs(ys, fs, n, sf) ; nfeval = nfeval + 1_ik
            RHt = wcap*(Ht - H0)/dt - sf%src_enth - gah*(fs%enth_atm - Ht)
            Rqt = wcap*(qt - q0)/dt - sf%src_vap  - gaw*(fs%shv_atm  - qt)
            if (RHt*RHt + Rqt*Rqt <= (1.0_wp - 1.0e-4_wp*lam)*rn0) exit          ! Armijo
            lam = 0.5_wp*lam
         end do
         Hk = Ht ; qk = qt ; R_H = RHt ; R_q = Rqt
         if (nfeval >= FEVAL_CAP) exit
      end do

      !----- authoritative final eval + flux-form commit (conservation holds for ANY sf). -----------!
      ys%cas_enthalpy = Hk ; ys%cas_shv = qk
      call surface_derivs(ys, fs, n, sf) ; nfeval = nfeval + 1_ik
      enth1 = (wcap*H0 + dt*(sf%src_enth + gah*fs%enth_atm)) / (wcap + dt*gah)
      shv1  = (wcap*q0 + dt*(sf%src_vap  + gaw*fs%shv_atm )) / (wcap + dt*gaw)
   end subroutine newton_surface_solve

   !----- 2x2 numerical Jacobian of (R_H, R_q) w.r.t. (H, q) by forward-differencing surface_derivs. --!
   subroutine jac_surface(Hk, qk, co2, fs, sf, n, wcap, gah, gaw, dt, J11, J12, J21, J22, nfeval)
      real(wp),               intent(in)    :: Hk, qk, co2, wcap, gah, gaw, dt
      type(surface_frozen_t), intent(in)    :: fs
      type(surface_tend_t),   intent(in)    :: sf         ! base eval at (Hk,qk)
      integer(ik),            intent(in)    :: n
      real(wp),               intent(out)   :: J11, J12, J21, J22
      integer(ik),            intent(inout) :: nfeval
      type(surface_state_t) :: ys
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq, dse_dH, dsv_dH, dse_dq, dsv_dq
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      ys%cas_co2 = co2
      ys%cas_enthalpy = Hk + dH ; ys%cas_shv = qk
      call surface_derivs(ys, fs, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dH = (sfp%src_enth - sf%src_enth)/dH ; dsv_dH = (sfp%src_vap - sf%src_vap)/dH
      ys%cas_enthalpy = Hk ; ys%cas_shv = qk + dq
      call surface_derivs(ys, fs, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dq = (sfp%src_enth - sf%src_enth)/dq ; dsv_dq = (sfp%src_vap - sf%src_vap)/dq
      J11 = wcap/dt + gah - dse_dH ; J12 =              - dse_dq
      J21 =              - dsv_dH  ; J22 = wcap/dt + gaw - dsv_dq
   end subroutine jac_surface
   !----- copy the prognostic state (used to seed the RK combination). --------------------!
   pure subroutine state_init(y, n, nsl, ys)
      type(column_state_t), intent(in)  :: y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: ys
      ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
      ys%soil_energy  = y%soil_energy  ; ys%theta   = y%theta
      allocate(ys%psi(N_HYDRO, n))
      ys%psi(:, 1:n) = y%psi(:, 1:n)
   end subroutine state_init
   !---------------------------------------------------------------------------------------!
   ! ark2_column_step -- one 2nd-order L-stable IMEX step via the ARS(2,2,2) additive Runge-Kutta      !
   ! (Ascher-Ruuth-Spiteri 1997, Appl.Numer.Math. 25:151; identical gamma in Giraldo et al. 2013     !
   ! "ARK2"). Stiffly accurate on both tableaux, so the two ESDIRK stages map onto imex_euler's       !
   ! "last BE solve = committed state" structure with dt -> gamma*dt. The biotic CO2 source is folded !
   ! IMPLICIT (stays in the CO2 BE numerator), so f_E == 0 and the scheme reduces to a clean 2-solve  !
   ! ESDIRK2. PLANT HYDRAULICS IS OPERATOR-SPLIT OUT of the tableau (solve_plant_water is an EXACT     !
   ! matrix exponential, NOT a backward-Euler stage -- putting it in the ESDIRK accumulation silently  !
   ! drops the order and overshoots psi): psi is frozen through the stages, then advanced once over    !
   ! the full dt, and excluded from the embedded error. y_err = (Y3-base3)-(Y2-y_n) is the free        !
   ! embedded 1st-order estimate for the adaptive controller (2 solves/step vs step-doubling's 3).     !
   !---------------------------------------------------------------------------------------!
   subroutine ark2_column_step(y, fro, n, nsl, dt, y_out, y_err, niter, bf, hyd_nsub, hyd_nonconv)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      integer(ik), optional, intent(in)  :: niter
      type(column_bflux_t), optional, intent(out) :: bf   !< b-weighted boundary-flux AMOUNTS over dt (ledger)
      integer(ik), optional, intent(out) :: hyd_nsub, hyd_nonconv  !< section 5.3 work counters (pass-through)
      real(wp), parameter :: GAMMA = 0.2928932188134524_wp   ! 1 - 1/sqrt(2)
      real(wp), parameter :: BETA  = 2.4142135623730951_wp   ! (1-gamma)/gamma = 1 + sqrt(2)
      type(column_state_t) :: Y2, base3, Y3
      type(stage_bflux_t)  :: bf2, bf3
      integer(ik) :: np
      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      !----- Stage 2: gamma*dt BE stage from y_n (CAS+soil only; psi frozen -- it is split out). -----!
      call column_be_stage(y, fro, n, nsl, GAMMA*dt, Y2, niter=np, bf=bf2)
      !----- Stage 3: extrapolated base. The BETA=2.414 extrapolation can overshoot BOTH the vG theta   !
      !      range AND the CAS enthalpy into a wild temperature where qsat(T) overflows to NaN; clamp     !
      !      both to physical ranges so the stage stays FINITE. This only bites on a genuinely oversized  !
      !      step (which the adaptive controller then rejects and shrinks normally) -- in range it is an  !
      !      identity, so no accuracy cost. Without the CAS clamp a big transient poisons the whole march !
      !      with NaN. --------------------------------------------------------------------------------!
      call state_extrap(y, BETA, Y2, n, nsl, base3)
      call clamp_theta(base3, fro, nsl)
      call clamp_cas(base3)
      call column_be_stage(base3, fro, n, nsl, GAMMA*dt, Y3, niter=np, bf=bf3)
      call state_init(Y3, n, nsl, y_out)
      !----- operator-split hydraulics: exact 2x2 over the FULL dt from y_n, endpoint transp. -------!
      call advance_hydraulics_full(y, fro, n, nsl, dt, y_out, nsub_out=hyd_nsub, nonconv_out=hyd_nonconv)
      !----- embedded 1st-order error estimate (psi is split out -> zeroed). -------------------------!
      call state_err_diff(Y3, base3, Y2, y, n, nsl, y_err)
      !----- b-weighted boundary-flux amounts: b^I = (0, 1-gamma, gamma) -> exact telescoping.         !
      !      (Water closure is exact only when clamp_theta is inactive; it barely moves over gamma*dt.) !
      if (present(bf)) call bflux_bweight(bf, bf2, bf3, dt, GAMMA)
   end subroutine ark2_column_step

   !----- ledger helpers: b-weight two stage RATE structs into accumulated AMOUNTS over dt (weights   !
   !      b^I = (1-gamma, gamma) times dt); zero an accumulator; add one substep's amounts. -----------!
   pure subroutine bflux_bweight(acc, s2, s3, dt, gam)
      type(column_bflux_t), intent(out) :: acc
      type(stage_bflux_t),  intent(in)  :: s2, s3
      real(wp),             intent(in)  :: dt, gam
      real(wp) :: b2, b3
      b2 = (1.0_wp - gam) * dt ; b3 = gam * dt
      acc%cas_enth_in   = b2*s2%cas_enth_in   + b3*s3%cas_enth_in
      acc%cas_enth_out  = b2*s2%cas_enth_out  + b3*s3%cas_enth_out
      acc%cas_vap_in    = b2*s2%cas_vap_in    + b3*s3%cas_vap_in
      acc%cas_vap_out   = b2*s2%cas_vap_out   + b3*s3%cas_vap_out
      acc%cas_co2_in    = b2*s2%cas_co2_in    + b3*s3%cas_co2_in
      acc%cas_co2_out   = b2*s2%cas_co2_out   + b3*s3%cas_co2_out
      acc%soil_enth_in  = b2*s2%soil_enth_in  + b3*s3%soil_enth_in
      acc%soil_enth_out = b2*s2%soil_enth_out + b3*s3%soil_enth_out
      acc%soil_wat_in   = b2*s2%soil_wat_in   + b3*s3%soil_wat_in
      acc%soil_wat_out  = b2*s2%soil_wat_out  + b3*s3%soil_wat_out
      acc%whole_enth_in = b2*s2%whole_enth_in + b3*s3%whole_enth_in
      acc%whole_enth_out= b2*s2%whole_enth_out+ b3*s3%whole_enth_out
      acc%whole_wat_in  = b2*s2%whole_wat_in  + b3*s3%whole_wat_in
      acc%whole_wat_out = b2*s2%whole_wat_out + b3*s3%whole_wat_out
   end subroutine bflux_bweight

   pure subroutine bflux_zero(acc)
      type(column_bflux_t), intent(out) :: acc
      acc = column_bflux_t()
   end subroutine bflux_zero

   pure subroutine bflux_add(acc, s)
      type(column_bflux_t), intent(inout) :: acc
      type(column_bflux_t), intent(in)    :: s
      acc%cas_enth_in   = acc%cas_enth_in   + s%cas_enth_in
      acc%cas_enth_out  = acc%cas_enth_out  + s%cas_enth_out
      acc%cas_vap_in    = acc%cas_vap_in    + s%cas_vap_in
      acc%cas_vap_out   = acc%cas_vap_out   + s%cas_vap_out
      acc%cas_co2_in    = acc%cas_co2_in    + s%cas_co2_in
      acc%cas_co2_out   = acc%cas_co2_out   + s%cas_co2_out
      acc%soil_enth_in  = acc%soil_enth_in  + s%soil_enth_in
      acc%soil_enth_out = acc%soil_enth_out + s%soil_enth_out
      acc%soil_wat_in   = acc%soil_wat_in   + s%soil_wat_in
      acc%soil_wat_out  = acc%soil_wat_out  + s%soil_wat_out
      acc%whole_enth_in = acc%whole_enth_in + s%whole_enth_in
      acc%whole_enth_out= acc%whole_enth_out+ s%whole_enth_out
      acc%whole_wat_in  = acc%whole_wat_in  + s%whole_wat_in
      acc%whole_wat_out = acc%whole_wat_out + s%whole_wat_out
   end subroutine bflux_add

   !----- out = (1-b)*y + b*Y2  (the ARS stage-3 extrapolation base). --------------------------!
   pure subroutine state_extrap(y, b, Y2, n, nsl, out)
      type(column_state_t), intent(in)  :: y, Y2
      real(wp),             intent(in)  :: b
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: out
      real(wp)    :: a
      integer(ik) :: k, i
      a = 1.0_wp - b
      out%cas_enthalpy = a*y%cas_enthalpy + b*Y2%cas_enthalpy
      out%cas_shv      = a*y%cas_shv      + b*Y2%cas_shv
      out%cas_co2      = a*y%cas_co2      + b*Y2%cas_co2
      out%soil_energy = y%soil_energy ; out%theta = y%theta
      do k = 1_ik, nsl
         out%soil_energy(k) = a*y%soil_energy(k) + b*Y2%soil_energy(k)
         out%theta(k)       = a*y%theta(k)       + b*Y2%theta(k)
      end do
      allocate(out%psi(N_HYDRO, n))
      do i = 1_ik, n
         out%psi(:, i) = a*y%psi(:, i) + b*Y2%psi(:, i)     ! == y%psi (psi frozen in the stages)
      end do
   end subroutine state_extrap

   !----- clamp the extrapolated CAS enthalpy + humidity into a wide PHYSICAL range so a BETA=2.414   !
   !      overshoot cannot drive cas_temp_of_enthalpy to a wild T where qsat(T) overflows to NaN. Only  !
   !      active on a pathological overshoot (then the step is rejected); an in-range base3 is untouched.!
   pure subroutine clamp_cas(s)
      type(column_state_t), intent(inout) :: s
      real(wp) :: t, shv_c
      real(wp), parameter :: T_LO = 180.0_wp, T_HI = 350.0_wp, SHV_LO = 1.0e-8_wp, SHV_HI = 0.06_wp
      shv_c = min(max(s%cas_shv, SHV_LO), SHV_HI)
      t     = cas_temp_of_enthalpy(s%cas_enthalpy, shv_c)
      t     = min(max(t, T_LO), T_HI)
      s%cas_shv      = shv_c
      s%cas_enthalpy = cas_enthalpy_of_temp(t, shv_c)
   end subroutine clamp_cas

   !----- clamp the extrapolated theta into [theta_res, theta_sat] (van Genuchten domain). -----!
   pure subroutine clamp_theta(s, fro, nsl)
      type(column_state_t),  intent(inout) :: s
      type(column_frozen_t), intent(in)    :: fro
      integer(ik),           intent(in)    :: nsl
      integer(ik) :: k
      do k = 1_ik, nsl
         s%theta(k) = min(max(s%theta(k), fro%soil%theta_res(k)), fro%soil%theta_sat(k))
      end do
   end subroutine clamp_theta

   !----- err = (Y3 - base3) - (Y2 - y)  (the embedded 2nd-1st order difference); psi zeroed. ---!
   pure subroutine state_err_diff(Y3, base3, Y2, y, n, nsl, err)
      type(column_state_t), intent(in)  :: Y3, base3, Y2, y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: err
      integer(ik) :: k
      err%cas_enthalpy = (Y3%cas_enthalpy - base3%cas_enthalpy) - (Y2%cas_enthalpy - y%cas_enthalpy)
      err%cas_shv      = (Y3%cas_shv      - base3%cas_shv)      - (Y2%cas_shv      - y%cas_shv)
      err%cas_co2      = (Y3%cas_co2      - base3%cas_co2)      - (Y2%cas_co2      - y%cas_co2)
      err%soil_energy = 0.0_wp ; err%theta = 0.0_wp
      do k = 1_ik, nsl
         err%soil_energy(k) = (Y3%soil_energy(k) - base3%soil_energy(k)) - (Y2%soil_energy(k) - y%soil_energy(k))
         err%theta(k)       = (Y3%theta(k)       - base3%theta(k))       - (Y2%theta(k)       - y%theta(k))
      end do
      allocate(err%psi(N_HYDRO, n))
      err%psi(:, 1:n) = 0.0_wp
   end subroutine state_err_diff

   !----- operator-split plant hydraulics: exact 2x2 matrix-exp over the FULL dt from y%psi, driven   !
   !      by the ENDPOINT transpiration demand (surface_derivs at the committed CAS). ---------------!
   subroutine advance_hydraulics_full(y, fro, n, nsl, dt, y_out, nsub_out, nonconv_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: fro
      integer(ik),           intent(in)    :: n, nsl
      real(wp),              intent(in)    :: dt
      type(column_state_t),  intent(inout) :: y_out
      !----- section 5.3 WORK counters (optional so the RK4 oracle's call is unchanged): hydraulics    !
      !      sub-steps summed over cohorts, and cohorts whose solve did not converge. ----------------!
      integer(ik), optional, intent(out)   :: nsub_out, nonconv_out
      type(surface_state_t)  :: ys
      type(surface_frozen_t) :: fs
      type(surface_tend_t)   :: sf
      type(hydro_env_t)      :: henv
      type(hydro_flux_t)     :: hfx
      real(wp)    :: tg, fl, psi_i(N_HYDRO)
      integer(ik) :: i
      if (present(nsub_out))    nsub_out    = 0_ik
      if (present(nonconv_out)) nonconv_out = 0_ik
      call uext_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        fro%therm%soil_dry_heat_capacity(1), tg, fl)
      fs = fro%surf ; fs%t_ground = tg
      ys%cas_enthalpy = y_out%cas_enthalpy ; ys%cas_shv = y_out%cas_shv ; ys%cas_co2 = y_out%cas_co2
      call surface_derivs(ys, fs, n, sf)
      do i = 1_ik, n
         henv%transp     = sf%transp_c(i) * fro%surf%src_frac / max(fro%nplant(i), tiny_num)
         henv%soil_psi   = fro%soil_psi_root ; henv%rhizo_cond = fro%rhizo_cond
         henv%bleaf      = fro%bleaf(i) ; henv%bsap = fro%bsap(i) ; henv%broot = fro%broot(i)
         henv%sap_area   = fro%sap_area(i) ; henv%height = fro%height(i) ; henv%leaf_area = fro%leaf_area(i)
         psi_i = y%psi(:, i)
         call solve_plant_water(henv, fro%hydro_p, fro%hydro_o, dt, psi_i, hfx)
         y_out%psi(:, i) = psi_i
         if (present(nsub_out))    nsub_out    = nsub_out    + hfx%nsub
         if (present(nonconv_out)) then
            if (.not. hfx%converged) nonconv_out = nonconv_out + 1_ik
         end if
      end do
   end subroutine advance_hydraulics_full

   !----- out = a - b  (state difference; used to form the low-order embedded solution). --------!
   pure subroutine state_sub(a, b, n, nsl, out)
      type(column_state_t), intent(in)  :: a, b
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: out
      integer(ik) :: k
      out%cas_enthalpy = a%cas_enthalpy - b%cas_enthalpy
      out%cas_shv      = a%cas_shv      - b%cas_shv
      out%cas_co2      = a%cas_co2      - b%cas_co2
      out%soil_energy = a%soil_energy ; out%theta = a%theta
      do k = 1_ik, nsl
         out%soil_energy(k) = a%soil_energy(k) - b%soil_energy(k)
         out%theta(k)       = a%theta(k)       - b%theta(k)
      end do
      allocate(out%psi(N_HYDRO, n))
      out%psi(:, 1:n) = a%psi(:, 1:n) - b%psi(:, 1:n)
   end subroutine state_sub

   !---------------------------------------------------------------------------------------!
   ! adaptive_ark_march -- integrate to t_end with the ARK2 embedded error estimate driving the       !
   ! step controller (2 solves/step, no step-doubling). y_err (already the 2nd-1st order difference)   !
   ! is the local error; the WRMS of it vs tolerance drives accept/reject via adaptive_step_update     !
   ! (p=1 embedded -> exponent -1/2). Reports the step + reject count.                                 !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_ark_march(y0, fro, n, nsl, t_end, ec, dt_init, y_out, nsteps, nrej, niter, acc, &
                                 hyd_nsub, hyd_nonconv)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec       !< tolerances + controller + strictness (meds_fast_control)
      type(column_state_t),  intent(out) :: y_out
      integer(ik),           intent(out) :: nsteps, nrej
      integer(ik), optional, intent(in)  :: niter    !< coupled leaf<->CAS Newton cap (default 8)
      integer(ik), optional, intent(out) :: hyd_nsub, hyd_nonconv   !< section 5.3 work counters (summed)
      type(column_bflux_t), optional, intent(out) :: acc  !< accumulated boundary-flux amounts (ledger)

      type(column_state_t) :: y, y_new, y_err, y_lo
      type(column_bflux_t) :: bfsub
      real(wp)             :: t, dt, err, err_prev, fac, dt_floor
      integer(ik)          :: np, hns, hnc, hns_tot, hnc_tot

      np = 8_ik ; if (present(niter)) np = max(1_ik, niter)
      hns_tot = 0_ik ; hnc_tot = 0_ik
      if (present(acc)) call bflux_zero(acc)
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
      do
         if (t >= t_end - tiny_num) exit
         dt = min(dt, t_end - t)
         call ark2_column_step(y, fro, n, nsl, dt, y_new, y_err, niter=np, bf=bfsub,                   &
                               hyd_nsub=hns, hyd_nonconv=hnc)
         !----- work accounting: count EVERY attempt, not just accepted ones -- a rejected step has    !
         !      already paid for its stages, and hiding that would flatter an over-rejecting scheme.   !
         hns_tot = hns_tot + hns ; hnc_tot = hnc_tot + hnc
         call state_sub(y_new, y_err, n, nsl, y_lo)               ! the 1st-order embedded solution
         !----- per-group WRMS(y_err). with_psi=.false.: state_err_diff zeroes err%psi (psi rides an     !
         !      operator-split exponential map outside the tableau), so y_lo%psi == y_new%psi exactly    !
         !      and those 2n terms are structurally zero. Counting them divided the norm by ~1.4 and     !
         !      ran the march looser than its stated tolerance -- see state_wrms_grouped's header. ------!
         err = state_wrms_grouped(y_new, y_lo, y, n, nsl, ec%tols, with_psi=.false.)
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
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         if (nsteps + nrej > 4096_ik) exit                        ! hard backstop (should never trigger)
      end do
      call state_init(y, n, nsl, y_out)
      if (present(hyd_nsub))    hyd_nsub    = hns_tot
      if (present(hyd_nonconv)) hyd_nonconv = hnc_tot
   end subroutine adaptive_ark_march
   !=======================================================================================!
   !  INTEG_ARK path: the coupled IMEX-ARK fast step (docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md). Shares the   !
   !  split's frozen pre-pass (build_column_frozen), packs the state into the pure column vector,     !
   !  advances one dt_fast with the ARK stepper, then unpacks. PARTIAL precip>0 guard-lift: the ARK   !
   !  now carries the split's soil-boundary water-enthalpy advection (rain/runoff/drainage liquid      !
   !  enthalpy, in column_be_stage) and persists the scratch hydrology's ponding/aquifer/water-table   !
   !  (column_state_t still doesn't advance them prognostically -> a lagged operator split, so the      !
   !  whole-WATER budget closes only to the split-error tolerance, not machine). STILL restricted to   !
   !  free-drain + no Zeng-Decker: those bottom BCs need prognostic aquifer/z_wt in the state vector.  !
   !=======================================================================================!
   subroutine column_fast_step_ark(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,  &
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
      type(column_state_t)   :: y, y_out, ycur, ytmp, yerr
      type(surface_state_t)  :: ys
      type(surface_frozen_t) :: fs
      type(surface_tend_t)   :: sf
      type(column_bflux_t)   :: acc, bfsub
      real(wp)    :: tg, fl, dt0, wcap, ccap, enth0, shv0, co20, enth1, shv1, co21, e_soil0, e_soil1, w_soil0, w_soil1
      real(wp)    :: w_surface0
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, k, isub, nsub, nsteps, nrej, hns, hnc, hns1, hnc1
      logical     :: halt_budgets     !< §5.1: hard-stop on a non-closing budget (full column + debug only)

      n = coh%n ; nsl = ccfg%soil%n_active

      !----- §5.1: a reduced column freezes a store while its fluxes still act on the neighbours, so it  !
      !      cannot conserve by construction -- suppress the HARD stops (soft n_fail counters still run). !
      halt_budgets = ccfg%energy%debug_error .and. mask_is_full(ccfg%mask)

      !----- bottom-BC guard (see header): free-drain + no Zeng-Decker only. precip>0 is now supported  !
      !      (partial guard-lift); the aquifer/water-table bottom BCs still need prognostic state.       !
      if (ccfg%hydro%zeng_decker .or. ccfg%hydro%bottom_bc /= SOIL_BC_FREE_DRAIN)                 &
         error stop 'column_fast_step_ark: INTEG_ARK requires a free-drain bottom BC (no aquifer/Zeng-Decker yet)'
      w_surface0 = bio%soil_w%w_surface

      call build_column_frozen(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, n, nsl, &
                               fro, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)

      !----- advance one dt_fast: adaptive (embedded-error) or GPU-warp-uniform fixed substeps. ----!
      if (cfg%ark_adaptive) then
         dt0 = dt_fast ; if (cfg%ark_dt_init > tiny_num) dt0 = min(cfg%ark_dt_init, dt_fast)
         !----- The UNIFIED error-control bundle (§8c Layer 1): build_error_control seeds every tolerance !
         !      group from the setting that governs it today (and honours the [fast].rtol_all master      !
         !      dial), plus the controller + strictness. Defaults (CTRL_I, CTRL_L1, rtol_all unset)       !
         !      reproduce the legacy march byte-for-byte. ------------------------------------------------!
         ec = build_error_control(cfg)
         call adaptive_ark_march(y, fro, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                                 niter=cfg%ark_niter, acc=acc, hyd_nsub=hns, hyd_nonconv=hnc)
      else
         nsub = max(1_ik, cfg%ark_fixed_substep) ; nrej = 0_ik ; ycur = y ; call bflux_zero(acc)
         hns = 0_ik ; hnc = 0_ik
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, fro, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr,          &
                                  niter=cfg%ark_niter, bf=bfsub, hyd_nsub=hns1, hyd_nonconv=hnc1)
            hns = hns + hns1 ; hnc = hnc + hnc1
            call bflux_add(acc, bfsub)
            ycur = ytmp
         end do
         y_out = ycur ; nsteps = nsub
      end if
      !----- section 5.3 WORK counters: record what the march actually cost. -------------------------!
      budg%integ_nsteps = nsteps ; budg%integ_nrej = nrej
      budg%hydro_nsub = hns    ; budg%hydro_nonconv = hnc

      !----- SOIL WATER is operator-split out: the ESDIRK stages passed theta through unchanged (=theta^n); !
      !      commit the AUTHORITATIVE end-of-step theta from the scratch column_hydrology_flux HERE, once,  !
      !      so a single consistent theta feeds the state commit, the soil_temp read-off, and BOTH the      !
      !      soil_water and whole_water storage terms (w_soil1 below). ------------------------------------!
      y_out%theta(1:nsl) = fro%theta1(1:nsl)

      !----- §5.1 PROCESS MASK. The mask must mean the same thing under every scheme, so it is applied  !
      !      at the ARK's single state-commit point: a masked-off component is restored to state^n (y),  !
      !      leaving the ODE one dimension smaller while its couplings still acted during the march.     !
      !      This mirrors the split path's freeze exactly. mask%veg_energy needs no case here -- the ARK !
      !      error-stops on prognostic leaf/wood, so no vegetation energy store exists on this path. ----!
      if (.not. ccfg%mask%cas_energy) y_out%cas_enthalpy        = y%cas_enthalpy
      if (.not. ccfg%mask%cas_vapour) y_out%cas_shv             = y%cas_shv
      if (.not. ccfg%mask%cas_co2)    y_out%cas_co2             = y%cas_co2
      if (.not. ccfg%mask%soil_heat)  y_out%soil_energy(1:nsl)  = y%soil_energy(1:nsl)
      if (.not. ccfg%mask%soil_water) y_out%theta(1:nsl)        = y%theta(1:nsl)
      if (.not. ccfg%mask%hydraulics) y_out%psi(:, 1:n)         = y%psi(:, 1:n)

      !----- unpack into bio + re-derive the diagnostic soil temperatures + leaf temperatures. -----!
      bio%cas%can_enthalpy = y_out%cas_enthalpy ; bio%cas%can_shv = y_out%cas_shv ; bio%cas%can_co2 = y_out%cas_co2
      bio%cas%can_temp = cas_temp_of_enthalpy(y_out%cas_enthalpy, y_out%cas_shv)
      bio%soil_e%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      bio%soil_w%theta(1:nsl)       = y_out%theta(1:nsl)
      bio%psi(:, 1:n)               = y_out%psi(:, 1:n)
      !----- persist the scratch hydrology's ponding/aquifer/water-table (lagged operator split).       !
      !      §5.1: these are part of the SOIL-WATER store, so they must obey the same freeze as theta   !
      !      -- otherwise mask%soil_water=.false. means something different on this path than on the    !
      !      split path (which restores the whole soil_column_t), and the two schemes are no longer     !
      !      running the same reduced system. They are the ONLY writes to bio%soil_w besides theta, and !
      !      the hydrology ran on soil_w_scratch, so skipping them leaves the store at state^n. --------!
      if (ccfg%mask%soil_water) then
         bio%soil_w%w_surface = fro%w_surface1
         bio%soil_w%w_aquifer = fro%w_aquifer1
         bio%soil_w%z_wt      = fro%z_wt1
      end if
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

      !----- WHOLE-COLUMN CONSERVATION LEDGER: close the same 7 budgets the split closes, using the     !
      !      b-weighted boundary-flux AMOUNTS accumulated over the substeps (acc). The flux-form CAS    !
      !      commits + the energy_resid=0 soil-heat column make the identity exact -> machine-precision !
      !      closure for ENERGY (incl. the frozen rain/runoff/drainage advection, a fixed source). dt=1  !
      !      because acc holds AMOUNTS, not rates. Whole-WATER carries the lagged ponding split, so it   !
      !      closes only to the operator-split tolerance. ---------------------------------------------!
      wcap = fro%surf%wcap ; ccap = fro%surf%ccap
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv ; co21 = y_out%cas_co2
      e_soil0 = 0.0_wp ; e_soil1 = 0.0_wp ; w_soil0 = 0.0_wp ; w_soil1 = 0.0_wp
      do k = 1_ik, nsl
         e_soil0 = e_soil0 + y%soil_energy(k)     * ccfg%soil%dz(k)
         e_soil1 = e_soil1 + y_out%soil_energy(k) * ccfg%soil%dz(k)
         w_soil0 = w_soil0 + y%theta(k)     * ccfg%soil%dz(k) * rho_h2o
         w_soil1 = w_soil1 + y_out%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      !----- L2/debug_error mode (ccfg%energy%debug_error) promotes a non-closing budget from a       !
      !      silently-counted n_fail to a hard `error stop` -- the enforced half of the conservation   !
      !      check (plan MEDS_NUMERICS_SCOPING.md sec 4/QW2), mirroring the split path; off by         !
      !      default so production behaviour is unchanged. Each check reuses budg%*%resid, which        !
      !      budget_accumulate just set as a side effect. --------------------------------------------!
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                             1.0_wp, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_check_stop(budg%cas_energy%resid, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp,        &
                             'cas_energy (ark)', halt_budgets)
      call budget_accumulate(budg%cas_water,  wcap*shv0,  wcap*shv1,  acc%cas_vap_in,  acc%cas_vap_out,  &
                             1.0_wp, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_check_stop(budg%cas_water%resid, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp,      &
                             1.0e-10_wp, 'cas_water (ark)', halt_budgets)
      call budget_accumulate(budg%cas_co2,    ccap*co20,  ccap*co21,  acc%cas_co2_in,  acc%cas_co2_out,  &
                             1.0_wp, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      call budget_check_stop(budg%cas_co2%resid, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp,             &
                             'cas_co2 (ark)', halt_budgets)
      call budget_accumulate(budg%soil_energy, e_soil0, e_soil1, acc%soil_enth_in, acc%soil_enth_out,    &
                             1.0_wp, abs(e_soil1) + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp)
      call budget_check_stop(budg%soil_energy%resid, abs(e_soil1) + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp,  &
                             'soil_energy (ark)', halt_budgets)
      !----- SOIL WATER (fully frozen now): storage theta^n -> theta1 (w_soil0 -> w_soil1, both from the    !
      !      scratch solve), inflow q_top*rho, outflow drainage + realized uptake -- all from the frozen    !
      !      hflux, which closed its OWN mass budget to machine precision inside column_hydrology_flux. -----!
      call budget_accumulate(budg%soil_water,  w_soil0, w_soil1,                                        &
                             fro%q_top*rho_h2o*dt_fast, (fro%drainage + fro%uptake)*dt_fast,            &
                             1.0_wp, max(w_soil1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
      call budget_check_stop(budg%soil_water%resid, max(w_soil1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp,    &
                             'soil_water (ark)', halt_budgets)
      !----- whole-WATER: precip IN; drainage + runoff + CAS-vapour OUT; ponding in the store. The soil +  !
      !      ponding + drainage/runoff/precip terms are frozen fast-step amounts; the CAS-vapour exchange   !
      !      gaw*(shv-shv_atm) is the ARK-accumulated part (acc%whole_wat_out). Unlike the SPLIT (one       !
      !      transp value feeds BOTH the soil sink and the CAS source, so it closes to ~machine), the ARK   !
      !      RE-EVALUATES transpiration per ESDIRK stage as the CAS VPD evolves, while the committed soil   !
      !      theta lost the FROZEN scratch uptake_total. That internal transp<->uptake flux therefore does  !
      !      NOT cancel to machine: the whole-water residual is the intra-step transpiration-demand swing,  !
      !      bounded by the transpiration flux over the step. Scale the tolerance to that lag (all OTHER 6  !
      !      budgets, incl. soil_water, still close to machine). ------------------------------------------!
      call budget_accumulate(budg%whole_water, w_soil0 + wcap*shv0 + w_surface0,                        &
                             w_soil1 + wcap*shv1 + fro%w_surface1,                                      &
                             acc%whole_wat_in + forc%precip*dt_fast,                                    &
                             acc%whole_wat_out + (fro%runoff_surf + fro%drainage)*dt_fast,              &
                             1.0_wp, max(w_soil1 + wcap*shv1 + fro%w_surface1, 1.0_wp), 1.0e-6_wp,      &
                             max(1.0e-3_wp, abs(fro%uptake)*dt_fast))
      call budget_check_stop(budg%whole_water%resid, max(w_soil1 + wcap*shv1 + fro%w_surface1, 1.0_wp), &
                             1.0e-6_wp, max(1.0e-3_wp, abs(fro%uptake)*dt_fast),                    &
                             'whole_water (ark)', halt_budgets)
      call budget_accumulate(budg%whole_energy, e_soil0 + wcap*enth0, e_soil1 + wcap*enth1, acc%whole_enth_in, &
                             acc%whole_enth_out, 1.0_wp, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp)
      call budget_check_stop(budg%whole_energy%resid, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp, &
                             'whole_energy (ark)', halt_budgets)

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_ark
   !----- The SHARED pre-pass (once per sub-step; ED2 freezes gs/hydraulics per DTLSM): refreshes the !
   !      aerodynamics + CAS-derived scalars from the current state, then computes LEAF gas exchange   !
   !      (GPP/gs/Rd), the FROZEN per-cohort leaf-energy coefficients h_coeff_f/g_tr_f, stem/root       !
   !      maintenance respiration, the NEE assembly, and the CAS capacities/atm-exchange conductances.  !
   !      ONE authority for both integrators -- this is what keeps split/ARK GPP bit-for-bit -- called   !
   !      from column_fast_step (meds_fast_split, which Picard-iterates the leaf<->CAS balance from      !
   !      here) and from build_column_frozen below (which freezes these as explicit ARK macro-step       !
   !      inputs). `bio` is intent(in): callers that need the CAS temperature persisted (the split)      !
   !      write bio%cas%can_temp = tcas themselves right after the call.                                 !
   subroutine column_prepass(cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,                       &
                             tcas, qcas, press, rho, t_ground, h_coeff_f, g_tr_f,                       &
                             wcap, ccap, gah, gaw, gac, nee_biotic,                                     &
                             gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(in)    :: bio
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budg
      real(wp),                intent(out)   :: tcas, qcas, press, rho, t_ground
      real(wp),                intent(out)   :: h_coeff_f(:), g_tr_f(:)
      real(wp),                intent(out)   :: wcap, ccap, gah, gaw, gac, nee_biotic
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)

      !----- Bare-array batch I/O for the per-cohort physiology kernels (MEDS_NUMERICS_SCOPING.md).   !
      real(wp) :: par_arr(coh%n), vpd_arr(coh%n), gb_arr(coh%n), rho_mol_arr(coh%n), psi_leaf_arr(coh%n)
      real(wp) :: a_gross_arr(coh%n), gs_arr(coh%n), rd_arr(coh%n)
      real(wp) :: stem_resp_arr(coh%n), root_resp_arr(coh%n)
      real(wp) :: e_air, gsw_ms, can_dmol
      real(wp) :: gpp, ra_leaf, ra_stem, ra_root, rh, soil_temp_root, theta_mean
      real(wp) :: xi(n_soil_pool), a_mat(n_soil_pool, n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      integer(ik) :: i, k, n, nsl

      n = coh%n ; nsl = ccfg%soil%n_active

      !----- aerodynamics from the current CAS state. -------------------------------------------!
      tcas = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air ; t_ground = bio%soil_e%soil_temp(1)
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = t_ground
      call aero_bottom_to_top(ccfg%aero, aenv, ageom, n, coh, bio%leaf_temp, aero)

      !----- Root-weighted soil temperature + column-mean moisture (root / heterotrophic resp). !
      soil_temp_root = 0.0_wp ; theta_mean = 0.0_wp
      do k = 1_ik, nsl
         soil_temp_root = soil_temp_root + bio%soil_e%soil_temp(k) * ccfg%soil%root_frac(k)
         theta_mean     = theta_mean     + bio%soil_w%theta(k) * ccfg%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-ccfg%soil%soil_layer_z(nsl+1_ik), tiny_num)

      !----- LEAF gas exchange (GPP/gs/Rd), frozen leaf-energy coefficients, stem+root maint. resp. --!
      !      BARE-ARRAY batch seam (MEDS_NUMERICS_SCOPING.md): (1) assemble the per-cohort leaf-env      !
      !      arrays, (2) call the three physiology kernels over the WHOLE patch at once, (3) accumulate   !
      !      the patch totals + frozen leaf-energy coefficients. The accumulation keeps the SAME          !
      !      i=1..n order as the old inline loop, so gpp/ra_* and every per-cohort output are             !
      !      bit-for-bit identical (verified vs a git-stash baseline). ---------------------------------!
      gpp = 0.0_wp ; ra_leaf = 0.0_wp ; ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      e_air = qcas * press / (0.622_wp + 0.378_wp * qcas)          ! loop-invariant (was recomputed each i)
      do i = 1_ik, n
         rho_mol_arr(i)  = press / (r_gas * bio%leaf_temp(i))
         par_arr(i)      = forc%abs_par(i) / max(coh%lai(i), 0.1_wp) * forc%par_per_w
         vpd_arr(i)      = max(sat_vapor_pressure(bio%leaf_temp(i)) - e_air, 0.0_wp)
         gb_arr(i)       = aero%leaf_gbw(i) * rho_mol_arr(i)
         psi_leaf_arr(i) = bio%psi(NODE_LEAF, i)                   ! gather (psi is node-major => strided)
      end do
      call leaf_gas_exchange_batch(n, par_arr, bio%leaf_temp(1:n), vpd_arr, bio%cas%can_co2, press, &
                                   psi_leaf_arr, gb_arr, cfg, coh%pft(1:n),                          &
                                   coh%vcmax25(1:n), coh%rd25(1:n), a_gross_arr, gs_arr, rd_arr)
      !----- Elemental (§11): the array actuals drive the element-wise broadcast; `ccfg%wood`/       !
      !      `ccfg%root` (scalar PODs) and the patch-uniform `soil_temp_root` broadcast. -------------!
      call stem_maintenance_respiration(bio%wood_temp(1:n), coh%dbh(1:n), coh%height(1:n),           &
                                   coh%wai(1:n), coh%nplant(1:n), ccfg%wood, stem_resp_arr(1:n))
      call fine_root_maintenance_respiration(soil_temp_root, coh%broot(1:n), ccfg%root, root_resp_arr(1:n))
      do i = 1_ik, n
         gsw_ms  = gs_arr(i) / max(rho_mol_arr(i), tiny_num)
         gpp     = gpp     + a_gross_arr(i) * coh%leaf_area(i) * coh%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = a_gross_arr(i) * coh%leaf_area(i)
         ra_leaf = ra_leaf + rd_arr(i)      * coh%leaf_area(i) * coh%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = rd_arr(i) * coh%leaf_area(i)
         h_coeff_f(i) = sensible_heat_coeff(ccfg%veg_thermal%effarea_heat * coh%lai(i), aero%leaf_gbh(i), rho, cp_air)
         g_tr_f(i)    = leaf_transp_coeff(ccfg%veg_thermal%effarea_transp, coh%lai(i), aero%leaf_gbw(i), gsw_ms)
         ra_stem = ra_stem + stem_resp_arr(i) * coh%nplant(i)
         ra_root = ra_root + root_resp_arr(i) * coh%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = stem_resp_arr(i)
         if (present(root_resp_coh)) root_resp_coh(i) = root_resp_arr(i)
      end do

      !----- NEE = autotrophic (leaf Rd + stem + root) + heterotrophic Rh - GPP. Rh is EITHER the    !
      !      OLD constant-pool scalar form (soil_carbon_on = .false., bit-identical to before Part   !
      !      II) OR the matrix form over bio%soil_carbon -- the FROZEN per-patch pool held constant    !
      !      across today's sub-steps (B2, MEDS_SLOW_DYNAMICS_DESIGN.md Part II section 9): the day's   !
      !      total fast Rh then equals the daily soil_carbon_step's pool debit BY CONSTRUCTION, since   !
      !      both read the same frozen pool + the same per-pool env scalar xi (accumulated into         !
      !      budg%xi_step below for the caller to integrate into xi_int). ------------------------------!
      if (cfg%soil_carbon_on) then
         call assemble_env_scalar(t_ground, soil_temp_root, theta_mean, ccfg%soil%theta_res(1),      &
                                  ccfg%soil%theta_sat(1), bio%soil_carbon, cfg%soil_carbon, xi)
         call assemble_transfer_matrix(bio%soil_carbon, cfg%soil_carbon, a_mat, k_diag, er)
         rh = heterotrophic_respiration_matrix(a_mat, k_diag, xi, bio%soil_carbon)
         budg%xi_step = xi ; budg%rh_matrix_step = rh
      else
         rh = heterotrophic_respiration_flux(ccfg%fast_soil_carbon, soil_temp_root, theta_mean,      &
                                             ccfg%soil%theta_res(1), ccfg%soil%theta_sat(1), ccfg%co2)
      end if
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budg%gpp_last = gpp ; budg%nee_last = nee_biotic

      !----- CAS capacities + atm-exchange conductances (frozen across passes / the ARK macro-step). --!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      wcap = rho      * bio%cas%can_depth
      ccap = can_dmol * bio%cas%can_depth
      gah  = rho      * aero%ustar * aero%temp1
      gaw  = rho      * aero%ustar * aero%temp2
      gac  = can_dmol * aero%ustar * aero%temp2
   end subroutine column_prepass
   !----- Build the frozen ARK inputs: the shared column_prepass above (leaf gas exchange /            !
   !      respiration / CAS caps / aero, bit-identical to the split) + this integrator's own per-cohort !
   !      geometry/radiation/wood packing + the frozen hydrology BCs, into a column_frozen_t; also      !
   !      packs the prognostic state into a column_state_t.                                             !
   subroutine build_column_frozen(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, n, nsl, &
                                  fro, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(in)    :: bio
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budg
      integer(ik),             intent(in)    :: n, nsl
      type(column_frozen_t),   intent(out)   :: fro
      type(column_state_t),    intent(out)   :: y
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)

      type(chydro_forcing_t) :: hforc ; type(chydro_flux_t) :: hflux
      type(soil_column_t)    :: soil_w_scratch
      type(surface_state_t)  :: ys ; type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, nee_biotic, wcap, ccap, gah, gaw, gac
      integer(ik) :: i

      allocate(fro%surf%h_coeff_f(n), fro%surf%g_tr_f(n), fro%surf%abs_sw(n), fro%surf%abs_lw(n), fro%surf%lai(n))
      allocate(fro%surf%h_coeff_w(n), fro%surf%abs_sw_wood(n), fro%surf%abs_lw_wood(n), fro%surf%wai(n))
      allocate(fro%psi_e(nsl), fro%nplant(n), fro%bleaf(n), fro%bsap(n), fro%broot(n),            &
               fro%sap_area(n), fro%height(n), fro%leaf_area(n))
      allocate(y%psi(N_HYDRO, n))

      !----- the SHARED pre-pass (column_prepass above): leaf gas exchange / respiration / CAS caps /   !
      !      aero -- writes directly into the frozen struct's h_coeff_f/g_tr_f arrays. ------------------!
      call column_prepass(cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,                       &
                          tcas, qcas, press, rho, t_ground, fro%surf%h_coeff_f, fro%surf%g_tr_f,      &
                          wcap, ccap, gah, gaw, gac, nee_biotic,                                      &
                          gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)

      !----- per-cohort geometry + radiation + WOOD frozen inputs the ARK path needs (not shared with   !
      !      the split, which reads coh%/forc% directly instead of packing a frozen struct). -----------!
      do i = 1_ik, n
         fro%surf%lai(i)    = coh%lai(i)
         fro%surf%abs_sw(i) = forc%abs_sw(i) ; fro%surf%abs_lw(i) = forc%abs_lw(i)
         !----- WOOD frozen inputs: real diagnostic values, or ZERO when wood is not diagnostic (so   !
         !      surface_derivs' wood branch is a no-op; prognostic wood is operator-split in P2). -----!
         if (ccfg%wood_energy_model == WOODEN_DIAGNOSTIC) then
            fro%surf%wai(i)         = coh%wai(i)
            fro%surf%h_coeff_w(i)   = sensible_heat_coeff(pi * coh%wai(i), aero%wood_gbh(i), rho, cp_air)
            fro%surf%abs_sw_wood(i) = forc%abs_sw_wood(i)
            fro%surf%abs_lw_wood(i) = forc%abs_lw_wood(i)
         else
            fro%surf%wai(i) = 0.0_wp ; fro%surf%h_coeff_w(i) = 0.0_wp
            fro%surf%abs_sw_wood(i) = 0.0_wp ; fro%surf%abs_lw_wood(i) = 0.0_wp
         end if
         fro%nplant(i)   = coh%nplant(i)  ; fro%bleaf(i)    = coh%bleaf(i)  ; fro%bsap(i) = coh%bsap(i)
         fro%broot(i)    = coh%broot(i)   ; fro%sap_area(i) = coh%sap_area(i)
         fro%height(i)   = coh%height(i)  ; fro%leaf_area(i) = coh%leaf_area(i)
      end do

      !----- the rest of the frozen surface inputs: CAS caps/conductances from column_prepass + atm     !
      !      state + NEE. ---------------------------------------------------------------------------!
      fro%surf%leaf_emiss = ccfg%veg_thermal%leaf_emiss
      fro%surf%wcap = wcap ; fro%surf%ccap = ccap
      fro%surf%gah  = gah  ; fro%surf%gaw  = gaw ; fro%surf%gac = gac
      fro%surf%enth_atm = forc%enthalpy_atm ; fro%surf%shv_atm = forc%shv_atm ; fro%surf%co2_atm = forc%co2_atm
      fro%surf%nee_biotic = nee_biotic
      fro%surf%abs_sw_ground = forc%abs_sw_ground ; fro%surf%abs_lw_ground = forc%abs_lw_ground
      fro%surf%ggnet = aero%ggnet ; fro%surf%rho = rho ; fro%surf%press = press
      fro%surf%src_frac = 1.0_wp ; fro%surf%t_ground = t_ground

      !----- params + hydraulics BCs. -----------------------------------------------------------!
      fro%soil = ccfg%soil ; fro%therm = ccfg%soil_thermal ; fro%energy_opts = ccfg%energy
      fro%hydro_opts = ccfg%hydro ; fro%hydro_p = ccfg%hydro_p ; fro%hydro_o = ccfg%hydro_o
      fro%geothermal = 0.0_wp ; fro%rhizo_cond = ccfg%rhizo_cond ; fro%psi_e(1:nsl) = 0.0_wp

      !----- FROZEN hydrology BCs: total transp demand (surface_derivs @ state^n, src_frac=1), then a  !
      !      SCRATCH column_hydrology_flux for soil_evap / infiltration / psi_soil / uptake_total. ----!
      ys%cas_enthalpy = bio%cas%can_enthalpy ; ys%cas_shv = bio%cas%can_shv ; ys%cas_co2 = bio%cas%can_co2
      call surface_derivs(ys, fro%surf, n, sf0)
      hforc%precip_ground      = forc%precip
      hforc%root_uptake(1:nsl) = sf0%coh_transp * ccfg%soil%root_frac(1:nsl)
      hforc%t_ground           = t_ground ; hforc%q_air = qcas ; hforc%rho_air = rho
      hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)
      soil_w_scratch = bio%soil_w
      call column_hydrology_flux(soil_w_scratch, hforc, ccfg%soil, ccfg%hydro, dt_fast, hflux)
      budg%soil_nsub = hflux%nsub                 ! section 5.3 work counter (same seam on both schemes)
      fro%surf%soil_evap = hflux%soil_evap
      fro%q_top          = (hflux%infiltration - hflux%soil_evap) / rho_h2o
      fro%soil_psi_root  = root_weighted_psi(hflux%psi_soil, ccfg%soil%root_frac, nsl)
      if (sf0%coh_transp > tiny_num) fro%surf%src_frac = min(1.0_wp, hflux%uptake_total / sf0%coh_transp)

      !----- FROZEN boundary hydrology for the guard-lift: the rain/drainage/runoff water-enthalpy       !
      !      advection (state^n temps, matching the split) + the scratch's end-of-step ponding/aquifer/  !
      !      water-table (soil_w_scratch was advanced in place by column_hydrology_flux). ---------------!
      fro%infiltration = hflux%infiltration ; fro%drainage    = hflux%drainage
      fro%runoff_surf  = hflux%runoff_surf  ; fro%rain_temp   = tcas
      fro%uptake       = hflux%uptake_total
      fro%t_bot        = bio%soil_e%soil_temp(nsl)
      fro%w_surface1   = soil_w_scratch%w_surface
      fro%w_aquifer1   = soil_w_scratch%w_aquifer
      fro%z_wt1        = soil_w_scratch%z_wt
      !----- the AUTHORITATIVE committed soil moisture: soil_w_scratch was advanced IN PLACE by the robust  !
      !      column_hydrology_flux above, so its theta IS the end-of-step (relieved) soil water. -----------!
      allocate(fro%theta1(nsl))
      fro%theta1(1:nsl) = soil_w_scratch%theta(1:nsl)

      !----- pack the prognostic state. ---------------------------------------------------------!
      y%cas_enthalpy = bio%cas%can_enthalpy ; y%cas_shv = bio%cas%can_shv ; y%cas_co2 = bio%cas%can_co2
      y%soil_energy(1:nsl) = bio%soil_e%soil_energy(1:nsl)
      y%theta(1:nsl)       = bio%soil_w%theta(1:nsl)
      y%psi(:, 1:n)        = bio%psi(:, 1:n)
   end subroutine build_column_frozen
   !----- Solve canopy aerodynamics with the cohort order it CONTRACTS for -- BOTTOM(1)->TOP(n)  !
   !      -- from the height-DESCENDING column buffer. Only the wind cascade + the per-cohort       !
   !      boundary layers depend on order; the whole-canopy scalars (ustar/temp1/temp2/uh) do not.  !
   !      An ascending-height permutation `ord` reverses the per-cohort inputs; the per-cohort wind  !
   !      and leaf/wood conductance outputs are scattered back to gather order. Identity for n<=1,   !
   !      so single-cohort behaviour is bit-unchanged.                                               !
   subroutine aero_bottom_to_top(acfg, aenv, ageom, n, coh, leaf_temp, aero)
      type(aero_cfg_t),      intent(in)    :: acfg
      type(aero_env_t),      intent(in)    :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      integer(ik),           intent(in)    :: n
      type(column_cohort_t), intent(in)    :: coh
      real(wp),              intent(in)    :: leaf_temp(:)
      type(aero_out_t),      intent(inout) :: aero
      integer(ik) :: ord(n), k, j, imin
      real(wp)    :: hmin
      logical     :: used(n)
      real(wp)    :: h_bt(n), lai_bt(n), cr_bt(n), lt_bt(n), lw_bt(n), bd_bt(n)
      real(wp)    :: wind_bt(n), lgbh_bt(n), lgbw_bt(n), wgbh_bt(n), wgbw_bt(n)

      !----- ord(k) = gather index of the k-th cohort counting from the canopy BOTTOM. ----------!
      used = .false.
      do k = 1_ik, n
         imin = 0_ik ; hmin = huge(1.0_wp)
         do j = 1_ik, n
            if (.not. used(j) .and. coh%height(j) <= hmin) then ; hmin = coh%height(j) ; imin = j ; end if
         end do
         ord(k)    = imin ; used(imin) = .true.
         h_bt(k)   = coh%height(imin)     ; lai_bt(k) = coh%lai(imin)
         cr_bt(k)  = coh%crown(imin)      ; lt_bt(k)  = leaf_temp(imin)
         lw_bt(k)  = coh%leaf_width(imin) ; bd_bt(k)  = coh%branch_diam(imin)
      end do

      call canopy_aerodynamics(acfg, aenv, ageom, n, h_bt, lai_bt, cr_bt, lt_bt, lt_bt, lw_bt, bd_bt, aero)

      !----- aero%*(k) is now bottom->top; copy out, then scatter back to gather order. ----------!
      do k = 1_ik, n
         wind_bt(k) = aero%wind(k)     ; lgbh_bt(k) = aero%leaf_gbh(k) ; lgbw_bt(k) = aero%leaf_gbw(k)
         wgbh_bt(k) = aero%wood_gbh(k) ; wgbw_bt(k) = aero%wood_gbw(k)
      end do
      do k = 1_ik, n
         aero%wind(ord(k))     = wind_bt(k)
         aero%leaf_gbh(ord(k)) = lgbh_bt(k) ; aero%leaf_gbw(ord(k)) = lgbw_bt(k)
         aero%wood_gbh(ord(k)) = wgbh_bt(k) ; aero%wood_gbw(ord(k)) = wgbw_bt(k)
      end do
   end subroutine aero_bottom_to_top

end module meds_fast_ark
