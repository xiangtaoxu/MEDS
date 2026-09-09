!==========================================================================================!
! meds_fast_ark -- the production fast-loop integrator: an L-stable ESDIRK2 (ARS(2,2,2) tableau;   !
! the explicit part is empty, so despite the historical "IMEX-ARK" name this is a diagonally       !
! implicit scheme). Design: docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md.                                !
!                                                                                          !
! Hosts: the scheme entry column_fast_step_ark (called from meds_fast_step%column_fast_step); the  !
! frozen pre-pass build_column_frozen + column_prepass (leaf gas exchange, respiration, CAS         !
! capacities and conductances, the scratch soil-water solve, plant hydraulics); the stage solve    !
! column_be_stage with the 2x2 leaf<->CAS Newton newton_surface_solve; ark2_column_step and its    !
! embedded-error controller adaptive_ark_march; the operator-split plant water-mass update          !
! advance_water_mass_full (with the transpiration corrector) and the canopy-film update            !
! advance_surf_water_full; the column_state_t algebra (state_init/state_axpy/state_accum/          !
! state_extrap/state_err_diff/state_sub, clamp_cas/clamp_theta/clamp_soil_energy) and the           !
! boundary-flux ledger accumulators (bflux_*), which meds_fast_rk45 and the test-only oracle        !
! meds_fast_rk4_oracle import from here.                                                           !
!                                                                                          !
! Soil water is committed ONCE per dt_fast from the scratch column_hydrology_flux solve (the ARK   !
! stages pass theta through); the pond is carried in column_state_t but committed the same way.    !
! Whole-column water and energy close to round-off on every bottom BC (free-drain, bedrock,        !
! aquifer); see meds_budget_check and the ledgers at the end of column_fast_step_ark.              !
!==========================================================================================!
module meds_fast_ark
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, latent_heat_vap, rho_h2o, r_gas, pi, &
                                     tsupercool_liq, grav_head, cp_liq, t_3ple
   use meds_plant_hydraulics, only : rhizosphere_cond, solve_plant_water_batch
   use meds_core_diag_types,  only : CD_ANET, CD_AGROSS, CD_GSW, CD_GBW, CD_CI, CD_CS, CD_RD,      &
                                     CD_TRANSP, CD_BETA_STOM, CD_BETA_NONSTOM, CD_LEAF_TEMP,      &
                                     CD_WOOD_TEMP, CD_LEAF_VPD, CD_PSI_LEAF, CD_PSI_WOOD, CD_PLC, &
                                     CD_SAPFLOW, CD_ROOT_UPTAKE, CD_ABS_PAR, CD_ABS_SW,           &
                                     CD_ABS_LW, CD_WIND, CD_GPP_RATE, CD_LEAF_WATER, CD_WOOD_WATER
   use meds_hydr_lib, only : soil_hydr_cond_from_theta, soil_psi_from_theta, psi_from_water_content, &
                             water_content
   use meds_config,           only : meds_config_t, hydraulics_config_t,                          &
                                     INTEG_ARK, CTRL_L2_STRICT
   use meds_fast_control,     only : error_control_t, build_error_control, state_wrms_grouped,   &
                                     step_control_factor
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,        &
                                     energy_forcing_t, energy_opts_t, energy_flux_t,           &
                                     soil_column_t, soil_energy_column_t, chydro_forcing_t, chydro_flux_t, &
                                     SOIL_BC_AQUIFER, &
                                     snow_params_t, snow_env_t, snow_flux_t, snow_melt_t
   use meds_fast_time_derivs, only : surface_derivs, cas_conductances
   use meds_cas_biophysics,   only : cas_column_step_implicit, cas_column_t, cas_source_t
   use meds_numerics,         only : weighted_mean
   use meds_column_state_ops, only : state_init, state_axpy, state_accum, state_extrap, state_err_diff, &
                                     state_sub, bflux_zero, bflux_add, bflux_bweight, clamp_cas,        &
                                     clamp_theta, clamp_soil_energy, soil_water_store, soil_energy_store, &
                                     plant_water_store, canopy_film_store, deposit_condensate,           &
                                     clamp_canopy_film, unpack_column_state, diagnose_soil_temps
   use meds_fast_snow,        only : snow_stage_t, advance_snow_stage
   use meds_fast_types,       only : column_config_t, column_cohort_t, column_forcing_t,       &
                                     column_budget_t, alloc_column_cohort,                      &
                                     column_state_t, column_frozen_t, surface_state_t,          &
                                     surface_frozen_t, surface_tend_t, stage_bflux_t, column_bflux_t, &
                                     column_tend_t, mask_is_full
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_soil_energy,      only : soil_energy_step_implicit
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_vegetation_biophysics, only : veg_energy_balance,                                &
                                     sensible_heat_coeff, leaf_transp_coeff, leaf_film_coeff,     &
                                     intercept_canopy_layer
   use meds_soil_water,       only : column_hydrology_flux
   use meds_ground_biophysics, only : snow_energy_step, snow_base_conductance,                  &
                                     snow_accumulate, snow_drain_meltwater, snow_cover_fraction, &
                                     ground_surface_fluxes
   use meds_plant_interface,  only : leaf_gas_exchange_batch,                                  &
                                     stem_maintenance_respiration,                             &
                                     fine_root_maintenance_respiration, solve_plant_water_batch, &
                                     N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_soil_biogeochem,  only : heterotrophic_respiration_flux, heterotrophic_respiration_matrix, &
                                     assemble_env_scalar, assemble_transfer_matrix
   use meds_biogeochem_types, only : co2_opts_t, n_soil_pool
   use meds_therm_lib,           only : cas_molar_density, cas_temp_of_enthalpy, cas_enthalpy_of_temp, sat_specific_humidity, &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_liquid,  &
                                     sat_vapor_pressure, uext_to_temp, temp_to_uext, internal_energy_ice,      &
                                     temp_of_liquid_enthalpy
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok, budget_check_stop,  &
                                     budget_check, budget_energy_rate_floor,                    &
                                     budget_water_rate_floor, budget_co2_rate_floor
   implicit none
   private

   !----- `niter` on column_be_stage / ark2_column_step / adaptive_ark_march is a two-valued port:   !
   !      <= 1 selects the uncoupled single-BE pass, anything > 1 selects the coupled 2x2 Newton      !
   !      (whose own iteration cap is NEWT_MAX, not this number). Passing a named sentinel rather      !
   !      than a bare 8 stops it reading as a tunable count it never was -- see meds_config's          !
   !      ark_coupled. -------------------------------------------------------------------------------!
   integer(ik), parameter :: NEWT_COUPLED = 2_ik

   !----- Per-cohort plant-hydraulics sub-step count above which the solve is judged pathological     !
   !      rather than merely stiff (issue #104). Measured band: 1.0-1.2 ordinary, ~136 collapsed. ----!
   integer(ik), parameter :: HYDRO_NSUB_THRASH = 16_ik

   public :: column_fast_step_ark, aero_bottom_to_top, column_prepass, build_column_frozen
   public :: ark2_column_step, adaptive_ark_march
   public :: column_be_stage, advance_water_mass_full, advance_surf_water_full

   !=========================================================================================!
   ! TISSUE HEAT STORE -- ACTIVATION SWITCH. 0 = zero-inertia tissue; 1 = the store live.          !
   ! Currently 1: the store is ON.                                                              !
   !                                                                                          !
   ! Everything the store needs is BUILT AND VERIFIED: the exact exponential relaxation, real WAI + !
   ! sapwood allometry, the dry-wood/sapwood-water capacity split, and -- the hard part -- EXACT    !
   ! conservation on both schemes via the b-weighted tissue-temperature time integrals in           !
   ! column_bflux_t. With this scale at 0 every path reproduces the pre-store answers bit for bit,  !
   ! which is the property the whole design was built around ("diagnostic is the a_store -> 0 limit !
   ! of one formula, not a separate mode") and it is checked by the full suite passing at 0.        !
   !                                                                                          !
   ! WHY IT IS NOT ON YET. Turning it on flips four PHYSICS assertions on ARK that had only ever    !
   ! been exercised on the retired split path: daytime NEE goes net-release, and the leaf water     !
   ! potential comes out POSITIVE (+0.72 MPa) instead of under tension, i.e. transpiration is being !
   ! suppressed and leaf water accumulates. That is not a conservation failure -- every budget still !
   ! closes -- but it is unexplained, and the leaf store is far too small to explain it directly     !
   ! (cap_leaf = 2205 J/m2/K against h_coeff = 100 W/m2/K, so tau <~ 22 s and w_avg <~ 0.012).       !
   ! Turning the LEAF store on alone is measurably WORSE than both together, which is non-monotonic  !
   ! in the store size and therefore points at something other than the store's own inertia.         !
   !                                                                                          !
   ! Deliberately a source constant, not a TOML knob: this is an unfinished feature, not a supported !
   ! configuration choice, and it must not look like one. Flip to 1.0 to resume the investigation.   !
   !=========================================================================================!
   real(wp), parameter :: TISSUE_STORE_SCALE = 1.0_wp

   !----- Prognostic-wood constants (retained; the split twin that these mirrored is retired). --------!
   !      build the same store from the same biomass. ----------------------------------------------!
   real(wp), parameter :: C2B_WOOD            = 2.0_wp  !< carbon -> biomass (carbon fraction 0.5)
   real(wp), parameter :: WOOD_MOIST_FRAC_ARK = 1.0_wp  !< [kg water/kg dry] fresh-sapwood moisture (MVP)


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
   subroutine column_be_stage(y, frozen, n, nsl, dt, y_out, niter, bf, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter    !< 1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux RATES for the ARK ledger
      !----- this stage's OWN surface tendencies (incl. transp_c), for the caller to b-weight into      !
      !      the plant water-mass update (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2) -- the SAME surf_tend the      !
      !      stage's own bf/CAS-source already used, so the mass debit and the CAS credit agree. --------!
      type(surface_tend_t), optional, intent(out) :: sf_out

      type(surface_state_t)      :: y_stage
      type(surface_frozen_t)     :: fs
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      real(wp)    :: t_ground, fliq1, wmass1, wcap, ccap, gah, gaw, gac
      real(wp)    :: enth1, shv1, e_infil, e_drain, e_clip, e_floor, t_cas1, qloss_total
      real(wp)    :: co21, enth_unused, shv_unused
      integer(ik) :: k, np, nfeval
      logical     :: ok

      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call uext_to_temp(y%soil_energy(1), wmass1, frozen%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      wcap = frozen%surf%wcap ; ccap = frozen%surf%ccap
      gah  = frozen%surf%gah  ; gaw  = frozen%surf%gaw ; gac = frozen%surf%gac
      fs = frozen%surf                    ! copied only to override the conductances below

      !----- Re-solve the Monin-Obukhov surface layer at THIS STAGE's canopy-air state, so the    !
      !      ventilation the stage is charged for is the ventilation its own temperature earns.    !
      !      Everything downstream (the BE commit, the Newton, and the boundary-flux ledger) reads !
      !      the LOCAL gah/gaw/gac, so refreshing them here keeps "one flux, both sides"           !
      !      automatically -- the state update and the ledger cannot disagree.                     !
      !                                                                                          !
      !      Mirror into fs and CLEAR ITS mo_live: surface_derivs would otherwise re-solve the     !
      !      surface layer on every one of the Newton's residual evaluations (up to 24 per stage)  !
      !      to fill a CAS tendency this scheme does not even read -- it commits the CAS through   !
      !      its own backward-Euler denominator.  So ARK pays for exactly ONE solve per stage. ----!
      call cas_conductances(frozen%surf, y%cas_enthalpy, y%cas_shv, gah, gaw, gac)
      fs%gah = gah ; fs%gaw = gaw ; fs%gac = gac
      fs%mo_live = .false.

      !----- CAS enthalpy + humidity. np==1: the uncoupled single-BE-pass baseline. np>1: a DIRECT 2x2  !
      !      Newton solve of the coupled backward-Euler surface block (the arrowhead). The FINAL surf_tend     !
      !      drives the soil sinks (single-flux-per-interface).                                         !
      if (np <= 1_ik) then
         y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
         call surface_derivs(y_stage, fs, t_ground, n, surf_tend)
         call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, fs, wcap, ccap, gah, gaw, gac, dt, &
                             enth1, shv1, co21)
      else
         call newton_surface_solve(y, fs, t_ground, n, dt, wcap, gah, gaw, enth1, shv1, surf_tend, nfeval, ok)
      end if
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      if (np > 1_ik) call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, fs, wcap, ccap, gah, gaw, gac, dt, &
                                         enth_unused, shv_unused, co21)   ! CO2 rides the same box
      y_out%cas_co2      = co21
      if (present(sf_out)) sf_out = surf_tend

      !----- soil-heat column: implicit BE-Thomas (soil_energy_step_implicit). ---------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top = surf_tend%g_top ; eforc%geothermal = frozen%geothermal
      !----- Root heat sink = qloss_total (uptake's advected enthalpy, sec 2/6, P2), the SAME sink     !
      !      column_derivs uses (meds_fast_time_derivs.f90), distributed by the static root_share       !
      !      profile: the soil pays once for the water the roots extract, the leaf/wood side gains it   !
      !      via qwflux_wl/q_wood_net, and the leaf pays the full vapour enthalpy of what it transpires. !
      !      (The old coh_qsoil proxy charged the soil a second time for that vapour's liquid part;      !
      !      2026-09 review, item 1A #10.) qloss_total sums to 0 when the P2 wiring is unset. ----------!
      qloss_total = sum(frozen%qloss_frozen(1:n))
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         !----- root sink + the two UNFACED post-solve mass corrections, valued at each layer's own    !
         !      state^n temperature in build_column_frozen. Sign: a SINK is positive-out, so the clip   !
         !      (water leaving layer k for the pond) ADDS and the theta_res floor (water created in     !
         !      layer k) SUBTRACTS. Both are 0 unless the hydrology actually corrected that layer. -----!
         eforc%root_heat_sink(k) = qloss_total * frozen%root_share(k)                                     &
                                 + frozen%clip_enth(k) - frozen%floor_enth(k)
         !----- INTERIOR advective faces (was hardcoded 0). Down-positive hydrology -> up-positive      !
         !      energy, same flip the split path applies. Without this the boundary enthalpy below has  !
         !      no path between layer 1 and the rest of the column -- see column_frozen_t%w_flux_frozen.!
         eforc%w_flux(k)         = -frozen%w_flux_frozen(k)
      end do
      !----- boundary water-enthalpy advection. The TOP face is now a KERNEL term with the same upwind  !
      !      rule and time level as the interior faces (the #71 fix, ported from the split path), not   !
      !      an ad-hoc layer-1 source. The BOTTOM face stays an explicit driver term at frozen%t_bot,      !
      !      matching the split -- it was never mis-timed, and re-basing it would mismatch the ledger.  !
      !      There is deliberately NO runoff term: runoff leaves the PONDING store, which holds mass    !
      !      but no enthalpy, so it has nothing to remove from soil layer 1 (it removed ~1 MJ/kg the    !
      !      layer never received). root_heat_sink is a SINK, so q_src = -sink/dz: add an outflow. ----!
      eforc%w_flux_top  = -frozen%infiltration / rho_h2o
      eforc%t_water_top = frozen%t_infil     ! #78 item 4: out of the pond
      eforc%w_flux_bot  = 0.0_wp
      e_infil = frozen%infiltration * internal_energy_liquid(frozen%t_infil)
      e_drain = frozen%drainage     * internal_energy_liquid(frozen%t_bot)
      e_clip  = sum(frozen%clip_enth(1:nsl))
      e_floor = sum(frozen%floor_enth(1:nsl))
      eforc%root_heat_sink(nsl) = eforc%root_heat_sink(nsl) + e_drain
      call soil_energy_step_implicit(se, eforc, frozen%therm, frozen%soil, frozen%energy_opts, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- soil water is OPERATOR-SPLIT OUT of the ESDIRK stages: theta is PASSED THROUGH (held at the   !
      !      stage input = theta^n) and the AUTHORITATIVE end-of-step theta is committed once, from the     !
      !      scratch column_hydrology_flux (frozen%theta1), in column_fast_step_ark. Re-solving it here with a !
      !      relief-free single-BE Richards drifted to saturation over long wet runs (no ponding/runoff),   !
      !      then hung the next scratch solve; the robust ponding/runoff/free-drain solve is the SOLE       !
      !      soil-water authority now (the ED2 "single soil-water authority" principle). theta feeds only   !
      !      the t_ground diagnosis + the soil-energy thermal property above, both correctly at theta^n. ---!
      y_out%theta(1:nsl) = y%theta(1:nsl)

      !----- pond PASSED THROUGH. y_out is intent(out), so without this it would default-initialise !
      !      to 0 rather than carry state^n -- harmless today (nothing reads it in a stage and it is   !
      !      excluded from the error norm) but wrong the moment #93 Phase 1 gives it a stage RHS. ----!
      y_out%w_surface      = y%w_surface
      y_out%w_surface_enth = y%w_surface_enth

      !----- plant water MASS PASSED THROUGH (advanced by advance_water_mass_full, not here). ----!
      y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)

      !----- emit this stage's boundary-flux RATES for the conservation ledger (§2.3). The b-weight   !
      !      + cross-substep accumulation happens in ark2_column_step / adaptive_ark_march. Every      !
      !      quantity is the committed-state flux, so the accumulated amounts telescope to closure.    !
      if (present(bf)) then
         t_cas1 = cas_temp_of_enthalpy(enth1, shv1)         ! committed CAS temp for the dew liquid enthalpy
         associate (fs2 => frozen%surf)
            bf%cas_enth_in  = surf_tend%src_enth + gah*fs2%enth_atm    ; bf%cas_enth_out = gah*enth1
            bf%cas_vap_in   = surf_tend%src_vap  + gaw*fs2%shv_atm     ; bf%cas_vap_out  = gaw*shv1
            bf%cas_co2_in   = fs2%nee_biotic + gac*fs2%co2_atm  ; bf%cas_co2_out  = gac*y_out%cas_co2
            bf%soil_enth_in = surf_tend%g_top + frozen%geothermal + e_infil + e_floor
            bf%soil_enth_out= qloss_total * sum(frozen%soil%root_frac(1:nsl)) + e_drain + e_clip
            !----- soil water is out of the ARK: its storage delta + q_top/drainage/uptake fluxes are     !
            !      re-sourced once/step from the frozen hflux in column_fast_step_ark, so the per-stage    !
            !      bf carries ONLY the CAS-vapour exchange (drainage/runoff/precip are frozen fast-step).  !
            bf%soil_wat_in  = 0.0_wp                            ; bf%soil_wat_out = 0.0_wp
            !----- condensation (dew) leaves the CAS as liquid at Tcas -> a whole-column water + liquid-   !
            !      enthalpy OUTPUT (the CAS-side loss is already in src_vap/src_enth, so cas_water/energy   !
            !      close automatically). -----------------------------------------------------------------!
            !----- ground_rad is the snowfac-BLENDED radiative input (= abs_sw+abs_lw when bare). ---!
            !----- #78 item 4: e_infil (pond -> soil) and e_clip (soil -> pond) are now transfers    !
            !      between two TRACKED stores, so they telescope and must NOT be boundary terms.     !
            !      The boundary precip input and the runoff output are added once at the outer level. !
            bf%whole_enth_in= surf_tend%coh_rnet + fs2%ground_rad + e_floor
            !----- row 1b: surf_tend%cond's enthalpy is NO LONGER a boundary loss -- the condensate is        !
            !      deposited into soil layer 1 by the caller, carrying this same u_liq(t_cas1). ------!
            bf%whole_enth_out= gah*(enth1 - fs2%enth_atm) + e_drain
            bf%whole_wat_in = 0.0_wp                            ; bf%whole_wat_out = gaw*(shv1 - fs2%shv_atm)
            bf%whole_cond   = surf_tend%cond                     ! row 1b: deposited into a store, not lost
            bf%whole_cond_enth = surf_tend%cond_enth   ! EXACTLY what surface_derivs debited from the CAS (one number, both sides)
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
   ! never error stops (GPU-safe). Commits the CAS via the FLUX form so budgets close for ANY surf_tend.      !
   !---------------------------------------------------------------------------------------!
   subroutine newton_surface_solve(y, fs, t_ground, n, dt, wcap, gah, gaw, enth1, shv1, surf_tend, nfeval, ok)
      type(column_state_t),   intent(in)    :: y
      type(surface_frozen_t), intent(in)    :: fs
      real(wp),               intent(in)    :: t_ground
      real(wp) :: co2_unused
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, wcap, gah, gaw
      real(wp),               intent(out)   :: enth1, shv1
      type(surface_tend_t),   intent(out)   :: surf_tend
      integer(ik),            intent(out)   :: nfeval
      logical,                intent(out)   :: ok

      type(surface_state_t) :: y_stage
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt
      integer(ik) :: it, ls
      real(wp),    parameter :: RTOL_N = 1.0e-7_wp, ATOL_H = 5.0e1_wp, ATOL_Q = 1.0e-6_wp
      real(wp),    parameter :: DETEPS = 1.0e-30_wp
      integer(ik), parameter :: NEWT_MAX = 4_ik, LS_MAX = 6_ik, FEVAL_CAP = 24_ik

      H0 = y%cas_enthalpy ; q0 = y%cas_shv
      Hk = H0 ; qk = q0 ; nfeval = 0_ik ; ok = .false.
      y_stage%cas_co2 = y%cas_co2
      y_stage%cas_enthalpy = Hk ; y_stage%cas_shv = qk
      call surface_derivs(y_stage, fs, t_ground, n, surf_tend) ; nfeval = nfeval + 1_ik
      R_H = wcap*(Hk - H0)/dt - surf_tend%src_enth - gah*(fs%enth_atm - Hk)
      R_q = wcap*(qk - q0)/dt - surf_tend%src_vap  - gaw*(fs%shv_atm  - qk)

      do it = 1_ik, NEWT_MAX
         if ( abs(R_H)*dt/wcap <= ATOL_H + RTOL_N*abs(Hk) .and.                                  &
              abs(R_q)*dt/wcap <= ATOL_Q + RTOL_N*abs(qk) ) then
            ok = .true. ; exit
         end if
         call jac_surface(Hk, qk, y%cas_co2, fs, t_ground, surf_tend, n, wcap, gah, gaw, dt, J11, J12, J21, J22, nfeval)
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
            y_stage%cas_enthalpy = Ht ; y_stage%cas_shv = qt
            call surface_derivs(y_stage, fs, t_ground, n, surf_tend) ; nfeval = nfeval + 1_ik
            RHt = wcap*(Ht - H0)/dt - surf_tend%src_enth - gah*(fs%enth_atm - Ht)
            Rqt = wcap*(qt - q0)/dt - surf_tend%src_vap  - gaw*(fs%shv_atm  - qt)
            if (RHt*RHt + Rqt*Rqt <= (1.0_wp - 1.0e-4_wp*lam)*rn0) exit          ! Armijo
            lam = 0.5_wp*lam
         end do
         Hk = Ht ; qk = qt ; R_H = RHt ; R_q = Rqt
         if (nfeval >= FEVAL_CAP) exit
      end do

      !----- authoritative final eval + flux-form commit (conservation holds for ANY surf_tend). -----------!
      y_stage%cas_enthalpy = Hk ; y_stage%cas_shv = qk
      call surface_derivs(y_stage, fs, t_ground, n, surf_tend) ; nfeval = nfeval + 1_ik
      call cas_box_commit(H0, q0, 0.0_wp, surf_tend, fs, wcap, 1.0_wp, gah, gaw, 0.0_wp, dt, enth1, shv1, co2_unused)
   end subroutine newton_surface_solve

   !---------------------------------------------------------------------------------------!
   ! The backward-Euler canopy-air box commit, routed through the SHARED kernel                      !
   ! meds_cas_biophysics%cas_column_step_implicit (which was exported but had no caller while this    !
   ! module re-implemented its three formulas inline). One implementation, both schemes' box.         !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_box_commit(h0, q0, c0, surf_tend, fs, wcap, ccap, gah, gaw, gac, dt, h1, q1, c1)
      real(wp),               intent(in)  :: h0, q0, c0, wcap, ccap, gah, gaw, gac, dt
      type(surface_tend_t),   intent(in)  :: surf_tend
      type(surface_frozen_t), intent(in)  :: fs
      real(wp),               intent(out) :: h1, q1, c1
      type(cas_column_t) :: box
      type(cas_source_t) :: src
      box%air_mass_capacity        = wcap ; box%air_molar_capacity     = ccap
      box%atm_conductance_enthalpy = gah  ; box%atm_conductance_vapor  = gaw ; box%atm_conductance_co2 = gac
      box%atm_enthalpy             = fs%enth_atm ; box%atm_specific_humidity = fs%shv_atm ; box%atm_co2 = fs%co2_atm
      src%surface_enthalpy_source  = surf_tend%src_enth ; src%surface_vapor_source = surf_tend%src_vap
      src%biotic_co2_source        = fs%nee_biotic
      call cas_column_step_implicit(h0, q0, c0, src, box, dt, h1, q1, c1)
   end subroutine cas_box_commit

   !----- 2x2 numerical Jacobian of (R_H, R_q) w.r.t. (H, q) by forward-differencing surface_derivs. --!
   subroutine jac_surface(Hk, qk, co2, fs, t_ground, surf_tend, n, wcap, gah, gaw, dt, J11, J12, J21, J22, nfeval)
      real(wp),               intent(in)    :: Hk, qk, co2, wcap, gah, gaw, dt
      type(surface_frozen_t), intent(in)    :: fs
      real(wp),               intent(in)    :: t_ground
      type(surface_tend_t),   intent(in)    :: surf_tend         ! base eval at (Hk,qk)
      integer(ik),            intent(in)    :: n
      real(wp),               intent(out)   :: J11, J12, J21, J22
      integer(ik),            intent(inout) :: nfeval
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq, dse_dH, dsv_dH, dse_dq, dsv_dq
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      y_stage%cas_co2 = co2
      y_stage%cas_enthalpy = Hk + dH ; y_stage%cas_shv = qk
      call surface_derivs(y_stage, fs, t_ground, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dH = (sfp%src_enth - surf_tend%src_enth)/dH ; dsv_dH = (sfp%src_vap - surf_tend%src_vap)/dH
      y_stage%cas_enthalpy = Hk ; y_stage%cas_shv = qk + dq
      call surface_derivs(y_stage, fs, t_ground, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dq = (sfp%src_enth - surf_tend%src_enth)/dq ; dsv_dq = (sfp%src_vap - surf_tend%src_vap)/dq
      J11 = wcap/dt + gah - dse_dH ; J12 =              - dse_dq
      J21 =              - dsv_dH  ; J22 = wcap/dt + gaw - dsv_dq
   end subroutine jac_surface
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
   ! advance_water_mass_full -- operator-split plant water MASS over the FULL dt from y%*_water_mass,   !
   ! driven by the CALLER-SUPPLIED transp_c_bw -- the SAME per-cohort transpiration (already b-weighted  !
   ! across the stage(s) that make up this step, sec below) the CAS's own vapour balance used, so the     !
   ! mass debit and the CAS credit are consistent to within the tableau's own stage algebra (not a         !
   ! separate, later evaluation at a possibly-different point).                                             !
   !                                                                                          !
   ! THE UPDATE IS A CORRECTOR, not the bare Euler step it used to be. The reason is exact algebra, not     !
   ! a tolerance judgement. solve_plant_water defines its own reported flux FROM its own storage change,    !
   !     flux%sapflow = dw_l/dt + transp,   flux%root_uptake = (dw_l + dw_w)/dt + transp,                   !
   ! so `W + dt*(sapflow - transp)` reproduces the kernel's matrix-exponential endpoint EXACTLY -- to        !
   ! machine precision, at any dt -- PROVIDED the transp used to price the sapflow is the transp actually    !
   ! debited. The pre-pass priced sapflow_frozen at transp_pp (the state^n FULL demand); this routine        !
   ! debits transp_c_bw (the b-weighted REALISED demand). The gap between the two is therefore not a          !
   ! discretisation error of the mass ODE at all -- it is exactly dt*(transp_pp - transp_bw), a pure FLUX     !
   ! INCONSISTENCY that grows linearly in dt and never converges away within a step. Measured on the 3 h      !
   ! midday fixture it was the WHOLE of the psi_leaf dt-divergence: 1.45 MPa of error at dt_fast = 900 s.     !
   !                                                                                          !
   ! Re-solving the kernel here on the REALISED transpiration -- same Category-0 frozen conductances and      !
   ! soil potentials as the pre-pass, so this is one semi-discretisation and not a second linearisation      !
   ! point -- makes the credit and the debit the same number again. Measured effect at dt_fast = 900 s:       !
   ! psi_leaf error 1.454 -> 0.0046 MPa (314x), against a common dt = 5 s reference that BOTH versions        !
   ! agree on to 8e-4 MPa (i.e. this changes the discretisation, NOT the dt->0 limit). T_cas differs by        !
   ! 6e-5 K at dt = 25 s and 0.007 K at 900 s -- also vanishing with dt, as a consistent scheme must.          !
   ! (MEDS_ED2_RK45_DESIGN.md sec 1/4/5 described the pre-corrector Euler form.) --------------------------!
   subroutine advance_water_mass_full(y, frozen, n, nsl, dt, transp_c_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: n, nsl
      real(wp),              intent(in)    :: dt
      real(wp),              intent(in)    :: transp_c_bw(n)  !< [kg/m2 ground/s] per-cohort, b-weighted over the step
      type(column_state_t),  intent(inout) :: y_out
      real(wp)    :: transp_i
      integer(ik) :: i
      !----- CORRECTOR locals (only touched on the re-solve path). --------------------------------!
      real(wp)    :: transp_pp(n), psi_c(N_HYDRO, n)
      real(wp)    :: sapflow_c(n), uptake_c(n), uptake_layer_c(nsl, n)
      real(wp)    :: psi_leaf_c(n), psi_wood_c(n), plc_c(n)
      integer(ik) :: nsub_c(n)
      logical     :: conv_c(n)
      real(wp)    :: sap_use(n), upt_use(n)

      !----- THE CORRECTOR (see the header). The pre-pass solved the kernel on transp_pp = the       !
      !      state^n FULL demand; the step actually debits transp_c_bw, the b-weighted REALISED       !
      !      demand. Since the kernel defines sapflow = dw_l/dt + transp, an Euler step on the mass   !
      !      reproduces the kernel's own dw EXACTLY when the two transpirations agree -- and differs   !
      !      by exactly dt*(transp_pp - transp_bw) when they do not. That residual is a pure FLUX      !
      !      INCONSISTENCY (credit priced at one transpiration, debit taken at another), it grows      !
      !      linearly in dt, and it is the whole of the psi_leaf dt-divergence. Re-solving here on     !
      !      the realised transpiration -- same FROZEN conductances and soil potentials, so this is    !
      !      still one Category-0 semi-discretisation, not a second linearisation point -- makes the   !
      !      credit and the debit the same number again.                                              !
      !                                                                                               !
      !      Guarded on allocation: the RK4/IMEX-Euler oracle hand-builds its own `frozen` and never      !
      !      populates the hydraulics boundary inputs, so it keeps the frozen fluxes (same rule as     !
      !      surface_frozen_t%mo_live). --------------------------------------------------------------!
      sap_use(1:n) = frozen%sapflow_frozen(1:n)
      upt_use(1:n) = frozen%uptake_frozen(1:n)
      if (allocated(frozen%rhizo_cond) .and. allocated(frozen%psi_soil_pre) .and. dt > tiny_num) then
         do i = 1_ik, n
            transp_pp(i) = transp_c_bw(i) / max(frozen%nplant(i), tiny_num)
         end do
         psi_c(NODE_LEAF, 1:n) = psi_from_water_content(y%leaf_water_mass(1:n),                      &
              frozen%hydro_p%leaf_pi0, frozen%hydro_p%leaf_elastic_mod, frozen%hydro_p%leaf_apoplast_frac,    &
              frozen%hydro_p%leaf_water_sat, frozen%bleaf(1:n))
         psi_c(NODE_WOOD, 1:n) = psi_from_water_content(y%wood_water_mass(1:n),                      &
              frozen%hydro_p%wood_pi0, frozen%hydro_p%wood_elastic_mod, frozen%hydro_p%wood_apoplast_frac,    &
              frozen%hydro_p%wood_water_sat, frozen%bsap(1:n) + frozen%broot(1:n))
         call solve_plant_water_batch(n, nsl, transp_pp(1:n), frozen%bleaf(1:n), frozen%bsap(1:n),         &
              frozen%broot(1:n), frozen%sap_area(1:n), frozen%height(1:n), frozen%leaf_area(1:n),                &
              frozen%psi_soil_pre(1:nsl), frozen%soil%z_node(1:nsl), frozen%rhizo_cond(1:nsl, 1:n),           &
              frozen%hydro_p, frozen%hydro_o, dt, psi_c(:, 1:n), sapflow_c(1:n), uptake_c(1:n),            &
              uptake_layer_c(1:nsl, 1:n), psi_leaf_c(1:n), psi_wood_c(1:n), plc_c(1:n),              &
              nsub_c(1:n), conv_c(1:n))
         !----- TAKE THE CORRECTOR'S SAPFLOW ONLY; the wood<->soil interface KEEPS uptake_frozen, which  !
         !      is the SAME number the soil column already committed as its root sink (frozen%uptake,        !
         !      rescaled by the soil's own fwilt-limited supply). That keeps "one flux, both sides"       !
         !      across an interface spanning two OPERATOR-SPLIT solves, so the whole-column water         !
         !      ledger closes exactly.                                                                    !
         !                                                                                          !
         !      MEASURED CONSEQUENCE, and why this is the right asymmetry rather than a shortcut.         !
         !      uptake_frozen is priced at the pre-pass transpiration, so at long dt_fast the plant       !
         !      receives ~1.3% less water than it should: plant total 1.9139 vs 1.9392 kg/m2 at 900 s     !
         !      against a 5 s reference. That TOTAL error is fixed by the committed uptake and is         !
         !      INDEPENDENT of the root boundary condition -- a specified-flux (Neumann) root BC was      !
         !      built and measured, and it produced the IDENTICAL total (1.9139). All the BC chooses is   !
         !      WHERE the deficit is parked:                                                              !
         !        potential BC (this code): all of it into the wood -- psi_leaf 4.6e-3, psi_wood 1.74e-1  !
         !        specified-flux BC:        split by capacitance -- psi_leaf 1.61e-1, psi_wood 1.60e-1    !
         !      The leaf is the consequential node (it drives the stomatal water-stress feedback, hence   !
         !      transpiration and GPP), while psi_wood at these potentials is far from wood_psi50 and     !
         !      drives almost no PLC. Parking the error in the wood is therefore the better trade, and    !
         !      the Neumann variant was NOT landed. Closing the deficit itself requires re-pricing the    !
         !      SOIL side (the pre-pass would have to see the realised transpiration); no plant-side      !
         !      boundary condition can do it. See <scratchpad>/neumann_root_bc.patch.                     !
         sap_use(1:n) = sapflow_c(1:n)
      end if

      do i = 1_ik, n
         transp_i = transp_c_bw(i) / max(frozen%nplant(i), tiny_num)
         !----- KNOWN DEFERRED EDGE CASE: unlike psi (whose PV-curve capacitance self-limits as        !
         !      tissue dries, dw/dpsi -> 0 in the flaccid tail), the mass ODE is a plain linear Euler    !
         !      step with no such restoring force -- sapflow_frozen/uptake_frozen are the STATE-n        !
         !      average, but transp_i here is the ENDPOINT-refreshed demand, so a large intra-step       !
         !      CAS swing could in principle debit more water than is actually in storage. Floor at a    !
         !      tiny positive value (not 0, so a subsequent psi_from_water_content diagnosis never        !
         !      divides by an exact 0 rwc) rather than let it go negative -- a bookkept boundary term      !
         !      for this clamp is deferred (mirrors the P1 surf_overflow precedent, not yet needed        !
         !      here: unobserved in this pass's test scenarios, see MEDS_ED2_RK45_DESIGN.md P2 notes). ---!
         y_out%leaf_water_mass(i) = max(y%leaf_water_mass(i) + dt*(sap_use(i) - transp_i), tiny_num)
         y_out%wood_water_mass(i) = max(y%wood_water_mass(i)                                          &
                                   + dt*(upt_use(i) - sap_use(i)), tiny_num)
      end do
   end subroutine advance_water_mass_full


   !---------------------------------------------------------------------------------------!
   ! advance_surf_water_full -- operator-split canopy-SURFACE water (sec 3.4, P2c) over the FULL dt   !
   ! from y%*_surf_water, mirroring advance_water_mass_full exactly: frozen%intercept_leaf/wood is the     !
   ! FROZEN Act-1 inflow, film_evap_*_bw is the b-weighted per-stage outflow (the SAME number the        !
   ! caller's sf2/sf3 credited to the CAS, sec 1/3/4/5's "one flux, both sides" discipline). Floors at    !
   ! exactly 0.0 (not tiny_num like the mass version -- surf_water feeds no psi/rwc division downstream).  !
   ! CAPACITY clamping + overflow bookkeeping (mirroring the split path's post-hoc treatment, sec 9)         !
   ! happens LATER, in column_fast_step_ark, not here -- this is the raw, unclamped Euler step only. --------!
   pure subroutine advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt
      real(wp),               intent(in)    :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      type(column_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      !----- NOT floored here (unlike an earlier version of this routine): the caller's post-hoc         !
      !      overflow/deficit bookkeeping (column_fast_step_ark, mirroring column_fast_step_rk45's own    !
      !      treatment) needs to SEE a possibly-negative value to correctly credit the shortfall back to    !
      !      the whole-column ledger's outflow -- clamping prematurely here would silently discard it,        !
      !      fabricating mass exactly like an unclamped store corrupts the NEXT step's frozen interception.  !
      do i = 1_ik, n
         y_out%leaf_surf_water(i) = y%leaf_surf_water(i) + dt*(frozen%intercept_leaf(i) - film_evap_leaf_bw(i))
         y_out%wood_surf_water(i) = y%wood_surf_water(i) + dt*(frozen%intercept_wood(i) - film_evap_wood_bw(i))
      end do
   end subroutine advance_surf_water_full


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
   !  advances one dt_fast with the ARK stepper, then unpacks. PARTIAL precip>0 guard-lift: the ARK   !
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
      type(surface_frozen_t) :: fs
      type(surface_tend_t)   :: surf_tend
      type(column_bflux_t)   :: acc, bfsub
      real(wp)    :: tg, fl, dt0, wcap, ccap, enth0, shv0, co20, enth1, shv1, co21, e_soil0, e_soil1, w_soil0, w_soil1
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
      !      base. The ARK commits the scratch column_hydrology_flux theta verbatim, so it inherits all  !
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
      if (cfg%ark_adaptive) then
         !----- section 8e WARM START: seed from the step the controller converged to on the PREVIOUS      !
         !      dt_fast call for this patch. Column stiffness barely changes call to call, so the cold     !
         !      start was re-discovering the same step size every call and paying ~1 rejection for it.     !
         !      An explicit [fast].ark_dt_init still wins (it is an operator override); with neither, the  !
         !      first call for a patch cold-starts exactly as before. ---------------------------------!
         dt0 = dt_fast
         if (biophys%adapt_dt_last > tiny_num) dt0 = min(biophys%adapt_dt_last, dt_fast)
         if (cfg%ark_dt_init  > tiny_num)  dt0 = min(cfg%ark_dt_init,   dt_fast)
         !----- The UNIFIED error-control bundle (§8c Layer 1): build_error_control seeds every tolerance !
         !      group from the setting that governs it today (and honours the [fast].rtol_all master      !
         !      dial), plus the controller + strictness. Defaults (CTRL_I, CTRL_L1, rtol_all unset)       !
         !      reproduce the legacy march byte-for-byte. ------------------------------------------------!
         ec = build_error_control(cfg)
         call adaptive_ark_march(y, frozen, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                                 niter=merge(NEWT_COUPLED, 1_ik, cfg%ark_coupled), acc=acc,          &
                                 dt_warm_out=dt_warm_next,                                          &
                                 clamp_n=budget%clamp_stage_n)
         biophys%adapt_dt_last = dt_warm_next
      else
         nsub = max(1_ik, cfg%ark_fixed_substep) ; nrej = 0_ik ; ycur = y ; call bflux_zero(acc, n)
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, frozen, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr,          &
                                  niter=merge(NEWT_COUPLED, 1_ik, cfg%ark_coupled), bf=bfsub,       &
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
      !      commit the AUTHORITATIVE end-of-step theta from the scratch column_hydrology_flux HERE, once,  !
      !      so a single consistent theta feeds the state commit, the soil_temp read-off, and BOTH the      !
      !      soil_water and whole_water storage terms (w_soil1 below). ------------------------------------!
      y_out%theta(1:nsl) = frozen%theta1(1:nsl)
      !----- pond: same treatment as theta -- committed from the scratch solve, but THROUGH the state !
      !      vector so there is one authority. Phase 1 replaces these two lines with a stage RHS. -----!
      y_out%w_surface      = frozen%w_surface1
      y_out%w_surface_enth = frozen%w_surface_enth1     ! #78 item 4: paired with the mass

      !----- §5.1 PROCESS MASK. The mask must mean the same thing under every scheme, so it is applied  !
      !      at the ARK's single state-commit point: a masked-off component is restored to state^n (y),  !
      !      leaving the ODE one dimension smaller while its couplings still acted during the march.     !
      !      This mirrors the split path's freeze exactly. mask%veg_energy needs no case here -- the ARK !
      !      error-stops on prognostic leaf/wood, so no vegetation energy store exists on this path. ----!
      if (.not. col_config%mask%cas_energy) y_out%cas_enthalpy        = y%cas_enthalpy
      if (.not. col_config%mask%cas_vapour) y_out%cas_shv             = y%cas_shv
      if (.not. col_config%mask%cas_co2)    y_out%cas_co2             = y%cas_co2
      if (.not. col_config%mask%soil_heat)  y_out%soil_energy(1:nsl)  = y%soil_energy(1:nsl)
      if (.not. col_config%mask%soil_water) then
         y_out%theta(1:nsl)   = y%theta(1:nsl)
         y_out%w_surface      = y%w_surface
         y_out%w_surface_enth = y%w_surface_enth
      end if
      if (.not. col_config%mask%hydraulics) then
         y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
         y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)
         !----- Canopy-SURFACE water (sec 3.4, P2c) rides the SAME hydraulics mask entry as internal   !
         !      water mass (mirrors meds_fast_rk45.f90's own choice; a dedicated mask field is deferred). !
         y_out%leaf_surf_water(1:n) = y%leaf_surf_water(1:n)
         y_out%wood_surf_water(1:n) = y%wood_surf_water(1:n)
      end if

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
         call clamp_canopy_film(y_out, col_cohort%lai, col_cohort%wai, col_config%hydro%dewmx, n, surf_overflow, surf_deficit)
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
      call uext_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        col_config%soil_thermal%soil_dry_heat_capacity(1), tg, fl)
      y_stage%cas_enthalpy = y_out%cas_enthalpy ; y_stage%cas_shv = y_out%cas_shv ; y_stage%cas_co2 = y_out%cas_co2
      call surface_derivs(y_stage, frozen%surf, tg, n, surf_tend)
      !----- Commit the tissue temperatures the frozen store already produced. surf_tend%leaf_temp/wood_temp !
      !      ARE the dt_fast endpoints: surface_derivs solved the balance with a_leaf/a_wood = cap/dt   !
      !      relaxing from frozen%surf%t_leaf0/t_wood0, so the store is already inside the CAS solve and    !
      !      there is nothing to correct afterwards. An earlier attempt applied the store as a POST-     !
      !      COMMIT lump on the committed CAS; that is unstable, because the tissue and the canopy air   !
      !      have comparable heat capacities (cap_wood/cap_cas ~ 0.5 here) and the tissue relaxes far    !
      !      inside one step, so dumping its demand on an already-solved CAS drove a growing oscillation !
      !      (measured: 49 K CAS kicks per step, wood diverging to 247 K against a 331 K balance).       !
      !      §5.1 mask: freezing veg_energy holds the tissue at its entry temperature while the fluxes   !
      !      it drives into the CAS are kept. -----------------------------------------------------------!
      tissue_store0 = 0.0_wp ; tissue_store1 = 0.0_wp
      do i = 1_ik, n
         !----- Derive the capacity from the SAME a_store the kernel relaxed against, not by         !
         !      recomputing it from dry_hcap + wmass. The two agree by construction today, but only    !
         !      this form guarantees that zeroing a_leaf/a_wood zeroes the ledger's store term too --  !
         !      i.e. that "no capacity" is a clean no-op end to end rather than a state change the     !
         !      fluxes never paid for. --------------------------------------------------------------!
         cap_leaf_a(i) = frozen%surf%a_leaf(i) * dt_fast
         cap_wood_a(i) = frozen%surf%a_wood(i) * dt_fast
         tissue_store0 = tissue_store0 + cap_leaf_a(i) * frozen%surf%t_leaf0(i)                          &
                                       + cap_wood_a(i) * frozen%surf%t_wood0(i)
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
         biophys%leaf_temp(1:n) = frozen%surf%t_leaf0(1:n) ; biophys%wood_temp(1:n) = frozen%surf%t_wood0(1:n)
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

      wcap = frozen%surf%wcap ; ccap = frozen%surf%ccap
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      enth1 = y_out%cas_enthalpy ; shv1 = y_out%cas_shv ; co21 = y_out%cas_co2
      e_pond1 = frozen%w_surface_enth1
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
      !      unlike w_plant0/1 above). Valued at u_liq(rain_temp) = frozen%surf%film_u_ref, the liquid       !
      !      enthalpy the intercepted water arrived with; the tissue pays enthalpy_vapor - film_u_ref per  !
      !      kg it evaporates (surface_derivs), so this store closes exactly against the CAS credit. All   !
      !      zero when canopy_water_on is off. -------------------------------------------------------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%rain_temp)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%rain_temp)
      intercept_total = sum(frozen%intercept_leaf(1:n) + frozen%intercept_wood(1:n))
      !----- L2/debug_error mode (col_config%energy%debug_error) promotes a non-closing budget from a       !
      !      silently-counted n_fail to a hard `error stop` -- the enforced half of the conservation   !
      !      check (plan MEDS_NUMERICS_SCOPING.md sec 4/QW2), mirroring the split path; off by         !
      !      default so production behaviour is unchanged. Each check reuses budget%*%resid, which        !
      !      budget_accumulate just set as a side effect. --------------------------------------------!
      !----- Tolerances are FLUX-scaled (meds_budget_check header): rtol * gross boundary flux over  !
      !      the step plus a rate floor * dt_fast. Store-scaled tolerances let a ~1 W/m2 leak through. !
      call budget_check(budget%cas_energy, wcap*enth0, wcap*enth1, acc%cas_enth_in, acc%cas_enth_out,     &
                        dt_fast, budget_energy_rate_floor, 'cas_energy (ark)', halt_budgets)
      call budget_check(budget%cas_water,  wcap*shv0,  wcap*shv1,  acc%cas_vap_in,  acc%cas_vap_out,      &
                        dt_fast, budget_water_rate_floor, 'cas_water (ark)', halt_budgets)
      call budget_check(budget%cas_co2,    ccap*co20,  ccap*co21,  acc%cas_co2_in,  acc%cas_co2_out,      &
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
      !      hflux, which closed its OWN mass budget to machine precision inside column_hydrology_flux. -----!
      !----- ...and the paired MASS, for the same reason (see the soil_energy note just above). ---------!
      !----- ...plus the scratch's own two post-solve corrections: clipped water LEFT the soil for the  !
      !      pond, floored water was CREATED in it -- both are in the committed theta1. -----------------!
      call budget_check(budget%soil_water,  w_soil0, w_soil1,                                              &
                        (frozen%q_top*rho_h2o + frozen%floor_mass)*dt_fast + cond_dep_mass,                    &
                        (frozen%drainage + frozen%uptake + frozen%clip_mass)*dt_fast,                             &
                        dt_fast, budget_water_rate_floor, 'soil_water (ark)', halt_budgets)
      !----- whole-WATER: precip IN; drainage + runoff + CAS-vapour OUT; ponding + plant internal water !
      !      MASS in the store. The soil + ponding + drainage/runoff/precip terms are frozen fast-step    !
      !      amounts; the CAS-vapour exchange gaw*(shv-shv_atm) is the ARK-accumulated part (acc%           !
      !      whole_wat_out). The plant's OWN w_plant0/1 store term (above) now absorbs the transp<->        !
      !      uptake mismatch the OLD comment here described (the ARK re-evaluates transpiration per          !
      !      ESDIRK stage as the CAS VPD evolves, while the committed soil theta lost the FROZEN scratch      !
      !      uptake_total) -- exactly the closure split's P0 already has, so the tolerance inflation that     !
      !      used to paper over that gap is gone; this closes to machine precision like the other 6 now. ----!
      !----- SNOW joins the ledgers (C4): the pack is a real mass + energy STORE, its accumulated      !
      !      precip enthalpy is a boundary INPUT, and SNOWFALL is a boundary water input -- forc%snowf  !
      !      was absent from w_in here (split has always carried it), so with a pack the ledger saw     !
      !      mass appear with no source and leaked exactly snowf*dt every step. Sublimation already     !
      !      leaves via the CAS vapour term and meltwater already moved pack -> soil as a paired        !
      !      transfer, so neither needs a term. All zero without snow. --------------------------------!
      !----- frozen%floor_mass: the theta_res floor's water is created inside the column and enters as a   !
      !      boundary INPUT, exactly as its enthalpy (e_floor -> acc%whole_enth_in) already did; the      !
      !      mass half was missing (2026-09 review, item 1A #5). -------------------------------------!
      call budget_check(budget%whole_water,                                                                &
                        w_soil0 + wcap*shv0 + w_surface0 + w_plant0 + surf_water0 + frozen%surf%snow_swe0,  &
                        w_soil1 + wcap*shv1 + frozen%w_surface1 + w_plant1 + surf_water1 + frozen%surf%snow_swe1, &
                        acc%whole_wat_in + (forc%precip + forc%snowf + biophys%shed_water_rate                &
                                            + frozen%floor_mass)*dt_fast,                                    &
                        acc%whole_wat_out + (frozen%runoff_surf + frozen%drainage)*dt_fast                      &
                                          + surf_overflow - surf_deficit,                                 &
                        dt_fast, budget_water_rate_floor, 'whole_water (ark)', halt_budgets)
      call budget_check(budget%whole_energy,                                                               &
                             !----- No melt rebase any more (#78 item 4): the pack hands its meltwater to  !
                             !      the POND, not to soil layer 1, so e_soil0 no longer contains the melt   !
                             !      enthalpy and the pack/pond pair telescopes on its own. -------------!
                        e_soil0                           + wcap*enth0 + surf_enth0                     &
                        + frozen%surf%snow_enth0 + e_pond0 + tissue_store0,                               &
                        e_soil1 + wcap*enth1 + surf_enth1 + frozen%surf%snow_enth1 + e_pond1                &
                        + tissue_store1,                                                                 &
                        acc%whole_enth_in + intercept_total*dt_fast*internal_energy_liquid(frozen%rain_temp) &
                                          + frozen%surf%snow_acc_enth                                       &
                                          + (frozen%precip_ground - frozen%surf%snow_melt_rate)*dt_fast         &
                                            * internal_energy_liquid(frozen%t_precip),                      &
                        acc%whole_enth_out + (surf_overflow - surf_deficit)*internal_energy_liquid(frozen%rain_temp) &
                                           + frozen%runoff_enth*dt_fast,                                    &
                        dt_fast, budget_energy_rate_floor, 'whole_energy (ark)', halt_budgets)

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_ark
   !----- The SHARED pre-pass (once per sub-step; ED2 freezes gs/hydraulics per DTLSM): refreshes the !
   !      aerodynamics + CAS-derived scalars from the current state, then computes LEAF gas exchange   !
   !      (GPP/gs/Rd), the FROZEN per-cohort leaf-energy coefficients h_coeff_f/g_tr_f, stem/root       !
   !      maintenance respiration, the NEE assembly, and the CAS capacities/atm-exchange conductances.  !
   !      ONE authority for both integrators -- this is what keeps split/ARK GPP bit-for-bit -- called   !
   !      from column_fast_step_ark / column_fast_step_rk45 once per dt_fast (the retired split path     !
   !      here) and from build_column_frozen below (which freezes these as explicit ARK macro-step       !
   !      inputs). `biophys` is intent(in): callers that need the CAS temperature persisted (the split)      !
   !      write biophys%cas%can_temp = tcas themselves right after the call.                                 !
   subroutine column_prepass(cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,                       &
                             tcas, qcas, press, rho, t_ground, h_coeff_f, g_tr_f,                       &
                             wcap, ccap, gah, gaw, gac, nee_biotic,                                     &
                             gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: col_config
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(in)    :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      real(wp),                intent(out)   :: tcas, qcas, press, rho, t_ground
      real(wp),                intent(out)   :: h_coeff_f(:), g_tr_f(:)
      real(wp),                intent(out)   :: wcap, ccap, gah, gaw, gac, nee_biotic
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      !----- OPTIONAL per-cohort DIAGNOSTIC capture (MEDS_IO_V01_PLAN.md section 3.4). Present only    !
      !      when the run reports per-cohort ecophysiology; absent, the extra leaf_flux_t fields are    !
      !      never even requested from the batch kernel, so this costs nothing.  ---------------------!
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) INSTANTANEOUS values

      !----- Bare-array batch I/O for the per-cohort physiology kernels (MEDS_NUMERICS_SCOPING.md).   !
      real(wp) :: par_arr(col_cohort%n), vpd_arr(col_cohort%n), gb_arr(col_cohort%n), rho_mol_arr(col_cohort%n), &
            psi_leaf_arr(col_cohort%n)
      real(wp) :: dmax_psi_arr(col_cohort%n), dmax_psi_seed
      real(wp) :: a_gross_arr(col_cohort%n), gs_arr(col_cohort%n), rd_arr(col_cohort%n)
      real(wp) :: a_net_arr(col_cohort%n), ci_arr(col_cohort%n), cs_arr(col_cohort%n), transp_arr(col_cohort%n)
      real(wp) :: bstom_arr(col_cohort%n), bnstom_arr(col_cohort%n)
      real(wp) :: stem_resp_arr(col_cohort%n), root_resp_arr(col_cohort%n)
      real(wp) :: e_air, gsw_ms, can_dmol
      real(wp) :: gpp, ra_leaf, ra_stem, ra_root, rh, soil_temp_root, theta_mean
      real(wp) :: xi(n_soil_pool), a_mat(n_soil_pool, n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      integer(ik) :: i, k, n, nsl

      n = col_cohort%n ; nsl = col_config%soil%n_active

      !----- aerodynamics from the current CAS state. -------------------------------------------!
      tcas = cas_temp_of_enthalpy(biophys%cas%can_enthalpy, biophys%cas%can_shv)
      qcas = biophys%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air ; t_ground = biophys%soil_e%soil_temp(1)
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = biophys%cas%can_co2
      aenv%t_ground = t_ground
      call aero_bottom_to_top(col_config%aero, aenv, ageom, n, col_cohort, biophys%leaf_temp, aero)

      !----- Root-weighted soil temperature + column-mean moisture (root / heterotrophic resp). !
      soil_temp_root = weighted_mean(biophys%soil_e%soil_temp(1:nsl), col_config%soil%root_frac, nsl)
      theta_mean = 0.0_wp
      do k = 1_ik, nsl
         theta_mean     = theta_mean     + biophys%soil_w%theta(k) * col_config%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-col_config%soil%soil_layer_z(nsl+1_ik), tiny_num)

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
         rho_mol_arr(i)  = press / (r_gas * biophys%leaf_temp(i))
         par_arr(i)      = forc%abs_par(i) / max(col_cohort%lai(i), 0.1_wp) * forc%par_per_w
         vpd_arr(i)      = max(sat_vapor_pressure(biophys%leaf_temp(i)) - e_air, 0.0_wp)
         gb_arr(i)       = aero%leaf_gbw(i) * rho_mol_arr(i)
         !----- psi_leaf for gs stays FROZEN (Category-0, ED2-faithful, section 12.6/Appendix A):     !
         !      diagnose ONCE per dt_fast from the prognostic leaf_water_mass^n (MEDS_ED2_RK45_       !
         !      DESIGN.md sec 4) -- do NOT refresh this per stage. -------------------------------------!
         psi_leaf_arr(i) = psi_from_water_content(biophys%leaf_water_mass(i), col_config%hydro_p%leaf_pi0,      &
              col_config%hydro_p%leaf_elastic_mod, col_config%hydro_p%leaf_apoplast_frac,                          &
              col_config%hydro_p%leaf_water_sat, col_cohort%bleaf(i))
      end do
      !----- STOMATAL WATER STRESS (issue #95). beta_stomata = min(1, exp(sref*psi)) is driven by     !
      !      YESTERDAY's daily-maximum leaf water potential -- the model's predawn potential -- which  !
      !      is what leaf_env_t%psi is documented to carry. Until this landed, psi was an    !
      !      OPTIONAL argument this driver never passed, so it defaulted to 0 and beta_stomata was     !
      !      IDENTICALLY 1: there was no stomatal water stress in the fast loop at all, and a plant    !
      !      would transpire at full rate with an empty wood store.                                    !
      !                                                                                          !
      !      DMAX_PSI_LEAF_UNSET (positive, so unmistakable -- a real leaf potential is <= 0) means the  !
      !      cohort has no history yet: a recruit, or the first step of a run. Seed it from the        !
      !      SURFACE-LAYER soil potential so it starts at its patch's actual water status rather than  !
      !      at 0, which would read as fully turgid. -----------------------------------------------!
      dmax_psi_seed = grav_head * soil_psi_from_theta(col_config%soil%retention, biophys%soil_w%theta(1),        &
                    col_config%soil%theta_sat(1), col_config%soil%theta_res(1), col_config%soil%vg_alpha(1),           &
                    col_config%soil%vg_n(1))
      do i = 1_ik, n
         if (col_cohort%dmax_psi_leaf(i) > 0.0_wp) then
            dmax_psi_arr(i) = dmax_psi_seed                    ! UNSET sentinel -> seed from the soil
         else
            dmax_psi_arr(i) = col_cohort%dmax_psi_leaf(i)
         end if
      end do
      !----- NOTE ON THE `psi=` KEYWORD: despite the name, what is passed is `dmax_psi_leaf` --  !
      !      the cohort's own predawn (daily-max) LEAF potential, NOT a soil potential. The leaf       !
      !      kernel's dummy is called psi because the Sabot stomatal limb is conventionally keyed !
      !      on soil/predawn potential, and predawn leaf psi IS the plant's overnight equilibration    !
      !      with the soil -- so the two coincide in WET soil. They do NOT coincide under drought:     !
      !      tau_w = C_wood/rhizo is ~9 s at theta 0.25 but ~4.8 DAYS at theta 0.10, which is exactly  !
      !      the regime this feedback exists for. Real soil potential enters here only as the seed for !
      !      a cohort with no history (dmax_psi_seed, above). Renaming the kernel dummy is deferred    !
      !      because `psi` is a published Python keyword (meds.plant.leaf) -- see issue #99. -----!
      if (present(cdiag)) then
         call leaf_gas_exchange_batch(n, par_arr, biophys%leaf_temp(1:n), vpd_arr, biophys%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, cfg, col_cohort%pft(1:n),                          &
                                      col_cohort%vcmax25(1:n), col_cohort%rd25(1:n), a_gross_arr, gs_arr, rd_arr,     &
                                      psi=dmax_psi_arr(1:n),                                            &
                                      a_net=a_net_arr, ci=ci_arr, cs=cs_arr, transp=transp_arr,         &
                                      beta_stom=bstom_arr, beta_nonstom=bnstom_arr)
         do i = 1_ik, n
            cdiag(CD_ANET,         i) = a_net_arr(i)
            cdiag(CD_AGROSS,       i) = a_gross_arr(i)
            cdiag(CD_GSW,          i) = gs_arr(i)
            cdiag(CD_GBW,          i) = aero%leaf_gbw(i)
            cdiag(CD_CI,           i) = ci_arr(i)
            cdiag(CD_CS,           i) = cs_arr(i)
            cdiag(CD_RD,           i) = rd_arr(i)
            cdiag(CD_TRANSP,       i) = transp_arr(i)
            cdiag(CD_BETA_STOM,    i) = bstom_arr(i)
            cdiag(CD_BETA_NONSTOM, i) = bnstom_arr(i)
            cdiag(CD_LEAF_TEMP,    i) = biophys%leaf_temp(i)
            cdiag(CD_WOOD_TEMP,    i) = biophys%wood_temp(i)
            cdiag(CD_LEAF_VPD,     i) = vpd_arr(i)
            cdiag(CD_PSI_LEAF,     i) = psi_leaf_arr(i)
            cdiag(CD_ABS_PAR,      i) = forc%abs_par(i)
            cdiag(CD_ABS_SW,       i) = forc%abs_sw(i)
            cdiag(CD_ABS_LW,       i) = forc%abs_lw(i)
            cdiag(CD_WIND,         i) = aero%wind(i)
            cdiag(CD_LEAF_WATER,   i) = biophys%leaf_water_mass(i)
            cdiag(CD_WOOD_WATER,   i) = biophys%wood_water_mass(i)
         end do
      else
         call leaf_gas_exchange_batch(n, par_arr, biophys%leaf_temp(1:n), vpd_arr, biophys%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, cfg, col_cohort%pft(1:n),                          &
                                      col_cohort%vcmax25(1:n), col_cohort%rd25(1:n), a_gross_arr, gs_arr, rd_arr,     &
                                      psi=dmax_psi_arr(1:n))
      end if
      !----- Elemental (§11): the array actuals drive the element-wise broadcast; `col_config%wood`/       !
      !      `col_config%root` (scalar PODs) and the patch-uniform `soil_temp_root` broadcast. -------------!
      call stem_maintenance_respiration(biophys%wood_temp(1:n), col_cohort%dbh(1:n), col_cohort%height(1:n),           &
                                   col_cohort%wai(1:n), col_cohort%nplant(1:n), col_config%wood, stem_resp_arr(1:n))
      call fine_root_maintenance_respiration(soil_temp_root, col_cohort%broot(1:n), col_config%root, root_resp_arr(1:n))
      do i = 1_ik, n
         gsw_ms  = gs_arr(i) / max(rho_mol_arr(i), tiny_num)
         gpp     = gpp     + a_gross_arr(i) * col_cohort%leaf_area(i) * col_cohort%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = a_gross_arr(i) * col_cohort%leaf_area(i)
         if (present(cdiag))   cdiag(CD_GPP_RATE, i) = a_gross_arr(i) * col_cohort%leaf_area(i)
         ra_leaf = ra_leaf + rd_arr(i)      * col_cohort%leaf_area(i) * col_cohort%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = rd_arr(i) * col_cohort%leaf_area(i)
         h_coeff_f(i) = sensible_heat_coeff(col_config%veg_thermal%effarea_heat * col_cohort%lai(i), aero%leaf_gbh(i), rho, cp_air)
         g_tr_f(i)    = leaf_transp_coeff(col_config%veg_thermal%effarea_transp, col_cohort%lai(i), aero%leaf_gbw(i), gsw_ms)
         ra_stem = ra_stem + stem_resp_arr(i) * col_cohort%nplant(i)
         ra_root = ra_root + root_resp_arr(i) * col_cohort%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = stem_resp_arr(i)
         if (present(root_resp_coh)) root_resp_coh(i) = root_resp_arr(i)
      end do

      !----- NEE = autotrophic (leaf Rd + stem + root) + heterotrophic Rh - GPP. Rh is EITHER the    !
      !      OLD constant-pool scalar form (soil_carbon_on = .false., bit-identical to before Part   !
      !      II) OR the matrix form over biophys%soil_carbon -- the FROZEN per-patch pool held constant    !
      !      across today's sub-steps (B2, MEDS_SLOW_DYNAMICS_DESIGN.md Part II section 9): the day's   !
      !      total fast Rh then equals the daily soil_carbon_step's pool debit BY CONSTRUCTION, since   !
      !      both read the same frozen pool + the same per-pool env scalar xi (accumulated into         !
      !      budget%xi_step below for the caller to integrate into xi_int). ------------------------------!
      if (cfg%soil_carbon_on) then
         call assemble_env_scalar(t_ground, soil_temp_root, theta_mean, col_config%soil%theta_res(1),      &
                                  col_config%soil%theta_sat(1), biophys%soil_carbon, cfg%soil_carbon, xi)
         call assemble_transfer_matrix(biophys%soil_carbon, cfg%soil_carbon, a_mat, k_diag, er)
         rh = heterotrophic_respiration_matrix(a_mat, k_diag, xi, biophys%soil_carbon)
         budget%xi_step = xi ; budget%rh_matrix_step = rh
      else
         rh = heterotrophic_respiration_flux(col_config%fast_soil_carbon, soil_temp_root, theta_mean,      &
                                             col_config%soil%theta_res(1), col_config%soil%theta_sat(1), col_config%co2)
      end if
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budget%gpp_last = gpp ; budget%nee_last = nee_biotic

      !----- CAS capacities + atm-exchange conductances (frozen across passes / the ARK macro-step). --!
      can_dmol = cas_molar_density(rho, qcas)
      wcap = rho      * biophys%cas%can_depth
      ccap = can_dmol * biophys%cas%can_depth
      gah  = rho      * aero%ustar * aero%temp1
      gaw  = rho      * aero%ustar * aero%temp2
      gac  = can_dmol * aero%ustar * aero%temp2
   end subroutine column_prepass
   !----- Build the frozen ARK inputs: the shared column_prepass above (leaf gas exchange /            !
   !      respiration / CAS caps / aero, bit-identical to the split) + this integrator's own per-cohort !
   !      geometry/radiation/wood packing + the frozen hydrology BCs, into a column_frozen_t; also      !
   !      packs the prognostic state into a column_state_t.                                             !
   subroutine build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, n, nsl, &
                                  frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: col_config
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      !----- intent(inout) since C4: the shared snow stage ADVANCES the pack (biophys%snow) and hands its  !
      !      melt enthalpy to biophys%soil_e%soil_energy(1) as a paired transfer, exactly as on the split   !
      !      path. Every other use of biophys here is still read-only. ------------------------------------!
      type(patch_biophys_t),   intent(inout) :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      integer(ik),             intent(in)    :: n, nsl
      type(column_frozen_t),   intent(out)   :: frozen
      type(column_state_t),    intent(out)   :: y
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) per-cohort diagnostic capture

      type(snow_stage_t)     :: snow_st   !< shared pre-column snow stage (C4, issue #76)
      type(chydro_forcing_t) :: hforc ; type(chydro_flux_t) :: hflux
      type(soil_column_t)    :: soil_w_scratch
      type(surface_state_t)  :: y_stage ; type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, nee_biotic, wcap, ccap, gah, gaw, gac
      integer(ik) :: i, k
      !----- Act-1 hydraulics pre-pass scratch (MEDS_ED2_RK45_DESIGN.md sec 1/3/5, P2): mirrors        !
      !      the original pre-pass order exactly (hydraulics BEFORE the soil solve, psi                   !
      !      diagnosed from state^n theta, the plant's own aggregate REQUEST becomes the soil's          !
      !      root-sink forcing, a post-hoc rescale if the soil can't honour it in full). --------------!
      real(wp) :: psi_soil_pre(nsl), psi_scratch(N_HYDRO, n), transp_pp(n)
      real(wp) :: sapflow_b(n), root_uptake_b(n), root_uptake_layer_b(nsl, n)
      real(wp) :: psi_leaf_b(n), psi_wood_b(n), plc_b(n)   !< batch outputs (unused downstream, complete SoA API)
      real(wp) :: rhizo_cond_all(nsl, n), k_theta_layer(nsl), total_uptake_b, scale, share_tot
      real(wp) :: t_up_wl, soil_temp_root, u_liq_soil, u_liq_up
      real(wp) :: sapflow_gnd(n), uptake_gnd(n)
      integer(ik) :: nsub_b(n)
      logical     :: converged_b(n)
      !----- Canopy-SURFACE water pre-pass scratch (MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c). ------------!
      real(wp) :: rain_above, combined_w, pai_i, throughfall_i, drip_i, throughfall_total
      real(wp) :: avail_leaf, avail_wood

      allocate(frozen%surf%h_coeff_f(n), frozen%surf%g_tr_f(n), frozen%surf%abs_sw(n), frozen%surf%abs_lw(n), frozen%surf%lai(n))
      allocate(frozen%surf%h_coeff_w(n), frozen%surf%abs_sw_wood(n), frozen%surf%abs_lw_wood(n), frozen%surf%wai(n))
      allocate(frozen%wood_dry_hcap(n), frozen%wood_wmass(n))
      allocate(frozen%leaf_dry_hcap(n), frozen%leaf_wmass(n))
      allocate(frozen%surf%a_leaf(n), frozen%surf%a_wood(n), frozen%surf%t_leaf0(n), frozen%surf%t_wood0(n))
      allocate(frozen%surf%qwflux_wl(n), frozen%surf%q_wood_net(n))
      !----- ZERO the advective-enthalpy terms AT ALLOCATION. They are only given their real values    !
      !      further down (after the plant-hydraulics batch supplies sapflow/uptake), but the sf0       !
      !      surface_derivs evaluations ABOVE that point already read them as veg_energy_balance's    !
      !      q_extra. Without this they are read UNINITIALISED -- undefined behaviour that injected       !
      !      whatever garbage the freshly-allocated heap happened to hold straight into the leaf energy   !
      !      balance. That is the origin of the 30-yr run's hang: garbage q_extra -> |dt_temp| ~1e64..     !
      !      1e292 -> absurd transp_c -> the hydraulics returned a non-finite root_uptake -> the soil-     !
      !      water adaptive loop spun forever on err = NaN (neither t nor nsub can advance). It also       !
      !      explains why the failure moved between builds and even between runs of the SAME binary: the   !
      !      garbage depends on prior heap contents, not on the physics. 0 is exactly the documented       !
      !      default for q_extra ("absent/0 for every caller but the P2 advective-enthalpy pre-pass"). ----!
      frozen%surf%qwflux_wl(1:n)  = 0.0_wp
      frozen%surf%q_wood_net(1:n) = 0.0_wp
      allocate(frozen%surf%f_wet_c(n), frozen%surf%g_film_f(n), frozen%surf%g_film_w(n))
      allocate(frozen%root_share(nsl), frozen%nplant(n), frozen%bleaf(n), frozen%bsap(n), frozen%broot(n),            &
               frozen%sap_area(n), frozen%height(n), frozen%leaf_area(n))
      allocate(frozen%sapflow_frozen(n), frozen%uptake_frozen(n), frozen%qloss_frozen(n))
      allocate(frozen%psi_soil_pre(nsl), frozen%rhizo_cond(nsl, n))
      allocate(frozen%intercept_leaf(n), frozen%intercept_wood(n))
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water(1:n) = biophys%leaf_surf_water(1:n) ; y%wood_surf_water(1:n) = biophys%wood_surf_water(1:n)

      !----- the SHARED pre-pass (column_prepass above): leaf gas exchange / respiration / CAS caps /   !
      !      aero -- writes directly into the frozen struct's h_coeff_f/g_tr_f arrays. ------------------!
      call column_prepass(cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,                       &
                          tcas, qcas, press, rho, t_ground, frozen%surf%h_coeff_f, frozen%surf%g_tr_f,      &
                          wcap, ccap, gah, gaw, gac, nee_biotic,                                      &
                          gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      !----- SHARED SNOW STAGE (C4, issue #76). AFTER column_prepass: it needs the CAS state           !
      !      (tcas/qcas), air properties (rho/press) and the aerodynamic conductance (aero%ggnet) that !
      !      the pre-pass establishes -- placing it earlier reads those undefined and yields a NaN     !
      !      pack. Split calls it after its own column_prepass for exactly this reason. Still BEFORE   !
      !      the hydrology forcing below, so meltwater reaches infiltration this step and the melt     !
      !      enthalpy is inside the soil column the state^n snapshot takes. No-op without a pack.      !
      call advance_snow_stage(col_config, forc, aero, biophys, dt_fast, tcas, qcas, rho, press, snow_st)
      frozen%surf%snowfac     = snow_st%snowfac   ; frozen%surf%h_snow       = snow_st%h_snow
      frozen%surf%le_snow     = snow_st%le_snow   ; frozen%surf%g_base_snow  = snow_st%g_base
      frozen%surf%subl_rate   = snow_st%subl_rate ; frozen%surf%ground_rad   = snow_st%ground_rad
      frozen%surf%snow_swe0   = snow_st%swe0      ; frozen%surf%snow_swe1    = snow_st%swe1
      frozen%surf%snow_enth0  = snow_st%enth0     ; frozen%surf%snow_enth1   = snow_st%enth1
      frozen%surf%snow_acc_enth = snow_st%acc_enth ; frozen%surf%snow_melt_rate = snow_st%melt_rate

      !----- Canopy INTERCEPTION (sec 3.4, P2c): frozen ONCE per dt_fast, frozen once per dt_fast:       !
      !      own "2c. CANOPY INTERCEPTION" sweep. ONE combined leaf+wood bucket per cohort, top-to-       !
      !      bottom over col_cohort's native height-DESCENDING gather order (the SAME direction the split path's  !
      !      own i=1..n loop already assumes is top-first), e_canopy=0 (capture/capacity only -- film       !
      !      evaporation is the SEPARATE, per-stage flux surface_derivs computes below from the frozen      !
      !      f_wet_c/g_film_f/g_film_w this block also sets). Converts the one-shot bucket update into an    !
      !      EQUIVALENT frozen RATE (matching sapflow_frozen/uptake_frozen's own "one frozen number, no      !
      !      per-stage re-solve" convention): integrating this rate by explicit Euler over dt_fast exactly    !
      !      reproduces the split's one-shot commit. Gated behind canopy_water_on so the untouched default    !
      !      path stays byte-identical (throughfall_total defaults to the FULL precip + snowf sum, matching   !
      !      the pre-P2c hforc%precip_ground line below unconditionally).                                    !
      !                                                                                                      !
      !      PLACEMENT (E-6, MEDS_INTEGRATOR_PARITY.md [RETIRED] sec 3e): this block MUST run after                     !
      !      advance_snow_stage, because it branches on snow_st%exists -- the pack OWNS the surface and        !
      !      liquid-rain interception is mutually exclusive with it, exactly as on the split path              !
      !      (the pack owns the surface, so liquid interception is skipped). It used to sit above            !
      !      column_prepass on the (correct, but insufficient) grounds that it has no aero dependency, and     !
      !      so read snow_st%exists BEFORE its only writer ran. snow_stage_t default-initialises exists to     !
      !      .false., so the read was defined rather than undefined -- but it was ALWAYS .false., i.e. the     !
      !      guard never fired. Nothing between here and the first surface_derivs call reads f_wet_c or        !
      !      intercept_leaf/wood, so moving the block down is otherwise inert.                                !
      !                                                                                                       !
      !      SNOWFALL (E-4): throughfall_total carries forc%precip + forc%snowf, as the split path always      !
      !      has. precip_phase splits the met precipitation into rain and snow WITHOUT consulting snow_on,     !
      !      so with the snow store off (the default) a sub-freezing step still delivers forc%snowf > 0.       !
      !      Omitting it here dropped that water on the floor -- it never reached hforc%precip_ground, while   !
      !      the whole-column ledger's w_in counted it (column_fast_step_ark / _rk45), so ARK and RK45 leaked  !
      !      exactly snowf*dt per step and diverged from split all winter (measured: split retained 39% more   !
      !      soil water over a January month; the gap vanished with snow_on = .true.). Under a pack this term  !
      !      is unused -- snow_accumulate has already taken BOTH precip and snowf, so precip_ground becomes    !
      !      meltwater only (the branch below). ---------------------------------------------------------------!
      frozen%intercept_leaf = 0.0_wp ; frozen%intercept_wood = 0.0_wp
      frozen%surf%f_wet_c = 0.0_wp ; frozen%surf%g_film_f = 0.0_wp ; frozen%surf%g_film_w = 0.0_wp
      throughfall_total = forc%precip + forc%snowf
      if (col_config%canopy_water_on .and. .not. snow_st%exists) then
         rain_above = forc%precip + forc%snowf
         do i = 1_ik, n
            pai_i      = col_cohort%lai(i) + col_cohort%wai(i)
            combined_w = biophys%leaf_surf_water(i) + biophys%wood_surf_water(i)
            call intercept_canopy_layer(combined_w, rain_above, col_cohort%lai(i), col_cohort%wai(i), 0.0_wp, dt_fast, &
                                        col_config%hydro%dewmx, col_config%hydro%intercept_k, col_config%hydro%intercept_alpha, &
                                        throughfall_i, drip_i, frozen%surf%f_wet_c(i))
            if (pai_i > tiny_num) then
               frozen%intercept_leaf(i) = (combined_w*col_cohort%lai(i)/pai_i - biophys%leaf_surf_water(i)) / dt_fast
               frozen%intercept_wood(i) = (combined_w*col_cohort%wai(i)/pai_i - biophys%wood_surf_water(i)) / dt_fast
            else
               frozen%intercept_leaf(i) = -biophys%leaf_surf_water(i) / dt_fast
               frozen%intercept_wood(i) = -biophys%wood_surf_water(i) / dt_fast
            end if
            rain_above = throughfall_i   ! cascades to the next (shorter) cohort
         end do
         throughfall_total = rain_above   ! whatever survives the shortest (last) cohort
      end if

      !----- per-cohort geometry + radiation + WOOD frozen inputs the ARK path needs (not shared with   !
      !      the split, which reads col_cohort%/forc% directly instead of packing a frozen struct). -----------!
      do i = 1_ik, n
         frozen%surf%lai(i)    = col_cohort%lai(i)
         frozen%surf%abs_sw(i) = forc%abs_sw(i) ; frozen%surf%abs_lw(i) = forc%abs_lw(i)
         !----- WOOD frozen inputs: real diagnostic values, or ZERO when wood is not diagnostic (so   !
         !      surface_derivs' wood branch is a no-op; prognostic wood is operator-split in P2). -----!
         !----- PROGNOSTIC-WOOD frozen set (Phase 4). Built unconditionally -- it is simply unread when  !
         !      wood is diagnostic -- so the two wood authorities never both feed the CAS: whichever mode !
         !      is active, exactly one of {frozen%surf's diagnostic inputs, frozen%wood_*} is non-zero. -------!
         !----- The two halves of the wood heat capacity take DIFFERENT masses. DRY tissue is ALL   !
         !      the wood (col_cohort%bwood) -- heartwood is dead structure but it still stores sensible    !
         !      heat, and a sapwood fraction defined on the bole under-counts branch wood, which is  !
         !      thin enough to be thermally active throughout. INTERNAL WATER is the sapwood ring    !
         !      ONLY (col_cohort%bsap): heartwood is taken as dry, which is what makes it heartwood. bsap   !
         !      is also the HYDRAULIC capacitance -- same ring, two consumers, consistently. --------!
         frozen%wood_dry_hcap(i) = max(col_cohort%bwood(i) * col_cohort%nplant(i) * C2B_WOOD * col_config%veg_thermal%c_sapw, &
                                    col_config%veg_thermal%veg_hcap_min)
         frozen%wood_wmass(i)    = col_cohort%bsap(i)  * col_cohort%nplant(i) * C2B_WOOD * WOOD_MOIST_FRAC_ARK
         !----- LEAF capacity, same two-part construction: dry leaf tissue + the INTERNAL (symplast) !
         !      water the hydraulics actually carries, biophys%leaf_water_mass [kg/plant] -> per m2 ground.!
         !      Read from biophys, NOT from y: y%leaf_water_mass is not filled until the very end of this  !
         !      routine, so reading it here would be an uninitialised read -- the same trap the        !
         !      advective-enthalpy zeroing note above records. ---------------------------------------!
         !      The intercepted FILM is deliberately NOT here: it is a separate store with its own     !
         !      phase (Step B), because folding it into a temperature-based capacity would commit to   !
         !      a film that can never freeze. -------------------------------------------------------!
         frozen%leaf_dry_hcap(i) = max(col_cohort%bleaf(i) * col_cohort%nplant(i) * C2B_WOOD * col_config%veg_thermal%c_leaf, &
                                    col_config%veg_thermal%veg_hcap_min)
         frozen%leaf_wmass(i)    = max(biophys%leaf_water_mass(i), 0.0_wp) * col_cohort%nplant(i)
         !----- The wood's diagnostic (zero-inertia) inputs are now filled UNCONDITIONALLY. They used  !
         !      to be zeroed whenever wood was "prognostic", because the prognostic store was a wholly  !
         !      separate operator-split solve that owned the wood's radiation and sensible flux. With   !
         !      the post-commit store (apply_tissue_store) there is only ONE wood solve: the diagnostic !
         !      balance here, relaxed afterwards by its own heat capacity. Zeroing these would leave    !
         !      apply_tissue_store relaxing towards a wood temperature of tcas with no absorbed         !
         !      radiation, i.e. silently deleting the wood's energy budget. -----------------------------!
         frozen%surf%wai(i)         = col_cohort%wai(i)
         frozen%surf%h_coeff_w(i)   = sensible_heat_coeff(pi * col_cohort%wai(i), aero%wood_gbh(i), rho, cp_air)
         frozen%surf%abs_sw_wood(i) = forc%abs_sw_wood(i)
         frozen%surf%abs_lw_wood(i) = forc%abs_lw_wood(i)
         !----- TISSUE STORE, frozen for the whole dt_fast. a = cap/dt_fast; the relaxation origin is  !
         !      the START-of-step tissue temperature. Every stage evaluation therefore returns the      !
         !      same dt_fast-averaged flux and dt_fast-endpoint temperature -- the store is an algebraic !
         !      closure, not a tableau DOF. --------------------------------------------------------------!
         frozen%surf%a_leaf(i)  = TISSUE_STORE_SCALE                                                   &
                               * (frozen%leaf_dry_hcap(i) + frozen%leaf_wmass(i) * cp_liq) / dt_fast
         frozen%surf%a_wood(i)  = TISSUE_STORE_SCALE                                                   &
                               * (frozen%wood_dry_hcap(i) + frozen%wood_wmass(i) * cp_liq) / dt_fast
         frozen%surf%t_leaf0(i) = biophys%leaf_temp(i)
         frozen%surf%t_wood0(i) = biophys%wood_temp(i)
         frozen%nplant(i)   = col_cohort%nplant(i)
         frozen%bleaf(i)    = col_cohort%bleaf(i)
         frozen%bsap(i) = col_cohort%bsap(i)
         frozen%broot(i)    = col_cohort%broot(i)   ; frozen%sap_area(i) = col_cohort%sap_area(i)
         frozen%height(i)   = col_cohort%height(i)  ; frozen%leaf_area(i) = col_cohort%leaf_area(i)
         !----- Canopy-SURFACE water film-evap conductances (sec 3.4, P2c): need aero%leaf_gbw/wood_gbw,   !
         !      so these run HERE (after column_prepass's aero solve above), not in the interception        !
         !      block before it. Gated behind canopy_water_on (not just "harmless via f_wet_c=0") per the    !
         !      P1 nvfortran lesson (this doc's own "P1 implementation notes": gate new ledger/diagnostic     !
         !      arithmetic behind its own flag from the start, rather than relying on it telescoping to a      !
         !      no-op). ------------------------------------------------------------------------------------!
         if (col_config%canopy_water_on) then
            frozen%surf%g_film_f(i) = leaf_film_coeff(col_config%veg_thermal%effarea_evap, col_cohort%lai(i), aero%leaf_gbw(i))
            frozen%surf%g_film_w(i) = leaf_film_coeff(col_config%veg_thermal%effarea_evap, col_cohort%wai(i), aero%wood_gbw(i))
         end if
      end do

      !----- the rest of the frozen surface inputs: CAS caps/conductances from column_prepass + atm     !
      !      state + NEE. ---------------------------------------------------------------------------!
      frozen%surf%leaf_emiss = col_config%veg_thermal%leaf_emiss
      frozen%surf%wcap = wcap ; frozen%surf%ccap = ccap
      frozen%surf%gah  = gah  ; frozen%surf%gaw  = gaw ; frozen%surf%gac = gac
      !----- Everything refresh_cas_conductances needs that is NOT the live CAS state. The         !
      !      geometry pair (displacement, roughness) is taken from the pre-pass's OWN aero output   !
      !      rather than recomputed, so a stage re-solve starts from the identical surface the      !
      !      state^n solve used, and the two agree exactly when the state has not moved. ----------!
      frozen%surf%aero_cfg     = col_config%aero
      frozen%surf%mo_u_ref     = aenv%u_ref     ; frozen%surf%mo_zref      = aenv%zref
      frozen%surf%mo_displace  = aero%displace  ; frozen%surf%mo_rough     = aero%rough
      frozen%surf%mo_theta_atm = aenv%theta_atm ; frozen%surf%mo_shv_atm   = aenv%shv_atm
      frozen%surf%mo_rho       = rho
      !----- ...and declare the inputs live, which is what licenses a per-stage re-solve. A bundle  !
      !      that has NOT been through here (a unit-test fixture, the RK4 oracle) leaves this false !
      !      and its supplied gah/gaw/gac are used verbatim. ----------------------------------------!
      frozen%surf%mo_live      = .true.
      frozen%surf%enth_atm = forc%enthalpy_atm ; frozen%surf%shv_atm = forc%shv_atm ; frozen%surf%co2_atm = forc%co2_atm
      frozen%surf%nee_biotic = nee_biotic
      frozen%surf%abs_sw_ground = forc%abs_sw_ground ; frozen%surf%abs_lw_ground = forc%abs_lw_ground
      frozen%surf%ggnet = aero%ggnet ; frozen%surf%rho = rho ; frozen%surf%press = press

      !----- params + hydraulics BCs. -----------------------------------------------------------!
      frozen%soil = col_config%soil ; frozen%therm = col_config%soil_thermal ; frozen%energy_opts = col_config%energy
      frozen%hydro_opts = col_config%hydro
      frozen%surf%cas_condensation = cfg%cas_condensation      ! §8g scheme-asymmetry guard
      frozen%geothermal = 0.0_wp

      !----- Act 1 (MEDS_ED2_RK45_DESIGN.md sec 1/3/5, P2): plant hydraulics runs BEFORE the soil     !
      !      solve, using psi diagnosed from state^n theta and the FULL transpiration demand -- no     !
      !      supply pre-throttle (the plant's own leaf/wood water MASS storage buffers any step-to-     !
      !      step soil-supply/demand mismatch instead; frozen%surf%src_frac stays at its 1.0 default).      !
      !      Hydraulics runs BEFORE the soil solve so the soil sees the realized uptake. ------------------!
      y_stage%cas_enthalpy = biophys%cas%can_enthalpy
      y_stage%cas_shv = biophys%cas%can_shv
      y_stage%cas_co2 = biophys%cas%can_co2
      call surface_derivs(y_stage, frozen%surf, t_ground, n, sf0)
      !----- Canopy-SURFACE water (sec 3.4, P2c): rescale the frozen film-evap conductance -- like        !
      !      uptake_frozen's own soil-limiting rescale above -- so a WORST-CASE potential evaporation       !
      !      over the FULL dt_fast (sf0's state^n film_evap, using the FULL unscaled g_film_f/w just         !
      !      computed) cannot drain MORE than will be available (current store + this step's frozen           !
      !      interception). Without this, film_evap (a state^n-frozen conductance, oblivious to               !
      !      depletion mid-step) can overdraw the store, driving leaf/wood_surf_water NEGATIVE -- which         !
      !      corrupts the NEXT step's frozen interception rate (intercept_canopy_layer's own internal            !
      !      floor silently "fixes" a negative starting bucket, fabricating mass that was never really            !
      !      lost) rather than genuinely closing the ledger. Re-evaluates sf0 afterward so the SAME             !
      !      call's transp_c (feeding the hydraulics batch solve below) sees the final, scaled conductance.      !
      !      Zero-effect when canopy_water_on is off (film_evap is 0 there, scale stays 1). -------------------!
      if (col_config%canopy_water_on) then
         do i = 1_ik, n
            avail_leaf = biophys%leaf_surf_water(i) + max(0.0_wp, frozen%intercept_leaf(i)) * dt_fast
            if (sf0%film_evap_leaf(i) > tiny_num) frozen%surf%g_film_f(i) = frozen%surf%g_film_f(i)          &
                 * min(1.0_wp, avail_leaf / (sf0%film_evap_leaf(i)*dt_fast))
            avail_wood = biophys%wood_surf_water(i) + max(0.0_wp, frozen%intercept_wood(i)) * dt_fast
            if (sf0%film_evap_wood(i) > tiny_num) frozen%surf%g_film_w(i) = frozen%surf%g_film_w(i)          &
                 * min(1.0_wp, avail_wood / (sf0%film_evap_wood(i)*dt_fast))
         end do
         call surface_derivs(y_stage, frozen%surf, t_ground, n, sf0)
      end if
      !----- UNITS: grav_head converts soil_psi_from_theta's METRES of head to the MPa the hydraulics   !
      !      seam expects. -----------------------------------------------------------------------------!
      psi_soil_pre(1:nsl) = grav_head * soil_psi_from_theta(col_config%soil%retention, biophys%soil_w%theta(1:nsl), &
           col_config%soil%theta_sat(1:nsl), col_config%soil%theta_res(1:nsl), col_config%soil%vg_alpha(1:nsl),          &
           col_config%soil%vg_n(1:nsl))
      !----- Per-layer rhizosphere conductance, UNCONDITIONAL since Phase 1 retired multilayer_roots.  !
      !      K(theta) is layer-only, so it is hoisted out of the cohort loop (hot-path work now). ------!
      do k = 1_ik, nsl
         k_theta_layer(k) = soil_hydr_cond_from_theta(col_config%soil%retention, biophys%soil_w%theta(k),        &
              col_config%soil%theta_sat(k), col_config%soil%theta_res(k), col_config%soil%vg_alpha(k),                   &
              col_config%soil%vg_n(k), col_config%soil%ksat(k))
      end do
      do i = 1_ik, n
         do k = 1_ik, nsl
            rhizo_cond_all(k, i) = rhizosphere_cond(rho_h2o*k_theta_layer(k)/grav_head, col_cohort%broot(i),  &
                 col_config%specific_root_area, col_config%soil%root_frac(k), col_config%soil%dz(k), col_cohort%nplant(i))
         end do
      end do
      !----- Hand the kernel's own frozen boundary inputs to the post-stage corrector, so it can       !
      !      re-solve on the SAME Category-0 coefficients this pre-pass used. -----------------------!
      frozen%psi_soil_pre(1:nsl)        = psi_soil_pre(1:nsl)
      frozen%rhizo_cond(1:nsl, 1:n)     = rhizo_cond_all(1:nsl, 1:n)
      frozen%hydro_p                    = col_config%hydro_p
      frozen%hydro_o                    = col_config%hydro_o
      psi_scratch(NODE_LEAF, 1:n) = psi_from_water_content(biophys%leaf_water_mass(1:n),                   &
           col_config%hydro_p%leaf_pi0, col_config%hydro_p%leaf_elastic_mod, col_config%hydro_p%leaf_apoplast_frac,      &
           col_config%hydro_p%leaf_water_sat, col_cohort%bleaf(1:n))
      psi_scratch(NODE_WOOD, 1:n) = psi_from_water_content(biophys%wood_water_mass(1:n),                   &
           col_config%hydro_p%wood_pi0, col_config%hydro_p%wood_elastic_mod, col_config%hydro_p%wood_apoplast_frac,      &
           col_config%hydro_p%wood_water_sat, col_cohort%bsap(1:n) + col_cohort%broot(1:n))
      transp_pp(1:n) = sf0%transp_c(1:n) / max(col_cohort%nplant(1:n), tiny_num)   ! [kg/plant/s] FULL demand
      call solve_plant_water_batch(n, nsl, transp_pp(1:n), col_cohort%bleaf(1:n),                             &
                                   col_cohort%bsap(1:n), col_cohort%broot(1:n), col_cohort%sap_area(1:n), &
                                         col_cohort%height(1:n),   &
                                   col_cohort%leaf_area(1:n),                                                  &
                                   psi_soil_pre(1:nsl), col_config%soil%z_node(1:nsl), rhizo_cond_all(1:nsl, 1:n), &
                                   col_config%hydro_p, col_config%hydro_o, dt_fast, psi_scratch(:, 1:n),                &
                                   sapflow_b(1:n), root_uptake_b(1:n), root_uptake_layer_b(1:nsl, 1:n),  &
                                   psi_leaf_b(1:n), psi_wood_b(1:n), plc_b(1:n), nsub_b(1:n), converged_b(1:n))
      budget%hydro_nsub    = sum(nsub_b(1:n))            ! section 5.3 work counter (same seam as split)
      budget%hydro_nonconv = count(.not. converged_b(1:n))
      !----- THRASH check (issue #104, plan E1b/E3). solve_plant_water takes 1.0-1.2 sub-steps per     !
      !      cohort in every physically ordinary state and ~136 once a wood store is pinned on its     !
      !      water floor. The threshold sits an order of magnitude above the ordinary band and an       !
      !      order below the collapsed one, so it fires on the pathology and nothing else. Under L2     !
      !      this is fatal -- a run that silently goes 13x slower AND returns floored stores is not a   !
      !      run anyone wants to keep. Under L1/L0 it is counted and reported through the existing      !
      !      work-counter output path (work_hydro_thrash_site). ------------------------------------!
      budget%hydro_thrash = merge(1_ik, 0_ik, budget%hydro_nsub > n * HYDRO_NSUB_THRASH)
      if (budget%hydro_thrash == 1_ik .and. cfg%error_level == CTRL_L2_STRICT)                          &
         error stop 'column_fast_step_ark: plant-hydraulics sub-stepping is pathological (issue #104) &
                    &-- a tissue store has almost certainly collapsed onto its water floor. Inspect &
                    &wood_water_mass; see docs/dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md sec 5c(v).'
      !----- HR (root efflux) intentionally NOT enabled anywhere in this model -- floor the aggregate  !
      !      (the same floor the hydrology kernel applies). ---------------------------------------------!
      root_uptake_b(1:n) = max(root_uptake_b(1:n), 0.0_wp)
      root_uptake_layer_b(1:nsl, 1:n) = max(root_uptake_layer_b(1:nsl, 1:n), 0.0_wp)
      !----- Per-cohort HYDRAULIC diagnostics, captured from the SAME solve the physics commits.    !
      !      psi_wood / plc were computed on every dt_fast and discarded; plc in particular is the    !
      !      one number that says whether a cohort is losing conductance, which nothing downstream    !
      !      of here could previously see. Taken AFTER the HR floor, so the reported uptake is the    !
      !      realized (non-negative) one the soil sink actually used.  ------------------------------!
      if (present(cdiag)) then
         cdiag(CD_SAPFLOW,     1:n) = sapflow_b(1:n)
         cdiag(CD_ROOT_UPTAKE, 1:n) = root_uptake_b(1:n)
         cdiag(CD_PSI_WOOD,    1:n) = psi_wood_b(1:n)
         cdiag(CD_PLC,         1:n) = plc_b(1:n)
      end if
      total_uptake_b = sum(root_uptake_b(1:n) * col_cohort%nplant(1:n))
      !----- THIS dt_fast's per-layer sink shares (Phase 1) -- the identical construction the split path  !
      !      does inline, so the three schemes place the root mass AND heat sink in the same layers. -----!
      frozen%root_share(1:nsl) = 0.0_wp
      do i = 1_ik, n
         do k = 1_ik, nsl
            frozen%root_share(k) = frozen%root_share(k) + root_uptake_layer_b(k, i) * col_cohort%nplant(i)
         end do
      end do
      share_tot = sum(frozen%root_share(1:nsl))
      if (share_tot > tiny_num) then
         frozen%root_share(1:nsl) = frozen%root_share(1:nsl) / share_tot
      else
         frozen%root_share(1:nsl) = col_config%soil%root_frac(1:nsl)
      end if

      !----- FROZEN hydrology BCs: the plant's OWN aggregate uptake REQUEST becomes the soil's root-   !
      !      sink forcing (not the raw transpiration demand), then a SCRATCH column_hydrology_flux      !
      !      for soil_evap / infiltration / uptake_total (the soil's TRUE, possibly fwilt-limited        !
      !      realized supply). throughfall_total (sec 3.4, P2c) is what SURVIVES the canopy -- the         !
      !      raw forc%precip when canopy_water_on is off (unchanged), or precip minus what the             !
      !      interception sweep above caught otherwise; feeding the soil the UN-reduced forc%precip         !
      !      here while ALSO crediting the intercepted share to the canopy surface store would create        !
      !      water from nothing (double-counted at the whole-column boundary). ------------------------!
      !----- PRECIP ROUTING under a pack (C4). snow_accumulate has                                      !
      !      ALREADY taken forc%snowf AND forc%precip into the pack, so ONLY meltwater may reach the   !
      !      ground -- adding throughfall on top double-counts the precip at the boundary. -----------!
      if (snow_st%exists) then
         hforc%precip_ground   = snow_st%melt_rate + biophys%shed_water_rate
      else
         hforc%precip_ground   = throughfall_total + biophys%shed_water_rate
      end if
      hforc%snow_free_frac     = 1.0_wp - snow_st%snowfac
      hforc%root_uptake(1:nsl) = total_uptake_b * frozen%root_share(1:nsl)
      hforc%t_ground           = t_ground ; hforc%q_air = qcas ; hforc%rho_air = rho
      !----- The hydrology kernel owns the ponding store's ENTHALPY too (#78 item 4): it needs each     !
      !      layer's temperature to value the saturation clip, and the temperature of the water         !
      !      entering the pond. Under a pack that is the MELTWATER temperature, not frozen%rain_temp --     !
      !      rain_temp is pinned to tsupercool_liq so the ledger books no boundary input for melt. -----!
      hforc%soil_temp(1:nsl)   = biophys%soil_e%soil_temp(1:nsl)
      !----- Temperature that VALUES the ground inflow. Under a pack it is the meltwater's. On bare      !
      !      ground it is the EFFECTIVE liquid temperature of the rain + sub-threshold-snowfall mixture:  !
      !      rain arrives as liquid at the canopy-air temperature, snow as ICE at min(t_3ple, tair) --   !
      !      the same valuation snow_accumulate gives snowfall that does form a pack -- and the mixture   !
      !      enthalpy per kg is expressed through temp_of_liquid_enthalpy (exact inverse of              !
      !      internal_energy_liquid; below t_3ple it represents water that must still melt, which the    !
      !      pond/soil plateau then does with soil heat). Valuing the snow as liquid at tcas, as this     !
      !      used to, created the fusion enthalpy L_f per kg of sub-threshold snow at the boundary        !
      !      (ledger-consistent, physically wrong; 2026-09 review). -----------------------------------!
      hforc%t_precip = tcas
      if (snow_st%exists) then
         hforc%t_precip = snow_st%t_melt
      else if (forc%precip + forc%snowf > tiny_num) then
         hforc%t_precip = temp_of_liquid_enthalpy(                                                    &
              (forc%precip * internal_energy_liquid(tcas)                                            &
               + forc%snowf * internal_energy_ice(min(t_3ple, forc%tair))) / (forc%precip + forc%snowf))
      end if
      !----- Bare-soil aerodynamic resistance, AREA-weighted by the snow-free fraction set above. This !
      !      path used to pin snow_free_frac at 1.0 because it modelled no snow at all; C4's shared     !
      !      stage removed that limitation, so the weighting is real here now. ------------------------!
      hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)
      soil_w_scratch = biophys%soil_w
      call column_hydrology_flux(soil_w_scratch, hforc, col_config%soil, col_config%hydro, dt_fast, hflux)
      budget%soil_nsub = hflux%nsub                 ! section 5.3 work counter (same seam on both schemes)
      frozen%surf%soil_evap = hflux%soil_evap
      frozen%q_top          = (hflux%infiltration - hflux%soil_evap) / rho_h2o

      !----- Soil-limiting rescale (MEDS_ED2_RK45_DESIGN.md sec 3): scale down ONLY the credit         !
      !      applied to wood_water_mass (not sapflow, an internal wood->leaf transfer) so the whole-    !
      !      column ledger closes to the soil's TRUE realized supply -- scale==1 exactly in the         !
      !      common (non-limited) case. sapflow_frozen/uptake_frozen are what column_derivs' mass ODE    !
      !      reads directly; frozen%uptake (== hflux%uptake_total) is ALSO the soil-water tendency's root    !
      !      sink authority in column_derivs, so the two sides of the wood<->soil interface use the       !
      !      identical number by construction. ---------------------------------------------------------!
      scale = 1.0_wp
      if (total_uptake_b > tiny_num) scale = min(1.0_wp, hflux%uptake_total / total_uptake_b)
      frozen%sapflow_frozen(1:n) = sapflow_b(1:n)
      frozen%uptake_frozen(1:n)  = root_uptake_b(1:n) * scale

      !----- ADVECTIVE ENTHALPY (MEDS_ED2_RK45_DESIGN.md sec 2/6, P2 -- ED2's qwflux_wl/qloss): both     !
      !      the mass FLUX (sapflow_b/uptake_frozen, already frozen above) AND the upwind reference        !
      !      TEMPERATURE are frozen at state^n here -- computing the reference temperature from the        !
      !      CURRENT (this-stage) leaf/wood/soil temperature would be circular (surface_derivs hasn't       !
      !      run yet this call, and it's what PRODUCES those temperatures). Freezing both still closes      !
      !      the energy ledger exactly (the same "one number, both sides" principle sapflow_frozen/          !
      !      uptake_frozen already rely on for water) -- it trades some fidelity in the thermal upwind       !
      !      choice for tractability, a documented, bounded approximation (see the design doc P2 notes).     !
      !      HR is disabled project-wide (uptake floored >=0), so qloss's upwind is unconditionally the       !
      !      root-frac-weighted mean soil temperature (weighted_mean is a generic weighted-sum, not      !
      !      psi-specific, so it is reused verbatim for temperature here). -------------------------------!
      soil_temp_root = weighted_mean(biophys%soil_e%soil_temp(1:nsl), col_config%soil%root_frac, nsl)
      u_liq_soil = internal_energy_liquid(soil_temp_root)
      do i = 1_ik, n
         t_up_wl   = merge(biophys%wood_temp(i), biophys%leaf_temp(i), sapflow_b(i) >= 0.0_wp)
         u_liq_up  = internal_energy_liquid(t_up_wl)
         sapflow_gnd(i) = frozen%sapflow_frozen(i) * col_cohort%nplant(i)   ! [kg/m2 ground/s]
         uptake_gnd(i)  = frozen%uptake_frozen(i)  * col_cohort%nplant(i)   ! [kg/m2 ground/s]
         frozen%surf%qwflux_wl(i)  = sapflow_gnd(i) * u_liq_up
         frozen%qloss_frozen(i)    = uptake_gnd(i)  * u_liq_soil
         frozen%surf%q_wood_net(i) = frozen%qloss_frozen(i) - frozen%surf%qwflux_wl(i)
      end do

      !----- FROZEN boundary hydrology for the guard-lift: the rain/drainage/runoff water-enthalpy       !
      !      advection (state^n temps, matching the split) + the scratch's end-of-step ponding/aquifer/  !
      !      water-table (soil_w_scratch was advanced in place by column_hydrology_flux). ---------------!
      frozen%infiltration = hflux%infiltration ; frozen%drainage    = hflux%drainage
      frozen%precip_ground = hforc%precip_ground
      frozen%runoff_surf  = hflux%runoff_surf
      !----- rain_temp = tsupercool_liq under a pack makes internal_energy_liquid vanish, so meltwater !
      !      infiltrates its MASS at zero enthalpy -- the enthalpy already moved, paired, inside        !
      !      advance_snow_stage. Without this the melt energy is counted twice at soil layer 1. -------!
      frozen%rain_temp = hforc%t_precip                    ! one valuation for the boundary inflow AND the film
      if (snow_st%exists) frozen%rain_temp = tsupercool_liq
      frozen%surf%film_u_ref = internal_energy_liquid(frozen%rain_temp)   ! what the film is valued at (surface_derivs)
      !----- The infiltrating water comes OUT OF THE POND, so the soil top-face advection is        !
      !      referenced to the pond temperature the kernel just reported (#78 item 4). -----------!
      frozen%t_infil = hflux%t_infil
      frozen%runoff_enth = hflux%runoff_enth
      frozen%uptake       = hflux%uptake_total
      !----- Interior face fluxes + the post-solve mass corrections, from the SAME scratch solve. The   !
      !      clip/floor enthalpies are valued HERE, at each layer's state^n temperature, because that   !
      !      is the temperature the correction has to be neutral against -- and it is the only place    !
      !      the per-layer soil temperature is in scope. ---------------------------------------------!
      frozen%w_flux_frozen(1:nsl) = hflux%w_flux(1:nsl)
      do k = 1_ik, nsl
         frozen%clip_enth(k)  = hflux%clip_layer(k)  * internal_energy_liquid(biophys%soil_e%soil_temp(k))
         frozen%floor_enth(k) = hflux%floor_layer(k) * internal_energy_liquid(biophys%soil_e%soil_temp(k))
      end do
      frozen%clip_mass  = sum(hflux%clip_layer(1:nsl))
      frozen%floor_mass = sum(hflux%floor_layer(1:nsl))
      frozen%t_bot        = biophys%soil_e%soil_temp(nsl)
      frozen%w_surface1   = soil_w_scratch%w_surface
      frozen%w_surface_enth1 = soil_w_scratch%w_surface_enth
      frozen%t_precip        = hforc%t_precip
      !----- the AUTHORITATIVE committed soil moisture: soil_w_scratch was advanced IN PLACE by the robust  !
      !      column_hydrology_flux above, so its theta IS the end-of-step (relieved) soil water. -----------!
      allocate(frozen%theta1(nsl))
      frozen%theta1(1:nsl) = soil_w_scratch%theta(1:nsl)

      !----- pack the prognostic state: plant water MASS is now NATIVE (MEDS_ED2_RK45_DESIGN.md sec 4, !
      !      P2) -- a direct copy from the persisted state, no psi round-trip needed any more. ----------!
      y%cas_enthalpy = biophys%cas%can_enthalpy ; y%cas_shv = biophys%cas%can_shv ; y%cas_co2 = biophys%cas%can_co2
      y%soil_energy(1:nsl) = biophys%soil_e%soil_energy(1:nsl)
      y%theta(1:nsl)       = biophys%soil_w%theta(1:nsl)
      !----- pond packed onto the state vector (#93 Phase 0). Still committed from the scratch      !
      !      hydrology below, so this is carriage only -- no behaviour change. --------------------!
      y%w_surface          = biophys%soil_w%w_surface
      y%w_surface_enth     = biophys%soil_w%w_surface_enth
      y%leaf_water_mass(1:n) = biophys%leaf_water_mass(1:n)
      y%wood_water_mass(1:n) = biophys%wood_water_mass(1:n)
   end subroutine build_column_frozen
   !----- Solve canopy aerodynamics with the cohort order it CONTRACTS for -- BOTTOM(1)->TOP(n)  !
   !      -- from the height-DESCENDING column buffer. Only the wind cascade + the per-cohort       !
   !      boundary layers depend on order; the whole-canopy scalars (ustar/temp1/temp2/uh) do not.  !
   !      An ascending-height permutation `ord` reverses the per-cohort inputs; the per-cohort wind  !
   !      and leaf/wood conductance outputs are scattered back to gather order. Identity for n<=1,   !
   !      so single-cohort behaviour is bit-unchanged.                                               !
   subroutine aero_bottom_to_top(acfg, aenv, ageom, n, col_cohort, leaf_temp, aero)
      type(aero_cfg_t),      intent(in)    :: acfg
      type(aero_env_t),      intent(in)    :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      integer(ik),           intent(in)    :: n
      type(column_cohort_t), intent(in)    :: col_cohort
      real(wp),              intent(in)    :: leaf_temp(:)
      type(aero_out_t),      intent(inout) :: aero
      integer(ik) :: ord(n), k, j, imin
      real(wp)    :: hmin
      logical     :: used(n), descending
      real(wp)    :: h_bt(n), lai_bt(n), cr_bt(n), lt_bt(n), lw_bt(n), bd_bt(n)
      real(wp)    :: wind_bt(n), lgbh_bt(n), lgbw_bt(n), wgbh_bt(n), wgbw_bt(n)

      !----- ord(k) = gather index of the k-th cohort counting from the canopy BOTTOM.            !
      !                                                                                          !
      !      This used to be an O(n^2) selection sort, and it is the ONLY superlinear term in the !
      !      whole fast loop -- run once per dt_fast on ALL THREE schemes (split reaches it via   !
      !      column_prepass, ark/rk45 via build_column_frozen). Measured in isolation over a      !
      !      6-day run: 0.97 s at n = 2000, 3.9 s at n = 4000, 16.5 s at n = 8000. Clean n^2, and !
      !      beyond n ~ 4000 it dominates the fast loop outright.                                 !
      !                                                                                          !
      !      It is also REDUNDANT in the normal case: sort_cohorts leaves the cohort block        !
      !      height-DESCENDING and the fast-loop gather preserves that order, so bottom-to-top is !
      !      simply the reverse, ord(k) = n-k+1. Detect that in O(n) and take the reverse; fall   !
      !      back to the original sort otherwise, because unit tests construct cohorts in         !
      !      arbitrary order and this routine must stay correct for them.                          !
      !                                                                                          !
      !      TIE-BREAK, and why the two branches agree exactly: the sort's `<=` keeps the LAST    !
      !      index achieving the running minimum, so among equal heights it emits the largest     !
      !      index first. In a descending array equal heights are consecutive and the largest     !
      !      remaining index is always minimal, so the reverse produces the identical permutation !
      !      -- ties included. This is bit-identical, not merely equivalent. -------------------!
      descending = .true.
      do j = 1_ik, n - 1_ik
         if (col_cohort%height(j) < col_cohort%height(j+1_ik)) then ; descending = .false. ; exit ; end if
      end do

      used = .false.
      do k = 1_ik, n
         if (descending) then
            imin = n - k + 1_ik                                  ! O(n) fast path
         else
            imin = 0_ik ; hmin = huge(1.0_wp)                    ! O(n^2) fallback (unsorted input)
            do j = 1_ik, n
               if (.not. used(j) .and. col_cohort%height(j) <= hmin) then ; hmin = col_cohort%height(j) ; imin = j ; end if
            end do
         end if
         ord(k)    = imin ; used(imin) = .true.
         h_bt(k)   = col_cohort%height(imin)     ; lai_bt(k) = col_cohort%lai(imin)
         cr_bt(k)  = col_cohort%crown(imin)      ; lt_bt(k)  = leaf_temp(imin)
         lw_bt(k)  = col_cohort%leaf_width(imin) ; bd_bt(k)  = col_cohort%branch_diam(imin)
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
