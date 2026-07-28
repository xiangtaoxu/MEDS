!==========================================================================================!
! meds_fast_split -- the operator-split (+ Picard) fast-loop stepper. column_fast_step is the   !
! entry point the fast-tier orchestrator (meds_fast_dynamics) calls once per dt_fast sub-step;   !
! it is ALSO the two-scheme DISPATCH gate: on cfg%time_integrator=="ark" it redirects to          !
! meds_fast_ark%column_fast_step_ark (a peer module, called cross-module) BEFORE the split path   !
! below is entered, so the golden-anchored split body is byte-for-byte unentered on that branch.   !
!                                                                                          !
! The split path couples the fast biophysics processes into one sub-daily operator-split sweep    !
! (design Part II). The canopy-air-space (CAS) twins are the shared coupling reservoir: each        !
! dt_fast the aerodynamics kernel sets the conductances; the leaf and GROUND surface fluxes feed    !
! the CAS; the soil WATER column (Richards) and soil THERMAL column (heat) are advanced; and the    !
! three CAS twins (enthalpy / specific humidity / CO2) are advanced IMPLICITLY in the atmospheric-  !
! exchange term, gated by u*.                                                                       !
!                                                                                          !
! Coupled so far: aerodynamics -> {diagnostic leaf balance, ground balance} -> soil WATER column     !
! (infiltration / DSL soil-evap / root uptake / drainage) -> CAS enthalpy/vapour/CO2 (implicit) ->    !
! soil THERMAL column (implicit BE-Thomas), which reads the just-updated soil moisture. Carries the    !
! §3.5 fix (atm<->CAS conductance = profile-factored rho*ustar*temp1/temp2) AND the §3.6 GROUND-        !
! evaporation single-authority: the hydrology kernel's DSL/alpha_soil soil_evap is THE ground latent    !
! flux -- it drives the CAS vapour twin AND the ground energy balance's LE (no double-count).           !
!                                                                                          !
! Leaf temperature is DIAGNOSED from a linearized steady-state balance. CARBON is now REAL: leaf gas     !
! exchange gives GPP + stomatal conductance (driving transpiration) + leaf Rd; stem + fine-root           !
! maintenance respiration (autotrophic) and heterotrophic Rh assemble NEE = (Rd+stem+root) + Rh - GPP     !
! for the CAS CO2 twin. Plant HYDRAULICS is now REAL too: each cohort carries prognostic node water       !
! potentials psi(NODE_LEAF/WOOD) in patch_biophys_t; solve_plant_water advances them from the realized     !
! (supply-limited) transpiration demand and the root-weighted soil psi, and the updated psi_leaf feeds     !
! NEXT step's leaf gas exchange -- the soil -> plant -> stomata drought feedback (lagged one dt_fast).     !
!                                                                                          !
! WHOLE-COLUMN CONSERVATION -- verified by budg%whole_energy / budg%whole_water (Δ of ALL stores vs the     !
! true boundary fluxes; these CATCH cross-seam leaks the per-kernel budgets miss). Water-borne enthalpy    !
! is transported consistently: the CAS latent uses enthalpy_vapor(tl) (matching the CAS inverter + ground); !
! the soil sheds the transpiration water's liquid enthalpy via root_heat_sink; infiltration/drainage water  !
! carry internal_energy_liquid across the soil boundaries. INTER-LAYER advective heat: the hydrology kernel  !
! now EXPOSES the time-mean per-face Darcy flux (hflux%w_flux), and soil_energy_step_implicit can advect the liquid   !
! enthalpy on it. cfg%advect_soil_heat is ON by default and is a CONSERVATION REQUIREMENT, not an accuracy   !
! knob: the boundary faces inject/remove liquid enthalpy at layer 1 and layer nsl, so lumping the INTERIOR   !
! faces leaves those boundary terms with no path between them. Under saturation (water entering the top and  !
! being clipped out of a deeper layer) that mis-places the whole boundary flux and runs the surface layer    !
! away -- test_column_dynamics RUN 7 reaches 361 K with it forced OFF, 297 K with it ON. Nothing outside     !
! that test sets the flag, so the OFF path is not reachable in production.                                   !
!==========================================================================================!
module meds_fast_split
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, cp_air, latent_heat_vap, rho_h2o, pi,               &
                                     tsupercool_liq, grav_head
   use meds_plant_hydraulics, only : rhizosphere_cond
   use meds_hydr_lib, only : soil_hydr_cond_from_theta, soil_psi_from_theta, psi_from_water_content
   use meds_config,           only : meds_config_t, hydraulics_config_t,                          &
                                     SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED,               &
                                     INTEG_SPLIT, INTEG_ARK, INTEG_RK4
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
                                     surface_frozen_t, surface_tend_t, column_bflux_t,          &
                                     apply_hydraulics_config, mask_is_full
   use meds_fast_ark,         only : column_fast_step_ark, aero_bottom_to_top, column_prepass
   use meds_fast_rk45,        only : column_fast_step_rk45
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_soil_energy,      only : soil_energy_step_implicit
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_vegetation_biophysics, only : veg_energy_diagnostic, veg_energy_step_implicit,      &
                                     sensible_heat_coeff, lw_emission_slope, le_conductance_flux, &
                                     leaf_film_coeff, intercept_canopy_layer
   use meds_soil_water,       only : column_hydrology_flux
   use meds_ground_biophysics, only : snow_energy_step, snow_base_conductance,                  &
                                     snow_accumulate, snow_drain_meltwater, snow_cover_fraction, &
                                     ground_surface_fluxes
   use meds_plant_interface,  only : solve_plant_water_batch, N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_therm_lib,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_liquid,  &
                                     temp_to_uext
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok, budget_check_stop
   implicit none
   private

   public :: column_fast_step

contains

   !=======================================================================================!
   !  One fast (dt_fast) operator-split sweep for a single patch. Cohort arrays BOTTOM(1)->TOP. !
   !  Leaf gas exchange (real GPP + stomata + leaf Rd), stem/root maintenance respiration and    !
   !  heterotrophic Rh feed a physically-decomposed NEE = (Rd_leaf + stem + root) + Rh - GPP.    !
   !=======================================================================================!
   subroutine column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh, &
                               leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters,        &
                               le_flux, h_flux)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg          !< PFT traits for leaf gas exchange
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv         !< can_* fields refreshed from CAS state
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero         !< preallocated (alloc_aero_out)
      type(column_budget_t),   intent(inout) :: budg
      real(wp), optional,      intent(out)   :: gpp_coh(:)   !< [umol CO2/plant/s] per-cohort GROSS GPP (fast->slow)
      real(wp), optional,      intent(out)   :: leaf_resp_coh(:) !< [umol CO2/plant/s] leaf dark respiration (fast->slow)
      real(wp), optional,      intent(out)   :: stem_resp_coh(:) !< [umol CO2/plant/s] stem maintenance resp (fast->slow)
      real(wp), optional,      intent(out)   :: root_resp_coh(:) !< [umol CO2/plant/s] fine-root maint. resp (fast->slow)
      logical,     optional,   intent(out)   :: converged    !< Picard converged this sub-step (true for split)
      integer(ik), optional,   intent(out)   :: iters        !< outer-iteration count taken (1 for split)
      real(wp),    optional,   intent(out)   :: le_flux      !< [W/m2] CAS->atm latent-heat (ET) flux this sub-step
      real(wp),    optional,   intent(out)   :: h_flux       !< [W/m2] CAS->atm sensible-heat flux this sub-step

      type(chydro_forcing_t) :: hforc
      type(chydro_flux_t)    :: hflux
      type(energy_forcing_t) :: eforc
      type(energy_flux_t)    :: sflux
      !----- Plant hydraulics (MEDS_NUMERICS_SCOPING.md BB1 phase 2): solve_plant_water_batch takes    !
      !      bare per-cohort arrays (see meds_plant_hydraulics.f90), not the henv/hfx scalar bundle     !
      !      the old inline per-cohort loop used -- these are its inputs/outputs, all sized to the      !
      !      ACTIVE (coh%n, ccfg%soil%n_active) extent, independent of any oversized backing capacity.  !
      real(wp)               :: transp_pp(coh%n)                                !< [kg/plant/s] per-plant demand
      real(wp)               :: rhizo_cond_all(ccfg%soil%n_active, coh%n)       !< [kg/s/MPa] per-(layer,cohort)
      real(wp)               :: sapflow_b(coh%n), root_uptake_b(coh%n)          !< batch outputs (unused downstream
      real(wp)               :: root_uptake_layer_b(ccfg%soil%n_active, coh%n)  !< today except root_uptake_layer_b,
      real(wp)               :: psi_leaf_b(coh%n), psi_wood_b(coh%n), plc_b(coh%n)  !< kept as a complete SoA API)
      integer(ik)            :: nsub_b(coh%n)
      logical                :: converged_b(coh%n)
      !----- MEDS_ED2_RK45_DESIGN.md sec 4: leaf_water_mass/wood_water_mass are the PERSISTED state  !
      !      (psi is now a per-pass DIAGNOSTIC, not stored). psi_scratch is the batch call's inout    !
      !      buffer, freshly diagnosed from the persisted mass every Picard pass (the mass itself is  !
      !      NOT touched until the post-loop commit, so no separate state^n snapshot/reset is needed   !
      !      the way psi used to require -- bio%leaf_water_mass/wood_water_mass simply stay put).      !
      real(wp)               :: psi_scratch(N_HYDRO, coh%n)
      real(wp)               :: psi_soil_pre(ccfg%soil%n_active)               !< state^n soil psi (PRE-solve BC)
      real(wp)               :: total_uptake_b, scale                          !< aggregate root-uptake REQUEST + rescale
      real(wp)               :: w_plant0, w_plant1                             !< whole-water: plant internal storage term
      real(wp)               :: transp_c(coh%n)     !< [kg/m2 ground/s] per-cohort transpiration demand (automatic)
      real(wp)               :: h_coeff_f(coh%n), g_tr_f(coh%n), leaf_in(coh%n)   !< frozen coeffs + prev-iterate leaf temp
      real(wp)               :: t_emit(coh%n)      !< LW emission base (start leaf_temp; matches the RT tcan_bt, P3c)
      real(wp)               :: wood_emit(coh%n)   !< start-of-sub-step wood temp (prognostic-wood seed; Picard-correct)
      real(wp)    :: te
      type(soil_column_t)        :: soil_w_n        !< snapshot of the soil-water column at state^n (Picard reset)
      type(soil_energy_column_t) :: soil_e_n        !< snapshot of the soil thermal column at state^n (Picard reset)
      real(wp)    :: soil_psi_root, t_ground_dia, t_bot_dia, k_theta, sink_tot
      real(wp)    :: tcas, qcas, press, rho, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: le_ref, dtl, tl, transp_i, dh, drnet, transp_w
      real(wp)    :: dtw, lw_slope_w, h_coeff_w, twood, te_w    !< diagnostic WOOD balance (own store)
      real(wp)    :: wood_store0, wood_store1, dry_hcap_w, wmass_w, dbio_w   !< prognostic WOOD store
      real(wp)    :: leaf_store0, leaf_store1, cap_leaf, a_leaf, dbio_leaf   !< prognostic LEAF store (BE cap/dt term)
      type(leaf_energy_env_t)  :: wenv_e
      type(leaf_energy_flux_t) :: wflux
      real(wp), parameter :: C2B = 2.0_wp                      !< carbon->biomass (carbon fraction 0.5)
      real(wp), parameter :: WOOD_MOIST_FRAC = 1.0_wp          !< [kg water / kg dry] fresh-sapwood moisture (MVP)
      real(wp)    :: resid_T, tcas_in, qcas_in, tcas_new
      real(wp), parameter :: LAI_SLAVE_MIN = 1.0e-2_wp    !< [m2/m2] below this a cohort is slaved to tcas (Picard)
      integer(ik) :: iter, niter, niter_taken
      logical     :: picard, nconv
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet, coh_film_evap
      real(wp)    :: nee_biotic
      real(wp)    :: t_ground, t_bot, g_top, h_ground, le_ground, soil_evap, rain_temp
      real(wp)    :: gah, gaw, gac, wcap, ccap, src_enth, src_vap
      type(cas_source_t) :: cas_src
      type(cas_column_t) :: cas_col
      real(wp)    :: enth0, shv0, co20, enth1, shv1, co21
      !----- Snow store (opt-in, operator-split; §ccfg%snow_on). ------------------------------------!
      type(snow_env_t)  :: senv
      type(snow_flux_t) :: sfx
      type(snow_melt_t) :: smelt
      logical  :: snow_exists
      real(wp) :: snowfac_col, h_snow_s, le_snow_s, g_base_snow, subl_mass, melt_rate
      real(wp) :: snow_e0, snow_e1, swe0_s, swe1_s, snow_acc_enth, ground_rad_col, h_bare, le_soil
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, e_in, e_out, w_in, w_out
      !----- Enthalpy paired with the hydrology's UNFACED post-solve mass corrections (sec 3d'). -------!
      real(wp)    :: e_clip, e_floor, u_corr
      logical     :: halt_budgets     !< §5.1: hard-stop on a non-closing budget (full column + debug only)
      !----- Canopy-surface water (opt-in, operator-split; ccfg%canopy_water_on). Film-evap flux         !
      !      arrays (NOT scalars): each Picard pass overwrites them, and -- exactly like sapflow_b/        !
      !      root_uptake_b in the P0 hydraulics commit -- the LAST pass's values are what the post-loop    !
      !      surface-water commit below applies, so the CAS credit (inside the loop) and the store debit   !
      !      (after it) are GUARANTEED to be the same number. ------------------------------------------!
      real(wp) :: f_wet_c(coh%n)                       !< [-] combined wetted fraction, frozen once per dt_fast
      real(wp) :: film_evap(coh%n), film_evap_w(coh%n) !< [kg/m2/s] leaf/wood film evaporation (dew if negative)
      real(wp) :: le_slope_wet, le_ref_wet             !< leaf film-evap latent terms (scratch, reused per cohort)
      real(wp) :: le_slope_wet_w, le_ref_wet_w         !< wood film-evap latent terms (scratch, diagnostic branch)
      real(wp) :: rain_above, combined_w, pai_i, throughfall_i, drip_i, throughfall_total
      real(wp) :: intercepted_total, surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp) :: surf_overflow, leaf_cap_i, wood_cap_i
      !----- Named temporary for sum(array_expr): nvfortran's whole-program optimizer has a documented   !
      !      miscompilation trap around implicit array-expression temporaries at -O2/-fast (CLAUDE.md);   !
      !      binding to a named array first is the established, verified workaround in this codebase.     !
      real(wp) :: surf_water_tmp(coh%n)
      integer(ik) :: i, n, nsl, k

      n = coh%n ; nsl = ccfg%soil%n_active

      !----- §5.1 process mask: the closed-budget HARD STOPS assume the full column. A reduced column  !
      !      freezes a store while its fluxes still act on the neighbours, so it cannot conserve by     !
      !      construction -- the soft n_fail counters still tally (useful signal), only the halt is     !
      !      suppressed. All-on (the default) leaves this exactly as debug_error, so byte-identical.    !
      halt_budgets = ccfg%energy%debug_error .and. mask_is_full(ccfg%mask)

      !----- Prognostic wood (P2) and prognostic leaf (P3) both advance own stores on the SPLIT/PICARD    !
      !      path via the BE cap/dt term; the ARK/RK45 arrowhead couplings (P4) are deferred, so gate       !
      !      those (RK45 shares column_derivs' diagnostic-only leaf/wood path with ARK). ------------------!
      if (ccfg%wood_energy_model == WOODEN_PROGNOSTIC .and. cfg%time_integrator /= INTEG_SPLIT) &
         error stop 'column_fast_step: wood_energy_model="prognostic" under ARK/RK45 is deferred (P2 split-only)'
      if (ccfg%leaf_energy_model == LEAFEN_PROGNOSTIC .and. cfg%time_integrator /= INTEG_SPLIT) &
         error stop 'column_fast_step: leaf_energy_model="prognostic" under ARK/RK45 is deferred (P4 arrowhead)'
      !----- Canopy-surface water (interception/film-evap/dew, MEDS_ED2_RK45_DESIGN.md sec 3.4, P1+P2c)  !
      !      is now wired into the shared ARK/RK45 tableau too: column_state_t carries leaf_surf_water/     !
      !      wood_surf_water natively, build_column_frozen freezes interception + the wetted fraction +      !
      !      film-evap conductances once per dt_fast (mirroring every other Act-1 frozen quantity), and       !
      !      surface_derivs computes film evaporation per stage via veg_energy_diagnostic's existing wet-      !
      !      canopy extension. ARK's §8g TAU_COND condensation sink is a separate, pre-existing mechanism      !
      !      left untouched (not unified with the dew pathway here) -- avoid enabling both simultaneously.      !

      !----- TIME-INTEGRATOR dispatch (inserted BEFORE the first bio mutation, so the split path       !
      !      below is byte-for-byte unentered -- the golden anchor is preserved structurally). Both      !
      !      the coupled IMEX-ARK and the ED2-faithful adaptive Cash-Karp RK45 are opt-in                !
      !      ([fast].time_integrator="ark"|"rk45"); default is the split. --------------------------------!
      if (cfg%time_integrator == INTEG_ARK) then
         call column_fast_step_ark(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,   &
                                   gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
         if (present(le_flux)) le_flux = aenv%rho_air * aero%ustar * aero%temp2                    &
                                         * (bio%cas%can_shv - forc%shv_atm) * latent_heat_vap
         if (present(h_flux))  h_flux  = aenv%rho_air * aero%ustar * aero%temp2                    &
                                         * (bio%cas%can_temp - aenv%theta_atm) * cp_air
         return
      else if (cfg%time_integrator == INTEG_RK4) then
         !----- RK45 is FULLY EXPLICIT over the whole column (no implicit CAS box). At high LAI + cold  !
         !      the diagnostic leaf<->CAS coupling (stiff -- see this routine's own split-path comment   !
         !      below: "a single explicit split pass is marginally UNSTABLE ... a 2*dt oscillation")     !
         !      drives an oscillatory divergence the explicit stages cannot resolve even at the sub-step !
         !      floor: the state rails to the clamp bounds (CAS/soil pinned at 180 or 350 K, never       !
         !      physical) and -- worse -- the 5th/4th embedded solutions rail TOGETHER, so the error     !
         !      controller sees err~0 and commits the railed garbage, locking into a cold-dead attractor !
         !      (GPP->0) that never recovers (MEDS_ED2_RK45_DESIGN.md P6). HYBRID RESCUE: snapshot        !
         !      state^n, take the explicit RK45 step, and roll back + REDO this dt_fast on the stable     !
         !      implicit-CAS split path below (the validated default -- the 30-yr split run is healthy    !
         !      through exactly this dense-canopy-winter regime) when EITHER trigger fires:               !
         !        * stiff_bail -- the explicit march burned its whole work budget (RK45_WORK_CAP) without !
         !          resolving the step. Catches the runaway EARLY and cheaply; without it a persistently  !
         !          stiff spell (a cold December under a closed canopy) makes every sub-step thrash to    !
         !          hundreds of explicit sub-steps BEFORE being discarded, which is what turned a 6-min   !
         !          30-yr run into one that could not get past a single month.                            !
         !        * rk45_state_railed -- the step did finish but committed a clamp-pinned CAS/soil state. !
         !      RK45 keeps its fast explicit path everywhere it is stable; only genuinely-stiff steps     !
         !      pay the split cost. ----------------------------------------------------------------------!
         block
            type(patch_biophys_t) :: bio_save
            type(column_budget_t) :: budg_save
            logical               :: rk45_stiff
            bio_save = bio ; budg_save = budg
            call column_fast_step_rk45(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, &
                                       gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, &
                                       stiff_bail=rk45_stiff)
            if (.not. rk45_stiff .and. .not. rk45_state_railed(bio, ccfg%soil%n_active)) then
               if (present(le_flux)) le_flux = aenv%rho_air * aero%ustar * aero%temp2               &
                                               * (bio%cas%can_shv - forc%shv_atm) * latent_heat_vap
               if (present(h_flux))  h_flux  = aenv%rho_air * aero%ustar * aero%temp2               &
                                               * (bio%cas%can_temp - aenv%theta_atm) * cp_air
               return
            end if
            bio = bio_save ; budg = budg_save         ! discard the railed/bailed RK45 step (rollback)
            budg%integ_nrej   = budg%integ_nrej   + 1_ik  ! count the rescue as a rejected integrator step
            budg%rk45_rescue  = budg%rk45_rescue  + 1_ik  ! ...and as a P6 RK45->split rescue (diagnostic)
         end block
         !----- fall through to the stable implicit-CAS split path with state^n restored ----------------!
      end if

      picard = (cfg%integration_scheme == SCHEME_PICARD_COUPLED)
      niter  = 1_ik ; if (picard) niter = max(1_ik, ccfg%picard_max_iter)
      !----- Prognostic leaf (P3): a BE storage term (cap_leaf/dt) is added to the diagnostic leaf     !
      !      linearization below (numerator + denominator), so diagnostic (cap_leaf=0) stays bit-       !
      !      identical and the leaf<->CAS coupling is exactly the Picard iterate. Leaf thermal inertia  !
      !      is tiny at dt_fast ~ 900 s, so prognostic ~ diagnostic; the option exists to quantify it.   !
      !      The leaf<->CAS coupling MUST be solved implicitly: the leaf is stiff and its sensible +     !
      !      latent flux feeds the CAS, which feeds back to the leaf. A single explicit split pass       !
      !      (SCHEME_SPLIT_SEQUENTIAL, niter=1) is marginally UNSTABLE with the storage term (a 2*dt      !
      !      oscillation, ~1.7 K midday spikes on the census stand); the Picard iterate damps it to       !
      !      ~0.2 K. So require Picard for a prognostic leaf. (Wood has no transpiration feedback and is   !
      !      stable on the pure-split path, so P2 is not gated this way.)                                 !
      if (ccfg%leaf_energy_model == LEAFEN_PROGNOSTIC .and. .not. picard) &
         error stop 'column_fast_step: leaf_energy_model="prognostic" needs integration_scheme="picard" &
                    &(explicit leaf<->CAS split is unstable with leaf storage)'
      !----- The soil-water coupling selector: LAGGED and COUPLED both currently re-solve the soil    !
      !      water from state^n each Picard pass (required for conservation while the leaf demand      !
      !      iterates). A true frozen/lagged optimization (thermal-only, cheaper) is deferred (P3f).   !

      !----- Snapshot start-of-step SOIL stores (for the whole-column budgets). --------------!
      e_soil0 = 0.0_wp ; w_soil0 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil0 = e_soil0 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil0 = w_soil0 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      !----- Plant internal water (MEDS_ED2_RK45_DESIGN.md sec 3): now a REAL, separately-evolving   !
      !      store (mass, not psi), so it must join the whole_water ledger below or a nonzero storage !
      !      change would read as a leak. ------------------------------------------------------------!
      w_plant0 = sum(coh%nplant(1:n) * (bio%leaf_water_mass(1:n) + bio%wood_water_mass(1:n)))

      !----- The SHARED pre-pass (once per sub-step; ED2 freezes gs/hydraulics per DTLSM): refreshes  !
      !      the aerodynamics + CAS-derived scalars, then computes LEAF gas exchange (GPP/gs/Rd), the  !
      !      FROZEN per-cohort leaf-energy coefficients h_coeff_f/g_tr_f, stem/root maintenance         !
      !      respiration, NEE, and the CAS capacities/conductances -- ONE authority shared with          !
      !      build_column_frozen (meds_fast_ark%column_prepass), so split/ARK GPP stay bit-for-bit.      !
      !      These do NOT change across the Picard passes (they use the lagged, start-of-sub-step        !
      !      leaf_temp). column_prepass takes bio as intent(in), so the CAS-temperature persistence       !
      !      write is done here, right after the call.                                                    !
      call column_prepass(cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,                          &
                          tcas, qcas, press, rho, t_ground, h_coeff_f, g_tr_f,                          &
                          wcap, ccap, gah, gaw, gac, nee_biotic,                                        &
                          gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)
      bio%cas%can_temp = tcas
      t_bot = bio%soil_e%soil_temp(nsl) ; rain_temp = tcas
      t_emit(1:n)    = bio%leaf_temp(1:n)   ! start-of-sub-step leaf temp = the RT LW emission base (P3c)
      wood_emit(1:n) = bio%wood_temp(1:n)   ! start-of-sub-step wood temp = prognostic-wood seed (Picard-correct)

      !----- Frozen CAS-column params for the shared cas_column_step_implicit box kernel. ---------!
      cas_col%air_mass_capacity        = wcap
      cas_col%air_molar_capacity       = ccap
      cas_col%atm_conductance_enthalpy = gah
      cas_col%atm_conductance_vapor    = gaw
      cas_col%atm_conductance_co2      = gac
      cas_col%atm_enthalpy             = forc%enthalpy_atm
      cas_col%atm_specific_humidity    = forc%shv_atm
      cas_col%atm_co2                  = forc%co2_atm

      !----- 2b. SNOW store (opt-in, OPERATOR-SPLIT out of the Picard iterate). Accumulate snowfall +  !
      !      rain-on-snow, advance the snow-surface energy balance at the LAGGED CAS, and drain melt-    !
      !      water INTO the soil top as a PAIRED (mass, enthalpy) transfer -- done BEFORE the state^n     !
      !      snapshot so the melt enthalpy is inside soil_e_n. When snow is present it OWNS the surface:  !
      !      the Picard surface balance below uses the snow sensible/sublimation + throttled base         !
      !      conduction, and the ground latent (soil evaporation) is suppressed. MVP: binary-when-present !
      !      (the snowfac continuity ramp of the design is a follow-up). No-op when snow_on is false.     !
      snow_exists = .false. ; snowfac_col = 0.0_wp
      h_snow_s = 0.0_wp ; le_snow_s = 0.0_wp ; g_base_snow = 0.0_wp ; subl_mass = 0.0_wp ; melt_rate = 0.0_wp
      snow_acc_enth = 0.0_wp
      ground_rad_col = forc%abs_sw_ground + forc%abs_lw_ground     ! bare-ground radiation boundary input (snowfac=0)
      snow_e0 = bio%snow%snow_energy(1) ; swe0_s = bio%snow%swe(1)
      if (ccfg%snow_on) then
         call snow_accumulate(bio%snow, forc%snowf, forc%precip, forc%tair, dt_fast, ccfg%snow)
         snow_acc_enth = bio%snow%snow_energy(1) - snow_e0         ! precip enthalpy that entered the pack (boundary in)
         snow_exists = bio%snow%nlayer >= 1_ik                     ! pack present: accumulate took snow+rain (precip routing)
         if (snow_exists .and. bio%snow%swe(1) > ccfg%snow%tiny_snow_mass) then
            !----- SUB-COLUMN: snowfac fraction is snow, (1-snowfac) is bare soil (design §4f/§4g/§6). The  !
            !      snow store advances with its boundary exchange SCALED by snowfac (so a thin/patchy pack    !
            !      barely exchanges -> continuous + stable, no threshold cliff); its fluxes are already        !
            !      snowfac-weighted. The bare-soil (1-snowfac) fraction is blended in below.                   !
            snowfac_col   = snow_cover_fraction(bio%snow%swe(1), bio%snow%snow_depth(1), ccfg%snow)
            senv%abs_sw = forc%abs_sw_ground ; senv%abs_lw = forc%abs_lw_ground
            senv%can_temp = tcas ; senv%can_shv = qcas ; senv%ggnet = aero%ggnet
            senv%rho_air = rho ; senv%press = press
            senv%t_soil_top = bio%soil_e%soil_temp(1)
            senv%dz_soil_top = max(-ccfg%soil%z_node(1), tiny_num)      ! |z_node(1)| = top-node depth
            call snow_energy_step(bio%snow, senv, ccfg%snow, dt_fast, snowfac_col, sfx)
            call snow_drain_meltwater(bio%snow, ccfg%snow, smelt)
            h_snow_s = sfx%h_snow ; le_snow_s = sfx%le_snow ; g_base_snow = sfx%g_base   ! snowfac-weighted
            snowfac_col = sfx%snowfac                                     ! the clamped/actual fraction the kernel used
            !----- Ground radiation boundary in = snow's snowfac-weighted net + bare's (1-snowfac) share. --!
            ground_rad_col = sfx%rnet + (1.0_wp - snowfac_col) * (forc%abs_sw_ground + forc%abs_lw_ground)
            subl_mass = sfx%w_flux * dt_fast
            melt_rate = (smelt%melt_mass + smelt%dump_mass) / dt_fast     ! meltwater -> infiltration
            !----- paired enthalpy: snow store -> soil top (extensive J/m2 -> volumetric J/m3). ------!
            bio%soil_e%soil_energy(1) = bio%soil_e%soil_energy(1)                                    &
                                      + (smelt%melt_enth + smelt%dump_enth) / ccfg%soil%dz(1)
         end if
      end if
      snow_e1 = bio%snow%snow_energy(1) ; swe1_s = bio%snow%swe(1)
      !----- Meltwater already handed its enthalpy to the soil top via the PAIRED transfer above, so    !
      !      the hydrology must infiltrate its MASS with ZERO enthalpy (rain_temp = tsupercool_liq =>     !
      !      internal_energy_liquid = 0) -- otherwise the melt enthalpy is double-counted at soil layer 1. !
      if (snow_exists) rain_temp = tsupercool_liq

      !----- 2c. CANOPY INTERCEPTION (opt-in, MEDS_ED2_RK45_DESIGN.md sec 3.4, P1): a top-to-bottom      !
      !      sweep over the height-DESCENDING cohort order (coh's native gather order), ONE combined       !
      !      leaf+wood bucket per cohort (intercept_canopy_layer's own Beer-law screen + capacity already   !
      !      combine lai+sai this way), split back into the two persisted stores by area-index share         !
      !      afterward. Frozen ONCE per dt_fast (not re-run per Picard pass), mirroring how gs/hydraulics     !
      !      are frozen. e_canopy=0 here DELIBERATELY: this call does ONLY capture/capacity (no              !
      !      evaporation) -- the ACTUAL film evaporation is computed from the Picard-CONVERGED leaf/wood      !
      !      energy balance below and committed to this SAME store afterward (mirroring exactly how P0        !
      !      commits leaf_water_mass/wood_water_mass once, after the loop, from the converged flux). Using    !
      !      a DIFFERENT one-shot evaporation estimate here (to size q_intr's headroom) would credit the       !
      !      CAS with one number while debiting the store a DIFFERENT number -- exactly the "two               !
      !      definitions at one interface" pitfall sec 9 warns against, just for a new interface. SKIPPED      !
      !      when snow owns the surface (canopy interception of LIQUID rain is a separate, orthogonal          !
      !      process this design doc does not specify interacting with snow, so the two stay mutually          !
      !      exclusive, exactly like the precip_ground branch below). --------------------------------------!
      f_wet_c(1:n) = 0.0_wp
      throughfall_total = forc%precip + forc%snowf   ! default: no canopy store -> everything reaches the ground
      surf_water0 = 0.0_wp ; intercepted_total = 0.0_wp ; surf_enth0 = 0.0_wp
      !----- Whole ledger block (surf_water0/intercepted_total/surf_enth0) is gated behind the SAME flag !
      !      that gates the interception sweep itself: every term here is exactly 0.0 whenever the        !
      !      feature is off (bio%leaf_surf_water/wood_surf_water never leave their 0.0 birth seed), so     !
      !      skipping the arithmetic entirely (not just computing-then-discarding it) is behavior-        !
      !      preserving and avoids needless per-step FLOPs on the byte-identical default path. -----------!
      if (ccfg%canopy_water_on) then
         surf_water_tmp(1:n) = bio%leaf_surf_water(1:n) + bio%wood_surf_water(1:n)
         surf_water0 = sum(surf_water_tmp(1:n))
         if (.not. snow_exists) then
            rain_above = forc%precip + forc%snowf
            do i = 1_ik, n
               pai_i      = coh%lai(i) + coh%wai(i)
               combined_w = bio%leaf_surf_water(i) + bio%wood_surf_water(i)
               call intercept_canopy_layer(combined_w, rain_above, coh%lai(i), coh%wai(i), 0.0_wp, dt_fast, &
                                           ccfg%hydro%dewmx, ccfg%hydro%intercept_k, ccfg%hydro%intercept_alpha, &
                                           throughfall_i, drip_i, f_wet_c(i))
               !----- Split the updated combined bucket back by area-index share (intercept_canopy_layer  !
               !      combines lai+sai for the interception SCREEN, but leaf/wood are separate energy-     !
               !      balance stores -- a simpler treatment than solving two independent, coupled          !
               !      buckets). ------------------------------------------------------------------------!
               if (pai_i > tiny_num) then
                  bio%leaf_surf_water(i) = combined_w * coh%lai(i) / pai_i
                  bio%wood_surf_water(i) = combined_w * coh%wai(i) / pai_i
               else
                  bio%leaf_surf_water(i) = 0.0_wp ; bio%wood_surf_water(i) = 0.0_wp
               end if
               rain_above = throughfall_i   ! cascades to the next (shorter) cohort
            end do
            throughfall_total = rain_above   ! whatever survives the shortest (last) cohort
         end if
         !----- Enthalpy of the intercepted (not-yet-evaporated) rain, at rain_temp -- the e_in boundary  !
         !      term added to whole_energy below (sec 7b); surf_enth0/1 (the STORE side) are computed      !
         !      from this SAME fixed reference once the post-Picard-loop film-evap commit gives the TRUE   !
         !      end-of-step surface water (see the commit alongside the internal-water-mass one, below).  !
         intercepted_total = (forc%precip + forc%snowf) - throughfall_total
         surf_enth0 = surf_water0 * internal_energy_liquid(rain_temp)
      end if

      !----- Snapshot state^n once: the Picard passes re-solve the SAME backward-Euler steps FROM  !
      !      this base each iteration (only the source, re-evaluated at the iterate, changes). The   !
      !      soil-water column is prognostic and column_hydrology_flux ADVANCES it, so it must be     !
      !      reset to state^n before each re-solve or it double-steps. Plant hydraulics needs NO such  !
      !      snapshot/reset any more: bio%leaf_water_mass/wood_water_mass (the persisted state) are    !
      !      only WRITTEN once, after the Picard loop converges (below), so every pass's psi_scratch   !
      !      is diagnosed fresh from the same untouched state^n mass -- see MEDS_ED2_RK45_DESIGN.md    !
      !      sec 4/5. -----------------------------------------------------------------------------!
      enth0 = bio%cas%can_enthalpy ; shv0 = bio%cas%can_shv ; co20 = bio%cas%can_co2
      soil_w_n = bio%soil_w ; soil_e_n = bio%soil_e

      !======================================================================================!
      !  3. Outer PICARD fixed point over { leaf energy -> soil water/hydraulics -> CAS twins }.  !
      !     niter = 1 under SCHEME_SPLIT_SEQUENTIAL reproduces the operator-split sweep EXACTLY   !
      !     (one pass, no convergence test, soil water solved that pass). Under PICARD the block   !
      !     iterates at the current tcas/qcas until the store temperatures converge; the soil       !
      !     water/hydraulics are frozen after pass 1 EVERY pass regardless of                       !
      !     ccfg%soil_water_coupling -- that selector is RESERVED for the P3f re-solve-inside-        !
      !     Picard optimization and has NO behavioral effect yet (both SOILH2O_LAGGED and             !
      !     SOILH2O_COUPLED take this same frozen path today).                                        !
      !======================================================================================!
      soil_evap = 0.0_wp ; nconv = .false. ; resid_T = 0.0_wp
      e_clip = 0.0_wp ; e_floor = 0.0_wp        ! read by the sec 7b ledger; set per Picard pass below
      niter_taken = 0_ik
      do iter = 1_ik, niter
         niter_taken = iter
         tcas_in = tcas ; qcas_in = qcas ; leaf_in(1:n) = bio%leaf_temp(1:n)

         !----- 3a. Leaf energy balance (diagnostic) at the CURRENT tcas/qcas, frozen coeffs. --!
         qsat_c = sat_specific_humidity(tcas, press)
         dqdt   = sat_specific_humidity_temp_deriv(tcas, press)
         coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
         coh_film_evap = 0.0_wp
         wood_store0 = 0.0_wp ; wood_store1 = 0.0_wp ; leaf_store0 = 0.0_wp ; leaf_store1 = 0.0_wp
         do i = 1_ik, n
            if (picard .and. coh%lai(i) < LAI_SLAVE_MIN) then    ! near-zero LAI: slave to CAS, no exchange
               bio%leaf_temp(i) = tcas ; bio%wood_temp(i) = tcas ; transp_c(i) = 0.0_wp ; cycle
            end if
            !----- LW emission linearized around T_emit: tcas for SPLIT (reduces to the current    !
            !      form, so split stays bit-identical) or the start leaf_temp for PICARD (which the   !
            !      two-stream also emits at via tcan_bt, P3c) -> leaf emission consistent at leaf_temp. !
            te = tcas ; if (picard) te = t_emit(i)
            lw_slope = lw_emission_slope(ccfg%veg_thermal%leaf_emiss, te, coh%lai(i))
            le_slope = le_conductance_flux(rho, g_tr_f(i), dqdt)
            le_ref   = le_conductance_flux(rho, g_tr_f(i), qsat_c - qcas)
            !----- Prognostic leaf (P3): backward-Euler storage term a_leaf = cap_leaf/dt. The leaf     !
            !      relaxes from its start-of-sub-step temperature t_emit(i); a_leaf=0 (diagnostic) makes  !
            !      dtl the steady-state solve EXACTLY. cap_leaf is the leaf dry heat capacity floored by   !
            !      veg_hcap_min so it is > 0 (a_leaf finite) even for a near-massless cohort.              !
            cap_leaf = 0.0_wp ; a_leaf = 0.0_wp
            if (ccfg%leaf_energy_model == LEAFEN_PROGNOSTIC) then
               dbio_leaf = coh%bleaf(i) * coh%nplant(i) * C2B                        ! [kg dry leaf/m2]
               cap_leaf  = max(dbio_leaf * ccfg%veg_thermal%c_leaf, ccfg%veg_thermal%veg_hcap_min)
               a_leaf    = cap_leaf / dt_fast                                        ! [W/m2/K] storage conductance
            end if
            !----- Wetted-canopy film-evap latent terms (sec 3.4, P1): boundary-layer-only conductance    !
            !      (no stomatal resistance), same effarea_evap sidedness convention leaf_film_coeff uses.  !
            !      Harmless to compute even with canopy_water_on off: f_wet_c(i) stays 0 there, which       !
            !      makes veg_energy_diagnostic's wet pathway vanish identically (proven in its own          !
            !      docstring), so these two lines are pure no-ops on the byte-identical default path. ------!
            le_slope_wet = le_conductance_flux(rho, leaf_film_coeff(ccfg%veg_thermal%effarea_evap,          &
                 coh%lai(i), aero%leaf_gbw(i)), dqdt)
            le_ref_wet   = le_conductance_flux(rho, leaf_film_coeff(ccfg%veg_thermal%effarea_evap,          &
                 coh%lai(i), aero%leaf_gbw(i)), qsat_c - qcas)
            call veg_energy_diagnostic(forc%abs_sw(i), forc%abs_lw(i), h_coeff_f(i), le_slope,       &
                                       lw_slope, le_ref, tcas, te, a_leaf, t_emit(i),                &
                                       dtl, tl, transp_i, dh, drnet,                                  &
                                       f_wet_c(i), le_slope_wet, le_ref_wet, film_evap(i))
            !----- Clamp EVAPORATION (film_evap>0) to the water actually sitting there (bio%leaf_surf_     !
            !      water(i) is fixed across Picard passes -- interception already ran; the debit commits    !
            !      post-loop), so the CAS credit here and the post-loop store debit can never disagree        !
            !      (sec 9's "closure-killers at the store's bounds": clamp, don't silently over-apply).        !
            !      DEW (film_evap<0) is deliberately left UNCLAMPED here: clamping it post-hoc would leave     !
            !      dh/drnet/tl (already solved above, self-consistently, AGAINST the unclamped value) out       !
            !      of sync with a different, clamped film_evap credited to the CAS -- a WORSE inconsistency     !
            !      than the one being guarded against, and empirically a large one (confirmed: at f_wet=1       !
            !      exactly -- a just-saturated film -- clamping dew to 0 while leaving the wet-pathway           !
            !      conductance fully active elsewhere in the SAME solve stranded the leaf's entire latent        !
            !      capacity into sensible heat, a multi-hundred-W/m2 energy error, not a small one). KNOWN       !
            !      DEFERRED LIMITATION: bio%leaf_surf_water can therefore transiently exceed its capacity        !
            !      after a large dew event; intercept_canopy_layer's own w_max ceiling silently re-clamps it     !
            !      on a later step (a small, bounded water debit this pass does not yet bookkeep -- verified     !
            !      empirically to stay within the whole_water tolerance for the scenarios this pass tests).      !
            !      Harmless when the feature is off (film_evap(i) is already exactly 0 there). ------------------!
            film_evap(i) = min(film_evap(i), bio%leaf_surf_water(i) / dt_fast)
            bio%leaf_temp(i) = tl
            leaf_store0 = leaf_store0 + cap_leaf * t_emit(i)     ! [J/m2] leaf internal energy (0 K ref; differenced)
            leaf_store1 = leaf_store1 + cap_leaf * tl            ! diagnostic: cap_leaf=0 -> telescopes to 0
            transp_c(i) = transp_i                                                       ! per-cohort demand (hydraulics)
            coh_h      = coh_h      + dh
            !----- CAS latent: transpiration (vapour enthalpy) + film evaporation (ditto; dew if <0). ----!
            coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl) + film_evap(i) * enthalpy_vapor(tl)
            coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)   ! liquid enthalpy soil sheds
            coh_transp = coh_transp + transp_i
            !----- CAS MASS: film evaporation is a real vapour-mass source (dew if <0), separate from    !
            !      coh_transp (which only carries the stomatal/transpiration mass) -- feeds src_vap below. !
            coh_film_evap = coh_film_evap + film_evap(i)
            coh_rnet   = coh_rnet   + drnet
            !----- 3a'. Diagnostic WOOD energy balance (own store; own boundary layer + net LW, NO      !
            !      transpiration). Wood sensible + net-LW join coh_h / coh_rnet -> CAS + energy budget.   !
            !      A diagnostic wood has no storage, so absorbed = emitted + sensible-to-CAS -> the        !
            !      coh_rnet and coh_h wood terms are EQUAL (h_coeff_w*dtw) and telescope in the ledger.    !
            if (ccfg%wood_energy_model == WOODEN_DIAGNOSTIC) then
               h_coeff_w  = sensible_heat_coeff(pi * coh%wai(i), aero%wood_gbh(i), rho, cp_air)
               te_w       = tcas
               lw_slope_w = lw_emission_slope(ccfg%veg_thermal%leaf_emiss, te_w, coh%wai(i))
               !----- Wood film-evap latent terms (sec 3.4, P1): lifts the OLD "dry bark, no film evap"    !
               !      MVP simplification that used to apply here -- aero%wood_gbw was already computed by   !
               !      the aerodynamics kernel (mirroring leaf's), just never consumed. Harmless no-op when   !
               !      canopy_water_on is off (f_wet_c(i)=0 there, same proof as the leaf call above). --------!
               le_slope_wet_w = le_conductance_flux(rho, leaf_film_coeff(ccfg%veg_thermal%effarea_evap,      &
                    coh%wai(i), aero%wood_gbw(i)), dqdt)
               le_ref_wet_w   = le_conductance_flux(rho, leaf_film_coeff(ccfg%veg_thermal%effarea_evap,      &
                    coh%wai(i), aero%wood_gbw(i)), qsat_c - qcas)
               !----- Diagnostic WOOD = the le_slope = le_ref = 0 case of the same kernel (no transp). --!
               call veg_energy_diagnostic(forc%abs_sw_wood(i), forc%abs_lw_wood(i), h_coeff_w,       &
                                          0.0_wp, lw_slope_w, 0.0_wp, tcas, te_w, 0.0_wp, tcas,      &
                                          dtw, twood, transp_w, dh, drnet,                            &
                                          f_wet_c(i), le_slope_wet_w, le_ref_wet_w, film_evap_w(i))
               !----- Same evaporation-only clamp as the leaf call above (dew left unclamped -- see the    !
               !      leaf site's comment for why; bio%wood_surf_water(i) is fixed across Picard passes). --!
               film_evap_w(i) = min(film_evap_w(i), bio%wood_surf_water(i) / dt_fast)
               bio%wood_temp(i) = twood
               coh_h    = coh_h    + dh
               coh_qw   = coh_qw   + film_evap_w(i) * enthalpy_vapor(twood)   ! wood film-evap -> CAS
               coh_film_evap = coh_film_evap + film_evap_w(i)   ! CAS mass source (dew if <0), see leaf site above
               coh_rnet = coh_rnet + drnet
            else   ! WOODEN_PROGNOSTIC: advance the wood internal-energy store (operator-split, non-stiff). !
               !----- Canopy water is NOT wired into the prognostic-wood path (a separate, already opt-in  !
               !      feature via veg_surface_fluxes/veg_energy_step_implicit) in this pass -- out of       !
               !      scope, mirrors the diagnostic-only wiring for prognostic leaf too (§3a is diagnostic- !
               !      only). film_evap_w stays 0 here, so the surface-water commit below is a no-op for      !
               !      this cohort regardless of canopy_water_on. --------------------------------------------!
               film_evap_w(i) = 0.0_wp
               dbio_w     = coh%bsap(i) * coh%nplant(i) * C2B              ! [kg dry biomass/m2]
               dry_hcap_w = max(dbio_w * ccfg%veg_thermal%c_sapw, ccfg%veg_thermal%veg_hcap_min)  ! absolute floor > 0 (cap/=0)
               wmass_w    = dbio_w * WOOD_MOIST_FRAC                       ! [kg water/m2] fresh-sapwood water
               wenv_e%abs_sw = forc%abs_sw_wood(i) ; wenv_e%abs_lw = forc%abs_lw_wood(i)
               wenv_e%can_temp = tcas ; wenv_e%can_shv = qcas
               wenv_e%gbh = aero%wood_gbh(i) ; wenv_e%gbw = 0.0_wp    ! MVP wood = dry bark: NO film evap/dew
               wenv_e%gsw = 0.0_wp ; wenv_e%fs_open = 0.0_wp ; wenv_e%area_index = coh%wai(i)
               wenv_e%leaf_water = 0.0_wp ; wenv_e%wmass = wmass_w ; wenv_e%dry_hcap = dry_hcap_w
               wenv_e%rho_air = rho ; wenv_e%press = press
               twood = temp_to_uext(dry_hcap_w, wmass_w, wood_emit(i), 1.0_wp)  ! seed store from start-of-sub-step temp
               wood_store0 = wood_store0 + twood
               call veg_energy_step_implicit(twood, wenv_e, ccfg%veg_thermal, dt_fast, .false., wflux)
               !----- §5.1 process mask: hold the wood store at its entry enthalpy (the fluxes it drives  !
               !      into the CAS below are kept, so the canopy still sees a constant-temperature stem). !
               if (.not. ccfg%mask%veg_energy) twood = temp_to_uext(dry_hcap_w, wmass_w, wood_emit(i), 1.0_wp)
               wood_store1 = wood_store1 + twood                          ! store energy AFTER the BE step
               bio%wood_temp(i) = wflux%temp
               if (.not. ccfg%mask%veg_energy) bio%wood_temp(i) = wood_emit(i)
               coh_h    = coh_h    + wflux%h_flux                         ! wood sensible -> CAS
               coh_qw   = coh_qw   + wflux%qw_flux                        ! wood film-evap -> CAS (0 in MVP)
               coh_rnet = coh_rnet + (forc%abs_sw_wood(i) + forc%abs_lw_wood(i))   ! net wood radiation into the column
            end if
         end do

         !----- 3b. Plant HYDRAULICS + soil WATER column, RE-SOLVED FROM state^n each pass so the     !
         !          realized uptake stays consistent with the iterated leaf demand -- freezing it      !
         !          while the demand iterates would leak water/enthalpy. REORDERED (MEDS_ED2_RK45_     !
         !          DESIGN.md sec 3/5): hydraulics now runs BEFORE the soil solve, using psi diagnosed   !
         !          from state^n theta (not the Richards-solved end-of-step psi) and the FULL            !
         !          transpiration demand -- no supply pre-throttle, since the plant's own leaf/wood       !
         !          water storage now buffers any step-to-step soil-supply/demand mismatch (the           !
         !          instantaneous src_frac clamp this design replaces is gone). The resulting aggregate    !
         !          root-uptake REQUEST becomes the soil's root-sink forcing; the soil's own fwilt-        !
         !          limited sink (meds_soil_water) may honour less than requested (e.g. near wilting       !
         !          point), so a rescale is applied to the credited uptake below. pass 1 is already at     !
         !          state^n; later passes reset the prognostic soil water to the snapshot before           !
         !          re-advancing -- plant hydraulics needs NO such reset (bio%leaf_water_mass/              !
         !          wood_water_mass are read, never written, until the post-loop commit). --------------!
         if (iter > 1_ik) then
            bio%soil_w = soil_w_n
         end if
         !----- Soil psi at state^n (PRE-solve), from the theta this pass starts from -- the hydraulics   !
         !      boundary condition, replacing the old post-solve hflux%psi_soil now that hydraulics runs   !
         !      first. Elemental broadcast over the per-layer theta/van-Genuchten arrays. -----------------!
         psi_soil_pre(1:nsl) = soil_psi_from_theta(ccfg%soil%retention, bio%soil_w%theta(1:nsl),          &
              ccfg%soil%theta_sat(1:nsl), ccfg%soil%theta_res(1:nsl), ccfg%soil%vg_alpha(1:nsl),          &
              ccfg%soil%vg_n(1:nsl))
         soil_psi_root = root_weighted_psi(psi_soil_pre, ccfg%soil%root_frac, nsl)
         !----- Per-cohort plant hydraulics via the BARE-ARRAY batch seam (MEDS_NUMERICS_SCOPING.md BB1  !
         !      phase 2). Precompute the per-(layer,cohort) rhizosphere conductance first (identical      !
         !      k_theta/rhizosphere_cond calls to the old inline loop, just gathered into an array;        !
         !      order-independent since no accumulator is shared across (k,i) pairs; now evaluated at       !
         !      state^n theta, ahead of the soil solve), then hand the WHOLE patch to solve_plant_water_    !
         !      batch in one call. --------------------------------------------------------------------!
         if (ccfg%multilayer_roots) then
            !----- Couple to the per-layer soil column: per-layer psi_soil [MPa] at state^n (same for      !
            !       every cohort, passed once below), and a per-(layer,cohort) rhizosphere conductance       !
            !       (Katul; ksat converted m/s -> kg/m/s/MPa). --------------------------------------------!
            do i = 1_ik, n
               do k = 1_ik, nsl
                  !----- UNSATURATED K(theta) [m/s] (not ksat): dry layers are conductance-down-weighted. !
                  k_theta = soil_hydr_cond_from_theta(ccfg%soil%retention, bio%soil_w%theta(k),                     &
                       ccfg%soil%theta_sat(k), ccfg%soil%theta_res(k), ccfg%soil%vg_alpha(k),            &
                       ccfg%soil%vg_n(k), ccfg%soil%ksat(k))
                  rhizo_cond_all(k, i) = rhizosphere_cond(rho_h2o*k_theta/grav_head,                 &
                       coh%broot(i), ccfg%specific_root_area, ccfg%soil%root_frac(k), ccfg%soil%dz(k),   &
                       coh%nplant(i))
               end do
            end do
         end if
         !----- Diagnose the entry psi from the PERSISTED water mass (state^n, untouched this whole      !
         !      Picard loop): psi_from_water_content is the exact inverse of water_content, so this        !
         !      round-trips losslessly into solve_plant_water's own internal psi representation. -----------!
         psi_scratch(NODE_LEAF, 1:n) = psi_from_water_content(bio%leaf_water_mass(1:n),                   &
              ccfg%hydro_p%leaf_pi0, ccfg%hydro_p%leaf_elastic_mod, ccfg%hydro_p%leaf_apoplast_frac,      &
              ccfg%hydro_p%leaf_water_sat, coh%bleaf(1:n))
         psi_scratch(NODE_WOOD, 1:n) = psi_from_water_content(bio%wood_water_mass(1:n),                   &
              ccfg%hydro_p%wood_pi0, ccfg%hydro_p%wood_elastic_mod, ccfg%hydro_p%wood_apoplast_frac,      &
              ccfg%hydro_p%wood_water_sat, coh%bsap(1:n) + coh%broot(1:n))
         transp_pp(1:n) = transp_c(1:n) / max(coh%nplant(1:n), tiny_num)   ! [kg/plant/s] FULL demand, unthrottled
         call solve_plant_water_batch(n, nsl, ccfg%multilayer_roots, transp_pp(1:n), coh%bleaf(1:n),      &
                                      coh%bsap(1:n), coh%broot(1:n), coh%sap_area(1:n), coh%height(1:n),   &
                                      coh%leaf_area(1:n), soil_psi_root, ccfg%rhizo_cond,                   &
                                      psi_soil_pre(1:nsl), ccfg%soil%z_node(1:nsl), rhizo_cond_all(1:nsl, 1:n), &
                                      ccfg%hydro_p, ccfg%hydro_o, dt_fast, psi_scratch(:, 1:n),                  &
                                      sapflow_b(1:n), root_uptake_b(1:n), root_uptake_layer_b(1:nsl, 1:n),    &
                                      psi_leaf_b(1:n), psi_wood_b(1:n), plc_b(1:n), nsub_b(1:n), converged_b(1:n))
         !----- section 5.3 work counters: hydraulics sub-steps summed over cohorts + non-convergences. !
         budg%hydro_nsub    = sum(nsub_b(1:n))
         budg%hydro_nonconv = count(.not. converged_b(1:n))
         !----- Hydraulic redistribution (root EFFLUX, a negative root_uptake_b) is intentionally NOT    !
         !      enabled anywhere in this model (CLAUDE.md; solve_plant_water already floors it OUT of    !
         !      the per-layer multilayer_roots distribution) -- but the AGGREGATE flux%root_uptake =      !
         !      (dw_leaf+dw_wood)/dt + transp returned by solve_plant_water itself carries no such floor,  !
         !      and CAN go negative in a real transient (e.g. right after the PSI_INIT lazy-init seed,      !
         !      if the seeded plant psi is wetter than the ambient soil, the solver's own psi trajectory     !
         !      implies net efflux while it re-equilibrates). Floor the WHOLE aggregate here, at the same   !
         !      place the per-layer breakdown is already floored -- this is the interface this design's      !
         !      new soil<->plant coupling introduces, so it needs the SAME convention applied to it.         !
         !      (root_uptake_layer_b sums to root_uptake_b by construction, so flooring both elementwise     !
         !      preserves that identity: every layer shares the aggregate's sign, per solve_plant_water's    !
         !      own proportional-distribution formula.) --------------------------------------------------!
         root_uptake_b(1:n) = max(root_uptake_b(1:n), 0.0_wp)
         if (ccfg%multilayer_roots) root_uptake_layer_b(1:nsl, 1:n) = max(root_uptake_layer_b(1:nsl, 1:n), 0.0_wp)
         !----- Aggregate root-uptake REQUEST (ground-area units): the plant's OWN computed demand on     !
         !      the soil this step, INCLUDING any storage-refill term (solve_plant_water's root_uptake =   !
         !      (dw_leaf+dw_wood)/dt + transp) -- this, not the raw transpiration demand, is what the       !
         !      soil solve below is asked to supply. ----------------------------------------------------!
         total_uptake_b = sum(root_uptake_b(1:n) * coh%nplant(1:n))
         !----- Under snow the ground water BC changes: only MELTWATER infiltrates (rain went to the    !
         !      pack, snowfall accumulated), and soil evaporation is suppressed by the snow_free_frac     !
         !      area factor below because the snow surface owns the latent flux (sublimation) -- at full  !
         !      cover it is exactly 0, where the old r_aero inflation only tended to 0. Off snow:         !
         !      bare-ground precip/evap. -----------------------------------------------------------------!
         if (snow_exists) then
            hforc%precip_ground   = melt_rate                          ! pack took snow+rain; only meltwater infiltrates
         else
            !----- throughfall_total is the raw forc%precip+forc%snowf sum unless canopy_water_on         !
            !      actually ran the interception sweep above, in which case it is what survived the        !
            !      canopy screen -- reduces to the OLD line exactly when the feature is off. --------------!
            hforc%precip_ground   = throughfall_total
         end if
         !----- Leaf/root-turnover shed water (P4): a PATCH-level input (bio%shed_water_rate, frozen   !
         !      daily by meds_vegetation_dynamics/meds_fast_dynamics -- not atmospheric forcing, so it  !
         !      does NOT live on forc) reaches the ground directly, bypassing interception entirely --  !
         !      regardless of snow cover (a rare, secondary combination not worth its own snow-pack       !
         !      pathway). bio%shed_water_rate is exactly 0 whenever no cohort shed water this slow        !
         !      step, so this is a no-op on every existing test. ------------------------------------!
         hforc%precip_ground = hforc%precip_ground + bio%shed_water_rate
         !----- Ground evaporation is an AREA-weighted tile flux: the BARE-SOIL resistance (1/ggnet)      !
         !      with an explicit (1-snowfac) area factor. It used to be throttled by dividing r_aero by   !
         !      (1-snowfac) instead, which scales only the AERODYNAMIC leg of the series and so is        !
         !      equivalent ONLY when the dry-surface-layer resistance r_soil is 0: with r_soil > 0 that    !
         !      form gives (1-f)*rho*dq/(1/g + (1-f)*r_soil) against the tile's (1-f)*rho*dq/(1/g+r_soil), !
         !      always the larger, so it over-predicted bare-soil evaporation at PARTIAL cover -- worst    !
         !      over a dry surface, where r_soil dominates. The two agree at snowfac = 0 and 1, so the     !
         !      snow-free path (every golden anchor) is bit-identical. This also matters downstream: g_top !
         !      area-weights radiation and sensible heat by (1-snowfac) but subtracts le_soil WHOLE, which !
         !      is only self-consistent once le_soil is itself an area-integrated tile flux. -------------!
         hforc%r_aero         = 1.0_wp / max(aero%ggnet, tiny_num)
         hforc%snow_free_frac = 1.0_wp - snowfac_col
         !----- Per-layer soil root sink: distribute the plant's aggregate REQUEST (total_uptake_b, in    !
         !       place of the raw transp demand) by the previous step's actual per-layer uptake shares      !
         !       when coupled (so the soil dries where roots take water), else the static root-fraction      !
         !       profile. Both share sets sum to 1 => column-total water balance (and every budget) is       !
         !       unchanged; only the vertical distribution differs. -------------------------------------!
         if (ccfg%multilayer_roots .and. sum(bio%root_sink_share(1:nsl)) > tiny_num) then
            hforc%root_uptake(1:nsl) = total_uptake_b * bio%root_sink_share(1:nsl)
         else
            hforc%root_uptake(1:nsl) = total_uptake_b * ccfg%soil%root_frac(1:nsl)
         end if
         hforc%t_ground           = t_ground
         hforc%q_air              = qcas
         hforc%rho_air            = rho
         call column_hydrology_flux(bio%soil_w, hforc, ccfg%soil, ccfg%hydro, dt_fast, hflux)
         !----- §5.1 process mask: a FROZEN soil-water column keeps supplying psi_soil / uptake / soil    !
         !      evaporation to its neighbours, but theta itself does not advance. The Picard state^n      !
         !      snapshot IS the freeze source -- no extra buffer needed. ----------------------------------!
         if (.not. ccfg%mask%soil_water) bio%soil_w = soil_w_n
         budg%soil_nsub = hflux%nsub               ! section 5.3 work counter (same seam on both schemes)
         budg%integ_nsteps = 1_ik ; budg%integ_nrej = 0_ik   ! the split takes exactly one step per dt_fast
         soil_evap = hflux%soil_evap                                     ! §3.6: THE ground latent authority
         !----- Soil-limiting rescale (MEDS_ED2_RK45_DESIGN.md sec 3): when the soil cannot honour the    !
         !      full request (e.g. near wilting point, meds_soil_water's own fwilt ramp), hflux%           !
         !      uptake_total < total_uptake_b -- scale down ONLY the credit applied to wood_water_mass      !
         !      below (the mass update, not sapflow/transp, which are internal to the plant) so the         !
         !      whole-column water budget still closes to the soil's TRUE realized supply. In the common    !
         !      (non-limited) case scale == 1 EXACTLY (column_hydrology_flux's own conservative              !
         !      accumulation), so this is a no-op there. ------------------------------------------------!
         scale = 1.0_wp
         if (total_uptake_b > tiny_num) scale = min(1.0_wp, hflux%uptake_total / total_uptake_b)
         if (ccfg%multilayer_roots) bio%root_sink_share(1:nsl) = 0.0_wp   ! re-accumulate this step's uptake shares
         if (ccfg%multilayer_roots) then                   ! accumulate per-layer uptake -> next step's sink shares
            do i = 1_ik, n
               do k = 1_ik, nsl
                  bio%root_sink_share(k) = bio%root_sink_share(k)                                 &
                                         + max(root_uptake_layer_b(k, i), 0.0_wp) * coh%nplant(i)
               end do
            end do
         end if
         if (ccfg%multilayer_roots) then                      ! normalize to shares (sum = 1); all-zero => root_frac
            sink_tot = sum(bio%root_sink_share(1:nsl))
            if (sink_tot > tiny_num) bio%root_sink_share(1:nsl) = bio%root_sink_share(1:nsl) / sink_tot
         end if

         !----- 3c. GROUND surface = snowfac-BLENDED snow + (1-snowfac) bare soil (design §4f/§4g/§6). The !
         !      snow terms (h_snow_s / le_snow_s / g_base_snow, computed operator-split in 2b) are already   !
         !      snowfac-weighted; the bare-soil terms are scaled by (1-snowfac). soil_evap carries its OWN   !
         !      (1-snowfac) AREA factor (hforc%snow_free_frac above), so le_soil is already an area-         !
         !      integrated tile flux and is subtracted whole -- do NOT scale it again here. snowfac=0        !
         !      reduces to the bare-soil skin exactly. ---------------------------------------------------!
         call ground_surface_fluxes(t_ground, tcas, aero%ggnet, rho, soil_evap, h_bare, le_soil)
         h_ground  = h_snow_s + (1.0_wp - snowfac_col) * h_bare
         le_ground = le_snow_s + le_soil
         g_top     = g_base_snow + (1.0_wp - snowfac_col) * (forc%abs_sw_ground + forc%abs_lw_ground - h_bare) - le_soil

         !----- 3d. CAS three-twin update: IMPLICIT atm exchange, FROM the state^n snapshot. ----!
         src_enth = coh_h + coh_qw + h_ground + le_ground                 ! [W/m2]  sensible + latent
         src_vap  = coh_transp + coh_film_evap + soil_evap + subl_mass / dt_fast   ! + snow sublimation vapour (0 off snow)
         cas_src%surface_enthalpy_source = src_enth
         cas_src%surface_vapor_source    = src_vap
         cas_src%biotic_co2_source       = nee_biotic                     ! passive CO2 twin (frozen source)
         call cas_column_step_implicit(enth0, shv0, co20, cas_src, cas_col, dt_fast, enth1, shv1, co21)
         !----- §5.1 process mask: the three CAS twins are solved together but freeze INDEPENDENTLY, so    !
         !      "no canopy-air temperature dynamics" and "no CO2 dynamics" are separable axes. -------------!
         if (.not. ccfg%mask%cas_energy) enth1 = enth0
         if (.not. ccfg%mask%cas_vapour) shv1  = shv0
         if (.not. ccfg%mask%cas_co2)    co21  = co20
         tcas_new = cas_temp_of_enthalpy(enth1, shv1)

         !----- 3d'. SOIL THERMAL step (P3b): reset to state^n, apply the water-enthalpy boundary  !
         !          advection (at THIS pass's t_ground/t_bot, so the budget uses the same values),   !
         !          then the BE heat diffusion with g_top; DIAGNOSE the surface/bottom temps for the  !
         !          NEXT pass. On the converged exit t_ground/t_bot stay at the values the advection   !
         !          + budget used (so split@niter=1 is bit-identical to the old post-loop step).       !
         if (iter > 1_ik) bio%soil_e = soil_e_n
         !----- BOUNDARY water enthalpy is now the KERNEL's job (eforc%w_flux_top/bot below), not an   !
         !      ad-hoc pre-addition to soil_energy here. The old form evaluated the inflow on state^n   !
         !      at rain_temp while the kernel's interior outflow used the post-conduction T^{n+1}: with !
         !      each face term ~1300 W/m2 (internal_energy_liquid ~1.0 MJ/kg) and the physical signal    !
         !      their ~30 W/m2 difference, that time-level split was larger than what it was computing.  !
         !      runoff_surf is also gone: it leaves the PONDING store, never entering the soil, so       !
         !      debiting it from layer 1 removed energy the layer never received. -------------------------!
         eforc%g_top = g_top ; eforc%geothermal = 0.0_wp
         !----- Up-positive: infiltration and drainage are both DOWNWARD, hence negative. Evaporation    !
         !      is excluded from the top water flux -- its enthalpy leaves as vapour inside g_top. -------!
         eforc%w_flux_top  = -hflux%infiltration / rho_h2o
         eforc%t_water_top = rain_temp
         !----- The BOTTOM face stays an explicit driver term at t_bot (below), NOT a kernel face.     !
         !      Only the TOP face was mis-timed: it is the one paired against the interior advection    !
         !      inside the same layer. Moving drainage too would re-base it from t_bot onto the kernel's !
         !      post-conduction T^{n+1}, which the whole-column ledger (e_out, sec 7b) evaluates at      !
         !      t_bot -- a gratuitous mismatch in a term that was never wrong. ------------------------!
         eforc%w_flux_bot  = 0.0_wp
         !----- Drainage OUT at the bottom. ------------------------------------------------------------!
         !      There is deliberately NO runoff term here any more. The PONDING store (col%w_surface) is  !
         !      a mass buffer with no enthalpy state of its own, so MEDS books water crossing its         !
         !      boundary at the temperature of the store it LEFT: rain that ponds never deposits its      !
         !      enthalpy in the column (the ledger's e_in counts only the infiltrated share), and runoff  !
         !      out of the pond therefore has no enthalpy to remove. Debiting it from soil layer 1 took   !
         !      out ~1 MJ per kg the layer never received -- with w_pond_max = 5 kg/m2 that is up to      !
         !      ~16 K of spurious cooling in a 0.1 m layer, unbounded under sustained runoff. The         !
         !      whole-column ledger closed only because e_out repeated the same mistake; both go          !
         !      together (sec 7b). KNOWN APPROXIMATION: ponded water that sits for several steps and      !
         !      later infiltrates still enters at rain_temp rather than at a relaxed pond temperature.    !
         !      Bounded by w_pond_max; removing it needs a prognostic pond enthalpy. -------------------!
         bio%soil_e%soil_energy(nsl) = bio%soil_e%soil_energy(nsl)                                  &
              - hflux%drainage * internal_energy_liquid(t_bot) * dt_fast / ccfg%soil%dz(nsl)
         !----- POST-SOLVE mass corrections the hydrology applied with NO face (meds_soil_water's        !
         !      saturation clip and theta_res hard floor). soil_energy_step_implicit consumes the         !
         !      CORRECTED theta but advects enthalpy only on the faces, so without this the corrected     !
         !      mass arrives/departs with no enthalpy and the whole discrepancy lands in the diagnosed    !
         !      temperature -- ~3 K per kg/m2 in a 0.1 m layer, since internal_energy_liquid is ~1.0      !
         !      MJ/kg in absolute terms. Move each layer's enthalpy at ITS OWN temperature so the         !
         !      correction is temperature-NEUTRAL (a numerical fix-up must not heat or cool the soil);    !
         !      soil_temp is state^n's read-off, which is what the pre-step energy inverts to. The two    !
         !      GIVE-BACKs are NOT here -- they undo part of the root sink, whose enthalpy counterpart    !
         !      is the deferred leaf-temperature approximation below. ----------------------------------!
         e_clip = 0.0_wp ; e_floor = 0.0_wp
         do k = 1_ik, nsl
            u_corr = internal_energy_liquid(bio%soil_e%soil_temp(k))
            bio%soil_e%soil_energy(k) = bio%soil_e%soil_energy(k)                                   &
                 + (hflux%floor_layer(k) - hflux%clip_layer(k)) * u_corr * dt_fast / ccfg%soil%dz(k)
            e_clip  = e_clip  + hflux%clip_layer(k)  * u_corr     ! -> pond, i.e. OUT of the column
            e_floor = e_floor + hflux%floor_layer(k) * u_corr     ! created with the floored mass, IN
         end do
         eforc%soil_water(1:nsl) = bio%soil_w%theta(1:nsl)
         if (ccfg%advect_soil_heat) then
            eforc%w_flux(1:nsl) = -hflux%w_flux(1:nsl)      ! down-positive (hydro) -> up-positive (energy)
         else
            eforc%w_flux(1:nsl) = 0.0_wp                    ! interior advection lumped (validated baseline)
         end if
         !----- Shed transpiration-water enthalpy, distributed over the SAME per-layer profile the root   !
         !      MASS sink used above (root_sink_share when multilayer_roots is on, else root_frac). It    !
         !      used to be pinned to root_frac unconditionally, so with multilayer roots the soil shed    !
         !      energy from layers the water had not left. Both profiles are normalized to sum to 1, so   !
         !      the column total -- and every budget -- is unchanged; only the vertical placement moves,  !
         !      and the single-layer/default path is bit-identical. --------------------------------------!
         if (ccfg%multilayer_roots .and. sum(bio%root_sink_share(1:nsl)) > tiny_num) then
            eforc%root_heat_sink(1:nsl) = coh_qsoil * bio%root_sink_share(1:nsl)
         else
            eforc%root_heat_sink(1:nsl) = coh_qsoil * ccfg%soil%root_frac(1:nsl)
         end if
         !----- KNOWN DEFERRED APPROXIMATION (MEDS_ED2_RK45_DESIGN.md sec 2, "part 2"): the MAGNITUDE     !
         !      above is still coh_qsoil -- liquid enthalpy at LEAF temperature for the FULL transpiration !
         !      mass -- rather than the soil's realized hflux%uptake_total at each layer's OWN temperature. !
         !      Neither axis can be fixed here alone, and the reason is concrete: rescaling by             !
         !      uptake_total/total_uptake_b would leave (1-scale)*coh_qsoil unmatched in whole_energy,      !
         !      because the water the plant draws from storage instead moves wood_water_mass -- which is    !
         !      a WATER store only. The energy ledger's wood term is built from wmass_w = bsap*C2B*         !
         !      WOOD_MOIST_FRAC, a fixed structural fraction, so plant water storage carries no enthalpy    !
         !      to supply the difference. Re-basing onto soil temperature has the same problem from the     !
         !      other end. Both need the upwind-temperature qloss/qwflux_wl chain sec 2 describes (soil ->  !
         !      wood -> leaf, each at its own temperature), which is where they belong. -------------------!
         call soil_energy_step_implicit(bio%soil_e, eforc, ccfg%soil_thermal, ccfg%soil, ccfg%energy, dt_fast, sflux)
         !----- §5.1 process mask: hold the soil column at state^n, so its temperature acts as a CONSTANT   !
         !      lower boundary for the surface system. Restoring the whole store also restores the          !
         !      soil_temp/soil_fliq read-offs, which stay consistent with the frozen energy. ---------------!
         if (.not. ccfg%mask%soil_heat) bio%soil_e = soil_e_n
         t_ground_dia = bio%soil_e%soil_temp(1) ; t_bot_dia = bio%soil_e%soil_temp(nsl)

         !----- 3e. Convergence: inter-iterate temperature (CAS + leaf + ground) + CAS humidity. -!
         resid_T = abs(tcas_new - tcas_in)
         resid_T = max(resid_T, abs(t_ground_dia - t_ground))
         do i = 1_ik, n
            if (picard .and. coh%lai(i) < LAI_SLAVE_MIN) cycle    ! slaved cohorts excluded from the norm
            resid_T = max(resid_T, abs(bio%leaf_temp(i) - leaf_in(i)))
         end do
         if (.not. picard) then
            nconv = .true. ; exit                                 ! split: single pass is the answer
         else if (.not. ccfg%picard_fixed_iter .and. resid_T < ccfg%picard_tol_temp                &
                  .and. abs(shv1 - qcas_in) < ccfg%picard_tol_shv) then
            nconv = .true. ; exit                                 ! converged: keep this pass's t_ground/t_bot
         end if
         !----- Seed the NEXT pass only (under-relax the CAS seed; committed enth1/shv1 stay exact;  !
         !      t_ground/t_bot take the freshly diagnosed values). Skipped on the LAST pass of a      !
         !      non-converged / picard_fixed_iter run so the post-loop budget + the last advection    !
         !      reference the SAME t_ground/t_bot (else a small energy asymmetry on those exits).      !
         if (iter < niter) then
            tcas = ccfg%picard_relax * tcas_new + (1.0_wp - ccfg%picard_relax) * tcas_in
            qcas = ccfg%picard_relax * shv1     + (1.0_wp - ccfg%picard_relax) * qcas_in
            t_ground = t_ground_dia ; t_bot = t_bot_dia
         end if
      end do
      if (picard .and. ccfg%picard_fixed_iter) nconv = .true.    ! fixed-count run: accept the last iterate

      !----- Commit the CAS (enth1/shv1/co21 = exact BE box solution at convergence, from the       !
      !      shared cas_column_step_implicit kernel; co21 is frozen-input so any pass gives it). ----!
      bio%cas%can_enthalpy = enth1 ; bio%cas%can_shv = shv1 ; bio%cas%can_co2 = co21
      bio%cas%can_temp     = cas_temp_of_enthalpy(enth1, shv1)

      !----- Commit the plant internal water MASS (MEDS_ED2_RK45_DESIGN.md sec 1/4/5): the frozen-    !
      !      averaged flux Euler update, using the LAST Picard pass's sapflow_b/root_uptake_b/scale/    !
      !      transp_pp (arrays declared outside the loop, so they hold that pass's values here). This   !
      !      is algebraically IDENTICAL to converting the matrix-exponential's own end-state psi via     !
      !      water_content() when scale==1 (the common case; see solve_plant_water's own dw_l/dw_w       !
      !      bookkeeping) and is what makes the soil-limited edge case (scale<1) closure-consistent:      !
      !      sapflow (an internal wood->leaf transfer) is untouched by any soil-supply shortfall; only    !
      !      the wood<->soil boundary credit is rescaled to what the soil actually gave up. §5.1 process   !
      !      mask: frozen hydraulics still transpires/takes up water (sapflow/uptake already fed the CAS   !
      !      and the soil sink above) but the tissue mass stays put -- the stiffest column mode removed.   !
      if (ccfg%mask%hydraulics) then
         bio%leaf_water_mass(1:n) = bio%leaf_water_mass(1:n) + dt_fast*(sapflow_b(1:n) - transp_pp(1:n))
         bio%wood_water_mass(1:n) = bio%wood_water_mass(1:n)                                             &
                                   + dt_fast*(root_uptake_b(1:n)*scale - sapflow_b(1:n))
      end if
      !----- Commit the surface (interception film) water debit (sec 3.4, P1): the SAME LAST-pass         !
      !      film_evap/film_evap_w that fed coh_qw above, so the CAS credit and this store debit are        !
      !      the identical number by construction -- already clamped to <= the water on hand each pass,     !
      !      so this can never go negative (the max(0,...) floor is defensive only). Dew (film_evap<0) is   !
      !      NOT pre-clamped to the store's capacity (see the veg_energy_diagnostic call site above for      !
      !      why -- clamping it there desyncs dh/drnet from the credited flux, a WORSE inconsistency), so    !
      !      it can transiently push a store above its own dewmx*lai/wai capacity here; the excess is        !
      !      explicitly bookkept below (sec 9: "clamp, don't silently over-apply") rather than left for       !
      !      intercept_canopy_layer's own w_max ceiling to silently discard on a later step. ------------------!
      surf_overflow = 0.0_wp
      if (ccfg%canopy_water_on) then
         bio%leaf_surf_water(1:n) = max(0.0_wp, bio%leaf_surf_water(1:n) - film_evap(1:n)   * dt_fast)
         bio%wood_surf_water(1:n) = max(0.0_wp, bio%wood_surf_water(1:n) - film_evap_w(1:n) * dt_fast)
         do i = 1_ik, n
            leaf_cap_i = ccfg%hydro%dewmx * coh%lai(i) ; wood_cap_i = ccfg%hydro%dewmx * coh%wai(i)
            surf_overflow = surf_overflow + max(0.0_wp, bio%leaf_surf_water(i) - leaf_cap_i)          &
                                           + max(0.0_wp, bio%wood_surf_water(i) - wood_cap_i)
            bio%leaf_surf_water(i) = min(bio%leaf_surf_water(i), leaf_cap_i)
            bio%wood_surf_water(i) = min(bio%wood_surf_water(i), wood_cap_i)
         end do
      end if

      !----- Picard diagnostics + non-convergence contract (clamp = last iterate, never partial). !
      budg%picard_iters       = max(budg%picard_iters, niter_taken)
      budg%picard_worst_resid = max(budg%picard_worst_resid, resid_T)
      if (picard .and. .not. nconv) then
         budg%picard_nonconv = budg%picard_nonconv + 1_ik
         if (ccfg%energy%debug_error) error stop 'column_fast_step: Picard did not converge (debug_error)'
      end if
      if (present(converged)) converged = nconv
      if (present(iters))     iters     = niter_taken

      !----- The soil thermal step now runs INSIDE the Picard loop (§3d'); bio%soil_e, sflux, the   !
      !      converged g_top, and t_ground/t_bot (the values the last pass's advection used) are all  !
      !      final here, so the budgets below close against the consistent boundary fluxes.           !

      !----- 7. Per-kernel closed budgets (each closes by construction). L2/debug_error mode      !
      !      (ccfg%energy%debug_error) promotes a non-closing budget from a silently-counted       !
      !      n_fail to a hard `error stop` -- the enforced half of the conservation check (plan     !
      !      MEDS_NUMERICS_SCOPING.md sec 4/QW2); off by default so production behaviour is         !
      !      unchanged. Each check reuses the accumulator's own resid/scale/rtol/atol (budget_      !
      !      accumulate/track_resid just set %resid as a side effect), so the tolerance can never    !
      !      drift out of sync between the soft count and the hard gate. ------------------------!
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, src_enth + gah*forc%enthalpy_atm, &
                             gah*enth1, dt_fast, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_check_stop(budg%cas_energy%resid, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp,        &
                             'cas_energy (split)', halt_budgets)
      call budget_accumulate(budg%cas_water,  wcap*shv0, wcap*shv1, src_vap + gaw*forc%shv_atm,        &
                             gaw*shv1, dt_fast, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_check_stop(budg%cas_water%resid, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp,      &
                             1.0e-10_wp, 'cas_water (split)', halt_budgets)
      call budget_accumulate(budg%cas_co2,    ccap*co20, ccap*co21, nee_biotic + gac*forc%co2_atm,&
                             gac*co21, dt_fast, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      call budget_check_stop(budg%cas_co2%resid, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp,             &
                             'cas_co2 (split)', halt_budgets)
      call track_resid(budg%soil_energy, sflux%energy_resid, abs(g_top)*dt_fast + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp)
      call budget_check_stop(budg%soil_energy%resid, abs(g_top)*dt_fast + 1.0_wp, 1.0e-6_wp,       &
                             1.0e-3_wp, 'soil_energy (split)', halt_budgets)
      call track_resid(budg%soil_water,  hflux%mass_resid,   1.0_wp,             1.0e-6_wp, 1.0e-4_wp)
      call budget_check_stop(budg%soil_water%resid, 1.0_wp, 1.0e-6_wp, 1.0e-4_wp,                  &
                             'soil_water (split)', halt_budgets)

      !----- 7b. WHOLE-COLUMN budgets: Δ(all stores) vs the TRUE boundary fluxes (catches leaks). !
      e_soil1 = 0.0_wp ; w_soil1 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil1 = e_soil1 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil1 = w_soil1 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      !----- Plant internal water AFTER the mass commit above -- joins the stores (MEDS_ED2_RK45_       !
      !      DESIGN.md sec 3): uptake/transpiration are already boundary transfers through the soil/CAS   !
      !      terms below, so this only needs to add the plant's OWN storage change, not a new boundary.   !
      w_plant1 = sum(coh%nplant(1:n) * (bio%leaf_water_mass(1:n) + bio%wood_water_mass(1:n)))
      !----- Surface (interception film) water AFTER the film-evap commit above -- the TRUE end-of-step  !
      !      value (already ground-area-referenced, so no nplant factor, unlike w_plant1 just above).      !
      !      Its enthalpy at the SAME fixed rain_temp reference as surf_enth0 (sec 2c above). Gated behind !
      !      the same flag as surf_water0 above -- exactly 0.0 whenever the feature is off. ---------------!
      surf_water1 = 0.0_wp ; surf_enth1 = 0.0_wp
      if (ccfg%canopy_water_on) then
         surf_water_tmp(1:n) = bio%leaf_surf_water(1:n) + bio%wood_surf_water(1:n)
         surf_water1 = sum(surf_water_tmp(1:n))
         surf_enth1  = surf_water1 * internal_energy_liquid(rain_temp)
      end if
      !----- Water: precip (rain + snow) IN; drainage + runoff + atm-vapour OUT. The SNOW store (swe)  !
      !      joins the stores; sublimation + melt are INTERNAL transfers (snow<->CAS/soil) that          !
      !      telescope with the CAS-vapour / soil-water changes, so they are NOT boundary terms. The       !
      !      surface water store (surf_water0/1) joins too -- interception/film-evap are ALSO internal      !
      !      transfers (precip already counts the FULL amount in w_in; evaporation telescopes with the       !
      !      CAS-vapour change) -- EXCEPT the dew-capacity overflow (surf_overflow, sec 5.3/9 above): that    !
      !      water is real (it left the CAS as vapour and was credited to the surface store), but this MVP    !
      !      does not model where excess-beyond-capacity condensate physically goes (drip is a candidate      !
      !      for a future pass), so it is booked as leaving the column boundary here rather than silently     !
      !      vanishing into intercept_canopy_layer's own w_max ceiling on a later step. surf_overflow is       !
      !      exactly 0 whenever canopy_water_on is off, so this is a no-op on the byte-identical default path.!
      w_in  = forc%precip + forc%snowf + bio%shed_water_rate   ! P4: shed water (patch-level) is a boundary input too
      w_out = hflux%drainage + hflux%runoff_surf + gaw * (shv1 - forc%shv_atm) + surf_overflow / dt_fast
      call budget_accumulate(budg%whole_water, w_soil0 + wcap*shv0 + swe0_s + w_plant0 + surf_water0, &
                             w_soil1 + wcap*shv1 + swe1_s + w_plant1 + surf_water1,                    &
                             w_in, w_out, dt_fast, max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
      call budget_check_stop(budg%whole_water%resid, max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp,  &
                             1.0e-4_wp, 'whole_water (split)', halt_budgets)
      !----- Energy: net ground+cohort radiation + precip enthalpy IN; atm exchange + drainage/runoff   !
      !      enthalpy OUT. Under snow the ground radiation is the snow's emission-corrected net input     !
      !      (ground_rad_col = sfx%rnet) and the precip enthalpy that entered the pack (snow_acc_enth) is  !
      !      a boundary input; the snow store (snow_e0/1) joins the stores so sublimation/melt/base-       !
      !      conduction telescope with the CAS/soil changes. Off snow this reduces to the bare-soil ledger. !
      !      The INTERCEPTED share of precip (routed to the surface store instead of infiltration) needs    !
      !      its own e_in term at the same rain_temp; intercepted_total is 0 whenever canopy_water_on is     !
      !      off (throughfall_total defaults to the full precip+snowf sum there), so this line is a no-op    !
      !      on the byte-identical default path. -----------------------------------------------------------!
      !----- e_floor / e_clip pair with the post-solve mass corrections applied at sec 3d' above: clip    !
      !      water leaves the column for the (enthalpy-free) pond, floored water is created with the      !
      !      mass. runoff carries NO enthalpy term -- the pond never received any (see 3d'). -------------!
      e_in  = coh_rnet + ground_rad_col + snow_acc_enth / dt_fast                                      &
              + hflux%infiltration * internal_energy_liquid(rain_temp)                                &
              + intercepted_total * internal_energy_liquid(rain_temp)   ! (0 under snow: rain_temp=tsupercool_liq)
      e_in  = e_in + e_floor
      e_out = gah * (enth1 - forc%enthalpy_atm) + hflux%drainage * internal_energy_liquid(t_bot)      &
              + surf_overflow / dt_fast * internal_energy_liquid(rain_temp)   ! pairs with w_out's surf_overflow
      e_out = e_out + e_clip
      !----- Prognostic wood/leaf + snow + surface water are real energy STORES: add their deltas to the  !
      !      ledger. All telescope to 0 when inactive (stores unchanged), so the split golden anchor is     !
      !      preserved. KNOWN DEFERRED IMPRECISION (mirrors the P0 root_heat_sink note): surf_enth0/1 use    !
      !      ONE fixed reference temperature (rain_temp) rather than a real prognostic surface-water         !
      !      temperature, so the residual proportional to (leaf/wood dt_temp) is not exactly zero -- the     !
      !      same category of upwind-temperature approximation sec 2's qloss/qwflux_wl coupling is meant     !
      !      to eventually resolve properly, deferred with the rest of the energy-advection work. -----------!
      call budget_accumulate(budg%whole_energy,                                                          &
                             e_soil0 + wcap*enth0 + wood_store0 + leaf_store0 + snow_e0 + surf_enth0,      &
                             e_soil1 + wcap*enth1 + wood_store1 + leaf_store1 + snow_e1 + surf_enth1,      &
                             e_in, e_out, dt_fast, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp)
      call budget_check_stop(budg%whole_energy%resid, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp, &
                             'whole_energy (split)', halt_budgets)
      !----- ET diagnostic: the CAS->atm latent-heat flux (matches the whole_water vapour OUT term). --!
      if (present(le_flux)) le_flux = gaw * (shv1 - forc%shv_atm) * latent_heat_vap
      !----- Sensible-heat diagnostic: CAS->atm flux via the heat conductance gah. ------------------!
      if (present(h_flux))  h_flux  = gah * (bio%cas%can_temp - aenv%theta_atm) * cp_air
   end subroutine column_fast_step
   !----- Track a kernel's own closed-budget residual into a budget_t (worst + fail count). ---!
   pure subroutine track_resid(b, resid, scale, rtol, atol)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: resid, scale, rtol, atol
      b%n_check = b%n_check + 1_ik
      b%resid   = resid
      b%worst   = max(b%worst, abs(resid))
      if (.not. closure_ok(resid, scale, rtol, atol)) b%n_fail = b%n_fail + 1_ik
   end subroutine track_resid

   !----- A committed CAS or soil-layer temperature within 5 K of the vegetation/soil-energy clamp     !
   !      bounds ([180,350] K, meds_fast_ark's clamp_cas/clamp_soil_energy) is never a physical         !
   !      land-surface state -- it is the fingerprint of the explicit RK45 surface instability railing  !
   !      to the clamp (MEDS_ED2_RK45_DESIGN.md P6). The dispatcher uses it to trigger the split rescue. !
   !      The soil freeze/thaw plateau (t_3ple = 273.16 K) is well inside the window, so a genuinely     !
   !      freezing soil is NOT flagged; only the clamp-pinned 180/350 K rails are. --------------------!
   pure function rk45_state_railed(bio, nsl) result(railed)
      type(patch_biophys_t), intent(in) :: bio
      integer(ik),           intent(in) :: nsl
      logical :: railed
      real(wp), parameter :: T_LO = 185.0_wp, T_HI = 345.0_wp
      railed =      bio%cas%can_temp <= T_LO .or. bio%cas%can_temp >= T_HI                          &
               .or. any(bio%soil_e%soil_temp(1:nsl) <= T_LO)                                        &
               .or. any(bio%soil_e%soil_temp(1:nsl) >= T_HI)
   end function rk45_state_railed

end module meds_fast_split
