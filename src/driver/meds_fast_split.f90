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
! enthalpy on it -- an OPT-IN coupling (cfg%advect_soil_heat, default OFF).                                  !
!==========================================================================================!
module meds_fast_split
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, cp_air, latent_heat_vap, rho_h2o, pi,               &
                                     tsupercool_liq, grav_head
   use meds_plant_hydraulics, only : rhizosphere_cond
   use meds_hydr_lib, only : soil_hydr_cond_from_theta
   use meds_config,           only : meds_config_t, hydraulics_config_t,                          &
                                     SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED,               &
                                     INTEG_SPLIT, INTEG_ARK
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
                                     apply_hydraulics_config
   use meds_fast_ark,         only : column_fast_step_ark, aero_bottom_to_top, column_prepass
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_soil_energy,      only : soil_energy_step_implicit
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_vegetation_biophysics, only : veg_energy_diagnostic, veg_energy_step_implicit,      &
                                     sensible_heat_coeff, lw_emission_slope, le_conductance_flux
   use meds_soil_water,       only : column_hydrology_flux
   use meds_ground_biophysics, only : snow_energy_step, snow_base_conductance,                  &
                                     snow_accumulate, snow_drain_meltwater, snow_cover_fraction, &
                                     ground_surface_fluxes
   use meds_plant_interface,  only : hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     solve_plant_water, N_HYDRO
   use meds_therm_lib,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_liquid,  &
                                     temp_to_uext
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok
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
      type(hydro_env_t)      :: henv
      type(hydro_flux_t)     :: hfx
      real(wp)               :: transp_c(coh%n)     !< [kg/m2 ground/s] per-cohort transpiration demand (automatic)
      real(wp)               :: h_coeff_f(coh%n), g_tr_f(coh%n), leaf_in(coh%n)   !< frozen coeffs + prev-iterate leaf temp
      real(wp)               :: t_emit(coh%n)      !< LW emission base (start leaf_temp; matches the RT tcan_bt, P3c)
      real(wp)               :: wood_emit(coh%n)   !< start-of-sub-step wood temp (prognostic-wood seed; Picard-correct)
      real(wp)    :: te
      type(soil_column_t)        :: soil_w_n        !< snapshot of the soil-water column at state^n (Picard reset)
      type(soil_energy_column_t) :: soil_e_n        !< snapshot of the soil thermal column at state^n (Picard reset)
      real(wp), allocatable  :: psi_n(:,:)          !< snapshot of the per-cohort plant water potentials at state^n
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
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet
      real(wp)    :: nee_biotic
      real(wp)    :: t_ground, t_bot, g_top, h_ground, le_ground, soil_evap, rain_temp
      real(wp)    :: gah, gaw, gac, wcap, ccap, src_enth, src_vap, src_frac
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
      integer(ik) :: i, n, nsl, k

      n = coh%n ; nsl = ccfg%soil%n_active

      !----- Prognostic wood (P2) and prognostic leaf (P3) both advance own stores on the SPLIT/PICARD    !
      !      path via the BE cap/dt term; the ARK arrowhead couplings (P4) are deferred, so gate those.    !
      if (ccfg%wood_energy_model == WOODEN_PROGNOSTIC .and. cfg%time_integrator == INTEG_ARK) &
         error stop 'column_fast_step: wood_energy_model="prognostic" under ARK is deferred (P2 split-only)'
      if (ccfg%leaf_energy_model == LEAFEN_PROGNOSTIC .and. cfg%time_integrator == INTEG_ARK) &
         error stop 'column_fast_step: leaf_energy_model="prognostic" under ARK is deferred (P4 arrowhead)'

      !----- TIME-INTEGRATOR dispatch (inserted BEFORE the first bio mutation, so the split path       !
      !      below is byte-for-byte unentered -- the golden anchor is preserved structurally). The     !
      !      coupled IMEX-ARK path is opt-in ([fast].time_integrator="ark"); default is the split.     !
      if (cfg%time_integrator == INTEG_ARK) then
         call column_fast_step_ark(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,   &
                                   gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
         if (present(le_flux)) le_flux = aenv%rho_air * aero%ustar * aero%temp2                    &
                                         * (bio%cas%can_shv - forc%shv_atm) * latent_heat_vap
         if (present(h_flux))  h_flux  = aenv%rho_air * aero%ustar * aero%temp2                    &
                                         * (bio%cas%can_temp - aenv%theta_atm) * cp_air
         return
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

      !----- Snapshot state^n once: the Picard passes re-solve the SAME backward-Euler steps FROM  !
      !      this base each iteration (only the source, re-evaluated at the iterate, changes). The   !
      !      soil-water column + plant-psi are prognostic and column_hydrology_flux / solve_plant_water !
      !      ADVANCE them, so they must be reset to state^n before each re-solve or they double-step.  !
      enth0 = bio%cas%can_enthalpy ; shv0 = bio%cas%can_shv ; co20 = bio%cas%can_co2
      soil_w_n = bio%soil_w ; soil_e_n = bio%soil_e ; psi_n = bio%psi

      !======================================================================================!
      !  3. Outer PICARD fixed point over { leaf energy -> soil water/src_frac -> CAS twins }.   !
      !     niter = 1 under SCHEME_SPLIT_SEQUENTIAL reproduces the operator-split sweep EXACTLY   !
      !     (one pass, no convergence test, soil water solved that pass). Under PICARD the block   !
      !     iterates at the current tcas/qcas until the store temperatures converge; the soil       !
      !     water/src_frac/hydraulics are frozen after pass 1 EVERY pass regardless of               !
      !     ccfg%soil_water_coupling -- that selector is RESERVED for the P3f re-solve-inside-        !
      !     Picard optimization and has NO behavioral effect yet (both SOILH2O_LAGGED and             !
      !     SOILH2O_COUPLED take this same frozen path today).                                        !
      !======================================================================================!
      src_frac = 1.0_wp ; soil_evap = 0.0_wp ; nconv = .false. ; resid_T = 0.0_wp
      niter_taken = 0_ik
      do iter = 1_ik, niter
         niter_taken = iter
         tcas_in = tcas ; qcas_in = qcas ; leaf_in(1:n) = bio%leaf_temp(1:n)

         !----- 3a. Leaf energy balance (diagnostic) at the CURRENT tcas/qcas, frozen coeffs. --!
         qsat_c = sat_specific_humidity(tcas, press)
         dqdt   = sat_specific_humidity_temp_deriv(tcas, press)
         coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
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
            call veg_energy_diagnostic(forc%abs_sw(i), forc%abs_lw(i), h_coeff_f(i), le_slope,       &
                                       lw_slope, le_ref, tcas, te, a_leaf, t_emit(i),                &
                                       dtl, tl, transp_i, dh, drnet)
            bio%leaf_temp(i) = tl
            leaf_store0 = leaf_store0 + cap_leaf * t_emit(i)     ! [J/m2] leaf internal energy (0 K ref; differenced)
            leaf_store1 = leaf_store1 + cap_leaf * tl            ! diagnostic: cap_leaf=0 -> telescopes to 0
            transp_c(i) = transp_i                                                       ! per-cohort demand (hydraulics)
            coh_h      = coh_h      + dh
            coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl)                       ! CAS latent (vapour enthalpy)
            coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)   ! liquid enthalpy soil sheds
            coh_transp = coh_transp + transp_i
            coh_rnet   = coh_rnet   + drnet
            !----- 3a'. Diagnostic WOOD energy balance (own store; own boundary layer + net LW, NO      !
            !      transpiration). Wood sensible + net-LW join coh_h / coh_rnet -> CAS + energy budget.   !
            !      A diagnostic wood has no storage, so absorbed = emitted + sensible-to-CAS -> the        !
            !      coh_rnet and coh_h wood terms are EQUAL (h_coeff_w*dtw) and telescope in the ledger.    !
            if (ccfg%wood_energy_model == WOODEN_DIAGNOSTIC) then
               h_coeff_w  = sensible_heat_coeff(pi * coh%wai(i), aero%wood_gbh(i), rho, cp_air)
               te_w       = tcas
               lw_slope_w = lw_emission_slope(ccfg%veg_thermal%leaf_emiss, te_w, coh%wai(i))
               !----- Diagnostic WOOD = the le_slope = le_ref = 0 case of the same kernel (no transp). --!
               call veg_energy_diagnostic(forc%abs_sw_wood(i), forc%abs_lw_wood(i), h_coeff_w,       &
                                          0.0_wp, lw_slope_w, 0.0_wp, tcas, te_w, 0.0_wp, tcas,      &
                                          dtw, twood, transp_w, dh, drnet)
               bio%wood_temp(i) = twood
               coh_h    = coh_h    + dh
               coh_rnet = coh_rnet + drnet
            else   ! WOODEN_PROGNOSTIC: advance the wood internal-energy store (operator-split, non-stiff). !
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
               wood_store1 = wood_store1 + twood                          ! store energy AFTER the BE step
               bio%wood_temp(i) = wflux%temp
               coh_h    = coh_h    + wflux%h_flux                         ! wood sensible -> CAS
               coh_qw   = coh_qw   + wflux%qw_flux                        ! wood film-evap -> CAS (0 in MVP)
               coh_rnet = coh_rnet + (forc%abs_sw_wood(i) + forc%abs_lw_wood(i))   ! net wood radiation into the column
            end if
         end do

         !----- 3b. Soil WATER column + supply limiter + plant hydraulics, RE-SOLVED FROM state^n  !
         !          each pass so the realized transpiration (= root uptake) stays consistent with    !
         !          the iterated leaf demand -- freezing it while the demand iterates would leak     !
         !          water/enthalpy. pass 1 is already at state^n; later passes reset the prognostic   !
         !          water + psi to the snapshot before re-advancing.                                  !
         if (iter > 1_ik) then
            bio%soil_w = soil_w_n ; bio%psi = psi_n
         end if
         !----- Under snow the ground water BC changes: only MELTWATER infiltrates (rain went to the    !
         !      pack, snowfall accumulated), and soil evaporation is suppressed (huge r_aero) because     !
         !      the snow surface owns the latent flux (sublimation). Off snow: bare-ground precip/evap.   !
         if (snow_exists) then
            hforc%precip_ground   = melt_rate                          ! pack took snow+rain; only meltwater infiltrates
         else
            hforc%precip_ground   = forc%precip + forc%snowf           ! no pack: sub-threshold snow melts straight in (MVP)
         end if
         hforc%r_aero = 1.0_wp / max((1.0_wp - snowfac_col) * aero%ggnet, tiny_num)  ! only the (1-snowfac) bare fraction evaporates
         !----- Per-layer soil root sink: distribute the SAME total (coh_transp) by the previous step's  !
         !       actual per-layer uptake shares when coupled (so the soil dries where roots take water),  !
         !       else the static root-fraction profile. Both share sets sum to 1 => column-total water    !
         !       balance (and every budget) is unchanged; only the vertical distribution differs. --------!
         if (ccfg%multilayer_roots .and. sum(bio%root_sink_share(1:nsl)) > tiny_num) then
            hforc%root_uptake(1:nsl) = coh_transp * bio%root_sink_share(1:nsl)
         else
            hforc%root_uptake(1:nsl) = coh_transp * ccfg%soil%root_frac(1:nsl)
         end if
         hforc%t_ground           = t_ground
         hforc%q_air              = qcas
         hforc%rho_air            = rho
         call column_hydrology_flux(bio%soil_w, hforc, ccfg%soil, ccfg%hydro, dt_fast, hflux)
         soil_evap = hflux%soil_evap                                     ! §3.6: THE ground latent authority
         src_frac  = 1.0_wp
         if (coh_transp > tiny_num) src_frac = min(1.0_wp, hflux%uptake_total / coh_transp)
         soil_psi_root = root_weighted_psi(hflux%psi_soil, ccfg%soil%root_frac, nsl)
         if (ccfg%multilayer_roots) bio%root_sink_share(1:nsl) = 0.0_wp   ! re-accumulate this step's uptake shares
         do i = 1_ik, n
            henv%transp     = transp_c(i) * src_frac / max(coh%nplant(i), tiny_num)   ! [kg/plant/s]
            if (ccfg%multilayer_roots) then
               !----- Couple to the per-layer soil column: per-layer psi_soil [MPa] from the hydrology     !
               !       solve, and a per-layer rhizosphere conductance (Katul; ksat converted m/s ->         !
               !       kg/m/s/MPa). The plant solver aggregates by conductance and returns per-layer uptake. !
               henv%n_root_layer = nsl
               do k = 1_ik, nsl
                  henv%soil_psi_layer(k)   = hflux%psi_soil(k)
                  henv%root_z_layer(k)     = ccfg%soil%z_node(k)
                  !----- UNSATURATED K(theta) [m/s] (not ksat): dry layers are conductance-down-weighted. !
                  k_theta = soil_hydr_cond_from_theta(ccfg%soil%retention, bio%soil_w%theta(k),                     &
                       ccfg%soil%theta_sat(k), ccfg%soil%theta_res(k), ccfg%soil%vg_alpha(k),            &
                       ccfg%soil%vg_n(k), ccfg%soil%ksat(k))
                  henv%rhizo_cond_layer(k) = rhizosphere_cond(rho_h2o*k_theta/grav_head,                 &
                       coh%broot(i), ccfg%specific_root_area, ccfg%soil%root_frac(k), ccfg%soil%dz(k),   &
                       coh%nplant(i))
               end do
            else
               henv%n_root_layer = 0_ik                       ! scalar BC (bit-identical single-layer path)
               henv%soil_psi     = soil_psi_root
               henv%rhizo_cond   = ccfg%rhizo_cond
            end if
            henv%bleaf      = coh%bleaf(i) ; henv%bsap = coh%bsap(i) ; henv%broot = coh%broot(i)
            henv%sap_area   = coh%sap_area(i) ; henv%height = coh%height(i) ; henv%leaf_area = coh%leaf_area(i)
            call solve_plant_water(henv, ccfg%hydro_p, ccfg%hydro_o, dt_fast, bio%psi(:,i), hfx)
            if (ccfg%multilayer_roots) then                   ! accumulate per-layer uptake -> next step's sink shares
               do k = 1_ik, nsl
                  bio%root_sink_share(k) = bio%root_sink_share(k)                                 &
                                         + max(hfx%root_uptake_layer(k), 0.0_wp) * coh%nplant(i)
               end do
            end if
         end do
         if (ccfg%multilayer_roots) then                      ! normalize to shares (sum = 1); all-zero => root_frac
            sink_tot = sum(bio%root_sink_share(1:nsl))
            if (sink_tot > tiny_num) bio%root_sink_share(1:nsl) = bio%root_sink_share(1:nsl) / sink_tot
         end if
         coh_qw     = coh_qw     * src_frac      ! the CAS gains only the water the soil gave up
         coh_qsoil  = coh_qsoil  * src_frac
         coh_transp = coh_transp * src_frac

         !----- 3c. GROUND surface = snowfac-BLENDED snow + (1-snowfac) bare soil (design §4f/§4g/§6). The !
         !      snow terms (h_snow_s / le_snow_s / g_base_snow, computed operator-split in 2b) are already   !
         !      snowfac-weighted; the bare-soil terms are scaled by (1-snowfac). soil_evap is already        !
         !      (1-snowfac)-scaled by the r_aero above. snowfac=0 reduces to the bare-soil skin exactly.  --!
         call ground_surface_fluxes(t_ground, tcas, aero%ggnet, rho, soil_evap, h_bare, le_soil)
         h_ground  = h_snow_s + (1.0_wp - snowfac_col) * h_bare
         le_ground = le_snow_s + le_soil
         g_top     = g_base_snow + (1.0_wp - snowfac_col) * (forc%abs_sw_ground + forc%abs_lw_ground - h_bare) - le_soil

         !----- 3d. CAS three-twin update: IMPLICIT atm exchange, FROM the state^n snapshot. ----!
         src_enth = coh_h + coh_qw + h_ground + le_ground                 ! [W/m2]  sensible + latent
         src_vap  = coh_transp + soil_evap + subl_mass / dt_fast          ! + snow sublimation vapour (0 off snow)
         cas_src%surface_enthalpy_source = src_enth
         cas_src%surface_vapor_source    = src_vap
         cas_src%biotic_co2_source       = nee_biotic                     ! passive CO2 twin (frozen source)
         call cas_column_step_implicit(enth0, shv0, co20, cas_src, cas_col, dt_fast, enth1, shv1, co21)
         tcas_new = cas_temp_of_enthalpy(enth1, shv1)

         !----- 3d'. SOIL THERMAL step (P3b): reset to state^n, apply the water-enthalpy boundary  !
         !          advection (at THIS pass's t_ground/t_bot, so the budget uses the same values),   !
         !          then the BE heat diffusion with g_top; DIAGNOSE the surface/bottom temps for the  !
         !          NEXT pass. On the converged exit t_ground/t_bot stay at the values the advection   !
         !          + budget used (so split@niter=1 is bit-identical to the old post-loop step).       !
         if (iter > 1_ik) bio%soil_e = soil_e_n
         bio%soil_e%soil_energy(1)   = bio%soil_e%soil_energy(1)                                    &
              + (hflux%infiltration * internal_energy_liquid(rain_temp)                             &
                 - hflux%runoff_surf * internal_energy_liquid(t_ground)) * dt_fast / ccfg%soil%dz(1)
         bio%soil_e%soil_energy(nsl) = bio%soil_e%soil_energy(nsl)                                  &
              - hflux%drainage * internal_energy_liquid(t_bot) * dt_fast / ccfg%soil%dz(nsl)
         eforc%g_top = g_top ; eforc%geothermal = 0.0_wp
         eforc%soil_water(1:nsl) = bio%soil_w%theta(1:nsl)
         if (ccfg%advect_soil_heat) then
            eforc%w_flux(1:nsl) = -hflux%w_flux(1:nsl)      ! down-positive (hydro) -> up-positive (energy)
         else
            eforc%w_flux(1:nsl) = 0.0_wp                    ! interior advection lumped (validated baseline)
         end if
         eforc%root_heat_sink(1:nsl) = coh_qsoil * ccfg%soil%root_frac(1:nsl)   ! shed transpiration-water enthalpy
         call soil_energy_step_implicit(bio%soil_e, eforc, ccfg%soil_thermal, ccfg%soil, ccfg%energy, dt_fast, sflux)
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

      !----- 7. Per-kernel closed budgets (each closes by construction). ----------------------!
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, src_enth + gah*forc%enthalpy_atm, &
                             gah*enth1, dt_fast, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(budg%cas_water,  wcap*shv0, wcap*shv1, src_vap + gaw*forc%shv_atm,        &
                             gaw*shv1, dt_fast, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(budg%cas_co2,    ccap*co20, ccap*co21, nee_biotic + gac*forc%co2_atm,&
                             gac*co21, dt_fast, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      call track_resid(budg%soil_energy, sflux%energy_resid, abs(g_top)*dt_fast + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp)
      call track_resid(budg%soil_water,  hflux%mass_resid,   1.0_wp,             1.0e-6_wp, 1.0e-4_wp)

      !----- 7b. WHOLE-COLUMN budgets: Δ(all stores) vs the TRUE boundary fluxes (catches leaks). !
      e_soil1 = 0.0_wp ; w_soil1 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil1 = e_soil1 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil1 = w_soil1 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      !----- Water: precip (rain + snow) IN; drainage + runoff + atm-vapour OUT. The SNOW store (swe)  !
      !      joins the stores; sublimation + melt are INTERNAL transfers (snow<->CAS/soil) that          !
      !      telescope with the CAS-vapour / soil-water changes, so they are NOT boundary terms.  ------!
      w_in  = forc%precip + forc%snowf
      w_out = hflux%drainage + hflux%runoff_surf + gaw * (shv1 - forc%shv_atm)
      call budget_accumulate(budg%whole_water, w_soil0 + wcap*shv0 + swe0_s, w_soil1 + wcap*shv1 + swe1_s, &
                             w_in, w_out, dt_fast, max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
      !----- Energy: net ground+cohort radiation + precip enthalpy IN; atm exchange + drainage/runoff   !
      !      enthalpy OUT. Under snow the ground radiation is the snow's emission-corrected net input     !
      !      (ground_rad_col = sfx%rnet) and the precip enthalpy that entered the pack (snow_acc_enth) is  !
      !      a boundary input; the snow store (snow_e0/1) joins the stores so sublimation/melt/base-       !
      !      conduction telescope with the CAS/soil changes. Off snow this reduces to the bare-soil ledger.!
      e_in  = coh_rnet + ground_rad_col + snow_acc_enth / dt_fast                                      &
              + hflux%infiltration * internal_energy_liquid(rain_temp)   ! (0 under snow: rain_temp=tsupercool_liq)
      e_out = gah * (enth1 - forc%enthalpy_atm) + hflux%drainage * internal_energy_liquid(t_bot)      &
              + hflux%runoff_surf * internal_energy_liquid(t_ground)
      !----- Prognostic wood/leaf + snow are real energy STORES: add their deltas to the ledger. All      !
      !      telescope to 0 when inactive (stores unchanged), so the split golden anchor is preserved.   !
      call budget_accumulate(budg%whole_energy, e_soil0 + wcap*enth0 + wood_store0 + leaf_store0 + snow_e0, &
                             e_soil1 + wcap*enth1 + wood_store1 + leaf_store1 + snow_e1, e_in, e_out,       &
                             dt_fast, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp)
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

end module meds_fast_split
