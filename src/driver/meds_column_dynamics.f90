!==========================================================================================!
! meds_column_dynamics -- the fast-timescale (dt_fast) integrator that couples the fast          !
! biophysics processes into one sub-daily operator-split sweep (design Part II). The canopy-air-  !
! space (CAS) twins are the shared coupling reservoir: each dt_fast the aerodynamics kernel sets   !
! the conductances; the leaf and GROUND surface fluxes feed the CAS; the soil WATER column         !
! (Richards) and soil THERMAL column (heat) are advanced; and the three CAS twins (enthalpy /       !
! specific humidity / CO2) are advanced IMPLICITLY in the atmospheric-exchange term, gated by u*.    !
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
! Still prescribed (next layers): absorbed radiation, precip. Canopy interception (leaf film), the         !
! stepper hook, and cross-demography persistence remain to wire.                                            !
!                                                                                          !
! WHOLE-COLUMN CONSERVATION -- verified by budg%whole_energy / budg%whole_water (Δ of ALL stores vs the     !
! true boundary fluxes; these CATCH cross-seam leaks the per-kernel budgets miss). Water-borne enthalpy    !
! is transported consistently: the CAS latent uses enthalpy_vapor(tl) (matching the CAS inverter + ground); !
! the soil sheds the transpiration water's liquid enthalpy via root_heat_sink; infiltration/drainage water  !
! carry internal_energy_liquid across the soil boundaries. INTER-LAYER advective heat: the hydrology kernel  !
! now EXPOSES the time-mean per-face Darcy flux (hflux%w_flux), and soil_energy_step_implicit can advect the liquid   !
! enthalpy on it -- an OPT-IN coupling (cfg%advect_soil_heat, default OFF). It reconciles the soil moisture   !
! <-> energy coupling that the standalone kernels never exercised (their unit tests forced moisture constant):!
! with the internal-energy soil store referenced at tsupercool_liq (~193 K), a per-step Δθ carries a large    !
! liquid enthalpy, so turning the advection ON shifts the coupled soil/CAS temperatures NOTABLY. It still      !
! conserves the whole-column energy to machine precision and keeps soil temps bounded, but the MAGNITUDE      !
! awaits validation against a reference (ED2) with real met forcing -- hence gated OFF by default. A          !
! SUPPLY-limited leaf is not re-solved for the extra sensible (small drought residual). The stepper hook +    !
! cross-demography persistence remain to wire.                                                                !
!==========================================================================================!
module meds_column_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, stefan, latent_heat_vap, rho_h2o, r_gas, pi, &
                                     tsupercool_liq, grav_head
   use meds_plant_hydraulics, only : rhizosphere_cond
   use meds_hydr_lib, only : soil_hydr_cond_from_theta
   use meds_config,           only : meds_config_t, hydraulics_config_t,                          &
                                     SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED,               &
                                     INTEG_SPLIT, INTEG_ARK
   use meds_hydr_lib,      only : build_hydro_table
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,        &
                                     energy_forcing_t, energy_opts_t, energy_flux_t,           &
                                     soil_column_t, soil_energy_column_t, chydro_forcing_t, chydro_flux_t, &
                                     leaf_energy_env_t, leaf_energy_flux_t, SOIL_BC_FREE_DRAIN, &
                                     snow_params_t, snow_env_t, snow_flux_t, snow_melt_t
   use meds_fast_time_derivs, only : column_state_t, column_frozen_t, surface_state_t,         &
                                     surface_frozen_t, surface_tend_t, surface_derivs, column_bflux_t
   use meds_ark_stepper,      only : ark2_column_step, adaptive_ark_march, bflux_zero, bflux_add
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_soil_energy,      only : soil_energy_step_implicit
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_vegetation_biophysics, only : veg_energy_diagnostic, veg_energy_step_implicit
   use meds_soil_water,       only : column_hydrology_flux
   use meds_ground_biophysics, only : snow_energy_step, snow_base_conductance,                  &
                                     snow_accumulate, snow_drain_meltwater, snow_cover_fraction, &
                                     ground_surface_fluxes
   use meds_plant_interface,  only : leaf_env_t, leaf_flux_t, leaf_gas_exchange,               &
                                     wood_env_t, wood_params_t, wood_flux_t, stem_maintenance_respiration, &
                                     root_env_t, root_params_t, root_flux_t, fine_root_maintenance_respiration, &
                                     hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     solve_plant_water, N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_soil_biogeochem,  only : heterotrophic_respiration_flux
   use meds_biogeochem_types, only : co2_opts_t
   use meds_therm_lib,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_liquid,  &
                                     sat_vapor_pressure, uext_to_temp, temp_to_uext
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok
   implicit none
   private

   public :: column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: alloc_column_cohort, column_fast_step, aero_bottom_to_top
   public :: apply_hydraulics_config

   !----- Static per-run column configuration (built once; constant across dt_fast steps). ----!
   type :: column_config_t
      type(aero_cfg_t)            :: aero            !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal    !< leaf/wood thermal params
      type(soil_params_t)         :: soil           !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal   !< soil thermal texture
      type(energy_opts_t)         :: energy         !< soil-thermal solver options
      type(soil_opts_t)           :: hydro          !< soil-water (Richards) solver options
      type(wood_params_t)         :: wood           !< stem-respiration parameters
      type(root_params_t)         :: root           !< fine-root-respiration parameters
      type(co2_opts_t)            :: co2            !< heterotrophic-respiration options
      type(hydro_params_t)        :: hydro_p        !< plant-hydraulics parameters (PV curves, vulnerability)
      type(hydro_opts_t)          :: hydro_o        !< plant-hydraulics solver options
      logical                     :: multilayer_roots  = .false.  !< opt-in soil->plant per-layer root coupling
      real(wp)                    :: specific_root_area = 20.0_wp  !< [m2/kgC] SRA (rhizosphere conductance)
      real(wp)                    :: fast_soil_carbon = 5.0_wp   !< [kgC/m2] decomposable soil-C pool (prescribed, MVP)
      real(wp)                    :: rhizo_cond       = 5.0e-4_wp !< [kg/s/MPa] soil->root conductance (prescribed, MVP)
      logical                     :: advect_soil_heat = .false.  !< opt-in: advect liquid enthalpy on the interior
                                                                 !< per-layer Darcy flux (moisture<->energy coupling)
      !----- P3 coupled-surface (Picard) solver knobs; only consulted under SCHEME_PICARD_COUPLED. !
      integer(ik) :: picard_max_iter = 20_ik        !< outer-iteration cap
      real(wp)    :: picard_tol_temp = 1.0e-3_wp     !< [K]     temperature convergence tolerance
      real(wp)    :: picard_tol_shv  = 1.0e-6_wp     !< [kg/kg] CAS specific-humidity convergence tolerance
      real(wp)    :: picard_relax    = 0.5_wp        !< [-]     under-relaxation of the next-pass seed. The CAS<->ground
                                                     !<         sensible coupling gives the fixed-point map a slope ~ -1
                                                     !<         (oscillatory); 0.5 makes the relaxed map a contraction.
      logical     :: picard_fixed_iter = .false.     !< run a uniform pass count (GPU warp-uniform; no early exit)
      integer(ik) :: leaf_energy_model  = 0_ik       !< LEAFEN_DIAGNOSTIC (0) | LEAFEN_PROGNOSTIC (1)
      integer(ik) :: wood_energy_model  = 0_ik       !< WOODEN_DIAGNOSTIC (0) | WOODEN_PROGNOSTIC (1)
      integer(ik) :: soil_water_coupling = 0_ik      !< SOILH2O_LAGGED (0) | SOILH2O_COUPLED (1)
      logical           :: snow_on = .false.         !< opt-in temporary-surface-water / snow store (P0, split path)
      type(snow_params_t) :: snow                    !< snow parameters (density, albedo, thresholds, conductivity)
   end type column_config_t

   !----- Leaf/wood thermal model + soil-water-in-loop selector codes (P3). ------------------!
   integer(ik), parameter, public :: LEAFEN_DIAGNOSTIC = 0_ik  !< steady-state leaf (tl = tcas + dtl)
   integer(ik), parameter, public :: LEAFEN_PROGNOSTIC = 1_ik  !< prognostic leaf_energy via veg_energy_step_implicit (P3e)
   integer(ik), parameter, public :: WOODEN_DIAGNOSTIC = 0_ik  !< steady-state wood (own balance; tw = tcas + dtw)
   integer(ik), parameter, public :: WOODEN_PROGNOSTIC = 1_ik  !< prognostic wood_energy via veg_energy_step_implicit
   integer(ik), parameter, public :: SOILH2O_LAGGED    = 0_ik  !< soil water/hydraulics frozen per sub-step
   integer(ik), parameter, public :: SOILH2O_COUPLED   = 1_ik  !< soil water re-solved inside the Picard loop (P3f)

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into cfg%pft)
      real(wp),    allocatable :: lai(:), wai(:), height(:), crown(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   !< [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              !< [kgC/plant],[kgC/plant],[m2] (hydraulics)
      real(wp),    allocatable :: vcmax25(:), rd25(:)                         !< [umol/m2/s] per-cohort (plastic) capacities
   end type column_cohort_t

   !----- Prescribed per-step forcing the higher layers (RT, met) supply; photosynthesis/    !
   !      respiration/NEE are now computed from the plant kernels (no longer prescribed).     !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: precip        = 0.0_wp   !< [kg/m2/s] ground-reaching rainfall (interception deferred)
      real(wp)              :: snowf         = 0.0_wp   !< [kg/m2/s] frozen precip (snowfall; drives snow accumulation)
      real(wp)              :: tair          = 288.0_wp !< [K] reference-level air temp (frozen/rain-on-snow precip enthalpy)
      real(wp)              :: par_per_w     = 2.1_wp   !< [umol photon / (W absorbed)] absorbed->PAR-photon factor
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW (VIS+NIR) / net LW per cohort (leaf ENERGY)
      real(wp), allocatable :: abs_par(:)               !< [W/m2] INCIDENT-equiv PAR (VIS) per cohort; the leaf
                                                        !< model re-applies leaf_absorptance internally (PHOTOSYNTHESIS)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] absorbed SW / net LW per cohort (WOOD energy)
   end type column_forcing_t

   !----- The per-patch conservation budgets (one place; the driver accumulates the closed resids).!
   !      The per-kernel budgets close BY CONSTRUCTION; whole_energy/whole_water are the CROSS-      !
   !      seam column totals (Δ all stores vs the true boundary fluxes) that actually catch leaks.   !
   type :: column_budget_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
      real(wp)       :: gpp_last = 0.0_wp, nee_last = 0.0_wp   !< [umol/m2/s] last-step diagnostics
      !----- P3 Picard diagnostics (reporting only; not conserved state). --------------------!
      integer(ik)    :: picard_iters       = 0_ik    !< worst outer-iteration count over the sub-steps
      integer(ik)    :: picard_nonconv     = 0_ik    !< number of sub-steps that hit picard_max_iter unconverged
      real(wp)       :: picard_worst_resid = 0.0_wp  !< [K] worst residual temperature at exit
   end type column_budget_t

contains

   !----- Flatten the shared [hydraulics] config into the plant hydro_params_t + rhizosphere      !
   !       conductance, and build the vulnerability lookup table from wood_kexp. The single seam    !
   !       between cfg%hydraulics (shared, TOML-driven) and the fast loop's hydro_params_t (plant),  !
   !       mirroring how the leaf seam flattens the PFT photosynthesis traits. -------------------!
   subroutine apply_hydraulics_config(hcfg, hydro_p, rhizo_cond)
      type(hydraulics_config_t), intent(in)    :: hcfg
      type(hydro_params_t),      intent(inout) :: hydro_p
      real(wp),                  intent(out)   :: rhizo_cond
      hydro_p%leaf_pi0       = hcfg%leaf_pi0       ; hydro_p%leaf_elastic_mod       = hcfg%leaf_elastic_mod
      hydro_p%leaf_apoplast_frac        = hcfg%leaf_apoplast_frac        ; hydro_p%leaf_water_sat = hcfg%leaf_water_sat
      hydro_p%wood_pi0       = hcfg%wood_pi0       ; hydro_p%wood_elastic_mod       = hcfg%wood_elastic_mod
      hydro_p%wood_apoplast_frac        = hcfg%wood_apoplast_frac        ; hydro_p%wood_water_sat = hcfg%wood_water_sat
      hydro_p%wood_psi50     = hcfg%wood_psi50     ; hydro_p%wood_kexp      = hcfg%wood_kexp
      hydro_p%k_plant_max    = hcfg%k_plant_max    ; hydro_p%wood_kmax      = hcfg%wood_kmax
      hydro_p%vessel_curl    = hcfg%vessel_curl
      call build_hydro_table(hydro_p%vuln_table, hydro_p%wood_kexp)
      rhizo_cond = hcfg%rhizo_cond
   end subroutine apply_hydraulics_config

   !----- Allocate a column_cohort_t (the per-patch cohort SoA the fast loop consumes). ------!
   subroutine alloc_column_cohort(coh, n)
      type(column_cohort_t), intent(out) :: coh
      integer(ik),           intent(in)  :: n
      coh%n = n
      allocate(coh%pft(n), coh%lai(n), coh%wai(n), coh%height(n), coh%crown(n),                &
               coh%leaf_width(n), coh%branch_diam(n), coh%leaf_area(n), coh%nplant(n),         &
               coh%dbh(n), coh%broot(n), coh%bleaf(n), coh%bsap(n), coh%sap_area(n),           &
               coh%vcmax25(n), coh%rd25(n))
      coh%pft = 1_ik
      coh%lai = 0.0_wp ; coh%wai = 0.0_wp ; coh%height = 0.0_wp ; coh%crown = 1.0_wp
      coh%leaf_width = 0.04_wp ; coh%branch_diam = 0.02_wp
      coh%leaf_area = 0.0_wp ; coh%nplant = 0.0_wp ; coh%dbh = 0.0_wp ; coh%broot = 0.0_wp
      coh%bleaf = 0.0_wp ; coh%bsap = 0.0_wp ; coh%sap_area = 0.0_wp
      coh%vcmax25 = 0.0_wp ; coh%rd25 = 0.0_wp
   end subroutine alloc_column_cohort

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
      type(leaf_env_t)       :: lenv
      type(leaf_flux_t)      :: lf
      type(wood_env_t)       :: wenv
      type(wood_flux_t)      :: wf
      type(root_env_t)       :: renv
      type(root_flux_t)      :: rf
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
      real(wp)    :: soil_psi_root, tground_in, t_ground_dia, t_bot_dia, k_theta, sink_tot
      real(wp)    :: tcas, qcas, press, rho, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: le_ref, dtl, tl, transp_i, gsw_ms, e_air, rho_mol, dh, drnet, transp_w
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
      real(wp)    :: gpp, ra_leaf, ra_stem, ra_root, rh, nee_biotic, soil_temp_root, theta_mean
      real(wp)    :: t_ground, t_bot, g_top, h_ground, le_ground, soil_evap, rain_temp
      real(wp)    :: gah, gaw, gac, wcap, ccap, can_dmol, src_enth, src_vap, src_frac
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

      !----- 1. Refresh the aerodynamics env from the current CAS state, then solve. ---------!
      bio%cas%can_temp = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      tcas = bio%cas%can_temp ; qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air
      t_ground = bio%soil_e%soil_temp(1) ; t_bot = bio%soil_e%soil_temp(nsl) ; rain_temp = tcas
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = t_ground
      !      canopy_aerodynamics expects cohorts ordered BOTTOM(1)->TOP(n) (its wind cascade walks !
      !      top->bottom), but the column buffer is gathered height-DESCENDING (index 1 = top). Feed !
      !      it the bottom->top order and scatter the per-cohort conductances back to gather order.   !
      call aero_bottom_to_top(ccfg%aero, aenv, ageom, n, coh, bio%leaf_temp, aero)

      !----- Root-weighted soil temperature + column-mean moisture (root / heterotrophic resp). !
      soil_temp_root = 0.0_wp ; theta_mean = 0.0_wp
      do k = 1_ik, nsl
         soil_temp_root = soil_temp_root + bio%soil_e%soil_temp(k) * ccfg%soil%root_frac(k)
         theta_mean     = theta_mean     + bio%soil_w%theta(k) * ccfg%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-ccfg%soil%soil_layer_z(nsl+1_ik), tiny_num)

      !----- 2. PRE-PASS (once per sub-step; ED2 freezes gs/hydraulics per DTLSM): LEAF gas      !
      !         exchange (GPP/gs/Rd), the FROZEN per-cohort leaf-energy coefficients h_coeff_f/    !
      !         g_tr_f, and stem/root maintenance respiration. These do NOT change across the      !
      !         Picard passes (they use the lagged, start-of-sub-step leaf_temp).                  !
      gpp = 0.0_wp ; ra_leaf = 0.0_wp ; ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      do i = 1_ik, n
         rho_mol       = press / (r_gas * bio%leaf_temp(i))                     ! [mol/m3] molar air density
         e_air         = qcas * press / (0.622_wp + 0.378_wp * qcas)            ! [Pa] canopy-air vapour pressure
         lenv%par      = forc%abs_par(i) / max(coh%lai(i), 0.1_wp) * forc%par_per_w   ! absorbed PAR (VIS), not total SW
         lenv%leaf_temp = bio%leaf_temp(i)
         t_emit(i)      = bio%leaf_temp(i)     ! start-of-sub-step leaf temp = the RT LW emission base (P3c)
         wood_emit(i)   = bio%wood_temp(i)     ! start-of-sub-step wood temp = prognostic-wood seed (Picard-correct)
         lenv%vpd      = max(sat_vapor_pressure(bio%leaf_temp(i)) - e_air, 0.0_wp)
         lenv%ca       = bio%cas%can_co2
         lenv%pressure = press
         lenv%psi_leaf = bio%psi(NODE_LEAF, i)                                  ! lagged plant water status (hydraulics)
         lenv%gb       = aero%leaf_gbw(i) * rho_mol                             ! m/s -> mol H2O/m2/s
         call leaf_gas_exchange(lenv, cfg, coh%pft(i), lf, vcmax25=coh%vcmax25(i), rd25=coh%rd25(i))
         gsw_ms        = lf%gs / max(rho_mol, tiny_num)                         ! mol/m2/s -> m/s
         gpp           = gpp     + lf%A_gross * coh%leaf_area(i) * coh%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = lf%A_gross * coh%leaf_area(i)          ! [umol/plant/s] per-plant gross
         ra_leaf       = ra_leaf + lf%rd      * coh%leaf_area(i) * coh%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = lf%rd * coh%leaf_area(i)   ! [umol/plant/s] per-plant leaf Rd
         !----- frozen leaf-energy coefficients (sensible + transpiration series conductance). --!
         h_coeff_f(i) = ccfg%veg_thermal%effarea_heat * coh%lai(i) * aero%leaf_gbh(i) * rho * cp_air
         g_tr_f(i)    = 0.0_wp
         if (aero%leaf_gbw(i) + gsw_ms > tiny_num) then
            g_tr_f(i) = ccfg%veg_thermal%effarea_transp * coh%lai(i)                            &
                        * aero%leaf_gbw(i) * gsw_ms / (aero%leaf_gbw(i) + gsw_ms)
         end if
         !----- Autotrophic maintenance respiration: stem + fine root (per plant -> per m2). ---!
         wenv%wood_temp = bio%wood_temp(i) ; wenv%dbh = coh%dbh(i) ; wenv%height = coh%height(i)
         wenv%wai = coh%wai(i) ; wenv%nplant = coh%nplant(i)
         call stem_maintenance_respiration(wenv, ccfg%wood, wf)
         renv%soil_temp = soil_temp_root ; renv%broot = coh%broot(i)
         call fine_root_maintenance_respiration(renv, ccfg%root, rf)
         ra_stem = ra_stem + wf%stem_resp * coh%nplant(i)
         ra_root = ra_root + rf%root_resp * coh%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = wf%stem_resp   ! [umol/plant/s] already per-plant
         if (present(root_resp_coh)) root_resp_coh(i) = rf%root_resp   ! [umol/plant/s] already per-plant
      end do

      !----- NEE = autotrophic (leaf Rd + stem + root) + heterotrophic Rh - GPP (thermally      !
      !      passive: frozen across passes; the CO2 twin is solved once after convergence). ----!
      rh = heterotrophic_respiration_flux(ccfg%fast_soil_carbon, soil_temp_root, theta_mean,   &
                                          ccfg%soil%theta_res(1), ccfg%soil%theta_sat(1), ccfg%co2)
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budg%gpp_last = gpp ; budg%nee_last = nee_biotic

      !----- CAS capacities + atm-exchange conductances (frozen across passes, §3.5). ---------!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      wcap     = rho      * bio%cas%can_depth
      ccap     = can_dmol * bio%cas%can_depth
      gah      = rho      * aero%ustar * aero%temp1
      gaw      = rho      * aero%ustar * aero%temp2
      gac      = can_dmol * aero%ustar * aero%temp2
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
      !     water/src_frac/hydraulics are frozen after pass 1 (SOILH2O_LAGGED) unless coupled.      !
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
            lw_slope = 4.0_wp * ccfg%veg_thermal%leaf_emiss * stefan * te**3 * coh%lai(i)
            le_slope = latent_heat_vap * rho * g_tr_f(i) * dqdt
            le_ref   = latent_heat_vap * rho * g_tr_f(i) * (qsat_c - qcas)
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
               h_coeff_w  = pi * coh%wai(i) * aero%wood_gbh(i) * rho * cp_air
               te_w       = tcas
               lw_slope_w = 4.0_wp * ccfg%veg_thermal%leaf_emiss * stefan * te_w**3 * coh%wai(i)
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
         soil_psi_root = 0.0_wp
         do k = 1_ik, nsl
            soil_psi_root = soil_psi_root + hflux%psi_soil(k) * ccfg%soil%root_frac(k)
         end do
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
      integer(ik) :: n, nsl, k, isub, nsub, nsteps, nrej

      n = coh%n ; nsl = ccfg%soil%n_active

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
         call adaptive_ark_march(y, fro, n, nsl, dt_fast, cfg%ark_rtol, dt0, y_out, nsteps, nrej,   &
                                 niter=cfg%ark_niter, relax=cfg%ark_relax, acc=acc)
      else
         nsub = max(1_ik, cfg%ark_fixed_substep) ; nrej = 0_ik ; ycur = y ; call bflux_zero(acc)
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, fro, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr,          &
                                  niter=cfg%ark_niter, relax=cfg%ark_relax, bf=bfsub)
            call bflux_add(acc, bfsub)
            ycur = ytmp
         end do
         y_out = ycur ; nsteps = nsub
      end if

      !----- SOIL WATER is operator-split out: the ESDIRK stages passed theta through unchanged (=theta^n); !
      !      commit the AUTHORITATIVE end-of-step theta from the scratch column_hydrology_flux HERE, once,  !
      !      so a single consistent theta feeds the state commit, the soil_temp read-off, and BOTH the      !
      !      soil_water and whole_water storage terms (w_soil1 below). ------------------------------------!
      y_out%theta(1:nsl) = fro%theta1(1:nsl)

      !----- unpack into bio + re-derive the diagnostic soil temperatures + leaf temperatures. -----!
      bio%cas%can_enthalpy = y_out%cas_enthalpy ; bio%cas%can_shv = y_out%cas_shv ; bio%cas%can_co2 = y_out%cas_co2
      bio%cas%can_temp = cas_temp_of_enthalpy(y_out%cas_enthalpy, y_out%cas_shv)
      bio%soil_e%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      bio%soil_w%theta(1:nsl)       = y_out%theta(1:nsl)
      bio%psi(:, 1:n)               = y_out%psi(:, 1:n)
      !----- persist the scratch hydrology's ponding/aquifer/water-table (lagged operator split). ------!
      bio%soil_w%w_surface = fro%w_surface1
      bio%soil_w%w_aquifer = fro%w_aquifer1
      bio%soil_w%z_wt      = fro%z_wt1
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
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                             1.0_wp, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(budg%cas_water,  wcap*shv0,  wcap*shv1,  acc%cas_vap_in,  acc%cas_vap_out,  &
                             1.0_wp, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(budg%cas_co2,    ccap*co20,  ccap*co21,  acc%cas_co2_in,  acc%cas_co2_out,  &
                             1.0_wp, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      call budget_accumulate(budg%soil_energy, e_soil0, e_soil1, acc%soil_enth_in, acc%soil_enth_out,    &
                             1.0_wp, abs(e_soil1) + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp)
      !----- SOIL WATER (fully frozen now): storage theta^n -> theta1 (w_soil0 -> w_soil1, both from the    !
      !      scratch solve), inflow q_top*rho, outflow drainage + realized uptake -- all from the frozen    !
      !      hflux, which closed its OWN mass budget to machine precision inside column_hydrology_flux. -----!
      call budget_accumulate(budg%soil_water,  w_soil0, w_soil1,                                        &
                             fro%q_top*rho_h2o*dt_fast, (fro%drainage + fro%uptake)*dt_fast,            &
                             1.0_wp, max(w_soil1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
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
      call budget_accumulate(budg%whole_energy, e_soil0 + wcap*enth0, e_soil1 + wcap*enth1, acc%whole_enth_in, &
                             acc%whole_enth_out, 1.0_wp, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp)

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_ark

   !----- Build the frozen pre-pass (leaf gas exchange / respiration / CAS caps / aero) + the frozen  !
   !      hydrology BCs into a column_frozen_t, and pack the prognostic state into a column_state_t.   !
   !      The pre-pass loop is a VERBATIM copy of column_fast_step's :237-314 (writing struct fields   !
   !      instead of locals), so gpp_coh/resp are bit-identical to the split; the split's inline       !
   !      pre-pass is untouched (the golden anchor stays byte-for-byte).                               !
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

      type(leaf_env_t)       :: lenv ; type(leaf_flux_t) :: lf
      type(wood_env_t)       :: wenv ; type(wood_flux_t) :: wf
      type(root_env_t)       :: renv ; type(root_flux_t) :: rf
      type(chydro_forcing_t) :: hforc ; type(chydro_flux_t) :: hflux
      type(soil_column_t)    :: soil_w_scratch
      type(surface_state_t)  :: ys ; type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, rho_mol, e_air, gsw_ms, can_dmol
      real(wp) :: gpp, ra_leaf, ra_stem, ra_root, rh, nee_biotic, soil_temp_root, theta_mean
      integer(ik) :: i, k

      allocate(fro%surf%h_coeff_f(n), fro%surf%g_tr_f(n), fro%surf%abs_sw(n), fro%surf%abs_lw(n), fro%surf%lai(n))
      allocate(fro%surf%h_coeff_w(n), fro%surf%abs_sw_wood(n), fro%surf%abs_lw_wood(n), fro%surf%wai(n))
      allocate(fro%psi_e(nsl), fro%nplant(n), fro%bleaf(n), fro%bsap(n), fro%broot(n),            &
               fro%sap_area(n), fro%height(n), fro%leaf_area(n))
      allocate(y%psi(N_HYDRO, n))

      !----- aerodynamics from the current CAS state. -------------------------------------------!
      tcas = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air ; t_ground = bio%soil_e%soil_temp(1)
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = t_ground
      call aero_bottom_to_top(ccfg%aero, aenv, ageom, n, coh, bio%leaf_temp, aero)

      soil_temp_root = 0.0_wp ; theta_mean = 0.0_wp
      do k = 1_ik, nsl
         soil_temp_root = soil_temp_root + bio%soil_e%soil_temp(k) * ccfg%soil%root_frac(k)
         theta_mean     = theta_mean     + bio%soil_w%theta(k) * ccfg%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-ccfg%soil%soil_layer_z(nsl+1_ik), tiny_num)

      !----- leaf gas exchange + frozen leaf-energy coefficients + maintenance respiration. -----!
      gpp = 0.0_wp ; ra_leaf = 0.0_wp ; ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      do i = 1_ik, n
         rho_mol        = press / (r_gas * bio%leaf_temp(i))
         e_air          = qcas * press / (0.622_wp + 0.378_wp * qcas)
         lenv%par       = forc%abs_par(i) / max(coh%lai(i), 0.1_wp) * forc%par_per_w
         lenv%leaf_temp = bio%leaf_temp(i)
         lenv%vpd       = max(sat_vapor_pressure(bio%leaf_temp(i)) - e_air, 0.0_wp)
         lenv%ca        = bio%cas%can_co2 ; lenv%pressure = press
         lenv%psi_leaf  = bio%psi(NODE_LEAF, i)
         lenv%gb        = aero%leaf_gbw(i) * rho_mol
         call leaf_gas_exchange(lenv, cfg, coh%pft(i), lf, vcmax25=coh%vcmax25(i), rd25=coh%rd25(i))
         gsw_ms         = lf%gs / max(rho_mol, tiny_num)
         gpp            = gpp     + lf%A_gross * coh%leaf_area(i) * coh%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = lf%A_gross * coh%leaf_area(i)
         ra_leaf        = ra_leaf + lf%rd      * coh%leaf_area(i) * coh%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = lf%rd * coh%leaf_area(i)
         fro%surf%h_coeff_f(i) = ccfg%veg_thermal%effarea_heat * coh%lai(i) * aero%leaf_gbh(i) * rho * cp_air
         fro%surf%g_tr_f(i)    = 0.0_wp
         if (aero%leaf_gbw(i) + gsw_ms > tiny_num) then
            fro%surf%g_tr_f(i) = ccfg%veg_thermal%effarea_transp * coh%lai(i)                     &
                                 * aero%leaf_gbw(i) * gsw_ms / (aero%leaf_gbw(i) + gsw_ms)
         end if
         wenv%wood_temp = bio%wood_temp(i) ; wenv%dbh = coh%dbh(i) ; wenv%height = coh%height(i)
         wenv%wai = coh%wai(i) ; wenv%nplant = coh%nplant(i)
         call stem_maintenance_respiration(wenv, ccfg%wood, wf)
         renv%soil_temp = soil_temp_root ; renv%broot = coh%broot(i)
         call fine_root_maintenance_respiration(renv, ccfg%root, rf)
         ra_stem = ra_stem + wf%stem_resp * coh%nplant(i)
         ra_root = ra_root + rf%root_resp * coh%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = wf%stem_resp
         if (present(root_resp_coh)) root_resp_coh(i) = rf%root_resp
         !----- per-cohort geometry + radiation the ARK frozen inputs need. --------------------!
         fro%surf%lai(i)    = coh%lai(i)
         fro%surf%abs_sw(i) = forc%abs_sw(i) ; fro%surf%abs_lw(i) = forc%abs_lw(i)
         !----- WOOD frozen inputs: real diagnostic values, or ZERO when wood is not diagnostic (so   !
         !      surface_derivs' wood branch is a no-op; prognostic wood is operator-split in P2). -----!
         if (ccfg%wood_energy_model == WOODEN_DIAGNOSTIC) then
            fro%surf%wai(i)         = coh%wai(i)
            fro%surf%h_coeff_w(i)   = pi * coh%wai(i) * aero%wood_gbh(i) * rho * cp_air
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
      rh = heterotrophic_respiration_flux(ccfg%fast_soil_carbon, soil_temp_root, theta_mean,      &
                                          ccfg%soil%theta_res(1), ccfg%soil%theta_sat(1), ccfg%co2)
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budg%gpp_last = gpp ; budg%nee_last = nee_biotic

      !----- CAS capacities + atm conductances + the rest of the frozen surface inputs. ---------!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      fro%surf%leaf_emiss = ccfg%veg_thermal%leaf_emiss
      fro%surf%wcap = rho * bio%cas%can_depth ; fro%surf%ccap = can_dmol * bio%cas%can_depth
      fro%surf%gah  = rho * aero%ustar * aero%temp1
      fro%surf%gaw  = rho * aero%ustar * aero%temp2
      fro%surf%gac  = can_dmol * aero%ustar * aero%temp2
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
      fro%surf%soil_evap = hflux%soil_evap
      fro%q_top          = (hflux%infiltration - hflux%soil_evap) / rho_h2o
      fro%soil_psi_root  = 0.0_wp
      do k = 1_ik, nsl
         fro%soil_psi_root = fro%soil_psi_root + hflux%psi_soil(k) * ccfg%soil%root_frac(k)
      end do
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

   !----- Track a kernel's own closed-budget residual into a budget_t (worst + fail count). ---!
   pure subroutine track_resid(b, resid, scale, rtol, atol)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: resid, scale, rtol, atol
      b%n_check = b%n_check + 1_ik
      b%resid   = resid
      b%worst   = max(b%worst, abs(resid))
      if (.not. closure_ok(resid, scale, rtol, atol)) b%n_fail = b%n_fail + 1_ik
   end subroutine track_resid

end module meds_column_dynamics
