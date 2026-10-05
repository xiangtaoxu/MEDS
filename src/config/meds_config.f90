! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_config -- immutable run configuration, threaded read-only through the engine.       !
!                                                                                          !
! Holds the time-stepping mode, the fusion/fission tunables (the diameter & size-           !
! distribution analogues of ED2's LAI/light tolerances), and the PFT trait table.          !
! `derive_config` fills derived quantities; `validate_config` error-stops on inconsistent    !
! settings (e.g. recruits born below the termination size, which would churn forever).      !
!==========================================================================================!
module meds_config
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : yr_day, yr_sec, day_sec
   use meds_pft_params, only : pft_table_t, PATH_C3, PATH_C4, derive_pft_rates, derive_leaf_params, HYD_UNSET
   use meds_allometry,  only : set_allometry
   use meds_time,       only : meds_time_t, time_lt, time_valid, time_to_string,                &
                               whole_years_between
   use meds_temp_response, only : TRESP_ARRHENIUS, TRESP_PEAKED
   use meds_leaf_opts,     only : SM_LEUNING, SM_MEDLYN, SM_KATUL, COLIM_MIN, COLIM_QUADRATIC
   use meds_water_retention, only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, pv_psi_tlp
   use meds_column_params, only : n_soil_layer_max, soil_params_t, build_soil_hydr_params
   use meds_forcing_config, only : forcing_config_t, LW_SYNTHESIZE, METAVG_INSTANT, METAVG_CENTER,   &
                                   METAVG_END, MET_BACKEND_ED_ERA5LAND, SWPART_PASSTHROUGH,        &
                                   WIND_EXPOSURE_OPEN_TERRAIN, ARCHIVE_DT_SEC
   use meds_output_config,  only : output_config_t
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   use meds_biophysics_opts, only : ENERGY_BC_DIRICHLET
   use meds_biogeochem_opts, only : decomp_opts_t
   use meds_region_opts,   only : region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   implicit none
   private

   public :: meds_config_t, allometry_config_t, hydraulics_config_t, soil_column_config_t
   public :: derive_config, derive_parameters
   public :: MAX_RECYCLE_YEARS
   public :: validate_config, growth_window_steps
   public :: pft_leaf_psi_tlp, pft_stomata_psi_onset
   public :: forcing_config_t, output_config_t
   public :: decomp_opts_t
   public :: region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   public :: BK_SERIAL, BK_MULTICORE, BK_GPU
   public :: DIST_PRIMARY, DIST_TREEFALL
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART
   public :: INTEG_ARK, INTEG_RK45
   public :: LWP_CONTROL_LINEAR_DECLINE
   public :: HYD_CONDUCTANCE_WHOLE_PLANT, HYD_CONDUCTANCE_SEGMENT
   public :: CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT, CTRL_I, CTRL_PI

   !----- Time-step modes. ----------------------------------------------------------------!
   !----- Parallel backend labels (the actual backend is chosen at COMPILE time via the    !
   !      compiler's do-concurrent target; this is for reporting/reproducibility only).    !
   integer(ik), parameter :: BK_SERIAL    = 0_ik
   integer(ik), parameter :: BK_MULTICORE = 1_ik
   integer(ik), parameter :: BK_GPU       = 2_ik
   !----- Patch disturbance / land-use classes. -------------------------------------------!
   !----- Upper bound on the declared forcing recycle window, in whole calendar years. Purely a  !
   !      search bound for the whole-year test (no real met record spans centuries). -------------!
   integer(ik), parameter :: MAX_RECYCLE_YEARS = 200_ik

   integer(ik), parameter :: DIST_PRIMARY  = 1_ik   !< undisturbed / primary stand
   integer(ik), parameter :: DIST_TREEFALL = 2_ik   !< treefall-gap (age-0) patch
   !----- Initialization modes (selected by [init].init_mode). ----------------------------!
   integer(ik), parameter :: INIT_BARE    = 0_ik    !< near-bare ground
   integer(ik), parameter :: INIT_CENSUS  = 1_ik    !< from a cohort census CSV (init_census_file)
   integer(ik), parameter :: INIT_RESTART = 2_ik    !< restart from a state .nc file (init_restart_file)
   !----- The leaf-model selectors (SM_*, COLIM_*) are owned by meds_leaf_opts and the        !
   !      temperature-response selectors (TRESP_*) by meds_temp_response. Imported above for    !
   !      LOADING and VALIDATION only -- not re-exported, so a caller that switches on one       !
   !      imports it from the module that defines it (decision #8).                              !


   !----- Fast-loop TIME integrator ([fast].time_integrator). TWO schemes, and there is deliberately  !
   !      no third: docs/science/numerical_scheme.md section 3 records why there is no operator-split   !
   !      path, and there is no separate coupling-sweep selector because nothing is swept.              !
   !                                                                                          !
   !      NAMING, stated once so it stops propagating: INTEG_ARK is NOT an IMEX method. The biotic     !
   !      CO2 source is folded implicit, so the explicit tableau is empty (f_E == 0) and the scheme    !
   !      is a 2-solve ESDIRK2 with gamma = 1 - 1/sqrt(2) (the ARS(2,2,2) value). The config string    !
   !      stays "ark" for compatibility. ------------------------------------------------------------!
   !----- STOMATAL CONTROL AT LOW LEAF WATER POTENTIAL ([leaf_physiology].low_water_potential_control, !
   !      #332). Without one, a plant transpires with an empty internal store: beta_stomata scales g1 !
   !      only, so conductance falls to the residual g0 and never reaches zero (issue #95). This     !
   !      multiplies the conductance the stomatal model calculates, g0 included, by a factor that    !
   !      falls LINEARLY from 1 at the turgor-loss point psi_tlp to 0 at 2*psi_tlp. It replaced a    !
   !      hard shutdown at 2*psi_tlp, whose step no calibration could see past. -------------------!
   integer(ik), parameter :: LWP_CONTROL_LINEAR_DECLINE = 1_ik  !< "linear_decline" (the only option)
   integer(ik), parameter :: HYD_CONDUCTANCE_WHOLE_PLANT = 1_ik  !< [hydraulics] conductance = "whole_plant"
   integer(ik), parameter :: HYD_CONDUCTANCE_SEGMENT     = 2_ik  !< [hydraulics] conductance = "segment"
   !----- RESERVED: a "dynamic vapour pressure" control -- the substomatal air held at the Kelvin     !
   !      humidity exp(psi/(rho_w*Rv*T)) rather than saturated, so the transpiration gradient shrinks  !
   !      with psi and REVERSES into foliar water uptake once e_i < e_a. It was implemented, measured  !
   !      and REMOVED (see issue #96 for the pro/con and docs/science/leaf_gas_exchange.md): scaling   !
   !      the existing latent flux by a humidity factor is not sound, because that flux is linearised  !
   !      about T_cas. It needs the flux AND its temperature derivative rebuilt from                   !
   !      rh_leaf*qsat(T_leaf) inside the leaf energy balance. Deferred, and interesting mainly as a   !
   !      route to representing foliar uptake rather than as a stress arrestor. ---------------------!
   integer(ik), parameter :: INTEG_ARK   = 2_ik  !< ESDIRK2 coupled implicit column (DEFAULT)
   integer(ik), parameter :: INTEG_RK45   = 3_ik  !< adaptive Cash-Karp RK45, the ACCURACY BASELINE
   !----- Above this dt_fast, an rk45 run is warned that it lacks the ARK-only transpiration       !
   !      corrector (#160). Not a hard limit -- rk45 at a fine step is the intended use.            !
   real(wp),    parameter :: RK45_UNCORRECTED_DT_WARN = 300.0_wp   !< [s]

   !----- Fast-loop ERROR-CONTROL selectors (MEDS_NUMERICS_SCOPING.md goal (a); consumed by            !
   !      meds_fast_control). Strictness LEVEL ([fast].error_level): L0 fixed / L1 adaptive (default) / !
   !      L2 strict (adaptive + hard failure when a floor step can't meet tolerance). CONTROLLER        !
   !      ([fast].step_controller): I = integral (default, legacy) / PI = Gustafsson proportional-      !
   !      integral (damps step-size hunting). ---------------------------------------------------------!
   integer(ik), parameter :: CTRL_L0_FIXED    = 0_ik
   integer(ik), parameter :: CTRL_L1_ADAPTIVE = 1_ik
   integer(ik), parameter :: CTRL_L2_STRICT   = 2_ik
   integer(ik), parameter :: CTRL_I  = 1_ik
   integer(ik), parameter :: CTRL_PI = 2_ik

   !----- NO hard-coded defaults: every field is set by the config reader (presence-mapped) or  !
   !       derived (derive_config / derive_pft_rates). DERIVED fields are noted.  --------------!
   !----- Global (PFT-independent, ED2 iallom==3) allometry coefficients. Held on cfg so the run  !
   !       config is the complete record; installed into meds_allometry's protected module state    !
   !       by derive_parameters (mirroring how every other derived quantity is computed at load).    !
   type :: allometry_config_t
      real(wp) :: b1Ht = 0.0_wp, b2Ht = 0.0_wp        !< height <-> diameter intercept / slope
      real(wp) :: agb_c1 = 0.0_wp, agb_c2 = 0.0_wp    !< AGB scale / exponent (Chave-2014)
      real(wp) :: ca_b1 = 0.0_wp, ca_b2 = 0.0_wp      !< crown-area scale / exponent
      real(wp) :: lai_b1 = 0.0_wp, lai_b2 = 0.0_wp    !< per-stem leaf-area scale / exponent
      real(wp) :: light_ext = 0.0_wp                  !< Beer-Lambert extinction through overtopping LAI
   end type allometry_config_t

   !----- SOIL COLUMN: the site's geometry, texture and thermal properties. Consumed once per run  !
   !       by the fast-context builder, which turns them into soil_params_t + soil_thermal_params_t. !
   !                                                                                          !
   !       These were HARD-CODED literals in build_fast_context, with a comment promising this      !
   !       block. That mattered for three reasons: the source is supposed to define only true        !
   !       constants; the column could not be given a second site's texture without a recompile;     !
   !       and `depth` is the value docs/dev_plans note is shallower than the ~2.5 m annual thermal  !
   !       damping depth, so the fix for that was unreachable from a config file.                     !
   !                                                                                          !
   !       Every key is OPTIONAL and defaults to the literal it replaced, so a config without the     !
   !       block is bit-identical to the old behaviour. NOTE the distinction from [soil], which is    !
   !       the Richards SOLVER's options (tolerances, bottom BC, substepping): this block is the       !
   !       PHYSICAL column those options are solved over.                                              !
   type :: soil_column_config_t
      !----- Vertical grid. `grid_growth` = 0 gives a uniform grid; > 0 thickens layers with depth. -!
      integer(ik) :: n_layer     = 10_ik      !< [-]  active soil layers (<= n_soil_layer_max)
      real(wp)    :: depth       = 2.0_wp     !< [m]  total column depth (positive; z is negative down)
      real(wp)    :: grid_growth = 3.0_wp     !< [-]  layer-thickness growth factor with depth
      !----- Hydraulic texture (uniform over the column in this MVP; per-layer texture is future). --!
      integer(ik) :: retention   = SOIL_RETENTION_VG  !< van_genuchten | campbell
      real(wp)    :: theta_sat   = 0.43_wp    !< [m3/m3] porosity
      real(wp)    :: theta_res   = 0.078_wp   !< [m3/m3] residual water content
      real(wp)    :: ksat        = 2.89e-6_wp !< [m/s]   saturated hydraulic conductivity
      real(wp)    :: curve_par_a = 3.6_wp     !< [1/m] van Genuchten alpha, OR [m] Campbell psi_sat (< 0)
      real(wp)    :: curve_par_n = 1.56_wp    !< [-]   van Genuchten n (>1), OR [-] Campbell b
      real(wp)    :: psi_fc      = -3.37_wp   !< [m]   field-capacity matric head (derives theta_fc)
      !----- Thermal texture (uniform over the column). ---------------------------------------------!
      real(wp)    :: solid_conductivity = 3.0_wp    !< [W/m/K]  mineral-solid conductivity
      real(wp)    :: dry_conductivity   = 0.15_wp   !< [W/m/K]  dry-matrix conductivity
      real(wp)    :: dry_heat_capacity  = 2.0e6_wp  !< [J/m3/K] dry-matrix volumetric heat capacity
   end type soil_column_config_t

   !----- Plant-hydraulics parameters (PFT-uniform MVP). Consumed only by the opt-in fast loop;    !
   !       flattened into the plant hydro_params_t by the fast-context builder. Defaults are the     !
   !       former hardcoded fast-loop placeholders; a [hydraulics] TOML block overrides any field.    !
   type :: hydraulics_config_t
      !----- Pressure-volume (Bartlett/Tyree-Hammel), per tissue. --------------------------!
      real(wp) :: leaf_pi0 = -1.5_wp, leaf_elastic_mod = 12.0_wp, leaf_apoplast_frac = 0.30_wp   !< [MPa],[MPa],[-]
      real(wp) :: leaf_water_sat = 2.0_wp                                     !< [kg H2O/kgC] at saturation
      real(wp) :: wood_pi0 = -1.0_wp, wood_elastic_mod =  8.0_wp, wood_apoplast_frac = 0.20_wp
      real(wp) :: wood_water_sat = 1.0_wp
      !----- Xylem vulnerability + conductance. -------------------------------------------!
      real(wp) :: wood_psi50 = -2.0_wp   !< [MPa,<0] potential at 50% loss
      real(wp) :: wood_kexp  =  2.0_wp   !< [-]  vulnerability shape (a)
      real(wp) :: k_plant_max = 6.0e-4_wp !< [kg/s/MPa/m2_leaf] whole-plant conductance
      !----- How the plant's maximum internal conductance is set: `whole_plant` (default) takes      !
      !      k_plant_max per unit leaf area; `segment` takes the sapwood's specific conductivity over    !
      !      the path, wood_kmax * sapwood area / (height * vessel_curl). -----------------------------!
      integer(ik) :: conductance = HYD_CONDUCTANCE_WHOLE_PLANT
      real(wp) :: wood_kmax   = 8.0_wp    !< [kg/m/s/MPa] sapwood specific conductivity (conductance = segment)
      real(wp) :: vessel_curl = 1.5_wp    !< [-] tortuosity / path-length factor (conductance = segment)
      !----- The root profile, a plant trait: ED2's root_beta^(depth/root_depth), normalized over the   !
      !       soil column; layers below root_depth hold no roots. It sets the per-layer root boundary  !
      !       (uptake and rhizosphere conductance), the root-zone temperature of root respiration, and !
      !       the root-weighted soil state. The default beta = exp(-4) with a 2 m rooting depth is the  !
      !       exponential profile exp(-2 m^-1 * depth) that [soil_column].root_beta = 2 used to give.   !
      !       Hydraulic redistribution stays off -- per-layer efflux is floored at zero in both the     !
      !       plant solver and the soil sink (docs/ROADMAP.md section 7). -----------------------------!
      real(wp) :: root_beta          = 0.018315638888734179_wp  !< [-] ED2 root-profile decay (0,1): exp(-4)
      real(wp) :: root_depth         = 2.0_wp   !< [m]      maximum rooting depth
      real(wp) :: specific_root_area = 20.0_wp  !< [m2/kgC] fine-root absorbing area per unit root carbon
      !----- OPT-IN: couple the plant hydraulics to the per-layer soil column (feed per-layer psi_soil  !
      !       + rhizosphere conductance into the multi-layer root boundary) instead of a single root-    !
   end type hydraulics_config_t

   type :: meds_config_t
      !----- Time stepping (run bounded by start/end calendar dates). ----------------------!
      real(wp)          :: dt_slow               !< [s] slow-process timestep (user resolution; default 1 d)
      real(wp)          :: dt_years              !< DERIVED = dt_slow / yr_sec
      type(meds_time_t) :: start_time, end_time
      logical     :: demography_on               !< if .false. structure is frozen
      !----- Master SLOW-TIER freeze ([run].slow_on, DEFAULTED true): .false. skips the WHOLE slow   !
      !      tier (vegetation_dynamics -- growth/mortality/phenology/traits/demography -- and the      !
      !      future biogeochemistry step) every step, holding cohort/patch/soil-carbon state static    !
      !      while the fast loop still runs. Broader than demography_on, which only freezes the        !
      !      structural fuse/fiss/disturbance triggers within an otherwise-active slow loop.           !
      logical     :: slow_on = .true.
      !----- The SLOW-tier conservation ledger ([run].slow_ledger_on, DEFAULTED true). It snapshots  !
      !      the site store either side of every slow-step operator and reports what did not close   !
      !      (plan §10.2). ON by default for the same reason the fast loop's whole-column ledgers    !
      !      are: a conservation check nobody runs is a conservation check nobody has. It writes no  !
      !      state and changes no answer -- only an end-of-run table -- so the only cost of leaving  !
      !      it on is a handful of reductions per simulated day.  ---------------------------------!
      logical     :: slow_ledger_on = .true.
      !----- Host THREADS ([run].n_threads, DEFAULTED 1). A site run puts them on the fast loop's    !
      !      PATCH axis (§7 C2); a region run puts them on its POLYGONS, each polygon's patch loop   !
      !      then running on one (#183 R3). Patches, and polygons, are independent within a step, so !
      !      this is the one lever that costs no accuracy -- PROVIDED the answer does not move with  !
      !      the thread count, which is why the site-level reductions are staged per (sub-step,       !
      !      patch) and folded back in patch order (§7 C3). Default 1 so no existing result moves without opt-in, and so a build that       !
      !      happens to carry OpenMP flags (NVHPC MEDS_GPU=multicore puts -mp PUBLIC on               !
      !      meds_demography, which its dependents inherit) stays serial until asked. Has effect only !
      !      in an OpenMP build, which is the default; with -DMEDS_OPENMP=OFF the directives are      !
      !      comments and this is ignored.                                                            !
      integer(ik) :: n_threads = 1_ik
      !----- [run].mode: one site, or every selected cell of a box as its own polygon, with the box  !
      !      and the selection rules in [region] (MEDS_POLYGON_RUNTIME_PLAN.md §9).                  !
      integer(ik)         :: run_mode = RUN_MODE_SITE
      type(region_opts_t) :: region
      !----- Fast (sub-daily) biophysics loop. --------------------------------------------!
      logical     :: fast_biophysics_on          !< master gate for the fast biophysics loop
      real(wp)    :: dt_fast                      !< [s] fast biophysics timestep (nested within dt_slow)
      integer(ik) :: n_fast_per_slow              !< DERIVED = max(1, nint(dt_slow / dt_fast))
      real(wp)    :: snow_init_swe        = 0.0_wp     !< [kg/m2] initial snow water-equivalent seeded at run start
      real(wp)    :: snow_init_temp       = 270.0_wp   !< [K] initial snow temperature (for the seeded pack)
      logical     :: canopy_water_on      = .false.    !< opt-in canopy interception film + film-evap/dew
      !----- Fast-loop TIME integrator selector + ARK knobs ([fast], DEFAULTED reads). ----------------!
      !      every existing config + the golden anchor byte-identical). --------------------------------!
      integer(ik) :: time_integrator      = INTEG_ARK !< INTEG_ARK (default) | INTEG_RK45
      !----- The stomatal control at low leaf water potential (LWP_CONTROL_*, #332): the stomatal     !
      !      conductance times a factor falling linearly from 1 at psi_tlp to 0 at 2*psi_tlp, on the    !
      !      previous day's daily-max leaf potential. The thermodynamic route -- the substomatal air   !
      !      at its Kelvin humidity, so the driving gradient shrinks with psi with no threshold at all  !
      !      -- was built and removed (issue #96); once #96 is in, this control matters less. -------!
      integer(ik) :: low_water_potential_control = LWP_CONTROL_LINEAR_DECLINE !< [leaf_physiology]
      logical     :: ark_adaptive         = .true.      !< adaptive (embedded-error) vs fixed-substep march
      real(wp)    :: ark_rtol             = 1.0e-3_wp   !< adaptive relative tolerance (broadcast to all tol groups)
      !----- ONE master relative-accuracy dial for the WHOLE fast loop (§8c Layer 1): when > 0 it       !
      !      overrides every tolerance group's rtol -- the ARK march AND the nested soil-water /         !
      !      soil-energy / plant-hydraulics sub-solvers. 0 (default) => each keeps its own per-sub-       !
      !      solver value, i.e. byte-identical to before the unification. -------------------------------!
      real(wp)    :: rtol_all             = 0.0_wp      !< [-] 0 = unset; > 0 = the single accuracy target
      !----- The ABSOLUTE-tolerance companion to rtol_all. Every group's atol is physically scaled (J/kg,  !
      !      kg/kg, m3/m3, K, ...), so they cannot share one number the way rtol can; instead this is a     !
      !      dimensionless MULTIPLIER applied to the whole atol vector. It matters because the WRMS         !
      !      denominator is atol + rtol*|y|: tightening rtol_all ALONE saturates once atol dominates the    !
      !      denominator (measured: 1e-3 -> 1e-6 on the soil-water group only raises the error estimate     !
      !      ~4x, too little to force a substep), so the accuracy dial is only half-effective without it.   !
      !      1.0 (default) is an EXACT IEEE identity => byte-identical. -----------------------------------!
      real(wp)    :: atol_scale           = 1.0_wp      !< [-] multiplier on every group's atol (1 = unset)
      !----- WHERE IN THE SUB-STEP THE MET FORCING IS SAMPLED, as a fraction of dt_fast (§8f). The     !
      !      fast loop samples met at (isub - 1 + f)*dt_fast while column_prepass freezes every        !
      !      coefficient (gs/GPP/Rd/aerodynamics/radiation) on the STATE at t^n. f = 0.5 (the historic !
      !      default) is the better quadrature of the forcing alone, but it pairs t+dt/2 forcing with  !
      !      t^n state -- a FIRST-ORDER inconsistency in the frozen coefficients. Measured on the      !
      !      forced ERA5 census stand at dt_fast = 900 s, CAS-T RMSE traces a clean U in f with its    !
      !      minimum at f = 0 (forcing and state agreeing), worth ~2x. Kept at 0.5 by default so       !
      !      existing runs are unchanged; f = 0 is the state-consistent choice, and a true midpoint    !
      !      freeze (f = 0.5 WITH a state predictor) is the second-order version. -----------------!
      real(wp)    :: forcing_sample_frac  = 0.5_wp      !< [-] met sample point within the sub-step, in [0,1]
      !----- §8g SCHEME-ASYMMETRY GUARD: the CAS supersaturation (condensation) sink lives in
      !      surface_derivs, which ONLY the ARK stages reach -- so split-vs-ARK has been comparing two
      !      different models. Default .true. keeps the ARK unchanged; set .false. for a like-for-like
      !      scheme comparison. Whether the sink should also exist on the split path is a MODEL question.
      logical     :: cas_condensation     = .true.     !< apply the CAS supersaturation sink (ARK path)
      !----- PROCESS MASK (§5.1): which column processes actually EVOLVE. All true = the full column   !
      !      (default, byte-identical); flipping one off freezes that store so the driver integrates a   !
      !      REDUCED ODE. This is the process-complexity axis of the goal-(b) sweep. The mask type       !
      !      itself lives in meds_fast_types (with column_config_t); config carries plain logicals so     !
      !      shared/ does not gain a driver dependency. --------------------------------------------------!
      logical     :: mask_veg_energy = .true.   !< leaf + wood energy stores
      logical     :: mask_cas_energy = .true.   !< canopy-air-space enthalpy
      logical     :: mask_cas_vapour = .true.   !< canopy-air-space specific humidity
      logical     :: mask_cas_co2    = .true.   !< canopy-air-space CO2
      logical     :: mask_soil_heat  = .true.   !< soil thermal column
      logical     :: mask_soil_water = .true.   !< soil water column
      logical     :: mask_hydraulics = .true.   !< plant hydraulics (psi)
      integer(ik) :: step_controller      = CTRL_I      !< CTRL_I (default, legacy) | CTRL_PI (goal a; §9.3)
      integer(ik) :: error_level          = CTRL_L1_ADAPTIVE !< CTRL_L0_FIXED | CTRL_L1_ADAPTIVE (default) | CTRL_L2_STRICT
      real(wp)    :: ark_dt_init          = 0.0_wp      !< [s] initial adaptive substep (0 => dt_fast)
      integer(ik) :: ark_fixed_substep    = 4_ik        !< fixed substeps/dt_fast (GPU warp-uniform path)
      !----- The leaf<->CAS surface solve is EITHER the uncoupled single-BE pass OR the coupled 2x2   !
      !      Newton -- there is nothing in between. `ark_niter` was typed as an iteration cap but is    !
      !      only ever tested as `np <= 1` (column_be_stage), so every value > 1 behaved identically    !
      !      and the real cap is the NEWT_MAX = 4 parameter. It is a boolean, so it is spelled as one   !
      !      now (plan E4): `fast.ark_coupled`. The old `fast.ark_niter` is refused as retired         !
      !      (meds_config_keys), naming it. `ark_relax` was deleted -- it was vestigial on the Newton     !
      !      branch and read by nothing. --------------------------------------------------------------!
      logical     :: ark_coupled          = .true.      !< .false. = uncoupled single BE pass; .true. = 2x2 Newton
      !----- Sub-daily fast-loop diagnostic PROBE (opt-in; for the integrator/dt_fast evaluation): dumps !
      !      per-(patch,sub-step) CAS temp / GPP / ET / soil-top temp / leaf temp to a CSV. -------------!
      logical            :: fast_probe      = .false.
      character(len=256) :: fast_probe_file = 'fast_probe.csv'
      integer(ik) :: backend                     !< reporting only

      !----- Structural master switches. --------------------------------------------------!
      logical     :: do_cohort_fissfuse, do_patch_fissfuse, do_patch_disturbance

      !----- Cohort fusion / termination. -------------------------------------------------!
      integer(ik) :: max_cohort, n_cohort_fusion_iter
      real(wp)    :: cohort_size_tol_min, cohort_size_tol_max
      real(wp)    :: cohort_size_tol_mult        !< DERIVED (geometric multiplier)
      real(wp)    :: cohort_lai_cap, min_cohort_agb, negligible_nplant, split_eps
      logical     :: enable_cohort_fission

      !----- Vertical light profile for patch fusion (cumulative-LAI by height layer). ----!
      integer(ik)           :: n_height_layers
      real(wp), allocatable :: height_edges(:)   !< DERIVED (ascending interior edges [m])

      !----- Patch fusion / termination. --------------------------------------------------!
      integer(ik) :: max_patch, n_patch_fusion_iter
      real(wp)    :: patch_light_tol, patch_light_maxdev_factor, patch_diff_age_tol
      !----- The light-profile tolerance steps from patch_light_tol to this ceiling and no further:   !
      !      two patches more different than it stay apart even above max_patch, which makes the   !
      !      patch count a target rather than a hard limit.                                          !
      real(wp)    :: patch_light_tol_max
      real(wp)    :: patch_light_tol_mult        !< DERIVED (geometric multiplier)
      real(wp)    :: min_patch_area, patch_min_area_remain
      logical     :: enable_patch_fission

      !----- Patch disturbance, growth memory, recruitment, conservation. -----------------!
      real(wp) :: patch_disturbance_rate, disturbance_survive_height
      real(wp) :: growth_memory_days, min_recruit_size, conservation_tol

      !----- Initial conditions (init_mode: 0 bare | 1 census | 2 restart). ---------------!
      integer(ik)        :: init_mode
      character(len=256) :: init_restart_file, init_census_file
      !----- The soil state a run starts from when no state file restores one: every layer of   !
      !      every patch at this temperature and volumetric water content. The defaults are the   !
      !      constants the fast context carried before they were keys.  -------------------------!
      real(wp)           :: init_soil_temp  = 288.0_wp   !< [K]
      real(wp)           :: init_soil_theta = 0.30_wp    !< [m3/m3]
      !----- A restart's plastic leaf traits from THIS run's PFT file rather than the state file's   !
      !      (reacclimate_plant_traits): each cohort re-acclimated to the LAI above it with          !
      !      plasticity on, the PFT's top-of-canopy values with it off. Restart only. --------------!
      logical            :: init_reacclimate_traits = .false.

      !----- netCDF output. ---------------------------------------------------------------!
      character(len=256) :: state_output_dir, state_output_prefix
      logical     :: state_write_state
      integer(ik) :: state_interval_years_cfg

      !----- Parameter-config controls. ---------------------------------------------------!
      character(len=256) :: pft_config        !< path to the PFT config file (named in the main file)
      logical            :: override_derived  !< if .true., a [derived] block overwrites computed values

      !----- Leaf physiology: model selection (non-PFT). ----------------------------------!
      integer(ik) :: stomatal_model           !< SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form        !< TRESP_ARRHENIUS | TRESP_PEAKED
      integer(ik) :: colimitation             !< COLIM_MIN | COLIM_QUADRATIC
      logical     :: leaf_use_boundary_layer  !< if .true., couple via the leaf boundary layer (gb)
      !----- NON-STOMATAL (capacity) water-stress limb: a linear psi_leaf ramp downregulating       !
      !      Vcmax/Jmax/TPU (Sabot 2022 / Zhou 2013's beta_nonstomata). OFF by default (issue #47): !
      !      the term is rarely measured directly, its two parameters (wstress_psi_open/_close) are !
      !      weakly constrained, and it acts as a LINEAR AMPLIFIER on psi_leaf -- with the ramp of  !
      !      the shipped PFT file, slope 0.5 per MPa, so a 1 MPa error in psi_leaf becomes a 50%    !
      !      error in Vcmax. Before the transpiration corrector (#91), psi_leaf was not converged   !
      !      in dt_fast (daytime mean -0.23 MPa at 12.5 s vs -1.19 MPa at 900 s), and this limb     !
      !      turned that into a 33% GPP shift; with it off, daily GPP is dt_fast-independent to     !
      !      0.05%. The stomatal limb (beta_stomata, driven by psi_SOIL) is unaffected and stays on. !
      logical     :: leaf_wstress_nonstomatal  !< if .true., apply the psi_leaf capacity limb
      !----- Leaf physiology: shared biochemistry at 25 degC + Arrhenius/deactivation terms.-!
      real(wp) :: kc25, ko25, gstar25                   !< [Pa]    Michaelis constants + CO2 compensation point
      real(wp) :: ea_kc, ea_ko, ea_gstar                !< [J/mol] activation energies (Bernacchi et al. 2001)
      real(wp) :: ea_vcmax, ea_jmax, ea_rd              !< [J/mol] activation energies
      real(wp) :: hd_vcmax, hd_jmax, hd_rd              !< [J/mol] deactivation energies (peaked form)
      real(wp) :: ds_vcmax, ds_jmax, ds_rd              !< [J/mol/K] entropy terms (peaked form)
      !----- THERMAL ACCLIMATION (#176, Kattge & Knorr 2007). OFF by default, so the shipped        !
      !      behaviour is unchanged and `ds_vcmax`/`ds_jmax`/`jmax_vcmax_ratio` stay exactly the     !
      !      fixed values above. ON, the driver recomputes those three from a running-mean GROWTH    !
      !      temperature once per slow step, so a warm-grown and a cold-grown stand of the same PFT  !
      !      no longer share one temperature response. Only meaningful with temp_response = peaked:  !
      !      the Arrhenius form has no dS to shift, and validate_config says so rather than letting  !
      !      the flag be silently inert.  -------------------------------------------------------------!
      logical  :: leaf_thermal_acclimation = .false.
      real(wp) :: acclim_ds_vcmax_a = 668.39_wp, acclim_ds_vcmax_b = 1.07_wp   !< [J/mol/K], [J/mol/K2]
      real(wp) :: acclim_ds_jmax_a  = 659.70_wp, acclim_ds_jmax_b  = 0.75_wp
      real(wp) :: acclim_jv_a       = 2.59_wp,   acclim_jv_b       = 0.035_wp  !< [-], [1/K]
      real(wp) :: acclim_window_days = 30.0_wp   !< [day] growth-temperature running-mean window
      real(wp) :: o2_mol_frac                           !< [mol/mol] atmospheric O2 mole fraction
      real(wp) :: leaf_absorptance                      !< [--] leaf PAR absorptance (for electron transport)
      real(wp) :: phi_psii                              !< [--] low-light electron yield (J slope 0.5*phi_psii/photon)
      !< [kPa] the leaf-to-air VPD the Medlyn stomatal model uses at least: g1/sqrt(D) is undefined at D = 0.
      !< CLM5's value (PhotosynthesisMod floors the Medlyn VPD at 50 Pa).
      real(wp) :: medlyn_vpd_min = 0.05_wp

      !----- Carbon growth: the model's demographic growth is carbon-prognostic (wood_carbon is  !
      !       the size anchor, driven by NPP). gpp_ref is the stub GPP when the fast loop is off.  !
      real(wp)    :: gpp_ref                  !< [kgC/m2 leaf/yr] stub GPP per unit leaf area (carbon mode)

      !----- Light trait plasticity ([trait_dynamics]). OPT-IN: default .false. => cohort leaf traits !
      !       (sla/vcmax25/rd25/llspan) stay at their top-of-canopy PFT values (bit-identical to the   !
      !       static path). When ON, the slow-loop driver acclimates them to cumulative LAI above.     !
      logical     :: trait_plasticity_on = .false.

      !----- Meteorological forcing ([forcing]/[site]). OPT-IN: forcing_on default .false. (the   !
      !       whole [forcing] block is gated on it), so a config with no [forcing] block runs the   !
      !       constant-forcing MVP unchanged. Defaults are the Ithaca NY / ERA5-Land reference.     !
      type(forcing_config_t) :: forcing

      !----- Diagnostic-aggregation output ([output]). OPT-IN: enabled default .false. (a config    !
      !       with no [output] block emits no diagnostic stream). The per-variable overrides        !
      !       live in the optional meds_io_config.toml named by output%io_config (§6, MEDS_IO_DESIGN). !
      type(output_config_t) :: output

      !----- PFT traits. ------------------------------------------------------------------!
      type(pft_table_t) :: pft

      !----- Global allometry coefficients (populated by the loader; installed by derive_parameters). !
      type(allometry_config_t) :: allom

      !----- Plant-hydraulics parameters ([hydraulics], opt-in; defaults = MVP placeholders). ------!
      type(soil_column_config_t) :: soil_column  !< [soil_column] the physical soil column
      type(hydraulics_config_t) :: hydraulics

      !----- Fast-loop biophysics run-config ([soil]/[energy]/[snow]/[aerodynamics], all opt-in;    !
      !       defaults = meds_biophysics_opts placeholders). build_fast_context copies each verbatim !
      !       into the column config (col_config%soil_water_opts/energy/snow/aero); an absent block is a no-op. ------!
      type(soil_opts_t)   :: soil        !< [soil]         soil-water Richards solver opts (-> col_config%soil_water_opts)
      type(energy_opts_t) :: energy      !< [energy]       soil-thermal solver opts       (-> col_config%energy)
      type(snow_params_t) :: snow        !< [snow]         snow physical parameter table  (-> col_config%snow)
      type(aero_cfg_t)    :: aero        !< [aerodynamics] canopy-aerodynamics constants  (-> col_config%aero)

      !----- Slow soil-carbon matrix ([soil_carbon], ON by default; MEDS_SLOW_DYNAMICS_DESIGN.md   !
      !      Part II B0). soil_carbon_on gates the FEATURE. It does NOT gate whether the fields       !
      !      are required: every key is a DEFAULTED read (like [snow]), falling back to its             !
      !      ED2-verified in-type default, so the feature needs no TOML edits beyond this flag.        !
      !                                                                                          !
      !      DEFAULT .true. since 2026-09-11. Off is not a coarser soil-carbon model, it is NO soil    !
      !      carbon: litter is discarded at the slow step and patch_heterotrophic_respiration returns  !
      !      rh = 0, so Reco carries only its autotrophic limb. Measured over a year at Ithaca, off    !
      !      reports annual-mean NEE at -4.657 against -2.464 umol/m2/s -- an 89% stronger apparent    !
      !      sink -- with Rh identically zero against 0.833 kgC/m2/yr, for no wall-clock saving.       !
      !      It also kept two real defects (PRs #139, #140) out of every code path anyone ran, for as  !
      !      long as they existed, because no shipped config turned it on.  --------------------------!
      logical            :: soil_carbon_on = .true.
      type(decomp_opts_t) :: soil_carbon   !< [soil_carbon] decomposition selectors + rate parameters
      !----- Cold-start spin-up ([soil_carbon], consumed only when soil_carbon_on): zero-init        !
      !      (default) leaves every pool at 0, matching bare-ground philosophy; steady-state solves     !
      !      SASU (solve_soil_carbon_steady_state) from a scalar climatological environmental factor    !
      !      and a constant litter-input estimate (pools 5-7 get zero input, matching build_litter_    !
      !      input's own behavior -- only 1-4 are ever populated from real litter). ---------------------!
      logical  :: soil_carbon_spinup_steady      = .false.
      real(wp) :: soil_carbon_spinup_xi          = 1.0_wp   !< [-] climatological mean env. scalar (broadcast, all pools)
      real(wp) :: soil_carbon_spinup_labile_grnd = 0.0_wp   !< [kgC/m2/day] u_bar(1), steady-state litter estimate
      real(wp) :: soil_carbon_spinup_labile_soil = 0.0_wp   !< [kgC/m2/day] u_bar(2)
      real(wp) :: soil_carbon_spinup_struct_grnd = 0.0_wp   !< [kgC/m2/day] u_bar(3)
      real(wp) :: soil_carbon_spinup_struct_soil = 0.0_wp   !< [kgC/m2/day] u_bar(4)
   end type meds_config_t

contains

   !---------------------------------------------------------------------------------------!
   ! Compute the DERIVED configuration from the (already-loaded) primary parameters: the      !
   ! timestep dt, the geometric cohort-fusion tolerance multiplier, and the evenly-spaced      !
   ! height-layer edges (0 to the tallest PFT's hgt_max). The mortality-hazard parameters are   !
   ! derived separately (derive_pft_rates).                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine derive_config(cfg)
      type(meds_config_t), intent(inout) :: cfg
      integer(ik) :: i

      !----- Slow-process timestep in years (demography currency); source is now dt_slow. --!
      cfg%dt_years = cfg%dt_slow / yr_sec

      !----- Fast sub-steps nested within one slow step (>=1; guarded against dt_fast = 0). --!
      if (cfg%dt_fast > 0.0_wp) then
         cfg%n_fast_per_slow = max(1_ik, nint(cfg%dt_slow / cfg%dt_fast, ik))
      else
         cfg%n_fast_per_slow = 1_ik
      end if

      !----- Geometric tolerance growth from min to max over niter iterations. ------------!
      if (cfg%n_cohort_fusion_iter > 1_ik) then
         cfg%cohort_size_tol_mult = (cfg%cohort_size_tol_max / cfg%cohort_size_tol_min)             &
                               ** (1.0_wp / real(cfg%n_cohort_fusion_iter - 1_ik, wp))
      else
         cfg%cohort_size_tol_mult = 1.0_wp
      end if
      !----- The patch-fusion tolerance grows the same way, from patch_light_tol to its ceiling. ----!
      if (cfg%n_patch_fusion_iter > 1_ik .and. cfg%patch_light_tol > 0.0_wp) then
         cfg%patch_light_tol_mult = (cfg%patch_light_tol_max / cfg%patch_light_tol)                  &
                               ** (1.0_wp / real(cfg%n_patch_fusion_iter - 1_ik, wp))
      else
         cfg%patch_light_tol_mult = 1.0_wp
      end if

      !----- Evenly spaced height-layer edges from 0 to the tallest PFT's height cap. ------!
      if (allocated(cfg%height_edges)) deallocate(cfg%height_edges)
      allocate(cfg%height_edges(cfg%n_height_layers - 1_ik))
      do i = 1_ik, cfg%n_height_layers - 1_ik
         cfg%height_edges(i) = real(i, wp) * maxval(cfg%pft%hgt_max(1:cfg%pft%n))                &
                               / real(cfg%n_height_layers, wp)
      end do
   end subroutine derive_config

   !---------------------------------------------------------------------------------------!
   ! Install + derive EVERY parameter that is a function of the primary (loaded) config: the   !
   ! global allometry coefficients (into meds_allometry's protected state), the derived run      !
   ! scalars (derive_config), the wood-density mortality hazard (derive_pft_rates), and the       !
   ! leaf-capacity ratios (derive_leaf_params). The single consolidation point for both the        !
   ! production loader (load_meds_config) and the test builder (build_test_config), so the          !
   ! derivation sequence lives in ONE place. The four calls are mutually order-independent (each     !
   ! reads only primary loaded fields; set_allometry installs a global none of them consume).        !
   ! Callers apply any [derived] override and validate AFTER this returns.                            !
   !---------------------------------------------------------------------------------------!
   subroutine derive_parameters(cfg)
      type(meds_config_t), intent(inout) :: cfg
      call set_allometry(cfg%allom%b1Ht, cfg%allom%b2Ht, cfg%allom%agb_c1, cfg%allom%agb_c2,       &
                         cfg%allom%ca_b1, cfg%allom%ca_b2, cfg%allom%lai_b1, cfg%allom%lai_b2,      &
                         cfg%allom%light_ext)
      call derive_config(cfg)
      call derive_pft_rates(cfg%pft)
      call derive_leaf_params(cfg%pft)
   end subroutine derive_parameters

   !---------------------------------------------------------------------------------------!
   ! Number of time steps spanned by the growth-memory window (>=1): the size of the per-    !
   ! cohort moving-average ring buffer. Derived from the memory window [days] and the step.   !
   !---------------------------------------------------------------------------------------!
   pure integer(ik) function growth_window_steps(cfg) result(nw)
      type(meds_config_t), intent(in) :: cfg
      nw = max(1_ik, nint(cfg%growth_memory_days / (cfg%dt_years * yr_day), ik))
   end function growth_window_steps

   !----- A PFT's leaf turgor-loss point [MPa], from its own pressure-volume traits: the [pft]       !
   !      leaf_pi0 and leaf_elastic_mod where given, else the shared [hydraulics] ones. That is the   !
   !      rule apply_hydraulics_config uses for the plant-water solver, so the stomatal closure and   !
   !      the phenology's drought counter see the curve the solver stores leaf water on. ------------!
   pure real(wp) function pft_leaf_psi_tlp(cfg, ipft) result(psi_tlp)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      real(wp) :: pi0, elastic_mod
      pi0 = cfg%hydraulics%leaf_pi0 ; elastic_mod = cfg%hydraulics%leaf_elastic_mod
      if (allocated(cfg%pft%hyd_leaf_pi0)) then
         if (cfg%pft%hyd_leaf_pi0(ipft) > HYD_UNSET) pi0 = cfg%pft%hyd_leaf_pi0(ipft)
      end if
      if (allocated(cfg%pft%hyd_leaf_elastic_mod)) then
         if (cfg%pft%hyd_leaf_elastic_mod(ipft) > HYD_UNSET) elastic_mod = cfg%pft%hyd_leaf_elastic_mod(ipft)
      end if
      psi_tlp = pv_psi_tlp(pi0, elastic_mod)
   end function pft_leaf_psi_tlp

   !----- Where a PFT's stomatal water stress begins [MPa of predawn leaf potential]: the [pft]       !
   !      stomata_psi_onset where given, else half the PFT's turgor-loss point. Above the onset the    !
   !      stomata feel no stress, so the predawn potential a tall tree has in wet soil -- its height's  !
   !      gravity head, -0.34 MPa at 35 m -- does not read as drought. -------------------------------!
   pure real(wp) function pft_stomata_psi_onset(cfg, ipft) result(psi_onset)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      psi_onset = 0.5_wp * pft_leaf_psi_tlp(cfg, ipft)
      if (allocated(cfg%pft%stomata_psi_onset)) then
         if (cfg%pft%stomata_psi_onset(ipft) > HYD_UNSET) psi_onset = cfg%pft%stomata_psi_onset(ipft)
      end if
   end function pft_stomata_psi_onset

   !---------------------------------------------------------------------------------------!
   ! Validate a configuration; halt on a setting that would corrupt the run.               !
   !---------------------------------------------------------------------------------------!
   subroutine validate_config(cfg)
      type(meds_config_t), intent(in) :: cfg
      character(len=*), parameter :: tag = 'meds_config: '

      !----- [hydraulics]: the rooting traits set the root profile over the soil column, and the    !
      !      segment conductance needs a positive conductivity and path factor, shared and per PFT. -!
      if (cfg%hydraulics%root_beta <= 0.0_wp .or. cfg%hydraulics%root_beta >= 1.0_wp)             &
         error stop tag//'hydraulics.root_beta must lie in (0, 1)'
      if (cfg%hydraulics%root_depth <= 0.0_wp) error stop tag//'hydraulics.root_depth <= 0'
      if (cfg%hydraulics%conductance == HYD_CONDUCTANCE_SEGMENT) then
         if (cfg%hydraulics%wood_kmax <= 0.0_wp .or. cfg%hydraulics%vessel_curl <= 0.0_wp)         &
            error stop tag//'hydraulics.conductance = "segment" needs wood_kmax > 0 and vessel_curl > 0'
         if (allocated(cfg%pft%hyd_wood_kmax)) then
            if (any(cfg%pft%hyd_wood_kmax > HYD_UNSET .and. cfg%pft%hyd_wood_kmax <= 0.0_wp))       &
               error stop tag//'pft.wood_kmax must be > 0 with conductance = "segment"'
         end if
         if (allocated(cfg%pft%hyd_vessel_curl)) then
            if (any(cfg%pft%hyd_vessel_curl > HYD_UNSET .and. cfg%pft%hyd_vessel_curl <= 0.0_wp))   &
               error stop tag//'pft.vessel_curl must be > 0 with conductance = "segment"'
         end if
      end if
      !----- [soil] ground optics. --------------------------------------------------------------!
      if (cfg%soil%ground_albedo_vis < 0.0_wp .or. cfg%soil%ground_albedo_vis >= 1.0_wp .or.       &
          cfg%soil%ground_albedo_nir < 0.0_wp .or. cfg%soil%ground_albedo_nir >= 1.0_wp)           &
         error stop tag//'soil.ground_albedo_vis and ground_albedo_nir must lie in [0, 1)'
      if (cfg%soil%ground_emissivity <= 0.0_wp .or. cfg%soil%ground_emissivity > 1.0_wp)           &
         error stop tag//'soil.ground_emissivity must lie in (0, 1]'

      !----- [soil_column]. Every one of these produces a silently WRONG column rather than a     !
      !      crash: a layer count over the compile-time ceiling writes past the active region, a   !
      !      non-positive depth or a theta_sat <= theta_res divides by zero in the retention curve, !
      !      and a van Genuchten n <= 1 makes the curve's exponent negative. Fail loud instead. ----!
      associate (sc => cfg%soil_column)
         if (sc%n_layer < 1_ik .or. sc%n_layer > n_soil_layer_max)                                &
            error stop tag//'soil_column.n_layer outside 1..n_soil_layer_max'
         if (sc%depth <= 0.0_wp)          error stop tag//'soil_column.depth <= 0'
         if (sc%grid_growth < 0.0_wp)     error stop tag//'soil_column.grid_growth < 0'
         if (sc%theta_sat <= sc%theta_res) error stop tag//'soil_column.theta_sat <= theta_res'
         if (sc%theta_res < 0.0_wp)       error stop tag//'soil_column.theta_res < 0'
         if (cfg%init_soil_theta <= sc%theta_res .or. cfg%init_soil_theta > sc%theta_sat)             &
            error stop tag//'init.soil_theta must lie in (soil_column.theta_res, theta_sat]'
         if (sc%ksat <= 0.0_wp)           error stop tag//'soil_column.ksat <= 0'
         !----- curve_par_a and curve_par_n mean different things per family: van Genuchten's alpha  !
         !      [1/m] > 0 and n > 1, Campbell's air-entry suction psi_sat [m] < 0 and exponent b > 0.  !
         !      The other family's pair makes the curve raise a negative base to a fractional power. -!
         if (sc%retention == SOIL_RETENTION_VG) then
            if (sc%curve_par_a <= 0.0_wp)                                                          &
               error stop tag//'soil_column.curve_par_a (van Genuchten alpha, 1/m) must be > 0'
            if (sc%curve_par_n <= 1.0_wp)                                                          &
               error stop tag//'soil_column.curve_par_n must exceed 1 for van Genuchten'
         else
            if (sc%curve_par_a >= 0.0_wp)                                                          &
               error stop tag//'soil_column.curve_par_a (Campbell psi_sat, m) must be < 0'
            if (sc%curve_par_n <= 0.0_wp)                                                          &
               error stop tag//'soil_column.curve_par_n (Campbell b) must be > 0'
         end if
         if (sc%psi_fc >= 0.0_wp)         error stop tag//'soil_column.psi_fc must be negative (a suction head)'
         if (sc%solid_conductivity <= 0.0_wp .or. sc%dry_conductivity <= 0.0_wp)                  &
            error stop tag//'soil_column conductivities must be positive'
         if (sc%dry_heat_capacity <= 0.0_wp) error stop tag//'soil_column.dry_heat_capacity <= 0'
      end associate

      !----- [energy] bottom thermal BC (#145). The Dirichlet anchor needs a SITE temperature and a  !
      !      depth BELOW the bottom node; neither can be guessed. `deep_temp` is the mean annual soil !
      !      temperature below the damping depth (close to mean annual AIR temperature), and an error !
      !      in it is a steady flux into the column base, so a silent default would be a silent mean- !
      !      annual bias -- exactly the defect #145 exists to remove. `deep_depth` at or above the    !
      !      bottom NODE gives a non-positive conduction length, i.e. a division by zero or a sign    !
      !      flip that would pump heat the wrong way.  ------------------------------------------------!
      if (cfg%energy%bottom_bc == ENERGY_BC_DIRICHLET) then
         block
            type(soil_params_t) :: sp
            if (cfg%energy%deep_temp <= 0.0_wp)                                                  &
               error stop tag//'energy.bottom_bc = "dirichlet" requires energy.deep_temp [K] -- '// &
                          'the mean annual soil temperature below the damping depth. It is a site '//&
                          'property, like latitude; there is no safe default.'
            if (cfg%energy%deep_temp < 200.0_wp .or. cfg%energy%deep_temp > 330.0_wp)             &
               error stop tag//'energy.deep_temp outside 200..330 K -- it is an absolute temperature'
            !----- Build the real grid rather than re-deriving the bottom node here: the exponential !
            !      generator lives in ONE place and a second copy would drift from it silently. -----!
            associate (sc => cfg%soil_column)
               call build_soil_hydr_params(sc%n_layer, sc%retention, sc%depth, sc%grid_growth,   &
                    sc%theta_sat, sc%theta_res, sc%ksat, sc%curve_par_a, sc%curve_par_n,         &
                    cfg%hydraulics%root_beta, cfg%hydraulics%root_depth, sc%psi_fc, sp)
            end associate
            if (cfg%energy%deep_depth <= abs(sp%z_node(cfg%soil_column%n_layer)))                &
               error stop tag//'energy.deep_depth must lie BELOW the bottom soil node '//          &
                          '(it is measured down from the surface, and the anchor conducts to it)'
         end block
      end if

      !----- Co-limitation and hyperbola curvatures must lie in (0, 1]. A curvature of zero makes the  !
      !      smoothing quadratic degenerate and a value above one has no root in the physical branch,   !
      !      so either produces a silently wrong assimilation rather than a crash. This guard exists    !
      !      because #118's own fixture reached the solver with the two new C3 curvatures unset. -------!
      block
         integer(ik) :: pf
         do pf = 1_ik, cfg%pft%n
            if (cfg%pft%theta_j(pf)     <= 0.0_wp .or. cfg%pft%theta_j(pf)     > 1.0_wp)         &
               error stop tag//'pft.theta_j outside (0, 1]'
            if (cfg%pft%theta_cj_c3(pf) <= 0.0_wp .or. cfg%pft%theta_cj_c3(pf) > 1.0_wp)         &
               error stop tag//'pft.theta_cj_c3 outside (0, 1]'
            if (cfg%pft%theta_ip_c3(pf) <= 0.0_wp .or. cfg%pft%theta_ip_c3(pf) > 1.0_wp)         &
               error stop tag//'pft.theta_ip_c3 outside (0, 1]'
            if (cfg%pft%theta_cj_c4(pf) <= 0.0_wp .or. cfg%pft%theta_cj_c4(pf) > 1.0_wp)         &
               error stop tag//'pft.theta_cj_c4 outside (0, 1]'
            if (cfg%pft%theta_ic_c4(pf) <= 0.0_wp .or. cfg%pft%theta_ic_c4(pf) > 1.0_wp)         &
               error stop tag//'pft.theta_ic_c4 outside (0, 1]'
         end do
      end block

      !----- Canopy optics. reflect + transmit is the single-scatter albedo: at or above 1 the      !
      !      two-stream conserves or creates energy in a scattering layer and the solve stops        !
      !      meaning anything, so it is a hard error rather than a clamp. Emissivity outside (0,1]   !
      !      and a leaf-angle mean outside (0,90) are equally unphysical.  ---------------------------!
      block
         integer(ik) :: pf
         do pf = 1_ik, cfg%pft%n
            if (cfg%pft%leaf_reflect_vis(pf) + cfg%pft%leaf_transmit_vis(pf) >= 1.0_wp)           &
               error stop tag//'pft.leaf_reflect_vis + leaf_transmit_vis >= 1'
            if (cfg%pft%leaf_reflect_nir(pf) + cfg%pft%leaf_transmit_nir(pf) >= 1.0_wp)           &
               error stop tag//'pft.leaf_reflect_nir + leaf_transmit_nir >= 1'
            if (cfg%pft%wood_reflect_vis(pf) + cfg%pft%wood_transmit_vis(pf) >= 1.0_wp)           &
               error stop tag//'pft.wood_reflect_vis + wood_transmit_vis >= 1'
            if (cfg%pft%wood_reflect_nir(pf) + cfg%pft%wood_transmit_nir(pf) >= 1.0_wp)           &
               error stop tag//'pft.wood_reflect_nir + wood_transmit_nir >= 1'
            if (cfg%pft%leaf_emissivity(pf) <= 0.0_wp .or. cfg%pft%leaf_emissivity(pf) > 1.0_wp)  &
               error stop tag//'pft.leaf_emissivity outside (0,1]'
            if (cfg%pft%wood_emissivity(pf) <= 0.0_wp .or. cfg%pft%wood_emissivity(pf) > 1.0_wp)  &
               error stop tag//'pft.wood_emissivity outside (0,1]'
            if (cfg%pft%leaf_clumping(pf) <= 0.0_wp .or. cfg%pft%leaf_clumping(pf) > 1.0_wp)      &
               error stop tag//'pft.leaf_clumping outside (0,1]'
            if (cfg%pft%wood_clumping(pf) <= 0.0_wp .or. cfg%pft%wood_clumping(pf) > 1.0_wp)      &
               error stop tag//'pft.wood_clumping outside (0,1]'
            if (cfg%pft%leaf_angle_mean(pf) <= 0.0_wp .or. cfg%pft%leaf_angle_mean(pf) >= 90.0_wp) &
               error stop tag//'pft.leaf_angle_mean outside (0,90) degrees'
            if (cfg%pft%leaf_angle_std(pf) < 0.0_wp)                                              &
               error stop tag//'pft.leaf_angle_std < 0'
         end do
      end block

      if (cfg%pft%n < 1_ik)                          error stop tag//'empty PFT table'
      if (.not. time_lt(cfg%start_time, cfg%end_time)) error stop tag//'end_time must be after start_time'
      if (cfg%dt_slow <= 0.0_wp)                        error stop tag//'dt_slow <= 0'
      if (cfg%fast_biophysics_on) then
         if (cfg%dt_fast <= 0.0_wp)                     error stop tag//'dt_fast <= 0'
         !----- dt_fast USED TO BE A STABILITY PARAMETER. It is now an ACCURACY parameter, and the     !
         !      binding quantity has changed with it.                                                  !
         !                                                                                          !
         !      WAS: the CAS<->atmosphere conductance was frozen for the whole step, and that one lag  !
         !      drove a sustained period-2 canopy-air oscillation (peak-to-peak ~8 K at 900 s) that no  !
         !      conservation ledger detected. The surface layer is now re-solved at every stage        !
         !      (meds_fast_types%mo_live), which removes it: the freeze-cadence multiplier goes from   !
         !      Phi' = -23.2 to -0.14 at 900 s, and peak-to-peak from 7.7 K to 0.10 K.                 !
         !                                                                                          !
         !      IS: what still degrades with dt_fast is the CARBON and WATER budget, because the leaf  !
         !      gas-exchange pre-pass -- and with it the leaf water potential that sets stomatal       !
         !      conductance -- is still frozen at state^n. MEASURED on a high-LAI sunlit stand against !
         !      a 12.5 s reference:                                                                   !
         !                                                                                          !
         !        dt_fast     150 s     225 s     300 s     450 s     900 s                            !
         !        GPP        -3.8%     -8.4%    -12.6%    -19.8%    -33.1%                             !
         !        ET         -2.6%     -5.8%     -8.7%    -14.5%    -23.8%                             !
         !        CAS-T RMSE  0.02 K    0.04 K    0.06 K    0.09 K    0.16 K                           !
         !                                                                                          !
         !      Attributed: freezing the plant-hydraulics store (mask%hydraulics) removes the GPP      !
         !      bias entirely (-33.1% -> +0.9% at 900 s, and flat in dt_fast), so it is the LAGGED     !
         !      leaf-water-potential feedback, not a quadrature error. Until that lag gets the same    !
         !      per-stage treatment the conductances just got, a long dt_fast is cheap in wall-clock   !
         !      and expensive in carbon. -------------------------------------------------------------!
         !----- ...and the size of that bias depends ENTIRELY on whether the non-stomatal        !
         !      (capacity) water-stress limb is active. It is a linear ramp on Vcmax/Jmax/TPU in    !
         !      psi_leaf, so it AMPLIFIES the psi error into a carbon error. With it off (the        !
         !      default, issue #47) daily GPP is dt_fast-independent to 0.05% and ET to ~1%. The     !
         !      table below predates the transpiration corrector (#91), which removed most of the    !
         !      psi error it amplifies; it has not been re-measured (numerical_scheme.md 5a). ------!
         if (cfg%dt_fast > 225.0_wp .and. cfg%leaf_wstress_nonstomatal) then
            print '(a)', 'WARNING [meds_config]: dt_fast > 225 s WITH the non-stomatal water-stress'
            print '(a)', '  limb on ([leaf_physiology].wstress_nonstomatal) biases the CARBON budget.'
            print '(a)', '  That limb is a linear amplifier on any psi_leaf error. Measured on a high-LAI'
            print '(a)', '  sunlit stand vs a 12.5 s reference, before the transpiration corrector (#91):'
            print '(a)', '    dt_fast    150 s   300 s   450 s   900 s'
            print '(a)', '    GPP       -3.8%  -12.6%  -19.8%  -33.1%'
            print '(a)', '    ET        -2.6%   -8.7%  -14.5%  -23.8%'
            print '(a)', '  With the limb off those become -0.0% / -1.2% at 900 s. No conservation'
            print '(a)', '  budget detects either. Prefer dt_fast <= 225 s, or leave the limb off.'
         else if (cfg%dt_fast > 900.0_wp) then
            print '(a)', 'WARNING [meds_config]: dt_fast > 900 s is beyond the measured range. At'
            print '(a)', '  900 s the canopy air is stable and carbon and daily leaf water potential'
            print '(a)', '  converge; above it neither has been measured.'
         end if
         if (cfg%dt_fast > 1800.0_wp) then
            print '(a)', 'WARNING [meds_config]: dt_fast > 1800 s is outside the measured range entirely.'
         end if
         if (cfg%dt_fast > cfg%dt_slow)                 error stop tag//'dt_fast > dt_slow'
         if (abs(cfg%dt_slow / cfg%dt_fast - real(nint(cfg%dt_slow / cfg%dt_fast, ik), wp)) > 1.0e-6_wp) &
            error stop tag//'dt_slow must be an integer multiple of dt_fast'
         if (cfg%time_integrator /= INTEG_ARK .and. cfg%time_integrator /= INTEG_RK45)          &
            error stop tag//'time_integrator out of range'
         if (cfg%rtol_all < 0.0_wp)           error stop tag//'rtol_all < 0 (0 = unset)'
         if (cfg%atol_scale <= 0.0_wp)        error stop tag//'atol_scale <= 0 (1 = unset)'
         if (cfg%forcing_sample_frac < 0.0_wp .or. cfg%forcing_sample_frac > 1.0_wp)                 &
            error stop tag//'forcing_sample_frac outside [0,1]'
         if (cfg%step_controller /= CTRL_I .and. cfg%step_controller /= CTRL_PI)                &
            error stop tag//'step_controller out of range (I|PI)'
         if (cfg%error_level /= CTRL_L0_FIXED .and. cfg%error_level /= CTRL_L1_ADAPTIVE .and.    &
             cfg%error_level /= CTRL_L2_STRICT)                                                 &
            error stop tag//'error_level out of range (L0|L1|L2)'
         if (cfg%ark_rtol <= 0.0_wp)          error stop tag//'ark_rtol <= 0'
         if (cfg%ark_fixed_substep < 1_ik)    error stop tag//'ark_fixed_substep < 1'
         !----- RK45 AT PRODUCTION CADENCE CARRIES AN ERROR THE DEFAULT PATH DOES NOT, and nothing  !
         !      used to say so (#160). The transpiration corrector that fixed a ~1 MPa psi_leaf      !
         !      error -- a 314x improvement, PR #91 -- lives in advance_water_mass_full, which the   !
         !      ARK path calls and the RK45 path does not. A run that selects rk45 at a coarse       !
         !      dt_fast therefore gets a psi_leaf the default scheme would not produce, silently.    !
         !      WARN rather than stop: rk45 is the deliberate accuracy baseline and is exactly what  !
         !      you want at a fine dt_fast, where the uncorrected error is small. The threshold is   !
         !      the cadence at which psi_leaf WITHOUT the corrector is known not to converge          !
         !      (docs/science/numerical_scheme.md section 7 item 1: daytime-mean -0.23 MPa at 12.5 s  !
         !      against -1.19 MPa at 900 s, measured before #91).                                                                 !
         if (cfg%time_integrator == INTEG_RK45 .and. cfg%dt_fast > RK45_UNCORRECTED_DT_WARN) then
            write(*,'(a)')    ' meds_config: WARNING -- time_integrator = "rk45" at dt_fast > 300 s.'
            write(*,'(a,f8.1,a)') '   dt_fast = ', cfg%dt_fast, ' s.'
            write(*,'(a)')    '   The transpiration corrector (PR #91) is ARK-only, so this run carries a'
            write(*,'(a)')    '   psi_leaf error the default scheme does not -- up to ~1 MPa at 900 s. RK45'
            write(*,'(a)')    '   is the ACCURACY BASELINE and is meant for a fine dt_fast; anything keyed'
            write(*,'(a)')    '   to psi_leaf (hydraulic stress, phenology cues) inherits the error here.'
            write(*,'(a)')    '   Use time_integrator = "ark" for production, or reduce dt_fast.'
         end if
         !----- §7 C3: the sub-daily probe is the ONE fast-loop writer that is not thread-safe --      !
         !      it holds a `save`d unit and appends to a shared file, so under threads its rows would   !
         !      interleave in thread-arrival order and stop being a reproducible diagnostic. Serializing !
         !      the write would fix the corruption but not the ORDER, which is the property the probe    !
         !      exists for, so the two are mutually exclusive by construction rather than by luck. ------!
         if (cfg%fast_probe .and. cfg%n_threads > 1_ik)                                            &
            error stop tag//'fast_probe requires n_threads = 1 (the probe writes one shared, &
                           &order-significant CSV; see plan sec 7 C3)'
      end if
      if (cfg%n_threads < 1_ik)               error stop tag//'n_threads < 1'
      !----- The FAST output tier closes a window every fast_interval_steps sub-steps, counted over   !
      !      the run. A window must not straddle two slow steps: the stand can be restructured between !
      !      them, and a window's cohort and patch slots are fixed when it opens (#312 O10).  ---------!
      if (cfg%fast_biophysics_on .and. cfg%output%enabled) then
         if (cfg%output%fast_interval_steps < 1_ik)                                                &
            error stop tag//'output.fast_interval_steps < 1'
         if (mod(cfg%n_fast_per_slow, cfg%output%fast_interval_steps) /= 0_ik)                    &
            error stop tag//'output.fast_interval_steps must divide the fast steps in a slow step '// &
                            '(dt_slow / dt_fast), so no fast-tier window straddles two slow steps'
      end if
      !----- REGION MODE (MEDS_POLYGON_RUNTIME_PLAN.md §9). Every polygon reads its own cell of the   !
      !      ED_ERA5land archive, so a region needs the archive, live forcing and the fast loop. Until  !
      !      the ragged restart exists (R4) a region starts from bare ground and writes no            !
      !      checkpoints, and the one-file probe stays off (every polygon would write into it). ------!
      if (cfg%run_mode == RUN_MODE_REGION) then
         if (.not. (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on))                         &
            error stop tag//'[run].mode = "region" needs fast.fast_biophysics_on and forcing.forcing_on'
         if (cfg%forcing%backend /= MET_BACKEND_ED_ERA5LAND)                                       &
            error stop tag//'[run].mode = "region" needs forcing.format = "ED_ERA5land"'
         if (cfg%fast_probe)                                                                       &
            error stop tag//'[run].mode = "region" cannot write fast.fast_probe (one CSV per run)'
         if (cfg%init_mode /= INIT_BARE)                                                           &
            error stop tag//'[run].mode = "region" starts from bare ground (init.init_mode = 0) until R4'
         if (cfg%state_write_state)                                                                &
            error stop tag//'[run].mode = "region" writes no checkpoints (state.write_state = false) until R4'
         associate (b => cfg%region%box_nwse)
            if (.not. (b(1) > b(3) .and. b(1) <= 90.0_wp .and. b(3) >= -90.0_wp))               &
               error stop tag//'region.box_nwse needs -90 <= S < N <= 90 ([N, W, S, E])'
            if (any(abs(b([2, 4])) > 360.0_wp))                                                    &
               error stop tag//'region.box_nwse longitudes must lie within [-360, 360]'
         end associate
         if (cfg%region%land_fraction_min < 0.0_wp .or. cfg%region%land_fraction_min > 1.0_wp)     &
            error stop tag//'region.land_fraction_min must lie in [0, 1]'
         !----- A region loads a month's forcing once, before the month (region_step_month), so the  !
         !      recycle seam must fall on a month boundary: a window starting at 00:00 or 01:00 on the !
         !      1st. A site run takes a window starting at any record stamp. ------------------------!
         if (cfg%forcing%recycle .and. time_valid(cfg%forcing%recycle_start)) then
            associate (a => cfg%forcing%recycle_start)
               if (a%day /= 1_ik .or. a%hour > 1_ik .or. a%minute /= 0_ik .or. a%second /= 0_ik)    &
                  error stop tag//'[run].mode = "region" needs forcing.recycle_start at 00:00 or '//  &
                                  '01:00 on the 1st of a month'
            end associate
         end if
      end if
      !----- Forcing: the forcing's own heights must be physical, and an open-terrain wind must have  !
      !      been made above its exposure roughness and below its blending height. No reference height  !
      !      has to clear the canopy any more: the forcing is moved to each patch's canopy-air top.     !
      if (cfg%forcing%forcing_on) then
         if (cfg%forcing%tq_height <= 0.0_wp)   error stop tag//'forcing.tq_height must be > 0'
         if (cfg%forcing%wind_height <= 0.0_wp) error stop tag//'forcing.wind_height must be > 0'
         if (cfg%forcing%wind_exposure == WIND_EXPOSURE_OPEN_TERRAIN) then
            if (cfg%forcing%wind_exposure_z0 <= 0.0_wp .or.                                       &
                cfg%forcing%wind_exposure_z0 >= cfg%forcing%wind_height)                           &
               error stop tag//'forcing.wind_exposure_z0 must lie in (0, wind_height)'
            if (cfg%forcing%wind_blending_height < cfg%forcing%wind_height)                         &
               error stop tag//'forcing.wind_blending_height must be >= wind_height'
         end if
         !----- The terrain lapse rate: between isothermal and about the dry adiabat. --------------!
         if (any(cfg%forcing%lapse_rate_tair < 0.0_wp) .or. any(cfg%forcing%lapse_rate_tair > 0.01_wp)) &
            error stop tag//'site.lapse_rate_tair must lie in [0, 0.01] K/m'
         !----- lwdown_source = "synthesize" is DECLARED but NOT IMPLEMENTED: meds_met_driver      !
         !      always reads LWdown from the file. Accepting the value silently ran the file path    !
         !      under a name that promised Brutsaert/Idso clear-sky synthesis. Reject it until the    !
         !      synthesis exists (MEDS_FORCING_DESIGN.md section 5.7, docs/ROADMAP.md), because a      !
         !      forcing switch that does nothing is worse than one that is absent. -------------------!
         !----- #182: the synthesis is implemented, so the rejection is gone. What is validated     !
         !      instead is that its coefficient is physical -- a negative cloud term would make a    !
         !      cloudy sky emit LESS than a clear one.  --------------------------------------------!
         if (cfg%forcing%lwdown_source == LW_SYNTHESIZE .and. cfg%forcing%lw_cloud_a < 0.0_wp)     &
            error stop tag//'forcing.lw_cloud_a must be >= 0 (a cloudy sky emits MORE, not less)'
         !----- avg_convention: "instant" and "center" PARSE and then run the end-of-interval      !
         !      disaggregation path anyway, because only METAVG_BEGIN has a branch and everything   !
         !      else falls through to METAVG_END. Selecting either therefore got a scheme the user   !
         !      did not ask for, silently. Reject until each has its own branch (#185), on the same  !
         !      precedent as lwdown_source above. ----------------------------------------------------!
         if (cfg%forcing%avg_convention == METAVG_INSTANT)                                       &
            error stop tag//'forcing.avg_convention = "instant" is not implemented (it would '//  &
                            'silently run the "end" path); use "end" or "begin"'
         if (cfg%forcing%avg_convention == METAVG_CENTER)                                        &
            error stop tag//'forcing.avg_convention = "center" is not implemented (it would '//   &
                            'silently run the "end" path); use "end" or "begin"'
         !----- The ED_ERA5land archive (§14-15) is hourly, end-stamped and stores TOTAL shortwave;  !
         !      a config that says otherwise would mis-time or mis-partition every record. ---------!
         if (cfg%forcing%backend == MET_BACKEND_ED_ERA5LAND) then
            if (len_trim(cfg%forcing%data_path) == 0) error stop tag//'forcing.data_path is empty'
            if (cfg%run_mode /= RUN_MODE_REGION .and. cfg%forcing%max_distance_km <= 0.0_wp)      &
               error stop tag//'forcing.max_distance_km must be > 0'
            if (abs(cfg%forcing%dt_forcing - ARCHIVE_DT_SEC) > 0.5_wp)                             &
               error stop tag//'forcing.timestep must be 1 hour for format = "ED_ERA5land"'
            if (cfg%forcing%avg_convention /= METAVG_END)                                          &
               error stop tag//'forcing.avg_convention must be "end" for format = "ED_ERA5land"'
            if (cfg%forcing%sw_partition == SWPART_PASSTHROUGH)                                    &
               error stop tag//'forcing.sw_partition cannot be "passthrough" for format = "ED_ERA5land" '// &
                               '(the archive stores total shortwave)'
            !----- The reader loads an archive month, plus the record before it, before each step     !
            !      (R1, MEDS_POLYGON_RUNTIME_PLAN.md §4). A daily step from midnight reads exactly that; !
            !      a longer step, or one starting mid-day, can straddle two months.  --------------------!
            if (abs(cfg%dt_slow - 86400.0_wp) > 0.5_wp)                                             &
               error stop tag//'format = "ED_ERA5land" needs [run].dt_slow = "1d" (the reader loads a month at a time)'
            if (cfg%start_time%hour /= 0_ik .or. cfg%start_time%minute /= 0_ik .or.                 &
                cfg%start_time%second /= 0_ik)                                                     &
               error stop tag//'format = "ED_ERA5land" needs [run].start_time at 00:00:00'
         end if
         !----- V1 RECYCLE WINDOW: declared, never inferred, and required to be an exact whole      !
         !      number of calendar years. A window of any other length cannot be wrapped without     !
         !      drifting BOTH hour-of-day and day-of-year: the real ERA5-Land Ithaca file spans       !
         !      366 d 22 h, and wrapping on that span shifted the sub-daily shortwave phase by ~10 h   !
         !      and the season by ~48 d over a 29-yr run, while the daily MEAN stayed correct (the     !
         !      cosz reconstruction is mean-conserving) so nothing downstream ever complained. Reject  !
         !      here rather than drift silently. The anchor may sit anywhere in the calendar -- a       !
         !      mid-year window is fine -- only the SPAN must be whole years. -------------------------!
         if (cfg%forcing%recycle) then
            if (.not. time_valid(cfg%forcing%recycle_start))                                    &
               error stop tag//'forcing.recycle_start is not a valid date/time'
            if (.not. time_valid(cfg%forcing%recycle_end))                                      &
               error stop tag//'forcing.recycle_end is not a valid date/time'
            if (.not. time_lt(cfg%forcing%recycle_start, cfg%forcing%recycle_end))              &
               error stop tag//'forcing.recycle_end must be after forcing.recycle_start'
            if (whole_years_between(cfg%forcing%recycle_start, cfg%forcing%recycle_end,         &
                                    MAX_RECYCLE_YEARS) < 1_ik) then
               write(*,'(4a)') ' meds_config: forcing recycle window ',                          &
                  time_to_string(cfg%forcing%recycle_start), ' .. ',                             &
                  time_to_string(cfg%forcing%recycle_end)
               write(*,'(a)')  '   is not an exact whole number of calendar years. Recycling a'
               write(*,'(a)')  '   partial-year window drifts hour-of-day and day-of-year on every'
               write(*,'(a)')  '   wrap. Set recycle_end to the SAME month/day/time as recycle_start,'
               write(*,'(a)')  '   N years later (it is the exclusive end of the window).'
               error stop tag//'forcing recycle window is not a whole number of years'
            end if
         end if
      end if
      !----- Leaf phenology (UNCONDITIONAL now -- no more phenology_on gate; a config's per-PFT cue   !
      !      params are either a deliberate [phenology] override in the PFT file, or the             !
      !      alloc_pft_table literature defaults, so validate them unconditionally either way. The    !
      !      daily-mean temperature that drives the active TEMP/GDD cue comes from the fast loop's     !
      !      met forcing WHEN a calendar context is supplied (advance_leaf_phenology no-ops otherwise, !
      !      leaving the cue drives at their vanilla-evergreen fixed point) -- so there is no fast/     !
      !      forcing precondition to enforce here any more. ALL FIVE cue bits are wired since #150 --   !
      !      TEMP(1), WATER(2), HYDRO(4), PHOTO(8) and LIGHT(16) -- so no bit is rejected any more; the  !
      !      WATER/HYDRO/LIGHT drivers are the root-weighted available-water fraction, the cohort's      !
      !      published predawn leaf potential and the daily-mean incident shortwave. Each mask is a bit  !
      !      set in [0,31]; the rate scales must be non-negative (flush strictly positive).               !
      !                                                                                          !
      !      SELECTABLE IS NOT VALIDATED. The kernel is unit-tested for all four strategies and the      !
      !      four cue drivers each have a hand-computed unit test, but no MEDS run's leaf-area cycle     !
      !      has ever been scored against a phenology observation, under ANY strategy. See the release   !
      !      notes and docs/science/plant_phenology.md.  -------------------------------------------------!
      if (any(cfg%pft%pheno_flush_cue_mask(1:cfg%pft%n) < 0_ik .or.                           &
              cfg%pft%pheno_flush_cue_mask(1:cfg%pft%n) > 31_ik) .or.                         &
          any(cfg%pft%pheno_shed_cue_mask(1:cfg%pft%n) < 0_ik .or.                            &
              cfg%pft%pheno_shed_cue_mask(1:cfg%pft%n) > 31_ik))                              &
         error stop tag//'pheno_{flush,shed}_cue_mask out of range [0,31]'
      !----- CUE_WATER thresholds bracket a FRACTION in [0,1], and flush must sit above shed or the   !
      !      logistic pair is inverted -- the cohort would flush when dry and shed when wet. -----------!
      if (any(cfg%pft%pheno_water_on_threshold(1:cfg%pft%n) <=                                &
              cfg%pft%pheno_water_off_threshold(1:cfg%pft%n)))                                &
         error stop tag//'phenology.water_on_threshold must exceed water_off_threshold'
      if (any(cfg%pft%pheno_water_window(1:cfg%pft%n) <= 0.0_wp) .or.                         &
          any(cfg%pft%pheno_light_window(1:cfg%pft%n) <= 0.0_wp))                             &
         error stop tag//'phenology water/light running-mean windows must be > 0'
      if (any(cfg%pft%pheno_low_psi_threshold(1:cfg%pft%n) <= 0.0_wp) .or.                    &
          any(cfg%pft%pheno_high_psi_threshold(1:cfg%pft%n) <= 0.0_wp))                       &
         error stop tag//'phenology low/high_psi_threshold (days) must be > 0'
      if (any(cfg%pft%pheno_k_flush_max(1:cfg%pft%n) <= 0.0_wp))                              &
         error stop tag//'pheno_k_flush_max must be > 0'
      if (any(cfg%pft%pheno_k_shed_max(1:cfg%pft%n) < 0.0_wp))                                &
         error stop tag//'pheno_k_shed_max must be >= 0'
      if (cfg%cohort_size_tol_min <= 0.0_wp)            error stop tag//'cohort_size_tol_min <= 0'
      if (cfg%cohort_size_tol_max < cfg%cohort_size_tol_min) error stop tag//'cohort_size_tol_max < min'
      if (cfg%n_cohort_fusion_iter < 1_ik)                   error stop tag//'n_cohort_fusion_iter < 1'
      if (cfg%n_patch_fusion_iter < 1_ik)                   error stop tag//'n_patch_fusion_iter < 1'
      if (cfg%patch_light_tol <= 0.0_wp)                    error stop tag//'patch_light_tol <= 0'
      if (cfg%init_soil_temp < 233.0_wp .or. cfg%init_soil_temp > 333.0_wp)                          &
         error stop tag//'init.soil_temp outside 233-333 K'
      if (cfg%init_reacclimate_traits .and. cfg%init_mode /= INIT_RESTART)                          &
         error stop tag//'init.reacclimate_traits applies to a restart (init.init_mode = 2) only'
      if (cfg%patch_light_tol_max < cfg%patch_light_tol)   error stop tag//'patch_light_tol_max < patch_light_tol'
      if (cfg%n_height_layers < 2_ik)                error stop tag//'n_height_layers < 2'
      if (cfg%min_patch_area <= 0.0_wp)              error stop tag//'min_patch_area <= 0'
      if (cfg%cohort_lai_cap <= 0.0_wp)              error stop tag//'cohort_lai_cap <= 0'
      if (cfg%growth_memory_days <= 0.0_wp)          error stop tag//'growth_memory_days <= 0'
      if (cfg%patch_disturbance_rate < 0.0_wp)       error stop tag//'patch_disturbance_rate < 0'
      if (cfg%disturbance_survive_height <= 0.0_wp)  error stop tag//'disturbance_survive_height <= 0'
      if (any(cfg%pft%wood_density <= 0.0_wp))       error stop tag//'wood_density <= 0'
      !----- A recruit must survive its own birth: pool threshold must exceed the cull. ---!
      if (cfg%min_recruit_size <= cfg%negligible_nplant)                                   &
         error stop tag//'min_recruit_size must exceed negligible_nplant'

      !----- THERMAL ACCLIMATION (#176) shifts the PEAKED form's entropy term. The Arrhenius form   !
      !      has no dS, so the flag would be silently inert -- refuse rather than let a config       !
      !      believe acclimation is on when nothing acclimates.  --------------------------------------!
      if (cfg%leaf_thermal_acclimation .and. cfg%temp_response_form /= TRESP_PEAKED)              &
         error stop tag//'leaf_physiology.thermal_acclimation requires temp_response_form = '//    &
                    '"peaked" (the Arrhenius form has no entropy term to shift)'
      if (cfg%leaf_thermal_acclimation .and. cfg%acclim_window_days <= 0.0_wp)                    &
         error stop tag//'leaf_physiology.acclim_window_days must be > 0'

      !----- Leaf physiology: shared biochemistry scalars (Kc/Ko/Gamma* are denominators). -!
      if (cfg%kc25 <= 0.0_wp)            error stop tag//'kc25 <= 0'
      if (cfg%ko25 <= 0.0_wp)            error stop tag//'ko25 <= 0'
      if (cfg%gstar25 <= 0.0_wp)         error stop tag//'gstar25 <= 0'
      if (cfg%o2_mol_frac <= 0.0_wp)     error stop tag//'o2_mol_frac <= 0'
      if (cfg%leaf_absorptance <= 0.0_wp) error stop tag//'leaf_absorptance <= 0'
      if (cfg%phi_psii <= 0.0_wp)        error stop tag//'phi_psii <= 0'
      if (cfg%medlyn_vpd_min <= 0.0_wp)  error stop tag//'medlyn_vpd_min <= 0 [kPa]'
      !----- Leaf physiology: per-PFT traits. ---------------------------------------------!
      if (any(cfg%pft%photosynthetic_pathway /= PATH_C3 .and.                              &
              cfg%pft%photosynthetic_pathway /= PATH_C4)) error stop tag//'photosynthetic_pathway not in {1,2}'
      if (any(cfg%pft%vcmax25 <= 0.0_wp))            error stop tag//'vcmax25 <= 0'
      if (any(cfg%pft%jmax_vcmax_ratio <= 0.0_wp))   error stop tag//'jmax_vcmax_ratio <= 0'
      if (any(cfg%pft%tpu_vcmax_ratio <= 0.0_wp))    error stop tag//'tpu_vcmax_ratio <= 0'
      if (any(cfg%pft%rd_vcmax_ratio < 0.0_wp))      error stop tag//'rd_vcmax_ratio < 0'
      if (any(cfg%pft%stomatal_g0 < 0.0_wp))         error stop tag//'stomatal_g0 < 0'
      if (any(cfg%pft%stomatal_g1 < 0.0_wp))         error stop tag//'stomatal_g1 < 0'
      if (any(cfg%pft%stomatal_d0 <= 0.0_wp))        error stop tag//'stomatal_d0 <= 0 (Leuning divides by it)'
      if (any(cfg%pft%katul_lambda25 <= 0.0_wp))     error stop tag//'katul_lambda25 <= 0'
      if (any(cfg%pft%wstress_lambda_exp < 0.0_wp .or. cfg%pft%wstress_lambda_exp > 8.0_wp)) &
         error stop tag//'wstress_lambda_exp must be in [0,8]'
      if (any(cfg%pft%wstress_psi_open > 0.0_wp))    error stop tag//'wstress_psi_open must be <= 0'
      if (any(cfg%pft%wstress_psi_close >= cfg%pft%wstress_psi_open))                      &
         error stop tag//'wstress_psi_close must be below wstress_psi_open'
      if (any(cfg%pft%wstress_sref_stomata <= 0.0_wp))                                     &
         error stop tag//'wstress_sref_stomata must be > 0 (beta_stomata = exp(sref*(psi - psi_onset)))'
      if (any(cfg%pft%stomata_psi_onset > 0.0_wp))                                         &
         error stop tag//'stomata_psi_onset must be <= 0 (a predawn leaf water potential)'
      if (any(cfg%pft%leaf_surf_water_max < 0.0_wp) .or. any(cfg%pft%wood_surf_water_max < 0.0_wp)) &
         error stop tag//'leaf_surf_water_max and wood_surf_water_max must be >= 0 [kg/m2 of surface]'
      !----- C3 uses theta_j (the J hyperbola / co-limitation curvature); C4 does not. -----!
      if (any(cfg%pft%photosynthetic_pathway == PATH_C3 .and.                              &
              (cfg%pft%theta_j <= 0.0_wp .or. cfg%pft%theta_j >= 1.0_wp)))                 &
         error stop tag//'C3 theta_j must be in (0,1)'
      !----- C4-only constraints (PEPcase slope, light slope, the two co-limitation curvatures).-!
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and. cfg%pft%kp25 <= 0.0_wp))     &
         error stop tag//'C4 PFT needs kp25 > 0'
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and. cfg%pft%quantum_yield_c4 <= 0.0_wp)) &
         error stop tag//'C4 PFT needs quantum_yield_c4 > 0'
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and.                              &
              (cfg%pft%theta_cj_c4 <= 0.0_wp .or. cfg%pft%theta_cj_c4 >= 1.0_wp)))         &
         error stop tag//'C4 theta_cj_c4 must be in (0,1)'
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and.                              &
              (cfg%pft%theta_ic_c4 <= 0.0_wp .or. cfg%pft%theta_ic_c4 >= 1.0_wp)))         &
         error stop tag//'C4 theta_ic_c4 must be in (0,1)'
   end subroutine validate_config

end module meds_config
