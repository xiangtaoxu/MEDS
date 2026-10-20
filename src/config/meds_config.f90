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
   use meds_constants,  only : yr_day, yr_sec, day_sec, n_soil_layer_max
   use meds_pft_params, only : pft_table_t, PATH_C3, PATH_C4, derive_leaf_params, HYD_UNSET
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use meds_allometry_lib, only : HEIGHT_CHAVE
   use meds_time,       only : meds_time_t, time_lt, time_valid, time_to_string,                &
                               whole_years_between
   use meds_therm_lib,     only : TEMP_RESPONSE_ARRHENIUS, TEMP_RESPONSE_PEAKED
   use meds_leaf_opts,     only : SM_LEUNING, SM_MEDLYN, SM_KATUL, COLIM_MIN, COLIM_QUADRATIC
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   use meds_plant_hydr_lib,  only : tissue_psi_tlp
   use meds_soil_column_config, only : soil_column_config_t, soil_texture_layers_t, soil_params_t,       &
                                       soil_column_layers, soil_column_texture, soil_column_params,   &
                                       SOIL_LAYER_TEMPLATE
   use meds_forcing_config, only : forcing_config_t, LW_SYNTHESIZE, METAVG_INSTANT, METAVG_CENTER,   &
                                   METAVG_END, MET_BACKEND_ED_ERA5LAND, SWPART_PASSTHROUGH,        &
                                   WIND_EXPOSURE_OPEN_TERRAIN, ARCHIVE_DT_SEC
   use meds_output_config,  only : output_config_t
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   use meds_biophysics_opts, only : ENERGY_BC_DIRICHLET
   use meds_biogeochem_opts, only : decomp_opts_t
   use meds_region_opts,   only : region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   use meds_maps_opts,  only : maps_opts_t, soil_maps_on
   implicit none
   private

   public :: meds_config_t, site_spec_t, site_config
   public :: derive_config, derive_parameters
   public :: MAX_RECYCLE_YEARS
   public :: validate_config, growth_window_steps
   public :: pft_leaf_psi_tlp, pft_stomata_psi_onset, pft_wstress_psi_open, pft_wstress_psi_close
   public :: forcing_config_t, output_config_t
   public :: decomp_opts_t
   public :: region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   public :: maps_opts_t, soil_maps_on
   public :: BK_SERIAL, BK_MULTICORE, BK_GPU
   public :: DIST_PRIMARY, DIST_TREEFALL
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART
   public :: INTEG_ARK, INTEG_RK45
   public :: CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT, CTRL_I, CTRL_PI

   !----- Parallel backend labels (the actual backend is chosen at COMPILE time via the    !
   !      compiler's do-concurrent target; this is for reporting/reproducibility only).    !
   integer(ik), parameter :: BK_SERIAL    = 0_ik
   integer(ik), parameter :: BK_MULTICORE = 1_ik
   integer(ik), parameter :: BK_GPU       = 2_ik
   !----- Upper bound on the declared forcing recycle window, in whole calendar years. Purely a  !
   !      search bound for the whole-year test (no real met record spans centuries). -------------!
   integer(ik), parameter :: MAX_RECYCLE_YEARS = 200_ik

   !----- Patch disturbance / land-use classes. -------------------------------------------!
   integer(ik), parameter :: DIST_PRIMARY  = 1_ik   !< undisturbed / primary stand
   integer(ik), parameter :: DIST_TREEFALL = 2_ik   !< treefall-gap (age-0) patch
   !----- Initialization modes (selected by [init].init_mode). ----------------------------!
   integer(ik), parameter :: INIT_BARE    = 0_ik    !< near-bare ground
   integer(ik), parameter :: INIT_CENSUS  = 1_ik    !< from a cohort census CSV (init_census_file)
   integer(ik), parameter :: INIT_RESTART = 2_ik    !< restart from a state .nc file (init_restart_file)
   !----- The leaf-model selectors (SM_*, COLIM_*) are owned by meds_leaf_opts and the         !
   !      temperature-response selectors (TEMP_RESPONSE_*) by meds_therm_lib. Imported above   !
   !      for LOADING and VALIDATION only -- not re-exported, so a caller that switches on     !
   !      one imports it from the module that defines it. -------------------------------------!


   !----- Fast-loop TIME integrator ([fast].time_integrator). TWO schemes, and there is             !
   !      deliberately no third: docs/science/numerical_scheme.md §3 ("Why there is no              !
   !      operator-split scheme") records why, and there is no separate coupling-sweep selector     !
   !      because nothing is swept.                                                                 !
   !                                                                                                !
   !      NAMING, stated once so it stops propagating: INTEG_ARK is NOT an IMEX method. The         !
   !      biotic CO2 source is folded implicit, so the explicit tableau is empty (f_E == 0) and     !
   !      the scheme is a 2-solve ESDIRK2 with gamma = 1 - 1/sqrt(2) (the ARS(2,2,2) value;         !
   !      docs/science/numerical_scheme.md eq. (1)). The config string stays "ark". ----------------!
   !----- NOT IMPLEMENTED: a "dynamic vapour pressure" control -- the substomatal air held at       !
   !      the Kelvin humidity exp(psi/(rho_w*Rv*T)) rather than saturated, so the transpiration     !
   !      gradient shrinks with psi and REVERSES into foliar water uptake once e_i < e_a.           !
   !      Scaling the existing latent flux by a humidity factor is not sound, because that flux     !
   !      is linearised about T_cas: it needs the flux AND its temperature derivative rebuilt       !
   !      from rh_leaf*qsat(T_leaf) inside the leaf energy balance. It matters mainly as a          !
   !      route to foliar uptake, not as a stress control: the leaf's water stress is the fast      !
   !      capacity ramp and the slow stomatal factor (docs/science/leaf_gas_exchange.md eqs.        !
   !      (9) and (10)). See docs/ROADMAP.md. ------------------------------------------------------!
   integer(ik), parameter :: INTEG_ARK   = 2_ik  !< ESDIRK2 coupled implicit column (DEFAULT)
   integer(ik), parameter :: INTEG_RK45   = 3_ik  !< adaptive Cash-Karp RK45, the ACCURACY BASELINE

   !----- Fast-loop ERROR-CONTROL selectors (consumed by meds_fast_control;                         !
   !      docs/science/numerical_scheme.md eqs. (3)-(4)). Strictness LEVEL                          !
   !      ([fast].error_level): L0 fixed / L1 adaptive (default) / L2 strict (adaptive + hard       !
   !      failure when a floor step can't meet tolerance). CONTROLLER ([fast].step_controller):     !
   !      I = integral (default) / PI = Gustafsson proportional-integral (damps step-size           !
   !      hunting). --------------------------------------------------------------------------------!
   integer(ik), parameter :: CTRL_L0_FIXED    = 0_ik
   integer(ik), parameter :: CTRL_L1_ADAPTIVE = 1_ik
   integer(ik), parameter :: CTRL_L2_STRICT   = 2_ik
   integer(ik), parameter :: CTRL_I  = 1_ik
   integer(ik), parameter :: CTRL_PI = 2_ik

   !----- A field with no in-type default is set by the config reader (presence-mapped) or     !
   !      derived (derive_config / derive_leaf_params). DERIVED fields are noted. -------------!

   !----- One site of the polygon ([[sites]]; docs/configuration.md, Several sites): its share of !
   !      the polygon, its terrain, the census rows it reads and, when its table sets any of      !
   !      [soil]'s column keys, its own ground, [soil]'s values for each key it leaves out. ------!
   type :: site_spec_t
      real(wp)    :: area_fraction  = 0.0_wp      !< [-] its share of the polygon
      logical     :: has_elevation  = .false.     !< .false.: the polygon's elevation
      real(wp)    :: elevation_m    = 0.0_wp      !< [m]
      real(wp)    :: slope_deg      = 0.0_wp      !< [degrees] from horizontal
      real(wp)    :: aspect_deg     = 0.0_wp      !< [degrees] clockwise from north
      integer(ik) :: census_site_id = 0_ik        !< the census rows it reads; 0: not given
      logical     :: own_ground     = .false.     !< .false.: the polygon's ground
      type(soil_column_config_t) :: soil_column
      integer(ik) :: water_bottom_bc = 0_ik       !< SOIL_BC_*
      real(wp)    :: drain_factor    = 1.0_wp     !< [-]
      real(wp)    :: deep_temp       = 0.0_wp     !< [K]
      real(wp)    :: deep_depth      = 0.0_wp     !< [m]
   end type site_spec_t

   type :: meds_config_t
      !----- Time stepping (run bounded by start/end calendar dates). ----------------------!
      real(wp)          :: dt_slow               !< [s] slow-process timestep (user resolution; default 1 d)
      real(wp)          :: dt_years              !< DERIVED = dt_slow / yr_sec
      type(meds_time_t) :: start_time, end_time
      logical     :: demography_on               !< if .false. structure is frozen
      !----- Master SLOW-TIER freeze ([run].slow_on, DEFAULTED true): .false. skips the WHOLE        !
      !      slow tier (vegetation dynamics -- growth/mortality/phenology/traits/demography --       !
      !      and the soil-carbon step) every step, holding cohort/patch/soil-carbon state static     !
      !      while the fast loop still runs. Broader than demography_on, which only freezes the      !
      !      structural fuse/fiss/disturbance triggers within an otherwise-active slow loop. --------!
      logical     :: slow_on = .true.
      !----- The SLOW-tier conservation ledger ([run].slow_ledger_on, DEFAULTED true). It            !
      !      snapshots the site store either side of every slow-step operator and reports what       !
      !      did not close (docs/science/order_of_processes.md §4). ON by default for the same       !
      !      reason the fast loop's whole-column ledgers are: a conservation check nobody runs       !
      !      is a conservation check nobody has. It writes no state and changes no answer --         !
      !      only an end-of-run table -- so the only cost of leaving it on is a handful of           !
      !      reductions per simulated day. ----------------------------------------------------------!
      logical     :: slow_ledger_on = .true.
      !----- Host THREADS ([run].n_threads, DEFAULTED 1; docs/science/numerical_scheme.md §6a).      !
      !      A site run puts them on the fast loop's PATCH axis; a region run puts them on its       !
      !      POLYGONS, each polygon's patch loop then running on one. Patches, and polygons, are     !
      !      independent within a step, so this is the one lever that costs no accuracy --           !
      !      PROVIDED the answer does not move with the thread count, which is why the site-         !
      !      level reductions are staged per (sub-step, patch) and folded back in patch order.       !
      !      Default 1, so a build that happens to carry OpenMP flags (NVHPC MEDS_GPU=multicore      !
      !      puts -mp PUBLIC on meds_demography, which its dependents inherit) stays serial          !
      !      until asked. Has effect only in an OpenMP build, which is the default; with             !
      !      -DMEDS_OPENMP=OFF the directives are comments and this is ignored. ---------------------!
      integer(ik) :: n_threads = 1_ik
      !----- [run].mode: one site, or every selected cell of a box as its own polygon, with the     !
      !      box and the selection rules in [region] (docs/configuration.md, "Regional runs"). -----!
      integer(ik)         :: run_mode = RUN_MODE_SITE
      type(region_opts_t) :: region
      !----- [maps]: the static maps a run takes its soil and its Chave site index from, at         !
      !      each polygon's cell (meds_surface_maps; docs/configuration.md, "Input maps"). ---------!
      type(maps_opts_t) :: maps
      !----- Fast (sub-daily) biophysics loop. --------------------------------------------!
      logical     :: fast_biophysics_on          !< master gate for the fast biophysics loop
      real(wp)    :: dt_fast                      !< [s] fast biophysics timestep (nested within dt_slow)
      integer(ik) :: n_fast_per_slow              !< DERIVED = max(1, nint(dt_slow / dt_fast))
      real(wp)    :: snow_init_swe        = 0.0_wp     !< [kg/m2] initial snow water-equivalent seeded at run start
      real(wp)    :: snow_init_temp       = 270.0_wp   !< [K] initial snow temperature (for the seeded pack)
      logical     :: canopy_water_on      = .false.    !< opt-in canopy interception film + film-evap/dew
      !----- Fast-loop TIME integrator selector + ARK knobs ([fast], DEFAULTED reads;               !
      !      docs/science/numerical_scheme.md §3, §6). ---------------------------------------------!
      integer(ik) :: time_integrator      = INTEG_ARK !< INTEG_ARK (default) | INTEG_RK45
      logical     :: ark_adaptive         = .true.      !< adaptive (embedded-error) vs fixed-substep march
      real(wp)    :: ark_rtol             = 1.0e-3_wp   !< adaptive relative tolerance (broadcast to all tol groups)
      !----- ONE master relative-accuracy dial for the WHOLE fast loop: when > 0 it overrides        !
      !      every tolerance group's rtol -- the ARK march AND the nested soil-water / soil-         !
      !      energy / plant-hydraulics sub-solvers. 0 (default) => each keeps its own per-sub-       !
      !      solver value. --------------------------------------------------------------------------!
      real(wp)    :: rtol_all             = 0.0_wp      !< [-] 0 = unset; > 0 = the single accuracy target
      !----- The ABSOLUTE-tolerance companion to rtol_all. Every group's atol is physically scaled     !
      !      (J/kg, kg/kg, m3/m3, K, ...), so they cannot share one number the way rtol can;           !
      !      instead this is a dimensionless MULTIPLIER applied to the whole atol vector. It           !
      !      matters because the WRMS denominator is atol + rtol*|y|                                   !
      !      (docs/science/numerical_scheme.md eq. (3)): tightening rtol_all ALONE saturates once      !
      !      atol dominates the denominator (measured: 1e-3 -> 1e-6 on the soil-water group only       !
      !      raises the error estimate ~4x, too little to force a substep), so the accuracy dial       !
      !      is only half-effective without it. 1.0 (default) leaves every atol exactly as set. -------!
      real(wp)    :: atol_scale           = 1.0_wp      !< [-] multiplier on every group's atol (1 = unset)
      !----- WHERE IN THE SUB-STEP THE MET FORCING IS SAMPLED, as a fraction of dt_fast. The         !
      !      fast loop samples met at (isub - 1 + f)*dt_fast while column_prepass freezes every      !
      !      coefficient (gs/GPP/Rd/aerodynamics/radiation) on the STATE at t^n. f = 0.5 (the        !
      !      default) is the better quadrature of the forcing alone, but it pairs t+dt/2 forcing     !
      !      with t^n state -- a FIRST-ORDER inconsistency in the frozen coefficients; f = 0 is      !
      !      the state-consistent choice, and a true midpoint freeze (f = 0.5 WITH a state           !
      !      predictor) is the second-order version. It is a trade: f = 0 improves canopy-air        !
      !      temperature and degrades soil temperature (docs/science/numerical_scheme.md §7,         !
      !      item 5). -------------------------------------------------------------------------------!
      real(wp)    :: forcing_sample_frac  = 0.5_wp      !< [-] met sample point within the sub-step, in [0,1]
      !----- The CAS supersaturation (condensation) sink in surface_derivs, which both schemes      !
      !      reach. Default .true.; .false. drops the sink, for an experiment that wants the        !
      !      canopy air allowed above saturation. --------------------------------------------------!
      logical     :: cas_condensation     = .true.     !< apply the CAS supersaturation sink (ARK path)
      !----- PROCESS MASK: which column processes actually EVOLVE. All true = the full column        !
      !      (default); flipping one off freezes that store so the driver integrates a REDUCED       !
      !      ODE, for numerical studies of process complexity. The mask type itself lives in         !
      !      meds_fast_types (with column_config_t); config carries plain logicals so config/        !
      !      does not gain a driver dependency. -----------------------------------------------------!
      logical     :: mask_veg_energy = .true.   !< leaf + wood energy stores
      logical     :: mask_cas_energy = .true.   !< canopy-air-space enthalpy
      logical     :: mask_cas_vapour = .true.   !< canopy-air-space specific humidity
      logical     :: mask_cas_co2    = .true.   !< canopy-air-space CO2
      logical     :: mask_soil_heat  = .true.   !< soil thermal column
      logical     :: mask_soil_water = .true.   !< soil water column
      logical     :: mask_hydraulics = .true.   !< plant hydraulics (psi)
      integer(ik) :: step_controller      = CTRL_I      !< CTRL_I (default) | CTRL_PI
      integer(ik) :: error_level          = CTRL_L1_ADAPTIVE !< CTRL_L0_FIXED | CTRL_L1_ADAPTIVE (default) | CTRL_L2_STRICT
      real(wp)    :: ark_dt_init          = 0.0_wp      !< [s] initial adaptive substep (0 => dt_fast)
      integer(ik) :: ark_fixed_substep    = 4_ik        !< fixed substeps/dt_fast (GPU warp-uniform path)
      !----- The leaf<->CAS surface solve is EITHER the uncoupled single-BE pass OR the coupled     !
      !      2x2 Newton (its own iteration cap is NEWT_MAX) -- there is nothing in between, so      !
      !      [fast].ark_coupled is a boolean. ------------------------------------------------------!
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
      real(wp)    :: light_ext_coeff_fusefiss = 0.5_wp   !< [-] extinction of the light profiles fusion compares
      real(wp)    :: min_patch_area, patch_min_area_remain
      logical     :: enable_patch_fission

      !----- Patch disturbance, growth memory, recruitment, conservation. -----------------!
      real(wp) :: patch_disturbance_rate, disturbance_survive_height
      real(wp) :: growth_memory_days, min_recruit_size, conservation_tol

      !----- Initial conditions (init_mode: 0 bare | 1 census | 2 restart). ---------------!
      integer(ik)        :: init_mode
      character(len=256) :: init_restart_file, init_census_file
      !----- The soil state a run starts from when no state file restores one: every layer of    !
      !      every patch at this temperature and volumetric water content. ----------------------!
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

      !----- Leaf physiology: model selection (non-PFT). ----------------------------------!
      integer(ik) :: stomatal_model           !< SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form        !< TEMP_RESPONSE_ARRHENIUS | TEMP_RESPONSE_PEAKED
      integer(ik) :: colimitation             !< COLIM_MIN | COLIM_QUADRATIC
      logical     :: leaf_use_boundary_layer  !< if .true., couple via the leaf boundary layer (gb)
      !----- Leaf physiology: shared biochemistry at 25 degC + Arrhenius/deactivation terms.-!
      real(wp) :: kc25, ko25, gstar25                   !< [Pa]    Michaelis constants + CO2 compensation point
      real(wp) :: ea_kc, ea_ko, ea_gstar                !< [J/mol] activation energies (Bernacchi et al. 2001)
      real(wp) :: ea_vcmax, ea_jmax, ea_rd              !< [J/mol] activation energies
      real(wp) :: hd_vcmax, hd_jmax, hd_rd              !< [J/mol] deactivation energies (peaked form)
      real(wp) :: ds_vcmax, ds_jmax, ds_rd              !< [J/mol/K] entropy terms (peaked form)
      !----- THERMAL ACCLIMATION (Kattge & Knorr 2007; docs/science/leaf_gas_exchange.md eq.         !
      !      (4)). OFF by default: `ds_vcmax`/`ds_jmax`/`jmax_vcmax_ratio` stay the fixed values     !
      !      above. ON, the driver recomputes those three from a running-mean GROWTH temperature     !
      !      once per slow step, so a warm-grown and a cold-grown stand of the same PFT do not       !
      !      share one temperature response. Only meaningful with temp_response = peaked: the        !
      !      Arrhenius form has no dS to shift, and validate_config says so rather than letting      !
      !      the flag be silently inert. ------------------------------------------------------------!
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

      !----- Light trait plasticity ([trait_dynamics]; docs/science/plant_traits.md). OPT-IN:       !
      !      default .false. => cohort leaf traits (sla/vcmax25/rd25/leaf_turnover_rate) stay at   !
      !      their top-of-canopy PFT values. When ON, the slow-loop driver acclimates them to       !
      !      cumulative LAI above. -----------------------------------------------------------------!
      logical     :: trait_plasticity_on = .false.
      !----- Fine roots turn over with the leaves they feed: with plasticity on, a cohort's        !
      !       fine-root turnover is its PFT's times (its own leaf turnover rate / the top-of-canopy   !
      !       one), so a shaded cohort replaces its roots as slowly as its leaves. OPT-IN: default    !
      !       .false. docs/science/plant_traits.md §4. ------------------------------------------!
      logical     :: fineroot_turnover_plastic = .false.

      !----- Meteorological forcing ([forcing]/[site]; docs/science/forcing.md). OPT-IN:           !
      !      forcing_on default .false. (the whole [forcing] block is gated on it), so a           !
      !      config with no [forcing] block runs on constant forcing, the met_forcing_t            !
      !      reference climate. Defaults are the Ithaca NY / ERA5-Land reference. -----------------!
      type(forcing_config_t) :: forcing

      !----- Diagnostic-aggregation output ([output]). OPT-IN: enabled default .false. (a           !
      !      config with no [output] block emits no diagnostic stream). The per-variable            !
      !      overrides live in the optional meds_io_config.toml named by output%io_config           !
      !      (docs/science/diagnostics.md §5). -----------------------------------------------------!
      type(output_config_t) :: output

      !----- PFT traits. ------------------------------------------------------------------!
      type(pft_table_t) :: pft

      !----- The physical soil column ([soil]). ------------------------------------------------!
      type(soil_column_config_t) :: soil_column  !< [soil] the column's layers and texture
      !----- The polygon's sites ([[sites]]); not allocated: one site, the polygon, on [soil]. ---!
      type(site_spec_t), allocatable :: sites(:)

      !----- Fast-loop biophysics run-config ([soil]/[snow]/[aerodynamics] and the solver settings  !
      !      in [fast], every key optional, defaults in meds_biophysics_opts). build_fast_context   !
      !      copies each verbatim into the column config (col_config%soil_water_opts/energy/snow/   !
      !      aero). ---------------------------------------------------------------------------------!
      type(soil_opts_t)   :: soil        !< [soil]+[fast]  soil water and water-network opts (-> col_config%soil_water_opts)
      type(energy_opts_t) :: energy      !< [soil]+[fast]  soil-thermal solver opts       (-> col_config%energy)
      type(snow_params_t) :: snow        !< [snow]         snow physical parameter table  (-> col_config%snow)
      type(aero_cfg_t)    :: aero        !< [aerodynamics] canopy-aerodynamics constants  (-> col_config%aero)

      !----- Slow soil-carbon matrix ([soil_carbon], ON by default;                                  !
      !      docs/science/soil_carbon.md). soil_carbon_on gates the FEATURE. It does NOT gate        !
      !      whether the fields are required: every key is a DEFAULTED read (like [snow]),           !
      !      falling back to its ED2 in-type default, so the feature needs no TOML edits beyond      !
      !      this flag.                                                                              !
      !                                                                                              !
      !      Off is not a coarser soil-carbon model, it is NO soil carbon: litter is discarded       !
      !      at the slow step and patch_heterotrophic_respiration returns rh = 0, so Reco            !
      !      carries only its autotrophic limb. Measured over a year at Ithaca, off reports          !
      !      annual-mean NEE at -4.657 against -2.464 umol/m2/s -- an 89% stronger apparent sink     !
      !      -- with Rh identically zero against 0.833 kgC/m2/yr, for no wall-clock saving. ---------!
      logical            :: soil_carbon_on = .true.
      type(decomp_opts_t) :: soil_carbon   !< [soil_carbon] decomposition selectors + rate parameters
      !----- Cold-start spin-up ([soil_carbon], consumed only when soil_carbon_on): zero-init        !
      !      (default) leaves every pool at 0, matching bare-ground philosophy; steady-state         !
      !      solves SASU (solve_soil_carbon_steady_state; docs/science/soil_carbon.md eq. (7))       !
      !      from a scalar climatological environmental factor and a constant litter-input           !
      !      estimate (pools 5-7 get zero input, as in build_litter_input -- only 1-4 ever           !
      !      receive litter). -----------------------------------------------------------------------!
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
   ! height-layer edges (0 to the tallest PFT's hgt_max).                                       !
   !---------------------------------------------------------------------------------------!
   subroutine derive_config(cfg)
      type(meds_config_t), intent(inout) :: cfg
      integer(ik) :: i

      !----- Slow-process timestep in years (demography currency). ---------------------!
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
         cfg%height_edges(i) = real(i, wp) * maxval(cfg%pft%allom(1:cfg%pft%n)%hgt_max)          &
                               / real(cfg%n_height_layers, wp)
      end do
   end subroutine derive_config

   !---------------------------------------------------------------------------------------!
   ! Derive EVERY parameter that is a function of the primary (loaded) config: the derived run   !
   ! scalars (derive_config) and the leaf capacities from their ratios to Vcmax25              !
   ! (derive_leaf_params). The single consolidation point for both the production loader       !
   ! (load_meds_config) and the test builder (build_test_config), so the derivation sequence   !
   ! lives in ONE place. Callers validate AFTER this returns.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine derive_parameters(cfg)
      type(meds_config_t), intent(inout) :: cfg
      call derive_config(cfg)
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

   !----- A PFT's leaf turgor-loss point [MPa], from its own leaf pressure-volume curve, so the      !
   !      stomatal closure and the phenology's drought counter see the curve the solver stores leaf  !
   !      water on. -----------------------------------------------------------------------------------!
   pure real(wp) function pft_leaf_psi_tlp(cfg, ipft) result(psi_tlp)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      associate (leaf_curve => cfg%pft%hydro(ipft)%leaf_curve)
         psi_tlp = tissue_psi_tlp(leaf_curve%pi0, leaf_curve%elastic_mod)
      end associate
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

   !----- The current leaf potential [MPa] at which the fast water-stress limb (beta_nonstomata,  !
   !      on Vcmax/Jmax/TPU) starts and completes: the [pft] wstress_psi_open / wstress_psi_close   !
   !      where given, else the PFT's turgor-loss point and twice it. A leaf keeps its capacity    !
   !      while it holds turgor and loses it as it wilts. ------------------------------------------!
   pure real(wp) function pft_wstress_psi_open(cfg, ipft) result(psi_open)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      psi_open = pft_leaf_psi_tlp(cfg, ipft)
      if (allocated(cfg%pft%wstress_psi_open)) then
         if (cfg%pft%wstress_psi_open(ipft) > HYD_UNSET) psi_open = cfg%pft%wstress_psi_open(ipft)
      end if
   end function pft_wstress_psi_open

   pure real(wp) function pft_wstress_psi_close(cfg, ipft) result(psi_close)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      psi_close = 2.0_wp * pft_leaf_psi_tlp(cfg, ipft)
      if (allocated(cfg%pft%wstress_psi_close)) then
         if (cfg%pft%wstress_psi_close(ipft) > HYD_UNSET) psi_close = cfg%pft%wstress_psi_close(ipft)
      end if
   end function pft_wstress_psi_close

   !----- Site `is`'s configuration: the polygon's, with the site's ground when its table gives  !
   !      one. make_site_surface builds the site's surface from it. -----------------------------!
   function site_config(cfg, is) result(scfg)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: is
      type(meds_config_t)             :: scfg
      scfg = cfg
      if (.not. cfg%sites(is)%own_ground) return
      scfg%soil_column       = cfg%sites(is)%soil_column
      scfg%soil%bottom_bc    = cfg%sites(is)%water_bottom_bc
      scfg%soil%drain_factor = cfg%sites(is)%drain_factor
      scfg%energy%deep_temp  = cfg%sites(is)%deep_temp
      scfg%energy%deep_depth = cfg%sites(is)%deep_depth
   end function site_config

   !---------------------------------------------------------------------------------------!
   ! Validate a configuration; halt on a setting that would corrupt the run. Each site     !
   ! with its own ground is validated as the run with that ground.                         !
   !---------------------------------------------------------------------------------------!
   recursive subroutine validate_config(cfg)
      type(meds_config_t), intent(in) :: cfg
      character(len=*), parameter :: TAG = 'meds_config: '
      integer(ik) :: ipft

      !----- [[sites]]: shares of the polygon that sum to 1, a terrain that is a surface, the      !
      !      census rows of each site in a census run, and a valid ground for each site that has   !
      !      its own. A region takes one site per cell until a site map gives its sites. ----------!
      if (allocated(cfg%sites)) then
         sites: block
            type(meds_config_t) :: scfg
            integer(ik) :: is, ns
            ns = size(cfg%sites, kind=ik)
            if (cfg%run_mode == RUN_MODE_REGION)                                                     &
               error stop TAG//'a region refuses [[sites]]: each cell is one site until a site map gives its sites'
            if (any(cfg%sites%area_fraction <= 0.0_wp)) error stop TAG//'[[sites]].area_fraction must be > 0'
            if (abs(sum(cfg%sites%area_fraction) - 1.0_wp) > 1.0e-6_wp)                             &
               error stop TAG//'the [[sites]] area_fraction values must sum to 1'
            if (any(cfg%sites%slope_deg < 0.0_wp .or. cfg%sites%slope_deg >= 90.0_wp))              &
               error stop TAG//'[[sites]].slope must lie in [0, 90) degrees'
            if (any(cfg%sites%aspect_deg < 0.0_wp .or. cfg%sites%aspect_deg > 360.0_wp))            &
               error stop TAG//'[[sites]].aspect must lie in [0, 360] degrees'
            if (cfg%init_mode == INIT_CENSUS) then
               if (any(cfg%sites%census_site_id <= 0_ik))                                           &
                  error stop TAG//'a census run needs [[sites]].census_site_id > 0 for every site'
               do is = 2_ik, ns
                  if (any(cfg%sites(1:is - 1_ik)%census_site_id == cfg%sites(is)%census_site_id))   &
                     error stop TAG//'two [[sites]] read the same census_site_id'
               end do
            end if
            do is = 1_ik, ns
               if (.not. cfg%sites(is)%own_ground) cycle
               write(*,'(a,i0,a)') ' config: checking the ground of [[sites]] ', is, ''
               scfg = site_config(cfg, is)
               deallocate(scfg%sites)
               call validate_config(scfg)
            end do
         end block sites
      end if

      !----- [hydraulics]: the rooting traits set the root profile over the soil column, and the    !
      !      sapwood conductance needs a positive conductivity and path factor. ---------------------!
      associate (hydro => cfg%pft%hydro(1:cfg%pft%n))
         if (any(hydro%root_beta <= 0.0_wp .or. hydro%root_beta >= 1.0_wp))                         &
            error stop TAG//'hydraulics.root_beta must lie in (0, 1)'
         if (any(hydro%root_depth_b1 <= 0.0_wp)) error stop TAG//'hydraulics.root_depth_b1 must be > 0'
         if (any(hydro%root_depth_b2 < 0.0_wp))  error stop TAG//'hydraulics.root_depth_b2 must be >= 0'
         if (any(hydro%specific_root_length <= 0.0_wp))                                             &
            error stop TAG//'hydraulics.specific_root_length must be > 0'
         if (any(hydro%fine_root_radius <= 0.0_wp)) error stop TAG//'hydraulics.fine_root_radius must be > 0'
         if (any(hydro%vessel_curl <= 0.0_wp)) error stop TAG//'hydraulics.vessel_curl must be > 0'
         if (any(hydro%wood_kmax <= 0.0_wp))   error stop TAG//'hydraulics.wood_kmax must be > 0'
      end associate
      !----- [soil] ground optics. --------------------------------------------------------------!
      if (cfg%soil%ground_albedo_vis < 0.0_wp .or. cfg%soil%ground_albedo_vis >= 1.0_wp .or.       &
          cfg%soil%ground_albedo_nir < 0.0_wp .or. cfg%soil%ground_albedo_nir >= 1.0_wp)           &
         error stop TAG//'soil.ground_albedo_vis and ground_albedo_nir must lie in [0, 1)'
      if (cfg%soil%ground_emissivity <= 0.0_wp .or. cfg%soil%ground_emissivity > 1.0_wp)           &
         error stop TAG//'soil.ground_emissivity must lie in (0, 1]'

      !----- The [soil] column. Every one of these produces a silently WRONG column rather than a !
      !      crash: a layer count over the compile-time ceiling writes past the active region, a   !
      !      non-positive depth or a theta_sat <= theta_res divides by zero in the retention curve, !
      !      and a van Genuchten n <= 1 makes the curve's exponent negative. Fail loud instead.    !
      !      A texture per layer (from the soil maps) is checked layer by layer. ------------------!
      soil_column: block
      type(soil_texture_layers_t) :: tex
      real(wp)    :: bottoms(n_soil_layer_max)
      integer(ik) :: n
      associate (sc => cfg%soil_column)
         if (sc%depth <= 0.0_wp)          error stop TAG//'soil.depth <= 0'
         if (sc%n_layer_bottoms > 0_ik) then
            if (sc%layer_bottoms(1) <= 0.0_wp .or. any(sc%layer_bottoms(2:sc%n_layer_bottoms)              &
                <= sc%layer_bottoms(1:sc%n_layer_bottoms - 1)))                                     &
               error stop TAG//'soil.layer_bottoms must be positive depths [m] that increase'
            if (sc%layer_bottoms(sc%n_layer_bottoms) < sc%depth)                                    &
               error stop TAG//'soil.layer_bottoms must reach soil.depth'
         else if (sc%depth > SOIL_LAYER_TEMPLATE(size(SOIL_LAYER_TEMPLATE))) then
            error stop TAG//'soil.depth is below the deepest default layer (27 m): '//       &
                       'set soil.layer_bottoms'
         end if
         call soil_column_layers(sc, bottoms, n)
         if (sc%n_layer_values > 0_ik .and. sc%n_layer_values /= n)                                 &
            error stop TAG//'soil: the per-layer texture does not match the column''s layers'
         call soil_column_texture(sc, tex)
         if (any(tex%theta_sat(1:n) <= tex%theta_res(1:n))) error stop TAG//'soil.theta_sat <= theta_res'
         if (any(tex%theta_res(1:n) < 0.0_wp))       error stop TAG//'soil.theta_res < 0'
         !----- With a texture per layer, the run starts each layer at its field capacity, so      !
         !      init.soil_theta is not used (build_fast_context). ----------------------------------!
         if (sc%n_layer_values == 0_ik) then
            if (cfg%init_soil_theta <= sc%theta_res .or. cfg%init_soil_theta > sc%theta_sat)          &
               error stop TAG//'init.soil_theta must lie in (soil.theta_res, theta_sat]'
         end if
         if (any(tex%ksat(1:n) <= 0.0_wp))           error stop TAG//'soil.ksat <= 0'
         !----- curve_par_a and curve_par_n mean different things per family: van Genuchten's alpha  !
         !      [1/m] > 0 and n > 1, Campbell's air-entry suction psi_sat [m] < 0 and exponent b > 0.  !
         !      The other family's pair makes the curve raise a negative base to a fractional power. -!
         if (sc%retention == SOIL_RETENTION_VG) then
            if (any(tex%par_a(1:n) <= 0.0_wp))                                                    &
               error stop TAG//'soil.curve_par_a (van Genuchten alpha, 1/m) must be > 0'
            if (any(tex%par_n(1:n) <= 1.0_wp))                                                    &
               error stop TAG//'soil.curve_par_n must exceed 1 for van Genuchten'
         else
            if (any(tex%par_a(1:n) >= 0.0_wp))                                                    &
               error stop TAG//'soil.curve_par_a (Campbell psi_sat, m) must be < 0'
            if (any(tex%par_n(1:n) <= 0.0_wp))                                                    &
               error stop TAG//'soil.curve_par_n (Campbell b) must be > 0'
         end if
         if (sc%psi_fc >= 0.0_wp)         error stop TAG//'soil.psi_fc must be negative (a suction head)'
         if (any(tex%k_solid(1:n) <= 0.0_wp) .or. any(tex%k_dry(1:n) <= 0.0_wp))                  &
            error stop TAG//'soil conductivities must be positive'
         if (any(tex%dry_cvol(1:n) <= 0.0_wp)) error stop TAG//'soil.dry_heat_capacity <= 0'
      end associate
      end block soil_column

      !----- [soil] bottom thermal BC (docs/science/soil_biophysics.md, "The bottom thermal         !
      !      boundary"). The Dirichlet anchor needs a SITE temperature and a depth BELOW the        !
      !      bottom node; neither can be guessed. `deep_temp` is the mean annual soil               !
      !      temperature below the damping depth (close to mean annual AIR temperature), and an     !
      !      error in it is a steady flux into the column base, so a silent default would be a      !
      !      silent mean-annual bias -- the very bias the anchor exists to remove. `deep_depth`     !
      !      at or above the bottom NODE gives a non-positive conduction length, i.e. a             !
      !      division by zero or a sign flip that would pump heat the wrong way. With the soil      !
      !      maps, the anchor is each polygon's (apply_map_inputs), so a configuration the maps     !
      !      have not filled yet has none to check. ------------------------------------------------!
      if (cfg%energy%bottom_bc == ENERGY_BC_DIRICHLET .and.                                         &
          .not. (soil_maps_on(cfg%maps) .and. cfg%soil_column%n_layer_values == 0_ik)) then
         block
            type(soil_params_t) :: sp
            if (cfg%energy%deep_temp <= 0.0_wp)                                                  &
               error stop TAG//'soil.thermal_bottom_bc = "dirichlet" requires soil.deep_temp [K] -- '// &
                          'the mean annual soil temperature below the damping depth. It is a site '//&
                          'property, like latitude; there is no safe default.'
            if (cfg%energy%deep_temp < 200.0_wp .or. cfg%energy%deep_temp > 330.0_wp)             &
               error stop TAG//'soil.deep_temp outside 200..330 K -- it is an absolute temperature'
            !----- Build the real grid rather than re-deriving the bottom node here: the grid rule  !
            !      lives in ONE place and a second copy would drift from it silently. --------------!
            call soil_column_params(cfg%soil_column, sp)
            if (cfg%energy%deep_depth <= abs(sp%z_node(sp%n_active)))                            &
               error stop TAG//'soil.deep_depth must lie BELOW the bottom soil node '//          &
                          '(it is measured down from the surface, and the anchor conducts to it)'
         end block
      end if

      if (any(cfg%pft%nbg_init_density < 0.0_wp)) error stop TAG//'pft.nbg_init_density must be >= 0'

      !----- Co-limitation and hyperbola curvatures must lie in (0, 1]. A curvature of zero          !
      !      makes the smoothing quadratic degenerate and a value above one has no root in the       !
      !      physical branch, so either produces a silently wrong assimilation rather than a         !
      !      crash; an unset curvature in a test fixture reaches the solver as zero. ----------------!
      block
         integer(ik) :: pf
         do pf = 1_ik, cfg%pft%n
            if (cfg%pft%theta_j(pf)     <= 0.0_wp .or. cfg%pft%theta_j(pf)     > 1.0_wp)         &
               error stop TAG//'pft.theta_j outside (0, 1]'
            if (cfg%pft%theta_cj_c3(pf) <= 0.0_wp .or. cfg%pft%theta_cj_c3(pf) > 1.0_wp)         &
               error stop TAG//'pft.theta_cj_c3 outside (0, 1]'
            if (cfg%pft%theta_ip_c3(pf) <= 0.0_wp .or. cfg%pft%theta_ip_c3(pf) > 1.0_wp)         &
               error stop TAG//'pft.theta_ip_c3 outside (0, 1]'
            if (cfg%pft%theta_cj_c4(pf) <= 0.0_wp .or. cfg%pft%theta_cj_c4(pf) > 1.0_wp)         &
               error stop TAG//'pft.theta_cj_c4 outside (0, 1]'
            if (cfg%pft%theta_ic_c4(pf) <= 0.0_wp .or. cfg%pft%theta_ic_c4(pf) > 1.0_wp)         &
               error stop TAG//'pft.theta_ic_c4 outside (0, 1]'
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
               error stop TAG//'pft.leaf_reflect_vis + leaf_transmit_vis >= 1'
            if (cfg%pft%leaf_reflect_nir(pf) + cfg%pft%leaf_transmit_nir(pf) >= 1.0_wp)           &
               error stop TAG//'pft.leaf_reflect_nir + leaf_transmit_nir >= 1'
            if (cfg%pft%wood_reflect_vis(pf) + cfg%pft%wood_transmit_vis(pf) >= 1.0_wp)           &
               error stop TAG//'pft.wood_reflect_vis + wood_transmit_vis >= 1'
            if (cfg%pft%wood_reflect_nir(pf) + cfg%pft%wood_transmit_nir(pf) >= 1.0_wp)           &
               error stop TAG//'pft.wood_reflect_nir + wood_transmit_nir >= 1'
            if (cfg%pft%leaf_emissivity(pf) <= 0.0_wp .or. cfg%pft%leaf_emissivity(pf) > 1.0_wp)  &
               error stop TAG//'pft.leaf_emissivity outside (0,1]'
            if (cfg%pft%wood_emissivity(pf) <= 0.0_wp .or. cfg%pft%wood_emissivity(pf) > 1.0_wp)  &
               error stop TAG//'pft.wood_emissivity outside (0,1]'
            if (cfg%pft%leaf_clumping(pf) <= 0.0_wp .or. cfg%pft%leaf_clumping(pf) > 1.0_wp)      &
               error stop TAG//'pft.leaf_clumping outside (0,1]'
            if (cfg%pft%wood_clumping(pf) <= 0.0_wp .or. cfg%pft%wood_clumping(pf) > 1.0_wp)      &
               error stop TAG//'pft.wood_clumping outside (0,1]'
            if (cfg%pft%leaf_angle_mean(pf) <= 0.0_wp .or. cfg%pft%leaf_angle_mean(pf) >= 90.0_wp) &
               error stop TAG//'pft.leaf_angle_mean outside (0,90) degrees'
            if (cfg%pft%leaf_angle_std(pf) < 0.0_wp)                                              &
               error stop TAG//'pft.leaf_angle_std < 0'
         end do
      end block

      if (cfg%pft%n < 1_ik)                          error stop TAG//'empty PFT table'
      if (.not. time_lt(cfg%start_time, cfg%end_time)) error stop TAG//'end_time must be after start_time'
      if (cfg%dt_slow <= 0.0_wp)                        error stop TAG//'dt_slow <= 0'
      if (cfg%fast_biophysics_on) then
         if (cfg%dt_fast <= 0.0_wp)                     error stop TAG//'dt_fast <= 0'
         !----- dt_fast is an ACCURACY parameter, not a stability one                               !
         !      (docs/science/numerical_scheme.md §5a): the surface layer is re-solved at every     !
         !      stage, so the canopy air is stable at every dt_fast measured, and what a longer     !
         !      step costs is in the carbon and water budget, through the leaf gas exchange         !
         !      frozen at the step start. The default 900 s suits production, including spin-       !
         !      ups; sub-daily fidelity (flux-tower comparison, diel cycles, energy                 !
         !      partitioning) wants <= 150 s (§6). Above 900 s nothing has been measured, hence     !
         !      the warnings. ----------------------------------------------------------------------!
         if (cfg%dt_fast > 900.0_wp) then
            print '(a)', 'WARNING [meds_config]: dt_fast > 900 s is beyond the measured range. At'
            print '(a)', '  900 s the canopy air is stable and carbon and daily leaf water potential'
            print '(a)', '  converge; above it neither has been measured.'
         end if
         if (cfg%dt_fast > 1800.0_wp) then
            print '(a)', 'WARNING [meds_config]: dt_fast > 1800 s is outside the measured range entirely.'
         end if
         if (cfg%dt_fast > cfg%dt_slow)                 error stop TAG//'dt_fast > dt_slow'
         if (abs(cfg%dt_slow / cfg%dt_fast - real(nint(cfg%dt_slow / cfg%dt_fast, ik), wp)) > 1.0e-6_wp) &
            error stop TAG//'dt_slow must be an integer multiple of dt_fast'
         if (cfg%time_integrator /= INTEG_ARK .and. cfg%time_integrator /= INTEG_RK45)          &
            error stop TAG//'time_integrator out of range'
         if (cfg%rtol_all < 0.0_wp)           error stop TAG//'rtol_all < 0 (0 = unset)'
         if (cfg%atol_scale <= 0.0_wp)        error stop TAG//'atol_scale <= 0 (1 = unset)'
         if (cfg%forcing_sample_frac < 0.0_wp .or. cfg%forcing_sample_frac > 1.0_wp)                 &
            error stop TAG//'forcing_sample_frac outside [0,1]'
         if (cfg%step_controller /= CTRL_I .and. cfg%step_controller /= CTRL_PI)                &
            error stop TAG//'step_controller out of range (I|PI)'
         if (cfg%error_level /= CTRL_L0_FIXED .and. cfg%error_level /= CTRL_L1_ADAPTIVE .and.    &
             cfg%error_level /= CTRL_L2_STRICT)                                                 &
            error stop TAG//'error_level out of range (L0|L1|L2)'
         if (cfg%ark_rtol <= 0.0_wp)          error stop TAG//'ark_rtol <= 0'
         if (cfg%ark_fixed_substep < 1_ik)    error stop TAG//'ark_fixed_substep < 1'
         !----- The sub-daily probe is the ONE fast-loop writer that is not thread-safe -- it        !
         !      holds a `save`d unit and appends to a shared file, so under threads its rows         !
         !      would interleave in thread-arrival order and stop being a reproducible               !
         !      diagnostic. Serializing the write would fix the corruption but not the ORDER,        !
         !      which is the property the probe exists for, so the two are mutually exclusive by     !
         !      construction rather than by luck (docs/science/numerical_scheme.md §6a). ------------!
         if (cfg%fast_probe .and. cfg%n_threads > 1_ik)                                            &
            error stop TAG//'fast_probe requires n_threads = 1 (the probe writes one shared, &
                           &order-significant CSV)'
      end if
      if (cfg%n_threads < 1_ik)               error stop TAG//'n_threads < 1'
      !----- The FAST output tier closes a window every fast_interval_steps sub-steps, counted       !
      !      over the run. A window must not straddle two slow steps: the stand can be               !
      !      restructured between them, and a window's cohort and patch slots are fixed when it      !
      !      opens. ---------------------------------------------------------------------------------!
      if (cfg%fast_biophysics_on .and. cfg%output%enabled) then
         if (cfg%output%fast_interval_steps < 1_ik)                                                &
            error stop TAG//'output.fast_interval_steps < 1'
         if (mod(cfg%n_fast_per_slow, cfg%output%fast_interval_steps) /= 0_ik)                    &
            error stop TAG//'output.fast_interval_steps must divide the fast steps in a slow step '// &
                            '(dt_slow / dt_fast), so no fast-tier window straddles two slow steps'
      end if
      !----- REGION MODE (docs/configuration.md, "Regional runs"). Every polygon reads its own       !
      !      cell of the ED_ERA5land archive, so a region needs the archive, live forcing and        !
      !      the fast loop. A region starts from bare ground or from a region state file (one        !
      !      census does not describe many cells), and the one-file probe stays off (every           !
      !      polygon would write into it). ----------------------------------------------------------!
      if (cfg%run_mode == RUN_MODE_REGION) then
         if (.not. (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on))                         &
            error stop TAG//'[run].mode = "region" needs fast.fast_biophysics_on and forcing.forcing_on'
         if (cfg%forcing%backend /= MET_BACKEND_ED_ERA5LAND)                                       &
            error stop TAG//'[run].mode = "region" needs forcing.format = "ED_ERA5land"'
         if (cfg%fast_probe)                                                                       &
            error stop TAG//'[run].mode = "region" cannot write fast.fast_probe (one CSV per run)'
         if (cfg%init_mode == INIT_CENSUS)                                                         &
            error stop TAG//'[run].mode = "region" starts from bare ground or a region state file '// &
                            '(init.init_mode = 0 or 2), not a census'
         associate (b => cfg%region%box_nwse)
            if (.not. (b(1) > b(3) .and. b(1) <= 90.0_wp .and. b(3) >= -90.0_wp))               &
               error stop TAG//'region.box_nwse needs -90 <= S < N <= 90 ([N, W, S, E])'
            if (any(abs(b([2, 4])) > 360.0_wp))                                                    &
               error stop TAG//'region.box_nwse longitudes must lie within [-360, 360]'
         end associate
         if (cfg%region%land_fraction_min < 0.0_wp .or. cfg%region%land_fraction_min > 1.0_wp)     &
            error stop TAG//'region.land_fraction_min must lie in [0, 1]'
         if (len_trim(cfg%region%cell_mask) > 0 .and. len_trim(cfg%region%cell_mask_var) == 0)       &
            error stop TAG//'region.cell_mask_var must name the mask''s variable'
         if (.not. ieee_is_finite(cfg%region%cell_mask_min))                                       &
            error stop TAG//'region.cell_mask_min must be a finite number'
         !----- Every polygon takes its soil from its own cell, and so does its height              !
         !      allometry's site index: a region does not run one hand-set soil everywhere. --------!
         if (.not. soil_maps_on(cfg%maps))                                                       &
            error stop TAG//'[run].mode = "region" needs the soil maps: maps.soil_depth_map, '//    &
                            'soil_texture_map and soil_chemistry_map'
         if (any(cfg%pft%allom(1:cfg%pft%n)%height_form == HEIGHT_CHAVE)                          &
             .and. len_trim(cfg%maps%chave_e_map) == 0)                                          &
            error stop TAG//'[run].mode = "region" with allometry.height_allometry = "chave2014" '//   &
                            'needs maps.chave_e_map'
         !----- A region loads a month's forcing once, before the month (region_step_month), so the  !
         !      recycle seam must fall on a month boundary: a window starting at 00:00 or 01:00 on the !
         !      1st. A site run takes a window starting at any record stamp. ------------------------!
         if (cfg%forcing%recycle .and. time_valid(cfg%forcing%recycle_start)) then
            associate (a => cfg%forcing%recycle_start)
               if (a%day /= 1_ik .or. a%hour > 1_ik .or. a%minute /= 0_ik .or. a%second /= 0_ik)    &
                  error stop TAG//'[run].mode = "region" needs forcing.recycle_start at 00:00 or '//  &
                                  '01:00 on the 1st of a month'
            end associate
         end if
      end if
      !----- [maps]: the maps are on the ED_ERA5land grid, and a site reads them at its forcing    !
      !      cell, so they need that archive. A chave_e map serves the Chave height allometry only. -!
      if (soil_maps_on(cfg%maps) .or. len_trim(cfg%maps%chave_e_map) > 0) then
         if (.not. (cfg%forcing%forcing_on .and. cfg%forcing%backend == MET_BACKEND_ED_ERA5LAND))     &
            error stop TAG//'the [maps] files need forcing.format = "ED_ERA5land" (they are on its grid)'
      end if
      if (len_trim(cfg%maps%chave_e_map) > 0 .and. .not. any(cfg%pft%allom(1:cfg%pft%n)%height_form == HEIGHT_CHAVE)) &
         error stop TAG//'maps.chave_e_map needs a PFT with allometry.height_allometry = "chave2014"'
      !----- Forcing: the forcing's own heights must be physical, and an open-terrain wind must have  !
      !      been made above its exposure roughness and below its blending height. No reference height  !
      !      has to clear the canopy: the forcing is moved to each patch's canopy-air top.              !
      if (cfg%forcing%forcing_on) then
         if (cfg%forcing%tq_height <= 0.0_wp)   error stop TAG//'forcing.tq_height must be > 0'
         if (cfg%forcing%wind_height <= 0.0_wp) error stop TAG//'forcing.wind_height must be > 0'
         if (cfg%forcing%wind_exposure == WIND_EXPOSURE_OPEN_TERRAIN) then
            if (cfg%forcing%wind_exposure_z0 <= 0.0_wp .or.                                       &
                cfg%forcing%wind_exposure_z0 >= cfg%forcing%wind_height)                           &
               error stop TAG//'forcing.wind_exposure_z0 must lie in (0, wind_height)'
            if (cfg%forcing%wind_blending_height < cfg%forcing%wind_height)                         &
               error stop TAG//'forcing.wind_blending_height must be >= wind_height'
         end if
         !----- The terrain lapse rate: between isothermal and about the dry adiabat. --------------!
         if (any(cfg%forcing%lapse_rate_tair < 0.0_wp) .or. any(cfg%forcing%lapse_rate_tair > 0.01_wp)) &
            error stop TAG//'site.lapse_rate_tair must lie in [0, 0.01] K/m'
         !----- Longwave synthesis (docs/science/forcing.md eq. (17)): its cloud coefficient       !
         !      must be physical -- a negative cloud term would make a cloudy sky emit LESS        !
         !      than a clear one. -----------------------------------------------------------------!
         if (cfg%forcing%lwdown_source == LW_SYNTHESIZE .and. cfg%forcing%lw_cloud_a < 0.0_wp)     &
            error stop TAG//'forcing.lw_cloud_a must be >= 0 (a cloudy sky emits MORE, not less)'
         !----- avg_convention: "instant" and "center" PARSE, but only METAVG_BEGIN and            !
         !      METAVG_END have a disaggregation branch, and everything else would fall            !
         !      through to METAVG_END -- a scheme the user did not ask for, silently. Each is      !
         !      refused until it has its own branch, because a forcing switch that does            !
         !      nothing is worse than one that is absent. -----------------------------------------!
         if (cfg%forcing%avg_convention == METAVG_INSTANT)                                       &
            error stop TAG//'forcing.avg_convention = "instant" is not implemented (it would '//  &
                            'silently run the "end" path); use "end" or "begin"'
         if (cfg%forcing%avg_convention == METAVG_CENTER)                                        &
            error stop TAG//'forcing.avg_convention = "center" is not implemented (it would '//   &
                            'silently run the "end" path); use "end" or "begin"'
         !----- The ED_ERA5land archive (docs/science/forcing.md §1) is hourly, end-stamped        !
         !      and stores TOTAL shortwave; a config that says otherwise would mis-time or         !
         !      mis-partition every record. -------------------------------------------------------!
         if (cfg%forcing%backend == MET_BACKEND_ED_ERA5LAND) then
            if (len_trim(cfg%forcing%data_path) == 0) error stop TAG//'forcing.data_path is empty'
            if (cfg%run_mode /= RUN_MODE_REGION .and. cfg%forcing%max_distance_km <= 0.0_wp)      &
               error stop TAG//'forcing.max_distance_km must be > 0'
            if (abs(cfg%forcing%dt_forcing - ARCHIVE_DT_SEC) > 0.5_wp)                             &
               error stop TAG//'forcing.timestep must be 1 hour for format = "ED_ERA5land"'
            if (cfg%forcing%avg_convention /= METAVG_END)                                          &
               error stop TAG//'forcing.avg_convention must be "end" for format = "ED_ERA5land"'
            if (cfg%forcing%sw_partition == SWPART_PASSTHROUGH)                                    &
               error stop TAG//'forcing.sw_partition cannot be "passthrough" for format = "ED_ERA5land" '// &
                               '(the archive stores total shortwave)'
            !----- The reader loads an archive month, plus the record before it, before each         !
            !      step. A daily step from midnight reads exactly that; a longer step, or one        !
            !      starting mid-day, can straddle two months. ---------------------------------------!
            if (abs(cfg%dt_slow - 86400.0_wp) > 0.5_wp)                                             &
               error stop TAG//'format = "ED_ERA5land" needs [run].dt_slow = "1d" (the reader loads a month at a time)'
            if (cfg%start_time%hour /= 0_ik .or. cfg%start_time%minute /= 0_ik .or.                 &
                cfg%start_time%second /= 0_ik)                                                     &
               error stop TAG//'format = "ED_ERA5land" needs [run].start_time at 00:00:00'
         end if
         !----- RECYCLE WINDOW (docs/science/forcing.md §9): declared, never inferred, and           !
         !      required to be an exact whole number of calendar years. A window of any other        !
         !      length cannot be wrapped without drifting BOTH hour-of-day and day-of-year:          !
         !      wrapping on the 366 d 22 h span of an ERA5-Land Ithaca file shifts the sub-daily     !
         !      shortwave phase by ~10 h and the season by ~48 d over a 29-yr run, while the         !
         !      daily MEAN stays correct (the cosz reconstruction is mean-conserving), so            !
         !      nothing downstream complains. Reject here rather than drift silently. The anchor     !
         !      may sit anywhere in the calendar -- a mid-year window is fine -- only the SPAN       !
         !      must be whole years. ----------------------------------------------------------------!
         if (cfg%forcing%recycle) then
            if (.not. time_valid(cfg%forcing%recycle_start))                                    &
               error stop TAG//'forcing.recycle_start is not a valid date/time'
            if (.not. time_valid(cfg%forcing%recycle_end))                                      &
               error stop TAG//'forcing.recycle_end is not a valid date/time'
            if (.not. time_lt(cfg%forcing%recycle_start, cfg%forcing%recycle_end))              &
               error stop TAG//'forcing.recycle_end must be after forcing.recycle_start'
            if (whole_years_between(cfg%forcing%recycle_start, cfg%forcing%recycle_end,         &
                                    MAX_RECYCLE_YEARS) < 1_ik) then
               write(*,'(4a)') ' meds_config: forcing recycle window ',                          &
                  time_to_string(cfg%forcing%recycle_start), ' .. ',                             &
                  time_to_string(cfg%forcing%recycle_end)
               write(*,'(a)')  '   is not an exact whole number of calendar years. Recycling a'
               write(*,'(a)')  '   partial-year window drifts hour-of-day and day-of-year on every'
               write(*,'(a)')  '   wrap. Set recycle_end to the SAME month/day/time as recycle_start,'
               write(*,'(a)')  '   N years later (it is the exclusive end of the window).'
               error stop TAG//'forcing recycle window is not a whole number of years'
            end if
         end if
      end if
      !----- Leaf phenology: every PFT's cue parameters are validated, from a [phenology] section  !
      !      or the alloc_pft_table defaults. Each mask is a set of TEMP(1), LIGHT(2), WATER(4).     !
      !      A sharpness of 0 would make its switch a constant 0.5, so it is refused.              !
      associate (p => cfg%pft, n => cfg%pft%n)
         if (any(p%pheno_flush_cue_mask(1:n) < 0_ik .or. p%pheno_flush_cue_mask(1:n) > 7_ik) .or.      &
             any(p%pheno_shed_cue_mask(1:n)  < 0_ik .or. p%pheno_shed_cue_mask(1:n)  > 7_ik))          &
            error stop TAG//'phenology flush_cue_mask / shed_cue_mask out of range [0,7]'
         if (any(p%pheno_flush_cue_timescale(1:n) <= 0.0_wp) .or.                                     &
             any(p%pheno_shed_cue_timescale(1:n)  <= 0.0_wp) .or.                                     &
             any(p%pheno_light_window(1:n)        <= 0.0_wp))                                         &
            error stop TAG//'phenology cue timescales and light_window must be > 0'
         if (any(p%pheno_par_min(1:n) < 0.0_wp))                                                      &
            error stop TAG//'phenology.par_min must be >= 0'
         if (any(p%pheno_flush_rate_max(1:n) <= 0.0_wp))                                              &
            error stop TAG//'phenology.flush_rate_max must be > 0'
         if (any(p%pheno_shed_rate_max(1:n) < 0.0_wp))                                                &
            error stop TAG//'phenology.shed_rate_max must be >= 0'
         if (any(p%pheno_min_leaf_cover(1:n) < 0.0_wp .or. p%pheno_min_leaf_cover(1:n) >= 1.0_wp) .or. &
             any(p%pheno_bare_leaf_cover(1:n) < 0.0_wp .or. p%pheno_bare_leaf_cover(1:n) >= 1.0_wp))   &
            error stop TAG//'phenology min_leaf_cover / bare_leaf_cover must be in [0,1)'
         if (any(p%pheno_flush_temp_sharpness(1:n) == 0.0_wp) .or.                                   &
             any(p%pheno_shed_temp_sharpness(1:n)  == 0.0_wp) .or.                                   &
             any(p%pheno_flush_light_sharpness(1:n) == 0.0_wp) .or.                                  &
             any(p%pheno_shed_light_sharpness(1:n)  == 0.0_wp) .or.                                  &
             any(p%pheno_flush_water_sharpness(1:n) == 0.0_wp) .or.                                  &
             any(p%pheno_shed_water_sharpness(1:n)  == 0.0_wp))                                      &
            error stop TAG//'a phenology sharpness is 0: its switch would be a constant 0.5'
      end associate
      if (cfg%cohort_size_tol_min <= 0.0_wp)            error stop TAG//'cohort_size_tol_min <= 0'
      if (cfg%cohort_size_tol_max < cfg%cohort_size_tol_min) error stop TAG//'cohort_size_tol_max < min'
      if (cfg%n_cohort_fusion_iter < 1_ik)                   error stop TAG//'n_cohort_fusion_iter < 1'
      if (cfg%n_patch_fusion_iter < 1_ik)                   error stop TAG//'n_patch_fusion_iter < 1'
      if (cfg%patch_light_tol <= 0.0_wp)                    error stop TAG//'patch_light_tol <= 0'
      if (cfg%init_soil_temp < 233.0_wp .or. cfg%init_soil_temp > 333.0_wp)                          &
         error stop TAG//'init.soil_temp outside 233-333 K'
      if (cfg%init_reacclimate_traits .and. cfg%init_mode /= INIT_RESTART)                          &
         error stop TAG//'init.reacclimate_traits applies to a restart (init.init_mode = 2) only'
      if (cfg%fineroot_turnover_plastic .and. .not. cfg%trait_plasticity_on)                          &
         error stop TAG//'trait_dynamics.fineroot_turnover_plastic needs trait_plasticity_on = true'
      if (cfg%patch_light_tol_max < cfg%patch_light_tol)   error stop TAG//'patch_light_tol_max < patch_light_tol'
      if (cfg%n_height_layers < 2_ik)                error stop TAG//'n_height_layers < 2'
      if (cfg%min_patch_area <= 0.0_wp)              error stop TAG//'min_patch_area <= 0'
      if (cfg%cohort_lai_cap <= 0.0_wp)              error stop TAG//'cohort_lai_cap <= 0'
      if (cfg%growth_memory_days <= 0.0_wp)          error stop TAG//'growth_memory_days <= 0'
      if (cfg%patch_disturbance_rate < 0.0_wp)       error stop TAG//'patch_disturbance_rate < 0'
      if (cfg%disturbance_survive_height <= 0.0_wp)  error stop TAG//'disturbance_survive_height <= 0'
      if (any(cfg%pft%wood_density <= 0.0_wp))       error stop TAG//'wood_density <= 0'
      !----- A leaf turnover rate of 0 would keep every leaf forever, and the plasticity slope    !
      !      (docs/science/plant_traits.md eq. (4)) and the fine roots' share of it divide by it.  !
      if (any(cfg%pft%leaf_turnover_rate_toc <= 0.0_wp)) error stop TAG//'leaf_turnover_rate_toc <= 0'
      if (any(cfg%pft%stem_resp_sapwood25 < 0.0_wp)) error stop TAG//'stem_resp_sapwood25 < 0'
      if (any(cfg%pft%max_relative_growth_rate < 0.0_wp))                                          &
         error stop TAG//'max_relative_growth_rate < 0'
      if (any(cfg%pft%max_absolute_growth_rate < 0.0_wp))                                          &
         error stop TAG//'max_absolute_growth_rate < 0'
      !----- A recruit must survive its own birth: pool threshold must exceed the cull. ---!
      if (cfg%min_recruit_size <= cfg%negligible_nplant)                                   &
         error stop TAG//'min_recruit_size must exceed negligible_nplant'

      !----- THERMAL ACCLIMATION shifts the PEAKED form's entropy term                              !
      !      (docs/science/leaf_gas_exchange.md eq. (4)). The Arrhenius form has no dS, so the      !
      !      flag would be silently inert -- refuse rather than let a config believe                !
      !      acclimation is on when nothing acclimates. --------------------------------------------!
      if (cfg%leaf_thermal_acclimation .and. cfg%temp_response_form /= TEMP_RESPONSE_PEAKED)              &
         error stop TAG//'leaf_physiology.thermal_acclimation requires temp_response_form = '//    &
                    '"peaked" (the Arrhenius form has no entropy term to shift)'
      if (cfg%leaf_thermal_acclimation .and. cfg%acclim_window_days <= 0.0_wp)                    &
         error stop TAG//'leaf_physiology.acclim_window_days must be > 0'

      !----- Leaf physiology: shared biochemistry scalars (Kc/Ko/Gamma* are denominators). -!
      if (cfg%kc25 <= 0.0_wp)            error stop TAG//'kc25 <= 0'
      if (cfg%ko25 <= 0.0_wp)            error stop TAG//'ko25 <= 0'
      if (cfg%gstar25 <= 0.0_wp)         error stop TAG//'gstar25 <= 0'
      if (cfg%o2_mol_frac <= 0.0_wp)     error stop TAG//'o2_mol_frac <= 0'
      if (cfg%leaf_absorptance <= 0.0_wp) error stop TAG//'leaf_absorptance <= 0'
      if (cfg%phi_psii <= 0.0_wp)        error stop TAG//'phi_psii <= 0'
      if (cfg%medlyn_vpd_min <= 0.0_wp)  error stop TAG//'medlyn_vpd_min <= 0 [kPa]'
      !----- Leaf physiology: per-PFT traits. ---------------------------------------------!
      if (any(cfg%pft%photosynthetic_pathway /= PATH_C3 .and.                              &
              cfg%pft%photosynthetic_pathway /= PATH_C4)) error stop TAG//'photosynthetic_pathway not in {1,2}'
      if (any(cfg%pft%vcmax25 <= 0.0_wp))            error stop TAG//'vcmax25 <= 0'
      if (any(cfg%pft%jmax_vcmax_ratio <= 0.0_wp))   error stop TAG//'jmax_vcmax_ratio <= 0'
      if (any(cfg%pft%tpu_vcmax_ratio <= 0.0_wp))    error stop TAG//'tpu_vcmax_ratio <= 0'
      if (any(cfg%pft%rd_vcmax_ratio < 0.0_wp))      error stop TAG//'rd_vcmax_ratio < 0'
      if (any(cfg%pft%stomatal_g0 < 0.0_wp))         error stop TAG//'stomatal_g0 < 0'
      if (any(cfg%pft%stomatal_g1 < 0.0_wp))         error stop TAG//'stomatal_g1 < 0'
      if (any(cfg%pft%stomatal_d0 <= 0.0_wp))        error stop TAG//'stomatal_d0 <= 0 (Leuning divides by it)'
      if (any(cfg%pft%katul_lambda25 <= 0.0_wp))     error stop TAG//'katul_lambda25 <= 0'
      if (any(cfg%pft%wstress_lambda_exp < 0.0_wp .or. cfg%pft%wstress_lambda_exp > 8.0_wp)) &
         error stop TAG//'wstress_lambda_exp must be in [0,8]'
      do ipft = 1_ik, cfg%pft%n
         if (pft_wstress_psi_open(cfg, ipft) > 0.0_wp) error stop TAG//'wstress_psi_open must be <= 0'
         if (pft_wstress_psi_close(cfg, ipft) >= pft_wstress_psi_open(cfg, ipft))           &
            error stop TAG//'wstress_psi_close must be below wstress_psi_open'
      end do
      if (any(cfg%pft%wstress_sref_stomata <= 0.0_wp))                                     &
         error stop TAG//'wstress_sref_stomata must be > 0 (beta_stomata = exp(sref*(psi - psi_onset)))'
      if (any(cfg%pft%stomata_psi_onset > 0.0_wp))                                         &
         error stop TAG//'stomata_psi_onset must be <= 0 (a predawn leaf water potential)'
      if (any(cfg%pft%leaf_surf_water_max < 0.0_wp) .or. any(cfg%pft%wood_surf_water_max < 0.0_wp)) &
         error stop TAG//'leaf_surf_water_max and wood_surf_water_max must be >= 0 [kg/m2 of surface]'
      !----- C3 uses theta_j (the J hyperbola / co-limitation curvature); C4 does not. -----!
      if (any(cfg%pft%photosynthetic_pathway == PATH_C3 .and.                              &
              (cfg%pft%theta_j <= 0.0_wp .or. cfg%pft%theta_j >= 1.0_wp)))                 &
         error stop TAG//'C3 theta_j must be in (0,1)'
      !----- C4-only constraints (PEPcase slope, light slope, the two co-limitation curvatures).-!
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and. cfg%pft%kp25 <= 0.0_wp))     &
         error stop TAG//'C4 PFT needs kp25 > 0'
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and. cfg%pft%quantum_yield_c4 <= 0.0_wp)) &
         error stop TAG//'C4 PFT needs quantum_yield_c4 > 0'
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and.                              &
              (cfg%pft%theta_cj_c4 <= 0.0_wp .or. cfg%pft%theta_cj_c4 >= 1.0_wp)))         &
         error stop TAG//'C4 theta_cj_c4 must be in (0,1)'
      if (any(cfg%pft%photosynthetic_pathway == PATH_C4 .and.                              &
              (cfg%pft%theta_ic_c4 <= 0.0_wp .or. cfg%pft%theta_ic_c4 >= 1.0_wp)))         &
         error stop TAG//'C4 theta_ic_c4 must be in (0,1)'
   end subroutine validate_config

end module meds_config
