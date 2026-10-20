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
   use meds_constants,  only : yr_day, yr_sec
   use meds_pft_params, only : pft_table_t, derive_leaf_params, HYD_UNSET
   use meds_time,       only : meds_time_t
   use meds_therm_lib,     only : TEMP_RESPONSE_PEAKED
   use meds_leaf_opts,     only : SM_MEDLYN, COLIM_QUADRATIC
   use meds_plant_hydr_lib,  only : tissue_psi_tlp
   use meds_soil_column_config, only : soil_column_config_t
   use meds_forcing_config, only : forcing_config_t
   use meds_output_config,  only : output_config_t
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   use meds_biogeochem_opts, only : decomp_opts_t
   use meds_region_opts,   only : region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   use meds_maps_opts,  only : maps_opts_t, soil_maps_on
   use meds_sites_opts, only : site_spec_t
   use meds_fast_opts,  only : INTEG_ARK, CTRL_I, CTRL_L1_ADAPTIVE
   implicit none
   private

   public :: meds_config_t, site_config
   public :: derive_config, derive_parameters
   public :: growth_window_steps
   public :: pft_leaf_psi_tlp, pft_stomata_psi_onset, pft_wstress_psi_open, pft_wstress_psi_close
   public :: forcing_config_t, output_config_t
   public :: decomp_opts_t
   public :: region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   public :: maps_opts_t, soil_maps_on
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART


   !----- Initialization modes ([init].init_mode). ------------------------------------------!
   integer(ik), parameter :: INIT_BARE    = 0_ik    !< near-bare ground
   integer(ik), parameter :: INIT_CENSUS  = 1_ik    !< from a cohort census CSV (init_census_file)
   integer(ik), parameter :: INIT_RESTART = 2_ik    !< restart from a state .nc file (init_restart_file)

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
      real(wp)    :: rtol                 = 1.0e-3_wp   !< [-] relative tolerance of every error-controlled store
      !----- The ABSOLUTE-tolerance companion to rtol. Every group's atol is physically scaled     !
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
      real(wp)    :: patch_light_tol, patch_light_maxdev_factor
      !----- The light-profile tolerance steps from patch_light_tol to this ceiling and no further:   !
      !      two patches more different than it stay apart even above max_patch, which makes the   !
      !      patch count a target rather than a hard limit.                                          !
      real(wp)    :: patch_light_tol_max
      real(wp)    :: patch_light_tol_mult        !< DERIVED (geometric multiplier)
      real(wp)    :: light_ext_coeff_fusefiss = 0.5_wp   !< [-] extinction of the light profiles fusion compares
      real(wp)    :: min_patch_area

      !----- Patch disturbance, growth memory, recruitment, conservation. -----------------!
      real(wp) :: patch_disturbance_rate, disturbance_survive_height
      real(wp) :: growth_memory_days, min_recruit_size, conservation_tol

      !----- Initial conditions (init_mode: 0 bare | 1 census | 2 restart). ---------------!
      integer(ik)        :: init_mode
      character(len=256) :: init_restart_file = '', init_census_file = ''
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
      real(wp)    :: gpp_ref = 0.0_wp         !< [kgC/m2 leaf/yr] stub GPP per unit leaf area (carbon mode)

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


end module meds_config
