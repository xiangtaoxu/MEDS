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
   use meds_pft_params, only : pft_table_t, PATH_C3, PATH_C4, derive_pft_rates, derive_leaf_params
   use meds_allometry,  only : set_allometry
   use meds_time,       only : meds_time_t, time_lt, time_valid, time_to_string,                &
                               whole_years_between
   use meds_temp_response, only : TRESP_ARRHENIUS, TRESP_PEAKED
   use meds_forcing_config, only : forcing_config_t
   use meds_output_config,  only : output_config_t
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   use meds_biogeochem_opts, only : decomp_opts_t
   implicit none
   private

   public :: meds_config_t, allometry_config_t, hydraulics_config_t, derive_config, derive_parameters
   public :: validate_config, growth_window_steps
   public :: forcing_config_t, output_config_t
   public :: decomp_opts_t
   public :: BK_SERIAL, BK_MULTICORE, BK_GPU
   public :: DIST_PRIMARY, DIST_TREEFALL
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART
   public :: SM_LEUNING, SM_MEDLYN, SM_KATUL
   public :: TRESP_ARRHENIUS, TRESP_PEAKED, COLIM_MIN, COLIM_QUADRATIC
   public :: INTEG_ARK, INTEG_RK4
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
   !----- Stomatal-conductance model (leaf physiology; [leaf_physiology].stomatal_model). -!
   integer(ik), parameter :: SM_LEUNING = 1_ik      !< Leuning (1995) BWB-VPD semi-empirical
   integer(ik), parameter :: SM_MEDLYN  = 2_ik      !< Medlyn et al. (2011) unified optimization (USO)
   integer(ik), parameter :: SM_KATUL   = 3_ik      !< Katul et al. (2010) analytical optimization
   !----- TRESP_ARRHENIUS / TRESP_PEAKED are owned by meds_temp_response, re-exported above. !
   !----- Co-limitation form combining the FvCB / C4 limitation rates. --------------------!
   integer(ik), parameter :: COLIM_MIN       = 1_ik !< sharp minimum
   integer(ik), parameter :: COLIM_QUADRATIC = 2_ik !< smoothed co-limitation quadratics


   !----- Fast-loop TIME integrator ([fast].time_integrator). TWO schemes; the operator-split third   !
   !      was RETIRED 2026-07-31 (it converged to a different limit and could not carry the coupled    !
   !      tissue heat store -- see meds_fast_step's header). The `integration_scheme` coupling-sweep   !
   !      selector went with it: it only ever chose the split's own Gauss-Seidel-vs-Picard sweep.      !
   !                                                                                          !
   !      NAMING, stated once so it stops propagating: INTEG_ARK is NOT an IMEX method. The biotic     !
   !      CO2 source is folded implicit, so the explicit tableau is empty (f_E == 0) and the scheme    !
   !      is a 2-solve ESDIRK2 with gamma = 1 - 1/sqrt(2) (the ARS(2,2,2) value). The config string    !
   !      stays "ark" for compatibility. ------------------------------------------------------------!
   integer(ik), parameter :: INTEG_ARK   = 2_ik  !< ESDIRK2 coupled implicit column (DEFAULT)
   integer(ik), parameter :: INTEG_RK4   = 3_ik  !< adaptive Cash-Karp RK45, the ACCURACY BASELINE

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
      real(wp) :: wood_kmax   = 8.0_wp    !< [kg/m/s/MPa] sapwood specific conductivity
      real(wp) :: vessel_curl = 1.5_wp    !< [-] tortuosity / path-length factor
      !----- Soil->root rhizosphere conductance (prescribed, single-layer MVP). ------------!
      !----- Multi-layer root distribution (MEDS_MULTILAYER_ROOTS_DESIGN). Consumed when the         !
      !       multi-layer root boundary is wired (per-layer soil state); the single-layer path         !
      !       ignores them, so defaults are inert. --------------------------------------------------!
      real(wp) :: root_beta          = 0.96_wp  !< [-]      ED2 root-profile decay (0,1); smaller => shallower
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
      !----- Fast (sub-daily) biophysics loop. --------------------------------------------!
      logical     :: fast_biophysics_on          !< master gate for the fast biophysics loop
      real(wp)    :: dt_fast                      !< [s] fast biophysics timestep (nested within dt_slow)
      integer(ik) :: n_fast_per_slow              !< DERIVED = max(1, nint(dt_slow / dt_fast))
      !----- P3 coupled-surface (Picard) solver knobs + option selectors ([fast], DEFAULTED reads, !
      !      solver tuning not physical params). ---------------------------------------------------!
      integer(ik) :: picard_max_iter     = 20_ik      !< outer-iteration cap
      real(wp)    :: picard_tol_temp      = 1.0e-3_wp  !< [K]     temperature convergence tolerance
      real(wp)    :: picard_tol_shv       = 1.0e-6_wp  !< [kg/kg] CAS humidity convergence tolerance
      real(wp)    :: picard_relax         = 0.5_wp     !< under-relaxation of the next-pass seed
      logical     :: picard_fixed_iter    = .false.    !< GPU warp-uniform fixed pass count (no early exit)
      integer(ik) :: leaf_energy_model    = 0_ik       !< 0 = diagnostic leaf | 1 = prognostic leaf_energy
      integer(ik) :: wood_energy_model    = 0_ik       !< 0 = diagnostic wood | 1 = prognostic wood (own store, never = leaf)
      real(wp)    :: snow_init_swe        = 0.0_wp     !< [kg/m2] initial snow water-equivalent seeded at run start
      real(wp)    :: snow_init_temp       = 270.0_wp   !< [K] initial snow temperature (for the seeded pack)
      logical     :: canopy_water_on      = .false.    !< opt-in canopy interception film + film-evap/dew (P1, split path)
      !----- Fast-loop TIME integrator selector + ARK knobs ([fast], DEFAULTED reads). ----------------!
      !      every existing config + the golden anchor byte-identical). --------------------------------!
      integer(ik) :: time_integrator      = INTEG_ARK !< INTEG_ARK (default) | INTEG_RK4
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
      integer(ik) :: ark_niter            = 8_ik        !< coupled leaf<->CAS Newton cap (>1 => coupled)
      real(wp)    :: ark_relax            = 0.6_wp      !< under-relaxation (vestigial on the Newton branch)
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
      real(wp)    :: min_patch_area, patch_min_area_remain
      logical     :: enable_patch_fission

      !----- Patch disturbance, growth memory, recruitment, conservation. -----------------!
      real(wp) :: patch_disturbance_rate, disturbance_survive_height
      real(wp) :: growth_memory_days, min_recruit_size, conservation_tol

      !----- Initial conditions (init_mode: 0 bare | 1 census | 2 restart). ---------------!
      integer(ik)        :: init_mode
      character(len=256) :: init_restart_file, init_census_file

      !----- netCDF output. ---------------------------------------------------------------!
      logical            :: io_write_output
      character(len=256) :: io_output_dir, io_output_prefix
      integer(ik) :: io_output_interval_years, io_cohort_max, io_patch_max
      logical     :: io_write_state
      integer(ik) :: io_state_interval_years

      !----- Parameter-config controls. ---------------------------------------------------!
      character(len=256) :: pft_config        !< path to the PFT config file (named in the main file)
      logical            :: override_derived  !< if .true., a [derived] block overwrites computed values

      !----- Leaf physiology: model selection (non-PFT). ----------------------------------!
      integer(ik) :: stomatal_model           !< SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form        !< TRESP_ARRHENIUS | TRESP_PEAKED
      integer(ik) :: colimitation             !< COLIM_MIN | COLIM_QUADRATIC
      logical     :: leaf_use_boundary_layer  !< if .true., couple via the leaf boundary layer (gb)
      !----- Leaf physiology: shared biochemistry at 25 degC + Arrhenius/deactivation terms.-!
      real(wp) :: kc25, ko25, gstar25                   !< [Pa]    Michaelis constants + CO2 compensation point
      real(wp) :: ea_kc, ea_ko, ea_gstar                !< [J/mol] activation energies (Bernacchi et al. 2001)
      real(wp) :: ea_vcmax, ea_jmax, ea_rd              !< [J/mol] activation energies
      real(wp) :: hd_vcmax, hd_jmax, hd_rd              !< [J/mol] deactivation energies (peaked form)
      real(wp) :: ds_vcmax, ds_jmax, ds_rd              !< [J/mol/K] entropy terms (peaked form)
      real(wp) :: o2_mol_frac                           !< [mol/mol] atmospheric O2 mole fraction
      real(wp) :: leaf_absorptance                      !< [--] leaf PAR absorptance (for electron transport)
      real(wp) :: phi_psii                              !< [--] PSII quantum yield (electrons/photon)

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
      !       with no [output] block runs the legacy [io] path unchanged). The per-variable overrides !
      !       live in the optional meds_io_config.toml named by output%io_config (§6, MEDS_IO_DESIGN). !
      type(output_config_t) :: output

      !----- PFT traits. ------------------------------------------------------------------!
      type(pft_table_t) :: pft

      !----- Global allometry coefficients (populated by the loader; installed by derive_parameters). !
      type(allometry_config_t) :: allom

      !----- Plant-hydraulics parameters ([hydraulics], opt-in; defaults = MVP placeholders). ------!
      type(hydraulics_config_t) :: hydraulics

      !----- Fast-loop biophysics run-config ([soil]/[energy]/[snow]/[aerodynamics], all opt-in;    !
      !       defaults = meds_biophysics_opts placeholders). build_fast_context copies each verbatim !
      !       into the column config (ccfg%hydro/energy/snow/aero); an absent block is a no-op. ------!
      type(soil_opts_t)   :: soil        !< [soil]         soil-water Richards solver opts (-> ccfg%hydro)
      type(energy_opts_t) :: energy      !< [energy]       soil-thermal solver opts       (-> ccfg%energy)
      type(snow_params_t) :: snow        !< [snow]         snow physical parameter table  (-> ccfg%snow)
      type(aero_cfg_t)    :: aero        !< [aerodynamics] canopy-aerodynamics constants  (-> ccfg%aero)

      !----- Slow soil-carbon matrix ([soil_carbon], opt-in; MEDS_SLOW_DYNAMICS_DESIGN.md Part II   !
      !      B0). soil_carbon_on (default .false.) gates the FEATURE (per-patch state alloc timing/  !
      !      spin-up now; the daily soil_carbon_step + fast-Rh reconciliation land at B2) -- it does   !
      !      NOT gate whether `soil_carbon`'s fields are required: every key is a DEFAULTED read       !
      !      (like [snow]), falling back to its ED2-verified in-type default, so turning the feature   !
      !      on needs no TOML edits beyond soil_carbon_on itself. --------------------------------------!
      logical            :: soil_carbon_on = .false.
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

   !---------------------------------------------------------------------------------------!
   ! Validate a configuration; halt on a setting that would corrupt the run.               !
   !---------------------------------------------------------------------------------------!
   subroutine validate_config(cfg)
      type(meds_config_t), intent(in) :: cfg
      character(len=*), parameter :: tag = 'meds_config: '

      if (cfg%pft%n < 1_ik)                          error stop tag//'empty PFT table'
      if (.not. time_lt(cfg%start_time, cfg%end_time)) error stop tag//'end_time must be after start_time'
      if (cfg%dt_slow <= 0.0_wp)                        error stop tag//'dt_slow <= 0'
      if (cfg%fast_biophysics_on) then
         if (cfg%dt_fast <= 0.0_wp)                     error stop tag//'dt_fast <= 0'
         !----- dt_fast IS A STABILITY PARAMETER, not just an accuracy one. Every surface coupling    !
         !      coefficient (aerodynamic conductances, leaf boundary layer, stomatal series            !
         !      conductance) is FROZEN for the whole dt_fast, while the canopy air it drives is a      !
         !      very low-capacity node: wcap*cp ~ 2.4e4 J/m2/K against surface fluxes of hundreds of   !
         !      W/m2, so 300 W/m2 over a 900 s step is an 11 K excursion. Above a threshold the lag    !
         !      turns into a sustained period-2 oscillation in canopy-air temperature.                 !
         !                                                                                          !
         !      MEASURED on a high-LAI sunlit stand (LAI 3, 500 W/m2), peak-to-peak over consecutive   !
         !      steps: 900 s -> ~8 K, 450 s -> ~8 K, 300 s -> ~5 K, 225 s -> ~4 K, 150 s -> ~1 K       !
         !      (trend only), 100 s -> smooth. Every conservation budget closes to ~1e-6 J THROUGHOUT, !
         !      so no ledger catches this -- hence the explicit warning here.                          !
         !                                                                                          !
         !      The threshold is stand-dependent (it scales with the flux-to-capacity ratio, so denser !
         !      canopies and stronger radiation push it lower). 300 s is a conservative alarm point,   !
         !      not a guarantee. See MEDS_VEG_ENERGY_INTEGRATION_PLAN.md sec 10. ---------------------!
         if (cfg%dt_fast > 300.0_wp) then
            print '(a)', 'WARNING [meds_config]: dt_fast > 300 s. The surface coupling coefficients'
            print '(a)', '  are frozen across the step, and above ~150-225 s that lag drives a'
            print '(a)', '  sustained period-2 canopy-air oscillation (~8 K at 900 s on a high-LAI'
            print '(a)', '  sunlit stand). Conservation budgets do NOT detect it. Prefer dt_fast <= 150 s.'
         end if
         if (cfg%dt_fast > cfg%dt_slow)                 error stop tag//'dt_fast > dt_slow'
         if (abs(cfg%dt_slow / cfg%dt_fast - real(nint(cfg%dt_slow / cfg%dt_fast, ik), wp)) > 1.0e-6_wp) &
            error stop tag//'dt_slow must be an integer multiple of dt_fast'
         if (cfg%time_integrator /= INTEG_ARK .and. cfg%time_integrator /= INTEG_RK4)          &
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
         if (cfg%ark_niter < 1_ik)            error stop tag//'ark_niter < 1'
      end if
      !----- Forcing: the reference height must clear every PFT canopy (ED2 aborts if zref<=hgt_max), !
      !      and the wind-profile roughness must be positive.                                          !
      if (cfg%forcing%forcing_on) then
         if (cfg%forcing%reference_height <= maxval(cfg%pft%hgt_max(1:cfg%pft%n)))               &
            error stop tag//'forcing reference_height must exceed every PFT hgt_max'
         if (cfg%forcing%wind_roughness_z0 <= 0.0_wp) error stop tag//'wind_roughness_z0 <= 0'
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
      !      forcing precondition to enforce here any more. P1-P2 wire the TEMP (bit 1) + PHOTO (bit 8) !
      !      cues only -- WATER (2), HYDRO (4) and LIGHT (16) are rejected in EITHER mask until their    !
      !      soil-water / dmax-leaf-psi / radiation drivers are threaded (P3). Each mask is a bit set    !
      !      in [0,31]; the rate scales must be non-negative (flush strictly positive). -----------------!
      if (any(cfg%pft%pheno_flush_cue_mask(1:cfg%pft%n) < 0_ik .or.                           &
              cfg%pft%pheno_flush_cue_mask(1:cfg%pft%n) > 31_ik) .or.                         &
          any(cfg%pft%pheno_shed_cue_mask(1:cfg%pft%n) < 0_ik .or.                            &
              cfg%pft%pheno_shed_cue_mask(1:cfg%pft%n) > 31_ik))                              &
         error stop tag//'pheno_{flush,shed}_cue_mask out of range [0,31]'
      !----- 22 = WATER(2) | HYDRO(4) | LIGHT(16): the not-yet-wired cue bits. ------------!
      if (any(iand(cfg%pft%pheno_flush_cue_mask(1:cfg%pft%n), 22_ik) /= 0_ik) .or.            &
          any(iand(cfg%pft%pheno_shed_cue_mask(1:cfg%pft%n),  22_ik) /= 0_ik))               &
         error stop tag//'pheno cue WATER(2)/HYDRO(4)/LIGHT(16) not yet wired (P1-P2: TEMP=1, PHOTO=8)'
      if (any(cfg%pft%pheno_k_flush_max(1:cfg%pft%n) <= 0.0_wp))                              &
         error stop tag//'pheno_k_flush_max must be > 0'
      if (any(cfg%pft%pheno_k_shed_max(1:cfg%pft%n) < 0.0_wp))                                &
         error stop tag//'pheno_k_shed_max must be >= 0'
      if (cfg%cohort_size_tol_min <= 0.0_wp)            error stop tag//'cohort_size_tol_min <= 0'
      if (cfg%cohort_size_tol_max < cfg%cohort_size_tol_min) error stop tag//'cohort_size_tol_max < min'
      if (cfg%n_cohort_fusion_iter < 1_ik)                   error stop tag//'n_cohort_fusion_iter < 1'
      if (cfg%n_patch_fusion_iter < 1_ik)                   error stop tag//'n_patch_fusion_iter < 1'
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

      !----- Leaf physiology: shared biochemistry scalars (Kc/Ko/Gamma* are denominators). -!
      if (cfg%kc25 <= 0.0_wp)            error stop tag//'kc25 <= 0'
      if (cfg%ko25 <= 0.0_wp)            error stop tag//'ko25 <= 0'
      if (cfg%gstar25 <= 0.0_wp)         error stop tag//'gstar25 <= 0'
      if (cfg%o2_mol_frac <= 0.0_wp)     error stop tag//'o2_mol_frac <= 0'
      if (cfg%leaf_absorptance <= 0.0_wp) error stop tag//'leaf_absorptance <= 0'
      if (cfg%phi_psii <= 0.0_wp)        error stop tag//'phi_psii <= 0'
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
         error stop tag//'wstress_sref_stomata must be > 0 (beta_stomata = exp(sref*psi_soil))'
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
