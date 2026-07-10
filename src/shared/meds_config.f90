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
   use meds_pft_params, only : pft_table_t, PATH_C3, PATH_C4
   use meds_time,       only : meds_time_t, time_lt
   use meds_temp_response, only : TRESP_ARRHENIUS, TRESP_PEAKED
   use meds_forcing_config, only : forcing_config_t
   use meds_output_config,  only : output_config_t
   implicit none
   private

   public :: meds_config_t, derive_config, validate_config, growth_window_steps
   public :: forcing_config_t, output_config_t
   public :: BK_SERIAL, BK_MULTICORE, BK_GPU
   public :: DIST_PRIMARY, DIST_TREEFALL
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART
   public :: SM_LEUNING, SM_MEDLYN, SM_KATUL
   public :: TRESP_ARRHENIUS, TRESP_PEAKED, COLIM_MIN, COLIM_QUADRATIC
   public :: GS_EMPIRICAL, GS_CARBON
   public :: SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED

   !----- Time-step modes. ----------------------------------------------------------------!
   !----- Parallel backend labels (the actual backend is chosen at COMPILE time via the    !
   !      compiler's do-concurrent target; this is for reporting/reproducibility only).    !
   integer(ik), parameter :: BK_SERIAL    = 0_ik
   integer(ik), parameter :: BK_MULTICORE = 1_ik
   integer(ik), parameter :: BK_GPU       = 2_ik
   !----- Patch disturbance / land-use classes. -------------------------------------------!
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

   !----- Growth source: which provider supplies the demographic rates. --------------------!
   integer(ik), parameter :: GS_EMPIRICAL = 1_ik    !< phenomenological growth/mortality/recruitment
   integer(ik), parameter :: GS_CARBON    = 2_ik    !< carbon-prognostic growth (wood_carbon is the size anchor)

   !----- Fast-loop coupling scheme ([fast].integration_scheme; how the driver couples the     !
   !      stores across a dt_fast step -- NOT a global-integrator switch). -------------------!
   integer(ik), parameter :: SCHEME_SPLIT_SEQUENTIAL = 1_ik  !< operator-split Gauss-Seidel sweep (default)
   integer(ik), parameter :: SCHEME_PICARD_COUPLED   = 2_ik  !< outer Picard leaf<->CAS fixed point (opt-in)

   !----- NO hard-coded defaults: every field is set by the config reader (presence-mapped) or  !
   !       derived (derive_config / derive_pft_rates). DERIVED fields are noted.  --------------!
   type :: meds_config_t
      !----- Time stepping (run bounded by start/end calendar dates). ----------------------!
      real(wp)          :: dt_slow               !< [s] slow-process timestep (user resolution; default 1 d)
      real(wp)          :: dt_years              !< DERIVED = dt_slow / yr_sec
      type(meds_time_t) :: start_time, end_time
      logical     :: demography_on               !< if .false. structure is frozen
      !----- Fast (sub-daily) biophysics loop. --------------------------------------------!
      logical     :: fast_biophysics_on          !< master gate for the fast biophysics loop
      real(wp)    :: dt_fast                      !< [s] fast biophysics timestep (nested within dt_slow)
      integer(ik) :: integration_scheme           !< SCHEME_SPLIT_SEQUENTIAL | SCHEME_PICARD_COUPLED
      integer(ik) :: n_fast_per_slow              !< DERIVED = max(1, nint(dt_slow / dt_fast))
      !----- P3 coupled-surface (Picard) solver knobs + option selectors ([fast], DEFAULTED reads, !
      !      solver tuning not physical params -- consumed only under SCHEME_PICARD_COUPLED). ------!
      integer(ik) :: picard_max_iter     = 20_ik      !< outer-iteration cap
      real(wp)    :: picard_tol_temp      = 1.0e-3_wp  !< [K]     temperature convergence tolerance
      real(wp)    :: picard_tol_shv       = 1.0e-6_wp  !< [kg/kg] CAS humidity convergence tolerance
      real(wp)    :: picard_relax         = 0.5_wp     !< under-relaxation of the next-pass seed
      logical     :: picard_fixed_iter    = .false.    !< GPU warp-uniform fixed pass count (no early exit)
      integer(ik) :: leaf_energy_model    = 0_ik       !< 0 = diagnostic leaf | 1 = prognostic leaf_energy
      integer(ik) :: soil_water_coupling  = 0_ik       !< 0 = soil water re-solved each pass (lagged=coupled for now)
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

      !----- Carbon-driven growth (opt-in). Default GS_EMPIRICAL keeps the phenomenological  !
      !       path (and a bit-identical spin-up); GS_CARBON drives wood_carbon from NPP.       !
      integer(ik) :: growth_source            !< GS_EMPIRICAL | GS_CARBON
      real(wp)    :: gpp_ref                  !< [kgC/m2 leaf/yr] stub GPP per unit leaf area (carbon mode)

      !----- Leaf phenology ([phenology]). OPT-IN: phenology_on default .false. (a config with no    !
      !       [phenology] block leaves the phenology status hard-wired ON, bit-identical to before).   !
      !       When ON, the slow-loop driver advances the per-cohort cue memory from the daily-mean     !
      !       temperature (needs the fast loop + met forcing) and the carbon leaf-flush gate obeys it.  !
      !       Per-PFT cue params live in the PFT trait table (cfg%pft%pheno_*).                         !
      logical     :: phenology_on = .false.

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
         if (cfg%dt_fast > cfg%dt_slow)                 error stop tag//'dt_fast > dt_slow'
         if (abs(cfg%dt_slow / cfg%dt_fast - real(nint(cfg%dt_slow / cfg%dt_fast, ik), wp)) > 1.0e-6_wp) &
            error stop tag//'dt_slow must be an integer multiple of dt_fast'
         if (cfg%integration_scheme /= SCHEME_SPLIT_SEQUENTIAL .and.                            &
             cfg%integration_scheme /= SCHEME_PICARD_COUPLED)                                   &
            error stop tag//'integration_scheme out of range'
      end if
      !----- Forcing: the reference height must clear every PFT canopy (ED2 aborts if zref<=hgt_max), !
      !      and the wind-profile roughness must be positive.                                          !
      if (cfg%forcing%forcing_on) then
         if (cfg%forcing%reference_height <= maxval(cfg%pft%hgt_max(1:cfg%pft%n)))               &
            error stop tag//'forcing reference_height must exceed every PFT hgt_max'
         if (cfg%forcing%wind_roughness_z0 <= 0.0_wp) error stop tag//'wind_roughness_z0 <= 0'
      end if
      !----- Leaf phenology (opt-in). The daily-mean temperature that drives the temperature/GDD cue  !
      !      comes from the fast loop's met forcing, so both must be on. v1 wires the TEMP (bit 1) and  !
      !      PHOTO (bit 8) cues only -- the WATER (bit 2) and HYDRO (bit 4) cue bits are rejected until  !
      !      their soil-water / leaf-psi drivers are threaded into the slow-loop env. cue_mask is a bit  !
      !      set in [0,15]; the tri-state deadband needs off_threshold < on_threshold.                  !
      if (cfg%phenology_on) then
         if (.not. cfg%fast_biophysics_on)                                                       &
            error stop tag//'phenology_on requires fast_biophysics_on (daily-mean temperature source)'
         if (.not. cfg%forcing%forcing_on)                                                       &
            error stop tag//'phenology_on requires forcing_on (met air temperature)'
         if (any(cfg%pft%pheno_cue_mask(1:cfg%pft%n) < 0_ik .or.                                 &
                 cfg%pft%pheno_cue_mask(1:cfg%pft%n) > 15_ik))                                   &
            error stop tag//'pheno_cue_mask out of range [0,15]'
         if (any(iand(cfg%pft%pheno_cue_mask(1:cfg%pft%n), 6_ik) /= 0_ik))                       &
            error stop tag//'pheno_cue_mask WATER(2)/HYDRO(4) cues not yet wired (v1: TEMP=1, PHOTO=8)'
         if (any(cfg%pft%pheno_off_threshold(1:cfg%pft%n) >=                                     &
                 cfg%pft%pheno_on_threshold(1:cfg%pft%n)))                                       &
            error stop tag//'pheno_off_threshold must be below pheno_on_threshold'
      end if
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
