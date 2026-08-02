!==========================================================================================!
! meds_plant_types -- ALL derived types of the plant-ecophysiology library, in ONE module.  !
!                                                                                          !
! Consolidated from the former per-domain meds_{leaf,hydro,pheno}_types modules (the plant   !
! library compiles as a whole, so there is no reason to keep them apart). The types are pure  !
! DATA -- no methods, no hidden state -- organized in clearly delimited sections:             !
!   * LEAF        -- leaf_env_t / leaf_flux_t / leaf_photo_params_t + limitation & pathway flags. !
!   * HYDRAULICS  -- hydro_env_t / hydro_params_t / hydro_opts_t / hydro_flux_t + topology flags.  !
!   * PHENOLOGY   -- pheno_env_t / pheno_params_t / pheno_state_t / pheno_out_t + cue/status flags. !
! (Respiration types wood_*/root_* will be added here when meds_plant_respiration lands.)      !
!==========================================================================================!
module meds_plant_types
   use meds_kinds,      only : wp, ik
   use meds_pft_params, only : PATH_C3, PATH_C4
   use meds_hydr_lib, only : hydro_table_t
   implicit none
   private

   !----- LEAF -----------------------------------------------------------------------------!
   public :: leaf_env_t, leaf_flux_t, leaf_photo_params_t
   public :: PATH_C3, PATH_C4                       ! re-export (pathway lives with the PFT traits)
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP
   !----- HYDRAULICS -----------------------------------------------------------------------!
   public :: hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t
   public :: N_HYDRO, NODE_LEAF, NODE_STEM, NODE_ROOT, NODE_WOOD, NROOT_MAX
   public :: HYDRO_NODES_2, HYDRO_NODES_3
   public :: HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE
   public :: HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT
   public :: HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED
   !----- PHENOLOGY ------------------------------------------------------------------------!
   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO, CUE_LIGHT
   !----- RESPIRATION (env/flux types removed: the kernels are now `elemental pure` over bare      !
   !      scalars, MEDS_NUMERICS_SCOPING.md §11; only the run-uniform trait PODs remain). ----------!
   public :: wood_params_t, root_params_t
   !----- CARBON ALLOCATION: the allocation kernel is now elemental over plain cohort scalars  !
   !       (meds_plant_carbon_allocation), so it needs NO derived types. Tissue turnover moved   !
   !       into the phenology section (baseline shed rate = degenerate phenology).               !

   !=======================================================================================!
   !     LEAF -- leaf-level gas-exchange interface seam.                                    !
   !=======================================================================================!
   !----- Limitation-regime flag of the binding term at the solution. ---------------------!
   integer(ik), parameter :: LIM_NONE    = 0_ik     !< degenerate (night / net <= 0)
   integer(ik), parameter :: LIM_RUBISCO = 1_ik     !< Rubisco-limited (C3 Ac) / Vcmax-limited (C4)
   integer(ik), parameter :: LIM_RUBP    = 2_ik     !< RuBP / light-limited (Aj)
   integer(ik), parameter :: LIM_PRODUCT = 3_ik     !< triose-phosphate-use limited (C3 Ap)
   integer(ik), parameter :: LIM_C4_PEP  = 4_ik     !< C4 PEPcase CO2 limitation (Ap)

   !----- Leaf-level environmental drivers (no canopy RT / energy balance / hydraulics). ---!
   type :: leaf_env_t
      real(wp) :: par        !< [umol photon/m2/s] incident PAR (leaf absorptance applied internally)
      real(wp) :: leaf_temp  !< [K]    leaf temperature
      real(wp) :: vpd        !< [Pa]   leaf-to-air vapour-pressure deficit
      real(wp) :: ca         !< [umol/mol] reference (canopy-air) CO2 mole fraction
      real(wp) :: pressure   !< [Pa]   air pressure
      real(wp) :: psi_leaf   !< [MPa]  leaf water potential (<= 0); drives beta_nonstomata (capacity)
      real(wp) :: gb         !< [mol H2O/m2/s] boundary-layer conductance (<= 0 => skip, Cs = Ca)
      !----- THE WATER-STATUS POTENTIAL DRIVING THE STOMATAL LIMB (<= 0; default 0 = well-watered).  !
      !                                                                                          !
      !      Deliberately just `psi`, not `psi_soil`. Sabot et al. (2022) write the stomatal limb as   !
      !      beta_stomata = min(1, exp(s_ref * psi)) with psi the SOIL (equivalently predawn) water     !
      !      potential -- the plant's water supply, as opposed to `psi_leaf` a few lines up, which is   !
      !      its instantaneous demand-side tension and drives the separate NON-stomatal (capacity)      !
      !      limb. This field is the supply term. Naming it `psi_soil` claimed more than the kernel     !
      !      knows, because what MEDS actually supplies is not soil potential:                          !
      !                                                                                          !
      !      MEDS passes `dmax_psi_leaf` -- the cohort's own PREDAWN (previous-day daily-maximum) LEAF  !
      !      potential. That is the right choice, and it is closer to Sabot's intent than a soil value  !
      !      would be. Predawn leaf potential is what the FIELD measurement behind these                !
      !      parameterizations actually is (a pre-sunrise leaf sample), and it is the plant's own       !
      !      integration over its whole rooted profile -- so it already contains rooting depth, per-PFT !
      !      vulnerability, and the ~0.3 MPa gravity head between a 30 m tree and a sapling, none of    !
      !      which any single soil layer's psi carries.                                                 !
      !                                                                                          !
      !      The two agree only in WET soil, where the plant re-equilibrates with the soil overnight.   !
      !      They do NOT agree under drought: the wood<->soil relaxation time tau_w = C_wood/rhizo is   !
      !      ~9 s at theta 0.25 but ~4.8 DAYS at theta 0.10, so a droughted cohort never equilibrates   !
      !      overnight and its predawn leaf potential stays well below the soil's. Drought is precisely !
      !      the regime this limb exists for, so the distinction is not academic -- and `s_ref` is      !
      !      calibrated against whichever quantity you believe is being passed. Hence the neutral name. !
      !                                                                                          !
      !      A caller with no predawn history (a recruit, or the first step of a run) is seeded from    !
      !      the surface-layer soil potential, which is the one place a genuine soil psi enters here.   !
      real(wp) :: psi = 0.0_wp  !< [MPa] supply-side water potential driving beta_stomata (Sabot 2022)
   end type leaf_env_t

   !----- Leaf-level fluxes returned by the solver. ---------------------------------------!
   type :: leaf_flux_t
      real(wp)    :: A_net        !< [umol CO2/m2 leaf/s] net assimilation (A_gross - Rd)
      real(wp)    :: A_gross      !< [umol CO2/m2 leaf/s] gross assimilation
      real(wp)    :: gs           !< [mol H2O/m2 leaf/s]  stomatal conductance to water vapour
      real(wp)    :: ci           !< [umol/mol]      intercellular CO2 mole fraction
      real(wp)    :: cs           !< [umol/mol]      leaf-surface CO2 mole fraction (= ca if no BL)
      real(wp)    :: transpiration!< [mol H2O/m2 leaf/s]  E = gs * VPD / pressure
      real(wp)    :: rd           !< [umol/m2 leaf/s]     leaf respiration used
      integer(ik) :: limitation   !< LIM_* of the binding term
      logical     :: converged    !< .true. if the Ci solve met tolerance
      !----- WATER-STRESS LIMBS, surfaced as DIAGNOSTICS (they do not feed back into anything     !
      !      here -- the kernel has already applied them). They exist because a stress closure      !
      !      whose two limbs are never observable is exactly the kind of thing that stays silently  !
      !      inert: env%psi defaulted to 0 and beta_stomata was identically 1 for the whole life    !
      !      of the feature until issue #95 caught it. Reporting them makes that failure mode        !
      !      visible in the output file instead of only in a code read.                              !
      real(wp)    :: beta_stomata    = 1.0_wp  !< [-] stomatal limb min(1, exp(sref*psi))  (1 = unstressed)
      real(wp)    :: beta_nonstomata = 1.0_wp  !< [-] capacity limb, the psi_leaf ramp on Vcmax/Jmax/TPU
   end type leaf_flux_t

   !----- Flat per-PFT parameter set (per-PFT traits + shared biochemistry), self-contained !
   !       so the solver never references meds_config. Filled by the interface from the config.!
   type :: leaf_photo_params_t
      integer(ik) :: pathway        !< PATH_C3 | PATH_C4
      !----- Per-PFT capacities at 25 degC and stomatal/water-stress traits. ---------------!
      real(wp) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(wp) :: g0, g1, d0, quantum_yield, theta_j, theta_cj, theta_ic
      real(wp) :: lambda25, psi_open, psi_close, lambda_psi_exp, sref_stomata
      !----- Leaf turgor-loss point [MPa], from pv_psi_tlp(leaf_pi0, leaf_elastic_mod). Below TWICE  !
      !      this (i.e. far past the point where the leaf has lost all turgor) the stomata are shut  !
      !      HARD -- see the psi_shut branch in solve_leaf_gas_exchange. -------------------------!
      real(wp) :: psi_tlp = -2.0_wp
      !----- ARREST_* selector; the clamp branch only fires for ARREST_GS_CLAMP. ---------------!
      integer(ik) :: stress_arrestor = 1_ik
      !----- Apply the NON-STOMATAL (capacity) water-stress limb? Default .false. -- see            !
      !      meds_config_t%leaf_wstress_nonstomatal for why (issue #47). The stomatal limb has no    !
      !      such switch: it is driven by psi and is the better-constrained of the two. --------!
      logical  :: wstress_nonstomatal = .false.
      !----- Shared biochemistry constants at 25 degC + activation/deactivation terms. -----!
      real(wp) :: kc25, ko25, gstar25
      real(wp) :: ea_kc, ea_ko, ea_gstar, ea_vcmax, ea_jmax, ea_rd
      real(wp) :: hd_vcmax, hd_jmax, hd_rd, ds_vcmax, ds_jmax, ds_rd
      real(wp) :: o2_mol_frac, absorptance, phi_psii
   end type leaf_photo_params_t

   !=======================================================================================!
   !     HYDRAULICS -- stateless per-individual water-transport interface seam.             !
   !=======================================================================================!
   !----- Fixed compile-time node dimension (max over supported topologies). The active count  !
   !      is a run option; the state array shape never changes, so it is GPU/SoA-friendly.      !
   integer(ik), parameter :: N_HYDRO   = 3_ik
   integer(ik), parameter :: NODE_LEAF = 1_ik   !< leaf water pool
   integer(ik), parameter :: NODE_STEM = 2_ik   !< stem (3-node) / lumped wood (2-node)
   integer(ik), parameter :: NODE_WOOD = 2_ik   !< alias of NODE_STEM in the 2-node run
   integer(ik), parameter :: NODE_ROOT = 3_ik   !< root (3-node only)

   !----- Fixed compile-time max soil/root layers for the multi-layer root boundary (§ MEDS_       !
   !      MULTILAYER_ROOTS_DESIGN). Value-type arrays of this length keep hydro_env_t/hydro_flux_t   !
   !      fixed-shape (GPU/SoA-safe). MUST be >= the soil column's n_soil_layer_max when coupled.    !
   integer(ik), parameter :: NROOT_MAX = 20_ik

   !----- Topology (number of active nodes). -----------------------------------------------!
   integer(ik), parameter :: HYDRO_NODES_2 = 2_ik   !< leaf + lumped wood (default)
   integer(ik), parameter :: HYDRO_NODES_3 = 3_ik   !< leaf + stem + root (opt-in)

   !----- Sub-step integrator. -------------------------------------------------------------!
   integer(ik), parameter :: HYDRO_SOLVER_EXPM = 1_ik  !< frozen-coefficient matrix exponential (default)
   integer(ik), parameter :: HYDRO_SOLVER_BE   = 2_ik  !< linearly-implicit backward Euler (3-node / stiff)

   !----- Internal-conductance parameterization. -------------------------------------------!
   integer(ik), parameter :: HYDRO_COND_KPLANT  = 1_ik !< leaf-area-specific whole-plant conductance (default)
   integer(ik), parameter :: HYDRO_COND_SEGMENT = 2_ik !< X16 stem allometry: kmax*sap_area/(height*curl)

   !----- Sub-step control. ----------------------------------------------------------------!
   integer(ik), parameter :: HYDRO_SUBSTEP_ADAPTIVE = 1_ik !< step-doubling error control (default)
   integer(ik), parameter :: HYDRO_SUBSTEP_FIXED    = 2_ik !< fixed n_sub equal steps (GPU lockstep)

   !----- Boundary conditions + per-plant geometry (all read-only; per plant). --------------!
   type :: hydro_env_t
      real(wp) :: transp     = 0.0_wp   !< [kg/s]  transpiration demand E (per plant)
      real(wp) :: soil_psi   = 0.0_wp   !< [MPa]   aggregated soil (rhizosphere) water potential
      real(wp) :: rhizo_cond = 0.0_wp   !< [kg/s/MPa] soil->root conductance (per plant, given BC)
      real(wp) :: bleaf      = 0.0_wp   !< [kgC]   leaf biomass  (sets leaf capacitance)
      real(wp) :: bsap       = 0.0_wp   !< [kgC]   sapwood biomass
      real(wp) :: broot      = 0.0_wp   !< [kgC]   fine-root biomass (bsap+broot set wood capacitance)
      real(wp) :: sap_area   = 0.0_wp   !< [m2]    sapwood cross-sectional area (segment cond. mode)
      real(wp) :: height     = 0.0_wp   !< [m]     plant height (gravity head + segment path length)
      real(wp) :: leaf_area  = 0.0_wp   !< [m2]    leaf area (scales whole-plant conductance)
      !----- Multi-layer root boundary (ED2-style; MEDS_MULTILAYER_ROOTS_DESIGN). When              !
      !       n_root_layer > 1 the solver aggregates these per-layer soil potentials + rhizosphere    !
      !       conductances into an effective (G_root, psi_soil_eff) at the wood node and distributes   !
      !       uptake back per layer; n_root_layer <= 1 uses the scalar soil_psi/rhizo_cond above       !
      !       (bit-identical single-BC path).                                                          !
      integer(ik) :: n_root_layer             = 0_ik    !< active root layers (<= 1 => scalar BC)
      real(wp) :: soil_psi_layer(NROOT_MAX)   = 0.0_wp  !< [MPa]       per-layer soil water potential
      real(wp) :: rhizo_cond_layer(NROOT_MAX) = 0.0_wp  !< [kg/s/MPa]  per-layer soil->root conductance
      real(wp) :: root_z_layer(NROOT_MAX)     = 0.0_wp  !< [m, <=0]    layer depth for gravity head
   end type hydro_env_t

   !----- Flat per-PFT hydraulic trait set (self-contained, filled by the seam from cfg%pft). !
   type :: hydro_params_t
      !----- Pressure-volume (Bartlett/Tyree-Hammel), per tissue. --------------------------!
      real(wp) :: leaf_pi0 = 0.0_wp, leaf_elastic_mod = 0.0_wp, leaf_apoplast_frac = 0.0_wp  !< [MPa],[MPa],[-]
      real(wp) :: wood_pi0 = 0.0_wp, wood_elastic_mod = 0.0_wp, wood_apoplast_frac = 0.0_wp
      real(wp) :: leaf_water_sat = 0.0_wp, wood_water_sat = 0.0_wp        !< [kg H2O / kgC] at saturation
      !----- Xylem vulnerability (loss of conductance). -----------------------------------!
      real(wp) :: wood_psi50 = 0.0_wp   !< [MPa, <0] potential at 50% loss
      real(wp) :: wood_kexp  = 0.0_wp   !< [-]  vulnerability shape (a)
      !----- Conductance parameterization. ------------------------------------------------!
      real(wp) :: k_plant_max = 0.0_wp  !< [kg/s/MPa/m2_leaf] whole-plant (HYDRO_COND_KPLANT)
      real(wp) :: wood_kmax   = 0.0_wp  !< [kg/m/s/MPa] sapwood specific conductivity (HYDRO_COND_SEGMENT)
      real(wp) :: vessel_curl = 1.0_wp  !< [-] tortuosity / path-length factor (HYDRO_COND_SEGMENT)
      !----- Precomputed Kirchhoff lookup table (built from wood_kexp by build_hydro_table). It is    !
      !       consulted on the hot path ONLY for wood_kexp not in {1,2} (the quadrature regime); for   !
      !       kexp in {1,2} the solver keeps the exact closed form, so it stays dormant there. --------!
      type(hydro_table_t) :: vuln_table
   end type hydro_params_t

   !----- Run selectors + numerical controls. ----------------------------------------------!
   type :: hydro_opts_t
      integer(ik) :: topology     = HYDRO_NODES_2
      integer(ik) :: solver       = HYDRO_SOLVER_EXPM
      integer(ik) :: cond_mode    = HYDRO_COND_KPLANT
      integer(ik) :: substep_mode = HYDRO_SUBSTEP_ADAPTIVE
      logical     :: gravity_on   = .true.
      real(wp)    :: rtol         = 1.0e-3_wp   !< [-]  relative sub-step tolerance
      real(wp)    :: atol         = 1.0e-3_wp   !< [MPa] absolute sub-step tolerance
      real(wp)    :: h_init       = 0.0_wp      !< [s]  initial sub-step (0 => start at dt)
      integer(ik) :: max_substep  = 200_ik      !< sub-step cap / fixed count
   end type hydro_opts_t

   !----- Outputs (per plant). -------------------------------------------------------------!
   type :: hydro_flux_t
      real(wp)    :: sapflow     = 0.0_wp   !< [kg/s]  wood->leaf sapflow (time-mean over dt)
      real(wp)    :: root_uptake = 0.0_wp   !< [kg/s]  soil->root uptake (time-mean; the budget term)
      real(wp)    :: root_uptake_layer(NROOT_MAX) = 0.0_wp !< [kg/s] per-layer uptake; sum = root_uptake
      real(wp)    :: psi_leaf    = 0.0_wp   !< [MPa]   leaf water potential at end of step
      real(wp)    :: psi_wood    = 0.0_wp   !< [MPa]   wood water potential at end of step
      real(wp)    :: plc         = 0.0_wp   !< [-]     plant loss of conductance (1 - retained)
      integer(ik) :: nsub        = 0_ik     !< sub-steps taken
      logical     :: converged   = .false.  !< .true. if fully integrated AND every sub-step met tolerance
   end type hydro_flux_t

   !=======================================================================================!
   !     PHENOLOGY -- pure SIGNAL kernel: env cues + traits -> two RELATIVE rate tendencies. !
   !     Emits leaf_flush_rate + leaf_shed_rate [1/day]; touches NO carbon, NO leaf/storage   !
   !     state, NO elongf. All leaf/storage carbon update lives in meds_plant_carbon_dynamics.!
   !     The kernel carries TWO governor accumulators (flush_drive, shed_drive) as its memory.!
   !     See docs/dev_plans/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md.                            !
   !=======================================================================================!
   !----- Cue-enable bits. flush_cue_mask and shed_cue_mask select which cues drive each side  !
   !      (min over flush cues, max over shed cues) -- the mechanism for the four target        !
   !      patterns (evergreen flushes on TEMP but never sheds; a leaf-exchanger flushes          !
   !      permissively but sheds on LIGHT).                                                       !
   integer(ik), parameter :: CUE_NONE  = 0_ik    !< no cues (permissive flush / no active shed)
   integer(ik), parameter :: CUE_TEMP  = 1_ik    !< temperature (GDD flush + autumn cold-drop shed)
   integer(ik), parameter :: CUE_WATER = 2_ik    !< soil-water running mean
   integer(ik), parameter :: CUE_HYDRO = 4_ik    !< daily-max leaf water potential (dmax_leaf_psi)
   integer(ik), parameter :: CUE_PHOTO = 8_ik    !< photoperiod (gates the temperature flush)
   integer(ik), parameter :: CUE_LIGHT = 16_ik   !< light: active shed rises with running-mean radiation

   !----- Raw daily environmental drivers (read-only; NO leaf status). ---------------------!
   type :: pheno_env_t
      real(wp)    :: temp_day      = 0.0_wp    !< [K]   daily-mean air/canopy temperature (thermal sums)
      real(wp)    :: soil_temp     = 0.0_wp    !< [K]   shallow-layer soil temperature (cold-drop trigger)
      real(wp)    :: avail_water   = 0.0_wp    !< [-] fraction OR [MPa] soil-water potential (CUE_WATER)
      real(wp)    :: dmax_leaf_psi = 0.0_wp    !< [MPa, <=0] DAILY-MAX leaf water potential (CUE_HYDRO)
      real(wp)    :: rad           = 0.0_wp    !< [W/m2] daily-mean radiation (CUE_LIGHT)
      real(wp)    :: daylength     = 12.0_wp   !< [h]   photoperiod (CUE_PHOTO; caller-supplied)
      integer(ik) :: doy           = 1_ik      !< [-]   day-of-year (thermal-sum season gating)
      logical     :: hemis_north   = .true.    !< northern hemisphere (season gating)
   end type pheno_env_t

   !----- The TWO governor accumulators (the prognostic memory) + cue sub-accumulators. ----!
   type :: pheno_state_t
      real(wp) :: flush_drive   = 1.0_wp   !< [-] smoothed flush permission in [0,1]; 0 => dormant
      real(wp) :: shed_drive    = 0.0_wp   !< [-] smoothed active-shed pressure in [0,1]; 0 => none
      real(wp) :: gdd           = 0.0_wp   !< [K day] growing-degree-day sum            (CUE_TEMP)
      real(wp) :: chill         = 0.0_wp   !< [day]   chilling-day count                (CUE_TEMP)
      real(wp) :: water_avg     = 0.0_wp   !< [-]|[MPa] running-mean available water    (CUE_WATER)
      real(wp) :: low_psi_days  = 0.0_wp   !< [day]   consecutive dry days (dmax<tlp)   (CUE_HYDRO)
      real(wp) :: high_psi_days = 0.0_wp   !< [day]   consecutive wet days (dmax>=.5tlp)(CUE_HYDRO)
      real(wp) :: light_avg     = 0.0_wp   !< [W/m2]  running-mean radiation            (CUE_LIGHT)
   end type pheno_state_t

   !----- Flat per-PFT trait set (self-contained; filled by the driver from cfg%pft). ------!
   type :: pheno_params_t
      !----- Selectors: the two cue masks + shared logistic sharpness. --------------------!
      integer(ik) :: flush_cue_mask = CUE_NONE   !< OR of CUE_* driving the flush side (min)
      integer(ik) :: shed_cue_mask  = CUE_NONE   !< OR of CUE_* driving the shed side (max)
      real(wp)    :: cue_sharpness   = 2.0_wp     !< [-] dimensionless logistic slope (large => ED2-sharp)
      !----- Per-cue transition WIDTHS (normalize each driver so cue_sharpness is a shared slope). !
      real(wp)    :: gdd_width       = 50.0_wp    !< [K day] GDD flush transition width
      real(wp)    :: daylen_width    = 1.0_wp     !< [h]     autumn daylength transition width
      real(wp)    :: soiltemp_width  = 2.0_wp     !< [K]     autumn soil-temperature transition width
      !----- Relative rate scales [1/day] + governor smoothing timescales [day]. -----------!
      real(wp)    :: k_flush_max = 0.06667_wp    !< [1/day] max relative flush rate (~full in 15 d)
      real(wp)    :: k_shed_max  = 0.05_wp       !< [1/day] max relative active-shed rate (~bare in 20 d)
      real(wp)    :: tau_flush   = 5.0_wp        !< [day]   flush governor low-pass timescale
      real(wp)    :: tau_shed    = 5.0_wp        !< [day]   shed governor low-pass timescale
      !----- Thermal (CUE_TEMP): GDD flush threshold a+b*exp(c*chill) + autumn cold drop. ---!
      real(wp)    :: gdd_base_temp       = 278.15_wp  !< [K] GDD accumulation base (5 degC)
      real(wp)    :: chill_base_temp     = 278.15_wp  !< [K] chilling-day base
      real(wp)    :: phen_a              = -68.0_wp   !< [K day] GDD threshold intercept  (Botta 2000)
      real(wp)    :: phen_b              = 638.0_wp   !< [K day] GDD threshold amplitude
      real(wp)    :: phen_c              = -0.01_wp   !< [1/day] chilling exponent (more chill => lower GDD need)
      real(wp)    :: cold_drop_daylength = 10.9_wp    !< [h] autumn short-day drop trigger  (White 1997)
      real(wp)    :: cold_drop_soiltemp1 = 284.3_wp   !< [K] cool-soil drop (with short days)
      real(wp)    :: cold_drop_soiltemp2 = 275.15_wp  !< [K] very-cold-soil drop (unconditional)
      !----- Water (CUE_WATER): running-mean ramp between off/on thresholds. ----------------!
      logical     :: water_use_potential = .false.    !< .false.: moisture fraction; .true.: soil-psi [MPa]
      real(wp)    :: water_off_threshold = 0.2_wp      !< available water at which shed = 1 / flush = 0
      real(wp)    :: water_on_threshold  = 0.5_wp      !< available water at which flush = 1 / shed = 0 (> off)
      real(wp)    :: water_window        = 10.0_wp     !< [day] running-mean window
      real(wp)    :: water_width         = 0.1_wp      !< transition width for the water logistics
      !----- Hydraulic (CUE_HYDRO): dmax_leaf_psi consecutive-day counters vs the TLP. ------!
      real(wp)    :: leaf_psi_tlp        = -2.0_wp     !< [MPa] turgor-loss point (Xu 2016)
      real(wp)    :: low_psi_threshold   = 10.0_wp     !< [day] dry days to full shed
      real(wp)    :: high_psi_threshold  = 10.0_wp     !< [day] wet days to full flush
      !----- Photoperiod (CUE_PHOTO): daylength logistic gate (multiplies the flush side). --!
      real(wp)    :: photo_crit          = 11.0_wp     !< [h] critical daylength
      real(wp)    :: photo_slope         = 2.0_wp      !< [1/h] daylength logistic slope
      !----- Light (CUE_LIGHT): active shed rises with running-mean radiation. --------------!
      real(wp)    :: light_on_threshold  = 200.0_wp    !< [W/m2] radiation at which the light shed onsets
      real(wp)    :: light_width         = 50.0_wp     !< [W/m2] transition width for the light shed logistic
      real(wp)    :: light_window        = 10.0_wp     !< [day] running-mean window (ED2 rad_avg)
   end type pheno_params_t

   !----- Outputs: the two RELATIVE rate tendencies + the governing shed cue (diagnostic). --!
   type :: pheno_out_t
      real(wp)    :: leaf_flush_rate = 0.0_wp    !< [1/day] relative flush tendency; 0 == dormancy
      real(wp)    :: leaf_shed_rate  = 0.0_wp    !< [1/day] relative active-shed tendency; 0 == no active shed
      integer(ik) :: cue_limiting    = CUE_NONE  !< strongest active shed cue (argmax over shed_cue_mask)
   end type pheno_out_t

   !=======================================================================================!
   !     RESPIRATION -- non-leaf MAINTENANCE respiration (stem + fine root). Per-plant       !
   !     fluxes [umol CO2 / plant / s]; x nplant -> per ground. All rates 25 degC-referenced. !
   !     (Leaf dark respiration Rd is computed in the leaf solver, returned as leaf_flux_t%rd.)!
   !=======================================================================================!
   !----- Woody-tissue (stem) maintenance respiration (ED2 Chambers surface-area form). ----!
   type :: wood_params_t
      logical  :: is_woody              = .true.       !< .false. (e.g. grass) => stem respiration is 0
      real(wp) :: stem_resp_factor25    = 0.0_wp       !< [umol CO2/m2 stem/s @25C] baseline (25C-based; ED2's
                                                       !< 15C Chambers value is converted once at parameter-init)
      real(wp) :: stem_resp_size_scaler = 0.0_wp       !< [1/cm]   DBH size effect (0 => flat; ED2 ~0.0041)
      real(wp) :: agf_bs                = 0.7_wp       !< [--]     aboveground fraction of structural biomass
      real(wp) :: ea                    = 46390.0_wp   !< [J/mol]   peaked-Arrhenius activation energy (leaf ea_rd)
      real(wp) :: hd                    = 200000.0_wp  !< [J/mol]   deactivation energy   (leaf hd_rd)
      real(wp) :: ds                    = 490.0_wp     !< [J/mol/K] entropy term          (leaf ds_rd)
   end type wood_params_t

   !----- Fine-root maintenance respiration (ED2 per-broot form; single effective soil T). --!
   type :: root_params_t
      real(wp) :: root_resp_factor25 = 0.0_wp      !< [umol CO2/kgC fine root/s @25C] (25C-based; = base_rate_per_N
                                                   !< * n_conc; ED2's 15C value converted once at parameter-init)
      real(wp) :: ea                 = 46390.0_wp  !< peaked-Arrhenius terms (default leaf ea_rd/hd_rd/ds_rd)
      real(wp) :: hd                 = 200000.0_wp
      real(wp) :: ds                 = 490.0_wp
   end type root_params_t


   !=======================================================================================!
   !     CARBON ALLOCATION -- see meds_plant_carbon_allocation. The allocation kernel and its !
   !     rate->amount helpers (leaf_shed_amount, flush_growth_cap) are ELEMENTAL over plain     !
   !     cohort scalars, so no derived types live here: the vegetation-dynamics driver hands    !
   !     the kernel raw per-cohort arrays (GPP, maintenance resp, storage, shed amounts, capped  !
   !     demands) and receives the per-pool NPP + growth respiration. Tissue turnover is a       !
   !     degenerate phenology (baseline shed rate) and lives in the PHENOLOGY section above.     !
   !=======================================================================================!

end module meds_plant_types
