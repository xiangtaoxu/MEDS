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
   implicit none
   private

   !----- LEAF -----------------------------------------------------------------------------!
   public :: leaf_env_t, leaf_flux_t, leaf_photo_params_t
   public :: PATH_C3, PATH_C4                       ! re-export (pathway lives with the PFT traits)
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP
   !----- HYDRAULICS -----------------------------------------------------------------------!
   public :: hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t
   public :: N_HYDRO, NODE_LEAF, NODE_STEM, NODE_ROOT, NODE_WOOD
   public :: HYDRO_NODES_2, HYDRO_NODES_3
   public :: HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE
   public :: HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT
   public :: HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED
   !----- PHENOLOGY ------------------------------------------------------------------------!
   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO
   public :: PHEN_ON, PHEN_DORMANT, PHEN_OFF
   !----- RESPIRATION ----------------------------------------------------------------------!
   public :: wood_env_t, wood_params_t, wood_flux_t
   public :: root_env_t, root_params_t, root_flux_t
   !----- CARBON DYNAMICS ------------------------------------------------------------------!
   public :: turnover_env_t, turnover_params_t, turnover_rates_t
   public :: carbon_gain_t, carbon_loss_t, carbon_demand_t, carbon_npp_t, carbon_env_t

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
      real(wp) :: psi_leaf   !< [MPa]  leaf water potential (<= 0); drives the water-stress factor
      real(wp) :: gb         !< [mol H2O/m2/s] boundary-layer conductance (<= 0 => skip, Cs = Ca)
   end type leaf_env_t

   !----- Leaf-level fluxes returned by the solver. ---------------------------------------!
   type :: leaf_flux_t
      real(wp)    :: a_net        !< [umol CO2/m2/s] net assimilation (A_gross - Rd)
      real(wp)    :: a_gross      !< [umol CO2/m2/s] gross assimilation
      real(wp)    :: gs           !< [mol H2O/m2/s]  stomatal conductance to water vapour
      real(wp)    :: ci           !< [umol/mol]      intercellular CO2 mole fraction
      real(wp)    :: cs           !< [umol/mol]      leaf-surface CO2 mole fraction (= ca if no BL)
      real(wp)    :: transpiration!< [mol H2O/m2/s]  E = gs * VPD / pressure
      real(wp)    :: rd           !< [umol/m2/s]     leaf respiration used
      integer(ik) :: limitation   !< LIM_* of the binding term
      logical     :: converged    !< .true. if the Ci solve met tolerance
   end type leaf_flux_t

   !----- Flat per-PFT parameter set (per-PFT traits + shared biochemistry), self-contained !
   !       so the solver never references meds_config. Filled by the interface from the config.!
   type :: leaf_photo_params_t
      integer(ik) :: pathway        !< PATH_C3 | PATH_C4
      !----- Per-PFT capacities at 25 degC and stomatal/water-stress traits. ---------------!
      real(wp) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(wp) :: g0, g1, d0, quantum_yield, theta_j, theta_cj, theta_ic
      real(wp) :: lambda25, psi_open, psi_close, lambda_psi_exp
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
   end type hydro_env_t

   !----- Flat per-PFT hydraulic trait set (self-contained, filled by the seam from cfg%pft). !
   type :: hydro_params_t
      !----- Pressure-volume (Bartlett/Tyree-Hammel), per tissue. --------------------------!
      real(wp) :: leaf_pi0 = 0.0_wp, leaf_eps = 0.0_wp, leaf_af = 0.0_wp  !< [MPa],[MPa],[-]
      real(wp) :: wood_pi0 = 0.0_wp, wood_eps = 0.0_wp, wood_af = 0.0_wp
      real(wp) :: leaf_water_sat = 0.0_wp, wood_water_sat = 0.0_wp        !< [kg H2O / kgC] at saturation
      !----- Xylem vulnerability (loss of conductance). -----------------------------------!
      real(wp) :: wood_psi50 = 0.0_wp   !< [MPa, <0] potential at 50% loss
      real(wp) :: wood_kexp  = 0.0_wp   !< [-]  vulnerability shape (a)
      !----- Conductance parameterization. ------------------------------------------------!
      real(wp) :: k_plant_max = 0.0_wp  !< [kg/s/MPa/m2_leaf] whole-plant (HYDRO_COND_KPLANT)
      real(wp) :: wood_kmax   = 0.0_wp  !< [kg/m/s/MPa] sapwood specific conductivity (HYDRO_COND_SEGMENT)
      real(wp) :: vessel_curl = 1.0_wp  !< [-] tortuosity / path-length factor (HYDRO_COND_SEGMENT)
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
      real(wp)    :: psi_leaf    = 0.0_wp   !< [MPa]   leaf water potential at end of step
      real(wp)    :: psi_wood    = 0.0_wp   !< [MPa]   wood water potential at end of step
      real(wp)    :: plc         = 0.0_wp   !< [-]     plant loss of conductance (1 - retained)
      integer(ik) :: nsub        = 0_ik     !< sub-steps taken
      logical     :: converged   = .false.  !< .true. if the sub-step cap was not hit
   end type hydro_flux_t

   !=======================================================================================!
   !     PHENOLOGY -- pure directional leaf-phenology interface seam.                       !
   !=======================================================================================!
   !----- Cue-enable mask bits (per-PFT, OR-combined). The ONLY strategy selector. ---------!
   integer(ik), parameter :: CUE_NONE  = 0_ik   !< no cues => evergreen (perpetually ON)
   integer(ik), parameter :: CUE_TEMP  = 1_ik   !< temperature (GDD/CDD cold-deciduous)
   integer(ik), parameter :: CUE_WATER = 2_ik   !< soil-water drought (running-mean available water)
   integer(ik), parameter :: CUE_HYDRO = 4_ik   !< leaf water potential (hydraulic)
   integer(ik), parameter :: CUE_PHOTO = 8_ik   !< photoperiod (daylength gate on the temperature cue)

   !----- Phenological status: the direction of leaf-display change. -----------------------!
   integer(ik), parameter :: PHEN_OFF     = -1_ik  !< unfavorable -- actively dropping leaves
   integer(ik), parameter :: PHEN_DORMANT =  0_ik  !< neutral deadband -- hold the current state
   integer(ik), parameter :: PHEN_ON      =  1_ik  !< favorable -- actively seeking leaf growth

   !----- Raw daily environmental drivers (read-only; all a caller boundary condition). ----!
   type :: pheno_env_t
      real(wp)    :: temp_day    = 0.0_wp    !< [K]   daily-mean air/canopy temperature (thermal sums)
      real(wp)    :: soil_temp   = 0.0_wp    !< [K]   shallow-layer soil temperature (cold-drop trigger)
      real(wp)    :: avail_water = 0.0_wp    !< [-] fraction OR [MPa] soil-water potential (CUE_WATER)
      real(wp)    :: psi_leaf    = 0.0_wp    !< [MPa, <=0] leaf water potential (CUE_HYDRO)
      real(wp)    :: daylength   = 12.0_wp   !< [h]   photoperiod (CUE_PHOTO; caller-supplied)
      integer(ik) :: doy         = 1_ik      !< [-]   day-of-year (thermal-sum season gating)
      logical     :: hemis_north = .true.    !< northern hemisphere (season gating)
   end type pheno_env_t

   !----- Cue accumulators: the prognostic phenological MEMORY (not leaf-mass state). ------!
   type :: pheno_state_t
      real(wp) :: gdd           = 0.0_wp   !< [K day] growing-degree-day sum          (CUE_TEMP)
      real(wp) :: chill         = 0.0_wp   !< [day]   chilling-day count              (CUE_TEMP)
      real(wp) :: water_avg     = 0.0_wp   !< [-]|[MPa] running-mean available water  (CUE_WATER)
      real(wp) :: low_psi_days  = 0.0_wp   !< [day]   consecutive days psi_leaf < psi_tlp   (CUE_HYDRO)
      real(wp) :: high_psi_days = 0.0_wp   !< [day]   consecutive days psi_leaf >= 0.5 psi_tlp (CUE_HYDRO)
   end type pheno_state_t

   !----- Flat per-PFT trait set (self-contained; filled by the seam from cfg%pft). --------!
   type :: pheno_params_t
      !----- Selector: the cue mask, the shared logistic sharpness, the on/off band. --------!
      integer(ik) :: cue_mask           = CUE_NONE   !< OR of CUE_* (evergreen = CUE_NONE)
      real(wp)    :: cue_sharpness       = 2.0_wp     !< [-] dimensionless logistic slope (large => ED2-sharp)
      real(wp)    :: phen_on_threshold   = 0.6_wp     !< favorability above which status = ON
      real(wp)    :: phen_off_threshold  = 0.4_wp     !< favorability below which status = OFF (< on => DORMANT band)
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
      real(wp)    :: water_off_threshold = 0.2_wp      !< available water at which favorability = 0
      real(wp)    :: water_on_threshold  = 0.5_wp      !< available water at which favorability = 1 (> off)
      real(wp)    :: water_window        = 10.0_wp     !< [day] running-mean window
      !----- Hydraulic (CUE_HYDRO): leaf-psi consecutive-day counters vs the TLP. -----------!
      real(wp)    :: leaf_psi_tlp        = -2.0_wp     !< [MPa] turgor-loss point (Xu 2016)
      real(wp)    :: low_psi_threshold   = 10.0_wp     !< [day] dry days to full unfavorable
      real(wp)    :: high_psi_threshold  = 10.0_wp     !< [day] wet days to full favorable
      !----- Photoperiod (CUE_PHOTO): daylength logistic gate. ------------------------------!
      real(wp)    :: photo_crit          = 11.0_wp     !< [h] critical daylength
      real(wp)    :: photo_slope         = 2.0_wp      !< [1/h] daylength logistic slope
   end type pheno_params_t

   !----- Outputs: the phenological status + the governing (most-limiting) cue. ------------!
   type :: pheno_out_t
      integer(ik) :: phenology_status = PHEN_DORMANT   !< PHEN_ON | PHEN_OFF | PHEN_DORMANT
      integer(ik) :: cue_limiting     = CUE_NONE       !< the CUE_* bit with the lowest favorability
   end type pheno_out_t

   !=======================================================================================!
   !     RESPIRATION -- non-leaf MAINTENANCE respiration (stem + fine root). Per-plant       !
   !     fluxes [umol CO2 / plant / s]; x nplant -> per ground. All rates 25 degC-referenced. !
   !     (Leaf dark respiration Rd is computed in the leaf solver, returned as leaf_flux_t%rd.)!
   !=======================================================================================!
   !----- Woody-tissue (stem) maintenance respiration (ED2 Chambers surface-area form). ----!
   type :: wood_env_t
      real(wp) :: wood_temp = 0.0_wp    !< [K]  woody-tissue temperature; drives the T-response
      real(wp) :: dbh       = 0.0_wp    !< [cm] stem diameter at breast height
      real(wp) :: height    = 0.0_wp    !< [m]  cohort height
      real(wp) :: wai       = 0.0_wp    !< [m2/m2 ground] wood area index (0 => cylinder-only stem area)
      real(wp) :: nplant    = 1.0_wp    !< [plant/m2] stem density (converts wai -> per-plant branch area)
      real(wp) :: t_acclim  = 0.0_wp    !< [K]  running-mean temperature for acclimation (RESERVED; unused v1)
   end type wood_env_t

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

   type :: wood_flux_t
      real(wp) :: stem_resp = 0.0_wp    !< [umol CO2 / plant / s] stem MAINTENANCE respiration (per plant)
   end type wood_flux_t

   !----- Fine-root maintenance respiration (ED2 per-broot form; single effective soil T). --!
   type :: root_env_t
      real(wp) :: soil_temp = 0.0_wp    !< [K]        effective (root-weighted mean) soil temperature
      real(wp) :: broot     = 0.0_wp    !< [kgC/plant] fine-root biomass (carbon); broot=0 => 0
      real(wp) :: t_acclim  = 0.0_wp    !< [K]        running-mean soil temperature for acclimation (RESERVED)
   end type root_env_t

   type :: root_params_t
      real(wp) :: root_resp_factor25 = 0.0_wp      !< [umol CO2/kgC fine root/s @25C] (25C-based; = base_rate_per_N
                                                   !< * n_conc; ED2's 15C value converted once at parameter-init)
      real(wp) :: ea                 = 46390.0_wp  !< peaked-Arrhenius terms (default leaf ea_rd/hd_rd/ds_rd)
      real(wp) :: hd                 = 200000.0_wp
      real(wp) :: ds                 = 490.0_wp
   end type root_params_t

   type :: root_flux_t
      real(wp) :: root_resp = 0.0_wp    !< [umol CO2 / plant / s] fine-root MAINTENANCE respiration (per plant)
   end type root_flux_t

   !=======================================================================================!
   !     CARBON DYNAMICS -- stateless per-cohort carbon budget + allocation. Every quantity  !
   !     is CARBON [kgC/plant]; all biomass<->carbon conversion is done ONCE at parameter      !
   !     initialization, so the kernels never convert. Two concerns handled by                 !
   !     meds_plant_carbon_dynamics: tissue turnover RATES (tissue_turnover_rates) and the      !
   !     allocation of the step's carbon gain among the pools (plant_carbon_allocation). Pools  !
   !     are carbon-explicit -- leaf_carbon / fineroot_carbon / wood_carbon / nonstructural --  !
   !     with wood_carbon the size anchor (dbh derived) and sapwood/heartwood diagnostic; see   !
   !     archive/MEDS_PLANT_CARBON_DYNAMICS_DESIGN.md.                                          !
   !=======================================================================================!
   !----- Turnover: driver, per-PFT baseline rates, and the returned (possibly cold-modified) !
   !      rates. Constant for now; the seam exists so a future light/tropical-phenology form    !
   !      can vary the leaf rate without touching allocation.                                   !
   type :: turnover_env_t
      real(wp) :: tissue_temp = 298.15_wp   !< [K]   tissue temperature (evergreen cold suppression)
   end type turnover_env_t

   type :: turnover_params_t
      real(wp) :: leaf_turnover_rate     = 0.0_wp   !< [1/yr] baseline leaf turnover
      real(wp) :: fineroot_turnover_rate = 0.0_wp   !< [1/yr] baseline fine-root turnover
      logical  :: evergreen              = .false.  !< .true. => cold-suppress turnover (ED2 evergreen form)
   end type turnover_params_t

   type :: turnover_rates_t
      real(wp) :: leaf     = 0.0_wp   !< [1/yr] effective leaf turnover rate
      real(wp) :: fineroot = 0.0_wp   !< [1/yr] effective fine-root turnover rate
   end type turnover_rates_t

   !----- Allocation inputs (gains / losses / demands) + outputs (net npp per pool). All are   !
   !      amounts integrated over ONE step [kgC/plant] (rate x pool x dt is done by the caller);  !
   !      the kernel is pure distribution arithmetic.                                             !
   type :: carbon_gain_t
      real(wp) :: net_carbon = 0.0_wp   !< [kgC/plant] NPP after growth respiration (may be < 0)
      real(wp) :: storage    = 0.0_wp   !< [kgC/plant] nonstructural carbon available to draw (>= 0)
   end type carbon_gain_t

   type :: carbon_loss_t
      real(wp) :: leaf     = 0.0_wp   !< [kgC/plant] leaf turnover this step (-> litter)
      real(wp) :: fineroot = 0.0_wp   !< [kgC/plant] fine-root turnover this step (-> litter)
   end type carbon_loss_t

   type :: carbon_demand_t
      real(wp) :: leaf     = 0.0_wp   !< [kgC/plant] deficit toward the leaf target (>= 0)
      real(wp) :: fineroot = 0.0_wp   !< [kgC/plant] deficit toward the fine-root target (>= 0)
      real(wp) :: storage  = 0.0_wp   !< [kgC/plant] deficit toward the storage target (>= 0)
      real(wp) :: wood     = 0.0_wp   !< [kgC/plant] structural-growth demand (residual sink; large => take all)
      real(wp) :: reproduction_fraction = 0.0_wp !< [--] fraction of the post-storage residual -> reproduction
   end type carbon_demand_t

   type :: carbon_npp_t
      real(wp) :: leaf          = 0.0_wp   !< [kgC/plant] NET leaf change (allocation - turnover; signed)
      real(wp) :: fineroot      = 0.0_wp   !< [kgC/plant] NET fine-root change (signed)
      real(wp) :: wood          = 0.0_wp   !< [kgC/plant] wood (structural) growth (>= 0)
      real(wp) :: nonstructural = 0.0_wp   !< [kgC/plant] NET storage change (refill - drawdown; signed)
      real(wp) :: repro         = 0.0_wp   !< [kgC/plant] carbon allocated to reproduction (>= 0; -> recruits)
      real(wp) :: deficit       = 0.0_wp   !< [kgC/plant] unpaid respiration after storage exhausted (>= 0)
      logical  :: starving      = .false.  !< .true. => storage could not cover the carbon debt
   end type carbon_npp_t

   !----- Per-cohort carbon state + drivers for the get_plant_flux_slow seam (the seam derives  !
   !      the turnover losses from the pools + params, then calls plant_carbon_allocation). -----!
   type :: carbon_env_t
      real(wp) :: net_carbon      = 0.0_wp     !< [kgC/plant] NPP after growth resp this step (may be < 0)
      real(wp) :: nonstructural   = 0.0_wp     !< [kgC/plant] current storage (available to draw)
      real(wp) :: leaf_carbon     = 0.0_wp     !< [kgC/plant] current leaf (for turnover)
      real(wp) :: fineroot_carbon = 0.0_wp     !< [kgC/plant] current fine root (for turnover)
      real(wp) :: tissue_temp     = 298.15_wp  !< [K] tissue temperature (evergreen turnover suppression)
      real(wp) :: dt_yr           = 0.0_wp     !< [yr] step length (turnover amount = rate*pool*dt)
      integer(ik) :: phenology_status = 0_ik   !< PHEN_ON | PHEN_OFF | PHEN_DORMANT
   end type carbon_env_t

end module meds_plant_types
