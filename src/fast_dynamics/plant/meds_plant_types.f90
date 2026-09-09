!==========================================================================================!
! meds_plant_types -- the derived types of the SUB-DAILY plant kernels: leaf gas exchange,   !
! hydraulics, and non-leaf maintenance respiration. Phenology's types live in meds_phenology_types.!
!                                                                                          !
! Consolidated from the former per-domain meds_{leaf,hydro,pheno}_types modules (the plant   !
! library compiles as a whole, so there is no reason to keep them apart). The types are pure  !
! DATA -- no methods, no hidden state -- organized in clearly delimited sections:             !
!   * LEAF        -- leaf_env_t / leaf_flux_t / leaf_photo_params_t + limitation & pathway flags. !
!   * HYDRAULICS  -- hydro_env_t / hydro_params_t / hydro_opts_t / hydro_flux_t + topology flags.  !
! (PHENOLOGY moved to meds_phenology_types when the plant library split by timescale.)              !
! (Respiration types wood_*/root_* will be added here when meds_plant_respiration lands.)      !
!==========================================================================================!
module meds_plant_types
   use meds_kinds,      only : wp, ik
   use meds_pft_params, only : PATH_C3, PATH_C4
   use meds_hydr_lib, only : hydro_table_t
   implicit none
   private

   !----- LEAF -----------------------------------------------------------------------------!
   public :: leaf_env_t, leaf_flux_t, leaf_photo_params_t, leaf_photo_table_t
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
   !----- RESPIRATION (env/flux types removed: the kernels are now `elemental pure` over bare      !
   !      scalars, MEDS_NUMERICS_SCOPING.md §11; only the run-uniform trait PODs remain). ----------!
   public :: wood_params_t, root_params_t
   public :: veg_thermal_params_t
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

   !----- leaf_photo_table_t -- the leaf-photosynthesis parameters of EVERY PFT, assembled ONCE at   !
   !      configuration time (meds_fast_config%build_leaf_photo_table) instead of re-flattened  !
   !      from the ~45 PFT trait arrays for every leaf on every fast step (2026-09 review, item 4   !
   !      #8). `pft(i)` carries PFT i's parameters with the TABLE's Vcmax25/Jmax25/TPU25/Rd25; a    !
   !      cohort's plastic capacities override those per leaf, and the two ratios are what scale   !
   !      Jmax25 and TPU25 with the overriding Vcmax25. The run-level solver selectors ride along  !
   !      so the batch kernel needs nothing else from the configuration. --------------------------!
   type :: leaf_photo_table_t
      integer(ik) :: n_pft = 0_ik
      type(leaf_photo_params_t), allocatable :: pft(:)              !< per-PFT parameters (table capacities)
      real(wp),                  allocatable :: jmax_vcmax_ratio(:) !< [--] Jmax25 / Vcmax25 per PFT
      real(wp),                  allocatable :: tpu_vcmax_ratio(:)  !< [--] TPU25 / Vcmax25 per PFT
      integer(ik) :: stomatal_model     = 0_ik    !< SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form = 0_ik    !< TRESP_ARRHENIUS | TRESP_PEAKED
      integer(ik) :: colimitation       = 0_ik    !< COLIM_MIN | COLIM_QUADRATIC
      logical     :: use_boundary_layer = .true.  !< couple through the leaf boundary layer (gb)
   end type leaf_photo_table_t

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
   !     RESPIRATION -- non-leaf MAINTENANCE respiration (stem + fine root). Per-plant       !
   !     fluxes [umol CO2 / plant / s]; x nplant -> per ground. All rates 25 degC-referenced. !
   !     (Leaf dark respiration Rd is computed in the leaf solver, returned as leaf_flux_t%rd.)!
   !=======================================================================================!
   !----- Woody-tissue (stem) maintenance respiration (ED2 Chambers surface-area form). ----!
   type :: wood_params_t
      real(wp) :: stem_resp_size_scaler = 0.0_wp       !< [1/cm]   DBH size effect (0 => flat; ED2 ~0.0041)
      !----- agf_bs is GONE (issue #128): the aboveground fraction is a per-PFT trait            !
      !      (`aboveground_frac`), gathered per cohort and passed to the kernel, not a run       !
      !      constant duplicated here.  ---------------------------------------------------------!
      real(wp) :: ea                    = 46390.0_wp   !< [J/mol]   peaked-Arrhenius activation energy (leaf ea_rd)
      real(wp) :: hd                    = 200000.0_wp  !< [J/mol]   deactivation energy   (leaf hd_rd)
      real(wp) :: ds                    = 490.0_wp     !< [J/mol/K] entropy term          (leaf ds_rd)
   end type wood_params_t

   !----- Fine-root maintenance respiration (ED2 per-broot form; single effective soil T). --!
   type :: root_params_t
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

   !----- (soil_energy_column_t + cas_state_t now live in meds_column_state_types; re-exported.) -!

   !----- soil_thermal_params_t (per-column soil thermal texture) is defined in (and re-exported !
   !      from) meds_column_state_types; the conductivity/heat-capacity kernels are in meds_therm_lib.!

   !----- Per-PFT vegetation thermal parameters. -------------------------------------------!
   type :: veg_thermal_params_t
      real(wp) :: leaf_emiss     = 0.95_wp                  !< [-] LW emissivity (Jacobian -8*eps*sigma*T^3 term)
      real(wp) :: effarea_heat   = 2.0_wp                   !< [-] sensible sidedness (both leaf sides)
      real(wp) :: effarea_evap   = 1.0_wp                   !< [-] film-evaporation sidedness
      real(wp) :: effarea_transp = 1.0_wp                   !< [-] transpiration sidedness (per PFT)
      real(wp) :: veg_hcap_min   = 20.0_wp                  !< [J/m2/K] resolvability floor
      real(wp) :: c_leaf = 3200.0_wp, c_sapw = 2700.0_wp    !< [J/kg/K] tissue specific heats
   end type veg_thermal_params_t

end module meds_plant_types
