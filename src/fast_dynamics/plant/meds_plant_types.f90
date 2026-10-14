! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_types -- the derived types of the SUB-DAILY plant kernels: leaf gas exchange,   !
! hydraulics, and non-leaf maintenance respiration. Phenology's types live in meds_phenology_types.!
!                                                                                          !
! Consolidated from the former per-domain meds_{leaf,hydro,pheno}_types modules (the plant   !
! library compiles as a whole, so there is no reason to keep them apart). The types are pure  !
! DATA -- no methods, no hidden state -- organized in clearly delimited sections:             !
!   * LEAF        -- leaf_env_t / leaf_flux_t / leaf_photo_params_t + limitation & pathway flags. !
!   * HYDRAULICS  -- hydro_params_t / hydro_params_table_t (per-PFT traits), hydro_opts_t.          !
! (PHENOLOGY moved to meds_phenology_types when the plant library split by timescale.)              !
! (Respiration types wood_*/root_* will be added here when meds_plant_respiration lands.)      !
!==========================================================================================!
module meds_plant_types
   use meds_kinds,      only : wp, ik
   use meds_pft_params, only : PATH_C3, PATH_C4
   use meds_hydr_lib, only : hydro_table_t
   use meds_water_retention, only : water_curve_t
   implicit none
   private

   !----- LEAF -----------------------------------------------------------------------------!
   public :: leaf_env_t, leaf_flux_t, leaf_photo_params_t, leaf_photo_table_t
   public :: PATH_C3, PATH_C4                       ! re-export (pathway lives with the PFT traits)
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP
   !----- HYDRAULICS -----------------------------------------------------------------------!
   public :: hydro_params_t, hydro_opts_t, hydro_params_table_t
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
      real(wp) :: psi_leaf   !< [MPa]  current leaf water potential (<= 0); drives beta_nonstomata (fast)
      real(wp) :: gb         !< [mol H2O/m2/s] boundary-layer conductance (<= 0 => skip, Cs = Ca)
      !----- THE WATER-STATUS POTENTIAL DRIVING THE STOMATAL LIMB (<= 0; default 0 = well-watered).  !
      !                                                                                          !
      !      Deliberately just `psi`, not `psi_soil`. Sabot et al. (2022, Eq. 5) write the stomatal     !
      !      limb as beta_stomata = 1 while the soil is at field capacity and exp(s_ref * psi) below,   !
      !      with psi the SOIL (equivalently predawn) water                                              !
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
      !      which any single soil layer's psi carries. Because that head is there even in wet soil,    !
      !      the stress begins only below an onset (leaf_photo_params_t%psi_onset), which stands in for !
      !      Sabot's field-capacity threshold.                                                          !
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
      real(wp)    :: beta_stomata    = 1.0_wp  !< [-] stomatal limb min(1, exp(sref*(psi-psi_onset))) (1 = unstressed)
      real(wp)    :: beta_nonstomata = 1.0_wp  !< [-] fast limb, the psi_leaf ramp on Vcmax/Jmax/TPU
   end type leaf_flux_t

   !----- Flat per-PFT parameter set (per-PFT traits + shared biochemistry), self-contained !
   !       so the solver never references meds_config. Filled by the interface from the config.!
   type :: leaf_photo_params_t
      integer(ik) :: pathway        !< PATH_C3 | PATH_C4
      !----- Per-PFT capacities at 25 degC and stomatal/water-stress traits. ---------------!
      real(wp) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(wp) :: g0, g1, d0, quantum_yield
      !----- theta_j is the ELECTRON-TRANSPORT hyperbola curvature; the theta_*_c3/c4 pairs are the  !
      !      CO-LIMITATION curvatures, one pair per pathway. The names carry the pathway because C3   !
      !      used to borrow theta_j for both of its smoothings, which is #118: a parameter named for  !
      !      one process doing the work of another, and costing ~29 % of assimilation while doing it. !
      real(wp) :: theta_j
      real(wp) :: theta_cj_c3, theta_ip_c3
      real(wp) :: theta_cj_c4, theta_ic_c4
      !----- psi_open/psi_close [MPa]: the current leaf potential at which the fast limb          !
      !      beta_nonstomata is 1 and 0. A MEDS run sets them per PFT (pft_wstress_psi_open/_close:   !
      !      the turgor-loss point and twice it unless the PFT file gives them). ------------------!
      real(wp) :: lambda25, psi_open, psi_close, lambda_psi_exp, sref_stomata
      !----- Predawn leaf potential [MPa] below which the stomatal water stress begins: beta_stomata  !
      !      is 1 above it and exp(sref_stomata*(psi - psi_onset)) below. A MEDS run sets it per PFT   !
      !      (pft_stomata_psi_onset: half the turgor-loss point unless the PFT file gives it); the 0    !
      !      here, stress from any negative potential, is what a caller building the record by hand   !
      !      gets. --------------------------------------------------------------------------------!
      real(wp) :: psi_onset = 0.0_wp
      !----- Shared biochemistry constants at 25 degC + activation/deactivation terms. -----!
      real(wp) :: kc25, ko25, gstar25
      real(wp) :: ea_kc, ea_ko, ea_gstar, ea_vcmax, ea_jmax, ea_rd
      real(wp) :: hd_vcmax, hd_jmax, hd_rd, ds_vcmax, ds_jmax, ds_rd
      real(wp) :: o2_mol_frac, absorptance, phi_psii
      real(wp) :: medlyn_vpd_min = 0.05_wp   !< [kPa] the VPD the Medlyn model uses at least ([leaf_physiology])
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
   !=======================================================================================!
   !     HYDRAULICS -- the per-PFT plant-water traits and the one run switch. How water moves  !
   !     is the soil-plant water solver's (meds_spac_solver).                                 !
   !=======================================================================================!
   !----- Flat per-PFT hydraulic trait set (self-contained, filled by the seam from cfg%pft). !
   type :: hydro_params_t
      !----- Pressure-volume (Bartlett/Tyree-Hammel), one curve per tissue. ------------------!
      type(water_curve_t) :: leaf_curve, wood_curve
      !----- Xylem vulnerability (loss of conductance). -----------------------------------!
      real(wp) :: wood_psi50 = 0.0_wp   !< [MPa, <0] potential at 50% loss
      real(wp) :: wood_kexp  = 0.0_wp   !< [-]  vulnerability shape (a)
      !----- Conductance: wood_kmax * sapwood area / (height * vessel_curl). -------------------!
      real(wp) :: wood_kmax   = 0.0_wp  !< [kg/m/s/MPa] sapwood specific conductivity
      real(wp) :: vessel_curl = 1.0_wp  !< [-] tortuosity / path-length factor
      !----- Roots: the depth allometry, the profile within it, and the absorbing length. ----------!
      real(wp) :: root_beta            = 0.1_wp       !< [-] (1 - beta^(d/D)) / (1 - beta) above depth d
      real(wp) :: root_depth_b1        = 1.1140580_wp !< [m] rooting depth = b1 * height^b2
      real(wp) :: root_depth_b2        = 0.4223014_wp !< [-]
      real(wp) :: specific_root_length = 2.0e4_wp     !< [m/kgC] fine-root length per unit fine-root carbon
      real(wp) :: fine_root_radius     = 2.5e-4_wp    !< [m] fine-root radius
      !----- Precomputed Kirchhoff lookup table (built from wood_kexp by build_hydro_table). It is    !
      !       consulted on the hot path ONLY for wood_kexp not in {1,2} (the quadrature regime); for   !
      !       kexp in {1,2} the solver keeps the exact closed form, so it stays dormant there. --------!
      type(hydro_table_t) :: vuln_table
   end type hydro_params_t

   !----- PER-PFT hydraulic parameter table (#179), the hydraulic twin of leaf_photo_table_t.     !
   !      Built ONCE per run (each entry carries its own precomputed Kirchhoff lookup, which is    !
   !      what makes per-PFT vulnerability shapes affordable on the hot path) and indexed by the   !
   !      cohort's PFT where the parameters are consumed. Until now one struct served every PFT,   !
   !      so wood density was the only axis on which PFTs could differ hydraulically.              !
   type :: hydro_params_table_t
      integer(ik) :: n_pft = 0_ik
      type(hydro_params_t), allocatable :: pft(:)
   end type hydro_params_table_t

   !----- Run switch: whether the leaves sit a hydrostatic head above the wood. --------------!
   type :: hydro_opts_t
      logical :: gravity_on = .true.
   end type hydro_opts_t


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
   !     CARBON ALLOCATION -- see meds_plant_carbon_allocation. The allocation kernel is      !
   !     ELEMENTAL over plain cohort scalars, so no derived types live here: the vegetation-     !
   !     dynamics driver hands the kernel raw per-cohort arrays (GPP, maintenance resp, storage, !
   !     shed amounts, capped demands) and receives the per-pool NPP + growth respiration. Leaf  !
   !     loss and the flush cap come from meds_phenology::leaf_turnover_step.                    !
   !=======================================================================================!

   !----- (soil_energy_column_t + cas_state_t now live in meds_column_state_types; re-exported.) -!

   !----- soil_thermal_params_t (per-column soil thermal texture) is defined in (and re-exported !
   !      from) meds_column_state_types; the conductivity/heat-capacity kernels are in meds_therm_lib.!

   !----- Per-PFT vegetation thermal parameters. -------------------------------------------!
   type :: veg_thermal_params_t
      real(wp) :: effarea_heat   = 2.0_wp                   !< [-] sensible sidedness (both leaf sides)
      real(wp) :: effarea_evap   = 1.0_wp                   !< [-] film-evaporation sidedness
      real(wp) :: effarea_transp = 1.0_wp                   !< [-] transpiration sidedness (per PFT)
      real(wp) :: veg_hcap_min   = 20.0_wp                  !< [J/m2/K] resolvability floor
      real(wp) :: c_leaf = 3200.0_wp, c_sapw = 2700.0_wp    !< [J/kg/K] tissue specific heats
   end type veg_thermal_params_t

end module meds_plant_types
