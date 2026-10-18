! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_types -- the derived types of the SUB-DAILY plant kernels: leaf gas exchange, !
! hydraulics, non-leaf maintenance respiration and the vegetation thermal parameters.      !
! Phenology's types live in meds_phenology_types.                                          !
!                                                                                          !
! The types are pure DATA -- no methods, no hidden state -- organized in clearly delimited !
! sections:                                                                                !
!   * LEAF        -- leaf_env_t / leaf_flux_t / leaf_photo_params_t / leaf_photo_table_t + !
!                    limitation & pathway flags.                                           !
!   * HYDRAULICS  -- hydro_params_t / hydro_params_table_t (per-PFT traits), hydro_opts_t. !
!   * RESPIRATION -- wood_params_t / root_params_t (run-uniform trait PODs).               !
!   * THERMAL     -- veg_thermal_params_t.                                                 !
!==========================================================================================!
module meds_plant_types
   use meds_kinds,      only : wp, ik
   use meds_pft_params, only : PATH_C3, PATH_C4
   use meds_hydr_lib, only : hydro_table_t
   use meds_water_retention, only : water_curve_t
   use meds_column_params, only : TISSUE_C_LEAF, TISSUE_C_WOOD_DRY, TISSUE_HCAP_MIN
   implicit none
   private

   !----- LEAF -----------------------------------------------------------------------------!
   public :: leaf_env_t, leaf_flux_t, leaf_photo_params_t, leaf_photo_table_t
   public :: PATH_C3, PATH_C4                       ! re-export (pathway lives with the PFT traits)
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP
   !----- HYDRAULICS -----------------------------------------------------------------------!
   public :: hydro_params_t, hydro_opts_t, hydro_params_table_t
   !----- RESPIRATION: the kernels are `elemental pure` over bare scalars (MEDS_NUMERICS_SCOPING.md  !
   !      §11), so only the run-uniform trait PODs live here. ---------------------------------------!
   public :: wood_params_t, root_params_t
   !----- THERMAL ----------------------------------------------------------------------------!
   public :: veg_thermal_params_t

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
      !----- THE WATER-STATUS POTENTIAL DRIVING THE STOMATAL LIMB (<= 0; default 0 = well-watered),  !
      !      docs/science/leaf_gas_exchange.md §4.2. Deliberately just `psi`, not `psi_soil`: it is  !
      !      the plant's water SUPPLY term (Sabot et al. 2022, Eq. 5 write it as the soil or predawn !
      !      potential), as opposed to `psi_leaf` above, the instantaneous demand-side tension that  !
      !      drives the separate NON-stomatal (capacity) limb. MEDS passes `dmax_psi_leaf`, the      !
      !      cohort's own PREDAWN (previous-day daily-maximum) LEAF potential: the quantity the      !
      !      field measurement behind s_ref is, and the plant's own integration over its rooted      !
      !      profile (rooting depth, per-PFT vulnerability, gravity head). It equals the soil        !
      !      potential only in WET soil; under drought a cohort never re-equilibrates overnight, so  !
      !      `s_ref` is calibrated against whichever quantity is passed -- hence the neutral name. A !
      !      caller with no predawn history (a recruit, or the first step of a run) is seeded from   !
      !      the surface-layer soil potential. The stress begins only below                          !
      !      leaf_photo_params_t%psi_onset, which stands in for Sabot's field-capacity threshold     !
      !      because the gravity head is there even in wet soil.                                     !
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
      !----- WATER-STRESS LIMBS, surfaced as DIAGNOSTICS (they do not feed back into anything here   !
      !      -- the kernel has already applied them). A stress closure whose two limbs are never     !
      !      observable can stay silently inert (env%psi left at its default 0 makes beta_stomata    !
      !      identically 1); reporting them makes that failure visible in the output file instead of !
      !      only in a code read.                                                                    !
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
      !----- theta_j is the ELECTRON-TRANSPORT hyperbola curvature; the theta_*_c3/c4 pairs are the   !
      !      CO-LIMITATION curvatures, one pair per pathway, and the names carry the pathway because  !
      !      the two kinds are different quantities (docs/science/leaf_gas_exchange.md §2, "The       !
      !      co-limitation curvatures are not theta_J"). ---------------------------------------------!
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

   !----- leaf_photo_table_t -- the leaf-photosynthesis parameters of EVERY PFT, assembled ONCE at  !
   !      configuration time (meds_fast_config%build_leaf_photo_table), not re-flattened from the   !
   !      ~45 PFT trait arrays for every leaf on every fast step. `pft(i)` carries PFT i's          !
   !      parameters with the TABLE's Vcmax25/Jmax25/TPU25/Rd25; a cohort's plastic capacities      !
   !      override those per leaf, and the two ratios are what scale Jmax25 and TPU25 with the      !
   !      overriding Vcmax25. The run-level solver selectors ride along so the batch kernel needs   !
   !      nothing else from the configuration. -----------------------------------------------------!
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
      !----- Conductance: wood_kmax * sapwood area / (height * vessel_curl),                        !
      !      docs/science/plant_hydraulics.md eq. (10). --------------------------------------------!
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

   !----- PER-PFT hydraulic parameter table, the hydraulic twin of leaf_photo_table_t. Built ONCE !
   !      per run (each entry carries its own precomputed Kirchhoff lookup, which is what makes   !
   !      per-PFT vulnerability shapes affordable on the hot path) and indexed by the cohort's    !
   !      PFT where the parameters are consumed. -------------------------------------------------!
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
      real(wp) :: ea                    = 46390.0_wp   !< [J/mol]   peaked-Arrhenius activation energy (leaf ea_rd)
      real(wp) :: hd                    = 200000.0_wp  !< [J/mol]   deactivation energy   (leaf hd_rd)
      real(wp) :: ds                    = 490.0_wp     !< [J/mol/K] entropy term          (leaf ds_rd)
   end type wood_params_t

   !----- Fine-root maintenance respiration (ED2 per-broot form; the temperature response is summed !
   !      over the root layers, root_zone_temp_scale). ---------------------------------------------!
   type :: root_params_t
      real(wp) :: ea                 = 46390.0_wp  !< peaked-Arrhenius terms (default leaf ea_rd/hd_rd/ds_rd)
      real(wp) :: hd                 = 200000.0_wp
      real(wp) :: ds                 = 490.0_wp
   end type root_params_t


   !=======================================================================================!
   !     CARBON ALLOCATION -- see meds_plant_carbon_allocation. The allocation kernel is   !
   !     ELEMENTAL over plain cohort scalars, so no derived types live here: the           !
   !     vegetation-dynamics driver hands the kernel raw per-cohort arrays (GPP,           !
   !     maintenance resp, storage, shed amounts, capped demands) and receives the         !
   !     per-pool NPP + growth respiration. Leaf loss and the flush cap come from          !
   !     meds_phenology::leaf_turnover_step.                                               !
   !=======================================================================================!

   !----- Per-PFT vegetation thermal parameters (docs/science/vegetation_energy_dynamics.md eqs.  !
   !      (1) and (8)). --------------------------------------------------------------------------!
   type :: veg_thermal_params_t
      real(wp) :: effarea_heat   = 2.0_wp                   !< [-] sensible sidedness (both leaf sides)
      real(wp) :: effarea_evap   = 1.0_wp                   !< [-] film-evaporation sidedness
      real(wp) :: effarea_transp = 1.0_wp                   !< [-] transpiration sidedness (per PFT)
      real(wp) :: veg_hcap_min   = TISSUE_HCAP_MIN          !< [J/m2 tissue/K] least dry heat capacity
      real(wp) :: c_leaf     = TISSUE_C_LEAF                !< [J/kg/K] dry leaf specific heat
      real(wp) :: c_wood_dry = TISSUE_C_WOOD_DRY            !< [J/kg/K] dry wood specific heat (ED2's, at 15 C)
   end type veg_thermal_params_t

end module meds_plant_types
