! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_types -- the data types of the sub-daily plant kernels (no methods, no state):!
!   * leaf        -- leaf_env_t / leaf_flux_t / leaf_photo_args_t / leaf_photo_table_t and !
!                    the limitation flags (docs/science/leaf_gas_exchange.md);             !
!   * hydraulics  -- hydro_args_t / hydro_args_table_t (per-PFT traits).                   !
! Phenology's types live in meds_phenology_types.                                          !
!==========================================================================================!
module meds_plant_types
   use meds_kinds,          only : wp, ik
   use meds_plant_hydr_lib, only : hydro_table_t, plant_hydraulics_t
   implicit none
   private

   !----- Leaf -----------------------------------------------------------------------------!
   public :: leaf_env_t, leaf_flux_t, leaf_photo_args_t, leaf_photo_table_t
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP
   !----- Hydraulics -----------------------------------------------------------------------!
   public :: hydro_args_t, hydro_args_table_t

   !=======================================================================================!
   !     Leaf -- the leaf-level gas-exchange interface.                                     !
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
      !----- The supply-side potential driving the slow stomatal control (default 0, well-watered): the !
      !      fast loop passes the cohort's predawn (previous-day maximum) leaf potential, dmax_psi_leaf. !
      !      Named neutrally: it equals the soil's only in wet soil (leaf_gas_exchange.md §4.2). --------!
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
      !----- The two water-stress limbs the kernel applied, reported so an inert control shows in the  !
      !      output (env%psi left at 0 makes beta_stomata identically 1). ---------------------------!
      real(wp)    :: beta_stomata    = 1.0_wp  !< [-] stomatal limb min(1, exp(sref*(psi-psi_onset))) (1 = unstressed)
      real(wp)    :: beta_nonstomata = 1.0_wp  !< [-] fast limb, the psi_leaf ramp on Vcmax/Jmax/TPU
   end type leaf_flux_t

   !----- Flat per-PFT parameter set (per-PFT traits + shared biochemistry), self-contained !
   !       so the solver never references meds_config. Filled by the interface from the config.!
   type :: leaf_photo_args_t
      integer(ik) :: pathway        !< PATH_C3 | PATH_C4 (meds_pft_params)
      !----- Per-PFT capacities at 25 degC and stomatal/water-stress traits. ---------------!
      real(wp) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(wp) :: g0, g1, d0, quantum_yield
      !----- theta_j is the electron-transport curvature; the theta_*_c3/c4 pairs are the co-limitation !
      !      curvatures, one pair per pathway (docs/science/leaf_gas_exchange.md §2). -----------------!
      real(wp) :: theta_j
      real(wp) :: theta_cj_c3, theta_ip_c3
      real(wp) :: theta_cj_c4, theta_ic_c4
      !----- psi_open/psi_close [MPa]: the leaf potential at which beta_nonstomata is 1 and 0, per PFT   !
      !      (pft_wstress_psi_open/_close: the turgor-loss point and twice it unless the PFT file sets them). !
      real(wp) :: lambda25, psi_open, psi_close, lambda_psi_exp, sref_stomata
      !----- Predawn potential [MPa] below which the stomatal stress begins, eq. (10); a run sets it per !
      !      PFT (pft_stomata_psi_onset); the 0 here is for a record built by hand. --------------------!
      real(wp) :: psi_onset = 0.0_wp
      !----- Shared biochemistry constants at 25 degC + activation/deactivation terms. -----!
      real(wp) :: kc25, ko25, gstar25
      real(wp) :: ea_kc, ea_ko, ea_gstar, ea_vcmax, ea_jmax, ea_rd
      real(wp) :: hd_vcmax, hd_jmax, hd_rd, ds_vcmax, ds_jmax, ds_rd
      real(wp) :: o2_mol_frac, absorptance, phi_psii
      real(wp) :: medlyn_vpd_min = 0.05_wp   !< [kPa] the VPD the Medlyn model uses at least ([leaf_physiology])
   end type leaf_photo_args_t

   !----- leaf_photo_table_t -- every PFT's leaf parameters, built once per run by                  !
   !      build_leaf_photo_table (meds_fast_config). A cohort's plastic capacities override the     !
   !      table's, Jmax25 and TPU25 scaling with Vcmax25 by the two ratios; the solver selectors   !
   !      ride along. ------------------------------------------------------------------------------!
   type :: leaf_photo_table_t
      integer(ik) :: n_pft = 0_ik
      type(leaf_photo_args_t),   allocatable :: pft(:)              !< per-PFT parameters (table capacities)
      real(wp),                  allocatable :: jmax_vcmax_ratio(:) !< [--] Jmax25 / Vcmax25 per PFT
      real(wp),                  allocatable :: tpu_vcmax_ratio(:)  !< [--] TPU25 / Vcmax25 per PFT
      integer(ik) :: stomatal_model     = 0_ik    !< SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form = 0_ik    !< TEMP_RESPONSE_ARRHENIUS | TEMP_RESPONSE_PEAKED
      integer(ik) :: colimitation       = 0_ik    !< COLIM_MIN | COLIM_QUADRATIC
      logical     :: use_boundary_layer = .true.  !< couple through the leaf boundary layer (gb)
   end type leaf_photo_table_t

   !=======================================================================================!
   !     Hydraulics -- the per-PFT plant-water traits. How water moves                         !
   !     is the soil-plant water solver's (meds_spac_solver).                                 !
   !=======================================================================================!
   !----- A PFT's hydraulic traits (cfg%pft%hydro) and the Kirchhoff table its wood_kexp sets. !
   type, extends(plant_hydraulics_t) :: hydro_args_t
      !----- Kirchhoff lookup from wood_kexp (build_hydro_table), read only for wood_kexp outside {1,2}, !
      !      where no closed form exists (docs/science/plant_hydraulics.md §3). ----------------------!
      type(hydro_table_t) :: vuln_table
      real(wp) :: rwc_close   !< [-] leaf relative water content at the PFT's stomatal closing potential
   end type hydro_args_t

   !----- Per-PFT hydraulic parameter table, built once per run; each entry carries its own       !
   !      Kirchhoff lookup, so each PFT keeps its vulnerability shape. ---------------------------!
   type :: hydro_args_table_t
      integer(ik) :: n_pft = 0_ik
      type(hydro_args_t),   allocatable :: pft(:)
   end type hydro_args_table_t

end module meds_plant_types
