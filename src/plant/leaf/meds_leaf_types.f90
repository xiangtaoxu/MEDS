!==========================================================================================!
! meds_leaf_types -- the data structures of the leaf-physiology interface seam.            !
!                                                                                          !
! Three public types describe the leaf-level gas-exchange calculation, plus the limitation- !
! regime and photosynthetic-pathway flags:                                                  !
!   * leaf_env_t          -- the EXPLICIT environmental drivers (PAR, leaf temperature, VPD, !
!                            reference CO2, pressure, leaf water potential, boundary layer).  !
!   * leaf_flux_t         -- the computed fluxes (net/gross assimilation, stomatal           !
!                            conductance, intercellular CO2, transpiration, the binding       !
!                            limitation, convergence).                                        !
!   * leaf_photo_params_t -- a flat, self-contained parameter set for ONE PFT (per-PFT traits !
!                            plus the shared biochemistry constants), so the solver needs no   !
!                            reference back into meds_config.                                  !
!                                                                                          !
! These are pure DATA: no methods, no hidden state. The sealed driver meds_leaf_physiology%   !
! leaf_gas_exchange consumes a leaf_env_t and (via the config) a PFT index and returns a       !
! leaf_flux_t -- the leaf-level analogue of the demography data-rate seam.                     !
!==========================================================================================!
module meds_leaf_types
   use meds_kinds,      only : wp, ik
   use meds_pft_params, only : PATH_C3, PATH_C4
   implicit none
   private

   public :: leaf_env_t, leaf_flux_t, leaf_photo_params_t
   public :: PATH_C3, PATH_C4                       ! re-export (pathway lives with the PFT traits)
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP

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

end module meds_leaf_types
