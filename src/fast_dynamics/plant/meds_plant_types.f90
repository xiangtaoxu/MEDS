! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_types -- the data types of the sub-daily leaf kernel (no methods, no state):  !
! leaf_env_t, leaf_flux_t and the limitation flags (docs/science/leaf_gas_exchange.md).    !
! The per-PFT tables the plant kernels read (leaf_photo_table_t, hydraulics_table_t) are   !
! the PFT table's (meds_pft_params); phenology's types live in meds_phenology_types.       !
!==========================================================================================!
module meds_plant_types
   use meds_kinds,          only : wp, ik
   implicit none
   private

   public :: leaf_env_t, leaf_flux_t
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP

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
      real(wp)    :: transpiration!< [mol H2O/m2 leaf/s]  E = g VPD / pressure: g = gs, or gs and gb in series
      real(wp)    :: rd           !< [umol/m2 leaf/s]     leaf respiration used
      integer(ik) :: limitation   !< LIM_* of the binding term
      logical     :: converged    !< .true. if the Ci solve met tolerance
      !----- The two water-stress limbs the kernel applied, reported so an inert control shows in the  !
      !      output (env%psi left at 0 makes beta_stomata identically 1). ---------------------------!
      real(wp)    :: beta_stomata    = 1.0_wp  !< [-] stomatal limb min(1, exp(sref*(psi-psi_onset))) (1 = unstressed)
      real(wp)    :: beta_nonstomata = 1.0_wp  !< [-] fast limb, the psi_leaf ramp on Vcmax/Jmax/TPU
   end type leaf_flux_t

end module meds_plant_types
