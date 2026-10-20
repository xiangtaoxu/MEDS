! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_leaf_opts -- the leaf-model selector codes, as a low-level config leaf.             !
!                                                                                          !
! The sibling of meds_biophysics_opts / meds_biogeochem_opts, for the same reason: these   !
! codes are read by both the meds_config aggregator (which loads and validates them from   !
! TOML) and the leaf kernel (which switches on them). Owning them here keeps the kernel    !
! meds_leaf_gas_exchange free of meds_config.                                              !
!                                                                                          !
! They do not sit in meds_plant_types, beside the leaf_photo_table_t fields that hold      !
! them: meds_plant_types imports the PFT trait table, so it is above the config layer, and !
! meds_config importing it would be a cycle.                                               !
!==========================================================================================!
module meds_leaf_opts
   use meds_kinds, only : ik
   implicit none
   private

   public :: SM_LEUNING, SM_MEDLYN, SM_KATUL
   public :: COLIM_MIN, COLIM_QUADRATIC

   !----- Stomatal-conductance model ([leaf_physiology].stomatal_model): docs/science/leaf_gas_exchange.md §3. ---!
   integer(ik), parameter :: SM_LEUNING = 1_ik      !< Leuning (1995) BWB-VPD semi-empirical
   integer(ik), parameter :: SM_MEDLYN  = 2_ik      !< Medlyn et al. (2011) unified optimization (USO)
   integer(ik), parameter :: SM_KATUL   = 3_ik      !< Katul et al. (2010) analytical optimization

   !----- Co-limitation of the FvCB / C4 rates ([leaf_physiology].colimitation): docs/science/leaf_gas_exchange.md §2. -!
   integer(ik), parameter :: COLIM_MIN       = 1_ik !< sharp minimum
   integer(ik), parameter :: COLIM_QUADRATIC = 2_ik !< smoothed co-limitation quadratics

   !----- tresp_arrhenius / tresp_peaked are the same KIND of selector but are owned by      !
   !      meds_therm_lib, which defines the responses they select. Import them there.         !
end module meds_leaf_opts
