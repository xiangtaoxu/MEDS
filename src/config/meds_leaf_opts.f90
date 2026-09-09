!==========================================================================================!
! meds_leaf_opts -- the LEAF-model selector codes, as a low-level config leaf.               !
!                                                                                          !
! The sibling of meds_biophysics_opts / meds_biogeochem_opts, and it exists for the same     !
! reason: these codes are read by BOTH the meds_config aggregator (which loads and validates  !
! them from TOML) and the leaf kernel (which switches on them). Owning them here is what lets  !
! the kernel stay config-free -- `meds_leaf_gas_exchange` imported `meds_config` for these five !
! integers alone, which is the last of the structure plan's decision #10.                       !
!                                                                                          !
! NOTE the plan filed them in `meds_plant_types`, next to the leaf_photo_table_t fields that    !
! hold them. That is not reachable: meds_config must see them too, and meds_plant_types sits     !
! ABOVE the config layer (it reads the PFT trait table), so the import would be a cycle. A       !
! config leaf is the same intent -- "beside the fields, not in the aggregator" -- one layer down. !
!==========================================================================================!
module meds_leaf_opts
   use meds_kinds, only : ik
   implicit none
   private

   public :: SM_LEUNING, SM_MEDLYN, SM_KATUL
   public :: COLIM_MIN, COLIM_QUADRATIC

   !----- Stomatal-conductance model ([leaf_physiology].stomatal_model). ------------------!
   integer(ik), parameter :: SM_LEUNING = 1_ik      !< Leuning (1995) BWB-VPD semi-empirical
   integer(ik), parameter :: SM_MEDLYN  = 2_ik      !< Medlyn et al. (2011) unified optimization (USO)
   integer(ik), parameter :: SM_KATUL   = 3_ik      !< Katul et al. (2010) analytical optimization

   !----- Co-limitation form combining the FvCB / C4 limitation rates. -------------------!
   integer(ik), parameter :: COLIM_MIN       = 1_ik !< sharp minimum
   integer(ik), parameter :: COLIM_QUADRATIC = 2_ik !< smoothed co-limitation quadratics

   !----- TRESP_ARRHENIUS / TRESP_PEAKED are the same KIND of selector but are owned by      !
   !      meds_temp_response, which defines the responses they select. Import them there.     !
end module meds_leaf_opts
