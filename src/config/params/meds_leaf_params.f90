! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_leaf_params -- the [leaf_physiology] selector codes, as a config leaf: meds_config  !
! reads and validates them, and the leaf kernel (meds_leaf_gas_exchange) switches on them  !
! without importing meds_config.                                                           !
!==========================================================================================!
module meds_leaf_params
   use meds_kinds, only : ik
   implicit none
   private

   public :: SM_LEUNING, SM_MEDLYN, SM_KATUL
   public :: COLIM_MIN, COLIM_QUADRATIC

   !----- [leaf_physiology].stomatal_model (docs/science/leaf_gas_exchange.md §3). --------!
   integer(ik), parameter :: SM_LEUNING = 1_ik      !< Leuning (1995) BWB-VPD semi-empirical
   integer(ik), parameter :: SM_MEDLYN  = 2_ik      !< Medlyn et al. (2011) unified optimization (USO)
   integer(ik), parameter :: SM_KATUL   = 3_ik      !< Katul et al. (2010) analytical optimization

   !----- [leaf_physiology].colimitation (docs/science/leaf_gas_exchange.md §2). ----------!
   integer(ik), parameter :: COLIM_MIN       = 1_ik !< sharp minimum
   integer(ik), parameter :: COLIM_QUADRATIC = 2_ik !< smoothed co-limitation quadratics

   !----- The temp_response_form codes, TEMP_RESPONSE_*, are in meds_therm_lib with the    !
   !      responses they select. ----------------------------------------------------------!
end module meds_leaf_params
