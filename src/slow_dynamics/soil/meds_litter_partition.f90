! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_litter_partition -- dead plant carbon to litter (docs/science/soil_carbon.md §5).   !
!   * necromass_to_litter -- the split, on plain scalars, so the demography operators'     !
!       kills can also call it and add its six outputs straight onto a soil_carbon_t.      !
!   * add_necromass       -- the slow driver's form: adds the split to a litter_input_t and !
!       records the carbon by tissue, so the litter diagnostics and soil input are one sum. !
!==========================================================================================!
module meds_litter_partition
   use meds_kinds,              only : wp
   use meds_patch_state_vars,   only : litter_input_t
   implicit none
   private

   public :: necromass_to_litter, add_necromass

contains

   !=======================================================================================!
   !  Split per-area dead carbon [kgC/m2] by tissue: canopy (leaves and reproduction) all  !
   !  above ground and fine roots all below, both labile by f_labile_leaf; wood and storage !
   !  above/below by aboveground_frac, wood labile by f_labile_stem and storage by         !
   !  f_labile_leaf. Only the structural streams carry lignin. Pass 0 for a tissue that    !
   !  does not apply (the turnover call has no wood or storage).                           !
   !=======================================================================================!
   pure subroutine necromass_to_litter(canopy_c, root_c, wood_c, storage_c,                          &
                                       f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac,   &
                                       labile_grnd, labile_soil, struct_grnd, struct_soil,             &
                                       lignin_grnd, lignin_soil)
      real(wp), intent(in)  :: canopy_c, root_c, wood_c, storage_c
      real(wp), intent(in)  :: f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac
      real(wp), intent(out) :: labile_grnd, labile_soil, struct_grnd, struct_soil
      real(wp), intent(out) :: lignin_grnd, lignin_soil
      real(wp) :: canopy_lab, canopy_str, root_lab, root_str, woody_lab, woody_str

      canopy_lab = canopy_c * f_labile_leaf
      canopy_str = canopy_c * (1.0_wp - f_labile_leaf)
      root_lab   = root_c * f_labile_leaf
      root_str   = root_c * (1.0_wp - f_labile_leaf)
      woody_lab  = wood_c * f_labile_stem + storage_c * f_labile_leaf
      woody_str  = wood_c * (1.0_wp - f_labile_stem) + storage_c * (1.0_wp - f_labile_leaf)

      labile_grnd = canopy_lab + woody_lab * aboveground_frac
      labile_soil = root_lab   + woody_lab * (1.0_wp - aboveground_frac)
      struct_grnd = canopy_str + woody_str * aboveground_frac
      struct_soil = root_str   + woody_str * (1.0_wp - aboveground_frac)
      lignin_grnd = lignin_frac * struct_grnd
      lignin_soil = lignin_frac * struct_soil
   end subroutine necromass_to_litter

   !=======================================================================================!
   !  Add one group of dead tissue [kgC/m2] to a patch's litter accumulator: the pool      !
   !  streams (reproduction going with the canopy) and the same carbon by tissue.          !
   !=======================================================================================!
   pure subroutine add_necromass(lit, leaf_c, root_c, wood_c, storage_c, repro_c,                   &
                                 f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac)
      type(litter_input_t), intent(inout) :: lit
      real(wp),             intent(in)    :: leaf_c, root_c, wood_c, storage_c, repro_c
      real(wp),             intent(in)    :: f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac
      real(wp) :: lab_g, lab_s, str_g, str_s, lig_g, lig_s

      call necromass_to_litter(leaf_c + repro_c, root_c, wood_c, storage_c,                         &
                               f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac,         &
                               lab_g, lab_s, str_g, str_s, lig_g, lig_s)
      lit%labile_grnd = lit%labile_grnd + lab_g
      lit%labile_soil = lit%labile_soil + lab_s
      lit%struct_grnd = lit%struct_grnd + str_g
      lit%struct_soil = lit%struct_soil + str_s
      lit%lignin_grnd = lit%lignin_grnd + lig_g
      lit%lignin_soil = lit%lignin_soil + lig_s
      lit%leaf_c      = lit%leaf_c      + leaf_c
      lit%fineroot_c  = lit%fineroot_c  + root_c
      lit%wood_c      = lit%wood_c      + wood_c
      lit%storage_c   = lit%storage_c   + storage_c
      lit%repro_c     = lit%repro_c     + repro_c
   end subroutine add_necromass

end module meds_litter_partition
