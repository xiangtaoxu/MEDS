! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_litter_partition -- the necromass -> litter-destination split.                        !
!                                                                                          !
! This is a biogeochemical LAW (it partitions dead plant carbon into the CENTURY litter pools) !
! that used to live in the column-STATE module for one reason only: src/core could not link   !
! biogeochemistry, and the state module was the one place both sides could see. It is filed    !
! where it belongs now; the DAG carries the edge instead of the law being displaced.           !
!                                                                                          !
! necromass_to_litter takes plain scalars, so it stays callable from the demographic engine's  !
! termination / disturbance kills (which add the six outputs straight onto a soil_carbon_t).   !
! add_necromass is the driver's form: it adds the same split to a litter_input_t and records    !
! the carbon by the tissue it came from, so the litter diagnostics and the soil input are one   !
! accumulation (#399).                                                                          !
!==========================================================================================!
module meds_litter_partition
   use meds_kinds,              only : wp
   use meds_column_state_types, only : litter_input_t
   implicit none
   private

   public :: necromass_to_litter, add_necromass

contains

   !=======================================================================================!
   !  NECROMASS -> LITTER-DESTINATION split (MEDS_SLOW_DYNAMICS_DESIGN.md Part II, B1). Every input !
   !  is a per-area carbon AMOUNT [kgC/m2], already converted from the per-plant cohort pools via    !
   !  nplant; every PFT input is a plain scalar trait. Each tissue goes where it is (#399):          !
   !    - canopy (leaves, and the flowers, fruit and seed of reproduction): all ABOVE ground,        !
   !      labile/structural by f_labile_leaf;                                                         !
   !    - fine roots: all BELOW ground, by f_labile_leaf (no separate root trait);                    !
   !    - wood (stem + coarse roots) and storage (non-structural carbon, held in the sapwood and the   !
   !      coarse roots): above/below by aboveground_frac, the woody split; wood labile/structural by   !
   !      f_labile_stem, storage by f_labile_leaf.                                                    !
   !  Only the STRUCTURAL streams carry lignin (struct_lignin_frac of each). Pass 0 for a channel     !
   !  that does not apply (the turnover call has no wood or storage).                                 !
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
   !  Add one group of dead plant tissue [kgC/m2] to a patch's litter accumulator: the pool    !
   !  streams the soil reads (necromass_to_litter, reproduction going with the canopy), and the !
   !  same carbon by the tissue it came from, which the litter diagnostics report.             !
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
