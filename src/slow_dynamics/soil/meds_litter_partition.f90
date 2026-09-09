!==========================================================================================!
! meds_litter_partition -- the necromass -> litter-destination split.                        !
!                                                                                          !
! This is a biogeochemical LAW (it partitions dead plant carbon into the CENTURY litter pools) !
! that used to live in the column-STATE module for one reason only: src/core could not link   !
! biogeochemistry, and the state module was the one place both sides could see. It is filed    !
! where it belongs now; the DAG carries the edge instead of the law being displaced.           !
!                                                                                          !
! Every argument is a plain scalar, so it has no derived-type dependency and stays callable    !
! from the driver (which routes the result through litter_input_t) and from the demographic    !
! engine's termination / disturbance kills (which add the six outputs straight onto a           !
! soil_carbon_t).                                                                               !
!==========================================================================================!
module meds_litter_partition
   use meds_kinds, only : wp
   implicit none
   private

   public :: necromass_to_litter

contains

   !=======================================================================================!
   !  NECROMASS -> LITTER-DESTINATION split (MEDS_SLOW_DYNAMICS_DESIGN.md Part II, B1). Every input !
   !  is a plain per-area carbon AMOUNT [kgC/m2] (leaf/storage necromass, fine-root necromass, wood  !
   !  necromass) already converted from the per-plant cohort pools via nplant; every PFT input is a   !
   !  plain scalar trait, so this has NO derived-type dependency (callable from BOTH the driver,      !
   !  which routes the result through litter_input_t/build_litter_input, and the core engine's         !
   !  termination/disturbance kills, which cannot link biogeochemistry -- it adds these six outputs    !
   !  straight onto a soil_carbon_t's carbon+lignin fields). Leaf necromass (+ storage, bundled since  !
   !  both are canopy-associated labile-eligible pools) splits labile/structural by f_labile_leaf,      !
   !  then each of those splits above/below by aboveground_frac; fine-root necromass reuses            !
   !  f_labile_leaf (no separate root trait) but lands ENTIRELY below-ground (no agf split -- fine      !
   !  roots are definitionally below-ground); wood/CWD necromass splits labile/structural by            !
   !  f_labile_stem, then above/below by aboveground_frac. Only the STRUCTURAL streams carry lignin     !
   !  (struct_lignin_frac of each). Pass 0 for any necromass channel that does not apply (e.g. the       !
   !  turnover/shed call site has no wood or storage component).                                         !
   !=======================================================================================!
   pure subroutine necromass_to_litter(leaf_c, root_c, wood_c, storage_c,                            &
                                       f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac,   &
                                       labile_grnd, labile_soil, struct_grnd, struct_soil,             &
                                       lignin_grnd, lignin_soil)
      real(wp), intent(in)  :: leaf_c, root_c, wood_c, storage_c
      real(wp), intent(in)  :: f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac
      real(wp), intent(out) :: labile_grnd, labile_soil, struct_grnd, struct_soil
      real(wp), intent(out) :: lignin_grnd, lignin_soil
      real(wp) :: canopy_c, canopy_lab, canopy_str, root_lab, root_str, wood_lab, wood_str

      canopy_c   = leaf_c + storage_c
      canopy_lab = canopy_c * f_labile_leaf
      canopy_str = canopy_c * (1.0_wp - f_labile_leaf)
      root_lab   = root_c * f_labile_leaf
      root_str   = root_c * (1.0_wp - f_labile_leaf)
      wood_lab   = wood_c * f_labile_stem
      wood_str   = wood_c * (1.0_wp - f_labile_stem)

      labile_grnd = (canopy_lab + wood_lab) * aboveground_frac
      labile_soil = (canopy_lab + wood_lab) * (1.0_wp - aboveground_frac) + root_lab
      struct_grnd = (canopy_str + wood_str) * aboveground_frac
      struct_soil = (canopy_str + wood_str) * (1.0_wp - aboveground_frac) + root_str
      lignin_grnd = lignin_frac * struct_grnd
      lignin_soil = lignin_frac * struct_soil
   end subroutine necromass_to_litter

end module meds_litter_partition
