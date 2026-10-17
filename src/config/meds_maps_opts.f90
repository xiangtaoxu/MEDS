! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_maps_opts -- the [maps] block, as a low-level config leaf: the static maps a run  !
! reads its soil and its height allometry's site index from (MEDS_POLYGON_RUNTIME_PLAN.md §14, !
! RP3). Each map is on the forcing grid; meds_surface_maps reads it at the run's cells.        !
!                                                                                          !
! A region run needs the three soil maps: every polygon takes its own soil column from its    !
! cell. A site run may use them too, at its forcing cell, or set [soil_column] by hand -- not  !
! both. A key a map provides is refused when the map is set, so a run never holds two values   !
! for one quantity.                                                                            !
!==========================================================================================!
module meds_maps_opts
   implicit none
   private

   public :: maps_opts_t, soil_maps_on, MAPS_PATH_LEN

   integer, parameter :: MAPS_PATH_LEN = 512

   type :: maps_opts_t
      !----- The soil: the column's depth and its bottom boundary (soil_depth, bc_class), its   !
      !      deep temperature (deep_temp); each SoilGrids layer's texture (sand, silt, clay,     !
      !      bulk_density, coarse_fragments) and chemistry (organic_carbon, cec, ph_h2o). -------!
      character(len=MAPS_PATH_LEN) :: soil_depth_map     = ''
      character(len=MAPS_PATH_LEN) :: soil_texture_map   = ''
      character(len=MAPS_PATH_LEN) :: soil_chemistry_map = ''
      !----- Chave et al. (2014)'s site stress index E (chave_e), for height_allometry = "chave2014". !
      character(len=MAPS_PATH_LEN) :: chave_e_map        = ''
   end type maps_opts_t

contains

   !----- Whether the run takes its soil from the maps (the loader sets all three or none). ----!
   pure logical function soil_maps_on(s)
      type(maps_opts_t), intent(in) :: s
      soil_maps_on = len_trim(s%soil_depth_map) > 0
   end function soil_maps_on

end module meds_maps_opts
