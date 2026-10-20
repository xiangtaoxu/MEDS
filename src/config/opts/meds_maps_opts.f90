! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_maps_opts -- the [maps] block, as a config leaf: the static maps a run reads its    !
! soil and its height allometry's site index from (docs/configuration.md, "Input maps").   !
! Each map is on the forcing grid; meds_surface_maps reads it at the run's cells. A region !
! needs the three soil maps; a site may use them or set [soil] by hand, and a key a map    !
! provides is refused while the map is set.                                                !
!==========================================================================================!
module meds_maps_opts
   implicit none
   private

   public :: maps_opts_t, soil_maps_on, MAPS_PATH_LEN

   integer, parameter :: MAPS_PATH_LEN = 512

   type :: maps_opts_t
      !----- The three soil maps (docs/configuration.md, "Input maps"). -------------------!
      character(len=MAPS_PATH_LEN) :: soil_depth_map     = ''
      character(len=MAPS_PATH_LEN) :: soil_texture_map   = ''
      character(len=MAPS_PATH_LEN) :: soil_chemistry_map = ''
      !----- Chave et al. (2014)'s stress index E, for height_allometry = "chave2014". ----!
      character(len=MAPS_PATH_LEN) :: chave_e_map        = ''
   end type maps_opts_t

contains

   !----- Whether the soil comes from the maps (the loader sets all three or none). -------!
   pure logical function soil_maps_on(s)
      type(maps_opts_t), intent(in) :: s
      soil_maps_on = len_trim(s%soil_depth_map) > 0
   end function soil_maps_on

end module meds_maps_opts
