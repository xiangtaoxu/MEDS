! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_region_opts -- the run mode and the [region] block, as a low-level config leaf           !
! (MEDS_POLYGON_RUNTIME_PLAN.md §9).                                                             !
!                                                                                          !
! A site run ([run].mode = "site", the default) simulates one location from the [site] block. A  !
! region run ([run].mode = "region") simulates every selected forcing cell of a contiguous box    !
! as its own polygon, each at its cell's centre, and writes region files with a `polygon`        !
! dimension. Scattered site networks are separate site runs, not a region (§8.1).                 !
!==========================================================================================!
module meds_region_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION, MAX_DETAIL_POLYGONS

   !----- [run].mode. ----------------------------------------------------------------------!
   integer(ik), parameter :: RUN_MODE_SITE   = 1_ik   !< one location, the [site] block
   integer(ik), parameter :: RUN_MODE_REGION = 2_ik   !< every selected cell of [region].box_nwse

   !----- The most polygons that may also write full single-site output. Each one writes its own  !
   !      file set, cohort and fast tiers included, so the list is meant to stay short. ----------!
   integer(ik), parameter :: MAX_DETAIL_POLYGONS = 64_ik

   type :: region_opts_t
      !----- [N, W, S, E] in degrees. The box may cross 0 or 180 degrees (W > E then). ---------!
      real(wp)    :: box_nwse(4) = 0.0_wp
      !----- Cells whose static land fraction is below this are not simulated: the forcing's     !
      !      valid mask includes large lakes (the whole Caspian) and fractional coastal cells.  ----!
      real(wp)    :: land_fraction_min = 0.5_wp
      !----- Polygon ids (the cell's row-major index on the forcing grid) that also write ordinary !
      !      single-site files, with cohort, patch and fast-tier output. ----------------------------!
      integer(ik) :: n_detail = 0_ik
      integer(ik) :: detail_polygons(MAX_DETAIL_POLYGONS) = -1_ik
   end type region_opts_t

end module meds_region_opts
