! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_region_opts -- the run mode and the [region] block, as a low-level config leaf           !
! (MEDS_POLYGON_RUNTIME_PLAN.md §9).                                                             !
!                                                                                          !
! A site run ([run].mode = "site", the default) simulates one location from the [site] block. A  !
! region run ([run].mode = "region") simulates every selected forcing cell of a contiguous box    !
! as its own polygon, each at its cell's centre, and writes region files with a `polygon`        !
! dimension, every variable and tier in them. Scattered site networks are separate site runs, not !
! a region (§8.1).                                                                                 !
!==========================================================================================!
module meds_region_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION

   !----- [run].mode. ----------------------------------------------------------------------!
   integer(ik), parameter :: RUN_MODE_SITE   = 1_ik   !< one location, the [site] block
   integer(ik), parameter :: RUN_MODE_REGION = 2_ik   !< every selected cell of [region].box_nwse

   type :: region_opts_t
      !----- [N, W, S, E] in degrees. The box may cross 0 or 180 degrees (W > E then). ---------!
      real(wp)    :: box_nwse(4) = 0.0_wp
      !----- Cells whose static land fraction is below this are not simulated: the forcing's     !
      !      valid mask includes large lakes (the whole Caspian) and fractional coastal cells.  ----!
      real(wp)    :: land_fraction_min = 0.5_wp
      !----- The cell mask (MEDS_POLYGON_RUNTIME_PLAN.md §14 RP2): a map on the forcing grid; a cell of  !
      !      the box runs where `cell_mask_var` is at least `cell_mask_min`, and is skipped where it is  !
      !      below or the map holds no value. '' runs every selected cell of the box. --------------------!
      character(len=512) :: cell_mask     = ''
      character(len=64)  :: cell_mask_var = 'mask'
      real(wp)           :: cell_mask_min = 0.5_wp
   end type region_opts_t

end module meds_region_opts
