! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_hydraulics -- the plant-water facts the rest of the fast loop needs about     !
! each cohort (docs/science/plant_hydraulics.md): the hydrostatic head from its wood to    !
! its crown, its leaves' water potential, and its roots -- how deep they reach and how its !
! fine roots spread over the soil layers -- and the limit that keeps the soil-plant water  !
! solver's step solvable when a leaf runs dry (supply_limited_transpiration).              !
!                                                                                          !
! How the plant's water MOVES -- uptake and sapflow -- and how much a tissue holds are in  !
! meds_plant_hydr_lib, and meds_spac_solver solves them together with the soil's           !
! (docs/science/soil_plant_water.md).                                                      !
!==========================================================================================!
module meds_plant_hydraulics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : grav_head, tiny_num
   use meds_plant_types, only : hydro_params_t, hydro_opts_t
   use meds_plant_hydr_lib, only : root_fraction_profile, tissue_water_to_psi
   use meds_numerics,       only : smooth_max0
   implicit none
   private

   public :: leaf_node_psi, crown_gravity_head
   public :: cohort_root_depth, cohort_root_profile
   public :: supply_limited_transpiration

   !----- [-] Width, in relative water content, over which smooth_max0 rounds the leaves'  !
   !      supply limit (docs/science/soil_plant_water.md eq. (3)). ------------------------!
   real(wp), parameter :: TRANSP_CLOSE_WIDTH = 2.0e-2_wp

contains

   !----- The water-potential drop from the wood node to the crown [MPa]: the hydrostatic head of  !
   !      the cohort's height, or zero with gravity switched off. ----------------------------------!
   pure real(wp) function crown_gravity_head(height, o) result(grav)
      real(wp),           intent(in) :: height   !< [m] cohort height
      type(hydro_opts_t), intent(in) :: o
      if (o%gravity_on) then ; grav = grav_head*height ; else ; grav = 0.0_wp ; end if
   end function crown_gravity_head

   !----- The leaf node's water potential [MPa], read from the leaf water and its pressure-volume    !
   !      curve (docs/science/plant_hydraulics.md eqs. (4)-(5)). A cohort with no leaves -- a        !
   !      deciduous canopy shed to bare -- has an empty leaf store and no curve to read: the reading !
   !      would be meaningless (the relative water content floors at 1e-4, about -15,000 MPa). Its   !
   !      leaf node instead sits at the wood potential less the gravity head to the crown, the value !
   !      the leaf node relaxes to with no transpiration and the one Xu et al. (2016) give a         !
   !      leafless cohort. So a bare cohort's predawn potential, which the stomata and the phenology !
   !      water cue read, follows its stem and roots. -----------------------------------------------!
   pure real(wp) function leaf_node_psi(leaf_water, bleaf, wood_water, bwood, height, p, o) result(psi_l)
      real(wp),             intent(in) :: leaf_water  !< [kg H2O/plant] leaf water
      real(wp),             intent(in) :: bleaf       !< [kgC/plant] leaf carbon
      real(wp),             intent(in) :: wood_water  !< [kg H2O/plant] wood (sapwood + fine-root) water
      real(wp),             intent(in) :: bwood       !< [kgC/plant] sapwood + fine-root carbon
      real(wp),             intent(in) :: height      !< [m] cohort height
      type(hydro_params_t), intent(in) :: p
      type(hydro_opts_t),   intent(in) :: o
      if (p%leaf_curve%water_sat*bleaf > tiny_num) then
         psi_l = tissue_water_to_psi(leaf_water, p%leaf_curve, bleaf)
      else
         psi_l = tissue_water_to_psi(wood_water, p%wood_curve, bwood) - crown_gravity_head(height, o)
      end if
   end function leaf_node_psi

   !----- A cohort's rooting depth [m], docs/science/plant_hydraulics.md eq. (12): root_depth_b1 *     !
   !      height^root_depth_b2 (ED2's height allometry, IALLOM 1: 5 m for a 35 m tree), capped at the  !
   !      soil column. --------------------------------------------------------------------------------!
   pure real(wp) function cohort_root_depth(height, hp, soil_depth) result(depth)
      real(wp),             intent(in) :: height      !< [m] cohort height
      type(hydro_params_t), intent(in) :: hp          !< the cohort's PFT
      real(wp),             intent(in) :: soil_depth  !< [m] depth of the soil column's bottom
      depth = min(hp%root_depth_b1 * max(height, 0.0_wp)**hp%root_depth_b2, soil_depth)
   end function cohort_root_depth

   !----- A cohort's fine roots in each soil layer [-], summing to one                               !
   !      (docs/science/plant_hydraulics.md eq. (13)): ED2's exponential profile within its own      !
   !      rooting depth D, (1 - beta^(d/D)) / (1 - beta) above depth d. Layers below D hold none. ---!
   pure subroutine cohort_root_profile(height, hp, z_top, z_bot, frac)
      real(wp),             intent(in)  :: height      !< [m] cohort height
      type(hydro_params_t), intent(in)  :: hp          !< the cohort's PFT
      real(wp),             intent(in)  :: z_top(:)    !< [m] depth of each layer's top (>= 0)
      real(wp),             intent(in)  :: z_bot(:)    !< [m] depth of each layer's bottom (> z_top)
      real(wp),             intent(out) :: frac(:)     !< [-] the cohort's fine roots in each layer
      real(wp) :: depth
      depth = cohort_root_depth(height, hp, z_bot(size(z_bot)))
      frac = root_fraction_profile(hp%root_beta, depth, z_top, z_bot) / (1.0_wp - hp%root_beta)
   end subroutine cohort_root_profile

   !----- The soil-plant water solver's limit on a drying leaf's transpiration [kg/m2/s]. !
   !      The stomata set the demand e_demand at the step's start and the solver holds it  !
   !      through the step; while the leaf's symplastic relative water content r stays     !
   !      above r_close, the content at the closing potential, the leaf transpires the     !
   !      demand, and below it the demand scaled by r/r_close, so the implicit step always !
   !      has a solution and the leaf never empties. The corner is rounded over            !
   !      TRANSP_CLOSE_WIDTH (docs/science/soil_plant_water.md eq. (3)). ------------------!
   elemental subroutine supply_limited_transpiration(e_demand, r, r_close, e, de_dr)
      real(wp), intent(in)  :: e_demand, r, r_close
      real(wp), intent(out) :: e, de_dr
      real(wp) :: excess, dexcess, share
      call smooth_max0(r - r_close, TRANSP_CLOSE_WIDTH, excess, dexcess)
      share = (r - excess) / r_close
      e     = e_demand * share
      de_dr = e_demand * (1.0_wp - dexcess) / r_close
   end subroutine supply_limited_transpiration

end module meds_plant_hydraulics
