! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_canopy_water -- the canopy's interception film: one cohort's capacity-limited       !
! bucket with a Beer interception fraction and the wetted fraction it returns              !
! (docs/science/soil_biophysics.md, "Canopy interception"). The film's evaporation is the  !
! tissue energy balance's wet pathway (meds_plant_energy).                                 !
!==========================================================================================!
module meds_canopy_water
   use meds_kinds,     only : wp
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: intercept_canopy_layer

contains

   !---------------------------------------------------------------------------------------!
   ! One cohort's canopy interception (docs/science/soil_biophysics.md, "Canopy            !
   ! interception"): a capacity-limited film bucket [kg/m2 ground] with a Beer interception!
   ! fraction, returning the wetted fraction sigma_w. The caller sweeps height-sorted      !
   ! cohorts top to bottom, each throughfall the next rain_above, so it is not elemental.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine intercept_canopy_layer(leaf_water, rain_above, lai, sai, e_canopy, dt,     &
                                          leaf_max, wood_max, k_int, alpha_pi, throughfall, drip, sigma_w)
      real(wp), intent(inout) :: leaf_water
      real(wp), intent(in)    :: rain_above, lai, sai, e_canopy, dt, k_int, alpha_pi
      real(wp), intent(in)    :: leaf_max, wood_max   !< [kg/m2 leaf], [kg/m2 wood] film capacities
      real(wp), intent(out)   :: throughfall, drip, sigma_w
      real(wp) :: pai, f_pi, w_max, q_grab, room, q_intr, w_new, over
      pai    = lai + sai
      f_pi   = alpha_pi * (1.0_wp - exp(-k_int * pai))
      w_max  = leaf_max * lai + wood_max * sai
      q_grab = f_pi * rain_above
      room   = max(0.0_wp, w_max - leaf_water) / dt
      q_intr = min(q_grab, room + e_canopy)                      ! bounded by capacity + evap headroom
      !----- A film above capacity (leaf area shed under a full film) drips its excess to the     !
      !      layer below, so the column ledger keeps the water. -------------------------------------!
      w_new = leaf_water + (q_intr - e_canopy) * dt
      over  = max(0.0_wp, w_new - w_max)
      leaf_water  = min(max(w_new, 0.0_wp), w_max)
      drip        = max(0.0_wp, q_grab - q_intr) + over / dt
      throughfall = (rain_above - q_grab) + drip                 ! gap-throughfall + drip -> below
      if (w_max > tiny_num) then
         sigma_w = min(1.0_wp, (leaf_water / w_max) ** (2.0_wp / 3.0_wp))
      else
         sigma_w = 0.0_wp
      end if
   end subroutine intercept_canopy_layer

end module meds_canopy_water
