! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_canopy_water -- the canopy's interception film (docs/science/soil_biophysics.md,    !
! "Canopy interception"): one cohort's capacity-limited bucket with a Beer interception    !
! fraction (intercept_canopy_layer); the step's capture, swept down the height-sorted      !
! cohorts (canopy_interception); each film's evaporation held to what it holds             !
! (limit_film_evaporation); and the stepped films kept within their capacities             !
! (clamp_canopy_film). The film's evaporation is the tissue energy balance's wet pathway   !
! (meds_plant_energy).                                                                     !
!==========================================================================================!
module meds_canopy_water
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: canopy_water_out_t
   public :: intercept_canopy_layer, canopy_interception, limit_film_evaporation, clamp_canopy_film

   !----- The interception film over one step: the wetted fraction and film-evaporation            !
   !      conductances, the liquid enthalpy the film is valued at, and the interception rates,     !
   !      held for the step. All zero when canopy_water_on is off. --------------------------------!
   type :: canopy_water_out_t
      real(wp), allocatable :: g_film_leaf(:), g_film_w(:)   !< [m/s] film-evaporation conductance, leaf/wood
      real(wp), allocatable :: f_wet_c(:)                    !< [-]   combined wetted fraction (sigma_w)
      !----- The film's liquid enthalpy, internal_energy_liquid(t_film_valuation). The tissue pays     !
      !      enthalpy_vapor(T) minus this per kg it evaporates, so film, tissue and CAS close exactly. -!
      real(wp) :: film_liquid_enthalpy = 0.0_wp              !< [J/kg]
      !----- Interception rates (capture and capacity only): integrating them by explicit Euler over  !
      !      the step reproduces the bucket's one-shot update exactly. ------------------------------!
      real(wp), allocatable :: intercept_leaf(:), intercept_wood(:)   !< [kg/m2 ground/s] (ncoh)
   end type canopy_water_out_t

contains

   !---------------------------------------------------------------------------------------!
   ! One cohort's canopy interception (docs/science/soil_biophysics.md, "Canopy            !
   ! interception"): a capacity-limited film bucket [kg/m2 ground] with a Beer             !
   ! interception fraction, returning the wetted fraction sigma_w. The caller sweeps       !
   ! height-sorted cohorts top to bottom, each throughfall the next rain_above, so it is   !
   ! not elemental.                                                                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine intercept_canopy_layer(leaf_water, rain_above, lai, sai, dt, leaf_max, wood_max, k_int,     &
                                          alpha_pi, throughfall, drip, sigma_w)
      real(wp), intent(inout) :: leaf_water             !< [kg/m2 ground] the film
      real(wp), intent(in)    :: rain_above             !< [kg/m2/s] precipitation reaching the cohort
      real(wp), intent(in)    :: lai, sai               !< [m2/m2] leaf and stem area index
      real(wp), intent(in)    :: dt                     !< [s]
      real(wp), intent(in)    :: leaf_max, wood_max     !< [kg/m2 leaf], [kg/m2 wood] film capacities
      real(wp), intent(in)    :: k_int, alpha_pi        !< [-] the interception fraction's extinction and scale
      real(wp), intent(out)   :: throughfall, drip      !< [kg/m2/s] to the layer below
      real(wp), intent(out)   :: sigma_w                !< [-] the wetted fraction
      real(wp) :: pai, f_pi, w_max, q_grab, room, q_intr, w_new, over
      pai    = lai + sai
      f_pi   = alpha_pi * (1.0_wp - exp(-k_int * pai))
      w_max  = leaf_max * lai + wood_max * sai
      q_grab = f_pi * rain_above
      room   = max(0.0_wp, w_max - leaf_water) / dt
      q_intr = min(q_grab, room)                                 ! bounded by the room left
      !----- A film above capacity (leaf area shed under a full film) drips its excess to the     !
      !      layer below, so the column ledger keeps the water. -------------------------------------!
      w_new = leaf_water + q_intr * dt
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

   !---------------------------------------------------------------------------------------!
   ! canopy_interception -- the step's capture by the canopy's films                       !
   ! (docs/science/soil_biophysics.md, "Canopy interception"): one leaf-and-wood bucket    !
   ! per cohort, tallest first, each cohort's throughfall the next one's rain. Capture     !
   ! only: the films' evaporation is the tissue energy balance's. Each cohort's film is    !
   ! shared between leaf and wood by capacity, so neither is left above its own, and held  !
   ! as the rate whose Euler step over dt is the bucket's update. Returns what reaches the !
   ! ground.                                                                               !
   !---------------------------------------------------------------------------------------!
   pure subroutine canopy_interception(n, rain, leaf_water, wood_water, lai, wai, leaf_max, wood_max, k_int,     &
                                       alpha_pi, dt, intercept_leaf, intercept_wood, f_wet, throughfall)
      integer(ik), intent(in)    :: n
      real(wp),    intent(in)    :: rain                             !< [kg/m2/s] rain and snow above the canopy
      real(wp),    intent(in)    :: leaf_water(:), wood_water(:)     !< [kg/m2 ground] the films at the start
      real(wp),    intent(in)    :: lai(:), wai(:)                   !< [m2/m2] leaf and wood area index
      real(wp),    intent(in)    :: leaf_max(:), wood_max(:)         !< [kg/m2 leaf], [kg/m2 wood] capacities
      real(wp),    intent(in)    :: k_int, alpha_pi                  !< [-] intercept_canopy_layer's
      real(wp),    intent(in)    :: dt                               !< [s]
      real(wp),    intent(inout) :: intercept_leaf(:), intercept_wood(:)   !< [kg/m2 ground/s]
      real(wp),    intent(inout) :: f_wet(:)                         !< [-] combined wetted fraction
      real(wp),    intent(out)   :: throughfall                      !< [kg/m2/s] what reaches the ground
      real(wp)    :: rain_above, combined_w, leaf_cap_i, wood_cap_i, throughfall_i, drip_i
      integer(ik) :: i
      rain_above = rain
      do i = 1_ik, n
         combined_w = leaf_water(i) + wood_water(i)
         leaf_cap_i = leaf_max(i) * lai(i)
         wood_cap_i = wood_max(i) * wai(i)
         call intercept_canopy_layer(combined_w, rain_above, lai(i), wai(i), dt, leaf_max(i), wood_max(i), k_int,  &
                                     alpha_pi, throughfall_i, drip_i, f_wet(i))
         if (leaf_cap_i + wood_cap_i > tiny_num) then
            intercept_leaf(i) = (combined_w*leaf_cap_i/(leaf_cap_i + wood_cap_i) - leaf_water(i)) / dt
            intercept_wood(i) = (combined_w*wood_cap_i/(leaf_cap_i + wood_cap_i) - wood_water(i)) / dt
         else
            intercept_leaf(i) = -leaf_water(i) / dt
            intercept_wood(i) = -wood_water(i) / dt
         end if
         rain_above = throughfall_i   ! cascades to the next (shorter) cohort
      end do
      throughfall = rain_above        ! whatever survives the shortest (last) cohort
   end subroutine canopy_interception

   !---------------------------------------------------------------------------------------!
   ! limit_film_evaporation -- each film's evaporation conductance scaled so that its      !
   ! evaporation over dt, at the rate film_evap the step starts with, cannot exceed what   !
   ! the film holds plus the step's capture. Dew (film_evap <= 0) is left alone.           !
   !---------------------------------------------------------------------------------------!
   pure subroutine limit_film_evaporation(n, dt, leaf_water, wood_water, film_evap_leaf, film_evap_wood, film)
      integer(ik),              intent(in)    :: n
      real(wp),                 intent(in)    :: dt                                !< [s]
      real(wp),                 intent(in)    :: leaf_water(:), wood_water(:)      !< [kg/m2 ground] the films
      real(wp),                 intent(in)    :: film_evap_leaf(:), film_evap_wood(:)   !< [kg/m2 ground/s]
      type(canopy_water_out_t), intent(inout) :: film
      real(wp)    :: avail_leaf, avail_wood
      integer(ik) :: i
      do i = 1_ik, n
         avail_leaf = leaf_water(i) + max(0.0_wp, film%intercept_leaf(i)) * dt
         if (film_evap_leaf(i) > tiny_num) film%g_film_leaf(i) = film%g_film_leaf(i)                          &
              * min(1.0_wp, avail_leaf / (film_evap_leaf(i)*dt))
         avail_wood = wood_water(i) + max(0.0_wp, film%intercept_wood(i)) * dt
         if (film_evap_wood(i) > tiny_num) film%g_film_w(i) = film%g_film_w(i)                                &
              * min(1.0_wp, avail_wood / (film_evap_wood(i)*dt))
      end do
   end subroutine limit_film_evaporation

   !---------------------------------------------------------------------------------------!
   ! clamp_canopy_film -- each cohort's stepped leaf and wood films capped at leaf_max*LAI !
   ! and wood_max*WAI and floored at 0. The excess cut (overflow) and the shortfall filled !
   ! (deficit) are returned so the caller books them in the whole-column ledger.           !
   !---------------------------------------------------------------------------------------!
   pure subroutine clamp_canopy_film(leaf_water, wood_water, lai, wai, leaf_max, wood_max, n, surf_overflow,      &
                                     surf_deficit)
      real(wp),    intent(inout) :: leaf_water(:), wood_water(:)   !< [kg/m2 ground] the films
      real(wp),    intent(in)    :: lai(:), wai(:)
      real(wp),    intent(in)    :: leaf_max(:), wood_max(:)       !< [kg/m2 leaf], [kg/m2 wood] per cohort
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: surf_overflow, surf_deficit    !< [kg/m2] accumulated
      real(wp)    :: leaf_cap_i, wood_cap_i
      integer(ik) :: i
      do i = 1_ik, n
         leaf_cap_i = leaf_max(i) * lai(i)
         wood_cap_i = wood_max(i) * wai(i)
         surf_overflow = surf_overflow + max(0.0_wp, leaf_water(i) - leaf_cap_i)                             &
                                       + max(0.0_wp, wood_water(i) - wood_cap_i)
         surf_deficit  = surf_deficit  + max(0.0_wp, -leaf_water(i))                                         &
                                       + max(0.0_wp, -wood_water(i))
         leaf_water(i) = min(max(leaf_water(i), 0.0_wp), leaf_cap_i)
         wood_water(i) = min(max(wood_water(i), 0.0_wp), wood_cap_i)
      end do
   end subroutine clamp_canopy_film

end module meds_canopy_water
