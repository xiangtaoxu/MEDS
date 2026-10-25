! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_sink_limitation -- the most wood a cohort can build in a step, whatever its carbon: !
! the meristems divide and expand cells at a bounded rate (sink limitation; Körner 2015).  !
! growth_sink_limitation returns the tighter of two per-PFT limits, each off at zero:      !
!   * relative: wood carbon grows by at most max_relative_growth_rate of itself a year;    !
!   * absolute: the diameter grows by at most max_absolute_growth_rate x                   !
!     dbh^max_absolute_growth_exponent a year (dbh in cm), a budget over the growth window. !
! plant_carbon_allocation builds wood up to the limit and exudes the rest from the roots   !
! (docs/science/plant_carbon_allocation.md §2, "Sink limitation of wood growth", eq. (4)). !
!==========================================================================================!
module meds_sink_limitation
   use meds_kinds,         only : wp
   use meds_allometry_lib, only : allometry_t, size_to_wood_carbon, dbh_to_height
   implicit none
   private

   public :: growth_sink_limitation

   real(wp), parameter :: DBH_FLOOR = 1.0e-3_wp   !< [cm] keeps dbh^exponent defined for any exponent

contains

   !---------------------------------------------------------------------------------------!
   ! The sink-limited wood growth of one step [kgC/plant]. The absolute limit is the wood  !
   ! carbon, on the model's own allometry, between today's diameter and today's plus the   !
   ! window's budget less grown_in_window, so a step it binds brings the window's growth   !
   ! to exactly the cap; window_yr = dt_yr with grown_in_window = 0 is a per-step cap.     !
   ! With both limits off the result is huge(1.0_wp): no limit.                            !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_sink_limitation(allom, wood_carbon, dbh, wood_density,                &
                                                  max_relative_growth_rate,                          &
                                                  max_absolute_growth_rate, max_absolute_growth_exponent, &
                                                  dt_yr, window_yr, grown_in_window) result(growth_max)
      type(allometry_t), intent(in) :: allom            !< the cohort's PFT's allometry
      real(wp), intent(in) :: wood_carbon               !< [kgC/plant] wood carbon at the step's start
      real(wp), intent(in) :: dbh                       !< [cm] stem diameter at the step's start
      real(wp), intent(in) :: wood_density              !< [g/cm3] the cohort's wood density
      real(wp), intent(in) :: max_relative_growth_rate  !< [1/yr] the PFT's maximum relative growth rate of wood
      real(wp), intent(in) :: max_absolute_growth_rate  !< [cm/yr] the PFT's maximum diameter growth of a 1 cm stem
      real(wp), intent(in) :: max_absolute_growth_exponent !< [--] how that maximum scales with dbh [cm]
      real(wp), intent(in) :: dt_yr                     !< [yr] step length
      real(wp), intent(in) :: window_yr                 !< [yr] the growth window this step closes (>= dt_yr)
      real(wp), intent(in) :: grown_in_window           !< [cm] diameter already grown in it, before this step
      real(wp)             :: growth_max                !< [kgC/plant] the most wood the step can build
      real(wp) :: dt, dbh_room, dbh_max, wood_now, wood_max

      dt = max(dt_yr, 0.0_wp)
      growth_max = huge(1.0_wp)
      if (max_relative_growth_rate > 0.0_wp)                                                        &
         growth_max = min(growth_max, max_relative_growth_rate * max(wood_carbon, 0.0_wp) * dt)
      if (max_absolute_growth_rate > 0.0_wp) then
         dbh_room   = max_absolute_growth_rate * max(dbh, DBH_FLOOR) ** max_absolute_growth_exponent       &
                    * max(window_yr, dt) - max(grown_in_window, 0.0_wp)
         dbh_max    = dbh + max(dbh_room, 0.0_wp)
         wood_now   = size_to_wood_carbon(allom, dbh,     dbh_to_height(allom, dbh),     wood_density)
         wood_max   = size_to_wood_carbon(allom, dbh_max, dbh_to_height(allom, dbh_max), wood_density)
         growth_max = min(growth_max, max(wood_max - wood_now, 0.0_wp))
      end if
   end function growth_sink_limitation

end module meds_sink_limitation
