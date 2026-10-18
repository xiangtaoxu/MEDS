! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_sink_limitation -- how much new tissue a plant can make in a step, whatever its carbon.  !
!                                                                                          !
! Carbon supply is not the only limit on growth: the meristems divide and expand cells at a     !
! bounded rate, however much carbon the leaves fix (sink limitation; Körner 2015). This module   !
! holds that limit on wood growth, in two forms, each a per-PFT trait and each off at zero:      !
!                                                                                          !
!   * RELATIVE (`pft.max_relative_growth_rate`, 1/yr): wood carbon grows by at most this        !
!     fraction of itself a year.                                                            !
!   * ABSOLUTE (`pft.max_absolute_growth_rate`, cm/yr): the stem diameter grows by at most this  !
!     much a year, times dbh^`pft.max_absolute_growth_exponent` (dbh in cm; exponent 0 by        !
!     default, a fixed width of new wood). The cambium lays down a bounded width of wood, so the  !
!     most wood a step can build scales with the stem's surface rather than its mass; at BCI the  !
!     upper quantiles of census diameter growth rise about as dbh^0.5 at a given light.          !
!     The cap is a RATE OVER THE GROWTH WINDOW (growth_memory_days, the one mortality averages     !
!     growth over), not over each step: a stem that grew nothing through a dry season may make up  !
!     the year's diameter in the wet season. The caller passes what the stem has already grown in  !
!     the window, from the cohort's growth ring buffer.                                          !
!                                                                                          !
! With both set, the tighter one applies. This is where a more mechanistic sink would go -- one  !
! that responds to temperature or water.                                                     !
!                                                                                          !
! plant_carbon_allocation builds wood up to the limit returned here; the carbon it cannot use    !
! leaves the plant as root exudate, charged the same growth respiration as tissue              !
! (docs/science/plant_carbon_allocation.md §2 "Sink limitation of wood growth", eq. (4)).      !
!==========================================================================================!
module meds_sink_limitation
   use meds_kinds,     only : wp
   use meds_allometry, only : allometry_t, size2wood_carbon, dbh_to_height
   implicit none
   private

   public :: growth_sink_limitation

   real(wp), parameter :: DBH_FLOOR = 1.0e-3_wp   !< [cm] keeps dbh^exponent defined for any exponent

contains

   !---------------------------------------------------------------------------------------!
   ! The sink-limited wood growth of one step [kgC/plant]: the most wood carbon the plant can  !
   ! build in dt_yr years. The relative limit is max_relative_growth_rate x wood_carbon x dt.   !
   ! The absolute limit is a budget over the growth window: the diameter may grow by            !
   ! max_absolute_growth_rate x dbh^max_absolute_growth_exponent x window_yr, less what it has   !
   ! already grown in the window (grown_in_window), and the limit is the wood carbon between the  !
   ! diameter now and that diameter, on the model's own allometry, so a step it binds brings the  !
   ! window's growth to exactly the cap. window_yr = dt_yr and grown_in_window = 0 is a per-step  !
   ! cap. A rate of zero or less turns its limit off; with both off there is no limit (huge),    !
   ! which is the default.                                                                     !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_sink_limitation(allom, wood_carbon, dbh, wood_density, hgt_max,      &
                                                  aboveground_frac, max_relative_growth_rate,        &
                                                  max_absolute_growth_rate, max_absolute_growth_exponent, &
                                                  dt_yr, window_yr, grown_in_window) result(growth_max)
      type(allometry_t), intent(in) :: allom            !< the site's allometry coefficients
      real(wp), intent(in) :: wood_carbon               !< [kgC/plant] wood carbon at the step's start
      real(wp), intent(in) :: dbh                       !< [cm] stem diameter at the step's start
      real(wp), intent(in) :: wood_density              !< [g/cm3] the cohort's wood density
      real(wp), intent(in) :: hgt_max                   !< [m] the PFT's height asymptote
      real(wp), intent(in) :: aboveground_frac          !< [--] aboveground share of wood carbon
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
         wood_now   = size2wood_carbon(allom, dbh,     dbh_to_height(allom, dbh,     hgt_max), wood_density,      &
                                       aboveground_frac)
         wood_max   = size2wood_carbon(allom, dbh_max, dbh_to_height(allom, dbh_max, hgt_max), wood_density,      &
                                       aboveground_frac)
         growth_max = min(growth_max, max(wood_max - wood_now, 0.0_wp))
      end if
   end function growth_sink_limitation

end module meds_sink_limitation
