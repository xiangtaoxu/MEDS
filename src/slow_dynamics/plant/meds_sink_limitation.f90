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
!                                                                                          !
! With both set, the tighter one applies. This is where a more mechanistic sink would go -- one  !
! that responds to temperature or water.                                                     !
!                                                                                          !
! plant_carbon_allocation builds wood up to the limit returned here; the carbon it cannot use    !
! leaves the plant as root exudate, charged the same growth respiration as tissue               !
! (docs/science/plant_carbon_allocation.md).                                                    !
!==========================================================================================!
module meds_sink_limitation
   use meds_kinds,     only : wp
   use meds_allometry, only : size2wood_carbon, dbh_to_height
   implicit none
   private

   public :: growth_sink_limitation

   real(wp), parameter :: DBH_FLOOR = 1.0e-3_wp   !< [cm] keeps dbh^exponent defined for any exponent

contains

   !---------------------------------------------------------------------------------------!
   ! The sink-limited wood growth of one step [kgC/plant]: the most wood carbon the plant can  !
   ! build in dt_yr years. The relative limit is max_relative_growth_rate x wood_carbon x dt;   !
   ! the absolute limit is the wood carbon between the stem's diameter now and that diameter     !
   ! plus max_absolute_growth_rate x dbh^max_absolute_growth_exponent x dt, on the model's own   !
   ! allometry, so a step it binds grows the diameter by exactly that much. A rate of zero or    !
   ! less turns its limit off; with both off there is no limit (huge), which is the default.     !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_sink_limitation(wood_carbon, dbh, wood_density, hgt_max,             &
                                                  aboveground_frac, max_relative_growth_rate,        &
                                                  max_absolute_growth_rate, max_absolute_growth_exponent, &
                                                  dt_yr) result(growth_max)
      real(wp), intent(in) :: wood_carbon               !< [kgC/plant] wood carbon at the step's start
      real(wp), intent(in) :: dbh                       !< [cm] stem diameter at the step's start
      real(wp), intent(in) :: wood_density              !< [g/cm3] the cohort's wood density
      real(wp), intent(in) :: hgt_max                   !< [m] the PFT's height asymptote
      real(wp), intent(in) :: aboveground_frac          !< [--] aboveground share of wood carbon
      real(wp), intent(in) :: max_relative_growth_rate  !< [1/yr] the PFT's maximum relative growth rate of wood
      real(wp), intent(in) :: max_absolute_growth_rate  !< [cm/yr] the PFT's maximum diameter growth of a 1 cm stem
      real(wp), intent(in) :: max_absolute_growth_exponent !< [--] how that maximum scales with dbh [cm]
      real(wp), intent(in) :: dt_yr                     !< [yr] step length
      real(wp)             :: growth_max                !< [kgC/plant] the most wood the step can build
      real(wp) :: dt, dbh_max, wood_now, wood_max

      dt = max(dt_yr, 0.0_wp)
      growth_max = huge(1.0_wp)
      if (max_relative_growth_rate > 0.0_wp)                                                        &
         growth_max = min(growth_max, max_relative_growth_rate * max(wood_carbon, 0.0_wp) * dt)
      if (max_absolute_growth_rate > 0.0_wp) then
         dbh_max    = dbh + max_absolute_growth_rate * max(dbh, DBH_FLOOR) ** max_absolute_growth_exponent * dt
         wood_now   = size2wood_carbon(dbh,     dbh_to_height(dbh,     hgt_max), wood_density, aboveground_frac)
         wood_max   = size2wood_carbon(dbh_max, dbh_to_height(dbh_max, hgt_max), wood_density, aboveground_frac)
         growth_max = min(growth_max, max(wood_max - wood_now, 0.0_wp))
      end if
   end function growth_sink_limitation

end module meds_sink_limitation
