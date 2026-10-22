! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_site_terrain -- a site's forcing within its polygon (docs/science/forcing.md §8).   !
!                                                                                          !
! 1. The sun on a slope: level ground sees cos z; at a mid-latitude winter noon a slope     !
!    facing the equator sees the sun higher than level ground does and one facing the pole  !
!    lower; a slope facing east sees more of the morning sun than of the afternoon's.       !
! 2. The beam on a slope: level ground keeps its beam, no beam stays none, and the beam     !
!    normal to the sun is held at the solar constant.                                       !
! 3. A site above its polygon is colder by the lapse rate times the height, at a lower      !
!    pressure and the same relative humidity, under less longwave; a site at the polygon's  !
!    elevation on level ground keeps the polygon's forcing exactly, and a run without the   !
!    lapse keeps its air and scales only the beam.                                          !
!==========================================================================================!
program test_site_terrain
   use meds_test_assert,     only : check_true, test_report
   use meds_kinds,           only : wp
   use meds_time,            only : meds_time_t, seconds_into_day, solar_cosz_utc, slope_incidence_cosz
   use meds_therm_lib,       only : sat_vapor_pressure, air_density
   use meds_forcing_params,  only : forcing_params_t
   use meds_met_types,   only : met_forcing_t
   use meds_radiation_lib,   only : slope_beam_factor
   use meds_met,             only : met_at_site
   implicit none

   real(wp), parameter :: LAT = 45.0_wp, LON = 0.0_wp, S0 = 1361.0_wp
   type(meds_time_t)      :: noon, morning, afternoon
   type(forcing_params_t) :: f
   type(met_forcing_t)    :: met, site
   real(wp) :: cosz, flat, south, north, east_am, east_pm

   !=== 1. The sun on a slope (45 N on the winter solstice; at 0 E, solar noon is ~12 UTC). ====!
   noon      = meds_time_t(year=2021, month=12, day=21, hour=12)
   morning   = meds_time_t(year=2021, month=12, day=21, hour=10)
   afternoon = meds_time_t(year=2021, month=12, day=21, hour=14)
   cosz  = solar_cosz_utc(noon, seconds_into_day(noon), LAT, LON)
   flat  = slope_incidence_cosz(noon, seconds_into_day(noon), LAT, LON, 0.0_wp, 0.0_wp)
   south = slope_incidence_cosz(noon, seconds_into_day(noon), LAT, LON, 30.0_wp, 180.0_wp)
   north = slope_incidence_cosz(noon, seconds_into_day(noon), LAT, LON, 30.0_wp, 0.0_wp)
   call check_true('level ground sees the sun at cos z', abs(flat - cosz) < 1.0e-12_wp)
   call check_true('a slope facing the equator sees the winter noon sun higher than level ground', south > cosz)
   call check_true('a slope facing the pole sees it lower', north < cosz)
   east_am = slope_incidence_cosz(morning, seconds_into_day(morning), LAT, LON, 30.0_wp, 90.0_wp)
   east_pm = slope_incidence_cosz(afternoon, seconds_into_day(afternoon), LAT, LON, 30.0_wp, 90.0_wp)
   call check_true('a slope facing east sees more of the morning sun than of the afternoon''s', east_am > east_pm)

   !=== 2. The beam on a slope. =================================================================!
   call check_true('level ground keeps its beam', slope_beam_factor(300.0_wp, 0.5_wp, 0.5_wp) == 1.0_wp)
   call check_true('no beam has nothing to scale', slope_beam_factor(0.0_wp, 0.5_wp, 0.7_wp) == 1.0_wp)
   call check_true('a slope facing the sun takes cos i / cos z',                                  &
                   abs(slope_beam_factor(300.0_wp, 0.5_wp, 0.7_wp) - 1.4_wp) < 1.0e-12_wp)
   call check_true('the beam normal to a low sun is held at the solar constant',                  &
                   abs(slope_beam_factor(100.0_wp, 0.01_wp, 0.5_wp) - S0 * 0.5_wp / 100.0_wp) < 1.0e-9_wp)

   !=== 3. A site above its polygon. ===========================================================!
   f%apply_elevation_lapse = .true.
   f%lapse_rate_tair = 0.0065_wp
   met%tair_k = 290.0_wp
   met%psurf_pa = 95000.0_wp
   met%qair = 0.010_wp
   met%lwdown = 350.0_wp
   met%par_beam = 200.0_wp
   met%nir_beam = 220.0_wp
   met%par_diffuse = 60.0_wp
   met%rho_air = air_density(met%tair_k, met%psurf_pa, met%qair)
   site = met_at_site(met, f, 1, 500.0_wp, 1.0_wp)
   call check_true('500 m up is 500 m x the lapse rate colder', abs(site%tair_k - (290.0_wp - 3.25_wp)) < 1.0e-12_wp)
   call check_true('and at a lower pressure', site%psurf_pa < met%psurf_pa)
   call check_true('at the same relative humidity',                                               &
                   abs(rh(site) - rh(met)) < 1.0e-9_wp)
   call check_true('under less longwave', site%lwdown < met%lwdown)
   call check_true('with its air density from its own air',                                       &
                   site%rho_air == air_density(site%tair_k, site%psurf_pa, site%qair))
   call check_true('the light is the polygon''s on level ground',                                  &
                   site%par_beam == met%par_beam .and. site%par_diffuse == met%par_diffuse)
   site = met_at_site(met, f, 1, 0.0_wp, 1.0_wp)
   call check_true('a site at the polygon''s elevation on level ground keeps its forcing exactly',  &
                   site%tair_k == met%tair_k .and. site%psurf_pa == met%psurf_pa .and.                &
                   site%qair == met%qair .and. site%lwdown == met%lwdown .and.                        &
                   site%par_beam == met%par_beam .and. site%nir_beam == met%nir_beam .and.            &
                   site%rho_air == met%rho_air)
   f%apply_elevation_lapse = .false.
   site = met_at_site(met, f, 1, 500.0_wp, 1.5_wp)
   call check_true('without the lapse the air is the polygon''s',                                  &
                   site%tair_k == met%tair_k .and. site%psurf_pa == met%psurf_pa)
   call check_true('and the slope scales the beam alone',                                          &
                   site%par_beam == 1.5_wp * met%par_beam .and. site%nir_beam == 1.5_wp * met%nir_beam &
                   .and. site%par_diffuse == met%par_diffuse)

   call test_report('test_site_terrain')

contains

   !----- Relative humidity of a forcing sample: its vapour pressure over saturation. ----------!
   real(wp) function rh(m)
      type(met_forcing_t), intent(in) :: m
      rh = m%qair * m%psurf_pa / (0.622_wp + 0.378_wp * m%qair) / sat_vapor_pressure(m%tair_k)
   end function rh

end program test_site_terrain
