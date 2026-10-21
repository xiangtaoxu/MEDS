! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_forcing_kernels -- the pure and elemental meteorological-forcing math                       !
! (docs/science/forcing.md §3, §5 and §7): the energy-form wind, the interval-mean-conserving       !
! shortwave disaggregation, humidity from dewpoint or RH and the rain/snow split. Uses only         !
! meds_share, so the saturation curve, the solar geometry (meds_time) and the radiation             !
! (meds_radiation_lib) are the model's own. Device-safe.                                            !
!==========================================================================================!
module meds_forcing_kernels
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : t_3ple, tiny_num
   use meds_therm_lib,      only : sat_vapor_pressure
   use meds_time,           only : meds_time_t, mean_solar_cosz
   use meds_radiation_lib,  only : COSZ_MIN
   implicit none
   private

   public :: interpolate_wind_energy
   public :: cosz_reconstruct_factor, disaggregate_sw
   public :: dewpoint_to_specific_humidity, rh_to_specific_humidity
   public :: precip_phase

   real(wp), parameter :: COSZ_BAR_MIN   = 1.0e-3_wp   !< [-] window-mean-cosz floor for a valid reconstruction
   real(wp), parameter :: PHASE_BAND_K   = 1.0_wp      !< [K] half-width of the rain/snow phase-transition band


contains

   !----- Wind interpolates as its square, so its kinetic energy varies linearly              !
   !      (docs/science/forcing.md eq. (1)). No floor here: the aerodynamics floors the speed  !
   !      at aerodynamics.ubmin where it uses it. --------------------------------------------!
   elemental function interpolate_wind_energy(prev, next, w_next) result(u)
      real(wp), intent(in) :: prev, next, w_next
      real(wp) :: u
      u = sqrt(max((1.0_wp - w_next) * prev * prev + w_next * next * next, 0.0_wp))
   end function interpolate_wind_energy

   !=======================================================================================!
   !  SHORTWAVE DISAGGREGATION (docs/science/forcing.md eqs. (5)-(6)): the factor 1/<cosz>_win  !
   !  that turns an interval mean F_avg into F(now) = F_avg*cosz(now)*factor, whose mean over the !
   !  same window is F_avg. It is the reciprocal of the mean cosine (night sub-samples as 0), not  !
   !  the mean secant, which is biased high at low sun. A fully night window returns 0.           !
   !=======================================================================================!
   pure function cosz_reconstruct_factor(t, win_start_sec, dt_sub, dt_win,                      &
                                         latitude_deg, longitude_deg) result(factor)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: win_start_sec, dt_sub, dt_win
      real(wp),          intent(in) :: latitude_deg, longitude_deg
      real(wp)    :: factor, cosz_bar
      cosz_bar = mean_solar_cosz(t, win_start_sec, dt_sub, dt_win, latitude_deg, longitude_deg)
      if (cosz_bar > COSZ_BAR_MIN) then
         factor = 1.0_wp / cosz_bar
      else
         factor = 0.0_wp                                  ! fully-night window; F_avg is 0 too
      end if
   end function cosz_reconstruct_factor

   !----- Instantaneous SW stream from its interval mean: F_now = F_avg*cosz(now)*factor. ------!
   !      Night (cosz(now) <= COSZ_MIN) -> 0 automatically; factor already 0 for a night window. !
   elemental function disaggregate_sw(f_avg, cosz_now, factor) result(f_now)
      real(wp), intent(in) :: f_avg, cosz_now, factor
      real(wp) :: f_now
      if (cosz_now > COSZ_MIN) then
         f_now = f_avg * cosz_now * factor
      else
         f_now = 0.0_wp
      end if
   end function disaggregate_sw

   !=======================================================================================!
   !  HUMIDITY (docs/science/forcing.md eq. (11)). Both forms use meds_therm_lib%sat_vapor_pressure !
   !  over liquid water (Bolton 1980), the curve of the rest of MEDS and of the flux-tower         !
   !  formatter, so the reader adds no saturation curve of its own.                               !
   !=======================================================================================!
   !----- Specific humidity from dewpoint: the actual e is e_sat at the dewpoint. -------------!
   elemental function dewpoint_to_specific_humidity(td_k, p_pa) result(q)
      real(wp), intent(in) :: td_k, p_pa
      real(wp) :: q, e
      e = sat_vapor_pressure(td_k)
      q = 0.622_wp * e / max(p_pa - 0.378_wp * e, tiny_num)
   end function dewpoint_to_specific_humidity

   !----- Specific humidity from relative humidity (fraction) + temperature + pressure. ------!
   elemental function rh_to_specific_humidity(rh, t_k, p_pa) result(q)
      real(wp), intent(in) :: rh, t_k, p_pa
      real(wp) :: q, e
      e = min(1.0_wp, max(0.0_wp, rh)) * sat_vapor_pressure(t_k)
      q = 0.622_wp * e / max(p_pa - 0.378_wp * e, tiny_num)
   end function rh_to_specific_humidity

   !=======================================================================================!
   !  PRECIP PHASE (docs/science/forcing.md eq. (12)): a linear rain/snow ramp across a 2 K band  !
   !  centred on the triple point. Mass-conserving: rainf + snowfall = total.                     !
   !=======================================================================================!
   elemental subroutine precip_phase(precip_total, tair_k, rainf, snowfall)
      real(wp), intent(in)  :: precip_total, tair_k
      real(wp), intent(out) :: rainf, snowfall
      real(wp) :: frac_liq
      frac_liq = min(1.0_wp, max(0.0_wp, (tair_k - (t_3ple - PHASE_BAND_K)) / (2.0_wp * PHASE_BAND_K)))
      rainf = frac_liq * precip_total
      snowfall = precip_total - rainf
   end subroutine precip_phase

end module meds_forcing_kernels
