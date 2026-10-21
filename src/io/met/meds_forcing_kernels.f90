! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_forcing_kernels -- the pure and elemental meteorological-forcing math                       !
! (docs/science/forcing.md §3 and §5): the energy-form wind and the interval-mean-conserving        !
! shortwave disaggregation, on the solar geometry of meds_time. Device-safe.                        !
!==========================================================================================!
module meds_forcing_kernels
   use meds_kinds,          only : wp, ik
   use meds_time,           only : meds_time_t, mean_solar_cosz
   use meds_radiation_lib,  only : COSZ_MIN
   implicit none
   private

   public :: interpolate_wind_energy
   public :: cosz_reconstruct_factor, disaggregate_sw

   real(wp), parameter :: COSZ_BAR_MIN   = 1.0e-3_wp   !< [-] window-mean-cosz floor for a valid reconstruction


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

end module meds_forcing_kernels
