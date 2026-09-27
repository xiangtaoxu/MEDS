! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_lapse_rate -- every VERTICAL correction of the forcing, in one place: the terrain lapse     !
! from the forcing cell's elevation to the site's, and the move from the forcing's own heights to  !
! the top of each patch's canopy air space. Pure/elemental kernels over plain scalars, so the       !
! forcing library still links only meds_shared; the callers (the reader at ingest, the fast loop     !
! per patch) supply the heights.                                                                    !
!==========================================================================================!
module meds_lapse_rate
   use meds_kinds,     only : wp
   use meds_constants, only : grav, r_dry
   implicit none
   private

   public :: wind_log_profile, lapse_air_temperature, lapse_pressure

contains

   !----- Neutral-log wind from the measurement height to the model reference height. The factor !
   !      is independent of u (commutes with the energy-form interpolation); degenerate z0 -> no-op. !
   elemental function wind_log_profile(u_meas, z_meas, z_ref, z0) result(u_ref)
      real(wp), intent(in) :: u_meas, z_meas, z_ref, z0
      real(wp) :: u_ref
      if (z0 <= 0.0_wp .or. z_meas <= z0 .or. z_ref <= z0) then
         u_ref = u_meas                                       ! degenerate: leave the wind untouched
      else
         u_ref = u_meas * log(z_ref / z0) / log(z_meas / z0)
      end if
   end function wind_log_profile

   !----- Linear environmental lapse of air temperature (ED2 calc_met_lapse). dz = site - grid,  !
   !      gamma > 0 cools with height, so a higher site is colder.                                  !
   elemental function lapse_air_temperature(tair_grid, dz, gamma) result(tair_site)
      real(wp), intent(in) :: tair_grid, dz, gamma
      real(wp) :: tair_site
      tair_site = tair_grid - gamma * dz
   end function lapse_air_temperature

   !----- Hydrostatic hypsometric pressure CONSISTENT with the same linear T(z): dP/dz=-Pg/(Rd T),  !
   !      T(z)=T_grid-gamma*z -> P_site = P_grid*(T_site/T_grid)^(g/(Rd*gamma)); isothermal limit    !
   !      (gamma -> 0) is the barometric exp(-g*dz/(Rd*T)). Keeps P and T ideal-gas-consistent.      !
   elemental function lapse_pressure(psurf_grid, tair_grid, dz, gamma) result(psurf_site)
      real(wp), intent(in) :: psurf_grid, tair_grid, dz, gamma
      real(wp) :: psurf_site, tair_site
      real(wp), parameter :: LAPSE_GAMMA_MIN = 1.0e-6_wp
      if (abs(gamma) > LAPSE_GAMMA_MIN) then
         tair_site  = tair_grid - gamma * dz
         psurf_site = psurf_grid * (tair_site / tair_grid) ** (grav / (r_dry * gamma))
      else
         psurf_site = psurf_grid * exp(-grav * dz / (r_dry * tair_grid))
      end if
   end function lapse_pressure

end module meds_lapse_rate
