! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_lapse_rate -- every VERTICAL correction of the forcing, in one place (docs/science/        !
! forcing.md §8). Two steps, in this order:                                                       !
!                                                                                          !
!   1. TERRAIN, per record at ingest (the reader): from the forcing cell's elevation to the site's. !
!      Temperature by the month's environmental lapse rate, pressure hydrostatically with the same  !
!      linear T(z), specific humidity at constant RELATIVE humidity, and file longwave by the ratio !
!      of clear-sky emission eps*T^4 (NLDAS, Cosgrove et al. 2003; ERA5-Land's own terrain step    !
!      also holds relative humidity).                                                              !
!   2. HEIGHT, per patch and sub-step (the fast loop): from the forcing's own heights to the top of !
!      the patch's canopy air space, z_top = can_depth, which grows with the stand. Neutral surface !
!      layer: potential temperature and specific humidity are conserved, and the wind follows the  !
!      patch's own log profile (the one its aerodynamics starts from). An open-terrain wind        !
!      diagnostic (ERA5's 10 m wind, made from a 40 m blending height with z0 = 0.03 m) is first   !
!      returned to its blending height. Pressure stays at the ground, where the canopy air, ground !
!      and leaves use it.                                                                          !
!                                                                                          !
! Pure/elemental kernels over plain scalars, so the forcing library still links only meds_shared    !
! and the config leaf; the callers supply the heights.                                             !
!==========================================================================================!
module meds_lapse_rate
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : grav, r_dry, cp_air
   use meds_therm_lib,       only : sat_vapor_pressure, air_density
   use meds_forcing_config,  only : forcing_config_t, HEIGHT_ABOVE_ZERO_PLANE, WIND_EXPOSURE_OPEN_TERRAIN
   use meds_forcing_kernels, only : clear_sky_emissivity, rh_to_specific_humidity
   use meds_forcing_types,   only : met_forcing_t
   implicit none
   private

   public :: lapse_air_temperature, lapse_pressure
   public :: monthly_lapse_rate, lapse_specific_humidity, lapse_longwave
   public :: cas_top_wind_factor, cas_top_air_temperature, met_to_cas_top

contains

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

   !=======================================================================================!
   !  1. TERRAIN (per record, at ingest).                                                         !
   !=======================================================================================!
   !----- The environmental lapse rate for a calendar month, from the twelve configured values. -!
   pure real(wp) function monthly_lapse_rate(gamma_month, month) result(gamma)
      real(wp),    intent(in) :: gamma_month(12)   !< [K/m] January .. December
      integer(ik), intent(in) :: month             !< 1 .. 12
      gamma = gamma_month(min(12_ik, max(1_ik, month)))
   end function monthly_lapse_rate

   !----- Specific humidity at the site, holding RELATIVE humidity across the lapse. Holding q   !
   !      instead would dry a site below its cell by ~6 % RH per K of warming (20 % over 500 m).   !
   elemental function lapse_specific_humidity(q_grid, t_grid, p_grid, t_site, p_site) result(q_site)
      real(wp), intent(in) :: q_grid, t_grid, p_grid, t_site, p_site
      real(wp) :: q_site, e, rh
      e      = max(q_grid, 0.0_wp) * p_grid / (0.622_wp + 0.378_wp * max(q_grid, 0.0_wp))
      rh     = e / sat_vapor_pressure(t_grid)
      q_site = rh_to_specific_humidity(rh, t_site, p_site)
   end function lapse_specific_humidity

   !----- File longwave at the site: scaled by the clear-sky emission eps*T^4 at the site over the   !
   !      cell (NLDAS), with the same clear-sky emissivity the longwave synthesis uses.             !
   pure function lapse_longwave(lw_grid, form, t_grid, q_grid, p_grid, t_site, q_site, p_site) result(lw_site)
      real(wp),    intent(in) :: lw_grid                   !< [W/m2] downward longwave at the cell
      integer(ik), intent(in) :: form                      !< LW_CLEAR_BRUTSAERT | LW_CLEAR_IDSO
      real(wp),    intent(in) :: t_grid, q_grid, p_grid, t_site, q_site, p_site
      real(wp) :: lw_site
      lw_site = lw_grid * (clear_sky_emissivity(form, t_site, q_site, p_site) * t_site**4)      &
                        / (clear_sky_emissivity(form, t_grid, q_grid, p_grid) * t_grid**4)
   end function lapse_longwave

   !=======================================================================================!
   !  2. HEIGHT (per patch and sub-step): the forcing's own heights -> the canopy-air top.        !
   !=======================================================================================!
   !----- The factor taking the forcing's wind to z_top over a patch of roughness z0 and zero-plane  !
   !      displacement d. An open-terrain diagnostic is first returned to its blending height z_b,   !
   !      then the patch's neutral log profile carries it from z_b to z_top:                        !
   !          u(z_top) = u_m ln(z_b/z0e)/ln(z_m/z0e) * ln(h_top/z0)/ln(h_b/z0)                      !
   !      with h the height above d (a reanalysis's heights already are; a tower's are above the    !
   !      ground), floored at 2 z0 like the aerodynamics' own reference height.                     !
   pure real(wp) function cas_top_wind_factor(f, z_top, displace, rough) result(factor)
      type(forcing_config_t), intent(in) :: f
      real(wp),               intent(in) :: z_top      !< [m] canopy-air top, above the ground
      real(wp),               intent(in) :: displace   !< [m] zero-plane displacement
      real(wp),               intent(in) :: rough      !< [m] roughness length (> 0)
      real(wp) :: z_b, h_b, h_top
      factor = 1.0_wp ; z_b = f%wind_height
      if (f%wind_exposure == WIND_EXPOSURE_OPEN_TERRAIN) then
         factor = log(f%wind_blending_height / f%wind_exposure_z0) / log(f%wind_height / f%wind_exposure_z0)
         z_b    = f%wind_blending_height
      end if
      h_b = z_b
      if (f%height_above /= HEIGHT_ABOVE_ZERO_PLANE) h_b = z_b - displace
      h_b    = max(h_b, 2.0_wp * rough)
      h_top  = max(z_top - displace, 2.0_wp * rough)
      factor = factor * log(h_top / rough) / log(h_b / rough)
   end function cas_top_wind_factor

   !----- Air temperature at z_top, conserving potential temperature from the forcing's height (a  !
   !      reanalysis's is above d): T(z_top) = T - (g/cp)(z_top - z_T). The aerodynamics' theta =   !
   !      T + (g/cp) z uses the same constants, so its theta is exactly the forcing's.               !
   pure real(wp) function cas_top_air_temperature(f, tair, z_top, displace) result(t_top)
      type(forcing_config_t), intent(in) :: f
      real(wp),               intent(in) :: tair       !< [K] at the forcing's temperature height
      real(wp),               intent(in) :: z_top      !< [m] canopy-air top, above the ground
      real(wp),               intent(in) :: displace   !< [m] zero-plane displacement
      real(wp) :: z_t
      z_t = f%tq_height
      if (f%height_above == HEIGHT_ABOVE_ZERO_PLANE) z_t = displace + f%tq_height
      t_top = tair - (grav / cp_air) * (z_top - z_t)
   end function cas_top_air_temperature

   !----- A forcing record at a patch's canopy-air top: wind (and its vector) by the patch's factor, !
   !      temperature along the dry adiabat, air density re-derived. Humidity is conserved and the   !
   !      radiation, rain and pressure are unchanged.                                                 !
   pure function met_to_cas_top(met, f, z_top, displace, rough) result(top)
      type(met_forcing_t),    intent(in) :: met
      type(forcing_config_t), intent(in) :: f
      real(wp),               intent(in) :: z_top, displace, rough
      type(met_forcing_t) :: top
      real(wp) :: factor
      top    = met
      factor = cas_top_wind_factor(f, z_top, displace, rough)
      top%wind = factor * met%wind
      if (met%has_wind_vector) then
         top%wind_u = factor * met%wind_u
         top%wind_v = factor * met%wind_v
      end if
      top%tair_k  = cas_top_air_temperature(f, met%tair_k, z_top, displace)
      top%rho_air = air_density(top%tair_k, top%psurf_pa, top%qair)
   end function met_to_cas_top

end module meds_lapse_rate
