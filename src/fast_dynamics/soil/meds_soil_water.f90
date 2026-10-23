! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_water -- how water moves through and out of the soil column, each flux with    !
! its derivatives for the water network (docs/science/soil_biophysics.md eqs. (1)-(4)):    !
!   * soil_water_flux         -- Darcy flow between two layers, or the pond and the top    !
!                                layer, eqs. (3)-(4);                                      !
!   * ground_evaporation_rate -- evaporation from the top layer, on the terms              !
!                                ground_evaporation_setup sets once per water step;        !
!   * bottom_flux             -- drainage through the bottom face, or a head-driven aquifer.!
! A layer's water relations are in meds_soil_hydr_lib, root uptake in meds_plant_hydr_lib. !
!==========================================================================================!
module meds_soil_water
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav, grav_head, r_wv, tiny_num, p_std, safe_exp, t_kelvin
   use meds_soil_params,      only : soil_profile_t, soil_params_t, SOIL_BC_AQUIFER
   use meds_soil_hydr_lib,    only : soil_head_to_theta
   use meds_numerics,         only : smooth_max0
   use meds_therm_lib,        only : sat_specific_humidity
   implicit none
   private

   public :: soil_water_flux, ground_evaporation_rate, bottom_flux
   public :: ground_evap_env_t, ground_evaporation_setup

   !----- Widths of the two rounded corners in ground_evaporation_rate: [kg/kg] where the surface's  !
   !      humidity meets the air's (evaporation stops; there is no dew), and [m3/m3] where the top   !
   !      layer dries below the moisture at which the dry surface layer starts to form. -----------!
   real(wp), parameter :: EVAP_HUMIDITY_WIDTH = 1.0e-6_wp
   real(wp), parameter :: DSL_THETA_WIDTH     = 1.0e-3_wp
   !----- CLM5 eq 5.78: the 'air dry' matric head [m], where the soil counts as dry (CLM5's 1e7 mm). !
   real(wp), parameter :: PSI_AIR_DRY = -1.0e4_wp

   !----- The boundary ground_evaporation_rate reads: the ground skin and the canopy air over it,    !
   !      and the terms that depend only on them and the top layer's curve. -------------------------!
   type :: ground_evap_env_t
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temperature
      !----- A frozen surface evaporates over ice (docs/science/snow_biophysics.md eq. (1)). ---------!
      real(wp) :: ground_fliq = 1.0_wp                    !< [-] top-layer liquid fraction
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance of the bare-soil tile
                                                          !<       (1/ggnet), not inflated for snow
      real(wp) :: snow_free_frac = 1.0_wp                 !< [-] snow-free area fraction (1-snowfac) that evaporates,
                                                          !<     an area-weighted tile flux (snow_biophysics.md)
      !----- Set by ground_evaporation_setup from the fields above and the top layer's curve. -------!
      real(wp) :: rv_t            = 0.0_wp   !< [J/kg] R_v T of the ground skin, the scale of Philip's humidity
      real(wp) :: q_sat           = 0.0_wp   !< [kg/kg] saturation specific humidity at the ground skin
      real(wp) :: theta_init      = 0.0_wp   !< [m3/m3] top-layer water content below which the dry layer forms
      real(wp) :: theta_span      = 0.0_wp   !< [m3/m3] theta_init less the air-dry water content
      real(wp) :: dsl_diffusivity = 0.0_wp   !< [m2/s] vapour diffusivity through the dry surface layer
   end type ground_evap_env_t

contains

   !----- The Darcy flux of soil water [kg/m2/s] from node a to node b (two layers, or the pond and  !
   !      the top layer), with its derivatives (docs/science/soil_biophysics.md eqs. (3)-(4)): the     !
   !      thickness-weighted geometric-mean face conductivity times the potential difference [MPa]  !
   !      plus the elevation drop z_drop [m], over the node distance [m]. ----------------------------!
   elemental subroutine soil_water_flux(k_a, k_b, w_b, psi_a, psi_b, z_drop, dist, q,               &
                                        dq_dpsi_a, dq_dpsi_b, dq_dk_a, dq_dk_b)
      real(wp), intent(in)  :: k_a, k_b, w_b, psi_a, psi_b, z_drop, dist
      real(wp), intent(out) :: q, dq_dpsi_a, dq_dpsi_b, dq_dk_a, dq_dk_b
      real(wp) :: k_face, drive
      k_face = k_a ** (1.0_wp - w_b) * k_b ** w_b
      drive  = ((psi_a - psi_b) / grav_head + z_drop) / dist
      q      = rho_h2o * k_face * drive
      dq_dpsi_a =  rho_h2o * k_face / (grav_head * dist)
      dq_dpsi_b = -dq_dpsi_a
      dq_dk_a   = (1.0_wp - w_b) * q / max(k_a, tiny_num)
      dq_dk_b   = w_b * q / max(k_b, tiny_num)
   end subroutine soil_water_flux

   !---------------------------------------------------------------------------------------!
   ! The bottom-face Darcy flux [m/s, down-positive] (docs/science/soil_biophysics.md,     !
   ! "Boundary conditions, and what the aquifer option actually means"): drainage          !
   ! q = drain_factor * K(theta_n), or the aquifer boundary q = K_bot (psi_n/Delta + 1),   !
   ! Delta = dz(n)/2, with K_bot upstream: K(theta_n) downward, K_sat upward.              !
   !---------------------------------------------------------------------------------------!
   pure subroutine bottom_flux(profile, soil, n, psi_n, k_n, qbot, dq_dpsi, dq_dk)
      type(soil_profile_t), intent(in)  :: profile
      type(soil_params_t),  intent(in)  :: soil
      integer(ik),          intent(in)  :: n
      real(wp),             intent(in)  :: psi_n    !< [m]   matric head of the deepest layer
      real(wp),             intent(in)  :: k_n      !< [m/s] K(theta_n) of the deepest layer
      real(wp),             intent(out) :: qbot     !< [m/s] down-positive
      real(wp), optional,   intent(out) :: dq_dpsi  !< [1/s] d(qbot)/d(psi_n), for the implicit solve
      real(wp), optional,   intent(out) :: dq_dk    !< [-]   d(qbot)/d(k_n), for the water network
      real(wp) :: delta, grad, k_bot
      qbot = 0.0_wp
      if (present(dq_dpsi)) dq_dpsi = 0.0_wp
      if (present(dq_dk))   dq_dk   = 0.0_wp
      if (soil%water_bottom_bc == SOIL_BC_AQUIFER) then
         delta = 0.5_wp * profile%dz(n)                                  ! z_node(n) - z_bottom > 0
         grad  = psi_n / delta + 1.0_wp
         if (grad >= 0.0_wp) then
            k_bot = k_n                                                 ! downward: upstream is layer n
            if (present(dq_dk)) dq_dk = grad
         else
            k_bot = profile%ksat(n)                                      ! upward: upstream is saturated
         end if
         qbot = k_bot * grad
         if (present(dq_dpsi)) dq_dpsi = k_bot / delta
      else                                                              ! drainage, q = f K(theta_n)
         qbot = soil%drain_factor * k_n
         if (present(dq_dk)) dq_dk = soil%drain_factor
      end if
   end subroutine bottom_flux

   !----- The terms of ground_evaporation_rate that depend only on its boundary and the top layer's  !
   !      curve: R_v T and the saturation humidity at the ground skin (over ice when frozen), and    !
   !      the dry surface layer's water-content span and vapour diffusivity. The water network holds !
   !      the boundary over its step, so water_step sets them once, before its first evaluation. ----!
   pure subroutine ground_evaporation_setup(profile, soil, forcing)
      type(soil_profile_t),    intent(in)    :: profile
      type(soil_params_t),     intent(in)    :: soil
      type(ground_evap_env_t), intent(inout) :: forcing
      real(wp) :: theta_air, dvap, phi, phi_air, tau
      forcing%rv_t  = r_wv * forcing%t_ground
      forcing%q_sat = sat_specific_humidity(forcing%t_ground, p_std, forcing%ground_fliq)
      !----- The dry surface layer, referenced to the air-dry water content. ---------------------!
      theta_air = soil_head_to_theta(profile%retention, PSI_AIR_DRY, profile%theta_sat(1), profile%theta_res(1),  &
                                     profile%curve_par_a(1), profile%curve_par_n(1))
      forcing%theta_init = soil%dsl_theta_init * profile%theta_sat(1)
      dvap    = 2.12e-5_wp * (forcing%t_ground / t_kelvin) ** 1.75_wp
      phi     = profile%theta_sat(1)
      phi_air = max(phi - theta_air, tiny_num)
      tau     = phi_air ** (10.0_wp / 3.0_wp) / max(phi * phi, tiny_num)
      forcing%theta_span      = max(forcing%theta_init - theta_air, tiny_num)
      forcing%dsl_diffusivity = max(dvap * tau, tiny_num)
   end subroutine ground_evaporation_setup

   !---------------------------------------------------------------------------------------!
   ! ground_evaporation_rate -- ground evaporation [kg/m2/s] from the top layer at water content        !
   ! theta1 [m3/m3] and matric head psi1 [m], with its derivatives (docs/science/soil_biophysics.md,    !
   ! "The soil water column"): Philip's pore-space humidity (over ice when frozen) through the          !
   ! Swenson-Lawrence dry-surface-layer resistance (CLM5 eqs. 5.77-5.80) in series with r_aero, from    !
   ! the snow-free share only (docs/science/snow_biophysics.md). Both corners are rounded for Newton.   !
   ! The boundary's constant terms are ground_evaporation_setup's.                                      !
   !---------------------------------------------------------------------------------------!
   pure subroutine ground_evaporation_rate(theta1, psi1, forcing, soil, e_soil, de_dtheta, de_dpsi)
      real(wp),               intent(in)  :: theta1      !< [m3/m3] top-layer water content
      real(wp),               intent(in)  :: psi1        !< [m]     top-layer matric head
      type(ground_evap_env_t), intent(in) :: forcing
      type(soil_params_t),    intent(in)  :: soil
      real(wp),               intent(out) :: e_soil      !< [kg/m2/s]
      real(wp),               intent(out) :: de_dtheta   !< [kg/m2/s per m3/m3]
      real(wp),               intent(out) :: de_dpsi     !< [kg/m2/s per m of head]
      real(wp) :: alpha_soil, dalpha, drive, ddrive, dry, ddry, r_soil, dr_dtheta, resist
      alpha_soil = safe_exp(psi1 * grav / forcing%rv_t)
      dalpha     = alpha_soil * grav / forcing%rv_t
      call smooth_max0(alpha_soil * forcing%q_sat - forcing%q_air, EVAP_HUMIDITY_WIDTH, drive, ddrive)
      !----- The dry surface layer, from theta_init down to the air-dry water content. ----------!
      call smooth_max0(forcing%theta_init - theta1, DSL_THETA_WIDTH, dry, ddry)
      r_soil    = soil%dsl_dmax * dry / forcing%theta_span / forcing%dsl_diffusivity
      dr_dtheta = -soil%dsl_dmax * ddry / forcing%theta_span / forcing%dsl_diffusivity
      resist    = forcing%r_aero + r_soil
      e_soil    = forcing%snow_free_frac * forcing%rho_air * drive / resist
      de_dpsi   = forcing%snow_free_frac * forcing%rho_air * ddrive * forcing%q_sat * dalpha / resist
      de_dtheta = -e_soil * dr_dtheta / resist
   end subroutine ground_evaporation_rate

end module meds_soil_water
