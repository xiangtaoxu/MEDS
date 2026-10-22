! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_water -- how water moves through and out of the soil column, each flux with    !
! its derivatives for meds_spac_solver (docs/science/soil_biophysics.md eqs. (1)-(4)):     !
!   * soil_water_flux         -- Darcy flow between two layers, or the pond and the top    !
!                                layer, eqs. (3)-(4);                                      !
!   * ground_evaporation_rate -- evaporation from the top layer;                           !
!   * bottom_flux             -- drainage through the bottom face, or a head-driven aquifer.!
! A layer's water relations are in meds_soil_hydr_lib, root uptake in meds_plant_hydr_lib. !
!==========================================================================================!
module meds_soil_water
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav, grav_head, r_wv, tiny_num, p_std, safe_exp
   use meds_soil_types, only : chydro_forcing_t
   use meds_soil_params,      only : soil_profile_t, soil_params_t, SOIL_BC_AQUIFER
   use meds_soil_hydr_lib,    only : soil_head_to_theta
   use meds_numerics,         only : smooth_max0
   use meds_therm_lib,        only : sat_specific_humidity
   implicit none
   private

   public :: soil_water_flux, ground_evaporation_rate, bottom_flux

   !----- Widths of the two rounded corners in ground_evaporation_rate: [kg/kg] where the surface's  !
   !      humidity meets the air's (evaporation stops; there is no dew), and [m3/m3] where the top   !
   !      layer dries below the moisture at which the dry surface layer starts to form. -----------!
   real(wp), parameter :: EVAP_HUMIDITY_WIDTH = 1.0e-6_wp
   real(wp), parameter :: DSL_THETA_WIDTH     = 1.0e-3_wp
   !----- CLM5 eq 5.78: the 'air dry' matric head [m], where the soil counts as dry (CLM5's 1e7 mm). !
   real(wp), parameter :: PSI_AIR_DRY = -1.0e4_wp

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
      type(soil_profile_t), intent(in) :: profile
      type(soil_params_t), intent(in)  :: soil
      integer(ik),         intent(in)  :: n
      real(wp),            intent(in)  :: psi_n     !< [m]   matric head of the deepest layer
      real(wp),            intent(in)  :: k_n       !< [m/s] K(theta_n) of the deepest layer
      real(wp),            intent(out) :: qbot      !< [m/s] down-positive
      real(wp), optional,  intent(out) :: dq_dpsi   !< [1/s] d(qbot)/d(psi_n), for the implicit solve
      real(wp), optional,  intent(out) :: dq_dk     !< [-]   d(qbot)/d(k_n), for the water network
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

   !---------------------------------------------------------------------------------------!
   ! ground_evaporation_rate -- ground evaporation [kg/m2/s] from the top layer at water content        !
   ! theta1 [m3/m3] and matric head psi1 [m], with its derivatives (docs/science/soil_biophysics.md,    !
   ! "The soil water column"): Philip's pore-space humidity (over ice when frozen) through the          !
   ! Swenson-Lawrence dry-surface-layer resistance (CLM5 eqs. 5.77-5.80) in series with r_aero, from    !
   ! the snow-free share only (docs/science/snow_biophysics.md). Both corners are rounded for Newton.   !
   !---------------------------------------------------------------------------------------!
   pure subroutine ground_evaporation_rate(theta1, psi1, profile, forcing, soil, e_soil, de_dtheta, de_dpsi)
      real(wp),               intent(in)  :: theta1      !< [m3/m3] top-layer water content
      real(wp),               intent(in)  :: psi1        !< [m]     top-layer matric head
      type(soil_profile_t),   intent(in)  :: profile
      type(chydro_forcing_t), intent(in)  :: forcing
      type(soil_params_t),    intent(in)  :: soil
      real(wp),               intent(out) :: e_soil      !< [kg/m2/s]
      real(wp),               intent(out) :: de_dtheta   !< [kg/m2/s per m3/m3]
      real(wp),               intent(out) :: de_dpsi     !< [kg/m2/s per m of head]
      real(wp) :: alpha_soil, dalpha, q_sat, drive, ddrive, theta_air, theta_init, dry, ddry
      real(wp) :: dvap, phi, phi_air, tau, r_soil, dr_dtheta, resist
      alpha_soil = safe_exp(psi1 * grav / (r_wv * forcing%t_ground))
      dalpha     = alpha_soil * grav / (r_wv * forcing%t_ground)
      q_sat      = sat_specific_humidity(forcing%t_ground, p_std, forcing%ground_fliq)
      call smooth_max0(alpha_soil * q_sat - forcing%q_air, EVAP_HUMIDITY_WIDTH, drive, ddrive)
      !----- The dry surface layer, referenced to the air-dry water content. ---------------------!
      theta_air  = soil_head_to_theta(profile%retention, PSI_AIR_DRY, profile%theta_sat(1), profile%theta_res(1),  &
                                       profile%curve_par_a(1), profile%curve_par_n(1))
      theta_init = soil%dsl_theta_init * profile%theta_sat(1)
      call smooth_max0(theta_init - theta1, DSL_THETA_WIDTH, dry, ddry)
      dvap    = 2.12e-5_wp * (forcing%t_ground / 273.15_wp) ** 1.75_wp
      phi     = profile%theta_sat(1)
      phi_air = max(phi - theta_air, tiny_num)
      tau     = phi_air ** (10.0_wp / 3.0_wp) / max(phi * phi, tiny_num)
      r_soil    = soil%dsl_dmax * dry / max(theta_init - theta_air, tiny_num) / max(dvap * tau, tiny_num)
      dr_dtheta = -soil%dsl_dmax * ddry / max(theta_init - theta_air, tiny_num) / max(dvap * tau, tiny_num)
      resist    = forcing%r_aero + r_soil
      e_soil    = forcing%snow_free_frac * forcing%rho_air * drive / resist
      de_dpsi   = forcing%snow_free_frac * forcing%rho_air * ddrive * q_sat * dalpha / resist
      de_dtheta = -e_soil * dr_dtheta / resist
   end subroutine ground_evaporation_rate

end module meds_soil_water
