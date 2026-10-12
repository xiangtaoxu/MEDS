! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_water -- the soil's water BOUNDARY processes, in the form the soil-plant water       !
! network takes them (MEDS_SPAC_WATER_NETWORK_PLAN.md §4.2): functions of the top and bottom     !
! layers' state that return their flux and its derivatives.                                     !
!   * ground_evaporation_rate -- evaporation from the top layer: Philip's pore-space humidity and !
!                                the Swenson-Lawrence dry-surface-layer resistance;              !
!   * bottom_flux             -- drainage through the bottom face: free drainage, a sealed       !
!                                bedrock, or a head-driven aquifer.                               !
! The column's water itself -- storage, Darcy flow between layers and from the pond, root uptake  !
! -- is solved with the plants' in meds_water_network. Canopy interception lives in                !
! meds_plant_biophysics (a per-cohort film).                                                       !
!==========================================================================================!
module meds_soil_water
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav, r_wv, tiny_num, p_std, safe_exp
   use meds_soil_types, only : chydro_forcing_t
   use meds_column_params, only : soil_params_t, curve_a, curve_n
   use meds_biophysics_opts, only : soil_opts_t, SOIL_BC_BEDROCK, SOIL_BC_AQUIFER
   use meds_water_retention,  only : soil_theta_from_psi
   use meds_water_paths,      only : smooth_max0
   use meds_therm_lib,        only : sat_specific_humidity
   implicit none
   private

   public :: ground_evaporation_rate, bottom_flux

   !----- Widths of the two rounded corners in ground_evaporation_rate: [kg/kg] where the surface's  !
   !      humidity meets the air's (evaporation stops; there is no dew), and [m3/m3] where the top   !
   !      layer dries below the moisture at which the dry surface layer starts to form. -----------!
   real(wp), parameter :: EVAP_HUMIDITY_WIDTH = 1.0e-6_wp
   real(wp), parameter :: DSL_THETA_WIDTH     = 1.0e-3_wp
   !----- CLM5 eq 5.78: the 'air dry' matric head [m], where the soil counts as dry (CLM5's 1e7 mm). !
   real(wp), parameter :: PSI_AIR_DRY = -1.0e4_wp

contains

   !---------------------------------------------------------------------------------------!
   ! BOTTOM-FACE Darcy flux [m/s, DOWN-positive]. One authority for all three bottom BCs.    !
   !                                                                                        !
   !   free_drain : unit gradient, q = K(theta_n).  The DEEP-water-table limit.              !
   !   bedrock    : sealed, q = 0.                                                          !
   !   aquifer    : head-driven against a SATURATED zone at the column base -- psi_wt = 0 at !
   !                z_wt = z_bottom, so the soil column IS the unsaturated zone above the    !
   !                water table (MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md Phase 0):            !
   !                    q = K_bot * (psi_n / Delta + 1),   Delta = z_n - z_bottom = dz(n)/2  !
   !                It reduces to free drainage when psi_n -> 0, and REVERSES (upward,       !
   !                capillary rise) once |psi_n| > Delta.  This is a pure BOUNDARY CONDITION:!
   !                there is no aquifer store, no baseflow bucket and no water-table state.  !
   !                                                                                        !
   ! K_bot is UPSTREAM-weighted on the head gradient (which is K-independent, so there is no !
   ! circularity): K(theta_n) for downward flow, K_sat for upward flow OUT of the saturated  !
   ! zone.  Using K(theta_n) upward would silently under-predict capillary rise, which is    !
   ! the entire point of the boundary.                                                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine bottom_flux(params, opts, n, psi_n, k_n, qbot, dq_dpsi, dq_dk)
      type(soil_params_t), intent(in)  :: params
      type(soil_opts_t),   intent(in)  :: opts
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
      select case (opts%bottom_bc)
      case (SOIL_BC_BEDROCK)
         qbot = 0.0_wp
      case (SOIL_BC_AQUIFER)
         delta = 0.5_wp * params%dz(n)                                  ! z_node(n) - z_bottom > 0
         grad  = psi_n / delta + 1.0_wp
         if (grad >= 0.0_wp) then
            k_bot = k_n                                                 ! downward: upstream is layer n
            if (present(dq_dk)) dq_dk = grad
         else
            k_bot = params%ksat(n)                                      ! upward: upstream is saturated
         end if
         qbot = k_bot * grad
         if (present(dq_dpsi)) dq_dpsi = k_bot / delta
      case default                                                      ! SOIL_BC_FREE_DRAIN
         qbot = k_n
         if (present(dq_dk)) dq_dk = 1.0_wp
      end select
   end subroutine bottom_flux

   !---------------------------------------------------------------------------------------!
   ! ground_evaporation_rate -- ground evaporation [kg/m2/s] from the top layer at water content   !
   ! theta1 [m3/m3] and matric head psi1 [m], with its derivatives, for the soil-plant water         !
   ! network, which takes it at the end-of-step state.                                                !
   !   * The surface's humidity is Philip's pore-space humidity, alpha = exp(psi g / (R_v T)), times  !
   !     saturation over ice when the top layer is frozen (#89: the liquid curve overstates e_sat by   !
   !     10 % at -10 C and 34 % at -30 C, and a frozen surface sublimes). As the layer dries its        !
   !     potential falls without bound and alpha with it, so evaporation stops by itself before the    !
   !     layer empties; there is no cap at the water above residual.                                     !
   !   * The dry-surface-layer resistance (Swenson and Lawrence; CLM5 eq 5.77-5.80) grows as the      !
   !     layer dries below theta_init. Its thickness and tortuosity are referenced to the 'air dry'    !
   !     water content theta_air, a texture constant at PSI_AIR_DRY: referenced to the bulk theta1     !
   !     instead, tortuosity falls faster than the layer thickens and a wetter soil evaporates LESS.   !
   !     Tortuosity is Millington-Quirk rather than CLM5's phi_air^2 (phi_air/phi)^(3/b): b is a        !
   !     Clapp-Hornberger parameter and the shipped curve is van Genuchten; with phi_air constant the  !
   !     two differ only by a constant factor, so the shape of the moisture response is unaffected.    !
   !   * Only the snow-free share of the ground evaporates, at the bare-soil resistance              !
   !     (chydro_forcing_t%snow_free_frac): an area-weighted tile flux, not an inflated r_aero.         !
   !   * The two corners -- evaporation stops where the surface's humidity meets the air's (there is   !
   !     no dew), and the dry surface layer starts below theta_init -- are rounded over                 !
   !     EVAP_HUMIDITY_WIDTH and DSL_THETA_WIDTH, so Newton sees them.                                   !
   !---------------------------------------------------------------------------------------!
   pure subroutine ground_evaporation_rate(theta1, psi1, params, forcing, opts, e_soil, de_dtheta, de_dpsi)
      real(wp),               intent(in)  :: theta1      !< [m3/m3] top-layer water content
      real(wp),               intent(in)  :: psi1        !< [m]     top-layer matric head
      type(soil_params_t),    intent(in)  :: params
      type(chydro_forcing_t), intent(in)  :: forcing
      type(soil_opts_t),      intent(in)  :: opts
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
      theta_air  = soil_theta_from_psi_l(params%retention, params, 1_ik, PSI_AIR_DRY)
      theta_init = opts%dsl_theta_init * params%theta_sat(1)
      call smooth_max0(theta_init - theta1, DSL_THETA_WIDTH, dry, ddry)
      dvap    = 2.12e-5_wp * (forcing%t_ground / 273.15_wp) ** 1.75_wp
      phi     = params%theta_sat(1)
      phi_air = max(phi - theta_air, tiny_num)
      tau     = phi_air ** (10.0_wp / 3.0_wp) / max(phi * phi, tiny_num)
      r_soil    = opts%dsl_dmax * dry / max(theta_init - theta_air, tiny_num) / max(dvap * tau, tiny_num)
      dr_dtheta = -opts%dsl_dmax * ddry / max(theta_init - theta_air, tiny_num) / max(dvap * tau, tiny_num)
      resist    = forcing%r_aero + r_soil
      e_soil    = forcing%snow_free_frac * forcing%rho_air * drive / resist
      de_dpsi   = forcing%snow_free_frac * forcing%rho_air * ddrive * q_sat * dalpha / resist
      de_dtheta = -e_soil * dr_dtheta / resist
   end subroutine ground_evaporation_rate

   !----- theta(psi) with the layer-k curve params (thin wrapper). ------------------------!
   pure function soil_theta_from_psi_l(rc, params, k, psi) result(theta)
      integer(ik),         intent(in) :: rc, k
      type(soil_params_t), intent(in) :: params
      real(wp),            intent(in) :: psi
      real(wp)                        :: theta
      theta = soil_theta_from_psi(rc, psi, params%theta_sat(k), params%theta_res(k),             &
                                curve_a(params,k), curve_n(params,k))
   end function soil_theta_from_psi_l

end module meds_soil_water
