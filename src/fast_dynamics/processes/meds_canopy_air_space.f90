! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_canopy_air_space -- the canopy-air-space (CAS) box                                  !
! (docs/science/canopy_air_space_biophysics.md): one well-mixed air layer with three       !
! prognostic twins (specific enthalpy, specific humidity, CO2 mixing ratio) sharing one    !
! capacity and the atm<->CAS conductance. cas_box_time_deriv is the explicit right-hand !
! side of eq. (1), which rk45 integrates; cas_box_step_implicit is the backward-Euler   !
! update of eq. (2), the box commit of each ESDIRK2 stage. The integrator assembles the    !
! sources (cas_env_t) and the capacities, conductances and atmospheric state            !
! (cas_env_t).                                                                             !
!==========================================================================================!
module meds_canopy_air_space
   use meds_kinds,            only : wp, ik
   implicit none
   private

   public :: cas_env_t
   public :: cas_box_time_deriv, cas_box_step_implicit

   !----- What the canopy-air box reads, per patch: its capacities, the atm<->CAS conductances and the   !
   !      atmosphere's state, held for the step, and the summed sources into it, set at each         !
   !      evaluation. Any adjustment (the supersaturation sink) is folded into the sources before the  !
   !      box update, so both kernels share it. ---------------------------------------------------!
   type :: cas_env_t
      real(wp) :: air_mass_capacity  = 0.0_wp   !< [kg/m2]  CAS air mass per ground area (enthalpy+vapour twins)
      real(wp) :: air_molar_capacity = 0.0_wp   !< [mol/m2] CAS dry-air molar capacity (CO2 twin)
      real(wp) :: atm_conductance_enthalpy = 0.0_wp  !< [kg/m2/s]  atm<->CAS enthalpy conductance (rho*ustar*heat_transfer_factor)
      real(wp) :: atm_conductance_vapor    = 0.0_wp  !< [kg/m2/s]  atm<->CAS vapour conductance
      real(wp) :: atm_conductance_co2      = 0.0_wp  !< [mol/m2/s] atm<->CAS CO2 conductance
      real(wp) :: atm_enthalpy         = 0.0_wp !< [J/kg]      reference-level specific enthalpy
      real(wp) :: atm_specific_humidity= 0.0_wp !< [kg/kg]     reference-level specific humidity
      real(wp) :: atm_co2              = 400.0_wp!< [umol/mol]  free-atmosphere CO2
      real(wp) :: surface_enthalpy_source = 0.0_wp   !< [W/m2]      sensible + latent into the CAS
      real(wp) :: surface_vapor_source    = 0.0_wp   !< [kg/m2/s]   vapour into the CAS
      real(wp) :: biotic_co2_source       = 0.0_wp   !< [umol/m2/s] respiration - GPP (Reco - GPP)
   end type cas_env_t

contains

   !----- The explicit CAS-box right-hand side d(twin)/dt at the current state                      !
   !      (docs/science/canopy_air_space_biophysics.md eq. (1)), which rk45 integrates. -----------!
   pure subroutine cas_box_time_deriv(cas_enthalpy, cas_shv, cas_co2, env,                &
                                         d_enthalpy, d_shv, d_co2)
      real(wp),           intent(in)  :: cas_enthalpy, cas_shv, cas_co2
      type(cas_env_t),    intent(in)  :: env
      real(wp),           intent(out) :: d_enthalpy, d_shv, d_co2
      d_enthalpy = (env%surface_enthalpy_source                                             &
                    + env%atm_conductance_enthalpy * (env%atm_enthalpy - cas_enthalpy))  &
                   / env%air_mass_capacity
      d_shv      = (env%surface_vapor_source                                                &
                    + env%atm_conductance_vapor * (env%atm_specific_humidity - cas_shv))  &
                   / env%air_mass_capacity
      d_co2      = (env%biotic_co2_source                                                    &
                    + env%atm_conductance_co2 * (env%atm_co2 - cas_co2))                  &
                   / env%air_molar_capacity
   end subroutine cas_box_time_deriv

   !----- One backward-Euler CAS-box step over dt, implicit in the atmospheric exchange (L-stable;  !
   !      docs/science/canopy_air_space_biophysics.md eq. (2)): the box commit of an ESDIRK2        !
   !      stage. -----------------------------------------------------------------------------------!
   pure subroutine cas_box_step_implicit(cas_enthalpy, cas_shv, cas_co2, env, dt,         &
                                            enthalpy_new, shv_new, co2_new)
      real(wp),           intent(in)  :: cas_enthalpy, cas_shv, cas_co2
      type(cas_env_t),    intent(in)  :: env
      real(wp),           intent(in)  :: dt
      real(wp),           intent(out) :: enthalpy_new, shv_new, co2_new
      enthalpy_new = (env%air_mass_capacity * cas_enthalpy                                  &
                      + dt * (env%surface_enthalpy_source                                   &
                              + env%atm_conductance_enthalpy * env%atm_enthalpy))        &
                     / (env%air_mass_capacity + dt * env%atm_conductance_enthalpy)
      shv_new      = (env%air_mass_capacity * cas_shv                                       &
                      + dt * (env%surface_vapor_source                                      &
                              + env%atm_conductance_vapor * env%atm_specific_humidity))  &
                     / (env%air_mass_capacity + dt * env%atm_conductance_vapor)
      co2_new      = (env%air_molar_capacity * cas_co2                                      &
                      + dt * (env%biotic_co2_source                                         &
                              + env%atm_conductance_co2 * env%atm_co2))                  &
                     / (env%air_molar_capacity + dt * env%atm_conductance_co2)
   end subroutine cas_box_step_implicit

end module meds_canopy_air_space
