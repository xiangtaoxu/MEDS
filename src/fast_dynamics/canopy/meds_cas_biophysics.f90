! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_cas_biophysics -- the CANOPY-AIR-SPACE (CAS) balance                                !
! (docs/science/canopy_air_space_biophysics.md): the well-mixed sub-canopy air box and its !
! three prognostic twins (specific enthalpy, specific humidity, CO2 mixing ratio), all     !
! sharing one capacity and the ustar-based atm<->CAS conductance. A fast diffusion/venting !
! exchange, so it lives with the biophysics kernels (not the slow soil-carbon pools).      !
!                                                                                          !
! The box is advanced by the two-form kernel pair the two integrators share:               !
! cas_column_time_deriv (the explicit RHS of eq. (1), which the rk45 stepper integrates)   !
! and cas_column_step_implicit (the backward-Euler update of eq. (2), the box commit of    !
! each ark stage). The driver assembles the summed surface + biotic sources (cas_source_t) !
! and the capacities / conductances / atmospheric BCs (cas_column_t). (Heterotrophic soil  !
! respiration -- a carbon-decomposition process -- lives in meds_soil_biogeochem.)         !
!==========================================================================================!
module meds_cas_biophysics
   use meds_kinds,            only : wp, ik
   implicit none
   private

   public :: cas_column_t, cas_source_t
   public :: cas_column_time_deriv, cas_column_step_implicit

   !----- Frozen canopy-air-space column params (per patch, per substep). A single well-mixed     !
   !      layer (see docs/ROADMAP.md), shaped so a MULTI-LAYER canopy-air column (K-theory        !
   !      Thomas, n>1) is the general case and this is the n=1 degenerate one                     !
   !      (MEDS_COLUMN_CO2_BALANCE_DESIGN.md §3.5).                                               !
   type :: cas_column_t
      real(wp) :: air_mass_capacity  = 0.0_wp   !< [kg/m2]  CAS air mass per ground area (enthalpy+vapour twins)
      real(wp) :: air_molar_capacity = 0.0_wp   !< [mol/m2] CAS dry-air molar capacity (CO2 twin)
      real(wp) :: atm_conductance_enthalpy = 0.0_wp  !< [kg/m2/s]  atm<->CAS enthalpy conductance (rho*ustar*temp1)
      real(wp) :: atm_conductance_vapor    = 0.0_wp  !< [kg/m2/s]  atm<->CAS vapour conductance
      real(wp) :: atm_conductance_co2      = 0.0_wp  !< [mol/m2/s] atm<->CAS CO2 conductance
      real(wp) :: atm_enthalpy         = 0.0_wp !< [J/kg]      reference-level specific enthalpy
      real(wp) :: atm_specific_humidity= 0.0_wp !< [kg/kg]     reference-level specific humidity
      real(wp) :: atm_co2              = 400.0_wp!< [umol/mol]  free-atmosphere CO2
   end type cas_column_t

   !----- Summed surface sources into the CAS, assembled by the driver from the surface fluxes.    !
   !      Any source adjustment (the CAS supersaturation sink, `[fast].cas_condensation`) is       !
   !      folded in BEFORE the box update, so the box kernels below are scheme-shared.             !
   type :: cas_source_t
      real(wp) :: surface_enthalpy_source = 0.0_wp   !< [W/m2]      sensible + latent into the CAS
      real(wp) :: surface_vapor_source    = 0.0_wp   !< [kg/m2/s]   vapour into the CAS
      real(wp) :: biotic_co2_source       = 0.0_wp   !< [umol/m2/s] respiration - GPP (Reco - GPP)
   end type cas_source_t

contains

   !---------------------------------------------------------------------------------------!
   ! cas_column_time_deriv -- the EXPLICIT CAS-box RHS d(twin)/dt at the current state,    !
   ! docs/science/canopy_air_space_biophysics.md eq. (1): the right-hand side the rk45     !
   ! stepper integrates (surface_derivs). Pure, commits nothing. The surface + biotic      !
   ! sources arrive already assembled.                                                     !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_column_time_deriv(cas_enthalpy, cas_shv, cas_co2, source, column,       &
                                         d_enthalpy, d_shv, d_co2)
      real(wp),           intent(in)  :: cas_enthalpy, cas_shv, cas_co2
      type(cas_source_t), intent(in)  :: source
      type(cas_column_t), intent(in)  :: column
      real(wp),           intent(out) :: d_enthalpy, d_shv, d_co2
      d_enthalpy = (source%surface_enthalpy_source                                             &
                    + column%atm_conductance_enthalpy * (column%atm_enthalpy - cas_enthalpy))  &
                   / column%air_mass_capacity
      d_shv      = (source%surface_vapor_source                                                &
                    + column%atm_conductance_vapor * (column%atm_specific_humidity - cas_shv))  &
                   / column%air_mass_capacity
      d_co2      = (source%biotic_co2_source                                                    &
                    + column%atm_conductance_co2 * (column%atm_co2 - cas_co2))                  &
                   / column%air_molar_capacity
   end subroutine cas_column_time_deriv

   !---------------------------------------------------------------------------------------!
   ! cas_column_step_implicit -- one backward-Euler CAS-box advance over dt, implicit in   !
   ! the atm exchange (L-stable), docs/science/canopy_air_space_biophysics.md eq. (2): the !
   ! box commit of each ark stage (meds_fast_be_stage). The surface + biotic sources       !
   ! arrive already assembled. Returns the updated twins; nothing else committed.          !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_column_step_implicit(cas_enthalpy, cas_shv, cas_co2, source, column, dt, &
                                            enthalpy_new, shv_new, co2_new)
      real(wp),           intent(in)  :: cas_enthalpy, cas_shv, cas_co2
      type(cas_source_t), intent(in)  :: source
      type(cas_column_t), intent(in)  :: column
      real(wp),           intent(in)  :: dt
      real(wp),           intent(out) :: enthalpy_new, shv_new, co2_new
      enthalpy_new = (column%air_mass_capacity * cas_enthalpy                                  &
                      + dt * (source%surface_enthalpy_source                                   &
                              + column%atm_conductance_enthalpy * column%atm_enthalpy))        &
                     / (column%air_mass_capacity + dt * column%atm_conductance_enthalpy)
      shv_new      = (column%air_mass_capacity * cas_shv                                       &
                      + dt * (source%surface_vapor_source                                      &
                              + column%atm_conductance_vapor * column%atm_specific_humidity))  &
                     / (column%air_mass_capacity + dt * column%atm_conductance_vapor)
      co2_new      = (column%air_molar_capacity * cas_co2                                      &
                      + dt * (source%biotic_co2_source                                         &
                              + column%atm_conductance_co2 * column%atm_co2))                  &
                     / (column%air_molar_capacity + dt * column%atm_conductance_co2)
   end subroutine cas_column_step_implicit

end module meds_cas_biophysics
