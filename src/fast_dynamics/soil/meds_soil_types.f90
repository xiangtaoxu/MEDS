! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_types -- the argument records (data only) of the ground-column kernels: the    !
! boundary ground evaporation reads, the soil thermal column's forcing and outputs, and    !
! the snow store's (docs/science/soil_biophysics.md, docs/science/snow_biophysics.md).     !
! The prognostic stores are in meds_patch_state_vars; the soil's parameters and options   !
! in meds_soil_params, the snow's parameters in meds_biophysics_params.                    !
!==========================================================================================!
module meds_soil_types
   use meds_kinds,             only : wp, ik
   use meds_constants,      only : n_soil_layer_max
   implicit none
   private

   public :: chydro_forcing_t
   public :: energy_forcing_t, energy_flux_t
   public :: snow_env_t, snow_flux_t, snow_melt_t

   !=======================================================================================!
   !  Soil-column water: the boundary ground evaporation reads. The soil's profile and     !
   !  parameters are soil_profile_t and soil_params_t (meds_soil_params).                  !
   !=======================================================================================!

   !----- The surface boundary ground_evaporation_rate reads (read-only). ---------------------------!
   type :: chydro_forcing_t
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temperature
      !----- A frozen surface evaporates over ice (docs/science/snow_biophysics.md eq. (1)). ---------!
      real(wp) :: ground_fliq = 1.0_wp                    !< [-] top-layer liquid fraction
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance of the bare-soil tile
                                                          !<       (1/ggnet), not inflated for snow
      real(wp) :: snow_free_frac = 1.0_wp                 !< [-] snow-free area fraction (1-snowfac) that evaporates,
                                                          !<     an area-weighted tile flux (snow_biophysics.md)
   end type chydro_forcing_t

   !----- Soil-column thermal boundary conditions (read-only). ------------------------------!
   type :: energy_forcing_t
      real(wp) :: g_top      = 0.0_wp                       !< [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !< [W/m2] prescribed bottom flux, positive up (held at 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !< [m3/m3] theta from hydrology (kappa, C_eff)
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !< [m/s]  inter-layer water flux, upward-positive like hf
                                                            !<        (assemble_soil_energy_forcing flips the network's)
      real(wp) :: root_heat_sink(n_soil_layer_max) = 0.0_wp !< [W/m2] enthalpy leaving each layer with water other
                                                            !<        than across a face (root uptake, drainage)
      !----- The top face's water, applied by the same upwind rule and time level as the interior      !
      !      faces (docs/science/soil_biophysics.md, "Water-thermal coupling"); no bottom face term. ---!
      real(wp) :: w_flux_top = 0.0_wp                       !< [m/s] water across the top face (<0 = infiltration in)
      real(wp) :: t_water_top = 0.0_wp                      !< [K] temperature of the infiltrating pond water
   end type energy_forcing_t

   !----- Soil-column energy outputs + diagnostics. ----------------------------------------!
   type :: energy_flux_t
      real(wp) :: bottom_heat  = 0.0_wp                     !< [W/m2] heat conducted out through the bottom face
      real(wp) :: energy_resid = 0.0_wp                     !< [J/m2] closed-budget residual (~0)
      real(wp) :: heat_capacity = 0.0_wp                    !< [J/m2/K] the column's, at the step's start
   end type energy_flux_t

   !=======================================================================================!
   !  Snow / temporary-surface-water types (meds_ground_biophysics;                        !
   !  docs/science/snow_biophysics.md); snow_params_t lives in meds_biophysics_params.     !
   !=======================================================================================!

   !----- Surface boundary conditions for the snow-surface energy balance (read-only, forced). --!
   type :: snow_env_t
      real(wp) :: abs_sw   = 0.0_wp, abs_lw = 0.0_wp    !< [W/m2] absorbed SW (snow albedo), NET LW at the snow surface
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp !< [K],[kg/kg] CAS state (forced sibling)
      real(wp) :: ggnet    = 0.0_wp                     !< [m/s] ground<->CAS conductance (heat = vapour)
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp !< [kg/m3],[Pa] CAS air
      real(wp) :: t_soil_top = 273.16_wp                !< [K] top soil-node temperature (snow-base conduction target)
      real(wp) :: k_soil_top = 1.5_wp                   !< [W/m/K] top soil thermal conductivity (k_face geom-mean)
      real(wp) :: dz_soil_top = 0.05_wp                 !< [m] top soil node depth |z_node(1)| (interface spacing)
   end type snow_env_t

   !----- Snow-surface energy-balance outputs (all positive upward; per unit ground area). ------!
   type :: snow_flux_t
      real(wp) :: t_surf  = 273.16_wp, fliq = 0.0_wp    !< [K],[-] diagnosed snow-surface temperature + liquid fraction
      real(wp) :: h_snow  = 0.0_wp, le_snow = 0.0_wp    !< [W/m2] sensible + latent (sublimation/evap) to the CAS
      real(wp) :: w_flux  = 0.0_wp                      !< [kg/m2/s] vapour mass to CAS (>0 sublimation, <0 deposition)
      real(wp) :: g_base  = 0.0_wp                      !< [W/m2] conduction snow-base -> soil top (the new soil top BC)
      real(wp) :: rnet    = 0.0_wp                      !< [W/m2] net all-wave radiation absorbed by the snow surface
      real(wp) :: snowfac = 0.0_wp                      !< [-] Niu-Yang07 snow-cover / burial fraction
      real(wp) :: energy_resid = 0.0_wp                 !< [J/m2] closed-budget residual (~0 by construction)
   end type snow_flux_t

   !----- Snow-mass outputs from accumulation + melt drainage (a store->store handoff to the pond). -!
   type :: snow_melt_t
      real(wp) :: melt_mass = 0.0_wp                    !< [kg/m2] free liquid drained this step -> the pond
      real(wp) :: melt_enth = 0.0_wp                    !< [J/m2]  matching enthalpy (internal_energy_liquid(t_3ple))
      real(wp) :: dump_mass = 0.0_wp                    !< [kg/m2] full-melt-out residual water dumped to the pond
      real(wp) :: dump_enth = 0.0_wp                    !< [J/m2]  full-melt-out residual energy dumped to the pond
      logical  :: melted_out = .false.                  !< .true. when the layer vanished (nlayer -> 0)
      real(wp) :: mass_resid = 0.0_wp                   !< [kg/m2] closed-budget residual (~0 by construction)
   end type snow_melt_t

end module meds_soil_types
