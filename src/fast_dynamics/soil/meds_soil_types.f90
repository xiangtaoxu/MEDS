! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_types -- the argument records of the ground-column kernels: soil hydrology,    !
! soil thermal energy, and the snow store (docs/science/soil_biophysics.md,                !
! docs/science/snow_biophysics.md). Pure DATA, no methods, no hidden state.                !
!                                                                                          !
! Each domain follows the same shape: a *_forcing_t / *_env_t the kernel reads, and a      !
! *_flux_t it returns, with the PROGNOSTIC stores themselves in meds_patch_state_vars and  !
! the run-config selectors in meds_biophysics_params.                                      !
!==========================================================================================!
module meds_soil_types
   use meds_kinds,             only : wp, ik
   use meds_constants,      only : n_soil_layer_max, n_snow_layer_max
   implicit none
   private

   public :: chydro_forcing_t
   public :: energy_forcing_t, energy_flux_t
   public :: snow_env_t, snow_flux_t, snow_melt_t

   !=======================================================================================!
   !  Soil-column water: the boundary the water network's evaporation reads. The soil's       !
   !  parameters are soil_params_t (meds_soil_params), its options soil_opts_t                !
   !  (meds_biophysics_params), its retention curves meds_soil_hydr_lib.                       !
   !=======================================================================================!

   !----- The soil column's surface boundary for the water network's evaporation (read-only). -----!
   type :: chydro_forcing_t
      real(wp) :: precip_ground = 0.0_wp                  !< [kg/m2/s] water reaching the ground (post interception)
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temperature
      !----- Liquid fraction of the top soil layer, so ground evaporation from a FROZEN surface       !
      !      saturates over ICE rather than over liquid (docs/science/snow_biophysics.md eq. (1)).    !
      !      Defaults to 1 (all liquid). -------------------------------------------------------------!
      real(wp) :: ground_fliq = 1.0_wp                    !< [-] top-layer liquid fraction
      !----- The temperature that values precip_ground as it joins the pond. Defaulted to a PHYSICAL   !
      !      temperature, not 0: internal_energy_liquid is referenced to tsupercool_liq (~57 K), so a 0 K !
      !      default would value the pond's inflow at a large NEGATIVE enthalpy. ----------------------!
      real(wp) :: t_pond_inflow = 298.15_wp               !< [K]
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance of the BARE-SOIL tile
                                                          !<       (1/ggnet) -- NOT snow-inflated; see below.
      real(wp) :: snow_free_frac = 1.0_wp                 !< [-] snow-free AREA fraction (1-snowfac). Ground
                                                          !<     evaporation is an AREA-weighted tile flux:
                                                          !<     E = snow_free_frac * rho*(q_g-q_air)/(r_aero+r_soil).
                                                          !<     Throttling partial cover by inflating r_aero
                                                          !<     instead is NOT equivalent whenever the dry-surface
                                                          !<     -layer resistance r_soil > 0 -- it divides only the
                                                          !<     aerodynamic leg by (1-snowfac), so it over-predicts
                                                          !<     bare-soil evaporation at intermediate cover, worst
                                                          !<     where a dry surface makes r_soil dominant. The two
                                                          !<     agree exactly at snowfac = 0 and snowfac = 1.
   end type chydro_forcing_t

   !----- Soil-column thermal boundary conditions (read-only). ------------------------------!
   type :: energy_forcing_t
      real(wp) :: g_top      = 0.0_wp                       !< [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !< [W/m2] bottom flux (default 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !< [m3/m3] theta from hydrology (kappa, C_eff)
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !< [m/s]  inter-layer water flux for advective heat,
                                                            !<        UPWARD-positive (matches the hf face convention;
                                                            !<        the caller flips the DOWNWARD-positive hydrology
                                                            !<        flux -- see assemble_soil_energy_forcing:
                                                            !<        eforc%w_flux = -face_flux). Populating it with
                                                            !<        the raw downward flux would reverse the advection.
      real(wp) :: root_heat_sink(n_soil_layer_max) = 0.0_wp !< [W/m2] enthalpy removed with root uptake
      !----- BOUNDARY water fluxes, same UPWARD-positive convention as w_flux above, with the          !
      !      temperature of water arriving from OUTSIDE the column (docs/science/soil_biophysics.md    !
      !      eq. (10)). These close the water-enthalpy advection at the boundary faces so ALL of it -- !
      !      top, interior, bottom -- is applied by ONE upwind rule at ONE time level. Adding the      !
      !      top/bottom terms to soil_energy directly before the step would evaluate the inflow on     !
      !      state^n while the interior outflow uses the post-conduction T^{n+1}. With                 !
      !      internal_energy_liquid ~1.0 MJ/kg (the tsupercool_liq datum), each face term reaches      !
      !      ~1300 W/m2 for a few mm/h of percolation, so the physical signal is the small DIFFERENCE  !
      !      of two very large numbers -- and a 5 K time-level inconsistency in one of them is worth   !
      !      ~30 W/m2, i.e. the entire signal. --------------------------------------------------------!
      real(wp) :: w_flux_top = 0.0_wp                       !< [m/s] water across the TOP face (<0 = infiltration in)
      real(wp) :: w_flux_bot = 0.0_wp                       !< [m/s] water across the BOTTOM face (<0 = drainage out)
      real(wp) :: t_water_top = 0.0_wp                      !< [K] temperature of water entering from above
      real(wp) :: t_water_bot = 0.0_wp                      !< [K] temperature of water entering from below
   end type energy_forcing_t

   !----- Soil-column energy outputs + diagnostics. ----------------------------------------!
   type :: energy_flux_t
      real(wp) :: ground_heat  = 0.0_wp                     !< [W/m2] conductive flux into layer 1
      real(wp) :: bottom_heat  = 0.0_wp                     !< [W/m2] advective+geothermal bottom loss
      real(wp) :: energy_resid = 0.0_wp                     !< [J/m2] closed-budget residual (~0)
      integer(ik) :: nsub = 0_ik
      logical  :: converged = .true.
   end type energy_flux_t

   !=======================================================================================!
   !  Canopy-air-space CO2 balance: the prognostic third twin is carried in cas_state_t    !
   !  and advanced by meds_cas_biophysics's cas_column_* kernels. Heterotrophic            !
   !  respiration is decomposition and lives in meds_soil_biogeochem, with ONE production  !
   !  authority: the CENTURY matrix.                                                       !
   !=======================================================================================!

   !=======================================================================================!
   !  SNOW / temporary-surface-water types (meds_ground_biophysics;                        !
   !  docs/science/snow_biophysics.md). STATELESS per-store kernels; forcing arrives as    !
   !  value types. (snow_params_t -- the physical parameter table -- lives in              !
   !  meds_biophysics_params.)                                                             !
   !=======================================================================================!

   !----- Surface boundary conditions for the snow-surface energy balance (read-only, FORCED). --!
   type :: snow_env_t
      real(wp) :: abs_sw   = 0.0_wp, abs_lw = 0.0_wp    !< [W/m2] absorbed SW (snow albedo), NET LW at the snow surface
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp !< [K],[kg/kg] CAS state (forced sibling)
      real(wp) :: ggnet    = 0.0_wp                     !< [m/s] ground<->CAS conductance (heat = vapour)
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp !< [kg/m3],[Pa] CAS air
      real(wp) :: t_soil_top = 273.16_wp                !< [K] top soil-node temperature (snow-base conduction target)
      real(wp) :: k_soil_top = 1.5_wp                   !< [W/m/K] top soil thermal conductivity (k_face geom-mean)
      real(wp) :: dz_soil_top = 0.05_wp                 !< [m] top soil node depth |z_node(1)| (interface spacing)
   end type snow_env_t

   !----- Snow-surface energy-balance outputs (all positive UPWARD; per unit ground area). ------!
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
