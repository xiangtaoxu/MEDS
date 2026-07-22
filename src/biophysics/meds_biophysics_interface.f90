!==========================================================================================!
! meds_biophysics_interface -- the public façade of the fast-loop biophysics library (the         !
! biophysics analogue of meds_plant_interface / meds_core_interface): it RE-EXPORTS the surface-   !
! subsystem seams + their shared types so a black-box caller can `use` this one module.            !
!                                                                                          !
! ORCHESTRATION is NOT here -- it lives in the drivers (meds_column_dynamics operator-split,        !
! meds_ark_stepper / meds_fast_time_derivs IMEX-ARK), which sequence these kernels + the per-patch    !
! state. This is a CONVENIENCE façade, not a sealed wall; two access patterns coexist: BLACK-BOX      !
! callers (e.g. meds_fast_dynamics) `use` this module for the seams, while the WHITE-BOX fast-loop     !
! integrators `use` the kernel modules directly (they need the `*_time_deriv` RHS at each ARK stage,   !
! which a per-call seam cannot expose). The shared optical-property + soil constitutive kernels stay   !
! in shared/functions (meds_optics_lib / meds_hydr_lib / meds_therm_lib), imported where needed.        !
!==========================================================================================!
module meds_biophysics_interface
   !----- Canopy radiative transfer (optics assembly + two-stream solver + seam). ----------!
   use meds_canopy_radiation, only : canopy_radiation, derive_rad_optics, blend_cohort_optics,  &
                                     ground_optics, solve_band,                                 &
                                     rad_pft_optics_t, rad_forcing_t, rad_flux_t, surface_state_t
   !----- Canopy aerodynamics. -------------------------------------------------------------!
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   !----- Soil thermal column (two forms). -------------------------------------------------!
   use meds_soil_energy,      only : soil_energy_step_implicit, soil_energy_time_deriv
   !----- Soil water column (seam + two forms). --------------------------------------------!
   use meds_soil_water,       only : column_hydrology_flux, soil_water_step_implicit,           &
                                     soil_water_time_deriv
   !----- Vegetation surface (leaf/wood energy + interception). ----------------------------!
   use meds_vegetation_biophysics, only : veg_energy_step_implicit, intercept_canopy_layer
   !----- Ground surface (bare skin + snow store). -----------------------------------------!
   use meds_ground_biophysics, only : ground_surface_fluxes, snow_cover_fraction,              &
                                     snow_accumulate, snow_drain_meltwater, snow_surface_fluxes, &
                                     snow_base_conductance, snow_energy_step
   !----- Canopy air space (enthalpy/humidity/CO2 twins) -- the two-form box kernels. ---------!
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t,                               &
                                     cas_column_time_deriv, cas_column_step_implicit
   implicit none
   private

   !----- Re-exported public surface (the seams + their re-exported types). ----------------!
   public :: canopy_radiation, derive_rad_optics, blend_cohort_optics, ground_optics, solve_band
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t, surface_state_t
   public :: canopy_aerodynamics
   public :: soil_energy_step_implicit, soil_energy_time_deriv
   public :: column_hydrology_flux, soil_water_step_implicit, soil_water_time_deriv
   public :: veg_energy_step_implicit, intercept_canopy_layer
   public :: ground_surface_fluxes, snow_cover_fraction, snow_accumulate, snow_drain_meltwater
   public :: snow_surface_fluxes, snow_base_conductance, snow_energy_step
   public :: cas_column_t, cas_source_t, cas_column_time_deriv, cas_column_step_implicit

end module meds_biophysics_interface
