! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_column_config -- the soil column a run is set on: the column keys of [soil]    !
! and the two bundles the soil kernels read, soil_params_t (the layers, the hydraulic      !
! texture and the field capacity) and soil_thermal_params_t (the thermal texture).         !
! meds_config builds them, once per column (soil_column_params). A config leaf (src/share  !
! only): the soil kernels import it without meds_config.                                   !
!==========================================================================================!
module meds_soil_column_config
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : n_soil_layer_max
   use meds_soil_hydr_lib,    only : SOIL_RETENTION_VG
   implicit none
   private

   public :: soil_column_config_t, soil_texture_layers_t, soil_params_t, soil_thermal_params_t

   !----- The column keys of [soil] (docs/science/soil_biophysics.md), each optional at    !
   !      the default below. [soil]'s surface and bottom boundaries are read into          !
   !      soil_opts_t and energy_opts_t (meds_biophysics_opts). ---------------------------!
   type :: soil_column_config_t
      !----- Layers: layer_bottoms or SOIL_LAYER_TEMPLATE, cut at depth. ------------------!
      real(wp)    :: depth           = 2.0_wp    !< [m]  total column depth (positive; z is negative down)
      integer(ik) :: n_layer_bottoms = 0_ik      !< [-]  entries of layer_bottoms; 0 = the template
      real(wp)    :: layer_bottoms(n_soil_layer_max) = 0.0_wp   !< [m] layer bottoms, positive down
      !----- Hydraulic texture, one for the column (or per layer, below). -----------------!
      integer(ik) :: retention   = SOIL_RETENTION_VG  !< van_genuchten | campbell
      real(wp)    :: theta_sat   = 0.43_wp    !< [m3/m3] porosity
      real(wp)    :: theta_res   = 0.078_wp   !< [m3/m3] residual water content
      real(wp)    :: ksat        = 2.89e-6_wp !< [m/s]   saturated hydraulic conductivity
      real(wp)    :: curve_par_a = 3.6_wp     !< [1/m] van Genuchten alpha, or [m] Campbell psi_sat (< 0)
      real(wp)    :: curve_par_n = 1.56_wp    !< [-]   van Genuchten n (>1), or [-] Campbell b
      real(wp)    :: psi_fc      = -3.37_wp   !< [m]   field-capacity matric head (derives theta_fc)
      !----- Thermal texture, one for the column (or per layer, below). -------------------!
      real(wp)    :: solid_conductivity = 3.0_wp    !< [W/m/K]  mineral-solid conductivity
      real(wp)    :: dry_conductivity   = 0.15_wp   !< [W/m/K]  dry-matrix conductivity
      real(wp)    :: dry_heat_capacity  = 2.0e6_wp  !< [J/m3/K] dry-matrix volumetric heat capacity
      !----- Per-layer texture from the soil maps (meds_surface_maps); no config key. -----!
      integer(ik) :: n_layer_values = 0_ik      !< [-] the column's layer count; 0: the one texture above
      real(wp)    :: layer_theta_sat(n_soil_layer_max)          = 0.0_wp
      real(wp)    :: layer_theta_res(n_soil_layer_max)          = 0.0_wp
      real(wp)    :: layer_ksat(n_soil_layer_max)               = 0.0_wp
      real(wp)    :: layer_curve_par_a(n_soil_layer_max)        = 0.0_wp
      real(wp)    :: layer_curve_par_n(n_soil_layer_max)        = 0.0_wp
      real(wp)    :: layer_solid_conductivity(n_soil_layer_max) = 0.0_wp
      real(wp)    :: layer_dry_conductivity(n_soil_layer_max)   = 0.0_wp
      real(wp)    :: layer_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp
   end type soil_column_config_t

   !----- The column's texture, one value per layer: the per-layer arrays, or the one      !
   !      texture repeated (soil_column_texture). -----------------------------------------!
   type :: soil_texture_layers_t
      real(wp) :: theta_sat(n_soil_layer_max), theta_res(n_soil_layer_max), ksat(n_soil_layer_max)
      real(wp) :: par_a(n_soil_layer_max), par_n(n_soil_layer_max)
      real(wp) :: k_solid(n_soil_layer_max), k_dry(n_soil_layer_max), dry_cvol(n_soil_layer_max)
   end type soil_texture_layers_t

   !----- Per-column geometry and texture, assembled once per site: elevations z <= 0      !
   !      below ground, dz and dz_node > 0. Fixed size (n_soil_layer_max), so the          !
   !      hydrology kernel stays allocatable-free and GPU-eligible. -----------------------!
   type :: soil_params_t
      integer(ik) :: n_active  = n_soil_layer_max            !< active layer count (<= n_soil_layer_max)
      integer(ik) :: retention = SOIL_RETENTION_VG           !< curve family
      real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp  !< [m] interface elevations (<= 0)
      real(wp) :: z_node(n_soil_layer_max)  = 0.0_wp         !< [m] node (mid) elevations (<= 0)
      real(wp) :: dz(n_soil_layer_max)      = 0.0_wp         !< [m] layer thickness (> 0)
      real(wp) :: dz_node(n_soil_layer_max) = 0.0_wp         !< [m] internode spacing (> 0)
      real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp       !< [m3/m3] porosity
      real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp       !< [m3/m3] residual water content
      real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp       !< [m/s] saturated conductivity
      !----- One pair for both curve families (docs/science/soil_biophysics.md). ----------!
      real(wp) :: curve_par_a(n_soil_layer_max) = 0.0_wp     !< [1/m] VG alpha, or [m] Campbell psi_sat
      real(wp) :: curve_par_n(n_soil_layer_max) = 0.0_wp     !< [-] VG n, or Campbell b
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
   end type soil_params_t

   !----- Per-column soil thermal texture (geometry and porosity are in soil_params_t). ---!
   type :: soil_thermal_params_t
      integer(ik) :: nzg_active = n_soil_layer_max
      real(wp) :: soil_solid_conductivity(n_soil_layer_max) = 0.0_wp   !< [W/m/K] kappa_solid
      real(wp) :: soil_dry_conductivity(n_soil_layer_max)   = 0.0_wp   !< [W/m/K] kappa_dry
      real(wp) :: soil_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   !< [J/m3/K] dry-matrix vol. heat cap
   end type soil_thermal_params_t

end module meds_soil_column_config
