! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_column_config -- the soil column a run is set on: the column keys of [soil],  !
! the rule that cuts the column into layers, and the parameters derived from it once per   !
! column.                                                                                  !
!                                                                                          !
! soil_column_params turns the block into the two bundles the soil kernels read:          !
! soil_params_t (the layers, the hydraulic texture and the field capacity) and             !
! soil_thermal_params_t (the thermal texture). They are derived from the config and never  !
! integrated -- parameters, not state -- so they live with the config, as the PFT table's  !
! derived rates do (meds_pft_params). A config leaf: it uses only src/shared/, so the soil !
! kernels import it without the meds_config aggregator.                                    !
!==========================================================================================!
module meds_soil_column_config
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : n_soil_layer_max, grav_head
   use meds_water_retention,  only : soil_theta_from_psi, soil_psi_from_theta, SOIL_RETENTION_VG
   implicit none
   private

   public :: soil_column_config_t, soil_texture_layers_t, soil_params_t, soil_thermal_params_t
   public :: soil_column_layers, soil_column_texture, soil_column_params
   public :: soil_column_start_theta, top_soil_psi
   public :: SOIL_LAYER_TEMPLATE, soil_layer_grid
   public :: build_soil_hydr_params, build_soil_therm_params
   public :: build_soil_hydr_params_layers, build_soil_therm_params_layers

   !----- The soil layers every column cuts at its own depth                               !
   !      (docs/science/soil_biophysics.md, "The soil layers"): bottoms [m, down] at       !
   !      SoilGrids 2.0's depths to 1 m (0-5, 5-15, 15-30, 30-60, 60-100 cm), its          !
   !      100-200 cm interval split at 1.5 m, then thickening to 27 m (each layer at       !
   !      most twice the one above). Every column's top layer is 5 cm, a shallow soil      !
   !      gets few layers rather than thin ones, per-layer soil maps match the layers      !
   !      one to one, and every column shares one depth axis. -----------------------------!
   real(wp), parameter :: SOIL_LAYER_TEMPLATE(15) = [0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp,     &
                          1.50_wp, 2.00_wp, 3.00_wp, 4.50_wp, 6.50_wp, 9.00_wp, 12.0_wp, 16.0_wp, 21.0_wp, 27.0_wp]
   !----- A remainder below the cut thinner than this share of the list layer it           !
   !      falls in joins the layer above, so no column ends in a sliver. ----------------- !
   real(wp), parameter :: LAYER_MERGE_SHARE = 0.5_wp
   real(wp), parameter :: GRID_TOL = 1.0e-6_wp     !< [m] a depth this close to a list bottom ends on it

   !----- The SOIL COLUMN ([soil]): the site's geometry, texture and thermal properties            !
   !      (docs/science/soil_biophysics.md). soil_column_params turns them into soil_params_t +     !
   !      soil_thermal_params_t. Every key is OPTIONAL and defaults to the value below. [soil] also  !
   !      holds the ground's surface and its two bottom boundaries, which are read into the solver  !
   !      options (soil_opts_t, energy_opts_t in meds_biophysics_opts). ----------------------------!
   type :: soil_column_config_t
      !----- Vertical grid: SOIL_LAYER_TEMPLATE, or the `layer_bottoms` list when          !
      !      set, cut at `depth` (soil_column_layers). ----------------------------------- !
      real(wp)    :: depth           = 2.0_wp    !< [m]  total column depth (positive; z is negative down)
      integer(ik) :: n_layer_bottoms = 0_ik      !< [-]  entries of layer_bottoms; 0 = the template
      real(wp)    :: layer_bottoms(n_soil_layer_max) = 0.0_wp   !< [m] layer bottoms, positive down
      !----- Hydraulic texture, one for the whole column (or per layer: n_layer_values below). ------!
      integer(ik) :: retention   = SOIL_RETENTION_VG  !< van_genuchten | campbell
      real(wp)    :: theta_sat   = 0.43_wp    !< [m3/m3] porosity
      real(wp)    :: theta_res   = 0.078_wp   !< [m3/m3] residual water content
      real(wp)    :: ksat        = 2.89e-6_wp !< [m/s]   saturated hydraulic conductivity
      real(wp)    :: curve_par_a = 3.6_wp     !< [1/m] van Genuchten alpha, OR [m] Campbell psi_sat (< 0)
      real(wp)    :: curve_par_n = 1.56_wp    !< [-]   van Genuchten n (>1), OR [-] Campbell b
      real(wp)    :: psi_fc      = -3.37_wp   !< [m]   field-capacity matric head (derives theta_fc)
      !----- Thermal texture (uniform over the column). ---------------------------------------------!
      real(wp)    :: solid_conductivity = 3.0_wp    !< [W/m/K]  mineral-solid conductivity
      real(wp)    :: dry_conductivity   = 0.15_wp   !< [W/m/K]  dry-matrix conductivity
      real(wp)    :: dry_heat_capacity  = 2.0e6_wp  !< [J/m3/K] dry-matrix volumetric heat capacity
      !----- A texture per layer, from the soil maps ([maps], meds_surface_maps), in place of the  !
      !      one texture above: n_layer_values is the column's layer count, and each array holds    !
      !      one value per layer. 0 = the one texture. No config key sets these. -------------------!
      integer(ik) :: n_layer_values = 0_ik
      real(wp)    :: layer_theta_sat(n_soil_layer_max)          = 0.0_wp
      real(wp)    :: layer_theta_res(n_soil_layer_max)          = 0.0_wp
      real(wp)    :: layer_ksat(n_soil_layer_max)               = 0.0_wp
      real(wp)    :: layer_curve_par_a(n_soil_layer_max)        = 0.0_wp
      real(wp)    :: layer_curve_par_n(n_soil_layer_max)        = 0.0_wp
      real(wp)    :: layer_solid_conductivity(n_soil_layer_max) = 0.0_wp
      real(wp)    :: layer_dry_conductivity(n_soil_layer_max)   = 0.0_wp
      real(wp)    :: layer_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp
   end type soil_column_config_t

   !----- The column's texture, one value per layer: the per-layer arrays, or the one texture   !
   !      repeated (soil_column_texture). --------------------------------------------------------!
   type :: soil_texture_layers_t
      real(wp) :: theta_sat(n_soil_layer_max), theta_res(n_soil_layer_max), ksat(n_soil_layer_max)
      real(wp) :: par_a(n_soil_layer_max), par_n(n_soil_layer_max)
      real(wp) :: k_solid(n_soil_layer_max), k_dry(n_soil_layer_max), dry_cvol(n_soil_layer_max)
   end type soil_texture_layers_t

   !----- Per-column geometry + texture (assembled once per site; ED2 negative-z convention:   !
   !      interface elevations z <= 0 below ground; dz, dz_node are positive magnitudes). Fixed- !
   !      size (n_soil_layer_max) so the hydrology kernel stays allocatable-free and GPU-eligible.!
   type :: soil_params_t
      integer(ik) :: n_active  = n_soil_layer_max            !< active layer count (<= n_soil_layer_max)
      integer(ik) :: retention = SOIL_RETENTION_VG           !< curve family
      real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp  !< [m] interface elevations (<= 0, ED2 slz)
      real(wp) :: z_node(n_soil_layer_max)  = 0.0_wp         !< [m] node (mid) elevations (<= 0)
      real(wp) :: dz(n_soil_layer_max)      = 0.0_wp         !< [m] layer thickness (> 0)
      real(wp) :: dz_node(n_soil_layer_max) = 0.0_wp         !< [m] internode spacing (> 0)
      real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp       !< [m3/m3] porosity
      real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp       !< [m3/m3] residual water content
      real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp       !< [m/s] saturated conductivity
      !----- The retention curve's two parameters, whose meaning the family sets: van Genuchten's !
      !      alpha [1/m] and n (> 1), or Campbell's air-entry head psi_sat [m, < 0] and exponent b !
      !      (docs/science/soil_biophysics.md, "The soil water column"). One pair for both        !
      !      families, so every reader -- the water network, the predawn seed, the diagnostics -- !
      !      reads the curve the column was built with. ------------------------------------------!
      real(wp) :: curve_par_a(n_soil_layer_max) = 0.0_wp     !< [1/m] VG alpha, or [m] Campbell psi_sat
      real(wp) :: curve_par_n(n_soil_layer_max) = 0.0_wp     !< [-] VG n, or Campbell b
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
   end type soil_params_t

   !----- Per-column soil THERMAL texture (geometry + porosity arrive via soil_params_t). ----!
   type :: soil_thermal_params_t
      integer(ik) :: nzg_active = n_soil_layer_max
      real(wp) :: soil_solid_conductivity(n_soil_layer_max) = 0.0_wp   !< [W/m/K] kappa_solid
      real(wp) :: soil_dry_conductivity(n_soil_layer_max)   = 0.0_wp   !< [W/m/K] kappa_dry
      real(wp) :: soil_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   !< [J/m3/K] dry-matrix vol. heat cap
   end type soil_thermal_params_t

contains

   !---------------------------------------------------------------------------------------!
   ! The column's layers: the bottoms [m, positive down] of `layer_bottoms` when the       !
   ! config sets it, else SOIL_LAYER_TEMPLATE, cut at `depth` (soil_layer_grid).           !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_column_layers(sc, bottoms, n_layer)
      type(soil_column_config_t), intent(in)  :: sc
      real(wp),                   intent(out) :: bottoms(n_soil_layer_max)
      integer(ik),                intent(out) :: n_layer
      if (sc%n_layer_bottoms > 0_ik) then
         call soil_layer_grid(sc%depth, sc%layer_bottoms, sc%n_layer_bottoms, bottoms, n_layer)
      else
         call soil_layer_grid(sc%depth, SOIL_LAYER_TEMPLATE, size(SOIL_LAYER_TEMPLATE, kind=ik),   &
                              bottoms, n_layer)
      end if
   end subroutine soil_column_layers

   !----- The column's texture per layer: the map's per-layer values, or the one texture. -------!
   pure subroutine soil_column_texture(sc, tex)
      type(soil_column_config_t),  intent(in)  :: sc
      type(soil_texture_layers_t), intent(out) :: tex
      if (sc%n_layer_values > 0_ik) then
         tex%theta_sat = sc%layer_theta_sat ; tex%theta_res = sc%layer_theta_res ; tex%ksat = sc%layer_ksat
         tex%par_a     = sc%layer_curve_par_a ; tex%par_n   = sc%layer_curve_par_n
         tex%k_solid   = sc%layer_solid_conductivity ; tex%k_dry = sc%layer_dry_conductivity
         tex%dry_cvol  = sc%layer_dry_heat_capacity
      else
         tex%theta_sat = sc%theta_sat ; tex%theta_res = sc%theta_res ; tex%ksat = sc%ksat
         tex%par_a     = sc%curve_par_a ; tex%par_n   = sc%curve_par_n
         tex%k_solid   = sc%solid_conductivity ; tex%k_dry = sc%dry_conductivity
         tex%dry_cvol  = sc%dry_heat_capacity
      end if
   end subroutine soil_column_texture

   !---------------------------------------------------------------------------------------!
   ! The soil column [soil] describes: its layers (soil_column_layers),                    !
   ! each layer's retention curve and Ks, and, when asked, each layer's thermal texture.   !
   ! The fast loop, the slow ledger and validate_config all build the column here.         !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_column_params(sc, soil, therm)
      type(soil_column_config_t),            intent(in)  :: sc
      type(soil_params_t),                   intent(out) :: soil
      type(soil_thermal_params_t), optional, intent(out) :: therm
      type(soil_texture_layers_t) :: tex
      real(wp)    :: bottoms(n_soil_layer_max)
      integer(ik) :: n_layer
      call soil_column_layers(sc, bottoms, n_layer)
      call soil_column_texture(sc, tex)
      call build_soil_hydr_params_layers(n_layer, sc%retention, bottoms, tex%theta_sat, tex%theta_res,  &
                                         tex%ksat, tex%par_a, tex%par_n, sc%psi_fc, soil)
      if (present(therm)) call build_soil_therm_params_layers(n_layer, tex%k_solid, tex%k_dry,        &
                                                              tex%dry_cvol, therm)
   end subroutine soil_column_params

   !---------------------------------------------------------------------------------------!
   ! The water each layer starts the run with [m3/m3]: [init].soil_theta, or each layer's   !
   ! field capacity for a column with a texture per layer (the soil maps).                  !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_column_start_theta(sc, soil, init_soil_theta, theta)
      type(soil_column_config_t), intent(in)  :: sc
      type(soil_params_t),        intent(in)  :: soil
      real(wp),                   intent(in)  :: init_soil_theta   !< [m3/m3] [init].soil_theta
      real(wp),                   intent(out) :: theta(n_soil_layer_max)
      if (sc%n_layer_values > 0_ik) then
         theta = soil%theta_fc
      else
         theta = init_soil_theta
      end if
   end subroutine soil_column_start_theta

   !---------------------------------------------------------------------------------------!
   ! The water potential [MPa] of the column's top layer at water content `theta_top`: the  !
   ! predawn potential of a plant that has relaxed to the soil overnight. It seeds the      !
   ! stomatal predawn potential and a new cohort's tissue water.                           !
   !---------------------------------------------------------------------------------------!
   pure real(wp) function top_soil_psi(soil, theta_top) result(psi)
      type(soil_params_t), intent(in) :: soil
      real(wp),            intent(in) :: theta_top   !< [m3/m3]
      psi = grav_head * soil_psi_from_theta(soil%retention, theta_top, soil%theta_sat(1), soil%theta_res(1), &
                                            soil%curve_par_a(1), soil%curve_par_n(1))
   end function top_soil_psi

   !=======================================================================================!
   !  Per-column soil PARAMETER assemblers. Both are `pure` -- they take plain scalar      !
   !  texture/geometry inputs and fill a parameter struct; they depend on NO soil STATE    !
   !  (theta/energy), so they are state-free constructors, not part of the fast loop.      !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! A column's layer bottoms [m, positive down]: `list` (SOIL_LAYER_TEMPLATE or a         !
   ! config list) cut at `depth`. The last layer ends at the depth; a remainder thinner    !
   ! than LAYER_MERGE_SHARE of the list layer it falls in joins the layer above. A list    !
   ! that ends above the depth gets one more layer down to it (the config refuses one).    !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_layer_grid(depth, list, n_list, bottoms, n_layer)
      real(wp),    intent(in)  :: depth, list(:)
      integer(ik), intent(in)  :: n_list
      real(wp),    intent(out) :: bottoms(n_soil_layer_max)
      integer(ik), intent(out) :: n_layer
      integer(ik) :: k
      real(wp)    :: top
      bottoms = 0.0_wp
      n_layer = 0_ik
      top     = 0.0_wp
      do k = 1_ik, n_list
         if (list(k) >= depth - GRID_TOL) exit                 ! list layer k holds the cut
         n_layer = n_layer + 1_ik
         bottoms(n_layer) = list(k)
         top = list(k)
      end do
      if (k <= n_list .and. n_layer > 0_ik) then
         if (depth - top < LAYER_MERGE_SHARE * (list(k) - top)) n_layer = n_layer - 1_ik
      end if
      n_layer = n_layer + 1_ik
      bottoms(n_layer) = depth
   end subroutine soil_layer_grid

   !---------------------------------------------------------------------------------------!
   ! Assemble a per-column soil_params_t: the ED2 negative-z geometry from the layer       !
   ! bottoms (soil_layer_grid), each layer's own texture, and the DERIVED threshold        !
   ! theta_fc (from each layer's retention curve). The arrays hold one value per active     !
   ! layer: from a soil map's layers, or one value broadcast (build_soil_hydr_params). The  !
   ! roots are the cohorts' own (meds_plant_hydraulics::cohort_root_profile), not the       !
   ! column's.                                                                              !
   !---------------------------------------------------------------------------------------!
   pure subroutine build_soil_hydr_params_layers(n_active, retention, layer_bottoms, theta_sat,       &
                                                 theta_res, ksat, par_a, par_n, psi_fc_m, params)
      integer(ik),         intent(in)  :: n_active, retention
      real(wp),            intent(in)  :: layer_bottoms(:)    !< [m] bottoms of the active layers, positive down
      real(wp),            intent(in)  :: theta_sat(:), theta_res(:), ksat(:), par_a(:), par_n(:)
      real(wp),            intent(in)  :: psi_fc_m
      type(soil_params_t), intent(out) :: params
      integer(ik) :: k

      params%n_active  = n_active
      params%retention = retention

      !----- Interface elevations soil_layer_z(k) <= 0 (k=1 surface, deeper = more negative). -!
      params%soil_layer_z(1)              = 0.0_wp
      params%soil_layer_z(2:n_active + 1) = -layer_bottoms(1:n_active)

      !----- Thicknesses, node elevations, internode spacings (dz, dz_node > 0). ----------!
      do k = 1_ik, n_active
         params%dz(k)     = params%soil_layer_z(k) - params%soil_layer_z(k+1)
         params%z_node(k) = 0.5_wp * (params%soil_layer_z(k) + params%soil_layer_z(k+1))
      end do
      do k = 1_ik, n_active - 1_ik
         params%dz_node(k) = params%z_node(k) - params%z_node(k+1)
      end do
      params%dz_node(n_active) = params%dz(n_active)      ! unused at the bottom face; kept > 0

      !----- Each layer's texture. ---------------------------------------------------------!
      params%theta_sat(1:n_active) = theta_sat(1:n_active)
      params%theta_res(1:n_active) = theta_res(1:n_active)
      params%ksat(1:n_active)      = ksat(1:n_active)
      params%curve_par_a(1:n_active) = par_a(1:n_active)
      params%curve_par_n(1:n_active) = par_n(1:n_active)

      !----- Derived threshold (from each layer's retention curve). ----------------------!
      params%theta_fc(1:n_active) = soil_theta_from_psi(retention, psi_fc_m, theta_sat(1:n_active),    &
                                        theta_res(1:n_active), par_a(1:n_active), par_n(1:n_active))
   end subroutine build_soil_hydr_params_layers

   !----- One texture for every layer: the per-layer build with each value broadcast. --------!
   pure subroutine build_soil_hydr_params(n_active, retention, layer_bottoms, theta_sat, theta_res,  &
                                ksat, par_a, par_n, psi_fc_m, params)
      integer(ik),         intent(in)  :: n_active, retention
      real(wp),            intent(in)  :: layer_bottoms(:)    !< [m] bottoms of the active layers, positive down
      real(wp),            intent(in)  :: theta_sat, theta_res, ksat, par_a, par_n, psi_fc_m
      type(soil_params_t), intent(out) :: params
      real(wp) :: ts(n_active), tr(n_active), ks(n_active), pa(n_active), pn(n_active)
      ts = theta_sat ; tr = theta_res ; ks = ksat ; pa = par_a ; pn = par_n
      call build_soil_hydr_params_layers(n_active, retention, layer_bottoms, ts, tr, ks, pa, pn,        &
                                         psi_fc_m, params)
   end subroutine build_soil_hydr_params

   !----- The column's thermal texture, one value per active layer. ---------------------------!
   pure subroutine build_soil_therm_params_layers(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                 intent(in)  :: n_active
      real(wp),                    intent(in)  :: k_solid(:), k_dry(:), dry_cvol(:)
      type(soil_thermal_params_t), intent(out) :: therm
      therm%nzg_active = n_active
      therm%soil_solid_conductivity(1:n_active) = k_solid(1:n_active)
      therm%soil_dry_conductivity(1:n_active)   = k_dry(1:n_active)
      therm%soil_dry_heat_capacity(1:n_active)  = dry_cvol(1:n_active)
   end subroutine build_soil_therm_params_layers

   !----- One thermal texture for every layer. ------------------------------------------------!
   pure subroutine build_soil_therm_params(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                 intent(in)  :: n_active
      real(wp),                    intent(in)  :: k_solid, k_dry, dry_cvol
      type(soil_thermal_params_t), intent(out) :: therm
      real(wp) :: ks(n_active), kd(n_active), cv(n_active)
      ks = k_solid ; kd = k_dry ; cv = dry_cvol
      call build_soil_therm_params_layers(n_active, ks, kd, cv, therm)
   end subroutine build_soil_therm_params

end module meds_soil_column_config
