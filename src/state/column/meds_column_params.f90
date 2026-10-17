! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_column_params -- everything about a column that is NOT state: the compile-time      !
! dimensions, the fresh-cohort initial values, and the static per-column parameter bundles !
! (geometry, texture, thermal properties) with their `pure` assemblers.                    !
!                                                                                          !
! These describe the same stores meds_column_state_types holds, but they are fixed or      !
! derived once per column and never integrated, which is the line src/README.md's rule 3   !
! draws (parameters are not state). Three kinds sit here, and they are worth telling       !
! apart:                                                                                   !
!   * DIMENSIONS (n_soil_layer_max, n_snow_layer_max) -- compile-time ceilings, not        !
!     settings. They are what keep the columns allocatable-free and GPU-eligible.          !
!   * BIRTH VALUES (LEAF_TEMP_INIT, PSI_INIT) -- the state a fresh cohort starts from.     !
!   * DERIVED BUNDLES (soil_params_t, soil_thermal_params_t) -- assembled per column at    !
!     run time from scalar geometry and texture inputs. NOT config: they stand to          !
!     `[soil_column]` as leaf_photo_table_t stands to the PFT table, which is why they     !
!     live beside the stores they describe rather than in src/config.                      !
!                                                                                          !
! The builders are state-free constructors -- they read plain scalar texture/geometry      !
! inputs and the shared retention curve, never soil state -- so this module stays below    !
! the kernels.                                                                             !
!==========================================================================================!
module meds_column_params
   use meds_kinds,            only : wp, ik
   use meds_water_retention,  only : soil_theta_from_psi, SOIL_RETENTION_VG
   implicit none
   private

   public :: n_soil_layer_max, n_snow_layer_max, LEAF_TEMP_INIT, PSI_INIT
   public :: C2B_WOOD, TISSUE_C_LEAF, TISSUE_C_WOOD_DRY, TISSUE_HCAP_MIN
   public :: soil_params_t, soil_thermal_params_t
   public :: build_soil_hydr_params, build_soil_therm_params
   public :: build_soil_hydr_params_layers, build_soil_therm_params_layers
   public :: SOIL_LAYER_TEMPLATE, soil_layer_grid

   integer(ik), parameter :: n_soil_layer_max = 20_ik      !< compile-time soil-column-depth ceiling
   integer(ik), parameter :: n_snow_layer_max = 1_ik        !< one bulk layer (multi-layer snow: see docs/ROADMAP.md)

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

   !----- Per-cohort fast state carried on cohort_block (rides the cohort lockstep): a fresh      !
   !      cohort starts at a mild tension + a neutral leaf temperature (both relax within one fast !
   !      step).                                                                                   !
   real(wp),    parameter :: LEAF_TEMP_INIT = 288.15_wp     !< [K]   fresh-cohort leaf temperature
   real(wp),    parameter :: PSI_INIT       = -0.1_wp       !< [MPa] fresh-cohort node water potential

   !----- A tissue's dry heat capacity (docs/science/vegetation_energy_dynamics.md eq. (8)). One set  !
   !      for the fast loop (veg_thermal_params_t starts from it, build_column_frozen) and the slow   !
   !      side (cohort_tissue_heat_capacity: the fusion weights and the slow ledger), so both weigh  !
   !      a tissue's heat alike. Hard-coded parameters, not config keys. --------------------------!
   real(wp),    parameter :: C2B_WOOD          = 2.0_wp    !< [-] carbon -> dry biomass (carbon fraction 0.5)
   real(wp),    parameter :: TISSUE_C_LEAF     = 3200.0_wp !< [J/kg/K] dry leaf specific heat
   real(wp),    parameter :: TISSUE_C_WOOD_DRY = 1217.0_wp !< [J/kg/K] dry wood specific heat (ED2's, at 15 C)
   real(wp),    parameter :: TISSUE_HCAP_MIN   = 20.0_wp   !< [J/m2 tissue/K] least dry heat capacity

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

end module meds_column_params
