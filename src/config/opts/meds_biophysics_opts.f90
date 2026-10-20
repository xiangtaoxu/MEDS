! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biophysics_opts -- the option bundles of the fast biophysics loop, from the [soil], !
! [fast], [snow] and [aerodynamics] blocks: the soil-water and soil-energy solver settings !
! and the snow and aerodynamics parameters. A config leaf (meds_kinds only), so the        !
! device-eligible kernels import it without meds_config. The SOIL_RETENTION_* selectors    !
! are in meds_soil_hydr_lib, the per-column soil descriptors in meds_soil_column_config.   !
!==========================================================================================!
module meds_biophysics_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   !----- Soil-water bottom boundary selectors. -------------------------------------------!
   public :: SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   !----- Soil-energy solver selectors. ---------------------------------------------------!
   public :: ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET
   public :: ANNUAL_DAMPING_DEPTH
   !----- The option / parameter bundles. -------------------------------------------------!
   public :: soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t

   !----- Soil-water bottom boundary: drainage at a fraction f of the deepest layer's      !
   !      conductivity, q = f K(theta_N), or an aquifer (docs/science/soil_biophysics.md,  !
   !      "Boundary conditions"). ---------------------------------------------------------!
   integer(ik), parameter :: SOIL_BC_DRAIN   = 1_ik     !< drainage q = f K(theta_N)  (default)
   integer(ik), parameter :: SOIL_BC_AQUIFER = 2_ik     !< water table at the column base (wet sites only)

   !----- Soil-energy bottom boundary (docs/science/soil_biophysics.md). ------------------!
   integer(ik), parameter :: ENERGY_BC_GEOTHERMAL   = 1_ik   !< bottom: prescribed flux, held at zero (Neumann)
   integer(ik), parameter :: ENERGY_BC_DIRICHLET     = 2_ik   !< bottom: conduction to a fixed deep temperature
   !----- Annual damping depth d = sqrt(2*alpha/omega) of a mid-latitude mineral soil at   !
   !      theta = 0.3. The Dirichlet anchor reflects the annual wave least d/sqrt(2) below !
   !      the bottom node (docs/science/soil_biophysics.md eq. (13)): deep_depth's         !
   !      default, and a soil map's column. -----------------------------------------------!
   real(wp), parameter :: ANNUAL_DAMPING_DEPTH = 1.973_wp   !< [m]

   !----- Soil water: the bottom boundary, the water network's tolerances, the pond, the   !
   !      canopy interception, and the ground's evaporation and optics. -------------------!
   type :: soil_opts_t
      integer(ik) :: water_bottom_bc = SOIL_BC_DRAIN    !< SOIL_BC_DRAIN | SOIL_BC_AQUIFER
      real(wp)    :: drain_factor = 1.0_wp              !< [-] f in [0, 1] of q = f K(theta_N): free_drain 1, bedrock 0
      !----- Water-network sub-step tolerances (docs/science/soil_plant_water.md §3). -----!
      real(wp)    :: rtol   = 1.0e-3_wp                 !< [-] [fast].rtol, which the driver copies here
      real(wp)    :: water_atol = 1.0e-4_wp             !< [-] [fast].water_atol
      real(wp)    :: w_pond_max = 5.0_wp                 !< [kg/m2] ponding capacity before surface runoff
      real(wp)    :: intercept_alpha = 1.0_wp            !< Beer interception efficiency ([fast])
      real(wp)    :: intercept_k     = 0.5_wp            !< [1/PAI] Beer extinction ([fast])
      real(wp)    :: dsl_dmax       = 0.015_wp           !< [m] DSL max thickness (soil evap)
      real(wp)    :: dsl_theta_init = 0.8_wp             !< DSL initiation theta_1/theta_sat
      logical     :: debug_error = .false.               !< [fast]: a budget breach stops the run
      !----- Bare-soil optics, blended with snow where there is snow. ---------------------!
      real(wp)    :: ground_albedo_vis = 0.15_wp        !< [-] soil albedo, visible
      real(wp)    :: ground_albedo_nir = 0.30_wp        !< [-] soil albedo, near infrared
      real(wp)    :: ground_emissivity = 0.95_wp        !< [-] soil longwave emissivity
   end type soil_opts_t

   !----- Soil energy: the bottom boundary and the debug budget check. The backward-Euler  !
   !      solve is unconditionally stable and the outer march controls its accuracy, so    !
   !      there are no sub-step controls (docs/science/numerical_scheme.md eq. (3)). ------!
   type :: energy_opts_t
      !----- The bottom boundary (docs/science/soil_biophysics.md eqs. (12)-(13)). --------!
      integer(ik) :: thermal_bottom_bc = ENERGY_BC_GEOTHERMAL   !< ENERGY_BC_GEOTHERMAL | ENERGY_BC_DIRICHLET
      real(wp)    :: deep_temp   = -1.0_wp        !< [K] anchor temperature; < 0: not set (dirichlet requires it)
      real(wp)    :: deep_depth  = 3.145_wp       !< [m] anchor plane depth: the default column's optimum, eq. (13)
      real(wp)    :: soil_energy_budget_atol = 1.0e-2_wp   !< [K] debug-halt budget threshold ([fast]); not a step tolerance
      logical     :: debug_error = .false.        !< [fast]: a budget breach stops the run
   end type energy_opts_t

   !----- The [snow] block (docs/science/snow_biophysics.md): one layer of fixed density   !
   !      and conductivity. ---------------------------------------------------------------!
   type :: snow_params_t
      real(wp) :: rho_snow          = 250.0_wp    !< [kg/m3] fixed seasonal-mean bulk density (depth = swe/rho)
      real(wp) :: k_snow            = 0.15_wp     !< [W/m/K] fixed snow thermal conductivity
      real(wp) :: albedo_vis_fresh  = 0.85_wp     !< [-] fresh-snow VIS albedo
      real(wp) :: albedo_nir_fresh  = 0.65_wp     !< [-] fresh-snow NIR albedo
      real(wp) :: albedo_vis_aged   = 0.518_wp    !< [-] aged VIS endpoint
      real(wp) :: albedo_nir_aged   = 0.435_wp    !< [-] aged NIR endpoint
      real(wp) :: snow_emiss        = 0.97_wp     !< [-] thermal (LW) emissivity
      real(wp) :: liquid_holding_frac = 0.10_wp   !< [-] percolation holding capacity (LEAF-3 1:9)
      real(wp) :: min_new_snow_mass = 1.0e-3_wp   !< [kg/m2] first-flake creation threshold
      real(wp) :: tiny_snow_mass    = 1.0e-6_wp   !< [kg/m2] melt-out threshold (dump residual to soil)
      real(wp) :: ny07_a            = 2.5_wp      !< [-] Niu-Yang07 snow-cover-fraction shape factor
      real(wp) :: ny07_m            = 1.0_wp      !< [-] Niu-Yang07 snow-cover-fraction exponent
   end type snow_params_t

   !----- The [aerodynamics] block (docs/science/canopy_aerodynamics.md); the von Karman   !
   !      constant is in meds_constants. --------------------------------------------------!
   type :: aero_cfg_t
      real(wp) :: z0m_ratio = 0.13_wp               !< z0m / canopy height
      real(wp) :: d_ratio   = 0.63_wp               !< displacement / canopy height
      real(wp) :: snow_rough = 0.001_wp             !< [m] snow-surface roughness (blend floor)
      !----- Monin-Obukhov (CLM5). --------------------------------------------------------!
      real(wp)    :: zeta_m = 1.574_wp, zeta_t = 0.465_wp   !< momentum / scalar range transitions
      real(wp)    :: zeta_max_stable = 0.5_wp       !< cap on the stable first guess of zeta
      integer(ik) :: n_iter_mo = 4_ik               !< fixed MO iterations (GPU-uniform, no root-find)
      real(wp)    :: wc = 0.5_wp                     !< [m/s] convective velocity scale
      !----- Floors. ----------------------------------------------------------------------!
      real(wp) :: ustmin = 0.10_wp, ubmin = 0.65_wp, ugbmin = 0.25_wp   !< [m/s] u*, bulk-wind and in-canopy wind
      real(wp) :: gbhmos_min = 1.0e-9_wp            !< [m/s] min boundary-layer conductance
      real(wp) :: min_canopy_depth = 5.0_wp         !< [m] CAS depth floor
      real(wp) :: canopy_freeboard  = 5.0_wp        !< [m] CAS depth above the tallest cohort
      !----- Ground conductance (CLM4 dense canopy). --------------------------------------!
      real(wp) :: cs_dense = 0.004_wp, gamma_g = 0.5_wp
      !----- Boundary-layer heat:vapour ratio + air properties (linear-in-T). -------------!
      real(wp) :: gbh_2_gbw = 1.075_wp
      real(wp) :: kin_visc0 = 1.33e-5_wp, dkin_visc = 0.007_wp   !< kinematic viscosity [m2/s] @ t_ref_air
      real(wp) :: th_diff0  = 1.89e-5_wp, dth_diff  = 0.007_wp   !< thermal diffusivity [m2/s] @ t_ref_air
      real(wp) :: t_ref_air = 293.15_wp             !< [K] 20 C reference for the linear-in-T props
      !----- Nusselt coefficients -- flat plate (leaf). -----------------------------------!
      real(wp) :: aflat_lami = 0.60_wp, nflat_lami = 0.50_wp
      real(wp) :: aflat_turb = 0.032_wp, nflat_turb = 0.80_wp
      real(wp) :: bflat_lami = 0.50_wp, mflat_lami = 0.25_wp
      real(wp) :: bflat_turb = 0.19_wp, mflat_turb = 1.0_wp/3.0_wp
      !----- Nusselt coefficients -- cylinder (wood). -------------------------------------!
      real(wp) :: ocyli_lami = 0.32_wp, acyli_lami = 0.51_wp, ncyli_lami = 0.52_wp
      real(wp) :: ocyli_turb = 0.0_wp,  acyli_turb = 0.24_wp, ncyli_turb = 0.60_wp
      real(wp) :: bcyli_lami = 0.48_wp, mcyli_lami = 0.25_wp
      real(wp) :: bcyli_turb = 0.09_wp, mcyli_turb = 1.0_wp/3.0_wp
   end type aero_cfg_t

end module meds_biophysics_opts
