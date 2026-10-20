! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biophysics_opts -- the run-configuration bundles of the fast (sub-daily) biophysics !
! loop: the solver selectors and tolerances (soil water, soil energy) and the physical     !
! parameter tables (snow, aerodynamics), filled from the [soil], [fast], [snow] and        !
! [aerodynamics] blocks, with the defaults below for standalone and test use. They are run !
! configuration, not prognostic state and not per-column descriptors, so they live in a    !
! low-level config leaf (`meds_kinds` only, like meds_pft_params) that the device-eligible !
! biophysics kernels can `use` without depending on the meds_config aggregator. The        !
! kernels that switch on them and the loader that reads them import this module directly.  !
!                                                                                          !
! (The SOIL_RETENTION_* selectors stay in meds_soil_hydr_lib with the retention curves;    !
! the per-column soil_params_t / soil_thermal_params_t descriptors are in                  !
! meds_soil_column_config.)                                                                !
!==========================================================================================!
module meds_biophysics_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   !----- Soil-water bottom boundary selectors. -------------------------------------------!
   public :: SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   !----- Soil-energy solver selectors. ----------------------------------------------------!
   public :: ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET
   public :: ANNUAL_DAMPING_DEPTH
   !----- The option / parameter bundles. --------------------------------------------------!
   public :: soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t

   !----- The column drains through its base at a fraction f of the deepest layer's conductivity,  !
   !      q = f K(theta_N) (free drainage f = 1, a leaky base 0 < f < 1, sealed bedrock f = 0), or it !
   !      sits on an aquifer (docs/science/soil_biophysics.md §Boundary conditions). ----------------!
   integer(ik), parameter :: SOIL_BC_DRAIN   = 1_ik     !< drainage q = f K(theta_N)  (default)
   integer(ik), parameter :: SOIL_BC_AQUIFER = 2_ik     !< water table at the column base (wet sites only)

   integer(ik), parameter :: ENERGY_BC_GEOTHERMAL   = 1_ik   !< bottom: prescribed flux (Neumann)
   integer(ik), parameter :: ENERGY_BC_DIRICHLET     = 2_ik   !< bottom: conduction to a fixed deep temperature
   !----- The annual damping depth d = sqrt(2*alpha/omega) of a mid-latitude mineral soil at theta = 0.3. !
   !      The Dirichlet anchor reflects the annual wave least at d/sqrt(2) below the bottom node: the   !
   !      default deep_depth below, and the one a soil map's column takes (meds_surface_maps). --------!
   real(wp), parameter :: ANNUAL_DAMPING_DEPTH = 1.973_wp   !< [m]

   !----- Soil water: the bottom boundary, the soil-plant water network's tolerances, the pond, the !
   !      canopy interception, the ground's evaporation and optics (NOT the whole config). ----!
   type :: soil_opts_t
      !----- Drainage (default) | aquifer, with drainage's factor f: [soil].water_bottom_bc =        !
      !      "free_drain" is f = 1, "bedrock" f = 0, and a soil map's leaky base its f_leak.  NOTE:   !
      !      `aquifer` is the WET-SITE boundary -- the water table is DEFINED as the column base, so  !
      !      the bottom face is head-driven against a saturated zone and supplies water upward       !
      !      without limit.  Right for riparian, floodplain and wetland stands; WRONG for a typical  !
      !      upland one, where free drainage is the correct choice.  It is not a drop-in "more        !
      !      realistic drainage" option. ---------------------------------------------------------!
      integer(ik) :: bottom_bc    = SOIL_BC_DRAIN
      real(wp)    :: drain_factor = 1.0_wp              !< [-] f in [0, 1] of q = f K(theta_N)
      !----- The water network's sub-step tolerances ([fast].water_rtol and water_atol;            !
      !      docs/science/soil_plant_water.md §3): a store's error is weighed against atol of its    !
      !      scale -- water content for soil layers and the pond, relative water content for plant   !
      !      tissues -- plus rtol of its water. ------------------------------------------------------!
      real(wp)    :: rtol   = 1.0e-3_wp
      real(wp)    :: atol   = 1.0e-4_wp
      real(wp)    :: w_pond_max = 5.0_wp                 !< [kg/m2] ponding capacity before surface runoff
      real(wp)    :: intercept_alpha = 1.0_wp            !< Beer interception efficiency ([fast])
      real(wp)    :: intercept_k     = 0.5_wp            !< [1/PAI] Beer extinction ([fast])
      real(wp)    :: dsl_dmax       = 0.015_wp           !< [m] DSL max thickness (soil evap)
      real(wp)    :: dsl_theta_init = 0.8_wp             !< DSL initiation theta_1/theta_sat
      logical     :: debug_error = .false.
      !----- Ground optics under the canopy: the bare soil's shortwave albedo per band and its        !
      !      longwave emissivity, as the canopy radiation solver sees the ground (blended with snow   !
      !      where there is snow). ------------------------------------------------------------------!
      real(wp)    :: ground_albedo_vis = 0.15_wp        !< [-] soil albedo, visible
      real(wp)    :: ground_albedo_nir = 0.30_wp        !< [-] soil albedo, near infrared
      real(wp)    :: ground_emissivity = 0.95_wp        !< [-] soil longwave emissivity
   end type soil_opts_t

   !----- Soil energy: the bottom boundary and the debug budget check (NOT the whole config).      !
   !      The soil-thermal solve is one backward-Euler step, unconditionally stable, and the       !
   !      integrator's adaptive step already controls the soil internal energy it carries          !
   !      (tolerance group GRP_SE), so there are no sub-step controls here. -----------------------!
   type :: energy_opts_t
      !----- geothermal (default) | dirichlet (docs/science/soil_biophysics.md, "The bottom          !
      !      thermal boundary"; the discretisation is in soil_heat_be_solve). `geothermal`           !
      !      prescribes the bottom FLUX (held at zero), which makes the base adiabatic and           !
      !      REFLECTS the annual temperature wave; `dirichlet` conducts to a fixed deep              !
      !      temperature `deep_temp` at depth `deep_depth`, which absorbs most of it. ---------------!
      integer(ik) :: bottom_bc    = ENERGY_BC_GEOTHERMAL
      !----- Dirichlet anchor (read only when bottom_bc == ENERGY_BC_DIRICHLET).  `deep_temp` is a     !
      !      SITE property -- the mean annual soil temperature below the annual damping depth, close    !
      !      to the mean annual air temperature -- so it has no defensible default and the loader        !
      !      REQUIRES it.  An error in it is a steady flux `kappa/l_deep * error` into the column base,  !
      !      i.e. a mean-annual bias of roughly the same size at depth, so guessing it is not safe.      !
      real(wp)    :: deep_temp   = -1.0_wp        !< [K] anchor temperature; < 0 means "not set"
      !----- Depth of the anchor PLANE below the surface. The bottom node conducts to it over        !
      !      l_deep = deep_depth - |z_node(n)|, and a resistive termination reflects the annual      !
      !      wave least at l_deep = d/sqrt(2), d the annual damping depth                            !
      !      (docs/science/soil_biophysics.md eq. (13)). The default, 1.75 + 1.973/sqrt(2) =         !
      !      3.145 m, is that optimum for the default 2 m column in a mid-latitude mineral soil      !
      !      at theta = 0.3; change it with the column depth or the soil's thermal texture. ---------!
      real(wp)    :: deep_depth  = 3.145_wp       !< [m] depth (positive) of the anchor plane
      !----- [K] budget-closure threshold for the debug halt, NOT an integrator tolerance              !
      !      ([fast].soil_energy_budget_atol). ---------------------------------------------------------!
      real(wp)    :: atol        = 1.0e-2_wp
      logical     :: debug_error = .false.
   end type energy_opts_t

   !----- Snow / temporary-surface-water physical parameters, the [snow] block             !
   !      (docs/science/snow_biophysics.md). One layer, with a fixed density and           !
   !      conductivity (multi-layer snow: see docs/ROADMAP.md). ---------------------------!
   type :: snow_params_t
      real(wp) :: rho_snow          = 250.0_wp    !< [kg/m3] fixed SEASONAL-mean bulk density (depth = swe/rho)
      real(wp) :: k_snow            = 0.15_wp     !< [W/m/K] fixed snow thermal conductivity
      real(wp) :: albedo_vis_fresh  = 0.85_wp     !< [-] fresh-snow VIS albedo
      real(wp) :: albedo_nir_fresh  = 0.65_wp     !< [-] fresh-snow NIR albedo
      real(wp) :: albedo_vis_aged   = 0.518_wp    !< [-] aged VIS endpoint (ED2)
      real(wp) :: albedo_nir_aged   = 0.435_wp    !< [-] aged NIR endpoint (ED2)
      real(wp) :: snow_emiss        = 0.97_wp     !< [-] thermal (LW) emissivity
      real(wp) :: liquid_holding_frac = 0.10_wp   !< [-] percolation holding capacity (LEAF-3 1:9)
      real(wp) :: snow_stab_thresh  = 3.0_wp      !< [kg/m2] active-integration floor; thinner packs stay passive
      real(wp) :: min_new_snow_mass = 1.0e-3_wp   !< [kg/m2] first-flake creation threshold (ED2 tiny_sfcwater_mass)
      real(wp) :: tiny_snow_mass    = 1.0e-6_wp   !< [kg/m2] melt-out threshold (dump residual to soil)
      real(wp) :: ny07_a            = 2.5_wp      !< [-] Niu-Yang07 snow-cover-fraction shape factor
      real(wp) :: ny07_m            = 1.0_wp      !< [-] Niu-Yang07 snow-cover-fraction exponent
      real(wp) :: z0_snow           = 0.001_wp    !< [m] snow surface roughness (aero burial floor)
   end type snow_params_t

   !----- Canopy-aerodynamics run constants, the [aerodynamics] block                       !
   !      (docs/science/canopy_aerodynamics.md); defaults are the ED2/CLM5 reference        !
   !      values. (The von Karman constant lives in meds_constants.) -----------------------!
   type :: aero_cfg_t
      real(wp) :: z0m_ratio = 0.13_wp               !< z0m / canopy height (ED2 vh2vr)
      real(wp) :: d_ratio   = 0.63_wp               !< displacement / canopy height (ED2 vh2dh)
      real(wp) :: snow_rough = 0.001_wp             !< [m] snow-surface roughness (blend floor)
      !----- Monin-Obukhov (CLM5). --------------------------------------------------------!
      real(wp)    :: zeta_m = 1.574_wp, zeta_t = 0.465_wp   !< momentum / scalar range transitions
      real(wp)    :: zeta_max_stable = 0.5_wp       !< CLM5 zetamaxstable (init-guess cap)
      integer(ik) :: n_iter_mo = 4_ik               !< fixed MO iterations (GPU-uniform, no root-find)
      real(wp)    :: wc = 0.5_wp                     !< convective velocity scale (MoninObukIni)
      !----- Floors. ----------------------------------------------------------------------!
      real(wp) :: ustmin = 0.10_wp, ubmin = 0.65_wp, ugbmin = 0.25_wp
      real(wp) :: gbhmos_min = 1.0e-9_wp            !< [m/s] min boundary-layer conductance
      real(wp) :: min_canopy_depth = 5.0_wp         !< [m] CAS depth floor
      !----- FREEBOARD above the tallest cohort. The canopy air space is not the canopy volume:   !
      !      it is the well-mixed layer the canopy exchanges with, which extends above the crowns. !
      !      can_depth = tallest cohort height + freeboard (floored by min_canopy_depth).          !
      real(wp) :: canopy_freeboard  = 5.0_wp        !< [m] CAS depth above the tallest cohort
      !----- Ground conductance (CLM4-like dense-canopy). ---------------------------------!
      real(wp) :: cs_dense = 0.004_wp, gamma_g = 0.5_wp
      !----- Boundary-layer heat:vapour ratio + air properties (linear-in-T). --------------!
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
