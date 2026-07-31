!==========================================================================================!
! meds_biophysics_opts -- the run-configuration bundles of the fast (sub-daily) biophysics loop: !
! solver selectors + tolerances (soil water / soil energy) and the physical parameter tables      !
! (snow, aerodynamics). These are NOT prognostic state and NOT per-column descriptors -- they are  !
! run-config, filled once (defaults here for standalone/test use; the tunable subset from the       !
! [soil]/[energy]/[snow]/[aerodynamics] TOML blocks when the wiring lands). So they live in the      !
! shared config layer, a LOW-LEVEL leaf (`meds_kinds` only, like meds_pft_params) that the sealed    !
! device-eligible biophysics kernels can `use` WITHOUT depending on the heavyweight meds_config       !
! aggregator. meds_biophysics_types re-exports every name below so the kernels + callers keep         !
! `use meds_biophysics_types, only : soil_opts_t` unchanged.                                          !
!                                                                                          !
! (The constitutive SOIL_RETENTION_* selectors stay in meds_hydr_lib with the retention curves; the   !
! per-column soil_params_t / soil_thermal_params_t descriptors stay in meds_column_state_types beside  !
! the prognostic soil columns they describe and travel with.)                                          !
!==========================================================================================!
module meds_biophysics_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   !----- Soil-water solver selectors. -----------------------------------------------------!
   public :: SOIL_SOLVER_BE
   public :: SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER, SOIL_BC_BEDROCK
   public :: SOIL_LIN_FROZEN, SOIL_LIN_PICARD
   public :: SOIL_SUBSTEP_ADAPTIVE, SOIL_SUBSTEP_FIXED
   !----- Soil-energy solver selectors. ----------------------------------------------------!
   public :: ENERGY_SOLVER_BE, ENERGY_BC_GEOTHERMAL
   public :: ENERGY_PHASE_OFF, ENERGY_PHASE_ON, ENERGY_SUBSTEP_ADAPTIVE
   !----- The option / parameter bundles. --------------------------------------------------!
   public :: soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t

   integer(ik), parameter :: SOIL_SOLVER_BE = 1_ik         !< linearly-implicit backward-Euler (only solver)

   integer(ik), parameter :: SOIL_BC_FREE_DRAIN = 1_ik     !< unit-gradient q = K(theta_N)  (MVP default)
   integer(ik), parameter :: SOIL_BC_AQUIFER    = 2_ik     !< SIMTOP aquifer / water table (P2)
   integer(ik), parameter :: SOIL_BC_BEDROCK    = 3_ik     !< no-flow  q = 0

   integer(ik), parameter :: SOIL_LIN_FROZEN = 1_ik        !< frozen-coefficient single solve (MVP)
   integer(ik), parameter :: SOIL_LIN_PICARD = 2_ik        !< Celia modified-Picard (P2)

   integer(ik), parameter :: SOIL_SUBSTEP_ADAPTIVE = 1_ik  !< adaptive step-doubling
   integer(ik), parameter :: SOIL_SUBSTEP_FIXED    = 2_ik  !< fixed count (GPU warp-uniform)

   integer(ik), parameter :: ENERGY_SOLVER_BE       = 1_ik   !< implicit backward-Euler (only solver)
   integer(ik), parameter :: ENERGY_BC_GEOTHERMAL   = 1_ik   !< bottom: zero/geothermal flux
   integer(ik), parameter :: ENERGY_PHASE_OFF = 0_ik, ENERGY_PHASE_ON = 1_ik     !< freeze/thaw plateau (P1 off)
   integer(ik), parameter :: ENERGY_SUBSTEP_ADAPTIVE = 1_ik

   !----- Soil-water: pre-extracted solver selectors + tolerances (NOT the whole config). --!
   type :: soil_opts_t
      integer(ik) :: solver    = SOIL_SOLVER_BE
      !----- free_drain (default) | bedrock | aquifer.  NOTE: `aquifer` is the WET-SITE boundary --   !
      !      the water table is DEFINED as the column base, so the bottom face is head-driven against  !
      !      a saturated zone and supplies water upward without limit.  Right for riparian, floodplain !
      !      and wetland stands; WRONG for a typical upland one, where free_drain is the correct        !
      !      choice.  It is not a drop-in "more realistic drainage" option.  ------------------------!
      integer(ik) :: bottom_bc = SOIL_BC_FREE_DRAIN
      integer(ik) :: linearize = SOIL_LIN_FROZEN
      integer(ik) :: substep   = SOIL_SUBSTEP_ADAPTIVE
      real(wp)    :: rtol   = 1.0e-3_wp
      real(wp)    :: atol   = 1.0e-4_wp
      real(wp)    :: h_init = 900.0_wp                   !< [s] initial substep
      integer(ik) :: max_substep = 200_ik
      integer(ik) :: max_picard  = 5_ik
      real(wp)    :: w_pond_max = 5.0_wp                 !< [kg/m2] ponding capacity before surface runoff
      real(wp)    :: dewmx      = 0.1_wp                 !< [kg/m2 per PAI] canopy storage capacity
      real(wp)    :: intercept_alpha = 1.0_wp            !< Beer interception efficiency
      real(wp)    :: intercept_k     = 0.5_wp            !< [1/PAI] Beer extinction
      real(wp)    :: dsl_dmax       = 0.015_wp           !< [m] DSL max thickness (soil evap)
      real(wp)    :: dsl_theta_init = 0.8_wp             !< DSL initiation theta_1/theta_sat
      real(wp)    :: psi_wilt = -152.96_wp               !< [m] wilting head (~ -1.5 MPa) for the sink f_wilt
      real(wp)    :: psi_open = -3.37_wp                 !< [m] onset head (~ -0.033 MPa) for the f_wilt ramp
      logical     :: debug_error = .false.
   end type soil_opts_t

   !----- Soil-energy: solver selectors + tolerances (NOT the whole config). ----------------!
   type :: energy_opts_t
      integer(ik) :: soil_solver  = ENERGY_SOLVER_BE
      integer(ik) :: bottom_bc    = ENERGY_BC_GEOTHERMAL
      integer(ik) :: phase_change = ENERGY_PHASE_OFF
      integer(ik) :: substep      = ENERGY_SUBSTEP_ADAPTIVE
      real(wp)    :: rtol = 1.0e-3_wp, atol = 1.0e-2_wp     !< atol in [K]
      real(wp)    :: h_init = 900.0_wp
      integer(ik) :: max_substep = 200_ik
      logical     :: debug_error = .false.
   end type energy_opts_t

   !----- Snow / temporary-surface-water physical parameters (MEDS_SNOW_DESIGN.md P0). ------!
   !      Defaults for standalone/test use; the tunable subset is filled from the [snow] TOML  !
   !      block when the wiring lands (like soil_opts_t / energy_opts_t).                       !
   type :: snow_params_t
      real(wp) :: rho_snow          = 250.0_wp    !< [kg/m3] P0 fixed SEASONAL-mean bulk density (depth = swe/rho)
      real(wp) :: k_snow            = 0.15_wp     !< [W/m/K] P0 fixed snow thermal conductivity
      real(wp) :: albedo_vis_fresh  = 0.85_wp     !< [-] fresh-snow VIS albedo
      real(wp) :: albedo_nir_fresh  = 0.65_wp     !< [-] fresh-snow NIR albedo
      real(wp) :: albedo_vis_aged   = 0.518_wp    !< [-] aged VIS endpoint (ED2)
      real(wp) :: albedo_nir_aged   = 0.435_wp    !< [-] aged NIR endpoint (ED2)
      real(wp) :: snow_emiss        = 0.97_wp     !< [-] thermal (LW) emissivity
      real(wp) :: liquid_holding_frac = 0.10_wp   !< [-] percolation holding capacity (LEAF-3 1:9)
      real(wp) :: snow_stab_thresh  = 3.0_wp      !< [kg/m2] active-integration floor; thinner packs stay passive (MVP)
      real(wp) :: min_new_snow_mass = 1.0e-3_wp   !< [kg/m2] first-flake creation threshold (ED2 tiny_sfcwater_mass)
      real(wp) :: tiny_snow_mass    = 1.0e-6_wp   !< [kg/m2] melt-out threshold (dump residual to soil)
      real(wp) :: ny07_a            = 2.5_wp      !< [-] Niu-Yang07 snow-cover-fraction shape factor
      real(wp) :: ny07_m            = 1.0_wp      !< [-] Niu-Yang07 snow-cover-fraction exponent
      real(wp) :: z0_snow           = 0.001_wp    !< [m] snow surface roughness (aero burial floor)
   end type snow_params_t

   !----- Canopy-aerodynamics run constants (device-constant; defaults = ED2/CLM5 reference     !
   !      values). The tunable subset is filled from the [aerodynamics] TOML block when the       !
   !      wiring lands, like soil_opts_t. (The von Karman constant lives in meds_constants.)       !
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
