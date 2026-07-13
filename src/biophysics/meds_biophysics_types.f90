!==========================================================================================!
! meds_biophysics_types -- shared derived types of the biophysics domain (fast, sub-daily     !
! physical processes over the ecosystem column). Aggregates FOUR now-landed process families:  !
! canopy radiative transfer, soil-column hydrology, energy balance, and canopy aerodynamics.    !
!                                                                                          !
! Pure DATA (no methods, no hidden state), the RT analogue of meds_plant_types. RT block:      !
!   * rad_pft_optics_t -- the PRECOMPUTED, MU-INDEPENDENT per-PFT optics table (single-scatter !
!                         albedo omega and the leaf-angle scattering asymmetry g per band and   !
!                         tissue, plus the leaf-angle distribution lidf and its 2nd moment bf).   !
!                         Filled ONCE by meds_optics%derive_rad_optics.                            !
!   * rad_forcing_t    -- the per-band incident fluxes (beam & diffuse), the solar cosine, and    !
!                         the ground reflectance/emission (from meds_optics%ground_optics).          !
!                         Absolute W/m2 throughout -- no ED2-style normalize-to-one.                  !
!   * rad_flux_t       -- the returned per-cohort absorbed radiation (leaf & wood, per band) plus    !
!                         the patch-level albedo and below-canopy transmission diagnostics.           !
!                                                                                          !
! Bands are configurable; the default set is {VIS, NIR, LW}. Every band carries a thermal-      !
! emission term, identically zero for VIS/NIR (has_emission = .false.).                          !
!==========================================================================================!
module meds_biophysics_types
   use meds_kinds,             only : wp, ik
   use meds_thermo,            only : cas_enthalpy_of_temp
   use meds_column_state_types, only : n_soil_layer_max, cas_state_t, soil_column_t,           &
                                       soil_energy_column_t, N_HYDRO_NODE
   implicit none
   private

   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t
   public :: alloc_rad_pft_optics, alloc_rad_forcing, alloc_rad_flux

   !----- Soil-column hydrology (see meds_column_hydrology / meds_soil_parameters). ----!
   public :: n_soil_layer_max
   public :: SOIL_SOLVER_BE
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   public :: SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER, SOIL_BC_BEDROCK, SOIL_BC_SLOPE
   public :: SOIL_LIN_FROZEN, SOIL_LIN_PICARD
   public :: SOIL_SUBSTEP_ADAPTIVE, SOIL_SUBSTEP_FIXED
   public :: soil_column_t, chydro_forcing_t, soil_params_t, soil_opts_t, chydro_flux_t

   !----- Default three-band layout (indices into the band dimension). --------------------!
   integer(ik), parameter :: RAD_VIS = 1_ik   !< visible / PAR (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_NIR = 2_ik   !< near-infrared (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_LW  = 3_ik   !< thermal / longwave (diffuse + emission, no beam)
   integer(ik), parameter :: N_RAD_BAND_DEFAULT = 3_ik

   !---------------------------------------------------------------------------------------!
   ! Precomputed per-PFT optics. omega = rho + tau (single-scatter albedo); g = bf*(rho-tau)/  !
   ! omega is the SCOPE leaf-angle asymmetry, so the diffuse backscatter is beta = 0.5*(1+g)   !
   ! and the direct-beam upscatter is beta0 = 0.5*(1 + g/k), k = G(mu)/mu. Leaf and wood are   !
   ! stored separately and blended per cohort by clumping-corrected area. `lidf` (per PFT) is    !
   ! kept so the mu-dependent G(mu) can be evaluated each timestep.                              !
   !---------------------------------------------------------------------------------------!
   type :: rad_pft_optics_t
      integer(ik) :: n_pft  = 0_ik
      integer(ik) :: n_band = 0_ik
      real(wp), allocatable :: omega_leaf(:,:)     !< (band,pft) single-scatter albedo, leaves
      real(wp), allocatable :: omega_wood(:,:)     !< (band,pft) single-scatter albedo, wood
      real(wp), allocatable :: g_leaf(:,:)         !< (band,pft) bf*(rho-tau)/omega, leaves
      real(wp), allocatable :: g_wood(:,:)         !< (band,pft) bf*(rho-tau)/omega, wood
      real(wp), allocatable :: clumping_leaf(:)    !< (pft) leaf clumping factor (0,1]
      real(wp), allocatable :: clumping_wood(:)    !< (pft) wood clumping factor (0,1]
      real(wp), allocatable :: lidf(:,:)           !< (class,pft) leaf-angle distribution weights
      real(wp), allocatable :: bf(:)               !< (pft) <cos^2(theta_leaf)>
      logical,  allocatable :: has_beam(:)         !< (band) band has a collimated (solar) beam
      logical,  allocatable :: has_emission(:)     !< (band) band has thermal emission (LW)
   end type rad_pft_optics_t

   !---------------------------------------------------------------------------------------!
   ! Per-band incident forcing + ground boundary (absolute W/m2). For emission bands the      !
   ! incident diffuse is the atmospheric downwelling (rlong) and grnd_emiss is the surface      !
   ! thermal source; for shortwave bands grnd_emiss = 0 and grnd_refl is the albedo.            !
   !---------------------------------------------------------------------------------------!
   type :: rad_forcing_t
      integer(ik) :: n_band = 0_ik
      real(wp)    :: cosz   = 1.0_wp               !< cosine of solar zenith / incidence (floored > 0)
      real(wp), allocatable :: incid_beam(:)       !< (band) [W/m2] direct-beam incident at canopy top
      real(wp), allocatable :: incid_diff(:)       !< (band) [W/m2] diffuse incident at canopy top
      real(wp), allocatable :: grnd_refl(:)        !< (band) ground reflectance (albedo, or 1-emiss)
      real(wp), allocatable :: grnd_emiss(:)       !< (band) [W/m2] ground thermal emission (0 for SW)
   end type rad_forcing_t

   !---------------------------------------------------------------------------------------!
   ! Returned fluxes. Per cohort, per band: radiation absorbed by leaves and by wood [W/m2 of   !
   ! ground]. Patch level: albedo (upward/incident) and the below-canopy downwelling fluxes.     !
   !---------------------------------------------------------------------------------------!
   type :: rad_flux_t
      integer(ik) :: n_band = 0_ik, n_coh = 0_ik
      real(wp), allocatable :: abs_leaf(:,:)       !< (band,coh) [W/m2] absorbed by leaves
      real(wp), allocatable :: abs_wood(:,:)       !< (band,coh) [W/m2] absorbed by wood
      real(wp), allocatable :: albedo(:)           !< (band) canopy+ground albedo (SW) / upward frac
      real(wp), allocatable :: dn_ground(:)        !< (band) [W/m2] downwelling below canopy (to ground)
      real(wp), allocatable :: up_ground(:)        !< (band) [W/m2] upwelling from ground into canopy
   end type rad_flux_t

   !=======================================================================================!
   !  Soil-column hydrology types + selector codes (meds_column_hydrology, design §4).  !
   !  Fixed-size (n_soil_layer_max) so the kernel stays allocatable-free and GPU-eligible.   !
   !  ED2 negative-z convention: elevation z <= 0 below ground; dz, dz_node are positive       !
   !  magnitudes. The SOIL_* codes live here for the standalone P0/P1 cut and migrate to        !
   !  meds_config when the TOML string->enum mapping lands at P3.                                !
   !=======================================================================================!
   !----- n_soil_layer_max + the prognostic column-state types now live in meds_column_state_types !
   !      (src/shared) so the state hub can own them; re-exported below for the fast kernels. -------!

   integer(ik), parameter :: SOIL_SOLVER_BE = 1_ik         !< linearly-implicit backward-Euler (only solver)

   integer(ik), parameter :: SOIL_RETENTION_VG       = 1_ik  !< van Genuchten-Mualem (default)
   integer(ik), parameter :: SOIL_RETENTION_CAMPBELL = 2_ik  !< Campbell / Clapp-Hornberger (option)

   integer(ik), parameter :: SOIL_BC_FREE_DRAIN = 1_ik     !< unit-gradient q = K(theta_N)  (MVP default)
   integer(ik), parameter :: SOIL_BC_AQUIFER    = 2_ik     !< SIMTOP aquifer / water table (P2)
   integer(ik), parameter :: SOIL_BC_BEDROCK    = 3_ik     !< no-flow  q = 0
   integer(ik), parameter :: SOIL_BC_SLOPE      = 4_ik     !< slope-reduced drainage (P2)

   integer(ik), parameter :: SOIL_LIN_FROZEN = 1_ik        !< frozen-coefficient single solve (MVP)
   integer(ik), parameter :: SOIL_LIN_PICARD = 2_ik        !< Celia modified-Picard (P2)

   integer(ik), parameter :: SOIL_SUBSTEP_ADAPTIVE = 1_ik  !< adaptive step-doubling
   integer(ik), parameter :: SOIL_SUBSTEP_FIXED    = 2_ik  !< fixed count (GPU warp-uniform)

   !----- Soil-column boundary conditions (read-only). ------------------------------------!
   type :: chydro_forcing_t
      real(wp) :: precip_ground = 0.0_wp                  !< [kg/m2/s] ground-reaching liquid (post interception)
      real(wp) :: root_uptake(n_soil_layer_max) = 0.0_wp  !< [kg/m2/s] per-layer transpiration DEMAND (x nplant)
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temp (FORCED = T_air until soil energy)
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance (soil evap series)
   end type chydro_forcing_t

   !----- Per-column geometry + texture (assembled once per site; ED2 negative-z). ---------!
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
      real(wp) :: vg_alpha(n_soil_layer_max)  = 0.0_wp       !< [1/m] van Genuchten inverse air-entry
      real(wp) :: vg_n(n_soil_layer_max)      = 0.0_wp       !< [-] van Genuchten pore-size index (> 1)
      real(wp) :: psi_sat(n_soil_layer_max)   = 0.0_wp       !< [m] Campbell air-entry potential (option)
      real(wp) :: b_camp(n_soil_layer_max)    = 0.0_wp       !< [-] Campbell exponent (option)
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
      real(wp) :: theta_wp(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] wilting point (DERIVED)
      real(wp) :: root_frac(n_soil_layer_max) = 0.0_wp       !< [-] normalized root fraction (sum = 1)
   end type soil_params_t

   !----- Pre-extracted solver selectors + tolerances (NOT the whole config). --------------!
   type :: soil_opts_t
      integer(ik) :: solver    = SOIL_SOLVER_BE
      integer(ik) :: bottom_bc = SOIL_BC_FREE_DRAIN
      integer(ik) :: linearize = SOIL_LIN_FROZEN
      integer(ik) :: substep   = SOIL_SUBSTEP_ADAPTIVE
      logical     :: zeng_decker = .false.               !< P2 only; MVP plain-gravity flux
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
      real(wp)    :: f_drai      = 2.5_wp                !< [1/m] SIMTOP baseflow decay (SOIL_BC_AQUIFER, P2)
      real(wp)    :: q_drai_max  = 5.5e-6_wp             !< [kg/m2/s] max baseflow (SOIL_BC_AQUIFER, P2)
      real(wp)    :: f_over      = 0.5_wp                !< [1/m] saturated-fraction decay (Dunne runoff, P2)
      real(wp)    :: f_max       = 0.4_wp                !< [-] max saturated fraction (Dunne runoff, P2)
      logical     :: debug_error = .false.
   end type soil_opts_t

   !----- Soil-column outputs + diagnostics. ----------------------------------------------!
   type :: chydro_flux_t
      real(wp) :: infiltration   = 0.0_wp                !< [kg/m2/s] top-face infiltration
      real(wp) :: drainage       = 0.0_wp                !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf    = 0.0_wp                !< [kg/m2/s] surface runoff
      real(wp) :: soil_evap      = 0.0_wp                !< [kg/m2/s] ground evaporation
      real(wp) :: uptake_total   = 0.0_wp                !< [kg/m2/s] realized root uptake (after theta_wp cap)
      real(wp) :: uptake_deficit = 0.0_wp                !< [kg/m2/s] capped (unmet) sink
      real(wp) :: clip_excess    = 0.0_wp                !< [kg/m2/s] theta-clip water routed to ponding
      real(wp) :: psi_soil(n_soil_layer_max) = 0.0_wp    !< [MPa] per-layer matric potential (EXPORTED to hydraulics)
      real(wp) :: w_flux(n_soil_layer_max)   = 0.0_wp    !< [m/s] time-mean DOWNWARD Darcy flux BELOW node k (k=1..n-1);
                                                         !<       interior interfaces only (EXPORTED for advective heat)
      real(wp) :: mass_resid     = 0.0_wp                !< [kg/m2] closed-budget residual (~0)
      integer(ik) :: nsub = 0_ik                         !< sub-steps taken
      logical  :: converged = .true.                     !< .false. on any cap-hit
   end type chydro_flux_t

   !=======================================================================================!
   !  Energy-balance types + selector codes (meds_column_energy -- whole soil-veg-air column, !
   !  design 5). Prognostic INTERNAL ENERGY / enthalpy (phase-safe); temperature diagnosed.  !
   !  Reuse the negative-z n_soil_layer_max grid + meds_soil_solver Thomas sweep.             !
   !=======================================================================================!
   integer(ik), parameter :: ENERGY_SOLVER_BE       = 1_ik   !< implicit backward-Euler (only solver)
   integer(ik), parameter :: ENERGY_BC_GEOTHERMAL   = 1_ik   !< bottom: zero/geothermal flux
   integer(ik), parameter :: ENERGY_BC_PRESCRIBED_T = 2_ik   !< bottom: prescribed deep temperature
   integer(ik), parameter :: ENERGY_PHASE_OFF = 0_ik, ENERGY_PHASE_ON = 1_ik     !< freeze/thaw plateau (P1 off)
   integer(ik), parameter :: ENERGY_SUBSTEP_ADAPTIVE = 1_ik, ENERGY_SUBSTEP_FIXED = 2_ik

   public :: ENERGY_SOLVER_BE, ENERGY_BC_GEOTHERMAL, ENERGY_BC_PRESCRIBED_T
   public :: ENERGY_PHASE_OFF, ENERGY_PHASE_ON, ENERGY_SUBSTEP_ADAPTIVE, ENERGY_SUBSTEP_FIXED
   public :: soil_energy_column_t, cas_state_t, soil_thermal_params_t, veg_thermal_params_t
   public :: energy_forcing_t, cas_atm_forcing_t, energy_flux_t
   public :: leaf_energy_env_t, leaf_energy_flux_t, energy_opts_t

   !----- (soil_energy_column_t + cas_state_t now live in meds_column_state_types; re-exported.) -!

   !----- Per-column soil THERMAL texture (geometry+porosity come via soil_params_t). -------!
   type :: soil_thermal_params_t
      integer(ik) :: nzg_active = n_soil_layer_max
      real(wp) :: soil_solid_conductivity(n_soil_layer_max) = 0.0_wp   !< [W/m/K] kappa_solid
      real(wp) :: soil_dry_conductivity(n_soil_layer_max)   = 0.0_wp   !< [W/m/K] kappa_dry
      real(wp) :: soil_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   !< [J/m3/K] dry-matrix vol. heat cap
   end type soil_thermal_params_t

   !----- Per-PFT vegetation thermal parameters. -------------------------------------------!
   type :: veg_thermal_params_t
      real(wp) :: leaf_emiss     = 0.95_wp                  !< [-] LW emissivity (Jacobian -8*eps*sigma*T^3 term)
      real(wp) :: effarea_heat   = 2.0_wp                   !< [-] sensible sidedness (both leaf sides)
      real(wp) :: effarea_evap   = 1.0_wp                   !< [-] film-evaporation sidedness
      real(wp) :: effarea_transp = 1.0_wp                   !< [-] transpiration sidedness (per PFT)
      real(wp) :: veg_hcap_min   = 20.0_wp                  !< [J/m2/K] resolvability floor
      real(wp) :: c_leaf = 3200.0_wp, c_sapw = 2700.0_wp    !< [J/kg/K] tissue specific heats
      real(wp) :: c_dead = 2300.0_wp, c_bark = 2000.0_wp
   end type veg_thermal_params_t

   !----- Soil-column thermal boundary conditions (read-only). ------------------------------!
   type :: energy_forcing_t
      real(wp) :: g_top      = 0.0_wp                       !< [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !< [W/m2] bottom flux (default 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !< [m3/m3] theta from hydrology (kappa, C_eff)
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !< [m/s]  inter-layer water flux for advective heat,
                                                            !<        UPWARD-positive (matches the hf face convention;
                                                            !<        the caller flips the DOWNWARD-positive hydrology
                                                            !<        flux -- see meds_column_dynamics.f90:
                                                            !<        eforc%w_flux = -hflux%w_flux). Populating it with
                                                            !<        the raw downward flux would reverse the advection.
      real(wp) :: root_heat_sink(n_soil_layer_max) = 0.0_wp !< [W/m2] enthalpy removed with root uptake
   end type energy_forcing_t

   !----- Atmospheric forcing feeding the canopy air space (read-only). ---------------------!
   type :: cas_atm_forcing_t
      real(wp) :: ustar        = 0.0_wp                     !< [m/s]  friction velocity
      real(wp) :: enthalpy_atm = 0.0_wp                     !< [J/kg] reference-level specific enthalpy
      real(wp) :: w_flux_ac    = 0.0_wp                     !< [kg/m2/s] atm<->CAS water-vapour mass flux
      real(wp) :: co2_atm      = 400.0_wp                   !< [umol/mol] reference-level (free-atmosphere) CO2
      real(wp) :: rho_air      = 1.2_wp                     !< [kg/m3] air density
   end type cas_atm_forcing_t

   !----- Soil-column energy outputs + diagnostics. ----------------------------------------!
   type :: energy_flux_t
      real(wp) :: ground_heat  = 0.0_wp                     !< [W/m2] conductive flux into layer 1
      real(wp) :: bottom_heat  = 0.0_wp                     !< [W/m2] advective+geothermal bottom loss
      real(wp) :: energy_resid = 0.0_wp                     !< [J/m2] closed-budget residual (~0)
      integer(ik) :: nsub = 0_ik
      logical  :: converged = .true.
   end type energy_flux_t

   !----- Per-cohort surface BCs (leaf OR wood). -------------------------------------------!
   type :: leaf_energy_env_t
      real(wp) :: abs_sw = 0.0_wp, abs_lw = 0.0_wp          !< [W/m2] absorbed SW, NET LW (from canopy RT)
      real(wp) :: can_temp = 298.15_wp, can_shv = 0.0_wp    !< [K],[kg/kg] CAS state (FORCED sibling)
      real(wp) :: gbh = 0.0_wp, gbw = 0.0_wp                !< [m/s] boundary-layer heat/vapour conductance
      real(wp) :: gsw = 0.0_wp, fs_open = 1.0_wp            !< [m/s] stomatal (leaf only), open fraction
      real(wp) :: area_index = 0.0_wp                       !< [m2/m2] LAI (leaf) or WAI (wood)
      real(wp) :: leaf_water = 0.0_wp, wmass = 0.0_wp       !< [kg/m2] film (sigma_w) and total water (heat cap)
      real(wp) :: dry_hcap = 0.0_wp                         !< [J/m2/K] tissue heat capacity
      real(wp) :: rho_air = 1.2_wp, press = 101325.0_wp     !< [kg/m3],[Pa] CAS air
   end type leaf_energy_env_t

   !----- Per-cohort energy outputs. -------------------------------------------------------!
   type :: leaf_energy_flux_t
      real(wp) :: temp = 298.15_wp, fliq = 1.0_wp           !< [K],[-] diagnosed store temperature
      real(wp) :: h_flux = 0.0_wp, qw_flux = 0.0_wp, q_transp = 0.0_wp   !< [W/m2] sensible, film-evap, transp
      real(wp) :: w_flux = 0.0_wp, transp = 0.0_wp          !< [kg/m2/s] mass twins (-> CAS can_shv)
      real(wp) :: energy_resid = 0.0_wp                     !< [J/m2] closed-budget residual (~0)
   end type leaf_energy_flux_t

   !----- Solver selectors + tolerances (NOT the whole config). -----------------------------!
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

   !=======================================================================================!
   !  Canopy-aerodynamics types (meds_canopy_aerodynamics; design Part I). STATELESS: the      !
   !  kernel is a pure function of (free-atm forcing, canopy-air-space state, canopy geometry). !
   !  Config carries the ED2/CLM literature constants (defaults here for standalone/test use;    !
   !  the tunable subset is filled from the [aerodynamics] TOML block at the aux layer, like       !
   !  soil_opts_t). The CLM Monin-Obukhov surface layer + ED2 Nusselt boundary layers + ED2         !
   !  per-cohort wind extinction + CLM-style ground conductance.                                    !
   !=======================================================================================!
   public :: aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   public :: patch_biophys_t, alloc_patch_biophys

   !----- Run constants (device-constant; defaults = ED2/CLM5 reference values). -------------!
   type :: aero_cfg_t
      real(wp) :: vonk = 0.4_wp                     !< von Karman
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

   !----- Per-patch forcing + canopy-air-space state (read-only). ---------------------------!
   type :: aero_env_t
      real(wp) :: u_ref     = 2.0_wp                !< [m/s]      wind at reference height
      real(wp) :: zref      = 30.0_wp               !< [m]        reference (measurement) height
      real(wp) :: theta_atm = 298.15_wp             !< [K]        potential temp at zref
      real(wp) :: shv_atm   = 0.010_wp              !< [kg/kg]    specific humidity at zref
      real(wp) :: co2_atm   = 400.0_wp              !< [umol/mol] free-atmosphere CO2
      real(wp) :: press     = 101325.0_wp           !< [Pa]
      real(wp) :: rho_air   = 1.2_wp                !< [kg/m3]    (diagnostic passthrough)
      real(wp) :: can_theta = 298.15_wp             !< [K]        CAS potential temp
      real(wp) :: can_temp  = 298.15_wp             !< [K]        CAS actual temp (buoyancy Grashof)
      real(wp) :: can_shv   = 0.010_wp              !< [kg/kg]    CAS specific humidity
      real(wp) :: can_co2   = 400.0_wp              !< [umol/mol] CAS CO2
      real(wp) :: t_ground  = 298.15_wp             !< [K]        ground skin temp (ground-conductance stability)
   end type aero_env_t

   !----- Per-patch canopy geometry. -------------------------------------------------------!
   type :: aero_geom_t
      real(wp) :: veg_height   = 20.0_wp            !< [m]  canopy top height
      real(wp) :: opencan_frac = 0.0_wp             !< [-]  open-sky fraction (0 = closed canopy)
      real(wp) :: snowfac      = 0.0_wp             !< [-]  snow burial fraction of the canopy
   end type aero_geom_t

   !----- Outputs: per-patch scalars + per-cohort arrays (caller-owned; written in place). ---!
   type :: aero_out_t
      integer(ik) :: n_coh = 0_ik
      real(wp) :: ustar = 0.0_wp, tstar = 0.0_wp, qstar = 0.0_wp, cstar = 0.0_wp   !< [m/s],[K],[kg/kg],[umol/mol]
      real(wp) :: temp1 = 0.0_wp, temp2 = 0.0_wp    !< scalar profile factors (gah=rho*ustar*temp1, gaw=..temp2)
      real(wp) :: zeta = 0.0_wp, rib = 0.0_wp, obu = 0.0_wp   !< stability diagnostics
      real(wp) :: ggbare = 0.0_wp, ggveg = 0.0_wp, ggnet = 0.0_wp   !< [m/s] ground conductances (r_aero = 1/ggnet)
      real(wp) :: rough = 0.0_wp, displace = 0.0_wp, can_depth = 0.0_wp   !< [m]
      real(wp) :: uh = 0.0_wp                        !< [m/s] canopy-top wind (diagnostic)
      real(wp), allocatable :: wind(:)               !< [m/s] per-cohort in-canopy wind
      real(wp), allocatable :: leaf_gbh(:), leaf_gbw(:)   !< [m/s] leaf boundary-layer heat/vapour conductance
      real(wp), allocatable :: wood_gbh(:), wood_gbw(:)   !< [m/s] wood boundary-layer heat/vapour conductance
   end type aero_out_t

   !----- Per-patch fast biophysics STATE (prognostic; carried between fast steps). The self- -!
   !      contained MVP block used by meds_column_dynamics; the eventual per-cohort/per-patch    !
   !      state threaded through the demographic SoA lockstep reorder is the fast<->slow step.    !
   type :: patch_biophys_t
      type(cas_state_t)          :: cas               !< canopy-air-space twins (enthalpy/shv/co2)
      type(soil_energy_column_t) :: soil_e            !< soil thermal column (internal energy; temp diagnosed)
      type(soil_column_t)        :: soil_w            !< soil water column (theta; psi_soil diagnosed)
      real(wp), allocatable      :: leaf_temp(:)      !< [K] per-cohort leaf temperature
      real(wp), allocatable      :: wood_temp(:)      !< [K] per-cohort wood/branch temperature (own store)
      real(wp), allocatable      :: psi(:,:)          !< [MPa] plant water potential (N_HYDRO=3 nodes, cohort)
   end type patch_biophys_t

contains

   subroutine alloc_rad_pft_optics(opt, n_band, n_pft, n_class)
      type(rad_pft_optics_t), intent(out) :: opt
      integer(ik),            intent(in)  :: n_band, n_pft, n_class
      opt%n_band = n_band ; opt%n_pft = n_pft
      allocate(opt%omega_leaf(n_band, n_pft), opt%omega_wood(n_band, n_pft))
      allocate(opt%g_leaf(n_band, n_pft),     opt%g_wood(n_band, n_pft))
      allocate(opt%clumping_leaf(n_pft),      opt%clumping_wood(n_pft))
      allocate(opt%lidf(n_class, n_pft),      opt%bf(n_pft))
      allocate(opt%has_beam(n_band),          opt%has_emission(n_band))
      opt%omega_leaf = 0.0_wp ; opt%omega_wood = 0.0_wp
      opt%g_leaf = 0.0_wp ; opt%g_wood = 0.0_wp
      opt%clumping_leaf = 1.0_wp ; opt%clumping_wood = 1.0_wp
      opt%lidf = 0.0_wp ; opt%bf = 0.0_wp
      opt%has_beam = .false. ; opt%has_emission = .false.
   end subroutine alloc_rad_pft_optics

   subroutine alloc_rad_forcing(f, n_band)
      type(rad_forcing_t), intent(out) :: f
      integer(ik),         intent(in)  :: n_band
      f%n_band = n_band
      allocate(f%incid_beam(n_band), f%incid_diff(n_band), f%grnd_refl(n_band), f%grnd_emiss(n_band))
      f%incid_beam = 0.0_wp ; f%incid_diff = 0.0_wp ; f%grnd_refl = 0.0_wp ; f%grnd_emiss = 0.0_wp
   end subroutine alloc_rad_forcing

   subroutine alloc_rad_flux(flux, n_band, n_coh)
      type(rad_flux_t), intent(out) :: flux
      integer(ik),      intent(in)  :: n_band, n_coh
      flux%n_band = n_band ; flux%n_coh = n_coh
      allocate(flux%abs_leaf(n_band, n_coh), flux%abs_wood(n_band, n_coh))
      allocate(flux%albedo(n_band), flux%dn_ground(n_band), flux%up_ground(n_band))
      flux%abs_leaf = 0.0_wp ; flux%abs_wood = 0.0_wp
      flux%albedo = 0.0_wp ; flux%dn_ground = 0.0_wp ; flux%up_ground = 0.0_wp
   end subroutine alloc_rad_flux

   !----- Allocate the per-cohort output arrays of an aero_out_t (the kernel writes in place). !
   subroutine alloc_aero_out(out, n_coh)
      type(aero_out_t), intent(out) :: out
      integer(ik),      intent(in)  :: n_coh
      out%n_coh = n_coh
      allocate(out%wind(n_coh), out%leaf_gbh(n_coh), out%leaf_gbw(n_coh),                      &
               out%wood_gbh(n_coh), out%wood_gbw(n_coh))
      out%wind = 0.0_wp
      out%leaf_gbh = 0.0_wp ; out%leaf_gbw = 0.0_wp
      out%wood_gbh = 0.0_wp ; out%wood_gbw = 0.0_wp
   end subroutine alloc_aero_out

   !----- Allocate + seed a patch_biophys_t from an initial CAS temperature (mirrors the other !
   !      alloc_* helpers; seeds can_enthalpy via the shared thermo inverter). ----------------!
   subroutine alloc_patch_biophys(bio, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(patch_biophys_t), intent(out) :: bio
      integer(ik),           intent(in)  :: n_coh
      real(wp),              intent(in)  :: can_temp0, can_shv0, can_co2, leaf_temp0
      allocate(bio%leaf_temp(n_coh), bio%wood_temp(n_coh), bio%psi(N_HYDRO_NODE, n_coh))  ! first dim = N_HYDRO_NODE (leaf/wood/root)
      bio%leaf_temp        = leaf_temp0
      bio%wood_temp        = leaf_temp0
      bio%psi              = -0.1_wp                            ! mild initial tension; hydraulics relaxes it
      bio%cas%can_temp     = can_temp0
      bio%cas%can_shv      = can_shv0
      bio%cas%can_co2      = can_co2
      bio%cas%can_enthalpy = cas_enthalpy_of_temp(can_temp0, can_shv0)
   end subroutine alloc_patch_biophys

end module meds_biophysics_types
