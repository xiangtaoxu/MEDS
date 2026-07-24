!==========================================================================================!
! meds_biophysics_types -- shared derived types of the biophysics domain (fast, sub-daily     !
! physical processes over the ecosystem column). Aggregates FOUR now-landed process families:  !
! canopy radiative transfer, soil-column hydrology, energy balance, and canopy aerodynamics.    !
!                                                                                          !
! Pure DATA (no methods, no hidden state), the RT analogue of meds_plant_types. RT block:      !
!   * rad_pft_optics_t -- the PRECOMPUTED, MU-INDEPENDENT per-PFT optics table (single-scatter !
!                         albedo omega and the leaf-angle scattering asymmetry g per band and   !
!                         tissue, plus the leaf-angle distribution lidf and its 2nd moment bf).   !
!                         Filled ONCE by meds_canopy_radiation%derive_rad_optics.                            !
!   * rad_forcing_t    -- the per-band incident fluxes (beam & diffuse), the solar cosine, and    !
!                         the ground reflectance/emission (from meds_canopy_radiation%ground_optics).          !
!                         Absolute W/m2 throughout -- no ED2-style normalize-to-one.                  !
!   * rad_flux_t       -- the returned per-cohort absorbed radiation (leaf & wood, per band) plus    !
!                         the patch-level albedo and below-canopy transmission diagnostics.           !
!                                                                                          !
! Bands are configurable; the default set is {VIS, NIR, LW}. Every band carries a thermal-      !
! emission term, identically zero for VIS/NIR (has_emission = .false.).                          !
!==========================================================================================!
module meds_biophysics_types
   use meds_kinds,             only : wp, ik
   use meds_therm_lib,            only : cas_enthalpy_of_temp
   !----- Soil constitutive curves live in meds_hydr_lib (retention family + SOIL_RETENTION_*) !
   !      and meds_therm_lib (thermal properties); the per-column PARAMETER types + their pure        !
   !      builders live in meds_column_state_types beside the prognostic soil columns. Re-exported !
   !      below so biophysics kernels + callers keep `use meds_biophysics_types, only : soil_params_t`.!
   use meds_column_state_types, only : n_soil_layer_max, cas_state_t, soil_column_t,           &
                                       soil_energy_column_t, snow_column_t,                     &
                                       soil_params_t, soil_thermal_params_t, soil_carbon_t
   use meds_hydr_lib,       only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   !----- Fast-loop run-config bundles (solver selectors + snow/aero parameters) live in the    !
   !      shared config layer (a low-level leaf, not the meds_config aggregator); re-exported    !
   !      below so callers keep `use meds_biophysics_types, only : soil_opts_t`. ----------------!
   use meds_biophysics_opts,    only : SOIL_SOLVER_BE, SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER,      &
                                       SOIL_BC_BEDROCK, SOIL_LIN_FROZEN, SOIL_LIN_PICARD,        &
                                       SOIL_SUBSTEP_ADAPTIVE, SOIL_SUBSTEP_FIXED,                &
                                       ENERGY_SOLVER_BE, ENERGY_BC_GEOTHERMAL, ENERGY_PHASE_OFF, &
                                       ENERGY_PHASE_ON, ENERGY_SUBSTEP_ADAPTIVE,                 &
                                       soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   implicit none
   private

   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t, surface_state_t
   public :: alloc_rad_pft_optics, alloc_rad_forcing, alloc_rad_flux

   !----- Soil-column hydrology (see meds_soil_water / meds_hydr_lib). ----!
   public :: n_soil_layer_max
   public :: SOIL_SOLVER_BE
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   public :: SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER, SOIL_BC_BEDROCK
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

   !---------------------------------------------------------------------------------------!
   ! Ground / surface optical state -- the two-stream lower boundary (bare-soil placeholder;   !
   ! only the soil fields are consulted now, the rest are reserved for the full surface model). !
   !---------------------------------------------------------------------------------------!
   type :: surface_state_t
      integer(ik)           :: n_band = 0_ik
      real(wp), allocatable :: soil_albedo(:)    !< (band) shortwave soil albedo; unused for emission bands
      real(wp)              :: soil_emiss = 0.96_wp   !< thermal emissivity of the ground
      real(wp)              :: soil_temp  = 298.0_wp  !< [K] ground (skin) temperature
   end type surface_state_t

   !=======================================================================================!
   !  Soil-column hydrology types + selector codes (meds_soil_water, design §4).  !
   !  Fixed-size (n_soil_layer_max) so the kernel stays allocatable-free and GPU-eligible.   !
   !  ED2 negative-z convention: elevation z <= 0 below ground; dz, dz_node are positive       !
   !  magnitudes.                                                                               !
   !=======================================================================================!
   !----- n_soil_layer_max + the prognostic column-state types live in meds_column_state_types    !
   !      (src/shared); the SOIL_* solver selectors + soil_opts_t live in meds_biophysics_opts     !
   !      (shared/config); the constitutive SOIL_RETENTION_* live in meds_hydr_lib. All re-exported !
   !      below so the fast kernels + callers keep `use meds_biophysics_types` unchanged. ----------!

   !----- Soil-column boundary conditions (read-only). ------------------------------------!
   type :: chydro_forcing_t
      real(wp) :: precip_ground = 0.0_wp                  !< [kg/m2/s] ground-reaching liquid (post interception)
      real(wp) :: root_uptake(n_soil_layer_max) = 0.0_wp  !< [kg/m2/s] per-layer transpiration DEMAND (x nplant)
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temp (FORCED = T_air until soil energy)
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance (soil evap series)
   end type chydro_forcing_t

   !----- soil_params_t (per-column geometry + texture) is defined in (and re-exported from)    !
   !      meds_column_state_types, beside the prognostic soil columns it describes. -------------!

   !----- soil_opts_t (soil-water solver selectors + tolerances) lives in meds_biophysics_opts    !
   !      (shared/config); re-exported above.                                                     !

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
   !  Energy-balance types + selector codes (feed meds_soil_energy / meds_vegetation_biophysics / meds_ground_biophysics / meds_cas_biophysics, !
   !  design 5). Prognostic INTERNAL ENERGY / enthalpy (phase-safe); temperature diagnosed.  !
   !  Reuse the negative-z n_soil_layer_max grid + meds_soil_solver Thomas sweep.             !
   !=======================================================================================!
   !----- The ENERGY_* solver selectors + energy_opts_t live in meds_biophysics_opts (shared/     !
   !      config); re-exported below.                                                             !

   public :: ENERGY_SOLVER_BE, ENERGY_BC_GEOTHERMAL
   public :: ENERGY_PHASE_OFF, ENERGY_PHASE_ON, ENERGY_SUBSTEP_ADAPTIVE
   public :: soil_energy_column_t, cas_state_t, soil_thermal_params_t, veg_thermal_params_t
   public :: energy_forcing_t, energy_flux_t
   public :: leaf_energy_env_t, leaf_energy_flux_t, energy_opts_t
   public :: snow_params_t, snow_env_t, snow_flux_t, snow_melt_t

   !----- (soil_energy_column_t + cas_state_t now live in meds_column_state_types; re-exported.) -!

   !----- soil_thermal_params_t (per-column soil thermal texture) is defined in (and re-exported !
   !      from) meds_column_state_types; the conductivity/heat-capacity kernels are in meds_therm_lib.!

   !----- Per-PFT vegetation thermal parameters. -------------------------------------------!
   type :: veg_thermal_params_t
      real(wp) :: leaf_emiss     = 0.95_wp                  !< [-] LW emissivity (Jacobian -8*eps*sigma*T^3 term)
      real(wp) :: effarea_heat   = 2.0_wp                   !< [-] sensible sidedness (both leaf sides)
      real(wp) :: effarea_evap   = 1.0_wp                   !< [-] film-evaporation sidedness
      real(wp) :: effarea_transp = 1.0_wp                   !< [-] transpiration sidedness (per PFT)
      real(wp) :: veg_hcap_min   = 20.0_wp                  !< [J/m2/K] resolvability floor
      real(wp) :: c_leaf = 3200.0_wp, c_sapw = 2700.0_wp    !< [J/kg/K] tissue specific heats
   end type veg_thermal_params_t

   !----- Soil-column thermal boundary conditions (read-only). ------------------------------!
   type :: energy_forcing_t
      real(wp) :: g_top      = 0.0_wp                       !< [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !< [W/m2] bottom flux (default 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !< [m3/m3] theta from hydrology (kappa, C_eff)
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !< [m/s]  inter-layer water flux for advective heat,
                                                            !<        UPWARD-positive (matches the hf face convention;
                                                            !<        the caller flips the DOWNWARD-positive hydrology
                                                            !<        flux -- see meds_fast_split.f90:
                                                            !<        eforc%w_flux = -hflux%w_flux). Populating it with
                                                            !<        the raw downward flux would reverse the advection.
      real(wp) :: root_heat_sink(n_soil_layer_max) = 0.0_wp !< [W/m2] enthalpy removed with root uptake
   end type energy_forcing_t

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

   !=======================================================================================!
   !  Canopy-air-space CO2 balance: the prognostic third twin is carried in cas_state_t and       !
   !  advanced by meds_cas_biophysics's cas_column_* kernels; the Rh selectors (HR_*) + co2_opts_t !
   !  / damm_params_t live in meds_biogeochem_types (heterotrophic respiration is decomposition).  !
   !=======================================================================================!

   !=======================================================================================!
   !  SNOW / temporary-surface-water types (meds_ground_biophysics; design MEDS_SNOW_DESIGN.md P0). !
   !  STATELESS per-store kernels; forcing arrives as value types. (snow_params_t -- the physical  !
   !  parameter table -- lives in meds_biophysics_opts, shared/config; re-exported above.)          !
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

   !----- Snow-mass outputs from accumulation + melt drainage (a store->store handoff to soil). -!
   type :: snow_melt_t
      real(wp) :: melt_mass = 0.0_wp                    !< [kg/m2] free liquid drained this step -> soil infiltration
      real(wp) :: melt_enth = 0.0_wp                    !< [J/m2]  matching enthalpy (internal_energy_liquid(t_3ple))
      real(wp) :: dump_mass = 0.0_wp                    !< [kg/m2] full-melt-out residual water dumped to the soil
      real(wp) :: dump_enth = 0.0_wp                    !< [J/m2]  full-melt-out residual energy dumped to the soil
      logical  :: melted_out = .false.                  !< .true. when the layer vanished (nlayer -> 0)
      real(wp) :: mass_resid = 0.0_wp                   !< [kg/m2] closed-budget residual (~0 by construction)
   end type snow_melt_t

   !=======================================================================================!
   !  Canopy-aerodynamics types (meds_canopy_aerodynamics; design Part I). STATELESS: the       !
   !  kernel is a pure function of (free-atm forcing, canopy-air-space state, canopy geometry).  !
   !  (aero_cfg_t -- the ED2/CLM literature run constants -- lives in meds_biophysics_opts,        !
   !  shared/config; re-exported above. The env/geom/out I/O types stay here.)                     !
   !=======================================================================================!
   public :: aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out, ensure_aero_out_capacity
   public :: patch_biophys_t, alloc_patch_biophys, ensure_patch_biophys_capacity

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
   !      contained MVP block used by meds_fast_split/meds_fast_ark; the eventual per-cohort/per-patch !
   !      state threaded through the demographic SoA lockstep reorder is the fast<->slow step.    !
   type :: patch_biophys_t
      type(cas_state_t)          :: cas               !< canopy-air-space twins (enthalpy/shv/co2)
      type(soil_energy_column_t) :: soil_e            !< soil thermal column (internal energy; temp diagnosed)
      type(soil_column_t)        :: soil_w            !< soil water column (theta; psi_soil diagnosed)
      type(snow_column_t)        :: snow              !< temporary-surface-water / snow store (swe + energy)
      !----- FROZEN slow soil-carbon pool (B2, MEDS_SLOW_DYNAMICS_DESIGN.md Part II): a read-only  !
      !      snapshot of site%patch%soil_carbon(ip), seeded ONCE at the top of the slow step and     !
      !      held constant across the day's fast sub-steps -- the fast loop's heterotrophic Rh        !
      !      respires against THIS frozen copy (never mutated here; the daily soil_carbon_step is     !
      !      the sole writer of the real site-level pool). Zero (soil_carbon_t's own default) when    !
      !      [soil_carbon].soil_carbon_on = .false., which reduces heterotrophic_respiration_matrix    !
      !      to Rh=0 -- so the OLD constant-pool scalar path is used instead in that case (gated in    !
      !      column_prepass on cfg%soil_carbon_on, not on this field being populated). ------------------!
      type(soil_carbon_t)        :: soil_carbon
      real(wp), allocatable      :: leaf_temp(:)      !< [K] per-cohort leaf temperature
      real(wp), allocatable      :: wood_temp(:)      !< [K] per-cohort wood/branch temperature (own store)
      !----- Internal (xylem/symplast) water mass [kg/plant] -- the prognostic hydraulic state;    !
      !      psi is diagnosed from it wherever needed (psi_from_water_content), never persisted.    !
      !      MEDS_ED2_RK45_DESIGN.md sec 4. -----------------------------------------------------!
      real(wp), allocatable      :: leaf_water_mass(:) !< [kg/plant] internal leaf water
      real(wp), allocatable      :: wood_water_mass(:) !< [kg/plant] internal wood water
      !----- Surface (interception film) water [kg/m2 ground] -- DISTINCT store from the internal      !
      !      water above; MEDS_ED2_RK45_DESIGN.md sec 3.4. Already ground-area-referenced (unlike the   !
      !      per-plant internal water), so fusion SUMS it (meds_core_cohort_fusefiss.f90). ------------!
      real(wp), allocatable      :: leaf_surf_water(:) !< [kg/m2 ground] leaf interception film
      real(wp), allocatable      :: wood_surf_water(:) !< [kg/m2 ground] wood interception film
      !----- Lagged per-layer root-uptake SHARES (sum = 1), from the previous fast step's multi-layer  !
      !       plant solve; the soil sink distributes coh_transp by these (vs static root_frac) so the   !
      !       soil dries where roots actually took water. Default 0 => root_frac fallback (first step /  !
      !       single-layer). Fixed-shape (no allocation); only used when [hydraulics].multilayer_roots. -!
      real(wp)                   :: root_sink_share(n_soil_layer_max) = 0.0_wp
      !----- WARM START for the adaptive march (MEDS_NUMERICS_SCOPING.md section 8e). The controller     !
      !      spends real work discovering the admissible step size, and that size is a property of the    !
      !      column's stiffness, which barely changes from one dt_fast to the next. Cold-starting each    !
      !      call at the full dt_fast threw that away and paid ~1 rejected step per call (measured        !
      !      1.65/1.08/0.39/0.014 rejections per call at dt_fast = 1800/900/450/225 s). Carrying the      !
      !      last controller proposal across calls is the standard ODE-solver warm restart. 0 = no        !
      !      history yet (first call / non-adaptive path) => cold start, i.e. the old behaviour. ---------!
      real(wp)                   :: adapt_dt_last = 0.0_wp   !< [s] last accepted controller step proposal
   end type patch_biophys_t

contains

   subroutine alloc_rad_pft_optics(optics, n_band, n_pft, n_class)
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik),            intent(in)  :: n_band, n_pft, n_class
      optics%n_band = n_band
      optics%n_pft  = n_pft
      allocate(optics%omega_leaf(n_band, n_pft), optics%omega_wood(n_band, n_pft))
      allocate(optics%g_leaf(n_band, n_pft),     optics%g_wood(n_band, n_pft))
      allocate(optics%clumping_leaf(n_pft),      optics%clumping_wood(n_pft))
      allocate(optics%lidf(n_class, n_pft),      optics%bf(n_pft))
      allocate(optics%has_beam(n_band),          optics%has_emission(n_band))
      optics%omega_leaf = 0.0_wp
      optics%omega_wood = 0.0_wp
      optics%g_leaf = 0.0_wp
      optics%g_wood = 0.0_wp
      optics%clumping_leaf = 1.0_wp
      optics%clumping_wood = 1.0_wp
      optics%lidf = 0.0_wp
      optics%bf = 0.0_wp
      optics%has_beam = .false.
      optics%has_emission = .false.
   end subroutine alloc_rad_pft_optics

   subroutine alloc_rad_forcing(f, n_band)
      type(rad_forcing_t), intent(out) :: f
      integer(ik),         intent(in)  :: n_band
      f%n_band = n_band
      allocate(f%incid_beam(n_band), f%incid_diff(n_band), f%grnd_refl(n_band), f%grnd_emiss(n_band))
      f%incid_beam = 0.0_wp
      f%incid_diff = 0.0_wp
      f%grnd_refl  = 0.0_wp
      f%grnd_emiss = 0.0_wp
   end subroutine alloc_rad_forcing

   subroutine alloc_rad_flux(flux, n_band, n_coh)
      type(rad_flux_t), intent(out) :: flux
      integer(ik),      intent(in)  :: n_band, n_coh
      flux%n_band = n_band
      flux%n_coh  = n_coh
      allocate(flux%abs_leaf(n_band, n_coh), flux%abs_wood(n_band, n_coh))
      allocate(flux%albedo(n_band), flux%dn_ground(n_band), flux%up_ground(n_band))
      flux%abs_leaf = 0.0_wp
      flux%abs_wood = 0.0_wp
      flux%albedo = 0.0_wp
      flux%dn_ground = 0.0_wp
      flux%up_ground = 0.0_wp
   end subroutine alloc_rad_flux

   !----- Allocate the per-cohort output arrays of an aero_out_t (the kernel writes in place). !
   subroutine alloc_aero_out(out, n_coh)
      type(aero_out_t), intent(out) :: out
      integer(ik),      intent(in)  :: n_coh
      out%n_coh = n_coh
      allocate(out%wind(n_coh), out%leaf_gbh(n_coh), out%leaf_gbw(n_coh),                      &
               out%wood_gbh(n_coh), out%wood_gbw(n_coh))
      out%wind = 0.0_wp
      out%leaf_gbh = 0.0_wp
      out%leaf_gbw = 0.0_wp
      out%wood_gbh = 0.0_wp
      out%wood_gbw = 0.0_wp
   end subroutine alloc_aero_out

   !----- Grow-only capacity check for aero_out_t (mirrors ensure_column_cohort_capacity,       !
   !      MEDS_NUMERICS_SCOPING.md BB1 phase 1): aero is a per-call scratch (canopy_aerodynamics  !
   !      overwrites every active element each call), so reusing over-sized capacity is safe.    !
   subroutine ensure_aero_out_capacity(out, n_coh)
      type(aero_out_t), intent(inout) :: out
      integer(ik),       intent(in)    :: n_coh
      if (.not. allocated(out%wind)) then
         call alloc_aero_out(out, n_coh)
      else if (size(out%wind) < n_coh) then
         call alloc_aero_out(out, n_coh)
      else
         out%n_coh = n_coh
      end if
   end subroutine ensure_aero_out_capacity

   !----- Allocate + seed a patch_biophys_t from an initial CAS temperature (mirrors the other !
   !      alloc_* helpers; seeds can_enthalpy via the shared thermo inverter). ----------------!
   subroutine alloc_patch_biophys(bio, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(patch_biophys_t), intent(out) :: bio
      integer(ik),           intent(in)  :: n_coh
      real(wp),              intent(in)  :: can_temp0, can_shv0, can_co2, leaf_temp0
      allocate(bio%leaf_temp(n_coh), bio%wood_temp(n_coh))
      allocate(bio%leaf_water_mass(n_coh), bio%wood_water_mass(n_coh))
      allocate(bio%leaf_surf_water(n_coh), bio%wood_surf_water(n_coh))
      bio%leaf_temp        = leaf_temp0
      bio%wood_temp        = leaf_temp0
      bio%leaf_water_mass  = 0.0_wp    ! scratch seed only -- always discarded by the next real gather
      bio%wood_water_mass  = 0.0_wp    ! (mirrors leaf_temp/wood_temp's own scratch-seed discipline)
      bio%leaf_surf_water  = 0.0_wp    ! ditto
      bio%wood_surf_water  = 0.0_wp
      bio%cas%can_temp     = can_temp0
      bio%cas%can_shv      = can_shv0
      bio%cas%can_co2      = can_co2
      bio%cas%can_enthalpy = cas_enthalpy_of_temp(can_temp0, can_shv0)
   end subroutine alloc_patch_biophys

   !---------------------------------------------------------------------------------------!
   ! Grow-only capacity check for the per-cohort arrays of patch_biophys_t (mirrors            !
   ! ensure_column_cohort_capacity, MEDS_NUMERICS_SCOPING.md BB1 phase 1). Does NOT touch        !
   ! bio%cas/soil_e/soil_w/snow/soil_carbon -- every caller overwrites those with the site's      !
   ! persisted per-patch reservoirs (site%patch%cas(ip) etc.) immediately after allocating, so    !
   ! their alloc_patch_biophys seed values are always discarded; only leaf_temp/wood_temp/         !
   ! leaf_water_mass/wood_water_mass/leaf_surf_water/wood_surf_water need their CAPACITY ensured     !
   ! here (the caller's gather loop fills indices 1..n_coh).                                          !
   !---------------------------------------------------------------------------------------!
   subroutine ensure_patch_biophys_capacity(bio, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(patch_biophys_t), intent(inout) :: bio
      integer(ik),            intent(in)    :: n_coh
      real(wp),               intent(in)    :: can_temp0, can_shv0, can_co2, leaf_temp0
      if (.not. allocated(bio%leaf_temp)) then
         call alloc_patch_biophys(bio, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      else if (size(bio%leaf_temp) < n_coh) then
         call alloc_patch_biophys(bio, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      end if
   end subroutine ensure_patch_biophys_capacity

end module meds_biophysics_types
