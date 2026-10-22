! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_types -- the shared types of the fast loop: the working bundles the dispatch   !
! hands to either integrator (the run and site parts of the configuration, column_config_t !
! and column_site_t; column_cohort_t, column_forcing_t, column_budget_t, patch_biophys_t)  !
! and the column state, frozen record, tendency and boundary-flux types the integrators    !
! advance (column_state_t, column_frozen_t, column_tend_t, stage_bflux_t, column_bflux_t). !
! A leaf module, so meds_fast_ark, meds_fast_rk45 and the shared right-hand side import it !
! rather than each other (docs/science/numerical_scheme.md §2-§3,                          !
! docs/science/column_biophysics.md).                                                      !
!==========================================================================================!
module meds_fast_types
   use meds_kinds, only : wp, ik
   use meds_constants, only : t_3ple, n_soil_layer_max
   use meds_plant_types, only : hydro_args_table_t, leaf_photo_table_t
   use meds_soil_params, only : soil_params_t, soil_thermal_params_t, soil_opts_t, energy_opts_t
   use meds_patch_state_vars,   only : cas_state_t, soil_water_column_t, soil_energy_column_t, snow_column_t, soil_carbon_t
   use meds_therm_lib, only : cas_temp_to_enthalpy
   use meds_biophysics_params, only : aero_cfg_t, snow_params_t
   use meds_biogeochem_types, only : N_SOIL_POOL
   use meds_biogeochem_params, only : decomp_opts_t
   use meds_budget_check, only : budget_t
   use meds_spac_solver, only : spac_solver_t
   use meds_fast_params, only : INTEG_ARK, CTRL_I
   use meds_state_defaults,    only : UNSET_NEGATIVE
   implicit none
   private

   public :: column_config_t, column_site_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W, GRP_THETA, N_TOL_GROUP
   public :: tol_set_t, error_control_t, integrator_opts_t
   public :: ensure_column_cohort_capacity
   public :: surface_state_t, surface_tend_t
   public :: patch_biophys_t, ensure_patch_biophys_capacity
   public :: snow_stage_t
   public :: cas_boundary_t, tissue_coefficients_t, canopy_film_capacity_t, ground_boundary_t
   public :: soil_hydrology_t, root_zone_t, plant_water_t, ark_stages_t, fast_workspace_t
   public :: column_state_t, column_frozen_t, column_tend_t
   public :: stage_bflux_t, column_bflux_t
   public :: alloc_column_cohort, alloc_patch_biophys   ! used only by the tests

   !----- The tolerance groups, one per kind of state, each scaling its entries in the error norm    !
   !      (docs/science/numerical_scheme.md eq. (3)). GRP_THETA also sets the soil-plant water       !
   !      solver's own sub-step control, so one tolerance set governs both. -------------------------!
   integer(ik), parameter :: GRP_ENTH    = 1_ik   !< CAS specific enthalpy   [J/kg]
   integer(ik), parameter :: GRP_SHV     = 2_ik   !< CAS specific humidity   [kg/kg]
   integer(ik), parameter :: GRP_CO2     = 3_ik   !< CAS CO2 mole fraction   [umol/mol]
   integer(ik), parameter :: GRP_SE      = 4_ik   !< soil internal energy    [J/m3]
   integer(ik), parameter :: GRP_LEAF_W  = 5_ik   !< leaf internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_WOOD_W  = 6_ik   !< wood internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_THETA   = 7_ik   !< soil moisture           [m3/m3] (soil-water sub-solver)
   !----- No soil-temperature group: the state carries soil energy (GRP_SE). ----------------------!
   integer(ik), parameter :: N_TOL_GROUP = 7_ik

   !----- Per-field default absolute tolerances: the group defaults unless config overrides them. ---!
   real(wp), parameter :: ATOL_ENTH_DEF   = 5.0e1_wp    !< [J/kg]      (~0.05 K in enthalpy)
   real(wp), parameter :: ATOL_SHV_DEF    = 1.0e-6_wp   !< [kg/kg]
   real(wp), parameter :: ATOL_CO2_DEF    = 1.0e-1_wp   !< [umol/mol]
   real(wp), parameter :: ATOL_SE_DEF     = 1.0e3_wp    !< [J/m3]
   real(wp), parameter :: ATOL_LEAF_W_DEF = 1.0e-4_wp   !< [kg/plant]
   real(wp), parameter :: ATOL_WOOD_W_DEF = 1.0e-4_wp   !< [kg/plant]
   real(wp), parameter :: ATOL_THETA_DEF  = 1.0e-4_wp   !< [m3/m3] (== soil_opts_t's own default)
   !----- Default PI gains for a first-order embedded pair (Gustafsson 1988; Soderlind), eq. (4). --!
   real(wp), parameter :: PI_ALPHA_DEF = 0.35_wp
   real(wp), parameter :: PI_BETA_DEF  = 0.20_wp

   !----- Per-group (rtol, atol). The WRMS normalizes state group g by atol(g) + rtol(g)*|y|. ---------!
   type :: tol_set_t
      real(wp) :: rtol(N_TOL_GROUP) = 1.0e-3_wp
      real(wp) :: atol(N_TOL_GROUP) = [ATOL_ENTH_DEF, ATOL_SHV_DEF, ATOL_CO2_DEF, ATOL_SE_DEF,   &
                                       ATOL_LEAF_W_DEF, ATOL_WOOD_W_DEF, ATOL_THETA_DEF]
   end type tol_set_t

   !----- What an adaptive march runs under: controller, step-factor limits, tolerances. ----------!
   type :: error_control_t
      integer(ik)      :: controller = CTRL_I
      real(wp)         :: safety     = 0.9_wp
      real(wp)         :: fmin       = 0.2_wp
      real(wp)         :: fmax       = 5.0_wp
      real(wp)         :: pi_alpha   = PI_ALPHA_DEF
      real(wp)         :: pi_beta    = PI_BETA_DEF
      !----- p of docs/science/numerical_scheme.md eq. (4), the embedded estimate's order: 1 ARK, 4 RK45. !
      integer(ik)      :: p_order    = 1_ik
      type(tol_set_t)  :: tols
   end type error_control_t

   !----- The integrator's configuration, built once per run (meds_fast_config%build_integrator_opts) !
   !      and carried on column_config_t. ---------------------------------------------------------!
   type :: integrator_opts_t
      integer(ik) :: scheme           = INTEG_ARK  !< INTEG_ARK | INTEG_RK45
      logical     :: adaptive         = .true.     !< ARK: adaptive sub-stepping (else fixed_substeps)
      integer(ik) :: fixed_substeps   = 1_ik       !< ARK, adaptive = .false.: equal sub-steps per dt_fast
      type(error_control_t) :: error_control       !< tolerances + controller
   end type integrator_opts_t

   !----- The run part of the column's configuration, built once per run (build_fast_context) and  !
   !      shared by every site: what the run's configuration and the PFT table set. ---------------!
   type :: column_config_t
      type(aero_cfg_t)            :: aero            !< aerodynamics constants
      !----- Stem and fine-root maintenance respiration's peaked temperature response: leaf Rd's,     !
      !      [leaf_physiology] ea_rd, hd_rd, ds_rd (docs/science/plant_respiration.md §4). ----------!
      real(wp)                    :: ea_rd, hd_rd, ds_rd
      !----- Per-PFT plant hydraulics, built by apply_hydraulics_config. There is no PFT-uniform copy,
      !      so a caller cannot run every PFT on the first PFT's hydraulics by mistake.
      type(hydro_args_table_t)    :: hydraulics_table   !< per-PFT PV curves, vulnerability, conductance
      type(integrator_opts_t)     :: integrator     !< the fast-loop integrator's configuration (built once per run)
      !----- Canopy-surface water, the interception film and its evaporation and dew                 !
      !      (docs/science/soil_biophysics.md §Canopy interception); opt-in. RK45 integrates the     !
      !      films in its stages, ARK advances them once per step. Off, the wetted fraction stays 0  !
      !      and every film term vanishes. ----------------------------------------------------------!
      logical                     :: canopy_water_on  = .false.
      !----- The films' capacities per PFT [kg/m2 leaf], [kg/m2 wood]: a cohort holds              !
      !      leaf_max*LAI + wood_max*WAI. Built by apply_canopy_film_config. ----------------------!
      real(wp), allocatable       :: leaf_surf_water_max(:), wood_surf_water_max(:)
      type(snow_params_t) :: snow                    !< snow parameters (density, albedo, thresholds, conductivity)
      type(decomp_opts_t) :: soil_carbon             !< [soil_carbon]: the switch and the CENTURY parameters (Rh)
      !----- The ARK's Newton iteration cap is NEWT_MAX in meds_fast_be_stage, not a config field. ---!
   end type column_config_t

   !----- The site part: what differs from site to site, made for each site every slow step       !
   !      (site_column_configs): its soil column, the water and heat solvers' options with its     !
   !      bottom boundary, and its leaf table, acclimated to its own air. ------------------------!
   type :: column_site_t
      type(soil_params_t)         :: soil            !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal    !< soil thermal texture
      type(soil_opts_t)           :: soil_water_opts !< soil water and water-network options, the bottom boundary
      type(energy_opts_t)         :: energy          !< soil-thermal solver options, the deep temperature
      type(leaf_photo_table_t)    :: leaf_photo      !< per-PFT leaf-photosynthesis parameters
   end type column_site_t

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into the PFT table)
      real(wp),    allocatable :: lai(:), wai(:), height(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: aboveground_frac(:)          !< [--] gathered per-PFT (stem respiration)
      logical,     allocatable :: is_woody(:)                   !< gathered per-PFT (stem respiration off for grass)
      real(wp),    allocatable :: stem_resp_factor25(:)         !< [umol CO2/m2 stem/s @25C] gathered per-PFT
      real(wp),    allocatable :: stem_resp_sapwood25(:)        !< [umol CO2/m3 sapwood/s @25C] gathered per-PFT
      real(wp),    allocatable :: sap_volume(:)                 !< [m3/plant] sapwood volume (stem respiration)
      real(wp),    allocatable :: root_resp_factor25(:)         !< [umol CO2/kgC root/s @25C] gathered per-PFT
      real(wp),    allocatable :: leaf_emissivity(:), wood_emissivity(:)   !< [-] gathered per-PFT (tissue longwave)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   !< [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              !< [kgC/plant],[kgC/plant],[m2] (hydraulics)
      !----- Total wood carbon, not interchangeable with bsap: the sapwood sets the hydraulic      !
      !      capacitance, while all the wood, thin branches included, is the one thermal store    !
      !      (docs/science/vegetation_energy_dynamics.md §The store is ON). ------------------------!
      real(wp),    allocatable :: bwood(:)                                    !< [kgC/plant] total wood (thermal store)
      real(wp),    allocatable :: vcmax25(:), rd25(:)                         !< [umol/m2/s] per-cohort (plastic) capacities
      !----- Yesterday's daily-maximum leaf water potential, the predawn potential behind the slow  !
      !      stomatal control (docs/science/leaf_gas_exchange.md §4.2). column_prepass always passes !
      !      it, since leaf_gas_exchange_batch reads an absent `psi` as no stress. -----------------!
      real(wp),    allocatable :: dmax_psi_leaf(:)                        !< [MPa] <= 0; UNSET_NEGATIVE => seed from soil
   end type column_cohort_t

   !----- Prescribed per-step forcing the higher layers (RT, met) supply. ---------------------!
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: rainfall        = 0.0_wp   !< [kg/m2/s] reference-level rainfall (intercepted downstream)
      real(wp)              :: snowfall         = 0.0_wp   !< [kg/m2/s] frozen rainfall (snowfall; drives snow accumulation)
      real(wp)              :: air_temp          = 288.0_wp !< [K] reference-level air temp (frozen/rain-on-snow rainfall enthalpy)
      real(wp)              :: par_per_w     = 2.1_wp   !< [umol photon / (W absorbed)] absorbed->PAR-photon factor
      !----- Top-of-canopy radiative fluxes per band, kept from the two-stream solve for comparison with a
      !      radiometer or a satellite product (docs/science/canopy_radiation_transfer.md §Surface radiative
      !      fluxes on the diagnostic path). Upwelling longwave includes the surface's own emission.
      real(wp) :: sw_in_vis = 0.0_wp, sw_in_nir = 0.0_wp   !< [W/m2] incident SW at canopy top, per band
      real(wp) :: sw_up_vis = 0.0_wp, sw_up_nir = 0.0_wp   !< [W/m2] upwelling SW leaving the canopy top
      real(wp) :: lw_up     = 0.0_wp                        !< [W/m2] upwelling LW leaving the canopy top
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW (VIS+NIR) / net LW per cohort (leaf energy)
      real(wp), allocatable :: abs_par(:)               !< [W/m2] incident-equivalent PAR (VIS) per cohort; the leaf
                                                        !< model applies leaf_absorptance itself (photosynthesis)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] absorbed SW / net LW per cohort (wood energy)
      real(wp), allocatable :: par_top(:)               !< [W/m2] PAR (VIS) reaching the cohort's top: the
                                                        !< phenology PAR cue's driver (canopy light gradient)
   end type column_forcing_t

   !----- The per-patch conservation budgets. The per-kernel budgets close by construction;         !
   !      whole_energy and whole_water set every store against the column's true boundary fluxes,  !
   !      which is what catches a leak across a seam (docs/science/numerical_scheme.md §5b). -------!
   type :: column_budget_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
      !----- Per-layer face mass closure, summed over layers: |the water a layer's theta moved - the   !
      !      water its faces were charged for|. It sees the vertical misplacement the whole-column     !
      !      ledgers cannot (docs/science/numerical_scheme.md §5b). A pre-formed residual, so only     !
      !      resid/worst/abs_sum/n_check/n_fail are set; n_fail counts unconverged water solves. ------!
      type(budget_t) :: soil_face_mass
      real(wp)       :: gpp_last = 0.0_wp, nee_last = 0.0_wp   !< [umol/m2/s] last-step diagnostics
      !----- Net CAS -> atmosphere export over the step, b-weighted like the ledgers: the reported H, LE. !
      real(wp)       :: atm_heat_export = 0.0_wp   !< [J/m2]  sensible heat over this dt_fast
      real(wp)       :: atm_vap_export  = 0.0_wp   !< [kg/m2] over this dt_fast
      !----- Water leaving the soil column over the step, as the water network moved it: the pond's   !
      !      spill (surface runoff) and the bottom flux (drainage; negative when an aquifer feeds it). !
      real(wp)       :: runoff_surf_export = 0.0_wp   !< [kg/m2] over this dt_fast
      real(wp)       :: drainage_export    = 0.0_wp   !< [kg/m2] over this dt_fast
      !----- Soil-carbon Rh diagnostics (docs/science/soil_carbon.md §7), set by column_prepass when    !
      !      soil carbon is on (else 0). The fast driver accumulates them into the patch's day          !
      !      integrals xi_accum and rh_fast_accum. ----------------------------------------------------!
      real(wp)       :: xi_step(N_SOIL_POOL) = 0.0_wp   !< [-] this sub-step's per-pool env scalar
      real(wp)       :: rh_matrix_step       = 0.0_wp   !< [kgC/m2/day] this sub-step's matrix Rh
      !----- Work counters, the cost axis of a numerics comparison: set per sub-step by the stepper,  !
      !      accumulated site-wide by the fast driver. -----------------------------------------------!
      integer(ik)    :: integ_nsteps   = 0_ik   !< accepted integrator sub-steps this dt_fast
      integer(ik)    :: integ_nrej   = 0_ik   !< rejected integrator steps this dt_fast
      !----- The soil-plant water solver's work this dt_fast (docs/science/soil_plant_water.md §6). ---!
      integer(ik)    :: water_nsub   = 0_ik   !< accepted sub-steps
      integer(ik)    :: water_newton = 0_ik   !< Newton iterations
      integer(ik)    :: water_reject = 0_ik   !< rejected sub-steps (too large an error, or Newton failed)
      integer(ik)    :: water_nfail  = 0_ik   !< sub-steps taken at the shortest length although unconverged
      !----- dt_fast steps whose RK45 step bailed or committed a railed state and were redone on ARK  !
      !      (docs/science/numerical_scheme.md §3, `rk45`); a persistently high count flags a stiff   !
      !      regime. ---------------------------------------------------------------------------------!
      integer(ik)    :: rk45_rescue  = 0_ik   !< dt_fast steps rescued RK45->ARK this sub-step (0 on ARK)
      !----- Stage-clamp activations: how often clamp_cas bounded a throwaway stage input's canopy air !
      !      so the right-hand side stays evaluable. A rising count means the trial step is too big;   !
      !      neither scheme clamps the committed state. ----------------------------------------------!
      integer(ik)    :: clamp_stage_n  = 0_ik   !< stage-input clamp activations (layers/stores clamped)
      !----- How far theta lay outside its constitutive domain at an RK45 stage-1 input, which no     !
      !      clamp sees; telemetry, not a correction (docs/science/numerical_scheme.md §5b). 0 on ARK. !
      real(wp)       :: theta_ood_max  = 0.0_wp !< [m3/m3] max excursion outside [theta_res, theta_sat]
   end type column_budget_t

   !----- The prognostic CAS surface state advanced by the fast loop. ---------------------------!
   type :: surface_state_t
      real(wp) :: cas_enthalpy = 0.0_wp    !< [J/kg]      canopy-air specific enthalpy
      real(wp) :: cas_shv      = 0.0_wp    !< [kg/kg]     canopy-air specific humidity
      real(wp) :: cas_co2      = 0.0_wp    !< [umol/mol]  canopy-air CO2 mixing ratio
   end type surface_state_t

   !=====================================================================================!
   ! The frozen record, by physical content (docs/science/numerical_scheme.md §2): what the     !
   ! pre-pass (build_column_frozen) evaluates once at state^n and every stage of either         !
   ! integrator reads. A kernel takes only the pieces it reads (surface_derivs: cas, tissue,    !
   ! film, ground, snow). The CAS <-> atmosphere conductances are re-solved per stage.          !
   !=====================================================================================!

   !----- CAS <-> atmosphere boundary: capacities, bulk conductances, the surface-layer inputs a  !
   !      per-stage re-solve needs, and the reference-level state. Plain scalars, so a stage can  !
   !      copy it and override the conductances. ------------------------------------------------!
   type :: cas_boundary_t
      real(wp) :: cas_mass_capacity          = 0.0_wp      !< [kg/m2]   CAS mass capacity  -> enthalpy & vapour
      real(wp) :: cas_molar_capacity          = 0.0_wp      !< [mol/m2]  CAS molar capacity -> CO2
      real(wp) :: g_atm_heat           = 0.0_wp      !< [kg/m2/s] CAS<->atm enthalpy conductance
      real(wp) :: g_atm_vapour           = 0.0_wp      !< [kg/m2/s] CAS<->atm vapour   conductance
      real(wp) :: g_atm_co2           = 0.0_wp      !< [mol/m2/s]CAS<->atm CO2      conductance
      !----- g_atm_* above are the state^n values, which seed the march and the diagnostics; a stage  !
      !      re-solves them at its own canopy-air state (docs/science/numerical_scheme.md §2, item 3; !
      !      docs/science/canopy_aerodynamics.md §2). mo_live says the mo_* inputs below are set: it  !
      !      is false for a hand-built record, which supplies its own conductances (a re-solve from    !
      !      unset roughness and wind would not converge), true from build_column_frozen, and false   !
      !      again in ARK's stage copy once refreshed, so the Newton's repeated evaluations reuse it.  !
      logical :: mo_live = .false.            !< the mo_* inputs are populated => re-solve per evaluation
      type(aero_cfg_t) :: aero_cfg            !< MO constants (copy; POD)
      real(wp) :: mo_u_ref      = 0.0_wp      !< [m/s]  wind at the reference height
      real(wp) :: mo_zref       = 0.0_wp      !< [m]    reference height
      real(wp) :: mo_displace   = 0.0_wp      !< [m]    displacement height  (canopy geometry, frozen)
      real(wp) :: mo_rough      = 0.0_wp      !< [m]    roughness length     (canopy geometry, frozen)
      real(wp) :: mo_theta_atm  = 0.0_wp      !< [K]    potential temperature referenced to zref (= the air temperature there)
      real(wp) :: mo_shv_atm    = 0.0_wp      !< [kg/kg] specific humidity at zref (the AERO reference,
                                              !<        which need not equal shv_atm below)
      real(wp) :: mo_rho        = 0.0_wp      !< [kg/m3] air density
      real(wp) :: enthalpy_atm      = 0.0_wp      !< [J/kg]    reference-level specific enthalpy
      real(wp) :: shv_atm       = 0.0_wp      !< [kg/kg]   reference-level specific humidity
      real(wp) :: co2_atm       = 400.0_wp    !< [umol/mol]free-atmosphere CO2
      real(wp) :: nee_biotic    = 0.0_wp      !< [umol/m2/s] frozen biotic CO2 source (Ra+Rh-GPP)
      real(wp) :: rho           = 0.0_wp      !< [kg/m3]   canopy-air density
      real(wp) :: press         = 0.0_wp      !< [Pa]      canopy-air pressure
   end type cas_boundary_t

   !----- Per-cohort leaf and wood energy-balance coefficients and the tissue heat store, frozen at  !
   !      state^n. Each stage relaxes a tissue from t_leaf0/t_wood0 against hcap_per_dt = cap/dt_fast, !
   !      an algebraic closure rather than a tableau unknown, so it needs no error group and no Newton !
   !      (docs/science/vegetation_energy_dynamics.md §One formulation: the tissue relaxes exactly). -!
   type :: tissue_coefficients_t
      real(wp), allocatable :: h_coeff_leaf(:)   !< [W/m2/K]  frozen sensible coefficient
      real(wp), allocatable :: g_transp_leaf(:)      !< [m/s]     frozen leaf transpiration series conductance
      real(wp), allocatable :: abs_sw(:)      !< [W/m2]    absorbed shortwave (frozen source)
      real(wp), allocatable :: abs_lw(:)      !< [W/m2]    net longwave at the emission base (frozen source)
      real(wp), allocatable :: lai(:)         !< [m2/m2]   cohort leaf area index
      real(wp), allocatable :: h_coeff_w(:)   !< [W/m2/K]  frozen WOOD sensible coefficient (pi*wai*wood_gbh*rho*cp)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] frozen absorbed SW / net LW on wood
      real(wp), allocatable :: wai(:)         !< [m2/m2]   cohort wood area index
      real(wp), allocatable :: leaf_hcap_per_dt(:), wood_hcap_per_dt(:)   !< [W/m2/K] cap/dt_fast
      real(wp), allocatable :: t_leaf0(:), t_wood0(:) !< [K]      start-of-step tissue temperatures
      !----- Heat the moving plant water brings each tissue, and the enthalpy each tissue's own water  !
      !      gains as its mass changes, both at start-of-step temperatures (docs/science/soil_plant_water.md !
      !      eqs. (4)-(5)). ----------------------------------------------------------------------------!
      real(wp), allocatable :: qwflux_wl(:)   !< [W/m2 ground] into the leaf, from sapflow (wood->leaf)
      real(wp), allocatable :: q_wood_net(:)  !< [W/m2 ground] into the wood, from root uptake less sapflow
      real(wp), allocatable :: leaf_water_enth(:)  !< [W/m2 ground] the leaf's water, eq. (5)
      real(wp), allocatable :: wood_water_enth(:)  !< [W/m2 ground] the wood's water, eq. (5)
      !----- Each cohort's transpiration over the step [kg/m2 ground/s], from the water network         !
      !      (docs/science/soil_plant_water.md §4). With transp_held the leaf energy balance holds it, so  !
      !      the canopy air gains exactly the water the leaves lose; without (the pre-pass's own call, a  !
      !      hand-built record) the demand follows the stomata and the canopy air's humidity. ----------!
      logical :: transp_held = .false.
      real(wp), allocatable :: transp_net(:)
      !----- Per-cohort longwave emissivities, the ones the radiation solver used, so each tissue's     !
      !      emission slope 4*eps*sigma*T^3 matches the longwave it was given. ----------------------!
      real(wp), allocatable :: leaf_emiss(:)  !< [-]       leaf LW emissivity
      real(wp), allocatable :: wood_emiss(:)  !< [-]       wood LW emissivity
      !----- The heat-capacity inputs behind leaf_hcap_per_dt/wood_hcap_per_dt = (dry hcap + water !
      !      mass*cp_liq)/dt (docs/science/vegetation_energy_dynamics.md eq. (8)). ----------------!
      real(wp), allocatable :: wood_dry_hcap(:)   !< [J/m2/K]  dry wood heat capacity
      real(wp), allocatable :: wood_wmass(:)      !< [kg/m2]   the wood store's water (sapwood and fine roots)
      real(wp), allocatable :: leaf_dry_hcap(:)   !< [J/m2/K]  dry leaf heat capacity
      real(wp), allocatable :: leaf_wmass(:)      !< [kg/m2]   internal (symplast) leaf water mass
   end type tissue_coefficients_t

   !----- Canopy interception film: the wetted fraction and film-evaporation conductances, the     !
   !      liquid enthalpy the film is valued at, and the interception rates, frozen per dt_fast. All !
   !      zero when canopy_water_on is off; only the state-dependent terms are evaluated per stage. !
   type :: canopy_film_capacity_t
      real(wp), allocatable :: g_film_leaf(:), g_film_w(:)   !< [m/s] frozen film-evap conductance, leaf/wood
      real(wp), allocatable :: f_wet_c(:)                 !< [-]   frozen combined wetted fraction (sigma_w)
      !----- The film's liquid enthalpy, internal_energy_liquid(t_film_valuation). The tissue pays     !
      !      enthalpy_vapor(T) minus this per kg it evaporates, so film, tissue and CAS close exactly. -!
      real(wp) :: film_liquid_enthalpy    = 0.0_wp      !< [J/kg]
      !----- Frozen interception rates (capture/capacity only, e_canopy = 0): integrating them by    !
      !      explicit Euler over dt_fast reproduces the one-shot bucket commit exactly. -------------!
      real(wp), allocatable :: intercept_leaf(:), intercept_wood(:)   !< [kg/m2 ground/s] (ncoh)
   end type canopy_film_capacity_t

   !----- Bare-ground boundary: radiation reaching the ground, the ground <-> CAS conductance and   !
   !      the ground evaporation the frozen hydrology authority committed to. ----------------------!
   type :: ground_boundary_t
      real(wp) :: abs_sw_ground = 0.0_wp      !< [W/m2]    shortwave reaching the ground (frozen source)
      real(wp) :: abs_lw_ground = 0.0_wp      !< [W/m2]    net longwave at the ground (frozen source)
      real(wp) :: ggnet         = 0.0_wp      !< [m/s]     ground<->CAS aerodynamic conductance
      real(wp) :: soil_evap     = 0.0_wp      !< [kg/m2/s] ground latent flux (frozen hydrology authority)
   end type ground_boundary_t

   !----- Surface-block tendencies and the diagnostics the ledgers and the plant-water update use   !
   !      (transp_c is each cohort's leaf water loss). ------------------------------------------!
   type :: surface_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp     !< [J/kg/s]     dH/dt
      real(wp) :: d_cas_shv      = 0.0_wp     !< [kg/kg/s]    dq/dt
      real(wp) :: d_cas_co2      = 0.0_wp     !< [umol/mol/s] dC/dt
      real(wp) :: src_enth       = 0.0_wp     !< [W/m2]     summed surface enthalpy source into the CAS
      real(wp) :: src_vap        = 0.0_wp     !< [kg/m2/s]  summed surface vapour source into the CAS
      real(wp) :: g_top          = 0.0_wp     !< [W/m2]     net energy into the soil-top store
      real(wp) :: h_ground       = 0.0_wp     !< [W/m2]     ground sensible flux to the CAS
      real(wp) :: le_ground      = 0.0_wp     !< [W/m2]     ground latent flux to the CAS
      real(wp) :: coh_rnet       = 0.0_wp     !< [W/m2]     net radiation absorbed by the canopy
      real(wp) :: cond           = 0.0_wp     !< [kg/m2/s]  smooth condensation sink (dew) draining CAS supersat
      real(wp) :: cond_enth      = 0.0_wp     !< [W/m2]     liquid enthalpy the sink debited from the CAS
      real(wp), allocatable :: leaf_temp(:)   !< [K]        diagnosed per-cohort leaf temperature
      real(wp), allocatable :: wood_temp(:)   !< [K]        diagnosed per-cohort wood temperature
      real(wp), allocatable :: transp_c(:)    !< [kg/m2/s]  per-cohort transpiration demand
      !----- Per-cohort film evaporation (dew if negative), the films' analogue of transp_c. Zero    !
      !      when canopy_water_on is off. -----------------------------------------------------------!
      real(wp), allocatable :: film_evap_leaf(:), film_evap_wood(:)   !< [kg/m2 ground/s]
   end type surface_tend_t

   !----- The prognostic column state the integrators advance: a flat copy of the patch's stores   !
   !      that the schemes combine and measure as one vector. Plant water is carried as mass; psi  !
   !      is diagnosed from it once per dt_fast for the pre-pass (docs/science/soil_plant_water.md). !
   type :: column_state_t
      real(wp) :: cas_enthalpy = 0.0_wp                    !< [J/kg]
      real(wp) :: cas_shv      = 0.0_wp                    !< [kg/kg]
      real(wp) :: cas_co2      = 0.0_wp                    !< [umol/mol]
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp   !< [J/m3]   per soil layer
      real(wp) :: theta(n_soil_layer_max)       = 0.0_wp   !< [m3/m3]  per soil layer
      real(wp), allocatable :: leaf_water_mass(:) !< [kg/plant] internal leaf water (ncoh)
      real(wp), allocatable :: wood_water_mass(:) !< [kg/plant] internal wood water (ncoh)
      !----- The pond: carried through the stages unchanged and committed from the water network.  !
      !      Like the films it carries no error group, so the controller scores no frozen quantity. !
      real(wp) :: w_surface      = 0.0_wp         !< [kg/m2] ponded surface water
      real(wp) :: w_surface_energy = 0.0_wp         !< [J/m2]  enthalpy of the ponded water
      !----- The interception films, distinct from the internal water above: a film evaporates with !
      !      no stomatal resistance. Per m2 of ground, unlike the per-plant fields above. ----------!
      real(wp), allocatable :: leaf_surf_water(:) !< [kg/m2 ground] canopy interception film on leaf
      real(wp), allocatable :: wood_surf_water(:) !< [kg/m2 ground] canopy interception film on wood
   end type column_state_t

   !----- Soil water over this dt_fast, from the soil-plant water solver (meds_spac_solver, run in    !
   !      build_column_frozen): its boundary flows, end-of-step stores, the temperatures that value   !
   !      the water crossing each boundary, and its interior faces. Every scheme commits its theta and !
   !      advects soil enthalpy on these same faces (docs/science/numerical_scheme.md §5b). ----------!
   type :: soil_hydrology_t
      real(wp) :: q_top         = 0.0_wp          !< [m/s]     net top-face water flux (infiltration - exfiltration - evaporation)
      real(wp) :: infiltration  = 0.0_wp          !< [kg/m2/s] pond -> top layer
      real(wp) :: drainage      = 0.0_wp          !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf   = 0.0_wp          !< [kg/m2/s] surface runoff
      real(wp) :: ground_inflow_enth = 0.0_wp     !< [W/m2]    the enthalpy it brings across the column's boundary
      real(wp) :: t_infil       = 0.0_wp          !< [K]       temperature of the infiltrating (pond) water
      real(wp) :: w_surface_energy1 = 0.0_wp        !< [J/m2]    end-of-step pond ENTHALPY
      real(wp) :: runoff_enth   = 0.0_wp          !< [W/m2]    enthalpy leaving with surface runoff
      real(wp) :: t_film_valuation     = 0.0_wp          !< [K]       valuation T of intercepted water (tsupercool_liq under a pack)
      !----- Water a supersaturated top layer pushed up into the pond, and the enthalpy it took, valued at !
      !      the layer's start-of-step temperature (the pond received the same number). ----------------!
      real(wp) :: exfiltration  = 0.0_wp          !< [kg/m2/s] top layer -> pond
      real(wp) :: exfil_enth    = 0.0_wp          !< [W/m2]    enthalpy leaving the top layer with it
      real(wp) :: t_bot         = 0.0_wp          !< [K]       bottom-layer soil temperature @ state^n
      real(wp) :: w_surface1    = 0.0_wp          !< [kg/m2]   end-of-step ponded surface water
      real(wp) :: w_flux_frozen(n_soil_layer_max) = 0.0_wp  !< [m/s]   DOWNWARD interior face flux, k=1..nsl-1
      !----- The network's end-of-step soil moisture, which every scheme commits. -------------------!
      real(wp), allocatable :: theta1(:)          !< [m3/m3]   committed post-step soil moisture (per layer)
      !----- (theta1 - theta^n) / dt_fast: the steady rate at which the stages move soil water, so the   !
      !      water a layer holds keeps pace with the enthalpy its faces carry in (see column_be_stage).  !
      real(wp) :: theta_rate(n_soil_layer_max) = 0.0_wp   !< [1/s]
   end type soil_hydrology_t

   !----- Root zone over this dt_fast, from the water network: the patch's root uptake, its shares by  !
   !      layer (sum 1; the root profile when the roots took nothing), and the liquid enthalpy it      !
   !      carries out of the soil, by cohort (what each wood is credited) and by layer, at each        !
   !      layer's start-of-step temperature. -----------------------------------------------------------!
   type :: root_zone_t
      real(wp) :: uptake        = 0.0_wp          !< [kg/m2/s] the patch's root uptake
      real(wp), allocatable :: root_share(:)      !< [-]       per-layer root-sink shares (sum = 1)
      real(wp), allocatable :: qloss_frozen(:)    !< [W/m2 ground] (ncoh) advected enthalpy leaving the soil with uptake
      !----- By layer [W/m2], so the soil loses heat where it loses water (docs/science/soil_plant_water.md §5). !
      real(wp), allocatable :: root_heat(:)       !< (nsl)
   end type root_zone_t

   !----- Plant water over this dt_fast, from the water network, per plant and constant across the  !
   !      stages (sapflow wood -> leaf, uptake soil -> wood). With the held transpiration (transp_net) !
   !      they move each cohort's water exactly to the network's end-of-step values. ----------------!
   type :: plant_water_t
      real(wp), allocatable :: nplant(:)                            !< [plants/m2] (ncoh)
      real(wp), allocatable :: sapflow_frozen(:), uptake_frozen(:)   !< [kg/plant/s] (ncoh)
   end type plant_water_t

   !----- The frozen outcome of the step's snow advance (meds_fast_frozen%advance_snow_stage). With  !
   !      no pack every field is 0 or .false., which reduces the consumers exactly to the snow-free   !
   !      form, so neither scheme needs a snow branch. ----------------------------------------------!
   type :: snow_stage_t
      logical  :: exists     = .false.   !< a pack is present (drives rainfall routing + t_film_valuation)
      real(wp) :: snowfac    = 0.0_wp    !< [-]        Niu-Yang cover fraction actually used
      real(wp) :: h_snow     = 0.0_wp    !< [W/m2]     snowfac-weighted sensible flux to the CAS
      real(wp) :: le_snow    = 0.0_wp    !< [W/m2]     snowfac-weighted latent (sublimation) flux
      real(wp) :: g_base     = 0.0_wp    !< [W/m2]     throttled base conduction into the soil top
      real(wp) :: subl_rate  = 0.0_wp    !< [kg/m2/s]  sublimation vapour source for the CAS
      real(wp) :: melt_rate  = 0.0_wp    !< [kg/m2/s]  meltwater to the ponding store (see t_melt)
      real(wp) :: ground_rad = 0.0_wp    !< [W/m2]     blended ground radiative input for the ledgers
      real(wp) :: acc_enth   = 0.0_wp    !< [J/m2]     rainfall enthalpy that entered the pack (boundary in)
      real(wp) :: swe0       = 0.0_wp    !< [kg/m2]    pack mass BEFORE the stage (ledger store term)
      real(wp) :: swe1       = 0.0_wp    !< [kg/m2]    pack mass AFTER  the stage (ledger store term)
      real(wp) :: enth0      = 0.0_wp    !< [J/m2]     pack internal energy BEFORE (ledger store term)
      real(wp) :: enth1      = 0.0_wp    !< [J/m2]     pack internal energy AFTER  (ledger store term)
      real(wp) :: melt_enth  = 0.0_wp    !< [J/m2] enthalpy the meltwater carries from the pack to the pond
      !----- The temperature with u_liq(T)*melt mass = melt_enth, handed to the water network as the  !
      !      pond inflow's, so the pond receives exactly melt_enth; t_3ple when nothing melts. --------!
      real(wp) :: t_melt     = t_3ple    !< [K] effective temperature of the meltwater
   end type snow_stage_t

   !----- The container: everything held constant over one dt_fast, by physical content. -----------!
   type :: column_frozen_t
      type(cas_boundary_t)         :: cas          !< CAS <-> atmosphere boundary
      type(tissue_coefficients_t)  :: tissue       !< per-cohort leaf/wood energy coefficients + heat store
      type(canopy_film_capacity_t) :: film         !< canopy interception film
      type(ground_boundary_t)      :: ground       !< bare-ground boundary
      type(snow_stage_t)           :: snow         !< the snow stage's outcome (all zero without a pack)
      type(soil_hydrology_t)       :: hydrology    !< the water network's soil flows over the step
      type(root_zone_t)            :: roots        !< the network's root uptake and its placement by layer
      type(plant_water_t)          :: plant        !< the network's sapflow and uptake per plant
   end type column_frozen_t

   !----- Storage one thread's fast step reuses from step to step: the soil-plant water solver and   !
   !      the ARK march's states and stages, filled in place and reallocated only when the cohort     !
   !      count changes. The stages are their own record so the march passes them beside its states. !
   type :: ark_stages_t
      type(column_state_t) :: y2, base3, y3            !< the ARS(2,2,2) stages
      type(surface_tend_t) :: sf2, sf3                 !< the stages' surface tendencies
   end type ark_stages_t

   type :: fast_workspace_t
      type(spac_solver_t)   :: spac                          !< the soil-plant water solver
      type(column_state_t)  :: y_march, y_new, y_err, y_lo   !< the ARK march's state and trial step
      type(ark_stages_t)    :: stages                        !< the ARK step's stages
   end type fast_workspace_t

   !----- The whole-column tendency vector + diagnostics. ---------------------------------------!
   type :: column_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp
      real(wp) :: d_cas_shv      = 0.0_wp
      real(wp) :: d_cas_co2      = 0.0_wp
      real(wp) :: dedt(n_soil_layer_max)   = 0.0_wp   !< [W/m3] dsoil_energy/dt
      real(wp) :: dtheta_dt(n_soil_layer_max) = 0.0_wp!< [1/s]  dtheta/dt
      real(wp), allocatable :: d_leaf_water_mass(:)   !< [kg/plant/s] frozen_sapflow - transp(stage)
      real(wp), allocatable :: d_wood_water_mass(:)   !< [kg/plant/s] frozen_uptake  - frozen_sapflow
      !----- The films' analogue of d_leaf_water_mass/d_wood_water_mass. ---------------------------!
      real(wp), allocatable :: d_leaf_surf_water(:)   !< [kg/m2 ground/s] frozen_intercept_leaf - film_evap
      real(wp), allocatable :: d_wood_surf_water(:)   !< [kg/m2 ground/s] frozen_intercept_wood - film_evap
      real(wp) :: g_top = 0.0_wp, drainage_rate = 0.0_wp, uptake_rate = 0.0_wp
      real(wp) :: g_bottom = 0.0_wp   !< [W/m2] heat across the soil column's bottom face, positive UP
      !----- The interior faces this stage's soil-energy equation advected enthalpy on, taken from what !
      !      assemble_soil_energy_forcing was handed, not from the water solve, so RK45 can check them  !
      !      against the committed theta change (docs/science/numerical_scheme.md §5b). ---------------!
      real(wp) :: soil_face(n_soil_layer_max) = 0.0_wp   !< [m/s] DOWNWARD interior face, k=1..nsl-1
      real(wp), allocatable :: leaf_temp(:)
   end type column_tend_t

   !----- ARK conservation ledger: per-stage boundary-flux rates (from column_be_stage), and their    !
   !      b-weighted amounts summed over the accepted sub-steps (column_bflux_t). Since Y3 - y =      !
   !      (1-gamma)*h*K2 + gamma*h*K3 (docs/science/numerical_scheme.md eq. (1)), the amounts         !
   !      telescope against the committed store change to machine precision. ----------------------!
   type :: stage_bflux_t                                    !< per-stage RATES
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp    !< [W/m2]
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp    !< [kg/m2/s]
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp    !< [umol/m2/s]
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp  !< [W/m2]
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp  !< [kg/m2/s]
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp!< [W/m2]
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp!< [kg/m2/s]
      real(wp) :: whole_cond    = 0.0_wp                         !< [kg/m2/s] condensate
      real(wp) :: whole_cond_enth = 0.0_wp                       !< [W/m2] its liquid enthalpy at the stage CAS temperature
      !----- Net CAS -> atmosphere export at the ledger's own stage states and conductances: the      !
      !      reported H and LE (meds_fast_step%atm_fluxes; docs/science/numerical_scheme.md §5b). -----!
      real(wp) :: atm_heat_out = 0.0_wp                          !< [W/m2]    g_atm_heat*cp_air*(T_cas - theta_atm)
      real(wp) :: atm_vap_out  = 0.0_wp                          !< [kg/m2/s] g_atm_vapour*(q_cas - q_atm)
   end type stage_bflux_t

   type :: column_bflux_t                                  !< accumulated AMOUNTS (J/m2, kg/m2, umol/m2)
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp
      !----- Condensate, apart from whole_wat_out: dew and fog land inside the column, and the caller  !
      !      deposits them into soil layer 1 (docs/science/numerical_scheme.md §7, item 4). ---------!
      real(wp) :: whole_cond    = 0.0_wp   !< [kg/m2] condensed vapour over the step (>= 0)
      !----- ...and the liquid enthalpy it left the CAS with, b-weighted at each stage's own CAS     !
      !      temperature, as surface_derivs debited it; the deposit carries this, not u_liq(T_end). -!
      real(wp) :: whole_cond_enth = 0.0_wp !< [J/m2]
      real(wp) :: atm_heat_out = 0.0_wp    !< [J/m2]  net CAS -> atmosphere sensible-heat export (reported H)
      real(wp) :: atm_vap_out  = 0.0_wp    !< [kg/m2] net CAS -> atmosphere vapour export (reported LE / L_v)
      !----- Per-cohort tissue-temperature time integrals [K*s], b-weighted over the stages and summed !
      !      over the accepted sub-steps. A tissue's energy is committed from this time average, the   !
      !      complement of what the canopy air received (docs/science/numerical_scheme.md §5b). -------!
      real(wp), allocatable :: tissue_leaf_int(:), tissue_wood_int(:)   !< [K*s]
      !----- Each cohort's transpiration over the step [kg/m2 ground], b-weighted the same way: the  !
      !      water its leaves were debited and the canopy air credited, which the diagnostics report. !
      real(wp), allocatable :: transp_coh(:)
   end type column_bflux_t

   !----- A patch's fast-loop state, the working copy the integrators advance: gathered from the   !
   !      persisted patch and cohort state by meds_fast_dynamics and written back after. ----------!
   type :: patch_biophys_t
      type(cas_state_t)          :: cas               !< canopy-air-space twins (enthalpy/shv/co2)
      type(soil_energy_column_t) :: soil_energy            !< soil thermal column (internal energy; temp diagnosed)
      type(soil_water_column_t)  :: soil_water            !< soil water column (theta; psi_soil diagnosed)
      type(snow_column_t)        :: snow              !< temporary-surface-water / snow store (swe + energy)
      !----- Read-only snapshot of the patch's soil-carbon pools (docs/science/soil_carbon.md §7),    !
      !      seeded once per slow step and held over the day's fast sub-steps: the fast Rh respires   !
      !      against it, and the daily soil_carbon_step alone writes the real pools. 0 when off. ----!
      type(soil_carbon_t)        :: soil_carbon
      !----- The day's leaf/root-turnover shed water, a read-only snapshot like soil_carbon, added to  !
      !      the ground's water input every sub-step. ---------------------------------------------!
      real(wp)                   :: shed_water_rate = 0.0_wp  !< [kg/m2 ground/s]
      real(wp)                   :: shed_water_enth = 0.0_wp  !< [W/m2 ground] its enthalpy, from the tissues it left
      !----- The slow tier's CO2 owing (growth respiration and its corrections,                    !
      !      docs/science/plant_respiration.md §6), seeded and held the same way. -----------------!
      real(wp)                   :: slow_co2_rate   = 0.0_wp  !< [umol/m2 ground/s]
      real(wp), allocatable      :: leaf_energy(:)    !< [J/m2 ground] per-cohort leaf internal energy
      real(wp), allocatable      :: wood_energy(:)    !< [J/m2 ground] per-cohort wood internal energy
      real(wp), allocatable      :: leaf_temp(:)      !< [K] derived from leaf_energy
      real(wp), allocatable      :: wood_temp(:)      !< [K] derived from wood_energy
      !----- Internal water mass, the prognostic hydraulic state; psi is diagnosed from it         !
      !      (tissue_water_to_psi) and never persisted (docs/science/soil_plant_water.md §2). ------!
      real(wp), allocatable      :: leaf_water_mass(:) !< [kg/plant] internal leaf water
      real(wp), allocatable      :: wood_water_mass(:) !< [kg/plant] internal wood water
      !----- The interception films, per m2 of ground here; the cohort holds them per plant. ------!
      real(wp), allocatable      :: leaf_surf_water(:) !< [kg/m2 ground] leaf interception film
      real(wp), allocatable      :: wood_surf_water(:) !< [kg/m2 ground] wood interception film
      !----- Warm start for the adaptive march (docs/science/numerical_scheme.md §3, `ark`): the step  !
      !      that seeds the next dt_fast's first sub-step; 0 = no history, start from dt_fast. --------!
      real(wp)                   :: adapt_dt_last = 0.0_wp   !< [s] last accepted sub-step not cut by the step's end
   end type patch_biophys_t

contains

   !----- Allocate a column_cohort_t (the per-patch cohort SoA the fast loop consumes). ------!
   subroutine alloc_column_cohort(col_cohort, n)
      type(column_cohort_t), intent(out) :: col_cohort
      integer(ik),           intent(in)  :: n
      col_cohort%n = n
      allocate(col_cohort%pft(n), col_cohort%lai(n), col_cohort%wai(n), col_cohort%height(n),                                     &
               col_cohort%leaf_width(n), col_cohort%branch_diam(n), col_cohort%aboveground_frac(n),                   &
               col_cohort%is_woody(n), col_cohort%stem_resp_factor25(n), col_cohort%root_resp_factor25(n),            &
               col_cohort%stem_resp_sapwood25(n), col_cohort%sap_volume(n),                                         &
               col_cohort%leaf_emissivity(n), col_cohort%wood_emissivity(n),                                         &
               col_cohort%leaf_area(n), col_cohort%nplant(n),                                                          &
               col_cohort%dbh(n), col_cohort%broot(n), col_cohort%bleaf(n), col_cohort%bsap(n), col_cohort%sap_area(n),           &
               col_cohort%bwood(n), col_cohort%vcmax25(n), col_cohort%rd25(n), col_cohort%dmax_psi_leaf(n))
      col_cohort%pft = 1_ik
      col_cohort%lai = 0.0_wp
      col_cohort%wai = 0.0_wp
      col_cohort%height = 0.0_wp
      col_cohort%leaf_width = 0.04_wp
      col_cohort%branch_diam = 0.02_wp
      col_cohort%aboveground_frac = 0.7_wp
      col_cohort%is_woody = .true.
      col_cohort%stem_resp_factor25 = 0.0_wp
      col_cohort%root_resp_factor25 = 0.0_wp
      col_cohort%stem_resp_sapwood25 = 0.0_wp
      col_cohort%sap_volume = 0.0_wp
      col_cohort%leaf_emissivity = 1.0_wp
      col_cohort%wood_emissivity = 1.0_wp
      col_cohort%leaf_area = 0.0_wp
      col_cohort%nplant = 0.0_wp
      col_cohort%dbh = 0.0_wp
      col_cohort%broot = 0.0_wp
      col_cohort%bleaf = 0.0_wp
      col_cohort%bsap = 0.0_wp
      col_cohort%sap_area = 0.0_wp
      col_cohort%bwood = 0.0_wp
      col_cohort%vcmax25 = 0.0_wp
      col_cohort%rd25 = 0.0_wp
      !----- Unset, not 0: a 0 would read as fully turgid and disable the stomatal stress limb for a  !
      !      caller that forgets to fill it; the sentinel makes column_prepass seed from the soil. ----!
      col_cohort%dmax_psi_leaf = UNSET_NEGATIVE
   end subroutine alloc_column_cohort

   !----- Grow-only capacity: reallocate only when the arrays hold fewer than n cohorts, else just   !
   !      set the active count. Every reader loops to col_cohort%n and the gather overwrites 1..n,   !
   !      so the caller sizes it once per fast_dynamics call, to the site's largest patch. ---------!
   subroutine ensure_column_cohort_capacity(col_cohort, n)
      type(column_cohort_t), intent(inout) :: col_cohort
      integer(ik),            intent(in)    :: n
      if (.not. allocated(col_cohort%pft)) then
         call alloc_column_cohort(col_cohort, n)
      else if (size(col_cohort%pft) < n) then
         call alloc_column_cohort(col_cohort, n)
      else
         col_cohort%n = n
      end if
   end subroutine ensure_column_cohort_capacity

   !----- Allocate a patch_biophys_t and seed it from an initial canopy-air state. --------------!
   subroutine alloc_patch_biophys(biophys, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(patch_biophys_t), intent(out) :: biophys
      integer(ik),           intent(in)  :: n_coh
      real(wp),              intent(in)  :: can_temp0, can_shv0, can_co2, leaf_temp0
      allocate(biophys%leaf_energy(n_coh), biophys%wood_energy(n_coh))
      allocate(biophys%leaf_temp(n_coh), biophys%wood_temp(n_coh))
      allocate(biophys%leaf_water_mass(n_coh), biophys%wood_water_mass(n_coh))
      allocate(biophys%leaf_surf_water(n_coh), biophys%wood_surf_water(n_coh))
      biophys%leaf_energy      = 0.0_wp    ! scratch seed: the gather sets it, a column fixture with tissue_temp_to_energy
      biophys%wood_energy      = 0.0_wp
      biophys%leaf_temp        = leaf_temp0
      biophys%wood_temp        = leaf_temp0
      biophys%leaf_water_mass  = 0.0_wp    ! scratch seeds: the next gather overwrites them
      biophys%wood_water_mass  = 0.0_wp
      biophys%leaf_surf_water  = 0.0_wp
      biophys%wood_surf_water  = 0.0_wp
      biophys%cas%can_temp     = can_temp0
      biophys%cas%can_shv      = can_shv0
      biophys%cas%can_co2      = can_co2
      biophys%cas%can_enthalpy = cas_temp_to_enthalpy(can_temp0, can_shv0)
   end subroutine alloc_patch_biophys

   !----- Grow-only capacity for patch_biophys_t's per-cohort arrays, as ensure_column_cohort_capacity. !
   !      The patch stores (cas, soil, snow, soil carbon) need no seed: every caller overwrites them   !
   !      with the patch's persisted state right after. -----------------------------------------------!
   subroutine ensure_patch_biophys_capacity(biophys, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(patch_biophys_t), intent(inout) :: biophys
      integer(ik),            intent(in)    :: n_coh
      real(wp),               intent(in)    :: can_temp0, can_shv0, can_co2, leaf_temp0
      if (.not. allocated(biophys%leaf_temp)) then
         call alloc_patch_biophys(biophys, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      else if (size(biophys%leaf_temp) < n_coh) then
         call alloc_patch_biophys(biophys, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      end if
   end subroutine ensure_patch_biophys_capacity

end module meds_fast_types
