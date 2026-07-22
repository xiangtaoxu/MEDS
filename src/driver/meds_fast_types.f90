!==========================================================================================!
! meds_fast_types -- the shared driver-scope TYPES of the fast-loop: the four working-buffer     !
! bundles both integrators pass through the dispatch seam (column_config_t/column_cohort_t/       !
! column_forcing_t/column_budget_t, from the former meds_column_dynamics) and the ARK POD state /  !
! frozen-input / tendency / boundary-flux-ledger types the whole-column RHS advances (surface_*_t, !
! column_state_t, column_frozen_t, column_tend_t, stage_bflux_t, column_bflux_t, from the former    !
! meds_column_derivs). Extracting them here is the prerequisite plumbing that lets meds_fast_split  !
! (the operator-split + Picard stepper) and meds_fast_ark (the IMEX-ARK stepper) be separate         !
! modules without a dynamics<->ark<->derivs cycle -- both link this leaf instead of each other's      !
! types.                                                                                            !
!                                                                                          !
! These are DRIVER-scope bundles, not kernel I/O contracts and not persistent state:                !
!   * column_config_t COMPOSES the biophysics/plant param+opts types (aero_cfg_t, soil_opts_t, ...) !
!     wholesale -- it is the aggregation seam meds_fast_dynamics%build_fast_context fills, not a      !
!     duplicate of them.                                                                             !
!   * column_state_t is a FLAT re-packing of the same prognostic quantities the persistent per-store  !
!     structs in meds_column_state_types hold (cas_state_t/soil_column_t/soil_energy_column_t) -- a    !
!     deliberate representation choice (the ARK needs a contiguous vector for state_axpy/state_wrms/   !
!     tableau linear combinations), not a duplication to unify.                                        !
!==========================================================================================!
module meds_fast_types
   use meds_kinds,            only : wp, ik
   use meds_biophysics_types, only : n_soil_layer_max, aero_cfg_t, veg_thermal_params_t,        &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,         &
                                     energy_opts_t, snow_params_t
   use meds_plant_types,      only : wood_params_t, root_params_t, hydro_params_t, hydro_opts_t
   use meds_biogeochem_types, only : co2_opts_t
   use meds_budget_check,     only : budget_t
   use meds_config,           only : hydraulics_config_t
   use meds_hydr_lib,         only : build_hydro_table
   implicit none
   private

   public :: LEAFEN_DIAGNOSTIC, LEAFEN_PROGNOSTIC, WOODEN_DIAGNOSTIC, WOODEN_PROGNOSTIC
   public :: SOILH2O_LAGGED, SOILH2O_COUPLED
   public :: column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: alloc_column_cohort, apply_hydraulics_config
   public :: surface_state_t, surface_frozen_t, surface_tend_t
   public :: column_state_t, column_frozen_t, column_tend_t
   public :: stage_bflux_t, column_bflux_t

   !----- Leaf/wood thermal model + soil-water-in-loop selector codes (P3). ------------------!
   integer(ik), parameter :: LEAFEN_DIAGNOSTIC = 0_ik  !< steady-state leaf (tl = tcas + dtl)
   integer(ik), parameter :: LEAFEN_PROGNOSTIC = 1_ik  !< prognostic leaf_energy via veg_energy_step_implicit (P3e)
   integer(ik), parameter :: WOODEN_DIAGNOSTIC = 0_ik  !< steady-state wood (own balance; tw = tcas + dtw)
   integer(ik), parameter :: WOODEN_PROGNOSTIC = 1_ik  !< prognostic wood_energy via veg_energy_step_implicit
   !----- RESERVED for the P3f re-solve-inside-Picard optimization; NOT YET WIRED -- both values   !
   !      take the identical frozen-after-pass-1 path in column_fast_step today (no behavioral      !
   !      branch exists on this selector; see the note in the Picard loop header there). ------------!
   integer(ik), parameter :: SOILH2O_LAGGED    = 0_ik  !< soil water/hydraulics frozen per sub-step
   integer(ik), parameter :: SOILH2O_COUPLED   = 1_ik  !< RESERVED (P3f): would re-solve inside the Picard loop

   !----- Static per-run column configuration (built once; constant across dt_fast steps). ----!
   type :: column_config_t
      type(aero_cfg_t)            :: aero            !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal    !< leaf/wood thermal params
      type(soil_params_t)         :: soil           !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal   !< soil thermal texture
      type(energy_opts_t)         :: energy         !< soil-thermal solver options
      type(soil_opts_t)           :: hydro          !< soil-water (Richards) solver options
      type(wood_params_t)         :: wood           !< stem-respiration parameters
      type(root_params_t)         :: root           !< fine-root-respiration parameters
      type(co2_opts_t)            :: co2            !< heterotrophic-respiration options
      type(hydro_params_t)        :: hydro_p        !< plant-hydraulics parameters (PV curves, vulnerability)
      type(hydro_opts_t)          :: hydro_o        !< plant-hydraulics solver options
      logical                     :: multilayer_roots  = .false.  !< opt-in soil->plant per-layer root coupling
      real(wp)                    :: specific_root_area = 20.0_wp  !< [m2/kgC] SRA (rhizosphere conductance)
      real(wp)                    :: fast_soil_carbon = 5.0_wp   !< [kgC/m2] decomposable soil-C pool (prescribed, MVP)
      real(wp)                    :: rhizo_cond       = 5.0e-4_wp !< [kg/s/MPa] soil->root conductance (prescribed, MVP)
      logical                     :: advect_soil_heat = .false.  !< opt-in: advect liquid enthalpy on the interior
                                                                 !< per-layer Darcy flux (moisture<->energy coupling)
      !----- P3 coupled-surface (Picard) solver knobs; only consulted under SCHEME_PICARD_COUPLED. !
      integer(ik) :: picard_max_iter = 20_ik        !< outer-iteration cap
      real(wp)    :: picard_tol_temp = 1.0e-3_wp     !< [K]     temperature convergence tolerance
      real(wp)    :: picard_tol_shv  = 1.0e-6_wp     !< [kg/kg] CAS specific-humidity convergence tolerance
      real(wp)    :: picard_relax    = 0.5_wp        !< [-]     under-relaxation of the next-pass seed. The CAS<->ground
                                                     !<         sensible coupling gives the fixed-point map a slope ~ -1
                                                     !<         (oscillatory); 0.5 makes the relaxed map a contraction.
      logical     :: picard_fixed_iter = .false.     !< run a uniform pass count (GPU warp-uniform; no early exit)
      integer(ik) :: leaf_energy_model  = 0_ik       !< LEAFEN_DIAGNOSTIC (0) | LEAFEN_PROGNOSTIC (1)
      integer(ik) :: wood_energy_model  = 0_ik       !< WOODEN_DIAGNOSTIC (0) | WOODEN_PROGNOSTIC (1)
      integer(ik) :: soil_water_coupling = 0_ik      !< SOILH2O_LAGGED (0) | SOILH2O_COUPLED (1); RESERVED, no
                                                     !< behavioral effect yet (P3f) -- see the selector comment above
      logical           :: snow_on = .false.         !< opt-in temporary-surface-water / snow store (P0, split path)
      type(snow_params_t) :: snow                    !< snow parameters (density, albedo, thresholds, conductivity)
   end type column_config_t

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into cfg%pft)
      real(wp),    allocatable :: lai(:), wai(:), height(:), crown(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   !< [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              !< [kgC/plant],[kgC/plant],[m2] (hydraulics)
      real(wp),    allocatable :: vcmax25(:), rd25(:)                         !< [umol/m2/s] per-cohort (plastic) capacities
   end type column_cohort_t

   !----- Prescribed per-step forcing the higher layers (RT, met) supply; photosynthesis/    !
   !      respiration/NEE are now computed from the plant kernels (no longer prescribed).     !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: precip        = 0.0_wp   !< [kg/m2/s] ground-reaching rainfall (interception deferred)
      real(wp)              :: snowf         = 0.0_wp   !< [kg/m2/s] frozen precip (snowfall; drives snow accumulation)
      real(wp)              :: tair          = 288.0_wp !< [K] reference-level air temp (frozen/rain-on-snow precip enthalpy)
      real(wp)              :: par_per_w     = 2.1_wp   !< [umol photon / (W absorbed)] absorbed->PAR-photon factor
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW (VIS+NIR) / net LW per cohort (leaf ENERGY)
      real(wp), allocatable :: abs_par(:)               !< [W/m2] INCIDENT-equiv PAR (VIS) per cohort; the leaf
                                                        !< model re-applies leaf_absorptance internally (PHOTOSYNTHESIS)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] absorbed SW / net LW per cohort (WOOD energy)
   end type column_forcing_t

   !----- The per-patch conservation budgets (one place; the driver accumulates the closed resids).!
   !      The per-kernel budgets close BY CONSTRUCTION; whole_energy/whole_water are the CROSS-      !
   !      seam column totals (Δ all stores vs the true boundary fluxes) that actually catch leaks.   !
   type :: column_budget_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
      real(wp)       :: gpp_last = 0.0_wp, nee_last = 0.0_wp   !< [umol/m2/s] last-step diagnostics
      !----- P3 Picard diagnostics (reporting only; not conserved state). --------------------!
      integer(ik)    :: picard_iters       = 0_ik    !< worst outer-iteration count over the sub-steps
      integer(ik)    :: picard_nonconv     = 0_ik    !< number of sub-steps that hit picard_max_iter unconverged
      real(wp)       :: picard_worst_resid = 0.0_wp  !< [K] worst residual temperature at exit
   end type column_budget_t

   !----- The prognostic CAS surface state advanced by the fast loop. ---------------------------!
   type :: surface_state_t
      real(wp) :: cas_enthalpy = 0.0_wp    !< [J/kg]      canopy-air specific enthalpy
      real(wp) :: cas_shv      = 0.0_wp    !< [kg/kg]     canopy-air specific humidity
      real(wp) :: cas_co2      = 0.0_wp    !< [umol/mol]  canopy-air CO2 mixing ratio
   end type surface_state_t

   !----- Frozen-per-substep inputs to the surface block (pre-pass coefficients, aerodynamic       !
   !      capacities/conductances, atmospheric BCs, lagged radiation + ground-latent forcing,       !
   !      and the soil-water supply fraction). Mirrors what column_fast_step freezes once per        !
   !      sub-step (meds_fast_split.f90).                                                            !
   type :: surface_frozen_t
      real(wp), allocatable :: h_coeff_f(:)   !< [W/m2/K]  frozen sensible coefficient
      real(wp), allocatable :: g_tr_f(:)      !< [m/s]     frozen leaf transpiration series conductance
      real(wp), allocatable :: abs_sw(:)      !< [W/m2]    absorbed shortwave (frozen source)
      real(wp), allocatable :: abs_lw(:)      !< [W/m2]    net longwave at the emission base (frozen source)
      real(wp), allocatable :: lai(:)         !< [m2/m2]   cohort leaf area index
      real(wp), allocatable :: h_coeff_w(:)   !< [W/m2/K]  frozen WOOD sensible coefficient (pi*wai*wood_gbh*rho*cp)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] frozen absorbed SW / net LW on wood
      real(wp), allocatable :: wai(:)         !< [m2/m2]   cohort wood area index
      real(wp) :: leaf_emiss    = 0.95_wp     !< [-]       leaf LW emissivity
      real(wp) :: wcap          = 0.0_wp      !< [kg/m2]   CAS mass capacity  -> enthalpy & vapour
      real(wp) :: ccap          = 0.0_wp      !< [mol/m2]  CAS molar capacity -> CO2
      real(wp) :: gah           = 0.0_wp      !< [kg/m2/s] CAS<->atm enthalpy conductance
      real(wp) :: gaw           = 0.0_wp      !< [kg/m2/s] CAS<->atm vapour   conductance
      real(wp) :: gac           = 0.0_wp      !< [mol/m2/s]CAS<->atm CO2      conductance
      real(wp) :: enth_atm      = 0.0_wp      !< [J/kg]    reference-level specific enthalpy
      real(wp) :: shv_atm       = 0.0_wp      !< [kg/kg]   reference-level specific humidity
      real(wp) :: co2_atm       = 400.0_wp    !< [umol/mol]free-atmosphere CO2
      real(wp) :: nee_biotic    = 0.0_wp      !< [umol/m2/s] frozen biotic CO2 source (Ra+Rh-GPP)
      real(wp) :: abs_sw_ground = 0.0_wp      !< [W/m2]    shortwave reaching the ground (frozen source)
      real(wp) :: abs_lw_ground = 0.0_wp      !< [W/m2]    net longwave at the ground (frozen source)
      real(wp) :: ggnet         = 0.0_wp      !< [m/s]     ground<->CAS aerodynamic conductance
      real(wp) :: soil_evap     = 0.0_wp      !< [kg/m2/s] ground latent flux (frozen hydrology authority)
      real(wp) :: rho           = 0.0_wp      !< [kg/m3]   canopy-air density
      real(wp) :: press         = 0.0_wp      !< [Pa]      canopy-air pressure
      real(wp) :: src_frac      = 1.0_wp      !< [-]       soil-water supply fraction (uptake / demand)
      real(wp) :: t_ground      = 0.0_wp      !< [K]       soil-top temperature (diagnosed from the state in column_derivs)
   end type surface_frozen_t

   !----- Surface-block tendencies + the diagnostics the ARK ledger and the soil/hydraulics         !
   !      tendencies consume (coh_qsoil -> soil-heat sink; coh_transp -> soil-water sink; transp_c   !
   !      -> per-cohort hydraulic demand). ------------------------------------------------------!
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
      real(wp) :: coh_qsoil      = 0.0_wp     !< [W/m2]     liquid enthalpy the soil sheds (post src_frac)
      real(wp) :: coh_transp     = 0.0_wp     !< [kg/m2/s]  total realized transpiration (post src_frac)
      real(wp) :: cond           = 0.0_wp     !< [kg/m2/s]  smooth condensation sink (dew) draining CAS supersat
      real(wp), allocatable :: leaf_temp(:)   !< [K]        diagnosed per-cohort leaf temperature
      real(wp), allocatable :: wood_temp(:)   !< [K]        diagnosed per-cohort wood temperature
      real(wp), allocatable :: transp_c(:)    !< [kg/m2/s]  per-cohort transpiration DEMAND (pre src_frac)
   end type surface_tend_t

   !----- The full prognostic column state advanced per dt_fast. --------------------------------!
   type :: column_state_t
      real(wp) :: cas_enthalpy = 0.0_wp                    !< [J/kg]
      real(wp) :: cas_shv      = 0.0_wp                    !< [kg/kg]
      real(wp) :: cas_co2      = 0.0_wp                    !< [umol/mol]
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp   !< [J/m3]   per soil layer
      real(wp) :: theta(n_soil_layer_max)       = 0.0_wp   !< [m3/m3]  per soil layer
      real(wp), allocatable :: psi(:,:)                    !< [MPa]    (N_HYDRO, ncoh) leaf/wood water potentials
   end type column_state_t

   !----- Frozen inputs for the whole column: the surface pre-pass + the soil/hydraulics params +   !
   !      the frozen hydrology surface BCs + per-cohort geometry the hydraulics kernel needs.        !
   type :: column_frozen_t
      type(surface_frozen_t)      :: surf         !< the surface-block frozen inputs (t_ground overwritten per call)
      type(soil_params_t)         :: soil         !< soil geometry + texture (dz, root_frac, ...)
      type(soil_thermal_params_t) :: therm        !< soil thermal texture
      type(energy_opts_t)         :: energy_opts  !< soil-thermal options (phase change)
      type(soil_opts_t)           :: hydro_opts   !< soil-water (Richards) options
      type(hydro_params_t)        :: hydro_p      !< plant-hydraulics parameters
      type(hydro_opts_t)          :: hydro_o      !< plant-hydraulics solver options
      real(wp) :: geothermal    = 0.0_wp          !< [W/m2]    bottom heat flux BC
      real(wp) :: q_top         = 0.0_wp          !< [m/s]     Richards top water flux (infiltration - evaporation)
      real(wp) :: soil_psi_root = 0.0_wp          !< [MPa]     root-zone soil water potential (hydraulics BC)
      real(wp) :: rhizo_cond    = 0.0_wp          !< [kg/s/MPa]soil->root conductance (hydraulics BC)
      !----- frozen boundary hydrology for the precip>0 guard-lift: the throughfall/drainage/runoff    !
      !      water carries internal_energy_liquid across the soil boundaries (matches the split's       !
      !      :436-439,518-520 advection), and the scratch column_hydrology_flux's end-of-step ponding/  !
      !      aquifer/water-table is persisted (column_state_t does NOT carry these surface stores). ----!
      real(wp) :: infiltration  = 0.0_wp          !< [kg/m2/s] throughfall reaching the soil top face
      real(wp) :: drainage      = 0.0_wp          !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf   = 0.0_wp          !< [kg/m2/s] surface runoff
      real(wp) :: rain_temp     = 0.0_wp          !< [K]       rain temperature (CAS temp @ state^n)
      real(wp) :: t_bot         = 0.0_wp          !< [K]       bottom-layer soil temperature @ state^n
      real(wp) :: w_surface1    = 0.0_wp          !< [kg/m2]   end-of-step ponded surface water
      real(wp) :: w_aquifer1    = 0.0_wp          !< [kg/m2]   end-of-step aquifer store
      real(wp) :: z_wt1         = 0.0_wp          !< [m]       end-of-step water-table elevation
      real(wp) :: uptake        = 0.0_wp          !< [kg/m2/s] realized root uptake (soil_wat_out ledger term)
      !----- the AUTHORITATIVE end-of-step soil moisture from the scratch column_hydrology_flux (the robust  !
      !      ponding/runoff/free-drain Richards solve). The ARK COMMITS this instead of re-solving theta in   !
      !      the ESDIRK stages (soil water is fully operator-split out; see column_fast_step_ark).            !
      real(wp), allocatable :: theta1(:)          !< [m3/m3]   committed post-step soil moisture (per layer)
      real(wp), allocatable :: psi_e(:)           !< [m]       Zeng-Decker equilibrium potential per layer (frozen)
      !----- per-cohort geometry the hydraulics kernel reads (frozen over the step). ------------!
      real(wp), allocatable :: nplant(:), bleaf(:), bsap(:), broot(:), sap_area(:), height(:), leaf_area(:)
   end type column_frozen_t

   !----- The whole-column tendency vector + diagnostics. ---------------------------------------!
   type :: column_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp
      real(wp) :: d_cas_shv      = 0.0_wp
      real(wp) :: d_cas_co2      = 0.0_wp
      real(wp) :: dedt(n_soil_layer_max)   = 0.0_wp   !< [W/m3] dsoil_energy/dt
      real(wp) :: dtheta_dt(n_soil_layer_max) = 0.0_wp!< [1/s]  dtheta/dt
      real(wp), allocatable :: dpsi_dt(:,:)           !< [MPa/s] (N_HYDRO, ncoh)
      real(wp) :: g_top = 0.0_wp, drainage_rate = 0.0_wp, uptake_rate = 0.0_wp
      real(wp), allocatable :: leaf_temp(:)
   end type column_tend_t

   !----- ARK conservation ledger: per-stage boundary-flux RATES (emitted by column_be_stage) and     !
   !      the b-weighted, cross-substep-ACCUMULATED amounts (the time-integral of the true boundary    !
   !      fluxes). Because Y3 - y = (1-gamma)*h*K2 + gamma*h*K3 exactly (BETA*gamma = 1-gamma), the     !
   !      accumulated in/out amounts telescope against the committed store change to machine precision  !
   !      for the flux-form CAS twins + the (energy_resid=0) soil-heat column -- the ARK path can then  !
   !      close the same 7 budgets the split closes. Reflects the CURRENT inert ARK (no soil-boundary   !
   !      water-enthalpy advection); the deferred precip>0 guard-lift adds those terms.                 !
   type :: stage_bflux_t                                    !< per-stage RATES
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp    !< [W/m2]
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp    !< [kg/m2/s]
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp    !< [umol/m2/s]
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp  !< [W/m2]
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp  !< [kg/m2/s]
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp!< [W/m2]
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp!< [kg/m2/s]
   end type stage_bflux_t

   type :: column_bflux_t                                  !< accumulated AMOUNTS (J/m2, kg/m2, umol/m2)
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp
   end type column_bflux_t

contains

   !----- Allocate a column_cohort_t (the per-patch cohort SoA the fast loop consumes). ------!
   subroutine alloc_column_cohort(coh, n)
      type(column_cohort_t), intent(out) :: coh
      integer(ik),           intent(in)  :: n
      coh%n = n
      allocate(coh%pft(n), coh%lai(n), coh%wai(n), coh%height(n), coh%crown(n),                &
               coh%leaf_width(n), coh%branch_diam(n), coh%leaf_area(n), coh%nplant(n),         &
               coh%dbh(n), coh%broot(n), coh%bleaf(n), coh%bsap(n), coh%sap_area(n),           &
               coh%vcmax25(n), coh%rd25(n))
      coh%pft = 1_ik
      coh%lai = 0.0_wp ; coh%wai = 0.0_wp ; coh%height = 0.0_wp ; coh%crown = 1.0_wp
      coh%leaf_width = 0.04_wp ; coh%branch_diam = 0.02_wp
      coh%leaf_area = 0.0_wp ; coh%nplant = 0.0_wp ; coh%dbh = 0.0_wp ; coh%broot = 0.0_wp
      coh%bleaf = 0.0_wp ; coh%bsap = 0.0_wp ; coh%sap_area = 0.0_wp
      coh%vcmax25 = 0.0_wp ; coh%rd25 = 0.0_wp
   end subroutine alloc_column_cohort

   !----- Flatten the shared [hydraulics] config into the plant hydro_params_t + rhizosphere      !
   !       conductance, and build the vulnerability lookup table from wood_kexp. The single seam    !
   !       between cfg%hydraulics (shared, TOML-driven) and the fast loop's hydro_params_t (plant),  !
   !       mirroring how the leaf seam flattens the PFT photosynthesis traits. -------------------!
   subroutine apply_hydraulics_config(hcfg, hydro_p, rhizo_cond)
      type(hydraulics_config_t), intent(in)    :: hcfg
      type(hydro_params_t),      intent(inout) :: hydro_p
      real(wp),                  intent(out)   :: rhizo_cond
      hydro_p%leaf_pi0       = hcfg%leaf_pi0       ; hydro_p%leaf_elastic_mod       = hcfg%leaf_elastic_mod
      hydro_p%leaf_apoplast_frac        = hcfg%leaf_apoplast_frac        ; hydro_p%leaf_water_sat = hcfg%leaf_water_sat
      hydro_p%wood_pi0       = hcfg%wood_pi0       ; hydro_p%wood_elastic_mod       = hcfg%wood_elastic_mod
      hydro_p%wood_apoplast_frac        = hcfg%wood_apoplast_frac        ; hydro_p%wood_water_sat = hcfg%wood_water_sat
      hydro_p%wood_psi50     = hcfg%wood_psi50     ; hydro_p%wood_kexp      = hcfg%wood_kexp
      hydro_p%k_plant_max    = hcfg%k_plant_max    ; hydro_p%wood_kmax      = hcfg%wood_kmax
      hydro_p%vessel_curl    = hcfg%vessel_curl
      call build_hydro_table(hydro_p%vuln_table, hydro_p%wood_kexp)
      rhizo_cond = hcfg%rhizo_cond
   end subroutine apply_hydraulics_config

end module meds_fast_types
