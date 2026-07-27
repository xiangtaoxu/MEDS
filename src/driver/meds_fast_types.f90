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
   use meds_biogeochem_types, only : co2_opts_t, n_soil_pool
   use meds_budget_check,     only : budget_t
   use meds_config,           only : hydraulics_config_t
   use meds_hydr_lib,         only : build_hydro_table
   implicit none
   private

   public :: LEAFEN_DIAGNOSTIC, LEAFEN_PROGNOSTIC, WOODEN_DIAGNOSTIC, WOODEN_PROGNOSTIC
   public :: SOILH2O_LAGGED, SOILH2O_COUPLED
   public :: column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: process_mask_t, mask_is_full
   public :: alloc_column_cohort, ensure_column_cohort_capacity, apply_hydraulics_config
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
   !----- The uniform PROCESS MASK (MEDS_NUMERICS_SCOPING.md §5.1). One switch set that every scheme    !
   !      honors, so the SAME driver can run a REDUCED column ODE: turning off energy leaves a water +   !
   !      CO2 column, turning off all but hydraulics isolates the stiffest mode. This is the process-     !
   !      complexity axis of the goal-(b) sweep, and it is scientifically useful on its own (attribute    !
   !      column behaviour to a process).                                                                 !
   !                                                                                          !
   !      SEMANTICS: a masked-OFF process is FROZEN -- its store does not evolve, so the ODE loses that   !
   !      dimension while the process still supplies its couplings to the others as a CONSTANT. The       !
   !      kernel is still invoked and its store restored afterwards (rather than the call being skipped)  !
   !      so no downstream consumer is ever handed an unset flux; the cost of a reduced system is         !
   !      therefore NOT lower, which matters when reading harness WORK metrics.                            !
   !                                                                                          !
   !      CONSERVATION: a reduced column is deliberately NOT closed -- freezing a store while its fluxes  !
   !      still act on its neighbours breaks the ledger by construction. mask_is_full() reports whether    !
   !      the budget halts are meaningful, and the driver suppresses them when they are not.               !
   !      All-true (the default) is the full column, so every existing path is byte-identical.             !
   type :: process_mask_t
      logical :: veg_energy = .true.   !< leaf + wood energy stores (prognostic modes)
      logical :: cas_energy = .true.   !< canopy-air-space enthalpy (temperature)
      logical :: cas_vapour = .true.   !< canopy-air-space specific humidity
      logical :: cas_co2    = .true.   !< canopy-air-space CO2 twin
      logical :: soil_heat  = .true.   !< soil thermal column
      logical :: soil_water = .true.   !< soil water (Richards) column
      logical :: hydraulics = .true.   !< plant hydraulics (psi)
   end type process_mask_t

   type :: column_config_t
      type(process_mask_t)        :: mask            !< process-complexity mask (all on = full column)
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
      !----- Advect liquid enthalpy on the interior per-layer Darcy flux. NOT optional -- it is a       !
      !      conservation requirement, not a tunable, and it was previously unreachable: the flag       !
      !      defaulted .false. and nothing anywhere set it. Water crossing a layer face carries its     !
      !      enthalpy with it; without this the top-face infiltration term deposits the whole of the    !
      !      infiltrating water's enthalpy in layer 1 while the water solver carries its MASS on down   !
      !      the column, so layer 1 keeps energy for water it no longer holds and the layers below      !
      !      gain mass carrying none (uext_to_temp then divides by a larger heat capacity). Measured    !
      !      at Ithaca before the fix: +1.23 K/h in layer 1 during infiltration hours against           !
      !      -0.075 K/h otherwise, layers 2-3 cooling up to 3.8 K/h at the same time, and soil-surface  !
      !      maxima landing at midnight after cloudy days. -----------------------------------------------!
      logical                     :: advect_soil_heat = .true.   !< advect liquid enthalpy on the interior
                                                                 !< per-layer Darcy flux (moisture<->energy coupling)
      !----- Canopy-surface water: interception film + film-evap/dew (MEDS_ED2_RK45_DESIGN.md sec 3.4, !
      !      P1) -- opt-in (default off, so existing configs are unchanged); SPLIT PATH ONLY for now,   !
      !      mirroring how snow (ccfg%snow_on) and prognostic leaf/wood energy both landed split-first  !
      !      with ARK support deferred (column_fast_step error-stops if this is on under INTEG_ARK). ---!
      logical                     :: canopy_water_on  = .false.
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
      !----- Soil-carbon matrix Rh diagnostics (B2, MEDS_SLOW_DYNAMICS_DESIGN.md Part II): filled     !
      !      by column_prepass ONLY when cfg%soil_carbon_on (else left at 0, matching the OLD          !
      !      constant-pool scalar path that runs instead). xi_step is this sub-step's per-pool          !
      !      environmental decomposition scalar [-]; the caller (fast_dynamics) accumulates             !
      !      xi_step*dt_fast_days into the per-patch day-integral xi_int, and rh_matrix_step*dt_fast_    !
      !      days into rh_fast_accum (design section 9's audit-only rh_seam_gap check). -----------------!
      real(wp)       :: xi_step(n_soil_pool) = 0.0_wp   !< [-] this sub-step's per-pool env scalar
      real(wp)       :: rh_matrix_step       = 0.0_wp   !< [kgC/m2/day] this sub-step's matrix Rh
      !----- P3 Picard diagnostics (reporting only; not conserved state). --------------------!
      integer(ik)    :: picard_iters       = 0_ik    !< worst outer-iteration count over the sub-steps
      integer(ik)    :: picard_nonconv     = 0_ik    !< number of sub-steps that hit picard_max_iter unconverged
      real(wp)       :: picard_worst_resid = 0.0_wp  !< [K] worst residual temperature at exit
      !----- WORK counters (MEDS_NUMERICS_SCOPING.md section 5.3). These are the COST axis of the      !
      !      benchmark: without them a sweep can report accuracy but not accuracy-per-unit-work, and   !
      !      wall-clock alone is too coarse and too machine-dependent to rank schemes. Every one of     !
      !      these was already computed somewhere and then discarded. Set PER SUB-STEP by the stepper;  !
      !      the fast driver accumulates them site-wide. -------------------------------------------!
      integer(ik)    :: integ_nsteps   = 0_ik   !< accepted ARK sub-steps this dt_fast (1 on the split path)
      integer(ik)    :: integ_nrej   = 0_ik   !< rejected integrator steps this dt_fast (0 on the split path)
      integer(ik)    :: soil_nsub    = 0_ik   !< soil-water Richards solver sub-steps
      integer(ik)    :: hydro_nsub   = 0_ik   !< plant-hydraulics sub-steps, summed over cohorts
      integer(ik)    :: hydro_nonconv = 0_ik  !< cohorts whose hydraulics solve did not converge
      !----- P6 (MEDS_ED2_RK45_DESIGN.md): count of sub-steps where the explicit RK45 step committed a   !
      !      railed (clamp-pinned, unphysical) CAS/soil state and the dispatcher rolled back + redid the   !
      !      step on the stable implicit-CAS split path. Rare (a handful over a healthy 30-yr run); a       !
      !      persistently-high value flags a genuinely stiff regime RK45 is degrading to split for. --------!
      integer(ik)    :: rk45_rescue  = 0_ik   !< dt_fast steps rescued RK45->split this sub-step (0 on split/ARK)
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
      !----- ADVECTIVE ENTHALPY (MEDS_ED2_RK45_DESIGN.md sec 2/6, P2): the water crossing the         !
      !      wood<->leaf and soil<->wood interfaces carries its own thermal energy (ED2's qwflux_wl/    !
      !      qloss) -- frozen (mass flux AND upwind reference temperature both fixed at state^n,         !
      !      build_column_frozen) per-cohort source terms folded directly into the diagnostic leaf/       !
      !      wood energy balance (surface_derivs), exactly like abs_sw/abs_sw_wood. Zero when unset       !
      !      (every existing caller/test fixture), so this is a no-op unless a caller populates it. ------!
      real(wp), allocatable :: qwflux_wl(:)   !< [W/m2 ground] sapflow's advected enthalpy INTO the leaf (wood->leaf)
      real(wp), allocatable :: q_wood_net(:)  !< [W/m2 ground] net advected enthalpy INTO wood (qloss - qwflux_wl)
      !----- Canopy-SURFACE water energy coupling (MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c): the wetted     !
      !      fraction is FROZEN once per dt_fast (the Act-1 pre-pass's intercept_canopy_layer sweep),      !
      !      mirroring every other frozen quantity in this tableau; only the STATE-dependent terms         !
      !      (dqdt, qsat_c-qcas) are re-evaluated per stage, exactly like the dry/stomatal g_tr_f pathway.  !
      !      g_film_f/w are the frozen boundary-layer-only (no stomatal resistance) film-evap conductances  !
      !      (leaf_film_coeff's result, sec 3.4/P1) feeding veg_energy_diagnostic's le_slope_wet/le_ref_wet !
      !      arguments; f_wet_c is its sigma_w output. All zero when canopy_water_on is off, so this is a    !
      !      no-op unless build_column_frozen populates it (mirrors qwflux_wl/q_wood_net above). ------------!
      real(wp), allocatable :: g_film_f(:), g_film_w(:)   !< [m/s] frozen film-evap conductance, leaf/wood
      real(wp), allocatable :: f_wet_c(:)                 !< [-]   frozen combined wetted fraction (sigma_w)
      real(wp) :: leaf_emiss    = 0.95_wp     !< [-]       leaf LW emissivity
      real(wp) :: wcap          = 0.0_wp      !< [kg/m2]   CAS mass capacity  -> enthalpy & vapour
      real(wp) :: ccap          = 0.0_wp      !< [mol/m2]  CAS molar capacity -> CO2
      real(wp) :: gah           = 0.0_wp      !< [kg/m2/s] CAS<->atm enthalpy conductance
      !----- SCHEME-ASYMMETRY GUARD (§8g). surface_derivs applies a smooth CAS supersaturation
      !      (condensation) sink, and surface_derivs is reached ONLY from the ARK stages -- the split
      !      stepper never calls it. So the two "schemes" have been integrating DIFFERENT MODELS, and
      !      every split-vs-ARK comparison conflated a physics term with a numerical method. This switch
      !      makes the term controllable so a like-for-like comparison is possible; .true. (default)
      !      preserves the historic ARK behaviour exactly. Whether the sink belongs on BOTH paths is a
      !      model question, deliberately left open here.
      logical  :: cas_condensation = .true.  !< apply the CAS supersaturation sink (ARK path only today)
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
      !----- Canopy-SURFACE water (sec 3.4, P2c): per-cohort film evaporation (dew if negative),        !
      !      exposed for column_derivs' surf-water ODE and ARK's per-stage b-weighted commit, mirroring   !
      !      transp_c's own role for the internal-water mass ODE. Zero when canopy_water_on is off. -------!
      real(wp), allocatable :: film_evap_leaf(:), film_evap_wood(:)   !< [kg/m2 ground/s]
   end type surface_tend_t

   !----- The full prognostic column state advanced per dt_fast. Plant hydraulics is represented   !
   !      NATIVELY as internal water MASS (MEDS_ED2_RK45_DESIGN.md sec 4, P2) -- not psi -- because   !
   !      mass's ODE is non-stiff (its inflow is a FROZEN constant and its outflow moves at the CAS    !
   !      timescale, sec 6), so it rides the same explicit stage machinery as CAS/soil with no        !
   !      operator split; psi is purely diagnostic (psi_from_water_content), read once per macro-step   !
   !      for the frozen gs pre-pass (column_prepass) and never advanced here. --------------------!
   type :: column_state_t
      real(wp) :: cas_enthalpy = 0.0_wp                    !< [J/kg]
      real(wp) :: cas_shv      = 0.0_wp                    !< [kg/kg]
      real(wp) :: cas_co2      = 0.0_wp                    !< [umol/mol]
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp   !< [J/m3]   per soil layer
      real(wp) :: theta(n_soil_layer_max)       = 0.0_wp   !< [m3/m3]  per soil layer
      real(wp), allocatable :: leaf_water_mass(:) !< [kg/plant] internal leaf water (ncoh)
      real(wp), allocatable :: wood_water_mass(:) !< [kg/plant] internal wood water (ncoh)
      !----- Canopy-SURFACE water (MEDS_ED2_RK45_DESIGN.md sec 3.4, P1+P2c): interception film on the   !
      !      leaf/wood boundary layer, DISTINCT from the internal (xylem/symplast) water above -- the     !
      !      surface film evaporates with no stomatal resistance, the internal water feeds transpiration   !
      !      through stomata. [kg/m2 GROUND] (already area-referenced, unlike the per-plant fields above,  !
      !      matching bio%leaf_surf_water/wood_surf_water's own convention from the split-path P1 landing). !
      real(wp), allocatable :: leaf_surf_water(:) !< [kg/m2 ground] canopy interception film on leaf
      real(wp), allocatable :: wood_surf_water(:) !< [kg/m2 ground] canopy interception film on wood
   end type column_state_t

   !----- Frozen inputs for the whole column: the surface pre-pass + the soil/hydraulics params +   !
   !      the frozen hydrology surface BCs + per-cohort geometry the hydraulics kernel needs.        !
   type :: column_frozen_t
      type(surface_frozen_t)      :: surf         !< the surface-block frozen inputs (t_ground overwritten per call)
      type(soil_params_t)         :: soil         !< soil geometry + texture (dz, root_frac, ...)
      type(soil_thermal_params_t) :: therm        !< soil thermal texture
      type(energy_opts_t)         :: energy_opts  !< soil-thermal options (phase change)
      type(soil_opts_t)           :: hydro_opts   !< soil-water (Richards) options
      real(wp) :: geothermal    = 0.0_wp          !< [W/m2]    bottom heat flux BC
      real(wp) :: q_top         = 0.0_wp          !< [m/s]     Richards top water flux (infiltration - evaporation)
      real(wp) :: soil_psi_root = 0.0_wp          !< [MPa]     root-zone soil water potential (hydraulics BC;
                                                  !<           DIAGNOSED from state^n theta, sec 3/5 -- the
                                                  !<           Act-1 pre-pass runs hydraulics BEFORE the soil
                                                  !<           solve, so this is no longer the scratch solve's
                                                  !<           own post-solve psi_soil)
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
      !----- realized root uptake [kg/m2/s]: the Act-1 pre-pass's plant-side REQUEST (total_uptake_b,     !
      !      sec 3), rescaled by the soil's OWN fwilt-limited supply (scale = uptake/requested <= 1) --   !
      !      the SAME number both the soil-water tendency's root sink (column_derivs) and the per-cohort   !
      !      uptake_frozen below are built from, so the wood<->soil interface closes to the soil's TRUE     !
      !      realized supply (mirrors the split path's own treatment, sec 3/5). ------------------------!
      real(wp) :: uptake        = 0.0_wp          !< [kg/m2/s] realized (post-rescale) aggregate root uptake
      !----- the AUTHORITATIVE end-of-step soil moisture from the scratch column_hydrology_flux (the robust  !
      !      ponding/runoff/free-drain Richards solve). The ARK COMMITS this instead of re-solving theta in   !
      !      the ESDIRK stages (soil water is fully operator-split out; see column_fast_step_ark).            !
      real(wp), allocatable :: theta1(:)          !< [m3/m3]   committed post-step soil moisture (per layer)
      real(wp), allocatable :: psi_e(:)           !< [m]       Zeng-Decker equilibrium potential per layer (frozen)
      !----- per-cohort geometry the hydraulics kernel reads (frozen over the step). ------------!
      real(wp), allocatable :: nplant(:), bleaf(:), bsap(:), broot(:), sap_area(:), height(:), leaf_area(:)
      !----- FROZEN plant-hydraulics fluxes (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2): the Act-1 pre-pass's  !
      !      time-averaged solve_plant_water output (per plant), held CONSTANT across every sub-stage of    !
      !      the macro-step -- sapflow_frozen is the wood->leaf transfer; uptake_frozen is the soil->wood    !
      !      transfer, ALREADY floored >=0 (no hydraulic redistribution, matching the project-wide           !
      !      convention) and rescaled by `scale` (uptake/requested) so sum(uptake_frozen*nplant) == uptake    !
      !      above EXACTLY -- one number used on both sides of the wood<->soil interface, the closure         !
      !      principle sec 3.2 specifies. column_derivs' mass ODE reads these directly; no PV-curve/          !
      !      conductance evaluation is needed per stage any more (that algebra lives ONLY in the pre-pass).    !
      real(wp), allocatable :: sapflow_frozen(:), uptake_frozen(:)   !< [kg/plant/s] (ncoh)
      !----- FROZEN advective enthalpy leaving the soil via root uptake (MEDS_ED2_RK45_DESIGN.md sec   !
      !      2/6, P2 -- ED2's qloss): uptake_frozen(i)*nplant(i) converted to per-ground-area, times      !
      !      the root-frac-weighted state^n soil temperature's liquid internal energy -- frozen ONCE      !
      !      in the Act-1 pre-pass alongside sapflow_frozen/uptake_frozen. column_derivs debits this        !
      !      from the soil-heat root_heat_sink (the same interface fro%surf%q_wood_net's wood credit        !
      !      pairs with, sec 2's qloss - qwflux_wl). -------------------------------------------------!
      real(wp), allocatable :: qloss_frozen(:)   !< [W/m2 ground] (ncoh)
      !----- FROZEN canopy interception (MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c): the Act-1 pre-pass's     !
      !      ONE height-sorted intercept_canopy_layer sweep (e_canopy=0, capture/capacity only), held     !
      !      CONSTANT across the whole macro-step -- mirrors sapflow_frozen/uptake_frozen's own            !
      !      "one frozen number, no per-stage re-solve" convention (sec 6 stability argument): re-running   !
      !      a capacity-limited bucket per RK/ESDIRK stage would need a stage-local dt, not dt_fast, and     !
      !      no precedent elsewhere in this tableau does that. Zero when canopy_water_on is off (every       !
      !      existing caller/test fixture), so this is a no-op unless build_column_frozen populates it. -----!
      real(wp), allocatable :: intercept_leaf(:), intercept_wood(:)   !< [kg/m2 ground/s] (ncoh)
   end type column_frozen_t

   !----- The whole-column tendency vector + diagnostics. ---------------------------------------!
   type :: column_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp
      real(wp) :: d_cas_shv      = 0.0_wp
      real(wp) :: d_cas_co2      = 0.0_wp
      real(wp) :: dedt(n_soil_layer_max)   = 0.0_wp   !< [W/m3] dsoil_energy/dt
      real(wp) :: dtheta_dt(n_soil_layer_max) = 0.0_wp!< [1/s]  dtheta/dt
      real(wp), allocatable :: d_leaf_water_mass(:)   !< [kg/plant/s] frozen_sapflow - transp(stage)
      real(wp), allocatable :: d_wood_water_mass(:)   !< [kg/plant/s] frozen_uptake  - frozen_sapflow
      !----- Canopy-SURFACE water (sec 3.4, P2c): frozen_intercept - film_evap(stage), the surface-film  !
      !      analogue of d_leaf_water_mass/d_wood_water_mass above. --------------------------------------!
      real(wp), allocatable :: d_leaf_surf_water(:)   !< [kg/m2 ground/s] frozen_intercept_leaf - film_evap
      real(wp), allocatable :: d_wood_surf_water(:)   !< [kg/m2 ground/s] frozen_intercept_wood - film_evap
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

   !----- Is the column the FULL system? Only then are the closed-budget halts meaningful (a frozen  !
   !      store still exchanges with its neighbours, so a reduced column cannot conserve). -----------!
   pure logical function mask_is_full(m)
      type(process_mask_t), intent(in) :: m
      mask_is_full = m%veg_energy .and. m%cas_energy .and. m%cas_vapour .and. m%cas_co2 .and.        &
                     m%soil_heat  .and. m%soil_water .and. m%hydraulics
   end function mask_is_full

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

   !---------------------------------------------------------------------------------------!
   ! Grow-only capacity check: reallocate ONLY when the current backing arrays are too small !
   ! for n cohorts; otherwise just update the active count coh%n and leave the (over-sized)    !
   ! capacity in place (MEDS_NUMERICS_SCOPING.md BB1 phase 1). Every downstream reader loops    !
   ! over coh%n (never size(coh%pft) -- verified: no caller does), and the per-patch gather      !
   ! that follows a call to this routine overwrites indices 1..n unconditionally, so reusing a   !
   ! larger patch's leftover capacity for a smaller patch is bit-identical: the caller sizes      !
   ! coh ONCE per fast_dynamics call (to the site-wide max cohort count) instead of once PER      !
   ! PATCH, cutting O(n_patch) heap allocations to O(1) per slow step.                            !
   !---------------------------------------------------------------------------------------!
   subroutine ensure_column_cohort_capacity(coh, n)
      type(column_cohort_t), intent(inout) :: coh
      integer(ik),            intent(in)    :: n
      if (.not. allocated(coh%pft)) then
         call alloc_column_cohort(coh, n)
      else if (size(coh%pft) < n) then
         call alloc_column_cohort(coh, n)
      else
         coh%n = n
      end if
   end subroutine ensure_column_cohort_capacity

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
