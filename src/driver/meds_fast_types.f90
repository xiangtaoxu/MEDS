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
      real(wp)                    :: specific_root_area = 20.0_wp  !< [m2/kgC] SRA (rhizosphere conductance)
      real(wp)                    :: fast_soil_carbon = 5.0_wp   !< [kgC/m2] decomposable soil-C pool (prescribed, MVP)
      !----- Canopy-surface water: interception film + film-evap/dew (MEDS_ED2_RK45_DESIGN.md sec 3.4, !
      !      P1) -- opt-in (default off, so existing configs are unchanged); SPLIT PATH ONLY for now,   !
      !      mirroring how snow (ccfg%snow_on) and prognostic leaf/wood energy both landed split-first  !
      !      with ARK support deferred (column_fast_step error-stops if this is on under INTEG_ARK). ---!
      logical                     :: canopy_water_on  = .false.
      !----- Inner-solve knobs (retained: ARK's newton_surface_solve uses the iteration cap). ----!
      type(snow_params_t) :: snow                    !< snow parameters (density, albedo, thresholds, conductivity)
      integer(ik) :: picard_max_iter = 20_ik        !< outer-iteration cap
      real(wp)    :: picard_tol_temp = 1.0e-3_wp     !< [K]     temperature convergence tolerance
      real(wp)    :: picard_tol_shv  = 1.0e-6_wp     !< [kg/kg] CAS specific-humidity convergence tolerance
      real(wp)    :: picard_relax    = 0.5_wp        !< [-]     under-relaxation of the next-pass seed. The CAS<->ground
                                                     !<         sensible coupling gives the fixed-point map a slope ~ -1
                                                     !<         (oscillatory); 0.5 makes the relaxed map a contraction.
      logical     :: picard_fixed_iter = .false.     !< run a uniform pass count (GPU warp-uniform; no early exit)
      integer(ik) :: leaf_energy_model  = 0_ik       !< LEAFEN_DIAGNOSTIC (0) | LEAFEN_PROGNOSTIC (1)
      integer(ik) :: wood_energy_model  = 0_ik       !< WOODEN_DIAGNOSTIC (0) | WOODEN_PROGNOSTIC (1)
   end type column_config_t

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into cfg%pft)
      real(wp),    allocatable :: lai(:), wai(:), height(:), crown(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   !< [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              !< [kgC/plant],[kgC/plant],[m2] (hydraulics)
      !----- TOTAL wood carbon, distinct from bsap and NOT interchangeable with it. bsap is the       !
      !      sapwood ring, which is the right quantity for the HYDRAULIC capacitance; the THERMAL     !
      !      store takes all the wood, because branch wood is thin enough to be thermally active      !
      !      throughout and a sapwood fraction defined on the bole systematically under-counts it     !
      !      (measured ~6x on a 70 cm cohort). One pool, one temperature -- a bole/branch partition   !
      !      would mean two of each and is deliberately not taken. ---------------------------------!
      real(wp),    allocatable :: bwood(:)                                    !< [kgC/plant] total wood (thermal store)
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
      !----- CLAMP activations (MEDS_INTEGRATOR_PARITY.md, Phase A). The stability clamps            !
      !      (clamp_theta / clamp_cas / clamp_soil_energy) are the one place where a scheme edits      !
      !      state outside the conservation ledger, and they are TRAJECTORY-dependent -- ifx and       !
      !      nvfortran do not fire them on the same steps, which is why a compiler-split test          !
      !      failure was the first symptom. Counting them separates two very different events:         !
      !                                                                                                !
      !        STAGE clamps bound a THROWAWAY stage input so the RHS stays evaluable. Harmless in      !
      !        itself; a rising count means "this dt is too big", and the step is normally rejected.   !
      !                                                                                                !
      !        COMMIT clamps edit the state that is actually kept. THESE fabricate mass/energy with    !
      !        no ledger term (clamp_theta moves theta with no mass debit; clamp_soil_energy then      !
      !        re-derives temperature at the new water mass), and at the sub-step FLOOR the accept     !
      !        branch takes the clamped state unconditionally -- the controller cannot veto it.        !
      !                                                                                                !
      !      Diluting the two into one counter would hide exactly the signal being sought, so they     !
      !      are separate. Magnitudes are accumulated for the two stores whose books they break;       !
      !      the CAS clamp is counted only (rk45_state_railed already covers CAS railing). ARK clamps  !
      !      its ARS stage-3 extrapolation base ONLY, so its commit counters stay 0 by construction.   !
      !                                                                                                !
      !      READ THE MAGNITUDE, NOT THE COUNT. Measured on test_column_rk45: the commit clamp fires    !
      !      on essentially every sub-step even in a well-behaved wet window (~1250 activations / 96    !
      !      steps), because the Richards solve routinely lands a whisker outside                       !
      !      [theta_res, theta_sat]. The count therefore barely separates a healthy column (~1250)      !
      !      from a saturated, railing one (~2840). The magnitudes separate them by SEVEN ORDERS        !
      !      (3e-5 vs 1.4e2 kg/m2). So the counts are cheap telemetry; clamp_mass/clamp_energy are      !
      !      the metric.                                                                                !
      !                                                                                                !
      !      Both magnitudes are GROSS: a running sum of |correction|, not a net ledger residual.       !
      !      Corrections of opposite sign do not cancel here (deliberately -- two large corrections     !
      !      that happen to cancel are still two moments when the state left its own model), so this    !
      !      is an upper bound on, not an estimate of, the resulting budget gap. -----------------------!
      integer(ik)    :: clamp_stage_n  = 0_ik   !< stage-input clamp activations (layers/stores clamped)
      integer(ik)    :: clamp_commit_n = 0_ik   !< COMMITTED-state clamp activations -- unbookkept
      real(wp)       :: clamp_mass     = 0.0_wp !< [kg/m2] sum |water| moved by a COMMIT clamp_theta
      real(wp)       :: clamp_energy   = 0.0_wp !< [J/m2]  sum |energy| moved by a COMMIT clamp_soil_energy
      !----- CONSTITUTIVE-DOMAIN excursion of theta (issue #78 item 2), max over the sub-steps of one   !
      !      dt_fast, in m3/m3. Complements clamp_stage_n rather than duplicating it, in two ways: it   !
      !      is a MAGNITUDE where that is a count, and it covers the RK45 stage-1 evaluation, which     !
      !      clamp_stage_n cannot see at all -- k1 reads the previous sub-step's committed state, and   !
      !      C1 deliberately stopped clamping what is committed, so no clamp fires there to be counted. !
      !                                                                                                !
      !      This is deliberately telemetry and not a correction. The constitutive kernels all clamp    !
      !      their own effective saturation (see test_rhs_domain_safety), so an excursion is harmless   !
      !      to evaluate -- an oversaturated cell is treated as exactly saturated, which is the right   !
      !      answer for one. What the number is FOR is noticing if that excursion ever stops being      !
      !      small: measured worst case is 7.9e-3 in theta (Se = 1.022) on a sealed 29 mm/h fixture,    !
      !      and identically 0 on a month-long forced Ithaca cell. Zero on the split and ARK paths,     !
      !      which have no explicit stages to overshoot in. --------------------------------------------!
      real(wp)       :: theta_ood_max  = 0.0_wp !< [m3/m3] max excursion outside [theta_res, theta_sat]
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
      !----- TISSUE HEAT STORE, frozen per dt_fast (Category 0). a_* = cap/dt_fast is the storage      !
      !      conductance veg_energy_diagnostic relaxes against; t_*0 is the start-of-step temperature   !
      !      it relaxes FROM. Both frozen for the whole fast step, like every other coefficient here,   !
      !      so each stage evaluation returns the SAME dt_fast-averaged flux and dt_fast-endpoint       !
      !      temperature. The store is deliberately NOT a tableau degree of freedom: it is an algebraic !
      !      closure evaluated at each stage, which is why it needs no new WRMS group, no arrowhead and !
      !      no Newton -- see MEDS_VEG_ENERGY_INTEGRATION_PLAN.md sec 2. --------------------------------!
      real(wp), allocatable :: a_leaf(:), a_wood(:)   !< [W/m2/K] cap/dt_fast
      real(wp), allocatable :: t_leaf0(:), t_wood0(:) !< [K]      start-of-step tissue temperatures
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
      !----- SHARED SNOW STAGE outputs (C4, issue #76). Frozen once per dt_fast by
      !      meds_fast_snow%advance_snow_stage and consumed by surface_derivs' ground blend.
      !      ALL DEFAULT TO ZERO, and the blend is written so zeros reduce it EXACTLY to the
      !      snow-free form -- that is what makes snow-off bit-identical structural rather than
      !      something each scheme has to re-verify.
      real(wp) :: snowfac    = 0.0_wp   !< [-]       snow cover fraction (0 = bare ground)
      real(wp) :: h_snow     = 0.0_wp   !< [W/m2]    snowfac-weighted sensible flux to the CAS
      real(wp) :: le_snow    = 0.0_wp   !< [W/m2]    snowfac-weighted latent (sublimation) flux
      real(wp) :: g_base_snow = 0.0_wp  !< [W/m2]    throttled base conduction into the soil top
      real(wp) :: subl_rate  = 0.0_wp   !< [kg/m2/s] sublimation vapour source for the CAS
      real(wp) :: ground_rad = 0.0_wp   !< [W/m2]    blended ground radiative input (= abs_sw+abs_lw when bare)
      !----- snow STORE + boundary terms the whole-column ledgers need (C4). All 0 without snow. ---!
      real(wp) :: snow_swe0  = 0.0_wp   !< [kg/m2]   pack mass BEFORE the stage
      real(wp) :: snow_swe1  = 0.0_wp   !< [kg/m2]   pack mass AFTER  the stage
      real(wp) :: snow_enth0 = 0.0_wp   !< [J/m2]    pack internal energy BEFORE
      real(wp) :: snow_enth1 = 0.0_wp   !< [J/m2]    pack internal energy AFTER
      real(wp) :: snow_acc_enth = 0.0_wp!< [J/m2]    precip enthalpy that entered the pack (boundary in)
      real(wp) :: snow_melt_enth= 0.0_wp!< [J/m2]    melt enthalpy pack -> POND (reported; no baseline rebase now)
      real(wp) :: snow_t_melt   = 0.0_wp!< [K]       temperature that values the meltwater (#78 item 4)
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
      !----- frozen boundary hydrology for the precip>0 guard-lift: the throughfall/drainage/runoff    !
      !      water carries internal_energy_liquid across the soil boundaries (matches the split's       !
      !      :436-439,518-520 advection), and the scratch column_hydrology_flux's end-of-step ponding/  !
      !      aquifer/water-table is persisted (column_state_t does NOT carry these surface stores). ----!
      real(wp) :: infiltration  = 0.0_wp          !< [kg/m2/s] throughfall reaching the soil top face
      real(wp) :: drainage      = 0.0_wp          !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf   = 0.0_wp          !< [kg/m2/s] surface runoff
      !----- ground water input the hydrology saw [kg/m2/s]. Needed by RK45 to rebuild its OWN     !
      !      ponding store from its own trajectory rather than inheriting the scratch solve's      !
      !      end-of-step pond (which already contains the SCRATCH clip -- see issue #75). ---------!
      real(wp) :: precip_ground = 0.0_wp          !< [kg/m2/s] water reaching the ground
      real(wp) :: t_infil     = 0.0_wp   !< [K] temperature of the infiltrating (pond) water, #78 item 4
      real(wp) :: w_surface_enth1 = 0.0_wp !< [J/m2] scratch solve's end-of-step pond ENTHALPY, #78 item 4
      real(wp) :: t_precip = 0.0_wp      !< [K] temperature of the water entering the pond, #78 item 4
      real(wp) :: runoff_enth = 0.0_wp   !< [W/m2] enthalpy leaving with surface runoff, #78 item 4
      real(wp) :: rain_temp     = 0.0_wp          !< [K]       rain temperature (CAS temp @ state^n)
      real(wp) :: t_bot         = 0.0_wp          !< [K]       bottom-layer soil temperature @ state^n
      real(wp) :: w_surface1    = 0.0_wp          !< [kg/m2]   end-of-step ponded surface water
      !----- realized root uptake [kg/m2/s]: the Act-1 pre-pass's plant-side REQUEST (total_uptake_b,     !
      !      sec 3), rescaled by the soil's OWN fwilt-limited supply (scale = uptake/requested <= 1) --   !
      !      the SAME number both the soil-water tendency's root sink (column_derivs) and the per-cohort   !
      !      uptake_frozen below are built from, so the wood<->soil interface closes to the soil's TRUE     !
      !      realized supply (mirrors the split path's own treatment, sec 3/5). ------------------------!
      real(wp) :: uptake        = 0.0_wp          !< [kg/m2/s] realized (post-rescale) aggregate root uptake
      !----- FROZEN interior face fluxes + post-solve mass corrections from the SAME Act-1 scratch      !
      !      column_hydrology_flux that produced theta1. Freezing them is not an approximation the ARK   !
      !      pays extra for: hflux%w_flux is itself the solver's TIME-MEAN face flux over the step, so    !
      !      this is exactly what the split path advects on (eforc%w_flux = -hflux%w_flux).               !
      !                                                                                                  !
      !      w_flux_frozen matters for CONSERVATION, not accuracy. The boundary faces put liquid          !
      !      enthalpy into layer 1 and take it out at layer nsl / wherever the clip fires; with the        !
      !      interior faces zeroed (as ARK and RK45 both had them) there is no path between them, so       !
      !      layer 1 accumulates the whole infiltration enthalpy while a deeper layer sheds it. The        !
      !      whole-column ledger still closes -- the error is purely vertical, which a column-vs-boundary  !
      !      sum cannot see. On the split path the same defect drives a soil surface to 361 K              !
      !      (test_column_dynamics RUN 7).                                                                !
      !                                                                                                  !
      !      THESE THREE ARE ARK-ONLY (issue #78 item 3). They are the right numbers for a scheme that     !
      !      commits the scratch solve's theta VERBATIM, which the ARK does (soil water is operator-split  !
      !      out of its stages). RK45 integrates its OWN theta, on which the scratch's faces move a        !
      !      different amount of water and the scratch's clip mass never moves at all, so column_derivs    !
      !      now takes both from the stage's own soil_water_time_deriv. Using these there cost ~2.6e6      !
      !      J/m2/step of vertical enthalpy misplacement (soil surface 345 K) against ~2.6e6 J/m2/step of  !
      !      spurious clip cooling -- two defects of matched magnitude and opposite sign, which is why     !
      !      each hid the other and why removing either one alone made the fixture worse. ----------------!
      real(wp) :: w_flux_frozen(n_soil_layer_max) = 0.0_wp  !< [m/s]   DOWNWARD interior face flux, k=1..nsl-1
      !----- Enthalpy paired with the hydrology's UNFACED post-solve mass corrections, already valued  !
      !      at each layer's own state^n temperature (so the correction is temperature-NEUTRAL) and     !
      !      carried as [W/m2] to join the root_heat_sink column the stages already assemble. ---------!
      real(wp) :: clip_enth(n_soil_layer_max)  = 0.0_wp     !< [W/m2] enthalpy leaving layer k with clipped water
      real(wp) :: floor_enth(n_soil_layer_max) = 0.0_wp     !< [W/m2] enthalpy created with theta_res-floored water
      !----- the AUTHORITATIVE end-of-step soil moisture from the scratch column_hydrology_flux (the robust  !
      !      ponding/runoff/free-drain Richards solve). The ARK COMMITS this instead of re-solving theta in   !
      !      the ESDIRK stages (soil water is fully operator-split out; see column_fast_step_ark).            !
      real(wp), allocatable :: theta1(:)          !< [m3/m3]   committed post-step soil moisture (per layer)
      !----- Per-layer root-sink placement (Phase 1, MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md): THIS       !
      !      dt_fast's realized per-layer uptake shares, normalized to 1, built in build_column_frozen   !
      !      from solve_plant_water_batch's own breakdown. Placed identically to the split path (which   !
      !      builds the same array inline), so all three schemes put the root MASS sink and the root      !
      !      HEAT sink in the same layers by construction. Falls back to the static root_frac profile     !
      !      when no layer supplies anything. ---------------------------------------------------------!
      real(wp), allocatable :: root_share(:)      !< [-]       per-layer root-sink shares (sum = 1)
      !----- PROGNOSTIC WOOD (Phase 4). Populated by build_column_frozen ONLY when                    !
      !      wood_energy_model == WOODEN_PROGNOSTIC, in which case fro%surf's DIAGNOSTIC wood inputs   !
      !      are zeroed instead, so surface_derivs contributes no wood term and there is exactly one   !
      !      wood authority per run. Advanced operator-split by advance_wood_energy_full at the        !
      !      COMMITTED CAS endpoint -- wood is stiff (measured 55-199 s vs dt_fast = 1800 s, see       !
      !      test_wood_stiffness_spread), so it rides the L-stable veg_energy_step_implicit kernel     !
      !      rather than either tableau. -------------------------------------------------------------!
      real(wp), allocatable :: wood_dry_hcap(:)   !< [J/m2/K]  dry sapwood heat capacity (floored)
      real(wp), allocatable :: wood_wmass(:)      !< [kg/m2]   fresh-sapwood water mass
      !----- LEAF twin of the two above, for the post-commit tissue store (veg_store_correction).    !
      !      The leaf's own tau is ~12.5 s -- far inside a dt_fast -- so its correction is nearly a   !
      !      no-op at production cadence; it is carried anyway so the leaf and wood stores are ONE    !
      !      mechanism rather than a wood special case, and so the store stays right as dt_fast       !
      !      shrinks (the sub-daily probe and any future short-step run). --------------------------!
      real(wp), allocatable :: leaf_dry_hcap(:)   !< [J/m2/K]  dry leaf heat capacity (floored)
      real(wp), allocatable :: leaf_wmass(:)      !< [kg/m2]   internal (symplast) leaf water mass
      real(wp), allocatable :: wood_gbh(:)        !< [m/s]     wood boundary-layer conductance
      real(wp), allocatable :: wood_abs_sw(:)     !< [W/m2]    absorbed SW on wood
      real(wp), allocatable :: wood_abs_lw(:)     !< [W/m2]    net LW on wood
      real(wp), allocatable :: wood_area(:)       !< [m2/m2]   cohort wood area index
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
      real(wp) :: whole_cond    = 0.0_wp                         !< [kg/m2/s] condensate (row 1b)
   end type stage_bflux_t

   type :: column_bflux_t                                  !< accumulated AMOUNTS (J/m2, kg/m2, umol/m2)
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp
      !----- CONDENSATE, tracked SEPARATELY from whole_wat_out (row 1b). Dew/fog is not a boundary     !
      !      loss -- it lands on a surface inside the column -- but it used to be summed into           !
      !      whole_wat_out with the atmospheric vapour flux, where it could not be told apart or        !
      !      redirected. Carrying it in its own slot is what lets the caller deposit it into a store.   !
      real(wp) :: whole_cond    = 0.0_wp   !< [kg/m2] condensed vapour over the step (>= 0)
      !----- TISSUE-TEMPERATURE TIME INTEGRALS [K*s], per cohort, b-weighted across stages and summed !
      !      over accepted sub-steps. These are what make the tissue store conserve EXACTLY on an      !
      !      adaptive scheme, and they are also the physically right answer rather than merely the     !
      !      conserving one.                                                                          !
      !                                                                                          !
      !      The frozen-store kernel's own balance is cap*(T_end - T_0)/dt_fast = numer - denom*dT_avg,!
      !      so the store's gain RATE at a stage is a_store*(sf%leaf_temp - t_leaf0) -- a plain        !
      !      function of what surface_derivs already returns. The canopy air receives the b-weighted   !
      !      denom*dT_avg, so the store must gain the b-weighted complement, i.e. the store's energy   !
      !      is set by the TIME INTEGRAL of the tissue temperature, not by its final-stage value.      !
      !      Committing T = integral/dt_fast makes the state and the ledger the SAME number.           !
      !                                                                                          !
      !      With a CAS that is constant over the step every stage returns the same temperature and    !
      !      the integral collapses to it, so this degenerates correctly. -----------------------------!
      real(wp), allocatable :: tissue_leaf_int(:), tissue_wood_int(:)   !< [K*s]
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
               coh%bwood(n), coh%vcmax25(n), coh%rd25(n))
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
   subroutine apply_hydraulics_config(hcfg, hydro_p)
      type(hydraulics_config_t), intent(in)    :: hcfg
      type(hydro_params_t),      intent(inout) :: hydro_p
      hydro_p%leaf_pi0       = hcfg%leaf_pi0       ; hydro_p%leaf_elastic_mod       = hcfg%leaf_elastic_mod
      hydro_p%leaf_apoplast_frac        = hcfg%leaf_apoplast_frac        ; hydro_p%leaf_water_sat = hcfg%leaf_water_sat
      hydro_p%wood_pi0       = hcfg%wood_pi0       ; hydro_p%wood_elastic_mod       = hcfg%wood_elastic_mod
      hydro_p%wood_apoplast_frac        = hcfg%wood_apoplast_frac        ; hydro_p%wood_water_sat = hcfg%wood_water_sat
      hydro_p%wood_psi50     = hcfg%wood_psi50     ; hydro_p%wood_kexp      = hcfg%wood_kexp
      hydro_p%k_plant_max    = hcfg%k_plant_max    ; hydro_p%wood_kmax      = hcfg%wood_kmax
      hydro_p%vessel_curl    = hcfg%vessel_curl
      call build_hydro_table(hydro_p%vuln_table, hydro_p%wood_kexp)
   end subroutine apply_hydraulics_config

end module meds_fast_types
