!==========================================================================================!
! meds_fast_types -- the shared driver-scope TYPES of the fast-loop: the four working-buffer     !
! bundles both integrators pass through the dispatch seam (column_config_t/column_cohort_t/       !
! column_forcing_t/column_budget_t, from the former meds_column_dynamics) and the ARK POD state /  !
! frozen-input / tendency / boundary-flux-ledger types the whole-column RHS advances (surface_*_t, !
! column_state_t, column_frozen_t, column_tend_t, stage_bflux_t, column_bflux_t, from the former    !
! meds_column_derivs). Extracting them here is the prerequisite plumbing that lets meds_fast_ark    !
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
   use meds_kinds, only : wp, ik
   use meds_plant_types, only : veg_thermal_params_t
   use meds_column_params, only : n_soil_layer_max, soil_params_t, soil_thermal_params_t
   use meds_column_state_types, only : cas_state_t, soil_column_t, soil_energy_column_t, snow_column_t, soil_carbon_t
   use meds_therm_lib, only : cas_enthalpy_of_temp
   use meds_biophysics_opts, only : aero_cfg_t, soil_opts_t, energy_opts_t, snow_params_t
   use meds_plant_types, only : wood_params_t, root_params_t, hydro_params_t, hydro_opts_t, leaf_photo_table_t
   use meds_biogeochem_types, only : co2_opts_t, n_soil_pool
   use meds_budget_check, only : budget_t
   use meds_config, only : hydraulics_config_t, INTEG_ARK, CTRL_L1_ADAPTIVE, CTRL_I
   use meds_hydr_lib, only : build_hydro_table
   use meds_site_state_types, only : DMAX_PSI_LEAF_UNSET
   implicit none
   private

   public :: column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W, GRP_THETA, GRP_SOIL_T, N_TOL_GROUP
   public :: tol_set_t, error_control_t, integrator_opts_t
   public :: process_mask_t, mask_is_full
   public :: alloc_column_cohort, ensure_column_cohort_capacity, apply_hydraulics_config
   public :: surface_state_t, surface_tend_t
   public :: patch_biophys_t, alloc_patch_biophys, ensure_patch_biophys_capacity
   public :: snow_stage_t
   public :: cas_boundary_t, tissue_coefficients_t, canopy_film_capacity_t, ground_boundary_t
   public :: soil_hydrology_t, root_zone_t, plant_water_t, column_params_t
   public :: column_state_t, column_frozen_t, column_tend_t
   public :: stage_bflux_t, column_bflux_t

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

   !----- The tolerance GROUPS -- one per physical field class in the fast-loop state. Groups 1-6 are   !
   !      the INTEGRATED state (what the embedded-error WRMS measures); groups 7-8 belong to the nested   !
   !      SUB-SOLVERS (soil-water Richards on theta, soil-energy on temperature) that the driver drives   !
   !      from this same set -- so ONE tolerance source governs the whole hierarchy (§8c Layer 1). -----!
   integer(ik), parameter :: GRP_ENTH    = 1_ik   !< CAS specific enthalpy   [J/kg]
   integer(ik), parameter :: GRP_SHV     = 2_ik   !< CAS specific humidity   [kg/kg]
   integer(ik), parameter :: GRP_CO2     = 3_ik   !< CAS CO2 mole fraction   [umol/mol]
   integer(ik), parameter :: GRP_SE      = 4_ik   !< soil internal energy    [J/m3]
   integer(ik), parameter :: GRP_LEAF_W  = 5_ik   !< leaf internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_WOOD_W  = 6_ik   !< wood internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_THETA   = 7_ik   !< soil moisture           [m3/m3] (soil-water sub-solver)
   integer(ik), parameter :: GRP_SOIL_T  = 8_ik   !< soil temperature        [K]     (soil-energy sub-solver)
   integer(ik), parameter :: N_TOL_GROUP = 8_ik

   !----- Historical per-field absolute tolerances, used as the group defaults so every path is       !
   !      byte-identical unless overridden. ----------------------------------------------------------!
   real(wp), parameter :: ATOL_ENTH_DEF   = 5.0e1_wp    !< [J/kg]      (~0.05 K in enthalpy)
   real(wp), parameter :: ATOL_SHV_DEF    = 1.0e-6_wp   !< [kg/kg]
   real(wp), parameter :: ATOL_CO2_DEF    = 1.0e-1_wp   !< [umol/mol]
   real(wp), parameter :: ATOL_SE_DEF     = 1.0e3_wp    !< [J/m3]
   real(wp), parameter :: ATOL_LEAF_W_DEF = 1.0e-4_wp   !< [kg/plant]
   real(wp), parameter :: ATOL_WOOD_W_DEF = 1.0e-4_wp   !< [kg/plant]
   real(wp), parameter :: ATOL_THETA_DEF  = 1.0e-4_wp   !< [m3/m3] (== soil_opts_t's own default)
   real(wp), parameter :: ATOL_SOIL_T_DEF = 1.0e-2_wp   !< [K]     (== energy_opts_t's own default)
   !----- Default PI gains for a 1st-order embedded pair (Gustafsson 1988 / Soderlind): a = 0.7/2,   !
   !      b = 0.4/2. fac = safety*err^-a*err_prev^b; b = 0 recovers a pure I-controller. -----------!
   real(wp), parameter :: PI_ALPHA_DEF = 0.35_wp
   real(wp), parameter :: PI_BETA_DEF  = 0.20_wp

   !----- Per-group (rtol, atol). The WRMS normalizes state group g by atol(g) + rtol(g)*|y|. ---------!
   type :: tol_set_t
      real(wp) :: rtol(N_TOL_GROUP) = 1.0e-3_wp
      real(wp) :: atol(N_TOL_GROUP) = [ATOL_ENTH_DEF, ATOL_SHV_DEF, ATOL_CO2_DEF, ATOL_SE_DEF,   &
                                       ATOL_LEAF_W_DEF, ATOL_WOOD_W_DEF, ATOL_THETA_DEF, ATOL_SOIL_T_DEF]
   end type tol_set_t

   !----- The bundle threaded into an adaptive march: strictness + controller + step-clamp knobs +     !
   !      PI gains + the tolerance set. -----------------------------------------------------------------!
   type :: error_control_t
      integer(ik)      :: level      = CTRL_L1_ADAPTIVE
      integer(ik)      :: controller = CTRL_I
      real(wp)         :: safety     = 0.9_wp
      real(wp)         :: fmin       = 0.2_wp
      real(wp)         :: fmax       = 5.0_wp
      real(wp)         :: pi_alpha   = PI_ALPHA_DEF
      real(wp)         :: pi_beta    = PI_BETA_DEF
      !----- Embedded-pair LOWER order (MEDS_ED2_RK45_DESIGN.md sec 6): default 1 matches ARK's       !
      !      ARS(2,2,2) 1st-order embedded estimate; Cash-Karp's RK45 sets this to 4 so the           !
      !      I-controller uses the correct -1/5 exponent instead of silently reusing ARK's -1/2. ------!
      integer(ik)      :: p_order    = 1_ik
      type(tol_set_t)  :: tols
   end type error_control_t

   !----- Everything the fast-loop INTEGRATOR is configured by, in one record built once per run     !
   !      (meds_fast_control%build_integrator_opts) and carried on column_config_t, so the schemes    !
   !      read one named record instead of eight loose fields of the run configuration (2026-09       !
   !      review, decisions after items 4-6). ---------------------------------------------------------!
   type :: integrator_opts_t
      integer(ik) :: scheme           = INTEG_ARK  !< INTEG_ARK | INTEG_RK45
      logical     :: adaptive         = .true.     !< ARK: adaptive sub-stepping (else fixed_substeps)
      real(wp)    :: dt_init          = 0.0_wp     !< [s] first sub-step (<= 0: warm start / dt_fast)
      logical     :: coupled_newton   = .true.     !< ARK: coupled leaf<->CAS Newton (else uncoupled BE)
      integer(ik) :: fixed_substeps   = 1_ik       !< ARK, adaptive = .false.: equal sub-steps per dt_fast
      logical     :: cas_condensation = .true.     !< apply the CAS supersaturation sink (both schemes)
      type(error_control_t) :: error_control       !< tolerances + controller + strictness
   end type integrator_opts_t

   type :: column_config_t
      type(process_mask_t)        :: mask            !< process-complexity mask (all on = full column)
      type(aero_cfg_t)            :: aero            !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal    !< leaf/wood thermal params
      type(soil_params_t)         :: soil           !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal   !< soil thermal texture
      type(energy_opts_t)         :: energy         !< soil-thermal solver options
      type(soil_opts_t)           :: soil_water_opts !< soil-water (Richards) solver options
      type(wood_params_t)         :: wood           !< stem-respiration parameters
      type(root_params_t)         :: root           !< fine-root-respiration parameters
      type(hydro_params_t)        :: hydraulics_params  !< plant-hydraulics parameters (PV curves, vulnerability)
      type(hydro_opts_t)          :: hydraulics_opts    !< plant-hydraulics solver options
      type(leaf_photo_table_t)    :: leaf_photo     !< per-PFT leaf-photosynthesis parameters (built once per run)
      type(integrator_opts_t)     :: integrator     !< the fast-loop integrator's configuration (built once per run)
      real(wp)                    :: specific_root_area = 20.0_wp  !< [m2/kgC] SRA (rhizosphere conductance)
      !----- Canopy-surface water: interception film + film-evap/dew (MEDS_ED2_RK45_DESIGN.md sec 3.4, !
      !      P1) -- opt-in (default off, so existing configs are unchanged); SPLIT PATH ONLY for now,   !
      !      mirroring how snow (col_config%snow_on) and prognostic leaf/wood energy both landed split-first  !
      !      with ARK support deferred (column_fast_step error-stops if this is on under INTEG_ARK). ---!
      logical                     :: canopy_water_on  = .false.
      type(snow_params_t) :: snow                    !< snow parameters (density, albedo, thresholds, conductivity)
      !----- The picard_* mirrors that used to sit here were DELETED (plan E4): fill_ctx copied five   !
      !      config fields into them every slow step and nothing ever read them back. The comment      !
      !      claiming "ARK's newton_surface_solve uses the iteration cap" was false -- that cap is the  !
      !      NEWT_MAX parameter in meds_fast_ark. ------------------------------------------------------!
   end type column_config_t

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into cfg%pft)
      real(wp),    allocatable :: lai(:), wai(:), height(:), crown(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: aboveground_frac(:)          !< [--] gathered per-PFT (stem respiration)
      logical,     allocatable :: is_woody(:)                   !< gathered per-PFT (stem respiration off for grass)
      real(wp),    allocatable :: stem_resp_factor25(:)         !< [umol CO2/m2 stem/s @25C] gathered per-PFT
      real(wp),    allocatable :: root_resp_factor25(:)         !< [umol CO2/kgC root/s @25C] gathered per-PFT
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
      !----- Yesterday's daily-max leaf water potential, gathered from site%cohort (issue #95). It    !
      !      drives beta_stomata in leaf_gas_exchange, which was previously inert: env%psi (formerly psi_soil) is an     !
      !      OPTIONAL argument to leaf_gas_exchange_batch and this driver never passed it, so         !
      !      env%psi defaulted to 0 and beta_stomata was identically 1. -----------------------!
      real(wp),    allocatable :: dmax_psi_leaf(:)                        !< [MPa] <= 0; DMAX_PSI_LEAF_UNSET => seed from soil
   end type column_cohort_t

   !----- Prescribed per-step forcing the higher layers (RT, met) supply; photosynthesis/    !
   !      respiration/NEE are now computed from the plant kernels (no longer prescribed).     !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: rainfall        = 0.0_wp
      !< [kg/m2/s] met rainfall at the reference level (interception is applied downstream)
      real(wp)              :: snowfall         = 0.0_wp   !< [kg/m2/s] frozen rainfall (snowfall; drives snow accumulation)
      real(wp)              :: air_temp          = 288.0_wp !< [K] reference-level air temp (frozen/rain-on-snow rainfall enthalpy)
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
      !----- The step's NET CAS -> atmosphere export, b-weighted over the accepted march exactly as  !
      !      the conservation ledgers are: what atm_fluxes reports as LE and H. ---------------------!
      real(wp)       :: atm_heat_export = 0.0_wp   !< [J/m2]  sensible heat over this dt_fast
      real(wp)       :: atm_vap_export  = 0.0_wp   !< [kg/m2] over this dt_fast
      !----- Soil-carbon matrix Rh diagnostics (B2, MEDS_SLOW_DYNAMICS_DESIGN.md Part II): filled     !
      !      by column_prepass ONLY when cfg%soil_carbon_on (else left at 0, matching the OLD          !
      !      constant-pool scalar path that runs instead). xi_step is this sub-step's per-pool          !
      !      environmental decomposition scalar [-]; the caller (fast_dynamics) accumulates             !
      !      xi_step*dt_fast_days into the per-patch day-integral xi_int, and rh_matrix_step*dt_fast_    !
      !      days into rh_fast_accum (design section 9's audit-only rh_seam_gap check). -----------------!
      real(wp)       :: xi_step(n_soil_pool) = 0.0_wp   !< [-] this sub-step's per-pool env scalar
      real(wp)       :: rh_matrix_step       = 0.0_wp   !< [kgC/m2/day] this sub-step's matrix Rh
      !----- WORK counters (MEDS_NUMERICS_SCOPING.md section 5.3). These are the COST axis of the      !
      !      benchmark: without them a sweep can report accuracy but not accuracy-per-unit-work, and   !
      !      wall-clock alone is too coarse and too machine-dependent to rank schemes. Every one of     !
      !      these was already computed somewhere and then discarded. Set PER SUB-STEP by the stepper;  !
      !      the fast driver accumulates them site-wide. -------------------------------------------!
      integer(ik)    :: integ_nsteps   = 0_ik   !< accepted integrator sub-steps this dt_fast
      integer(ik)    :: integ_nrej   = 0_ik   !< rejected integrator steps this dt_fast
      integer(ik)    :: soil_nsub    = 0_ik   !< soil-water Richards solver sub-steps
      integer(ik)    :: hydro_nsub   = 0_ik   !< plant-hydraulics sub-steps, summed over cohorts
      !----- Set when this step's hydraulics sub-stepping crossed HYDRO_NSUB_THRASH per cohort       !
      !      (issue #104). Per step, like every other counter here; the driver area-weights it. -------!
      integer(ik)    :: hydro_thrash       = 0_ik    !< 1 = pathological hydraulics sub-stepping this step
      integer(ik)    :: hydro_nonconv = 0_ik  !< cohorts whose hydraulics solve did not converge
      !----- P6 (MEDS_ED2_RK45_DESIGN.md): count of sub-steps where the explicit RK45 step committed a   !
      !      railed (clamp-pinned, unphysical) CAS/soil state and the dispatcher rolled back + redid the   !
      !      step on the implicit-CAS ARK path. Rare (a handful over a healthy 30-yr run); a                !
      !      persistently-high value flags a genuinely stiff regime RK45 is degrading to ARK for. ----------!
      integer(ik)    :: rk45_rescue  = 0_ik   !< dt_fast steps rescued RK45->ARK this sub-step (0 on ARK)
      !----- CLAMP activations (MEDS_INTEGRATOR_PARITY.md [RETIRED], Phase A). The stability clamps            !
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

   !=====================================================================================!
   ! THE FROZEN RECORD, decomposed by physical content (2026-09 review, item 4 #1).           !
   !                                                                                          !
   ! Everything below is held constant across one dt_fast: the pre-pass (build_column_frozen)  !
   ! evaluates it once at state^n and every stage of either integrator reads it. "Frozen" is a  !
   ! statement about LIFETIME, so it is carried by the container column_frozen_t and by         !
   ! intent(in) at the call sites -- the pieces are named for what they describe, and a piece   !
   ! that a later scheme refreshes per stage (cas_boundary_t under a live surface-layer solve)  !
   ! keeps its name. A kernel takes only the pieces it reads: surface_derivs takes cas, tissue, !
   ! film, ground and snow; the soil-water and soil-heat tendencies take hydrology, roots and    !
   ! params; the hydraulics corrector takes plant, roots and params.                            !
   !=====================================================================================!

   !----- CAS <-> atmosphere boundary: capacities, bulk conductances, the surface-layer inputs a  !
   !      live re-solve needs, and the reference-level state. Plain scalars, so a stage can copy  !
   !      it to override the conductances at its own canopy-air state without touching an array. !
   type :: cas_boundary_t
      real(wp) :: cas_mass_capacity          = 0.0_wp      !< [kg/m2]   CAS mass capacity  -> enthalpy & vapour
      real(wp) :: cas_molar_capacity          = 0.0_wp      !< [mol/m2]  CAS molar capacity -> CO2
      real(wp) :: g_atm_heat           = 0.0_wp      !< [kg/m2/s] CAS<->atm enthalpy conductance
      real(wp) :: g_atm_vapour           = 0.0_wp      !< [kg/m2/s] CAS<->atm vapour   conductance
      real(wp) :: g_atm_co2           = 0.0_wp      !< [mol/m2/s]CAS<->atm CO2      conductance
      !=====================================================================================!
      ! THE CAS<->ATM CONDUCTANCES ARE RE-SOLVED AT EVERY STAGE, not frozen per dt_fast          !
      ! (MEDS_PRODUCTION_INTEGRATOR_PLAN.md sec 1g/5-N2).  g_atm_heat/g_atm_vapour/g_atm_co2 above are still WRITTEN    !
      ! by the state^n pre-pass -- they seed the march and every diagnostic that wants a          !
      ! representative value -- but a stage evaluates them at its OWN canopy-air state.           !
      !                                                                                          !
      ! WHY THIS ONE COEFFICIENT AND NOT THE REST.  Measured: holding g_atm_heat at its unperturbed      !
      ! value takes the freeze-cadence map's multiplier from Phi' = -23.2 to +0.80 at            !
      ! dt_fast = 900 s, while g_transp_leaf, h_coeff_leaf, abs_lw and f_wet_c each contribute <= 0.7%.     !
      ! The mechanism is the Monin-Obukhov feedback -- a warmer canopy air is a more unstable     !
      ! surface layer, which vents it harder (d ln g_atm_heat/dT ~ 2.2 /K) -- so lagging it by a whole   !
      ! dt_fast is what made the canopy air oscillate.  It is also one of the CHEAPEST things in  !
      ! the pre-pass (canopy_aerodynamics is 2% of it; leaf gas exchange, which stays frozen, is  !
      ! 89% at 30 cohorts), so this costs ~3% of a sub-step and pays for itself in fewer          !
      ! sub-steps.  ED2 does the same thing: canopy_turbulence8 at EVERY RK stage.                !
      !                                                                                          !
      ! There is deliberately NO SWITCH.  The frozen alternative is unstable at the production    !
      ! cadence on four of five stand heights measured (Phi' down to -8.2 at dt_fast = 150 s), so !
      ! it is known-wrong physics, not a supported configuration -- the same reasoning that       !
      ! deleted snow_on and the with_theta norm switch.                                           !
      !                                                                                          !
      ! mo_live is NOT that switch.  It records whether the mo_* INPUTS below are populated and     !
      ! the g_atm_heat/g_atm_vapour/g_atm_co2 above are therefore stale for any state other than the one they were       !
      ! solved at.  Three consistent uses, none of them a user choice:                             !
      !   * default .false. -- "use g_atm_heat/g_atm_vapour/g_atm_co2 exactly as given".  A hand-built record (a unit    !
      !     test, the RK4 oracle) supplies its own conductances and never populates the mo_*       !
      !     inputs, so this MUST be the default: re-solving from zeroed roughness and wind does    !
      !     not converge.  (It hung test_column_derivs when the default was the other way round.)  !
      !   * build_column_frozen sets .true. -- the inputs are populated, so re-solve per stage.    !
      !   * ARK's stage-local copy sets it back to .false. after refreshing ONCE, because          !
      !     surface_derivs is called up to 24 times per stage by the Newton to fill a CAS tendency !
      !     ARK never reads.  RK45 leaves it .true.: its every RHS evaluation IS a stage.          !
      !=====================================================================================!
      logical :: mo_live = .false.            !< the mo_* inputs are populated => re-solve per evaluation
      type(aero_cfg_t) :: aero_cfg            !< MO constants (copy; POD)
      real(wp) :: mo_u_ref      = 0.0_wp      !< [m/s]  wind at the reference height
      real(wp) :: mo_zref       = 0.0_wp      !< [m]    reference height
      real(wp) :: mo_displace   = 0.0_wp      !< [m]    displacement height  (canopy geometry, frozen)
      real(wp) :: mo_rough      = 0.0_wp      !< [m]    roughness length     (canopy geometry, frozen)
      real(wp) :: mo_theta_atm  = 0.0_wp      !< [K]    potential temperature at zref
      real(wp) :: mo_shv_atm    = 0.0_wp      !< [kg/kg] specific humidity at zref (the AERO reference,
                                              !<        which need not equal shv_atm below)
      real(wp) :: mo_rho        = 0.0_wp      !< [kg/m3] air density
      real(wp) :: enthalpy_atm      = 0.0_wp      !< [J/kg]    reference-level specific enthalpy
      real(wp) :: shv_atm       = 0.0_wp      !< [kg/kg]   reference-level specific humidity
      real(wp) :: co2_atm       = 400.0_wp    !< [umol/mol]free-atmosphere CO2
      real(wp) :: nee_biotic    = 0.0_wp      !< [umol/m2/s] frozen biotic CO2 source (Ra+Rh-GPP)
      real(wp) :: rho           = 0.0_wp      !< [kg/m3]   canopy-air density
      real(wp) :: press         = 0.0_wp      !< [Pa]      canopy-air pressure
      !----- SCHEME-ASYMMETRY GUARD (§8g). surface_derivs applies a smooth CAS supersaturation
      !      (condensation) sink. This switch makes the term controllable so a like-for-like
      !      comparison between schemes is possible; .true. (default) preserves the historic ARK
      !      behaviour exactly. Whether the sink belongs on BOTH paths is a model question,
      !      deliberately left open here.
      logical  :: cas_condensation = .true.  !< apply the CAS supersaturation sink (both schemes)
   end type cas_boundary_t

   !----- Per-cohort LEAF and WOOD energy-balance coefficients and the tissue heat store, frozen at !
   !      state^n. a_* = cap/dt_fast is the storage conductance veg_energy_balance relaxes against; !
   !      t_*0 is the start-of-step temperature it relaxes FROM. Both frozen for the whole fast    !
   !      step, so each stage evaluation returns the SAME dt_fast-averaged flux and dt_fast-endpoint !
   !      temperature. The store is deliberately NOT a tableau degree of freedom: it is an          !
   !      algebraic closure evaluated at each stage, which is why it needs no new WRMS group, no    !
   !      arrowhead and no Newton -- see MEDS_VEG_ENERGY_INTEGRATION_PLAN.md sec 2. The advective  !
   !      enthalpy terms (qwflux_wl, q_wood_net; ED2's qwflux_wl/qloss) are the water crossing the  !
   !      wood<->leaf and soil<->wood interfaces carrying its own thermal energy, frozen at state^n. !
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
      real(wp), allocatable :: qwflux_wl(:)   !< [W/m2 ground] sapflow's advected enthalpy INTO the leaf (wood->leaf)
      real(wp), allocatable :: q_wood_net(:)  !< [W/m2 ground] net advected enthalpy INTO wood (qloss - qwflux_wl)
      real(wp) :: leaf_emiss    = 0.95_wp     !< [-]       leaf LW emissivity
      !----- heat-capacity inputs behind leaf_hcap_per_dt/wood_hcap_per_dt: a_* = (dry hcap + water mass*cp_liq)/dt. !
      real(wp), allocatable :: wood_dry_hcap(:)   !< [J/m2/K]  dry sapwood heat capacity (floored)
      real(wp), allocatable :: wood_wmass(:)      !< [kg/m2]   fresh-sapwood water mass
      real(wp), allocatable :: leaf_dry_hcap(:)   !< [J/m2/K]  dry leaf heat capacity (floored)
      real(wp), allocatable :: leaf_wmass(:)      !< [kg/m2]   internal (symplast) leaf water mass
   end type tissue_coefficients_t

   !----- Canopy interception FILM: the wetted fraction and film-evaporation conductances frozen   !
   !      once per dt_fast (the pre-pass's intercept_canopy_layer sweep), the liquid enthalpy the   !
   !      film is valued at, and the frozen interception rates the film state integrates. All zero  !
   !      when canopy_water_on is off, so the film is a no-op unless build_column_frozen populates  !
   !      it. Only the STATE-dependent terms (dqdt, qsat_c - qcas) are re-evaluated per stage.      !
   type :: canopy_film_capacity_t
      real(wp), allocatable :: g_film_leaf(:), g_film_w(:)   !< [m/s] frozen film-evap conductance, leaf/wood
      real(wp), allocatable :: f_wet_c(:)                 !< [-]   frozen combined wetted fraction (sigma_w)
      !----- Liquid enthalpy the film is valued at (= internal_energy_liquid(t_film_valuation), the        !
      !      temperature intercepted water arrives with; 0 under a pack). The tissue pays            !
      !      enthalpy_vapor(T) - film_liquid_enthalpy per kg of film it evaporates, so film store + tissue +   !
      !      CAS close exactly (see surface_derivs). ------------------------------------------------!
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

   !----- Surface-block tendencies + the diagnostics the ARK ledger and the soil/hydraulics         !
   !      tendencies consume (coh_transp -> soil-water sink; transp_c   !
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
      real(wp) :: coh_transp     = 0.0_wp     !< [kg/m2/s]  total realized transpiration
      real(wp) :: cond           = 0.0_wp     !< [kg/m2/s]  smooth condensation sink (dew) draining CAS supersat
      real(wp) :: cond_enth      = 0.0_wp
      !< [W/m2]     the liquid enthalpy that sink debited from the CAS (one number, both sides)
      real(wp), allocatable :: leaf_temp(:)   !< [K]        diagnosed per-cohort leaf temperature
      real(wp), allocatable :: wood_temp(:)   !< [K]        diagnosed per-cohort wood temperature
      real(wp), allocatable :: transp_c(:)    !< [kg/m2/s]  per-cohort transpiration demand
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
      !      matching biophys%leaf_surf_water/wood_surf_water's own convention from the split-path P1 landing). !
      !----- SURFACE (ponding) store. In this phase it is PASSED THROUGH the stages and committed  !
      !      from the scratch hydrology solve, exactly as theta is (see column_fast_step_ark), so    !
      !      carrying it here changes nothing yet. It lives on the state vector so it can become     !
      !      PROGNOSTIC without a second lockstep sweep over every combinator -- issue #93 Phase 1.  !
      !      Until it has a stage RHS it must stay OUT of state_wrms_grouped, like the canopy-surface !
      !      films below, or the error controller would score a frozen quantity. -------------------!
      real(wp) :: w_surface      = 0.0_wp         !< [kg/m2] ponded surface water
      real(wp) :: w_surface_enth = 0.0_wp         !< [J/m2]  enthalpy of the ponded water
      real(wp), allocatable :: leaf_surf_water(:) !< [kg/m2 ground] canopy interception film on leaf
      real(wp), allocatable :: wood_surf_water(:) !< [kg/m2 ground] canopy interception film on wood
   end type column_state_t

   !----- Frozen inputs for the whole column: the surface pre-pass + the soil/hydraulics params +   !
   !      the frozen hydrology surface BCs + per-cohort geometry the hydraulics kernel needs.        !
   !----- The operator-split SOIL-WATER solve's outcome for this dt_fast: the scratch                !
   !      advance_soil_water_column's boundary fluxes, its end-of-step stores, the temperatures that      !
   !      value the water crossing each boundary, and its interior faces + post-solve corrections.    !
   !                                                                                                  !
   !      THE FACES AND CORRECTIONS ARE ARK-ONLY (issue #78 item 3). They are the right numbers for a  !
   !      scheme that commits the scratch solve's theta VERBATIM, which the ARK does. RK45 integrates  !
   !      its OWN theta, on which the scratch's faces move a different amount of water and the         !
   !      scratch's clip mass never moves at all, so column_derivs takes both from the stage's own     !
   !      soil_water_time_deriv. Using these there cost ~2.6e6 J/m2/step of vertical enthalpy          !
   !      misplacement (soil surface 345 K) against ~2.6e6 J/m2/step of spurious clip cooling -- two   !
   !      defects of matched magnitude and opposite sign, which is why each hid the other. The whole-  !
   !      column ledger cannot see a purely VERTICAL misplacement; only the faces' provenance protects  !
   !      against it. w_flux_frozen is the solver's TIME-MEAN face flux over the step, so it is exactly !
   !      what the split path advected on (eforc%w_flux = -hflux%w_flux).                              !
   type :: soil_hydrology_t
      real(wp) :: geothermal    = 0.0_wp          !< [W/m2]    bottom heat flux BC
      real(wp) :: q_top         = 0.0_wp          !< [m/s]     Richards top water flux (infiltration - evaporation)
      real(wp) :: infiltration  = 0.0_wp          !< [kg/m2/s] throughfall reaching the soil top face
      real(wp) :: drainage      = 0.0_wp          !< [kg/m2/s] bottom-face drainage
      !----- The scratch solve's two post-solve MASS corrections, summed over layers (their per-layer  !
      !      enthalpies are clip_enth/floor_enth below). The ARK commits the scratch theta verbatim, so !
      !      these are water that really left (clip -> pond) or was created (theta_res floor) in the    !
      !      committed state, and the ledgers must book them (2026-09 review, item 1A #5). ------------!
      real(wp) :: clip_mass     = 0.0_wp          !< [kg/m2/s] saturation-clip water leaving the soil for the pond
      real(wp) :: floor_mass    = 0.0_wp          !< [kg/m2/s] theta_res-floor water created in the soil
      real(wp) :: runoff_surf   = 0.0_wp          !< [kg/m2/s] surface runoff
      real(wp) :: precip_ground = 0.0_wp          !< [kg/m2/s] water reaching the ground (RK45 rebuilds its OWN pond from it)
      real(wp) :: t_infil       = 0.0_wp          !< [K]       temperature of the infiltrating (pond) water, #78 item 4
      real(wp) :: w_surface_enth1 = 0.0_wp        !< [J/m2]    scratch solve's end-of-step pond ENTHALPY, #78 item 4
      real(wp) :: t_pond_inflow      = 0.0_wp          !< [K]       temperature of the water entering the pond, #78 item 4
      real(wp) :: runoff_enth   = 0.0_wp          !< [W/m2]    enthalpy leaving with surface runoff, #78 item 4
      real(wp) :: t_film_valuation     = 0.0_wp          !< [K]       valuation T of intercepted water (tsupercool_liq under a pack)
      real(wp) :: t_bot         = 0.0_wp          !< [K]       bottom-layer soil temperature @ state^n
      real(wp) :: w_surface1    = 0.0_wp          !< [kg/m2]   end-of-step ponded surface water
      real(wp) :: w_flux_frozen(n_soil_layer_max) = 0.0_wp  !< [m/s]   DOWNWARD interior face flux, k=1..nsl-1
      !----- Enthalpy paired with the hydrology's UNFACED post-solve mass corrections, already valued  !
      !      at each layer's own state^n temperature (so the correction is temperature-NEUTRAL) and     !
      !      carried as [W/m2] to join the root_heat_sink column the stages already assemble. ---------!
      real(wp) :: clip_enth(n_soil_layer_max)  = 0.0_wp     !< [W/m2] enthalpy leaving layer k with clipped water
      real(wp) :: floor_enth(n_soil_layer_max) = 0.0_wp     !< [W/m2] enthalpy created with theta_res-floored water
      !----- the AUTHORITATIVE end-of-step soil moisture from the scratch advance_soil_water_column (the robust  !
      !      ponding/runoff/free-drain Richards solve). The ARK COMMITS this instead of re-solving theta in   !
      !      the ESDIRK stages (soil water is fully operator-split out; see column_fast_step_ark).            !
      real(wp), allocatable :: theta1(:)          !< [m3/m3]   committed post-step soil moisture (per layer)
   end type soil_hydrology_t

   !----- ROOT ZONE: the realized aggregate uptake, where it is placed, and the soil-side hydraulic  !
   !      boundary the corrector re-solves against. uptake is the pre-pass's plant-side REQUEST      !
   !      rescaled by the soil's OWN fwilt-limited supply -- the SAME number both the soil-water     !
   !      tendency's root sink and the per-cohort uptake_frozen are built from, so the wood<->soil   !
   !      interface closes to the soil's TRUE realized supply. root_share is THIS dt_fast's realized  !
   !      per-layer uptake shares (sum = 1), built from solve_plant_water_batch's own breakdown so   !
   !      the root MASS sink and the root HEAT sink land in the same layers by construction; it      !
   !      falls back to the static root_frac profile when no layer supplies anything. qloss_frozen  !
   !      is ED2's qloss: the liquid enthalpy the uptake carries out of the soil, per cohort.         !
   type :: root_zone_t
      real(wp) :: uptake        = 0.0_wp          !< [kg/m2/s] realized (post-rescale) aggregate root uptake
      real(wp), allocatable :: root_share(:)      !< [-]       per-layer root-sink shares (sum = 1)
      real(wp), allocatable :: psi_soil_pre(:)    !< [MPa] per-layer soil water potential @ state^n (nsl)
      real(wp), allocatable :: rhizo_cond(:,:)    !< [kg/plant/s/MPa] rhizosphere conductance (nsl, ncoh)
      real(wp), allocatable :: qloss_frozen(:)    !< [W/m2 ground] (ncoh) advected enthalpy leaving the soil with uptake
   end type root_zone_t

   !----- PLANT WATER: the pre-pass's time-averaged solve_plant_water output (per plant), held      !
   !      CONSTANT across every sub-stage of the macro-step -- sapflow_frozen is the wood->leaf       !
   !      transfer; uptake_frozen is the soil->wood transfer, floored >= 0 (no hydraulic              !
   !      redistribution) and rescaled so sum(uptake_frozen*nplant) == roots%uptake EXACTLY -- and    !
   !      the cohort geometry the hydraulics kernel reads when the post-stage corrector re-solves it   !
   !      with the REALISED b-weighted transpiration. column_derivs' mass ODE reads the two fluxes     !
   !      directly; no PV-curve/conductance evaluation is needed per stage.                           !
   type :: plant_water_t
      real(wp), allocatable :: nplant(:), bleaf(:), bsap(:), broot(:), sap_area(:), height(:), leaf_area(:)
      real(wp), allocatable :: sapflow_frozen(:), uptake_frozen(:)   !< [kg/plant/s] (ncoh)
   end type plant_water_t

   !----- Parameter records the stages read, COPIED from column_config_t once per dt_fast. They are  !
   !      here only because the march signatures (both schemes and the RK4 oracle) carry the frozen  !
   !      record and not the column configuration; passing them instead of copying them is the       !
   !      remaining step of the decomposition (2026-09 review, decisions after items 4-6). ----------!
   type :: column_params_t
      type(soil_params_t)         :: soil         !< soil geometry + texture (dz, root_frac, ...)
      type(soil_thermal_params_t) :: therm        !< soil thermal texture
      type(energy_opts_t)         :: energy_opts  !< soil-thermal options (phase change)
      type(soil_opts_t)           :: hydro_opts   !< soil-water (Richards) options
      type(hydro_params_t)        :: hydraulics_params      !< PV curves + vulnerability (for the corrector)
      type(hydro_opts_t)          :: hydraulics_opts      !< hydraulics kernel solver options (for the corrector)
   end type column_params_t

   !----- THE CONTAINER: everything held constant over one dt_fast, by physical content. -----------!
   !----- The frozen outcome of one pre-column snow advance. Every field is 0/.false. when snow is  !
   !      off or no pack exists, and the consumers are written so that those values reduce their     !
   !      arithmetic EXACTLY to the pre-C4 snow-free form -- which is what makes "snow-off            !
   !      bit-identical" a structural property rather than something to re-verify per scheme. -------!
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
      !----- enthalpy the melt transfer moved pack -> soil layer 1. Needed by any caller whose soil    !
      !      baseline is snapshotted AFTER this stage runs: that snapshot already contains the melt    !
      !      energy while enth0 still contains it too, so the pair double-counts it by exactly this    !
      !      amount. Split snapshots BEFORE the stage and needs no correction. ---------------------!
      real(wp) :: melt_enth  = 0.0_wp    !< [J/m2] melt enthalpy leaving the pack with the meltwater
      !----- Temperature that VALUES the meltwater, i.e. the T with u_liq(T)*melt_mass == melt_enth.     !
      !      The caller hands this to the hydrology kernel as chydro_forcing_t%t_pond_inflow so the pond      !
      !      receives exactly melt_enth when it receives melt_rate*dt of mass -- one number, both        !
      !      sides. Falls back to t_3ple when there is no melt mass to value. ------------------------!
      real(wp) :: t_melt     = 0.0_wp    !< [K] effective temperature of the meltwater
   end type snow_stage_t

   type :: column_frozen_t
      type(cas_boundary_t)         :: cas          !< CAS <-> atmosphere boundary
      type(tissue_coefficients_t)  :: tissue       !< per-cohort leaf/wood energy coefficients + heat store
      type(canopy_film_capacity_t) :: film         !< canopy interception film
      type(ground_boundary_t)      :: ground       !< bare-ground boundary
      type(snow_stage_t)           :: snow         !< the snow stage's outcome (all zero without a pack)
      type(soil_hydrology_t)       :: hydrology    !< the scratch soil-water solve's outcome
      type(root_zone_t)            :: roots        !< realized uptake, its placement, the rhizosphere boundary
      type(plant_water_t)          :: plant        !< frozen sapflow/uptake + cohort geometry for the corrector
      type(column_params_t)        :: params       !< parameter copies (see column_params_t)
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
   !      water-enthalpy advection); the deferred rainfall>0 guard-lift adds those terms.                 !
   type :: stage_bflux_t                                    !< per-stage RATES
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp    !< [W/m2]
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp    !< [kg/m2/s]
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp    !< [umol/m2/s]
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp  !< [W/m2]
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp  !< [kg/m2/s]
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp!< [W/m2]
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp!< [kg/m2/s]
      real(wp) :: whole_cond    = 0.0_wp                         !< [kg/m2/s] condensate (row 1b)
      real(wp) :: whole_cond_enth = 0.0_wp                       !< [W/m2] its liquid enthalpy at the stage CAS temperature
      !----- NET CAS -> atmosphere export, the two turbulent fluxes the run REPORTS as H and LE       !
      !      (meds_fast_step%atm_fluxes): the SAME stage states and conductances as the ledger terms  !
      !      above, so the reported vapour export IS the conserved one (2026-09 review, item 4 #11).  !
      !      The sensible flux is spelled out (cp_air, dry-air basis) rather than taken as "enthalpy   !
      !      export minus L_v*vapour": the CAS enthalpy values vapour at cp_vap*(T - tsupercool_vap)  !
      !      ~ 3.4 MJ/kg, so that difference would carry the water's LIQUID-datum enthalpy (~40% of   !
      !      LE) inside H. ---------------------------------------------------------------------------!
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
      !----- CONDENSATE, tracked SEPARATELY from whole_wat_out (row 1b). Dew/fog is not a boundary     !
      !      loss -- it lands on a surface inside the column -- but it used to be summed into           !
      !      whole_wat_out with the atmospheric vapour flux, where it could not be told apart or        !
      !      redirected. Carrying it in its own slot is what lets the caller deposit it into a store.   !
      real(wp) :: whole_cond    = 0.0_wp   !< [kg/m2] condensed vapour over the step (>= 0)
      !----- ...and the liquid enthalpy it left the CAS with, b-weighted at each stage's OWN CAS      !
      !      temperature -- the SAME number surface_derivs debited from the CAS. The deposit into    !
      !      soil layer 1 must carry this, not cond*u_liq(T_end): valuing the deposit at the         !
      !      end-of-step temperature while the debit ran per stage left sum b_i*cond_i*(u(T_end) -   !
      !      u(T_i)) unbooked on every dew step (2026-09 review, item 1A #6). ---------------------!
      real(wp) :: whole_cond_enth = 0.0_wp !< [J/m2]
      real(wp) :: atm_heat_out = 0.0_wp    !< [J/m2]  net CAS -> atmosphere sensible-heat export (reported H)
      real(wp) :: atm_vap_out  = 0.0_wp    !< [kg/m2] net CAS -> atmosphere vapour export (reported LE / L_v)
      !----- TISSUE-TEMPERATURE TIME INTEGRALS [K*s], per cohort, b-weighted across stages and summed !
      !      over accepted sub-steps. These are what make the tissue store conserve EXACTLY on an      !
      !      adaptive scheme, and they are also the physically right answer rather than merely the     !
      !      conserving one.                                                                          !
      !                                                                                          !
      !      The frozen-store kernel's own balance is cap*(T_end - T_0)/dt_fast = numer - denom*dT_avg,!
      !      so the store's gain RATE at a stage is store_hcap_per_dt*(surf_tend%leaf_temp - t_leaf0) -- a plain        !
      !      function of what surface_derivs already returns. The canopy air receives the b-weighted   !
      !      denom*dT_avg, so the store must gain the b-weighted complement, i.e. the store's energy   !
      !      is set by the TIME INTEGRAL of the tissue temperature, not by its final-stage value.      !
      !      Committing T = integral/dt_fast makes the state and the ledger the SAME number.           !
      !                                                                                          !
      !      With a CAS that is constant over the step every stage returns the same temperature and    !
      !      the integral collapses to it, so this degenerates correctly. -----------------------------!
      real(wp), allocatable :: tissue_leaf_int(:), tissue_wood_int(:)   !< [K*s]
   end type column_bflux_t


   !----- Per-patch fast biophysics STATE (prognostic; carried between fast steps). The self- -!
   !      contained MVP block used by meds_fast_ark/meds_fast_rk45; the eventual per-cohort/per-patch !
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
      !----- FROZEN daily leaf/root-turnover shed-water rate (P4, MEDS_ED2_RK45_DESIGN.md): a       !
      !      read-only snapshot of site%patch%shed_water_rate(ip), seeded ONCE at the top of the      !
      !      slow step and held constant across the day's fast sub-steps, exactly like soil_carbon     !
      !      just above -- the fast loop adds it to its ground-water input every sub-step (never        !
      !      mutated here; meds_vegetation_dynamics is the sole writer of the real site-level rate). ---!
      real(wp)                   :: shed_water_rate = 0.0_wp  !< [kg/m2 ground/s]
      !----- The slow tier's CO2 owing, same lifetime and same read-only discipline as the shed   !
      !      water above: seeded ONCE per patch at the top of the fast window, constant across    !
      !      every sub-step, never written here.  -------------------------------------------!
      real(wp)                   :: slow_co2_rate   = 0.0_wp  !< [umol/m2 ground/s]
      real(wp), allocatable      :: leaf_temp(:)      !< [K] per-cohort leaf temperature
      real(wp), allocatable      :: wood_temp(:)      !< [K] per-cohort wood/branch temperature (own store)
      !----- Internal (xylem/symplast) water mass [kg/plant] -- the prognostic hydraulic state;    !
      !      psi is diagnosed from it wherever needed (psi_from_water_content), never persisted.    !
      !      MEDS_ED2_RK45_DESIGN.md sec 4. -----------------------------------------------------!
      real(wp), allocatable      :: leaf_water_mass(:) !< [kg/plant] internal leaf water
      real(wp), allocatable      :: wood_water_mass(:) !< [kg/plant] internal wood water
      !----- Surface (interception film) water [kg/m2 ground] -- DISTINCT store from the internal      !
      !      water above; MEDS_ED2_RK45_DESIGN.md sec 3.4. Already ground-area-referenced (unlike the   !
      !      per-plant internal water), so fusion SUMS it (meds_demography_cohort_fusefiss.f90). ------------!
      real(wp), allocatable      :: leaf_surf_water(:) !< [kg/m2 ground] leaf interception film
      real(wp), allocatable      :: wood_surf_water(:) !< [kg/m2 ground] wood interception film
      !----- Lagged per-layer root-uptake SHARES (sum = 1), from the previous fast step's multi-layer  !
      !       plant solve; the soil sink distributes coh_transp by these (vs static root_frac) so the   !
      !       soil dries where roots actually took water. Default 0 => root_frac fallback (first step /  !
      !       single-layer). RETIRED with the multilayer_roots flag (Phase 1): the per-layer shares are
      !       now built from THIS step's realized uptake, inline on the split path and on column_frozen_t
      !       for ARK/RK45, so nothing lags through patch state any more. ---------------------------!
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

   !----- Is the column the FULL system? Only then are the closed-budget halts meaningful (a frozen  !
   !      store still exchanges with its neighbours, so a reduced column cannot conserve). -----------!
   pure logical function mask_is_full(m)
      type(process_mask_t), intent(in) :: m
      mask_is_full = m%veg_energy .and. m%cas_energy .and. m%cas_vapour .and. m%cas_co2 .and.        &
                     m%soil_heat  .and. m%soil_water .and. m%hydraulics
   end function mask_is_full

   !----- Allocate a column_cohort_t (the per-patch cohort SoA the fast loop consumes). ------!
   subroutine alloc_column_cohort(col_cohort, n)
      type(column_cohort_t), intent(out) :: col_cohort
      integer(ik),           intent(in)  :: n
      col_cohort%n = n
      allocate(col_cohort%pft(n), col_cohort%lai(n), col_cohort%wai(n), col_cohort%height(n), col_cohort%crown(n),                &
               col_cohort%leaf_width(n), col_cohort%branch_diam(n), col_cohort%aboveground_frac(n),                   &
               col_cohort%is_woody(n), col_cohort%stem_resp_factor25(n), col_cohort%root_resp_factor25(n),            &
               col_cohort%leaf_area(n), col_cohort%nplant(n),                                                          &
               col_cohort%dbh(n), col_cohort%broot(n), col_cohort%bleaf(n), col_cohort%bsap(n), col_cohort%sap_area(n),           &
               col_cohort%bwood(n), col_cohort%vcmax25(n), col_cohort%rd25(n), col_cohort%dmax_psi_leaf(n))
      col_cohort%pft = 1_ik
      col_cohort%lai = 0.0_wp ; col_cohort%wai = 0.0_wp ; col_cohort%height = 0.0_wp ; col_cohort%crown = 1.0_wp
      col_cohort%leaf_width = 0.04_wp ; col_cohort%branch_diam = 0.02_wp ; col_cohort%aboveground_frac = 0.7_wp
      col_cohort%is_woody = .true.
      col_cohort%stem_resp_factor25 = 0.0_wp ; col_cohort%root_resp_factor25 = 0.0_wp
      col_cohort%leaf_area = 0.0_wp ; col_cohort%nplant = 0.0_wp ; col_cohort%dbh = 0.0_wp ; col_cohort%broot = 0.0_wp
      col_cohort%bleaf = 0.0_wp ; col_cohort%bsap = 0.0_wp ; col_cohort%sap_area = 0.0_wp
      col_cohort%bwood = 0.0_wp                     ! was ALLOCATED and never initialized
      col_cohort%vcmax25 = 0.0_wp ; col_cohort%rd25 = 0.0_wp
      !----- UNSET, not 0. A 0 here would read as FULLY TURGID and silently disable the stomatal    !
      !      stress limb for any caller that forgets to fill it; the sentinel makes column_prepass    !
      !      seed from the soil instead, which is the safe default. --------------------------------!
      col_cohort%dmax_psi_leaf = DMAX_PSI_LEAF_UNSET
   end subroutine alloc_column_cohort

   !---------------------------------------------------------------------------------------!
   ! Grow-only capacity check: reallocate ONLY when the current backing arrays are too small !
   ! for n cohorts; otherwise just update the active count col_cohort%n and leave the (over-sized)    !
   ! capacity in place (MEDS_NUMERICS_SCOPING.md BB1 phase 1). Every downstream reader loops    !
   ! over col_cohort%n (never size(col_cohort%pft) -- verified: no caller does), and the per-patch gather      !
   ! that follows a call to this routine overwrites indices 1..n unconditionally, so reusing a   !
   ! larger patch's leftover capacity for a smaller patch is bit-identical: the caller sizes      !
   ! col_cohort ONCE per fast_dynamics call (to the site-wide max cohort count) instead of once PER      !
   ! PATCH, cutting O(n_patch) heap allocations to O(1) per slow step.                            !
   !---------------------------------------------------------------------------------------!
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

   !----- Flatten the shared [hydraulics] config into the plant hydro_params_t + rhizosphere      !
   !       conductance, and build the vulnerability lookup table from wood_kexp. The single seam    !
   !       between cfg%hydraulics (shared, TOML-driven) and the fast loop's hydro_params_t (plant),  !
   !       mirroring how the leaf seam flattens the PFT photosynthesis traits. -------------------!
   subroutine apply_hydraulics_config(hcfg, hydraulics_params)
      type(hydraulics_config_t), intent(in)    :: hcfg
      type(hydro_params_t),      intent(inout) :: hydraulics_params
      hydraulics_params%leaf_pi0       = hcfg%leaf_pi0       ; hydraulics_params%leaf_elastic_mod       = hcfg%leaf_elastic_mod
      hydraulics_params%leaf_apoplast_frac = hcfg%leaf_apoplast_frac
      hydraulics_params%leaf_water_sat     = hcfg%leaf_water_sat
      hydraulics_params%wood_pi0       = hcfg%wood_pi0       ; hydraulics_params%wood_elastic_mod       = hcfg%wood_elastic_mod
      hydraulics_params%wood_apoplast_frac = hcfg%wood_apoplast_frac
      hydraulics_params%wood_water_sat     = hcfg%wood_water_sat
      hydraulics_params%wood_psi50     = hcfg%wood_psi50     ; hydraulics_params%wood_kexp      = hcfg%wood_kexp
      hydraulics_params%k_plant_max    = hcfg%k_plant_max    ; hydraulics_params%wood_kmax      = hcfg%wood_kmax
      hydraulics_params%vessel_curl    = hcfg%vessel_curl
      call build_hydro_table(hydraulics_params%vuln_table, hydraulics_params%wood_kexp)
   end subroutine apply_hydraulics_config

   !----- Allocate + seed a patch_biophys_t from an initial CAS temperature (mirrors the other !
   !      alloc_* helpers; seeds can_enthalpy via the shared thermo inverter). ----------------!
   subroutine alloc_patch_biophys(biophys, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(patch_biophys_t), intent(out) :: biophys
      integer(ik),           intent(in)  :: n_coh
      real(wp),              intent(in)  :: can_temp0, can_shv0, can_co2, leaf_temp0
      allocate(biophys%leaf_temp(n_coh), biophys%wood_temp(n_coh))
      allocate(biophys%leaf_water_mass(n_coh), biophys%wood_water_mass(n_coh))
      allocate(biophys%leaf_surf_water(n_coh), biophys%wood_surf_water(n_coh))
      biophys%leaf_temp        = leaf_temp0
      biophys%wood_temp        = leaf_temp0
      biophys%leaf_water_mass  = 0.0_wp    ! scratch seed only -- always discarded by the next real gather
      biophys%wood_water_mass  = 0.0_wp    ! (mirrors leaf_temp/wood_temp's own scratch-seed discipline)
      biophys%leaf_surf_water  = 0.0_wp    ! ditto
      biophys%wood_surf_water  = 0.0_wp
      biophys%cas%can_temp     = can_temp0
      biophys%cas%can_shv      = can_shv0
      biophys%cas%can_co2      = can_co2
      biophys%cas%can_enthalpy = cas_enthalpy_of_temp(can_temp0, can_shv0)
   end subroutine alloc_patch_biophys

   !---------------------------------------------------------------------------------------!
   ! Grow-only capacity check for the per-cohort arrays of patch_biophys_t (mirrors            !
   ! ensure_column_cohort_capacity, MEDS_NUMERICS_SCOPING.md BB1 phase 1). Does NOT touch        !
   ! biophys%cas/soil_e/soil_w/snow/soil_carbon -- every caller overwrites those with the site's      !
   ! persisted per-patch reservoirs (site%patch%cas(ip) etc.) immediately after allocating, so    !
   ! their alloc_patch_biophys seed values are always discarded; only leaf_temp/wood_temp/         !
   ! leaf_water_mass/wood_water_mass/leaf_surf_water/wood_surf_water need their CAPACITY ensured     !
   ! here (the caller's gather loop fills indices 1..n_coh).                                          !
   !---------------------------------------------------------------------------------------!
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
