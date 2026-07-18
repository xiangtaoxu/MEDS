# MEDS adversarial code review — 2026-07-06

**Scope.** A full-tree adversarial review of `MEDS/src/**` (~10.4k LOC, 47 modules), run as a
multi-agent workflow: one *finder* per review unit (16 units covering biophysics, plant physiology,
demography, state, driver, shared/numerics, biogeochemistry, IO, allometry/init, and a whole-tree
organization sweep), followed by an *adversarial verifier* that re-read the actual code and tried to
**refute** each finding. Three defect classes were targeted: **physical-process bugs** (wrong
science/math/units), **performance** (numerical-solver underperformance, wasted work, GPU-offload
cost), and **code organization** (against the project's own rules — stateless process modules, the
`shared ← {allometry, plant} ← state ← demography ← aux` DAG, one-module-per-file, readable names,
≤132 columns, no hidden global mutable state).

**Verification status (read this first).** The workflow was interrupted by a session limit and
resumed once, then stopped by the user before the verify fan-out finished. Consequently:

- **7 units were fully adversarially verified** — `bio_rt`, `bio_hydro`, `bio_energy`, `bio_aero`,
  `plant_leaf`, `plant_hydraulics`, `dem_dynamics` (31 machine verdicts: **17 CONFIRMED, 8 PLAUSIBLE,
  6 REFUTED**). Refuted findings were dropped (listed in Appendix B).
- **9 units produced findings but were not machine-verified** — `plant_crp`, `dem_fusefiss`,
  `dem_state`, `allom_init`, `shared_num`, `driver`, `biogeochem`, `io`, `org_sweep`. For these, the
  **top-severity findings were spot-verified by hand against the source** for this report (marked
  *verified-by-hand*); the remainder are **REPORTED (unverified)** and should be confirmed before
  acting.

Each finding below carries a confidence tag: **CONFIRMED** (verified, real), **PLAUSIBLE** (real code
observation, impact depends on a stated assumption), **REPORTED** (from a finder, not yet
independently verified). Line numbers are 1-indexed anchors on the working tree.

---

## ✅ Resolution status (updated 2026-07-08) — FULLY ADDRESSED

Every finding in this review has been triaged and dispositioned. Each was independently
re-verified against the source (multi-agent adversarial triage + hand review) and then either
**FIXED** (behavior-preserving / bit-identical, validated on ifx + nvfortran multicore +
nvfortran GPU — 24/24 each; carbon-mode/GPU/IO paths validated end-to-end, incl. live netCDF on
both back ends), **DEFERRED** with a stated reason (would change numerics, or a large/architectural
rework), or found **NOT-REAL**. Full per-finding log in **Appendix C**.

| Section | Merged PR | Outcome |
|---|---|---|
| §2.1 Critical / high bugs | #29 | 8 addressed (7 fixed, 1 verified-correct) |
| §2.2 Medium bugs | #30 | 10 fixed, 7 not-real/deferred |
| §2.3 Low bugs | #31 | 16 fixed, 12 not-real/deferred (+3 regression tests) |
| §3 Performance / solver | #32 | 6 fixed, 6 deferred |
| §4 Code organization | #33 | 10 fixed, 6 deferred/not-real |

The default `empirical` + `fast_biophysics_on=false` run is **unchanged** by all of the above.

---

## 1. Executive summary

| Category | CONFIRMED | PLAUSIBLE | REPORTED (unverified) | Total distinct |
|---|---:|---:|---:|---:|
| Physical-process bugs | 12 | 7 | 18 | 37 |
| Performance / solver | 6 | 1 | 6 | 13 |
| Code organization | 8 | 0 | 11 | 19 |

**Top must-fix items (triage):**

1. **[CRITICAL · bug · verified-by-hand]** The fast biophysics loop is never invoked by the
   executable, so `growth_source="carbon"` + `fast_biophysics_on=true` grows silently on **zero
   GPP**. — `src/driver/meds_main.f90:124`, `src/driver/meds_stepper.f90:37`.
2. **[HIGH · bug · verified-by-hand]** Cohort fusion & fission overwrite the carbon-mode prognostic
   pools with empirical allometry, destroying prognostic carbon on every fuse/split. —
   `src/demography/meds_demography_fusefiss.f90:250,291` → `src/state/meds_demography_types.f90:429`.
3. **[HIGH · bug · verified-by-hand]** Autotrophic **maintenance respiration is omitted** from the
   carbon-growth balance (only growth respiration is deducted), so NPP is overestimated in carbon
   mode. — `src/driver/meds_vegetation_dynamics.f90:103`.
4. **[HIGH · bug · CONFIRMED]** C4 CO₂/PEP-limited rate `kp_eff` has **no temperature response** —
   every other C4 term is temp-scaled. — `src/plant/meds_leaf_gas_exchange.f90:227`.
5. **[HIGH · bug · verified-by-hand]** `toml_real_array` ignores `iostat` and **silently zero-fills
   unparseable array entries** while reporting full length — a malformed *required* PFT trait vector
   becomes zeros with no error. — `src/io/meds_toml.f90:191`.
6. **[MEDIUM · perf · CONFIRMED]** The daily `growth_step`/`mortality_step` map **~18 full cohort
   arrays host↔device every step** (8 write-only arrays as `tofrom`), migrating the whole cohort
   state each day. — `src/demography/meds_demography_dynamics.f90:75,78`.
7. **[MEDIUM · correctness-on-nvfortran · REPORTED]** Array-valued `cstr()` passed straight into
   netCDF calls hits the documented nvfortran miscompile trap (silent wrong values at `-O2`). —
   `src/io/meds_io.f90:64`.

**Cross-cutting theme.** The most consequential defects cluster at the **carbon-mode / fast-loop
seam** (items 1–3): the opt-in mechanistic-carbon pathway is wired at the module level but not
integrated end-to-end, and several routines silently assume the empirical default. The default
`empirical` + `fast_biophysics_on=false` run is unaffected by items 1–3. A second theme is
**silent failure on malformed/edge input** (config parsing, zero-guards, sentinel handling) — many
individually-low bugs that share the anti-pattern of swallowing an error instead of `error stop`.

---

## 2. Physical-process bugs

### 2.1 Critical / high

> **✅ Resolution:** all 8 findings addressed in **PR #29** (7 fixed, soil-advection verified correct — see Appendix B/C).

**[CRITICAL · verified-by-hand] Fast biophysics loop never runs in the executable → carbon+fast grows on zero GPP.**
`src/driver/meds_main.f90:124` calls `advance_one_step(site, cfg, is_new_month, is_new_year)` with no
`fast_ctx`. `meds_stepper.f90:37` runs the fast loop only `if (cfg%fast_biophysics_on .and.
present(fast_ctx))`, so `run_fast_biophysics` — the **sole** writer of `cohort%gpp_accum`
(`meds_fast_loop.f90:164`) — never executes from the real program. `gpp_accum` stays at its init value
0. Then `meds_vegetation_dynamics.f90:98–99` sets `a_carbon = cohort%gpp_accum(j)` when
`fast_biophysics_on`, so the flagship `growth_source="carbon"` + `fast_biophysics_on=true`
configuration produces `net_carbon ≈ 0` → **zero structural growth, with no error or warning**.
*Impact:* the headline mechanistic pathway is silently dead at the top level. *Fix:* construct a
`fast_context_t` in `meds_main` and pass it to `advance_one_step`; or hard-`error stop` when
`fast_biophysics_on` is set but no `fast_ctx` is wired, so the trap can't pass silently.

**[HIGH · verified-by-hand] Cohort fusion & fission destroy carbon-mode prognostic pools.**
`fuse_2_cohorts` (`meds_demography_fusefiss.f90:249–250`) re-derives `dbh` from the conserved
per-plant AGB via the **empirical** `agb_to_dbh`, then calls `set_cohort_size`, which overwrites all
four carbon pools from allometry (`meds_demography_types.f90:429–432`:
`leaf_carbon`, `wood_carbon`, `fineroot_carbon`, `nonstructural_carbon`). `split_cohorts` does the same
(`:291`). In carbon mode `wood_carbon` is the *prognostic size anchor* and the pools carry independent
state — the module even provides `set_cohort_size_from_carbon` (`:441`) for exactly this inverse — but
neither fuse nor split uses it, and the pools are never nplant-weight-merged the way `agb`/`gpp_accum`
are. *Impact:* every monthly fuse/split silently resets prognostic carbon to the empirical allometric
manifold, breaking carbon conservation and erasing pool history. *Fix:* in carbon mode, nplant-weight
the four pools across donor+survivor and call `set_cohort_size_from_carbon` (anchor on conserved
`wood_carbon`), mirroring the empirical branch.

**[HIGH · verified-by-hand] Maintenance respiration omitted from the carbon-growth balance.**
`meds_vegetation_dynamics.f90:103` computes `env%net_carbon = a_carbon -
growth_respiration(a_carbon, …)` — only **growth** respiration is deducted. Leaf/stem/fine-root
**maintenance** respiration (`meds_plant_respiration` exists but is never called on this seam) is
absent, so `net_carbon` = GPP − Rgrowth instead of GPP − Rgrowth − Rmaint. *Impact:* NPP (hence
growth) is systematically overestimated in carbon mode; the sign of the bias grows with standing
biomass. *Fix:* call the maintenance-respiration kernel per cohort and subtract it here (or subtract
it inside `gpp_accum` during the fast loop), consistent with the CAS-CO₂ respiration accounting.

**[HIGH · CONFIRMED] C4 CO₂-limited rate `kp_eff` has no temperature response.**
`src/plant/meds_leaf_gas_exchange.f90:227` builds `kp_eff = p%kp25 * pressure / p_std` (pressure
correction only); `assim_demand_c4` uses `ap = kp_eff * ci`. Lines 206–212 temperature-scale every
other C4 term (Vcmax, Vpmax, Rd). *Impact:* the PEP/CO₂-limited branch of C4 photosynthesis is frozen
at its 25 °C value → wrong low-CO₂/low-temperature C4 assimilation. *Fix:* apply the appropriate
Arrhenius/Q10 factor to `kp` (ED2/Collatz scale `kp` with temperature).

**[HIGH · verified-by-hand] `toml_real_array` silently zero-fills unparseable entries.**
`src/io/meds_toml.f90:190–197`: `tmp = 0.0_wp`; `read(s,*,iostat=ios) …` — **`ios` is never checked**;
`out(k)=tmp(k)` for all k and `nout = count_tokens(s)`. A malformed numeric token makes the
list-directed read stop, leaving trailing entries at 0.0, and `nout` still reports the full token
count. *Impact:* given the project's "every parameter required from TOML, no defaults" contract, a
typo'd required array (e.g. a per-PFT trait vector) becomes zeros that pass validation — a silent
scientific-wrongness vector. *Fix:* check `ios` and `error stop` with the key and offending text;
derive `nout` from the actual parsed count.

**[HIGH · REPORTED] Soil interior heat advection uses the downwind cell / wrong sign.**
`src/biophysics/meds_column_energy.f90:82` — the advective heat term for interior soil faces selects
the wrong upstream cell relative to the downward-positive `w_flux` convention. *Impact:* sensible heat
carried by infiltrating/draining water moves the wrong direction, biasing the soil thermal profile
(worse with strong fluxes). *Fix:* upwind on the sign of `w_flux` at each face; add a regression that
advects a temperature step through the column. *(Bio_energy was a verified unit, but this specific
finding was not among the machine verdicts — verify before acting.)*

**[HIGH · REPORTED] Dunne saturation-excess runoff spuriously large under the default free-drain BC.**
`src/biophysics/meds_column_hydrology.f90:122` — with the default free-drain bottom BC the diagnosed
water table `z_wt` is pinned at the shallow column bottom, so the `f_sat` saturated fraction (and thus
Dunne runoff) is systematically overstated. *Impact:* over-generates surface runoff, drying the
column and starving the root sink. *Fix:* decouple `f_sat` from a bottom-pinned `z_wt` under
free-drain, or diagnose `z_wt` from the actual saturated-from-below depth.

**[HIGH · CONFIRMED] atm↔CAS conductance omits the `temp1` profile factor.**
`src/biophysics/meds_column_energy.f90:267` sets `gatm = rho_air*ustar` with no profile-function
factor, and uses it for the implicit enthalpy exchange (`:269`). `cas_atm_forcing_t` carries only
`ustar` and has no `temp1` field. *Impact:* the canopy-air ↔ free-atmosphere exchange coefficient is
mis-scaled → biased CAS temperature/humidity and surface fluxes. *Fix:* carry and apply the
`temp1`-style profile factor (ED2 `rho*ustar*temp1`), as the CO₂ twin design already prescribes.

### 2.2 Medium

> **✅ Resolution:** 10 fixed, 7 not-real/deferred in **PR #30** — see Appendix C.

- **[CONFIRMED] Hard `theta_res` moisture floor injects water with no bookkeeping.**
  `meds_column_hydrology.f90:162` — `theta1(k)=max(theta1(k),theta_res(k))` can raise moisture with no
  matching debit; the give-back loop only compensates rooted layers. Breaks the machine-precision
  water budget for undershooting rootless layers. *Fix:* route the floor correction through the same
  `deficit` accounting used for the root sink.
- **[CONFIRMED] `energy_resid` is identically zero by construction.**
  `meds_column_energy.f90:91–110` — the conservative commit makes the residual telescope to exactly
  the boundary-flux definition, so the "conservation self-check" can never detect a solver error.
  *Fix:* form the residual from an independent energy inventory (pre/post internal-energy sum vs
  imposed boundary fluxes), not from the update identity.
- **[CONFIRMED] Transpiration uses stomatal conductance only, ignoring boundary-layer resistance.**
  `meds_leaf_gas_exchange.f90:350` — `flux%transpiration = gs*vpd/pressure`, dropping the `gb` series
  term that the CO₂ solve *does* include. Inconsistent water vs carbon conductance → transpiration
  overestimated when `gb` is finite. *Fix:* use the series `1/(1/gs+1/gb)` (or `gb`-aware total),
  consistent with the assimilation path; gated by `use_boundary_layer`.
- **[REPORTED] `include_pft` gates seed rain but not the reproduction recruit flux.**
  `meds_demography_rates.f90:98` — a PFT excluded from establishment can still spawn cohorts via the
  reproduction-carbon recruit path. *Fix:* apply the `include_pft` mask to both recruitment terms.
- **[REPORTED] `fast_soil_carbon` respired every step but never decremented.**
  `src/biogeochemistry/meds_column_co2.f90:203` — heterotrophic respiration draws from a pool that is
  never debited; `resid==0` masks the non-conservation. *Fix:* decrement the pool by the respired
  flux (this is the P0 prognostic-pool seam the biogeochem design formalizes).
- **[REPORTED] Unrecognized `hr_model` selector silently computes Q10.**
  `meds_column_co2.f90:118` — an unknown heterotrophic-respiration model falls through to Q10 instead
  of `error stop`. *Fix:* `error stop` on an unrecognized selector.
- **[REPORTED] `phenology_status` only gates the storage-flush half of leaf/fineroot growth.**
  `src/plant/meds_plant_carbon_dynamics.f90:107` — the NPP-funded half of leaf/root growth is not
  gated by phenology, so tissue can grow while phenology says "off". *Fix:* gate both allocation paths
  on `phenology_status`.
- **[REPORTED] `split_cohorts` renorm uses the uncapped AGB exponent.**
  `meds_demography_fusefiss.f90:271` — the ±eps daughter-diameter renormalization uses the uncapped
  `dbh→AGB` exponent, so conservation is violated for cohorts pinned at `hgt_max` (where the height
  allometry saturates). *Fix:* use the height-capped AGB relation in the renorm, or conserve on AGB
  directly.
- **[REPORTED] Actual temperatures fed to aerodynamics as potential temperatures.**
  `src/driver/meds_fast_loop.f90:220` — real (not potential) temperatures are passed into the M-O
  stability inputs, biasing the atm↔canopy stability parameter. *Fix:* convert to potential
  temperature before the stability calculation.
- **[REPORTED] PFT count inferred from `wood_density` has no `MAXPFT` bound.**
  `src/io/meds_config_io.f90:351` — an over-length trait array reads past the 64-element buffer. *Fix:*
  bound the inferred PFT count by `MAXPFT` and `error stop` on overflow.
- **[REPORTED] `init_from_census` unconditionally discards the first data row as a header.**
  `src/init/meds_init.f90:133` — a header-less inventory CSV silently loses its first cohort. *Fix:*
  detect the header (non-numeric first field) or make it an explicit config flag.
- **[REPORTED] `adaptive_step_update` divides by zero when `err==0`.**
  `src/shared/meds_numerics.f90:76` — a perfectly-converged step raises SIGFPE under Debug
  (`-fpe0`). *Fix:* clamp `err` to a tiny floor before forming the step ratio.
- **[REPORTED] `thomas_solve` divides by `b(1)`/`denom` with no guard.**
  `src/shared/meds_numerics.f90:37` — a singular/near-singular tridiagonal silently yields Inf/NaN
  that propagate into every hydrology/energy solve that calls it. *Fix:* guard the pivots and return a
  status flag the callers already have machinery to act on.
- **[PLAUSIBLE] Modified-Picard convergence tests a head increment against a moisture tolerance.**
  `meds_column_hydrology.f90:333` — `opts%atol` is reused for three dimensionally-incompatible roles
  (a volumetric-moisture error-norm floor, a matric-head Picard tolerance, and a mass tolerance), so
  the Picard stopping test compares Δψ [m head] against a θ [m³/m³] tolerance. *Impact:* convergence is
  effectively mis-scaled; whether it stops too early or too late depends on soil. *Fix:* separate
  `atol` into head/moisture/mass tolerances with correct units.
- **[PLAUSIBLE] Ground evaporation is potential (RH=100%) with no soil resistance.**
  `meds_column_energy.f90:243` — `ground_surface_balance` uses a saturated-surface humidity with only
  the aerodynamic conductance, no `α_soil`/DSL resistance — even though the hydrology module *has* a
  DSL formulation. *Impact:* bare-soil evaporation overestimated in dry conditions; inconsistent with
  the water side. *Fix:* apply the α/β soil-resistance (share the hydrology DSL).
- **[PLAUSIBLE] CAS humidity twin fully explicit while the enthalpy twin is implicit.**
  `meds_column_energy.f90:271` — the specific-humidity update is bare forward-Euler (no
  `1+dt·wci·gatm` denominator, no lower clamp) while enthalpy is implicit in the atm term. Asymmetric
  stability; humidity can undershoot. *Fix:* make the humidity twin implicit-in-atm like the enthalpy
  and CO₂ twins.
- **[PLAUSIBLE] `effective_growth` conflates negative carbon-fed growth with the 'unset' sentinel.**
  `meds_demography_rates.f90:261` tests `growth_avg < 0.0` to mean "unset", but the true sentinel is
  `GROWTH_AVG_UNSET = -1.0`. A genuinely shrinking cohort in carbon mode (negative growth_avg) is
  misread as unset. *Fix:* test against the explicit sentinel, not sign.

### 2.3 Low (guards, edge cases, minor physics)

> **✅ Resolution:** 16 fixed, 12 not-real/deferred in **PR #31** (+3 regression tests) — see Appendix C.

Compact list — all worth a guard/`error stop`, individually low impact:

| Finding | Location | Conf. |
|---|---|---|
| C4 `kp` also appears with pressure-only scaling in a second path | `meds_leaf_gas_exchange.f90:227` | CONFIRMED |
| `smaller_root` divides by θ with no θ=0 guard | `meds_leaf_gas_exchange.f90:122` | REPORTED |
| Katul optima below `g0` floored → breaks A–gs–Ci consistency | `meds_leaf_gas_exchange.f90:274` | PLAUSIBLE |
| TPU temperature response reuses Vcmax activation/deactivation params | `meds_leaf_gas_exchange.f90:211` | REPORTED |
| Integration finishing in exactly `max_substep` steps reported not-converged | `meds_plant_hydraulics.f90:290` | PLAUSIBLE |
| Step force-accepted at `h≤h_floor` with `err>1` still reports converged | `meds_plant_hydraulics.f90:283` | PLAUSIBLE |
| `rhizosphere_cond` units don't resolve to kg/s/MPa from documented inputs | `meds_plant_hydraulics.f90:143` | REPORTED |
| Soil-grid generator divides by zero when `grid_growth=0` → NaN depths | `meds_soil_parameters.f90:151` | PLAUSIBLE |
| `soil_thermal_cond` uses unfrozen log10 Kersten number for frozen soil | `meds_soil_thermal.f90:33` | CONFIRMED |
| Wood film-evap area (`effarea_evap·WAI`) inconsistent with wood sensible (`π·WAI`) | `meds_column_energy.f90:212` | REPORTED |
| `down0(ncoh+1)=incid_beam` set before `has_beam` guard → phantom top-layer interception | `meds_twostream_band.f90:56` | PLAUSIBLE |
| Thermal emission uses one blended `canopy_temp`, not separate leaf/wood temps | `meds_canopy_radiation.f90:76` | PLAUSIBLE |
| `beta_params_from_mean` ÷ variance, no lower floor → NaN LIDF at `std_deg=0` | `meds_optics.f90:135` | CONFIRMED |
| Fixed-substep hydrology path discards per-step convergence, always reports converged | `meds_column_hydrology.f90:238` | CONFIRMED |
| `rebuild_csr` no guard on `np==0` / `owner_patch∈[1,np]` → silent corruption in Release | `meds_demography_types.f90:460` | REPORTED |
| `carbon_flux_block` is a parallel per-cohort SoA no lockstep routine reorders (latent desync) | `meds_demography_types.f90:131` | REPORTED |
| `site_free` guards a 25-array deallocate on one unrelated array's status | `meds_demography_types.f90:152` | REPORTED |
| Persistent global ids/counters use `ik`=int32 though wide `lk`=int64 exists | `meds_demography_types.f90:122` | REPORTED |
| DAMM `sx_total` divides by `damm%depth_cm` with no positive guard | `meds_column_co2.f90:152` | REPORTED |
| `n_fast_per_slow` rounds `dt_slow/dt_fast` without enforcing exact divisibility | `meds_config.f90:154` | REPORTED |
| `quadratic_smaller_root` divides by `2θ` with no θ=0 guard | `meds_numerics.f90:61` | REPORTED |
| `toml_logical` fixed-position prefix compare mis-parses some boolean spellings | `meds_toml.f90:147` | REPORTED |
| GDD/chilling reset zeroes sums on cold out-of-season days (unlike ED2) | `meds_phenology.f90:125` | REPORTED |
| `req_dur` silently falls back to 1 day when a present duration string won't parse | `meds_config_io.f90:117` | REPORTED |
| Allometry coeffs `protected` with no default/installed-guard (indeterminate if used pre-`set_allometry`) | `meds_allometry.f90:39` | REPORTED |
| Census malformed rows abort in pass 1 but are silently `cycle`d in pass 2 | `meds_init.f90:135` | REPORTED |
| `fuse_2_patches` zero-combined-area early return still lets caller drop donor plants | `meds_demography_fusefiss.f90:458` | REPORTED |
| `has_nan` self-comparison may optimize to always-false; only checks dbh/nplant | `meds_demography_diagnostics.f90:102` | REPORTED |

---

## 3. Performance & numerical-solver issues

> **✅ Resolution:** 6 fixed (bit-identical, incl. the nvfortran `cstr()` correctness bug), 6 deferred in **PR #32** — see Appendix C.

**[MEDIUM · CONFIRMED] Whole cohort state migrates host↔device every daily step.**
`src/demography/meds_demography_dynamics.f90:75,78` — `growth_step`/`mortality_step` `map` ~18 full
cohort arrays each step, and 8 **write-only** derived arrays (`height, basal_area, agb, leaf_area,
leaf_carbon, fineroot_carbon, wood_carbon, nonstructural_carbon`) are mapped `tofrom`, doubling their
per-step transfer. On GPU the daily loop is migration-bound. *Fix:* mark the write-only arrays `map(from:)`;
wrap the daily loop in an `!$omp target data` region so the cohort SoA stays device-resident across
steps (this is the "reserved follow-up" already noted in `CLAUDE.md`).

**[MEDIUM · CONFIRMED] Hydraulics recomputes identical frozen coefficients for full- and half-step.**
`src/plant/meds_plant_hydraulics.f90:277` — the step-doubling error estimate recomputes the frozen
capacitance and the Kirchhoff Gauss–Legendre integral from the *same start state* for both the full
step and the first half step. *Fix:* compute the start-state coefficients once and reuse for both.

**[MEDIUM · REPORTED] Fusion/split loops re-allocate scratch and re-sort the whole site every pass.**
`meds_demography_fusefiss.f90:213` — each tolerance-relaxation iteration re-sorts and re-allocates.
*Fix:* hoist scratch out of the loop; sort once and maintain order incrementally.

**[MEDIUM · REPORTED] Allocation churn inside the per-patch / per-substep fast loop.**
`src/driver/meds_column_dynamics.f90:180` — allocatables created/destroyed every substep. *Fix:*
pre-allocate per-patch scratch in the `fast_context_t`.

**[MEDIUM · correctness · REPORTED] Array-valued `cstr()` result passed straight into netCDF calls.**
`src/io/meds_io.f90:64` — exactly the nvfortran array-temp miscompile trap `CLAUDE.md` warns about
(silent wrong values at `-O2`/`-fast`). *Fix:* bind `tmp = cstr(...)` to a named array first. (Flagged
as perf-unit but it is a portability/correctness bug on the NVHPC back end.)

**[CONFIRMED] Ross G-function & beam extinction recomputed per band though band-independent.**
`src/biophysics/meds_optics.f90:280` — `gfun_direct` (13-class Ross-G trig) is called once per band
inside the band loop with identical `cosz`. *Fix:* compute `gfun` once per cohort/PFT before the band
loop.

**[CONFIRMED] Harmonic face conductivity `kf` recomputed identically in two routines.**
`meds_column_energy.f90:80` — `soil_heat_be_step` and `soil_energy_flux` compute the same face
conductivity. *Fix:* compute once and pass.

**[Others]**
- Adaptive step-doubling does 3 implicit solves (6 `face_and_sink` sweeps) per accepted step for an
  L-stable integrator — `meds_column_hydrology.f90:252` (CONFIRMED, low): consider a cheaper embedded
  error estimate.
- Per-leaf `Ci` solved by ~19-iteration bisection where the diffusion regimes admit a closed-form
  root — `meds_leaf_gas_exchange.f90:258` (REPORTED, low–med).
- Temperature-dependent air properties recomputed per cohort in the boundary-layer loop —
  `meds_canopy_aerodynamics.f90:97` (PLAUSIBLE, low): hoist to per-patch.
- `cohort_reorder` materializes ~27 heap gather-temporaries per call; `cohort_compact` heap-allocates
  `perm` every call — `meds_demography_types.f90:318` (REPORTED, low).
- Pass-1 census patch discovery is O(N_rows·N_patches) with repeated whole-array realloc —
  `meds_init.f90:138` (REPORTED, low).

---

## 4. Code organization

> **✅ Resolution:** 10 fixed (bit-identical), 6 deferred/not-real in **PR #33** — see Appendix C.

Measured against the project's own rules (`CLAUDE.md` "Modernization guidelines").

**[MEDIUM · CONFIRMED] `meds_temp_response` pulls in all of `meds_config` for two integer selectors.**
`src/shared/meds_temp_response.f90:18` — a leaf-math kernel gains a whole-config prerequisite. *Fix:*
pass the two selectors as arguments (or a tiny options type), keeping the math module dependency-free.

**[MEDIUM · CONFIRMED] Public stateless energy kernels re-implemented inline in the coupler.**
`src/biophysics/meds_column_energy.f90:31` — `meds_column_dynamics` imports only `soil_energy_flux`
and re-derives `canopy_air_update`/`ground_surface_balance`/`veg_energy_balance` inline with different
formulations. Two sources of truth for the same physics. *Fix:* call the public seam kernels; delete
the inline copies.

**[MEDIUM · REPORTED] The "one centralized field list" is actually hand-duplicated across six routines.**
`src/state/meds_demography_types.f90:9` — the header advertises a single lockstep field list as the fix
for ED2's "forgot to reallocate" class, but the full 27-field per-cohort enumeration is copy-pasted in
`copy`/`reorder`/`compact`/`move_alloc`/`ensure_capacity`/`copy_cohort_slot`. Adding a field means
editing six lists (and `carbon_flux_block`/`gpp_accum` are already easy to miss — see §2.3). *Fix:*
generate the per-array operations from one macro/include list, or a single "for each field do X"
codegen, so the field set has one authority.

**[MEDIUM · REPORTED] `meds_column_dynamics` reimplements shared `sat_vapor_pressure` locally.**
`src/driver/meds_column_dynamics.f90:389` — a local `sat_e` duplicates `meds_thermo`'s
`sat_vapor_pressure` while the module already imports that function's derivative from `meds_thermo`.
*Fix:* use the shared function.

**[MEDIUM · REPORTED] Hydraulic node count `3` defined in three independent places.**
`src/biophysics/meds_biophysics_types.f90:440` hardcodes `psi(3,n_coh)`; `meds_demography_types.f90:83`
carries `N_HYDRO_NODE` vs `N_HYDRO`, linked only by a comment. *Fix:* one named constant in `shared`,
used everywhere.

**[Lower-severity organization]**
- Beer–Lambert `light_ext` lives in and is installed by `meds_allometry` (a structural-geometry
  module) — `meds_allometry.f90:47` (REPORTED): optical extinction belongs with RT/optics.
- Canopy-interception kernel lives in the soil-column hydrology module —
  `meds_column_hydrology.f90:48` (REPORTED): one-responsibility breach.
- `meds_optics` bundles leaf-angle + canopy + surface responsibilities — `meds_optics.f90:17`
  (REPORTED, low).
- Ring-buffer SMA eviction duplicated verbatim in `growth_step` and `apply_carbon_npp` —
  `meds_demography_dynamics.f90:156` (CONFIRMED, low): extract one helper.
- `meds_kinds` imported three times; fatal message names the old split module `meds_hydro_solver` —
  `meds_plant_hydraulics.f90:10,361` (CONFIRMED, low).
- Duplicate/unused `use` statements (module-merge artifacts) — `meds_leaf_gas_exchange.f90:9`
  (CONFIRMED, low).
- `meds_biophysics_types` aggregates four process domains; header comment stale —
  `meds_biophysics_types.f90:3` (CONFIRMED, low).
- `GROWTH_AVG_UNSET` value and its detection threshold live in different modules —
  `meds_demography_types.f90:34` / `meds_demography_rates.f90:28` (CONFIRMED coupling; note the
  "duplicated magic number" framing was **refuted** — they are a sentinel *value* −1.0 vs a *threshold*
  0.0, semantically distinct; still worth co-locating with a comment linking them).
- Module headers reference removed `build_config` / relocated `meds_plant_vital_rates` —
  `meds_config.f90:6` (REPORTED, low): stale docs.
- `column_co2_budget_t` / `damm_params_t` comments disagree with the code (loss2atm sign; whether
  `alpha_sx` absorbs depth) — `meds_biogeochem_types.f90:48,58` (REPORTED, low).
- Local variable named `target` shadows the Fortran `target` attribute keyword — `meds_init.f90:114`
  (REPORTED, low).
- `meds_io.f90:298` is 134 columns (>132 limit) — (REPORTED, low).

---

## Appendix A — coverage matrix

| Unit | Files | Machine-verified? | Notable findings |
|---|---|---|---|
| bio_rt | canopy_radiation, twostream_band, optics | ✅ | gfun per-band (perf), NaN LIDF, phantom beam |
| bio_hydro | column_hydrology, soil_parameters | ✅ | theta_res floor, Dunne runoff, Picard unit-mismatch |
| bio_energy | column_energy, soil_thermal, thermo | ✅ | soil-advection sign, gatm factor, resid-identically-zero |
| bio_aero | canopy_aerodynamics, biophysics_types | ✅ | um frozen pre-MO-iter, stale header |
| plant_leaf | leaf_gas_exchange | ✅ | **C4 kp temp**, transpiration gb, guards |
| plant_hydraulics | plant_hydraulics | ✅ | convergence-report bugs, frozen-coeff recompute |
| dem_dynamics | dynamics, rates, interface | ✅ | phantom-growth SMA, device migration, include_pft |
| plant_crp | carbon_dynamics, respiration, phenology | ⚠️ unverified | phenology gate, water-mode selector unimpl. |
| dem_fusefiss | fusefiss, diagnostics | ⚠️ unverified | **carbon-pool destruction** (hand-verified) |
| dem_state | demography_types (state) | ⚠️ unverified | field-list duplication, id int32 |
| allom_init | allometry, init | ⚠️ unverified | census header drop, light_ext placement |
| shared_num | numerics, temp_response, time, config, pft_params | ⚠️ unverified | thomas/adaptive guards, temp_response↔config |
| driver | main, stepper, fast_loop, column_dynamics, veg_dynamics | ⚠️ unverified | **fast-loop wiring, Rmaint** (hand-verified) |
| biogeochem | column_co2, biogeochem_types | ⚠️ unverified | pool non-conservation, silent Q10 fallback |
| io | io, netcdf_c, config_io, toml | ⚠️ unverified | **toml zero-fill** (hand-verified), nvfortran trap |
| org_sweep | whole tree | ⚠️ unverified | sat_e dup, hydro-node constant, 132-col |

To finish the machine verification of the ⚠️ units, resume the workflow
(`resumeFromRunId: wf_f1de6c1a-a4a`) — the finders replay from cache and only the verify/synth tail
re-runs.

## Appendix B — investigated and dismissed (REFUTED by adversarial verification)

- **Canopy interception over-reports evaporation** (`meds_column_hydrology.f90:60`) — capping
  `e_canopy` at `leaf_water/dt + q_grab` is the caller's documented contract (design §3c), not a bug.
- **`solve_plant_water` dt=0 divide-by-zero** (`meds_plant_hydraulics.f90:302`) — `dt=0` never reaches
  this path in the wiring.
- **Non-converged adaptive exit dilutes fluxes over the full `dt`** (`meds_plant_hydraulics.f90:290`) —
  the fluxes are formed over the actually-advanced interval; the dilution claim misreads the code.
- **`zeta` not clamped inside the M-O iteration → runaway-stable decoupling**
  (`meds_canopy_aerodynamics.f90:153`) — `zeta_max_stable=0.5` (CLM5) is applied; no runaway.
- **`apply_recruitment` `/mon_per_yr` mis-integration** (`meds_demography_dynamics.f90:283`) — correct
  for the only possible wiring (`apply_recruitment` runs strictly on the monthly `do_cohort_fissfuse`
  cadence).
- **Growth-unset "duplicated magic number"** (`meds_demography_rates.f90:28`) — a sentinel *value*
  (−1.0) and a detection *threshold* (0.0) are semantically distinct, not one convention duplicated.

---

*Generated from a 16-unit multi-agent adversarial review (31 machine verdicts) plus hand-verification
of the top-severity findings in the units whose verification was interrupted. Confidence tags are
conservative: treat REPORTED items as leads to confirm, CONFIRMED/verified-by-hand as actionable.*

---

## Appendix C — Resolution log (2026-07-08)

Every finding was re-verified against the source and dispositioned. **FIXED** = behavior-preserving
change landed and validated on ifx + nvfortran multicore + nvfortran GPU (24/24 each; carbon/GPU/IO
paths validated end-to-end). **DEFERRED** = real but would change numerics or needs a large/
architectural rework (reason stated). **NOT-REAL** = already correct or invariant-protected.

### §2.1 Critical / high — PR #29
- Fast biophysics loop never invoked → **FIXED** (construct `fast_context_t` in `meds_main`; hard `error stop` if `fast_biophysics_on` without a wired context).
- Cohort fusion/fission destroy carbon-mode pools → **FIXED** (carbon-aware fuse/split: nplant-weight the four pools, anchor on conserved `wood_carbon` via `set_cohort_size_from_carbon`).
- Maintenance respiration omitted from carbon growth → **FIXED** (`net_carbon = gpp − resp_maint − growth_resp`).
- C4 `kp_eff` no temperature response → **FIXED** (Arrhenius `temp_response` on `kp`, ED2/Collatz).
- `toml_real_array` silently zero-fills unparseable entries → **FIXED** (check `iostat`, `error stop` with key/token).
- atm↔CAS conductance omits the `temp1` profile factor → **FIXED** (threaded `temp1`/`temp2` through the energy + CO₂ twins).
- Dunne runoff spuriously large under free-drain → **FIXED** (gated to the SIMTOP-aquifer BC).
- Soil interior heat advection cell/sign → **verified CORRECT** (upwind on `w_flux` matches ED2 `qw_flux_g`; refuted, see Appendix B).

### §2.2 Medium — PR #30
- **FIXED (10):** `theta_res` floor routed through deficit accounting; Picard convergence on |Δθ| (dimensionally consistent); transpiration `gs`–`gb` **series** (+test); `effective_growth` tests the exact `GROWTH_AVG_UNSET` sentinel; `include_pft` gates **both** recruit paths; census header auto-detect (+test); unknown `hr_model` → NaN (not silent Q10); `adaptive_step_update` floors `err` before `err^(−½)`; actual→potential temperature into the aero solver; CAS humidity twin clamped ≥ 0.
- **NOT-REAL / by-design (7):** `thomas_solve` pivot (invariant-protected, pure); `energy_resid` (conservation-by-construction); ground-evap "potential" (test-only kernel; live uses the hydrology DSL); `fast_soil_carbon` (MVP prescribed pool); phenology NPP-funded growth (deliberate, tested); `split_cohorts` renorm (2nd-order in `split_eps`); `MAXPFT` bound (already fixed in #29).

### §2.3 Low — PR #31 (+3 regression tests)
- **FIXED — live correctness (5):** frozen-soil Kersten uses the linear `S_r` form, ice-aware blend (+test); Katul optimum-below-`g0` re-solve keeps A/gs/Ci consistent (+test); hydraulics `converged` false-negative/false-positive; fixed-substep convergence aggregation; `dt_slow` must be an integer multiple of `dt_fast`.
- **FIXED — latent guards (11):** `grid_growth=0` uniform-grid limit; no-beam phantom top-face beam; `beta_params_from_mean` variance floor (+test); DAMM `depth_cm` divide guard; two `smaller_root` θ-floors; `rebuild_csr` range guards; `fuse_2_patches` zero-area donor-cohort keep; `has_nan` → `ieee_is_nan` (+`agb`/`wood_carbon`); `rhizosphere_cond` unit annotation.
- **NOT-REAL / intended / deferred (12):** C4 `kp` "second path"; TPU reuses Vcmax T-response (faithful ED2); wood π/1 evap-area asymmetry (faithful ED2); `carbon_flux_block` reorder (transient input); `site_free` single guard (all-or-nothing alloc); int32 ids (overflow unreachable); `toml_logical` (already exact-match); `req_dur` / census pass-2 (already `error stop`); canopy thermal single-temp (MVP; P3); GDD/chill reset (intended ED2 divergence); allometry install-guard (init-order-safe).

### §3 Performance / solver — PR #32
- **FIXED — bit-identical (6):** **C1** `cstr()` nvfortran array-temp trap → typed `_f` wrappers binding to a named local (validated by live netCDF on **nvfortran** + ifx); **A1** write-only cohort arrays `map(from: arr(1:n))` (GPU transfer halved, tail-safe); **B1** hydraulics frozen-coefficient reuse; **B2** Ross-G computed once per cohort (not per band); **B3** `kf` reused across the two soil-energy routines; **C3** `transp_c` → automatic array.
- **DEFERRED (6):** B5 `Ci` closed form (not bit-identical); B6 hydrology embedded error estimator (not bit-identical); C2 fusefiss "sort once" (changes which cohorts fuse); A3 `cohort_reorder` temporaries (off the daily hot path, lockstep-critical); B4 per-cohort air properties (~6 flops); C4 census O(N_rows·N_patches) (one-time init).

### §4 Code organization — PR #33
- **FIXED — bit-identical (10):** **S1** reversed the `meds_temp_response`↔`meds_config` dependency edge (selectors owned by the kernel, re-exported); **S4** deleted the duplicate local `sat_e`, use shared `sat_vapor_pressure`; **S5** `bio%psi(3,…)` → shared `N_HYDRO_NODE`; **T1** `use meds_kinds` ×3 dedup + stale fatal message; **T2** leaf `use`-block dedup; **T3** stale `meds_biophysics_types` header; **T5** stale `build_config`/`meds_plant_vital_rates` doc refs; **T6** `column_co2_budget_t`/`damm_params_t` comment corrections; **T7** local `target` → `sel_site`; **T8** 134-col line wrapped.
- **DEFERRED / NOT-REAL (6):** **S2** inline energy kernels vs the seam kernels — **different formulations** (implicit-vs-explicit vapour twin + a CO₂ twin the kernel lacks; hydrology-authoritative vs self-derived ground latent; diagnostic-vs-prognostic leaf) → unifying is a **modeling change, not a refactor**; **S3** SoA field-list codegen (needs the C preprocessor on the lockstep-critical file + rank-1/2 X-macros across 3 compilers); **T4** `GROWTH_AVG_UNSET` (**not real** — already a single named-constant equality); **T9** `light_ext` / **T10** interception module moves (cross/strain the library DAG for a single consumer); **T11** SMA-eviction dedup (the `growth_step` copy runs in an `!$omp target` region → a shared helper needs `declare target`, conflicting with the arithmetic-only-kernel invariant).
