# Changelog

All notable changes to MEDS. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
MEDS does not yet promise semantic versioning, because v0.x is explicitly pre-benchmark.

**This file is the only place change history lives.** Source comments state present-tense
rationale; what changed and when goes here; what is deferred goes in
[`docs/ROADMAP.md`](docs/ROADMAP.md). See `CLAUDE.md` for the rule.

Entries cite the pull request that shipped them. Anything that moves a number states the
before and after.

---

## [Unreleased]

### Fixed

- **The tissue-water floor now reports the water it creates** (#148). `advance_water_mass_full`
  floors `leaf_water_mass` and `wood_water_mass` at a tiny positive value so the linear mass Euler
  step cannot go negative — and creates water doing so. The code's own comment called the case
  "unobserved in this pass's test scenarios", which was **a belief, not a measurement**: the floor
  fires per cohort per tissue, while the whole-column water ledger sums leaf + wood over all
  cohorts, so water created in one cohort's wood is indistinguishable from a redistribution between
  cohorts. The budget closed to ~4×10⁻¹² kg m⁻² with the floor entirely unmonitored. It is now
  reported through the existing `budget%clamp_mass` / `clamp_commit_n` commit-clamp channel, which
  already reduces to `site%work_clamp_mass` and already has an output variable — so it surfaces
  end-to-end with no new reporting surface, which is what the issue asked for. Measured on a forced
  fixture: 2 activations creating 2.92×10⁻² kg m⁻², and **zero on an ordinary step**. The claim in
  `ark2_column_step` that its commit counter "stays 0 by construction" was corrected — that claim
  was the whole of the defect.
- **`PD_DISTURB_AREA` is written and emitted** (#170). The slot was declared in the patch
  diagnostic block with no writer and no registry row, so the disturbed-area flux read as a **silent
  zero** rather than a missing variable — the harder failure to notice, and one no conservation
  check can see, because zero disturbed area is a perfectly conservative answer. Written on the
  donor patches inside `apply_patch_disturbance` *before* the gap is appended, which is where the
  diag slots still line up with the donors; writing after the append is the out-of-bounds trap this
  file has already paid for once. New output variable `disturb_area_site`.
- **A run selecting `time_integrator = "rk45"` above `dt_fast` = 300 s is now warned** (#160). The
  transpiration corrector that cut a ~1 MPa `psi_leaf` error by 314× (PR #91) lives in
  `advance_water_mass_full`, which ARK calls and RK45 does not, so an RK45 production run silently
  carried an error the default path does not. A warning rather than an error, because RK45 is the
  accuracy baseline and is meant for a fine step.
- **Γ\* now responds to `o2_mol_frac`, so the O₂ knob propagates to both places oxygen enters the
  C3 demand** (#117). The compensation point is set by Rubisco's CO₂/O₂ specificity and is
  proportional to the O₂ partial pressure, but it was computed with no O₂ dependence at all — so
  raising O₂ correctly inhibited carboxylation through `Kc(1 + O/Ko)` while leaving the entire
  photorespiratory penalty on `Aj` untouched. Scaled by `o2_mol_frac / 0.209`, the O₂ the shipped
  `gstar25` was measured at (Bernacchi et al. 2001). **The factor is exactly 1 at the default, so
  every shipped configuration is bit-identical** (10.18931 µmol m⁻² s⁻¹ before and after on the test
  fixture). Away from it, measured on a midday tropical leaf: halving O₂ raises A_net by 22.0 %
  against 8.2 % before, and 35 % O₂ cuts it by 23.7 % against 10.5 %.
- **Two roadmap items described code that no longer exists**, found by re-measuring every filed
  premise before scheduling it (#168, #166). `bsap` stopped being a placeholder in PR #125 —
  `set_cohort_wood_geometry` derives it from ED2's real `b1SA`/`b2SA` sapwood-area allometry. The
  old `0.10 * wood_carbon` placeholder made the wood thermal time constant **6.5–10× too short**
  across the whole size range (`f_sap` runs 1.00 at dbh ≤ 19 cm to 0.655 at 117 cm, against 0.10).
  `veg_energy_step_implicit` was deleted in PR #120, and `veg_energy_diagnostic` does not exist
  either; `veg_energy_balance` is the single closure. Both entries corrected in `docs/ROADMAP.md`
  and their four stale source/doc references removed.
- **`scripts/numerics_sweep.py` could not run a cross-scheme comparison.** The `--parity` preset
  pinned three config keys that no longer exist (`fast.integration_scheme`,
  `fast.leaf_energy_model`, `fast.wood_energy_model` — #199 named two), and the `SCHEMES` table
  still offered `split` and `picard`, which are now a hard error. The preset is removed rather than
  repointed: every difference it pinned has since been closed by making the schemes agree.
- **The C-API demography shim inlined the allometry** instead of calling `meds_allometry`, so a
  coefficient change updated the model and not the shim (#200). It now calls `dbh_to_height`,
  `dbh_to_agb`, `dbh_to_leaf_area`, `size2leaf_carbon` and `size2wood_carbon`.
  `examples/example_demography/empirical_laws.py` reads the four allometry coefficients from the
  `[allometry]` block of its shipped PFT config instead of hard-coding them; the values are
  unchanged, so the example's behaviour is unchanged.

### Documentation

- **The frozen-seam contract is written down** (#201):
  [`docs/dev_plans/MEDS_FROZEN_SEAM_CONTRACT.md`](docs/dev_plans/MEDS_FROZEN_SEAM_CONTRACT.md). The
  Λ = F·dt/S criterion and its case split (linear-in-store is scale-free and sound; a prescribed
  flux against a prognostic store is unsound, with the +30.7 µmol m⁻² s⁻¹ NEE scar to prove it), the
  four seams classified against it, why Λ is meaningless for a *rate* seam and what replaces it
  (debit-before-credit), and the arbitration rule — scale all demands by `min(1, S/D)` — for a store
  with several consumers, which is order-independent where per-process clamping is not.
- **The ED2 two-stream defects found during the port are folded into
  [`docs/ed2_comparison.md`](docs/ed2_comparison.md) §5a** (#6), with the MEDS ↔ ED2 RT structure
  mapping. All six are reported upstream. Two of them — the stale-PFT diffuse index and the missing
  clumping factor in the longwave split — are *impossible by construction* in MEDS, and that is why
  its RT is shaped the way it is.
- **`test/` stays flat, deliberately** (#193, structure-plan decision #13), recorded in
  `src/README.md` with the reason: the discipline that matters is the link line, not the directory.
- **The year-rollover `rh_seam_gap` residual is attributed** (#192). 8.370×10⁻⁴ kgC m⁻² at a year
  boundary is the annual patch cadence, not a leak: the fast loop accumulates against one patch
  composition and the daily step debits another, blended through a matrix nonlinear in the lignin
  fraction. Documented at the field, with the instruction not to widen a tolerance to absorb it.
- **The build files no longer describe GPU offload as the parallel path** (#194). Measured, the
  offload build runs **1.4× slower** than the CPU (49.8 s against 36.2 s), one kernel sits at 0.4 %
  occupancy, and the device treated as 20 slow cores is 26× slower than 4 CPU cores. Patch-axis CPU
  threading is the parallel path. `MEDS_GPU=gpu` is kept as a reproducible experiment.
  `CMakeLists.txt`, `docs/building.md`, `src/README.md` and `docs/ed2_comparison.md` corrected.
- **A v0.2.0 release plan** in [`docs/dev_plans/MEDS_V02_RELEASE_PLAN.md`](docs/dev_plans/MEDS_V02_RELEASE_PLAN.md):
  all 66 open issues triaged into six phases plus release, 48 in scope and 18 deferred to v0.3+,
  with twelve decisions recorded.
- **The documentation was reorganized against the restructured source tree.** Thirty design
  plans moved to `docs/dev_plans/archive/` with tombstones; eleven stay live or as reference.
  The README dropped from 310 lines to a reader's entry point, with building, configuration
  and post-processing moved to their own pages. A new [`src/README.md`](src/README.md)
  documents the source layout, the placement rules and the library graph.
  `CLAUDE.md` split into a short always-loaded file plus path-scoped rules. This file and
  `docs/ROADMAP.md` were created. Three missing science pages were written:
  `docs/science/soil_carbon.md`, `forcing.md`, `plant_respiration.md`.

### Changed

- **The IMEX-Euler oracle tier is retired** (#198). It was not an oracle: it returned no reference
  trajectory, and `imex_euler_column_step` was a two-line wrapper around `column_be_stage` plus
  `advance_water_mass_full` — the production scheme's own kernels — so its independence was in the
  tableau, not in the machinery it was supposed to check. All five of its consumers kept their
  coverage: four now call a test-local `be_euler_step` (the same two-line composition, living in the
  code that degrades the scheme), and `test_adaptive_march` was **ported to `adaptive_ark_march`**,
  so it now exercises the production controller instead of one that existed only to serve the tier
  (8 sub-steps at `rtol` 1e-3 against 25 at 1e-6). `meds_fast_rk4_oracle` holds one oracle and its
  name is accurate again.
- **The RK4 oracle's independence is now structural.** With the tier gone, the module no longer
  imports `meds_fast_be_stage` at all — no `column_be_stage`, no `newton_surface_solve`, no
  `advance_water_mass_full`. It sees the pure right-hand side and the state algebra and nothing
  else, which is what makes agreement between it and an implicit scheme rule out a shared-bug false
  pass. Before, the module imported the BE machinery for the tier's benefit while the oracle itself
  never touched it.
- **One implementation of each test assertion helper** (#191). Twenty-three local copies across
  nineteen files, consolidated behind generic interfaces so all ~1000 call sites compile unchanged;
  net −403 lines. The two families (fatal condition-first, accumulating name-first) are kept
  deliberately — they differ in failure behaviour, not just signature. A new `meds_test_assert`
  holds the assertions and depends on nothing but a kind, because seventeen tests deliberately link
  one narrow library each. Found and closed one coverage hole: `test_biogeochem_dynamics`' local
  helper error-stopped while the shared one accumulates, which turned four assertions into no-ops
  until a verdict call was added.

### Fixed

- **`soil_carbon_on` now actually defaults to on.** PR #144 flipped the in-type default to
  `.true.` but the TOML loader still passed `.false.` as its absent-key default, and
  `toml_logical` returns the supplied default when the key is absent. Every config that
  omitted the key, including the shipped `meds_config_main.toml`, therefore ran with soil
  carbon off. Off is not a coarser soil-carbon model, it is no soil carbon: litter is
  discarded and `rh = 0`.
- **`forcing.lwdown_source = "synthesize"` is rejected rather than silently ignored.** The
  value parsed to `LW_SYNTHESIZE` but the met driver always read LWdown from the file, so it
  selected the file path under a name promising Brutsaert/Idso synthesis. `validate_config`
  now stops on it until the synthesis exists.

### Changed

- Test coverage for the three untested fast-loop state combinators, with two silent
  exclusions named (#147).
- `soil_carbon_on` default flipped to `.true.` in the type (#144). Measured over a year at
  Ithaca, off reports annual-mean NEE at −4.657 against −2.464 µmol m⁻² s⁻¹, an 89 % stronger
  apparent sink, with Rh identically zero against 0.833 kgC m⁻² yr⁻¹, for no wall-clock saving.
- `[energy].phase_change` deleted; ice-aware soil thermal properties are unconditional (#143).
  The flag only ever gated ice-aware `κ_sat(f_liq)` and `C_eff(f_liq)`; the difference was
  ≤ 0.61 K and 35 692 substeps either way.

### Added

- **The full coupled model is drivable from Python** (`meds.model.Run`): open, step, finalize
  through a C-API shim (#139). The Python and executable paths differ only by libm-versus-libimf
  interposition: worst relative difference ~1×10⁻¹² over a simulated day.
- **One `libmeds.so`** built by scikit-build-core, with a mandatory ctest target per C-API shim
  so an ABI change is a build failure in a default build rather than a silent break in an
  optional one (#138).
- **A slow-loop conservation ledger**, with birth and death as paired transfers (#132).

### Fixed (slow-loop conservation, #132–#137)

The ledger closed a series of real carbon, energy and water leaks that a green test suite and
the per-store fast budgets had both missed:

- **Growth respiration was 23 % of growth carbon that was never exhaled.** The allocator's
  outputs now all have a destination (#133).
- Tissue thermal mass is an exchange, not an appearance (#134).
- The time-averaged diagnostics now close for energy and water everywhere (#135).
- Mortality is valued on what the applier actually removed, not on what was requested (#136).
- Reproduction carbon became a flow rather than a disappearance, which closed the ledger
  (#137). Side effect: the recruit pool now accrues every slow step instead of as a monthly
  lump, so the first cohorts appear about a month later and the early trajectory is offset.
  The offset decays as the stand fills: 70 % in AGB at year 2, 1.6 % at year 10, 0.1 % at
  year 40.
- **Soil respiration reached the atmosphere 964× too small** — a CENTURY Rh unit error, hidden
  because no shipped config turned soil carbon on. The seam check that would have caught it was
  computed but never reported (#140).
- An out-of-bounds litter read: a per-patch quantity indexed outside the patch block when
  disturbance added a patch mid-step. The ledger declared the same garbage it consumed, so
  conservation balanced on it (#139).

### Changed (source-tree restructure, #125–#127, #141)

The tree is now **timescale-first for processes** over a **state layer in two halves**. Moving a
file changes only CMake wiring, because Fortran `use` is by module name, so every step was
verified byte-identical on both back ends.

- Steps 1–6: state as a layer (`state/column`, `state/site`), processes timescale-first
  (`fast_dynamics/`, `slow_dynamics/`) (#125). `src/shared/` dissolved into `base/`,
  `functions/`, `config/`, `state/`.
- Steps 9–10: every façade module dissolved; the fast-loop argument vocabulary renamed so
  names tell the truth — `wcap`/`ccap` → `cas_mass_capacity`/`cas_molar_capacity`,
  `gah`/`gaw`/`gac` → `g_atm_heat`/`g_atm_vapour`/`g_atm_co2`,
  `uext_to_temp`/`temp_to_uext` → `internal_energy_to_temp`/`temp_to_internal_energy`, and the
  three different things called `hydro` split into `soil_water_opts` versus
  `hydraulics_params`/`hydraulics_opts` (#126). netCDF registry strings and TOML keys were
  deliberately not renamed, so output files and configs are unchanged.
- Steps 0 and 7: the fast-loop state vector (#127).
- The library name `demography` was renamed to `core` in July and back to `demography` here:
  `core` stopped carrying information once its state half became `state/site`.
- Structure plan §15 phased remainder, Phase 1 executed (#141).

### Changed (configuration, #129, #130)

- **The soil column comes from config, not from literals in the driver** (#129): the new
  `[soil_column]` block carries layer count, depth, grid growth, hydraulic texture, retention
  family, root profile and the three thermal properties, validated at load. `depth` is the knob
  for the known too-shallow-column defect (2.0 m against a ~2.5 m annual damping depth).
- **The fast driver's hard-coded parameters found their homes** (#130). `agf_bs` duplicated the
  per-PFT `aboveground_frac` with a hard-coded 0.7, so a run whose PFTs differed in allocation
  used their values everywhere except stem respiration (#128). Canopy optics — leaf and wood
  reflectance and transmittance per band, clumping, leaf-angle mean and standard deviation —
  became `[pft]` traits, so two PFTs can finally differ in how they intercept light (#131).
  Longwave is configured as emissivity, with reflectance derived as `1 − emissivity` and
  transmittance zero, because a leaf is opaque at thermal wavelengths.

### Fixed (2026-09 review, #119–#124)

- Conservation fixes: a closed whole-column energy ledger, RK45 and threading defects,
  film-water and seam conservation (#119). A 50-year spin-up's energy imbalance went from
  +6.75 MJ m⁻² to 3.6×10⁻¹¹ W m⁻².
- Dead code removed, shared column-state algebra extracted, `t_ground` made explicit, one CAS
  box kernel instead of two (#120).
- The fast-loop argument vocabulary renamed; over-long lines wrapped (#122).
- `column_prepass` split into the five processes it fused; a per-PFT leaf photosynthesis table;
  scalar signatures; shared assemblers; ledger-consistent latent and sensible heat (#123).
- The frozen work record decomposed by physical content; `integrator_opts_t`;
  `apply_process_mask` (#124).

---

## [0.1.0] — 2026-08-02

First tagged release. **Unbenchmarked**: no EDTS-equivalent regression suite has been run, no
site compared flux-for-flux, no output scored against observations. What is verified is
internal — the test suite on two compilers, per-step conservation ledgers, and thread-invariant
output.

### Added

- **The v0.1 diagnostic output layer** (#111): per-variable, per-timescale output control.
  ~208 registered variables across 8 groups and 7 axes (cohort, patch, site, soil, PFT,
  DBH size class, and the 2-D soil×patch slab), each switchable individually per timescale
  from TOML. Extensive quantities carry their own aggregation weight, so one registry line
  emits a cohort field's patch, site, PFT and size-class rollups. `meds_main --dump-io-config`
  lists every variable. Verified byte-identical at 1 versus 4 threads across all 75 files of a
  3-year run.
- **The ED2-to-MEDS comparison** for people who already run ED2 (#113, #115),
  [`docs/ed2_comparison.md`](docs/ed2_comparison.md).
- **Patch-axis threading** (#109): byte-identical output at 1, 2, 4 and 8 threads. 2.03× at
  4 threads against a measured 3.03× hardware ceiling. Three silent compiler traps were found
  and worked around: Intel's `-auto-scalar` default placing local arrays in static storage,
  nvfortran rejecting `BLOCK` inside a parallel region, and ifx building `private` copies of a
  derived type through a compiler-generated static mold.

### Changed

- **`dt_fast` became an accuracy parameter, not a stability one** (#90). The period-2
  canopy-air oscillation was traced to one frozen coefficient — the canopy-air-to-atmosphere
  conductance — and re-solving it at every integrator stage removed the oscillation while
  *reducing* integrator work. The production default went from 150 s to **900 s**. The
  non-stomatal water-stress limb was gated off by default in the same change.
- The legacy `[io]` diagnostic writer was retired at v0.1: an annual-cadence, instantaneous,
  hard-coded 21-variable schema that duplicated `[output]` and collided with it on the `-D-`
  filename prefix. `meds_io` now carries the state (restart) stream only.

### Fixed

- **The transpiration seam in the plant water update** (#91): a pure flux inconsistency,
  exactly `dt·(transp_pp − transp_bw)`, corrected by a transpiration corrector. Improved
  `psi_leaf` convergence 314×. ARK only; RK45 does not carry the corrector.
- Stomatal water-stress closure (#98, issue #95): `beta_stomata` was identically 1 until fixed.
  Restart persistence for the daily-maximum leaf water potential.
- Pathological hydraulics sub-stepping is now detected (#105): a collapsed, floored wood store
  burned 13× the wall clock, silently (issue #104).
- The adaptive warm start got per-patch storage (#108, issue #106) — it had been loop-carried,
  so patch 1 cold-started and patches 2..N inherited a neighbour's step size.

### Documentation

- **GPU offload evaluated against measurement and found not viable as scoped** (#110). The GPU
  build ran 1.4× *slower* than the CPU (49.8 s against 36.2 s), with one kernel at 0.4 %
  occupancy; the GPU treated as 20 slow cores was 26× slower than 4 CPU cores. The
  recommendation is to thread and vectorise the cohort axis on the CPU instead.

---

## Pre-0.1 development — 2026-06-23 to 2026-08-02

MEDS began on 2026-06-23. The sections below group the first six weeks by subsystem rather than
by date, because the work proceeded as a dozen parallel subsystem builds.

### Demographic core

- Source tree organized by process domain; per-PFT `hgt_max` (#3).
- The engine became **carbon-driven**: carbon pools on the cohort structure-of-arrays (#17),
  carbon-driven growth wired in (#19), the empirical vital-rate laws moved out to the Python
  example and the engine reduced to law-free apply-primitives (#16, #18, #43).
- The demography engine reorganized into `core` with a tendency-seam growth/mortality
  interface (#44, #45).
- **Run-model decision** (#63): the fast biophysics loop is always on, and phenology is
  unconditional. A slow-only, empirical-demography run is the Python C-API path.

### Fast-loop biophysics

- **Canopy radiative transfer** (#5): a faithful modernized ED2 two-stream (`icanrad = 2`) with
  a unified multi-band solver, SCOPE/4SAIL leaf-angle scattering over a Beta leaf-angle
  distribution, and an in-house block-tridiagonal solve. ED2 bugs found during the port are
  recorded in issue #6.
- **Plant hydraulics** (#9, #49): a stateless matrix-exponential network solver, shared
  constitutive curves (pressure-volume, Kirchhoff conductance, xylem vulnerability) extracted
  into `meds_hydr_lib`, and multi-layer root water uptake.
- **Leaf gas exchange** (#2, #46): FvCB C3 and Collatz C4 demand, Leuning / Medlyn / Katul
  stomatal models, a bracketed C_i solver, and Sabot two-limb water stress.
- **Column soil hydrology** (#22, #23): implicit backward-Euler Thomas Richards with Celia
  modified-Picard linearization, upstream-weighted conductivity, adaptive step-doubling,
  infiltration and ponding, and a free-drain / bedrock / aquifer bottom boundary.
- **Energy balance** (#24, #25): four stateless per-store thermal kernels carrying prognostic
  **internal energy rather than temperature**, so the freeze/thaw plateau is a read-off of the
  enthalpy inverter rather than a special case.
- **Canopy air space CO₂** (#26, #27): `can_co2` as the third prognostic CAS twin, with a DAMM
  heterotrophic-respiration option.
- **Fast-loop coupling capstone** (#28): the sub-daily loop coupled and owning its own state.
- **Prognostic leaf and wood energy** (#41) with a separate wood temperature; FAST diurnal
  diagnostics in the same change.
- **Snow and temporary surface water** (#42): stateless kernels and a conserving fast-loop
  coupling that closes whole-column mass and energy budgets to machine precision through
  accumulation, sublimation, melt into infiltration, and the snow-albedo ramp.

### Numerics

- **IMEX-ARK** (2026-07-11, e77f7f0, no PR number): an L-stable ESDIRK2 on the ARS(2,2,2)
  tableau with an arrowhead Newton surface solve. The config string stays `"ark"`; the explicit
  part is empty, so despite the name it is a diagonally implicit scheme.
- **Error-control infrastructure** (#65) and tolerance unification, a process mask, a sweep
  harness (#66).
- **ED2-faithful RK45** (#67, #68): an adaptive Cash-Karp 5(4) march, internal water carried as
  **mass** rather than potential, and a mass-conserving demographic seam.
- **Integrator parity** (#77, #80–#87): a long sequence making the schemes solve the same
  physics. Along the way: the reference used to score them was mis-timed; ARK and RK45 dropped
  snowfall all winter; `veg_coupling_floor` destroyed energy; the split scheme's soil-energy
  budget was a tautology; the pond held no enthalpy so water crossing into it shed its heat
  while the books still closed.
- **The operator-split integrator was retired** (#88). It converged to a different limit
  (~0.45 K in canopy-air temperature) that refinement never removed and nobody ever attributed,
  and it could not carry the coupled tissue heat store. `ark` became the default and `rk45` the
  accuracy baseline, with the RK45 stiff rescue redoing the step on `ark`. The tissue heat
  store was turned on in the same change, integrated by an exact exponential with two weights
  (endpoint and step-average) because the tissue ODE is linear under the frozen coefficients.

### Slow loop

- **Phenology**: a stateless signal module (#10), wired into the run loop (#39), then refactored
  to a **rate-based** signal-only kernel emitting two relative tendencies instead of a
  directional tri-state (#51).
- **Plant carbon**: pure carbon-dynamics kernels (#14), carbon-pool allometry and PFT traits
  (#15), then an **elemental growth-allocation kernel** following FATES PARTEH Hypothesis-1
  (#52), with growth respiration charged inside the kernel on realized growth.
- **Trait plasticity** (#54): light-driven acclimation of specific leaf area, V_cmax, dark
  respiration and leaf lifespan.
- **Non-leaf maintenance respiration** (#13): stem and fine-root, an ED2 Chambers-2004 port.
- **Slow soil-carbon biogeochemistry** (#35): ED2's CENTURY decomposition reorganized as the
  carbon matrix ODE `dX/dt = B·I + A·ξ·K·X`, with a 7-pool state, a lignin sub-tracer, an exact
  augmented matrix exponential for accelerated steps, and a SASU steady-state solve. Wired into
  the slow loop (#64): the fast loop's heterotrophic respiration respires the same frozen pool
  the daily step debits, so the day's total fast Rh equals the daily debit by construction.

### Forcing

- **Meteorological forcing** (#36): a single-site NetCDF reader over a multi-grid `(time, grid)`
  file, `pure`/`elemental` disaggregation kernels including an interval-mean-conserving
  shortwave reconstruction, the canopy-RT join, net longwave, Weiss-Norman band-specific
  shortwave, multi-year calendar recycling, nearest-grid matching, and wind-height and
  elevation lapse. MEDS never gap-fills: a missing or NaN required value is a hard error.
- **The recycle window became declared rather than inferred** (#69). The previous classifier
  accepted only Jan-1 00:00 files and silently fell back to an absolute-seconds span wrap
  otherwise. Real ERA5-Land records are stamped at the *end* of each interval, so their first
  record is 01:00:00 and they always took that fallback: the Ithaca file spans 366 d 22 h, and
  a 29-year run ended up reading late May at a ~10 h offset. Nothing caught it, because the
  cosz reconstruction is mean-conserving, so daily-mean shortwave stayed correct and the slow
  demography looked healthy.

### Output and I/O

- **Diagnostic aggregation and output subsystem** (#38, #40): the registry, the per-tier
  temporal integrators, and the netCDF serializer, written through the netCDF **C** library via
  `iso_c_binding` so the output layer builds under ifx and nvfortran (netCDF-Fortran's module
  format is gfortran-only).
- A **FAST (sub-daily) output tier** (2026-07-13, 01347bb, no PR number) for diurnal-cycle
  analysis.

### Code review and bug fixes

- **Adversarial code review, 2026-07-06**, over ~10.4k lines and 47 modules, targeting
  physical-process bugs, numerical defects and organization. All sections addressed across
  #29–#34: 8 critical/high physical-process bugs, then medium, then low-priority guards, then
  performance and solver issues, then organization.

### Tooling and infrastructure

- CMake with automatic Fortran module-dependency resolution — the deliberate fix for ED2's
  "run `make` six times" hack and per-platform `include.mk` files.
- **nvfortran portability trap documented** (#8, issue #7): never pass an array-valued function
  result straight into a call. nvfortran's whole-program optimizer miscompiles the temporary
  descriptor — silently wrong values at `-O2`, segfault at `-O0` — while
  `ifx -stand f18 -check all` tolerates it. A green ifx run is not sufficient.
- 3D visualization of forest structure (#4); the biophysics example (#70).
- Design plans relocated to `docs/dev_plans/` (#48); GitHub math rendering fixed in the science
  pages (#50).

---

[Unreleased]: https://github.com/xiangtaoxu/MEDS/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/xiangtaoxu/MEDS/releases/tag/v0.1.0
