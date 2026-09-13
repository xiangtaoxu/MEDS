# MEDS roadmap

**What this file is.** Every piece of work MEDS has deliberately deferred, in one place. It
replaces the "reserved follow-ups" sections and the scattered `deferred` / `MVP` / `placeholder`
comments that used to be the only record of intent — those were findable only by grepping, and
several of them had quietly become true while nobody updated the comment.

**The rule.** A source comment states present-tense rationale. What changed and when goes in
[`CHANGELOG.md`](../CHANGELOG.md). What is deferred goes here, with a pointer to the issue or the
design section that carries the detail.

**Status vocabulary.** *Planned* — intended, with a design. *Candidate* — worth doing, no design
yet. *Open question* — the decision itself is undecided. Nothing here is a commitment to a date.

---

## 1. Known defects and open issues

Tracked on GitHub; listed here so the roadmap is one read rather than two.

| # | Title | Note |
|---|---|---|
| [#1](https://github.com/xiangtaoxu/MEDS/issues/1) | Equal-height shading ambiguity: capped large trees and same-height recruits | Enhancement |
| [#6](https://github.com/xiangtaoxu/MEDS/issues/6) | ED2 two-stream RT bugs found during the port, and the MEDS ↔ ED2 mapping | Upstream record |
| [#7](https://github.com/xiangtaoxu/MEDS/issues/7) | nvfortran miscompiles array-valued function results passed as actual arguments | Standing toolchain rule, not a fix |
| [#47](https://github.com/xiangtaoxu/MEDS/issues/47) | Reconcile leaf water-stress (`beta_stomata` / `beta_nonstomata`) with ED2's Manzoni-style form | |
| [#74](https://github.com/xiangtaoxu/MEDS/issues/74) | Condensate is deposited into soil layer 1, not onto leaf/wood surface water | |
| [#89](https://github.com/xiangtaoxu/MEDS/issues/89) | Saturated vapour calculation | |
| [#96](https://github.com/xiangtaoxu/MEDS/issues/96) | Dynamic vapour pressure for leaf transpiration (Kelvin $`e_i`$) | Built, measured, removed; likely route to foliar water uptake |
| [#104](https://github.com/xiangtaoxu/MEDS/issues/104) | Plant hydraulics burns 13× wall clock on a collapsed (floored) wood store | Detector shipped (#105); the physics decision is open — see §4 |
| [#114](https://github.com/xiangtaoxu/MEDS/issues/114) | Comparison of ED2 and MEDS v0.1.0 | |
| [#117](https://github.com/xiangtaoxu/MEDS/issues/117) | $`\Gamma^*`$ does not respond to `o2_mol_frac`, so the O₂ knob only half-propagates | |
| [#118](https://github.com/xiangtaoxu/MEDS/issues/118) | C3 co-limitation reuses `theta_j` (the J-hyperbola curvature); C4 has its own. Costs ~30 % of A | |
| [#145](https://github.com/xiangtaoxu/MEDS/issues/145) | Soil thermal bottom boundary: an adiabatic base reflects the annual wave; no column depth in 2–3 m is converged | `[soil_column].depth` is the knob; 2→3 m moved the base layer by −8.5 K |
| [#146](https://github.com/xiangtaoxu/MEDS/issues/146) | Fast-integrator state vector: only a packed layout gives compile-time omission safety (1 207 field references) | |
| [#148](https://github.com/xiangtaoxu/MEDS/issues/148) | Tissue-water floor in `advance_water_mass_full` creates water with no ledger term | |

---

## 2. Phenology — finish the cue set

Source: `docs/dev_plans/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md` §9. The kernel already computes
all five cues; what is missing is the wiring from the fast loop, so three of the four strategies
the model claims to support cannot actually be selected.

- **P3 — thread the real drivers into the phenology kernel.** *Planned.* The slow driver
  currently hard-codes `avail_water = 0`, `dmax_leaf_psi = 0` and `rad = 0`, and uses
  `temp_day` as a proxy for soil temperature. The fast loop now produces all three. Needs:
  a soil-water running mean, a daily-**maximum** leaf water potential (a fast-loop daily
  reduction), a running-mean radiation, and the shallow soil-layer temperature in place of the
  air-temperature proxy.
- **P3 — the four cue-state columns.** *Planned.* `pheno_water_avg`, `pheno_low_psi_days`,
  `pheno_high_psi_days`, `pheno_light_avg` must become cohort structure-of-arrays columns (they
  are per-cohort state today only inside a routine, so they are re-zeroed every day). Adding a
  column means the lockstep reorder, every creation site, and the fusion blend.
- **P3 — lift the config rejection.** `validate_config` rejects the WATER, HYDRO and LIGHT cue
  bits. Lift each as its driver lands. Acceptance: the tropical drought-deciduous and
  light-driven leaf-exchanging strategies run from configuration alone.
- **P4 — `retained_carbon_fraction`.** *Planned.* The carbon trait for resorption on leaf shed,
  with the full-removal closure. Optional companion: `root_phen_factor`.
- **P5 — one solar declination.** *Candidate.* `solar_cosz` and `daylength` currently use two
  different declination formulae (Cooper 1969 and White 1997). Unify on one
  `solar_declination(doy)` with a deliberate golden re-baseline.
- **Open question** (§10 Q2b): whether flush should match shed for the light-driven
  leaf-exchanging strategy, rather than the present fixed high `k_flush_max`.

---

## 3. Soil biogeochemistry

Source: `docs/dev_plans/MEDS_BIOGEOCHEMISTRY_DESIGN.md` §7. Science page:
[`science/soil_carbon.md`](science/soil_carbon.md).

- **Make DAMM reachable, or delete it.** *Open question.* `heterotrophic_respiration_damm`
  exists and is tested, but has no TOML key and no caller. Either add the `hr_model` selector
  the design specifies, or remove the kernel. A tested-but-unreachable kernel is the worst of
  both.
- **The nitrogen twin.** *Planned.* Shaped in already: `n_cycle_on` is parsed and the N fields
  exist on `soil_carbon_t`, but no kernel reads the flag and the restart skips the fields.
  Needs the decomposition N limitation (`f_decomp`), the mineralization/immobilization flux, and
  N limitation of NPP on the plant side.
- **Vertically resolved soil-carbon pools.** *Candidate.* The `n = 1` layer is the degenerate
  case of the intended profile.
- **An explicit coarse woody debris pool.** *Candidate.* Today treefall necromass enters the
  structural pools directly.
- **Fire consumption of the ground pools.** *Candidate.* No fire anywhere in MEDS.

---

## 4. Fast-loop numerics

Source: `docs/dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §5–§8. Science page:
[`science/numerical_scheme.md`](science/numerical_scheme.md).

- **N5 — an adaptive freeze cadence with a real error estimator.** *Planned.* The last
  remaining efficiency item in the plan: decide per step how long the frozen coefficients stay
  valid, rather than freezing for exactly one `dt_fast`.
- **Fold soil water into the ARK tableau.** *Planned.* Tracked as issue #93 Phase 1; the pond is
  already on the state vector (Phase 0). Measured cost of in-stage soil water is +14–26 %, not
  the +5 % first estimated.
- **Warn when RK45 runs at production cadence.** *Planned.* The transpiration corrector that
  fixed a ~1 MPa `psi_leaf` error is ARK-only. RK45 at 900 s therefore carries an error the
  default path does not, and nothing says so.
- **The `rwc_floor` clamp artefact** (issue #104). *Open question.* A floored relative water
  content maps to a potential of about −10⁴ MPa, which is not a pressure any tissue reaches.
  The detector ships; whether to clamp the potential, arrest the solve, or kill the cohort is a
  physics decision.
- **E5 — the RK45 rescue snapshot.** *Candidate.* The rescue currently re-runs the step on ARK
  from the last accepted state.
- **`psi_leaf` is the one state that does not converge at 900 s.** *Open question.* Its error is
  inherited from the canopy air and amplified roughly 4×; the residual relocates to `psi_wood`
  through the frozen uptake seam. Every other state and flux converges.

Source: `docs/dev_plans/MEDS_NUMERICS_SCOPING.md`.

- **Mark BB2/BB3 as refuted** in that document's §7, which still lists them as committed. The
  GPU evaluation refuted them.
- **MB2 — soil-energy substepping.** *Candidate.* The adaptive knobs exist on `energy_opts_t`
  but are never read. Re-verify that the need is real before building.
- **Bare-array forms** (§11.3) for `cas_column_step_implicit`, `soil_energy_step_implicit`,
  `soil_carbon_step` and the snow kernels, so they match the device-eligible convention the
  other kernels follow.

---

## 5. Vegetation energy

Source: `docs/dev_plans/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §6–§7. Science page:
[`science/vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md).

- **A separate canopy film store with phase change.** *Planned.* Intercepted water currently
  has no independent thermal state and cannot freeze.
- **Retire `veg_energy_step_implicit`.** *Candidate.* The exact-exponential
  `veg_energy_diagnostic` is the one closure both paths share; the prognostic sibling is a
  second code path for the same physics.
- **The free-convection slope.** *Candidate.* The 1.25·h term is absent from the linearization.
- **Honest wood sizing.** *Planned.* `bsap` is still a placeholder in the vegetation driver,
  which makes the modelled wood time constant unrealistic even though the wood-area allometry is
  now ED2's real `b1WAI`/`b2WAI`.

---

## 6. Output and diagnostics

Source: `docs/dev_plans/MEDS_IO_V01_PLAN.md` §4.5–§4.6, §6; `MEDS_IO_DESIGN.md` §3.5. Science
page: [`science/diagnostics.md`](science/diagnostics.md).

- **Mortality carbon by pathway.** *Planned.* Background, cull and disturbance separately. The
  slow-loop ledger now values mortality per phase, which is the seam this needs.
- **Disturbance area flux.** *Planned.* `PD_DISTURB_AREA` is declared as an accumulator slot but
  has no writer and no registry row.
- **Per-band albedo and up-welling shortwave and longwave.** *Planned.* Needs a surface
  radiative-flux record.
- **Unify the FAST tier's extraction path.** *Planned.* `fast_sample_t`, `extract_fast_scalar`
  and `output_integrate_fast` were slated for deletion when the FAST tier moved onto the general
  registry, and were not deleted. The tier still has a bespoke staging path.
- **Rename `[io]` to `[state]`.** *Candidate.* The block now carries only restart settings; the
  diagnostic writer it was named for is gone. Several comments still say "legacy `[io]` path".
- **Variance output** (`AGG_MEANSQ`). *Candidate.* The aggregation code exists and has no
  consumer; it needs a companion output slot.
- **An evaluation notebook** against the Ithaca test bed, and a PFT / size-class plotter in
  `post_proc/`. *Candidate.*

---

## 7. Plant physiology

- **Thermal acclimation of photosynthesis and respiration.** *Planned.* A running-mean tissue
  temperature shifting the peaked-Arrhenius reference. The seams (`t_acclim` on the wood and
  root environment records) were shaped in and never filled. Sources:
  `MEDS_PLANT_TRAIT_DYNAMICS_DESIGN.md` §5, `MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md` §17.
- **Storage-pool respiration.** *Candidate.* Same source.
- **Vertically resolved root respiration.** *Candidate.* One bulk soil temperature drives fine-root
  maintenance respiration today.
- **Per-PFT hydraulic traits.** *Candidate.* The hydraulics parameters are PFT-uniform. Source:
  `MEDS_HYDRO_CURVE_EXTRACTION_DESIGN.md` §8.
- **Phase B — per-layer root nodes.** *Planned.* The plant hydraulic network resolves the root
  system as one node against a weighted soil boundary. Sources:
  `MEDS_MULTILAYER_ROOTS_DESIGN.md` §5, `MEDS_HYDRAULICS_DESIGN.md` §16.
- **Hydraulic redistribution.** *Planned.* Per-layer root efflux is floored at zero in both the
  plant solver and the soil sink, so uptake is non-negative by construction. Turning it on means
  deciding how the ledger treats water moving between layers through the plant. Same source.

---

## 8. Forcing

Source: `docs/dev_plans/MEDS_FORCING_DESIGN.md` §5.7, §8. Science page:
[`science/forcing.md`](science/forcing.md).

- **LWdown synthesis.** *Planned.* `lwdown_source = "synthesize"` is currently rejected by
  `validate_config` because the Brutsaert/Idso clear-sky synthesis does not exist. Sources
  lacking longwave cannot drive MEDS until it does.
- **The multi-polygon runtime.** *Candidate.* A grid → polygon → site state hierarchy, an array
  of readers, a polygon loop and MPI. Large and orthogonal to everything else.
  `nearest_grid_index` is the reusable atom, already built.
- **A transient or observed CO₂ stream.** *Planned.* Atmospheric CO₂ is a constant 420 ppm.

---

## 9. Snow

Source: `docs/dev_plans/MEDS_SNOW_DESIGN.md` §7. Science page:
[`science/snow_biophysics.md`](science/snow_biophysics.md).

- **P1 — multi-layer snow.** *Planned.* Compaction and densification, an aging albedo, and a
  density-dependent thermal conductivity. The single-layer store is the degenerate case.
- **P2 — canopy snow interception and unloading.** *Planned.* Snow currently reaches the ground
  through the canopy unimpeded.

---

## 10. Source structure and testing

Source: `docs/dev_plans/MEDS_CODE_STRUCTURE_DESIGN.md` §15.

- **Pass `column_params_t` through `column_config_t`** rather than copying it into the frozen
  record every step. *Planned.*
- **Per-layer face budget imbalance on the committed path**, per-cohort tissue residuals, and
  RK45 ledgers asserted after the rail decision. *Planned.*
- **Delete `column_cohort_t`** in favour of `cohort_fast_slice_t` / `patch_fast_slice_t` with a
  per-field policy table. *Planned.* 38 references across 11 files.
- **Consolidate the per-test `check` routines.** *Planned.* Fourteen of the 46 test files define
  their own.
- **A packed `column_state_t`** (issue #146). *Planned.* Only a packed layout makes field
  omission a compile-time error; there are 1 207 field references today.
- **Decide the year-rollover `seam[soil_carbon_rh]` residual** (§15.2): 8.37×10⁻⁴ kgC m⁻² at a
  year boundary. Either document it or mark it unmeasurable. *Open question.*
- **Whether `test/` should mirror the source tree** (decision #13). Never decided; `test/` is
  flat. *Open question.*
- **A frozen-seam contract note** (§15.4, Phase 3): the Λ criterion, the four-seam
  classification, and the arbitration rule. *Candidate*, marked optional in the plan.

---

## 11. Performance

Source: `docs/dev_plans/MEDS_GPU_EVALUATION.md` §12.

- **Correct the GPU claims** in `CMakeLists.txt` and the project documentation. The build
  comments still describe the offload as the parallel path; the measurement says otherwise.
  *Planned.*
- **Attack the allocator traffic.** *Planned.* About 24 % of fast-loop self time is allocator
  work in `build_column_frozen`.
- **Thread and vectorise the cohort axis on the CPU.** *Planned.* This is the evaluation's
  headline recommendation. Patch-axis threading already ships; the cohort axis is untouched.
- **A single-precision experiment** (`wp = real32`). *Candidate.* Measured 56× on the device for
  FP32, which is why it is worth knowing what MEDS actually needs.

---

## 12. Smaller residues

- **Retire the IMEX-Euler oracle tier**, or record the decision to keep it
  (`MEDS_DRIVER_REORG_DESIGN.md` §7). *Open question.*
- **Remove the dead `--parity` preset** in `scripts/numerics_sweep.py`: it pins two config keys
  that no longer exist, so the preset cannot run. *Planned.*
- **The demography C-API shim inlines allometry** rather than calling it, and the Python
  empirical laws hard-code two allometry coefficients. *Planned.*
