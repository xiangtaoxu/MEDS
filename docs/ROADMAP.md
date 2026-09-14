# MEDS roadmap

**What this file is.** Every piece of work MEDS has deliberately deferred, in one place. It
replaces the "reserved follow-ups" sections and the scattered `deferred` / `MVP` / `placeholder`
comments that used to be the only record of intent — those were findable only by grepping, and
several of them had quietly become true while nobody updated the comment.

**The rule.** A source comment states present-tense rationale. What changed and when goes in
[`CHANGELOG.md`](../CHANGELOG.md). What is deferred goes here, **and every item carries a GitHub
issue number** — so a deferred decision is findable, assignable and closable rather than a sentence
in a document nobody greps.

Items #150–#200 were filed on 2026-09-13 from this file. When you add an item here, open the issue
in the same change.

**Status vocabulary.** *Planned* — intended, with a design. *Candidate* — worth doing, no design
yet. *Open question* — the decision itself is undecided. Nothing here is a commitment to a date.

---

## 1. Known defects and open issues

The issues that predate this roadmap. Everything in sections 2–12 has an issue too (#150–#201,
filed 2026-09-13); these are listed separately because they were filed from measurement or from a
failure, not from a plan.

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

- **P3 — thread the real drivers into the phenology kernel.** *Planned.* [#150](https://github.com/xiangtaoxu/MEDS/issues/150) The slow driver
  currently hard-codes `avail_water = 0`, `dmax_leaf_psi = 0` and `rad = 0`, and uses
  `temp_day` as a proxy for soil temperature. The fast loop now produces all three. Needs:
  a soil-water running mean, a daily-**maximum** leaf water potential (a fast-loop daily
  reduction), a running-mean radiation, and the shallow soil-layer temperature in place of the
  air-temperature proxy.
- **P3 — the four cue-state columns.** *Planned.* [#150](https://github.com/xiangtaoxu/MEDS/issues/150) `pheno_water_avg`, `pheno_low_psi_days`,
  `pheno_high_psi_days`, `pheno_light_avg` must become cohort structure-of-arrays columns (they
  are per-cohort state today only inside a routine, so they are re-zeroed every day). Adding a
  column means the lockstep reorder, every creation site, and the fusion blend.
- **P3 — lift the config rejection.** [#150](https://github.com/xiangtaoxu/MEDS/issues/150) `validate_config` rejects the WATER, HYDRO and LIGHT cue
  bits. Lift each as its driver lands. Acceptance: the tropical drought-deciduous and
  light-driven leaf-exchanging strategies run from configuration alone.
- **P4 — `retained_carbon_fraction`.** *Planned.* [#151](https://github.com/xiangtaoxu/MEDS/issues/151) The carbon trait for resorption on leaf shed,
  with the full-removal closure. Optional companion: `root_phen_factor`.
- **P5 — one solar declination.** *Candidate.* [#152](https://github.com/xiangtaoxu/MEDS/issues/152) `solar_cosz` and `daylength` currently use two
  different declination formulae (Cooper 1969 and White 1997). Unify on one
  `solar_declination(doy)` with a deliberate golden re-baseline.
- **Open question** ([#150](https://github.com/xiangtaoxu/MEDS/issues/150)): whether flush should match shed for the light-driven
  leaf-exchanging strategy, rather than the present fixed high `k_flush_max`.

---

## 3. Soil biogeochemistry

Source: `docs/dev_plans/MEDS_BIOGEOCHEMISTRY_DESIGN.md` §7. Science page:
[`science/soil_carbon.md`](science/soil_carbon.md).

- **Make DAMM reachable, or delete it.** *Open question.* [#153](https://github.com/xiangtaoxu/MEDS/issues/153) `heterotrophic_respiration_damm`
  exists and is tested, but has no TOML key and no caller. Either add the `hr_model` selector
  the design specifies, or remove the kernel. A tested-but-unreachable kernel is the worst of
  both.
- **The nitrogen twin.** *Planned.* [#154](https://github.com/xiangtaoxu/MEDS/issues/154) Shaped in already: `n_cycle_on` is parsed and the N fields
  exist on `soil_carbon_t`, but no kernel reads the flag and the restart skips the fields.
  Needs the decomposition N limitation (`f_decomp`), the mineralization/immobilization flux, and
  N limitation of NPP on the plant side.
- **Vertically resolved soil-carbon pools.** *Candidate.* [#155](https://github.com/xiangtaoxu/MEDS/issues/155) The `n = 1` layer is the degenerate
  case of the intended profile.
- **An explicit coarse woody debris pool.** *Candidate.* [#156](https://github.com/xiangtaoxu/MEDS/issues/156) Today treefall necromass enters the
  structural pools directly.
- **Fire consumption of the ground pools.** *Candidate.* [#157](https://github.com/xiangtaoxu/MEDS/issues/157) No fire anywhere in MEDS.

---

## 4. Fast-loop numerics

Source: `docs/dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §5–§8. Science page:
[`science/numerical_scheme.md`](science/numerical_scheme.md).

- **N5 — an adaptive freeze cadence with a real error estimator.** *Planned.* [#158](https://github.com/xiangtaoxu/MEDS/issues/158) The last
  remaining efficiency item in the plan: decide per step how long the frozen coefficients stay
  valid, rather than freezing for exactly one `dt_fast`.
- **Fold soil water into the ARK tableau.** *Planned.* [#159](https://github.com/xiangtaoxu/MEDS/issues/159) — successor to
  #93, now closed; the pond is already on the state vector (Phase 0). Measured cost of in-stage soil water is +14–26 %, not
  the +5 % first estimated.
- **Warn when RK45 runs at production cadence.** *Planned.* [#160](https://github.com/xiangtaoxu/MEDS/issues/160) The transpiration corrector that
  fixed a ~1 MPa `psi_leaf` error is ARK-only. RK45 at 900 s therefore carries an error the
  default path does not, and nothing says so.
- **The `rwc_floor` clamp artefact** ([#104](https://github.com/xiangtaoxu/MEDS/issues/104)). *Open question.* A floored relative water
  content maps to a potential of about −10⁴ MPa, which is not a pressure any tissue reaches.
  The detector ships; whether to clamp the potential, arrest the solve, or kill the cohort is a
  physics decision. **Deferred to v0.3.0** (2026-09-13). Note that *arresting* is not the free
  option it looks: the collapsed store diagnoses ψ at about −1.5×10⁴ MPa against a soil at perhaps
  −2 MPa, so the cohort recovers today — that enormous artificial gradient IS the 13× cost — and
  removing uptake would make a transiently desiccated cohort permanently dead.
- **E5 — the RK45 rescue snapshot.** *Candidate.* [#161](https://github.com/xiangtaoxu/MEDS/issues/161) The rescue currently re-runs the step on ARK
  from the last accepted state.
- **`psi_leaf` is the one state that does not converge at 900 s.** *Open question.* [#162](https://github.com/xiangtaoxu/MEDS/issues/162) Its error is
  inherited from the canopy air and amplified roughly 4×; the residual relocates to `psi_wood`
  through the frozen uptake seam. Every other state and flux converges.

Source: `docs/dev_plans/MEDS_NUMERICS_SCOPING.md`.

- **Mark BB2/BB3 as refuted** in that document's §7, which still lists them as committed. The
  GPU evaluation refuted them. [#194](https://github.com/xiangtaoxu/MEDS/issues/194)
- **MB2 — soil-energy substepping.** *Candidate.* [#163](https://github.com/xiangtaoxu/MEDS/issues/163) The adaptive knobs exist on `energy_opts_t`
  but are never read. Re-verify that the need is real before building.
- **Bare-array forms** ([#164](https://github.com/xiangtaoxu/MEDS/issues/164), §11.3) for `cas_column_step_implicit`, `soil_energy_step_implicit`,
  `soil_carbon_step` and the snow kernels, so they match the device-eligible convention the
  other kernels follow.

---

## 5. Vegetation energy

Source: `docs/dev_plans/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §6–§7. Science page:
[`science/vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md).

- **A separate canopy film store with phase change.** *Planned.* [#165](https://github.com/xiangtaoxu/MEDS/issues/165) Intercepted water currently
  has no independent thermal state and cannot freeze.
- **Retire `veg_energy_step_implicit`.** *Done* ([#166](https://github.com/xiangtaoxu/MEDS/issues/166)). This entry was **stale**: the kernel was
  deleted in PR #120, and `veg_energy_diagnostic` does not exist either. `veg_energy_balance` is
  the single closure, diagnostic at `store_hcap_per_dt = 0` and prognostic above it.
- **The free-convection slope.** *Candidate.* [#167](https://github.com/xiangtaoxu/MEDS/issues/167) The 1.25·h term is absent from the linearization.
- **Honest wood sizing.** *Done, 2026-09-13* ([#168](https://github.com/xiangtaoxu/MEDS/issues/168)). This entry was **stale**: `bsap` stopped
  being a placeholder in PR #125. `set_cohort_wood_geometry` derives it from ED2's real
  `b1SA`/`b2SA` sapwood-area allometry. Measured, the old `0.10 * wood_carbon` placeholder made
  the wood thermal time constant 6.5–10× too short across the whole size range.

---

## 6. Output and diagnostics

Source: `docs/dev_plans/MEDS_IO_V01_PLAN.md` §4.5–§4.6, §6; `MEDS_IO_DESIGN.md` §3.5. Science
page: [`science/diagnostics.md`](science/diagnostics.md).

- **Mortality carbon by pathway.** *Planned.* [#169](https://github.com/xiangtaoxu/MEDS/issues/169) Background, cull and disturbance separately. The
  slow-loop ledger now values mortality per phase, which is the seam this needs.
- **Disturbance area flux.** *Planned.* [#170](https://github.com/xiangtaoxu/MEDS/issues/170) `PD_DISTURB_AREA` is declared as an accumulator slot but
  has no writer and no registry row.
- **Per-band albedo and up-welling shortwave and longwave.** *Planned.* [#171](https://github.com/xiangtaoxu/MEDS/issues/171) Needs a surface
  radiative-flux record.
- **Unify the FAST tier's extraction path.** *Planned.* [#172](https://github.com/xiangtaoxu/MEDS/issues/172) `fast_sample_t`, `extract_fast_scalar`
  and `output_integrate_fast` were slated for deletion when the FAST tier moved onto the general
  registry, and were not deleted. The tier still has a bespoke staging path.
- **Rename `[io]` to `[state]`.** *Candidate.* [#173](https://github.com/xiangtaoxu/MEDS/issues/173) The block now carries only restart settings; the
  diagnostic writer it was named for is gone. Several comments still say "legacy `[io]` path".
- **Variance output** (`AGG_MEANSQ`). *Candidate.* [#174](https://github.com/xiangtaoxu/MEDS/issues/174) The aggregation code exists and has no
  consumer; it needs a companion output slot.
- **An evaluation notebook** against the Ithaca test bed, and a PFT / size-class plotter in
  `post_proc/`. *Candidate.* [#175](https://github.com/xiangtaoxu/MEDS/issues/175)

---

## 7. Plant physiology

- **Thermal acclimation of photosynthesis and respiration.** *Planned.* [#176](https://github.com/xiangtaoxu/MEDS/issues/176) A running-mean tissue
  temperature shifting the peaked-Arrhenius reference. The seams (`t_acclim` on the wood and
  root environment records) were shaped in and never filled. Sources:
  `MEDS_PLANT_TRAIT_DYNAMICS_DESIGN.md` §5, `MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md` §17.
- **Storage-pool respiration.** *Candidate.* [#177](https://github.com/xiangtaoxu/MEDS/issues/177) Same source.
- **Vertically resolved root respiration.** *Candidate.* [#178](https://github.com/xiangtaoxu/MEDS/issues/178) One bulk soil temperature drives fine-root
  maintenance respiration today.
- **Per-PFT hydraulic traits.** *Candidate.* [#179](https://github.com/xiangtaoxu/MEDS/issues/179) The hydraulics parameters are PFT-uniform. Source:
  `MEDS_HYDRO_CURVE_EXTRACTION_DESIGN.md` §8.
- **Phase B — per-layer root nodes.** *Planned.* [#180](https://github.com/xiangtaoxu/MEDS/issues/180) The plant hydraulic network resolves the root
  system as one node against a weighted soil boundary. Sources:
  `MEDS_MULTILAYER_ROOTS_DESIGN.md` §5, `MEDS_HYDRAULICS_DESIGN.md` §16.
- **Hydraulic redistribution.** *Planned.* [#181](https://github.com/xiangtaoxu/MEDS/issues/181) Per-layer root efflux is floored at zero in both the
  plant solver and the soil sink, so uptake is non-negative by construction. Turning it on means
  deciding how the ledger treats water moving between layers through the plant. Same source.

---

## 8. Forcing

Source: `docs/dev_plans/MEDS_FORCING_DESIGN.md` §5.7, §8. Science page:
[`science/forcing.md`](science/forcing.md).

- **LWdown synthesis.** *Planned.* [#182](https://github.com/xiangtaoxu/MEDS/issues/182) `lwdown_source = "synthesize"` is currently rejected by
  `validate_config` because the Brutsaert/Idso clear-sky synthesis does not exist. Sources
  lacking longwave cannot drive MEDS until it does.
- **The multi-polygon runtime.** *Candidate.* [#183](https://github.com/xiangtaoxu/MEDS/issues/183) A grid → polygon → site state hierarchy, an array
  of readers, a polygon loop and MPI. Large and orthogonal to everything else.
  `nearest_grid_index` is the reusable atom, already built.
- **A transient or observed CO₂ stream.** *Planned.* [#184](https://github.com/xiangtaoxu/MEDS/issues/184)
- **Forcing-file global attributes are written but never read.** *Planned.* [#185](https://github.com/xiangtaoxu/MEDS/issues/185) The
  `sw_input_kind`, `timestep_seconds` and `avg_convention` attributes are produced by the prep
  script and ignored by the reader, so a file that disagrees with the config is undetected. Atmospheric CO₂ is a constant 420 ppm.

---

## 9. Snow

Source: `docs/dev_plans/MEDS_SNOW_DESIGN.md` §7. Science page:
[`science/snow_biophysics.md`](science/snow_biophysics.md).

- **P1 — multi-layer snow.** *Planned.* [#186](https://github.com/xiangtaoxu/MEDS/issues/186) Compaction and densification, an aging albedo, and a
  density-dependent thermal conductivity. The single-layer store is the degenerate case.
- **P2 — canopy snow interception and unloading.** *Planned.* [#187](https://github.com/xiangtaoxu/MEDS/issues/187) Snow currently reaches the ground
  through the canopy unimpeded.

---

## 10. Source structure and testing

Source: `docs/dev_plans/MEDS_CODE_STRUCTURE_DESIGN.md` §15.

- **Pass `column_params_t` through `column_config_t`** rather than copying it into the frozen
  record every step. *Planned.* [#188](https://github.com/xiangtaoxu/MEDS/issues/188)
- **Per-layer face budget imbalance on the committed path**, per-cohort tissue residuals, and
  RK45 ledgers asserted after the rail decision. *Planned.* [#189](https://github.com/xiangtaoxu/MEDS/issues/189)
- **Delete `column_cohort_t`** in favour of `cohort_fast_slice_t` / `patch_fast_slice_t` with a
  per-field policy table. *Deferred to v0.3.0, **paired with #146*** (2026-09-13).
  [#190](https://github.com/xiangtaoxu/MEDS/issues/190) 38 references across 11 files. Four of the
  five benefits the design claimed have since landed piecemeal: `column_cohort_init` gives the test
  fixtures allometric consistency, the three hard-coded constants are PFT parameters, the derived
  geometry is on the cohort block, and `reconcile_tissue_water_capacity` took the seed and clamp out
  of the gather. The fusion/scaling policy is already centralised in `fuse_cohort_fast_state` and
  `scale_cohort_ground_fields`. What is left is **completeness you cannot forget** — a table the
  blend iterates cannot omit a field a hand-written routine can — and that is #146's hazard class,
  which is why the two now travel together.
- **Consolidate the per-test `check` routines.** *Planned.* [#191](https://github.com/xiangtaoxu/MEDS/issues/191) Fourteen of the 46
  test files define their own.
- **A packed `column_state_t`** ([#146](https://github.com/xiangtaoxu/MEDS/issues/146)). *Deferred to v0.3.0, **paired with #190***.
  Only a packed layout makes field omission a compile-time error; there are 1 207 field references
  today. One packed, policy-carrying layout should serve the fast state vector and the cohort slice
  together — separately, each is a large refactor buying a fraction of one property.
- **Decide the year-rollover `seam[soil_carbon_rh]` residual** (§15.2): 8.37×10⁻⁴ kgC m⁻² at a
  year boundary. Either document it or mark it unmeasurable. *Open question.* [#192](https://github.com/xiangtaoxu/MEDS/issues/192)
- **Whether `test/` should mirror the source tree** (decision #13). Never decided; `test/` is
  flat. *Open question.* [#193](https://github.com/xiangtaoxu/MEDS/issues/193)
- **A frozen-seam contract note** (§15.4, Phase 3): the Λ criterion, the four-seam
  classification, and the arbitration rule. *Candidate*, marked optional in the plan. [#201](https://github.com/xiangtaoxu/MEDS/issues/201)

---

## 11. Performance

Source: `docs/dev_plans/MEDS_GPU_EVALUATION.md` §12.

- **Correct the GPU claims** [#194](https://github.com/xiangtaoxu/MEDS/issues/194) in `CMakeLists.txt` and the project documentation. The build
  comments still describe the offload as the parallel path; the measurement says otherwise.
  *Planned.*
- **Attack the allocator traffic.** *Planned.* [#195](https://github.com/xiangtaoxu/MEDS/issues/195) About 24 % of fast-loop self time is allocator
  work in `build_column_frozen`.
- **Thread and vectorise the cohort axis on the CPU.** *Planned.* [#196](https://github.com/xiangtaoxu/MEDS/issues/196) This is the evaluation's
  headline recommendation. Patch-axis threading already ships; the cohort axis is untouched.
- **A single-precision experiment** (`wp = real32`). *Candidate.* [#197](https://github.com/xiangtaoxu/MEDS/issues/197) Measured 56× on the device for
  FP32, which is why it is worth knowing what MEDS actually needs.

---

## 12. Smaller residues

- **Retire the IMEX-Euler oracle tier**, or record the decision to keep it
  (`archive/MEDS_DRIVER_REORG_DESIGN.md` §7). *Open question.* [#198](https://github.com/xiangtaoxu/MEDS/issues/198)
- **Remove the dead `--parity` preset** in `scripts/numerics_sweep.py`: it pins two config keys
  that no longer exist, so the preset cannot run. *Planned.* [#199](https://github.com/xiangtaoxu/MEDS/issues/199)
- **The demography C-API shim inlines allometry** rather than calling it, and the Python
  empirical laws hard-code two allometry coefficients. *Planned.* [#200](https://github.com/xiangtaoxu/MEDS/issues/200)
