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
| [#74](https://github.com/xiangtaoxu/MEDS/issues/74) | Condensate is deposited into soil layer 1, not onto leaf/wood surface water | |
| [#96](https://github.com/xiangtaoxu/MEDS/issues/96) | Dynamic vapour pressure for leaf transpiration (Kelvin $`e_i`$) | Built, measured, removed; likely route to foliar water uptake |
| [#104](https://github.com/xiangtaoxu/MEDS/issues/104) | Plant hydraulics burns 13× wall clock on a collapsed (floored) wood store | Detector shipped (#105); the physics decision is open — see §4 |
| [#254](https://github.com/xiangtaoxu/MEDS/issues/254) | Passive deep **thermal** layers below the hydrologically active column | The Dirichlet anchor shipped and cut the base-layer amplitude error from +82 % to −2 %, but a purely resistive termination cannot reflect less than 0.41 — closing the rest needs heat *capacity* below the column, i.e. a thermal grid that extends past the water grid |
| [#146](https://github.com/xiangtaoxu/MEDS/issues/146) | Fast-integrator state vector: only a packed layout gives compile-time omission safety (1 207 field references) | |

---

## 2. Phenology

Source: `docs/dev_plans/archive/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md`. Science page:
[`science/plant_phenology.md`](science/plant_phenology.md).

All five cues are wired and selectable as of v0.2.0, and the design plan is archived. What is
**not** done is validation: no MEDS leaf-area cycle has been scored against an observation, at any
site, under any strategy. That is not a roadmap item with a design — it is the benchmarking gap in
§1, of which this is one instance.

- **`root_phen_factor`** — *Candidate.* [#258](https://github.com/xiangtaoxu/MEDS/issues/258) The fine-root side of phenological shedding.
  Leaf resorption shipped in v0.2.0 (`retained_carbon_fraction`, default 0); the fine-root coupling
  did not. Phenology stays leaf-only by design, so if adopted it belongs in the carbon layer beside
  the resorption split.

---

## 3. Soil biogeochemistry

Source: `docs/dev_plans/MEDS_BIOGEOCHEMISTRY_DESIGN.md` §7. Science page:
[`science/soil_carbon.md`](science/soil_carbon.md).

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
- **The `rwc_floor` clamp artefact** ([#104](https://github.com/xiangtaoxu/MEDS/issues/104)). *Open question.* A floored relative water
  content maps to a potential of about −10⁴ MPa, which is not a pressure any tissue reaches.
  The detector ships; whether to clamp the potential, arrest the solve, or kill the cohort is a
  physics decision. **Deferred to v0.3.0** (2026-09-13). Note that *arresting* is not the free
  option it looks: the collapsed store diagnoses ψ at about −1.5×10⁴ MPa against a soil at perhaps
  −2 MPa, so the cohort recovers today — that enormous artificial gradient IS the 13× cost — and
  removing uptake would make a transiently desiccated cohort permanently dead.
- **`psi_leaf` is the one state that does not converge at 900 s.** *Open question.* [#162](https://github.com/xiangtaoxu/MEDS/issues/162) Its error is
  inherited from the canopy air and amplified roughly 4×; the residual relocates to `psi_wood`
  through the frozen uptake seam. Every other state and flux converges.

Source: `docs/dev_plans/MEDS_NUMERICS_SCOPING.md`.

- **Bare-array forms** ([#164](https://github.com/xiangtaoxu/MEDS/issues/164), §11.3) for `cas_column_step_implicit`, `soil_energy_step_implicit`,
  `soil_carbon_step` and the snow kernels, so they match the device-eligible convention the
  other kernels follow.

---

## 5. Vegetation energy

Source: `docs/dev_plans/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §6–§7. Science page:
[`science/vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md).

- **A separate canopy film store with phase change.** *Planned.* [#165](https://github.com/xiangtaoxu/MEDS/issues/165) Intercepted water currently
  has no independent thermal state and cannot freeze.
- **The free-convection slope.** *Deferred to v0.3.0, premise re-measured.*
  [#167](https://github.com/xiangtaoxu/MEDS/issues/167) The design note says the true sensible-heat
  slope is `1.25·h` and the solved `ΔT_leaf` is overstated ~20 % in calm conditions. Measured, it is
  not. `1.25` is the **pure free-convection** limit: `H ∝ ΔT^{1+m}` holds only for the Grashof part
  of the Nusselt number, so the real slope factor is `1 + m·f_free` with
  `f_free = Nu_free/(Nu_forced+Nu_free)`. At the shipped `leaf_width = 0.04 m`, and at the
  `ugbmin = 0.25 m/s` in-canopy wind **floor** — the most free-convection-favourable state the model
  can reach — that factor is **1.02–1.07** over the leaf-minus-CAS temperature range an Ithaca run
  actually produces, and only **1.15** at an extreme 15 K. It reaches 1.25 only at exactly zero
  wind, which `ugbmin` makes unreachable.

  The effect on `ΔT_leaf` is smaller again, because `h_coeff` is one of four terms in
  `veg_energy_balance`'s denominator (`h_coeff + le_slope + le_slope_wet + lw_slope`).

  And the fix is not the one the note describes. A flat factor on the denominator alone would break
  the kernel's exact conservation identity (`store_hcap_per_dt*(dt_end − dt_prev) + denom*dt_avg ==
  numer`). The correct form is a **Newton linearization about the frozen `ΔT₀`**, which changes the
  numerator, the denominator, the reported flux `dh`, *and* the exponential relaxation time constant
  `τ = cap/denom` together — in the kernel with the most delicate conservation invariant in the
  model, whose adaptive controller has already been broken once by a discontinuity here. That is not
  a proportionate trade for a few percent of one denominator term.

---

## 6. Output and diagnostics

Source: `docs/dev_plans/MEDS_IO_V01_PLAN.md` §4.5–§4.6, §6; `MEDS_IO_DESIGN.md` §3.5. Science
page: [`science/diagnostics.md`](science/diagnostics.md).

- **A spectrally resolved surface radiative record.** *Candidate.* [#255](https://github.com/xiangtaoxu/MEDS/issues/255)
  Follow-up to #171: Per-band incident and upwelling fluxes
  shipped in v0.2.0 on the VIS/NIR/LW three-band grid the two-stream solves. Comparing against a
  multispectral product (MODIS bands, Sentinel-2) needs finer bands, which is a change to the RT's
  band structure rather than to its output.
- **Remove the `[io]` deprecation shim.** *Scheduled, post-v0.2.x.* The block was renamed to
  `[state]` in v0.2.0 (#173) with `[io]` still loading behind one warning. Drop the
  `req_*_renamed` readers and the warning once users have had a minor release to migrate.

---

## 7. Plant physiology

- **PER-COHORT thermal acclimation.** *Candidate.* [#256](https://github.com/xiangtaoxu/MEDS/issues/256) Follow-up to #176:
  Photosynthetic acclimation shipped in v0.2.0 as a **site-level** growth temperature (Kattge &
  Knorr 2007), which is what that fit is calibrated on — the mean air temperature of the preceding
  weeks. Letting a shaded understory cohort acclimate differently from a sunlit canopy one needs a
  relation calibrated on **leaf** rather than air temperature; applying the air-temperature fit to a
  per-cohort leaf temperature would use it well outside its range. Acclimation of **respiration**
  (as opposed to photosynthetic capacity) is also still open.
- **Whether storage maintenance should be ON by default.** *Open, v0.3.0 question.* The mechanism
  shipped in v0.2.0 (#177) with `storage_turnover_rate` defaulting to **0**, which reproduces the
  earlier behaviour. ED2's temperate-broadleaf value of 0.6243 yr⁻¹ costs an Ithaca run 35 % of GPP
  and 43 % of AGB over five years, so switching the default on is a rebaseline decision, not a
  parameter tweak — and it lands with the same question #118 and #176 raised about presets that were
  never climate- or biome-specific.
- **Vertically resolved root respiration** — *shipped in v0.2.0* (#178). Kept here only as the
  pointer that the remaining vertical-resolution gap is **root biomass**, not temperature: the
  response is now summed per layer against the root profile, but `broot` itself is still a single
  per-plant pool distributed by the static `root_beta` profile rather than a prognostic per-layer
  one.
- **Per-PFT hydraulic traits are shipped** (#179, v0.2.0); what remains is that nothing has
  *calibrated* them. The thirteen traits are selectable per PFT and default to the shared
  `[hydraulics]` block, so a run that does not set them is unchanged — which means MEDS ships with
  hydraulically identical PFTs until somebody puts real trait values in. Wood density is the obvious
  axis to derive them from (denser wood ⇒ more negative `wood_psi50`, lower `wood_kmax`), as the
  Camac mortality coefficients already are.
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

- **A better cloud term for the LWdown synthesis.** *Candidate.* [#257](https://github.com/xiangtaoxu/MEDS/issues/257)
  Follow-up to #182: The Brutsaert/Idso synthesis shipped in
  v0.2.0, so a source lacking longwave can now drive MEDS. Its cloud correction is one empirical
  coefficient on the SW clearness index, and it holds the last daytime index through the night —
  adequate for a fallback, but the residual is real: driving Ithaca from synthesis leaves the soil
  surface 1.37 K cooler than the file's `strd`. A cloud-fraction formulation (Crawford & Duchon
  1999) or a nocturnal index carried from a longer window would close more of it.
- **The multi-polygon runtime.** *Candidate.* [#183](https://github.com/xiangtaoxu/MEDS/issues/183) A grid → polygon → site state hierarchy, an array
  of readers, a polygon loop and MPI. Large and orthogonal to everything else.
  `nearest_grid_index` is the reusable atom, already built.
- **A transient or observed CO₂ stream.** *Planned.* [#184](https://github.com/xiangtaoxu/MEDS/issues/184)

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
- **A packed `column_state_t`** ([#146](https://github.com/xiangtaoxu/MEDS/issues/146)). *Deferred to v0.3.0, **paired with #190***.
  Only a packed layout makes field omission a compile-time error; there are 1 207 field references
  today. One packed, policy-carrying layout should serve the fast state vector and the cohort slice
  together — separately, each is a large refactor buying a fraction of one property.

---

## 11. Performance

Source: `docs/dev_plans/MEDS_GPU_EVALUATION.md` §12.

- **Attack the allocator traffic.** *Planned.* [#195](https://github.com/xiangtaoxu/MEDS/issues/195) About 24 % of fast-loop self time is allocator
  work in `build_column_frozen`.
- **Thread and vectorise the cohort axis on the CPU.** *Planned.* [#196](https://github.com/xiangtaoxu/MEDS/issues/196) This is the evaluation's
  headline recommendation. Patch-axis threading already ships; the cohort axis is untouched.
- **A single-precision experiment** (`wp = real32`). *Candidate.* [#197](https://github.com/xiangtaoxu/MEDS/issues/197) Measured 56× on the device for
  FP32, which is why it is worth knowing what MEDS actually needs.

---

## 12. Smaller residues

*(Emptied by v0.2.0 — every item that was here shipped. New small residues go here as they are
filed.)*
