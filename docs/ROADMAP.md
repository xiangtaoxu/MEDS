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
| [#265](https://github.com/xiangtaoxu/MEDS/issues/265) | Sub-canopy conductance: MEDS ports ED2's non-default `icanturb = 4`, giving an 8–16× too-stiff ground resistance | Matters in gaps |
| [#268](https://github.com/xiangtaoxu/MEDS/issues/268) | No litter layer: no surface organic horizon for the ground energy balance or soil evaporation to act on | |
| [#269](https://github.com/xiangtaoxu/MEDS/issues/269) | The canopy air space is one well-mixed slab: 1.2 K warmer than the free air at midday, where a real sub-canopy is cooler and steadier | |

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

Source: `docs/dev_plans/archive/MEDS_BIOGEOCHEMISTRY_DESIGN.md` §7. Science page:
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

- **N5 — an adaptive freeze cadence with a real error estimator.** *Closed 2026-09-30; revisit with
  the numerical scheme.* [#158](https://github.com/xiangtaoxu/MEDS/issues/158) It was proposed when the
  stability limit forced `dt_fast = 150 s` for the whole day. N2a removed that limit and the shipped
  `dt_fast` is 900 s, so what is left is a bounded efficiency gain: the coefficient pre-pass is
  4.4–11 % of a step. Any cadence would have to stay site-wide, because the threaded patch loop shares
  the sub-step's forcing samples and output staging.
- **Fold soil water into the ARK tableau.** *Closed 2026-09-30; revisit with the numerical scheme.*
  [#159](https://github.com/xiangtaoxu/MEDS/issues/159) — successor to #93. In-stage soil water costs
  +14–26 % and removes a split error of 0.011–0.018 % of column water. On ARK the committed soil state
  and the fluxes its soil-heat stages use come from one solve, so the borrowed-flux defect class is
  absent there. The drought result behind closing #93 (2.3 % at 900 s) was measured under
  `ARREST_GS_CLAMP`, which #335 replaced; re-measure it in that round.
- **The `rwc_floor` clamp artefact** ([#104](https://github.com/xiangtaoxu/MEDS/issues/104)). *Planned* (decided 2026-09-30,
  `dev_plans/MEDS_EFFICIENCY_SWEEP_PLAN.md` Phase 3 and Appendix A): the wood's apoplastic water drains
  as its conduits embolise, as TFS-Hydro and SurEau treat it, so every wood mass has a finite
  potential. A floored relative water content maps today to a potential of about −10⁴ MPa, which is
  not a pressure any tissue reaches. Note that *arresting* is not the free
  option it looks: the collapsed store diagnoses ψ at about −1.5×10⁴ MPa against a soil at perhaps
  −2 MPa, so the cohort recovers today — that enormous artificial gradient IS the 13× cost — and
  removing uptake would make a transiently desiccated cohort permanently dead.
- **`psi_wood` at 900 s, and a restart's first day.** *Known limitation; revisit with the numerical
  scheme.* [#162](https://github.com/xiangtaoxu/MEDS/issues/162) `psi_leaf` converges at 900 s on ARK since the
  transpiration corrector (#91): re-measured 2026-09-30 over July on the established Ithaca stand,
  daily means at 900 s match a 12.5 s run to 0.001 MPa on every day after the first
  (`science/numerical_scheme.md` §5a). The first day after a restart still carries a start-up
  transient (−1.83 against −0.28 MPa), and `psi_wood` keeps 0.17 MPa at 900 s on the 3-hour midday
  probe through the frozen uptake seam.

---

## 5. Vegetation energy

Source: `docs/dev_plans/archive/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §6–§7. Science page:
[`science/vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md).

- **A separate canopy film store with phase change.** *Planned.* [#165](https://github.com/xiangtaoxu/MEDS/issues/165) Intercepted water currently
  has no independent thermal state and cannot freeze.
- **The free-convection slope.** *Closed 2026-09-30; revisit with the numerical scheme.*
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
- **A patch axis on the FAST tier.** *Candidate.* [#270](https://github.com/xiangtaoxu/MEDS/issues/270) Sub-daily output is site-mean only,
  although the fast staging already carries the patch dimension.
- **Within-step variances and a skin temperature.** *Candidate.* [#275](https://github.com/xiangtaoxu/MEDS/issues/275) v0.3.0 dropped
  `ground_temp_site` and relabelled the four variances as what they are, variances of end-of-step
  samples. What remains: a sum-of-squares row per variance in the patch block (`PD_*_SQ`) with an
  aggregation that consumes a mean and a mean square, and the skin temperature the ground balance
  already computes.
- **Fast-only rows in a slow-only run read 0.** *Candidate.* [#299](https://github.com/xiangtaoxu/MEDS/issues/299) v0.3.0 weights the patch block
  with the slow step when the fast loop is off, so the slow rows report their rates. The rows only
  the fast loop fills (fluxes, the forcing echo) still read 0 there, where `_FillValue` would say
  "not simulated"; that needs a per-variable flag in the registry.
- **One serializer for site and region files.** *Candidate.* [#312](https://github.com/xiangtaoxu/MEDS/issues/312) The region writer duplicates
  the site writer, and the forcing echo is enumerated in three places.

---

## 7. Plant physiology

- **PER-COHORT thermal acclimation.** *Candidate.* [#256](https://github.com/xiangtaoxu/MEDS/issues/256) Follow-up to #176:
  Photosynthetic acclimation shipped in v0.2.0 as a **site-level** growth temperature (Kattge &
  Knorr 2007), which is what that fit is calibrated on — the mean air temperature of the preceding
  weeks. Letting a shaded understory cohort acclimate differently from a sunlit canopy one needs a
  relation calibrated on **leaf** rather than air temperature; applying the air-temperature fit to a
  per-cohort leaf temperature would use it well outside its range. Acclimation of **respiration**
  (as opposed to photosynthetic capacity) is also still open.
- **Whether storage maintenance should be ON by default.** *Open question.* The mechanism
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

Source: `docs/dev_plans/archive/MEDS_FORCING_DESIGN.md` §5.7, §8. Science page:
[`science/forcing.md`](science/forcing.md).

- **A better cloud term for the LWdown synthesis.** *Candidate.* [#257](https://github.com/xiangtaoxu/MEDS/issues/257)
  Follow-up to #182: The Brutsaert/Idso synthesis shipped in
  v0.2.0, so a source lacking longwave can now drive MEDS. Its cloud correction is one empirical
  coefficient on the SW clearness index, and it holds the last daytime index through the night —
  adequate for a fallback, but the residual is real: driving Ithaca from synthesis leaves the soil
  surface 1.37 K cooler than the file's `strd`. A cloud-fraction formulation (Crawford & Duchon
  1999) or a nocturnal index carried from a longer window would close more of it.
- **Forcing adapters for NLDAS-3, Daymet and CHIRPS** (forcing phase F6). *Candidate.*
  [#302](https://github.com/xiangtaoxu/MEDS/issues/302) Deferred 2026-09-27; ERA5-Land stays the
  only product. Each would follow the ED_ERA5land pattern: a downloader, the same flat monthly
  archive, and a `met_source` entry (`MEDS_FORCING_DESIGN.md` §16). NLDAS-3 (~1 km, hourly) could
  drive MEDS alone but needs regional archives and its own chunk size (OD3); Daymet and CHIRPS are
  daily corrections to an hourly base.
- **ED_ERA5land archive years before June 2002.** *Planned.*
  [#303](https://github.com/xiangtaoxu/MEDS/issues/303) Deferred until a study needs them. GDEX
  starts in June 2002, so they come from the CDS, with the tools that already built 2026-04 and
  2026-06.
- **The multi-polygon runtime.** *Planned.* [#183](https://github.com/xiangtaoxu/MEDS/issues/183) Region runs exist since R2 (#289):
  `[run].mode = "region"` runs every selected ED_ERA5land cell of a `[region]` box as its own polygon, in
  one process, all sharing one forcing reader. What remains is in `MEDS_POLYGON_RUNTIME_PLAN.md` §10:
  - R3, threads: the fast loop over all patches of all polygons (plan §10.4; today the polygons
    are stepped one after another on one thread), with the runtime consolidation of [#310](https://github.com/xiangtaoxu/MEDS/issues/310);
  - R4, region checkpoints and restarts (a region writes none yet);
  - R5, failure isolation and batching by tiles;
  - R6, a C API and Python entry point, and an example.

  MPI is not planned: a large region runs as tiles in a job array (§8 of that plan).
- **A source interface for the met reader.** *Candidate.* [#311](https://github.com/xiangtaoxu/MEDS/issues/311) `meds_met_driver` carries three
  sources behind ten backend branches; the hourly cadence and the time-units parsing are repeated.

---

## 9. Snow

Source: `docs/dev_plans/archive/MEDS_SNOW_DESIGN.md` §7. Science page:
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
  per-field policy table. *Deferred, **paired with #146*** (2026-09-13).
  [#190](https://github.com/xiangtaoxu/MEDS/issues/190) 38 references across 11 files. Four of the
  five benefits the design claimed have since landed piecemeal: `column_cohort_init` gives the test
  fixtures allometric consistency, the three hard-coded constants are PFT parameters, the derived
  geometry is on the cohort block, and `reconcile_tissue_water_capacity` took the seed and clamp out
  of the gather. The fusion/scaling policy is already centralised in `fuse_cohort_fast_state` and
  `scale_cohort_ground_fields`. What is left is **completeness you cannot forget** — a table the
  blend iterates cannot omit a field a hand-written routine can — and that is #146's hazard class,
  which is why the two now travel together.
- **A packed `column_state_t`** ([#146](https://github.com/xiangtaoxu/MEDS/issues/146)). *Deferred, **paired with #190***.
  Only a packed layout makes field omission a compile-time error; there are 1 207 field references
  today. One packed, policy-carrying layout should serve the fast state vector and the cohort slice
  together — separately, each is a large refactor buying a fraction of one property.

---

## 11. Performance

Source: `docs/dev_plans/archive/MEDS_GPU_EVALUATION.md` §12.

- **Attack the allocator traffic.** *Planned.* [#195](https://github.com/xiangtaoxu/MEDS/issues/195) About 24 % of fast-loop self time is allocator
  work in `build_column_frozen`.
- **Thread and vectorise the cohort axis on the CPU.** *Planned.* [#196](https://github.com/xiangtaoxu/MEDS/issues/196) This is the evaluation's
  headline recommendation. Patch-axis threading already ships; the cohort axis is untouched.
- **A single-precision experiment** (`wp = real32`). *Decided: no.* [#197](https://github.com/xiangtaoxu/MEDS/issues/197) Measured 56× on the device for
  FP32, which is why it was worth asking what MEDS needs; the answer (2026-09-13) is that MEDS stays
  `real64` throughout. Kept here as the record; the issue closes with v0.3.0.

---

## 12. Smaller residues

*(Emptied by v0.2.0 — every item that was here shipped. New small residues go here as they are
filed.)*
