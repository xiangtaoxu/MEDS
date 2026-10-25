# Development plans

Design and implementation records: how each MEDS subsystem was built, and what is still intended.
They are working documents, not user-facing reference. For the science, read
[`docs/science/`](../science/); for the code layout, [`src/README.md`](../../src/README.md).

## The rule for this directory

A plan stays in `dev_plans/` only while it has open items, or a section that source comments cite by
number. Every other plan is in [`archive/`](archive/), with a tombstone at the top that says what
shipped, in which pull requests, and where the live description is.

When a plan's last open item ships, give it a tombstone and move it in the same pull request. What
changed goes in [`CHANGELOG.md`](../../CHANGELOG.md); what is deferred goes in
[`docs/ROADMAP.md`](../ROADMAP.md).

Section numbers are never renumbered, in either directory. Code cites a plan by section only where
the plan is the one written record of a design rationale; most citations use the bare file name, so
a move is safe and a renumbering is not.

---

## Live: open items

| Document | What is open |
|---|---|
| [`MEDS_BCI_CENSUS_INIT_PLAN.md`](MEDS_BCI_CENSUS_INIT_PLAN.md) | The census start of the BCI example is built. Open: the census-plot figure and the validation gates of §8 that need model runs. |
| [`MEDS_EFFICIENCY_SWEEP_PLAN.md`](MEDS_EFFICIENCY_SWEEP_PLAN.md) | The sweep's phases 0–3 shipped. Open: the runtime consolidation and region threading of Phase 5 (#310, #183, #400), the cohort axis (#196), and Phase 4, moved to a later numerics round. |
| [`MEDS_MULTI_SITE_PLAN.md`](MEDS_MULTI_SITE_PLAN.md) | Several sites per polygon, steps S0–S5, is built. Open: the diffuse light of a slope, a site map for a region, lateral water and seed dispersal between sites, and the example (#455). |
| [`MEDS_PHENOLOGY_SENESCENCE_PLAN.md`](MEDS_PHENOLOGY_SENESCENCE_PLAN.md) | The generic cue model is built. Open: senescence as an up-regulation of the PFT's baseline leaf turnover (#379). |
| [`MEDS_POLYGON_RUNTIME_PLAN.md`](MEDS_POLYGON_RUNTIME_PLAN.md) | Region runs are built (§14). Open: patch threads inside a region's polygons (#400), one run container and one loop (#310), reading the next month's forcing during the current step (#407). |
| [`MEDS_PRODUCTION_INTEGRATOR_PLAN.md`](MEDS_PRODUCTION_INTEGRATOR_PLAN.md) | The numerics roadmap's remainder and the record of what measurement refuted. Open: the `dt_fast` checks (#415) and the closed items to revisit with the numerical scheme (#158, #159). |
| [`MEDS_STATE_HIERARCHY_PLAN.md`](MEDS_STATE_HIERARCHY_PLAN.md) | The state layer by model hierarchy is built. Open: the follow-ups of §9 (the cohort file's field lists, folding the pool changes into the tendency bundle, the output-only accumulators on the polygon state). |

## Reference: no open items, cited from code or rules

| Document | Why it stays |
|---|---|
| [`MEDS_CODE_STRUCTURE_DESIGN.md`](MEDS_CODE_STRUCTURE_DESIGN.md) | The structure decisions behind `src/`'s layout; `src/c_api/meds_c_api_demography.f90` cites §7.6. Its phased remainder is closed, and [`src/README.md`](../../src/README.md) is the live description of the layout. |
| [`MEDS_EXAMPLE06_DRY_FOREST_PHENOLOGY_PLAN.md`](MEDS_EXAMPLE06_DRY_FOREST_PHENOLOGY_PLAN.md) | The PFT and phenology design of example 06; example 05's Palo Verde PFT file cites §6. |
| [`MEDS_FAST_CALIBRATION_BEST_PRACTICE.md`](MEDS_FAST_CALIBRATION_BEST_PRACTICE.md) | The calibration protocol, cited by section from `scripts/calibrate_fast/`. |
| [`MEDS_FAST_CALIBRATION_PLAN.md`](MEDS_FAST_CALIBRATION_PLAN.md) | The method and the argument for each fitted key's range, cited by `scripts/calibrate_fast/parameters.toml` (§4). |
| [`MEDS_FLUX_TOWER_FORCING_PLAN.md`](MEDS_FLUX_TOWER_FORCING_PLAN.md) | How the tower forcing conversions were found and validated; example 04's site file cites §2–§3. |
| [`MEDS_FROZEN_SEAM_CONTRACT.md`](MEDS_FROZEN_SEAM_CONTRACT.md) | The Λ criterion for freezing a store, the debit-before-credit rule and the arbitration rule for a shared store, written nowhere else. Moves to `archive/` once folded into the science pages. |
| [`MEDS_IO_DESIGN.md`](MEDS_IO_DESIGN.md) | The output registry, integrators and serializer as built; cited by section from the output modules and `.claude/rules/output.md`. |
| [`MEDS_IO_V01_PLAN.md`](MEDS_IO_V01_PLAN.md) | The diagnostic layer, cited by `.claude/rules/output.md` and `docs/science/diagnostics.md`. |

## Untracked

`MEDS_CAPI_REORG_DESIGN.md` is a local working document, excluded in `.gitignore` by the author. Its one
load-bearing decision (a single `libmeds.so`) is recorded in the archived structure plan. Nothing in a
clone needs it.

---

## `archive/`

Forty-seven documents whose work is done. Each opens with a tombstone: what shipped, in which pull
request, what changed name on the way in, and where the live description is. Three grades appear:

- **🗃️ ARCHIVED**: complete or superseded. Safe to read as history.
- **⚰️ RETIRED**: the conclusions are wrong. The tombstone says which sections and why. Four
  documents carry this: `MEDS_INTEGRATOR_PARITY.md` and `MEDS_INTEGRATOR_TEST.md` (scored against the
  retired operator-split integrator), `MEDS_PHENOLOGY_DESIGN.md` (its tri-state algorithm no longer
  exists) and `MEDS_COLUMN_HYDROLOGY_DESIGN.md` (one section describes code that was later deleted).
  `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` carries ⚰️ on §9–§11 and §14.
- **Status header is wrong**: `MEDS_REORG_DESIGN.md`, `MEDS_DRIVER_REORG_DESIGN.md` and
  `MEDS_CORE_MODULE_REORG_DESIGN.md` say "design-only" for work that merged within a day of being
  written. Their tombstones say so first.

`MEDS_CODE_STRUCTURE_MIGRATION_LOG.md` and `MEDS_SLOW_LOOP_CONSERVATION_LEDGER.md` were split out of
the structure plan ([`MEDS_CODE_STRUCTURE_DESIGN.md`](MEDS_CODE_STRUCTURE_DESIGN.md)) rather than being plans of
their own.

Line numbers and module paths in archived documents are from the tree at the time they were written;
the 2026-09 reorganizations moved almost every file since.
