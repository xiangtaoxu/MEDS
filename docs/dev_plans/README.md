# Development plans

Dated design and implementation records: the account of *how* each MEDS subsystem was built, and
what is still intended. These are working documents, not user-facing reference. For the science,
read [`docs/science/`](../science/); for the code layout, [`src/README.md`](../../src/README.md).

## The rule for this directory

**`dev_plans/` holds a document only while it is still doing work** — it has open items, or a
section that source comments cite by number. Everything else is in
[`archive/`](archive/) with a tombstone at the top saying what shipped, in which pull request, and
where the live description now lives.

That is a change of policy, made 2026-09-13. The previous rule kept every document here on the
grounds that a point-in-time record is *supposed* to go stale, and reserved `archive/` for
documents whose conclusions were actively misleading. The result was 41 files in one flat
directory, three of which carried "design-only" headers for work that had merged the same day, and
no way to tell from a listing which ones still mattered.

**When a plan's last open item ships, give it a tombstone and move it in the same pull request.**
What changed goes in [`CHANGELOG.md`](../../CHANGELOG.md); what is deferred goes in
[`docs/ROADMAP.md`](../ROADMAP.md).

**Section numbers are never renumbered**, in either directory. Roughly 170 source comments cite
these documents by section (`§7.6 #3`, `§10.2.11`, `phase P2`), and most cite by bare filename, so
a move is safe but a renumber is not.

---

## Live — open items still intended

| Document | What is still open |
|---|---|
| [`MEDS_CODE_STRUCTURE_DESIGN.md`](MEDS_CODE_STRUCTURE_DESIGN.md) | The structure decisions, the placement rules, and §15 — the phased remainder: `column_cohort_t` removal, the packed state vector, per-layer face budgets, test-support consolidation. |
| [`MEDS_PRODUCTION_INTEGRATOR_PLAN.md`](MEDS_PRODUCTION_INTEGRATOR_PLAN.md) | The active numerics roadmap: adaptive freeze cadence, soil water in the tableau, the RK45 production warning, the `rwc_floor` clamp artefact. Also the record of what was refuted by measurement. |
| [`MEDS_BIOGEOCHEMISTRY_DESIGN.md`](MEDS_BIOGEOCHEMISTRY_DESIGN.md) | P1 nitrogen, P1 DAMM (kernel exists, unreachable), P2 vertically resolved pools, fire, coarse woody debris. |
| [`MEDS_FORCING_DESIGN.md`](MEDS_FORCING_DESIGN.md) | LWdown synthesis, the multi-polygon runtime, a transient CO₂ stream. Also the reference for the forcing NetCDF format (§7.1) and the ERA5-Land de-accumulation recipe (§7.3). |
| [`MEDS_SNOW_DESIGN.md`](MEDS_SNOW_DESIGN.md) | P1 multi-layer snow with compaction and an aging albedo; P2 canopy snow interception. |
| [`MEDS_ERA5_PREP_PLAN.md`](MEDS_ERA5_PREP_PLAN.md) | Design only (2026-09-25), no code yet. `scripts/prepare_era5/`: the single-grid scripts moved and renamed, plus a continental (North America first) download, then a yearly MEDS forcing store, then site extraction. Phases P0–P6; decisions D1–D7 open. |
| [`MEDS_GPU_EVALUATION.md`](MEDS_GPU_EVALUATION.md) | Five of seven recommendations, including cohort-axis threading and the allocator traffic. The measurement itself is closed: GPU offload is not viable as scoped. |

## Reference — no open items, but cited by section from the code

| Document | Why it stays |
|---|---|
| [`MEDS_FROZEN_SEAM_CONTRACT.md`](MEDS_FROZEN_SEAM_CONTRACT.md) | The Λ criterion for when freezing a store is admissible, the four seams classified against it, the debit-before-credit rule for rate seams, and the arbitration rule for a shared store. Design note, no code (#201). |
| [`MEDS_NUMERICS_SCOPING.md`](MEDS_NUMERICS_SCOPING.md) | §5.1 process mask, §11 bare-array convention, §12.6 ED2 `DTLSM` catalogue, §8b–§8g measurements — cited from ~20 files. Its scheme roadmap is superseded; see its header. |
| [`MEDS_IO_DESIGN.md`](MEDS_IO_DESIGN.md) | §3–§6 describe the registry, integrators and serializer as built. Cited by section from CMake and 9 source files. |
| [`MEDS_IO_V01_PLAN.md`](MEDS_IO_V01_PLAN.md) | The v0.1 diagnostic layer, cited by section from 12 source and test files. Carries its own deferred list. |
| [`MEDS_VEG_ENERGY_INTEGRATION_PLAN.md`](MEDS_VEG_ENERGY_INTEGRATION_PLAN.md) | §1–§8 and §12–§13 are the exact-exponential tissue store. **§9–§11 and §14 are overturned** — read its correction header first. |

## Reviews

| Document | |
|---|---|
| [`MEDS_DOCS_REVIEW_2026-09-13.md`](MEDS_DOCS_REVIEW_2026-09-13.md) | The documentation audit this reorganization executed: per-file verdicts on all 41 plans, and the plans for the README, `src/README.md`, `CLAUDE.md` and the changelog. |

## Untracked

`MEDS_CAPI_REORG_DESIGN.md` is a local working document, excluded in `.gitignore` by the author's
choice. It is referenced once from the migration log; its one load-bearing decision (a single
`libmeds.so`) is decision #1 of the structure plan and shipped in PR #138. Nothing in a clone
needs it.

---

## `archive/`

Thirty-three documents whose work is done. Each opens with a tombstone: what shipped, in which pull
request, what changed name on the way in, and where the live description is. Three grades appear:

- **🗃️ ARCHIVED** — complete or superseded. Safe to read as history.
- **⚰️ RETIRED** — the conclusions are actively wrong. The tombstone says which sections and why.
  Four documents carry this: `MEDS_INTEGRATOR_PARITY.md` and `MEDS_INTEGRATOR_TEST.md` (both scored
  against the retired operator-split integrator, inside an oscillation no ledger detected),
  `MEDS_PHENOLOGY_DESIGN.md` (its tri-state algorithm no longer exists), and
  `MEDS_COLUMN_HYDROLOGY_DESIGN.md` (one of its "IMPLEMENTED" sections describes code that was
  later deleted).
- **Status header is wrong** — three documents say "design-only" for work that merged within a day
  of being written: `MEDS_REORG_DESIGN.md`, `MEDS_DRIVER_REORG_DESIGN.md` and
  `MEDS_CORE_MODULE_REORG_DESIGN.md`. Their tombstones say so first, because anyone reading the
  header alone gets the history backwards.

Two files in `archive/` were split out of the structure plan rather than being plans in their own
right: `MEDS_CODE_STRUCTURE_MIGRATION_LOG.md` (how the 2026-09 reorganization ran) and
`MEDS_SLOW_LOOP_CONSERVATION_LEDGER.md` (the measurement record behind PRs #132–#137, formerly
§10.2).

**Line numbers in archived documents are the pre-2026-09 tree.** Two reorganizations have moved
almost every file since most of them were written.
