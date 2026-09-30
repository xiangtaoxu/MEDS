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
| [`MEDS_BCI_CENSUS_INIT_PLAN.md`](MEDS_BCI_CENSUS_INIT_PLAN.md) | Written 2026-09-28. The BCI flux-tower example starts from the 2010 census of the BCI 50-ha plot instead of a 50-year spin-up, with one PFT for every tree and MEDS's default allometry: a ForestGEO-to-MEDS census tool that only maps trees to patches (`patch_id`, `patch_area`, `nplant`), a header-matched census reader with `patch_area`, MEDS's own fusion at census initialization, `[init]` keys for the initial soil state, and a steady-state soil-carbon start. Builds on the allometry fix of #321. Design only; phases P0–P3 open. |
| [`MEDS_FAST_CALIBRATION_PLAN.md`](MEDS_FAST_CALIBRATION_PLAN.md) | Written and revised 2026-09-29. A fast, crude calibration of the fast-timescale parameters against flux-tower data, with the vegetation structure frozen (`[run].slow_on = false`) and trait plasticity on: 10-day windows restarted from shared states, the Bowen-ratio closure correction, Levenberg–Marquardt with parallel finite-difference Jacobians and priors, the Laplace covariance for a longer-term calibration, the model changes it needs first (hourly SW and LW up, trait re-acclimation at restart, provenance, physics fixes), and gates. BCI first, as part of the example. Design only; the order of the physics fixes is open (§12). |
| [`MEDS_FLUX_TOWER_FORCING_PLAN.md`](MEDS_FLUX_TOWER_FORCING_PLAN.md) | Written 2026-09-28. Forcing from flux-tower data, with Barro Colorado Island as the worked example: the forcing-format changes it needs (UTC only, the measured humidity in the file, the `ED_default`/`ED_ERA5land` names, the rain-timing fix), a shared Python module for the forcing file, the tower tool with its validation gates, longwave gap filling two ways, the BCI example and its tests. Phases P0–P5 open. |
| [`MEDS_POLYGON_RUNTIME_PLAN.md`](MEDS_POLYGON_RUNTIME_PLAN.md) | Written 2026-09-26, revised 2026-09-27. Runs of a contiguous region in one process, without MPI: month-synchronous compute/IO split, shared forcing month buffer, polygon-dimension output and its performance, ragged restart, batching and job arrays for large regions; site networks run as separate processes. Phases R0–R6: R0 measured (§10.1), R1, the compute/I-O split, and R2, region runs (#289), implemented (§10.2–§10.3); R3–R6 open (threads — the fast loop over all patches of all polygons, §10.4 — region restarts, failure isolation and tiles, interfaces). The non-MPI part of ROADMAP #183. |
| [`MEDS_PRODUCTION_INTEGRATOR_PLAN.md`](MEDS_PRODUCTION_INTEGRATOR_PLAN.md) | The numerics roadmap's remainder: adaptive freeze cadence (#158), soil water in the tableau (#159), the `rwc_floor` clamp artefact (#104), and two items without a live issue — `psi_leaf` at 900 s (#162, closed without a comment) and E5. Archive it once those two are re-homed. Also the record of what was refuted by measurement. |

## Reference — no open items, but cited by section from the code

| Document | Why it stays |
|---|---|
| [`MEDS_CODE_STRUCTURE_DESIGN.md`](MEDS_CODE_STRUCTURE_DESIGN.md) | The structure decisions (§1), the straddlers (§5), the placement rules (§6) and §7.6 #3–#4, cited by number from ~20 files; `src/README.md` carries only some of them. §15, the phased remainder, is closed. |
| [`MEDS_FROZEN_SEAM_CONTRACT.md`](MEDS_FROZEN_SEAM_CONTRACT.md) | The Λ criterion for when freezing a store is admissible, the four seams classified against it, the debit-before-credit rule for rate seams, and the arbitration rule for a shared store — stated nowhere else. Design note, no code (#201); no source comment cites it by section, so it moves to `archive/` once folded into the science pages. |
| [`MEDS_NUMERICS_SCOPING.md`](MEDS_NUMERICS_SCOPING.md) | §5.1 process mask, §11 bare-array convention, §12.6 ED2 `DTLSM` catalogue, §8b–§8g measurements — cited from ~20 files. Its scheme roadmap is superseded; see its header. |
| [`MEDS_IO_DESIGN.md`](MEDS_IO_DESIGN.md) | §3–§6 describe the registry, integrators and serializer as built. Cited by section from CMake and 9 source files. |
| [`MEDS_IO_V01_PLAN.md`](MEDS_IO_V01_PLAN.md) | The v0.1 diagnostic layer, cited by section from 12 source and test files. Its deferred list shipped in v0.2.0. |

## Reviews

| Document | |
|---|---|
| [`MEDS_CODE_REVIEW_2026-09-27.md`](MEDS_CODE_REVIEW_2026-09-27.md) | The review of `beta` before the v0.3.0 merge: the ranked pre-merge set (three verified must-fix items), the open-issue triage, the documentation drift, the build and test results, the consolidation items deferred to R3, the phased fix plan (§9) and this directory's archive audit (§10). |

## Untracked

`MEDS_CAPI_REORG_DESIGN.md` is a local working document, excluded in `.gitignore` by the author's
choice. It is referenced once from the migration log; its one load-bearing decision (a single
`libmeds.so`) is decision #1 of the structure plan and shipped in PR #138. Nothing in a clone
needs it.

---

## `archive/`

Forty documents whose work is done. Each opens with a tombstone: what shipped, in which pull
request, what changed name on the way in, and where the live description is. Three grades appear:

- **🗃️ ARCHIVED** — complete or superseded. Safe to read as history.
- **⚰️ RETIRED** — the conclusions are actively wrong. The tombstone says which sections and why.
  Four documents carry this: `MEDS_INTEGRATOR_PARITY.md` and `MEDS_INTEGRATOR_TEST.md` (both scored
  against the retired operator-split integrator, inside an oscillation no ledger detected),
  `MEDS_PHENOLOGY_DESIGN.md` (its tri-state algorithm no longer exists), and
  `MEDS_COLUMN_HYDROLOGY_DESIGN.md` (one of its "IMPLEMENTED" sections describes code that was
  later deleted). `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` is archived with ⚰️ on §9–§11 and §14,
  overturned the day they were written.
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
