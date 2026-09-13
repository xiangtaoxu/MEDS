---
paths:
  - "src/state/**"
  - "src/slow_dynamics/**"
  - "src/init/**"
  - "test/test_containers.f90"
  - "test/test_fusion_cohort.f90"
  - "test/test_patch.f90"
  - "test/test_disturbance.f90"
  - "test/test_carbon_growth.f90"
  - "test/test_slow_ledger.f90"
  - "test/test_state_*.f90"
---

# State and demography

See [`src/README.md`](../../src/README.md) for the layering.

## The state layout

**A flat, site-wide structure of arrays.** Every cohort of the whole site lives in one contiguous
set of 1-D arrays, with patch membership as a CSR map (offset plus count, and an owner index). The
dominant daily kernels are a single unit-stride sweep.

**One centralized lockstep reorder** in `state/site/meds_site_state_types.f90`: the reorder,
compaction, slot copy, CSR rebuild and capacity routines, plus the size setter that fills one slot's
cached geometry. **When you add a per-cohort field, update these** — the single place that touches
every array. This is the structural fix for ED2's "forgot to reallocate one array" bug class.

Patch arrays have no single reorder routine. Their permute and pack sites are `sort_patches` and
`patch_compact` in `meds_demography_patch_fusefiss`; update both.

**Persistent identity.** Every cohort and patch gets a monotonic `global_id` at creation and carries
it in lockstep through every sort, fusion and compaction. Ids are never reused; fusion keeps the
survivor's; a split daughter, a recruit and a disturbance fragment each get a fresh one. Creation
sites that must stamp: the initial-community builders, recruitment, cohort split, and patch
disturbance.

**Cached geometry is per plant.** Leaf area, wood area, sapwood carbon and sapwood area are cached
on the cohort block; per-ground indices are formed at the point of use, so no stored value carries a
plant density that mortality can invalidate. Re-derive wherever size changes. A stale cache is
invisible in a conservation ledger — that is what `test_carbon_growth` exists to catch.

## Laws and operators do not touch

`demography/` holds both, and they stay apart: the operators take rate **arrays** as arguments and
never import `meds_demography_rates`. The slow driver is the one place a rate meets its application.
The Python `apply_rates` path, feeding externally computed rates through the same operators, is the
standing test.

**The engine never computes a rate.** The driver computes a per-cohort tendency bundle and the
engine applies it. The bundle is transient — deliberately not lockstep-reordered — which is correct
for its own consumer and wrong for anything that runs after a fuse-fission.

## Conservation

Every fuse and split asserts within tolerance or stops the run: cohort fusion and splitting conserve
total aboveground biomass and plant number, with diameter **re-derived** from the conserved per-plant
biomass rather than averaged; patch fusion conserves site-level plant number through area rescaling,
and patch area renormalizes to one. Patch **disturbance** conserves area but intentionally does not
conserve plant number — the killed canopy is the disturbance.

The slow loop carries a site-wide ledger on top of the per-store fast budgets. It found real leaks
that every fast budget had passed, including growth respiration that was 23 % of growth carbon and
was never exhaled. **A closed budget proves bookkeeping, not plausibility.**

## Order of operations

Every step: growth, then mortality, then patch ageing. Monthly: recruit, then cohort
fuse/terminate/split, then sort. Annually: disturbance, then patch restructuring. Disturbance and
patch restructuring are **independent triggers**. Disturbance integrates its yearly rate over the
one-year patch-dynamics interval, not over the per-step interval.

## Output must not read the tendency bundle

The output tick runs after the monthly fuse-fission, so the bundle's index `i` and the cohort array's
index `i` are different plants on exactly the boundary steps. Read the diagnostic block instead.
This mistake was made once and caught only by the thread-invariance test, because thread count
perturbs which cohorts fuse.
