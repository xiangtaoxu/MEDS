---
paths:
  - "src/io/**"
  - "src/config/meds_output_config.f90"
  - "test/test_output_*.f90"
  - "test/test_diagnostic_reduce.f90"
---

# The diagnostic output subsystem

User-facing description: [`docs/science/diagnostics.md`](../../docs/science/diagnostics.md). Design
records: `docs/dev_plans/MEDS_IO_DESIGN.md` and `MEDS_IO_V01_PLAN.md`, both cited by section number
from this code — do not renumber their sections.

## The five-stage wall

**derive → capture → reduce → integrate → serialize.**

1. **Derive** — pure closed-form quantities (leaf area index, conductances, water-use efficiency,
   soil potential, vapour-pressure deficit, size-class index). Calls the owning physics library
   rather than re-deriving.
2. **Capture** — the per-cohort and per-patch dt-weighted accumulators for everything the fast loop
   computes per sub-step and would otherwise discard.
3. **Reduce** — one weighted aggregation emitting the cohort → patch, site, PFT and size-class
   family.
4. **Integrate** — the per-variable, per-tier temporal folding.
5. **Serialize** — per-tier, per-time-chunk netCDF.

## The extensive/intensive contract is DATA, not code

A variable descriptor carries its weight, whether it is a mean, and its scale, declared once at
registration beside the units. An extensive per-plant quantity reduces as a weighted **sum** over
plant number; an intensive one as a weighted **mean**, and *which* weight is a physical statement —
basal area for diameter, leaf area for canopy temperatures and conductances, plant number for
demographic rates. This is the classic place a diagnostic goes silently wrong, so it is stated
rather than implied.

**Empty sets:** a mean over an empty patch, PFT or class is the fill value; a sum over one is a true
zero. Otherwise the identity that per-PFT and per-class sums equal the site total breaks the moment
a PFT goes locally extinct.

## Things that bite

- **The diagnostic blocks ride the cohort and patch lockstep**, and the obligation is discharged by
  **layout**: the fields are rows of one 2-D array, so each reorder or clear is a single whole-array
  statement that cannot omit a field. Adding a diagnostic costs zero edits to the reorder machinery.
- **Never read the transient tendency bundle from here.** Use the diagnostic block. The output tick
  runs after the monthly fuse-fission, so the two are different plants on boundary steps.
- **Source ids are range-partitioned by entity**, with a class dispatcher. This is load-bearing: the
  previous flat numbering had six literal collisions, invisible only because two switchboards
  consumed them separately — moving a variable between tiers would have silently written a
  different quantity.
- **Cohort and patch axes are forbidden on the annual stream**, because a window longer than a month
  straddles the disturbance restructuring.
- **The netCDF-free half is a separate CMake target from the serializer**, and it is an **explicit
  file list, not a glob**. Adding a diagnostic module there is a deliberate edit. That split is what
  keeps the stepper's edge free of a C dependency.
- **`[io]` is the restart stream only.** The legacy diagnostic writer under that block was retired
  at v0.1; it collided with the registry on the output filename prefix. A checkpoint is raw
  prognostic state at an instant, never a time average.

## Adding a variable

One `add_variable` line in the registry, with the weight, mean flag, scale, units and group. Then
check it appears in `meds_main --dump-io-config`. If it needs a new accumulator, add the row to the
diagnostic block — not a new standalone array.
