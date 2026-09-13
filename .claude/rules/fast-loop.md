---
paths:
  - "src/fast_dynamics/**"
  - "test/test_column_*.f90"
  - "test/test_fast_loop.f90"
  - "test/test_surface_energy.f90"
  - "test/test_snow.f90"
  - "test/test_aerodynamics.f90"
  - "test/test_canopy_radiation.f90"
---

# The fast loop

Sub-daily biophysics on the `dt_fast` tier. See [`src/fast_dynamics/README.md`](../../src/fast_dynamics/README.md)
for the layout and [`docs/science/numerical_scheme.md`](../../docs/science/numerical_scheme.md) for
the integration.

## The freeze, and what it costs

Each `dt_fast` step builds a **frozen work record** and marches against it. Leaf gas exchange and
leaf water potential are frozen across the step, so daily photosynthesis drifts with `dt_fast` when
the non-stomatal water-stress limb is on.

**`dt_fast` is an ACCURACY parameter, not a stability one.** It buys sub-daily detail. It used to be
stability-bound — a period-2 canopy-air oscillation above ~150–225 s that no conservation ledger
detected — and that was traced to a single frozen coefficient, the canopy-air-to-atmosphere
conductance, now refreshed at every integrator stage. Shipped configs use **900 s**; the config
warns above 225 s only when the non-stomatal limb is on.

## Two integrators, one dispatch

- **`ark`** (default) — an L-stable ESDIRK2 on the ARS(2,2,2) tableau, with a 2×2 Newton arrowhead
  for the leaf ↔ canopy-air surface solve. Despite the historical "IMEX-ARK" name the explicit part
  is empty, so it is diagonally implicit.
- **`rk45`** — adaptive Cash-Karp 5(4), the accuracy baseline.

`meds_fast_step` dispatches and owns the **RK45-to-ARK stiff rescue**. There is no operator-split
scheme; it was retired for converging to a different limit that refinement never removed.

**The transpiration corrector is ARK-only.** RK45 at production cadence therefore carries a leaf
water-potential error the default path does not.

## Rules that break things quietly

- **Kernels never see `site_t`.** `canopy/`, `plant/` and `soil/` link `state_column` and `config`
  only. That is what keeps them device-eligible and standalone-buildable. A routine that needs
  `site_t` is driver code.
- **Build a column cohort view through `meds_column_view`.** `copy_column_cohort` for production,
  `column_cohort_init` for tests and probes — the latter builds the block through the canonical
  birth path, so a fixture tree is on-allometry by construction. Hand-assembling one is how three
  column tests ended up running on trees that could not exist.
- **Declare a fast-loop field's fusion kind once**, in `fuse_cohort_fast_state`: intensive fields
  blend leaf-area-weighted, extensive ones by plant number, ground-referenced ones sum. Getting it
  wrong at a call site is invisible — the biomass assertion passes and both ledgers close.
- **Every store closes a budget every step**, asserted in the suite. But **a closed budget proves
  bookkeeping, not plausibility**: several real leaks sat behind ledgers that balanced, because the
  ledger declared the same quantity it consumed. Watch consecutive-step traces, not just residuals.
- **Do not borrow the scratch solve's flux while committing your own state.** That defect class
  produces a vertical-only error that whole-column ledgers are blind to.

## Threading

The patch axis is threaded, and output is **byte-identical at any thread count** — the suite asserts
it. Per-patch scratch is an explicit per-thread pool indexed by thread number, not an OpenMP
data-sharing clause, because ifx builds `private` copies of a derived type through a
compiler-generated static mold that every thread writes.

## The one fast/slow seam

`heterotrophic_respiration_matrix`, called from the pre-pass so the sub-daily respiration debits the
same CENTURY pool the daily step does. The residual is asserted and closes to machine precision.
Do not add a second cross-tier kernel call without the same treatment.
