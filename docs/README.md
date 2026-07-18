# MEDS documentation

Two kinds of docs live here:

- **[`science/`](science/)** — the model's **equations and their rationale** (reader-facing theory
  that tracks the code). Math renders natively on GitHub (LaTeX in `$…$` and ` ```math ` fences).
  - [`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) — photosynthetic demand, the three
    stomatal-conductance models (Leuning / Medlyn / Katul), two-limb water stress, and the coupled
    Cᵢ solver.
- **[`dev_plans/`](dev_plans/)** — dated **design & implementation plans** and code reviews: the
  historical record of *how* each subsystem was built (formerly the top-level `archive/`). These are
  working documents, not user-facing reference.

These pages are versioned with the source, so an equation change and its documentation land in the
same pull request. For the *what-each-routine-is* API level, read the doc-comments in the Fortran
source directly.

### Planned pages
`science/photosynthesis.md` (FvCB C3 / Collatz C4 demand), `science/plant_hydraulics.md`,
`science/canopy_radiation.md`, `science/soil_carbon.md` — added as each subsystem is documented.
