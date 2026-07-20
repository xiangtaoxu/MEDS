# MEDS documentation

Two kinds of docs live here:

- **[`science/`](science/)** — the model's **equations and their rationale** (reader-facing theory
  that tracks the code). Math renders natively on GitHub (LaTeX in `$…$` and ` ```math ` fences).
  - [`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) — photosynthetic demand, the three
    stomatal-conductance models (Leuning / Medlyn / Katul), two-limb water stress, and the coupled
    Cᵢ solver.
  - [`plant_hydraulics.md`](science/plant_hydraulics.md) — the plant water-transport ODE (PV /
    Kirchhoff curves, xylem vulnerability, the matrix-exponential solver).
  - [`plant_phenology.md`](science/plant_phenology.md) — the leaf-phenology signal kernel: cues →
    two governor drives → flush/shed rates (with the baseline-turnover floor).
  - [`plant_carbon_allocation.md`](science/plant_carbon_allocation.md) — the daily PARTEH-H1 carbon
    budget: net carbon → priority ladder → per-pool growth, with growth respiration on realized growth.
  - [`plant_traits.md`](science/plant_traits.md) — light-driven leaf-trait plasticity (SLA / Vcmax /
    Rd / leaf lifespan) and its carbon consequences.
- **[`dev_plans/`](dev_plans/)** — dated **design & implementation plans** and code reviews: the
  historical record of *how* each subsystem was built (formerly the top-level `archive/`). These are
  working documents, not user-facing reference.

These pages are versioned with the source, so an equation change and its documentation land in the
same pull request. For the *what-each-routine-is* API level, read the doc-comments in the Fortran
source directly.

### Planned pages
`science/canopy_radiation.md`, `science/soil_carbon.md`, `science/energy_water_budgets.md` — added as
each subsystem is documented.
