# MEDS documentation

**Coming from ED2?** Start with **[`ed2_comparison.md`](ed2_comparison.md)** — what carried over,
what changed, what is missing, and what a MEDS run costs to set up, keyed to the `ED2IN` namelist
options you already know.

**Want to change the code?** Start with [`../src/README.md`](../src/README.md) — the source layout,
the placement rules, and the library graph.

## What is here

- **[`science/`](science/)** — the model's **equations and their rationale**: reader-facing theory
  that tracks the code. Math renders natively on GitHub.
- **[`dev_plans/`](dev_plans/)** — dated **design and implementation records**, and
  [`dev_plans/archive/`](dev_plans/archive/) for the ones whose work is done. Working documents,
  not reference. [`dev_plans/README.md`](dev_plans/README.md) says which is which and why.
- **[`ROADMAP.md`](ROADMAP.md)** — everything MEDS has deliberately deferred, in one place.
- **[`configuration.md`](configuration.md)** — the two TOML files, the parameter philosophy, and
  the three ways a run can be initialized.
- **[`building.md`](building.md)** — compilers, dependencies, build types, and the parallel builds.
- **[`ed2_comparison.md`](ed2_comparison.md)** — MEDS against ED2, process by process.

The changelog is at the repository root: [`../CHANGELOG.md`](../CHANGELOG.md).

## The science pages

**Sub-daily, the fast loop**

| Page | Covers |
|---|---|
| [`column_biophysics.md`](science/column_biophysics.md) | The integrative hub: how the surface stores are coupled and advanced together. |
| [`numerical_scheme.md`](science/numerical_scheme.md) | How the fast loop is advanced in time: what is frozen each `dt_fast`, the two integrators, what differs between them, and which to use. |
| [`canopy_radiation_transfer.md`](science/canopy_radiation_transfer.md) | The ED2 two-stream solver, leaf-angle distributions, and the multi-band adding solution. |
| [`canopy_aerodynamics.md`](science/canopy_aerodynamics.md) | Monin-Obukhov surface layer, leaf and wood boundary layers, in-canopy wind. |
| [`canopy_air_space_biophysics.md`](science/canopy_air_space_biophysics.md) | The three prognostic canopy-air twins: enthalpy, humidity, CO₂. |
| [`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) | Photosynthetic demand, three stomatal-conductance models, two-limb water stress, the coupled Cᵢ solver. |
| [`plant_hydraulics.md`](science/plant_hydraulics.md) | The plant water-transport ODE: pressure-volume and Kirchhoff curves, xylem vulnerability, the matrix-exponential solver. |
| [`vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md) | Leaf and wood temperature as an exact exponential relaxation. |
| [`soil_biophysics.md`](science/soil_biophysics.md) | Soil water (implicit Richards) and soil heat, both as conserved quantities. |
| [`snow_biophysics.md`](science/snow_biophysics.md) | The snow and temporary-surface-water store. |
| [`forcing.md`](science/forcing.md) | Meteorological forcing: the file format, temporal interpolation, the interval-mean-conserving shortwave disaggregation, calendar recycling. |

**Daily to annual, the slow loop**

| Page | Covers |
|---|---|
| [`plant_phenology.md`](science/plant_phenology.md) | The leaf-phenology signal kernel: cues → two governor drives → flush and shed rates. |
| [`plant_carbon_allocation.md`](science/plant_carbon_allocation.md) | The daily PARTEH-H1 carbon budget: net carbon → priority ladder → per-pool growth. |
| [`plant_respiration.md`](science/plant_respiration.md) | The autotrophic budget: leaf, stem and fine-root maintenance, and growth respiration on realized growth. |
| [`plant_traits.md`](science/plant_traits.md) | Light-driven leaf-trait plasticity and its carbon consequences. |
| [`soil_carbon.md`](science/soil_carbon.md) | The CENTURY decomposition network as a carbon matrix ODE, and the fast/slow respiration seam. |

**Output**

| Page | Covers |
|---|---|
| [`diagnostics.md`](science/diagnostics.md) | The diagnostic subsystem: what is emitted, on which axes and timescales, and the extensive/intensive contract. |

These pages are versioned with the source, so an equation change and its documentation land in the
same pull request. For the *what-each-routine-is* API level, read the doc comments in the Fortran
source directly.

## Writing a science page

Math renders on GitHub, whose Markdown sanitizer runs *inside* `$…$` before MathJax — so display
equations go in ` ```math ` fenced blocks, never `$$…$$`, and inline math containing backslash
spacing or paired brace subscripts needs the dollar-backtick form. The full rule set is in
`.claude/rules/docs-math.md`.

There is **no local preview**. Validate structurally (no `$$` left, fences balanced, every
conflicting construct protected) and eyeball the rendered file on a branch before merging.
