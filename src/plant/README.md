# src/plant — plant ecophysiology (`libmeds_plant`)

One flat, self-contained **plant-ecophysiology** library. It links `meds_shared` **only** — no cohort/
patch state (`site_t`) — so it compiles and unit-tests standalone
(`cmake --build … --target meds_plant`), orthogonal to the demographic core. Every process is a
stateless per-individual kernel driven by an environment struct; none is wired into the demographic
stepper yet. Design: [`docs/dev_plans/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md`](../../docs/dev_plans/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md).

Structural **allometry** lives in `src/allometry/` (a shared foundation below `state`, not here — see
issue #11), and the **empirical vital rates** (growth/mortality/recruitment) live in `src/demography/`
(they are empirical demographic functions, not ecophysiology).

## Contents

One **`meds_plant_interface`** hosts the public seams for the whole module (the single door for
callers); each domain's math lives in a dedicated **compute** module behind it.

| File | Role | Contents |
|------|------|----------|
| `meds_plant_types` | types | ALL derived types (leaf / hydraulics / phenology), one module, sectioned |
| `meds_plant_interface` | **the seams** | `leaf_gas_exchange(env, cfg, ipft, flux)` (flattens `cfg%pft`), `plant_water_flux(...)`, `update_phenology(...)` — thin wrappers; re-exports the public types |
| `meds_leaf_gas_exchange` | leaf compute | FvCB C3 + Collatz C4 demand, Leuning / Medlyn / Katul stomata, the bracketed Ci solver (`solve_leaf_gas_exchange`) |
| `meds_plant_hydraulics` | hydraulics compute | pressure-volume (Bartlett/Tyree-Hammel), Kirchhoff conductance, matrix-exp sub-step solver |
| `meds_phenology` | phenology compute | the cue engine → directional status (`phenology_kernel`) |
| `meds_plant_respiration` | respiration compute | non-leaf maintenance respiration: `stem_maintenance_respiration` + `fine_root_maintenance_respiration` + `growth_respiration`; re-exported via `meds_plant_interface` |
| `meds_plant_capi` | Python C-API | → `libmeds_plant_c` (`-DMEDS_BUILD_PYLIB=ON`; GLOB `*_capi.f90`); calls the leaf compute kernels directly, not the seam |

The shared temperature response (`meds_temp_response`, Arrhenius / peaked deactivation) lives in
`meds_shared` so leaf, respiration, and any tissue reach it without a plant→plant library edge.

## Python

The `*_capi.f90` shims are compiled only into the optional shared library `libmeds_plant_c`, exposed
through process-oriented Python packages (`python/meds/`): `meds.plant.leaf` (leaf gas exchange, reproduces
Slot & Winter 2017 in `examples/example_leaf_gas_exchange/`) and `meds.plant.pheno` (the leaf-phenology
kernel — the four phenology strategies in `examples/example_phenology/`).
