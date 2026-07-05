# src/plant — plant ecophysiology (`libmeds_plant`)

One flat, self-contained **plant-ecophysiology** library. It links `meds_shared` **only** — no cohort/
patch state (`site_t`) — so it compiles and unit-tests standalone
(`cmake --build … --target meds_plant`), orthogonal to the demographic core. Every process is a
stateless per-individual kernel driven by an environment struct; none is wired into the demographic
stepper yet. Design: [`archive/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md`](../../archive/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md).

Structural **allometry** lives in `src/allometry/` (a shared foundation below `state`, not here — see
issue #11), and the **empirical vital rates** (growth/mortality/recruitment) live in `src/demography/`
(they are empirical demographic functions, not ecophysiology).

## Contents

| Domain | Seam | Modules |
|--------|------|---------|
| **Types** | — | `meds_plant_types` (ALL derived types: leaf / hydraulics / phenology; one module, sectioned) |
| **Leaf gas exchange** | `meds_leaf_physiology%leaf_gas_exchange(env, cfg, ipft, flux)` | `meds_leaf_photosynthesis` (FvCB C3 + Collatz C4), `meds_leaf_stomata` (Leuning / Medlyn / Katul), `meds_leaf_solver` (bracketed Ci root-find) |
| **Hydraulics** | `meds_plant_hydraulics` | `meds_hydro_conductance`, `meds_hydro_pv`, `meds_hydro_solver` |
| **Phenology** | `meds_plant_phenology%update_phenology(env, params, dt, state, out)` | `meds_pheno_engine` |
| **Respiration** *(planned)* | `meds_plant_respiration` | stem + fine-root maintenance respiration + growth respiration |
| **Python C-API** | — | `meds_plant_capi` → `libmeds_plant_c` (`-DMEDS_BUILD_PYLIB=ON`; GLOB `*_capi.f90`) |

The shared temperature response (`meds_temp_response`, Arrhenius / peaked deactivation) lives in
`meds_shared` so leaf, respiration, and any tissue reach it without a plant→plant library edge.

## Python

The C-API shim is compiled only into the optional shared library `libmeds_plant_c`, exposed through the
process-oriented `meds.leaf` Python package (`python/meds/`), which reproduces Slot & Winter (2017) in
`examples/example_leaf_physiology/`.
