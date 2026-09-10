# fast_dynamics/plant — sub-daily plant ecophysiology

The **sub-daily** plant kernels: leaf gas exchange, plant hydraulics, non-leaf maintenance
respiration, and the leaf/wood energy balance. No cohort/patch state (`site_t`), so they compile and
unit-test standalone (`cmake --build … --target meds_fast_kernels`), orthogonal to demography. The
**daily** plant kernels — phenology, carbon allocation, trait dynamics — live in
`src/slow_dynamics/plant/`; the vital-rate laws are demography and live in
`src/slow_dynamics/demography/`. Every process is a
stateless per-individual kernel driven by an environment struct; none is wired into the demographic
stepper yet. Design: [`docs/dev_plans/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md`](../../docs/dev_plans/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md).

Structural **allometry** lives in `src/allometry/` (a shared foundation below `state`, not here — see
issue #11), and the **empirical vital rates** (growth/mortality/recruitment) live in `src/demography/`
(they are empirical demographic functions, not ecophysiology).

## Contents

Each kernel module exposes its own seams. There is no façade: `meds_plant_interface` was deleted
in the step-9 normalization, because it mixed re-export with config-flattening logic and a facade
must be one or the other. Its assemblers are driver code and live in `fast_dynamics/driver/`.

| File | Role | Contents |
|------|------|----------|
| `meds_plant_types` | types | ALL derived types (leaf / hydraulics / phenology), one module, sectioned |
| `meds_leaf_gas_exchange` | **the leaf seam** | `solve_leaf_gas_exchange` (one leaf) and `leaf_gas_exchange_batch` (bare arrays over n leaves). Both take a `leaf_photo_table_t`, never a config handle. The façade `meds_plant_interface` was deleted in the step-9 normalization: its config-flattening assemblers are driver code and live in `fast_dynamics/driver/meds_fast_config` |
| `meds_leaf_gas_exchange` | leaf compute | FvCB C3 + Collatz C4 demand, Leuning / Medlyn / Katul stomata, the bracketed Ci solver (`solve_leaf_gas_exchange`) |
| `meds_plant_hydraulics` | hydraulics compute | pressure-volume (Bartlett/Tyree-Hammel), Kirchhoff conductance, matrix-exp sub-step solver |
| `meds_phenology` | phenology compute | the cue engine → directional status (`phenology_kernel`) |
| `meds_plant_respiration` | respiration compute | non-leaf maintenance respiration: `stem_maintenance_respiration` + `fine_root_maintenance_respiration` + `growth_respiration` |
| (C-API) | Python C-API | `src/capi/meds_capi_leaf.f90` → the single `libmeds.so` (`-DMEDS_BUILD_PYLIB=ON`); calls the leaf compute kernels directly, not the seam |

The shared temperature response (`meds_temp_response`, Arrhenius / peaked deactivation) lives in
`meds_shared` so leaf, respiration, and any tissue reach it without a plant→plant library edge.

## Python

The `src/capi/meds_capi_*.f90` shims go into the single optional `libmeds.so` AND into one mandatory
ctest target each, so they cannot rot unnoticed. Exposed
through process-oriented Python packages (`python/meds/`): `meds.plant.leaf` (leaf gas exchange, reproduces
Slot & Winter 2017 in `examples/example_leaf_gas_exchange/`) and `meds.plant.pheno` (the leaf-phenology
kernel — the four phenology strategies in `examples/example_phenology/`).
