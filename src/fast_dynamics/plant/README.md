# `fast_dynamics/plant/` — sub-daily plant ecophysiology

The **sub-daily** plant kernels: leaf gas exchange, plant hydraulics, non-leaf maintenance
respiration, and the leaf and wood energy balance. No cohort or patch state, so they compile and
unit-test standalone (`cmake --build … --target meds_fast_kernels`), orthogonal to demography.

**The daily plant kernels are next door** in `src/slow_dynamics/plant/` — phenology, carbon
allocation, trait plasticity. The vital-rate laws are demography, not ecophysiology, and live in
`src/slow_dynamics/demography/`. Structural allometry is a shared geometric constraint reachable by
state, demography, plant and config alike, so it lives in `src/functions/meds_allometry.f90`.

Every process here is a stateless per-individual kernel driven by an environment record, and all of
them are wired into the fast loop: the pre-pass calls gas exchange and the two maintenance-respiration
kernels each sub-step, and the integrator advances the hydraulic network and the tissue energy store.

## Contents

Each kernel module exposes its own seams; there is no façade.

| Module | Role |
|---|---|
| `meds_plant_types` | All the sub-daily derived types — leaf, hydraulics, tissue energy — in one sectioned module. Phenology's types are in `meds_phenology_types` under `slow_dynamics/plant/`. |
| `meds_leaf_gas_exchange` | **The leaf seam.** FvCB C3 and Collatz C4 demand, the electron-transport hyperbola, Leuning / Medlyn / Katul stomata, and the bracketed Cᵢ solver. `solve_leaf_gas_exchange` for one leaf, `leaf_gas_exchange_batch` over bare arrays. Both take a photosynthesis parameter table, never a config handle. |
| `meds_plant_hydraulics` | The coupled matrix-exponential network solve, over the tissue pressure-volume and Kirchhoff conductance curves in `meds_hydr_lib`. |
| `meds_plant_biophysics` | The leaf and wood tissue energy balance, and the canopy interception film they share. |
| `meds_plant_respiration` | Non-leaf maintenance respiration: `stem_maintenance_respiration` and `fine_root_maintenance_respiration`. Growth respiration is charged in the daily allocator, on realized growth. |

The shared temperature response (`meds_temp_response`, Arrhenius and peaked deactivation) lives in
`src/functions/`, so leaf, stem and root reach one code path rather than three drifting copies.

The config-flattening assemblers that turn `meds_config_t` into these kernels' option records are
driver code and live in `fast_dynamics/driver/meds_fast_config`.

## Python

The `src/capi/meds_capi_*.f90` shims go into the single optional `libmeds.so`
(`-DMEDS_BUILD_PYLIB=ON`) **and** into one mandatory ctest target each, so an ABI or signature
change is a build failure in a default build rather than a silent break in an optional one.

Exposed through process-oriented packages in `python/meds/`:

- `meds.plant.leaf` — leaf gas exchange. Reproduces Slot & Winter (2017) in
  [`examples/example_leaf_gas_exchange/`](../../../examples/example_leaf_gas_exchange/).
- `meds.plant.pheno` — the leaf-phenology kernel (its Fortran side is in `slow_dynamics/plant/`).
  The four phenology strategies are in
  [`examples/example_phenology/`](../../../examples/example_phenology/).

## Science

[`leaf_gas_exchange.md`](../../../docs/science/leaf_gas_exchange.md),
[`plant_hydraulics.md`](../../../docs/science/plant_hydraulics.md),
[`plant_respiration.md`](../../../docs/science/plant_respiration.md),
[`vegetation_energy_dynamics.md`](../../../docs/science/vegetation_energy_dynamics.md).
