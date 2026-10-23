# `fast_dynamics/plant/` — sub-daily plant ecophysiology

The **sub-daily** plant kernels: leaf gas exchange, plant hydraulics, non-leaf maintenance
respiration, and the leaf and wood energy balance. No cohort or patch state, so they compile and
unit-test standalone (`cmake --build … --target meds_fast_kernels`), orthogonal to demography.

**The daily plant kernels are next door** in `src/slow_dynamics/plant/` — phenology, carbon
allocation, trait plasticity. The vital-rate laws are demography, not ecophysiology, and live in
`src/slow_dynamics/demography/`. Structural allometry is a shared geometric constraint reachable by
state, demography, plant and config alike, so it lives in `src/share/functions/meds_allometry_lib.f90`.

Every process here is a stateless per-individual kernel driven by an environment record, and all of
them are wired into the fast loop: the frozen build calls gas exchange and the two maintenance-respiration
kernels each sub-step, the soil-plant water network (`water_step`, in `integrator/`) moves the plant
water, and the integrator advances the tissue energy.

## Contents

Each kernel module exposes its own seams; there is no façade.

| Module | Role |
|---|---|
| `meds_plant_types` | All the sub-daily derived types — leaf, hydraulics, respiration, tissue sidedness — in one sectioned module. Phenology's types are in `meds_phenology_types` under `slow_dynamics/plant/`. |
| `meds_leaf_gas_exchange` | **The leaf seam.** FvCB C3 and Collatz C4 demand, the electron-transport hyperbola, Leuning / Medlyn / Katul stomata, and the bracketed Cᵢ solver. `solve_leaf_gas_exchange` for one leaf, `leaf_gas_exchange_batch` over bare arrays. They take a PFT's photosynthesis parameters or the per-PFT table, never a config handle. |
| `meds_plant_hydraulics` | How water moves through a cohort's plants: the sapflow through the Kirchhoff mean of the vulnerability curve (each PFT's lookup built at load by `meds_pft_derive`), root water uptake and the leaves' supply limit on transpiration, each with its derivatives; and the plant-water facts the fast loop needs: the leaf node potential, the crown gravity head, the rooting depth and root profile. The tissue curves, the vulnerability curve and the conductances' laws are in `meds_plant_hydr_lib`; the soil-plant water network is solved in `integrator/` (`water_step`). |
| `meds_plant_biophysics` | The leaf and wood tissue energy balance, and the canopy interception film they share. |
| `meds_plant_respiration` | Non-leaf maintenance respiration: `stem_maintenance_respiration` and `fine_root_maintenance_respiration`. Growth respiration is charged in the daily allocator, on realized growth. |

The shared temperature response (`meds_therm_lib`, Arrhenius and peaked deactivation) lives in
`src/share/functions/`, so leaf, stem and root reach one code path rather than three drifting
copies.

The per-PFT tables these kernels read (the leaf parameters and the hydraulics) are derived from the
configuration once, at load, onto the PFT table (`config/meds_pft_derive`).

## Python

The `src/c_api/meds_c_api_*.f90` shims compile once, into the object library `meds_c_api`, which the
single optional `libmeds.so` (`-DMEDS_BUILD_PYLIB=ON`) **and** the mandatory C-API ctest targets
link, so an ABI or signature change is a build failure in a default build rather than a silent break
in an optional one.

Exposed through process-oriented packages in `python/meds/`:

- `meds.plant.leaf` — leaf gas exchange. Reproduces Slot & Winter (2017) in
  [`examples/example01_leaf_gas_exchange/`](../../../examples/example01_leaf_gas_exchange/).
- `meds.plant.pheno` — the leaf-phenology kernel (its Fortran side is in `slow_dynamics/plant/`).
  The four phenology strategies are in
  [`examples/example02_canopy_phenology/`](../../../examples/example02_canopy_phenology/).

## Science

[`leaf_gas_exchange.md`](../../../docs/science/leaf_gas_exchange.md),
[`plant_hydraulics.md`](../../../docs/science/plant_hydraulics.md),
[`plant_respiration.md`](../../../docs/science/plant_respiration.md),
[`vegetation_energy_dynamics.md`](../../../docs/science/vegetation_energy_dynamics.md).
