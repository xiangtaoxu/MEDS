# `fast_dynamics/` — the sub-daily half of the model

**Everything that runs on the `dt_fast` tier.** The kernels are mostly-stateless physical flux
calculators: device-eligible, netCDF-free (forcing arrives as passed-in value types, never a direct
`use netcdf`), and free of `polygon_t`, so they compile and unit-test standalone.

| folder | what it is |
|---|---|
| `processes/` | the **process modules**, one per process, named after its science page: canopy radiation, canopy aerodynamics, the canopy air space; leaf gas exchange, plant respiration, the plant's energy, the canopy's water, plant hydraulics; soil water, soil heat, the surface (snow and bare ground). See [`processes/README.md`](processes/README.md) |
| `integrator/` | the integrator, by stage: the records and the march state (`meds_fast_integ_types`), error control (`_control`), the right-hand sides (`_derivs`), the backward-Euler solves (`_backward_euler`), the two schemes (`_esdirk2`, `_rk45`), the frozen record (`_frozen`), the update and ledgers (`_update`) and the step itself (`_stepper`) |
| `driver/` | the loop that walks one slow step in `dt_fast` sub-steps over the patch axis (`meds_fast_dynamics`): the gather and scatter of each patch's state and the sums the slow step reads (`meds_fast_state`), the forcing (`meds_fast_forcing`) and what the output reads (`meds_fast_diag`) |

`integrator/` and `driver/` may see `polygon_t`; **`processes/` may not**, and that is what
keeps them device-eligible and standalone-buildable. The slow half is `src/slow_dynamics/`; the
split between them is by **timescale**, and within each half by domain. See
[`../README.md`](../README.md) for the tree-wide rules.

Each process module owns its argument records, beside the kernel that reads them, named by role
(`src/README.md`, rule 3): `*_env_t` what a kernel reads, `*_args_t` what is fixed for the slow
step, `*_out_t` what it returns. `rad_env_t`/`rad_out_t` are `meds_canopy_radiation`'s,
`aero_env_t`/`aero_args_t`/`aero_out_t` `meds_canopy_aerodynamics`', `cas_env_t`
`meds_canopy_air_space`'s, `leaf_env_t`/`leaf_out_t` `meds_leaf_gas_exchange`'s,
`soil_heat_env_t`/`soil_heat_out_t` `meds_soil_heat`'s, `surface_env_t`/`snow_out_t` `meds_surface`'s,
`ground_evap_env_t` `meds_soil_water`'s. The
`SOIL_*` / `ENERGY_*` selector codes live one layer down in
`config/params/meds_soil_params`, and the prognostic per-store column state — the canopy-air state, the
two soil columns, the snow store — lives in `state/meds_patch_state_vars`, with the rest of what a
patch holds, and the kernels take those types from there.

## The process families

Grouped **by surface subsystem**: one module per thermal or chemical store, with the
radiative-transfer pair on the side. There is no façade — each kernel module exposes its own seams,
so every symbol has exactly one legal spelling.

- **Canopy radiative transfer** — ED2 two-stream (`icanrad = 2`). The pure optical-property kernels
  (the Beta leaf-angle distribution, the Ross `G` function, the single-scatter pair) live in the
  shared `meds_optics_lib` in `src/share/functions/`. The optics assembly, the unified multi-band
  VIS / NIR / LW adding solver, and the public seam `canopy_radiation` are together in
  **`meds_canopy_radiation`**.
- **Canopy aerodynamics** — **`meds_canopy_aerodynamics`**: a CLM5 Monin-Obukhov surface layer, ED2
  Nusselt leaf and wood boundary layers, per-cohort in-canopy wind extinction, a CLM ground
  conductance, and the scalar-transfer factors that set the shared friction-velocity conductance for
  all three canopy-air twins.
- **Canopy air space** — **`meds_canopy_air_space`**: the three prognostic twins (specific enthalpy,
  specific humidity, molar CO₂), all advanced by one shared box kernel in two forms — a tendency for
  the explicit path and an implicit step, implicit in the exchange with the atmosphere. The
  integrator assembles the summed surface and biotic sources and the capacities, conductances and
  atmospheric boundary conditions.
- **Soil and plant water** — solved in `integrator/` (`water_step`): the pond, every
  soil layer and every cohort's wood and leaf water are stores in one network, joined by flow
  paths (soil faces, infiltration, root uptake, xylem, transpiration). Each `dt_fast` it is solved
  implicitly, with Newton iterations on each store's own unknown and adaptive sub-steps. The
  processes stay with their owners: **`meds_soil_water`** has the soil water flux between layers,
  the dry-surface-layer ground evaporation and the bottom boundary (drainage at a factor of free
  drainage, from 1 to 0 for bedrock, or an aquifer);
  **`meds_plant_hydraulics`** the plant's flows (sapflow through the Kirchhoff mean of the
  vulnerability curve, root uptake, the leaves' supply limit); `meds_plant_hydr_lib` the tissues'
  water relations and the plant's conductances, and `meds_soil_hydr_lib` the soil layers' and the
  pond's (both in `src/share/functions/`). Closes a machine-precision water budget.
- **Soil heat** — **`meds_soil_heat`**: the soil-heat store as prognostic **internal energy**,
  not temperature, so freeze/thaw is a read-off of the shared inverter. An implicit
  backward-Euler Thomas heat-diffusion solve, with an explicit tendency sibling.
- **Plant energy** — **`meds_plant_energy`**: the leaf and wood tissue energy balance, and the
  coefficients it is built from; the canopy interception film it evaporates is
  **`meds_canopy_water`**'s. The tissue relaxes **exactly** over the step, because under
  the frozen coefficients its ODE is linear: the kernel uses the closed form with two weights, an
  endpoint weight for the committed state and a step-average weight for every reported flux. Pairing
  them is what makes the balance close identically; using one for both does not. "Diagnostic" is the
  zero-heat-capacity limit of that one formula, which is why there is no leaf/wood energy-model
  selector.
- **Surface** — **`meds_surface`**: the bare-ground skin energy balance and the snow store —
  Niu-Yang cover fraction, snowfall and rain-on-snow
  accumulation, meltwater percolation, the snow-surface energy balance, and the snow-base to
  soil-top conductance. The two are mutually exclusive modes of one interface, blended by the cover
  fraction.

**Respiration** — **`meds_respiration`**: the plants' stem and fine-root maintenance
respiration, and the soil's heterotrophic respiration — the decomposition scalar, the CENTURY
transfer matrix and Rh. The daily step of the soil carbon pools
(`slow_dynamics/processes/meds_soil_carbon`) integrates the same matrix, so the sub-daily
respiration debits the same CENTURY pool the daily step does — the one documented fast/slow kernel
seam in the model.

## Shared constitutive kernels (in `src/share/functions/`)

The soil **material-property** kernels are stateless and `elemental`, grouped with the other
constitutive relations by physical quantity rather than by caller:

- **Retention curves** → `meds_soil_hydr_lib`: the water content, potential, conductivity and moisture
  capacity relations (van Genuchten by default, Campbell available) plus the `SOIL_RETENTION_*`
  selectors, and the pedotransfer function that sets them from the soil maps; the plant's
  pressure-volume and vulnerability curves are the matching `meds_plant_hydr_lib`.
- **Thermal properties** → `meds_therm_lib`: the ice-aware Johansen conductivity and the volumetric
  heat capacity — beside the moist-air psychrometrics, which are their thermal twin.
- **Per-column parameter bundles** → `meds_soil_params` in `src/config/params/`, and their
  `pure` builders → `meds_config` in `src/config/`, with the soil-column block they are derived from.
  The builders are state-free: they do not depend on water content or energy.

## Coupling

The stateless kernels are woven per sub-step by **`meds_fast_integ_stepper`**: the frozen record, the
water network's step, one adaptive march of **`meds_fast_integ_esdirk2`** (an L-stable ESDIRK2, the
default) or **`meds_fast_integ_rk45`** (adaptive Cash-Karp, the explicit alternative, redone on
ESDIRK2 when it cannot resolve a step), and the update. Every store closes a machine-precision
budget residual.

The integration story is in [`docs/science/column_biophysics.md`](../../docs/science/column_biophysics.md),
with per-store pages for the canopy air space, soil, vegetation energy and snow, and the time
integration itself in [`numerical_scheme.md`](../../docs/science/numerical_scheme.md).

**Tests.** One per process module, named after it: `test_canopy_aerodynamics`, `test_canopy_air_space`,
`test_canopy_radiation`, `test_canopy_water`, `test_leaf_gas_exchange`, `test_plant_energy`,
`test_plant_hydraulics`, `test_respiration`, `test_soil_heat`, `test_soil_water`, `test_surface`; the laws
in `src/share/functions/` in `test_therm_lib`, `test_optics_lib`, `test_plant_hydr_lib` and
`test_soil_hydr_lib`. The water network, and the soil column alone: `test_water_network`. The integrator,
one test per concern on the shared fixture `meds_test_fast_fixture`: `test_fast_integ_kernels` (the right-hand
sides, the backward-Euler stage, the state vector and the error control), `test_fast_integ_schemes` (ESDIRK2,
RK45 and the adaptive march) and `test_fast_integ_step` (one fast step end to end). The driver:
`test_fast_driver`.
