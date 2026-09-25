# `fast_dynamics/` — the sub-daily half of the model

**Everything that runs on the `dt_fast` tier.** The kernels are mostly-stateless physical flux
calculators: device-eligible, netCDF-free (forcing arrives as passed-in value types, never a direct
`use netcdf`), and free of `site_t`, so they compile and unit-test standalone.

| folder | what it is |
|---|---|
| `canopy/` | the **medium** — radiative transfer, aerodynamics, the canopy air space |
| `plant/` | the **organisms** — leaf gas exchange, hydraulics, maintenance respiration, tissue energy |
| `soil/` | the **ground** — soil water, soil energy, the ground skin, the snow store |
| `numerics/` | the integrator machinery: the state vector, the frozen work record, ARK and RK45, the pre-pass, error control |
| `driver/` | the loop that walks one slow step in `dt_fast` sub-steps over the patch axis |

`numerics/` and `driver/` may see `site_t`; **the three kernel folders may not**, and that is what
keeps them device-eligible and standalone-buildable. The slow half is `src/slow_dynamics/`; the
split between them is by **timescale**, and within each half by domain. See
[`../README.md`](../README.md) for the tree-wide rules.

Each domain owns its argument records: **`meds_canopy_types`** (radiative transfer and
aerodynamics), **`meds_soil_types`** (hydrology, thermal, snow), **`meds_plant_types`** (leaf,
hydraulics, tissue energy). The `SOIL_*` / `ENERGY_*` / `HR_*` selector codes live one layer down in
`config/meds_biophysics_opts`, and the prognostic per-store column state — the canopy-air state, the
two soil columns, the snow store — lives in `state/column`, so the demographic state hub can own it.

## The process families

Grouped **by surface subsystem**: one module per thermal or chemical store, with the
radiative-transfer pair on the side. There is no façade — each kernel module exposes its own seams,
so every symbol has exactly one legal spelling.

- **Canopy radiative transfer** — ED2 two-stream (`icanrad = 2`). The pure optical-property kernels
  (the Beta leaf-angle distribution, the Ross `G` function, the single-scatter pair) live in the
  shared `meds_optics_lib` in `src/shared/functions/`. The optics assembly, the unified multi-band
  VIS / NIR / LW adding solver, and the public seam `canopy_radiation` are together in
  **`meds_canopy_radiation`**.
- **Canopy aerodynamics** — **`meds_canopy_aerodynamics`**: a CLM5 Monin-Obukhov surface layer, ED2
  Nusselt leaf and wood boundary layers, per-cohort in-canopy wind extinction, a CLM ground
  conductance, and the scalar-transfer factors that set the shared friction-velocity conductance for
  all three canopy-air twins.
- **Canopy air space** — **`meds_cas_biophysics`**: the three prognostic twins (specific enthalpy,
  specific humidity, molar CO₂), all advanced by one shared box kernel in two forms — a tendency for
  the explicit path and an implicit step — both implicit in the exchange with the atmosphere. The
  driver assembles the summed surface and biotic sources and the capacities, conductances and
  atmospheric boundary conditions.
- **Soil water** — **`meds_soil_water`**: implicit backward-Euler Thomas Richards with Celia
  modified-Picard or frozen-coefficient linearization, upstream-weighted conductivity, adaptive
  substepping, conductivity-limited infiltration and ponding, dry-surface-layer evaporation, a
  ψ-limited root sink, and a free-drain / bedrock / aquifer bottom boundary. Closes a
  machine-precision water budget.
- **Soil thermal** — **`meds_soil_energy`**: the soil-heat store as prognostic **internal energy**,
  not temperature, so freeze/thaw is a read-off of the shared inverter. An implicit
  backward-Euler Thomas heat-diffusion solve, with an explicit tendency sibling.
- **Vegetation** — **`meds_plant_biophysics`**: the leaf and wood tissue energy balance and the
  canopy interception film they share. The tissue relaxes **exactly** over the step, because under
  the frozen coefficients its ODE is linear: the kernel uses the closed form with two weights, an
  endpoint weight for the committed state and a step-average weight for every reported flux. Pairing
  them is what makes the balance close identically; using one for both does not. "Diagnostic" is the
  zero-heat-capacity limit of that one formula, which is why there is no leaf/wood energy-model
  selector.
- **Ground and snow** — **`meds_ground_biophysics`**: the bare-ground skin energy balance and the
  full snow / temporary-surface-water store — Niu-Yang cover fraction, snowfall and rain-on-snow
  accumulation, meltwater percolation, the snow-surface energy balance, and the snow-base to
  soil-top conductance. The two are mutually exclusive modes of one interface, blended by the cover
  fraction.

**Heterotrophic respiration is not here.** It is a carbon-decomposition process and lives in
`slow_dynamics/soil/meds_soil_biogeochem`. The fast loop calls its matrix form so that the sub-daily
respiration debits the same CENTURY pool the daily step does — the one documented fast/slow kernel
seam in the model.

## Shared constitutive kernels (in `src/shared/functions/`)

The soil **material-property** kernels are stateless and `elemental`, grouped with the other
constitutive relations by physical quantity rather than by caller:

- **Retention curves** → `meds_hydr_lib`: the water content, potential, conductivity and moisture
  capacity relations (van Genuchten by default, Campbell available) plus the `SOIL_RETENTION_*`
  selectors — beside the plant pressure-volume and vulnerability curves, which are the same kind of
  object.
- **Thermal properties** → `meds_therm_lib`: the ice-aware Johansen conductivity and the volumetric
  heat capacity — beside the moist-air psychrometrics, which are their thermal twin.
- **Per-column parameter bundles and their `pure` builders** → `meds_column_params` in
  `state/column`, beside the prognostic columns they describe. The builders are state-free: they do
  not depend on water content or energy.

## Coupling

The stateless kernels are woven per sub-step by **`meds_fast_ark`** (an L-stable ESDIRK2, the
default) or **`meds_fast_rk45`** (adaptive Cash-Karp, the accuracy baseline), dispatched by
`meds_fast_step`, which also owns the RK45-to-ARK stiff rescue. Every store closes a
machine-precision budget residual.

The integration story is in [`docs/science/column_biophysics.md`](../../docs/science/column_biophysics.md),
with per-store pages for the canopy air space, soil, vegetation energy and snow, and the time
integration itself in [`numerical_scheme.md`](../../docs/science/numerical_scheme.md).

**Tests.** Individual kernels: `test_canopy_radiation`, `test_aerodynamics`, `test_column_hydrology`,
`test_column_energy`, `test_surface_energy`, `test_snow`, `test_column_co2`. The coupled loop:
`test_column_dynamics`, `test_column_derivs`, `test_column_ark`, `test_column_rk45`,
`test_fast_loop`.
