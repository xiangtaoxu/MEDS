# biophysics

**Fast, sub-daily, mostly-stateless physical flux calculators.** Device-eligible and netCDF-free
(forcing enters via passed-in value types, never a direct `use netcdf`). Sealed and orthogonal to the
demographic engine: the modules here link `meds_shared` only (no `site_t`), so they compile and test
standalone, exactly like `src/plant/`. The boundary against `biogeochemistry/` is **by domain, not by
timescale** — biophysics owns the fast energy / water / momentum / CO2 exchange physics; the **slow**
soil-carbon pools live in `biogeochemistry/`.

Shared derived types + `SOIL_*` / `ENERGY_*` / `HR_*` selector codes are consolidated in
**`meds_biophysics_types`**; the prognostic per-store column state (`cas_state_t`, the two soil
columns, the snow store) lives in `src/shared/state/meds_column_state_types` so the demographic state
hub can own it. The science pages are `docs/science/canopy_radiation_transfer.md`,
`canopy_aerodynamics.md`, and `column_biophysics.md`.

## Process families

The physical kernels are grouped **by surface subsystem** (one module per thermal/chemical store),
with the radiative-transfer pair on the side and a logic-free re-export façade
(**`meds_biophysics_interface`**, the analogue of `meds_plant_interface`) exposing every seam through
one `use`.

- **Canopy radiative transfer** — ED2 two-stream (`icanrad=2`). The **pure optical-property kernels**
  (leaf-angle Beta distribution, Ross `G(mu)`, `omega`/`g` `scatter_pair`, the `beta_*`/`leaf_bf`/
  `gfun_direct`/`leaf_class_angle` family) live in the shared library **`meds_optics_lib`**
  (`src/shared/functions/`). The **RT assembly** (`derive_rad_optics`, `blend_cohort_optics`,
  `ground_optics`), the unified multi-band (VIS/NIR/LW) O(N) adding solver (`solve_band`/`layer_rt`),
  and the public seam `canopy_radiation` all live together in **`meds_canopy_radiation`**. See
  `docs/science/canopy_radiation_transfer.md`.
- **Canopy aerodynamics** — **`meds_canopy_aerodynamics`**: CLM5 Monin-Obukhov surface layer, ED2
  Nusselt leaf/wood boundary layers, per-cohort in-canopy wind extinction, CLM ground conductance, and
  the `temp1`/`temp2` scalar-transfer factors that set the shared `ustar`-based conductance for all
  three CAS twins. See `docs/science/canopy_aerodynamics.md`.
- **Soil water** — **`meds_soil_water`**: implicit backward-Euler Thomas Richards
  (`column_hydrology_flux` / `soil_water_step_implicit` / `soil_water_advance`; Celia/frozen
  linearization, upstream K, Zeng-Decker equilibrium, adaptive substepping, infiltration/ponding, Dunne
  runoff, DSL soil evaporation, psi-limited root sink, free-drain/bedrock/aquifer BC); the
  `column_hydrology_flux` seam and `ground_evaporation` live here too. Closes a machine-precision water
  budget.
- **Soil thermal** — **`meds_soil_energy`**: the soil-heat store (`soil_energy_step_implicit`, its
  explicit sibling `soil_energy_time_deriv`, and the `soil_heat_be_solve` BE-Thomas heat-diffusion
  solve). Prognostic **internal energy** (not temperature), so freeze/thaw is a shared-inverter read-off.
- **Vegetation biophysics** — **`meds_vegetation_biophysics`**: the leaf/wood energy store
  (`veg_energy_step_implicit`, `veg_surface_fluxes`) plus per-cohort canopy interception
  (`intercept_canopy_layer`). Prognostic internal energy, same freeze/thaw read-off.
- **Ground biophysics** — **`meds_ground_biophysics`**: the ground-skin balance
  (`ground_surface_balance`) and the full snow / temporary-surface-water store (all `snow_*` kernels —
  Niu-Yang cover fraction, snowfall/rain-on-snow accumulation, meltwater percolation, snow-surface
  energy balance, and the snow-base → soil-top conductance).
- **Canopy-air-space (CAS) biophysics** — **`meds_cas_biophysics`**: the three prognostic CAS twins.
  Enthalpy and humidity via `canopy_air_update`; the molar CO2 twin `can_co2 [umol/mol]` via
  `canopy_air_co2_update` (implicit atmosphere exchange, closed budget, mirroring `canopy_air_update`);
  plus `aggregate_cohort_co2_fluxes` and the shared two-form CAS box (`cas_column_time_deriv` /
  `cas_column_step_implicit`, called by both integrators). The fast CO2 is a diffusion/venting exchange
  — hence biophysics. Heterotrophic soil **respiration** (`heterotrophic_respiration_flux`/`_damm`) is a
  carbon-decomposition process, so it lives in `biogeochemistry` (`meds_soil_biogeochem`); the driver is
  its single authority and passes the resulting CO2 source into the CAS box.

## Shared constitutive kernels (in `src/shared/`)

The soil **material-property** kernels are stateless, `elemental`, scalar-in, grouped with the other
constitutive relations by physical quantity (the soil analogue of the tissue curves):

- **Retention curves** → `meds_hydr_lib` (`shared/functions/`): `soil_theta_from_psi` /
  `soil_psi_from_theta` / `soil_hydr_cond_from_theta` / `soil_moist_cap_from_psi` (van Genuchten
  default + Campbell) + the `SOIL_RETENTION_*` selectors — beside the plant PV / vulnerability curves.
- **Thermal properties** → `meds_therm_lib` (`shared/functions/`): `soil_thermal_cond` (Johansen,
  ice-aware) / `soil_heat_cap_vol` — beside the moist-air psychrometrics (the thermal twin).
- **Per-column parameter bundles + their `pure` builders** → `meds_column_state_types`
  (`shared/state/`): `soil_params_t` + `build_soil_hydr_params`, `soil_thermal_params_t` +
  `build_soil_therm_params` — beside the prognostic soil columns they describe. The builders are
  state-free constructors (no `theta`/energy dependence).

`meds_biophysics_types` re-exports `soil_params_t` / `soil_thermal_params_t` / `SOIL_RETENTION_*` so the
biophysics kernels and callers keep `use meds_biophysics_types` unchanged.

## Coupling

The stateless kernels are woven per fast sub-step by the driver **`meds_column_dynamics`** (operator-
split or IMEX-ARK via `meds_column_derivs` / `meds_ark_stepper`), with leaf↔CAS Picard coupling; every
store closes a machine-precision budget residual. See `docs/science/column_biophysics.md` for the full
integration story. Individual kernels are exercised by `test/test_{canopy_radiation,aerodynamics,
column_hydrology,column_energy,surface_energy,snow,column_co2}.f90`; the coupled loop by
`test/test_{column_dynamics,column_derivs,picard_coupling,column_ark,fast_loop}.f90`.
