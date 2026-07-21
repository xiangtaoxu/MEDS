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

- **Canopy radiative transfer** — ED2 two-stream (`icanrad=2`): optics (leaf-angle Beta distribution,
  Ross `G(mu)`, `omega`/`g` scattering, ground optics) in **`meds_optics`**, the unified multi-band
  (VIS/NIR/LW) O(N) adding solve in **`meds_twostream_band`**, behind the seam
  **`meds_canopy_radiation`**. See `docs/science/canopy_radiation_transfer.md`.
- **Canopy aerodynamics** — **`meds_canopy_aerodynamics`**: CLM5 Monin-Obukhov surface layer, ED2
  Nusselt leaf/wood boundary layers, per-cohort in-canopy wind extinction, CLM ground conductance, and
  the `temp1`/`temp2` scalar-transfer factors that set the shared `ustar`-based conductance for all
  three CAS twins. See `docs/science/canopy_aerodynamics.md`.
- **Soil-column hydrology** — **`meds_column_hydrology`**: implicit backward-Euler Thomas Richards
  (Celia/frozen linearization, upstream K, Zeng-Decker equilibrium, adaptive substepping, infiltration/
  ponding, Dunne runoff, DSL soil evaporation, psi-limited root sink, free-drain/bedrock/aquifer BC),
  plus per-cohort interception; closes a machine-precision water budget.
- **Energy balance** — **`meds_column_energy`**: four stateless per-store kernels (leaf/wood
  `veg_energy_balance`, ground `ground_surface_balance`, canopy air space `canopy_air_update`, soil
  thermal column `soil_energy_flux`). Prognostic **internal energy / enthalpy** (not temperature), so
  freeze/thaw is a shared-inverter read-off.
- **Snow / temporary-surface-water** — **`meds_snow`** (the merged mass + energy sides): Niu-Yang
  cover fraction, snowfall/rain-on-snow accumulation, meltwater percolation, snow-surface energy
  balance, and the snow-base → soil-top conductance.
- **Canopy-air-space CO2 balance** — **`meds_column_co2`**: `can_co2 [umol/mol]` is the **third
  prognostic CAS twin** beside `can_enthalpy` / `can_shv`, advanced by `canopy_air_co2_update` (molar
  capacity, implicit atmosphere exchange, closed budget) exactly mirroring `canopy_air_update`; plus
  `aggregate_cohort_co2_fluxes`, `heterotrophic_respiration_flux` (Q10 / ED2 capped-exp × moisture, and
  the mechanistic `HR_DAMM` dual-Arrhenius), and the `column_co2_step` NEE/NEP assembler. This is a
  fast diffusion/venting exchange — hence biophysics, not biogeochemistry.

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
