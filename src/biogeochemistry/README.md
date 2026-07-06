# biogeochemistry

The **carbon / nutrient cycle** of the ecosystem column. The clean boundary against `biophysics/`
is **by domain, not by timescale**: `biophysics/` holds the fast energy / water / momentum physics;
`biogeochemistry/` holds carbon — spanning the **fast** canopy-air CO2 exchange *and* the **slow**
soil-organic-matter / litter pools. Like `biophysics/` it links `src/shared` **only** and keeps every
compute kernel stateless / `pure` / GPU-eligible (per-patch state, TOML config, and cross-store
coupling land at P3).

**Implemented — P0 column CO2 balance** (design `archive/MEDS_COLUMN_CO2_BALANCE_DESIGN.md`):
- `meds_biogeochem_types` — shared derived types + `HR_*` selector codes (`co2_opts_t`,
  `soil_carbon_t`, `cohort_co2_flux_t`, `column_co2_budget_t`).
- `meds_column_co2` — the canopy-air-space CO2 balance: `can_co2 [umol/mol]` is the **third
  prognostic CAS twin** beside `can_enthalpy` / `can_shv`, advanced by `canopy_air_co2_update`
  (molar capacity, implicit atmosphere exchange, closed budget) exactly mirroring
  `meds_column_energy%canopy_air_update`; plus `aggregate_cohort_co2_fluxes`
  (per-cohort → `[umol CO2 / m2 ground / s]`), `heterotrophic_respiration_flux` (MVP Q10 / ED2
  capped-exp × moisture on a frozen soil-C pool), and the `column_co2_step` assembler
  (NEE / NEP / loss-to-atmosphere + a machine-precision residual). Tested in `test/test_column_co2.f90`.

**Reserved follow-ups** (see the design doc): DAMM heterotrophic respiration (`HR_DAMM`, P1) behind a
generalizable soil-gas-flux seam for CH4/N2O; a multi-layer canopy air space (P2, the `n=1` box is the
degenerate case) sharing an eddy-diffusivity transport column with the energy/water twins; the slow
CENTURY soil-carbon pool dynamics; and P3 state + config + the fast-loop master-step coupling
(`ca = can_co2` feedback on photosynthesis).
