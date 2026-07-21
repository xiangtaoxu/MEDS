# biogeochemistry

The **slow soil-carbon / nutrient cycle** of the ecosystem column: the CENTURY-family
soil-organic-matter / litter pools, advanced daily. Like `biophysics/` it links `src/shared`
**only** and keeps every compute kernel stateless / `pure` / GPU-eligible (per-patch state, TOML
config, and cross-store coupling land at P3).

> **Note (module reorg):** the **fast** canopy-air-space CO2 exchange — `meds_column_co2` and its
> `co2_opts_t` / `cohort_co2_flux_t` / `column_co2_budget_t` / `damm_params_t` types — is a
> **sub-daily biophysical** process (turbulent diffusion / venting of the third CAS twin), so it now
> lives in `src/biophysics/` (types in `meds_biophysics_types`). `biogeochemistry/` is therefore the
> **slow** soil-carbon half only. The `heterotrophic_respiration_flux` kernel it re-uses for the
> fast/slow-seam reconciliation is imported from `meds_column_co2` (biophysics).

**Implemented — P0 slow soil-carbon matrix** (design `docs/dev_plans/MEDS_BIOGEOCHEMISTRY_DESIGN.md`):
- `meds_biogeochem_types` — shared derived types + selector codes: the SLOW `decomp_opts_t` /
  `litter_input_t` / `soilc_audit_t` / `soilc_diag_t`, the (7-pool + lignin + optional N)
  `soil_carbon_t`, and the `n_soil_pool` / `IP_*` / `DECOMP_STEP_*` / `DECOMP_SCHEME_*` parameters.
- `meds_soil_biogeochem` — the ED2-faithful CENTURY decomposition network organized as the carbon
  matrix ODE `dX/dt = B·I + A·xi·K·X`, as stateless `pure`-where-possible kernels: `assemble_env_scalar`
  (per-pool temperature × moisture × oxygen + lignin brake, matched to the fast Rh chemistry),
  `assemble_transfer_matrix` (scheme-0 3-active / scheme-5 5-active CENTURY topology; `a_jj=-1`,
  `er_j = 1-Σa_ij`), `build_litter_input`, the daily `soil_carbon_step` (forward **EULER** decrementing
  each donor by the fast loop's *accumulated* loss integral `xi_int`, plus an exact augmented **EXPM**
  matrix-exponential for large accelerated steps; carbon-mass + lignin-tracer audit), the respired
  complement `heterotrophic_respiration_matrix = -1ᵀA·xi·K·X`, the **SASU** `solve_soil_carbon_steady_state`
  (active-block `{K_j>0}` Gaussian solve — the full 7×7 is singular with inert pools), and
  `soil_carbon_diagnostics` (storage capacity/potential, residence time). Tested in
  `test/test_soil_biogeochem.f90` (mass closure, Rh complement, scalar placement, scheme topology, SASU,
  residence/capacity, litter+lignin, T/θ response, and the fast/slow-seam Jensen reconciliation).

**Reserved follow-ups** (see the design docs): DAMM heterotrophic respiration (`HR_DAMM`, P1) behind a
generalizable soil-gas-flux seam for CH4/N2O; a multi-layer canopy air space (P2, the `n=1` box is the
degenerate case) sharing an eddy-diffusivity transport column with the energy/water twins; the optional
nitrogen twin (P1, shaped-in: `n_cycle_on` fields present, `f_decomp ≡ 1` when off); vertically-resolved
soil-C pools (P2); and P3 state + config + the fast-loop master-step coupling — per-patch `soil_carbon_t`
in `state/`, `[soil_carbon]` TOML + `ed_params.f90` provenance check, netCDF restart serialization, and
the demography→litter→`soil_carbon_step`→CAS-Rh driver seam (single Rh authority, `xi_int` accumulator).
