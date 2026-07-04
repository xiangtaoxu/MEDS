# plant/hydraulics

Per-individual **plant hydraulics** — a stateless xylem water-transport / water-stress compute
library (the prognostic state, e.g. xylem water potential, lives in the cohort SoA under
`src/state/` and is passed by argument — the FATES `*Mem`/compute split). Links `meds_shared`
only, no `site_t`, so it compiles and tests standalone (`cmake --build … --target meds_hydraulics`),
exactly like `src/plant/leaf/` and `src/biophys/`.

Design: `archive/MEDS_HYDRAULICS_DESIGN.md`.

## Public seam

`meds_plant_hydraulics%plant_water_flux(env, params, opts, dt, psi, flux)` advances a cohort's node
water potentials `psi(N_HYDRO)` [MPa] over one fast step `dt`, given the boundary conditions
(`hydro_env_t`: transpiration demand, soil potential, rhizosphere conductance, per-plant
biomass/geometry), the per-PFT traits (`hydro_params_t`), and the run selectors (`hydro_opts_t`).

## Modules

- `meds_hydro_types`       — the interface types + node/solver/conductance enums + compile-time `N_HYDRO`.
- `meds_hydro_pv`          — nonlinear pressure–volume curve (Tyree & Hammel 1972 symplast in the
  Bartlett et al. 2012 form + Christoffersen et al. 2016 apoplast): `water_content`, `capacitance`,
  closed-form `psi<->RWC`, and the legacy-linear `pv_water_cap_from_traits` mapping.
- `meds_hydro_conductance` — vulnerability `plc_retained` and the **Kirchhoff flux potential** Φ
  (`flux_potential`, `phi_inverse`, `kirchhoff_edge`): `plc` enters only through `Φ = ∫plc dψ`,
  never pointwise (the single conductance law).
- `meds_hydro_solver`      — the coupled adaptive integrator `solve_plant_water` (2-node leaf+wood,
  frozen-coefficient matrix exponential + step-doubling; mass-closing ΔW-based fluxes).
- `meds_plant_hydraulics`  — the sealed public seam.

## Status

Implemented: the standalone 2-node (leaf + lumped wood) network with whole-plant conductance, the
nonlinear PV curve, the Kirchhoff conductance, the matrix-exponential adaptive solver, and
`test/test_plant_hydraulics.f90` (PV, Kirchhoff, conservation, steady-state Ohm's law, a fine
explicit-reference match, diurnal cycle, degenerate cohorts). Validated on ifx and nvfortran multicore.

Follow-ups (see the design doc): the 3-node (leaf/stem/root) backward-Euler path and the stem ψ(z)
diagnostic; the cohort-SoA `psi_node` state + demography wiring; the `meds_config` `[hydraulics]`
block + TOML traits + the config-driven seam; per-soil-layer fine-root nodes.
