# `slow_dynamics/soil/` — the daily soil carbon cycle

The **slow soil-carbon cycle** of the ecosystem column: the CENTURY-family soil-organic-matter and
litter pools, advanced once a day. Links the column state, the config leaves and the fast kernels'
`meds_respiration`, whose decomposition and respiration chemistry it integrates, and every compute
kernel is stateless, `pure` where possible, and device-eligible.

**The fast canopy-air CO₂ exchange is not here.** It is a sub-daily biophysical process —
turbulent diffusion and venting — and lives in `fast_dynamics/processes/meds_canopy_air_space`. This
folder is the slow carbon half: it advances the pools whose heterotrophic respiration the fast loop
charges, across the one documented fast/slow kernel seam.

## Modules

- **`meds_biogeochem_types`** — the shared derived types and selector codes: the decomposition
  options, the litter input record, the carbon and lignin audits, the traceability diagnostics, the
  seven-pool `soil_carbon_t` with its lignin sub-tracer, and the pool index and scheme parameters.
- **`meds_soil_biogeochem`** — ED2's CENTURY decomposition network, reorganized as the carbon matrix
  ODE `dX/dt = B·I + A·ξ·K·X`, as stateless kernels. ξ, A, K and Rh are `meds_respiration`'s
  (`decomposition_env_scalar`, `decomposition_transfer_matrix`, `heterotrophic_respiration`):
  - `build_litter_input` — where dead plant carbon enters;
  - `soil_carbon_step` — the daily advance, as forward Euler over the fast loop's *accumulated*
    environmental integral, or as an exact augmented matrix exponential for large accelerated steps;
  - `solve_soil_carbon_steady_state` — the SASU solve, over the active block only, because the full
    system is singular once a pool is inert;
  - `soil_carbon_diagnostics` — storage capacity and potential, and residence time.
- **`meds_litter_partition`** — the necromass-to-litter destination split. A biogeochemical law, not
  a demographic one, which is why it lives here and the demography modules call it.

## The seam that matters

The fast loop's heterotrophic respiration respires the **same frozen pool** the daily step debits,
through the one matrix both take from `meds_respiration`. That shared matrix is deliberate and it is the reason the day's
total sub-daily respiration equals the daily pool debit *by construction* rather than approximately.
The residual is asserted every step and closes to machine precision.

This is the one place in MEDS where a kernel is called from both timescale tiers. It is a documented
seam, not a folder.

## Science and tests

The equations, the pool table, the scalar placement and the seam are documented in
[`docs/science/soil_carbon.md`](../../../docs/science/soil_carbon.md); the design record is
[`docs/dev_plans/archive/MEDS_BIOGEOCHEMISTRY_DESIGN.md`](../../../docs/dev_plans/archive/MEDS_BIOGEOCHEMISTRY_DESIGN.md).

Tested in `test/test_soil_biogeochem.f90` (mass closure, the respiration complement, scalar
placement, scheme topology, the steady-state solve, residence and capacity, litter and lignin,
temperature and moisture response, and the fast/slow seam reconciliation) and
`test/test_biogeochem_dynamics.f90` (the daily driver).

## Not here yet

There is no nitrogen: the soil is carbon-only, and the nitrogen cycle is not built (it is planned
after the move to ED-MORE).
There is no DAMM decomposition-moisture alternative; an unwired kernel for it is kept on branch
`archive/damm-hr`.
Vertically resolved pools, nitrogen limitation of productivity, fire, and an explicit coarse woody
debris pool are all unbuilt. See [`docs/ROADMAP.md`](../../../docs/ROADMAP.md) §3.
