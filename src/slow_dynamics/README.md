# `slow_dynamics/` — the daily-to-annual half of the model

**Everything that runs once per slow step, or at a calendar boundary.** The process modules are the
laws: stateless, `pure` where possible, and free of `polygon_t`, so they compile and unit-test
standalone. The driver reads the stand, calls the laws, and writes the stand back.

| folder | what it is |
|---|---|
| `processes/` | the **process modules**, one per process, named after its science page. See [`processes/README.md`](processes/README.md) |
| `driver/` | the slow loop over the stand: its order, the reads and writes around each process, the one writer of the stand's carbon and structure, the slow diagnostic rows and the site ledger |

**`processes/`** (library `meds_slow_kernels`):

- `meds_plant_phenology` — the daily cues become the leaf flush and shed tendencies; the leaf
  turnover step.
- `meds_plant_carbon_allocation` — a cohort's carbon demand, the sink limit on wood growth, the
  priority ladder and growth respiration.
- `meds_plant_traits` — the light-plastic leaf traits, and a site's leaf table at its growth
  temperature.
- `meds_soil_carbon` — dead plant carbon into litter, the daily step of a patch's soil carbon pools,
  the steady-state solve and the pools' plausibility.
- `meds_demography` — the demographic laws: the growth rate, the mortality hazard, the decline to the
  density floor, the recruit efficiency and a cohort's share of the recruitment; and the decisions
  the update applies: whether a recruit pool founds a cohort, the share of the stand a treefall
  opens, and who survives it.

**`driver/`** (library `meds_slow`):

- `meds_slow_dynamics` — the order of the slow step and of the calendar boundary, and every
  slow-ledger mark. It loops over nothing.
- `meds_slow_state` — reads the stand, calls the processes and writes each process's own fields:
  the phenology memory, the plastic traits, the tissue water, the litter. What changes the stand's
  carbon or structure it hands on: the step's tendencies in the step record `slow_step_data_t`,
  a trait change's leaf and storage carbon, the soil step's new pools, and the founding and
  treefall decisions.
- `meds_slow_update` — the one place a slow step writes the stand's carbon and structure: the
  step's commit (the cohorts' pools and density, the trait change's carbon, the soil pools),
  recruitment and disturbance. It also sorts, fuses, splits and culls the cohorts, and sorts,
  fuses and terminates the patches.
- `meds_slow_diag` — the slow rows of the diagnostic blocks.
- `meds_slow_ledger` — the site conservation ledger across a slow step, with the run's worst
  soil-carbon seam checks.

`driver/` may see `polygon_t`; **`processes/` may not**, and ctest `kernels_patch_only` checks it.
The fast half is `src/fast_dynamics/`; the split between them is by **timescale**. See
[`../README.md`](../README.md) for the tree-wide rules.

Each process module owns its argument records, beside the kernel that reads them (`src/README.md`,
rule 3): `pheno_env_t`, `pheno_state_t` and `pheno_out_t` are
`meds_plant_phenology`'s, and `soil_carbon_out_t` is `meds_soil_carbon`'s. The
patch stores a slow step advances — the soil carbon pools (`soil_carbon_state_t`), the litter
(`litter_input_t`) and the decomposition integral (`decomp_xi_t`) — live in
`state/meds_patch_state_vars`. The phenology's cue bits (`CUE_*`) are in
`config/params/meds_pft_params`, beside the masks they code.

## Two libraries

- **`meds_slow_kernels`** (`processes/`) links `meds_state`, `meds_config` and `meds_fast_kernels`,
  for the seam below. It builds on its own (`cmake --build … --target meds_slow_kernels`).
- **`meds_slow`** (`driver/`) links `meds_slow_kernels`. The stepper and `meds_init` link it.

## The laws decide, the update applies

`meds_slow_state` calls the laws of `meds_demography`; `meds_slow_update` applies the rates and
decisions it is given and never uses them. `meds_slow_dynamics` is where the two meet
(`src/README.md`, rule 4). The C API's `apply_rates` path
feeds externally computed rates through the same update routines. That path is the standing test
that the separation holds.

## The one fast/slow seam

The soil's decomposition chemistry is `meds_respiration`'s, in `fast_dynamics/processes/`: the
transfer matrix, the decay rates and Rh. The fast loop respires the day's frozen copy of the pools
with it at every sub-step. The daily step (`meds_soil_carbon`) debits the real pools through the same
matrix. So the day's sub-daily respiration equals the daily pool debit by construction, not
approximately. `advance_soil_carbon` checks the gap every step, and the slow ledger keeps the run's
worst value (`slow_ledger_soil_check`); it closes to machine precision.

This is the one place in MEDS where a kernel is called from both timescales. It is a documented seam,
not a folder. The canopy air's sub-daily CO₂ exchange is not here: it is
`fast_dynamics/processes/meds_canopy_air_space`.

## Science and tests

The order of the slow step and the calendar boundary is
[`order_of_processes.md`](../../docs/science/order_of_processes.md) §4–§5. Each process has its page:
[`plant_phenology.md`](../../docs/science/plant_phenology.md),
[`plant_carbon_allocation.md`](../../docs/science/plant_carbon_allocation.md),
[`plant_traits.md`](../../docs/science/plant_traits.md),
[`soil_carbon.md`](../../docs/science/soil_carbon.md) and
[`demography.md`](../../docs/science/demography.md). The slow ledger is
[`conservation.md`](../../docs/science/conservation.md), and the slow rows are
[`diagnostics.md`](../../docs/science/diagnostics.md) §4.

**Tests**, one per module. The processes: `test_plant_phenology` (with its C-API shim),
`test_plant_carbon_allocation`, `test_plant_traits`, `test_soil_carbon`, `test_demography`. The
driver: `test_slow_update` (the operators, and the restructuring through the ledger),
`test_slow_step` (the step through `meds_slow_state` and `meds_slow_dynamics`) and
`test_slow_ledger`. Through the C API: `test_c_api_demography`.

## Not here yet

- **Nitrogen.** The soil is carbon-only, and the nitrogen cycle is not built. It is planned after the
  move to ED-MORE.
- **DAMM.** There is no DAMM decomposition-moisture alternative. An unwired kernel for it is kept on
  branch `archive/damm-hr`.
- **Vertically resolved pools, nitrogen limitation of productivity, fire, and an explicit coarse
  woody debris pool** are all unbuilt. See [`docs/ROADMAP.md`](../../docs/ROADMAP.md) §3.
