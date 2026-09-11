# MEDS source-tree structure — reorganization plan

**Status (2026-09-11): migration steps 0-10 are ALL MERGED** — steps 1-6 (PR #125), steps 9/10 and
part of 8 (PR #126), steps 0 and 7 (`refactor/fast-loop-state-vector`), step 8 (PR #138, see §14),
and §10.2 slow-loop conservation (PRs #132-#137). What remains is the §10.3/§10.1 review remainder
plus two items found afterwards; **§15 is the phased plan for all of it**, and §15.0 strikes the four
§10 entries that have since closed. §13 records what step 0 and step 7 changed
about this plan, including a verification gap that affects how §8's acceptance criterion should
be read. Every implemented step was verified on BOTH back ends (ifx 38/38 +
nvfortran 38/38 multicore) and **byte-identical** in all 75 netCDF outputs of a 3-year, 4-thread
reference run -- including the two module SPLITS (steps 2 and 4), which the plan expected to be
only data-identical at round-off. See §11 for what the implementation changed about this plan.
Written 2026-09-08, revised after review PRs #119-#124, naming pass 2026-09-09 (folder names settled
with the author: `fast_dynamics/`/`slow_dynamics/`, `demography/` reinstated for the operator half of
core, `main/`, and the vital-rate laws move to demography -- see #2, #4, #5, #6, #14, §4, §5).

The branch this plan was blocked on (`refactor/review-prepass-signatures-assemblers`, PR #123)
has merged; PR #124 (`refactor/review-frozen-decomposition`: `column_frozen_t` decomposed into
nine content-named pieces, `integrator_opts_t`, `apply_process_mask`) is open and is the new
step 0. Every number and line reference below was re-verified against the PR #124 tip; the
diagnoses that those PRs changed are marked **[revised]**. §10 folds the still-open items of
`MEDS_CODE_REVIEW_2026-09-08.md` into the migration order so the two plans are one sequence.

Supersedes two placements locked by `MEDS_REORG_DESIGN.md` (v6):
- *"column types: `meds_column_state_types` STAYS in `shared/state/` — it's the BOUNDARY type"*
- *"each domain kernel library links `meds_shared` ONLY"*

Both were correct **given** the two-shared-library Python split. That split is now retired by
decision **#1** below, which removes the constraint that forced state into `shared`.

Scope: file/folder placement, library layering, and the Python package surface. **No physics,
no numerics, no algorithm changes.** Steps 1–6 are `git mv` + CMake target edits + `use`-line
renames; the only behavioural change in the whole plan is decision **#12** (allometry defaults),
which is a bug fix the reorg makes reachable.

---

## 1. Decisions locked

| # | Decision |
|---|----------|
| **#1** | **ONE Python shared library.** `libmeds_plant_c` + `libmeds_c` → **`libmeds.so`**. Per-subsystem independence moves from the *link* level to the *call* level (verbs, not libraries). The per-domain STATIC libs all stay. |
| **#2** | **State becomes a LAYER, in two halves.** `src/state/column/` (per-patch reservoirs + params + fusion blends — the boundary types kernels legitimately need) and `src/state/site/` (cohort SoA, patch CSR, lockstep, `site_t`, diag blocks — drivers only). Kernels link `state/column` and stay `site_t`-free. Two folders, not two file prefixes in one folder: every library in the tree is a per-folder GLOB, and the folder boundary *is* the link guard. **Naming (2026-09-09):** `column` is kept although a column is exactly one patch's vertical profile — `patch/` was rejected because `patch_block` (the CSR container of all patches) lives in the *site* half, so a `patch/` folder would point readers at the wrong half. The `state/column` README states the synonym once. **Grid-ready:** a future `state/grid/` (`grid_t` = array of `site_t` + geolocation + the met-forcing handle) is a third folder, not a rework: `met_driver_t` is already passed beside `site` (the stepper takes both), `meds_forcing` links shared + netCDF only and never `site_t`, so `grid_t` can own the forcing with no DAG cycle. Do **not** create an empty `state/grid/` now. |
| **#3** | **`shared/` dissolves** into `base/` + `functions/` + `config/` + `state/`. The word "shared" stops being a place where things go when no other place fits. |
| **#4** | **Timescale-first top level for processes:** `src/fast_dynamics/` and `src/slow_dynamics/`, each with domain subfolders and its own `driver/`. |
| **#5** | **`src/core/` splits along the seam it already has.** `meds_core_state_types` + `meds_core_diag_types` → `state/site/` (memory structure). `meds_core_state_update` + both `*_fusefiss` → `slow_dynamics/demography/` (they are daily-cadence apply-operators, not a foundation); library `meds_core` → **`meds_demography`** (the offload flags follow it, as they follow `meds_core` today). **The name goes back to `demography`** (reversing `MEDS_CORE_MODULE_REORG_DESIGN.md`'s rename, not its content): it is what ecologists call birth/death/growth bookkeeping and it is already the name of the Python module (`meds.demography`) and the C shim (`meds_demography_capi.f90`) that wrap exactly these files. **`meds_plant_vital_rates` moves in too**, as `slow_dynamics/demography/meds_demography_rates.f90` (module `meds_demography_rates`): its three `elemental pure` laws (carbon → diameter growth, Camac additive mortality hazard, reproduction carbon → recruits) are demographic rates, depend only on `meds_allometry` (functions layer) and PFT traits, and are called only by the slow driver — so `demography/` needs nothing from `plant/`. `slow_dynamics/plant/` is then pure physiology. **What this costs:** today `meds_core` links shared only, so an operator physically cannot call a rate law; after the move that guard is a signature convention — **rule 8** below — kept honest by review and by the Python `apply_rates` path, which feeds externally computed rates through the same operators every build. Optional mechanical guard: a one-line ctest grepping the operator files for `use meds_demography_rates`. |
| **#6** | **`src/plant/` splits by timescale**, exactly at file granularity (verified caller-by-caller, §5). `meds_plant_types` splits with it. The fast half lives in `fast_dynamics/plant/` (not `vegetation/`) so the two same-named folders make the timescale split of the plant library visible and match the Python names; `meds_plant_vital_rates` goes to `demography/` per #5, not to `slow_dynamics/plant/`. |
| **#7** | **Fix the `surface_state_t` homonym.** Two unrelated types, same name, both in scope in `src/driver/`. The biophysics one → `ground_optics_state_t`. |
| **#8** | **`*_types` modules hold only types they define.** No re-export shims: `grep 'type :: foo_t'` must return exactly one hit. |
| **#9** | **`toml` + `config_io` leave `src/io/`** for `src/config/`. They are configuration *input*, not I/O of model results. |
| **#10** | **Kernels take options leaves, never `meds_config_t`.** **[revised]** True for every kernel since PR #123: `leaf_gas_exchange_batch` takes `leaf_photo_table_t`, built once per run. What remains is (a) `meds_leaf_gas_exchange` importing five *selector constants* (`COLIM_*`, `SM_*`) from `meds_config` — they belong next to the table fields that hold them, in `meds_plant_types`; and (b) the two `cfg`-taking assemblers (`leaf_photo_params_for_pft`, `build_leaf_photo_table`) plus the single-leaf wrapper `leaf_gas_exchange(env, cfg, ...)` sitting in `meds_plant_interface`, i.e. logic in a facade (rule 6). They move to the fast driver, which is their only caller. |
| **#11** | **Python: declare `netcdf4` as a runtime dependency; do NOT vendor it into the wheel.** (§7.3) |
| **#12** | **Give `meds_allometry`'s nine coefficients their pan-tropical defaults as initializers.** The only behavioural change in this plan; a latent bug the merge makes reachable. (§7.4) |
| **#13** | **`test/` mirrors the new tree in the same commit as each move.** |
| **#14** | **Renames are a SEPARATE, optional, last step.** The *moves* fix cognition; the renames mostly do not. Two settled exceptions: (a) the `core` word retires *to* something concrete — operator modules `meds_core_{state_update,cohort_fusefiss,patch_fusefiss,interface}` → `meds_demography_*`, types modules → `meds_site_state_types` / `meds_site_diag_types` (churn from the current tree: state_types 14 src + 10 test `use` sites, diag_types 9 + 1, interface 8 + 3, the three operators 6 + 4 together); (b) `meds_plant_vital_rates` → `meds_demography_rates` is done **with** its move in step 6, since the new name is the reason for the move and it touches only `meds_vegetation_dynamics` and `meds_plant_interface`. Note the old, deleted `meds_demography_rates` held the *empirical* laws now in the Python example; the reused name holds the *carbon-driven* laws — a README line, so git archaeology does not mislead. |

**Load-bearing invariants preserved:** the acyclic library DAG; kernels never see `site_t`
(hence stay OpenMP-`target` device-eligible); the 12 ctests that link a kernel library alone
keep working; `meds_plant` / `meds_biophysics` / `meds_biogeochemistry` keep building as
standalone targets; byte-identical output at any thread count; every conservation ledger.

---

## 2. Diagnosis — why the current tree reads awkwardly

Numbers are the PR #124 tip (75 modules in 77 files, ~27.7k LOC).

### D1. One concept, three homes
A soil column's water is **declared** in `shared/state/meds_column_state_types.f90`
(`soil_column_t`), **owned** in `core/meds_core_state_types.f90` (`patch_block%soil_w`), and
**mirrored** in `driver/meds_fast_types.f90` (`column_state_t%theta`). Nothing in the tree says
which of the three is the state. This is the root of the "states defined across modules" feeling.

### D2. `shared/state` is a workaround, not a concept
It exists because CMake enforces *"biophysics/plant/biogeochemistry link `meds_shared` ONLY"*.
`patch_block` needs `soil_column_t`; the soil kernels need `soil_column_t`; `shared` was the only
mutually visible node. So a **state** module got filed under **shared**, which by its own
docstring is *"the foundation, NOT tied to any process."* The fix is not a better folder name —
it is admitting that **boundary state is a layer** and giving it one (decision #2).

### D3. A genuine homonym collision
`surface_state_t` is **defined twice**, with unrelated meanings:
- `biophysics/meds_biophysics_types.f90:118` — ground optical / skin state (albedo, emissivity, skin T)
- `driver/meds_fast_types.f90:306` — the prognostic CAS vector (enthalpy, shv, co2)

Both are in scope in the same subsystem: `meds_fast_dynamics` imports the biophysics one while
its siblings `meds_fast_ark` / `_rk45` / `_time_derivs` import the driver one. Highest
confusion-per-line in the tree, and the cheapest thing on this list to fix.

### D4. `*_types` modules are grab-bags
Nine `*_types` modules; the domain ones are largely re-export shims.
`meds_biophysics_types` (598 lines) re-exports `soil_column_t`, `soil_energy_column_t`,
`cas_state_t`, `soil_params_t`, `soil_thermal_params_t`, `n_soil_layer_max`, `curve_a`/`curve_n`
from three other modules. `meds_biogeochem_types` re-exports `soil_carbon_t` the same way. So
`use meds_biophysics_types` hands you state you cannot find by grepping that file for `type ::`.

`meds_column_state_types` itself mixes four concerns: prognostic reservoirs + fusion blends;
**parameter** bundles (`soil_params_t`, `build_soil_hydr_params`); init constants (`PSI_INIT`,
`LEAF_TEMP_INIT`); and one process **law** (`necromass_to_litter`, a litter-partition rule that
lives there only because `src/core` "cannot link biogeochemistry").

### D5. Timescale is invisible
MEDS *is* an operator-split fast/slow model — that is its dominant organizing axis — but the tree
is organized by neither timescale nor a clean domain split. `src/driver/` is 7.8k lines of which
**6.4k is fast-loop machinery (eleven modules, now including `meds_fast_prepass` and
`meds_column_state_ops` from the review PRs), 0.9k is the slow loop (three modules) and 0.5k is
the application (`meds_main`, `meds_stepper`)**, while the fast *processes* it drives live three
folders away in `src/biophysics/`. Nothing says `meds_leaf_gas_exchange` is sub-daily and
`meds_plant_carbon_allocation` is daily. The review PRs made the fast/slow boundary sharper (the
prepass, the frozen record and the integrator options are all fast-only records) without moving
a single file — the tree still hides it.

### D6. `meds_config` is a god object
628 lines, one flat `meds_config_t`, `use`d by **22** modules. It aggregates `soil_opts_t` /
`energy_opts_t` / `snow_params_t` / `aero_cfg_t` / `decomp_opts_t`, which are then **re-copied**
into `column_config_t` in `meds_fast_types`: two config layers. **[revised]** The leaf kernel no
longer bypasses them — PR #123 gave it `leaf_photo_table_t`, built once per run into
`column_config_t%leaf_photo`, and PR #124 did the same for the integrator knobs
(`integrator_opts_t`, built by `build_integrator_opts`). So `column_config_t` is now the
fast loop's *complete* options record and `meds_config_t` is read by the fast loop only at
`build_fast_context` time. The kernel's remaining `meds_config` import is five selector
constants (see #10). The god-object complaint is therefore about the *slow* loop, io and init,
where `meds_config_t` is still passed whole.

### D7. `src/io/` is one folder producing four libraries
`meds_config_io`, `meds_netcdf_c`, `meds_io_prep`, `meds_io_stream` — and two of them are not
I/O of model results at all (TOML config *input*; site-level diagnostic *reductions*). One is a
hand-maintained explicit file list while its neighbours `GLOB`, which is an easy trap for a
newly added diagnostic module.

### D8. Inconsistent facades
Three `*_interface` modules, three contracts: `meds_core_interface` (33 lines, pure re-export),
`meds_biophysics_interface` (49, pure re-export), `meds_plant_interface` (209, re-export **plus
logic** — four procedures, two of which take `meds_config_t`). So there are always two legal
ways to `use` anything, and callers mix them — `meds_fast_ark` imports `solve_plant_water_batch`
from `meds_plant_hydraulics` (line 25) *and* from `meds_plant_interface` (line 64) in the same
module. The core facade is bypassed by 12 src modules and 10 tests (item 4 #13 of the review),
mostly for `site_t`, the allocators, `rebuild_csr`, the `DMAX_PSI_LEAF_*` sentinels and every
`meds_core_diag_types` constant — i.e. for the *state* half of core, which a verb facade was never
the right shape for. Decision #5 dissolves that complaint rather than fixing the facade (§10).

---

## 3. The lever: one Python library dissolves the constraint

Today:

| artifact | links | size | netCDF/HDF5 |
|---|---|---|---|
| `libmeds_plant_c.so` | `meds_plant` only | 59 KB | no |
| `libmeds_c.so` | `meds_aux` + `meds_io_stream` | 728 KB | yes |

The two-`.so` split is the **only** reason for the "each domain library links `shared` only"
rule, and that rule is what pushed column state into `shared`. One library over per-subsystem
*verbs* replaces the flat rule with a layered DAG:

```
base <- functions <- config <- state/column <- kernels <- state/site <- demography <- drivers <- main | capi
                                                       \      forcing, io          /
```

Splitting state in two (decision #2) is what preserves the property worth keeping: kernels link
`state/column`, never see `site_t`, stay device-eligible, and the 12 kernel-only ctests keep
working. `shared/state` was groping for exactly this distinction; it just had no layer to put it in.

---

## 4. Target tree

```
src/
├── base/        kinds, constants, numerics, time, budget_check          [libmeds_base]
├── functions/   allometry, therm_lib, hydr_lib, optics_lib, temp_response
├── config/      pft_params, *_opts, meds_config, toml, config_io        <- toml/config_io ex-io/
├── state/
│   ├── column/  soil_water / soil_energy / snow / cas / soil_carbon reservoirs,
│   │            blends, params + assemblers  (one patch, seen vertically) <- ex-shared/state, SPLIT
│   └── site/    cohort SoA, patch CSR + reservoir ownership, lockstep,
│                site_t, diag blocks, fast/slow slices                  <- ex-core (types half)
│   (grid/      future: grid_t = site_t(:) + geolocation + met handle — NOT created now)
├── fast_dynamics/                                            # sub-daily
│   ├── canopy/     canopy_radiation, canopy_aerodynamics, cas_biophysics   (the medium)
│   ├── plant/      vegetation_biophysics, leaf_gas_exchange, plant_hydraulics,
│   │               plant_maintenance_respiration   (+ leaf / hydro arg types) (the organisms)
│   ├── soil/       soil_water, soil_energy, ground_biophysics (incl. snow kernels)
│   ├── numerics/   fast_types, column_state_ops, fast_control, fast_time_derivs,
│   │               fast_ark, fast_rk45, fast_rk4_oracle, fast_snow, fast_step
│   └── driver/     fast_prepass, fast_dynamics
├── slow_dynamics/                                            # daily / annual
│   ├── plant/      phenology, carbon_allocation, trait_dynamics (+ pheno arg types)
│   ├── soil/       soil_biogeochem, biogeochem_types, necromass_to_litter
│   ├── demography/ demography_rates (ex-plant_vital_rates), state_update,
│   │               cohort_fusefiss, patch_fusefiss             <- ex-core (operator half) [libmeds_demography]
│   └── driver/     vegetation_dynamics, biogeochem_dynamics, slow_dynamics
├── forcing/     (unchanged)
├── io/          netcdf_c, diagnostic_kernels, diagnostic_reduce, output_*, meds_io
├── init/        (unchanged)
├── capi/        meds_capi_{leaf,hydraulics,phenology,demography,fast,site}.f90
└── main/        meds_stepper, meds_main
```

Folder-name rationale (settled 2026-09-09, so it is not re-opened per file):

- **`fast_dynamics/` / `slow_dynamics/`**, not `fast/`/`slow/` (adjective without a noun) and not
  `*_proc/` (each folder also holds numerics and a driver; "proc" reads as "procedure" in Fortran).
  The folder's namesake module (`meds_fast_dynamics`, `meds_slow_dynamics`) is its entry point,
  the same pattern as `plant/meds_plant_interface` today.
- **`canopy/` absorbs the canopy-air-space module** (95 lines, two procedures; aerodynamics computes
  the conductances it consumes). Reading: `canopy/` is the medium, `plant/` is the organisms.
- **`soil/` on both timescales**, not `ground/`: snow rides inside. The word "ground" survives where
  it means the soil-or-snow *surface* (`meds_ground_biophysics`, `ground_optics_state_t` from #7).
  The fast→slow Rh edge (S1) then reads `fast_dynamics/soil` → `slow_dynamics/soil`.
- **`demography/`**, not `structure/` (which was shorthand for "stand structure" and read as nothing).
  Its README states the two-part rule (rule 8) in the reader's vocabulary: the rate laws live in
  `meds_demography_rates`; the operators take rate *arrays* and never import them.
- **`main/`**, not `app/`: the program entry plus the one-step scheduler, and later the grid loop.

Three placements the review PRs settled, so they are not re-litigated by the moves:

- **`fast_dynamics/numerics/` is the home of the fast loop's work records.** `meds_fast_types` (761
  lines, 24 types) now holds `column_config_t` with its `leaf_photo` and `integrator` leaves,
  `column_frozen_t` and its nine content-named pieces (`cas_boundary_t`, `tissue_coefficients_t`,
  `canopy_film_capacity_t`, `ground_boundary_t`, `soil_hydrology_t`, `root_zone_t`,
  `plant_water_t`, `column_params_t`, plus `snow_stage_t` from `meds_fast_snow`), the state
  vector and its tendencies, the stage/column boundary-flux ledgers, and `column_cohort_t`
  (which §10 deletes). The review decided this file is **not** split and the pieces carry no
  `frozen` in their names — "frozen" is a lifetime statement made by the container and by
  `intent(in)`. The moves keep that: one file, one folder.
- **`patch_biophys_t` is a driver record filed in a kernel library.** It is the per-patch
  reservoir gather (`cas`, `soil_e`, `soil_w`, `snow`, the frozen soil-carbon snapshot,
  `adapt_dt_last`) used only by `src/driver` and the column tests, yet it is defined in
  `meds_biophysics_types`. It moves to `fast_dynamics/numerics/meds_fast_types` with the rest of the
  fast-loop records (rule 1: a type lives with whoever mutates it).
- **Assemblers from `meds_config_t` are driver code.** `build_leaf_photo_table`,
  `build_integrator_opts`, `build_tol_set`, `build_error_control` read the whole config and
  produce an options leaf. They belong in `fast_dynamics/driver` (their only caller is
  `build_fast_context`), not in a kernel facade (`meds_plant_interface`) or a numerics module
  (`meds_fast_control`, which then keeps only `state_wrms_grouped`/`step_control_factor`).

Two notes on decision #5: after `src/core/` splits, the word "core" stops carrying information
and the library name can simply retire (see #14 — not urgent). And the split is along a seam the
code already documents: *"the engine NEVER computes a rate — it APPLIES"*. Applying on a daily
cadence is a slow-loop process, not a foundation.

---

## 5. The three straddlers

Be explicit about these rather than letting them blur the tree.

**Verified caller map for `src/plant/`** — the timescale split is exact at file granularity, and
the kernels never call each other across it. Their *only* coupling is `meds_plant_types`, which
splits just as cleanly:

| module | called from | tier |
|---|---|---|
| `meds_leaf_gas_exchange` | `meds_fast_prepass` (`canopy_leaf_gas_exchange`) | fast |
| `meds_plant_hydraulics` | `meds_fast_ark` (`solve_plant_water_batch`) | fast |
| `meds_plant_respiration` | `meds_fast_prepass` (`canopy_maintenance_respiration`) | fast |
| `meds_phenology` | `meds_vegetation_dynamics` | slow |
| `meds_plant_carbon_allocation` (incl. `growth_respiration`) | `meds_vegetation_dynamics` | slow |
| `meds_plant_trait_dynamics` | `meds_vegetation_dynamics` | slow |
| `meds_plant_vital_rates` | `meds_vegetation_dynamics` | slow → `demography/` as `meds_demography_rates` (#5) |

`meds_plant_types` → `leaf_*` (including the PR #123 `leaf_photo_table_t` and, per #10, the
`COLIM_*`/`SM_*` selector constants) + `hydro_*` + `wood_params_t` + `root_params_t` to **fast**
(`wood_params_t`/`root_params_t` are read only by maintenance respiration), `pheno_*` to **slow**.
`meds_plant_capi` imports leaf, temp-response and phenology symbols, so after the split it links
both halves — fine for a capi shim, and one more reason it moves to `src/capi/` (§7.6).

**S1. `meds_soil_biogeochem` — do NOT split.** `soil_carbon_step` is daily;
`heterotrophic_respiration_matrix` is called every `dt_fast`. The fast Rh must respire the *same*
matrix the slow step debits, and that co-location is precisely what makes `rh_seam_gap` close to
machine precision. Keep the module whole in `slow_dynamics/soil/` and document the one deliberate
fast→slow edge.

**S2. `meds_plant_respiration` — nothing to split.** Maintenance is fast; growth respiration
already lives in `carbon_allocation`. Goes to `fast_dynamics/plant/` under a name that says
"maintenance".

**S2b. `meds_plant_vital_rates` — move to `demography/`, not `slow_dynamics/plant/`.** Its three laws
consume the wood NPP that `carbon_allocation` produces, so producer and consumer end up in different
folders. That is already the case in substance — `meds_vegetation_dynamics` is the only bridge between
them — and it is the honest seam: allocation is physiology, the hazard it drives is demography. The
module imports only `meds_kinds` and `meds_allometry`, so `demography/` still links nothing from
`plant/`.

**S3. `necromass_to_litter` — move it.** A litter-partition law living in the column *state*
module only because `src/core` cannot link biogeochemistry. Under the new DAG that constraint is
gone; it belongs in `slow_dynamics/soil/`.

---

## 6. Placement rules

So the next new file is not a judgement call:

1. A derived type lives with **whoever mutates it**. If two subsystems mutate it, it is boundary
   state → `state/column/`.
2. `*_types` modules hold **only types they define**. No re-export shims (decision #8).
3. **Parameters ≠ state.** `soil_params_t` + `build_soil_hydr_params` are derived-once-per-column
   parameters, not reservoirs — a separate module from the reservoirs they describe.
4. A kernel goes where its **caller's timescale** is. A kernel called from both timescales is a
   documented seam (§5), never a folder.
5. Kernels take **options leaves** (`soil_opts_t`, `aero_cfg_t`), never `meds_config_t`
   (decision #10).
6. One facade per library, and it is **either** pure re-export **or** it has logic — not both.
7. C-API verbs take **explicit parameter structs**, never a config handle, except for the
   coupled-model verbs that legitimately own a `Config`. This is rule 5 at the API level, and it
   is what keeps `meds.leaf` runnable with no TOML (§7.2).
8. **Demography operators take rate ARRAYS as arguments and never `use` a rate module.** The rate
   laws (`meds_demography_rates`) and the operators (`*_state_update`, `*_fusefiss`) share a folder
   and a library, but the only place a rate meets its application is the slow driver. This is the
   "engine never computes a rate — it APPLIES" invariant restated for the new tree; the Python
   `apply_rates` path (externally computed rates through the same operators) is its standing test.

---

## 7. Python package

### 7.1 Compiled artifact

One `libmeds.so` (~1 MB), linking netCDF + HDF5 — unavoidable, since the forcing reader and
restart writer are in it. All the per-domain STATIC CMake targets stay: fast incremental builds,
and 12 of the 40 ctests keep linking a kernel library alone. Only the *shared* library collapses.

### 7.2 Package tree

```
python/
├── pyproject.toml                  # scikit-build-core: compiles + bundles the .so into a wheel
└── meds/
    ├── __init__.py                 # cheap; NO dlopen (keep today's property)
    ├── _core/                      # the ONLY ctypes code in the package
    │   ├── _lib.py                 #   lazy dlopen: MEDS_LIB env -> bundled -> build*/
    │   ├── _sig.py                 #   one signature declaration per bind(c) verb
    │   ├── _structs.py             #   bind(c) struct mirrors
    │   └── libmeds.so              #   bundled by the wheel
    │
    ├── leaf.py                     # photosynthesis + stomata + Ci solver   (was meds.plant.leaf)
    ├── hydraulics.py               # plant water transport
    ├── radiation.py                # two-stream canopy RT
    ├── aerodynamics.py
    ├── soil.py                     # Richards / soil thermal
    ├── phenology.py                #                                       (was meds.plant.pheno)
    ├── allocation.py               # carbon allocation + growth respiration
    ├── soil_carbon.py              # CENTURY decomposition
    ├── allometry.py                # size <-> height <-> AGB  (gated on #12)
    ├── demography.py               # Site / apply_rates / fuse-fiss + the vital-rate laws  (1:1 with slow_dynamics/demography/)
    │
    ├── fast.py                     # driver verb: advance one dt_fast
    ├── slow.py                     # driver verb: advance one dt_slow
    ├── site.py                     # the coupled model (Config + Site handles)
    └── output.py                   # read MEDS netCDF output
```

The naming split mirrors decision #4 where it carries information and ignores it where it does not:

- **Kernels are named by domain** (`meds.leaf`, `meds.soil`). A kernel's timescale is a property
  of the kernel, not a choice the caller makes; a user asking for photosynthesis should not need
  to know it is a sub-daily-tier module.
- **Drivers are named by timescale** (`meds.fast`, `meds.slow`) — that *is* the choice the caller
  makes.

**Independence is preserved, and rests where it already actually rests — the call signature:**

```fortran
meds_leaf_solve(env_c, p_c, sm, tresp, colim, use_boundary_layer, flux_c)  ! no site, no config, no file
meds_config_load(path, path_len) -> handle                                 ! demography needs a TOML
```

`meds_leaf_solve` takes all 35 leaf parameters as an explicit struct — which is what makes
`examples/example_leaf_gas_exchange/reproduce_slot2017.py` possible ("the model lives in Fortran,
but no parameters are hard-coded there"). Merging the libraries does not touch that. Both existing
examples change **only their import line**:

```python
from meds import leaf           # was: from meds.plant import leaf
from meds import demography     # was: from meds.demography import Site
```

**New capability unlocked.** There is no C-API for the fast loop today — "run submodules
independently" currently means `{photosynthesis, phenology}` plus `{demography slow loop}` and
nothing else. A `meds_capi_fast` verb driving `column_fast_step` on a column handle would expose
the sub-daily integrator to Python, which is where most recent diagnostic work lives
(`scripts/numerics_sweep.py`, the `runs/ithaca_ark30` probes). Merging the libraries is a
**precondition**: that verb needs both the kernels *and* `site_t`, so under today's split it has
no library to live in.

### 7.3 Deployment cost (the real one)

The `.so` size is noise. What matters is the dependency closure, and netCDF-4 drags in a lot via
its remote-access support:

```
libnetcdf.so.19   1.7 MB
libhdf5.so.310    4.6 MB
+ libcurl, libssl, libcrypto, libxml2, libzip, libjpeg, libsz, libnghttp2, libssh2 ...
                  ~ 67 MB transitive closure (measured, this machine)
```

- **Declare `netcdf4` as a runtime dependency, do not vendor** → cost ≈ **zero** for the realistic
  audience, who already have netCDF installed to read the output.
  `[project.optional-dependencies].viz` already lists `netCDF4`.
- Vendoring into a self-contained wheel (`auditwheel repair`) → ~67 MB wheel. Don't.

And the cost is smaller than it first looks: **the project already cannot be configured without
netCDF.** `CMakeLists.txt` calls `find_package(netCDF CONFIG REQUIRED)` unconditionally, before
any target is declared, so a collaborator who wants only the leaf model already needs netCDF to
build `meds_plant_c` today. The 59 KB `.so` was netCDF-free *at link time*, never netCDF-free to
*build*. Merging adds a dependency nobody is currently escaping.

### 7.4 The cost that is not about import — decision #12

`shared/functions/meds_allometry.f90` declares nine `protected` module variables with **no
initializers**, written only by `set_allometry`, which is called from exactly one place
(`derive_parameters` in `meds_config`). Its own docstring says the pan-tropical values are
*"the canonical defaults shipped in the config, not baked into the source."*

Today that is unreachable from Python: `libmeds_plant_c` has no config loader and its capi does
not expose allometry. A merged package that exposes `meds.allometry` — which is wanted, since
`examples/example_demography/empirical_laws.py` currently reimplements the pan-tropical allometry
in numpy purely to have it — makes it a live footgun:

```python
import meds
meds.allometry.dbh_to_height(30.0, hgt_max=35.0)   # reads uninitialized module state
```

Zeros in practice on Linux (`.bss`), formally undefined, silently wrong either way. Also note
that one `.so` means one process now shares those coefficients between `meds.demography` and
`meds.allometry`, where two `.so`s gave two independent copies — arguably better (one process,
one allometry, no silent divergence), but a real semantic change.

Two fixes, in order:
1. **Now:** give the nine coefficients their pan-tropical values as initializers, so unconfigured
   use is *correct* rather than undefined. Config still overrides via `set_allometry`.
2. **Eventually:** make them an argument bundle instead of module state. The path exists — the
   offloaded growth kernel already takes them as scalar arguments ("it cannot read host module
   state on the device").

This is the last significant instance in `shared/` of the module-level-mutable-state pattern that
made the `-auto` / static-locals trap so expensive to find.

### 7.5 What is genuinely free

Checked, and there is no cost here: **no numerical change** (same objects, same kernels); **no GPU
change** (`MEDS_BUILD_PYLIB=ON` already forces `CMAKE_POSITION_INDEPENDENT_CODE` globally); **no
build-time change** (the static libs stay); the handle registry (`g_site` / `g_cfg` in
`meds_demography_capi.f90`) is already a single `save` array inside one `.so`.

**The whole bill: one deployment decision (#11) and one code fix (#12).**

### 7.6 Python-side work items

1. Swap the build backend `setuptools` → `scikit-build-core` in `python/pyproject.toml`, so
   `pip install python/` compiles and bundles `libmeds.so`. The existing comment there already
   flags this as the plan; one `.so` is what makes it tractable.
2. Move `src/plant/meds_plant_capi.f90` into `src/capi/`; delete the
   `GLOB src/plant/*_capi.f90` special case in CMake.
3. One capi file per subsystem, mirroring the Fortran tree: `meds_capi_leaf`,
   `meds_capi_hydraulics`, `meds_capi_phenology`, `meds_capi_demography`, `meds_capi_fast`,
   `meds_capi_site`.
4. **Keep compiling every `*_capi.f90` into a ctest target** (`test_capi_leaf`,
   `test_capi_demography`, ...). This is the #95 → #100 lesson already recorded in
   `CMakeLists.txt`: a component inserted mid-type in `leaf_photo_params_t` broke the C API while
   the whole suite stayed green, because the shim was compiled only by the optional pylib target.
   One `.so` means one ABI and one place for that to break — keep it a *build* failure.

---

## 8. Migration order

Each step leaves the suite green and touches no physics. Because Fortran `use` is by module name
and all `.mod` files land in one directory, steps 3/5/6 are **pure `git mv` + CMake target edits
with zero source changes** — the cheapest part of this, and the part that buys the most.

| # | Step | Churn | Kills |
|---|---|---|---|
| **0** **DONE** | Merge PR #124 (frozen decomposition, `integrator_opts_t`, `apply_process_mask`); then do review step 5 — delete `column_cohort_t`, fast/slow slices, PFT geometry params — **on the current tree, before any file moves** (§10.1) | content edit to core/io/driver | the last content change to `src/core` files, so every later move commit is a pure `git mv` |
| **1** **DONE** | Rename biophysics `surface_state_t` → `ground_optics_state_t` | ~6 sites | **D3** — highest ratio on this list |
| **2** **DONE** | Split `meds_column_state_types` → reservoirs / params / init-constants; move `necromass_to_litter` to `slow_dynamics/soil/` | ~27 `use` sites, mechanical | **D4** |
| **3** **DONE** | Introduce `src/state/{column,site}` as real layers; delete `shared/state`; relink kernels to `state/column` | CMake + moves | **D1 + D2** |
| **4** **DONE** | Drop the re-export blocks from `meds_biophysics_types` / `meds_biogeochem_types` | ~15 `use` lines | the "invisible state" half of **D4** |
| **5** **DONE** | Create `src/config/` (absorb `toml` + `config_io`); `src/io/` becomes netCDF + diagnostics only | moves only | **D7** |
| **6** **DONE** | Create `src/fast_dynamics/`, `src/slow_dynamics/`, `src/main/`; split `src/plant/` (§5) and `src/core/` (#5); `meds_core` → `meds_demography` target; the one move-with-rename: `plant/meds_plant_vital_rates.f90` → `slow_dynamics/demography/meds_demography_rates.f90` (module renamed, 2 `use` sites) | ~40 file moves, 2 `use`-line edits | **D5** |
| **7** **DONE** | Continue splitting `meds_fast_ark` (1581 lines). PR #120 already moved the state algebra to `meds_column_state_ops`; what non-ARK code still imports from it is exactly three symbols: `build_column_frozen` (RK45), `column_be_stage` and `advance_water_mass_full` (oracle). Move the pre-pass builder to `meds_fast_prepass` and the BE-stage/Newton machinery to its own module; `meds_fast_ark` keeps the tableau and the march. Fold in the deferred review item "pass `column_params_t` instead of copying it into the frozen record" — this is the one step that touches every march signature anyway (§10.3) | procedure moves between modules → data-identity criterion, not byte-identity | the `rk45 → ark` and `oracle → ark` edges, which are not about ARK |
| **8** **DONE** | Python: decisions #1, #11, #12 + §7.6 (§14) | small | the two real costs |
| **9** **DONE** | Facade normalization, now concrete (§10.4): `state/site` has **no** facade — drivers, io and tests import `site_t`, the allocators and the diag blocks from the state module directly, which is what 22 of them already do; `meds_core_interface` becomes `meds_demography_interface`, re-exporting the `slow_dynamics/demography` verbs only, or is deleted. `meds_plant_interface` loses its logic (§4 note 3) and becomes pure re-export like the other two. *Optional:* config decomposition for the slow loop/io (#10, D6), the `meds_core_*` → `meds_demography_*` / `meds_site_*` renames (#14) | larger | **D6, D8** |
| **10** **DONE** | Renames, last and byte-identical (§10.5): the review's remaining field and routine renames merge into decision #14's list | `sed -I -w` per group | the names that lie |

Steps 1–5 are worth doing **regardless** of whether the fast/slow tree (#4) is adopted.

**Acceptance criterion per step, learned in the review PRs.** A pure `git mv` of a file keeps the
module name, so the objects are unchanged and the July example (`examples/example_biophysics`,
`bi.toml`, 4 threads) must be **byte-identical** to main in all 31 output files — steps 3, 5, 6.
Splitting a module or moving a procedure between modules changes ifx code generation; the outputs
are then **data-identical at round-off** (hour-1 differences ~1e-11 relative, ~1e-6 after a month)
— steps 2, 7 and the §10 content edits that move procedures. Renames are byte-identical. Compare
netCDF *data*, not bytes, whenever a module boundary moved; headers carry a timestamp. And build
the NVHPC back end at every step that moves a procedure: nvfortran 25.11 rejected
`acc = column_bflux_t()` on an `intent(out)` dummy after PR #120 moved `bflux_zero`, a line it
had compiled for months in its old module.

### Non-goals / cautions

- **`test/` moves with each step** (#13). 38 flat test programs (21 `add_test` statements, some
  in loops) already do not say which subsystem they cover (`test_column_energy` vs
  `test_surface_energy` vs `test_snow`); if the mapping is deferred it rots immediately. The
  review's item 5 #17 rides along: 21 of them define their own `check`/`check_true`, three
  hand-build the same column fixture, and `meds_test_support` (`check`, `check_close`, `banner`,
  `build_test_config`) is used by none of the numerics tests. Consolidate into `meds_test_support`
  in the same commit that moves each test — the move already touches every line of the file list.
  Note `meds_testsupport` links `meds_core` today; after step 3 it links `state/site`.
- **Do not rename for its own sake** (#14).
- **Do not split `meds_soil_biogeochem`** (S1) — it would break the `rh_seam_gap` co-location.
- **Do not re-create a kernels-only `.so`** to dodge §7.3. A `-DMEDS_PYLIB_KERNELS_ONLY` variant
  is the escape hatch if deployment ever really bites, but adding it back re-creates exactly the
  split decision #1 removes.
- Verification per step is the existing suite on **both** back ends (ifx 38/38 + nvfortran
  multicore), per the CLAUDE.md portability rule — a green ifx run is not sufficient.

---

## 9. Relation to other plans

- `MEDS_REORG_DESIGN.md` (v6) — supersedes its two placement decisions quoted in the header;
  everything else in v6 stands (allometry → `shared/functions`, empirical laws → Python, carbon
  rates → `meds_plant_vital_rates` — now `meds_demography_rates` under `demography/`, #5 — `libmeds_c` opt-in).
- `MEDS_CORE_MODULE_REORG_DESIGN.md` — decision #5 splits the 4-file core it established; the
  4-file *content* is unchanged, only its two halves move to different layers. Its `demography` →
  `core` rename is **reversed** for the operator half (the author's 2026-09-09 decision); the state
  half takes the `site` name, so `core` means nothing afterwards.
- `MEDS_CAPI_REORG_DESIGN.md`, `MEDS_DRIVER_REORG_DESIGN.md` — §7.6 and step 7 extend these.
- `MEDS_CODE_REVIEW_2026-09-08.md` — PRs #119–#124 closed its items 1A, 2 and most of 4–6; §10
  below absorbs what it left open, so that document's "order of the remaining work" is superseded
  by this plan's migration table from step 0 onward.

---

## 10. Folding the 2026-09 review remainder into this plan

The review closed the fast-loop ledgers (energy to machine precision over a seasonal cycle), the
RK45/threading defects, film-water conservation across split/fuse/disturbance, the dead
selectors, the duplicated post-march blocks, the `column_prepass` fusion, the frozen-record
grab-bag and the integrator knobs. What it left open falls into five groups. Each is placed at
the migration step where it costs least, with the rule: **content edits before the file that
holds them moves; renames after everything.**

### 10.1 Review step 5 — `column_cohort_t` out, fast/slow slices in → **migration step 0**

This is the last content change to `src/core` and `src/io` files and the one that touches the
lockstep lists, so it goes **before** `src/core` splits (step 3/6). Done after the split it would
be the same edit spread over two folders and two commits' worth of `use` lines.

- **Delete `column_cohort_t`** (`meds_fast_types:176`, used by 6 src + 3 test files). It is a
  read-only view of demographic inputs plus derived geometry, not state; the per-thread scratch
  exists only because the gather writes into it. The hand-built test views are allometrically
  inconsistent and never set `bwood` (the wood heat capacity in those tests runs on unset memory
  floored to the minimum). The driver reads the patch's **contiguous CSR cohort section** of the
  cohort block directly — drivers may see `site_t`; kernels in `fast_dynamics/` still never do.
- **Derived geometry becomes cohort-block fields** refreshed after growth, next to the cached
  `height`/`basal_area`/`agb`/`leaf_area` that already live there: `lai`, `wai`, sapwood carbon
  and area, total wood carbon. **The three hard-coded constants** (`meds_fast_types:711-712`:
  leaf width 0.04 m, branch diameter 0.02 m, crown fraction 1.0) become **PFT parameters** in
  `config/pft_params`, gathered per cohort like `p_hgt_max`.
- **Fast/slow slices.** The 12 cohort + 8 patch fields the fast loop writes back are today
  repeated in four lockstep lists in core plus the restart writer (47 mentions of the two probe
  fields alone across those five files). They become `cohort_fast_slice_t` / `patch_fast_slice_t`
  components of the blocks in `state/site`, with **one** per-field policy table (extensive vs
  intensive vs ground-referenced) that fusion, fission, disturbance copy and the restart writer
  all consume. That table is where item 1B #1/#2 (film water is ground-referenced) was fixed by
  hand in PR #119; making it declarative is what stops the next such bug.
- **Seed and clamp leave the fast gather.** The lazy PSI_INIT tissue-water seed and
  `clamp_water_to_capacity` (item 1A #7, unbooked mass edit) move to a slow-loop
  `reconcile_tissue_water_capacity` in `slow_dynamics/driver`, run after allocation changes the capacity,
  with its discarded/added mass declared to the ledger of 10.2. This also retires the
  `leaf_water_mass <= 0` sentinel that PR #119 only half-fixed (item 1B #3).
- **Extensibility check** (the author's stated reason for the decision): after this, a new
  per-cohort *input* to column physics is one field on the cohort block, and a new *prognostic*
  cohort field is one slice component plus its policy row. The remaining ~12 enumeration sites
  for `column_state_t` (combinators, pack/unpack, mask, error norm) are the fast integrator's own
  problem and are addressed in 10.3.
- Verification: data-identical for the production path (the geometry values are unchanged, only
  their home moves); the three column tests change numerically because their fixtures become
  allometrically consistent — record the new golden values with the reason.

### 10.2 Slow-loop conservation (review item 1B #4–#10) → **after migration step 6, in `slow_dynamics/`**

Physics, so outside this plan's "no numerics" scope — but the *placement* is this plan's business, and
the order matters: these land after step 6 so they are written once into their final home instead of
being moved a week later. **Surveyed against the code 2026-09-09**; what that survey found is recorded
below, because several of the review's original items turned out to be larger, smaller or differently
shaped than the review's one-line description of them.

#### 10.2.0 Why none of this is currently visible

`budget_t` lives in `util/meds_budget_check`, **not** `base/` (an earlier draft of this section said
`base/`; the ledger reuses it from `util/`). It accumulates **per-fast-step flux residuals**, which
`fast_dynamics` merges up into the run-level `run_energy_budget` / `run_water_budget`. Nothing anywhere
compares a **store** across the slow step. A discontinuous jump between the end of one fast window and
the start of the next is therefore invisible **by construction**, not by oversight — which is why every
item below has survived a full review and a conservation-focused PR.

Two consequences for the design:

- The slow ledger is a **snapshot** ledger — store before, declared boundary terms, store after — not a
  flux accumulator. It is the fast loop's ledger turned inside out.
- It must be **site-level and area-weighted** (`Σ_p area_p × store_p`). Patch identity does not survive
  the step: `apply_patch_disturbance` creates a patch, `fuse_2_patches` destroys one, and
  `terminate_patches` renormalizes every area. A per-patch ledger has nothing to compare against.

#### 10.2.1 The stores

`meds_fast_ark`'s `whole_water` / `whole_energy` check already names the canonical list; the slow ledger
spans the same stores so the two tiers cannot disagree about what exists, and adds carbon, which the fast
whole-column ledger does not track.

| currency | stores (each area-weighted to the site) |
|---|---|
| water  | soil column; pond `w_surface`; snow `swe`; CAS vapour `cas_mass_capacity·shv`; tissue water (`leaf_water_mass`+`wood_water_mass`, ×nplant); interception films (already ground-referenced — **not** ×nplant) |
| energy | soil energy; pond enthalpy; snow enthalpy; CAS `cas_mass_capacity·enth`; tissue heat |
| carbon | live pools (leaf/fineroot/wood/nonstructural, ×nplant); the 7 CENTURY pools; CAS CO2 (`cas_molar_capacity·can_co2`); `recruit_pool` as carbon-in-transit |

#### 10.2.2 Carbon terms, largest first

1. **Growth respiration is destroyed, not exhaled** (this is 1B #5, and it is much the largest item in
   this section). `plant_carbon_allocation` charges `growth_resp` against the plant's carbon and its own
   header states the identity it must satisfy; nothing consumes it — `growth_resp` reaches only the
   `CS_GROWTH_RESP` diagnostic. Meanwhile `nee_biotic = ra_leaf + ra_stem + ra_root + rh − gpp` carries
   **maintenance** respiration only. With `growth_resp_factor = 0.3` that is `0.3/1.3` = **23 % of every
   unit of carbon entering growth**, debited from the plant and never reaching the atmosphere. Same
   defect class as the deleted soil-carbon fallback (#128), opposite sign: the CAS runs too low, so
   photosynthesis is under-fertilized and the site reports as a larger sink than its own carbon flow
   implies.
2. **The starvation `deficit` is the mirror image.** When storage cannot cover maintenance, the fast loop
   has *already* respired that carbon into the CAS, but no pool is debited. Carbon from nothing, per
   cohort, on every step a cohort starves.
3. **Recruitment endows ~5–6× what it debits** (1B #4). The debit is `min_cohort_carbon` =
   `dbh_to_agb(...)`, **AGB only**; `init_cohort` → `set_cohort_size` endows `wood_carbon = agb /
   aboveground_frac` **plus** leaf, fineroot and storage. Measured on the shipped 3-PFT table
   (`min_cohort_height` 2 m): total/agb = 6.43, 5.52, 5.18. `repro_carbon_efficiency = 1e-3` then
   over-corrects, so **net ~99.4 % of reproduction carbon vanishes**. Both halves need declaring: the
   establishment loss is physically real but it is necromass, not nothing, and the endowment/debit ratio
   should be 1.
4. **Silent floors.** `max(pool, 0)` appears in both `update_cohort_derivatives` and
   `update_cohort_states_kernel`; `nplant` is floored at `negligible_nplant`. Each is creation. The
   nplant floor also **double-counts**: `accumulate_mortality_litter` values litter on the *unfloored*
   `died_nplant`, so litter is produced for individuals the floor then refuses to kill.
5. **Litter is gated on `soil_carbon_on`** in four places (turnover, continuous mortality, cull,
   disturbance kill). With soil carbon off — the default — all necromass vanishes. Defensible as an
   *export*, but after #128 it is declared as one rather than left implicit.
6. **Operator-split offset.** Mortality litter is valued on pre-growth pools while the nplant decrement
   lands on post-growth pools; the error is `died_nplant × npp`, one-signed.
7. **External seed rain** (`seed_rain_recruits`, fully endowed recruits) is a legitimate boundary
   **import**. Declared, not eliminated.
8. **Reproduction is aliased, not leaked** — recorded here so it is not mistaken for a leak later.
   `npp_repro` is debited every slow step; `apply_recruitment` runs only under `is_new_month` and treats
   that one day's rate as the whole month. Unbiased for a steady rate, but reproduction tracks NPP
   seasonally, so it is a 12-sample estimator of a daily flux. A correctness question for the demography
   cadence, not a ledger term.

#### 10.2.3 Water terms

- `shed_turnover_water` → `patch%shed_water_rate` → the fast loop's `precip_ground` is **the model of
  what every other seam here should look like**: the mass leaves one store, a named variable carries it,
  and another store receives it. Already closed; the ledger just declares it.
- `reconcile_tissue_water_capacity` **already returns `seeded` and `discarded`**, and its header already
  says they are booked nowhere and become ledger terms when this lands. Free.
- **Cull and disturbance-kill water is discarded** (1B #7). `terminate_cohorts` routes carbon to litter
  and drops `leaf_water_mass`, `wood_water_mass` and both films; `apply_patch_disturbance` does the same
  for the killed canopy. A pure sink.

#### 10.2.4 Energy terms — and two that are bugs, not bookkeeping

- **`blend_cas` is wrong for patches of different height** (1B #9, and it needs fixing whether or not the
  ledger lands). `can_enthalpy`, `can_shv` and `can_co2` are *specific*; the extensive content is
  `area × depth × ρ × value`. The function area-weights the intensive values while *separately*
  area-weighting `can_depth`, so fusion drops the covariance term `w₁w₂(d₁−d₂)(h₁−h₂)`. Fusing a tall
  patch with a gap corrupts canopy-air energy, humidity **and** CO2 at once.
- **`cas_set_depth`'s open-volume term is computed nowhere.** The routine already carries the `de_open`
  hook and a header explaining that the jumps that matter are disturbance and fusion — a 20 m canopy
  becoming a 1 m gap in one step. `refresh_canopy_depth` calls it without `rho_air`/`de_open`, so the
  term is never formed.
- **Fusion weights both tissue temperatures by leaf area.** `fuse_cohort_fast_state` blends `leaf_temp`
  *and* `wood_temp` on leaf area. For leaves that roughly tracks heat capacity; for **wood** it does not
  — wood heat capacity follows sapwood+fine-root carbon and the water in it, which does not scale with
  leaf area. And `leaf_water_mass` is nplant-weighted while its temperature is leaf-area-weighted, so
  their product — the actual energy — is conserved by neither weighting.
- Shed water leaves tissue at tissue temperature and re-enters the ground at `t_film_valuation`
  (deliberate, and documented as such — still an energy term).
- Cull and disturbance-kill tissue heat, alongside the water above.

#### 10.2.5 Tolerances

`fuse_2_cohorts` asserts its carbon pools at `conservation_tol` (`1e-3` in the shipped configs). That
operation is exact algebra and should hold at round-off — `split_cohorts` already asserts its film water
at `1e-12`. A 0.1 % window on every fusion is eleven orders too loose. Per 1B #10: fuse/split assert at
~1e-12 relative, `conservation_tol` is **reserved for the site ledger**, where the declared terms
genuinely carry approximation, and `size_tol` is deleted.

#### 10.2.6 Implementation order

**The skeleton lands before any fix.** Snapshot the stores, declare the terms that are already honest
(shed water, reconcile seed/discard, seed rain, litter export), and route everything else into an
explicit per-currency `unattributed` bucket. That makes each gap *measurable* before deciding what it
deserves — the lesson from #128 and #129, where measurement inverted the expected ranking twice. The
prediction on record is that growth respiration dominates the carbon bucket by an order of magnitude and
that `blend_cas` dominates energy only in runs with active patch fusion; the numbers, not this
paragraph, decide the order of the fixes that follow.

Then, in their final homes:

- `slow_dynamics/driver/meds_slow_ledger` — `site_ledger_t`, the snapshot, and the daily/annual assert.
- `slow_dynamics/plant/`: recruitment carbon debit vs `init_cohort` endowment (1B #4); `growth_resp` and
  `deficit` routed into the CAS carbon balance (1B #5).
- `slow_dynamics/demography/`: mortality, cull and disturbance-kill hand tissue water, film water and
  tissue heat to a declared sink (1B #7); `blend_cas` depth blending and the `terminate_patches`
  survivor rescale fixed (1B #9); the fuse/split tolerances of 10.2.5 (1B #10).
- `slow_dynamics/soil/`: CENTURY transfer-matrix column conservation asserted; `rh_seam_gap` asserted in
  production; `audit%litter_in` includes cull/disturbance litter (1B #8, #6).
- Not structure at all, listed so it is not lost: `docs/ed2_comparison.md` soil-temperature statements
  predate the PR #119 root-heat-sink fix and must be re-run.

#### 10.2.7 What the skeleton measured (2026-09-09)

The skeleton shipped and ran. Numbers are **site totals per m², cumulative over a 1096-day (3-year)
Ithaca run** from a small establishing stand, ifx, byte-identical to `main` on all 75 outputs. Read
`residual` against `declared`, which is the gross boundary flux the phase did account for.

| phase | carbon [kgC] | declared | water [kg] | energy [J] | marks |
|---|---|---|---|---|---|
| allocate         | –          | 9.1e-4 | **+1.3e-12** | **+1.3e-6** | 1096 |
| grow + mortality | **−3.0e-3** (−3.4e-4 with soil C on) | 2.5e-3 (5.2e-3) | **−4.0e-4** | **−272** | 1096 |
| recruit          | **+7.7e-3** | 0 | 1e-13 | **+1.08e5** | 36 |
| cohort fuse/fiss | 3e-18 | 0 | 1e-13 | **−6.8e4** | 36 |
| disturbance      | 2e-18 | 0 | 5e-13 | **+3345** | 3 |
| patch fuse/term  | 2e-17 | 0 | 1e-13 | **−3345** | 3 |
| canopy depth     | **−9e-18** | 3.6e-3 | **+7e-13** | **+3.6e-7** | 1096 |
| soil carbon      | **+9.0e-19** | 3.5e-3 | – | – | 1096 |

**Three phases are verified closed**, which is what makes the rest trustworthy: the canopy-depth
open-volume declaration closes to 1e-6 J against **5.77e6 J** declared — eleven orders — so the
entrainment term is both large and now fully accounted; the CENTURY step closes to 9e-19 against
3.5e-3 declared once `rh_today` is declared; and the turnover-shed seam closes to 1e-12.

**The prediction in 10.2.6 was wrong, and this is the record of it.** It said growth respiration
would dominate the carbon bucket by an order of magnitude. It does not. **Recruitment creates
+7.7e-3 kgC/m², twenty-three times the growth phase's −3.4e-4**, and it does so in 36 events rather
than 1096. Two honest caveats before that is read as settled: this is a 3-year *establishing* stand,
where recruitment is proportionally far larger than in a mature forest, and the growth-phase figure
is a *net* of terms with opposite signs (growth respiration destroys, the starvation deficit and the
floors create). Both need a mature run and a decomposition before the fix order is set. What is not
in doubt is that the ranking was not the one predicted.

**With soil carbon off — the default — the necromass export is 2.7e-3 kgC/m²**, the difference
between the two grow-phase figures, and the largest single carbon term in a stock configuration.
That is item 5 of 10.2.2, now measured. With soil carbon on, the growth phase loses **6.4 % of the
carbon handed to it by the fast tier** (−3.4e-4 against 5.2e-3 declared).

**Two findings the survey did not predict:**

- **Recruits are born with tissue heat from nowhere: +1.08e5 J over 36 events.** `init_cohort` sets
  `leaf_temp`/`wood_temp` to `LEAF_TEMP_INIT` and `set_cohort_size` immediately gives the cohort a
  heat capacity, so a recruitment event creates sensible heat in proportion to the biomass it also
  creates. This is the energy twin of 10.2.2 item 3 and was not on the list.
- **Disturbance (+3345 J) and patch fusion (−3345 J) are exactly equal and opposite**, to every
  digit printed. The gap patch's canopy air is created by `blend_cas` and reabsorbed by it; the sign
  symmetry says the two are the same arithmetic run forwards and backwards, which is a useful
  constraint on any fix to `blend_cas`.

Cohort fusion's **−6.8e4 J** is 10.2.4's leaf-area-weighted temperature blend, measured; mortality's
**−4.0e-4 kg** is 1B #7's discarded tissue water, measured.

**One correction to this section's own design.** `shed_water_rate` is *not* a store, though 10.2.3
implies it is. It is a handoff whose lifetime is the FAST window, not the slow step: by the time the
next slow step opens, that water is already in the soil and the rate variable still holds it, so
carrying it as a store double-counts at every open — a steady one-signed ~8e-7 kg/step phantom leak
in the allocate phase, which is exactly how it was found. It is declared instead, the mirror of the
fast→slow carbon handover.

**Deferred from the skeleton, on purpose.** `slow_site_store` duplicates four store loops from
`meds_column_state_ops` (soil water, soil energy, plant water, canopy film) and the tissue
heat-capacity construction from `meds_fast_frozen`, because those live in `meds_fast` and `meds_slow`
must not link it. Moving them down to `state/column`, so both tiers value the stores with one piece
of code rather than two that agree today, is the natural next commit and belongs with the energy fix.

#### 10.2.8 The mature stand overturns 10.2.7's ranking (2026-09-09)

10.2.7's headline — that recruitment dominates the carbon bucket — was measured on a 3-year
*establishing* stand and 10.2.7 said in as many words that the stand flattered recruitment. It did.
Re-run from the `runs/ithaca_ark30` 2074 spin-up restart (114 cohorts, LAI 5.6, AGB 16.7 kgC/m²,
mean dbh 37 cm, 1 PFT, soil carbon on), same 3-year span, same ledger:

| | establishing | **mature** |
|---|---|---|
| grow+mortality carbon | −3.02e-3 kgC/m² | **−2.5102 kgC/m²** |
| …as % of declared | 119 % | **25.3 %** |
| recruit carbon | +7.67e-3 | +7.61e-3 |
| **grow : recruit** | **0.39 : 1** | **330 : 1** |

The establishing stand nets **0.85 gC/m²/yr**; the mature stand **837 gC/m²/yr lost in the growth
phase against ~2925 gC/m²/yr GPP**. Seed rain is a fixed 0.01 plant/m²/yr, so on an unproductive
stand it swamps everything and on a real one it is noise. **10.2.6's original prediction was right:
growth respiration leads, and the measured 25.3 % sits right beside the `g/(1+g)` = 23.1 % the
construction cost implies.** The lesson is not about growth respiration; it is that a conservation
ranking measured on a stand that is not growing measures the stand, not the model.

The soil-carbon phase closes to **−1.1e-13 against 3.41 kgC/m² declared** on a productive forest,
and the canopy-depth phase to 1e-6 J against 4.42e6 J. Both hold at the mature scale.

#### 10.2.9 Birth and death become paired transfers (2026-09-09)

Implementing the author's framing: a recruit **draws** what it arrives with from outside the model,
and a death **hands on** what it carried. Measured on the mature stand, before → after:

| term | before | after | |
|---|---|---|---|
| grow+mortality water | −1.2325 kg | **+6.0e-12** | closed |
| disturbance water | −0.8852 kg | **−3.8e-6** | closed |
| recruit energy | +9.54e4 J | **−7.7e-7** | closed |
| disturbance energy | −3.498e6 J | **+9.9e3** | 99.7 % |
| recruit carbon | +7.61e-3 | **+9.60e-4** | 87 % |
| cohort fuse energy | −8.38e4 J | −8.30e4 J | unchanged (weighting, not a transfer) |
| grow+mortality energy | +5.18e6 J | **+1.00e7 J** | *exposed*, see below |

**Why birth draws externally.** Recruitment stands in for everything between a seed and a 2 m
sapling — germination and the seedling's own photosynthesis, transpiration and energy balance —
and the model tracks no cohort below `min_cohort_height`. What a recruit arrives with was fixed and
absorbed by a size class that is not represented, so it is genuinely external. Drawing it from the
free atmosphere rather than the patch's canopy air is deliberate: crediting the CAS with seedling
uptake while representing none of the seedling's respiration, transpiration or shading would add one
term of a missing process and call it an improvement.

Only the part the model did **not** already pay for is declared. `recruit_pool` carries reproduction
carbon at `carbon_min` per plant and that much *is* debited from the parents, so the draw is
`endowment − carbon_min`, plus the baseline seed rain at its true entry point (the monthly pool
credit — which arrives whether or not anything is born that month, a distinction worth 2/3 of the
term). The residual 9.6e-4 that remains is the productivity-driven reproduction carbon: debited from
parents in the growth phase, re-created here at a quantity the `repro_carbon_efficiency / carbon_min`
conversion does not preserve. That is 10.2.2 item 3, deliberately left visible rather than declared
away.

**Birth temperature was the real energy bug.** `init_cohort` stamped `LEAF_TEMP_INIT = 288.15 K`, one
global constant, so a sapling appearing in an Ithaca January was born ~20 K warmer than the air it
stood in. Recruits now start at their patch's canopy-air temperature (guarded: an unstepped CAS falls
back to the constant). Declaring the old term would have been declaring an artifact.

**Death hands its water to the ground down the channel turnover shedding already uses** — one
verified path rather than a second mechanism — for all three death paths (background mortality, the
cull, the disturbance kill). The tissue **heat** leaves the thermal system with the necromass and is
reported, because the CENTURY pools it becomes carry no temperature; a litter thermal store is where
it would belong if one existed.

**What this exposed.** The growth phase's energy residual *rose*, from +5.18e6 to +1.00e7 J, because
mortality's heat was partly cancelling it. Growing biomass raises the tissue heat capacity at
constant temperature, so `cap × T` rises with no flux: **the model creates thermal mass out of
carbon**. That is the birth-side twin of the death-side term just fixed, it is now the largest energy
item in the ledger, and it was invisible until the offsetting term was removed.

**Still open after this**, in measured order: growth respiration and the growth phase's carbon
(−2.51 kgC/m², 25 % of the handover); the growth-side thermal mass (+1.00e7 J); cohort fusion's
leaf-area-weighted `wood_temp` blend (−8.3e4 J) and `blend_cas`; the item-3 reproduction conversion
(+9.6e-4).

#### 10.2.10 The allocator's outputs all get a destination (2026-09-09)

`plant_carbon_allocation` produces seven outputs. Four are pools the driver commits. **Three left the
plant and went nowhere**, and they are one sentence, not three problems: the allocator's outputs had
no destinations. All three are fixed together.

- **Growth respiration → the canopy air.** Charged against the plant, reaching only the
  `CS_GROWTH_RESP` diagnostic, while `nee_biotic` carried maintenance respiration alone. It is now
  handed down as a frozen daily rate on a new per-patch `slow_co2_rate`, added to `nee_biotic` by the
  prepass — the carbon twin of `shed_water_rate`, built to the same read-only, seeded-once discipline.
- **The starvation `deficit` → the same channel, opposite sign.** The fast loop has *already* exhaled
  the full maintenance respiration; when storage could not fund it, `deficit` is the part no pool paid
  for. Real maintenance respiration is substrate-limited, so the honest reading is that the fast loop
  **over-reported**, and an operator-split model corrects an over-report on the next step rather than
  rewriting the last one. Netting it into one signed channel keeps one sign convention.
- **The unestablished seed fraction → litter.** Reproduction carbon is debited in full from the parent
  and only `repro_carbon_efficiency` of it establishes. The remaining ~99.9 % is dead seed and dead
  seedling — **necromass, not nothing** (author's decision). It enters `necromass_to_litter` as
  `storage_c`, which pools with the canopy and splits on `f_labile_leaf`: seed tissue is labile and
  canopy-derived, which is what that argument means.

**Result on the mature stand.** The growth phase's carbon residual falls from **−2.5102 to −7.337e-4
kgC/m²**, a factor of **3420**, and the worst single step from 5.06e-3 to 1.81e-6. What remains is
0.006 % of the declared flux: the pool and `nplant` floors (item 4) and the pre/post-growth offset in
the mortality valuation (item 6). The soil-carbon phase still closes at 7e-15, now against 4.50
declared rather than 3.41 — the seed litter is real carbon arriving in real pools.

**The NEE shift is the missing respiration, to within 0.4 %.** Site NEE moves from −5.988 to −4.616
µmol/m²/s, **+1.372**. The growth respiration the ledger said was absent is ~1.55 kgC/m² over 3 years
= 0.517 kgC/m²/yr = **1.366 µmol/m²/s**. Two independent routes to the same number: the ledger's
residual before the fix, and the flux difference after it.

| | main | branch | |
|---|---|---|---|
| `nee_site` | −5.988 | **−4.616** | +22.9 % — the site is a much smaller sink |
| `soilc_total_site` | 1.775 | 2.181 | +22.8 % — seed necromass accumulates |
| `rh_site` | 0.01091 | 0.01555 | +42.6 % — and decomposes |
| `gpp_site` | 0.243774 | 0.243946 | **+0.070 % — the CO2 fertilization feedback** |
| `agb_site` / `lai_site` / `nplant_site` | — | — | +0.005…0.007 % over 3 yr |

The GPP rise is small but it is the feedback predicted when the soil-carbon fallback was deleted in
#128, running the other way: putting CO2 back into the canopy air raises photosynthesis. Structure
barely moves, which is right — the carbon *through* the plant is unchanged; only its fate after
leaving is.

**A test the ledgers cannot replace.** Neither ledger can catch a break in this channel: the slow
ledger declares the **handoff**, so it closes whether or not the fast loop ever picks the rate up, and
the fast CAS ledger closes around whatever `nee_biotic` it is handed. Only a test binds the two ends.
`test_slow_ledger` gains three differential assertions — growth respiration reaches the channel and a
zero construction cost empties it, the unestablished seed fraction reaches litter, and a **starving**
stand owes a **negative** flux. All three mutation-tested.

Two flaws in the first version of that test, worth recording because both made it pass while asserting
nothing: it called `vegetation_dynamics` twice on one site, and since that call *commits* growth the
two variants compared different forests; and its fixture cohort was below `min_reproduction_height`,
so no seed carbon existed to lose. Each variant now runs on a fresh stand, sized above the threshold
and given a leaf lifespan long enough that turnover does not consume the supply before reproduction is
reached.

**Still open**, in measured order: the growth-side thermal mass (+1.00e7 J — biomass growth raises the
tissue heat capacity at constant temperature, the birth-side twin of the death-side term 10.2.9
fixed); cohort fusion's leaf-area-weighted `wood_temp` blend (−8.3e4 J) and `blend_cas`; item 3's
remaining half, the recruit pool as a carbon quantity, together with the monthly-sampling aliasing
(+9.6e-4); the pool and `nplant` floors (−7.3e-4).

#### 10.2.11 Tissue thermal mass (2026-09-09)

A cohort's heat capacity is a function of its biomass and its density, so growing or dying changes
`cap × T` **with no flux at all**: the model makes thermal mass out of carbon, and unmakes it. Both
directions were undeclared and they partly cancelled, which is why the growth side stayed hidden
until 10.2.9 declared the death side and the residual *rose* from +5.18e6 to +1.00e7 J.

**Declared as one exchange, not two.** The whole change across the growth commit is one mechanism,
and splitting it into a growth part and a mortality part needs a cross-term convention the physics
does not supply (the `hcap_min` floor is not linear in density). Both `shed_mortality_water` and
`update_cohort_states` leave the temperatures alone, so the change across them is *purely* thermal
mass. 10.2.9's separate mortality-heat term is therefore withdrawn — subsumed, not repealed — and its
water term stays, because water going to the ground is a transfer and this is not.

**Why it is an exchange and not a leak.** New tissue is assembled from CO2 and water at the plant's
own temperature and arrives carrying the sensible heat of that mass. The model tracks no thermal
content for CO2 — the canopy air's capacity is dry air alone — so that heat genuinely crosses the
boundary of the modelled thermal system. A fuller treatment would give CO2 a heat capacity in the CAS
and carry enthalpy on root water uptake; both are far larger than this, and naming the approximation
here is better than burying it.

**Result: +1.00e7 → +2.06e-6 J against 8.00e6 declared.** Twelve orders. The growth phase now closes
on all three currencies to round-off except carbon's −7.3e-4 (the pool and `nplant` floors).

**Byte-identical**, because nothing here changes the model — `slow_tissue_heat` is a pure read, and
the only other change is deleting an output that is now computed elsewhere. Mutation-tested: with the
declaration removed the phase carries +5.20e6 J, which is the *net* of the two directions and exactly
the figure 10.2.7 reported before the death side was declared.

**Also folded in:** `slow_site_store` carried its own copy of the tissue heat-capacity formula from
before the shared one existed. It now calls `cohort_tissue_heat_capacity`, so the ledger, the
demography operators and the slow driver cannot drift on what a cohort's thermal mass is. Two of the
duplications 10.2.7 flagged are gone; the four store loops in `meds_column_state_ops` remain.

**Still open**, in measured order — and the character of the list has changed. What is left is no
longer *missing transfers* but **wrong averages**: cohort fusion's leaf-area-weighted `wood_temp`
blend (−8.3e4 J), `blend_cas`'s depth-blended intensive quantities (disturbance +9.9e3 J and patch
fusion −1.0e4 J, still equal and opposite), item 3's remaining half (+9.6e-4), and the pool and
`nplant` floors (−7.3e-4). Those need a fix to the weighting, not a new declaration.

#### 10.2.12 The averages (2026-09-09)

Three wrong averages, and the thermal-mass declaration that 10.2.11 established extended to the
structural phases. **Energy and water now close everywhere in the ledger; only two carbon terms
remain.**

**`blend_cas` weights on AIR MASS, not ground area.** Enthalpy, specific humidity and CO2 mixing
ratio are per kg of air, and the air mass is `area × depth`, so the conserving weight is `area ×
depth`. Weighting on area alone dropped the covariance `w₁w₂(d₁−d₂)(v₁−v₂)`, corrupting canopy-air
energy, humidity **and** CO2 together by the same relative error. It was exact only when the two
depths matched — precisely the case that needs no blend. `rho` cancels (site-uniform), and
`can_depth` stays area-weighted because that is what a depth is.

**Cohort fusion weights the tissue temperatures on HEAT CAPACITY.** They were leaf-area weighted,
which is roughly right for leaves — leaf capacity tracks leaf carbon, and leaf area tracks that
through `sla` — and not even approximately right for **wood**, whose capacity follows wood carbon
and the sapwood ring. Capacity weighting conserves `cap × T` exactly for the additive part, since
every carbon pool and both tissue waters are nplant-weighted and so add across the merge.

The diagnostic twins keep their leaf-area weights, and 10.2.4's claim that the two tables share a
weight is withdrawn: a per-leaf-area *diagnostic* really is leaf-area weighted; a prognostic
*temperature* is per unit heat capacity. Treating them as one kind is what made the wood wrong.

**`terminate_patches` MERGES the sliver instead of deleting it.** It used to drop any patch under
`min_patch_area` and renormalize the survivors' areas back to 1, which silently redistributed the
whole site: every conserved quantity changed by `(a/(1−a))·Σ_kept a_i X_i − a·X_dropped`, zero only
if the doomed patch held the survivors' mean — and a doomed patch is atypical by construction,
usually a fresh gap. It is now fused into the largest survivor with `fuse_2_patches`, which
conserves exactly and reuses the operator patch fusion already depends on. The renormalisation
becomes a round-off correction rather than a redistribution.

**And the structural phases now bracket their thermal mass**, as the growth commit does. Fusion and
fission re-derive the sapwood ring from the merged or perturbed diameter, and `sapwood_fraction` is
**nonlinear in dbh**, so the merged capacity is not the sum of the two even when every carbon pool
is. That is a thermal-mass change of exactly the kind growth makes. It is declared *only after* the
weighting was fixed — otherwise the declaration would have been hiding the wrong average rather than
accounting for what the right one leaves behind. 10.2.9's separate cull and disturbance-kill heat
reports are withdrawn into it, as 10.2.11's mortality report already was; their **water** reports
stay, because water reaching the ground is a transfer.

| phase, energy [J] | before | after |
|---|---|---|
| cohort fuse/fiss | −8.30e4 | **−1.30e-6** |
| disturbance | +1.01e4 | **+1.49e-7** |
| patch fuse/term | −1.00e4 | **+3.43e-7** |

Water at those phases closed too: patch fusion −2.2e-5 → 3.4e-13, disturbance −3.8e-6 → 1.1e-13,
and `patch fuse/term` carbon 8.3e-8 → −1.5e-11.

**What is left in the entire ledger is two carbon terms**, both known and both named here already:
the recruit pool's productivity-driven credit (+9.61e-4, item 3's remaining half, needing the pool
to become a carbon quantity and the monthly sampling to stop aliasing) and the growth phase's pool
and `nplant` floors (−7.34e-4, items 4 and 6). Energy and water close on every phase.

**A bug this nearly shipped with.** The first version of the sliver merge rebuilt the CSR map once
after the loop. `fuse_2_patches` reads the *receptor's* CSR slice to rescale its cohorts, and
`patch_fuse_pass` — the existing caller — rebuilds after **every** fusion for exactly that reason.
With two slivers the second merge would have read a stale slice. The test now uses two.

**And a blind fixture.** `test_patch`'s CAS-fusion assertion had both patches at the default 20 m
depth, where the area-weighted and mass-weighted answers agree exactly. It passed against the wrong
code and would have passed against the right one. It now uses 30 m and 10 m, checks the
mass-weighted value, and separately asserts the extensive content survives — the property, not the
formula. `test_fusion_cohort` was rewritten the same way: it asserted the old leaf-area formula, and
now asserts that leaf tissue energy is conserved **exactly** and wood to within the sapwood
re-derivation.

#### 10.2.13 Mortality is valued on what the applier removed (2026-09-09)

The live carbon store is `nplant × pool`, and its change over the growth commit decomposes exactly:

    n₁p₁ − n₀p₀  =  n₀(p₁ − p₀)  +  (n₁ − n₀)p₁

— **growth at the old density, mortality at the new pools.** Every declaration was on the other
side of both terms: `accumulate_mortality_litter` valued the litter on the pre-growth pools `p₀`,
and the fast→slow handover was taken after the commit, at `n₁`. Both now match the decomposition:
the mortality routines run **after** `update_cohort_states` on the committed pools, and the handover
is captured **before** it.

The density drop is taken as `n₀ − n₁` rather than recomputed from the Camac hazard, which fixes
item 4's other half for free: when the `negligible_nplant` floor stops a cohort dying, `n₀ − n₁` is
smaller than the hazard implies, so litter is no longer credited for individuals still standing.

**What this leaves is one term, not several.** The growth phase's carbon residual is now
**−9.170e-4** and the recruit phase's **+9.607e-4** — the same reproduction carbon, debited from
parents in one phase and re-created in the other, with nothing connecting them. They sum to
**+4.4e-5**, which is the seed rain, the pool floors and the monthly-sampling asymmetry. Energy and
water close on every phase; the *entire* remaining ledger is item 3.

**Two mistakes in getting here, both recorded because both were invisible to everything else.**

The first estimate had item 6's sign backwards. Re-basing the litter on `p₁` makes the declared
export *larger* (the pools grew), which pushes the residual positive; re-basing the handover on `n₀`
makes the declared import larger, which pushes it negative — and the handover term dominates. The
measured residual moved from −7.34e-4 to −9.17e-4, in the direction the arithmetic says once both
changes are counted rather than one.

The second was a real bug that survived a green suite: moving `shed_mortality_water` after the
commit left it calling `cohort_tissue_water`, which carries the cohort's *current* nplant — so the
loss was scaled by `n₁/n₀` instead of being `(n₀−n₁) × per-plant`. The growth phase's water went
from 1e-12 to **1.22e-4** and nothing but the ledger noticed. `test_slow_ledger` now asserts the
identity directly — every kg that leaves tissue arrives in the patch shed channel — and both that
bug and the no-routing case fail it.

#### 10.2.14 Reproduction carbon becomes a flow, and §10.2 closes (2026-09-09)

The last item, and it needed **no new state**. `recruit_pool` was credited inside
`apply_recruitment`, monthly, from whatever recruitment rate the driver had computed on that one
day, scaled up to stand for the month. Crediting it **every step instead** fixes both halves at
once:

- the **cadence**: a 12-point sample of a quantity the model computes 365 times a year, whose value
  depended on which days happened to be month boundaries, becomes the exact integral;
- the **carbon link**: `recruitment × dt_yr` is `n·npp_repro·efficiency / carbon_min`, so the pool
  valued at `carbon_min` is exactly the establishing share of the reproduction carbon the parents
  were debited for *in that same step*. Debit and credit land in the same phase, and the ledger sees
  a transfer rather than carbon vanishing in one phase and appearing in another.

The seed-rain declaration moves with it, from monthly to daily.

**And the gap inherits the seed bank.** `apply_patch_disturbance` zeroed `recruit_pool` on the new
gap while inheriting `soil_carbon` and `xi_accum` from the same donors. A treefall gap does not
sterilise the ground it opens: the carry-forward pool sits in the soil with the litter and the
CENTURY carbon. Zeroing it destroyed the pool's carbon on the disturbed fraction — the last
non-round-off term in the ledger, **−2.96e-6 kgC/m²** over three events.

| carbon [kgC/m²] | before 10.2.13 | after 10.2.13 | **now** |
|---|---|---|---|
| grow + mortality | −7.34e-4 | −9.170e-4 | **−1.98e-13** |
| recruit | +9.61e-4 | +9.607e-4 | **+1.30e-14** |
| disturbance | −4.02e-6 | −2.96e-6 | **−1.13e-11** |

**§10.2 is closed.** Every phase closes on every currency to round-off; the largest residual
anywhere in the ledger is 1.5e-11 against declared fluxes of order 1 to 12. The report now ends in
a **verdict**: each phase and currency is judged against the flux it declared plus a per-mark
absolute floor, and the run says so in one line. Reported, not fatal — the same choice the fast
loop's whole-column ledgers make, because a run that breaches this is telling you something and
stopping it half-way tells you less than finishing it.

**One correction to what 10.2.13 said was coming.** It called the monthly sampling a case where
"353 of 365 days' reproduction carbon is never sampled". That overstates it: the estimator is a
12-node rectangle rule targeting the annual total, not a wholesale loss, so the error was quadrature
rather than a leak — and the conservation break was the separate fact that the debit and the credit
were computed from different quantities. Both are fixed by the same one-line move, but they were
two faults, not one.

**The verdict's own tolerance had to be corrected on first use**, which is worth recording because
it is the same class of mistake the ledger keeps finding. The first version used a per-mark
absolute floor plus a relative test on the declared flux — and immediately flagged two phases.
It was right to: `disturbance` and `patch fuse/term` declare *no* carbon, so their whole allowance
was `3 × 1e-12`, while the structural operators permute, merge and renormalise a ~25 kgC/m² store
and cost ~1.3e-11 in arithmetic doing it. Round-off scales with the **store**, not with the number
of checks. The tolerance now carries a store-proportional round-off allowance of 1e-11 per mark.

That is *not* the store-relative tolerance `budget_check`'s header warns against: 1e-11 is a
round-off bound five orders tighter than the 1e-6 that let a sustained 1 W/m² leak hide, and the
smallest **real** term this ledger ever caught — the disturbance seed bank, 2.96e-6 — is still four
orders above it.

**Verification.** 42/42 on ifx and nvfortran. Not byte-identical, and cannot be — mortality water now
reaches the soil and recruits are born at a different temperature. Site integrals move by ≤4e-8
relative (`veg_carbon_site`, `gpp_site`, `nee_site`, `agb_site`, `nplant_site` all unchanged to six
digits); soil carbon −0.004 %, Rh +0.04 %. The per-cohort `dmax_psi_leaf` distribution moves by a
median of 6e-7 with a handful of cohorts — the recruits whose birth temperature changed — moving up
to 0.068 across a range of 0.35, which is the intended change and confined to them.

### 10.3 Fast-integrator leftovers → **migration step 7**

Step 7 is the one step that rewrites the march signatures on both schemes and the oracle, so
everything that needs those signatures rides with it:

- Pass `column_params_t` (the six parameter records) via `column_config_t` instead of copying it
  into `column_frozen_t` — the review deferred this precisely because of the signature churn.
- The **silent-omission matrix**: every `column_state_t` field is enumerated in ~12 places with
  no compile-time failure on omission (`zero_like` never allocates the films; `state_err_diff`
  leaves the pond default-initialised while `state_sub` subtracts it). Give the state vector one
  field-visitor (`state_visit`) or one contiguous packed layout, and make the combinators, the
  process mask, pack/unpack and the error norm all go through it. A regression test that adds a
  dummy field is the acceptance check.
- Item 1A (vi)/(vii): per-layer `budget_imbalance` over faces on the committed path, per-cohort
  tissue water/energy residuals, and RK45 ledgers asserted **after** the rail decision, not on a
  step the dispatcher may roll back.
- Item 5 helpers still duplicated, moved verbatim: `relieve_theta_bounds` (three implementations
  → `fast_dynamics/soil/soil_water`), `soil_layer_temp` (six inline copies → `base/therm_lib`),
  `seed_soil_column` and `seed_plant_water` (driver + 4 and + 3 tests → `init/` and the 10.1
  reconcile respectively).

### 10.4 Core facade (review item 4 #13) → **migration step 9, dissolved by decision #5**

The review found 12 src modules and 10 tests reaching past `meds_core_interface`. Listing what
they reach for settles the question: `site_t`, `site_alloc`/`site_free`, `cohort_ensure_capacity`,
`rebuild_csr`, `set_cohort_size`, `gather_pft_params`, `DMAX_PSI_LEAF_*`, `GROWTH_AVG_UNSET`, and
every `meds_core_diag_types` block and constant — all **state**, none of it a verb. A verb facade
was the wrong shape for that half of core, which is why it was bypassed. After decision #5 the
state half is the `state/site` layer and is imported directly by design (a state layer has no
facade; rule 6 applies to libraries of operations). `meds_core_interface` then either becomes
`meds_demography_interface` over the `slow_dynamics/demography` operators (`update_*`, `*_fuse_*`, `terminate_*`, `split_cohorts`,
`apply_recruitment`, `apply_patch_disturbance`, `sort_*`) or is deleted in favour of importing
those from `slow_dynamics/demography` directly. Either is consistent; deleting is fewer legal spellings.

### 10.5 Remaining renames (review step 6) → **migration step 10, merged into decision #14**

All byte-identical, `sed -I -w`, one commit per group, after every move so the paths in the
commit are final:

- Fields: `wcap/ccap` → `cas_mass_capacity/cas_molar_capacity`; `gah/gaw/gac` →
  `g_cas_atm_heat/vapour/co2`; `hydro` (soil options) vs `hydro_p/hydro_o` (plant) →
  `soil_water_opts` / `hydraulics_params/opts`; `h_coeff_f/g_tr_f/g_film_f` → `_leaf`;
  `a_leaf/a_wood/a_store` → `*_hcap_per_dt`; `enth_atm` → `enthalpy_atm`; `snowf/tair/precip` →
  `snowfall/air_temp/rainfall`.
- Routines whose names describe a mechanism they no longer have: `column_hydrology_flux`
  (commits state) → `advance_soil_water_column`; `uext_to_temp`/`temp_to_uext` (ED2 token) →
  `internal_energy_to_temp`/`temp_to_internal_energy`. `column_prepass` stays: since PR #123 it
  is an orchestrator over six named routines and "prepass" describes when it runs.
- The `t_precip`/`rain_temp`/`film_u_ref` question (item 4 #9) is a *naming* question, not a
  merge — `rain_temp` is the film valuation temperature pinned to `tsupercool_liq` under a pack,
  `t_precip` is the pond inflow temperature. Name them for what they are here.
- Keep the review's exclusions: netCDF registry strings and TOML keys are not renamed without a
  compatibility note; `agb`, `lai`, `wai`, `swe`, `ustar`, `ggnet`, `can_*`, `gbh/gbw/gsw`,
  `dbh`, `nplant`, `pft` stay.

### 10.6 What this changes in the migration table

Only the ends. Step 0 gains a content edit (10.1) and is no longer "wait"; step 7 gains the
signature-dependent leftovers (10.3); step 9 becomes concrete (10.4); step 10 is new (10.5); the
slow-loop physics (10.2) is scheduled after step 6 but is not a step of this plan. Steps 1–6
are unchanged and remain pure moves.


---

## 11. What the implementation changed about this plan (2026-09-09)

Steps 1-6 are done. Five things the plan got wrong or left implicit, recorded here so the remaining
steps are planned against the tree that exists.

### 11.1 `meds_fast_prepass` belongs in `numerics/`, not `driver/`

§4's tree filed it under `fast/driver/`. It cannot go there: `meds_fast_ark`, `meds_fast_rk45` and
the RK4 oracle all call `column_prepass` mid-stage, so `numerics -> driver` and `driver -> numerics`
both exist and the folder split makes that a hard CMake cycle. The plan placed it on the strength of
its name. `driver/` holds `meds_fast_dynamics` alone: the loop that walks a slow step in `dt_fast`
sub-steps over the patch axis. Everything the marches call is `numerics/`.

### 11.2 The plant-types split was load-bearing, not cosmetic

§5 noted that `meds_plant_types` "splits just as cleanly" as the kernels. It is stronger than that:
without the split the two halves of `src/plant/` import each other (`meds_phenology` needs the pheno
types; the facade needs phenology), which is the same cycle as 11.1. The pheno half is now
`meds_pheno_types` in `slow_dynamics/plant/`.

### 11.3 `meds_plant_interface` forces a fast -> slow kernel-library edge

The facade re-exports phenology and carbon allocation, so `meds_fast_kernels` links
`meds_slow_kernels`. Acyclic (nothing in `slow_dynamics/` imports a fast kernel) and harmless, but it
is an artefact of the facade, not of the physics. Step 9 removes it. The other fast -> slow edge,
`meds_fast -> meds_slow_kernels`, is the real one: the documented S1 Rh seam.

### 11.4 `necromass_to_litter` moved at step 2, and `meds_core` grew one link

§8 scheduled the move for step 2 but §5's S3 justified it "under the new DAG", which did not exist
until step 3. It moved at step 2 anyway, as its own module `meds_litter_partition` in the
biogeochemistry folder, with `meds_core` gaining a link to it. That edge is exactly what the plan's
DAG carries -- the demographic operators sit ABOVE the stateless kernels -- so paying for it early
cost nothing and kept step 6 a pure move.

### 11.5 `libmeds_plant_c` must survive step 6

Moving `meds_plant_capi.f90` into `src/capi/` silently dropped the `meds_plant_c` target, because the
old one globbed `src/plant/*_capi.f90`. `python/meds/plant/` dlopens `libmeds_plant_c.so`, so the
package breaks with no build error. The target is restored, globbing the new location. Merging the
two `.so` files is decision #1 and belongs to step 8 with the packaging work, not to a file move.

### 11.6 Step 0 is NOT a move, and is re-scoped out of the first PR

The author's expectation was that steps 0-6 are "mostly renaming and deleting unused code". That is
true of 1-6 and false of 0, which is the largest content change in the whole plan:

- `column_cohort_t` is not dead code. It is threaded through `column_fast_step`, `column_prepass`,
  both marches, the oracle and three kernels' signatures (6 src + 3 test files), and it carries 18
  allocatable array components.
- §10.1's replacement mechanism -- "the driver reads the patch's contiguous CSR cohort section
  directly" -- is under-specified for those 18 components. The kernels must stay `site_t`-free, so
  the caller must slice and the callee must take plain array dummies; Fortran has no way to hand a
  struct of array SECTIONS across that boundary without pointer components. Whether the answer is 18
  dummies, a pointer-component view, or keeping a gather struct that no longer computes anything is a
  design question this plan does not answer.
- It changes the **TOML schema** (leaf width, branch diameter and crown fraction become PFT
  parameters, and `load_meds_config` hard-errors on a missing key), so every shipped PFT file and
  every test config changes with it.
- It changes **numbers**: §10.1 already says the three column tests get new golden values because
  their fixtures become allometrically consistent, and moving the seed/clamp out of the fast gather
  into a slow-loop `reconcile_tissue_water_capacity` changes when a mass edit happens.

Bundling that with steps 1-6 would also destroy their verification. The acceptance criterion for a
pure move is byte-identical output; a PR that legitimately changes numbers cannot demonstrate it. The
ordering argument for putting step 0 first ("so every later move commit is a pure `git mv`") has
already been satisfied -- the moves are done and every one of them verified byte-identical -- so
step 0 now simply runs on the new tree, at the cost of spanning two folders instead of one.

**Recommendation:** step 0 as its own PR, after this one, with §10.1's five bullets split into at
least two commits (the geometry/PFT-parameter change, which is data-identical for the production
path, and the fast/slow slice policy table, which is the extensibility payload).

### 11.7 Incidental finding, not part of the reorg

`meds_config_main.toml` cannot be run with `[output].enabled = true`: the tier subsections write
`enabled = false ; interval_steps = 4 ; file_chunk = "day"` on one line, and `meds_toml` does not
accept `;` as a key separator, so it hard-errors on `output.fast.enabled`. Nothing reads those lines
while output is off, which is why it has never fired. Worth a one-line fix to the shipped config.


---

## 12. What implementing steps 9 and 10 changed about this plan (2026-09-09)

Facade normalization and the renames landed as one PR, because a facade deletion IS a naming
decision (which spelling survives) and because the field sweep touches ~540 occurrences: doing it
after step 0/7 would rewrite lines those steps had just authored. The plan scheduled renames last so
the PATHS would be final; paths became final at step 6, and the constraint then inverts.

Four corrections.

### 12.1 The leaf selectors cannot live in `meds_plant_types`

Decision #10's last item said `COLIM_*`/`SM_*` "belong next to the table fields that hold them, in
`meds_plant_types`". Not reachable: `meds_config` must see them to load and validate them, and
`meds_plant_types` sits ABOVE the config layer because it reads the PFT trait table. They went to
**`meds_leaf_opts`**, a low-level config leaf beside `meds_biophysics_opts` and `meds_biogeochem_opts`
— the same intent one layer down. The payoff is that `meds_leaf_gas_exchange` no longer imports
`meds_config` at all, which is what decision #10 was actually for.

### 12.2 `leaf_gas_exchange_batch` was never a facade concern

§4 note 3 listed the facade's contents as "four procedures, two of which take `meds_config_t`", and
decision #10 sent them all to the fast driver. Three of the four take `meds_config_t` and are driver
code; the fourth, `leaf_gas_exchange_batch`, takes a `leaf_photo_table_t` and no config, so it is the
bare-array leaf KERNEL and it belongs in `meds_leaf_gas_exchange` beside `solve_leaf_gas_exchange`.
The cfg-taking three, plus `build_tol_set`/`build_error_control`/`build_integrator_opts` from
`meds_fast_control`, are now **`fast_dynamics/driver/meds_fast_config`** — the one place
`meds_config_t` is read on the fast path.

### 12.3 Deleting the plant facade removed a library edge

Recorded in §11.3 as an artefact to be cleaned up later; it is cleaned up now.
`meds_fast_kernels` linked `meds_slow_kernels` only because the facade re-exported phenology and
allocation. With the facade gone the fast kernels import nothing from `slow_dynamics/`, and the edge
is deleted from CMake. Tests that genuinely span both tiers (the plant kernel loop, the CAS-CO2 and
soil-biogeochem tests, the plant C-API shim) now say so in their own link lines instead of riding it.

### 12.4 `g_cas_atm_*` was not worth its width

§10.5 proposed `gah`/`gaw`/`gac` → `g_cas_atm_heat/vapour/co2`. Naming both endpoints reads well in
prose and cost **32 over-length lines** in code, because these three appear together on argument
lists and in the CAS box solve. They took `g_atm_heat`/`g_atm_vapour`/`g_atm_co2`, matching the
`_atm` suffix the sibling fields in the same records already use for the atmospheric side
(`enthalpy_atm`, `shv_atm`, `co2_atm`). The CAS side is stated by the record that holds them.

Everything else in §10.5 landed as written, plus the three-way `hydro` split and the
`t_precip`/`rain_temp`/`film_u_ref` question, which §10.5 left open as a naming problem:
`t_pond_inflow`, `t_film_valuation`, `film_liquid_enthalpy`. `rain_temp` was the dangerous one — under
a snow pack it is deliberately NOT the rain's temperature but `tsupercool_liq`, so meltwater carries
zero enthalpy and the pack's latent heat is not double-counted.

### 12.5 Renames are verified by absence of string change

Every rename commit was checked mechanically for **new string literals**: none introduced one. That
is the evidence for §10.5's exclusion ("netCDF registry strings and TOML keys are not renamed without
a compatibility note") — output files and shipped configs are untouched. The 132-column limit holds
tree-wide; ~70 lines were rewrapped across the sweep.

### 12.6 `meds_biophysics_interface` stays

Rule 6 bans a facade that is re-export AND logic, not a facade as such. This one is 49 lines of pure
re-export with a single consumer that never bypasses it, so it is left alone. It is the only facade
left in the tree.


---

## 13. What implementing step 0 changed about this plan (2026-09-09)

### 13.1 The reference run had no light, and therefore no growth

The byte-identity harness used through PRs #125 and #126 ran with `[forcing].forcing_on = false`.
That means **GPP was identically zero for the whole run**: no cohort ever grew past the recruit
size (every cohort sat at dbh 0.4534 cm for three simulated years), and nothing that depends on a
size CHANGE could be exercised at all.

For steps 1-6, 9 and 10 that is not a defect in the conclusion -- those steps are file moves,
module splits and renames, and a deterministic run is a valid witness for them; the suite and the
conservation ledgers covered the rest. But it is a much narrower witness than "byte-identical in
all 75 outputs" sounds, and it could not have caught a stale-cache bug, which is exactly what
step 0 risks. **The harness now runs real ERA5-Land forcing over a recycled year** (built with
`scripts/prep_era5land_forcing.py`), so trees grow, and it compares against `main` in a git
worktree rather than against a saved snapshot.

Two lessons for the steps still open:
- Byte-identity of a run that exercises nothing is not evidence. Before trusting a comparison,
  perturb the thing you changed and confirm the harness notices. (Perturbing `sapwood_carbon` by
  2x moves 51 of 75 output files; that is what makes the null result meaningful.)
- A conservation ledger cannot see a stale cache. Ledgers closed to machine precision throughout,
  in both the correct and the deliberately-broken builds.

### 13.2 Cache the geometry PER PLANT, not per ground

§10.1 says the derived geometry "becomes cohort-block fields", listing `lai` and `wai`. Cache those
and every mortality step makes them stale, because `nplant` changes without any geometry changing.
`dbh_to_wai` is exactly linear in `nplant`, so the per-plant quantity exists: the block caches
`wood_area`, `sapwood_carbon` and `sapwood_area`, all per plant, and the per-ground index is formed
as `nplant*area` at the point of use, mirroring `LAI = nplant*leaf_area`. Nothing stored carries a
plant density.

The remaining staleness is real and is handled explicitly: `update_cohort_states` advances dbh,
basal_area and wood_carbon by their own tendencies without re-deriving geometry, so it re-derives
the wood cache afterwards. Measured: without that, 51 of 75 output files differ.

### 13.3 Step 0 is NOT byte-identical, and that is correct

Moving `f_sap*wood_carbon` from the driver to the state module changes ifx code generation. First
divergence is **1.4e-12 relative**, at the first output tick after cohorts appear, growing over
three simulated years into a sub-percent trajectory difference whose largest survivor is a discrete
integrator counter. §8's acceptance table already provides for this ("data-identical at round-off");
the point worth adding is that the correct evidence is **the onset**, not the endpoint: a chaotic
coupled model will turn one ULP into a percent given enough time, so quoting the final difference
says nothing about whether the change was faithful.

### 13.4 `column_cohort_t` keeps its shape; the FILLER is what was wrong

The author chose this over the pointer-view and bare-array alternatives, and the evidence supports
it. Once 0a removed the computation, the gather is 18 plain copies, and the extensibility argument
that motivated deleting the type does not actually favour deleting it: adding a per-cohort input
costs three adjacent edits either way. What §10.1 correctly identified is the *fixtures*:

- the hand-built test views were not allometrically consistent (dbh 20 cm with a 16 m height and a
  leaf area of 10 m2/plant against an allometric 134), and
- **`bwood` was allocated by `alloc_column_cohort` and initialized nowhere**, so the wood heat
  capacity in three column tests ran on uninitialized memory. That is a live defect, now fixed.

So there is one filler with two entry points (`meds_column_view`), and the fixture entry builds a
cohort block through the canonical birth path before gathering it, which makes a fixture tree
on-allometry by construction. The old fixtures also implied 3000 stems/ha of 20 cm trees; made
consistent at that density the stand has LAI 40 and intercepts all rain, which broke two tests on a
forest that cannot exist. They are 224 stems/ha now.

### 13.5 Still open in step 0

The fast/slow slice components and their per-field policy table (§10.1 bullet 3), and moving the
lazy PSI_INIT seed and `clamp_water_to_capacity` out of the fast gather into a slow-loop
`reconcile_tissue_water_capacity` (§10.1 bullet 5). Both are independent of what has landed.


### 13.6 Step 7, and three §10.3 items that no longer exist

`meds_fast_ark` is split (1580 -> 629 lines): `meds_fast_frozen` takes the frozen-record builder
that RK45 also uses, `meds_fast_be_stage` takes the implicit stage, its Newton and Jacobian, and
the water-mass and canopy-film advance that the RK4 oracle also uses. Nothing outside the module
imports `meds_fast_ark` now except the dispatcher and the RHS test, which is what §8 step 7 asked
for. The move was byte-identical, which the plan does not promise for procedure moves.

Three of §10.3's leftovers are stale and should be struck: `relieve_theta_bounds`,
`soil_layer_temp`, `seed_soil_column` and `seed_plant_water` do not exist anywhere in `src/` or
`test/` any more -- the review PRs removed them. §10.3's `zero_like` claim is stale in the same
way (it allocates the films).

What §10.3 leaves genuinely open, and why it is not done here:

- **Pass `column_params_t` via `column_config_t` instead of copying it into the frozen record.**
  Independent of the split and worth doing; it is a signature change across both marches and the
  oracle, and this PR is already large.
- **Item 1A (vi)/(vii) ledgers.** Per-layer and per-cohort residuals asserted after the rail
  decision. These are ledger work and belong with §10.2, not with a refactor.

The silent-omission matrix is addressed as far as Fortran allows: `test_state_combinators` fills
every `column_state_t` field with a distinct value and asserts each combinator field by field,
including the fields the embedded-error estimate deliberately excludes. A true completeness check
is impossible -- Fortran cannot enumerate a derived type's components -- so the acceptance check
§10.3 proposes ("a regression test that adds a dummy field") cannot be written either. What this
does catch is an omission in an EXISTING combinator, which is the failure that has happened.


### 13.7 A second naming pass, and decision #3 finished (2026-09-09)

Review of the step-0/7 branch produced six more corrections. Four are naming; two are real.

**`src/shared/` is gone.** Decision #3 said it dissolves into base + functions + config + state.
Config left at step 5 and state at step 3, but `base/`, `functions/` and `util/` kept the prefix, so
the folder survived as a wrapper around three folders that already had names. They are top level now.
The LIBRARY keeps the name `meds_shared`, which is earned: it is a fact about the link graph, not a
folder for leftovers.

**`meds_fast_snow` is gone, and it should never have existed.** Its two symbols belonged elsewhere:
`snow_stage_t` is one of the nine content-named pieces of the frozen record and the other eight live
in `meds_fast_types`, and `advance_snow_stage` has exactly one caller, in `meds_fast_frozen`. Keeping
them apart also created a backwards edge -- a types module importing from a process module to obtain
one of its own components. Snow was never a separate concern from the other fast processes; it just
had a separate file. §4's tree lists it as a numerics module; strike it.

Names, all mechanical:

| was | is | why |
|---|---|---|
| `meds_tissue_water` | `meds_fast_reconcile` | named for the act, so future state repairs have a home. NOT `check`: everything in it writes |
| `meds_column_gather` | `meds_column_view` | "gather" is HPC jargon for what is a copy |
| `gather_column_cohort` | `copy_column_cohort` | ditto |
| `column_cohort_fixture` | `column_cohort_init` | it initializes a cohort block through the birth path |
| `meds_column_reservoirs` | `meds_column_state_types` | matches `meds_site_state_types` |
| `meds_column_constants` | folded into `meds_column_params` | two modules per state half, named alike |

**Where the column parameters live, settled.** `soil_params_t` stays in `state/column` rather than
moving to `src/config`: it is DERIVED per column from the `[soil]` scalars, not loaded, and stands to
them exactly as `leaf_photo_table_t` stands to the PFT table -- which lives with its kernel. Recorded
in the module header so it is not re-litigated.

### 13.8 Finding: the soil column is hard-coded, and config never sees it

Raised by the question "could the column parameters move to config like the PFT params". They cannot
move to config because **they are not in config at all**. `build_fast_context` calls
`build_soil_hydr_params` with literals:

```fortran
call build_soil_hydr_params(NSL_MVP, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp, &
                       2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ctx%col_config%soil)
```

Ten layers, two metres deep, a loam texture and a saturated conductivity, with the wood and root
respiration factors and the prescribed soil-carbon pool immediately below. The comment above it
admits this and names the follow-up. CLAUDE.md meanwhile states that the source defines only true
constants and every parameter is required from TOML.

This is the same defect just fixed for leaf width, branch diameter and crown fraction, at much larger
scale, and it has a second consequence: the hardwired 2.0 m soil depth is the one
`MEDS_SOIL_BOTTOM_THERMAL_BC` flags as shallower than the annual damping depth, so the fix for that
defect is currently unreachable from a config file.

**RESOLVED for the column itself (2026-09-09, branch `feature/soil-column-config`).** The fourteen
geometry, texture and thermal literals are a `[soil_column]` TOML block -- deliberately separate from
`[soil]`, which is the Richards SOLVER's options: one is the ground, the other is how it is solved.
Validated at load, defaults reproduce the literals, and `test_soil_column_config` asserts the keys
reach the column rather than merely compiling.

**Still open:** the five respiration / Rh / prescribed-pool literals in the same block. `agf_bs` is
the sharp one -- it silently shadows the per-PFT `aboveground_frac`, the same physical quantity, so a
user who differentiates PFTs by allocation gets demography using their values and stem respiration
using 0.7 (issue #128). `stem_resp_factor25`, `root_resp_factor25` and `is_woody` are per-PFT in ED2
and want traits, not global keys, so that is a PFT-table change rather than a config-block one.


---

## 14. What implementing step 8 changed about this plan (2026-09-10)

Step 8 is done: one `libmeds.so`, scikit-build-core, `netcdf4` declared, the shims split by subsystem,
and a mandatory ctest target for each. What the implementation changed:

### 14.1 Both things §7.6 #4 warns about had ALREADY happened

The rule exists because of issue #95 → #100 (a component inserted mid-type in `leaf_photo_params_t`
broke the C API while the suite stayed green). Starting step 8 found **two more instances of the same
class, live on main**:

- **`meds_demography_capi.f90` did not compile.** PR #137 changed `apply_recruitment`'s signature and
  the shim still passed the old rate array. It sat on main for a day through a green 42/42 suite on
  both back ends, because the only thing that built it was the optional `-DMEDS_BUILD_PYLIB=ON`
  library. **A shim nothing compiles rots.**
- **`examples/example_demography/example_config_pft.toml` was 37 required keys behind the schema**, so
  the example could not load its own config. **A config nothing loads rots the same way.**

So the fix is broader than the plan's: every shim gets a ctest target (§7.6 #4 as written), *and*
`test_capi_demography` runs from the source directory against the shipped example config, giving that
config a build-time consumer too.

### 14.2 FOUR subsystems are exposed, not six, and the demography shim is not split

§7.6 #3 lists six shims (`leaf`, `hydraulics`, `phenology`, `demography`, `fast`, `site`).
`hydraulics` and `fast` expose **no `bind(c)` entry points at all**, so creating them would invent API
surface. And `site` and `demography` share a module-`save` handle registry (`g_site` / `g_cfg` /
`site_used` / `cfg_used` / `g_generation`) that every entry point indexes; splitting them would hoist
that into a third module both `use`, which buys nothing at ~290 lines and gives the registry two places
to be got wrong. **Three shims: `meds_capi_leaf`, `meds_capi_phenology`, `meds_capi_demography`.**

### 14.3 Decision #11 is necessary but NOT sufficient

"Declare `netcdf4`, do not vendor" is right, and it does not make the wheel load. `libmeds.so` has
`NEEDED` entries for `libnetcdf.so.19`, HDF5 and the Fortran runtime; the pip `netCDF4` wheel bundles
its own libnetcdf under a mangled soname *inside its own package*, so it cannot satisfy a plain
`libnetcdf.so.19` for a different library. Freshly installed into a bare venv, the library failed on
**ten unresolved sonames**.

The fix is RPATH, and it takes two parts, because `INSTALL_RPATH_USE_LINK_PATH` alone is not enough:
it captures the directories of libraries CMake was *told* to link (netCDF, HDF5) but not the Fortran
runtime, which the compiler driver adds implicitly. `CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES` is
exactly that set. With both, the wheel loads in a bare venv with no `LD_LIBRARY_PATH` and no
`source setvars.sh`.

**The consequence, stated plainly: the wheel is MACHINE-LOCAL.** Its RPATH points at this machine's
conda prefix and oneAPI install. That is the right trade for `pip install python/` from a source
checkout, which is how MEDS is used; a redistributable manylinux wheel means gfortran + `auditwheel`
and the vendoring §7.3 measured and rejected.

**A measurement trap worth recording:** `ldd` run from a shell that has sourced `setvars.sh` inherits
`LD_LIBRARY_PATH` and reports everything resolved. It said "0 not found" while `dlopen` failed on
`libimf.so`. Check with `env -u LD_LIBRARY_PATH ldd`, or by importing.

### 14.4 Two dependencies the plan never listed

`numpy` is imported at **module level** by `meds/demography/_site.py`, so `from meds.demography import
Site` fails without it — but it was listed only under the `plot` extra. A clean-venv install of the
wheel caught it immediately; nothing before could, because the leaf and phenology sub-modules import no
third-party package at all. `dependencies = ["netCDF4", "numpy"]`.

### 14.5 A missing CMake edge that only Ninja exposed

`meds_config` uses `n_soil_layer_max` from `meds_column_params` (since the `[soil_column]` block
landed) but linked `meds_shared` only. Every everyday build uses the Makefile generator, which happened
to order the two `.mod` files correctly; scikit-build-core uses **Ninja**, which scheduled them
concurrently and failed. The edge is acyclic (`state/column` links `meds_shared` only) and is now
declared. **A build that passes under one generator is not a build that passes.**

### 14.6 The examples now prefer an INSTALLED package

They inserted `python/` on `sys.path` unconditionally, which shadowed an installed `meds` — so a wheel
could never be exercised by its own examples. The insert is a fallback now (`try: import meds / except
ImportError:`), and all four Python example entry points run against the installed wheel with **no
environment variables**, producing output byte-identical to the source-tree runs.

---

## 15. The remaining open items, phased (2026-09-11)

Migration steps 0–10 are all merged. What is left is **not migration** — it is the review remainder
of §10.3 and §10.1, plus two items found after §10 was written that it could not have anticipated.
This section is the plan for that remainder. It is written after a day (PRs #139, #140) in which two
real defects were found, both of a class the whole test suite was blind to, and the ordering below is
a direct consequence of what those two cost to find.

### 15.0 What has already closed — strike these from §10

The §10 lists are stale in four places. Recording it here so nobody re-plans finished work:

| §10 item | status | where |
|---|---|---|
| §10.2 slow-loop conservation (review 1B #4–#10) | **DONE** | PRs #132–#137 |
| §10.3 item 5 helpers: `relieve_theta_bounds` ×3, `soil_layer_temp` ×6 | **DONE** — neither identifier exists anywhere in `src/` or `test/` any more | — |
| §10.4 core facade | **DISSOLVED**, as decision #5 predicted. There is no `meds_demography_interface` and no `meds_core_interface`; `slow_dynamics/demography/` holds four operator modules and nothing else | PR #126 |
| soil-audit: `worst_rh_seam_gap` is dead output | **DONE** | PR #140 |

`worst_rh_seam_gap` is worth a sentence because wiring it up was not cosmetic. It had been computed
inside `advance_biogeochem_dynamics` and discarded on every step of every run since it was written,
because `meds_stepper` never asked for it. Within minutes of reporting it, it caught a state-carryover
bug in the brand-new C-API run registry. **A diagnostic nobody reads is not a diagnostic.** That is
the argument for Phase 1 below.

### 15.1 Ordering principle

Three rules, in priority order:

1. **Prevention before cleanup.** Both defects found on 2026-09-10 were *field/state bookkeeping*
   failures — a per-patch quantity that did not ride the patch lockstep, and a run object whose
   accumulators were not reset on reuse. The silent-omission matrix (Phase 4) is the same class,
   industrialised: about a dozen field-enumerating routines in `meds_column_state_ops` plus per-scheme
   enumerations in
   `meds_fast_rk45` (21 field references) and `meds_fast_ark` (12), none of which fails to compile
   when a field is forgotten. It is the single largest remaining defect *generator* in the tree.
2. **Cheap-and-loud before expensive-and-quiet.** A check that already exists and is merely unread
   (Phase 1) costs a day and can find something immediately. A refactor that prevents future defects
   (Phase 4) costs a week and finds nothing today. Do the former first — not because it matters more,
   but because it changes what you know before you commit the week.
3. **Signature churn rides together.** This is §10.3's own rule and it still holds: Phase 5's items
   all touch the march signatures, so they go with Phase 4, which rewrites them.

Each phase below is independently mergeable and carries its own acceptance test. **A phase without a
falsifiable acceptance test is not ready to start.**

### 15.2 Phase 1 — make the checks that already exist speak

Four numbers that already exist and no *run* ever looks at. Each is a one-day item. (Precision
matters here: three of them ARE asserted by unit tests on synthetic pools. What none of them has is a
consumer in a 50-year run, which is the regime where they would have something to say.)

- **`Lambda`, the freeze number** — `dvec(j) = xi_int(j)*k_diag(j)`, the fraction of pool `j` the
  slow step withdraws, computed on line 321 of `meds_soil_biogeochem` and discarded as soon as the
  step uses it. `soilc_audit_t` has seven fields and no Λ. This is the most valuable of the four and
  the reason it was promoted here out of Phase 3: it is strictly more informative than the seam gap.
  **The seam gap catches a broken contract; Λ catches the approximation degrading before anything
  breaks.** Measured at **3.6e-3 per day** (fast_grnd, mature Ithaca stand, July), so a run where it
  climbs toward 1 is telling you `dt_slow`
  is too long for that store, well before a pool goes negative.
- **`audit%lignin_resid`** — the lignin passive-tracer balance. Asserted in `test_soil_biogeochem`
  on synthetic pools, never reported from a run. It is NOT a tautology: it is what keeps
  `0 <= L <= C` under the EXPM solver, where using the Euler survival fraction drives lignin
  negative. Report it beside Λ and the seam gap.
- **`audit%resid`** — `dC_pool - (litter_in - rh_today)`, with `rh_today` defined as
  `litter_in - dC_pool`. Substitute and it is identically zero for any inputs. `test_soil_biogeochem`
  asserts it twice, labelled "reporting invariant", and **the genuinely independent check sits four
  lines below it** — recomputing Rh from `er*xi_int*K*X0` rather than from the closure definition,
  with a comment saying in so many words that `audit%resid` cannot catch what it catches. **Delete
  the field and the two vacuous assertions**; a check that cannot fail teaches readers to discount
  the ones that can.
- **`audit%litter_in`** is `sum(u)` from `build_litter_input`, i.e. the turnover + continuous-mortality
  channel ONLY. Cull-termination and disturbance-kill necromass are added directly onto
  `patch%soil_carbon` by the demography operators, so the audit's "litter in" is not the patch's
  litter in. Widening it is not cheap (the operators bypass `u` by design), so **rename it
  `litter_in_matrix`** — silently wrong naming on an audit field is worse than no field.

Also in Phase 1, because it is the same kind of work: **the seam caveat measured in PR #140.**
`seam[soil_carbon_rh]` is machine-zero (1.1e-14 kgC/m²) over a 31-day July but reaches **8.370e-4
kgC/m² over the 50-year spin-up, at 2073-01-01 — a year rollover**, when the annual patch cadence
fires. So the "equal by construction" contract holds exactly *except* on days when patch structure
changes between the fast window and the slow step: the fast loop accumulates `rh_fast_accum` against
one patch composition, and the daily step debits a different one, with `blend_xi_accum` and
`blend_soil_carbon` area-weighting the two ends separately through a matrix that is nonlinear in the
lignin fraction. Decide between: (a) accept and document the exception, (b) mark the seam
unmeasurable on structural-change days rather than reporting a number that is expected to be nonzero.
Do NOT "fix" it by widening the tolerance.

**What Phase 1 deliberately does NOT add:** a debit-before-credit assertion on the two *rate* seams
(`shed_water_rate`, `slow_co2_rate`). Λ is meaningless for them — there is no store in the
denominator — and the invariant that makes them sound is that the slow tier removes the water/carbon
from its own store before handing the fast loop a rate to deliver. **The slow ledger already tests
exactly that**: both are declared as handoffs, and a credit without its debit shows up as a phase
residual. That is how the discarded mortality water was found in the first place. Adding a second
check over the same transfer would be redundant; §15.4 records the classification instead.

**Acceptance:** a deliberate sign error in the CENTURY lignin update is caught by a reported number,
not by reading the code; and `Lambda` appears in the run output with its measured value.

### 15.3 Phase 2 — the defaults that gate known-wrong physics

`feedback_delete_flags_that_gate_wrong_physics` already settled the principle: a config flag whose OFF
path is known-wrong physics should be deleted, not plumbed. Three flags are currently in that state,
and all three are defaulted to the wrong side.

| flag | default | the problem |
|---|---|---|
| ~~`energy.phase_change`~~ | — | **RETIRED 2026-09-11.** The premise was wrong (see the measurements below): the plateau was never gated, and the flag only selected liquid-only conductivity and heat capacity. No use case survived the measurement — ≤0.61 K, zero wall-clock cost, *identical* solver work — so the flag was deleted rather than re-defaulted, per `feedback_delete_flags_that_gate_wrong_physics`. A config still carrying the key is now a hard error rather than a silent no-op. |
| `soil_column.depth` | `2.0` m | Against a ~2.5 m annual damping depth, so the annual wave reflects off a zero-flux base. Measured 2→3 m: base layer **−8.5 K**. Reachable since the `[soil_column]` block landed; no shipped config sets it past 2.0. |
| `soil_carbon.soil_carbon_on` | `false` | Off is not a coarser soil model, it is *no* soil carbon. This default is what kept two defects (PR #139's out-of-bounds litter read, PR #140's 964× Rh units error) out of every code path anyone ran, for as long as they existed. |

This phase is **science, not structure**, so it needs the author's decision per flag rather than a
default recommendation from the plan. What the plan can say is that the current state — a documented
"this is wrong" next to a default that selects it — is the worst of the three options.

#### Measured, 2026-09-11 — one full year through a winter

Four variants, all restarting from the same 2074-01-01 spin-up state and running 2074-01-01 →
2075-01-01 at `dt_fast = 900 s`, run concurrently under identical load so the wall times compare.

| | base (as shipped) | `phase_change = on` | `depth = 3.0` | `soil_carbon_on = false` |
|---|---|---|---|---|
| wall clock | 53.8 s | 53.2 s | 53.4 s | 53.1 s |
| soil T layer 1, annual min | 271.10 K | 271.47 K | 270.5 K | 271.10 K |
| annual swing at −1.73 m | 18.55 K | 18.55 K | **13.56 K** | 18.55 K |
| Rh, annual total | 0.833 kgC/m² | 0.833 | 0.820 | **0** |
| NEE, annual mean | −2.464 µmol/m²/s | −2.465 | −2.502 | **−4.657** |

**None of the three costs anything measurable in wall clock.** That removes the usual argument for a
cheap default, and it is the single most useful thing the measurement produced.

**`phase_change`: the premise was wrong, and that is a finding.** The ed2_comparison claim this phase
was built on — "freeze/thaw plateau implemented but opt-in, default off" — is **not what the code
does**. `internal_energy_to_temp` is called unconditionally and always inverts through the phase
change, because the column is prognostic in internal energy and temperature is a read-off. **The
plateau is always on.** The flag forces `fliq_use = 1.0` at exactly two sites
(`meds_soil_energy.f90:104` and `:182`), i.e. it evaluates conductivity and heat capacity as if all
the water were liquid while the temperature still shows the plateau. Both settings produce the same
zero-curtain signature (102 vs 105 layer-days with partial `fliq` at 0 °C) and the total difference
over the year is **≤0.61 K**. The doc has been corrected. The remaining decision is small and
cosmetic by comparison: turn on ice-aware κ/C by default because it is free and more nearly right,
or leave it and stop calling it a phase-change switch. **Resolved: DELETED** (2026-09-11). With
the plateau unconditional anyway, no scenario wanted liquid-only κ/C — it was P1 staging, not an
option — and the measurement closed the last argument for keeping the branch reachable: *identical*
solver work, 35 692 soil substeps either way. Deleting it reproduces the old `on` setting
**bit-identically** (max |ΔT| = 0.0 K over the verification year) and differs from the old default
by 0.611 K. A config still carrying the key is now a hard error, not a silent no-op.

**Phase 2 is CLOSED (2026-09-11).** `phase_change` deleted; `soil_carbon_on` defaulted to true
(the flag stays — unlike `phase_change`'s, its off branch is a legitimate diagnostic configuration,
and it is how this section's own numbers were produced); `soil_column.depth` deferred to **issue
#145** by decision, because the measurement showed it is not a default change at all.

**`soil_column.depth`: the real one, and the reason it became an issue rather than a config edit.** Comparing the two columns **at matched physical depths**
rather than at their own base layers is what makes it clear, and the error grows monotonically
toward the boundary — 2.2 K at −0.64 m, 3.2 K at −0.90 m, 4.1 K at −1.25 m, **5.0 K at −1.73 m**.
The 2 m column overstates the annual swing at its own lower third by ~37 %. And **3 m is not the
answer either**: its own base layer still swings 13.1 K, so the wave is not damped there either. The
decision is therefore not "2 → 3 m" but "how deep, or does the adiabatic bottom BC need replacing".
Fitting an e-folding depth to each column's own amplitude profile makes it plainest: **5.01 m (2 m
column) and 3.98 m (3 m column) against a physical ~2.0–2.5 m** — both far too slow, the shallower
one worse. Filed as issue #145 with the reproduction; the recommendation there is a Dirichlet
temperature anchor at the base rather than a deeper column, following the pattern the *water* BC
already sets with `free_drain | bedrock | aquifer`.

**`soil_carbon_on`:** annual-mean NEE reads **−4.657 vs −2.464 µmol/m²/s** with it off — the stand
looks like an 89 % stronger sink — and Rh is identically zero against 0.833 kgC/m²/yr. This is the
default that kept PR #139's and PR #140's defects out of every path anyone ran.

**Acceptance:** for each flag, a measured number for what the correct setting costs (wall clock and
the headline diagnostics), and either a changed default or a recorded decision not to change it.

### 15.4 Phase 3 — the frozen-seam contract (design only, no code)

The slow→fast seam freezes a slow state across the whole slow step while the fast loop computes
fluxes from it. That pattern is everywhere in MEDS (`soil_carbon` as a frozen store; `xi_accum` as
the return accumulator; `shed_water_rate` and `slow_co2_rate` as frozen daily rates) and it has a
sharp, checkable soundness criterion that is currently written down nowhere. Write it down.

For a frozen store `S` consumed by flux `F` over `dt`, the governing number is
`Lambda = F*dt/S` — the fraction of the store the step withdraws — and the case split is what `F`
does as `S` goes to zero:

- **`F = k*S` (linear in the store):** `Lambda = k*dt`, **independent of S**. The freeze is
  *scale-free*; near-bare ground is no more dangerous than a mature soil. This is the CENTURY case:
  `dvec = xi_int*k_diag` measures **3.6e-3 per day**, well under 1 by any margin that matters, and it
  is why bare-ground spin-up works.
- **`F` prescribed, or otherwise not vanishing with `S`:** `Lambda` diverges as `S` goes to zero and
  the freeze is not merely inaccurate, it is **unsound**. MEDS has a scar from exactly this: the old
  `soil_carbon_on = false` branch respired a *prescribed* 5 kgC/m² pool that did not exist, held the
  canopy air 47 ppm above ambient and drove site NEE to +30.7 µmol/m²/s.

The note should carry: the criterion above; a classification of all four existing seams by it; the
rule that source seams are safe because the slow tier **debits before it credits** (currently a
convention in comments, not an assertion); `Lambda` as a per-seam reported diagnostic, which is the
generalisation of `rh_seam_gap` and strictly better than it — the seam gap catches a broken contract,
`Lambda` catches a *degrading approximation* before it breaks anything; and the note that
multi-consumer stores (soil water: transpiration, evaporation, drainage, runoff) need *arbitration*
(scale all demands by `min(1, S/D)`) rather than per-process clamping, which is order-dependent.

**The Λ instrument itself moved to Phase 1** (it is a number the code already computes, so it belongs
with the other unread ones). What stays here is the *reasoning*: why Λ is the right instrument for a
store seam and the wrong one for a rate seam, and how to classify the next seam somebody adds. That
makes this phase genuinely optional — the instruments exist whether or not the note gets written —
which is the right status for a note whose value is conceptual.

Explicitly **not** recommended: making the seam implicit. `SOIL_LIN_PICARD` exists, but paying an
outer iteration to fix a 3.6e-3 error is the wrong trade.

**Acceptance:** a design note in `docs/dev_plans/`. No code. It exists to make Phase 1's seam decision
and any future slow→fast coupling a lookup rather than a re-derivation.

### 15.5 Phase 4 — the silent-omission matrix (§10.3)

The largest remaining item, and the one most likely to prevent the next bug. Every `column_state_t`
field is enumerated independently in about a dozen routines (`state_init`, `state_axpy`, `state_accum`,
`state_extrap`, `state_err_diff`, `state_sub`, `zero_like`, the three clamps, `unpack_column_state`)
plus the per-scheme sites in `meds_fast_rk45` and `meds_fast_ark`, and **omitting a field anywhere
compiles clean**. Two known live instances are recorded in §10.3: `zero_like` never allocates the
films, and `state_err_diff` leaves the pond default-initialised while `state_sub` subtracts it.

Give the state vector one field-visitor (`state_visit`) or one contiguous packed layout, and route
the combinators, the process mask, pack/unpack and the error norm through it.

**Acceptance is already stated in §10.3 and should be held to literally: add a dummy field to
`column_state_t` and the build fails, or a test fails.** Not "a reviewer would notice".

### 15.6 Phase 5 — fast-integrator structural leftovers (rides Phase 4)

These go with Phase 4 because they touch the same signatures:

- **§10.3a** — pass `column_params_t` via `column_config_t` instead of copying it into
  `column_frozen_t` (`meds_fast_types.f90:639`). Deferred by the review precisely for this churn.
- **§10.3c / review item 1A (vi)–(vii)** — per-layer `budget_imbalance` over faces on the *committed*
  path, per-cohort tissue water/energy residuals, and **RK45 ledgers asserted after the rail decision,
  not on a step the dispatcher may roll back**. `rk45_state_railed` and the hybrid rescue already
  exist, so this is placement, not new machinery.
- **§10.1** — `column_cohort_t` out, fast/slow slice types in. Still entirely open;
  `column_cohort_t` is alive in `meds_fast_prepass` and elsewhere.

### 15.7 Phase 6 — test-support consolidation

**14 of 46 test files define their own `check`.** `test/meds_test_support.f90` has `check`,
`check_close` and `banner`; the duplicates predate it. Mechanical, no risk, and it makes every future
test's failure output uniform. Last because it prevents nothing — it is genuine hygiene, and should
not be allowed to displace Phases 1–4.

### 15.8 What is deliberately NOT on this list

- **Anything found by re-reading rather than re-measuring.** Both PR #140 corrections to this
  repository's own documentation (a Reco comparison that compared two spin-up states rather than two
  respiration limbs; an agreement figure quoted for "one July" that was measured over one *day*) were
  caught by re-running, not by proofreading. A plan item that cannot be measured is a note, not a
  phase.
- **Bit-identity between the Python driver and the executable.** Measured, explained (glibc `libm`
  interposes on Intel `libimf` inside a `dlopen`ed library; `LD_PRELOAD=libimf.so` reproduces the
  executable byte for byte), and not worth a link-time fix: `-static-intel` would change the RPATH
  story §14 settled, for a round-off difference.
- **The `dt_fast` accuracy work, the ground-conductance refresh, and the soil-water arbitration.**
  All real, all tracked in their own plans (`MEDS_NUMERICS_SCOPING.md`,
  `MEDS_PRODUCTION_INTEGRATOR`). This section is the *structure-plan* remainder and should not grow
  into a general backlog.
