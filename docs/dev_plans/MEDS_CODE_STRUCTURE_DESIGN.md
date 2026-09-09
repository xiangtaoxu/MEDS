# MEDS source-tree structure — reorganization plan

**Status:** **steps 1-6 MERGED (PR #125); steps 9, 10 and part of 8 MERGED (PR #126). Steps 0 and
most of 7 IMPLEMENTED** on `refactor/fast-loop-state-vector`. Remaining: §10.2 slow-loop
conservation (physics, always outside this plan's scope), the rest of §8 step 8 (Python
packaging), and two §10.3 items that need the ledger. §13 records what step 0 and step 7 changed
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
| **8** **PARTIAL** (#12 done; #11 + §7.6 packaging open) | Python: decisions #11, #12 + §7.6 | small | the two real costs |
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

Physics, so outside this plan's "no numerics" scope — but the *placement* is this plan's
business, and the order matters: these land after step 6 so they are written once into their
final home instead of being moved a week later.

- `slow_dynamics/driver/meds_slow_ledger` — `site_ledger_t` snapshot and daily/annual assert for carbon,
  water and energy across the slow step, with seam events (recruitment seed rain, tissue-water
  reconcile, mortality/cull water and heat, disturbance) as **declared** boundary terms. It
  reuses `budget_t` from `base/`. This is the fifth entry of the review's original fix order and
  the largest open physics item; nothing else in this section is checkable without it.
- `slow_dynamics/plant/`: recruitment carbon debit vs `init_cohort` endowment (1B #4); starvation
  `deficit` and `growth_resp` routed into the CAS carbon balance (1B #5).
- `slow_dynamics/demography/`: mortality, cull and disturbance-kill hand tissue water, film water and
  tissue heat to a declared sink instead of discarding them (1B #7); `terminate_patches` survivor
  rescale and `blend_cas` depth blending revisited (1B #9); fuse/split asserts at ~1e-12 relative,
  `conservation_tol` reserved for the site ledger, `size_tol` deleted (1B #10).
- `slow_dynamics/soil/`: CENTURY transfer-matrix column conservation asserted; `rh_seam_gap` asserted in
  production; `audit%litter_in` includes cull/disturbance litter (1B #8, #6).
- Not structure at all, listed so it is not lost: `docs/ed2_comparison.md` soil-temperature
  statements predate the PR #119 root-heat-sink fix and must be re-run.

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
`project_meds_soil_bottom_thermal_bc` flags as shallower than the annual damping depth, so the fix
for that defect is currently unreachable from a config file. A `[soil]` TOML block is its own piece of
work -- it changes the schema -- and is not part of this plan.
