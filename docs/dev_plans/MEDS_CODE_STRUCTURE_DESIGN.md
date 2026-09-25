# MEDS source-tree structure — decisions, rules, and what is left

> # ✅ LIVE — split and status-rewritten 2026-09-13.
>
> **Migration steps 0–10 are all merged** (PRs #125, #126, #127, #138, with §15 Phase 1 in #141,
> Phase 2 in #143 and Phase 4 in #147). The tree is timescale-first for processes over a state
> layer in two halves, and every step was verified byte-identical on both back ends.
>
> **This file was 1894 lines and is now the decisions, the rules and the open work.** The record of
> *how* the migration ran moved to `archive/MEDS_CODE_STRUCTURE_MIGRATION_LOG.md` (the diagnosis,
> the target-tree specification, the migration order, and the four "what implementing this changed"
> sections), and the 549-line slow-loop conservation measurement log that was §10.2 moved to
> `archive/MEDS_SLOW_LOOP_CONSERVATION_LEDGER.md`.
>
> **Section numbers are unchanged across all three files**, because roughly forty source comments
> cite them by number (`§7.6 #3`, `§10.2.11`, `§15.6`). The numbering is therefore not contiguous in
> any one file; that is deliberate.
>
> **For the tree itself, read [`src/README.md`](../../src/README.md)** — the live description of the
> layout, the placement rules and the library graph, kept beside the code it describes. §4 of this
> plan specified that tree before it was built and has since drifted (nine module renames, seven
> modules added); it is in the migration log as the specification, not as a description.
>
> **Open items are §15 below**, and are also indexed in [`docs/ROADMAP.md`](../ROADMAP.md) §10.

**Scope.** File and folder placement, library layering, and the Python package surface. **No
physics, no numerics, no algorithm changes.** Steps 1–6 were `git mv` plus CMake target edits plus
`use`-line renames; the only behavioural change in the whole plan was decision **#12** (allometry
defaults), which is a bug fix the reorganization made reachable.

**Supersedes two placements** locked by `archive/MEDS_REORG_DESIGN.md` (v6):

- *"column types: `meds_column_state_types` STAYS in `shared/state/` — it's the BOUNDARY type"*
- *"each domain kernel library links `meds_shared` ONLY"*

Both were correct **given** the two-shared-library Python split. That split is retired by decision
**#1** below, which removes the constraint that forced state into `shared`.

---

## 1. Decisions locked

| # | Decision |
|---|----------|
| **#1** | **ONE Python shared library.** `libmeds_plant_c` + `libmeds_c` → **`libmeds.so`**. Per-subsystem independence moves from the *link* level to the *call* level (verbs, not libraries). The per-domain STATIC libs all stay. |
| **#2** | **State becomes a LAYER, in two halves.** `src/state/column/` (per-patch reservoirs + params + fusion blends — the boundary types kernels legitimately need) and `src/state/site/` (cohort SoA, patch CSR, lockstep, `site_t`, diag blocks — drivers only). Kernels link `state/column` and stay `site_t`-free. Two folders, not two file prefixes in one folder: every library in the tree is a per-folder GLOB, and the folder boundary *is* the link guard. **Naming (2026-09-09):** `column` is kept although a column is exactly one patch's vertical profile — `patch/` was rejected because `patch_block` (the CSR container of all patches) lives in the *site* half, so a `patch/` folder would point readers at the wrong half. The `state/column` README states the synonym once. **Grid-ready:** a future `state/grid/` (`grid_t` = array of `site_t` + geolocation + the met-forcing handle) is a third folder, not a rework: `met_driver_t` is already passed beside `site` (the stepper takes both), `meds_forcing` links shared + netCDF only and never `site_t`, so `grid_t` can own the forcing with no DAG cycle. Do **not** create an empty `state/grid/` now. |
| **#3** | **`shared/` dissolves** into `base/` + `functions/` + `config/` + `state/`. The word "shared" stops being a place where things go when no other place fits. **Revised 2026-09-25:** `base/` + `functions/` + `util/` regroup as **`src/shared/`**, because together they *are* the `meds_shared` library, which was the only library whose sources spanned three top-level folders; folder = library again. `config/` and `state/` stay top-level. What made the old `shared/` a catch-all was the two-`.so` link rule that forced `config/` and `state/` into it, and decision **#1** retired that rule. The regrouped folder admits only what `use`s nothing outside it (no model state, no configuration, no external library), a test a reviewer can check against the `use` lines. |
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
| ~~`soil_carbon.soil_carbon_on`~~ | — | **DEFAULT FLIPPED TO `true` 2026-09-11.** Off is not a coarser soil model, it is *no* soil carbon. It cost nothing in wall clock, reported an 89 % stronger apparent sink, and kept two defects (PR #139's out-of-bounds litter read, PR #140's 964× Rh units error) out of every code path anyone ran for as long as they existed. The flag stays — off is a legitimate *diagnostic* configuration, unlike `phase_change`'s off branch — but it is no longer what a user gets by accident. |

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

#### Executed 2026-09-11 — and the acceptance criterion above is NOT achievable as written

Running the acceptance test *first*, before building anything, changed the phase:

1. **The hazard is real.** Adding two probe fields to `column_state_t` — a scalar *and* an
   allocatable — compiles clean and passes 45/45 with no combinator touching either.
2. **Both live instances §10.3 recorded are already fixed.** `zero_like` allocates all four
   arrays; `state_err_diff` assigns the pond explicitly with a comment. So this is prevention, not
   repair — which lowers its priority against §15.1's own ranking.
3. **Compile-time enforcement needs the packed layout, and that is disproportionate.** Fortran has
   no reflection. The one cheap trick that might have served — a positional structure constructor
   in a test as a field-count tripwire — **does not work**: ifx accepts a constructor with a
   trailing default-initialised component missing (tested). Real enforcement means
   `y%cas_enthalpy` → `y%v(I_CAS_ENTHALPY)` throughout: **1 207 field references** across `src/`
   and `test/`, reaching into the physics kernels. Filed rather than done, with that number
   attached.

**What was done instead (the proportionate version).** `test_state_combinators` already covered
five combinators. The gap was the three that produce or propagate a whole state and had no test —
`state_accum`, `state_extrap`, and **`unpack_column_state`, the commit path**. They are covered now.

**And the gap was hiding something worth naming.** None of the three had a defect, but two of them
— `state_accum` and `unpack_column_state` — write 9 of the 11 fields and silently omit the pond.
Both omissions are *deliberate* (the pond has no stage RHS and is committed from the scratch
hydrology solve), and neither said so. **A deliberate exclusion that is indistinguishable from a
forgotten field is the condition under which the next real omission hides** — which is exactly what
`state_err_diff` had already recognised and marked. Both sites now state the exclusion, and the
test asserts it: `state_accum` must leave the pond *unchanged* (not zeroed), `unpack_column_state`
must not commit it. All three new coverage areas were mutation-tested.

**Residual risk, stated honestly:** a field added to `column_state_t` and wired nowhere is still
invisible. It is also inert. The failure this now catches is the realistic one — an *existing*
field dropped from a combinator, which is what has actually happened here before.

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
