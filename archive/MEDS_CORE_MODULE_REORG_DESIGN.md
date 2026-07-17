# MEDS `core` module reorganization — implementation plan

**Status:** design-only (2026-07-16). Supersedes the naming decision in
`MEDS_REORG_DESIGN.md` (v6 "keep the name `demography`, NOT `core`") — the user has
now decided to **rename the demography engine to `core`**. Extends v6; does not undo the
v6 work already shipped (state merged into the engine library, `meds_ecosystem_state`,
allometry in `shared/functions`, empirical laws → Python, carbon rates → `meds_plant_vital_rates`).

This plan is the outcome of an evaluated 5-proposal reorg (P1–P5) plus four follow-up
decisions. It is behaviour-preserving except for **P1** (a numerically-identical-to-ULP
tendency reformulation) and the **P3** birth-field divergence fix.

---

## 1. Decisions locked

| # | Decision |
|---|----------|
| **Library rename** | `libmeds_demography` → **`libmeds_core`**; directory `src/demography/` → **`src/core/`**; all modules `meds_demography_*`/`meds_ecosystem_state` → **`meds_core_*`**. |
| **P1** | Fold the three cohort-mutators into a single **`update_cohort_states`** (in `core`) that APPLIES a **transient per-cohort tendency bundle** (`cohort_deriv_block`). The tendencies (all field derivatives + `dln_nplant_dt`) are COMPUTED in the driver `meds_vegetation_dynamics` (and the empirical path in the capi) — the old `growth_step`/`apply_growth` arithmetic + ring-buffer + mode branch move OUT of `core`. `core` stays mode-agnostic. |
| **P2** | Add a placeholder **`update_patch_states`** (in `core`) applying patch-level tendencies — patch aging now; the per-patch `soil_carbon_t` step is deferred to P3 (its CENTURY-kernel call will live in the driver, not `core`). No fast-loop consolidation. |
| **P3** | Extract **`init_cohort`** (in `core`) — the ~19-field birth block duplicated between `apply_recruitment` and `meds_init%add_cohort`. Resolves the `overtopping_lai`-at-birth divergence (zero it in the helper). |
| **P4** | Move the **pure output reductions** (`meds_demography_diagnostics`) UP into the io preparation library, renamed **`meds_output_diagnostics.f90`**. Keep `compute_overtopping_lai` in `core` (it writes a derived STATE field and feeds growth) — folded into `meds_core_state_update` as **`update_overtopping_lai`**. |
| **4-file core** | The structural core is **4 files** — `meds_core_state_types`, `meds_core_state_update`, `meds_core_cohort_fusefiss`, `meds_core_patch_fusefiss` — plus the public façade **`meds_core_interface`** (renamed from `meds_demography_interface`). |
| **io targets** | No subfolders. Rename CMake targets: `meds_output_core` → **`meds_io_prep`**; `meds_io` → **`meds_io_stream`**. `meds_netcdf_c` and `meds_config_io` unchanged (shared binding / config reader — neither is prep nor stream). |
| **`core` word** | Because the engine is now `core`, the io aggregation half is **`meds_io_prep`** (NOT `meds_io_core`) so "core" means exactly one thing. |

**Load-bearing invariants preserved:** the acyclic library DAG (`meds_core` links
`meds_shared` only); OpenMP-`target` device eligibility (bare-array + firstprivate-scalar
kernels); AGB/plant-number conservation on fuse/split; the single lockstep reorder machinery;
persistent `global_id`; the two-timescale wall (fast biophysics ⟂ slow `core`).

---

## 2. Naming map

### Library / directory / target
| Old | New |
|-----|-----|
| `src/demography/` (dir) | `src/core/` |
| `add_library(meds_demography …)` glob `src/demography/*.f90` | `add_library(meds_core …)` glob `src/core/*.f90` |
| `target_compile_options(meds_demography PUBLIC ${MEDS_MP_FLAGS})` | `…(meds_core …)` — offload flags follow the engine |
| target `meds_output_core` (`src/io/meds_output_{types,integrate,registry}.f90`) | target **`meds_io_prep`** (+ `meds_output_diagnostics.f90`) |
| target `meds_io` (`src/io/meds_io.f90`, `meds_output_{stream,manager}.f90`) | target **`meds_io_stream`** |

### Modules / files
| Old file (module) | New file (module) | Fate |
|-------------------|-------------------|------|
| `meds_ecosystem_state.f90` (`meds_ecosystem_state`) | `meds_core_state_types.f90` (`meds_core_state_types`) | rename + gains `init_cohort`, `cohort_deriv_block` |
| `meds_demography_dynamics.f90` (`meds_demography_dynamics`) | **dissolved** | `apply_recruitment`→cohort_fusefiss; `apply_patch_disturbance`→patch_fusefiss; `growth_step`/`mortality_step`/`apply_growth` arithmetic → driver; apply → `state_update` |
| `meds_demography_fusefiss.f90` (`meds_demography_fusefiss`) | **split** → `meds_core_cohort_fusefiss.f90` + `meds_core_patch_fusefiss.f90` | |
| `meds_competition.f90` (`meds_competition`) | **dissolved** → `update_overtopping_lai` in `meds_core_state_update` | |
| `meds_demography_diagnostics.f90` (`meds_demography_diagnostics`) | `src/io/meds_output_diagnostics.f90` (`meds_output_diagnostics`) | move to `meds_io_prep` |
| `meds_demography_interface.f90` (`meds_demography_interface`) | `meds_core_interface.f90` (`meds_core_interface`) | rename + updated re-exports |
| — (new) | `meds_core_state_update.f90` (`meds_core_state_update`) | `update_cohort_states`, `update_patch_states`, `update_overtopping_lai` |

---

## 3. The 4-file core + interface — detailed design

The core library DAG is acyclic **only** if all low-level primitives stay at the bottom:

```
meds_shared
   └── meds_core_state_types      (types + primitives + init + CSR bridge)
          ├── meds_core_state_update      (apply/update derived state)
          ├── meds_core_cohort_fusefiss   (cohort restructuring + recruit)
          └── meds_core_patch_fusefiss    (patch restructuring + disturbance)  ── may use cohort_fusefiss
                 └── meds_core_interface   (public façade, re-exports)
```

**The one rule that prevents a compile cycle:** `cohort_reorder`, `rebuild_csr`,
`copy_cohort_slot`, `cohort_compact`, `set_cohort_size*`, `assign_*_id` MUST stay in
`meds_core_state_types`. `rebuild_csr` calls `cohort_reorder`; if `cohort_reorder` were
lifted into `cohort_fusefiss`, and `cohort_fusefiss` calls back into state_types
(`assign_cohort_id`), the two modules would `use` each other. Keep every primitive at the bottom.

### 3.1 `meds_core_state_types` — types + primitives (bottom layer)
- **Types:** `cohort_block`, `patch_index`, `site_t`, `carbon_flux_block`, and the NEW
  **`cohort_deriv_block`** (transient tendency SoA; §4).
- **Sentinels:** `GROWTH_AVG_UNSET`, `PHENOLOGY_STATUS_INIT` (defined here today).
- **Lockstep + CSR primitives** (the "one place that touches every array"): `site_alloc`,
  `site_free`, `cohort_ensure_capacity`, `cohort_reorder`, `cohort_compact`,
  `copy_cohort_slot`, `rebuild_csr`, `patch_ensure_capacity`, `gather_pft_params`,
  `set_cohort_size`, `set_cohort_size_from_carbon`, `assign_cohort_id`, `assign_patch_id`.
- **NEW `init_cohort`** (P3, §6): `init_cohort(cohort, m, pft_table, ipft, owner_patch, nplant, dbh)`
  — per-slot field assignment + forward `set_cohort_size`, zeroing `overtopping_lai`. Does
  NOT ensure capacity / assign id / rebuild CSR (callers differ: single-append vs batch).

### 3.2 `meds_core_state_update` — apply/update derived state (middle layer)
Pure `state ← state ⊕ tendency` and derived-state refresh; `use`s state_types only.
- **`update_cohort_states(cohort, deriv, dt)`** (P1, §4): applies the tendency bundle —
  additive for geometry/pools, multiplicative-in-log for `nplant`; defensive nonneg floor on
  pools; density floor → mark-for-cull. Mode-agnostic, offloadable (pure array kernel).
- **`update_patch_states(patch, …, dt)`** (P2, §5): applies patch tendencies — `age += dt`
  now; a `soil_carbon` block later (its tendency computed driver-side, applied here).
- **`update_overtopping_lai(site)`** (P4): the folded `compute_overtopping_lai` — a
  height-sorted, per-patch top-down cumulative-LAI sweep that writes `cohort%overtopping_lai`.
  Must run AFTER `sort_cohorts` (sequenced by the driver).

### 3.3 `meds_core_cohort_fusefiss` — cohort restructuring (middle layer)
`use`s state_types. Holds: `sort_cohorts`, `new_fuse_cohorts`, `fuse_2_cohorts` (+ private
`fuse_pass`), `split_cohorts`, `terminate_cohorts`, `max_cohort_count`, and
**`apply_recruitment`** (moved from dynamics; spawns cohorts from `recruit_pool` via `init_cohort`).

### 3.4 `meds_core_patch_fusefiss` — patch restructuring (upper layer)
`use`s state_types AND `meds_core_cohort_fusefiss` (patch fusion cascades into `sort_cohorts`
— an allowed patch→cohort edge). Holds: `sort_patches`, `new_fuse_patches`, `fuse_2_patches`
(+ private `patch_fuse_pass`), `terminate_patches`, `patch_light_profile`, private
`patch_compact`, and **`apply_patch_disturbance`** (moved from dynamics; creates the age-0 gap).

### 3.5 `meds_core_interface` — public façade (top of core)
Re-exports the engine surface so one `use meds_core_interface` gives callers everything:
`site_t`, `carbon_flux_block`, `cohort_deriv_block`; `update_cohort_states`,
`update_patch_states`, `update_overtopping_lai`; `init_cohort`; `apply_recruitment`,
`apply_patch_disturbance`; `sort_cohorts`, `new_fuse_cohorts`, `split_cohorts`,
`terminate_cohorts`; `sort_patches`, `new_fuse_patches`, `terminate_patches`. (The header
note already records that `update_demography` was dissolved into the driver — keep it.)

---

## 4. P1 — tendency-based `update_cohort_states` (the one behavioural change)

**Contract.** `core` no longer computes any growth/mortality arithmetic. It exposes a pure
applier over a **transient** tendency SoA:

```
type :: cohort_deriv_block           ! NEW, in meds_core_state_types
   real(wp), allocatable :: d_dbh_dt(:), d_height_dt(:), d_basal_area_dt(:), &
                            d_agb_dt(:), d_leaf_area_dt(:),                   &
                            d_leaf_carbon_dt(:), d_fineroot_carbon_dt(:),     &
                            d_wood_carbon_dt(:), d_nonstructural_carbon_dt(:), &
                            dln_nplant_dt(:)      ! summed mortality hazard, log-space
end type
```

**Why a bundle (not new persistent fields).** The tendencies are produced and consumed within
the SAME slow step, BEFORE any sort/fuse/split. So `cohort_deriv_block` is allocated per step,
index-aligned with `cohort_block`, and discarded — it is **not** part of `site_t`, **not**
reordered, and therefore **does NOT touch the six lockstep sites**. This is what makes P1
cheap; storing `d_*_dt` as persistent `cohort_block` fields (taxing all six sites + two birth
sites) is the design to avoid.

**Why store ALL field tendencies (not just the primary).** Each derived field's tendency is
`(allom(primary_new) − old)/dt`. Applying every one reconstructs the on-allometry value
exactly, so `update_cohort_states` needs NO allometry and NO empirical/carbon branch — "which
variable is primary and which direction the allometry runs" stays entirely in the driver, the
detail outside `core`.

**`update_cohort_states(cohort, deriv, dt)`** (pure, offloadable SIMD loop):
```
dbh += d_dbh_dt*dt ; height += d_height_dt*dt ; basal_area += d_basal_area_dt*dt
agb += d_agb_dt*dt ; leaf_area += d_leaf_area_dt*dt
{leaf,fineroot,wood,nonstructural}_carbon += d_*_dt*dt   ! then max(.,0) physical floor
nplant *= exp(min(max(dln_nplant_dt*dt, lnexp_min), lnexp_max))  ! then negligible-floor cull mark
```

**Where the arithmetic goes (driver / capi).** The old `growth_step`/`apply_growth` bodies
become the derivative computer:
- **Carbon path** — `meds_vegetation_dynamics` (`meds_aux`): `compute_slow_derivatives`
  runs `carbon_growth`→pools, the `wood_to_dbh` flip (`set_cohort_size_from_carbon`) to BACK
  OUT `d_dbh_dt`/`d_height_dt`/…, computes the Camac mortality → `dln_nplant_dt`, and updates
  the persistent `growth_avg` ring buffer. Then calls `update_cohort_states`.
- **Empirical path** — the capi `meds_apply_rates`: from the Python-supplied per-cohort growth
  `[cm/yr]`, forward-derive the full tendency set via allometry (the old `growth_step` math),
  set `dln_nplant_dt` from the supplied mortality, update the ring buffer, then call
  `update_cohort_states`. `growth_step`/`mortality_step` are removed.

**Consequences to accept:**
- **Device kernel relocates** from `core` to the driver (`compute_slow_derivatives`).
  `meds_aux` inherits the `-mp` flags (PUBLIC on `meds_core`), so it still offloads;
  `update_cohort_states` (pure `+=`) is itself offloadable. GPU path preserved.
- **Ring buffer** (`growth_avg` SMA read by Camac) stays PERSISTENT in `cohort_block`, updated
  in the derivative computer using the realized increment. Ordering unchanged: mortality reads
  `growth_avg` BEFORE this step's growth updates it — keep the driver's current sequence.
- **`dbh_critical` cap** and **pool nonneg** are rate-computation concerns → in the derivative
  computer (cap the increment) with a defensive floor in the applier.
- **Golden shift risk:** `old + (new−old)/dt·dt ≠ new` to the ULP. The empirical C-API golden
  reproduction may move by rounding. **Verify the golden's tolerance before committing P1**; if
  it is exact-match, either re-baseline or keep the empirical applier storing target-state
  deltas that reconstruct exactly.

---

## 5. P2 — `update_patch_states` placeholder
- Today: applies `age += dt` (moved from the driver's inline loop) — pure, no new deps.
- Future (P3): when `soil_carbon_t` becomes a per-patch field on `patch_index` (threaded
  through `patch_alloc`/`patch_compact`/`sort_patches` + `blend_*` on fusion/disturbance), its
  CENTURY tendency `d_soil_carbon_dt` is computed **driver-side** (calling `meds_soil_biogeochem`)
  and `update_patch_states` APPLIES it — so `core` gains NO biogeochem edge. No fast-loop merge.

## 6. P3 — `init_cohort`
Factor the ~19-field birth block (currently duplicated at
`meds_demography_dynamics%apply_recruitment` and `meds_init%add_cohort`) into
`meds_core_state_types%init_cohort(cohort, m, pft_table, ipft, owner_patch, nplant, dbh)`:
sets `pft/owner_patch/nplant/dbh`, `growth_avg=GROWTH_AVG_UNSET`, growth/pheno init,
`phenology_status=PHENOLOGY_STATUS_INIT`, the 7 `p_*` params, `leaf_temp=LEAF_TEMP_INIT`,
`psi=PSI_INIT`, **`overtopping_lai=0`** (resolves the add_cohort-omits-it divergence), then
forward `set_cohort_size`. Leaves capacity-ensure / `assign_cohort_id` / `rebuild_csr` / sort
at the two call sites. Does NOT absorb the clone-based creators (`split_cohorts`,
`apply_patch_disturbance` keep `copy_cohort_slot`).

## 7. P4 — diagnostics move
- `meds_demography_diagnostics.f90` → `src/io/meds_output_diagnostics.f90` (module rename),
  compiled into **`meds_io_prep`**. Verified clean: no consumer lives in `meds_core` or the
  `meds_aux` driver; `meds_output_integrate` (already in the prep target) is a consumer, so the
  move co-locates it. The module only `use`s `meds_kinds`/`meds_constants`/state_types/
  `meds_column_state_types` — all ≤ `meds_core`, satisfiable by `meds_io_prep`.
- `compute_overtopping_lai` does NOT move up — it becomes `update_overtopping_lai` in
  `meds_core_state_update` (structural, `intent(inout)`, feeds growth).

---

## 8. CMake changes (`CMakeLists.txt`)
1. `file(GLOB DEMOG_SOURCES … src/demography/*.f90)` → `src/core/*.f90`; `add_library(meds_demography …)` → `add_library(meds_core …)`; `target_link_libraries(meds_demography PUBLIC meds_shared)` → `meds_core`.
2. `target_compile_options(meds_demography PUBLIC ${MEDS_MP_FLAGS})` → `meds_core`.
3. `meds_output_core` target → rename **`meds_io_prep`**; add `src/io/meds_output_diagnostics.f90`; `target_link_libraries(meds_io_prep PUBLIC meds_core)`.
4. `meds_io` target → rename **`meds_io_stream`**; `target_link_libraries(meds_io_stream PUBLIC meds_core meds_netcdf_c meds_io_prep)`; `set(MEDS_IO_LIB meds_io_stream)`.
5. Repoint every `meds_demography` link → `meds_core` (in: `meds_io_prep`, `meds_io_stream`, `meds_aux`, `meds_testsupport`, the `foreach` test block) and every `meds_output_core` link → `meds_io_prep`.
6. `meds_aux`, `meds_main`, `meds_c` reach `MEDS_IO_LIB` via the variable → no edit beyond step 4.

## 9. `use`-site update checklist (module renames)
| Rename | Sites | Notes |
|--------|-------|-------|
| `meds_ecosystem_state` → `meds_core_state_types` | **20** (capi, io×2, init, driver×3, 8 tests, + the core files) | mechanical |
| `meds_demography_interface` → `meds_core_interface` | **9** (io, capi, phenology_driver, main, stepper, vegetation_dynamics, 2 tests) | mechanical |
| `meds_demography_diagnostics` → `meds_output_diagnostics` | **8** (main, capi, meds_io_stream, output_integrate, 4 tests) | domain change |
| `meds_demography_fusefiss` → `meds_core_{cohort,patch}_fusefiss` | io, init, 4 tests | repoint each `only:` symbol to the module (cohort vs patch) that now owns it |
| `meds_demography_dynamics` → gone | io, init | repoint to `state_update` / `cohort_fusefiss` / `patch_fusefiss` |
| `meds_competition` → gone | (only interface today) | consumers get `update_overtopping_lai` via `meds_core_interface` |

## 10. Ordered implementation (each phase compiles + `ctest` green on ifx AND nvfortran)
- **Phase A — rename + reshape (behaviour-preserving, bit-identical).** Rename dir/library
  `demography`→`core`; `meds_ecosystem_state`→`meds_core_state_types`; split fusefiss into
  cohort/patch; move `apply_recruitment`→cohort_fusefiss, `apply_patch_disturbance`→patch_fusefiss;
  fold `compute_overtopping_lai`→`meds_core_state_update%update_overtopping_lai`; **temporarily**
  keep `growth_step`/`mortality_step`/`apply_growth` in `meds_core_state_update` (still
  compute+apply); move diagnostics→`meds_output_diagnostics` in `meds_io_prep`; rename
  `meds_io`→`meds_io_stream`; rename interface; update all CMake + `use` sites. Validate
  bit-identical vs a pre-reorg output.
- **Phase B — P3 `init_cohort`.** Factor the two birth blocks; near bit-identical (only the
  `overtopping_lai`-at-birth zeroing is new — safe, recomputed each step).
- **Phase C — P1 tendency redesign.** Add `cohort_deriv_block`; create `update_cohort_states`;
  move the growth/mortality arithmetic + ring buffer to `compute_slow_derivatives` (driver) and
  the empirical tendency build (capi); delete `growth_step`/`mortality_step`/`apply_growth`.
  **Check the empirical golden tolerance** (ULP shift). Validate carbon spinup + Python examples.
- **Phase D — P2 placeholder.** Add `update_patch_states` (patch age); wire the driver to call it.
- **Phase E — docs.** Rewrite the CLAUDE.md "Demographic core" / source-layout section to
  `core`; fix the stale `meds_demography_rates`/`apply_carbon_npp` notes.

## 11. Risks & validation
- **Module cycle** if a primitive escapes state_types → keep §3's rule; build `--target meds_core` standalone.
- **Offload** — build the nvfortran multicore + GPU back ends every phase (a green ifx run hides the array-temp/offload traps).
- **Golden** — Phase C is the only numeric change; gate on the empirical C-API golden + carbon `tc_split` anchors.
- **Conservation** — fuse/split AGB + plant-number asserts must stay green through Phases A–D.

---

## 12. Final file organization of the `core` module + exposed functions

**Library `libmeds_core`** — `src/core/*.f90`, links `meds_shared` only, carries the OpenMP
`-mp` offload flags. Five files:

### `meds_core_state_types.f90` — module `meds_core_state_types`
> Types + the single lockstep/CSR primitive layer + cohort birth.

Exposed:
- **Types:** `cohort_block`, `patch_index`, `site_t`, `carbon_flux_block`, `cohort_deriv_block`
- **Constants:** `GROWTH_AVG_UNSET`, `PHENOLOGY_STATUS_INIT`
- **Primitives:** `site_alloc`, `site_free`, `cohort_ensure_capacity`, `cohort_reorder`,
  `cohort_compact`, `copy_cohort_slot`, `rebuild_csr`, `patch_ensure_capacity`,
  `gather_pft_params`, `set_cohort_size`, `set_cohort_size_from_carbon`, `assign_cohort_id`,
  `assign_patch_id`
- **Birth:** `init_cohort`

### `meds_core_state_update.f90` — module `meds_core_state_update`
> Pure apply/refresh of (derived) state. `use`s state_types.

Exposed: `update_cohort_states`, `update_patch_states`, `update_overtopping_lai`

### `meds_core_cohort_fusefiss.f90` — module `meds_core_cohort_fusefiss`
> Cohort restructuring + recruitment. `use`s state_types.

Exposed: `sort_cohorts`, `new_fuse_cohorts`, `fuse_2_cohorts`, `split_cohorts`,
`terminate_cohorts`, `max_cohort_count`, `apply_recruitment`
(private: `fuse_pass`)

### `meds_core_patch_fusefiss.f90` — module `meds_core_patch_fusefiss`
> Patch restructuring + disturbance. `use`s state_types + cohort_fusefiss.

Exposed: `sort_patches`, `new_fuse_patches`, `fuse_2_patches`, `terminate_patches`,
`patch_light_profile`, `apply_patch_disturbance`
(private: `patch_fuse_pass`, `patch_compact`)

### `meds_core_interface.f90` — module `meds_core_interface`
> The one-`use` public façade of the engine (re-exports only).

Re-exports: `site_t`, `carbon_flux_block`, `cohort_deriv_block`; `init_cohort`;
`update_cohort_states`, `update_patch_states`, `update_overtopping_lai`; `apply_recruitment`,
`apply_patch_disturbance`; `sort_cohorts`, `new_fuse_cohorts`, `split_cohorts`,
`terminate_cohorts`; `sort_patches`, `new_fuse_patches`, `terminate_patches`

**Moved OUT of core:**
- `meds_output_diagnostics` (`src/io/`, target `meds_io_prep`) — the pure reductions
  (`total_agb`, `total_lai`, `total_nplant`, `total_basal_area`, `total_area`, `mean_dbh`,
  `total_gpp`, `total_npp`, `mean_can_temp`, `mean_soil_temp_top`, `total_et`,
  `site_soil_temp_column`, `site_soil_water_column`, `count_cohorts`, `has_nan`, `print_summary`).
- The growth/mortality **arithmetic** (old `growth_step`/`apply_growth`/`mortality_step`) →
  `meds_vegetation_dynamics%compute_slow_derivatives` (driver) + the capi empirical tendency build.
