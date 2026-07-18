# MEDS Codebase Reorganization — a Carbon-Driven Standalone Fortran Model, a Python Companion, and a Driver-Orchestrated `demography ⊥ plant` Engine

> **STATUS: DESIGN-ONLY** (2026-07-16, v6, not implemented, no branch). A cross-cutting reorganization that
> retires two zero-isolation library boundaries, moves the geometric allometry (incl. the carbon→structure
> composite) into `shared`, merges the cohort/patch state ontology into the `demography` library **as pure
> apply-primitives**, relocates the per-individual vital-rate kernels into `plant`, and **extracts the
> orchestration (`update_demography`) into the driver (`meds_vegetation_dynamics`)** — which is the single
> place plant and demography meet. This **keeps `demography ⊥ plant`**: the demography module only *applies*
> supplied changes (using shared allometry), while the driver *computes* the rates (calling the plant
> kernels). The coupled fast+slow loop is exposed through **one opt-in `libmeds_c.so`** so a **`meds` Python
> package** can drive, inspect, and override it — **without MEDS depending on Python to run** (the standalone
> Fortran model is the **carbon** path).
>
> Locked decisions: **(#1)** allometry→`shared/functions/` (+ fold in `carbon_to_structure`); **(#2)** merge
> the SoA into `libmeds_demography` (keep the name `demography`) + rename `meds_demography_types →
> meds_ecosystem_state`; **(#3)** split `carbon_vital_rates` — **rate kernels → `plant/meds_plant_vital_rates.f90`**,
> the **rate computation + orchestration → the driver**, and the **apply-primitives stay in `demography`**
> (`demography ⊥ plant` PRESERVED); empirical growth/recruitment → **Python** (`example_demography`);
> **(#4)** `libmeds_c` opt-in + CPU-only first; **(#5)** `apply_carbon_npp → apply_growth`, a demography
> apply-primitive. Retained: **(A)** `overtopping_lai` a stored per-cohort field; **(B)** fast+slow always
> coupled; `meds_column_state_types` stays in `shared/state/`; optional/orthogonal `meds_column_co2 →
> biophysics`. Empirical growth removed in **P1** (standalone Fortran demography is carbon-only from P1 —
> sequencing risk, §10). No implementation until approved.

Grounded on verified facts:

1. **`meds_allometry` uses only `meds_kinds` + `meds_constants`** → into `shared` is cycle-free; it already
   hosts the carbon↔structure helpers (`wood_to_dbh`, `size2wood_carbon`, `size2leaf_carbon`), so
   `carbon_to_structure` folds in cleanly.
2. **No target links `meds_state` without `meds_demography`** → merging the SoA into the library is free.
3. **The stateless-kernel wall** (`plant`/`biophysics`/`biogeochem` link `shared` only, no `site_t`) is
   preserved, AND **`demography ⊥ plant` is preserved** (§2.3): the demography apply-primitives never call a
   plant kernel; only the driver does.
4. **`meds_column_co2` does not use `meds_soil_biogeochem`** → §2's optional fast/slow split is clean.
5. **`growth_dbh_*`/`growth_lai_slope` appear only** in `meds_demography_rates` (empirical), `meds_config_io`,
   `meds_pft_params` — the carbon path never uses them, so they are cleanly removable.

---

## 1. Governing decisions [locked]

| # | Decision | Realization |
|---|----------|-------------|
| 1 | Allometry is a geometric CONSTRAINT relationship. | Move `meds_allometry.f90` → `src/shared/functions/`; **fold `carbon_to_structure` into it**; delete the `meds_allometry` library target. |
| 2 | The cohort/patch SoA is the state ONTOLOGY; keep the library named `demography`. | Merge `src/state/` into `libmeds_demography`; rename `meds_demography_types → meds_ecosystem_state`. `column_state_types` stays in `shared/state/`. |
| 3 | Vital rates are per-plant physiology; the engine only APPLIES; the driver ORCHESTRATES. | Rate kernels → `plant/meds_plant_vital_rates.f90` (pure). The rate computation + the old `update_demography` sequencing → **`meds_vegetation_dynamics` (driver)**. The **apply-primitives stay in `demography`** (`apply_growth`, `mortality_step`, `apply_recruitment`, `apply_patch_disturbance`, fuse/fiss) — all ⊥ plant. Empirical growth/recruitment + params → **Python**. Standalone Fortran = the **carbon** path. |
| 4 | MEDS ships as one Python package that can call functions in each module. | One opt-in `libmeds_c.so` + one `meds` package; CPU-only first. Python is a **companion**, not required. |
| 5 | Carbon growth application is a demography apply-primitive. | `apply_carbon_npp → apply_growth`: pools += npp → `carbon_to_structure` (shared) → geometry + the Camac ring buffer. No plant call. |
| A | `overtopping_lai` is a useful diagnostic. | Stored per-cohort field on `meds_ecosystem_state`, filled by the law-free `meds_competition` sweep; lockstep + output registry + C-API getter. |
| B | Fast + slow always run together. | The standalone Fortran executable runs both (carbon path); the C-API exposes both loops so Python can drive/inspect/override. |

---

## 2. Target architecture

### 2.1 Folders / compile units

```
src/shared/        libmeds_shared.a  (root; ONE library, cosmetic subfolders)
   base/       meds_kinds, meds_constants
   config/     meds_config, meds_forcing_config, meds_output_config, meds_pft_params   # (loses growth_dbh_*/growth_lai_slope)
   state/      meds_column_state_types                # boundary type; STAYS in shared
   util/       meds_time, meds_numerics, meds_budget_check
   functions/  meds_temp_response, meds_thermo, meds_allometry (+ carbon_to_structure)  # pure relationship kernels

src/plant/         libmeds_plant.a  (stateless kernels, shared only; standalone; backs meds.leaf)
   ... existing leaf/hydraulics/phenology/respiration/carbon kernels ...
   meds_plant_vital_rates    (NEW: carbon_growth_rate, camac_mortality, recruitment_contribution;
                              pure scalar-in/scalar-out, NO site_t; CALLED BY THE DRIVER, not the engine)

src/demography/    libmeds_demography.a  (MERGED src/state/ + apply-primitives; links shared only; ⊥ PLANT)
   meds_ecosystem_state      (RENAMED from meds_demography_types; SoA + lockstep + CSR; composes shared column types)
   meds_demography_dynamics  (apply-PRIMITIVES: growth_step[supplied rate, offloaded], apply_growth[carbon npp],
                              mortality_step[supplied hazard], apply_recruitment[supplied density],
                              apply_patch_disturbance)  # NONE call plant
   meds_demography_fusefiss  (sort/fuse/fission; carbon re-derivation → carbon_to_structure)
   meds_demography_diagnostics
   meds_competition          (compute_overtopping_lai sweep; law-free)
   -- REMOVED: meds_demography_interface/update_demography (orchestration → driver); meds_demography_rates

src/biophysics/    libmeds_biophysics.a   (stateless, shared only; -- OPTIONAL: gains meds_column_co2)
src/biogeochemistry/ libmeds_biogeochemistry.a  (-- OPTIONAL: reduces to slow soil carbon)
src/forcing/       libmeds_forcing.a
src/io/            libmeds_io.a / output_core / config_io / netcdf_c   # plant-FREE again (⊥ preserved)

src/driver/        libmeds_aux.a  (the RUNNABLE model + the ORCHESTRATOR; links demography + plant + others)
   meds_vegetation_dynamics  (THE ORCHESTRATOR — absorbs the old update_demography sequencing: assembles NPP
                              [get_plant_flux_slow], calls the plant vital-rate kernels for mortality/recruitment,
                              interleaves the demography apply-primitives, handles cadence + fuse/fission.
                              THE single plant↔demography meeting point.)
   meds_fast_loop, meds_column_dynamics, meds_stepper   (the FAST loop + cadence)
   meds_main                 (the STANDALONE executable — runs fast+slow, CARBON path)

src/capi/          libmeds_c.so  (NEW, opt-in -DMEDS_BUILD_PYLIB)
   meds_fastloop_capi (advance_fast), meds_slowloop_capi (advance_slow = the driver's carbon orchestration),
   meds_demography_capi (opaque site handle + getters + the apply-primitives for custom orchestration),
   meds_plant_capi (leaf/carbon kernels + the vital-rate kernels), meds_config_io_capi, meds_forcing_capi, meds_io_capi

python/meds/       companion package (NOT required); examples/example_demography/ = the EMPIRICAL example (numpy)
```

### 2.2 Target library DAG (acyclic)

```
shared ─┬─ demography ─┬─ output_core ─┐
        │              ├─ io ──────────┤   aux → demography, plant, biophysics, biogeochem, forcing, config_io, io
        │              └─ aux ─────────┤   capi → aux, demography, plant, io, forcing (opt-in)
        ├─ plant ───────── (aux, capi) ┤   demography ⊥ plant PRESERVED (siblings above shared)
        ├─ biophysics ─────────────────┤
        ├─ biogeochemistry ────────────┼─→ libmeds_c.so (opt-in)
        ├─ config_io ──────────────────┤
        └─ netcdf_c ── forcing ────────┘   meds_main → aux, io   (STANDALONE, no Python, carbon path)
```

### 2.3 `demography ⊥ plant` — preserved by the mechanism/policy split [#3]

The plant dependency lives in the *rate computation*, not the *apply*. So:

- **The demography apply-primitives never call plant.** `apply_growth` (carbon npp → geometry) uses
  `carbon_to_structure` in **shared** `meds_allometry` — not a plant kernel. `mortality_step`/`apply_recruitment`
  apply **supplied** arrays. Fuse/fission use shared allometry.
- **The driver (`meds_vegetation_dynamics`) computes the rates** — calling `get_plant_flux_slow` for NPP and the
  `meds_plant_vital_rates` kernels for the Camac hazard + recruitment — and **interleaves** them with the
  demography apply-primitives. It is the single plant↔demography meeting point.

The carbon slow step (driver) sequences cleanly (mortality's `growth_avg` dependency resolves because
`apply_growth` folds the current sample into `growth_avg`, so `camac_mortality(growth_avg)` needs no separate
`dbh_rate`):

```
1. per cohort: get_plant_flux_slow → npp            [plant]
2. demography.apply_growth(npp)                       [demography ⊥ plant]   # pools+=npp → carbon_to_structure → geometry + ring buffer
3. read growth_avg (updated)                          [getter]
4. per cohort: camac_mortality(growth_avg, traits)   [plant kernel]
5. demography.mortality_step(mortality)               [demography ⊥ plant]
6. per cohort: recruitment_contribution(...)         [plant kernel]
7. demography.apply_recruitment(recruitment)          [demography ⊥ plant]
8. patch ageing, fuse/fission per cadence            [demography]
```

Two benign consequences: the **order-of-operations invariant moves** from `update_demography` into
`meds_vegetation_dynamics` (appropriate — it is the orchestrator, ED2's `veg_dynamics_driver` pattern); and a
**one-sample `growth_avg` shift** in the carbon path (mortality reads the current-inclusive `growth_avg`),
acceptable since the carbon path is new (define the golden with this order).

**Optional orthogonal split (independent PR):** `meds_column_co2 → biophysics`; `biogeochemistry` reduces to
slow soil carbon; split `meds_biogeochem_types`.

---

## 3. The `demography` state layer — merge, rename, `overtopping_lai`

**Merge + rename [#2].** `src/state/meds_demography_types.f90` + the apply-primitives compile into
`libmeds_demography`; rename `meds_demography_types → meds_ecosystem_state` (~20 `use` sites).
`column_state_types` is the boundary type (biophysics + demography + io/init/driver) and stays in `shared/state/`.

**`overtopping_lai` stored field [A].** `real(wp), allocatable :: overtopping_lai(:)` on `cohort_block`; filled
by `meds_competition%compute_overtopping_lai(site)` — a **law-free geometric reduction** (stays in the engine,
NOT plant). Refreshed after every structural change + once per slow step; recomputed (not conserved) on fusion;
added to the lockstep (recommend a `meds_cohort_fields.inc` manifest), the output registry, and a C-API getter
(Python reads it as competition context for the empirical example).

**`carbon_to_structure` [#1/#2].** Folded into `meds_allometry` (shared/functions) beside `wood_to_dbh` /
`size2*carbon`; `set_cohort_size_from_carbon` deleted (its body IS `carbon_to_structure`).

---

## 4. `apply_growth` — the carbon growth apply-primitive [#5, comment 1]

`apply_carbon_npp → apply_growth`, a demography apply-primitive in `meds_demography_dynamics`. Per cohort: adds
`npp` to the pools, flips carbon → geometry once via `carbon_to_structure` (shared) → new
`dbh`/`height`/`agb`/`leaf_area`, derives the `dbh`-rate, and folds it into the Camac `growth_avg` ring buffer.
**No plant call, no mortality computation** (that is the driver's job, §2.3). The offloaded supplied-rate sibling
`growth_step` (used by the Python/empirical path) is unchanged. **Bonus:** this removes the duplicate
`wood_to_dbh(wood_carbon+npp_wood)` (today in both `carbon_vital_rates` and `apply_carbon_npp`) — one pass now.

---

## 5. Fortran conducts (carbon path); the empirical path is a Python example [#3, B]

**Standalone Fortran run (no Python) = the carbon path.** `meds_main` → `meds_stepper` → per step the fast loop
(`meds_fast_loop`) + `meds_vegetation_dynamics`, which orchestrates the carbon slow step (§2.3): NPP + plant
vital-rate kernels + demography apply-primitives + cadence.

**The empirical path is a Python example.** `examples/example_demography/` implements the phenomenological growth
+ recruitment laws + params (`growth_dbh_*`, `growth_lai_slope`) in numpy; it reads state (`dbh`, `height`,
`overtopping_lai`, `growth_avg` via getters), computes growth + recruitment (its own laws) + mortality (via the
Fortran `camac_mortality` kernel through the C-API), and **orchestrates the demography apply-primitives itself**
(`growth_step` for its supplied rate, `mortality_step`, `apply_recruitment`, fuse/fission). A demography-only
experiment. **No `meds/demography` package folder.**

**Coarse verbs:** `advance_fast(...)` (all fast sub-steps in Fortran); `advance_slow(...)` (the driver's carbon
orchestration); the demography apply-primitives + the plant vital-rate kernels (for Python-orchestrated
experiments); getters.

---

## 6. The C-API contract

**One `libmeds_c.so`, opt-in** (`-DMEDS_BUILD_PYLIB=ON`, default OFF), CPU-only first. PIC on every static lib.
Hand-written `iso_c_binding` shims (f2py can't handle the allocatable/derived-type design).

**`site_t` crosses as an opaque handle** — a module `save` registry; Python holds an `integer(c_int)` index.
Only flat arrays (copied) + scalars cross. Config + forcing are opaque handles too.

**The generation-counter guard.** The SoA reorders on every fuse/fission, so a cached Python array is positionally
invalid — `global_id` is the only stable key. Getters return data + `global_id[]` + a per-handle `generation`;
the apply-primitives error on a stale `generation`. Getters always copy.

**Exposed verbs:** lifecycle; `advance_fast`, `advance_slow`; read getters (incl. `overtopping_lai`, `growth_avg`,
`gpp_accum`, resp accumulators, `phenology_status`, `patch_csr`) with `global_id`+`generation`; the **demography
apply-primitives** (`growth_step`, `apply_growth`, `mortality_step`, `apply_recruitment`, `apply_patch_disturbance`,
fuse/fission) for Python orchestration; the **plant kernels** (leaf + carbon + `camac_mortality`/
`carbon_growth_rate`/`recruitment_contribution`); `io_*` output verbs.

---

## 7. Python package + the empirical example

`python/meds/` is a companion package over `libmeds_c.so` (a single `CDLL`; `meds.leaf` folds in). It exposes
`site.py` (handle + generation guard), the kernel wrappers (incl. the plant vital-rate kernels), the demography
apply-primitives, and an optional `run.py`. The empirical demography example lives in
`examples/example_demography/` (numpy growth/recruitment + the `camac_mortality` kernel via C-API, orchestrating
the apply-primitives). No `meds/demography` submodule.

---

## 8. Regression + reproducing the current demography example

- **P0 — capture the golden.** From the current code, run the empirical demography spinup; serialize the
  trajectory (keyed by `global_id`).
- **Empirical example → Python, to TOLERANCE.** `examples/example_demography/` reproduces the golden to ~1e-8:
  numpy growth/recruitment mirror the removed Fortran formulas; mortality is the same Fortran `camac_mortality`
  kernel via the C-API.
- **Fortran refactor → BIT-IDENTICAL on surviving paths.** The relinks/merge/rename + the mechanism/policy split
  are behavior-preserving for the **engine invariants** and (with the one-sample `growth_avg` note) the **carbon
  path** — fusefiss conservation, disturbance, sort/CSR, init, fast loop, a carbon spinup, `tc_split(54)=292.450065`.
  `test_rates` removed; `test_spinup` (empirical) → the Python example test. C-API round-trip test added.

---

## 9. Phased migration plan

- **P0 — Capture the empirical spinup golden.**
- **P1 — Relinks + merge/rename + subfolders + `overtopping_lai` + the mechanism/policy split.** `allometry →
  shared/functions/` (**fold in `carbon_to_structure`**); reorganize `shared` into
  `base/config/state/util/functions/`; merge `src/state/` into `libmeds_demography` (rename
  `meds_ecosystem_state`); add the stored `overtopping_lai` field + `meds_competition` sweep. **Split
  `meds_demography_rates`:** the rate kernels → `src/plant/meds_plant_vital_rates.f90`; **move the
  `update_demography` orchestration + the plant-kernel calls into `meds_vegetation_dynamics` (driver)**; keep the
  apply-primitives in `meds_demography_dynamics` (rename `apply_carbon_npp → apply_growth`); **delete the empirical
  growth/recruitment laws + `growth_dbh_*`/`growth_lai_slope`** (pft_params, config_io). `demography ⊥ plant`
  preserved (io/output_core stay plant-free). Standalone Fortran demography is **carbon-only** from here. Verify
  the carbon spinup + engine invariants.
- **P2 — Finish the carbon flip.** Route `apply_growth`'s geometry through `carbon_to_structure` in
  `meds_allometry`; delete `set_cohort_size_from_carbon`.
- **P3 — C-API + `libmeds_c`.** `src/capi/` with the opaque handle + generation guard, `advance_fast`,
  `advance_slow`, the apply-primitives, the plant vital-rate kernels, getters, config/forcing/io shims; behind
  `-DMEDS_BUILD_PYLIB`. Land read-only getters + a no-op step first.
- **P4 — Python companion + `example_demography`.** Unify `python/meds/`; implement the empirical
  growth/recruitment in `examples/example_demography/`; orchestrate the apply-primitives via the C-API and
  reproduce the P0 golden to tolerance; add the C-API round-trip test.

*(Optional parallel track: `meds_column_co2 → biophysics` + `meds_biogeochem_types` split.)*

---

## 10. Pros / cons / risks / deferred

**Pros.** Acyclic domain DAG; **`demography ⊥ plant` preserved** (io/output_core plant-free; demographic core
builds standalone); mechanism (demography apply-primitives) cleanly separated from policy (driver orchestration,
ED2-like); vital-rate kernels grouped with the per-plant physiology; the carbon rate calc dedups the
`dbh`-from-carbon computation; `plant` stays stateless/standalone (backs `meds.leaf`); the empirical growth curve
demoted to an editable Python example; the hot fast loop + fuse/fission + GPU offload stay in Fortran;
`overtopping_lai` a first-class diagnostic; `carbon_to_structure` consolidated in `meds_allometry`.

**Cons.** The stateful-engine C-API is a new hand-written marshalling subsystem (opt-in; CI needs a pylib job) —
now with more surface, since Python orchestrates via the granular apply-primitives rather than one
`update_demography` verb. The empirical path exists only in Python.

**Biggest risk (elevated by "remove empirical now").** From P1 the standalone Fortran demographic path is
carbon-only, so it depends on a functional carbon spinup (fast loop → GPP → NPP → growth). RT is joined and
forcing is wired, but a full carbon spinup is not yet a validated regression. **De-risk:** validate a Fortran
carbon-path spinup early in P1, keeping the P0 golden + the Python empirical reproduction as the behaviour anchor.
Secondary: the opaque-handle + resizing-SoA + Python-snapshot trap — generation counter + `global_id` keying +
copy-only getters.

**Deferred / open.** Engine module names kept `meds_demography_*`; `growth_step` naming (parallel to
`apply_growth`?); `meds_pft_params` in `shared/config/`; re-entrancy (process-global allometry coeffs/config);
Python-owned config; GPU pylib; the optional `column_co2 → biophysics` split.
