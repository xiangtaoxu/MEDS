# MEDS C-API Reorganization — Design

**Status:** design-only (no code changed). 2026-07-20.
**Scope:** restructure the ISO_C_BINDING layer (`src/capi/`, `src/plant/meds_plant_capi.f90`)
and the `meds` Python package + `example_demography` scripts. No numerics change; the
existing goldens are the regression oracle throughout.

---

## 1. Principles (the contract this design must satisfy)

1. **capi invents nothing.** A `*_capi.f90` file contains only thin `bind(c)` shims that
   forward to *pre-existing public Fortran API*. No orchestration, no allometry, no
   derived logic authored inside a shim. (This retires the current
   `compute_empirical_derivatives`, which exists only in the shim.)
2. **The container lives in Fortran** behind the opaque `site_t` handle (decided). Python
   never owns the cohort/patch SoA; it reads copies out and drives the engine through the
   handle. This keeps the marshalling boundary narrow and leaves fuse/fission — a deep,
   whole-record conservation merge — entirely inside Fortran.
3. **One consolidated shared library** (`libmeds.so`) backs every Python submodule, loaded
   once. No more per-domain `.so` files.
4. **capi files mirror the source folder/module they wrap.** One shim per domain, all in
   `src/capi/`, so folder ↔ library is 1:1 and the folder glob is correct.
5. **Orchestration lives in Python; laws live in Python.** The *empirical* vital-rate laws
   and their apply-sequence are a Python **example**. The *carbon* path is the Fortran
   **model**, exposed as a single call.

### 1a. The organizing asymmetry (carbon vs. empirical)

This is the spine of the design and resolves the apparent tension between Principle 1
("no orchestration in capi") and the fact that the fuse/fission *ordering* is engine
mechanism:

| Path | What it is | capi surface | Who orchestrates |
|---|---|---|---|
| **Carbon** (`advance_slow`) | the Fortran **model** | **one** thin binding → `vegetation_dynamics` | Fortran (it *is* the model) |
| **Empirical** (rates) | a Python **example** | thin bindings to the **law-free apply-primitives** | **Python** (laws + sequence) |

The carbon cadence stays one Fortran routine because it is the model. The empirical
cadence moves to Python because it is an example — and Python drives *published core
primitives*, so capi still invents nothing. The two paths never share a hand-rolled order,
so there is no drift risk between them.

---

## 2. Part A — the one enabling Fortran change

Today `meds_apply_rates` ([src/capi/meds_demography_capi.f90:104](../../src/capi/meds_demography_capi.f90))
does three things inside the shim: (i) the private `compute_empirical_derivatives`
forward-allometry→tendency map, (ii) the pure applier `update_cohort_states`, (iii) the
hand-rolled fuse/fission cadence. (i) and the cadence sequencing violate Principle 1.

**Change:** add exactly **one** public core apply-primitive and delete the shim-invented code.

- **New public primitive** `grow_cohorts_by_rate(site, cfg, growth, mortality)` in
  `meds_core_state_update` (exported via `meds_core_interface`). It absorbs the current
  `compute_empirical_derivatives` + the `growth_hist_pos`/`cohort_deriv_alloc` bookkeeping +
  `update_cohort_states`. It is **law-free** (rates are inputs) — a sibling of the existing
  `apply_recruitment`, which is also a rate-driven apply-primitive. Its allometry comes from
  the shared `meds_allometry` module (no re-inlined constants).
- **Delete** `compute_empirical_derivatives` and the hand-rolled sequence from the shim.
- **Everything else the empirical path needs is already public** in `meds_core_interface`:
  `apply_recruitment`, `update_patch_states`, `new_fuse_cohorts`, `terminate_cohorts`,
  `split_cohorts`, `sort_cohorts`, `apply_patch_disturbance`, `new_fuse_patches`,
  `terminate_patches`, `sort_patches`, `update_overtopping_lai`.

This is behavior-preserving: `grow_cohorts_by_rate` is a verbatim relocation of the shim's
current geometry math, and the Python-side sequence (Part D) reproduces the exact order at
[src/capi/meds_demography_capi.f90:130-154](../../src/capi/meds_demography_capi.f90). The
`empirical_spinup` golden must stay bit-identical.

> **Optional dedup (secondary):** the carbon driver `vegetation_dynamics` runs its own
> fuse/fission cadence. A shared `apply_structural_dynamics(site, cfg, flags)` could serve
> both, but the empirical path sorts *only* on the monthly/annual cadence while carbon
> re-sorts every step — so the sort placement differs. Treat the shared-cadence extraction
> as a later cleanup, not part of this reorg.

---

## 3. Part B — capi file decomposition

All shims move into `src/capi/`, one module per source domain. `meds_plant_capi.f90`
relocates from `src/plant/` unchanged. The monolithic `meds_demography_capi.f90` unpacks:

| File (`src/capi/`) | Module | Wraps (source domain) | `bind(c)` symbols |
|---|---|---|---|
| `meds_capi_registry.f90` | `meds_capi_registry` | **infra** — no `bind(c)` | owns `MAXH`, `g_site`, `g_cfg`, `*_used`, `g_generation`, `f_string`; helpers `claim_cfg_slot`, `claim_site_slot`, `release_site`, `bump_generation` |
| `meds_config_capi.f90` | `meds_config_capi` | `shared/config` | `meds_config_load`, `meds_config_dt_years`, `meds_config_n_pft` |
| `meds_site_capi.f90` | `meds_site_capi` | `core` (`meds_core_state_types`) | `meds_site_create`, `meds_site_free`, `meds_site_generation`, `meds_site_n_cohort`, `meds_site_n_patch`, `meds_site_get_real`, `meds_site_get_int` |
| `meds_init_capi.f90` | `meds_init_capi` | `init` | `meds_site_init_bare` |
| `meds_diagnostics_capi.f90` | `meds_diagnostics_capi` | `io` (`meds_output_diagnostics`) | `meds_site_total_agb`, `_lai`, `_nplant`, `_basal_area` |
| `meds_dynamics_capi.f90` | `meds_dynamics_capi` | `driver` | `meds_advance_slow` *(carbon path)* |
| `meds_demography_capi.f90` | `meds_demography_capi` | `core` (`meds_core_interface`) | the **law-free apply-primitives**: `meds_grow_cohorts_by_rate`, `meds_apply_recruitment`, `meds_age_patches`, `meds_new_fuse_cohorts`, `meds_terminate_cohorts`, `meds_split_cohorts`, `meds_sort_cohorts`, `meds_apply_patch_disturbance`, `meds_new_fuse_patches`, `meds_terminate_patches`, `meds_sort_patches`, `meds_update_overtopping_lai` |
| `meds_allometry_capi.f90` | `meds_allometry_capi` | `shared/functions` (`meds_allometry`) | stateless: `meds_dbh_to_height`, `meds_dbh_to_agb`, `meds_height_to_dbh`, `meds_min_cohort_carbon`, … (as needed by the empirical example) |
| `meds_plant_capi.f90` | `meds_plant_capi` | `plant` *(moved, unchanged)* | leaf gas exchange + phenology (existing) |

**Naming:** domain shims keep `meds_<domain>_capi`; the infra module flips to
`meds_capi_registry` to signal "plumbing, not a domain shim." The name
`meds_demography_capi` now genuinely means *the demographic apply-primitives* — the
old monolith's `advance_slow`, config, site, init, and diagnostics symbols have moved to
their own files.

**Dependency shape** — a clean fan-out; no shim depends on another shim:

```
meds_capi_registry ◄── config / site / init / diagnostics / dynamics / demography
meds_allometry_capi, meds_plant_capi ── independent (link shared/plant only)
```

The Fortran module-dependency scanner orders `registry.mod` automatically.

---

## 4. Part C — build (one library)

Move `meds_plant_capi.f90` into `src/capi/`, delete the separate `meds_plant_c` target and
the `PLANT_CAPI` glob/`REMOVE_ITEM` in the plant block, and collapse to one shared target
(folder ↔ library is now 1:1, so the folder glob is correct):

```cmake
if(MEDS_BUILD_PYLIB)
   file(GLOB_RECURSE CAPI_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/capi/*.f90)
   add_library(meds_c SHARED ${CAPI_SOURCES})
   target_link_libraries(meds_c PRIVATE meds_aux ${MEDS_IO_LIB})   # superset: covers plant/core/driver/io
   set_target_properties(meds_c PROPERTIES OUTPUT_NAME meds)       # -> libmeds.so
   meds_fortran_flags(meds_c)
endif()
```

- `-fPIC` for the static deps is already forced under `MEDS_BUILD_PYLIB`
  ([CMakeLists.txt:28-30](../../CMakeLists.txt)); no change.
- `meds_aux` already PUBLIC-links `meds_plant`, so the plant symbols resolve in the merged lib.
- `bind(c)` names stay globally unique (verified — no overlap across shims). **Standing rule:**
  every exported symbol is domain-prefixed (`meds_site_*`, `meds_config_*`, `meds_grow_*`, …)
  since the single `.so` flattens them into one symbol table.

**Cost accepted:** `meds.plant` no longer ships as a tiny stateless `.so`; importing it now
`dlopen`s a library that *contains* the demographic engine + netCDF. Nothing extra runs; the
file is just bigger. Correct trade for a single distributable `meds` wheel.

---

## 5. Part D — Python package changes

### 5.1 One loader

Replace the two `_ffi.py` loaders with a single package-level loader that `CDLL`s
`libmeds.so` once; submodules import the shared handle and register their own signatures.

```python
# python/meds/_ffi.py  (new, shared)
import ctypes, os, glob
from pathlib import Path

def _find_lib():
    if (p := os.environ.get("MEDS_LIB")): return p
    here = Path(__file__).resolve().parent
    cands = [here / "libmeds.so"] + [
        c for bd in ("build", "build-py", "build-pylib")
        for c in glob.glob(str(here.parents[1] / bd / "libmeds.so"))]
    for c in cands:
        if os.path.exists(c): return str(c)
    raise OSError("libmeds.so not found. Build -DMEDS_BUILD_PYLIB=ON or set MEDS_LIB.")

lib = ctypes.CDLL(_find_lib())

def bind(name, restype, argtypes):
    fn = getattr(lib, name); fn.restype = restype; fn.argtypes = argtypes
    return fn
```

`meds/plant/_ffi.py` and `meds/demography/_ffi.py` become thin: `from .._ffi import lib, bind`,
then register their symbols. No submodule loads its own `.so`.

### 5.2 `meds.demography` — primitive bindings + a Python cadence helper

The `Config`, `Site` handle lifecycle, getters, `advance_slow` (carbon), and the scalar/array
diagnostics are **unchanged** — the container stays in Fortran. The single change is that
`Site.apply_rates` (which bound the shim-orchestrated `meds_apply_rates`) is replaced by
`Site.apply_demography(...)`, a **Python** orchestration over the new primitive bindings:

```python
# python/meds/demography/_site.py  (replaces apply_rates)
def apply_demography(self, growth, mortality, recruitment, is_new_month, is_new_year):
    """Empirical apply-cadence, orchestrated in Python over the Fortran apply-primitives.

    Mirrors the demographic engine's structural order (formerly hand-rolled inside the
    capi shim). Laws (growth/mortality/recruitment) come from the caller; every step
    here is a published law-free core primitive driven through the handle.
    """
    s, c = self.handle, self.cfg.handle
    g = _as_c_double(growth); m = _as_c_double(mortality); r = _as_c_double_F(recruitment)
    cff = bool(is_new_month); pd = pf = bool(is_new_year)          # cadence gates

    lib.meds_grow_cohorts_by_rate(s, c, g, m)                      # was compute_empirical_derivatives + applier
    lib.meds_age_patches(s, c)
    if cff:
        lib.meds_apply_recruitment(s, c, r)
        lib.meds_new_fuse_cohorts(s, c); lib.meds_terminate_cohorts(s, c)
        lib.meds_split_cohorts(s, c);    lib.meds_sort_cohorts(s)
    if pd:
        lib.meds_apply_patch_disturbance(s, c)
    if pf:
        lib.meds_sort_patches(s);        lib.meds_new_fuse_patches(s, c)
        lib.meds_terminate_patches(s, c);lib.meds_new_fuse_cohorts(s, c)
        lib.meds_terminate_cohorts(s, c);lib.meds_sort_cohorts(s)
    lib.meds_update_overtopping_lai(s)
    self._bump_generation_local()   # or read lib.meds_site_generation(s)
```

This is the *visible architectural change*: the empirical cadence moved from the Fortran
shim into the Python package, driving published primitives. The carbon path
(`Site.advance_slow` → `vegetation_dynamics`) stays one Fortran call.

### 5.3 `meds.allometry` — new stateless submodule (removes duplication)

Expose the `meds_allometry_capi` kernels as a small submodule so the example stops
hardcoding constants (see 6.2).

---

## 6. How the `example_demography` scripts change

### 6.1 `run_carbon.py` — **only the run command + header**

The carbon path is untouched (`Site.advance_slow` still binds `vegetation_dynamics`). The
only edits:

- `MEDS_LIB=build-pylib/libmeds_c.so` → `MEDS_LIB=build-pylib/libmeds.so` (docstring + any
  wrapper script).
- Header comment: "via `libmeds_c`" → "via `libmeds`".

Driving logic, output, and the `meds_main` agreement are unchanged.

### 6.2 `empirical_laws.py` — **call real allometry, drop the hardcoded constants**

Currently duplicates the Fortran allometry
([empirical_laws.py:17-31](../../examples/example_demography/empirical_laws.py)):
`B1HT, B2HT = 1.139963, 0.564899`, `AGB_C1, AGB_C2 = …`, and local `dbh_to_height` /
`dbh_to_agb` / `height_to_dbh`. Replace these with `meds.allometry` kernel calls so there is
a single source of truth (kills silent drift if the Fortran allometry changes):

```python
from meds import allometry as al
# was: dbh_to_height(dbh, hgt_max)  with local B1HT/B2HT
h   = al.dbh_to_height(dbh, hgt_max)
agb = al.dbh_to_agb(new_dbh, al.dbh_to_height(new_dbh, hgt_max), rho)
```

The *laws themselves* (Camac mortality, intrinsic growth, competition suppression,
reproduction→recruits) stay here unchanged — they are the empirical example.

### 6.3 `empirical_spinup.py` — **one call site changes**

```python
# before
site.apply_rates(growth, mortality, recruitment, new_month, new_year)
# after
site.apply_demography(growth, mortality, recruitment, new_month, new_year)
```

Plus the `MEDS_LIB` path in the header. The read-state / rate-compute / golden-compare loop
is otherwise identical, and it must still reproduce
`test/golden/empirical_spinup_golden.csv` bit-for-bit (Part A is behavior-preserving).

### 6.4 `_cadence.py` — **unchanged.**

---

## 7. Verification

Build `-DMEDS_BUILD_PYLIB=ON` → `libmeds.so`, then:

1. `run_carbon.py` — matches the current output / `meds_main` (carbon path untouched).
2. `empirical_spinup.py` — reproduces `empirical_spinup_golden.csv` to the same tolerance as
   today. This is the load-bearing check: it proves the Python-orchestrated primitive
   sequence + `grow_cohorts_by_rate` are an exact relocation of the old shim body.
3. `meds.plant` leaf/phenology examples — unchanged symbols, still pass.

Because the C symbols and kernels are unchanged (only relocated/renamed) and the container
stays in Fortran, the goldens are a complete oracle.

---

## 8. Migration phases (behavior-preserving at every step)

- **P0 — Fortran primitive.** Add `grow_cohorts_by_rate` to core; have the *still-monolithic*
  shim call it in place of `compute_empirical_derivatives`; delete the private helper. Build,
  run `empirical_spinup` → golden identical. (Isolates the only real Fortran move.)
- **P1 — consolidate the library.** Move `meds_plant_capi.f90` into `src/capi/`, collapse to
  one `libmeds.so` target, delete `meds_plant_c`. Update both `_ffi.py` → one loader. Run all
  examples.
- **P2 — unpack the shim.** Split the monolith into the Part B files (registry first, then
  config/site/init/diagnostics/dynamics), each a pure relocation of existing symbols. Build,
  re-run examples (symbols unchanged → goldens identical).
- **P3 — empirical primitives + Python cadence.** Add the apply-primitive bindings
  (`meds_demography_capi`), replace `Site.apply_rates` with `Site.apply_demography`, update
  `empirical_spinup.py`. Golden identical.
- **P4 — allometry kernels.** Add `meds_allometry_capi` + `meds.allometry`; dedup
  `empirical_laws.py`. Golden identical (kernels are the same math the constants encoded).

Each phase is independently verifiable against the goldens; stop after any phase and the tree
still builds and passes.

---

## 9. Open decision (one)

**Empirical apply-cadence: Python-orchestrated primitives (this doc) vs. one Fortran call.**
This doc puts the empirical cadence in Python (Principle 1, and it makes the example fully
hackable). The alternative is a single public `apply_demographic_rates(site, cfg, rates,
flags)` in the driver that capi binds thinly — fewer boundary crossings and the cadence
lives in Fortran, at the cost of the empirical sequence no longer being Python-editable.
The recommendation is the Python-orchestrated version; flag if you'd prefer the single-call
form and Part B/D collapse `meds_demography_capi` back to one `meds_apply_rates` binding.
