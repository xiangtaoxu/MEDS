# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What MEDS is

MEDS (**Modular Ecosystem Dynamics Simulator**, also read as **Modern Ecosystem Demography
Simulator** to flag its descent from the ED / ED2 lineage) is a ground-up reimplementation, in
**Fortran 2018**, of the Ecosystem Demography model **ED2** (https://github.com/EDmodel/ED2). It
keeps ED2's science — a size- and age-structured (demographic) terrestrial biosphere model — while
replacing the legacy code structure with modular, testable, standards-conformant Fortran.

The goal is *not* a line-by-line translation. It is to re-express ED2's proven algorithms with
modern language features (modules + submodules, derived-type encapsulation, explicit interfaces,
`allocatable` over `pointer`, parameterized real kinds) so the model is easier to extend, test, and
reason about. When in doubt about *what* a process should compute, the ED2 source is the reference;
when deciding *how* to structure it, follow the modernization guidelines below — do not copy ED2's
structure.

### Reference implementation
The original ED2 tree is checked out as a **sibling directory**: `../ED2`. The core model lives in
`../ED2/ED/src` (driver, dynamics, init, io, memory, mpi, utils). Treat it as read-only reference.
Key entry points to read when porting a feature:
- `../ED2/ED/src/memory/ed_state_vars.F90` — the entire 34k-line state hierarchy (see below).
- `../ED2/ED/src/driver/ed_model.F90` — the main time-integration loop.
- `../ED2/ED/src/driver/edmain.F90` → `ed_driver.F90` — program entry and dispatch.

Background papers: Moorcroft et al. 2001 (Ecol. Monogr.); Medvigy et al. 2009 (JGR);
**Longo et al. 2019 (GMD 12:4309)** — the definitive ED-2.2 technical description (read the
Supporting Information). ED2 wiki: https://github.com/EDmodel/ED2/wiki.

## Repository status

**Implemented: the standalone demographic core** (`src/`, `app/`, `test/`, `CMakeLists.txt`). It runs
cohort & patch dynamics — individual-tree growth, mortality, recruitment, cohort/patch fusion/fission,
and treefall **patch disturbance** — driven by demographic *rates supplied from outside* the engine as
three plain arrays. Size follows the pan-tropical (ED2 `iallom==3`) allometry (`meds_allometry`); each
cohort carries **AGB (carbon)** and **leaf area**, and cohort fusion/fission conserve total AGB. There
is deliberately **no** radiative transfer or full biogeochemistry yet — the rates are the
phenomenological (structure-only) relationships in `meds_phenomenological_vital_rates` (light
competition through overtopping LAI), and the seam for a mechanistic replacement is the data-array
interface `update_demography`. See "Demographic core" below.

Toolchain on this machine (installed, but **off the default PATH** — activate before building):
- **Intel `ifx` 2026** — `source /opt/intel/oneapi/setvars.sh`. Strict-standards CPU compiler; runs the
  full test suite. (`scripts/install_ifx.sh` installs it elsewhere.)
- **NVIDIA `nvfortran` 25.11** (HPC SDK) — add `…/hpc_sdk/Linux_x86_64/25.11/compilers/bin` to PATH.
  This is the parallel/GPU path. The hot kernels (`meds_demography_dynamics`) carry explicit OpenMP
  `target` regions over plain arrays, so the build picks the device via the NVHPC `-mp` flag
  (`MEDS_GPU=gpu` → `-mp=gpu -gpu=mem:separate`; `MEDS_GPU=multicore` → `-mp`; no flag → serial). All
  three back ends are validated on the RTX 3050 Ti: ifx 7/7, nvfortran multicore 7/7, nvfortran GPU
  7/7 (CPU↔GPU results identical). **Do NOT use `-stdpar=gpu`**:
  it forces the global CUDA-managed allocator, whose deep-copy/finalize of the allocatable-component
  `site` type double-frees on the host. OpenMP `target` + `-gpu=mem:separate` keeps all state in
  normal host memory and moves only the mapped arrays.
- **CMake 4.3.4** (conda, on PATH) — recent enough for both IntelLLVM and NVHPC compiler IDs.
- `gfortran` is not installed (`scripts/install_gfortran.sh` adds it); the build supports it too.
- **Portability trap (nvfortran):** never pass an array-valued *function result* straight into a call
  (`call foo(bar(x))` where `bar` returns an array). nvfortran's whole-program optimizer miscompiles the
  temporary descriptor — **silently wrong values at `-O2`/`-fast`, segfault at `-O0`** — while
  `ifx -stand f18 -check all` tolerates it (only an `arg_temp_created` remark), so a green ifx suite
  hides it. Bind to a named array first (`tmp = bar(x); call foo(tmp)`); this also clears the ifx remark.
  Corollary: a green ifx run is **not** sufficient — build the nvfortran multicore back end on new
  modules too. (See issue #7; found porting `src/biophys`.)

## Build (CMake)

CMake auto-resolves Fortran `.mod` dependencies — this is the deliberate fix for ED2's "run `make` six
times" hack and per-platform `include.mk` files. Never reintroduce manual object lists or repeated builds.

netCDF output is ON by default, so the standard build needs the netCDF C library (its CMake config
also pulls HDF5 -> the build enables the C language for the IO target); `scripts/install_netcdf.sh`
installs it (conda or apt) and prints the `-DCMAKE_PREFIX_PATH` to pass. Add `-DMEDS_ENABLE_IO=OFF` for
a netCDF-free build (links a no-op stub I/O layer); that build, the engine, and the tests have no
external dependency. Release builds the C-binding IO layer cleanly (Debug `-check all` emits a harmless
arg_temp_created remark per netCDF call).

```bash
# Default build WITH netCDF (Intel ifx). Point at the netCDF-C install via CMAKE_PREFIX_PATH.
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake --build build-ifx -j
ctest --test-dir build-ifx --output-on-failure          # 8 tests
LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib ./build-ifx/meds_main meds_config_main.toml  # -> <prefix>-D-output.nc (+ -S-<ts>.nc if write_state)

# netCDF-free test/debug build (no netCDF/C needed; strict -check all). Quick compiler/engine checks.
cmake -S . -B build-debug -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug -DMEDS_ENABLE_IO=OFF
cmake --build build-debug -j && ctest --test-dir build-debug --output-on-failure

# Multicore CPU / GPU via OpenMP target (NVIDIA nvfortran); add -DMEDS_ENABLE_IO=OFF to skip netCDF
export PATH=…/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
cmake -S . -B build-mc  -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=multicore \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=gpu \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake --build build-gpu -j

# Run a single test by regex
ctest --test-dir build-debug -R fusion_cohort --output-on-failure
```

`MEDS_GPU` (`none|multicore|gpu`) only affects NVHPC builds; the `-mp` flags are `PUBLIC` on the
`meds_demography` target so `meds_main`/tests inherit the offload compile+link flags. ifx/gfortran ignore
it — the `!$omp` lines are comments without an OpenMP flag, so those builds stay serial. Per-compiler
Debug flags live in the `meds_fortran_flags()` function in `CMakeLists.txt` (ifx `-stand f18 -check all
-fpe0`; gfortran `-std=f2018 -fcheck=all -ffpe-trap`; nvfortran `-Mbounds -Ktrap=fp`).

## Demographic core

### Source layout & libraries
The source tree is organised by process domain over a **state/process wall** (an acyclic library
DAG `shared ← {allometry, plant} ← state ← demography ← aux ← main`, where `plant` links `shared`
ONLY and is orthogonal to `state`); moving a file changes only CMake wiring because Fortran `use` is
by module name and all `.mod`s share one directory. **The 2026-07-04 plant refactor** flattened
`src/plant/` into one ecophysiology library and moved the empirical vital rates into `demography`
(making it self-contained); design: `archive/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md`.
- **`src/shared/`** → `libmeds_shared.a` — the foundation, NOT tied to any process: `meds_kinds`
  (precision), `meds_constants`, `meds_pft_params` (the PFT trait table, incl. per-PFT `hgt_max`),
  `meds_time` (calendar + leap-year-aware Gregorian arithmetic), `meds_config` (run config), and
  `meds_temp_response` (Arrhenius / peaked deactivation — promoted here from the leaf module so leaf,
  respiration and any tissue share one code path without a plant→plant library edge). Root of the DAG.
- **`src/allometry/meds_allometry.f90`** → `libmeds_allometry.a` — pan-tropical (`iallom==3`)
  size↔height↔AGB↔leaf-area relations. A shared structural-geometry foundation used by BOTH `state`
  (cohort geometry caching / fusion via `set_cohort_size`/`agb_to_dbh`) and the plant ecophysiology
  library, so it is its OWN library BELOW `state` — it cannot live in `libmeds_plant` without making
  the demographic core depend on all of ecophysiology (see issue #11). `hgt_max` is a per-PFT argument
  to `dbh_to_height`/`agb_to_dbh`.
- **`src/state/meds_demography_types.f90`** → `libmeds_state.a` — the cohort/patch STATE hub: the
  flat site-wide Structure-of-Arrays (`cohort_block` + patch CSR) and the ONE centralized lockstep
  `cohort_reorder`/`rebuild_csr`/`set_cohort_size` machinery. Links `shared` + `allometry`.
- **`src/demography/`** → `libmeds_demography.a` — a SELF-CONTAINED demography model (links `state` +
  `allometry`; NO plant-ecophysiology dependency, so it runs the spin-up on its own). The PROCESS
  engine — `meds_demography_interface` (the data-rate seam `update_demography`), `meds_demography_dynamics`
  (the OpenMP-target `growth_step`/`mortality_step` kernels + treefall + recruitment),
  `meds_demography_structure` (sort + cohort/patch fuse/fission), `meds_demography_diagnostics` — PLUS
  the phenomenological empirical RATE provider moved here from the former `plant/vitality/`:
  per-individual formulas `meds_growth`/`meds_mortality`/`meds_recruitment` behind the
  `meds_phenomenological_vital_rates` orchestrator (the `vital_rates` seam — the drop-in point for a
  future mechanistic provider).
- **`src/plant/`** → `libmeds_plant.a` — ONE flat, self-contained plant-ECOPHYSIOLOGY library (links
  `meds_shared` only; NO `site_t`, compiles standalone via `cmake --build … --target meds_plant`). All
  derived types are consolidated in **`meds_plant_types`**. It holds: **leaf gas exchange** — the seam
  `meds_leaf_physiology%leaf_gas_exchange(env, cfg, ipft, flux)` over `meds_leaf_photosynthesis`
  (FvCB C3 + Collatz C4), `meds_leaf_stomata` (Leuning / Medlyn / Katul), `meds_leaf_solver` (bracketed
  Ci root-find); **hydraulics** (`meds_plant_hydraulics` + `meds_hydro_*`); **phenology**
  (`meds_plant_phenology` + `meds_pheno_engine`); and (to come) `meds_plant_respiration`. The optional
  Python C-API (`meds_plant_capi.f90`, `-DMEDS_BUILD_PYLIB=ON` → `libmeds_plant_c`, GLOB
  `src/plant/*_capi.f90`) is compiled only into the shared lib, exposed through the `meds.leaf` Python
  package (reproduces Slot & Winter 2017 in `examples/example_leaf_physiology/`). NOT yet wired into
  the demographic stepper.
- **`src/biophys/`** → `libmeds_biophys.a` — self-contained canopy radiative transfer (ED2 two-stream
  `icanrad=2`), links `shared` only; a sibling stateless-kernel library to `plant`. `src/biogeochem/`
  and `src/utils/` remain empty placeholders.
- **`src/driver/`, `src/init/`** → all part of `libmeds_aux.a` — the top-level utilities that wire the
  process modules together: `meds_stepper` (the master stepper, `src/driver`; seed of a future
  all-process **master loop**, ED2-`ed_model` analogue) and `meds_init` (`src/init` — the
  initial-community builders: `init_bare_ground`, `add_cohort`, and `init_from_census`). The `meds_aux`
  target globs `src/driver/*.f90` + `src/init/*.f90`, EXCLUDING `src/driver/meds_main.f90`, and links
  `meds_demography` + `meds_vitality` + `meds_config_io`.
- **`src/driver/meds_main.f90`** → the executable `meds_main`, the single entry point (read config →
  build community → run → save output → exit; merged from the former `app/meds_demo` + `app/meds_io_demo`).
  Excluded from `meds_aux` (it is a PROGRAM) and built as the `meds_main` target, linking `meds_aux` +
  the I/O library (real or stub, below).
- **`src/io/`** — netCDF I/O via the netCDF **C** library through `iso_c_binding` (`meds_netcdf_c`
  bindings + `meds_io`, `libmeds_io.a`), built **by default**. `meds_io` has two output streams and a
  restart reader: a DIAGNOSTIC timeseries (`io_create`/`io_write_snapshot`/`io_close` ->
  `<prefix>-D-output.nc`, with derived diagnostics), instantaneous STATE checkpoints
  (`io_write_state` -> `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, prognostic state only — no diagnostics, cached
  geometry re-derived on read), and `io_read_state` (restart). `enable_language(C)`
  fires here so netCDF-C's config can resolve HDF5/Threads. netCDF-Fortran is unavailable for
  ifx/nvfortran here (its `.mod` is gfortran-only); the C API is compiler-agnostic. Point CMake at the
  netCDF-C install with `-DCMAKE_PREFIX_PATH=<prefix>` (here the
  conda env `~/miniforge3/envs/common`). When IO is OFF, `meds_io_stub.f90` provides a no-op module of
  the SAME name (`meds_io`) so `meds_main` always builds and runs; keep the two interfaces in sync. The
  core engine/tests never reference netCDF either way. Also here: the always-on TOML config reader
  (`meds_toml` + `meds_config_io`, `libmeds_config_io.a`, no external deps), which also writes the
  per-PFT parameter table to `<output_dir>/<prefix>_pft_parameters.csv` (`write_pft_params_csv`,
  netCDF-free, one row per PFT — a provenance record of the run's trait values).
  - **No hard-coded model parameters.** The source defines only true constants (`meds_constants`,
    plus `meds_allometry`'s coefficients which are `protected` module vars *installed at load* via
    `set_allometry`). Every parameter is REQUIRED from TOML across TWO files — a MAIN file (all non-PFT
    settings, names the PFT file via `[init].pft_config`) and a PFT file (`[pft]` traits + `[camac]`
    mortality coefficients + `[allometry]` coefficients). `load_meds_config` builds a **presence map**
    while reading and `error stop`s listing every missing key; a missing file is also a hard error.
    There is no `build_config`/defaults: derived quantities come from `derive_config` + `derive_pft_rates`
    (overridable via `[options].override_derived` + a `[derived]` block). Tests get a complete config
    from `build_test_config()` in `test/meds_test_support.f90` (the only place "default" values live in
    code). The offloaded `growth_step` takes the allometry coefficients as **firstprivate scalar args**
    (it can't read host module vars on the device); host allometry functions read them from the module.
- Build order via link deps: `meds_shared → meds_demography → meds_aux` (+ optional `meds_io`). The
  demography engine compiles standalone (`cmake --build <dir> --target meds_demography`).

### Invariants to build on when extending the engine
- **State = flat site-wide Structure-of-Arrays** (`meds_demography_types`): all cohorts of the whole
  site in one contiguous set of 1-D arrays (`cohort_block`), patch membership as a CSR map
  (`cohort_offset`/`cohort_count` + `owner_patch`). The dominant daily kernels are a single
  unit-stride sweep.
- **OpenMP `target` over plain arrays for the hot kernels** (`meds_demography_dynamics`: `growth_step`,
  `mortality_step`) — they take bare arrays (no `site_t`, no derived types), so the `map` clauses are
  clean and the host keeps all state in normal memory. Keep them arithmetic-only (intrinsics only).
  All restructuring (sort/fuse/split/terminate/recruit) and the rate evaluation are **host-only**.
- **Rates arrive as plain DATA** through `update_demography` (`meds_demography_interface`): three arrays
  — per-cohort growth `[cm/yr]`, per-cohort total mortality `[1/yr]`, per-(PFT,patch) recruitment —
  plus `dt_yr` and the structural triggers. The engine never computes a rate; the empirical rate layer
  (in `src/demography/`) does: `meds_phenomenological_vital_rates%vital_rates` produces all three from
  structure alone — **growth** = intrinsic (a capped log-linear function of dbh) × competition
  suppression (`exp(growth_lai_slope·overtopping LAI)`) × reproductive-allocation suppression;
  **mortality** = the Camac-2018 additive hazard `mort_gamma + mort_alpha·exp(−mort_beta·growth_avg)`,
  where `growth_avg` is the cohort's tracked simple moving-average growth (window `growth_memory_days`,
  a prognostic per-cohort field); and **recruitment** = baseline `seed_rain_recruits` plus a
  reproduction flux (the carbon diverted from growth to reproduction, per PFT and patch, converted to
  min-size recruits). The engine applies recruitment via `apply_recruitment` (pool + cohort spawn).
  There is
  **no environmental driver** (no temperature) in the current model. A mechanistic module produces the
  same arrays with no engine change. (There is no `rate_provider_t` — data arrays, not a class seam, so
  the engine stays callable from Python.)
- **Conserved invariants** (every fuse/split asserts within 1%, else `error stop`): cohort fusion/split
  conserve **total aboveground biomass (carbon)** and plant number — `dbh` is *re-derived* from the
  conserved per-plant AGB (`agb_to_dbh`), never averaged; patch fusion conserves **site-level plant
  number** via area-fraction rescaling, and patch area always renormalizes to 1. Patch **disturbance**
  conserves area (donors shed a fraction into a new age-0 gap) but intentionally does NOT conserve
  plant number — the killed canopy is the disturbance.
- **One centralized lockstep reorder** (`meds_demography_types`: `cohort_reorder`/`cohort_compact`/
  `copy_cohort_slot`/`rebuild_csr`/`cohort_ensure_capacity`/`move_alloc_block`, plus `set_cohort_size`
  which fills the cached height/BA/AGB/leaf-area of one slot). When you add a per-cohort field, update
  *these* — the single place that touches every array (the fix for ED2's "forgot to reallocate" class).
  (Patch arrays have no single reorder routine; their permute/pack sites are `sort_patches` and
  `patch_compact` in `meds_demography_structure` — update both when adding a per-patch field.)
- **Persistent identity** (`global_id` on `cohort_block` and `patch_index`, monotonic `next_*_id`
  counters on `site_t`): every cohort/patch is stamped at creation via `assign_cohort_id`/
  `assign_patch_id` and carries that id, in lockstep with all other fields, through every
  sort/fusion/compaction; ids are never reused. Fusion keeps the *survivor's* id (the absorbed one
  disappears); a split daughter, a recruit, and a disturbance-gap fragment each get a *fresh* id. This
  lets an external reader (e.g. `post_proc/plot_forest_structure.py`) track one cohort/patch across
  output records until it fuses away or is culled. Creation sites that must stamp ids: `add_cohort`/
  `init_bare_ground` (`meds_init`), `apply_recruitment`, `split_cohorts`, `apply_patch_disturbance`.
- **Order of operations** (`meds_stepper%advance_one_step` in `src/driver/`; cadence flags from the
  caller's calendar): every step `growth → mortality → patch-age`; monthly (`do_cohort_fissfuse`)
  `recruit → cohort fuse/terminate/split → sort`; annual `disturbance` (`do_patch_disturbance`) then
  patch restructuring (`do_patch_fissfuse`: `patch sort/fuse/terminate → cohort consolidate`) — the two
  are INDEPENDENT triggers. Disturbance integrates its yearly rate over the 1-year patch-dynamics
  interval, NOT the per-step `dt_yr`. Same routine serves daily/weekly/monthly steps.
- Cohort fusion/fission keys on **height & LAI** (size anchor is height): similarity is
  `|Δheight| < hgt_max·tol` (the cohort's per-PFT height cap) with geometric tolerance relaxation up to `max_cohort`, and a cohort
  splits (or is blocked from fusing) when its LAI exceeds `cohort_lai_cap`. **Patch fusion** compares a
  cumulative-LAI **light profile** (ED2 scheme, `patch_light_profile`): per-height-layer
  `light = exp(-light_ext·LAI_above)`, fused when the mean AND max layer light difference are within
  `patch_light_tol` / `patch_light_maxdev_factor`, same disturbance class only.
- **PFT contrast is wood density** (`meds_pft_params`). Low ρ ⇒ higher mortality hazard (the
  Camac-2018 `mort_gamma`/`mort_alpha`/`mort_beta` are a power law in ρ, `param_0·(ρ/0.6)^exp`, always
  positive). The growth,
  competition and reproduction parameters (`growth_dbh_*`, `growth_lai_slope`,
  `reproduction_investment_fraction`, `repro_carbon_efficiency`, `seed_rain_recruits`) are per-PFT but
  default to uniform values; the two height thresholds (`min_cohort_height` = recruit/birth size,
  `min_reproduction_height` = reproductive maturity) are shared. Add traits here; derive ρ-dependent
  rates in `derive_pft_rates`.

### Reserved follow-ups (not yet implemented)
A top-level master loop over all processes; an `!$omp target data` region keeping the cohort arrays
device-resident across the daily loop (cuts the per-step map overhead that currently makes the GPU
spin-up migration-bound); a `bind(c)` C-API + shared library for Python (`ctypes`/`cffi`) for the
DEMOGRAPHIC engine — `f2py` will not handle the derived-type/allocatable design, and the data-array
interface (not a Fortran class) is the intended foreign-call layer (the PLANT module already has this:
`src/plant/meds_plant_capi.f90` + `-DMEDS_BUILD_PYLIB=ON` → `libmeds_plant_c`, exposed through
the `meds.leaf` Python package (`python/meds/leaf`, a clean ctypes-free API installed with
`pip install -e python/`), exercised by `examples/example_leaf_physiology/reproduce_slot2017.py`);
and **coupling the leaf-physiology module into the demographic
growth** — `src/plant/` exists as a standalone plant-ecophysiology library (FvCB C3 + Collatz
C4, Leuning/Medlyn/Katul stomata, Arrhenius/peaked temperature response), but wiring its assimilation
into the rate provider needs the still-missing canopy radiative transfer (per-cohort absorbed PAR),
leaf energy balance (leaf temperature), a meteorological forcing source, and plant hydraulics
(`psi_leaf`). (Done since the first cut: a single `meds_main` entry point
(`src/driver`); netCDF output — `src/io`, on by default (`-DMEDS_ENABLE_IO=OFF` for a netCDF-free
build) — split into a diagnostic timeseries and instantaneous STATE checkpoints; initialization from a
cohort census or a restart state file; recruitment moved to the physiology rate layer; temperature
(environmental control) removed; the cumulative-LAI light patch-fusion; the wood-density PFT axis;
treefall disturbance; a real leap-year calendar (`meds_time`) — the run is bounded by `start_time`/
`end_time` dates (not a year count), `meds_main` walks the Gregorian calendar to set the monthly/annual
cadence, both output streams stamp the simulated date, and state files are named `-S-<sim-date>.nc`.)

## Architecture (inherited from ED2)

The science is organized around one central idea: **ecosystem state is a nested demographic
hierarchy**, and the model advances it on **two coupled timescales**. Understanding these two things
is the key to the whole codebase.

### 1. The demographic state hierarchy
State is a tree of nested, ragged arrays. In ED2 this is the monolithic `ed_state_vars.F90`
(`edtype` → `polygontype` → `sitetype` → `patchtype`, ~34k lines). In MEDS this becomes **one focused
module per level** with explicit allocate/deallocate/copy routines (ideally type-bound):

- **grid** — bookkeeping container linking polygons.
- **polygon** — a location sharing one meteorological forcing.
- **site** — shares a soil column and ground hydrology.
- **patch** — shares a disturbance history and stand age (a successional cohort of land).
- **cohort** — plants of one PFT in one height bin (the demographic atom).

Height and age are continuous but **binned**; "identical" means same bin. Bins are defined
dynamically, which is why **fusion/fission** (merging/splitting cohorts and patches to control
resolution) is a first-class operation — see `../ED2/ED/src/utils/fuse_fiss_utils.f90` and
`fusion_fission_coms.f90`. Plant functional types (**PFTs**) parameterize cohort behavior; PFT traits
live in `pft_coms.f90`.

### 2. The two-timescale integration loop
`ed_model.F90` drives the run. Per polygon:

- **Fast loop (sub-hourly, `DTLSM`)** — biophysics: radiation (two-stream), canopy air space,
  photosynthesis + stomatal conductance (Farquhar–Leuning and Katul variants), respiration,
  energy/water/CO₂ budgets, soil/leaf thermodynamics, plant hydraulics. These are stiff ODEs solved
  by a **selectable integrator**: Euler / RK4 / Heun / hybrid (`integration_scheme`; see
  `rk4_*`, `euler_driver.f90`, `heun_driver.f90`, `hybrid_driver.f90`, `bdf2_solver.f90`).
- **Slow loop (daily/monthly/yearly)** — vegetation dynamics (`veg_dynamics_driver`): carbon
  allocation and structural growth, phenology, mortality, reproduction/recruitment, and disturbance
  (treefall, fire, land use / forestry). This is where the demographic structure actually changes,
  followed by fusion/fission to keep cohort/patch counts bounded.

Process families and their ED2 source (under `../ED2/ED/src/dynamics` unless noted):
- Photosynthesis: `farq_leuning.f90`, `farq_katul.f90`, `photosyn_driv.f90`
- Radiation: `twostream_rad.f90`, `multiple_scatter.f90`, `radiate_driver.f90`
- Canopy/turbulence: `canopy_struct_dynamics.f90`
- Phenology: `phenology_driv.f90`, `phenology_aux.f90`
- Growth/allocation: `growth_balive.f90`, `structural_growth.f90`, `stem_resp_driv.f90`
- Mortality / reproduction: `mortality.f90`, `reproduction.f90`
- Disturbance: `disturbance.f90`, `fire.f90`, `forestry.f90`
- Soil C/N + hydrology: `soil_respiration.f90`, `plant_hydro.f90`, `lsm_hyd.f90`, `decomp_coms.f90`

### 3. Parameters, I/O, and parallelism
- **Parameters** — ED2 holds tunable parameters/constants in many global `*_coms` modules
  (`consts_coms.F90`, `pft_coms.f90`, `soil_coms.F90`, `physiology_coms.f90`, …). In MEDS these
  become **derived-type configuration objects passed explicitly** (see modernization guidelines) so
  state is not hidden global mutable data.
- **I/O** — ED2 reads an `ED2IN` Fortran namelist plus optional XML config, and writes HDF5 history
  (restart) and analysis output (`../ED2/ED/src/io`). Keep I/O behind a dedicated layer; isolate the
  HDF5/netCDF dependency so the science modules never touch file formats directly.
- **Parallelism** — MPI distributes polygons across ranks (master/slave, guarded by `-DRAMS_MPI`),
  with OpenMP threading inside a rank. Both are optional and must build cleanly when disabled.

## Modernization guidelines (ED2 → MEDS)

These are the rules that make MEDS "modern and modular." They override ED2's structure wherever they
conflict.

- **One responsibility per module; one module per file.** Never reproduce a 34k-line file. Put type
  definitions in `*_types` modules and heavy procedures in **submodules** to cut compile coupling.
- **`implicit none` everywhere; explicit interfaces via modules.** Drop ED2's `USE_INTERF` flag that
  *bypassed* interface checking — the compiler should check every call.
- **Parameterized kinds.** A single `shared`/`kinds` module defines `sp`, `dp` (and integer kinds);
  use `real(dp)` etc. Never write bare `real*8`, `kind=8`, or `1.0d0` literals scattered in code.
- **No hidden global mutable state.** Replace mutable `*_coms` globals with derived-type config/state
  passed as arguments. This is what makes routines unit-testable, thread-safe, and reentrant.
- **`allocatable` over `pointer`** for ownership; use `move_alloc` and automatic deallocation instead
  of manual bookkeeping. Reserve `pointer` for genuine linked structures.
- **`pure`/`elemental`** for leaf math kernels (photosynthesis, allometry, thermodynamics) — they are
  the easiest to test and the most reusable.
- **Errors via `error stop` / status codes**, not ED2's `fatal_error` halt pattern, so failures are
  catchable in tests.
- **≤132 columns, free-form, lowercase keywords.** ED2 is intentionally 132-column compliant; keep
  that — do not rely on `-ffree-line-length-none` to mask over-long lines.
- **Readability over succinctness when performance is equal**, especially for ecologically/physically
  meaningful names. Spell out `site`, `cohort`, `patch`, `dbh_critical`, `hgt_max` (the site-level
  type is `site_t`, the instance `site`; containers are `site%cohort` / `site%patch`); keep terse
  *loop indices* (`i`, `ip`, `pf`/`ipft`) and the established domain tokens `dbh`, `nplant`, `pft`.
  Inside `associate`, alias to the full word (`cohort => site%cohort`, `patch => site%patch`,
  `pft => cfg%pft`), not single letters.
- **Test as you port.** Each ported kernel gets a unit test (CTest target); validate whole-model
  behavior against ED2 outputs (the analog of ED2's `EDTS` regression suite). A port is not done
  until it reproduces the reference within tolerance.

## Git conventions

Track Fortran source (`*.f90 *.F90`), CMake files (`CMakeLists.txt`, `*.cmake`), namelists, and docs.
Build artifacts (`*.o *.mod *.smod *.a`, the `build*/` trees) are ignored — see `.gitignore`. Line
endings are normalized to LF via `.gitattributes`. Keep generated HDF5/netCDF output and large input
datasets out of the repo.
