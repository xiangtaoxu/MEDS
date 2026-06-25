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
cohort & patch dynamics — individual-tree growth, mortality, recruitment, and cohort/patch
fusion/fission — driven by demographic *rates supplied from outside* the engine. There is deliberately
**no** radiative transfer, biophysics, or carbon balance yet (those rates are currently the test
empirical relationships in `meds_rates_empirical`; the seam for mechanistic rates is `rate_provider_t`).
See "Demographic core" below.

Toolchain on this machine (installed, but **off the default PATH** — activate before building):
- **Intel `ifx` 2026** — `source /opt/intel/oneapi/setvars.sh`. Strict-standards CPU compiler; runs the
  full test suite. (`scripts/install_ifx.sh` installs it elsewhere.)
- **NVIDIA `nvfortran` 25.11** (HPC SDK) — add `…/hpc_sdk/Linux_x86_64/25.11/compilers/bin` to PATH.
  This is the parallel/GPU path. The hot kernels (`meds_demography_kernels`) carry explicit OpenMP
  `target` regions over plain arrays, so the build picks the device via the NVHPC `-mp` flag
  (`MEDS_GPU=gpu` → `-mp=gpu -gpu=mem:separate`; `MEDS_GPU=multicore` → `-mp`; no flag → serial). All
  three back ends are validated on the RTX 3050 Ti: ifx 5/5, nvfortran multicore 5/5, nvfortran GPU
  (offload codegen confirmed via `-Minfo=mp`; CPU↔GPU results identical). **Do NOT use `-stdpar=gpu`**:
  it forces the global CUDA-managed allocator, whose deep-copy/finalize of the allocatable-component
  `site` type double-frees on the host. OpenMP `target` + `-gpu=mem:separate` keeps all state in
  normal host memory and moves only the mapped arrays.
- **CMake 4.3.4** (conda, on PATH) — recent enough for both IntelLLVM and NVHPC compiler IDs.
- `gfortran` is not installed (`scripts/install_gfortran.sh` adds it); the build supports it too.

## Build (CMake)

CMake auto-resolves Fortran `.mod` dependencies — this is the deliberate fix for ED2's "run `make` six
times" hack and per-platform `include.mk` files. Never reintroduce manual object lists or repeated builds.

```bash
# CPU, strict checks (Intel ifx) — runs the test suite
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug
cmake --build build-ifx -j
ctest --test-dir build-ifx --output-on-failure          # 5 tests
./build-ifx/meds_demo [years]                            # spin-up demo (default 60 yr)

# Multicore CPU via OpenMP target on host (NVIDIA nvfortran)
export PATH=…/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
cmake -S . -B build-mc -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=multicore

# GPU offload of the OpenMP-target kernels (mem:separate => no managed memory)
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=gpu
cmake --build build-gpu -j

# Run a single test by regex
ctest --test-dir build-ifx -R fusion_cohort --output-on-failure
```

`MEDS_GPU` (`none|multicore|gpu`) only affects NVHPC builds; the `-mp` flags are `PUBLIC` on the
`meds_demography` target so the demo/tests inherit the offload compile+link flags. ifx/gfortran ignore
it — the `!$omp` lines are comments without an OpenMP flag, so those builds stay serial. Per-compiler
Debug flags live in the `meds_fortran_flags()` function in `CMakeLists.txt` (ifx `-stand f18 -check all
-fpe0`; gfortran `-std=f2018 -fcheck=all -ffpe-trap`; nvfortran `-Mbounds -Ktrap=fp`).

## Demographic core

### Source layout & libraries
- **`src/shared_util/`** → `libmeds_shared.a` — cross-cutting, NOT demography-specific: `meds_kinds`
  (precision; every future module needs it), `meds_constants`, `meds_pft_params` (allometry + PFT
  trait table), `meds_config`. Closed dependency set (no reference back into demography).
- **`src/demography/`** → `libmeds_demography.a` — the demographic mechanism (building-block
  operations), links `meds_shared`. Modules: `meds_demography_interface` (the rate seam),
  `meds_demography_types`, `meds_demography_kernels`, `meds_sort`, `meds_cohort_dynamics`,
  `meds_patch_dynamics`, `meds_recruitment`. Generic-named modules (`meds_sort`, …) may gain the
  `meds_demography_` prefix later.
- **`src/` root (interim `libmeds_aux.a`)** — not part of the demography mechanism, home TBD:
  `meds_step` (the stepper; seed of a future all-process **master loop**, ED2-`ed_model` analogue),
  `meds_diagnostics` (future dedicated IO/diagnostics module), `meds_setup` (community init helpers),
  `meds_rates_empirical` (the example provider). The `meds_aux` target is an explicit placeholder.
- Build order via link deps: `meds_shared → meds_demography → meds_aux`. The demography engine
  compiles standalone (`cmake --build <dir> --target meds_demography`).

### Invariants to build on when extending the engine
- **State = flat site-wide Structure-of-Arrays** (`meds_demography_types`): all cohorts of the whole
  site in one contiguous set of 1-D arrays (`cohort_block`), patch membership as a CSR map
  (`cohort_offset`/`cohort_count` + `owner_patch`). The dominant daily kernels are a single
  unit-stride sweep.
- **OpenMP `target` over plain arrays for the hot kernels** (`meds_demography_kernels`: `growth_step`,
  `mortality_step`) — they take bare arrays (no `site`, no derived types), so the `map` clauses are
  clean and the host keeps all state in normal memory. Keep them arithmetic-only (intrinsics only).
  All restructuring (sort/fuse/split/terminate/recruit) and the rate evaluation are **host-only**.
- **Rates come from `rate_provider_t`** (`meds_demography_interface`) — an abstract type with deferred
  `pure` scalar bindings (kept a clean Fortran seam; do not overload it with C/Python concerns).
  `meds_rates_empirical` is the test implementation (growth = size×competition, mortality = f(growth),
  constant recruitment). A mechanistic provider extends the same type with no engine change.
- **Conserved invariants** (every fuse/split asserts within 1%, else `error stop`): cohort fusion/split
  conserve **total basal area** (the diameter analogue of ED2's biomass) and plant number — `dbh` is
  *re-derived* `= sqrt(basarea/pio4)`, never averaged; patch fusion conserves **site-level plant number**
  via area-fraction rescaling, and patch area always renormalizes to 1.
- **One centralized lockstep reorder** (`meds_demography_types`: `cohort_reorder`/`cohort_compact`/
  `copy_cohort_slot`/`rebuild_csr`/`cohort_ensure_capacity`/`move_alloc_block`). When you add a
  per-cohort field, update *these* — the single place that touches every array (the fix for ED2's
  "forgot to reallocate an array" bug class).
- **Order of operations** (`meds_step%advance_one_step` — currently at `src/` root; cadence flags from
  the caller's calendar): daily `competition → growth → mortality accumulate`; monthly `patch-age →
  apply mortality → recruit → cohort fuse/terminate/split → sort`; annual `patch fuse → patch
  terminate → cohort consolidate`. Same routine serves daily (default) and monthly time-steps
  (`TS_DAILY`/`TS_MONTHLY` in `meds_config`).
- Fusion/fission thresholds are **diameter & size-distribution** based (no LAI): cohort similarity is
  `|ΔDBH| < dbh_crit·tol ∧ |Δhite| < hgt_max·tol` with geometric tolerance relaxation up to
  `maxcohort`; patches compare normalized cumulative-basal-area profiles by DBH class.

### Reserved follow-ups (not yet implemented)
Patch-level disturbance (`meds_disturbance` + a `disturbance_if` deferred binding); a dedicated
IO/diagnostics module; a top-level master loop over all processes; and a `bind(c)` C-API + shared
library for Python (`ctypes`/`cffi`) — `f2py` will not handle the derived-type/allocatable/polymorphic
design, and the abstract rate interface is not the foreign-call layer.

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
- **Test as you port.** Each ported kernel gets a unit test (CTest target); validate whole-model
  behavior against ED2 outputs (the analog of ED2's `EDTS` regression suite). A port is not done
  until it reproduces the reference within tolerance.

## Git conventions

Track Fortran source (`*.f90 *.F90`), CMake files (`CMakeLists.txt`, `*.cmake`), namelists, and docs.
Build artifacts (`*.o *.mod *.smod *.a`, the `build*/` trees) are ignored — see `.gitignore`. Line
endings are normalized to LF via `.gitattributes`. Keep generated HDF5/netCDF output and large input
datasets out of the repo.
