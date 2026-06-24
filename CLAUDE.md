# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What MEDS is

MEDS (Modern and modular Ecosystem Demography Simulator) is a ground-up reimplementation, in
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

**Greenfield.** This repository currently contains documentation and git configuration only — there
is no Fortran source or build tree yet. The "Build" and "Layout" sections below describe the *target*
structure to grow into, not files that exist today. Keep this section honest: update it as real code
lands.

No Fortran compiler is installed in the current environment (`gfortran`/`ifort`/`ifx` all absent).
Install one (and HDF5/netCDF) before expecting builds to run. Two supported toolchains:
- **gfortran** + `libhdf5-dev` — ED2's reference platform; the open-source default. Run
  `./scripts/install_gfortran.sh` (add `--hdf5` to also install HDF5 dev files).
- **Intel `ifx`** (oneAPI LLVM Fortran) — run `./scripts/install_ifx.sh`, which checks for `ifx`
  and offers to install it via Intel's APT repo (Debian/Ubuntu/WSL). After install, activate it with
  `source /opt/intel/oneapi/setvars.sh`. **`ifx` needs CMake ≥ 3.20** (the system CMake here is
  3.16, too old to recognize the IntelLLVM compiler — install a newer one, e.g. `conda install
  -c conda-forge cmake`).

## Build (CMake — target workflow)

MEDS uses **CMake**. This is a deliberate fix for ED2's biggest build pain: ED2 hand-maintains object
lists and platform `include.mk.<platform>` files and has to run `make` six times to resolve Fortran
module ordering. CMake's Fortran support tracks `.mod` dependencies automatically — never reintroduce
manual dependency lists or repeated-build hacks.

```bash
# Configure (Debug = strict checks; Release = optimized science runs)
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j

# Optimized build
cmake -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j

# Run the test suite (CTest)
ctest --test-dir build --output-on-failure

# Run a single test by name (regex)
ctest --test-dir build -R <test_name> --output-on-failure
```

Select the compiler at configure time with `-DCMAKE_Fortran_COMPILER=gfortran` (or `ifx`); for
`ifx`, `source /opt/intel/oneapi/setvars.sh` first.

Build-type conventions (carry over ED2's intent: strict-for-dev, fast-for-science):
- **Debug** — gfortran: `-g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -std=f2018
  -Wall`; ifx: `-g -O0 -check all -traceback -fpe0 -warn all -stand f18`. Use while developing; the
  FPE trap catches the NaN/Inf bugs that silently corrupt ED2 runs.
- **Release** — gfortran `-O2`/`-O3`; ifx `-O2`/`-O3 -xHost`. No checks; only for validated code.
Toggle MPI and OpenMP via cache options (e.g. `-DMEDS_USE_MPI=ON -DMEDS_USE_OPENMP=ON`), not by
editing flags inline. Discover HDF5/netCDF/MPI with `find_package`, never hard-coded `-I`/`-L` paths.

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
