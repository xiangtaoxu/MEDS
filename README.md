# MEDS — Modern and modular Ecosystem Demography Simulator

MEDS is a ground-up reimplementation in **Fortran 2018** of the Ecosystem Demography model
[ED2](https://github.com/EDmodel/ED2). It preserves ED2's size- and age-structured (demographic)
representation of terrestrial ecosystems — hydrology, land-surface biophysics, vegetation dynamics,
and soil biogeochemistry — while replacing the legacy code structure with modular, testable,
standards-conformant Fortran.

## Status

**Standalone demographic core implemented.** MEDS currently simulates cohort & patch dynamics —
individual-tree growth, mortality, recruitment, and cohort/patch fusion/fission — driven by
demographic rates supplied from *outside* the engine. There is deliberately no radiative transfer,
biophysics, or carbon balance yet: the rates are supplied by a swappable provider, and the only
provider today is a set of **test empirical relationships** (growth = size × competition; mortality =
f(growth); a simple recruitment schedule). A future mechanistic provider plugs into the same
interface with no engine change.

Highlights:
- Runs at a user-defined timestep (default **daily**, optionally **monthly**).
- Fusion/fission use **diameter & size-distribution** thresholds (no LAI), conserving basal area
  (cohorts) and site-level plant number (patches).
- Parallel by construction: hot kernels are standard-Fortran `do concurrent`, which **nvfortran**
  targets to **multicore CPU** (`-stdpar=multicore`) or **GPU** (`-stdpar=gpu`); CPU and GPU results
  match bit-for-bit.
- Builds clean and passes its CTest suite under **ifx** and **nvfortran**.

## Design goals

- **Modern Fortran 2018** — modules, derived-type encapsulation, explicit interfaces, `do concurrent`,
  `allocatable` ownership, parameterized real kinds, `pure`/`elemental` kernels.
- **Modular** — one responsibility per module, no hidden global mutable state, rates and parameters
  passed explicitly (the abstract `rate_provider_t` is the only rate seam).
- **Testable** — unit tests for conservation, container integrity, rate math, and full spin-up.
- **CMake-based build** — automatic Fortran module-dependency resolution; no hand-maintained object
  lists or repeated-build hacks.

## Building

Requires a Fortran 2018 compiler and CMake ≥ 3.20. Compilers may need activation first (Intel:
`source /opt/intel/oneapi/setvars.sh`; NVIDIA: put the HPC SDK `compilers/bin` on `PATH`).

```bash
# CPU, strict checks (Intel ifx) + run the tests + spin-up demo
cmake -S . -B build -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/meds_demo            # 60-year demographic spin-up (optional: meds_demo <years>)

# Multicore / GPU via standard do concurrent (NVIDIA nvfortran)
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_STDPAR=gpu
cmake --build build-gpu -j   # use -DMEDS_STDPAR=multicore for CPU threads
```

### Installing a compiler

On Debian/Ubuntu/WSL, helper scripts check whether a compiler is already available and offer to
install it. Each prompts before making changes (`-y` to skip the prompt, `-h` for help):

```bash
./scripts/install_gfortran.sh          # GNU gfortran (ED2's reference toolchain)
./scripts/install_gfortran.sh --hdf5   # gfortran + HDF5 dev files (libhdf5-dev)
./scripts/install_ifx.sh               # Intel ifx via the oneAPI APT repository
```

## Scientific reference

- Moorcroft, Hurtt & Pacala (2001), *Ecol. Monogr.* — the original ED formulation.
- Medvigy et al. (2009), *JGR Biogeosciences* — ED2.
- Longo et al. (2019), *Geosci. Model Dev.* 12:4309 — the ED-2.2 technical description.

See [CLAUDE.md](CLAUDE.md) for the architecture overview and contribution conventions.
