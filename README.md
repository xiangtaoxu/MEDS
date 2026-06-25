# MEDS — Modular Ecosystem Dynamics Simulator

**MEDS** stands for **Modular Ecosystem Dynamics Simulator** — and can equally be read as **Modern
Ecosystem Demography Simulator**, a nod to its descent from the ED / ED2 lineage.

MEDS is a ground-up reimplementation in **Fortran 2018** of the Ecosystem Demography model
[ED2](https://github.com/EDmodel/ED2). It preserves ED2's size- and age-structured (demographic)
representation of terrestrial ecosystems — hydrology, land-surface biophysics, vegetation dynamics,
and soil biogeochemistry — while replacing the legacy code structure with modular, testable,
standards-conformant Fortran.

## Status

**Standalone demographic core implemented.** MEDS currently simulates cohort & patch dynamics —
individual-tree growth, mortality, recruitment, cohort/patch fusion/fission, and treefall **patch
disturbance** — driven by demographic rates supplied from *outside* the engine as three plain arrays.
Size follows the **pan-tropical (ED2 `iallom==3`) allometry**, and each cohort carries **aboveground
biomass (carbon)** and **leaf area**. There is deliberately no radiative transfer or full
biogeochemistry yet: the rates come from a set of **test empirical relationships** — light competition
through overtopping LAI (growth = max-relative-growth × light × dbh; mortality = baseline +
shade-driven). A future mechanistic module produces the same arrays with no engine change.

Highlights:
- Runs at a user-defined timestep (default **daily**, optionally **weekly/monthly**).
- **Wood density** is the single PFT trade-off axis: low density ⇒ fast growth but high mortality
  (especially under shade), high density ⇒ slow but tolerant — the classic ED growth–survival trade-off.
- Cohort fusion/fission key on **height & LAI**, conserving **total aboveground biomass (carbon)**;
  patch fusion compares an ED2-style **cumulative-LAI light profile**; treefall disturbance opens
  **age-0 gaps**, giving the site a successional patch age structure.
- Parallel by construction: the hot kernels carry explicit **OpenMP `target`** regions over plain
  arrays, which **nvfortran** offloads to **multicore CPU** (`-DMEDS_GPU=multicore` → `-mp`) or the
  **GPU** (`-DMEDS_GPU=gpu` → `-mp=gpu`); CPU and GPU results match.
- Optional **netCDF output** (`-DMEDS_ENABLE_IO=ON`) writes the full cohort/patch/site state over time
  via the netCDF C library — no netCDF-Fortran dependency, so it works under both ifx and nvfortran.
- Builds clean and passes its CTest suite under **ifx** and **nvfortran**.

## Design goals

- **Modern Fortran 2018** — modules, derived-type encapsulation, explicit interfaces, OpenMP-target
  array kernels, `allocatable` ownership, parameterized real kinds, `pure`/`elemental` helpers.
- **Modular** — one responsibility per module, no hidden global mutable state, rates and parameters
  passed explicitly as data (the `update_demography` array interface is the only rate seam).
- **Testable** — unit tests for allometry round-trips, carbon/area conservation, container integrity,
  rate math, disturbance, and full spin-up.
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

# Multicore / GPU via OpenMP target (NVIDIA nvfortran)
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=gpu
cmake --build build-gpu -j   # use -DMEDS_GPU=multicore for CPU threads
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
