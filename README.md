# MEDS — Modern and modular Ecosystem Demography Simulator

MEDS is a ground-up reimplementation in **Fortran 2018** of the Ecosystem Demography model
[ED2](https://github.com/EDmodel/ED2). It preserves ED2's size- and age-structured (demographic)
representation of terrestrial ecosystems — hydrology, land-surface biophysics, vegetation dynamics,
and soil biogeochemistry — while replacing the legacy code structure with modular, testable,
standards-conformant Fortran.

## Status

Greenfield: documentation and project scaffolding only. No Fortran source or build tree yet.

## Design goals

- **Modern Fortran 2018** — modules + submodules, derived-type encapsulation, explicit interfaces,
  `allocatable` ownership, parameterized real kinds.
- **Modular** — one responsibility per module, no hidden global mutable state, parameters passed
  explicitly as configuration objects.
- **Testable** — unit tests for math kernels plus whole-model regression against ED2.
- **CMake-based build** — automatic Fortran module-dependency resolution; no hand-maintained object
  lists or repeated-build hacks.

## Building (target workflow)

Requires a Fortran 2018 compiler (e.g. `gfortran`), CMake ≥ 3.16, and HDF5/netCDF.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Scientific reference

- Moorcroft, Hurtt & Pacala (2001), *Ecol. Monogr.* — the original ED formulation.
- Medvigy et al. (2009), *JGR Biogeosciences* — ED2.
- Longo et al. (2019), *Geosci. Model Dev.* 12:4309 — the ED-2.2 technical description.

See [CLAUDE.md](CLAUDE.md) for the architecture overview and contribution conventions.
