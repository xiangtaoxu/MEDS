# Building MEDS

**Requirements:** a Fortran 2018 compiler, CMake ≥ 3.20, and the netCDF **C** library.

`meds_main` is the single entry point: it reads a config, runs the simulation, writes netCDF output,
and exits.

## Quick start

```bash
./scripts/install_netcdf.sh                    # if you do not have netCDF; prints the prefix to use

cmake -S . -B build -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build -j
ctest --test-dir build --output-on-failure

./build/meds_main examples/example_demography/example_config_main.toml
```

If the run cannot find `libnetcdf` or the compiler runtime, put both on the loader path:
`LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`. Note the `:$LD_LIBRARY_PATH` — replacing the
variable rather than extending it drops the Fortran runtime and every test fails to load.

## netCDF is a hard dependency

There is no netCDF-free build. Every configure needs `-DCMAKE_PREFIX_PATH=<prefix>` pointing at the
netCDF-C install, and because the netCDF CMake config pulls in HDF5, the build enables the C
language too.

**netCDF-Fortran is deliberately not used.** MEDS writes netCDF through the C API via
`iso_c_binding`, because netCDF-Fortran's compiled module format is gfortran-only and would not
load under ifx or nvfortran. The C ABI is compiler-agnostic, so the same `libnetcdf` links
identically under all three compilers.

## Compilers

| Compiler | Role | Activation |
|---|---|---|
| Intel `ifx` | The everyday compiler. Strict standards checking; runs the full suite. | `source /opt/intel/oneapi/setvars.sh` |
| NVIDIA `nvfortran` | The second back end, and the host-multicore path. (GPU offload builds but is slower than the CPU — see below.) | put the HPC SDK `compilers/bin` on `PATH` |
| GNU `gfortran` | Supported; ED2's reference toolchain. | usually already on `PATH` |

**A green ifx run is not sufficient.** Build the nvfortran multicore back end on new modules too.
Three portability traps have each bitten once, and each was invisible to ifx:

- **Never pass an array-valued function result straight into a call.** nvfortran's whole-program
  optimizer miscompiles the temporary descriptor — silently wrong values at `-O2`, a segfault at
  `-O0` — while `ifx -stand f18 -check all` emits only a remark. Bind to a named array first.
  (Issue #7.)
- **nvfortran rejects a `BLOCK` construct** anywhere inside a parallel region.
- **ifx builds `private` and `firstprivate` copies of a derived type through a compiler-generated
  static mold** that every thread writes. This is why the fast loop's per-patch scratch is an
  explicit per-thread pool indexed by thread number rather than an OpenMP data-sharing clause.

## Build types

```bash
# Release -- production runs and the netCDF layer.
cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$CONDA_PREFIX

# Debug -- strict checking (-stand f18 -check all -fpe0). Use it for engine work.
cmake -S . -B build-debug -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
```

Debug emits one harmless `arg_temp_created` remark per netCDF call, so prefer Release when the
output layer is what you are exercising. Per-compiler flags live in the `meds_fortran_flags()`
function in `CMakeLists.txt`.

## Parallel builds

**Host threading over the patch axis** — opt in at both build and run time:

```bash
cmake -S . -B build-omp -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DMEDS_OPENMP=ON -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
# then set [run].n_threads in the TOML (default 1).
```

Output is **byte-identical at any thread count**, and the test suite asserts it.

`-DMEDS_OPENMP=ON` does two things, and the second is load-bearing: it puts the OpenMP flag on the
fast-loop target, **and it adds the per-compiler "all locals on the stack" flag** (`-auto`,
`-frecursive`, `-Mrecursive`) to *every* target. Intel Fortran defaults to `-auto-scalar`, which
places local arrays and derived types in static storage shared by every thread; without that flag
the kernels race and return plausible, silently thread-count-dependent numbers.

**OpenMP `target` offload** (NVHPC only):

```bash
cmake -S . -B build-mc  -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release \
      -DMEDS_GPU=multicore -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release \
      -DMEDS_GPU=gpu -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
```

`MEDS_GPU` (`none|multicore|gpu`) affects NVHPC builds only; ifx and gfortran ignore it, because
the `!$omp` lines are comments without an OpenMP flag.

**On the GPU, honestly:** offload was measured in 2026-08 and is **not viable for the current
model**. The GPU build ran 1.4× slower than the CPU (49.8 s against 36.2 s), with one kernel at
0.4 % occupancy. `MEDS_GPU=multicore` is useful; `MEDS_GPU=gpu` builds and runs correctly and is
kept so the offload path does not rot, but it is not a speedup. The measurement and what to do
instead are in [`dev_plans/MEDS_GPU_EVALUATION.md`](dev_plans/MEDS_GPU_EVALUATION.md).

**Do not use `-stdpar=gpu`.** It forces the global CUDA managed allocator, whose deep copy and
finalize of the allocatable-component site type double-frees on the host. OpenMP `target` with
`-gpu=mem:separate` keeps all state in normal host memory and moves only the mapped arrays.

## Installing dependencies

Helper scripts check whether a dependency is already present and prompt before changing anything
(`-y` to skip the prompt, `-h` for help). They target Debian, Ubuntu and WSL.

```bash
./scripts/install_netcdf.sh            # netCDF C library, via conda-forge or apt
./scripts/install_ifx.sh               # Intel ifx via the oneAPI APT repository
./scripts/install_gfortran.sh          # GNU gfortran
./scripts/install_gfortran.sh --hdf5   # ... plus HDF5 dev files
```

`install_netcdf.sh` detects an existing netCDF via `nc-config` and prints the
`-DCMAKE_PREFIX_PATH` to pass. conda-forge is recommended because it ships the CMake config.

## Python and post-processing

The post-processing scripts need Python with **numpy**, **matplotlib**, **netCDF4** and **pillow**.
The conda environment in [`../environment.yml`](../environment.yml) provides them:

```bash
mamba env create -f environment.yml    # or: conda env create -f environment.yml
conda activate meds
```

To drive the model from Python, `pip install python/` compiles the Fortran and bundles
`libmeds.so` beside the package. See [`../python/README.md`](../python/README.md).

## Running the tests

```bash
ctest --test-dir build-ifx --output-on-failure          # all of them
ctest --test-dir build-debug -R fusion_cohort --output-on-failure   # one, by regex
```

The suite covers allometry round-trips, carbon and area conservation, container integrity, the
individual physical kernels, the coupled column on both integrators, the output round-trip, the
conservation ledgers, and a full spin-up. Every C-API shim is compiled by a mandatory test target,
so an ABI change is a build failure in a default build rather than a silent break in an optional
one.
