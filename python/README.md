# `meds` — Python interface to MEDS

The eventual pip/conda-installable Python front end to MEDS. Today it ships one submodule,
[`meds.leaf`](meds/leaf/__init__.py) — a clean, ctypes-free wrapper over the compiled Fortran
leaf-physiology model (FvCB C3 / Collatz C4, Leuning / Medlyn / Katul stomata, the coupled
A–gs–Ci solver). Future submodules (`meds.demography`, …) attach here as their Fortran C-APIs land.

```python
import meds.leaf as leaf
params = leaf.c3_params(vcmax25=60.0, jmax25=108.0)          # or leaf.c4_params(...)
flux = leaf.gas_exchange(par=1500.0, leaf_temp=298.15,       # leaf_temp in KELVIN
                         vpd=1000.0, ca=400.0, params=params,
                         stomata=leaf.Stomata.MEDLYN)
print(flux.a_net, flux.gs, flux.ci, flux.limitation, flux.converged)
```

## Dev install (rapid iteration)

The package is pure-Python for now: the compiled leaf library `libmeds_plant_c` is built **once** by
the top-level CMake and located at runtime, so an editable install means Python edits are live with
no rebuild — you only re-run `cmake --build` when the *Fortran* changes.

```bash
# 1. Build the shared library once (from the repo root; netCDF is required now -> pass its prefix):
cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DMEDS_BUILD_PYLIB=ON -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build-py --target meds_plant_c            # -> build-py/libmeds_plant_c.so

# 2. Editable install of the Python package (offline; this conda env's setuptools needs the flag):
SETUPTOOLS_USE_DISTUTILS=stdlib pip install -e python/ --no-build-isolation

# 3. Use it (the Fortran runtime must be on LD_LIBRARY_PATH):
source /opt/intel/oneapi/setvars.sh                   # ifx build -> Intel runtime
python -m meds.leaf                                    # round-trip self-test
pytest python/tests                                    # API smoke tests
```

`meds.leaf` finds the `.so` via (in order) the `MEDS_PLANT_LIB` env var, a copy beside the package
(a future bundled wheel), then a CMake build dir in the source tree — so the editable install just
works from `build-py/`.

## Road to a distributable package

This dev layout is deliberately the skeleton of the shipped package, so the next rungs are additive:

1. **Now — editable install** (this file): fastest edit/test loop; API still moving.
2. **Build-on-install** — switch `[build-system]` to `scikit-build-core`, which runs CMake during
   `pip install .` and bundles the `.so` into the wheel. Needs a Fortran compiler at install time.
3. **Portable wheels** — build `libmeds_plant_c` with **gfortran** (drops the Intel-runtime
   dependency; the top-level CMake already supports the GNU compiler) and add `cibuildwheel` +
   `auditwheel`, so end users get `pip install meds` with no compiler. Best for ecology users.
