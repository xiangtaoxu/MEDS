# `meds` — Python interface to MEDS

The Python front end to MEDS: clean, ctypes-free wrappers over the compiled Fortran model. **One
shared library, `libmeds.so`, backs every sub-module** (structure-plan decision #1) — per-subsystem
independence lives at the CALL level, not the link level, so there is one binary, one ABI, and one
place for it to break.

**`meds.plant`** — the plant-ecophysiology model, with two submodules:

- [`meds.plant.leaf`](meds/plant/leaf.py) — leaf gas exchange (FvCB C3 / Collatz C4, Leuning / Medlyn /
  Katul stomata, the coupled A–gs–Ci solver);
- [`meds.plant.pheno`](meds/plant/pheno.py) — leaf phenology (the two per-day flush / shed rate
  tendencies + the four strategy presets; see `examples/example_phenology/`).

**`meds.demography`** — [`Config` and `Site`](meds/demography/_site.py): load a TOML config, build a
site, and step the carbon slow loop or feed it externally computed rates (`apply_rates`). This is the
same Fortran engine `meds_main` runs, driven through an opaque handle.

Future submodules (`meds.plant.hydraulics`, `meds.fast`, …) attach as their Fortran C-APIs land.

```python
import meds.plant.leaf as leaf
params = leaf.c3_params(vcmax25=60.0, jmax25=108.0)          # or leaf.c4_params(...)
flux = leaf.gas_exchange(par=1500.0, leaf_temp=298.15,       # leaf_temp in KELVIN
                         vpd=1000.0, ca=400.0, params=params,
                         stomata=leaf.Stomata.MEDLYN)
print(flux.A_net, flux.gs, flux.ci, flux.limitation, flux.converged)

import meds.plant.pheno as pheno
ph = pheno.Phenology(pheno.temperate_deciduous())            # a stateful phenology driver
out = ph.step(temp_day=290.0, soil_temp=290.0, daylength=13.0, doy=150)
print(out.leaf_flush_rate, out.leaf_shed_rate)               # [1/day] flush / shed tendencies
```

## Install

`pip install python/` **compiles the Fortran and bundles the library**: scikit-build-core drives the
top-level CMake, and the wheel ships `meds/libmeds.so` beside the Python modules.

```bash
source /opt/intel/oneapi/setvars.sh                    # so CMake finds ifx
CMAKE_PREFIX_PATH=$CONDA_PREFIX pip install python/    # netCDF prefix; builds + bundles libmeds.so
python -m meds.plant                                    # round-trip self-test
```

**The installed library needs no `LD_LIBRARY_PATH`.** Its RPATH is written at build time from the
directories it actually linked — netCDF, HDF5 and the Fortran runtime — so it loads in a bare venv.
The corollary is that **the wheel is machine-local**: those paths are this machine's conda prefix and
oneAPI install. That is the right trade for `pip install` from a source checkout, which is how MEDS
is used; it is not a redistributable manylinux wheel (see the road map below).

### Dev loop (edit Python without recompiling)

```bash
# 1. Build the library once from the repo root:
cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DMEDS_BUILD_PYLIB=ON -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build-py --target meds_py                 # -> build-py/libmeds.so

# 2. Run against the source tree (the Fortran runtime must be on LD_LIBRARY_PATH here):
source /opt/intel/oneapi/setvars.sh
PYTHONPATH=python python -m meds.plant
PYTHONPATH=python pytest python/tests
```

`meds._libmeds` finds the library via (in order) the `MEDS_LIB` env var, a copy beside the package
(how a wheel ships it), then a CMake build dir in the source tree — so both paths above work with no
configuration. The retired `MEDS_PLANT_LIB` is recognised only to tell you it has been replaced.

## Road to a distributable package

This dev layout is deliberately the skeleton of the shipped package, so the next rungs are additive:

1. ~~**Editable install**~~ — was the starting point; the dev loop above replaces it.
2. **Build-on-install** — **DONE**: `[build-system]` is `scikit-build-core`, so `pip install python/`
   runs CMake and bundles the `.so`. Needs a Fortran compiler at install time.
3. **Portable wheels** — build with **gfortran** (drops the Intel-runtime dependency; the top-level
   CMake already supports GNU) and add `cibuildwheel` + `auditwheel`, so end users get
   `pip install meds` with no compiler. That is what turns today's machine-local wheel into a
   redistributable one, and it means vendoring netCDF's closure — the ~67 MB cost the structure plan
   measured (§7.3) and deferred. Best for ecology users.
