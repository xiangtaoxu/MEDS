"""Locate and load ``libmeds.so`` — the ONE compiled backend behind every sub-package.

There used to be two shared libraries (``libmeds_plant_c`` for leaf + phenology, ``libmeds_c``
for the site driver) and two copies of this search logic, one per sub-package. Structure-plan
decision #1 merged them: per-subsystem independence lives at the CALL level now (``meds.plant.leaf``,
``meds.plant.pheno``, ``meds.demography`` are verbs over one backend), not at the link level. So the
search lives here once, and a sub-package that wants the library asks for it.

Search order:
  1. ``$MEDS_LIB`` — an explicit path, which always wins.
  2. Next to this package, i.e. ``meds/libmeds.so`` — how a built wheel ships it.
  3. A CMake build directory in the source tree, for an editable install.

The library needs the Fortran runtime on ``LD_LIBRARY_PATH`` (``source /opt/intel/oneapi/setvars.sh``
for ifx) and, because it links the model's netCDF I/O, a netCDF runtime — which is why ``netcdf4``
is a hard dependency of the wheel rather than an extra.
"""
import ctypes
import os
from pathlib import Path

_SONAME = "libmeds.so"
#----- Build dirs an editable install might be paired with, most specific first. ----------------#
_BUILD_DIRS = ("build-py", "build-pylib", "build", "build-ifx")
#----- Retired env vars, still recognised so an old command line says what to do rather than -----#
#      failing with a path that no longer exists.                                                 #
_LEGACY_ENV = ("MEDS_PLANT_LIB",)

_LIB = None


def _candidates():
    here = Path(__file__).resolve()                  # .../python/meds/_libmeds.py
    cands = [here.parent / _SONAME]                  # bundled beside the package (wheel)
    if len(here.parents) > 2:
        repo_root = here.parents[2]                  # .../python/meds/ -> repo root
        cands += [repo_root / b / _SONAME for b in _BUILD_DIRS]
    return cands


def library_path():
    """Absolute path to libmeds.so, or a FileNotFoundError that says how to build it."""
    override = os.environ.get("MEDS_LIB")
    if override:
        if not os.path.exists(override):
            raise FileNotFoundError(f"MEDS_LIB points at a missing file: {override}")
        return override

    legacy = [v for v in _LEGACY_ENV if os.environ.get(v)]
    for cand in _candidates():
        if cand.exists():
            return str(cand)

    hint = ""
    if legacy:
        hint = (f"\n({', '.join(legacy)} is set, but the two shared libraries were merged into "
                f"one {_SONAME}; use MEDS_LIB instead.)")
    raise FileNotFoundError(
        f"{_SONAME} not found. Build it with:\n"
        "  cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON\n"
        "  cmake --build build-py --target meds_py\n"
        "then put the Fortran runtime on LD_LIBRARY_PATH, or set MEDS_LIB to the .so path."
        + hint
        + f"\nLooked in: {[str(c) for c in _candidates()]}")


def lib():
    """The loaded CDLL, cached. Every sub-package's _ffi goes through this."""
    global _LIB
    if _LIB is None:
        _LIB = ctypes.CDLL(library_path())
    return _LIB
