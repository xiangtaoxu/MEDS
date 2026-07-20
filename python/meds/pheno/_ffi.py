"""Internal ctypes bridge to libmeds_plant_c (phenology). NOT part of the public API.

The struct field order below MUST match the bind(c) mirror types in src/plant/meds_pheno_capi.f90
exactly. The shared library is the same libmeds_plant_c that meds.leaf uses (both *_capi.f90 shims
are compiled into it); see meds.leaf._ffi for the one-time CMake build.
"""
import os
import ctypes
from ctypes import c_double, c_int, byref
from pathlib import Path

#----- Field orders — must mirror meds_pheno_capi.f90. --------------------------------------#
_ENV_FIELDS = [
    ("temp_day", c_double), ("soil_temp", c_double), ("avail_water", c_double),
    ("dmax_leaf_psi", c_double), ("rad", c_double), ("daylength", c_double),
    ("doy", c_int), ("hemis_north", c_int),
]
_PARAM_FIELDS = [
    ("flush_cue_mask", c_int), ("shed_cue_mask", c_int),
    ("cue_sharpness", c_double), ("k_flush_max", c_double), ("k_shed_max", c_double),
    ("tau_flush", c_double), ("tau_shed", c_double),
    ("gdd_base_temp", c_double), ("chill_base_temp", c_double),
    ("phen_a", c_double), ("phen_b", c_double), ("phen_c", c_double),
    ("cold_drop_daylength", c_double), ("cold_drop_soiltemp1", c_double),
    ("cold_drop_soiltemp2", c_double),
    ("water_use_potential", c_int),
    ("water_off_threshold", c_double), ("water_on_threshold", c_double),
    ("water_window", c_double), ("water_width", c_double),
    ("leaf_psi_tlp", c_double), ("low_psi_threshold", c_double), ("high_psi_threshold", c_double),
    ("photo_crit", c_double), ("photo_slope", c_double),
    ("light_on_threshold", c_double), ("light_width", c_double), ("light_window", c_double),
]
_STATE_FIELDS = [
    ("flush_drive", c_double), ("shed_drive", c_double), ("gdd", c_double), ("chill", c_double),
    ("water_avg", c_double), ("low_psi_days", c_double), ("high_psi_days", c_double),
    ("light_avg", c_double),
]
_OUT_FIELDS = [
    ("leaf_flush_rate", c_double), ("leaf_shed_rate", c_double), ("cue_limiting", c_int),
]

# Field names that are integers on the C side (Python may pass bool/int; we coerce to int).
_INT_NAMES = {"doy", "hemis_north", "flush_cue_mask", "shed_cue_mask", "water_use_potential",
              "cue_limiting"}


class _EnvC(ctypes.Structure):
    _fields_ = _ENV_FIELDS


class _ParamsC(ctypes.Structure):
    _fields_ = _PARAM_FIELDS


class _StateC(ctypes.Structure):
    _fields_ = _STATE_FIELDS


class _OutC(ctypes.Structure):
    _fields_ = _OUT_FIELDS


ENV_NAMES = tuple(n for n, _ in _ENV_FIELDS)
PARAM_NAMES = tuple(n for n, _ in _PARAM_FIELDS)
STATE_NAMES = tuple(n for n, _ in _STATE_FIELDS)
OUT_NAMES = tuple(n for n, _ in _OUT_FIELDS)


def _find_lib():
    """Locate libmeds_plant_c.so; raise a helpful error if the CMake build hasn't run."""
    override = os.environ.get("MEDS_PLANT_LIB")
    if override:
        return override
    here = Path(__file__).resolve()
    candidates = [here.parent / "libmeds_plant_c.so"]           # bundled beside the package (future wheel)
    if len(here.parents) > 3:                                   # editable install: repo root is parents[3]
        repo_root = here.parents[3]                             # .../python/meds/pheno/_ffi.py -> repo root
        for build_dir in ("build-pylib", "build-py", "build"):
            candidates.append(repo_root / build_dir / "libmeds_plant_c.so")
    for cand in candidates:
        if cand.exists():
            return str(cand)
    raise FileNotFoundError(
        "libmeds_plant_c.so not found. Build it with:\n"
        "  cmake -S . -B build-pylib -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON "
        "-DCMAKE_PREFIX_PATH=<netcdf-prefix>\n"
        "  cmake --build build-pylib --target meds_plant_c\n"
        "then put the Intel/gfortran runtime on LD_LIBRARY_PATH, or set MEDS_PLANT_LIB to the .so path.\n"
        f"Looked in: {[str(c) for c in candidates]}")


_LIB = None


def _lib():
    """Load libmeds_plant_c once and cache it (with argtypes set)."""
    global _LIB
    if _LIB is None:
        lib = ctypes.CDLL(_find_lib())
        lib.meds_phenology_step.restype = None
        lib.meds_phenology_step.argtypes = [
            ctypes.POINTER(_EnvC), ctypes.POINTER(_ParamsC), c_double,
            ctypes.POINTER(_StateC), ctypes.POINTER(_OutC),
        ]
        _LIB = lib
    return _LIB


def _coerce(name, value):
    return int(value) if name in _INT_NAMES else float(value)


def _make(struct_cls, names, src):
    return struct_cls(**{n: _coerce(n, src[n]) for n in names})


def step(env, params, state, dt=1.0):
    """Advance the phenology kernel ONE step.

    `env`, `params`, `state` are plain dicts with the field names in _ENV_FIELDS / _PARAM_FIELDS /
    _STATE_FIELDS. Returns (out_dict, new_state_dict): out has leaf_flush_rate / leaf_shed_rate /
    cue_limiting; new_state is the advanced accumulators (feed it back next call).
    """
    env_c = _make(_EnvC, ENV_NAMES, env)
    params_c = _make(_ParamsC, PARAM_NAMES, params)
    state_c = _make(_StateC, STATE_NAMES, state)
    out_c = _OutC()
    _lib().meds_phenology_step(byref(env_c), byref(params_c), float(dt), byref(state_c), byref(out_c))
    out = {"leaf_flush_rate": out_c.leaf_flush_rate,
           "leaf_shed_rate": out_c.leaf_shed_rate,
           "cue_limiting": int(out_c.cue_limiting)}
    new_state = {n: getattr(state_c, n) for n in STATE_NAMES}
    return out, new_state
