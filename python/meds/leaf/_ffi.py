"""Internal ctypes bridge to libmeds_plant_c. NOT part of the public API.

This is the only module that touches ctypes; everything user-facing lives in `meds.leaf`
(dataclasses + enums). The struct field order below MUST match the bind(c) mirror types in
src/plant/meds_plant_capi.f90 exactly.

The shared library is built by CMake (it is NOT bundled in this pure-Python dev package yet):

    cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON -DMEDS_ENABLE_IO=OFF
    cmake --build build-py --target meds_plant_c

At import time the library is located from (in order): the MEDS_PLANT_LIB env var, a copy sitting
next to this package (a future bundled wheel), then a CMake build dir in the source tree.
"""
import os
import ctypes
from ctypes import c_double, c_int, byref, POINTER
from pathlib import Path

#----- Field orders — must mirror meds_plant_capi.f90. --------------------------------------#
_ENV_FIELDS = ("par", "leaf_temp", "vpd", "ca", "pressure", "psi_leaf", "gb")
_FLUX_REALS = ("A_net", "A_gross", "gs", "ci", "cs", "transpiration", "rd")
PARAM_FIELDS = (
    "vcmax25", "jmax25", "tpu25", "rd25", "kp25",
    "g0", "g1", "d0", "quantum_yield", "theta_j", "theta_cj", "theta_ic",
    "lambda25", "psi_open", "psi_close", "lambda_psi_exp",
    "kc25", "ko25", "gstar25",
    "ea_kc", "ea_ko", "ea_gstar", "ea_vcmax", "ea_jmax", "ea_rd",
    "hd_vcmax", "hd_jmax", "hd_rd", "ds_vcmax", "ds_jmax", "ds_rd",
    "o2_mol_frac", "absorptance", "phi_psii",
)


class _EnvC(ctypes.Structure):
    _fields_ = [(n, c_double) for n in _ENV_FIELDS]


class _FluxC(ctypes.Structure):
    _fields_ = ([(n, c_double) for n in _FLUX_REALS]
                + [("limitation", c_int), ("converged", c_int)])


class _ParamsC(ctypes.Structure):
    _fields_ = [("pathway", c_int)] + [(n, c_double) for n in PARAM_FIELDS]


_C3_DEMAND_FIELDS = ("A_gross", "Ac", "Aj", "Ap")


class _C3DemandC(ctypes.Structure):
    _fields_ = [(n, c_double) for n in _C3_DEMAND_FIELDS]


def _find_lib():
    """Locate libmeds_plant_c.so; raise a helpful error if the CMake build hasn't run."""
    override = os.environ.get("MEDS_PLANT_LIB")
    if override:
        return override
    here = Path(__file__).resolve()
    candidates = [here.parent / "libmeds_plant_c.so"]           # bundled beside the package (future wheel)
    if len(here.parents) > 3:                                  # editable install: repo root is parents[3]
        repo_root = here.parents[3]                            # .../python/meds/leaf/_ffi.py -> repo root
        for build_dir in ("build-py", "build-pylib", "build"):
            candidates.append(repo_root / build_dir / "libmeds_plant_c.so")
    for cand in candidates:
        if cand.exists():
            return str(cand)
    raise FileNotFoundError(
        "libmeds_plant_c.so not found. Build it with:\n"
        "  cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON -DMEDS_ENABLE_IO=OFF\n"
        "  cmake --build build-py --target meds_plant_c\n"
        "then put the Intel/gfortran runtime on LD_LIBRARY_PATH, or set MEDS_PLANT_LIB to the .so path.\n"
        f"Looked in: {[str(c) for c in candidates]}")


_LIB = None


def _lib():
    """Load libmeds_plant_c once and cache it (with argtypes/restypes set)."""
    global _LIB
    if _LIB is None:
        lib = ctypes.CDLL(_find_lib())
        lib.meds_leaf_solve.restype = None
        lib.meds_leaf_solve.argtypes = [POINTER(_EnvC), POINTER(_ParamsC),
                                        c_int, c_int, c_int, c_int, POINTER(_FluxC)]
        lib.meds_assimilation_demand_c3.restype = None
        lib.meds_assimilation_demand_c3.argtypes = [c_double] * 8 + [c_int, c_double, POINTER(_C3DemandC)]
        lib.meds_electron_transport_j.restype = c_double
        lib.meds_electron_transport_j.argtypes = [c_double] * 5
        lib.meds_peaked_arrhenius.restype = c_double
        lib.meds_peaked_arrhenius.argtypes = [c_double] * 5
        lib.meds_arrhenius.restype = c_double
        lib.meds_arrhenius.argtypes = [c_double] * 3
        _LIB = lib
    return _LIB


def solve(env, params, stomata, temp_response, colimitation, boundary_layer):
    """Call the coupled A-gs-Ci solver. `env`/`params` are plain dicts; returns a plain dict."""
    env_c = _EnvC(**{k: float(env[k]) for k in _ENV_FIELDS})
    params_c = _ParamsC(pathway=int(params["pathway"]),
                        **{k: float(params[k]) for k in PARAM_FIELDS})
    flux_c = _FluxC()
    _lib().meds_leaf_solve(byref(env_c), byref(params_c),
                           int(stomata), int(temp_response), int(colimitation),
                           1 if boundary_layer else 0, byref(flux_c))
    out = {n: getattr(flux_c, n) for n in _FLUX_REALS}
    out["limitation"] = flux_c.limitation
    out["converged"] = bool(flux_c.converged)
    return out


def assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colimitation, theta):
    """Raw C3 FvCB demand at a prescribed Ci (mole-fraction kinetics, no T-scaling). Plain dict."""
    dem_c = _C3DemandC()
    _lib().meds_assimilation_demand_c3(float(ci), float(vcmax), float(j), float(tpu), float(gstar),
                                float(kc), float(ko), float(o2), int(colimitation), float(theta),
                                byref(dem_c))
    return {n: getattr(dem_c, n) for n in _C3_DEMAND_FIELDS}


def electron_transport_j(par, absorptance, phi_psii, jmax, theta):
    """Electron-transport rate J from Jmax and PAR (non-rectangular hyperbola)."""
    return _lib().meds_electron_transport_j(par, absorptance, phi_psii, jmax, theta)


def peaked(k25, ea, hd, ds, t_leaf):
    return _lib().meds_peaked_arrhenius(k25, ea, hd, ds, t_leaf)


def arrhenius(k25, ea, t_leaf):
    return _lib().meds_arrhenius(k25, ea, t_leaf)
