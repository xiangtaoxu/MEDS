"""ctypes backend for the MEDS demography C-API (libmeds_c.so).

Loads the shared library built with ``-DMEDS_BUILD_PYLIB=ON`` and declares the
argument/return signatures of every ``bind(c)`` entry point. Set ``MEDS_LIB`` to
point at the .so explicitly; otherwise the usual build dirs are searched.
"""
import ctypes
import glob
import os

c_int, c_long, c_double, c_char_p = (
    ctypes.c_int, ctypes.c_long, ctypes.c_double, ctypes.c_char_p,
)
_dptr = ctypes.POINTER(c_double)
_iptr = ctypes.POINTER(c_int)


def _find_lib():
    env = os.environ.get("MEDS_LIB")
    if env and os.path.exists(env):
        return env
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    hits = glob.glob(os.path.join(root, "build*", "libmeds_c.so"))
    if hits:
        return hits[0]
    raise OSError(
        "libmeds_c.so not found. Build with -DMEDS_BUILD_PYLIB=ON or set MEDS_LIB."
    )


lib = ctypes.CDLL(_find_lib())


def _sig(name, restype, argtypes):
    fn = getattr(lib, name)
    fn.restype = restype
    fn.argtypes = argtypes
    return fn


# lifecycle
_sig("meds_config_load",   c_int,    [c_char_p, c_int])
_sig("meds_config_dt_years", c_double, [c_int])
_sig("meds_config_n_pft",  c_int,    [c_int])
_sig("meds_site_create",   c_int,    [])
_sig("meds_site_init_bare", None,    [c_int, c_int, c_int])
_sig("meds_site_n_patch",  c_int,    [c_int])
_sig("meds_advance_slow",  None,     [c_int, c_int, c_int, c_int])
_sig("meds_apply_rates",   None,     [c_int, c_int, _dptr, _dptr, _dptr, c_int, c_int])
_sig("meds_site_free",     None,     [c_int])
# scalar reads
_sig("meds_site_generation",   c_long,   [c_int])
_sig("meds_site_n_cohort",     c_int,    [c_int])
_sig("meds_site_total_agb",    c_double, [c_int])
_sig("meds_site_total_lai",    c_double, [c_int])
_sig("meds_site_total_nplant", c_double, [c_int])
_sig("meds_site_total_basal_area", c_double, [c_int])
# array copy-out
_sig("meds_site_get_real", None, [c_int, c_int, _dptr])
_sig("meds_site_get_int",  None, [c_int, c_int, _iptr])
