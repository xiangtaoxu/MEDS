"""High-level wrapper over the MEDS demography C-API.

``Config`` loads a run configuration from TOML; ``Site`` is an opaque handle to a
live cohort/patch state that Python drives through the Fortran carbon slow loop.
Array getters COPY out of the Fortran SoA into fresh numpy arrays (never alias),
and every getter is paired with the site's ``generation`` counter — bumped on
every ``advance_slow`` (which reorders the SoA on fuse/fission), so a snapshot
cached across a step is detectably stale; ``global_id`` is the only stable key.
"""
import ctypes

import numpy as np

from ._ffi import lib

_REAL = {"dbh": 0, "height": 1, "nplant": 2, "agb": 3, "leaf_area": 4,
         "overtopping_lai": 5, "growth_avg": 6, "wood_carbon": 7}
_INT = {"pft": 0, "owner_patch": 1, "global_id": 2}


class Config:
    def __init__(self, path):
        b = str(path).encode("utf-8")
        self.handle = lib.meds_config_load(b, len(b))
        if self.handle < 0:
            raise RuntimeError("MEDS config registry full")

    @property
    def dt_years(self):
        return lib.meds_config_dt_years(self.handle)


class Site:
    def __init__(self, cfg, n_patch=1):
        self.cfg = cfg
        self.handle = lib.meds_site_create()
        if self.handle < 0:
            raise RuntimeError("MEDS site registry full")
        lib.meds_site_init_bare(self.handle, cfg.handle, int(n_patch))

    def advance_slow(self, is_new_month, is_new_year):
        lib.meds_advance_slow(self.handle, self.cfg.handle,
                              int(bool(is_new_month)), int(bool(is_new_year)))

    # ---- scalar site diagnostics -------------------------------------------------
    @property
    def generation(self):
        return lib.meds_site_generation(self.handle)

    @property
    def n_cohort(self):
        return lib.meds_site_n_cohort(self.handle)

    @property
    def total_agb(self):
        return lib.meds_site_total_agb(self.handle)

    @property
    def total_lai(self):
        return lib.meds_site_total_lai(self.handle)

    @property
    def total_nplant(self):
        return lib.meds_site_total_nplant(self.handle)

    @property
    def total_basal_area(self):
        return lib.meds_site_total_basal_area(self.handle)

    # ---- per-cohort copy-out (fresh numpy arrays) --------------------------------
    def get(self, field):
        n = self.n_cohort
        if field in _REAL:
            buf = (ctypes.c_double * max(n, 1))()
            lib.meds_site_get_real(self.handle, _REAL[field], buf)
            return np.frombuffer(buf, dtype=np.float64, count=n).copy()
        if field in _INT:
            buf = (ctypes.c_int * max(n, 1))()
            lib.meds_site_get_int(self.handle, _INT[field], buf)
            return np.frombuffer(buf, dtype=np.int32, count=n).copy()
        raise KeyError(f"unknown field {field!r}; "
                       f"real={sorted(_REAL)} int={sorted(_INT)}")

    def snapshot(self):
        """A generation-stamped dict of the current stand (global_id-keyed rows)."""
        return {
            "generation": self.generation,
            "global_id": self.get("global_id"),
            "pft": self.get("pft"),
            "dbh": self.get("dbh"),
            "height": self.get("height"),
            "nplant": self.get("nplant"),
            "agb": self.get("agb"),
        }

    def free(self):
        lib.meds_site_free(self.handle)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.free()
