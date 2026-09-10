"""High-level wrapper over the MEDS full-coupled-model C-API.

``Run`` is an opaque handle to a live simulation -- the same thing the ``meds_main``
executable holds -- with the time loop handed back to Python:

    with Run("meds_config_july.toml") as run:
        for step in run:
            if step.is_new_year:
                print(step.date, run.total_agb, run.soil_carbon)

``driver_step`` is the identical call the executable makes, so a Python-driven run
and a ``meds_main`` run of the same config produce byte-identical netCDF. Nothing is
re-implemented here; this is the executable's loop, in Python.

Array getters COPY out of the Fortran SoA into fresh numpy arrays (they never alias),
and the SoA is reordered by cohort/patch fusion and fission at the monthly cadence, so
a positional snapshot cached across a step is stale -- ``global_id`` is the only key
that survives a step.
"""
import ctypes
from datetime import date

import numpy as np

from ._ffi import lib

#----- Field ids match meds_capi_demography's, so the two sub-packages read the same way. -----#
_REAL = {"dbh": 0, "height": 1, "nplant": 2, "agb": 3, "leaf_area": 4,
         "overtopping_lai": 5, "growth_avg": 6, "wood_carbon": 7}
_INT = {"pft": 0, "owner_patch": 1, "global_id": 2}

#----- driver_step status codes (meds_driver's DRIVER_* parameters). --------------------------#
OK, FINISHED, ERR_NAN, ERR_AREA, ERR_SOILC = 0, 1, 2, 3, 4


class StepInfo:
    """What one slow step advanced to. Yielded by iterating a Run."""

    __slots__ = ("date", "istep", "iyear", "is_new_month", "is_new_year")

    def __init__(self, d, istep, iyear, is_new_month, is_new_year):
        self.date, self.istep, self.iyear = d, istep, iyear
        self.is_new_month, self.is_new_year = is_new_month, is_new_year

    def __repr__(self):
        return f"<StepInfo {self.date} step={self.istep}>"


class Run:
    """A live MEDS simulation driven from Python.

    Parameters
    ----------
    config : str or Path
        The run TOML -- the same file ``meds_main`` takes on its command line.
    verbose : bool
        Keep the Fortran progress lines (init/forcing/summary/ledger) on stdout.
        They go to the process's stdout, not through Python's ``sys.stdout``.
    """

    def __init__(self, config, verbose=True):
        b = str(config).encode("utf-8")
        self.handle = lib.meds_run_open(b, len(b), 1 if verbose else 0)
        if self.handle == -1:
            raise RuntimeError("MEDS run registry full (max 4 concurrent runs)")
        if self.handle < 0:
            raise RuntimeError(f"MEDS could not open a run from {config!r}")
        self._finalized = False

    #----- Lifecycle ---------------------------------------------------------------------#
    def step(self):
        """Advance ONE slow step. Returns a StepInfo, or None once the run is over."""
        if self.is_done:
            return None
        prev = self.date
        status = lib.meds_run_step(self.handle)
        if status == ERR_NAN:
            #----- The Fortran driver RETURNS these rather than `error stop`ping, which is
            #      the only reason a bad state is an exception here instead of a dead
            #      interpreter. See the meds_driver module header.
            raise RuntimeError(f"MEDS: NaN detected in state at {self.date}")
        if status == ERR_SOILC:
            #----- A CENTURY pool went negative or past the divergence ceiling. The Fortran
            #      side has already printed which pool and patch; the slow ledger's report
            #      names the PHASE that did it.
            raise RuntimeError(
                f"MEDS: physically impossible soil-carbon pool at {self.date} "
                "(see the printed pool/patch and the slow-ledger IMPOSSIBLE STORE line)")
        if status == FINISHED:
            return None
        if status < 0:
            raise RuntimeError(f"MEDS: invalid run handle {self.handle}")
        now = self.date
        return StepInfo(now, self.istep, self.iyear,
                        (now.year, now.month) != (prev.year, prev.month),
                        now.year != prev.year)

    def run_to_end(self, callback=None):
        """Step until end_time. ``callback(run, step)`` is called after each step."""
        n = 0
        for info in self:
            n += 1
            if callback is not None:
                callback(self, info)
        return n

    def finalize(self):
        """Terminal checkpoint, conservation reports, close the netCDF streams.

        The state stays readable afterwards -- freeing it is ``close()``.
        """
        if self._finalized:
            return 0
        status = lib.meds_run_finalize(self.handle)
        self._finalized = True
        if status == ERR_AREA:
            raise RuntimeError("MEDS: site area not conserved")
        return status

    def close(self):
        """Finalize if needed, then release the site's Fortran allocations."""
        if self.handle is None:
            return
        self.finalize()
        lib.meds_run_free(self.handle)
        self.handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def __iter__(self):
        while True:
            info = self.step()
            if info is None:
                return
            yield info

    #----- Where the calendar is ---------------------------------------------------------#
    @property
    def is_done(self):
        return bool(lib.meds_run_is_done(self.handle))

    @property
    def date(self):
        return date(lib.meds_run_year(self.handle),
                    lib.meds_run_month(self.handle),
                    lib.meds_run_day(self.handle))

    @property
    def istep(self):
        return lib.meds_run_istep(self.handle)

    @property
    def iyear(self):
        return lib.meds_run_iyear(self.handle)

    #----- Site aggregates: the same reducers the executable's summary line prints. -------#
    @property
    def n_patch(self):
        return lib.meds_run_n_patch(self.handle)

    @property
    def n_cohort(self):
        return lib.meds_run_n_cohort(self.handle)

    @property
    def total_agb(self):
        """Above-ground biomass [kgC/m2], area-weighted over patches."""
        return lib.meds_run_total_agb(self.handle)

    @property
    def total_lai(self):
        """Leaf area index [m2/m2]."""
        return lib.meds_run_total_lai(self.handle)

    @property
    def total_nplant(self):
        """Stem density [plants/m2]."""
        return lib.meds_run_total_nplant(self.handle)

    @property
    def total_basal_area(self):
        """Basal area [cm2/m2]."""
        return lib.meds_run_total_basal_area(self.handle)

    @property
    def soil_carbon(self):
        """All seven CENTURY pools [kgC/m2], area-weighted.

        Identically zero unless ``[soil_carbon].soil_carbon_on`` -- with the feature
        off, litter is DISCARDED rather than stored, so this stays flat at zero while
        the stand grows. That is worth watching from the loop.
        """
        return lib.meds_run_soil_carbon(self.handle)

    #----- Per-cohort state, copied out ---------------------------------------------------#
    def cohorts(self, *fields):
        """Copy per-cohort fields into numpy arrays: ``run.cohorts("dbh", "pft")``.

        With no arguments, returns every field. Valid names are the keys of
        ``Run.real_fields`` and ``Run.int_fields``.
        """
        if not fields:
            fields = tuple(_REAL) + tuple(_INT)
        n = self.n_cohort
        out = {}
        for f in fields:
            if f in _REAL:
                buf = np.zeros(max(n, 1), dtype=np.float64)
                lib.meds_run_get_real(self.handle, _REAL[f],
                                      buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            elif f in _INT:
                buf = np.zeros(max(n, 1), dtype=np.int32)
                lib.meds_run_get_int(self.handle, _INT[f],
                                     buf.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
            else:
                raise KeyError(f"unknown cohort field {f!r}; "
                               f"have {sorted(set(_REAL) | set(_INT))}")
            out[f] = buf[:n]
        return out

    real_fields = tuple(_REAL)
    int_fields = tuple(_INT)

    def __repr__(self):
        if self.handle is None:
            return "<Run closed>"
        return (f"<Run {self.date} step={self.istep} "
                f"cohorts={self.n_cohort} patches={self.n_patch}>")
