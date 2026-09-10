"""meds.model — drive the FULL coupled MEDS model from Python.

The other sub-packages expose pieces: ``meds.plant.leaf`` and ``meds.plant.pheno`` are
stateless kernels, and ``meds.demography`` drives the slow loop alone. This one drives
what the ``meds_main`` executable drives — coupled fast biophysics + slow demography,
live met forcing, netCDF output, both conservation ledgers — with the time loop in
Python:

    from meds.model import Run

    with Run("meds_config_july.toml") as run:
        for step in run:
            if step.is_new_year:
                print(step.date, run.total_agb, run.soil_carbon)

``Run.step`` calls the same ``driver_step`` the executable calls, so a Python-driven run
and a ``meds_main`` run of the same config produce byte-identical output. Importing this
sub-module loads the shared library; ``import meds`` alone does not.
"""
from ._run import Run, StepInfo

__all__ = ["Run", "StepInfo"]
