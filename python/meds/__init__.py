"""MEDS — Modular Ecosystem Dynamics Simulator (Python interface).

Four sub-modules, all backed by ONE compiled library (``libmeds.so``, see ``meds._libmeds``):

    meds.plant       — plant ecophysiology kernels: leaf gas exchange (``meds.plant.leaf``) and
                       the leaf-phenology signal kernel (``meds.plant.pheno``). Stateless.
    meds.demography  — the demographic carbon slow loop alone, behind an opaque site handle:
                       no biophysics, no forcing, no output streams.
    meds.model       — the FULL coupled model, exactly what the ``meds_main`` executable runs:
                       fast biophysics + slow demography, live met forcing, netCDF output and
                       both conservation ledgers, with the time loop handed to Python.

``meds.model`` is the whole model; the other two are pieces of it exposed for their own sake.

Importing `meds` is cheap and does NOT load any compiled library; each sub-module loads the
shared library lazily on first import, so `import meds` works even without it built.
"""
__version__ = "0.1.0"
__all__ = ["plant", "demography", "model"]
