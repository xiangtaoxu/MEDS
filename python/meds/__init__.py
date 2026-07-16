"""MEDS — Modular Ecosystem Dynamics Simulator (Python interface).

The eventual home of a pip/conda-installable Python front end to MEDS. Today it exposes the first
submodule:

    meds.leaf        — leaf-level photosynthesis + stomatal conductance (FvCB C3 / Collatz C4).
    meds.demography  — drive the demographic carbon slow loop via libmeds_c (opaque site handle).

Importing `meds` is cheap and does NOT load any compiled library; each submodule loads its shared
library lazily (meds.leaf -> libmeds_plant_c, meds.demography -> libmeds_c) on first import, so
`import meds` works even without the shared libraries built.
"""
__version__ = "0.1.0"
__all__ = ["leaf", "demography"]
