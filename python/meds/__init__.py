"""MEDS — Modular Ecosystem Dynamics Simulator (Python interface).

The eventual home of a pip/conda-installable Python front end to MEDS. Today it exposes the first
submodule:

    meds.leaf   — leaf-level photosynthesis + stomatal conductance (FvCB C3 / Collatz C4).

Future submodules (meds.demography, ...) will attach here as their Fortran C-APIs land. Importing
`meds` is cheap and does NOT load any compiled library; the leaf submodule loads libmeds_leaf_c
lazily on first use, so `import meds` works even without the shared library built.
"""
__version__ = "0.1.0"
__all__ = ["leaf"]
