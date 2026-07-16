#!/usr/bin/env python3
"""example_demography — drive the MEDS carbon slow loop from Python via libmeds_c.

Loads the shipped config, builds a bare-ground site, and steps the demographic
carbon loop (the Fortran ``vegetation_dynamics`` orchestration), printing the
annual stand trajectory read back through the opaque-handle C-API. Because
``advance_slow`` runs the SAME Fortran engine as the standalone executable, this
reproduces the ``meds_main`` carbon run for the same config — the Python package
is a companion, not a reimplementation.

Run from the MEDS repo root (the config's PFT path is relative):

    MEDS_LIB=build-pylib/libmeds_c.so \\
    LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib:$LD_LIBRARY_PATH \\
    PYTHONPATH=python python3 examples/example_demography/run_carbon.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _cadence import slow_steps                                  # noqa: E402

from meds.demography import Config, Site                         # noqa: E402

CONFIG = sys.argv[1] if len(sys.argv) > 1 else "meds_config_main.toml"
N_YEARS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
N_PATCH = 4


def main():
    cfg = Config(CONFIG)
    print(f"# MEDS carbon spin-up via libmeds_c  (dt={cfg.dt_years*365:.3g} d, "
          f"{N_YEARS} yr, {N_PATCH} patches)")
    print(f"{'year':>4} {'gen':>6} {'n_cohort':>8} {'total_agb':>14} "
          f"{'total_lai':>12} {'total_nplant':>13}")

    with Site(cfg, n_patch=N_PATCH) as site:
        year = 0
        for istep, new_month, new_year in slow_steps(cfg.dt_years, N_YEARS):
            site.advance_slow(new_month, new_year)
            if new_year:
                year += 1
                print(f"{year:>4} {site.generation:>6} {site.n_cohort:>8} "
                      f"{site.total_agb:>14.6e} {site.total_lai:>12.4e} "
                      f"{site.total_nplant:>13.4e}")

        dbh = site.get("dbh")
        rng = f"[{dbh.min():.3g}, {dbh.max():.3g}]" if dbh.size else "[]"
        print(f"# final stand: {site.n_cohort} cohorts, dbh range {rng} cm, "
              f"generation {site.generation}")


if __name__ == "__main__":
    main()
