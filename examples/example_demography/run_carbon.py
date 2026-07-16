#!/usr/bin/env python3
"""example_demography — drive the MEDS carbon slow loop from Python via libmeds_c.

Loads the shipped config, builds a bare-ground site, and steps the demographic
carbon loop (the Fortran ``vegetation_dynamics`` orchestration) with the standard
daily/monthly/annual cadence, printing the annual stand trajectory read back
through the opaque-handle C-API. Because ``advance_slow`` runs the SAME Fortran
engine as the standalone executable, this reproduces the ``meds_main`` carbon run
for the same config — the Python package is a companion, not a reimplementation.

Run from the MEDS repo root (the config's PFT path is relative):

    MEDS_LIB=build-pylib/libmeds_c.so \\
    LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib:$LD_LIBRARY_PATH \\
    PYTHONPATH=python python3 examples/example_demography/run_carbon.py
"""
import sys

from meds.demography import Config, Site

CONFIG = sys.argv[1] if len(sys.argv) > 1 else "meds_config_main.toml"
N_YEARS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
N_PATCH = 4
NDAY = 365


def main():
    cfg = Config(CONFIG)
    step_days = max(1, round(cfg.dt_years * NDAY))
    nsteps = (N_YEARS * NDAY) // step_days

    print(f"# MEDS carbon spin-up via libmeds_c  (dt={cfg.dt_years*NDAY:.3g} d, "
          f"{N_YEARS} yr, {N_PATCH} patches)")
    print(f"{'year':>4} {'gen':>6} {'n_cohort':>8} {'total_agb':>14} "
          f"{'total_lai':>12} {'total_nplant':>13}")

    with Site(cfg, n_patch=N_PATCH) as site:
        yday = prev_month = year = 0
        for istep in range(1, nsteps + 1):
            yday += step_days
            new_year = yday > NDAY
            if new_year:
                yday -= NDAY
            month = min(12, (yday - 1) * 12 // NDAY + 1)
            new_month = (month != prev_month) or new_year or istep == 1
            if istep == 1:
                new_year = True
            prev_month = month

            site.advance_slow(new_month, new_year)

            if new_year:
                year += 1
                print(f"{year:>4} {site.generation:>6} {site.n_cohort:>8} "
                      f"{site.total_agb:>14.6e} {site.total_lai:>12.4e} "
                      f"{site.total_nplant:>13.4e}")

        # A generation-guarded snapshot of the final stand (global_id-keyed).
        snap = site.get("dbh")
        print(f"# final stand: {site.n_cohort} cohorts, "
              f"dbh range [{snap.min():.3g}, {snap.max():.3g}] cm, "
              f"generation {site.generation}")


if __name__ == "__main__":
    main()
