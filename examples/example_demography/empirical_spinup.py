#!/usr/bin/env python3
"""example_demography — the EMPIRICAL demography spin-up, driven from Python.

The phenomenological vital-rate laws live HERE, in the example (``empirical_laws``),
not in the ``meds`` package — the empirical laws are themselves an experiment
"example", removed from the Fortran core in the reorg (the Fortran model is the
carbon path). This script reimplements them in numpy and drives the refactored
Fortran engine's law-free apply-primitives through the C-API: each step reads the
stand state, computes growth/mortality/recruitment, and calls ``Site.apply_rates``.

It reproduces the golden (``test/golden/empirical_spinup_golden.csv``) and reports the
agreement. That file was RECAPTURED 2026-09-10 from this driver: it was originally taken from
the original Fortran empirical model, which the reorg deleted, and PR #137 then changed the
recruit-pool cadence (it accrues every step now rather than as a monthly lump), which shifted
the first cohorts about a month later and put the run 70% off the old reference in year 2,
decaying to 3% by year 7. Pass ``--emit-golden`` to rewrite it from the current model.

Expected: cohort counts + ``total_agb`` + ``total_nplant`` track the golden closely;
``total_lai`` diverges once cohorts restructure -- the *designed* consequence of the
carbon-only fusion/fission cleanup (carbon restructuring conserves the carbon pools;
the old empirical path re-derived leaf area from the perturbed dbh). See reorg doc S8.

Run from the MEDS repo root:

    MEDS_LIB=build-pylib/libmeds_c.so \\
    LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib:$LD_LIBRARY_PATH \\
    PYTHONPATH=python python3 examples/example_demography/empirical_spinup.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _cadence import slow_steps                                  # noqa: E402
from empirical_laws import empirical_rates                       # noqa: E402

from meds.demography import Config, Site                         # noqa: E402

#----- `--emit-golden` REWRITES the reference CSV from this run. It exists because the golden
#      had no reproducible recapture path: it was taken by hand from the original Fortran
#      empirical model, that model was then deleted in the reorg, and when PR #137 changed the
#      recruit-pool cadence there was no documented way to refresh it. Now there is one.
_ARGV = [a for a in sys.argv[1:] if a != "--emit-golden"]
EMIT_GOLDEN = "--emit-golden" in sys.argv[1:]
CONFIG = _ARGV[0] if _ARGV else "meds_config_main.toml"
GOLDEN = "test/golden/empirical_spinup_golden.csv"
N_YEARS = 40
N_PATCH = 4

_FIELDS = ("dbh", "height", "overtopping_lai", "growth_avg", "agb", "nplant",
           "pft", "owner_patch")


def read_state(site):
    return {k: site.get(k) for k in _FIELDS}


def load_golden():
    g = {}
    with open(GOLDEN) as f:
        next(f)
        for line in f:
            p = line.split(",")
            g[int(p[0])] = dict(n=int(p[1]), agb=float(p[2]), lai=float(p[3]),
                                nplant=float(p[4]))
    return g


def main():
    cfg = Config(CONFIG)
    rows = []
    with Site(cfg, n_patch=N_PATCH) as site:
        year = 0
        for istep, new_month, new_year in slow_steps(cfg.dt_years, N_YEARS):
            state = read_state(site)
            growth, mortality, recruitment = empirical_rates(state, site.n_patch)
            site.apply_rates(growth, mortality, recruitment, new_month, new_year)
            if new_year:
                year += 1
                rows.append(dict(year=year, n=site.n_cohort, agb=site.total_agb,
                                 lai=site.total_lai, nplant=site.total_nplant,
                                 ba=site.total_basal_area))

    if EMIT_GOLDEN:
        with open(GOLDEN, "w") as f:
            f.write("year,n_cohort,total_agb,total_lai,total_nplant,total_basal_area\n")
            for r in rows:
                f.write(f"{r['year']},{r['n']},{r['agb']:>24.16E},{r['lai']:>24.16E},"
                        f"{r['nplant']:>24.16E},{r['ba']:>24.16E}\n")
        print(f"# wrote {GOLDEN} ({len(rows)} years) from THIS run")

    golden = load_golden()
    print("# Python empirical spin-up vs Fortran golden (test/golden)")
    print(f"{'yr':>3} {'n_py':>5} {'n_g':>5} {'agb_py':>13} {'agb_gold':>13} "
          f"{'agb_relerr':>11} {'lai_relerr':>11} {'nplant_relerr':>13}")
    max_agb = max_np = 0.0
    for r in rows:
        g = golden.get(r["year"])
        if not g:
            continue
        ea = abs(r["agb"] - g["agb"]) / max(abs(g["agb"]), 1e-30)
        el = abs(r["lai"] - g["lai"]) / max(abs(g["lai"]), 1e-30)
        en = abs(r["nplant"] - g["nplant"]) / max(abs(g["nplant"]), 1e-30)
        max_agb = max(max_agb, ea)
        max_np = max(max_np, en)
        print(f"{r['year']:>3} {r['n']:>5} {g['n']:>5} {r['agb']:>13.6e} "
              f"{g['agb']:>13.6e} {ea:>11.2e} {el:>11.2e} {en:>13.2e}")

    print(f"# max rel-err vs golden:  total_agb={max_agb:.2e}  total_nplant={max_np:.2e}")
    print("# (total_agb/total_nplant are conserved through carbon restructuring; total_lai")
    print("#  diverges once cohorts fuse/split -- the designed carbon-only cleanup. See")
    print("#  docs/dev_plans/MEDS_REORG_DESIGN.md S8.)")


if __name__ == "__main__":
    main()
