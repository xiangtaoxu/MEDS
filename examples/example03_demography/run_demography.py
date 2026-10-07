#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Drive the MEDS demography engine with the census-trained vital rates.

    python run_demography.py --start census --years 300   # from the 1985 census
    python run_demography.py --start bare --years 300     # from near-bare ground

Each year it records stems and basal area by PFT and size class in output/stand_<start>.csv;
``--write-nc PATH`` also writes the stand itself, one record per year, for the post_proc/ renders.
Run from the repository root (the config's paths are relative to it).
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _cadence import slow_steps                                  # noqa: E402
from census_laws import CensusLaws                               # noqa: E402

from meds.demography import Config, Site                         # noqa: E402

FIELDS = ("dbh", "height", "nplant", "pft", "owner_patch")
N_PATCH_BARE = 4
FIRST_YEAR = {"census": 1985, "bare": 0}
SIZE_CLASSES = [1, 2, 5, 10, 20, 50, 100, 1000]   # dbh class edges [cm], as in census_stand.csv

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--start", choices=("census", "bare"), required=True)
ap.add_argument("--years", type=int, default=300)
ap.add_argument("--config", default=os.path.join(HERE, "example_config_main.toml"))
ap.add_argument("--write-nc", metavar="PATH")
args = ap.parse_args()

cfg = Config(args.config)
laws = CensusLaws(args.config)
writer = None
if args.write_nc:
    from _write_nc import StandWriter
    writer = StandWriter(args.write_nc)
rows = []


def record(site, year, show):
    """Stems [1/ha] and basal area [m2/ha] by PFT and size class, per hectare of the site."""
    dbh, pft = site.get("dbh"), site.get("pft")
    stems = site.get("nplant") * site.get_patch("area")[site.get("owner_patch").astype(int) - 1] * 1e4
    ba = stems * np.pi / 4 * (dbh / 100) ** 2
    size = np.digitize(dbh, SIZE_CLASSES[1:-1])
    for p in (1, 2, 3):
        for k in range(len(SIZE_CLASSES) - 1):
            m = (pft == p) & (size == k)
            rows.append((year, p, SIZE_CLASSES[k], stems[m].sum(), ba[m].sum()))
    if show:
        by_pft = " ".join(f"{ba[pft == p].sum():5.1f}" for p in (1, 2, 3))
        print(f"{year:5d} {site.n_cohort:7d} {site.n_patch:6d} {stems.sum():9.0f} "
              f"{stems[dbh >= 10].sum():7.0f} {ba.sum():7.1f}   {by_pft}")


t0 = time.time()
with Site(cfg, n_patch=N_PATCH_BARE, census=(args.start == "census")) as site:
    print(" year cohorts patches  stems/ha  >=10 cm  BA m2/ha  by PFT 1 2 3")
    year = FIRST_YEAR[args.start]
    record(site, year, True)
    for _, new_month, new_year in slow_steps(cfg.dt_years, args.years):
        state = {k: site.get(k) for k in FIELDS}
        growth, mortality, recruitment = laws.rates(state, site.n_patch)
        site.apply_rates(growth, mortality, recruitment, new_month, new_year)
        if new_year:
            year += 1
            record(site, year, (year - FIRST_YEAR[args.start]) % max(1, args.years // 20) == 0)
            if writer is not None:
                writer.sample(site, year)
    if writer is not None:
        print(f"wrote {writer.write()} ({len(writer.records)} yearly records)")

out = os.path.join(HERE, "output", f"stand_{args.start}.csv")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as fh:
    fh.write("year,pft,dbh_class,stems,basal_area\n")
    fh.writelines(f"{y},{p},{c},{n:.6g},{b:.6g}\n" for y, p, c, n, b in rows)
print(f"wrote {out}; {args.years} years in {time.time() - t0:.0f} s")
