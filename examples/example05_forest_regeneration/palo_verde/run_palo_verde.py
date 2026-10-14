#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MEDS example 05, Palo Verde -- a tropical dry forest regrowing from bare ground for 150 years.

Steps:
  1. cut Palo Verde's ERA5-Land cell, 2003-2022, from the archive into data/
     (scripts/prepare_era5/make_forcing_file.py); kept once made
  2. the PFTs' check: the 2009 inventory, July 2008 - December 2010 (meds_config_inventory.toml)
  3. the regrowth: 150 years from bare ground, June 1871 - May 2021 (meds_config_regrowth.toml)
  4. draw the regrowth (plot_palo_verde.py -> palo_verde.png)

Usage:
  python run_palo_verde.py --era5-archive /path/to/ED_ERA5land   # everything
  python run_palo_verde.py --skip-inventory                      # the regrowth only
  python run_palo_verde.py --plot-only

Needs numpy, netCDF4, matplotlib and a built meds_main (--meds-main, default ../../../build-ifx/meds_main).
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
FORCING = os.path.join(HERE, "data", "pv_era5land_2003-2022.nc")
CUT_TOOL = os.path.join(ROOT, "scripts", "prepare_era5", "make_forcing_file.py")
SITE = (10.35, -85.35)             # [deg N, deg E] the Palo Verde plots, as meds_config_*.toml [site]


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log is None:
        subprocess.check_call(cmd, cwd=HERE)
        return
    with open(log, "w") as fh:
        if subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT):
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"),
                    help="the processed ED_ERA5land archive, to cut the forcing from (step 1)")
    ap.add_argument("--meds-main", default=os.environ.get("MEDS_MAIN", os.path.join(ROOT, "build-ifx", "meds_main")))
    ap.add_argument("--skip-inventory", action="store_true", help="skip the 2009 inventory check (step 2)")
    ap.add_argument("--plot-only", action="store_true", help="only redraw the figure")
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(HERE, "output"), exist_ok=True)

    if not args.plot_only:
        if not os.path.exists(FORCING):
            if not args.era5_archive:
                raise SystemExit(f"ERROR: no {FORCING}; pass --era5-archive to cut it from ERA5-Land")
            run([sys.executable, CUT_TOOL, "--data-path", args.era5_archive, "--start", "2003-01-01",
                 "--end", "2022-12-31", "--lat", str(SITE[0]), "--lon", str(SITE[1]), "--out", FORCING])
        if not args.skip_inventory:
            run([args.meds_main, "meds_config_inventory.toml"], log=os.path.join(HERE, "output", "inventory.log"))
        run([args.meds_main, "meds_config_regrowth.toml"], log=os.path.join(HERE, "output", "regrowth.log"))
    run([sys.executable, os.path.join(HERE, "plot_palo_verde.py")])


if __name__ == "__main__":
    main()
