#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MEDS example_flux_tower_bci -- forcing from a flux tower, end to end, at Barro Colorado Island.

Steps (each skipped when its product already exists, unless --force):
  1. fetch the tower data into data/ (Zenodo 6456527, CC0; checked by md5)
  2. build the forcing file with scripts/prepare_flux_tower/make_tower_forcing.py -- the longwave
     filled by the model's synthesis regressed onto the tower, and, when data/bci_era5land.nc exists,
     a second file filled from ERA5-Land
  3. score the longwave fills on hidden observations (compare_longwave_fill.py)
  4. draw the forcing's fill flags (plot_forcing.py)
  5. spin up 50 years from bare ground on the recycled tower years (meds_config_spinup.toml)
  6. run the five tower years with hourly output (meds_config_eval.toml)
  7. compare with the tower in local time (plot_evaluation.py)

Usage:
  python run_example.py                          # everything
  python run_example.py --copy-from ~/BCI_flux   # take the data from a local copy
  python run_example.py --forcing-only           # steps 1-4: no model run
  python run_example.py --meds-main ../../build-ifx/meds_main

Needs numpy, pandas, netCDF4, matplotlib, and a built meds_main.
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
TOOLS = os.path.join(ROOT, "scripts", "prepare_flux_tower")
DATA = os.path.join(HERE, "data")
OUTPUT = os.path.join(HERE, "output")
STATE = os.path.join(OUTPUT, "spinup-S-20120801000000.nc")


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log:
        with open(log, "w") as fh:
            status = subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
        if status:
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")
    else:
        subprocess.check_call(cmd, cwd=HERE)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--copy-from", help="a folder holding BCI_v5.1.csv and README.txt")
    ap.add_argument("--meds-main", default=os.environ.get("MEDS_MAIN", os.path.join(ROOT, "build-ifx", "meds_main")))
    ap.add_argument("--forcing-only", action="store_true", help="stop after the forcing and its figures")
    ap.add_argument("--force", action="store_true", help="redo every step")
    args = ap.parse_args(argv)
    py = sys.executable
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(OUTPUT, exist_ok=True)

    run([py, "fetch_bci_data.py"] + (["--copy-from", args.copy_from] if args.copy_from else []))
    synth = os.path.join(DATA, "bci_forcing_lw-synth.nc")
    if args.force or not os.path.exists(synth):
        run([py, os.path.join(TOOLS, "make_tower_forcing.py"), "--site", "bci_site.toml", "--out", synth,
             "--lw-fill", "synth"])
    era5 = os.path.join(DATA, "bci_era5land.nc")
    compare = [py, os.path.join(TOOLS, "compare_longwave_fill.py"), "--site", "bci_site.toml",
               "--out", os.path.join(DATA, "lw_comparison.json"), "--figure", os.path.join(HERE, "lw_comparison.png")]
    if os.path.exists(era5):
        lw_era5 = os.path.join(DATA, "bci_forcing_lw-era5.nc")
        if args.force or not os.path.exists(lw_era5):
            run([py, os.path.join(TOOLS, "make_tower_forcing.py"), "--site", "bci_site.toml", "--out", lw_era5,
                 "--lw-fill", "era5", "--era5-file", era5])
        compare += ["--era5-file", era5]
    else:
        print(f"note: no {os.path.relpath(era5, HERE)}; the ERA5-Land longwave fill is skipped (README, 'Longwave')")
    run(compare)
    run([py, "plot_forcing.py", "--forcing", synth])
    if args.forcing_only:
        return

    if not os.path.exists(args.meds_main):
        raise SystemExit(f"ERROR: no meds_main at {args.meds_main}; build MEDS or pass --meds-main")
    if args.force or not os.path.exists(STATE):
        run([args.meds_main, "meds_config_spinup.toml"], log=os.path.join(OUTPUT, "spinup.log"))
    run([args.meds_main, "meds_config_eval.toml"], log=os.path.join(OUTPUT, "eval.log"))
    run([py, "plot_evaluation.py"])


if __name__ == "__main__":
    main()
