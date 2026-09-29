#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MEDS example_flux_tower_bci -- forcing from a flux tower, end to end, at Barro Colorado Island.

Steps (each skipped when its product already exists, unless --force):
  1. fetch the tower data into data/ (Zenodo 6456527, CC0; checked by md5)
  2. build the forcing file data/bci_forcing.nc with scripts/prepare_flux_tower/make_tower_forcing.py,
     the longwave filled by the model's synthesis regressed onto the tower
  3. score that longwave fill on hidden observations (compare_longwave_fill.py)
  4. draw the forcing's fill flags (plot_forcing.py)
  5. build the census file from the 2010 BCI census (bci_census.toml, scripts/prepare_census)
  6. run the five tower years from the census, with hourly output (meds_config_eval.toml)
  7. compare with the tower in local time: mean diurnal and seasonal cycles (plot_evaluation.py)

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
CENSUS_TOOL = os.path.join(ROOT, "scripts", "prepare_census", "make_census.py")
DATA = os.path.join(HERE, "data")
OUTPUT = os.path.join(HERE, "output")
CENSUS = os.path.join(DATA, "bci_census2010_meds.csv")


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
    forcing = os.path.join(DATA, "bci_forcing.nc")
    if args.force or not os.path.exists(forcing):
        run([py, os.path.join(TOOLS, "make_tower_forcing.py"), "--site", "bci_site.toml", "--out", forcing])
    run([py, os.path.join(TOOLS, "compare_longwave_fill.py"), "--site", "bci_site.toml",
         "--out", os.path.join(DATA, "lw_comparison.json"), "--figure", os.path.join(HERE, "lw_comparison.png")])
    run([py, "plot_forcing.py", "--forcing", forcing])
    if args.forcing_only:
        return

    if not os.path.exists(args.meds_main):
        raise SystemExit(f"ERROR: no meds_main at {args.meds_main}; build MEDS or pass --meds-main")
    if args.force or not os.path.exists(CENSUS):
        run([py, CENSUS_TOOL, "--declaration", "bci_census.toml", "--out", CENSUS])
    run([args.meds_main, "meds_config_eval.toml"], log=os.path.join(OUTPUT, "eval.log"))
    run([py, "plot_evaluation.py"])


if __name__ == "__main__":
    main()
