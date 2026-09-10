#!/usr/bin/env python3
"""MEDS example_biophysics -- run the FULL coupled model from Python, then build the figures.

This example used to be a bash script that exec'd ``meds_main`` twice and then ran three
plotting scripts over the netCDF it left behind. It now drives the same model in-process
through ``meds.model.Run``: Python owns the time loop, so the stand can be inspected AS IT
GROWS rather than only through the output files.

``Run.step`` calls the identical ``driver_step`` the executable calls -- nothing here
re-implements any physics. See "Reproducibility" in the README for the measured comparison
against a ``meds_main`` run of the same configs.

Two stages, as before:

  1. a 50-year spin-up from bare ground (2024-07-01 -> 2074-07-01), whose only product is the
     restart state file. Skipped automatically when that file is already present.
  2. one July at hourly resolution, restarting from it, which writes the FAST output tier the
     figures are built from.

What the Python loop adds is the fourth figure: the spin-up TRAJECTORY (AGB, LAI, stem density
and soil carbon against year) is collected in memory from the running model, one sample per
simulated year. The Fortran path cannot produce it without turning on an annual netCDF stream
and reading it back.

Usage
-----
    python run_example.py                 # spin-up (if needed) + July + all four figures
    python run_example.py --replot        # skip the model, rebuild the figures from existing output
    python run_example.py --force-spinup  # re-run stage 1 even if the state file exists
    python run_example.py --stage 2       # run only stage 2

Requires the `meds` package (``pip install python/`` from the repo root) plus matplotlib for
the figures.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPINUP_CFG = "meds_config_spinup.toml"
JULY_CFG = "meds_config_july.toml"
STATE = "spinup-S-20740701000000.nc"
FORCING = HERE.parents[1] / "data" / "forcing" / "ithaca_forcing.nc"


#----- The Fortran driver writes its progress straight to fd 1 while Python's stdout is
#      block-buffered whenever the output is redirected. Without flushing at every boundary the
#      two streams come out in the wrong order -- all the Fortran lines, then all the Python ones.
def _require_meds():
    """Import meds.model with an error that says what to do, not just what failed."""
    try:
        from meds.model import Run
    except ImportError as exc:
        sys.exit(
            f"error: cannot import meds.model ({exc}).\n"
            "  Install the package (it compiles and bundles libmeds.so):\n"
            "    pip install ../../python\n"
            "  or point MEDS_LIB at an existing build and put ../../python on PYTHONPATH."
        )
    return Run


def stage1_spinup(Run, force=False):
    """50-year spin-up from bare ground. Returns the annual trajectory."""
    if Path(STATE).exists() and not force:
        print(f"==> stage 1: {STATE} already present -- skipping the 50-year spin-up")
        print("    (--force-spinup to re-run it)", flush=True)
        return None

    print("==> stage 1: 50-year spin-up from bare ground (2024-07-01 -> 2074-07-01)")
    print("    ~9 min on 4 threads. Python samples the stand once per simulated year.",
          flush=True)
    #----- The trajectory is the reason this loop is in Python. Each entry is read straight off
    #      the live model through the C-API -- no netCDF, no parsing of the Fortran log.
    traj = {k: [] for k in ("year", "agb", "lai", "nplant", "soil_carbon", "n_cohort")}
    with Run(SPINUP_CFG) as run:
        for step in run:
            if not step.is_new_year:
                continue
            traj["year"].append(step.date.year)
            traj["agb"].append(run.total_agb)
            traj["lai"].append(run.total_lai)
            traj["nplant"].append(run.total_nplant)
            traj["soil_carbon"].append(run.soil_carbon)
            traj["n_cohort"].append(run.n_cohort)
            if step.date.year % 5 == 0:
                print(f"    {step.date}  cohorts={run.n_cohort:4d}  LAI={run.total_lai:6.3f}"
                      f"  AGB={run.total_agb:7.3f} kgC/m2  soilC={run.soil_carbon:7.3f} kgC/m2",
                      flush=True)
    return traj


def stage2_july(Run):
    """One July at hourly resolution, restarting from the spin-up state."""
    if not Path(STATE).exists():
        sys.exit(f"error: {STATE} not found -- run stage 1 first (drop --stage 2).")
    print("\n==> stage 2: July 2074 at hourly resolution (restart from the spin-up state)",
          flush=True)
    Path("out").mkdir(exist_ok=True)
    with Run(JULY_CFG) as run:
        for step in run:
            #----- Daily, because a single July is 31 steps. The tallest cohort's height is
            #      pulled from the SoA copy-out; it is the cohort plot_biophysics.py tracks.
            heights = run.cohorts("height")["height"]
            tallest = heights.max() if heights.size else float("nan")
            print(f"    {step.date}  cohorts={run.n_cohort:4d}  tallest={tallest:5.2f} m"
                  f"  LAI={run.total_lai:6.3f}", flush=True)


def plot_trajectory(traj, path="spinup_trajectory.png"):
    """The figure only the Python driver can make: 50 years of stand development, in memory."""
    if not traj:
        print(f"==> {path}: no trajectory collected this run (stage 1 was skipped) -- keeping"
              " the existing figure")
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.5), sharex=True)
    panels = [
        ("agb", "above-ground biomass", "kgC m$^{-2}$", "tab:green"),
        ("lai", "leaf area index", "m$^2$ m$^{-2}$", "tab:olive"),
        ("nplant", "stem density", "plants m$^{-2}$", "tab:brown"),
        ("soil_carbon", "soil carbon (7 CENTURY pools)", "kgC m$^{-2}$", "tab:gray"),
    ]
    for ax, (key, title, unit, colour) in zip(axes.flat, panels):
        ax.plot(traj["year"], traj[key], color=colour, lw=1.8)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(unit, fontsize=9)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("year")
    fig.suptitle("MEDS 50-year spin-up at Ithaca NY -- sampled from the running model in Python",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"==> wrote {path}")


def build_figures():
    """The three netCDF figures. They read the FAST output tier, so they are unchanged."""
    for script in ("plot_biophysics.py", "plot_carbon.py", "plot_soil.py"):
        print(f"==> {script}", flush=True)
        subprocess.run([sys.executable, script], check=True, cwd=HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replot", action="store_true",
                    help="skip the model; rebuild the figures from existing output")
    ap.add_argument("--force-spinup", action="store_true",
                    help="re-run stage 1 even if the state file is present")
    ap.add_argument("--stage", choices=("1", "2", "both"), default="both")
    args = ap.parse_args()

    os.chdir(HERE)          # every path in the configs is relative to this directory

    if args.replot:
        build_figures()
        return

    if not FORCING.exists():
        sys.exit(f"error: forcing file {FORCING} not found.\n"
                 "  It is not tracked in git (netCDF files are ignored). Build it with:\n"
                 "    python ../../scripts/download_era5land.py     # needs a CDS API key\n"
                 "    python ../../scripts/prep_era5land_forcing.py")

    Run = _require_meds()
    traj = None
    if args.stage in ("1", "both"):
        traj = stage1_spinup(Run, force=args.force_spinup)
    if args.stage in ("2", "both"):
        stage2_july(Run)

    print()
    plot_trajectory(traj)
    build_figures()


if __name__ == "__main__":
    main()
