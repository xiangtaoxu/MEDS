#!/usr/bin/env python3
"""Highlight figure for MEDS' fast (sub-daily) biophysics.

Reads the hourly FAST-tier output written by ``meds_config_july.toml`` and plots the four
temperatures that define the canopy energy balance, over one July at 1 h resolution:

  * air            -- above-canopy air temperature, straight from the met forcing (the boundary
                      condition; every other curve is something MEDS solved for)
  * canopy air     -- the prognostic canopy-air-space (CAS) temperature
  * leaf (tallest) -- leaf temperature of the tallest cohort, i.e. the sunlit upper canopy
  * soil surface   -- top soil-layer temperature

The point of the figure is that the three solved temperatures separate from the forcing in
*different* directions and with *different* phase: sunlit leaves run above air by day and below it
at night (radiative coupling plus transpirational cooling), the canopy air space sits between leaf
and soil, and the soil surface is damped and lagged by its heat capacity. Reproducing that
structure from a met file is the whole job of the fast loop.

Usage
-----
    python plot_biophysics.py                 # writes biophysics_july.png
    python plot_biophysics.py --show          # also open an interactive window

Requires numpy, netCDF4 and matplotlib.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

try:
    import netCDF4 as nc
except ImportError:  # pragma: no cover
    sys.exit("error: netCDF4 is required (conda install netcdf4 / pip install netCDF4)")

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
PATTERN = os.path.join(HERE, "out", "july-F-2074*.nc")
OUTPNG = os.path.join(HERE, "biophysics_july.png")

# One colour per store, chosen so the three solved temperatures read as a family against the
# forcing: the forcing is a neutral slate, the canopy pair is warm (leaf hottest), soil is earthy.
C_AIR = "#54677a"
C_CAS = "#1f8a9b"
C_LEAF = "#d0492b"
C_SOIL = "#8a6a3d"
C_SW = "#e8b52f"


def load(pattern: str) -> dict[str, np.ndarray]:
    """Concatenate the per-day FAST files into flat hourly arrays."""
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(
            f"error: no output found matching {pattern}\n"
            "  Run the model first:  ./run_example.sh"
        )

    cols: dict[str, list] = {k: [] for k in
                             ("air", "cas", "leaf", "soil", "sw", "day", "hour")}
    for path in files:
        with nc.Dataset(path) as d:
            v = d.variables
            get = lambda name: np.asarray(v[name][:]).ravel()

            cols["air"] += list(get("air_temp_fast"))
            cols["cas"] += list(get("cas_temp_site"))
            cols["soil"] += list(get("soil_temp_top_site"))
            cols["sw"] += list(get("sw_in_fast"))
            cols["hour"] += list(get("hour"))
            cols["day"] += list(get("day"))

            # Leaf temperature of the TALLEST cohort, resolved per record: cohort composition
            # changes as the run proceeds, so "cohort 1" is not a stable identity. Slots beyond
            # n_cohort hold fill values, hence the explicit mask rather than a bare argmax.
            leaf = np.asarray(v["leaf_temp_cohort_fast"][:])
            height = np.asarray(v["height_cohort_fast"][:])
            ncoh = np.asarray(v["n_cohort"][:]).ravel().astype(int)
            leaf = np.ma.filled(leaf, np.nan)
            height = np.ma.filled(height, np.nan)
            for t in range(leaf.shape[0]):
                n = max(int(ncoh[t]), 0)
                h = height[t, :n]
                valid = np.isfinite(h)
                if n == 0 or not valid.any():
                    cols["leaf"].append(np.nan)
                else:
                    idx = np.arange(n)[valid][np.nanargmax(h[valid])]
                    cols["leaf"].append(float(leaf[t, idx]))

    out = {k: np.asarray(vals, dtype=float) for k, vals in cols.items()}
    out["t"] = out["day"] + out["hour"] / 24.0  # day-of-month, fractional
    return out


def diel(values: np.ndarray, hour: np.ndarray):
    """Mean and standard deviation by hour of day, NaN-safe."""
    mean = np.full(24, np.nan)
    sd = np.full(24, np.nan)
    for h in range(24):
        sel = values[hour == h]
        sel = sel[np.isfinite(sel)]
        if sel.size:
            mean[h], sd[h] = sel.mean(), sel.std()
    return mean, sd


def celsius(kelvin: np.ndarray) -> np.ndarray:
    return kelvin - 273.15


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    ap.add_argument("--dpi", type=int, default=200, help="output resolution (default 200)")
    args = ap.parse_args()

    d = load(PATTERN)
    n_hours = d["t"].size
    print(f"loaded {n_hours} hourly records "
          f"(day {d['day'].min():.0f}-{d['day'].max():.0f} of July 2074)")

    series = [
        ("Air (forcing)", d["air"], C_AIR, 1.6, "-"),
        ("Canopy air space", d["cas"], C_CAS, 1.3, "-"),
        ("Leaf, tallest cohort", d["leaf"], C_LEAF, 1.3, "-"),
        ("Soil surface", d["soil"], C_SOIL, 1.3, "-"),
    ]

    plt.rcParams.update({
        "font.size": 9,
        "axes.edgecolor": "#3a4147",
        "axes.labelcolor": "#22282c",
        "text.color": "#22282c",
        "xtick.color": "#3a4147",
        "ytick.color": "#3a4147",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(11.0, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92], width_ratios=[1.0, 0.46],
                          hspace=0.34, wspace=0.16,
                          left=0.065, right=0.985, top=0.90, bottom=0.085)
    ax_ts = fig.add_subplot(gs[0, :])
    ax_diel = fig.add_subplot(gs[1, 0])
    ax_dev = fig.add_subplot(gs[1, 1])

    # ---- (1) the full month, hourly ------------------------------------------------------------
    # Shortwave behind the temperatures, as the driver of everything above it.
    ax_sw = ax_ts.twinx()
    ax_sw.fill_between(d["t"], 0.0, d["sw"], color=C_SW, alpha=0.16, linewidth=0, zorder=0)
    ax_sw.set_ylim(0, np.nanmax(d["sw"]) * 3.1)   # squash it into the lower third
    ax_sw.set_yticks([])
    ax_sw.spines[:].set_visible(False)

    for label, vals, colour, lw, ls in series:
        ax_ts.plot(d["t"], celsius(vals), color=colour, lw=lw, ls=ls, label=label,
                   solid_joinstyle="round", zorder=3)

    ax_ts.set_ylabel("temperature  (°C)")
    ax_ts.set_xlabel("day of July")
    ax_ts.set_xlim(d["t"].min(), d["t"].max())
    ax_ts.xaxis.set_major_locator(MultipleLocator(5))
    ax_ts.xaxis.set_minor_locator(MultipleLocator(1))
    ax_ts.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_ts.set_title("Hourly canopy energy balance — July, year 50 of an Ithaca NY simulation",
                    loc="left", fontsize=11.5, pad=8)
    leg = ax_ts.legend(loc="upper left", ncol=4, frameon=False, fontsize=8.6,
                       borderaxespad=0.2, columnspacing=1.4, handlelength=1.8)
    for line in leg.get_lines():
        line.set_linewidth(2.0)
    ax_ts.text(1.0, 1.012, "shaded: incoming shortwave", transform=ax_ts.transAxes,
               ha="right", va="bottom", fontsize=8.0, color="#8a7420")

    # ---- (2) mean diel cycle, with spread ------------------------------------------------------
    hours = np.arange(24)
    for label, vals, colour, lw, ls in series:
        mean, sd = diel(vals, d["hour"])
        ax_diel.fill_between(hours, celsius(mean - sd), celsius(mean + sd),
                             color=colour, alpha=0.13, linewidth=0)
        ax_diel.plot(hours, celsius(mean), color=colour, lw=1.9, ls=ls, zorder=3)

    ax_diel.set_xlabel("hour (UTC)")
    ax_diel.set_ylabel("temperature  (°C)")
    ax_diel.set_xlim(0, 23)
    ax_diel.xaxis.set_major_locator(MultipleLocator(6))
    ax_diel.xaxis.set_minor_locator(MultipleLocator(1))
    ax_diel.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_diel.set_title("Mean diel cycle  (band: ±1 sd across the month)",
                      loc="left", fontsize=9.5, pad=6)

    # ---- (3) departure from the forcing --------------------------------------------------------
    # The actual claim of the figure: each store separates from the driving air temperature by a
    # different amount and at a different time of day. Plotting the differences makes that legible
    # in a way three near-parallel absolute curves never can.
    air_mean, _ = diel(d["air"], d["hour"])
    for label, vals, colour, lw, ls in series[1:]:
        mean, _ = diel(vals, d["hour"])
        ax_dev.plot(hours, mean - air_mean, color=colour, lw=1.9, ls=ls, zorder=3)
    ax_dev.axhline(0.0, color=C_AIR, lw=1.4, zorder=2)
    ax_dev.text(0.6, 0.12, "air", color=C_AIR, fontsize=8.2, va="bottom")

    ax_dev.set_xlabel("hour (UTC)")
    ax_dev.set_ylabel("departure from air  (K)")
    ax_dev.set_xlim(0, 23)
    ax_dev.xaxis.set_major_locator(MultipleLocator(6))
    ax_dev.xaxis.set_minor_locator(MultipleLocator(1))
    ax_dev.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_dev.set_title("Departure from the forcing", loc="left", fontsize=9.5, pad=6)

    fig.savefig(OUTPNG, dpi=args.dpi, facecolor="white")
    print(f"wrote {OUTPNG}")

    # A short numeric summary -- useful on its own, and it keeps the caption honest.
    print("\nJuly means (°C) and diel range:")
    for label, vals, _, _, _ in series:
        v = vals[np.isfinite(vals)]
        mean, _ = diel(vals, d["hour"])
        print(f"  {label:<22s} mean {celsius(v.mean()):6.2f}   "
              f"diel amplitude {np.nanmax(mean) - np.nanmin(mean):5.2f} K   "
              f"[{celsius(v.min()):6.2f}, {celsius(v.max()):6.2f}]")
    lead = np.nanmax(diel(d['leaf'], d['hour'])[0] - air_mean)
    drop = np.nanmin(diel(d['leaf'], d['hour'])[0] - air_mean)
    print(f"\n  tallest-cohort leaf vs air: up to {lead:+.2f} K by day, {drop:+.2f} K at night")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
