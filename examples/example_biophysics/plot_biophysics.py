#!/usr/bin/env python3
"""Highlight figure for MEDS' fast (sub-daily) biophysics, for ONE closed-canopy patch.

Plots the four temperatures that define the canopy energy balance over one July at 1 h
resolution:

  * air                 -- above-canopy air temperature, straight from the met forcing (the
                           boundary condition; every other curve is something MEDS solved for)
  * canopy air          -- the prognostic canopy-air-space (CAS) temperature
  * leaf (tallest)      -- leaf temperature of the tallest cohort, i.e. the sunlit upper canopy
  * surface soil layer  -- temperature of the TOP SOIL LAYER (node at ~1.8 cm, 0-4 cm thick).
                           This is a soil temperature, not a skin or litter temperature.

ONE PATCH, NOT THE SITE MEAN. The stand carries a disturbance gap alongside its closed canopy --
patch LAI here spans 0.6 to 5.4 -- and a site mean across that is an average of a shaded forest
floor and a sunlit clearing, which is not a state any part of the forest is in. The gap alone
contributes ~70% of the site-mean soil warm anomaly from ~20% of the area. So the figure selects
the patch with the highest leaf area index and plots only that, which is the case the eye is
being invited to read.

Per-patch sub-daily temperatures come from the opt-in ``[fast].fast_probe`` CSV, because the FAST
netCDF tier is staged as a site mean. The daily tier supplies ``lai_patch`` (which patch to pick)
and ``cohort_offset``/``cohort_count`` (which cohorts are in it, for the leaf curve).

The point of the figure is that the three solved temperatures separate from the forcing in
*different* directions and with *different* phase: sunlit leaves run above air by day and below it
at night (radiative coupling plus transpirational cooling), the canopy air space sits between leaf
and soil, and the surface soil layer is damped and lagged by its heat capacity -- under a closed
canopy it tracks the DAILY MEAN air temperature, running below air by day and above it at night.
Reproducing that structure from a met file is the whole job of the fast loop.

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
DAILY = os.path.join(HERE, "out", "july-D-207407.nc")
PROBE = os.path.join(HERE, "out", "fast_probe.csv")
OUTPNG = os.path.join(HERE, "biophysics_july.png")

# One colour per store, chosen so the three solved temperatures read as a family against the
# forcing: the forcing is a neutral slate, the canopy pair is warm (leaf hottest), soil is earthy.
C_AIR = "#54677a"
C_CAS = "#1f8a9b"
C_LEAF = "#d0492b"
C_SOIL = "#8a6a3d"
C_SW = "#e8b52f"


def pick_patch(daily: str):
    """Return (patch_index, mean LAI, mean area, per-day cohort slice) for the closed-canopy patch.

    Selected on mean July leaf area index, not on area: the question the figure answers is what a
    SHADED forest floor does, so the most-shaded patch is the right one even if it is not the
    largest. Cohort ranges are read per day because fusion and fission can move them within the
    month; `cohort_offset` is 1-based (Fortran), so it is shifted here once and only here.
    """
    if not os.path.exists(daily):
        sys.exit(f"error: {daily} not found -- the daily tier must be on "
                 "([output.daily] enabled = true) so the patch can be identified")
    with nc.Dataset(daily) as d:
        lai = np.ma.filled(np.asarray(d.variables["lai_patch"][:]), np.nan)
        area = np.ma.filled(np.asarray(d.variables["area_patch"][:]), np.nan)
        off = np.ma.filled(np.asarray(d.variables["cohort_offset"][:]), 0).astype(int)
        cnt = np.ma.filled(np.asarray(d.variables["cohort_count"][:]), 0).astype(int)
        day = np.asarray(d.variables["day"][:]).ravel().astype(int)
    lai = np.where(lai > 1e30, np.nan, lai)
    area = np.where(area > 1e30, np.nan, area)
    idx = int(np.nanargmax(np.nanmean(lai, axis=0)))
    rng = {int(day[t]): (off[t, idx] - 1, off[t, idx] - 1 + cnt[t, idx]) for t in range(lai.shape[0])}
    return idx, float(np.nanmean(lai[:, idx])), float(np.nanmean(area[:, idx])), rng


def load(pattern: str, patch: int, coh_range: dict) -> dict[str, np.ndarray]:
    """Hourly series for ONE patch: air and shortwave from the FAST tier, canopy air and surface
    soil from the per-patch probe, leaf from the FAST cohort slab restricted to this patch.

    The FAST tier's scalars are site means and are deliberately NOT used for the solved
    temperatures. Its per-cohort slabs are written by GLOBAL cohort slot, which is what makes the
    daily `cohort_offset`/`cohort_count` usable to pick this patch's cohorts out of them.
    """
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"error: no output found matching {pattern}\n"
                 "  Run the model first:  python run_example.py")

    cols: dict[str, list] = {k: [] for k in ("air", "leaf", "sw", "day", "hour")}
    for path in files:
        with nc.Dataset(path) as d:
            v = d.variables
            get = lambda name: np.asarray(v[name][:]).ravel()
            cols["air"] += list(get("air_temp_fast"))
            cols["sw"] += list(get("sw_in_fast"))
            cols["hour"] += list(get("hour"))
            cols["day"] += list(get("day"))

            # Leaf temperature of the tallest cohort IN THIS PATCH, resolved per record. Cohort
            # composition changes as the run proceeds, so "cohort 1" is not a stable identity, and
            # slots outside the patch (or beyond n_cohort) hold other patches' trees or fill
            # values -- hence the explicit per-day slice plus a finite mask, not a bare argmax.
            leaf = np.ma.filled(np.asarray(v["leaf_temp_cohort_fast"][:]), np.nan)
            height = np.ma.filled(np.asarray(v["height_cohort_fast"][:]), np.nan)
            days = get("day").astype(int)
            for t in range(leaf.shape[0]):
                lo, hi = coh_range.get(int(days[t]), (0, 0))
                hi = min(hi, leaf.shape[1])
                h = height[t, lo:hi]
                valid = np.isfinite(h)
                if hi <= lo or not valid.any():
                    cols["leaf"].append(np.nan)
                else:
                    j = np.arange(lo, hi)[valid][np.nanargmax(h[valid])]
                    cols["leaf"].append(float(leaf[t, j]))

    out = {k: np.asarray(vals, dtype=float) for k, vals in cols.items()}

    #----- Canopy air and surface soil come from the probe, which is the only per-patch sub-daily
    #      record. It samples every fast sub-step, so it is averaged onto the FAST tier's hourly
    #      records rather than subsampled -- those records are themselves 4-sub-step means, and
    #      mixing an instantaneous series with an averaged one would fake a phase difference.
    if not os.path.exists(PROBE):
        sys.exit(f"error: {PROBE} not found -- set [fast].fast_probe = true "
                 "(it is the only per-patch sub-daily output)")
    raw = np.genfromtxt(PROBE, delimiter=",", names=True, dtype=None, encoding="utf-8")
    sel = raw["patch"] == patch + 1                      # probe patch ids are 1-based
    stamp = raw["datetime"][sel]
    key = np.array([f"{t[8:10]}{t[11:13]}" for t in stamp])   # DDHH
    for name, col in (("cas", "cas_temp_K"), ("soil", "soil_temp_top_K")):
        vals = raw[col][sel]
        table = {}
        for k, x in zip(key, vals):
            table.setdefault(k, []).append(x)
        table = {k: float(np.mean(x)) for k, x in table.items()}
        out[name] = np.array([table.get(f"{int(dd):02d}{int(hh):02d}", np.nan)
                              for dd, hh in zip(out["day"], out["hour"])])

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

    patch, lai, area, coh_range = pick_patch(DAILY)
    d = load(PATTERN, patch, coh_range)
    n_hours = d["t"].size
    print(f"loaded {n_hours} hourly records "
          f"(day {d['day'].min():.0f}-{d['day'].max():.0f} of July 2074)")
    print(f"closed-canopy patch: index {patch}, LAI {lai:.2f}, area fraction {area:.3f}")

    series = [
        ("Air (forcing)", d["air"], C_AIR, 1.6, "-"),
        ("Canopy air space", d["cas"], C_CAS, 1.3, "-"),
        ("Leaf, tallest cohort", d["leaf"], C_LEAF, 1.3, "-"),
        ("Surface soil layer", d["soil"], C_SOIL, 1.3, "-"),
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
    #----- Title on its own line, with the patch caveat and the shortwave note sharing the line
    #      below it. The caveat is not decoration: the figure is one patch, not the stand.
    ax_ts.set_title("Hourly canopy energy balance — July, year 50 of an Ithaca NY simulation",
                    loc="left", fontsize=11.5, pad=22)
    ax_ts.text(0.0, 1.012, f"closed-canopy patch only:  LAI {lai:.1f},  {area * 100:.0f}% of stand area",
               transform=ax_ts.transAxes, ha="left", va="bottom", fontsize=8.0, color="#5c6b78")
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
    #----- The surface soil layer's signature under a closed canopy: it sits at the DAILY MEAN air
    #      temperature, below air by day and above it at night. Printed because it is the claim a
    #      reader is most likely to want to check against their own site.
    soil_mean, _ = diel(d["soil"], d["hour"])
    dev = soil_mean - air_mean
    print(f"  surface soil layer vs air:  {np.nanmin(dev):+.2f} K by day, {np.nanmax(dev):+.2f} K at night, "
          f"{np.nanmean(dev):+.2f} K on the daily mean")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
