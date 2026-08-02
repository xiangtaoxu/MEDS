#!/usr/bin/env python3
"""Carbon-flux figure for MEDS' fast (sub-daily) loop.

Reads the same hourly FAST-tier output as ``plot_biophysics.py`` and plots the sub-daily carbon
cycle: the three site fluxes, the canopy-air CO2 they act on, and the tallest cohort's share of
canopy photosynthesis.

  * GPP   -- gross primary productivity
  * NPP   -- GPP net of autotrophic MAINTENANCE respiration (growth respiration is charged in the
             slow allocator, so it is deliberately not subtracted twice)
  * Reco  -- ecosystem respiration, autotrophic + heterotrophic
  * NEE   -- net ecosystem exchange, POSITIVE TO THE ATMOSPHERE (so negative = the stand is a sink)

The claim of the figure is that the model produces a coherent diel carbon cycle rather than a
plausible-looking GPP curve alone: NEE flips sign twice a day, the canopy-air CO2 draws down while
it does, and the two are the same number seen from either side -- the CAS box is driven by exactly
the NEE plotted here, so a disagreement between panels 1 and 3 would be a real inconsistency and
not a plotting artefact.

Usage
-----
    python plot_carbon.py                 # writes carbon_july.png
    python plot_carbon.py --show          # also open an interactive window

Requires numpy, netCDF4 and matplotlib (the `meds` conda environment).
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
    sys.exit("error: netCDF4 is required (conda activate meds)")

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
PATTERN = os.path.join(HERE, "out", "july-F-2074*.nc")
OUTPNG = os.path.join(HERE, "carbon_july.png")

# Categorical hues taken in FIXED SLOT ORDER from the validated default palette (slots 1-4), so a
# series keeps its hue no matter which subset a panel draws. Never cycled, never generated.
C_GPP = "#2a78d6"   # slot 1  blue
C_NPP = "#eb6834"   # slot 2  orange
C_RECO = "#1baf7a"  # slot 3  aqua
C_NEE = "#4a3aa7"   # slot 7  violet -- NEE is the headline series and wants the strongest ink
C_CO2 = "#8a6a3d"   # the CAS store, matched to the soil/earthy family of the energy figure
INK = "#22282c"
INK_MUTED = "#6b7378"


def load(pattern: str) -> dict[str, np.ndarray]:
    """Concatenate the per-day FAST files into flat hourly arrays."""
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"error: no output found matching {pattern}\n  Run the model first:  ./run_example.sh")

    keys = ("gpp", "npp", "reco", "nee", "co2", "day", "hour", "gpp_tall", "sw")
    cols: dict[str, list] = {k: [] for k in keys}
    for path in files:
        with nc.Dataset(path) as d:
            v = d.variables
            get = lambda name: np.asarray(v[name][:]).ravel()

            cols["gpp"] += list(get("gpp_rate_fast"))
            cols["npp"] += list(get("npp_rate_fast"))
            cols["reco"] += list(get("reco_fast"))
            cols["nee"] += list(get("nee_fast"))
            cols["co2"] += list(get("cas_co2_fast"))
            cols["sw"] += list(get("sw_in_fast"))
            cols["hour"] += list(get("hour"))
            cols["day"] += list(get("day"))

            # GPP of the TALLEST cohort, resolved per record. Cohort composition changes as the run
            # proceeds, so "cohort 1" is not a stable identity; slots past n_cohort hold fill.
            # Units differ from the site flux on purpose: the per-cohort variable is per PLANT
            # (umol/plant/s), so it is converted to a per-ground-area contribution below by the
            # caller, not here -- there is no nplant in the FAST stream, so panel 4 reports the
            # per-plant rate itself rather than inventing a conversion.
            gpp_c = np.ma.filled(np.asarray(v["gpp_cohort_fast"][:]), np.nan)
            height = np.ma.filled(np.asarray(v["height_cohort_fast"][:]), np.nan)
            ncoh = np.asarray(v["n_cohort"][:]).ravel().astype(int)
            for t in range(gpp_c.shape[0]):
                n = max(int(ncoh[t]), 0)
                h = height[t, :n]
                ok = np.isfinite(h)
                cols["gpp_tall"].append(float(gpp_c[t, np.arange(n)[ok][np.nanargmax(h[ok])]])
                                        if n and ok.any() else np.nan)

    out = {k: np.asarray(vals, dtype=float) for k, vals in cols.items()}
    out["t"] = out["day"] + out["hour"] / 24.0
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


def style():
    plt.rcParams.update({
        "font.size": 9,
        "axes.edgecolor": "#3a4147",
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": "#3a4147",
        "ytick.color": "#3a4147",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    ap.add_argument("--dpi", type=int, default=200, help="output resolution (default 200)")
    args = ap.parse_args()

    d = load(PATTERN)
    print(f"loaded {d['t'].size} hourly records "
          f"(day {d['day'].min():.0f}-{d['day'].max():.0f} of July 2074)")

    style()
    fig = plt.figure(figsize=(11.0, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.85, 0.85], width_ratios=[1.0, 0.52],
                          hspace=0.46, wspace=0.20,
                          left=0.065, right=0.985, top=0.915, bottom=0.075)
    ax_ts = fig.add_subplot(gs[0, :])
    ax_diel = fig.add_subplot(gs[1, 0])
    ax_co2 = fig.add_subplot(gs[1, 1])
    ax_cum = fig.add_subplot(gs[2, 0])
    ax_tall = fig.add_subplot(gs[2, 1])

    hours = np.arange(24)

    # ---- (1) the full month, hourly -----------------------------------------------------------!
    # GPP is drawn as a filled area because it is bounded at zero and is the envelope the other
    # two live inside; Reco and NEE are lines on the same single axis (never a second y-scale).
    ax_ts.fill_between(d["t"], 0.0, d["gpp"], color=C_GPP, alpha=0.20, linewidth=0, zorder=1)
    ax_ts.plot(d["t"], d["gpp"], color=C_GPP, lw=1.2, zorder=3, label="GPP")
    ax_ts.plot(d["t"], d["reco"], color=C_RECO, lw=1.2, zorder=3, label="Reco")
    ax_ts.plot(d["t"], d["nee"], color=C_NEE, lw=1.2, zorder=4, label="NEE (+ to atmosphere)")
    ax_ts.axhline(0.0, color=INK_MUTED, lw=1.0, zorder=2)

    ax_ts.set_ylabel("CO$_2$ flux  (µmol m$^{-2}$ s$^{-1}$)")
    ax_ts.set_xlabel("day of July")
    ax_ts.set_xlim(d["t"].min(), d["t"].max())
    ax_ts.xaxis.set_major_locator(MultipleLocator(5))
    ax_ts.xaxis.set_minor_locator(MultipleLocator(1))
    ax_ts.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_ts.set_title("Hourly carbon fluxes — July, year 50 of an Ithaca NY simulation",
                    loc="left", fontsize=11.5, pad=8)
    leg = ax_ts.legend(loc="upper left", ncol=3, frameon=False, fontsize=8.6,
                       borderaxespad=0.2, columnspacing=1.6, handlelength=1.8)
    for line in leg.get_lines():
        line.set_linewidth(2.0)

    # ---- (2) mean diel cycle, all four fluxes -------------------------------------------------!
    for label, key, colour in (("GPP", "gpp", C_GPP), ("NPP", "npp", C_NPP),
                               ("Reco", "reco", C_RECO), ("NEE", "nee", C_NEE)):
        mean, sd = diel(d[key], d["hour"])
        ax_diel.fill_between(hours, mean - sd, mean + sd, color=colour, alpha=0.12, linewidth=0)
        ax_diel.plot(hours, mean, color=colour, lw=1.9, zorder=3)
        # Direct labels at the curve's own peak: four series is exactly the ceiling for labelling
        # every one, and it removes the legend-to-curve lookup entirely.
        h_at = int(np.nanargmax(np.abs(mean)))
        # Label ABOVE the curve always: below-the-curve placement clipped NEE against the axis,
        # and a label that runs off the frame is worse than one that sits a little close.
        ax_diel.annotate(label, (h_at, mean[h_at]), textcoords="offset points",
                         xytext=(4, 6), color=colour, fontsize=8.4)
    ax_diel.axhline(0.0, color=INK_MUTED, lw=1.0, zorder=2)
    ax_diel.set_xlabel("hour (UTC)")
    ax_diel.set_ylabel("CO$_2$ flux  (µmol m$^{-2}$ s$^{-1}$)")
    ax_diel.set_xlim(0, 23)
    ax_diel.xaxis.set_major_locator(MultipleLocator(6))
    ax_diel.xaxis.set_minor_locator(MultipleLocator(1))
    ax_diel.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_diel.margins(y=0.16)                     # headroom for the direct labels
    ax_diel.set_title("Mean diel cycle  (band: ±1 sd across the month)",
                      loc="left", fontsize=9.5, pad=6)

    # ---- (3) canopy-air CO2 -- its own panel, NOT a second axis on panel 2 ---------------------!
    # Different unit, so a different chart. This is the state the fluxes above act on: the CAS box
    # is driven by exactly the NEE plotted there, so the drawdown here and the sign flip there are
    # two views of one number.
    co2_mean, co2_sd = diel(d["co2"], d["hour"])
    ax_co2.fill_between(hours, co2_mean - co2_sd, co2_mean + co2_sd,
                        color=C_CO2, alpha=0.14, linewidth=0)
    ax_co2.plot(hours, co2_mean, color=C_CO2, lw=1.9, zorder=3)
    ax_co2.axhline(400.0, color=INK_MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax_co2.text(0.4, 400.0, "free atmosphere", color=INK_MUTED, fontsize=7.8,
                va="bottom", ha="left")
    ax_co2.set_xlabel("hour (UTC)")
    ax_co2.set_ylabel("canopy-air CO$_2$  (µmol mol$^{-1}$)")
    ax_co2.set_xlim(0, 23)
    ax_co2.xaxis.set_major_locator(MultipleLocator(6))
    ax_co2.xaxis.set_minor_locator(MultipleLocator(1))
    ax_co2.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_co2.set_title("Canopy-air CO$_2$ — the state the fluxes act on",
                     loc="left", fontsize=9.5, pad=6)

    # ---- (4) cumulative net carbon over the month ----------------------------------------------!
    # The hourly panel shows a stand that is a sink by day and a source by night; only the running
    # integral says which one wins. umol/m2/s -> gC/m2 over an hour: 3600 s * 12.011 g/mol * 1e-6.
    dt_h = 3600.0 * 12.011e-6
    nee_ok = np.where(np.isfinite(d["nee"]), d["nee"], 0.0)
    cum_nee = np.cumsum(-nee_ok) * dt_h                       # sign-flip: positive = net UPTAKE
    cum_gpp = np.cumsum(np.where(np.isfinite(d["gpp"]), d["gpp"], 0.0)) * dt_h
    cum_reco = np.cumsum(np.where(np.isfinite(d["reco"]), d["reco"], 0.0)) * dt_h
    ax_cum.plot(d["t"], cum_gpp, color=C_GPP, lw=1.8, label="GPP")
    ax_cum.plot(d["t"], cum_reco, color=C_RECO, lw=1.8, label="Reco")
    ax_cum.plot(d["t"], cum_nee, color=C_NEE, lw=2.1, label="net uptake (−NEE)")
    ax_cum.axhline(0.0, color=INK_MUTED, lw=1.0, zorder=2)
    ax_cum.set_xlabel("day of July")
    ax_cum.set_ylabel("cumulative C  (gC m$^{-2}$)")
    ax_cum.set_xlim(d["t"].min(), d["t"].max())
    ax_cum.xaxis.set_major_locator(MultipleLocator(5))
    ax_cum.xaxis.set_minor_locator(MultipleLocator(1))
    ax_cum.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_cum.set_title("Cumulative carbon over the month", loc="left", fontsize=9.5, pad=6)
    leg2 = ax_cum.legend(loc="upper left", ncol=3, frameon=False, fontsize=8.2,
                         borderaxespad=0.2, columnspacing=1.4, handlelength=1.6)
    for line in leg2.get_lines():
        line.set_linewidth(2.0)

    # ---- (5) the tallest cohort's own photosynthesis --------------------------------------------!
    # Per PLANT, not per ground area -- the FAST stream carries no nplant, so converting would mean
    # inventing a number. The point is the shape: one sunlit crown, resolved individually, which is
    # what a cohort model has that a big-leaf model does not.
    tall_mean, tall_sd = diel(d["gpp_tall"], d["hour"])
    ax_tall.fill_between(hours, tall_mean - tall_sd, tall_mean + tall_sd,
                         color=C_GPP, alpha=0.14, linewidth=0)
    ax_tall.plot(hours, tall_mean, color=C_GPP, lw=1.9, zorder=3)
    ax_tall.set_xlabel("hour (UTC)")
    ax_tall.set_ylabel("GPP  (µmol plant$^{-1}$ s$^{-1}$)")
    ax_tall.set_xlim(0, 23)
    ax_tall.xaxis.set_major_locator(MultipleLocator(6))
    ax_tall.xaxis.set_minor_locator(MultipleLocator(1))
    ax_tall.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_tall.set_title("Tallest cohort — per-plant GPP", loc="left", fontsize=9.5, pad=6)

    fig.savefig(OUTPNG, dpi=args.dpi, facecolor="white")
    print(f"wrote {OUTPNG}")

    # A short numeric summary -- keeps the caption honest and doubles as a sanity check.
    gpp_m, _ = diel(d["gpp"], d["hour"])
    nee_m, _ = diel(d["nee"], d["hour"])
    co2_m, _ = diel(d["co2"], d["hour"])
    print("\nJuly diel summary (umol/m2/s unless noted):")
    print(f"  GPP        peak {np.nanmax(gpp_m):6.2f}   night {np.nanmin(gpp_m):6.2f}")
    print(f"  Reco       mean {np.nanmean(diel(d['reco'], d['hour'])[0]):6.2f}")
    print(f"  NEE        min  {np.nanmin(nee_m):6.2f} (uptake)   max {np.nanmax(nee_m):6.2f} (release)")
    print(f"  CAS CO2    min  {np.nanmin(co2_m):6.1f}   max {np.nanmax(co2_m):6.1f} umol/mol")
    print(f"\n  month totals: GPP {cum_gpp[-1]:7.1f}  Reco {cum_reco[-1]:7.1f}  "
          f"net uptake {cum_nee[-1]:7.1f} gC/m2")
    sink_hours = int(np.sum(nee_ok < 0))
    print(f"  the stand is a net sink in {sink_hours} of {nee_ok.size} hours "
          f"({100*sink_hours/nee_ok.size:.0f}%)")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
