#!/usr/bin/env python3
"""Soil-moisture figure for MEDS' fast (sub-daily) loop.

Reads the same hourly FAST-tier output as ``plot_biophysics.py`` and shows the soil water column
over one July: depth on the vertical axis, time on the horizontal, moisture in colour.

  * panel 1 -- the whole month as a depth-time field. This is where the vertical structure of a
               drydown is legible: the surface layers respond to every rain event and to the diel
               evaporative cycle, while the deep layers move slowly and monotonically, and the
               boundary between the two migrates downward as the month dries.
  * panel 2 -- four depths as ordinary lines, because a heat map communicates PATTERN but hides
               MAGNITUDE. Any claim about how much water actually moved has to be read here.
  * panel 3 -- the mean profile with its monthly envelope: how far each layer travelled.
  * panel 4 -- the soil TEMPERATURE field on the same axes, for comparison. Heat and water are
               driven through the same surface but with very different penetration depths, and
               putting them on identical axes is the cleanest way to see that.

Only the ACTIVE soil layers are drawn. The output carries all `n_soil_layer_max` slots; the
inactive tail holds fill, and the layer depths come from the file's own ``soil_z`` coordinate
rather than being re-derived from the run configuration.

Usage
-----
    python plot_soil.py                 # writes soil_july.png
    python plot_soil.py --show          # also open an interactive window

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
OUTPNG = os.path.join(HERE, "soil_july.png")

# SEQUENTIAL encoding: one hue, light -> dark, because moisture is a magnitude with no meaningful
# midpoint. Never a rainbow. Water gets blue, heat gets a warm single hue -- two separate ramps for
# two separate quantities, each with its own colour bar, rather than one shared scale.
CMAP_WATER = "Blues"
CMAP_TEMP = "YlOrBr"
INK = "#22282c"
INK_MUTED = "#6b7378"
# Depth traces: categorical hues in FIXED SLOT ORDER from the validated default palette. Four
# depths is exactly the ceiling for direct-labelling every one, so no legend box is needed.
TRACE_COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")


def load(pattern: str):
    """Concatenate the per-day FAST files into (time, layer) fields plus the depth coordinate."""
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"error: no output found matching {pattern}\n  Run the model first:  ./run_example.sh")

    water, temp, day, hour = [], [], [], []
    z = None
    for path in files:
        with nc.Dataset(path) as d:
            v = d.variables
            water.append(np.ma.filled(np.asarray(v["soil_water_site_fast"][:]), np.nan))
            temp.append(np.ma.filled(np.asarray(v["soil_temp_site_fast"][:]), np.nan))
            day += list(np.asarray(v["day"][:]).ravel())
            hour += list(np.asarray(v["hour"][:]).ravel())
            if z is None:
                # The file's OWN depth coordinate. Re-deriving the vertical grid from the run
                # configuration would silently go wrong the moment a run changed soil_depth,
                # layer count or the geometric growth factor.
                z = np.asarray(v["soil_z"][:]).ravel()

    water = np.concatenate(water, axis=0)
    temp = np.concatenate(temp, axis=0)
    t = np.asarray(day, dtype=float) + np.asarray(hour, dtype=float) / 24.0

    # ACTIVE layers only: the inactive tail is exactly zero in the depth coordinate and carries no
    # water. Masking on the coordinate rather than on the data avoids mistaking a genuinely dry
    # layer for an inactive one.
    active = z < 0.0
    return {
        "t": t,
        "hour": np.asarray(hour, dtype=float),
        "z": z[active],
        "water": water[:, active],
        "temp": temp[:, active],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    ap.add_argument("--dpi", type=int, default=200, help="output resolution (default 200)")
    args = ap.parse_args()

    d = load(PATTERN)
    nt, nz = d["water"].shape
    print(f"loaded {nt} hourly records x {nz} active soil layers "
          f"(surface {d['z'][0]:.3f} m to {d['z'][-1]:.3f} m)")

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

    fig = plt.figure(figsize=(11.0, 7.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 0.85, 0.85], width_ratios=[1.0, 0.42],
                          hspace=0.42, wspace=0.22,
                          left=0.068, right=0.955, top=0.915, bottom=0.075)
    ax_field = fig.add_subplot(gs[0, :])
    ax_traces = fig.add_subplot(gs[1, 0])
    ax_prof = fig.add_subplot(gs[1, 1])
    ax_temp = fig.add_subplot(gs[2, :])

    # Layer INTERFACES for pcolormesh, so each cell spans its real layer thickness rather than
    # being centred on the node with a made-up half-width. Interfaces are midpoints between nodes,
    # with the top clamped to the surface and the bottom extrapolated by the last half-thickness.
    z = d["z"]
    edges = np.empty(z.size + 1)
    edges[1:-1] = 0.5 * (z[:-1] + z[1:])
    edges[0] = 0.0
    edges[-1] = z[-1] - (edges[-2] - z[-1])

    t = d["t"]
    t_edges = np.empty(t.size + 1)
    t_edges[1:-1] = 0.5 * (t[:-1] + t[1:])
    t_edges[0] = t[0] - 0.5 / 24.0
    t_edges[-1] = t[-1] + 0.5 / 24.0

    # ---- (1) the depth-time moisture field -----------------------------------------------------!
    # ROBUST colour range (2nd-98th percentile). The raw range is set by a handful of saturated
    # surface cells during rain, which pushes the entire bulk of the field -- the drydown that the
    # panel exists to show -- into the bottom fifth of the ramp. Clipping is stated on the colour
    # bar rather than left for a reader to infer, and the traces in panel 2 carry the unclipped
    # values, so nothing is hidden by it.
    vmin, vmax = np.nanpercentile(d["water"], [2.0, 98.0])
    m = ax_field.pcolormesh(t_edges, edges, d["water"].T, cmap=CMAP_WATER,
                            vmin=vmin, vmax=vmax, shading="flat", rasterized=True)
    cb = fig.colorbar(m, ax=ax_field, pad=0.012, aspect=26, extend="both")
    cb.set_label("volumetric soil moisture  (m$^3$ m$^{-3}$)\n"
                 "colour clipped to the 2nd–98th percentile", fontsize=8.2)
    cb.outline.set_visible(False)
    ax_field.set_ylabel("depth  (m)")
    ax_field.set_xlabel("day of July")
    ax_field.set_xlim(t_edges[0], t_edges[-1])
    ax_field.xaxis.set_major_locator(MultipleLocator(5))
    ax_field.xaxis.set_minor_locator(MultipleLocator(1))
    ax_field.set_title("Soil moisture through July — year 50 of an Ithaca NY simulation",
                       loc="left", fontsize=11.5, pad=8)

    # ---- (2) four depths as lines --------------------------------------------------------------!
    # A heat map shows pattern and hides magnitude; these are the numbers.
    picks = [0, max(1, nz // 4), max(2, nz // 2), nz - 1]
    picks = sorted(set(min(p, nz - 1) for p in picks))
    for colour, k in zip(TRACE_COLOURS, picks):
        ax_traces.plot(t, d["water"][:, k], color=colour, lw=1.5, zorder=3,
                       label=f"{-z[k]:.2f} m")
    # A legend rather than outside-the-axes direct labels: the labels ran into the neighbouring
    # profile panel's axis, and a label that collides with another chart is worse than a legend.
    leg = ax_traces.legend(loc="upper left", ncol=4, frameon=False, fontsize=8.2,
                           borderaxespad=0.2, columnspacing=1.2, handlelength=1.5,
                           title="depth below surface", title_fontsize=8.0)
    leg._legend_box.align = "left"
    for line in leg.get_lines():
        line.set_linewidth(2.0)
    ax_traces.set_xlabel("day of July")
    ax_traces.set_ylabel("moisture  (m$^3$ m$^{-3}$)")
    ax_traces.set_xlim(t.min(), t.max())
    ax_traces.xaxis.set_major_locator(MultipleLocator(5))
    ax_traces.xaxis.set_minor_locator(MultipleLocator(1))
    ax_traces.grid(axis="y", color="#000000", alpha=0.06, lw=0.8)
    ax_traces.margins(y=0.22)                   # headroom for the legend
    ax_traces.set_title("Selected depths — the magnitudes the field above cannot show",
                        loc="left", fontsize=9.5, pad=6)

    # ---- (3) the mean profile and its monthly envelope -----------------------------------------!
    wmean = np.nanmean(d["water"], axis=0)
    wmin = np.nanmin(d["water"], axis=0)
    wmax = np.nanmax(d["water"], axis=0)
    ax_prof.fill_betweenx(z, wmin, wmax, color=TRACE_COLOURS[0], alpha=0.16, linewidth=0)
    ax_prof.plot(wmean, z, color=TRACE_COLOURS[0], lw=1.9, zorder=3)
    ax_prof.plot(wmin, z, color=TRACE_COLOURS[0], lw=0.8, alpha=0.55, zorder=3)
    ax_prof.plot(wmax, z, color=TRACE_COLOURS[0], lw=0.8, alpha=0.55, zorder=3)
    ax_prof.set_xlabel("moisture  (m$^3$ m$^{-3}$)")
    ax_prof.set_ylabel("depth  (m)")
    ax_prof.grid(axis="x", color="#000000", alpha=0.06, lw=0.8)
    ax_prof.set_title("Monthly profile  (band: min–max)", loc="left", fontsize=9.5, pad=6)

    # ---- (4) the temperature field on identical axes -------------------------------------------!
    mt = ax_temp.pcolormesh(t_edges, edges, d["temp"].T - 273.15, cmap=CMAP_TEMP,
                            shading="flat", rasterized=True)
    cbt = fig.colorbar(mt, ax=ax_temp, pad=0.012, aspect=26)
    cbt.set_label("soil temperature  (°C)")
    cbt.outline.set_visible(False)
    ax_temp.set_ylabel("depth  (m)")
    ax_temp.set_xlabel("day of July")
    ax_temp.set_xlim(t_edges[0], t_edges[-1])
    ax_temp.xaxis.set_major_locator(MultipleLocator(5))
    ax_temp.xaxis.set_minor_locator(MultipleLocator(1))
    ax_temp.set_title("Soil temperature, same axes — heat penetrates a different distance than water",
                      loc="left", fontsize=9.5, pad=6)

    fig.savefig(OUTPNG, dpi=args.dpi, facecolor="white")
    print(f"wrote {OUTPNG}")

    # Numeric summary -- the claims the figure makes, as numbers.
    print("\nJuly soil-moisture summary:")
    print(f"  {'depth (m)':>10s}  {'mean':>7s}  {'min':>7s}  {'max':>7s}  {'range':>7s}")
    for k in range(nz):
        print(f"  {-z[k]:10.3f}  {wmean[k]:7.4f}  {wmin[k]:7.4f}  {wmax[k]:7.4f}  "
              f"{wmax[k]-wmin[k]:7.4f}")
    surf, deep = 0, nz - 1
    print(f"\n  surface layer travels {wmax[surf]-wmin[surf]:.4f} m3/m3 over the month; "
          f"the deepest travels {wmax[deep]-wmin[deep]:.4f}")
    dry = d["water"][-1] - d["water"][0]
    print(f"  net change over July: {dry[surf]:+.4f} at the surface, {dry[deep]:+.4f} at depth")
    tmean = np.nanmean(d["temp"], axis=0) - 273.15
    trange = np.nanmax(d["temp"], axis=0) - np.nanmin(d["temp"], axis=0)
    print(f"  soil temperature: {tmean[surf]:.2f} degC at the surface "
          f"(range {trange[surf]:.2f} K), {tmean[deep]:.2f} degC at depth "
          f"(range {trange[deep]:.2f} K)")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
