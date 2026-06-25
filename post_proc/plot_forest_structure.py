#!/usr/bin/env python3
"""Animate MEDS forest structure as a canopy-layer stand profile.

MEDS is spatially implicit: within a patch a cohort's canopy is a thin horizontal disk that
covers the WHOLE patch area at the cohort's height (ED2's flat-canopy assumption). This script
renders that picture directly. For each output record it draws a pseudo-spatial cross-section
of the site:

  * The x-axis tiles the patches left -> right by DESCENDING age (oldest stand on the left,
    youngest treefall gap on the right); each patch's width is proportional to its area, so the
    row of patches is a faithful area cross-section of the unit site.
  * Each COHORT is a thin horizontal RECTANGLE spanning its patch's full width (the canopy disk
    seen edge-on), vertically centred at the cohort's height; the rectangle's THICKNESS is
    proportional to the cohort's leaf area index (LAI = nplant * leaf_area). Colour encodes PFT.
    Stacking these by height gives the patch's vertical canopy structure, and the patch row gives
    the horizontal (age-mosaic) heterogeneity.
  * Patches are keyed to their persistent `global_patch_id`, so a surviving patch keeps a stable
    slot (oldest-first ordering, ties broken by id) as the mosaic reshuffles frame to frame.

The frames are written to an animated GIF.

Usage:
    python plot_forest_structure.py STATE.nc [-o OUT.gif] [--fps N] [--lai-scale M] [--min-thick T]

Requires the `meds` conda environment (numpy, matplotlib, netCDF4, Pillow). The netCDF must
carry `global_patch_id` and per-cohort `leaf_area` (MEDS writes these from -DMEDS_ENABLE_IO=ON).
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")                       # headless: write a GIF, no display needed
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Patch
from netCDF4 import Dataset

DIST_LABEL = {1: "primary", 2: "treefall"}


def patch_bands(area, age, gid, area_floor=2.0e-3):
    """Left edge and width on a [0,1] axis for each patch, ordered oldest (left) -> youngest.

    Width is proportional to area (a small floor keeps slivers visible); ties in age break by
    global id so the ordering is stable across frames. Returns dict: local patch index -> (left,
    width), plus the draw order (list of local indices, left to right)."""
    order = sorted(range(len(area)), key=lambda ip: (-age[ip], gid[ip]))
    w = np.array([area[ip] + area_floor for ip in order], dtype=float)
    w /= w.sum()
    lefts = np.concatenate([[0.0], np.cumsum(w)[:-1]])
    band = {ip: (float(lefts[s]), float(w[s])) for s, ip in enumerate(order)}
    return band, order


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ncfile", help="MEDS netCDF output file")
    ap.add_argument("-o", "--out", default=None, help="output GIF (default: <ncfile>_forest.gif)")
    ap.add_argument("--fps", type=float, default=5.0, help="frames per second (default 5)")
    ap.add_argument("--lai-scale", type=float, default=2.0,
                    help="canopy thickness in metres per unit cohort LAI (default 2)")
    ap.add_argument("--min-thick", type=float, default=0.2,
                    help="minimum drawn canopy thickness in metres (default 0.2)")
    ap.add_argument("--stride", type=int, default=1,
                    help="plot every Nth record (default 1; the final record is always included)")
    ap.add_argument("--dpi", type=int, default=100, help="output resolution (default 100)")
    args = ap.parse_args()
    out = args.out or (args.ncfile.rsplit(".", 1)[0] + "_forest.gif")

    ds = Dataset(args.ncfile)
    t       = ds.variables["time"][:]
    ncoh    = ds.variables["n_cohort"][:].astype(int)
    npat    = ds.variables["n_patch"][:].astype(int)
    c_pft   = ds.variables["pft"][:]
    c_hgt   = ds.variables["height"][:]
    c_nplant = ds.variables["nplant"][:]
    c_la    = ds.variables["leaf_area"][:]                        # per-plant leaf area [m2]
    c_owner = ds.variables["owner_patch"][:].astype(int)         # 1-based local patch index
    p_area  = ds.variables["patch_area"][:]
    p_age   = ds.variables["patch_age"][:]
    p_dist  = ds.variables["dist_type"][:].astype(int)
    p_gid   = ds.variables["global_patch_id"][:].astype(int)
    nt = ds.dimensions["time"].size

    # Reduce over each record's VALID prefix [:nc] only -- the cohort dimension has a fill tail
    # (unused slots hold the netCDF _FillValue ~1e37), which must not leak into the limits.
    npft, ytop = 1, 5.0
    for r in range(nt):
        nc = int(ncoh[r])
        if nc > 0:
            npft = max(npft, int(np.max(c_pft[r, :nc])))
            lai = np.asarray(c_nplant[r, :nc] * c_la[r, :nc], dtype=float)
            depth = np.maximum(args.min_thick, lai * args.lai_scale)
            ytop = max(ytop, float(np.max(np.asarray(c_hgt[r, :nc], dtype=float) + depth / 2.0)))
    ytop *= 1.05
    pft_colors = plt.cm.YlGn(np.linspace(0.45, 0.95, npft))      # pale -> deep green by PFT
    legend_handles = [Patch(facecolor=pft_colors[i], edgecolor="0.3", label=f"PFT {i + 1}")
                      for i in range(npft)]

    fig, ax = plt.subplots(figsize=(12.0, 6.0))
    fig.subplots_adjust(left=0.06, right=0.84, top=0.92, bottom=0.09)   # room for an outside legend

    def draw(r):
        ax.clear()
        nc, npt = int(ncoh[r]), int(npat[r])
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, ytop)
        ax.set_ylabel("height (m)")
        ax.set_xlabel("patches: oldest (left) -> youngest (right); width $\\propto$ area")
        ax.set_xticks([])
        ax.set_title(f"MEDS forest structure — year {t[r]:.0f}   "
                     f"({nc} cohorts, {npt} patches)", fontsize=12)

        if npt == 0:
            return
        area = p_area[r, :npt]; age = p_age[r, :npt]
        dist = p_dist[r, :npt]; pgid = p_gid[r, :npt]
        band, order = patch_bands(area, age, pgid)

        #----- Patch backgrounds, separators, and age/disturbance labels. -----------------#
        for s, ip in enumerate(order):
            left, w = band[ip]
            ax.axvspan(left, left + w, color="0.5", alpha=0.05 if s % 2 else 0.11, zorder=0)
            if s > 0:
                ax.axvline(left, color="0.7", lw=0.6, zorder=0.5)
            tag = DIST_LABEL.get(int(dist[ip]), "?")[0].upper()
            if w > 0.035:                                        # skip text on slivers it can't fit
                ax.text(left + w / 2.0, ytop * 0.98, f"{tag}\n{age[ip]:.0f} yr",
                        ha="center", va="top", fontsize=7, color="0.35")
        ax.axhline(0.0, color="0.3", lw=1.0, zorder=1)            # ground line

        if nc == 0:
            return
        #----- Each cohort: a patch-spanning canopy disk (rectangle) at its height, thickness #
        #      proportional to its LAI. Drawn shortest-first so taller layers read in front.  #
        owner = c_owner[r, :nc] - 1                              # 0-based local patch index
        hs  = np.asarray(c_hgt[r, :nc], dtype=float)
        lai = np.asarray(c_nplant[r, :nc] * c_la[r, :nc], dtype=float)
        pft = np.asarray(c_pft[r, :nc], dtype=int)
        depth = np.maximum(args.min_thick, lai * args.lai_scale)
        lefts = np.array([band[o][0] for o in owner])
        widths = np.array([band[o][1] for o in owner])
        zsort = np.argsort(hs)
        ax.barh(hs[zsort], widths[zsort], left=lefts[zsort], height=depth[zsort],
                align="center", color=pft_colors[pft[zsort] - 1], alpha=0.7,
                edgecolor="0.25", linewidth=0.3, zorder=2)

        ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  fontsize=8, framealpha=0.9,
                  title="canopy disk = one cohort\n(spans its patch; thickness\n$\\propto$ LAI,"
                  " colour = PFT)\n\npatch tag:\n P primary  T treefall", title_fontsize=7)

    frames = sorted(set(range(0, nt, max(1, args.stride))) | {nt - 1}) if nt else [0]
    ani = FuncAnimation(fig, draw, frames=frames, interval=1000.0 / args.fps)
    ani.save(out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    print(f"wrote {out}  ({len(frames)} frames of {nt} records)")
    ds.close()


if __name__ == "__main__":
    main()
