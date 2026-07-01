#!/usr/bin/env python3
"""Plot the MEDS reproduction of Slot & Winter (2017, Plant Cell Environ. 40:3055) as ONE figure.

Layout: a large A-Ci demand-curve panel (F. insipida) on the LEFT, and five stacked leaf-temperature
response panels (Vcmax, Jmax, gs, Anet, Rlight) on the RIGHT, species distinguished by COLOUR.

Reads CSVs written by examples/example_leaf_physiology/reproduce_slot2017.py:
  <prefix>_aci.csv        -- ci, ac, aj, anet  (NET rates for F. insipida; ac=RuBP carboxylation,
                             aj=RuBP regeneration, anet=limiting rate = min(ac,aj))
  <prefix>_<species>.csv  -- tleaf_c, vcmax, jmax, rlight, gs, anet

Requires numpy + matplotlib.
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPECIES = ["F_insipida", "L_speciosa", "C_longifolium", "G_madruno"]
TITLES  = ["F. insipida", "L. speciosa", "C. longifolium", "G. madruno"]
COLORS  = ["tab:green", "tab:blue", "tab:orange", "tab:red"]


def _zero_crossing(x, y):
    """First x where y crosses zero (linear interpolation), or None."""
    for k in range(len(x) - 1):
        if y[k] <= 0.0 <= y[k + 1]:
            return x[k] - y[k] * (x[k + 1] - x[k]) / (y[k + 1] - y[k])
    return None


def _crossover(x, a, b):
    """First x where curves a and b cross (a - b changes sign), or None."""
    d = a - b
    for k in range(len(x) - 1):
        if d[k] == 0.0 or d[k] * d[k + 1] < 0.0:
            step = (d[k + 1] - d[k])
            frac = -d[k] / step if step != 0.0 else 0.0
            return x[k] + frac * (x[k + 1] - x[k])
    return None


def plot_combined(prefix, out):
    """One figure: left A-Ci panel (F. insipida) + five stacked temperature panels (colour = species)."""
    fig = plt.figure(figsize=(14, 8))
    grid = fig.add_gridspec(5, 2, width_ratios=[1.55, 1.0], hspace=0.16, wspace=0.24,
                            left=0.06, right=0.985, top=0.90, bottom=0.09)

    #----- Left: A-Ci demand curve for F. insipida (NET rates). --------------------------------#
    axc = fig.add_subplot(grid[:, 0])
    d = np.genfromtxt(f"{prefix}_aci.csv", delimiter=",", names=True)
    ci = d["ci"]
    axc.plot(ci, d["ac"],   color="tab:red",  lw=1.7, label=r"RuBP carboxylation (A$_c$)")
    axc.plot(ci, d["aj"],   color="tab:blue", lw=1.7, label=r"RuBP regeneration (A$_j$)")
    axc.plot(ci, d["anet"], color="black",    lw=2.6, label=r"limiting rate (A$_{net}$)")
    axc.axhline(0.0, color="0.6", lw=0.8)
    xo = _crossover(ci, d["ac"], d["aj"])
    if xo is not None:
        axc.plot([xo], [np.interp(xo, ci, d["anet"])], marker="*", color="purple", ms=16, zorder=5)
    cp = _zero_crossing(ci, d["anet"])
    if cp is not None:
        axc.plot([cp], [0.0], "o", color="black", ms=4)
        axc.annotate(rf"$\Gamma$ = {cp:.0f}", (cp, 0.0), textcoords="offset points",
                     xytext=(8, -14), fontsize=9)
    axc.set_xlabel(r"intercellular CO$_2$   C$_i$ [$\mu$mol mol$^{-1}$]", fontsize=11)
    axc.set_ylabel(r"net photosynthesis  A [$\mu$mol m$^{-2}$ s$^{-1}$]", fontsize=11)
    axc.set_title("(a)  A–C$_i$ demand curve — F. insipida", fontsize=12, loc="left")
    axc.grid(alpha=0.3)
    axc.legend(fontsize=9.5, loc="lower right")
    axc.text(0.035, 0.96, "V$_{cmax}$ = 161\nJ$_{max}$ = 238\n(corrected, 25 °C)",
             transform=axc.transAxes, va="top", fontsize=9.5,
             bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    #----- Right: five stacked leaf-temperature response panels (colour = species). ------------#
    panels = [("vcmax",  r"V$_{cmax}$" + "\n[$\\mu$mol m$^{-2}$s$^{-1}$]", (0, 400)),
              ("jmax",   r"J$_{max}$"  + "\n[$\\mu$mol m$^{-2}$s$^{-1}$]", (0, 260)),
              ("gs",     r"g$_s$"      + "\n[mol m$^{-2}$s$^{-1}$]",       (0, 0.6)),
              ("anet",   r"A$_{net}$"  + "\n[$\\mu$mol m$^{-2}$s$^{-1}$]", (0, 26)),
              ("rlight", r"R$_{light}$"+ "\n[$\\mu$mol m$^{-2}$s$^{-1}$]", (0, 1.6))]
    data = {sp: np.genfromtxt(f"{prefix}_{sp}.csv", delimiter=",", names=True) for sp in SPECIES}
    tags = ["(b)", "(c)", "(d)", "(e)", "(f)"]
    axes = []
    for i, (col, ylab, ylim) in enumerate(panels):
        ax = fig.add_subplot(grid[i, 1])
        for j, sp in enumerate(SPECIES):
            ax.plot(data[sp]["tleaf_c"], data[sp][col], color=COLORS[j], lw=1.8, label=TITLES[j])
        ax.set_ylim(*ylim)
        ax.set_xlim(25, 42)
        ax.grid(alpha=0.3)
        ax.set_ylabel(ylab, fontsize=8.5)
        ax.text(0.015, 0.9, tags[i], transform=ax.transAxes, va="top", fontsize=10)
        if i < 4:
            ax.tick_params(labelbottom=False)
        axes.append(ax)
    axes[-1].set_xlabel("leaf temperature [°C]", fontsize=11)
    axes[0].legend(fontsize=7.5, loc="upper left")

    fig.suptitle("Slot & Winter (2017) reproduced with the MEDS leaf model:   (a) A–C$_i$ demand curve "
                 "(F. insipida);   (b–f) leaf-temperature responses (colour = species)", fontsize=12.5)
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="CSV prefix (e.g. .../slot2017)")
    ap.add_argument("-o", "--out", required=True, help="output PNG")
    args = ap.parse_args()
    plot_combined(args.prefix, args.out)


if __name__ == "__main__":
    main()
