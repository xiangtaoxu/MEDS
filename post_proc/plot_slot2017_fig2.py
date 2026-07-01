#!/usr/bin/env python3
"""Reproduce Figure 2 of Slot & Winter (2017, Plant Cell Environ. 40:3055) with the MEDS
leaf-physiology model.

Reads the per-species CSVs written by reproduce_slot2017.py
(<prefix>_<species>.csv, columns tleaf_c, vcmax_m4, vcmax_m3, jmax_m4, jmax_m3, rlight, gs, anet)
and renders the 5x4 grid of leaf-temperature responses (lines only, no observed points):

  row a  VCMax  -- the model's peaked temperature response from Table 2 (solid = 4-parameter fit,
                  dashed = Hd fixed at 200 kJ/mol). EXACT reproduction of the paper's fitted lines.
  row b  JMax   -- same, from Table 2.
  row c  gs     -- OUTPUT of the coupled A-gs-Ci solver (Medlyn stomata).
  row d  ANet   -- OUTPUT of the coupled solver at PAR=1500, Ca=400 ppm with a rising leaf VPD.
  row e  RLight -- the model's Rd(T).

Usage:
    python plot_slot2017_fig2.py PREFIX [-o OUT.png]

Requires the `meds` conda environment (numpy, matplotlib).
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPECIES = ["F_insipida", "L_speciosa", "C_longifolium", "G_madruno"]
TITLES = ["F. insipida", "L. speciosa", "C. longifolium", "G. madruno"]


def plot(prefix, out):
    """Render the 5x4 Figure-2 grid from the per-species CSVs at `prefix`_<species>.csv."""
    fig, ax = plt.subplots(5, 4, figsize=(12, 13), sharex=True)

    for j, sp in enumerate(SPECIES):
        d = np.genfromtxt(f"{prefix}_{sp}.csv", delimiter=",", names=True)
        t = d["tleaf_c"]
        #----- a: VCMax (green); solid = 4-param, dashed = Hd=200. -------------------------#
        ax[0, j].plot(t, d["vcmax_m4"], color="tab:green", lw=2)
        ax[0, j].plot(t, d["vcmax_m3"], color="tab:green", lw=1.2, ls="--")
        #----- b: JMax (blue). ------------------------------------------------------------#
        ax[1, j].plot(t, d["jmax_m4"], color="tab:blue", lw=2)
        ax[1, j].plot(t, d["jmax_m3"], color="tab:blue", lw=1.2, ls="--")
        #----- c: gs (red), model output. -------------------------------------------------#
        ax[2, j].plot(t, d["gs"], color="tab:red", lw=2)
        #----- d: ANet (black), model output. ---------------------------------------------#
        ax[3, j].plot(t, d["anet"], color="black", lw=2)
        #----- e: RLight (brown), model Rd(T). --------------------------------------------#
        ax[4, j].plot(t, d["rlight"], color="saddlebrown", lw=2)
        ax[0, j].set_title(TITLES[j], style="italic", fontsize=12)

    #----- Row y-limits (paper Fig 2 ranges) and labels. ----------------------------------#
    ylims = [(0, 400), (0, 260), (0, 0.6), (0, 26), (0, 1.6)]
    ylabels = [r"V$_{cmax}$ [$\mu$mol m$^{-2}$ s$^{-1}$]",
               r"J$_{max}$ [$\mu$mol m$^{-2}$ s$^{-1}$]",
               r"g$_s$ [mol m$^{-2}$ s$^{-1}$]",
               r"A$_{net}$ [$\mu$mol m$^{-2}$ s$^{-1}$]",
               r"R$_{light}$ [$\mu$mol m$^{-2}$ s$^{-1}$]"]
    for r in range(5):
        for j in range(4):
            ax[r, j].set_ylim(*ylims[r])
            ax[r, j].grid(alpha=0.3)
        ax[r, 0].set_ylabel(ylabels[r], fontsize=10)
    for j in range(4):
        ax[4, j].set_xlabel("leaf temperature [°C]")
    ax[0, 0].set_xlim(25, 42)

    #----- Annotate the two fit variants on the VCMax row. --------------------------------#
    ax[0, 3].plot([], [], color="tab:green", lw=2, label="4-parameter fit")
    ax[0, 3].plot([], [], color="tab:green", lw=1.2, ls="--", label="H$_d$ = 200 kJ/mol")
    ax[0, 3].legend(fontsize=7, loc="upper left")

    fig.suptitle("Slot & Winter (2017) Fig. 2 reproduced with the MEDS leaf-physiology model\n"
                 "rows a,b: model peaked T-response from Table 2 (exact);  "
                 "rows c,d: coupled-solver output (PAR=1500, Ca=400);  row e: model R$_d$(T)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="CSV prefix (e.g. .../slot2017)")
    ap.add_argument("-o", "--out", default=None, help="output PNG (default: <prefix>_fig2.png)")
    args = ap.parse_args()
    plot(args.prefix, args.out or (args.prefix + "_fig2.png"))


if __name__ == "__main__":
    main()
