#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""One figure: the three census-trained laws (top) and what the engine does with them (bottom).

    python plot_demography.py      # reads census_stand.csv and output/stand_*.csv, writes demography.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census_laws import RateTable                                # noqa: E402

PFT_COLOR = {1: "#008300", 2: "#2a78d6", 3: "#e87ba4"}   # the ED colours post_proc renders with
PFT_NAME = {1: "PFT 1, light wood", 2: "PFT 2", 3: "PFT 3, dense wood"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "legend.frameon": False})

growth = RateTable("growth", "dbh", "bal", log_x=True)
mortality = RateTable("mortality", "dbh", "growth", log_x=True)
recruitment = RateTable("recruitment", "ba_tot", "ba_pft")
cen = pd.read_csv(os.path.join(HERE, "census_stand.csv"))
mc = pd.read_csv(os.path.join(HERE, "output", "stand_census.csv"))
mb = pd.read_csv(os.path.join(HERE, "output", "stand_bare.csv"))

fig, ax = plt.subplots(2, 3, figsize=(11, 6.6), layout="constrained")
dbh = np.geomspace(1, 120, 200)

# (a) growth: open (BAL 0) and shaded (BAL 40 m2/ha)
a = ax[0, 0]
for p in (1, 2, 3):
    for bal, ls in ((0.0, "-"), (40.0, "--")):
        a.plot(dbh, growth(p - 1, dbh, bal), ls, color=PFT_COLOR[p], lw=1.6,
               label=PFT_NAME[p] if bal == 0 else None)
a.set(xscale="log", xlabel="dbh [cm]", ylabel="diameter growth [cm/yr]",
      title="a  Growth (solid: open, dashed: BAL 40 m²/ha)")
a.legend(loc="upper left")

# (b) mortality at the growth each neighbourhood gives
a = ax[0, 1]
for p in (1, 2, 3):
    for bal, ls in ((0.0, "-"), (40.0, "--")):
        a.plot(dbh, 100 * mortality(p - 1, dbh, growth(p - 1, dbh, bal)), ls, color=PFT_COLOR[p], lw=1.6)
    a.annotate(f"PFT {p}", (dbh[0], 100 * mortality(p - 1, dbh[:1], growth(p - 1, dbh[:1], 40.0))[0]),
               xytext=(4, 0), textcoords="offset points", color=INK, va="center", fontsize=8)
a.set(xscale="log", xlabel="dbh [cm]", ylabel="death rate [%/yr]",
      title="b  Mortality at the predicted growth", ylim=(0, None))

# (c) recruitment against the patch's basal area, each PFT holding a third of it
a = ax[0, 2]
ba = np.linspace(0, 80, 161)
for p in (1, 2, 3):
    a.plot(ba, 1e4 * recruitment(p - 1, ba, ba / 3), color=PFT_COLOR[p], lw=1.6)
    a.annotate(f"PFT {p}", (ba[-1], 1e4 * recruitment(p - 1, ba[-1:], ba[-1:] / 3)[0]),
               xytext=(-4, 6), textcoords="offset points", ha="right", color=INK, fontsize=8)
a.set(xlabel="patch basal area [m²/ha]", ylabel="recruits ≥ 1 cm [1/ha/yr]",
      title="c  Recruitment", ylim=(0, None))


def by_pft(t):
    return t.groupby(["year", "pft"]).basal_area.sum().unstack("pft")


# (d) the census start against the next five censuses, then on alone
a = ax[1, 0]
m, c = by_pft(mc), by_pft(cen)
for p in (1, 2, 3):
    a.plot(m.index, m[p], color=PFT_COLOR[p], lw=1.6)
    a.plot(c.index, c[p], "o", ms=5, mfc=PFT_COLOR[p], mec="white", mew=1.0)
    a.annotate(f"PFT {p}", (m.index[-1], m[p].iloc[-1]), xytext=(-4, 5), textcoords="offset points",
               ha="right", color=INK, fontsize=8)
a.set(xlabel="year", ylabel="basal area [m²/ha]", ylim=(0, None),
      title="d  From the 1985 census (dots: censuses)")

# (e) size distribution in 2010, model against census
a = ax[1, 1]
edges = sorted(cen.dbh_class.unique())
x = np.arange(len(edges))
labels = [f"{lo}–{hi}" for lo, hi in zip(edges[:-1], edges[1:])] + [f"≥{edges[-1]}"]
for t, style, name in ((cen[cen.year == 2010], dict(marker="o", ls="none", ms=7, mfc="white", mec=INK, mew=1.4),
                        "census 2010"),
                       (mc[mc.year == 2010], dict(marker="", ls="-", lw=1.8, color=INK), "model 2010"),
                       (mc[mc.year == mc.year.max()], dict(marker="", ls="--", lw=1.8, color=MUTED),
                        f"model {mc.year.max()}")):
    a.plot(x, t.groupby("dbh_class").stems.sum().reindex(edges).values, label=name, **style)
a.set(yscale="log", xticks=x, xticklabels=labels, xlabel="dbh class [cm]", ylabel="stems [1/ha]",
      title="e  Size distribution")
a.legend(loc="upper right")

# (f) the near-bare-ground spin-up, with the census level of each PFT on the right
a = ax[1, 2]
m = by_pft(mb)
for p in (1, 2, 3):
    a.plot(m.index - m.index[0], m[p], color=PFT_COLOR[p], lw=1.6)
    a.annotate(f"PFT {p}", (m.index[-1] - m.index[0], m[p].iloc[-1]), xytext=(-4, 5),
               textcoords="offset points", ha="right", color=INK, fontsize=8)
level = c.mean().mean()
a.axhline(level, color=MUTED, lw=1.0, ls=":")
a.annotate("census, each PFT", (0, level), xytext=(2, 4), textcoords="offset points", color=MUTED, fontsize=8)
a.set(xlabel="years from bare ground", ylabel="basal area [m²/ha]", ylim=(0, None),
      title="f  From near-bare ground")

out = os.path.join(HERE, "demography.png")
fig.savefig(out, dpi=150)
print(f"wrote {out}")
