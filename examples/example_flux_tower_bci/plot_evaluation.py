#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_evaluation.py -- MEDS at Barro Colorado Island against the tower it was driven by.

The model's hourly records are stamped at the START of their hour in UTC; the tower's half hours
at their start in Panama time (UTC-5). This moves the model to local time -- the post-processing
step that UTC-only forcing leaves to the user -- and averages the tower's two half hours of each
model hour. Turbulent fluxes are compared only where the tower's FLAG says they were measured.

Panels: the mean diurnal cycle (local time) of net radiation, latent and sensible heat and GPP,
over the evaluation years; and the monthly means of latent heat and GPP, for the months with at
least ten days of measured hours (the tower's dry-season records are sparse).

Usage: python plot_evaluation.py [--out evaluation.png]
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
UTC_OFFSET_H = -5.0
MIN_MONTH_HOURS = 240          # a monthly mean needs at least ten days of hours the tower measured
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
MODEL, TOWER = "#2a78d6", INK
PAIRS = [("rnet_fast", "Rnet", "net radiation", "W m⁻²", False),
         ("le_flux_fast", "LE", "latent heat", "W m⁻²", True),
         ("h_flux_fast", "H", "sensible heat", "W m⁻²", True),
         ("gpp_rate_fast", "gpp", "GPP", "µmol m⁻² s⁻¹", True)]


def read_model(pattern):
    """The hourly records as a DataFrame on LOCAL time (the start of each hour)."""
    frames = []
    for path in sorted(glob.glob(pattern)):
        with Dataset(path) as ds:
            cols = {v: np.asarray(ds[v][:], dtype=float).squeeze() for v, *_ in PAIRS}
            stamp = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:],
                                        hour=ds["hour"][:], minute=ds["minute"][:]))
        frames.append(pd.DataFrame(cols, index=stamp + pd.Timedelta(hours=UTC_OFFSET_H)))
    if not frames:
        raise SystemExit(f"ERROR: no model output matches {pattern}; run the evaluation stage first")
    model = pd.concat(frames).sort_index()
    return model.where(model.abs() < 1e30)


def read_tower(path):
    """The tower's half hours averaged to local clock hours; turbulent fluxes only where measured."""
    tower = pd.read_csv(path, parse_dates=["date"], index_col="date")
    for _, col, _, _, turbulent in PAIRS:
        if turbulent:
            tower[col] = tower[col].where(tower["FLAG"] == 1)
    hourly = tower[[col for _, col, *_ in PAIRS]].resample("1h").mean()
    counts = tower[[col for _, col, *_ in PAIRS]].resample("1h").count()
    return hourly.where(counts == 2)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval-F-*.nc"))
    ap.add_argument("--tower", default=os.path.join(HERE, "data", "BCI_v5.1.csv"))
    ap.add_argument("--out", default=os.path.join(HERE, "evaluation.png"))
    args = ap.parse_args(argv)
    model = read_model(args.model)
    tower = read_tower(args.tower).reindex(model.index)
    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(2, 4, figsize=(15, 6.6))
    hour = model.index.hour
    for i, (mv, tv, name, units, _) in enumerate(PAIRS):
        both = model[mv].notna() & tower[tv].notna()
        m = model.loc[both, mv].groupby(hour[both]).mean()
        t = tower.loc[both, tv].groupby(hour[both]).mean()
        ax[0, i].plot(t.index + 0.5, t.values, color=TOWER, lw=2, label="tower")
        ax[0, i].plot(m.index + 0.5, m.values, color=MODEL, lw=1.8, label="MEDS")
        ax[0, i].set(title=f"{name} ({units})", xlabel=f"local time (UTC{UTC_OFFSET_H:+g})", xticks=range(0, 25, 6))
        bias = float((model.loc[both, mv] - tower.loc[both, tv]).mean())
        ax[0, i].text(0.02, 0.95, f"bias {bias:+.1f}", transform=ax[0, i].transAxes, va="top", color=MUTED, fontsize=9)
    ax[0, 0].legend(frameon=False, fontsize=9)
    for j, (mv, tv, name, units) in enumerate([("le_flux_fast", "LE", "latent heat", "W m⁻²"),
                                               ("gpp_rate_fast", "gpp", "GPP", "µmol m⁻² s⁻¹")]):
        both = model[mv].notna() & tower[tv].notna()
        enough = both.resample("MS").sum() >= MIN_MONTH_HOURS      # the tower's measured hours vary a lot
        mm = model.loc[both, mv].resample("MS").mean().where(enough)
        tm = tower.loc[both, tv].resample("MS").mean().where(enough)
        a = plt.subplot(2, 2, 3 + j)
        a.plot(tm.index, tm.values, color=TOWER, lw=2, label="tower")
        a.plot(mm.index, mm.values, color=MODEL, lw=1.8, label="MEDS")
        a.set(title=f"monthly mean {name} ({units})")
    ax[1, 0].set_visible(False) ; ax[1, 1].set_visible(False) ; ax[1, 2].set_visible(False) ; ax[1, 3].set_visible(False)
    fig.suptitle("MEDS driven by the Barro Colorado Island tower, against the same tower "
                 "(turbulent fluxes where FLAG = 1)", color=INK)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
