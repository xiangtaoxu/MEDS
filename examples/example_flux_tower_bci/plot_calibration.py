#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_calibration.py -- the fast-parameter calibration at BCI in one figure.

Left: how far the fit moved each fitted key from its default, in prior standard deviations of the
transformed parameter, with the posterior's +-1 sigma; a bar near 0 with a short whisker is a key
the tower confirmed, a long whisker one it could not inform. Right: each target's RMSE, in units
of its observation error, on the validation windows, which the fit never saw, for the default and
the calibrated set. One colour per variant of the fit (calibration.toml [variants]).

Usage: python plot_calibration.py [--fits calibration/fit_*.json] [--out calibration.png]
"""
import argparse
import glob
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts", "calibrate_fast"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))   # the registry reads TOML with meds.config
from registry import SIGMA_U, load_registry  # noqa: E402

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
COLOURS = ["#2a78d6", "#d4661c", "#3a9a5b"]
TARGET_NAMES = {"albedo": "albedo", "lw_up": "LW up", "rnet": "Rnet", "le": "LE", "h": "H", "ef": "EF",
                "gpp": "GPP", "nee_night": "night\nNEE", "ustar": "u*"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fits", default=os.path.join(HERE, "calibration", "fit_*.json"))
    ap.add_argument("--registry", default=os.path.join(HERE, "..", "..", "scripts", "calibrate_fast", "parameters.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "calibration.png"))
    args = ap.parse_args(argv)
    fits = {}
    for path in sorted(glob.glob(args.fits)):
        with open(path) as fh:
            rep = json.load(fh)
        fits[rep.get("variant") or os.path.basename(path)] = rep
    if not fits:
        raise SystemExit(f"no fit results match {args.fits}")
    registry = {p.name: p for p in load_registry(args.registry, "interception_on")}
    keys = []
    for rep in fits.values():
        keys += [k for k in rep["fitted"] if k not in keys]

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, max(4.5, 0.32 * len(keys) + 1.5)),
                                 gridspec_kw={"width_ratios": [1.15, 1.0]})
    n = len(fits)
    y = np.arange(len(keys))
    for i, (variant, rep) in enumerate(fits.items()):
        shift, err = [], []
        for k in keys:
            if k in rep["fitted"]:
                p = registry[k]
                p.default = rep["default"][k]
                shift.append(float((p.to_u(rep["map"][k]) - p.to_u(rep["default"][k])) / SIGMA_U))
                err.append(rep["sigma_ratio"][k])
            else:
                shift.append(np.nan)
                err.append(np.nan)
        off = (i - (n - 1) / 2) * 0.32
        a1.errorbar(shift, y + off, xerr=err, fmt="o", color=COLOURS[i], ms=4, lw=1.2, capsize=2,
                    label=variant.replace("_", " "))
    a1.axvline(0, color=MUTED, lw=0.8)
    a1.set(yticks=y, yticklabels=keys, xlabel="MAP - default, in prior sigma (whisker: posterior sigma)")
    a1.invert_yaxis()
    a1.legend(frameon=False, fontsize=9, loc="lower right")
    a1.set_title("fitted keys", color=INK, fontsize=10)

    targets = [t for t in TARGET_NAMES if any(t in r.get("scores_val", {}).get("default", {}) for r in fits.values())]
    x = np.arange(len(targets))
    width = 0.8 / (n + 1)
    first = next(iter(fits.values()))
    a2.bar(x - 0.4 + width / 2, [first["scores_val"]["default"][t]["nrmse"] for t in targets], width,
           color=GRID, edgecolor=MUTED, label="default")
    for i, (variant, rep) in enumerate(fits.items()):
        a2.bar(x - 0.4 + width * (i + 1.5), [rep["scores_val"]["map"][t]["nrmse"] for t in targets], width,
               color=COLOURS[i], label=f"calibrated, {variant.replace('_', ' ')}")
    top = max(first["scores_val"]["default"][t]["nrmse"] for t in targets)
    a2.set(xticks=x, xticklabels=[TARGET_NAMES[t] for t in targets], ylim=(0, 1.35 * top),
           ylabel="RMSE / observation error, validation windows")
    a2.legend(frameon=False, fontsize=9)
    a2.set_title("the validation windows (never fitted)", color=INK, fontsize=10)
    fig.suptitle("MEDS fast-parameter calibration at Barro Colorado Island (scripts/calibrate_fast)",
                 color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")
    for variant, rep in fits.items():
        print(f"{variant}: Phi {rep['cost']['default']:.4g} -> {rep['cost']['map']:.4g} (calibration windows); "
              f"validation {rep.get('cost_val', {}).get('default', math.nan):.4g} -> "
              f"{rep.get('cost_val', {}).get('map', math.nan):.4g}")


if __name__ == "__main__":
    main()
