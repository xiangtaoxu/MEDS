#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The four leaf-phenology strategies of MEDS, from one kernel over four synthetic climates.

Each strategy is a preset of `meds.plant.pheno` (its cue masks and rates) driven by a synthetic
daily climate. The kernel returns the flush and shed tendencies; `pheno.leaf_step` turns them into
relative LAI and leaf litter as the coupled model's daily leaf update does.

Writes phenology_patterns.png next to this script and prints a summary table. Build libmeds.so first
(see the README).
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import meds  # noqa: F401  -- an installed package (pip install python/)
except ImportError:
    sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))  # or the source tree
import meds.plant.pheno as pheno                                   # noqa: E402

YEAR = 365
N_YEARS = 3                       # one spin-up year, then two plotted years
LAT_TEMPERATE, LAT_TROPICAL = 44.0, 5.0


def seasonal(doy, mean, amplitude, peak_doy):
    return mean + amplitude * math.cos(2.0 * math.pi * (doy - peak_doy) / YEAR)


#----- Synthetic daily climates: each returns the kernel's inputs and the driver to plot. -------#
def temperate(doy):
    """Air temperature from -3 C in winter to 23 C in mid-July; it also stands in for soil."""
    t = seasonal(doy, 283.15, 13.0, 200)
    return dict(temp_day=t, soil_temp=t, daylength=pheno.daylength(LAT_TEMPERATE, doy)), t - 273.15


def tropical_dry(doy):
    """Predawn leaf water potential from -0.3 MPa in the wet season to -2.5 MPa around day 220."""
    psi = seasonal(doy, -1.4, -1.1, 220)
    return dict(temp_day=300.0, soil_temp=300.0, daylength=pheno.daylength(LAT_TROPICAL, doy),
                dmax_leaf_psi=psi), psi


def tropical_light(doy):
    """Daily mean shortwave from 150 W/m2 in the cloudy season to 530 W/m2 around day 220."""
    rad = seasonal(doy, 340.0, 190.0, 220)
    return dict(temp_day=300.0, soil_temp=300.0, daylength=pheno.daylength(LAT_TROPICAL, doy),
                rad=rad), rad


#----- (title, preset, climate, baseline leaf turnover [1/yr], driver). The deciduous canopy
#      turns over little within its season (Harvard Forest litter baskets: ~0.03 per year). ---#
PATTERNS = [
    ("1. Temperate deciduous", pheno.temperate_deciduous(), temperate,
     0.03, "air temperature [°C]"),
    ("2. Temperate evergreen", pheno.temperate_evergreen(flush_cue_mask=pheno.Cue.TEMP), temperate,
     1.0 / 3.0, "air temperature [°C]"),
    ("3. Tropical drought-deciduous", pheno.drought_deciduous(), tropical_dry,
     1.0 / 1.5, "predawn leaf ψ [MPa]"),
    ("4. Tropical light-driven leaf exchange", pheno.light_exchanging(), tropical_light,
     1.0 / 2.0, "shortwave [W m⁻²]"),
]


def simulate(params, climate, turnover):
    """Daily relative LAI, flush and shed tendencies (relative to their maxima), litter and driver
    over the plotted years."""
    ph, lai, rows = pheno.Phenology(params), 1.0, []
    for day in range(N_YEARS * YEAR):
        env, driver = climate(day % YEAR + 1)
        out = ph.step(doy=day % YEAR + 1, **env)
        lai, litter = pheno.leaf_step(lai, out.leaf_flush_rate, out.leaf_shed_rate,
                                      baseline_turnover=turnover / YEAR)
        rows.append((lai, out.leaf_flush_rate / params.k_flush_max,
                     out.leaf_shed_rate / params.k_shed_max, litter, driver))
    return dict(zip(("lai", "flush", "shed", "litter", "driver"), np.array(rows[YEAR:]).T))


def plot(results, out):
    days = np.arange(2 * YEAR)
    month = (days // (YEAR / 12)).astype(int)
    fig, axes = plt.subplots(len(results), 1, figsize=(10.5, 11.5), sharex=True)
    for ax, (title, label, r) in zip(axes, results):
        axr = ax.twinx()                                       # monthly litter behind the lines
        axr.bar(np.arange(24) + 0.5, np.bincount(month, weights=r["litter"]), width=0.75,
                color="#8D6E63", alpha=0.35, label="leaf litter (per month)")
        axr.set_ylim(bottom=0.0)
        axr.set_ylabel("leaf litter\n[canopies / month]", color="#5D4037", fontsize=9)
        axr.tick_params(axis="y", labelsize=8, colors="#5D4037")
        ax.set_zorder(axr.get_zorder() + 1)
        ax.patch.set_visible(False)
        x = days / (YEAR / 12)
        ax.plot(x, r["lai"], color="#2E7D32", lw=2.6, label="relative LAI")
        ax.plot(x, r["flush"], color="#1565C0", lw=1.6, label="flush tendency")
        ax.plot(x, r["shed"], color="#C62828", lw=1.6, ls="--", label="shed tendency")
        ax.set_ylim(-0.03, 1.08)
        ax.set_ylabel("relative [0–1]")
        ax.set_title(f"{title}   ·   driver: {label}", loc="left", fontsize=11, fontweight="bold")
        ax.grid(alpha=0.25)
        ax.axvline(12, color="0.7", lw=0.8, ls=":")
    handles = axes[0].get_legend_handles_labels()
    bars = axr.get_legend_handles_labels()
    axes[0].legend(handles[0] + bars[0], handles[1] + bars[1], loc="upper right", fontsize=8,
                   framealpha=0.9, ncol=2)
    axes[-1].set_xlabel("month (two years after a one-year spin-up)")
    axes[-1].set_xticks(range(0, 25, 3))
    fig.suptitle("Four leaf-phenology strategies from one MEDS kernel", fontsize=13,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


def main():
    results = [(title, label, simulate(params, climate, turnover))
               for title, params, climate, turnover, label in PATTERNS]
    plot(results, os.path.join(HERE, "phenology_patterns.png"))
    print(f"\n{'strategy':40s} {'LAI min':>8s} {'LAI max':>8s} {'litter/yr':>10s}")
    for title, _, r in results:
        print(f"{title:40s} {r['lai'].min():8.2f} {r['lai'].max():8.2f} "
              f"{r['litter'].sum() / 2:10.2f}")


if __name__ == "__main__":
    main()
