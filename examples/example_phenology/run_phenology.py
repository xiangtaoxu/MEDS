#!/usr/bin/env python3
"""example_phenology -- drive the MEDS leaf-phenology kernel over synthetic climates.

For each of the FOUR phenological strategies the Fortran kernel (meds_phenology.f90) supports, this
builds an artificial DAILY environmental time series, calls the real kernel once per day through the
`meds.pheno` C-API (which advances the two governor drives), and shows the resulting behaviour:

    1. temperate deciduous            -- flush on spring warmth (GDD), shed on the autumn cold drop
    2. temperate evergreen            -- same seasonal flush, but NO active shed: the canopy persists
    3. tropical drought-deciduous     -- flush/shed keyed on the daily-max leaf water potential
    4. tropical light-driven exchange -- flush stays high while shed RISES with radiation (leaves are
                                         turned over fast but the canopy stays ~full)

For every strategy it plots three RELATIVE quantities in [0,1]:
    (1) relative LAI (canopy fullness),        integrated IN PYTHON from the two rates,
    (2) relative leaf-flush rate = leaf_flush_rate / k_flush_max  (the flush governor drive),
    (3) relative leaf-shed  rate = leaf_shed_rate  / k_shed_max   (the shed  governor drive),
with the driving environmental variable on a secondary axis for context.

The kernel is SIGNAL-only (it never touches leaf mass); the relative LAI is integrated here with
`meds.pheno.integrate_lai`, a compact relative-unit analogue of the Fortran carbon leaf update.

Run (needs the compiled libmeds_plant_c on the search path + the Intel/gfortran runtime):

    cmake -S . -B build-pylib -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON \
          -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
    cmake --build build-pylib --target meds_plant_c
    source /opt/intel/oneapi/setvars.sh
    LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib:$LD_LIBRARY_PATH \
        python examples/example_phenology/run_phenology.py
"""
import math
import sys
from pathlib import Path

# Make `import meds.pheno` work straight from the repo without `pip install -e python/`.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import meds.pheno as pheno   # noqa: E402

YEAR = 365
LAT_TEMPERATE = 44.0    # deg N (temperate: strong day-length + temperature seasonality)
LAT_TROPICAL = 5.0      # deg N (tropical: ~12 h day, weak temperature seasonality)


def daylength(lat_deg, doy):
    """Day length [h] (White et al. 1997), matching meds_time.daylength (incl. the polar fix)."""
    latr = math.radians(lat_deg)
    decl = math.radians(-23.44) * math.cos(2.0 * math.pi / 365.0 * (doy + 9.0))
    arg = -math.tan(latr) * math.tan(decl)
    if arg >= 1.0:
        return 0.0
    if arg <= -1.0:
        return 24.0
    return 24.0 / math.pi * math.acos(arg)


#===========================================================================================#
#  Synthetic daily environments. Each returns (env_kwargs, driver_value_for_plotting).        #
#===========================================================================================#
def temperate_env(day, doy):
    """Ithaca-like seasonal air temperature + real day length; summer peak ~ doy 200."""
    temp = 283.15 + 13.0 * math.sin(2.0 * math.pi * (doy - 109) / YEAR)   # ~270 K winter .. ~296 K summer
    return dict(temp_day=temp, soil_temp=temp - 1.0, daylength=daylength(LAT_TEMPERATE, doy),
                doy=doy, hemis_north=True), temp - 273.15                  # driver plotted in degC


def drought_env(day, doy):
    """Warm tropics; daily-max leaf psi swings between a wet (~-0.3 MPa) and a dry (~-2.5 MPa) season."""
    psi = -1.4 - 1.1 * math.cos(2.0 * math.pi * (doy - 220) / YEAR)       # driest ~ doy 220
    return dict(temp_day=300.0, soil_temp=300.0, daylength=daylength(LAT_TROPICAL, doy),
                doy=doy, dmax_leaf_psi=psi), psi                          # driver = leaf psi [MPa]


def light_env(day, doy):
    """Warm tropics; radiation peaks in the high-light (dry) season, dips in the cloudy (wet) season."""
    rad = 340.0 + 190.0 * math.cos(2.0 * math.pi * (doy - 220) / YEAR)    # ~150 .. ~530 W/m2
    return dict(temp_day=300.0, soil_temp=300.0, daylength=daylength(LAT_TROPICAL, doy),
                doy=doy, rad=rad), rad                                     # driver = radiation [W/m2]


#----- The four strategies: (title, Params, env fn, baseline leaf turnover [1/day], driver label). --#
def leaf_turnover_per_day(years):
    """Baseline (replaceable) leaf turnover expressed per day from a leaf lifespan in years."""
    return 1.0 / (years * YEAR)


PATTERNS = [
    ("1. Temperate deciduous",
     pheno.temperate_deciduous(),
     temperate_env, leaf_turnover_per_day(1.0), "air T (degC)"),
    ("2. Temperate evergreen",
     pheno.temperate_evergreen(flush_cue_mask=pheno.Cue.TEMP),   # seasonal flush, but shed stays OFF
     temperate_env, leaf_turnover_per_day(3.0), "air T (degC)"),
    ("3. Tropical drought-deciduous",
     pheno.drought_deciduous(),
     drought_env, leaf_turnover_per_day(1.5), "leaf psi (MPa)"),
    ("4. Tropical light-driven leaf-exchanging",
     pheno.light_exchanging(),
     light_env, leaf_turnover_per_day(2.0), "radiation (W/m2)"),
]

N_YEARS = 3        # 1-year spin-up + 2 plotted years
SPINUP = 1 * YEAR


def simulate(params, env_fn, baseline_turnover):
    """Run one strategy day by day; return dict of time series (relative LAI + the two rates + driver)."""
    ph = pheno.Phenology(params)
    lai = 1.0
    rec = dict(lai=[], flush_rel=[], shed_rel=[], driver=[])
    kf, ks = params.k_flush_max, params.k_shed_max
    for day in range(N_YEARS * YEAR):
        doy = day % YEAR + 1
        env, driver = env_fn(day, doy)
        out = ph.step(dt=1.0, **env)
        lai = pheno.integrate_lai(lai, out.leaf_flush_rate, out.leaf_shed_rate,
                                  dt=1.0, baseline_turnover=baseline_turnover)
        if day >= SPINUP:                                   # keep only the post-spin-up window
            rec["lai"].append(lai)
            rec["flush_rel"].append(out.leaf_flush_rate / kf)   # = flush_drive in [0,1]
            rec["shed_rel"].append(out.leaf_shed_rate / ks if ks > 0 else 0.0)  # = shed_drive
            rec["driver"].append(driver)
    return rec


def main():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("This example needs matplotlib: pip install matplotlib")

    results = [(title, simulate(params, env_fn, base), drv_label)
               for (title, params, env_fn, base, drv_label) in PATTERNS]

    months = [i / (YEAR / 12.0) for i in range(2 * YEAR)]   # x-axis in months over the 2 plotted years

    fig, axes = plt.subplots(len(PATTERNS), 1, figsize=(10.5, 11.0), sharex=True)
    for ax, (title, rec, drv_label) in zip(axes, results):
        ax.plot(months, rec["lai"], color="#2E7D32", lw=2.6, label="relative LAI (canopy fullness)")
        ax.plot(months, rec["flush_rel"], color="#1565C0", lw=1.6, label="relative flush rate")
        ax.plot(months, rec["shed_rel"], color="#C62828", lw=1.6, ls="--", label="relative shed rate")
        ax.set_ylim(-0.03, 1.08)
        ax.set_ylabel("relative [0-1]")
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
        ax.grid(alpha=0.25)
        for x in (12,):                                     # mark the year boundary
            ax.axvline(x, color="0.7", lw=0.8, ls=":")
        # driver on a secondary axis (light grey for context)
        axd = ax.twinx()
        axd.plot(months, rec["driver"], color="0.55", lw=1.0, alpha=0.7)
        axd.set_ylabel(drv_label, color="0.45", fontsize=9)
        axd.tick_params(axis="y", labelsize=8, colors="0.45")

    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.9, ncol=3)
    axes[-1].set_xlabel("month (two years after a one-year spin-up)")
    axes[-1].set_xticks(range(0, 25, 3))
    fig.suptitle("MEDS leaf phenology: four strategies from the same kernel\n"
                 "(relative LAI integrated in Python from the kernel's two per-day rates)",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out_png = Path(__file__).resolve().parent / "phenology_patterns.png"
    fig.savefig(out_png, dpi=130)
    print(f"wrote {out_png}")

    #----- A compact numeric summary (min/max LAI, mean rates) per strategy. -----------------#
    print(f"\n{'strategy':40s} {'LAI min':>8s} {'LAI max':>8s} {'flush~':>8s} {'shed~':>8s}")
    for title, rec, _ in results:
        n = len(rec["lai"])
        print(f"{title:40s} {min(rec['lai']):8.3f} {max(rec['lai']):8.3f} "
              f"{sum(rec['flush_rel'])/n:8.3f} {sum(rec['shed_rel'])/n:8.3f}")


if __name__ == "__main__":
    main()
