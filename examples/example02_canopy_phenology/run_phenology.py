#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 02: one MEDS phenology kernel at a deciduous and an evergreen forest.

Harvard Forest (deciduous broadleaf, Massachusetts) and Hyytiala (Scots pine, Finland) run the same
kernel with the same cues: degree-day sums of the daily air temperature, and the day length. There
is no water stress. Only the parameter values differ; the pine keeps its canopy because its
senescence stops at a leaf-cover floor (min_leaf_cover).

The kernel and the leaf rule are the compiled Fortran of the coupled model, reached through
meds.plant.pheno, with carbon never limiting the flush.

Usage (after fetch_phenology_data.py; build libmeds.so first, see the README):
  python run_phenology.py           # run with fitted_parameters.json: print scores, write figures
  python run_phenology.py --fit     # refit both sites first (differential evolution, all cores)
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import meds  # noqa: F401  -- an installed package (pip install python/)
except ImportError:
    sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))  # or the source tree
import meds.plant.pheno as pheno                                   # noqa: E402

DATA = os.path.join(HERE, "data")
FIT_FILE = os.path.join(HERE, "fitted_parameters.json")

#----- Both sites: flush and senescence on temperature and day length. Everything not fitted keeps
#      its meds.plant.pheno default: the 5 C warmth base, the other sharpnesses, 5-day smoothing.
SHARED = dict(flush_cue_mask=pheno.Cue.TEMP | pheno.Cue.LIGHT,
              shed_cue_mask=pheno.Cue.TEMP | pheno.Cue.LIGHT)


def data(name):
    return os.path.join(DATA, name)


def crossing_day(series, first_doy, last_doy):
    """Per year, the (interpolated) day of year the series first rises through 0.5."""
    out = {}
    for year, g in series.groupby(series.index.year):
        g = g[(g.index.dayofyear >= first_doy) & (g.index.dayofyear <= last_doy)]
        k = np.flatnonzero(np.diff(np.sign(g.values - 0.5)) > 0)
        if len(k):
            i = k[0]
            x0, x1 = g.index[i].dayofyear, g.index[i + 1].dayofyear
            out[year] = x0 + (0.5 - g.values[i]) * (x1 - x0) / (g.values[i + 1] - g.values[i])
    return pd.Series(out)


def timing_skill(obs, mod):
    j = obs.index.intersection(mod.index)
    return float(np.corrcoef(obs[j], mod[j])[0, 1]), float(np.sqrt(np.mean((obs[j] - mod[j]) ** 2)))


class Site:
    """A site's daily drivers and observations; simulate() runs the kernel over them."""
    key = title = ""
    lat = cover0 = 0.0
    dates = tair = None

    def free(self):
        """The fitted parameters and their bounds. The flush day length stays an hour below the
        longest day, so the gate opens in summer. Its sharpness is fitted too: a deciduous canopy
        needs a step (a gradual gate is still partly open in October, when the warmth sum is high,
        and refills the canopy while it senesces), a pine a gradual flush window."""
        longest = pheno.daylength(self.lat, 172)
        return {"flush_degree_days": (5.0, 1000.0),       # [K day] warmth centre
                "shed_base_temp": (278.0, 300.0),         # [K]
                "shed_degree_days": (5.0, 600.0),         # [K day] cold centre
                "flush_rate_max": (1 / 60, 1 / 3),        # [1/day]
                "shed_rate_max": (0.0005, 1 / 3),         # [1/day]
                "leaf_turnover_rate": (0.0, 0.6),         # [1/yr]
                "flush_light_threshold": (8.0, longest - 1.0),   # [h]
                "flush_light_sharpness": (0.5, 8.0),             # [1/h]
                "shed_light_threshold": (8.0, longest)}          # [h]

    def simulate(self, values):
        ph = pheno.Phenology(pheno.Params(**SHARED, **values), leaf_cover=self.cover0)
        days = [ph.step(temp_day=t, daylength=pheno.daylength(self.lat, d), doy=d)
                for t, d in zip(self.tair, self.dates.dayofyear)]
        return pd.DataFrame({"cover": [x.leaf_cover for x in days],
                             "litter": [x.senescence + x.background for x in days],
                             "flush": [x.leaf_flush_tendency for x in days],
                             "shed": [x.leaf_shed_tendency for x in days]}, index=self.dates)


class HarvardForest(Site):
    key, title = "harvard_forest", "Harvard Forest (deciduous broadleaf, 42.5° N)"
    lat, cover0 = 42.538, 0.0
    first, last = 2002, 2023                          # 2002 grows the first canopy (spin-up)
    #----- Leaf litter per year as a share of the canopy. The baskets give the timing within each
    #      year, not the amount: a deciduous canopy is built once a year and every leaf falls, so
    #      one canopy. Without it the fit can keep flushing in October, when the warmth sum is still
    #      high, and drop nearly two canopies a year while matching every observed timing.
    ANNUAL_LEAF_FALL = 1.0
    CONIFERS = {"hemlock", "w.pine", "r.pine", "spruce"}
    NOT_LEAVES = ("twig", "bark", "fruit", "flower", "bud", "acorn", "cone", "seed", "non.leaf")

    def __init__(self):
        y0, y1 = self.first, self.last
        met = pd.read_csv(data("hf001-06-daily-m.csv"), parse_dates=["date"], usecols=["date", "airt"])
        met = met.set_index("date").loc[f"{y0}":f"{y1}"].asfreq("D").interpolate()
        self.dates, self.tair = met.index, met.airt.values + 273.15
        #----- MODIS LAI as canopy fullness: a 3-point running median, scaled each year between
        #      its winter level (the evergreen understorey and hemlock) and its summer top.
        lai = pd.read_csv(data("harvard_modis_lai.csv"), parse_dates=["date"]).set_index("date").lai
        lai = lai.dropna().rolling(3, center=True).median().dropna().loc[f"{y0 + 1}":f"{y1}"]
        rel = []
        for _, g in lai.groupby(lai.index.year):
            doy = g.index.dayofyear
            winter = g[(doy <= 90) | (doy >= 330)].median()
            summer = g[(doy >= 170) & (doy <= 240)].quantile(0.9)
            rel.append(((g - winter) / (summer - winter)).clip(-0.1, 1.2))
        self.modis = pd.concat(rel)
        #----- HF003: the mean fraction of leaves fallen on the tagged trees, by date. ---------#
        fall = pd.read_csv(data("hf003-04-fall.csv"), parse_dates=["date"]).dropna(subset=["lfall"])
        fall = fall[(fall.date.dt.year > y0) & (fall.date.dt.year <= y1)]
        self.fall = fall.groupby("date").lfall.mean().clip(0, 100) / 100.0
        #----- HF069: broadleaf litter in the EMS baskets, cumulative within each leaf year. ----#
        lit = pd.read_csv(data("hf069-05-litter.csv"), parse_dates=["date"])
        lit = lit[(lit.site == "ems") & (lit.year > y0) & (lit.year < y1)]
        spp = lit.spp.str.lower()
        lit = lit[~spp.isin(self.CONIFERS) & ~spp.str.contains("|".join(self.NOT_LEAVES))]
        b = lit.groupby(["year", "date"]).wt_g.sum().reset_index()
        b["frac"] = b.groupby("year").wt_g.cumsum() / b.groupby("year").wt_g.transform("sum")
        self.basket = b

    def fallen(self, sim):
        """Fraction of the year's peak canopy lost, from August on (as HF003 counts it)."""
        peak = sim.cover.groupby(sim.index.year).transform("max")
        return (1.0 - sim.cover / peak.clip(lower=1e-6)).where(sim.index.dayofyear >= 213, 0.0)

    def basket_fraction(self, sim):
        """Cumulative fraction of each leaf year's model litter (from 1 April) at the basket dates."""
        cum, out = sim.litter.cumsum(), []
        for year, g in self.basket.groupby("year"):
            start, end = cum.asof(pd.Timestamp(f"{year}-04-01")), cum.asof(g.date.max())
            out.append((cum.reindex(g.date, method="nearest").values - start) / max(end - start, 1e-9))
        return np.concatenate(out)

    def errors(self, sim):
        return {"MODIS LAI": sim.cover.reindex(self.modis.index).values - self.modis.values,
                "leaf fall": self.fallen(sim).reindex(self.fall.index).values - self.fall.values,
                "baskets": self.basket_fraction(sim) - self.basket.frac.values}

    def annual(self, sim):
        return float(sim.litter.loc[f"{self.first + 1}":].sum() / (self.last - self.first))

    def loss(self, sim):
        return (sum(np.nanmean(e ** 2) for e in self.errors(sim).values())
                + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2)

    def timing(self, sim):
        spring = (crossing_day(self.modis, 60, 200), crossing_day(sim.cover, 60, 200))
        autumn = (crossing_day(self.fall, 213, 366), crossing_day(self.fallen(sim), 213, 366))
        return spring, autumn

    def scores(self, sim):
        out = {f"RMSE {k}": float(np.sqrt(np.nanmean(e ** 2))) for k, e in self.errors(sim).items()}
        (s_obs, s_mod), (a_obs, a_mod) = self.timing(sim)
        out["spring half-green r"], out["spring half-green RMSE [d]"] = timing_skill(s_obs, s_mod)
        out["autumn half-fallen r"], out["autumn half-fallen RMSE [d]"] = timing_skill(a_obs, a_mod)
        out["leaf litter [canopies/yr]"] = self.annual(sim)
        return out

    def plot(self, sim, path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(11, 7.5), constrained_layout=True)
        gs = fig.add_gridspec(2, 2)
        ax = fig.add_subplot(gs[0, :])
        win = slice("2014", "2018")
        ax.plot(sim.loc[win].index, sim.cover.loc[win], color="tab:green", lw=2, label="model leaf cover")
        ax.plot(self.modis.loc[win].index, self.modis.loc[win], "k.", ms=4, label="MODIS LAI (scaled)")
        ax.plot(sim.loc[win].index, sim.flush.loc[win], color="tab:blue", lw=0.8, label="flush tendency")
        ax.plot(sim.loc[win].index, sim.shed.loc[win], color="tab:red", lw=0.8, ls="--",
                label="senescence tendency")
        ax.set(ylabel="fraction", title=f"(a) {self.title}, 2014–2018", ylim=(-0.12, 1.25))
        ax.legend(ncol=4, loc="upper center", fontsize=8)

        ax = fig.add_subplot(gs[1, 0])
        (s_obs, s_mod), (a_obs, a_mod) = self.timing(sim)
        ax.plot(s_obs.index, s_obs, "o", color="tab:green", ms=4, label="half green: MODIS")
        ax.plot(s_mod.index, s_mod, "-", color="tab:green", label="half green: model")
        ax.plot(a_obs.index, a_obs, "o", color="tab:brown", ms=4, label="half fallen: HF003")
        ax.plot(a_mod.index, a_mod, "-", color="tab:brown", label="half fallen: model")
        ax.set(xlabel="year", ylabel="day of year", title="(b) spring and autumn, year by year")
        ax.legend(fontsize=8)

        ax = fig.add_subplot(gs[1, 1])
        fallen = self.fallen(sim)
        for _, g in fallen.groupby(fallen.index.year):
            g = g[(g.index.dayofyear >= 230) & (g.index.dayofyear <= 330)]
            ax.plot(g.index.dayofyear, g.values, color="tab:brown", lw=0.6, alpha=0.5)
        ax.plot(self.fall.index.dayofyear, self.fall.values, "k.", ms=3)
        ax.plot([], [], color="tab:brown", label="model, each year")
        ax.plot([], [], "k.", label="HF003, tagged trees")
        ax.set(xlabel="day of year", ylabel="fraction of leaves fallen", xlim=(230, 330),
               title="(c) autumn leaf fall")
        ax.legend(fontsize=8)
        fig.savefig(path, dpi=130)


class Hyytiala(Site):
    key, title = "hyytiala", "Hyytiälä (Scots pine, 61.8° N)"
    lat, cover0 = 61.8475, 1.0
    start, end = "2018-01-01", "2024-07-31"
    #----- Needle fall per year as a share of the canopy. The traps give the timing, not the share:
    #      southern-Finnish Scots pine keeps 3.4-4.2 needle cohorts and needles live about three
    #      years (Pensa & Jalkanen 1999, Silva Fennica 33:654), so about 0.3 of it falls each year.
    ANNUAL_NEEDLE_FALL = 0.30

    def __init__(self):
        met = pd.read_csv(data("hyytiala_fluxnet_dd.csv"), usecols=["TIMESTAMP", "TA_F"], na_values=[-9999])
        met["date"] = pd.to_datetime(met.TIMESTAMP.astype(str), format="%Y%m%d")
        tair = met.set_index("date").TA_F.loc[self.start:self.end].interpolate()
        self.dates, self.tair = tair.index, tair.values + 273.15
        lit = self.needle_litter(data("hyytiala_ancillary.csv"))
        first = self.dates[0] + pd.Timedelta(days=200)              # after a 200-day spin-up
        self.litter = lit[(lit.start >= first) & (lit.end <= self.dates[-1])].reset_index(drop=True)
        rate = (self.litter.mass / self.litter.days).values
        self.share = rate / rate.mean()          # needle fall rate per interval / its mean

    def free(self):
        return {**super().free(), "min_leaf_cover": (0.4, 0.95)}

    @staticmethod
    def needle_litter(path):
        """ICOS ancillary litter: the plot-mean Scots pine foliage litter, per collection interval."""
        a = pd.read_csv(path, dtype=str)
        a = a[a.VARIABLE_GROUP == "GRP_LITTER"]
        g = a.pivot_table(index="GROUP_ID", columns="VARIABLE", values="DATAVALUE", aggfunc="first")
        g = g[(g.LITTER_DEBRIS_STATISTIC == "Mean") & (g.LITTER_DEBRIS_ORGAN == "Foliage")]
        g = g[~g.LITTER_DEBRIS_APPROACH.fillna("").str.contains("Yearly")]
        g = g[g.LITTER_DEBRIS_SPP.fillna("").str.startswith("Pinus sylvestris")]
        end = pd.to_datetime(g.LITTER_DEBRIS_DATE.fillna(g.LITTER_DEBRIS_DATE_END), format="%Y%m%d")
        start = pd.to_datetime(g.LITTER_DEBRIS_DATE_START, format="%Y%m%d", errors="coerce")
        out = pd.DataFrame({"end": end.values, "start": start.values,
                            "mass": g.LITTER_DEBRIS.astype(float).values * 1000.0})   # [g DM m-2]
        out = out.sort_values("end").drop_duplicates("end").reset_index(drop=True)
        out["start"] = out.start.fillna(out.end.shift(1))
        out["days"] = (out.end - out.start).dt.days
        return out.dropna(subset=["start"])

    def model_share(self, sim):
        cum, day = sim.litter.cumsum(), pd.Timedelta(days=1)
        rate = ((cum.reindex(self.litter.end - day, method="nearest").values
                 - cum.reindex(self.litter.start - day, method="nearest").values) / self.litter.days.values)
        return rate / rate.mean()

    def annual(self, sim):
        return float(sim.litter.loc["2019":"2023"].sum() / 5.0)

    def loss(self, sim):
        rmse = np.sqrt(np.mean((self.model_share(sim) - self.share) ** 2))
        return rmse ** 2 + 4.0 * (self.annual(sim) - self.ANNUAL_NEEDLE_FALL) ** 2

    def scores(self, sim):
        mod = self.model_share(sim)
        return {"RMSE needle-fall share": float(np.sqrt(np.mean((mod - self.share) ** 2))),
                "needle-fall share r": float(np.corrcoef(mod, self.share)[0, 1]),
                "needle fall [canopies/yr]": self.annual(sim),
                "lowest leaf cover": float(sim.cover.min())}

    def plot(self, sim, path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 1, figsize=(11, 8.5), constrained_layout=True)
        ax = axs[0]
        ax.plot(sim.index, sim.cover, color="tab:green", lw=2, label="model leaf cover")
        ax.plot(sim.index, sim.flush, color="tab:blue", lw=0.8, label="flush tendency")
        ax.plot(sim.index, sim.shed, color="tab:red", lw=0.8, ls="--", label="senescence tendency")
        ax.set(ylabel="fraction", ylim=(-0.05, 1.1), title=f"(a) {self.title}")
        ax.legend(ncol=3, loc="lower left", fontsize=8)

        ax = axs[1]
        mod, lit = self.model_share(sim), self.litter
        for k in range(len(lit)):
            span = [lit.start[k], lit.end[k]]
            ax.plot(span, [self.share[k]] * 2, color="k", lw=2, label="traps (ICOS)" if k == 0 else None)
            ax.plot(span, [mod[k]] * 2, color="tab:green", lw=2, label="model" if k == 0 else None)
        ax.set(ylabel="needle fall rate / mean", title="(b) needle fall, collection by collection")
        ax.legend(fontsize=8)

        ax = axs[2]
        month = lit.end.dt.month.values
        obs_m = pd.Series(self.share).groupby(month).mean()
        mod_m = pd.Series(mod).groupby(month).mean()
        ax.bar(obs_m.index - 0.2, obs_m.values, 0.4, color="k", label="traps (ICOS)")
        ax.bar(mod_m.index + 0.2, mod_m.reindex(obs_m.index).values, 0.4, color="tab:green", label="model")
        ax.set(xlabel="month the collection ended", ylabel="needle fall rate / mean",
               xticks=range(1, 13), title="(c) mean by month")
        ax.legend(fontsize=8)
        fig.savefig(path, dpi=130)


SITES = {cls.key: cls for cls in (HarvardForest, Hyytiala)}
_LOADED = {}


def site(key):
    """Each site's data, loaded once per process (also in each fitting worker)."""
    if key not in _LOADED:
        _LOADED[key] = SITES[key]()
    return _LOADED[key]


def _loss(x, key, names):
    s = site(key)
    value = s.loss(s.simulate(dict(zip(names, x))))
    return float(value) if np.isfinite(value) else 1e3


def fit(key, workers):
    """Best of three differential-evolution seeds over the site's free parameters."""
    from scipy.optimize import differential_evolution
    free = site(key).free()
    names, best = list(free), None
    for seed in (1, 2, 3):
        r = differential_evolution(_loss, [free[n] for n in names], args=(key, names), popsize=15,
                                   maxiter=120, tol=1e-6, seed=seed, polish=False, workers=workers,
                                   updating="deferred")
        print(f"  {key} seed {seed}: loss {r.fun:.5f}", flush=True)
        if best is None or r.fun < best.fun:
            best = r
    return {n: float(f"{v:.4g}") for n, v in zip(names, best.x)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fit", action="store_true", help="refit both sites and rewrite " + os.path.basename(FIT_FILE))
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="processes for the fit")
    args = ap.parse_args(argv)
    if args.fit:
        fitted = {key: fit(key, args.workers) for key in SITES}
        with open(FIT_FILE, "w") as fh:
            json.dump(fitted, fh, indent=2)
            fh.write("\n")
    with open(FIT_FILE) as fh:
        fitted = json.load(fh)
    for key in SITES:
        s = site(key)
        sim = s.simulate(fitted[key])
        print(f"\n{s.title}")
        for name, value in fitted[key].items():
            print(f"  {name:24s} {value:g}")
        for name, value in s.scores(sim).items():
            print(f"  {name:30s} {value:.3f}")
        s.plot(sim, os.path.join(HERE, f"{key}.png"))
    print(f"\nwrote {', '.join(k + '.png' for k in SITES)}")


if __name__ == "__main__":
    main()
