#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 02: one MEDS phenology kernel, three cues, four forests.

  Harvard Forest  deciduous broadleaf: warmth and the hours of light at a low par_min (the photoperiod)
  Hyytiala        Scots pine (evergreen): warmth and the bright hours; needle fall stops at a floor
  BCI             a light-driven leaf exchanger: senescence on many bright hours in the dry season
  Palo Verde      a drought-deciduous forest: predawn leaf water potential and the hours of light

Every site's light is ERA5-Land through MEDS's own forcing reader: the PAR at the canopy top every
900 s, counted as hours a day above the site's par_min (drivers/era5_par_hours_*.csv.gz, made by
make_era5_par_hours.py). The kernel and the leaf rule are the compiled Fortran of the coupled model,
reached through meds.plant.pheno, with carbon never limiting the flush. The water potentials are
stand-ins: at BCI the tower's soil water content through a retention curve, at Palo Verde a MEDS
run's canopy predawn leaf water potential (drivers/palo_verde_predawn_psi.csv).

Usage (after fetch_phenology_data.py; build libmeds.so first, see the README):
  python run_phenology.py                  # run fitted_parameters.json: print scores, write figures
  python run_phenology.py --fit            # refit every site first
  python run_phenology.py --fit bci        # refit one site
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
DRIVERS = os.path.join(HERE, "drivers")
FIT_FILE = os.path.join(HERE, "fitted_parameters.json")
Cue = pheno.Cue


def data(name):
    return os.path.join(DATA, name)


def driver(name):
    return os.path.join(DRIVERS, name)


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


def plot_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


class LightHours:
    """The hours a day of PAR above any par_min. Each column of drivers/era5_par_hours_<site>.csv.gz
    counts the 900 s steps of a UTC day above one value of a grid; between grid values the hours are
    interpolated in log par_min."""

    def __init__(self, key, dates):
        t = pd.read_csv(driver(f"era5_par_hours_{key}.csv.gz"), comment="#", parse_dates=["date"])
        t = t.set_index("date").reindex(dates)
        if t.isna().any().any():
            raise ValueError(f"{key}: the hours of light do not cover {dates[0]:%Y-%m-%d} - {dates[-1]:%Y-%m-%d}")
        self.log_grid = np.log([float(c.rsplit("_", 1)[1]) for c in t.columns])
        self.hours = t.values * 0.25

    def __call__(self, par_min):
        x = np.clip(np.log(max(par_min, 1e-9)), self.log_grid[0], self.log_grid[-1])
        i = min(int(np.searchsorted(self.log_grid, x, side="right")) - 1, len(self.log_grid) - 2)
        w = (x - self.log_grid[i]) / (self.log_grid[i + 1] - self.log_grid[i])
        return (1.0 - w) * self.hours[:, i] + w * self.hours[:, i + 1]


def collection_amounts(litter, rounds):
    """Model litter [canopy fraction] over each collection (start inclusive, end exclusive)."""
    cum = np.r_[0.0, np.cumsum(litter.values)]
    return cum[litter.index.get_indexer(rounds.end)] - cum[litter.index.get_indexer(rounds.start)]


def cumulative_errors(mod, obs, years):
    """Each calendar year's cumulative fraction of its fall at the collection ends: model minus traps.
    A single big collection cannot dominate it, and the traps' unit cancels."""
    out = []
    for y in sorted(set(years)):
        k = years == y
        cm, co = np.cumsum(mod[k]), np.cumsum(obs[k])
        if cm[-1] <= 0.0:                         # a trial with no litter that year
            return np.ones(1)
        out.append(cm / cm[-1] - co / co[-1])
    return np.concatenate(out)


class Site:
    """A site's daily drivers and observations; simulate() runs the kernel over them."""
    key = title = ""
    lat = cover0 = 0.0
    BASE = {}                                     # the cue masks
    FIXED = {}                                    # parameters held, not fitted
    dates = tair = psi = light = None

    def params(self, values):
        return pheno.Params(**{**self.BASE, **self.FIXED, **values})

    def simulate(self, values):
        p = self.params(values)
        ph = pheno.Phenology(p, leaf_cover=self.cover0)
        n = len(self.dates)
        uses_light = (int(p.flush_cue_mask) | int(p.shed_cue_mask)) & int(Cue.LIGHT)
        hours = self.light(p.par_min) if uses_light else np.zeros(n)
        tair = self.tair if self.tair is not None else np.full(n, 298.15)
        psi = self.psi if self.psi is not None else np.zeros(n)
        days = [ph.step(temp_day=t, par_hours=h, predawn_leaf_psi=w, doy=d)
                for t, h, w, d in zip(tair, hours, psi, self.dates.dayofyear)]
        return pd.DataFrame({"cover": [x.leaf_cover for x in days],
                             "litter": [x.senescence + x.background for x in days],
                             "flush": [x.leaf_flush_tendency for x in days],
                             "shed": [x.leaf_shed_tendency for x in days]}, index=self.dates)

    def light_mean(self, values):
        """The kernel's running mean of the hours of light, for the figures."""
        p = self.params(values)
        return pd.Series(self.light(p.par_min), index=self.dates).ewm(alpha=1.0 / max(p.light_window, 1.0)).mean()


class LitterSite(Site):
    """A site scored on litter traps: each calendar year's cumulative fall, and the canopies a year."""
    rounds = None                                 # DataFrame: start, end, days, mass (any unit)
    SCORED_YEARS = ()
    ANNUAL_LEAF_FALL = 1.0

    def annual(self, sim):
        return float(sim.litter.groupby(sim.index.year).sum().reindex(list(self.SCORED_YEARS)).mean())

    def litter_errors(self, sim):
        r = self.rounds[self.rounds.end.dt.year.isin(self.SCORED_YEARS)]
        return cumulative_errors(collection_amounts(sim.litter, r), r.mass.values, r.end.dt.year.values)

    def shares(self, sim):
        """Litter rate per collection over its mean: model and traps."""
        mod = collection_amounts(sim.litter, self.rounds) / self.rounds.days.values
        obs = (self.rounds.mass / self.rounds.days).values
        with np.errstate(invalid="ignore", divide="ignore"):
            return mod / mod.mean(), obs / obs.mean()

    def litter_scores(self, sim):
        mod, obs = self.shares(sim)
        return {"RMSE cumulative litter": float(np.sqrt(np.mean(self.litter_errors(sim) ** 2))),
                "litter share r": float(np.corrcoef(mod, obs)[0, 1]),
                "leaf litter [canopies/yr]": self.annual(sim)}


def plot_collections(ax_t, ax_m, site, sim, what, obs_label):
    """Collection-by-collection litter rates over their mean, and their mean by month."""
    mod, obs = site.shares(sim)
    r = site.rounds
    for k in range(len(r)):
        span = [r.start.iloc[k], r.end.iloc[k]]
        ax_t.plot(span, [obs[k]] * 2, color="k", lw=2, label=obs_label if k == 0 else None)
        ax_t.plot(span, [mod[k]] * 2, color="tab:green", lw=2, label="model" if k == 0 else None)
    ax_t.set(ylabel=f"{what} fall rate / mean")
    ax_t.legend(fontsize=8)
    month = r.end.dt.month.values
    obs_m, mod_m = pd.Series(obs).groupby(month).mean(), pd.Series(mod).groupby(month).mean()
    ax_m.bar(obs_m.index - 0.2, obs_m.values, 0.4, color="k", label=obs_label)
    ax_m.bar(mod_m.index + 0.2, mod_m.reindex(obs_m.index).values, 0.4, color="tab:green", label="model")
    ax_m.set(xlabel="month the collection ended", ylabel=f"{what} fall rate / mean", xticks=range(1, 13))
    ax_m.legend(fontsize=8)


def temperate_free():
    """The temperature and rate parameters fitted at Harvard Forest."""
    return {"flush_degree_days": (5.0, 1000.0),       # [K day] warmth centre
            "shed_base_temp": (278.0, 300.0),         # [K]
            "shed_degree_days": (5.0, 600.0),         # [K day] cold centre
            "flush_rate_max": (1 / 60, 1 / 3),        # [1/day]
            "shed_rate_max": (0.0005, 1 / 3),         # [1/day]
            "leaf_turnover_rate": (0.0, 0.6)}         # [1/yr]


class HarvardForest(Site):
    key, title = "harvard_forest", "Harvard Forest (deciduous broadleaf, 42.5° N)"
    lat, cover0 = 42.538, 0.0
    first, last = 2003, 2023                          # 2003 grows the first canopy (spin-up)
    BASE = dict(flush_cue_mask=Cue.TEMP | Cue.LIGHT, shed_cue_mask=Cue.TEMP | Cue.LIGHT)
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
        self.light = LightHours(self.key, self.dates)
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

    def free(self):
        """Fitted parameters and bounds. par_min stays low (the hours of light ~ the day length);
        the flush gate's sharpness is fitted (a gradual gate is still partly open in October and
        refills the canopy while it senesces)."""
        longest = pheno.daylength(self.lat, 172)
        return {**temperate_free(),
                "log10_par_min": (0.0, 2.0),                     # 1-100 umol/m2/s
                "flush_light_hours": (8.0, longest - 1.0),       # [h/day]
                "flush_light_sharpness": (0.5, 8.0),             # [1/h]
                "shed_light_hours": (6.0, longest),              # [h/day]
                "light_window": (1.0, 15.0)}                     # [day]

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

    def plot(self, sim, values, path):
        plt = plot_setup()
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
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
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


class Hyytiala(LitterSite):
    key, title = "hyytiala", "Hyytiälä (Scots pine, 61.8° N)"
    lat, cover0 = 61.8475, 1.0
    start, end = "2018-01-01", "2024-07-31"
    BASE = dict(flush_cue_mask=Cue.TEMP | Cue.LIGHT, shed_cue_mask=Cue.TEMP | Cue.LIGHT)
    #----- Held, not fitted: the needle-fall timing does not constrain the flush (its loss moves
    #      by less than the fit's noise across their whole ranges). -----------------------------#
    FIXED = dict(flush_degree_days=88.0, flush_rate_max=0.2866)
    SCORED_YEARS = (2019, 2020, 2021, 2022, 2023)
    #----- Needle fall per year as a share of the canopy. The traps give the timing, not the share:
    #      southern-Finnish Scots pine keeps 3.4-4.2 needle cohorts and needles live about three
    #      years (Pensa & Jalkanen 1999, Silva Fennica 33:654), so about 0.3 of it falls each year.
    ANNUAL_LEAF_FALL = 0.30

    def __init__(self):
        met = pd.read_csv(data("hyytiala_fluxnet_dd.csv"), usecols=["TIMESTAMP", "TA_F"], na_values=[-9999])
        met["date"] = pd.to_datetime(met.TIMESTAMP.astype(str), format="%Y%m%d")
        tair = met.set_index("date").loc[self.start:self.end].TA_F.interpolate()
        self.dates, self.tair = tair.index, tair.values + 273.15
        self.light = LightHours(self.key, self.dates)
        lit = self.needle_litter(data("hyytiala_ancillary.csv"))
        first = self.dates[0] + pd.Timedelta(days=200)              # after a 200-day spin-up
        self.rounds = lit[(lit.start >= first) & (lit.end <= self.dates[-1])].reset_index(drop=True)

    def free(self):
        """Bounded so that no single cloudy spell can trigger needle fall: a light window of at
        least a week and gentle light switches."""
        return {"shed_base_temp": (278.0, 300.0), "shed_degree_days": (5.0, 600.0),
                "shed_rate_max": (0.0005, 1 / 3), "leaf_turnover_rate": (0.0, 0.6),
                "min_leaf_cover": (0.4, 0.95),
                "log10_par_min": (1.7, 2.78),                    # 50-600 umol/m2/s
                "flush_light_hours": (0.0, 16.0), "flush_light_sharpness": (0.5, 2.0),
                "shed_light_hours": (0.0, 16.0), "shed_light_sharpness": (-2.0, -0.25),
                "light_window": (7.0, 60.0)}

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

    def loss(self, sim):
        return float(np.mean(self.litter_errors(sim) ** 2)
                     + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2)

    def scores(self, sim):
        out = self.litter_scores(sim)
        out["lowest leaf cover"] = float(sim.cover.min())
        return out

    def plot(self, sim, values, path):
        plt = plot_setup()
        fig, axs = plt.subplots(3, 1, figsize=(11, 9), constrained_layout=True)
        ax = axs[0]
        ax.plot(sim.index, sim.cover, color="tab:green", lw=2, label="model leaf cover")
        ax.plot(sim.index, sim.flush, color="tab:blue", lw=0.8, label="flush tendency")
        ax.plot(sim.index, sim.shed, color="tab:red", lw=0.8, ls="--", label="senescence tendency")
        ax.set(ylabel="fraction", ylim=(-0.05, 1.1), title=f"(a) {self.title}")
        ax.legend(ncol=3, loc="lower left", fontsize=8)
        ax2 = ax.twinx()
        ax2.plot(self.dates, self.light_mean(values), color="tab:orange", lw=1)
        ax2.set_ylabel(f"hours above {self.params(values).par_min:.0f} µmol m$^{{-2}}$ s$^{{-1}}$",
                       color="tab:orange")
        plot_collections(axs[1], axs[2], self, sim, "needle", "traps (ICOS)")
        axs[1].set_title("(b) needle fall, collection by collection")
        axs[2].set_title("(c) mean by month")
        fig.savefig(path, dpi=130)


class BCI(LitterSite):
    """A light-driven leaf exchanger: it flushes while its water threshold is never reached and
    senesces on many bright hours, so its leaves turn over in the dry season and the canopy stays
    full."""
    key, title = "bci", "Barro Colorado Island (moist tropical forest, 9.2° N)"
    lat, cover0 = 9.15, 1.0
    start, end = "2012-07-03", "2017-08-31"
    BASE = dict(flush_cue_mask=Cue.WATER, shed_cue_mask=Cue.LIGHT | Cue.WATER)
    #----- A water threshold far below the surrogate's range (-0.05 to -0.9 MPa): the exchanger
    #      never sheds for drought, and its water switch keeps it flushing. --------------------#
    FIXED = dict(leaf_psi_tlp=-1.5, flush_water_sum=3.0, flush_water_sharpness=6.667,
                 shed_water_sum=1.0, shed_water_sharpness=20.0)
    SCORED_YEARS = (2013, 2014, 2015, 2016)
    #----- Leaf litter per year as a share of the canopy: the leaves falling at BCI each year have
    #      about the canopy's area, 7.3 m2 per m2 of ground (Leigh 1999, ORNL NPP data set BRR).
    ANNUAL_LEAF_FALL = 1.0

    def __init__(self):
        tower = pd.read_csv(data("BCI_v5.1.csv"), parse_dates=["date"], na_values=["NaN"]).set_index("date")
        daily = tower[["tair", "SWC"]].resample("D").mean().interpolate()
        self.soil_psi_fit(daily.SWC)
        daily = daily.loc[self.start:self.end]
        self.dates = daily.index
        self.tair = daily.tair.values + 273.15
        self.psi = self.soil_psi(daily.SWC.values)
        self.light = LightHours(self.key, self.dates)
        #----- GLiMP control plots: mean fine litter per trap and collection. -------------------#
        lit = pd.read_csv(data("glimp_litterfall.csv"), encoding="utf-8-sig")
        lit = lit[lit.treat == "CT"]
        lit["date"] = pd.to_datetime(lit.date, dayfirst=True)
        c = lit.groupby("date").agg(mass=("mass.g.trap", "mean"), days=("days.accum", "median"))
        c["end"] = c.index
        c["start"] = c.end - pd.to_timedelta(c.days, unit="D")
        first = self.dates[0] + pd.Timedelta(days=180)              # after a half-year spin-up
        self.rounds = c[(c.start >= first) & (c.end <= self.dates[-1])].reset_index(drop=True)

    def soil_psi_fit(self, swc):
        """A Campbell retention curve, psi = -exp(a) theta_g^b, fitted to the soil water content
        and potential Kupers et al. (2019) sampled together; the tower's volumetric water content
        converts to gravimetric by its ratio to the plot's on their four sampling periods."""
        pairs = pd.concat([pd.read_csv(data(f), sep="\t") for f in
                           ("kupers_soil_moisture_mapping.txt", "kupers_soil_moisture_small_scale.txt")])
        pairs = pairs[(pairs.swp < -0.05) & (pairs.swc > 5.0)]      # the WP4C reads ~0 when wet
        self.campbell_b, self.campbell_a = np.polyfit(np.log(pairs.swc / 100.0), np.log(-pairs.swp), 1)
        plot = pd.read_csv(data("kupers_soil_moisture_mapping.txt"), sep="\t")
        plot["date"] = pd.to_datetime(plot.date, format="%m/%d/%Y")
        periods = plot.groupby("period").agg(date=("date", "median"), swc=("swc", "median"))
        tower = [swc.loc[t - pd.Timedelta(days=3): t + pd.Timedelta(days=3)].mean() for t in periods.date]
        self.bulk_density = float(np.median(np.array(tower) / (periods.swc.values / 100.0)))

    def soil_psi(self, swc):
        theta_g = swc / self.bulk_density
        return -np.exp(self.campbell_a + self.campbell_b * np.log(theta_g))

    def free(self):
        """Rates, the evergreen floor and the light trigger. par_min is held to bright light:
        free, it drifts to the photoperiod, a calendar with no year-to-year signal."""
        return {"flush_rate_max": (1 / 60, 1 / 3), "shed_rate_max": (0.0005, 1 / 3),
                "leaf_turnover_rate": (0.0, 3.0), "min_leaf_cover": (0.8, 0.95),
                "log10_par_min": (2.0, 3.08),                    # 100-1200 umol/m2/s
                "shed_light_hours": (0.0, 13.0),                 # [h/day]
                "shed_light_sharpness": (0.1, 4.0),              # > 0: many bright hours trigger senescence
                "light_window": (2.0, 60.0)}                     # [day]

    def loss(self, sim):
        return float(np.mean(self.litter_errors(sim) ** 2)
                     + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2)

    def scores(self, sim):
        out = self.litter_scores(sim)
        out["lowest leaf cover"] = float(sim.cover.min())
        mod, obs = self.shares(sim)
        days, dry = self.rounds.days.values, self.rounds.end.dt.month.isin([1, 2, 3, 4]).values
        out["Jan-Apr share of the leaf fall, model"] = float((mod * days)[dry].sum() / (mod * days).sum())
        out["Jan-Apr share of the leaf fall, traps"] = float((obs * days)[dry].sum() / (obs * days).sum())
        return out

    def plot(self, sim, values, path):
        plt = plot_setup()
        p = self.params(values)
        fig, axs = plt.subplots(4, 1, figsize=(11, 11), constrained_layout=True)
        ax = axs[0]
        ax.plot(self.dates, self.psi, color="tab:blue", lw=1, label="soil water potential (surrogate)")
        ax.set(ylabel="MPa", title=f"(a) {self.title}: drivers", ylim=(-1.0, 0.05))
        ax2 = ax.twinx()
        ax2.plot(self.dates, self.light_mean(values), color="tab:orange", lw=1,
                 label=f"hours above {p.par_min:.0f} µmol m$^{{-2}}$ s$^{{-1}}$, running mean")
        ax2.axhline(p.shed_light_hours, color="tab:orange", ls=":", lw=1)
        ax2.set_ylabel("h day$^{-1}$")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=8, ncol=2)
        ax = axs[1]
        ax.plot(sim.index, sim.cover, color="tab:green", lw=2, label="leaf cover")
        ax.plot(sim.index, sim.shed, color="tab:red", lw=0.8, ls="--", label="senescence tendency")
        ax.set(ylabel="fraction", ylim=(-0.05, 1.1), title="(b) the light exchanger")
        ax.legend(ncol=2, loc="lower left", fontsize=8)
        plot_collections(axs[2], axs[3], self, sim, "litter", "traps (GLiMP)")
        axs[2].set_title("(c) litter, collection by collection")
        axs[3].set_title("(d) mean by month")
        fig.savefig(path, dpi=130)


class PaloVerde(LitterSite):
    """A drought-deciduous forest: flushing and senescence on the predawn leaf water potential and
    the hours of light. The water potential is hypothetical: a MEDS run's canopy predawn value for
    an evergreen stand at the site, so the fit asks what cues map that driver onto the observed
    canopy, not how this forest's own water status moves."""
    key, title = "palo_verde", "Palo Verde (seasonally dry tropical forest, 10.4° N)"
    lat, cover0 = 10.35, 1.0
    start, end = "2008-01-01", "2013-12-31"
    BASE = dict(flush_cue_mask=Cue.WATER | Cue.LIGHT, shed_cue_mask=Cue.WATER | Cue.LIGHT)
    SCORED_YEARS = (2009, 2010, 2011, 2012, 2013)
    #----- About one canopy a year: 0.56 kg m-2 yr-1 of leaf litter over a canopy of LAI 5.5 at
    #      80-100 g m-2 of leaf is 1.0-1.3 canopies. Without the term the fit churns leaves. -----#
    ANNUAL_LEAF_FALL = 1.0

    def __init__(self):
        d = pd.read_csv(driver("palo_verde_predawn_psi.csv"), comment="#", parse_dates=["date"]).set_index("date")
        d = d.loc[self.start:self.end]
        self.dates, self.psi = d.index, d.predawn_leaf_psi.values
        self.light = LightHours(self.key, self.dates)
        #----- MODIS LAI: the monthly maximum composite (Xu et al. 2016), over its 95th percentile. -#
        m = pd.read_csv(driver("palo_verde_modis_lai.csv"), comment="#", parse_dates=["date"]).set_index("date")
        lai = m.median_mainalg_clear.fillna(m.median_all_5x5).resample("MS").max()
        self.lai_full = float(lai.quantile(0.95))
        self.lai_rel = (lai / self.lai_full).clip(upper=1.1).loc[f"{self.SCORED_YEARS[0]}":f"{self.SCORED_YEARS[-1]}"]
        #----- Litter traps: g per 0.25 m2 trap, about monthly (drivers/palo_verde_leaf_litter.csv). -#
        r = pd.read_csv(driver("palo_verde_leaf_litter.csv"), comment="#", parse_dates=["start", "end"])
        r["days"] = (r.end - r.start).dt.days
        r["mass"] = r.leaf_g_m2_day * r.days
        first = self.dates[0] + pd.Timedelta(days=150)              # after a five-month spin-up
        self.rounds = r[(r.start >= first) & (r.end <= self.dates[-1])].reset_index(drop=True)

    def free(self):
        """Water and light on both sides, with both light signs free."""
        return {"leaf_psi_tlp": (-3.0, -0.3),                    # [MPa]
                "flush_water_sum": (0.5, 60.0), "flush_water_sharpness": (0.5, 20.0),
                "shed_water_sum": (0.5, 60.0), "shed_water_sharpness": (0.5, 20.0),
                "flush_rate_max": (1 / 60, 1 / 3), "shed_rate_max": (0.0005, 1 / 3),
                "leaf_turnover_rate": (0.0, 1.0), "min_leaf_cover": (0.0, 0.5),
                "log10_par_min": (0.0, 3.08),                    # 1-1200 umol/m2/s
                "flush_light_hours": (0.0, 13.0), "flush_light_sharpness": (-8.0, 8.0),
                "shed_light_hours": (0.0, 13.0), "shed_light_sharpness": (-8.0, 8.0),
                "light_window": (1.0, 60.0)}

    def lai_errors(self, sim):
        cm = sim.cover.resample("MS").mean().reindex(self.lai_rel.index)
        return cm.values - self.lai_rel.values

    def loss(self, sim):
        return float(np.nanmean(self.lai_errors(sim) ** 2) + np.mean(self.litter_errors(sim) ** 2)
                     + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2)

    def leafless(self, sim):
        """Per year, the days the canopy is below half cover."""
        return {y: int((g < 0.5).sum()) for y, g in sim.cover.groupby(sim.index.year)}

    def scores(self, sim):
        out = {"RMSE relative LAI": float(np.sqrt(np.nanmean(self.lai_errors(sim) ** 2)))}
        out.update(self.litter_scores(sim))
        out["lowest leaf cover"] = float(sim.cover.min())
        return out

    def plot(self, sim, values, path):
        plt = plot_setup()
        p = self.params(values)
        fig, axs = plt.subplots(4, 1, figsize=(11, 12), constrained_layout=True)
        ax = axs[0]
        ax.plot(self.dates, self.psi, color="tab:blue", lw=1, label="canopy predawn leaf ψ (MEDS run)")
        ax.axhline(p.leaf_psi_tlp, color="tab:blue", ls=":", lw=1)
        ax.set(ylabel="MPa", title=f"(a) {self.title}: drivers")
        ax2 = ax.twinx()
        ax2.plot(self.dates, self.light_mean(values), color="tab:orange", lw=1,
                 label=f"hours above {p.par_min:.0f} µmol m$^{{-2}}$ s$^{{-1}}$, running mean")
        ax2.set_ylabel("h day$^{-1}$")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=8)
        ax = axs[1]
        ax.plot(self.lai_rel.index + pd.Timedelta(days=14), self.lai_rel.values, "ko", ms=4,
                label=f"MODIS LAI / {self.lai_full:.1f} (monthly maximum)")
        ax.plot(sim.index, sim.cover, color="tab:green", lw=1.5, label="model leaf cover")
        ax.plot(sim.index, sim.shed, color="tab:red", lw=0.8, ls="--", label="senescence tendency")
        ax.set(ylabel="fraction", ylim=(-0.05, 1.15), title="(b) the canopy",
               xlim=(pd.Timestamp("2009-01-01"), pd.Timestamp("2013-12-31")))
        ax.legend(ncol=3, loc="lower left", fontsize=8)
        plot_collections(axs[2], axs[3], self, sim, "leaf", "traps (Xu et al. 2016)")
        axs[2].set_title("(c) leaf fall, collection by collection")
        axs[3].set_title("(d) mean by month")
        fig.savefig(path, dpi=130)


SITES = {cls.key: cls for cls in (HarvardForest, Hyytiala, BCI, PaloVerde)}
_LOADED = {}


def site(key):
    """Each site's data, loaded once per process (also in each fitting worker)."""
    if key not in _LOADED:
        _LOADED[key] = SITES[key]()
    return _LOADED[key]


def as_values(names, x):
    """Fitted vector -> parameter values; par_min is fitted as its log10."""
    v = dict(zip(names, map(float, x)))
    if "log10_par_min" in v:
        v["par_min"] = 10.0 ** v.pop("log10_par_min")
    return v


def _loss(x, key, names):
    s = site(key)
    try:
        value = s.loss(s.simulate(as_values(names, x)))
    except (ValueError, FloatingPointError):
        return 1e3
    return float(value) if np.isfinite(value) else 1e3


def fit(key, workers, seeds=4):
    """Differential evolution from several seeds, each polished by a bounded Nelder-Mead; the best
    wins, and the spread of the seeds is reported (a wide spread means the optimum is not unique)."""
    from scipy.optimize import differential_evolution, minimize
    free = site(key).free()
    names, bounds, runs = list(free), list(free.values()), []
    for seed in range(1, seeds + 1):
        r = differential_evolution(_loss, bounds, args=(key, names), popsize=25, maxiter=250, tol=1e-8,
                                   seed=seed, polish=False, workers=workers, updating="deferred")
        q = minimize(_loss, r.x, args=(key, names), method="Nelder-Mead", bounds=bounds,
                     options=dict(maxiter=2000, xatol=1e-5, fatol=1e-9))
        x, v = (q.x, q.fun) if q.fun < r.fun else (r.x, r.fun)
        runs.append((float(v), x))
        print(f"  {key} seed {seed}: {r.fun:.5f} -> polished {v:.5f}", flush=True)
    runs.sort(key=lambda t: t[0])
    print(f"  {key}: best {runs[0][0]:.5f}; seeds {', '.join(f'{v:.5f}' for v, _ in runs)}", flush=True)
    return {n: float(f"{v:.6g}") for n, v in as_values(names, runs[0][1]).items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fit", nargs="*", metavar="SITE",
                    help="refit these sites (all if none given) and rewrite " + os.path.basename(FIT_FILE))
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="processes for the fit")
    args = ap.parse_args(argv)
    with open(FIT_FILE) as fh:
        values = json.load(fh)
    if args.fit is not None:
        for key in args.fit or list(SITES):
            values[key] = fit(key, args.workers)
            with open(FIT_FILE) as fh:                       # other sites may have been refitted meanwhile
                current = json.load(fh)
            current[key] = values[key]
            with open(FIT_FILE, "w") as fh:
                json.dump(current, fh, indent=2)
                fh.write("\n")
            values = current
    for key in SITES:
        if key not in values:
            continue
        s = site(key)
        sim = s.simulate(values[key])
        print(f"\n{s.title}")
        for name, value in {**s.FIXED, **values[key]}.items():
            print(f"  {name:26s} {value:g}")
        for name, value in s.scores(sim).items():
            print(f"  {name:46s} {value:.3f}")
        if key == "bci":
            print(f"  soil water potential = -exp({s.campbell_a:.2f}) (SWC / {s.bulk_density:.3f})^{s.campbell_b:.2f} MPa")
        if key == "palo_verde":
            print("  days below half cover: " + ", ".join(f"{y} {n}" for y, n in s.leafless(sim).items()))
        s.plot(sim, values[key], os.path.join(HERE, f"{key}.png"))
    print("\nwrote " + ", ".join(f"{k}.png" for k in SITES if k in values))


if __name__ == "__main__":
    main()
