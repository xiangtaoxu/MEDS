#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 02: one MEDS phenology kernel at three forests and four leaf habits.

  Harvard Forest  deciduous broadleaf: warmth and day length
  Hyytiala        Scots pine (evergreen): warmth and PAR, senescence stopping at a leaf-cover floor
  BCI             two tropical species under one climate, both cued on predawn water potential and
                  PAR: a drought-deciduous species and a light-driven leaf exchanger

The kernel and the leaf rule are the compiled Fortran of the coupled model, reached through
meds.plant.pheno, with carbon never limiting the flush. At BCI the predawn water potential is a
surrogate: the tower's soil water content through a retention curve fitted to paired soil samples.

Usage (after fetch_phenology_data.py; build libmeds.so first, see the README):
  python run_phenology.py           # run with fitted_parameters.json: print scores, write figures
  python run_phenology.py --fit     # refit the fitted parameters first (differential evolution)
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
Cue = pheno.Cue


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


def collection_share(litter, coll):
    """Model litter rate over each trap collection, over its mean (the traps' unit is relative)."""
    cum, day = litter.cumsum(), pd.Timedelta(days=1)
    rate = ((cum.reindex(coll.index - day, method="nearest").values
             - cum.reindex(coll.start - day, method="nearest").values) / coll.days.values)
    with np.errstate(invalid="ignore", divide="ignore"):    # a trial with no litter scores NaN
        return rate / rate.mean()


def plot_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


class Site:
    """A site's daily drivers and observations; simulate() runs the kernel over them."""
    key = title = ""
    lat = cover0 = 0.0
    BASE = {}                                   # the cue masks and any fixed parameters
    dates = tair = par = psi = None

    def daylengths(self):
        if not hasattr(self, "_daylen"):
            self._daylen = np.array([pheno.daylength(self.lat, d) for d in self.dates.dayofyear])
        return self._daylen

    def simulate(self, values, base=None):
        ph = pheno.Phenology(pheno.Params(**(self.BASE if base is None else base), **values),
                             leaf_cover=self.cover0)
        n = len(self.dates)
        par = self.par if self.par is not None else np.zeros(n)
        psi = self.psi if self.psi is not None else np.zeros(n)
        days = [ph.step(temp_day=t, daylength=dl, par=p, predawn_leaf_psi=w, doy=d)
                for t, dl, p, w, d in zip(self.tair, self.daylengths(), par, psi, self.dates.dayofyear)]
        return pd.DataFrame({"cover": [x.leaf_cover for x in days],
                             "litter": [x.senescence + x.background for x in days],
                             "flush": [x.leaf_flush_tendency for x in days],
                             "shed": [x.leaf_shed_tendency for x in days]}, index=self.dates)

    def temperate_free(self):
        """The temperature and rate parameters fitted at both temperate sites."""
        return {"flush_degree_days": (5.0, 1000.0),       # [K day] warmth centre
                "shed_base_temp": (278.0, 300.0),         # [K]
                "shed_degree_days": (5.0, 600.0),         # [K day] cold centre
                "flush_rate_max": (1 / 60, 1 / 3),        # [1/day]
                "shed_rate_max": (0.0005, 1 / 3),         # [1/day]
                "leaf_turnover_rate": (0.0, 0.6)}         # [1/yr]


class HarvardForest(Site):
    key, title = "harvard_forest", "Harvard Forest (deciduous broadleaf, 42.5° N)"
    lat, cover0 = 42.538, 0.0
    first, last = 2002, 2023                          # 2002 grows the first canopy (spin-up)
    BASE = dict(flush_cue_mask=Cue.TEMP | Cue.DAYLENGTH, shed_cue_mask=Cue.TEMP | Cue.DAYLENGTH)
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

    def free(self):
        """Fitted parameters and bounds. The flush day length stays an hour below the longest
        day, so the gate opens in summer; its sharpness is fitted too (a gradual gate is still
        partly open in October and refills the canopy while it senesces)."""
        longest = pheno.daylength(self.lat, 172)
        return {**self.temperate_free(),
                "flush_daylength_threshold": (8.0, longest - 1.0),   # [h]
                "flush_daylength_sharpness": (0.5, 8.0),             # [1/h]
                "shed_daylength_threshold": (8.0, longest)}          # [h]

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
    BASE = dict(flush_cue_mask=Cue.TEMP | Cue.PAR, shed_cue_mask=Cue.TEMP | Cue.PAR)
    #----- Needle fall per year as a share of the canopy. The traps give the timing, not the share:
    #      southern-Finnish Scots pine keeps 3.4-4.2 needle cohorts and needles live about three
    #      years (Pensa & Jalkanen 1999, Silva Fennica 33:654), so about 0.3 of it falls each year.
    ANNUAL_NEEDLE_FALL = 0.30

    def __init__(self):
        met = pd.read_csv(data("hyytiala_fluxnet_dd.csv"), usecols=["TIMESTAMP", "TA_F", "PPFD_IN", "SW_IN_F"],
                          na_values=[-9999])
        met["date"] = pd.to_datetime(met.TIMESTAMP.astype(str), format="%Y%m%d")
        met = met.set_index("date").loc[self.start:self.end]
        #----- PAR [umol/m2/s, daily mean]: the measured PPFD, its gaps filled from the shortwave
        #      by their ratio over the days with both. ------------------------------------------#
        ratio = float((met.PPFD_IN / met.SW_IN_F).median())
        par = met.PPFD_IN.where(met.PPFD_IN.notna(), ratio * met.SW_IN_F).interpolate()
        tair = met.TA_F.interpolate()
        self.dates, self.tair, self.par = tair.index, tair.values + 273.15, par.values
        lit = self.needle_litter(data("hyytiala_ancillary.csv"))
        first = self.dates[0] + pd.Timedelta(days=200)              # after a 200-day spin-up
        lit = lit[(lit.start >= first) & (lit.end <= self.dates[-1])].set_index("end")
        self.litter = lit
        rate = (lit.mass / lit.days).values
        self.share = rate / rate.mean()          # needle fall rate per interval / its mean

    def free(self):
        return {**self.temperate_free(),
                "min_leaf_cover": (0.4, 0.95),
                "flush_par_threshold": (20.0, 700.0),     # [umol/m2/s]
                "flush_par_sharpness": (0.002, 0.5),      # [m2 s/umol]
                "shed_par_threshold": (20.0, 700.0),      # [umol/m2/s]
                "shed_par_sharpness": (-0.5, -0.002),     # < 0: dim light triggers senescence
                "par_window": (2.0, 60.0)}                # [day]

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

    def annual(self, sim):
        return float(sim.litter.loc["2019":"2023"].sum() / 5.0)

    def loss(self, sim):
        rmse = np.sqrt(np.mean((collection_share(sim.litter, self.litter) - self.share) ** 2))
        return rmse ** 2 + 4.0 * (self.annual(sim) - self.ANNUAL_NEEDLE_FALL) ** 2

    def scores(self, sim):
        mod = collection_share(sim.litter, self.litter)
        return {"RMSE needle-fall share": float(np.sqrt(np.mean((mod - self.share) ** 2))),
                "needle-fall share r": float(np.corrcoef(mod, self.share)[0, 1]),
                "needle fall [canopies/yr]": self.annual(sim),
                "lowest leaf cover": float(sim.cover.min())}

    def plot(self, sim, path):
        plt = plot_setup()
        fig, axs = plt.subplots(3, 1, figsize=(11, 8.5), constrained_layout=True)
        ax = axs[0]
        ax.plot(sim.index, sim.cover, color="tab:green", lw=2, label="model leaf cover")
        ax.plot(sim.index, sim.flush, color="tab:blue", lw=0.8, label="flush tendency")
        ax.plot(sim.index, sim.shed, color="tab:red", lw=0.8, ls="--", label="senescence tendency")
        ax.set(ylabel="fraction", ylim=(-0.05, 1.1), title=f"(a) {self.title}")
        ax.legend(ncol=3, loc="lower left", fontsize=8)
        plot_collections(axs[1], axs[2], self.litter, self.share, collection_share(sim.litter, self.litter),
                         "needle", "traps (ICOS)")
        fig.savefig(path, dpi=130)


def plot_collections(ax_t, ax_m, coll, obs, mod, what, obs_label):
    """Collection-by-collection litter rates over their mean, and their mean by month."""
    for k in range(len(coll)):
        span = [coll.start.iloc[k], coll.index[k]]
        ax_t.plot(span, [obs[k]] * 2, color="k", lw=2, label=obs_label if k == 0 else None)
        ax_t.plot(span, [mod[k]] * 2, color="tab:green", lw=2, label="model" if k == 0 else None)
    ax_t.set(ylabel=f"{what} fall rate / mean", title=f"(b) {what} fall, collection by collection")
    ax_t.legend(fontsize=8)
    month = coll.index.month
    obs_m = pd.Series(obs).groupby(month).mean()
    mod_m = pd.Series(mod).groupby(month).mean()
    ax_m.bar(obs_m.index - 0.2, obs_m.values, 0.4, color="k", label=obs_label)
    ax_m.bar(mod_m.index + 0.2, mod_m.reindex(obs_m.index).values, 0.4, color="tab:green", label="model")
    ax_m.set(xlabel="month the collection ended", ylabel=f"{what} fall rate / mean",
             xticks=range(1, 13), title="(c) mean by month")
    ax_m.legend(fontsize=8)


class BCI(Site):
    """Two species under one climate. Both flush on water and senesce on water or PAR; they
    differ in their water threshold (the predawn potential they shed at) and PAR sensitivity."""
    key, title = "bci", "Barro Colorado Island (moist tropical forest, 9.2° N)"
    lat, cover0 = 9.15, 1.0
    start, end = "2012-07-03", "2017-08-31"
    BASE = dict(flush_cue_mask=Cue.WATER, shed_cue_mask=Cue.PAR | Cue.WATER)
    FITTED = "light_exchanger"
    #----- Leaf litter per year as a share of the canopy: the leaves falling at BCI each year have
    #      about the canopy's area, 7.3 m2 per m2 of ground (Leigh 1999, ORNL NPP data set BRR).
    ANNUAL_LEAF_FALL = 1.0

    def __init__(self):
        tower = pd.read_csv(data("BCI_v5.1.csv"), parse_dates=["date"], na_values=["NaN"]).set_index("date")
        #----- PAR [umol/m2/s]: the tower's shortwave times the PAR/shortwave ratio over the
        #      half-hours with both (its PAR sensor covers only part of the record). --------------#
        both = tower[["Rs", "Par_tot"]].dropna()
        both = both[both.Rs > 50.0]
        self.par_per_sw = float((both.Par_tot / both.Rs).median())
        daily = tower[["tair", "Rs", "SWC"]].resample("D").mean().interpolate()
        self.soil_psi_fit(daily.SWC)
        daily = daily.loc[self.start:self.end]
        self.dates = daily.index
        self.tair = daily.tair.values + 273.15
        self.par = self.par_per_sw * daily.Rs.values
        self.psi = self.soil_psi(daily.SWC.values)
        #----- GLiMP control plots: mean fine litter per trap and collection. -------------------#
        lit = pd.read_csv(data("glimp_litterfall.csv"), encoding="utf-8-sig")
        lit = lit[lit.treat == "CT"]
        lit["date"] = pd.to_datetime(lit.date, dayfirst=True)
        c = lit.groupby("date").agg(mass=("mass.g.trap", "mean"), days=("days.accum", "median"))
        c["start"] = c.index - pd.to_timedelta(c.days, unit="D")
        first = self.dates[0] + pd.Timedelta(days=180)              # after a half-year spin-up
        self.litter = c[(c.start >= first) & (c.index <= self.dates[-1])]
        rate = (self.litter.mass / self.litter.days).values
        self.share = rate / rate.mean()

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
        """The light exchanger: rates, its evergreen floor and its PAR trigger."""
        return {"flush_rate_max": (1 / 60, 1 / 3), "shed_rate_max": (0.0005, 1 / 3),
                "leaf_turnover_rate": (0.0, 3.0), "min_leaf_cover": (0.8, 0.95),
                "shed_par_threshold": (100.0, 900.0),    # [umol/m2/s]
                "shed_par_sharpness": (0.002, 0.5),      # > 0: bright light triggers senescence
                "par_window": (2.0, 60.0)}               # [day]

    def annual(self, sim):
        return float(sim.litter.loc["2013":"2016"].sum() / 4.0)

    def loss(self, sim):
        rmse = np.sqrt(np.mean((collection_share(sim.litter, self.litter) - self.share) ** 2))
        return rmse ** 2 + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2

    def scores(self, sim):
        mod = collection_share(sim.litter, self.litter)
        return {"RMSE litter share": float(np.sqrt(np.mean((mod - self.share) ** 2))),
                "litter share r": float(np.corrcoef(mod, self.share)[0, 1]),
                "leaf litter [canopies/yr]": self.annual(sim),
                "lowest leaf cover": float(sim.cover.min())}

    def leafless(self, sim):
        """Per dry season, the days the deciduous canopy is below half cover, and the dates."""
        out = {}
        for year, g in sim.cover.groupby(sim.index.year):
            low = g[g < 0.5]
            if len(low):
                out[year] = (len(low), low.index[0], low.index[-1])
        return out

    def plot(self, sims, path):
        plt = plot_setup()
        fig, axs = plt.subplots(4, 1, figsize=(11, 11), constrained_layout=True)
        ax = axs[0]
        ax.plot(self.dates, self.psi, color="tab:blue", lw=1, label="soil water potential (surrogate)")
        for name, colour in (("drought_deciduous", "tab:brown"), ("light_exchanger", "tab:green")):
            ax.axhline(self.values[name]["leaf_psi_tlp"], color=colour, ls=":", lw=1)
        ax.set(ylabel="MPa", title=f"(a) {self.title}: drivers", ylim=(-1.0, 0.05))
        ax2 = ax.twinx()
        parm = pd.Series(self.par, index=self.dates).ewm(alpha=1.0 / self.values["light_exchanger"]["par_window"]).mean()
        ax2.plot(self.dates, parm, color="tab:orange", lw=1, label="PAR, running mean")
        ax2.axhline(self.values["light_exchanger"]["shed_par_threshold"], color="tab:orange", ls=":", lw=1)
        ax2.set_ylabel("umol m$^{-2}$ s$^{-1}$")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=8, ncol=2)

        ax = axs[1]
        for name, colour, label in (("drought_deciduous", "tab:brown", "drought-deciduous"),
                                    ("light_exchanger", "tab:green", "light exchanger")):
            ax.plot(sims[name].index, sims[name].cover, color=colour, lw=2, label=f"{label}: leaf cover")
            ax.plot(sims[name].index, sims[name].shed, color=colour, lw=0.8, ls="--",
                    label=f"{label}: senescence tendency")
        ax.set(ylabel="fraction", ylim=(-0.05, 1.1), title="(b) two species, one climate")
        ax.legend(ncol=2, loc="lower left", fontsize=8)
        plot_collections(axs[2], axs[3], self.litter, self.share,
                         collection_share(sims["light_exchanger"].litter, self.litter), "litter", "traps (GLiMP)")
        axs[2].set_title("(c) light exchanger: litter, collection by collection")
        axs[3].set_title("(d) light exchanger: mean by month")
        fig.savefig(path, dpi=130)


SITES = {cls.key: cls for cls in (HarvardForest, Hyytiala, BCI)}
_LOADED = {}


def site(key):
    """Each site's data, loaded once per process (also in each fitting worker)."""
    if key not in _LOADED:
        _LOADED[key] = SITES[key]()
    return _LOADED[key]


def species_values(key, fitted, species=None):
    v = fitted[key]
    return v[species] if species else v


def _loss(x, key, names, fixed):
    s = site(key)
    value = s.loss(s.simulate({**fixed, **dict(zip(names, x))}))
    return float(value) if np.isfinite(value) else 1e3


def fit(key, workers, fixed):
    """Best of three differential-evolution seeds over the site's free parameters."""
    from scipy.optimize import differential_evolution
    free = site(key).free()
    names, best = list(free), None
    fixed = {k: v for k, v in fixed.items() if k not in free}
    for seed in (1, 2, 3):
        r = differential_evolution(_loss, [free[n] for n in names], args=(key, names, fixed), popsize=15,
                                   maxiter=120, tol=1e-6, seed=seed, polish=False, workers=workers,
                                   updating="deferred")
        print(f"  {key} seed {seed}: loss {r.fun:.5f}", flush=True)
        if best is None or r.fun < best.fun:
            best = r
    return {**fixed, **{n: float(f"{v:.4g}") for n, v in zip(names, best.x)}}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fit", action="store_true", help="refit and rewrite " + os.path.basename(FIT_FILE))
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="processes for the fit")
    args = ap.parse_args(argv)
    with open(FIT_FILE) as fh:
        values = json.load(fh)
    if args.fit:
        values["harvard_forest"] = fit("harvard_forest", args.workers, {})
        values["hyytiala"] = fit("hyytiala", args.workers, {})
        values["bci"][BCI.FITTED] = fit("bci", args.workers, values["bci"][BCI.FITTED])
        with open(FIT_FILE, "w") as fh:
            json.dump(values, fh, indent=2)
            fh.write("\n")
    for key in ("harvard_forest", "hyytiala"):
        s = site(key)
        sim = s.simulate(values[key])
        print(f"\n{s.title}")
        for name, value in values[key].items():
            print(f"  {name:26s} {value:g}")
        for name, value in s.scores(sim).items():
            print(f"  {name:30s} {value:.3f}")
        s.plot(sim, os.path.join(HERE, f"{key}.png"))
    s = site("bci")
    s.values = values["bci"]
    sims = {name: s.simulate(v) for name, v in values["bci"].items()}
    print(f"\n{s.title}: soil water potential = -exp({s.campbell_a:.2f}) (SWC / {s.bulk_density:.3f})"
          f"^{s.campbell_b:.2f} MPa; PAR = {s.par_per_sw:.3f} x shortwave")
    for name, v in values["bci"].items():
        print(f"  {name}:")
        for k, value in v.items():
            print(f"    {k:26s} {value:g}")
    for k, value in s.scores(sims["light_exchanger"]).items():
        print(f"  light exchanger {k:30s} {value:.3f}")
    for year, (n, a, b) in s.leafless(sims["drought_deciduous"]).items():
        print(f"  drought-deciduous {year}: below half cover {n} days, {a:%d %b} - {b:%d %b}")
    s.plot(sims, os.path.join(HERE, "bci.png"))
    print("\nwrote harvard_forest.png, hyytiala.png, bci.png")


if __name__ == "__main__":
    main()
