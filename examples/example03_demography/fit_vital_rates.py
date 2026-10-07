#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Fit the three vital-rate forests to the census tables and tabulate them for the MEDS driver.

    python fit_vital_rates.py          # reads data/ (prepare_census.py), writes vital_rates/

  growth       dbh growth [cm/yr]           from dbh, BAL, PFT
  mortality    death rate [1/yr]            from dbh, PREDICTED growth, PFT
  recruitment  new stems >= 1 cm [1/m2/yr]  from the quadrat's basal area, the PFT's share of it, PFT

Mortality learns from out-of-fold predicted growth, not measured growth: a MEDS cohort carries the
mean growth of its trees, so the law has to say how fast trees die that are EXPECTED to grow this
fast; this is also how the neighbourhood reaches mortality. Skill is cross-validated on 1-ha blocks,
so that neighbouring trees never sit on both sides of a split. Each forest is then evaluated on a
regular grid and written as a table; runs read the tables, not the forests. Needs scikit-learn;
about two minutes on 40 cores.
"""
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "vital_rates")
MIN_LEAF = {"growth": 300, "mortality": 1000, "recruitment": 50}   # trees (quadrats) per leaf
N_FOLD = 5


def forest(name):
    return RandomForestRegressor(n_estimators=200, min_samples_leaf=MIN_LEAF[name], max_features=1.0,
                                 max_samples=0.5, n_jobs=-1, random_state=0)


def class_mean_r2(obs, pred, df):
    """R2 of the means by PFT, size class and BAL class -- what a cohort carries."""
    keys = [df.pft.values, np.digitize(df.dbh.values, [2, 5, 10, 20, 40, 80]),
            np.digitize(df.bal.values, [5, 10, 20, 30, 45, 60, 80])]
    t = pd.DataFrame({"o": obs, "p": pred}).groupby(keys).agg(o=("o", "mean"), p=("p", "mean"),
                                                              n=("o", "size"))
    t = t[t.n >= 50]
    return r2_score(t.o, t.p, sample_weight=t.n)


g = pd.read_csv(os.path.join(DATA, "growth.csv.gz"))
s = pd.read_csv(os.path.join(DATA, "survival.csv.gz"))
r = pd.read_csv(os.path.join(DATA, "recruits.csv.gz"))
folds = GroupKFold(N_FOLD)

# growth, predicting each held-out block's growth rows (skill) and survival rows (mortality's input)
XG = ["dbh", "bal", "pft"]
g_cv = np.zeros(len(g))
s["growth"] = np.nan
for train, test in folds.split(g, groups=g.block):
    model = forest("growth").fit(g.loc[train, XG], g.g.iloc[train])
    g_cv[test] = model.predict(g.loc[test, XG])
    held = s.block.isin(g.block.iloc[test].unique())
    s.loc[held, "growth"] = model.predict(s.loc[held, XG])
print(f"growth       R2 {r2_score(g.g, g_cv):.2f} for trees, {class_mean_r2(g.g, g_cv, g):.3f} for class means")
growth = forest("growth").fit(g[XG], g.g)

XM = ["dbh", "growth", "pft"]
dt = s.dt.mean()
m_cv = np.zeros(len(s))
for train, test in folds.split(s, groups=s.block):
    m_cv[test] = forest("mortality").fit(s.loc[train, XM], s.dead.iloc[train]).predict(s.loc[test, XM])
print(f"mortality    AUC {roc_auc_score(s.dead, m_cv):.2f} for trees, R2 {class_mean_r2(s.dead, m_cv, s):.3f} "
      f"for class means")
mortality = forest("mortality").fit(s[XM], s.dead)

XR = ["ba_tot", "ba_pft", "pft"]
r_cv = np.zeros(len(r))
for train, test in folds.split(r, groups=r.block):
    r_cv[test] = forest("recruitment").fit(r.loc[train, XR], r.rate.iloc[train]).predict(r.loc[test, XR])
print(f"recruitment  R2 {r2_score(r.rate, r_cv):.2f} for quadrats; plot total {r.rate.mean() * 3e4:.0f} "
      f"observed, {r_cv.mean() * 3e4:.0f} predicted [stems/ha/yr]")
recruitment = forest("recruitment").fit(r[XR], r.rate)


def table(model, names, axes, value):
    """The forest on the grid spanned by ``axes``; the driver reads between the nodes."""
    mesh = np.meshgrid(*axes, indexing="ij")
    t = pd.DataFrame({n: a.ravel() for n, a in zip(names, mesh)})
    t[value] = model.predict(t[names])
    return t


pft = np.array([1, 2, 3])
dbh = np.round(np.geomspace(1.0, 150.0, 41), 4)
bal = np.arange(0.0, 112.5, 2.5)
gro = np.round(np.linspace(0.0, np.quantile(s.growth, 0.995), 41), 5)
ba = np.arange(0.0, 112.5, 5.0)
tables = {
    "growth": table(growth, XG, [dbh, bal, pft], "growth"),
    "mortality": table(mortality, XM, [dbh, gro, pft], "p_dead"),
    "recruitment": table(recruitment, XR, [ba, ba, pft], "recruitment"),
}
# the forest gives the chance of dying within a census interval; the engine wants a rate per year
m = tables["mortality"]
m["mortality"] = -np.log(1.0 - m.pop("p_dead").clip(upper=0.999)) / dt   # clamp-ok: certain death
os.makedirs(OUT, exist_ok=True)
for name, t in tables.items():
    cols = ["pft"] + [c for c in t.columns if c not in ("pft", name)] + [name]
    t[cols].sort_values(cols[:3]).to_csv(os.path.join(OUT, f"{name}.csv"), index=False, float_format="%.6g")
print(f"wrote {OUT}/growth.csv, mortality.csv, recruitment.csv")
