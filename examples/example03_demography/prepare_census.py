#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""BCI 50-ha plot censuses -> what the example fits, starts from and compares with.

    python prepare_census.py CENSUS_DIR

CENSUS_DIR holds the ForestGEO tree tables of censuses 2-7 as bci_1985.csv ... bci_2010.csv (one row
per tree; dbh in mm). Writes, in data/ (not committed):

  growth.csv.gz      one row per tree and interval: dbh growth [cm/yr] of a tree alive at both ends
  survival.csv.gz    one row per tree and interval: alive at the start, dead (1) or alive (0) at the end
  recruits.csv.gz    one row per quadrat, PFT and interval: new stems >= 1 cm [plants/m2/yr]
  bci_1985_census.csv  the 1985 stand as a MEDS census: a patch per quadrat, nplant per m2

and, committed, census_stand.csv: stems and basal area by PFT and size class in every census.

A tree is described by what it is at the start of an interval: its dbh, its PFT (from its species'
wood density) and the basal area of LARGER trees in its 20 m quadrat (BAL), the census twin of a MEDS
patch. A MEDS tree is one stem, so the census is read the same way: a tree's dbh is its main stem's;
a tree alive with no main stem to measure has lost its size, which is a death, and when a stem is
measured again it is a recruit; a main stem that breaks or changes is a change of the tree's dbh.
Only increments that are not comparable (the measuring height moved) or not credible (more than
75 mm/yr) are dropped. The 1982 census is not used: it rounded small stems to 5 mm.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS = sys.argv[1]
OUT = os.path.join(HERE, "data")
YEARS = (1985, 1990, 1995, 2000, 2005, 2010)
PLOT_AREA = 50.0                # [ha]
QUAD = 20.0                     # quadrat side [m]: the MEDS patch the census stands in for
BLOCK = 100.0                   # cross-validation block side [m]
RHO_BREAKS = (0.42, 0.58)       # [g/cm3] terciles of 1985 basal area: PFT 1 below, 2 between, 3 above
SIZE_CLASSES = [1, 2, 5, 10, 20, 50, 100, 1000]   # dbh class edges [cm] for census_stand.csv

wood_density = pd.read_csv(os.path.join(HERE, "bci_wood_density.csv")).set_index("sp").rho


def load(year):
    c = pd.read_csv(os.path.join(CENSUS, f"bci_{year}.csv"), low_memory=False,
                    usecols=["sp", "gx", "gy", "dbh", "hom", "date", "status"])
    c["dbh"] = c.dbh / 10.0                                   # mm -> cm
    c["pft"] = np.digitize(c.sp.map(wood_density), RHO_BREAKS) + 1
    c["quad"] = (np.floor(c.gx / QUAD) * 100 + np.floor(c.gy / QUAD)).astype("Int64")
    c["block"] = (np.floor(c.gx / BLOCK) * 10 + np.floor(c.gy / BLOCK)).astype("Int64")
    c["live"] = (c.status == "A") & c.dbh.notna() & c.quad.notna()
    return c


def basal_area_m2ha(dbh_cm, area_m2):
    return np.pi / 4 * (dbh_cm / 100) ** 2 / area_m2 * 1e4


def basal_area_of_larger(c):
    """BAL [m2/ha] of each live tree: the basal area of strictly larger trees in its quadrat."""
    a = c[c.live][["quad", "dbh"]].copy()
    a["ba"] = basal_area_m2ha(a.dbh, QUAD * QUAD)
    a = a.sort_values(["quad", "dbh"], ascending=[True, False])
    a["above"] = a.groupby("quad").ba.cumsum() - a.ba
    return a.groupby(["quad", "dbh"]).above.transform("min").reindex(c.index)   # ties: none above


os.makedirs(OUT, exist_ok=True)
census = {y: load(y) for y in YEARS}
for y in YEARS:
    census[y]["bal"] = basal_area_of_larger(census[y])

growth, survival, recruits = [], [], []
for y0, y1 in zip(YEARS[:-1], YEARS[1:]):
    c0, c1 = census[y0], census[y1]
    dt = (c1.date - c0.date) / 365.25
    start = pd.DataFrame({"interval": y0, "block": c0.block, "pft": c0.pft, "dbh": c0.dbh,
                          "bal": c0.bal, "dt": dt})

    inc_mm = (c1.dbh - c0.dbh) * 10
    ok = (c0.live & (c1.status == "A") & c1.dbh.notna()
          & (np.abs(c1.hom - c0.hom) <= 0.05 * c0.hom) & (inc_mm / dt <= 75))
    growth.append(start[ok].assign(g=(inc_mm / 10 / dt)[ok]))

    ok = c0.live & c1.status.isin(["A", "D"])                # trees not found are left out
    survival.append(start[ok].assign(dead=((c1.status == "D") | c1.dbh.isna())[ok].astype(int)))

    q_dt = (c1.date[c1.status == "A"].groupby(c1.quad).mean()
            - c0.date[c0.status == "A"].groupby(c0.quad).mean()) / 365.25
    live = c0[c0.live]
    ba = basal_area_m2ha(live.dbh, QUAD * QUAD)
    ba_tot = ba.groupby(live.quad).sum()
    ba_pft = ba.groupby([live.quad, live.pft]).sum().unstack(fill_value=0.0)
    new = (((c0.status == "P") | ((c0.status == "A") & c0.dbh.isna()))
           & (c1.status == "A") & c1.dbh.notna())
    n_new = new.groupby([c1.quad, c1.pft]).sum().unstack(fill_value=0)
    quads = ba_tot.index
    block = c0.groupby("quad").block.first().reindex(quads).values
    for p in (1, 2, 3):
        r = pd.DataFrame({"interval": y0, "block": block, "pft": p, "ba_tot": ba_tot.values,
                          "ba_pft": ba_pft.get(p, 0.0 * ba_tot).reindex(quads, fill_value=0.0).values,
                          "n_new": n_new.get(p, 0 * ba_tot).reindex(quads, fill_value=0).values,
                          "dt": q_dt.reindex(quads).values})
        recruits.append(r.assign(rate=r.n_new / (QUAD * QUAD) / r.dt))
    print(f"{y0}-{y1}: {len(growth[-1])} growth rows, {len(survival[-1])} survival rows "
          f"({survival[-1].dead.mean() * 100:.1f} % die), {int(new.sum())} recruits")

for name, parts in (("growth", growth), ("survival", survival), ("recruits", recruits)):
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(OUT, f"{name}.csv.gz"), index=False,
                                               float_format="%.6g")

# the 1985 stand as a MEDS census: one row per quadrat, PFT and dbh; MEDS fuses it when it reads it
live = census[1985][census[1985].live]
rows = live.groupby(["quad", "pft", "dbh"]).size().rename("n").reset_index()
rows["patch_id"] = rows.quad.astype(int) + 1
rows["patch_area"] = QUAD * QUAD
rows["nplant"] = rows.n / (QUAD * QUAD)
rows[["patch_id", "patch_area", "pft", "dbh", "nplant"]].to_csv(
    os.path.join(OUT, "bci_1985_census.csv"), index=False, float_format="%.6g")

# every census's stand by PFT and size class, for the comparison with the model
stand = []
for y in YEARS:
    live = census[y][census[y].live]
    k = np.digitize(live.dbh, SIZE_CLASSES[1:-1])
    for p in (1, 2, 3):
        for j in range(len(SIZE_CLASSES) - 1):
            m = (live.pft == p) & (k == j)
            stand.append((y, p, SIZE_CLASSES[j], m.sum() / PLOT_AREA,
                          basal_area_m2ha(live.dbh[m], PLOT_AREA * 1e4).sum()))
pd.DataFrame(stand, columns=["year", "pft", "dbh_class", "stems", "basal_area"]).to_csv(
    os.path.join(HERE, "census_stand.csv"), index=False, float_format="%.6g")
print(f"wrote data/ and census_stand.csv")
