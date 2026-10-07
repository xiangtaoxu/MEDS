# SPDX-License-Identifier: Apache-2.0
"""Census-trained vital rates for the demography engine.

The three laws are tables that fit_vital_rates.py tabulated from random forests trained on the BCI
censuses (vital_rates/*.csv). Each step this module describes every cohort the way the census
described every tree -- its dbh, its PFT and the basal area of LARGER trees in its own patch (BAL) --
looks its rates up, and returns the three arrays ``Site.apply_rates`` takes:

  growth       [cm/yr]        per cohort, from dbh, BAL and PFT
  mortality    [1/yr]         per cohort, from dbh, the growth above and PFT
  recruitment  [plants/m2/yr] per PFT and patch, from the patch's basal area and the PFT's share of it

The census counts every death, including the canopy trees that fall and open a gap. The engine kills
those through treefall disturbance (canopy cohorts on the disturbed area die), so that rate is taken
off the mortality of cohorts tall enough to die in a gap.
"""
import os

try:
    import tomllib                      # py3.11+
except ModuleNotFoundError:             # pragma: no cover -- py3.10 and older
    import tomli as tomllib

import numpy as np
import pandas as pd

TABLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vital_rates")


def _locate(grid, x):
    """Lower node and weight of x on a sorted grid. The weight is held to [0, 1], so beyond the grid
    the edge value holds -- what the forest itself does outside its data."""
    i = np.clip(np.searchsorted(grid, x, side="right") - 1, 0, len(grid) - 2)
    w = np.clip((x - grid[i]) / (grid[i + 1] - grid[i]), 0.0, 1.0)
    return i, w


class RateTable:
    """One rate on a regular (x, y) grid per PFT, read bilinearly; ``log_x`` grids dbh in log space."""

    def __init__(self, name, x, y, log_x=False):
        t = pd.read_csv(os.path.join(TABLES, f"{name}.csv")).sort_values(["pft", x, y])
        self.x = np.unique(t[x].to_numpy())
        self.y = np.unique(t[y].to_numpy())
        self.log_x = log_x
        if log_x:
            self.x = np.log(self.x)
        self.v = t[name].to_numpy().reshape(t.pft.nunique(), len(self.x), len(self.y))

    def __call__(self, pft, x, y):
        """The rate of 0-based ``pft`` at (x, y); the arguments broadcast against each other."""
        ix, wx = _locate(self.x, np.log(x) if self.log_x else x)
        iy, wy = _locate(self.y, y)
        v = self.v
        return ((1 - wx) * (1 - wy) * v[pft, ix, iy] + wx * (1 - wy) * v[pft, ix + 1, iy]
                + (1 - wx) * wy * v[pft, ix, iy + 1] + wx * wy * v[pft, ix + 1, iy + 1])


def basal_area_of_larger(dbh, nplant, patch):
    """[m2/ha] the basal area of larger trees that the average member of each cohort has: all of the
    larger cohorts in its own patch, and half of the trees of its own size -- within a cohort, the
    average tree has half of its peers above it. Cohorts of equal dbh (recruits of different PFTs
    born the same month) count as one size."""
    ba = nplant * np.pi / 4 * dbh ** 2                 # cm2 per m2 of patch == m2/ha
    o = np.lexsort((-dbh, patch))                      # by patch, largest first
    p, d, b = patch[o], dbh[o], ba[o]
    above = np.cumsum(b) - b                           # everything earlier in the sorted order
    new_patch = np.r_[True, p[1:] != p[:-1]]
    above -= np.maximum.accumulate(np.where(new_patch, above, 0.0))   # ... less the earlier patches
    new_size = new_patch | np.r_[True, d[1:] != d[:-1]]
    first = np.maximum.accumulate(np.where(new_size, np.arange(len(p)), 0))
    same_size = np.add.reduceat(b, np.flatnonzero(new_size))[np.cumsum(new_size) - 1]
    out = np.empty_like(above)
    out[o] = above[first] + 0.5 * same_size
    return out


class CensusLaws:
    def __init__(self, main_config):
        with open(main_config, "rb") as fh:
            disturbance = tomllib.load(fh)["disturbance"]
        self.treefall_rate = disturbance["patch_disturbance_rate"]
        self.treefall_height = disturbance["disturbance_survive_height"]
        self.growth = RateTable("growth", "dbh", "bal", log_x=True)
        self.mortality = RateTable("mortality", "dbh", "growth", log_x=True)
        self.recruitment = RateTable("recruitment", "ba_tot", "ba_pft")
        self.n_pft = self.growth.v.shape[0]

    def rates(self, state, n_patch):
        """(growth[n], mortality[n], recruitment[n_pft, n_patch]) for the stand in ``state``."""
        dbh, nplant = state["dbh"], state["nplant"]
        pft = state["pft"].astype(int) - 1
        patch = state["owner_patch"].astype(int) - 1

        ba_pft = np.zeros((self.n_pft, max(n_patch, 1)))
        np.add.at(ba_pft, (pft, patch), nplant * np.pi / 4 * dbh ** 2)
        recruitment = self.recruitment(np.arange(self.n_pft)[:, None], ba_pft.sum(0)[None, :], ba_pft)
        if len(dbh) == 0:
            return np.zeros(0), np.zeros(0), recruitment

        growth = self.growth(pft, dbh, basal_area_of_larger(dbh, nplant, patch))
        mortality = self.mortality(pft, dbh, growth)
        tall = state["height"] >= self.treefall_height
        mortality[tall] = np.maximum(mortality[tall] - self.treefall_rate, 0.0)
        return growth, mortality, recruitment
