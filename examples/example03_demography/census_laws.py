# SPDX-License-Identifier: Apache-2.0
"""Census-trained vital rates for the demography engine.

fit_vital_rates.py fits three laws per PFT to the BCI censuses and writes their coefficients to
vital_rates.json. Each step this module describes every cohort the way the census described every
tree -- its dbh D, its PFT and the basal area of LARGER trees in its own patch (BAL) -- and returns
the three arrays ``Site.apply_rates`` takes:

  growth       [cm/yr]        [g_min + (g_max - g_min) / (1 + exp(-k (lnD - lnD0)))] exp(-b BAL)
  mortality    [1/yr]         gamma + alpha exp(-beta growth)          (Camac et al. 2018)
  recruitment  [plants/m2/yr] exp(c0 + c1 BA + c2 BA_pft), per PFT and patch

Recruitment's basal areas beyond the census's range take the range's edge. The census counts every death, including the
canopy trees that fall and open a gap. The engine kills those through treefall disturbance (canopy
cohorts on the disturbed area die), so that rate is taken off the mortality of cohorts tall enough to
die in a gap.
"""
import json
import os

try:
    import tomllib                      # py3.11+
except ModuleNotFoundError:             # pragma: no cover -- py3.10 and older
    import tomli as tomllib

import numpy as np

COEFFICIENTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vital_rates.json")


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
        with open(COEFFICIENTS) as fh:
            laws = json.load(fh)
        self.g = np.array([laws["growth"][k] for k in ("g_min", "g_max", "D0", "k", "b")]).T   # [pft, 5]
        self.camac = np.array([laws["mortality"][k] for k in ("gamma", "alpha", "beta")]).T
        self.c = np.array(laws["recruitment"]["coefficients"])                 # [pft, 3]
        self.range = laws["range"]
        self.n_pft = len(self.g)

    def growth(self, pft, dbh, bal):
        """[cm/yr] for 0-based ``pft``; the arguments broadcast against each other."""
        t = self.g[pft]
        size = 1.0 / (1.0 + np.exp(-t[..., 3] * (np.log(dbh) - np.log(t[..., 2]))))
        return (t[..., 0] + (t[..., 1] - t[..., 0]) * size) * np.exp(-t[..., 4] * bal)

    def mortality(self, pft, growth):
        """[1/yr] at the given growth [cm/yr]."""
        t = self.camac[pft]
        return t[..., 0] + t[..., 1] * np.exp(-t[..., 2] * growth)

    def recruitment(self, pft, ba, ba_pft):
        """[plants/m2/yr] in a patch of basal area ``ba`` [m2/ha], ``ba_pft`` of it the PFT's own."""
        ba, ba_pft = np.clip(ba, *self.range["ba"]), np.clip(ba_pft, *self.range["ba"])
        c = self.c[pft]
        return np.exp(c[..., 0] + c[..., 1] * ba + c[..., 2] * ba_pft)

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
        mortality = self.mortality(pft, growth)
        tall = state["height"] >= self.treefall_height
        mortality[tall] = np.maximum(mortality[tall] - self.treefall_rate, 0.0)
        return growth, mortality, recruitment
