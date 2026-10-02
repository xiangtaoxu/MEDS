# SPDX-License-Identifier: Apache-2.0
"""The staged fit (MEDS_FAST_CALIBRATION_REVISION_PLAN.md §7): the kernel models of the optics and
photosynthesis stages, the water stage's search, and the photosynthesis stage's gate G8.

Within the fast loop the coupling is mostly one-directional: optics -> light per leaf ->
photosynthesis -> stomata -> energy balance -> leaf temperature (a weak feedback). So a stage fits
its own keys with the upstream results fixed and the downstream drivers taken from a full run.

  optics          the two-stream alone (meds.canopy) over each window's stand and shortwave, against
                  the albedo: seconds per evaluation, no model run.
  photosynthesis  a canopy of leaf solves (meds.canopy) over each window's hourly per-cohort drivers,
                  against GPP; outer passes re-run the windows with the new keys and refresh the
                  drivers (the leaf-temperature feedback).
  energy, polish  the coupled fast loop: fit.Model on the 10-day windows.
  water           frozen seasonal runs, searched without derivatives (a grid, then a quadratic).

Both kernel models are ANCHORED to the full run they take their drivers from: their value for a
parameter set is the kernel's, times, per hour, the model's over the kernel's at the drivers'
parameters. At the anchor they reproduce the model exactly; elsewhere the kernel carries the
response. The anchor removes what the kernel cannot see -- the hourly output averages four 15-min
steps, and the kernel evaluates the hour's mean drivers once (gate G8 measures that error, before
the anchor).

Both expose `residuals(thetas, windows)` like fit.Model, so fit.lm and fit.jacobian drive them.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np

import residuals as res
import trials as T


def _patches(stand: dict):
    own = stand["stand_owner_patch"].astype(int)
    return [(p, np.flatnonzero(own == p)) for p in np.unique(own)]


class KernelBase:
    """What both kernel models share: the parameter list, a config writer, and the windows'
    drivers and full-run series."""

    def __init__(self, params, windows, specs, drivers: dict, series: dict, make_config, log=print):
        self.params = params
        self.windows = windows
        self.specs = specs                      # window name -> WindowSpec (this stage's target only)
        self.drivers = drivers                  # window name -> trials.load_drivers(...)
        self.series = series                    # window name -> the driver run's hourly DataFrame
        self.make_config = make_config          # theta -> path of a main TOML with those parameters
        self.log = log
        self.ratio = {}                         # window name -> per-hour anchor ratio
        self.n_evals = 0

    def set_anchor(self, theta_anchor):
        """The per-hour ratio of the full run's value to the kernel's at the drivers' parameters."""
        cfg = self.make_config(theta_anchor)
        for w in self.windows:
            k = self.kernel(cfg, w.name)
            m = self.model_value(w.name)
            self.ratio[w.name] = np.where(np.abs(k) > 1e-9, m / np.where(np.abs(k) > 1e-9, k, 1.0), 1.0)

    def residuals(self, thetas: list, windows=None) -> list:
        windows = self.windows if windows is None else windows
        out = []
        for th in thetas:
            cfg = self.make_config(th)
            parts = []
            for w in windows:
                df = self.series[w.name].copy()
                df[self.column] = self.kernel(cfg, w.name) * self.ratio.get(w.name, 1.0)
                parts.append(res.residual(self.specs[w.name], df))
            self.n_evals += 1
            out.append(np.concatenate(parts) if parts else np.zeros(0))
        return out


class OpticsModel(KernelBase):
    """Stage 1: the site's upwelling shortwave from the two-stream alone, against the albedo."""
    column = "sw_up_fast"

    def model_value(self, name):
        return self.series[name]["sw_up_fast"].to_numpy()

    def kernel(self, cfg, name):
        from meds.canopy import Canopy
        d = self.drivers[name]
        c = Canopy(cfg)
        lai = d["stand_nplant"] * d["stand_leaf_area"]
        up = np.zeros(len(d["index"]))
        for p, sel in _patches(d):
            r = c.radiation(d["cosz_fast"], d["par_beam_fast"], d["par_diffuse_fast"], d["nir_beam_fast"],
                            d["nir_diffuse_fast"], d["stand_pft"].astype(int)[sel], d["stand_height"][sel],
                            lai[sel], d["stand_wai"][sel])
            up += d["stand_patch_area"][p - 1] * (r["up_vis"] + r["up_nir"])
        return up


class PhotoModel(KernelBase):
    """Stage 2: the site's GPP from a canopy of leaf solves over the hourly per-cohort drivers,
    against the GPP target's hours (only those hours are computed)."""
    column = "gpp_rate_fast"

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.rows = {}
        for w in self.windows:
            d = self.drivers[w.name]
            spec = self.specs[w.name]
            hours = np.unique(np.concatenate([t.hours for t in spec.targets])) if spec.targets else np.zeros(0, int)
            lai = d["stand_nplant"] * d["stand_leaf_area"]
            hh, jj = np.meshgrid(hours, np.flatnonzero(lai > 0), indexing="ij")
            self.rows[w.name] = (hours, hh.ravel(), jj.ravel(),
                                 d["stand_patch_area"][d["stand_owner_patch"].astype(int) - 1] * lai)

    def model_value(self, name):
        return self.series[name]["gpp_rate_fast"].to_numpy()

    def kernel(self, cfg, name):
        from meds.canopy import Canopy
        d = self.drivers[name]
        hours, h, j, w = self.rows[name]
        c = Canopy(cfg)
        pft = d["stand_pft"].astype(int)
        vc, rd = c.plastic_traits(pft, d["stand_overtopping_lai"])
        a = c.leaf(pft[j], vc[j], rd[j], par=d["gx_par_cohort_fast"][h, j],
                   leaf_temp=d["gx_leaf_temp_cohort_fast"][h, j], vpd=d["gx_vpd_cohort_fast"][h, j],
                   ca=d["gx_ca_cohort_fast"][h, j], pressure=d["gx_pressure_cohort_fast"][h, j],
                   psi_leaf=d["gx_psi_leaf_cohort_fast"][h, j], gb=d["gx_gb_cohort_fast"][h, j],
                   psi=d["gx_psi_predawn_cohort_fast"][h, j])
        gpp = self.model_value(name).astype(float).copy()       # hours outside the target keep the model's
        acc = np.zeros(len(gpp))
        np.add.at(acc, h, a["A_gross"] * w[j])
        gpp[hours] = acc[hours]
        return gpp

    def g8(self, theta_anchor) -> dict:
        """Gate G8: the kernel, BEFORE the anchor, against the full run's hourly GPP on the target's
        hours (GPP > 1 umol m-2 s-1): median and largest relative difference, and the mean bias."""
        cfg = self.make_config(theta_anchor)
        rel = []
        for w in self.windows:
            hours = self.rows[w.name][0]
            k, m = self.kernel(cfg, w.name)[hours], self.model_value(w.name)[hours]
            ok = m > 1.0
            rel.append(k[ok] / m[ok] - 1.0)
        rel = np.concatenate(rel) if rel else np.zeros(0)
        if not len(rel):
            return {"pass": False, "note": "no GPP hours"}
        med = float(np.median(np.abs(rel)))
        return {"pass": med <= 0.01, "median_abs_rel": med, "max_abs_rel": float(np.max(np.abs(rel))),
                "mean_rel": float(np.mean(rel)), "hours": int(len(rel))}


def config_writer(base, params, overrides, root: Path):
    """theta -> the path of a main TOML (and its PFT file) with those parameters, for meds.canopy."""
    cache = {}

    def make(theta):
        cfg = T.with_params(base, params, theta, overrides)
        key = T.digest(cfg.main, cfg.pft)[:16]
        if key not in cache:
            d = root / key
            if not (d / "main.toml").exists():
                d.mkdir(parents=True, exist_ok=True)
                cfg.write(d)
            cache[key] = str(d / "main.toml")
        return cache[key]
    return make


def grid_search(prob, u_start, rounds=2, width=1.0, log=print, label="water"):
    """A derivative-free search for keys whose response can be rough (the water stage): a 3^k grid
    around the current point (width prior sds in u), a quadratic fitted to the grid's objective and
    its minimum tried inside the grid's box, the best point kept; then a grid half as wide around it.
    Returns the best u, its objective, and the history."""
    k = len(u_start)
    sig = prob.sigma
    centre = np.asarray(u_start, dtype=float)
    r0 = prob.data([centre])[0]
    if r0 is None:
        raise RuntimeError(f"{label}: the starting point's trials failed")
    best_u, best_c = centre, prob.cost(centre, r0)
    hist = [{"round": 0, "u": centre.tolist(), "cost": best_c}]
    log(f"{label} start: Phi = {best_c:.6g}")
    for rnd in range(1, rounds + 1):
        offs = [np.array(o) for o in itertools.product((-1.0, 0.0, 1.0), repeat=k) if any(o)]
        pts = [centre + o * width * sig for o in offs]
        R = prob.data(pts)
        rows, ys = [], []
        for o, pt, rv in zip(offs, pts, R):
            if rv is None:
                continue
            c = prob.cost(pt, rv)
            ys.append(c)
            rows.append(o)
            if c < best_c:
                best_u, best_c = pt, c
        ys.append(prob.cost(centre, prob.data([centre])[0]))
        rows.append(np.zeros(k))
        #----- the quadratic c + g.x + x'Hx/2 in grid units x, by least squares
        X = []
        for o in rows:
            quad = [o[a] * o[b] * (0.5 if a == b else 1.0) for a in range(k) for b in range(a, k)]
            X.append(np.concatenate([[1.0], o, quad]))
        coef, *_ = np.linalg.lstsq(np.array(X), np.array(ys), rcond=None)
        g = coef[1:1 + k]
        H = np.zeros((k, k))
        idx = 1 + k
        for a in range(k):
            for b in range(a, k):
                H[a, b] = H[b, a] = coef[idx]
                idx += 1
        tried = None
        if np.all(np.linalg.eigvalsh(H) > 0):
            x = np.clip(-np.linalg.solve(H, g), -1.0, 1.0)
            tried = centre + x * width * sig
            rv = prob.data([tried])[0]
            if rv is not None:
                c = prob.cost(tried, rv)
                if c < best_c:
                    best_u, best_c = tried, c
        hist.append({"round": rnd, "width": width, "best_u": np.asarray(best_u).tolist(), "cost": best_c,
                     "quadratic_min": None if tried is None else tried.tolist()})
        log(f"{label} round {rnd} (width {width:g} sd): Phi = {best_c:.6g}")
        centre = np.asarray(best_u, dtype=float)
        width /= 2.0
    return best_u, best_c, hist
