# SPDX-License-Identifier: Apache-2.0
"""The estimation (MEDS_FAST_CALIBRATION_PLAN.md §6): Levenberg-Marquardt on the stacked, weighted
residuals of every target in every window, with the Jacobian from central finite differences run
in parallel, Gaussian priors in the transformed space, the screening from the first Jacobian, and
the Laplace covariance at the MAP with a linearity check.

The objective in u (free parameters only) is

    Phi(u) = || r_data(u) ||^2 + || (u - u_prior) / sigma_prior ||^2

The data residuals come from `Model.residuals`, which runs one trial per (point, window) through
the worker pool; every trial of a batch runs at once.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import residuals as res
import trials
from pool import Task
from registry import H_U, SIGMA_U


@dataclass
class Model:
    """Runs candidate parameter sets over a set of windows and returns their data residuals."""
    params: list                     # every registry parameter (free and fixed)
    windows: list                    # trials.Window
    specs: dict                      # window name -> residuals.WindowSpec
    states: dict                     # window name -> state file
    base_main: dict
    base_pft: dict
    overrides: dict
    exe: str
    pool: object
    root: Path
    utc_offset_h: float
    growth_resp: dict | None = None
    timeout: float = 900.0
    keep_netcdf: bool = False
    seconds: list = field(default_factory=list)
    seconds_ok: list = field(default_factory=list)
    n_trials: int = 0
    n_failed: int = 0
    log: object = print

    def _timeout(self) -> float:
        """A runaway trial (the hydraulic cost cliff) is cut off; a slow one under a full node is
        not. The median is of completed trials only, and the floor is `self.timeout`, because a
        trial's time depends on how many share its node: the first trials of a fit run on an idle
        node, several times faster than a full one."""
        ok = self.seconds_ok[-2000:]
        if len(ok) < 8:
            return self.timeout
        return max(3.0 * float(np.median(ok)), self.timeout)

    def residuals(self, thetas: list, windows=None) -> list:
        """One stacked data-residual vector per parameter set (theta over self.params), or None
        where any of its trials failed."""
        windows = self.windows if windows is None else windows
        dirs = [[trials.build_trial(self.base_main, self.base_pft, self.params, th, w,
                                    self.states[w.name], self.root, self.overrides)
                 for w in windows] for th in thetas]
        todo, seen = [], set()
        for tds in dirs:
            for td in tds:
                if td not in seen and not (td / "series.npz").exists():
                    seen.add(td)
                    todo.append(Task(td.name, trials.command(self.exe, td), str(td), str(td / "run.log"),
                                     self._timeout()))
        status = self.pool.run(todo) if todo else {}
        for t in todo:
            self.seconds.append(status[t.id][1])
            if status[t.id][0] == "ok":
                self.seconds_ok.append(status[t.id][1])
        self.n_trials += len(todo)
        out = []
        for th, tds in zip(thetas, dirs):
            parts = []
            for w, td in zip(windows, tds):
                try:
                    if (td / "series.npz").exists():
                        df = trials.load_series(td)
                    else:
                        st = status.get(td.name, ("missing", 0.0))[0]
                        if st != "ok":
                            raise trials.TrialError(f"{td.name}: {st}")
                        df = trials.finish(td, self.params, th, self.utc_offset_h, self.keep_netcdf)
                    parts.append(res.residual(self.specs[w.name], df, self.growth_resp))
                except (trials.TrialError, ValueError) as e:
                    if "parameter record check failed" in str(e):
                        raise                          # a harness bug, not a model failure
                    self.n_failed += 1
                    self.log(f"  trial failed: {e}".splitlines()[0])
                    parts = None
                    break
            out.append(None if parts is None else np.concatenate(parts))
        return out


@dataclass
class Problem:
    """The free parameters, their prior, and the map u -> theta over every parameter."""
    model: Model
    free: list                       # indices into model.params
    theta_fixed: np.ndarray          # theta of every parameter; free ones are overwritten

    @property
    def params(self):
        return [self.model.params[i] for i in self.free]

    @property
    def u_prior(self) -> np.ndarray:
        return np.array([p.u0 for p in self.params])

    def theta(self, u) -> np.ndarray:
        th = self.theta_fixed.copy()
        for k, i in enumerate(self.free):
            th[i] = float(self.model.params[i].to_theta(u[k]))
        return th

    def prior_r(self, u) -> np.ndarray:
        return (np.asarray(u) - self.u_prior) / SIGMA_U

    def cost(self, u, r_data) -> float:
        return float(r_data @ r_data + self.prior_r(u) @ self.prior_r(u))

    def data(self, us: list, windows=None) -> list:
        return self.model.residuals([self.theta(u) for u in us], windows)


def jacobian(prob: Problem, u, r0, h=H_U, log=print):
    """Central differences in u, every column's two trials in one batch. Returns the data
    Jacobian, the per-column smoothness ratio (the + and - one-sided slopes' agreement, 1 for a
    smooth key) and the columns that failed twice (set to zero)."""
    k = len(u)
    J = np.zeros((len(r0), k))
    smooth = np.full(k, np.nan)
    pending = list(range(k))
    steps = np.full(k, h)
    for attempt in range(2):
        pts = []
        for j in pending:
            e = np.zeros(k)
            e[j] = steps[j]
            pts += [u + e, u - e]
        R = prob.data(pts)
        retry = []
        for n, j in enumerate(pending):
            rp, rm = R[2 * n], R[2 * n + 1]
            if rp is None or rm is None:
                retry.append(j)
                continue
            J[:, j] = (rp - rm) / (2.0 * steps[j])
            dp, dm = rp - r0, r0 - rm
            den = float(dm @ dm)
            smooth[j] = float(dp @ dm) / den if den > 0 else (1.0 if float(dp @ dp) == 0 else np.inf)
        pending = retry
        steps[retry] = steps[retry] / 2.0
        if not pending:
            break
        log(f"  Jacobian: {len(pending)} column(s) failed; retrying at h/2")
    return J, smooth, pending


def lm(prob: Problem, u_start, max_iter=15, rtol=1e-3, lam0=1e-2, max_step=2.0, log=print,
       label="", jac0=None):
    """Levenberg-Marquardt from u_start. Each iteration tries three damping values at once."""
    u = np.array(u_start, dtype=float)
    r = prob.data([u])[0]
    if r is None:
        raise RuntimeError(f"{label}: the starting point's trials failed")
    cost = prob.cost(u, r)
    J, smooth, failed = jac0 if jac0 is not None else jacobian(prob, u, r, log=log)
    lam = lam0
    hist = [{"iter": 0, "cost": cost, "lambda": lam, "u": u.tolist()}]
    log(f"{label} start: Phi = {cost:.6g}")
    for it in range(1, max_iter + 1):
        Jf = np.vstack([J, np.eye(len(u)) / SIGMA_U])
        rf = np.concatenate([r, prob.prior_r(u)])
        A, g = Jf.T @ Jf, Jf.T @ rf
        D = np.diag(np.maximum(np.diag(A), 1e-12))
        accepted = False
        for _ in range(4):
            lams = [lam / 10.0, lam, lam * 10.0]
            cands = []
            for l in lams:
                d = np.linalg.solve(A + l * D, -g)
                s = np.max(np.abs(d))
                cands.append(u + (d * max_step / s if s > max_step else d))
            R = prob.data(cands)
            costs = [prob.cost(c, rc) if rc is not None else np.inf for c, rc in zip(cands, R)]
            b = int(np.argmin(costs))
            if costs[b] < cost:
                rel = (cost - costs[b]) / cost
                u, r, cost, lam = cands[b], R[b], costs[b], lams[b]
                accepted = True
                break
            lam *= 100.0
        hist.append({"iter": it, "cost": cost, "lambda": lam, "u": u.tolist(), "accepted": accepted})
        if not accepted:
            log(f"{label} iter {it}: no descent at any damping; stopping")
            break
        log(f"{label} iter {it}: Phi = {cost:.6g} (drop {100 * rel:.2f} %, lambda {lam:.3g})")
        J, smooth, failed = jacobian(prob, u, r, log=log)
        if rel < rtol:
            break
    return {"u": u, "r": r, "cost": cost, "J": J, "smooth": smooth, "failed": failed, "history": hist}


def posterior(prob: Problem, J, r, specs, weighted=True):
    """Laplace covariance in u: (J' W J + Sigma_prior^-1)^-1, W the effective-sample-size weights."""
    w = res.ess_weights(specs, r) if weighted else np.ones_like(r)
    Jw = J * np.sqrt(w)[:, None]
    H = Jw.T @ Jw + np.eye(J.shape[1]) / SIGMA_U ** 2
    cov_u = np.linalg.inv(H)
    return cov_u, w


def theta_cov(prob: Problem, u, cov_u):
    dth = np.array([p.dtheta_du(uk) for p, uk in zip(prob.params, u)])
    return cov_u * np.outer(dth, dth)


def correlation(cov):
    s = np.sqrt(np.diag(cov))
    return cov / np.outer(s, s)


def screening(prob: Problem, J, r0, smooth, failed, specs, max_free=20, min_teach=0.05):
    """The first Jacobian's screening (§6.3): sensitivity per target, identifiability as each key's
    posterior-to-prior sigma ratio, collinear pairs, dead columns and rough keys. Returns the keys
    to fit and a report."""
    params = prob.params
    cov_u, _ = posterior(prob, J, r0, specs)
    ratio = np.sqrt(np.diag(cov_u)) / SIGMA_U
    corr = correlation(cov_u)
    sl = []
    i = 0
    for spec in specs:
        for t in spec.targets:
            sl.append((t.name, slice(i, i + len(t.obs))))
            i += len(t.obs)
    sens = {}
    for j, p in enumerate(params):
        per = {}
        for name, s in sl:
            col = J[s, j] * SIGMA_U
            acc = per.setdefault(name, [0.0, 0])
            acc[0] += float(col @ col)
            acc[1] += len(col)
        sens[p.name] = {k: math.sqrt(v[0] / v[1]) for k, v in per.items()}
    dead = [p.name for j, p in enumerate(params) if not np.any(J[:, j])]
    rough = [p.name for j, p in enumerate(params)
             if np.isfinite(smooth[j]) and (smooth[j] < 1 / 3 or smooth[j] > 3) and p.name not in dead]
    rough += [params[j].name for j in failed]
    candidates = [j for j, p in enumerate(params)
                  if p.name not in dead and p.name not in rough and ratio[j] < 1.0 - min_teach]
    candidates.sort(key=lambda j: ratio[j])
    keep = sorted(candidates[:max_free])
    pairs = [(params[a].name, params[b].name, float(corr[a, b]))
             for a in range(len(params)) for b in range(a + 1, len(params)) if abs(corr[a, b]) > 0.9]
    sv = np.linalg.svd(J * SIGMA_U, compute_uv=False)
    report = {"keys": [p.name for p in params],
              "sigma_ratio": {p.name: float(ratio[j]) for j, p in enumerate(params)},
              "smoothness": {p.name: float(smooth[j]) for j, p in enumerate(params)},
              "sensitivity": sens, "dead": dead, "rough": rough, "collinear": pairs,
              "singular_values": sv.tolist(),
              "fitted": [params[j].name for j in keep],
              "fixed": [p.name for j, p in enumerate(params) if j not in keep]}
    return keep, report


def linearity(prob: Problem, u, cov_u, J, n_dir=3, log=print):
    """Along the leading principal directions of the posterior, +-1 sigma: the actual change in the
    fit's objective Phi against the quadratic prediction d'(J'J + P)d from the same Jacobian. A
    ratio outside [0.5, 2] marks the covariance as local only (§6.4)."""
    vals, vecs = np.linalg.eigh(cov_u)
    order = np.argsort(vals)[::-1][:n_dir]
    r0 = prob.data([u])[0]
    base = prob.cost(u, r0)
    H = J.T @ J + np.eye(len(u)) / SIGMA_U ** 2
    pts, pred = [], []
    for k in order:
        d = math.sqrt(vals[k]) * vecs[:, k]
        pts += [u + d, u - d]
        pred.append(float(d @ H @ d))
    R = prob.data(pts)
    out = []
    for n, k in enumerate(order):
        acts = [prob.cost(pts[2 * n + s], R[2 * n + s]) - base if R[2 * n + s] is not None else None
                for s in (0, 1)]
        ratios = [a / pred[n] if a is not None and pred[n] > 0 else None for a in acts]
        out.append({"direction": int(k), "sigma_u": float(math.sqrt(vals[k])), "predicted": pred[n],
                    "delta_phi_plus": acts[0], "delta_phi_minus": acts[1], "ratio_plus": ratios[0],
                    "ratio_minus": ratios[1],
                    "local_only": any(q is None or not (0.5 <= q <= 2.0) for q in ratios),
                    "loadings": {p.name: float(vecs[j, k]) for j, p in enumerate(prob.params)}})
    return out


def line_search(prob: Problem, u, j, grid=(-2.0, -1.0, 1.0, 2.0)):
    """1-D search for a rough key j, others held: the best of u_j + g * sigma_prior."""
    pts = []
    for g in grid:
        v = u.copy()
        v[j] = u[j] + g * SIGMA_U
        pts.append(v)
    R = prob.data(pts)
    base = prob.data([u])[0]
    best_u, best_c = u, prob.cost(u, base)
    for v, rv in zip(pts, R):
        if rv is not None and prob.cost(v, rv) < best_c:
            best_u, best_c = v, prob.cost(v, rv)
    return best_u, best_c


def bound_pushers(prob: Problem, u, J, r, specs, frac=0.05):
    """Parameters within `frac` of a bound, with the target whose gradient pushes them there."""
    out = []
    sl = []
    i = 0
    for spec in specs:
        for t in spec.targets:
            sl.append((t.name, slice(i, i + len(t.obs))))
            i += len(t.obs)
    for j, p in enumerate(prob.params):
        th = float(p.to_theta(u[j]))
        pos = (th - p.lo) / (p.hi - p.lo)
        if pos < frac or pos > 1.0 - frac:
            push = {}
            for name, s in sl:
                push[name] = push.get(name, 0.0) + float(J[s, j] @ r[s])
            # the descent direction for key j is -gradient; the target pushing toward the near
            # bound is the one whose contribution has the sign that moves u that way
            toward_hi = pos > 0.5
            pusher = min(push, key=lambda k: push[k]) if toward_hi else max(push, key=lambda k: push[k])
            out.append({"key": p.name, "theta": th, "range": [p.lo, p.hi], "pushed_by": pusher})
    return out


def save_json(path, obj):
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, Path):
            return str(o)
        raise TypeError(type(o).__name__)
    Path(path).write_text(json.dumps(obj, indent=1, default=conv))
