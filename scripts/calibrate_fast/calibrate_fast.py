#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""calibrate_fast.py -- fit MEDS's fast (sub-daily) parameters to a flux tower, with the stand frozen
at its initial structure (docs/dev_plans/MEDS_FAST_CALIBRATION_PLAN.md).

Every trial is a 10-day frozen run restarted from a shared state at its window's start; the fit is
Levenberg-Marquardt on the stacked residuals of every target in every window, with a parallel
central-difference Jacobian, Gaussian priors, and the Laplace covariance at the MAP.

Commands (each reads the site declaration, e.g. examples/example_flux_tower_bci/calibration.toml):
  select-windows  pick the 10-day windows with the best tower coverage in each declared season
  growth-resp     the stand's monthly growth-respiration climatology from a slow-loop run
  check           one trial per window at the default: runs, the parameter record, repeatability
  fit             screening, the fit from several starts, the state refresh, covariance,
                  linearity, validation, the gates and the calibrated configs
  worker          a node's worker for --pool queue (started inside the Slurm allocation)
  smoke           the CTest smoke test: one Jacobian column on a short window, a repeated trial

Usage:
  calibrate_fast.py fit --site calibration.toml --variant interception_off --work runs/cal \
      --meds-main build-ifx/meds_main --pool local --workers 40
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fit as F                      # noqa: E402
import residuals as R                # noqa: E402
import states as S                   # noqa: E402
import tomlio                        # noqa: E402
import tower as TW                   # noqa: E402
import trials as T                   # noqa: E402
from pool import make_pool, worker   # noqa: E402
from registry import SIGMA_U, load_registry, resolve_defaults   # noqa: E402


# ------------------------------------------------------------------------------------------------
# The site declaration
# ------------------------------------------------------------------------------------------------
class Site:
    def __init__(self, path, variant=None):
        self.path = Path(path).resolve()
        self.dir = self.path.parent
        self.decl = tomlio.load(self.path)
        d = self.decl
        b = d["base"]
        self.main_path = self._p(b["main"])
        self.pft_path = self._p(b["pft"])
        self.base_main = T.absolutize(tomlio.load(self.main_path), self.main_path.parent)
        self.base_pft = tomlio.load(self.pft_path)
        self.registry_path = self._p(b["registry"])
        self.growth_resp_path = self._p(b["growth_resp"]) if b.get("growth_resp") else None
        self.variant = variant
        self.overrides = dict(d.get("overrides", {}))
        if variant is not None:
            if variant not in d.get("variants", {}):
                raise SystemExit(f"unknown variant '{variant}'; declared: {list(d.get('variants', {}))}")
            self.overrides.update(d["variants"][variant])
        tw = d["tower"]
        self.tower = TW.TowerSpec(path=str(self._p(tw["path"])), columns=dict(tw["columns"]),
                                  time_column=tw.get("time_column", "date"),
                                  flag_column=tw.get("flag_column", "FLAG"),
                                  flag_good=tw.get("flag_good", 1),
                                  utc_offset_h=float(tw.get("utc_offset_h", 0.0)),
                                  forcing=str(self._p(tw["forcing"])) if tw.get("forcing") else None,
                                  forcing_qc=tuple(tw.get("forcing_qc", ("LWdown_qc", "Wind_qc"))),
                                  forcing_qc_val=tuple(tw.get("forcing_qc_val", ("Wind_qc",))),
                                  forcing_grid=int(tw.get("forcing_grid", 1)),
                                  closure_days=int(tw.get("closure_days", 31)),
                                  daytime_sw=float(tw.get("daytime_sw", 10.0)))
        w = d["windows"]
        self.days = int(w.get("days", 10))
        self.skip_hours = int(w.get("skip_hours", 0))
        self.chains = {k: dt.datetime.fromisoformat(v) for k, v in w["chains"].items()}
        self.windows = [T.Window(e["name"], dt.datetime.fromisoformat(e["start"]),
                                 int(e.get("days", self.days)), e["role"], e["chain"])
                        for e in w.get("list", [])]
        self.targets = d["targets"]
        self.fitcfg = d.get("fit", {})

    def _p(self, rel) -> Path:
        return (self.dir / rel).resolve()

    def params(self, record=None):
        ps = load_registry(self.registry_path, self.variant)
        only = self.fitcfg.get("keys")
        if only:
            ps = [p for p in ps if p.name in only]
        return resolve_defaults(ps, self.base_main, self.base_pft, record)


def log_to(path):
    fh = open(path, "a")
    lock = threading.Lock()

    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        with lock:
            print(line, flush=True)
            fh.write(line + "\n")
            fh.flush()
    return log


def observations(site: Site):
    """The tower's hours, and the observed-forcing masks for calibration and validation windows."""
    obs = TW.load_tower(site.tower)
    fok = {"cal": TW.forcing_observed(site.tower, obs.index),
           "val": TW.forcing_observed(site.tower, obs.index, site.tower.forcing_qc_val)}
    return obs, fok


def specs_for(site: Site, windows, obs, fok, log=print):
    out = {}
    for w in windows:
        idx = R.window_index(w.start, w.days, site.tower.utc_offset_h)
        out[w.name] = R.build_spec(w.name, idx, obs, fok[w.role], site.targets, site.skip_hours,
                                   site.tower.daytime_sw)
        if out[w.name].n == 0:
            log(f"WARNING: window {w.name} has no usable hours; it contributes nothing")
    return out


def run_chains(site, params, theta, root, pool, exe, log, windows=None):
    windows = site.windows if windows is None else windows
    states, errs, threads = {}, [], []

    def one(name):
        ws = [w for w in windows if w.chain == name]
        if not ws:
            return
        try:
            log(f"chain {name}: {len(ws)} windows from {site.chains[name]:%Y-%m-%d}")
            states.update(S.run_chain(name, site.chains[name], ws, site.base_main, site.base_pft,
                                      params, theta, root, pool, exe, site.overrides))
        except Exception as e:           # noqa: BLE001 -- reported below
            errs.append(e)
    for name in site.chains:
        th = threading.Thread(target=one, args=(name,))
        th.start()
        threads.append(th)
    for th in threads:
        th.join()
    if errs:
        raise errs[0]
    return states


def base_record(site, work, pool, exe, log):
    """Run one short trial at the base configs to read the parameter record: the defaults of keys
    the base TOML leaves to the model."""
    rdir = work / "record"
    rec_csv = rdir / "out" / f"{T.PREFIX}_parameters.csv"
    if not rec_csv.exists():
        rdir.mkdir(parents=True, exist_ok=True)
        (rdir / "out").mkdir(exist_ok=True)
        main, pft = tomlio.clone(site.base_main), tomlio.clone(site.base_pft)
        for k, v in site.overrides.items():
            tomlio.deep_set(main, k, v)
        (rdir / "pft.toml").write_text(tomlio.dumps(pft))
        w = site.windows[0]
        for k, v in {"run.start_time": T.stamp(w.start), "run.end_time": T.stamp(w.start + dt.timedelta(days=1)),
                     "run.slow_on": False, "run.n_threads": 1, "init.pft_config": str(rdir / "pft.toml"),
                     "state.write_state": False, "output.enabled": True, "output.dir": str(rdir / "out"),
                     "output.prefix": T.PREFIX, "output.daily.enabled": False, "output.fast.enabled": False,
                     "output.monthly.enabled": False, "output.annual.enabled": False}.items():
            tomlio.deep_set(main, k, v)
        tomlio.write(rdir / "main.toml", main)
        from pool import Task
        log("reading the parameter record from a one-day base run")
        pool.run([Task("record", [exe, str(rdir / "main.toml")], str(rdir), str(rdir / "run.log"), 3600)])
    return T.read_record(rec_csv, rdir / "main.toml", rdir / "pft.toml")


# ------------------------------------------------------------------------------------------------
# select-windows
# ------------------------------------------------------------------------------------------------
def cmd_select_windows(args):
    site = Site(args.site)
    obs, fok = observations(site)
    decl = site.decl["windows"].get("seasons", [])
    if not decl:
        raise SystemExit("declare [[windows.seasons]] (name, from, to, role, chain) to select from")
    day_ok = pd.DataFrame({
        "flux": (obs["h"].notna() & obs["le"].notna()),
        "lw": fok["cal"] & obs["lw_up"].notna(),
    }).astype(float).resample("1D").mean()
    rows = []
    taken = []
    for s in decl:
        a, b = pd.Timestamp(s["from"]), pd.Timestamp(s["to"])
        best, best_score = None, -1.0
        for start in pd.date_range(a, b - pd.Timedelta(days=site.days), freq="1D"):
            span = day_ok.loc[start:start + pd.Timedelta(days=site.days - 1)]
            if len(span) < site.days:
                continue
            need_lw = s.get("role", "cal") == "cal"
            score = float(span["flux"].mean()) * (float(span["lw"].mean()) if need_lw else 1.0)
            clash = any(abs((start - t).days) < site.days for t in taken)
            if not clash and score > best_score:
                best, best_score = start, score
        if best is None:
            print(f"# {s['name']}: no complete {site.days}-day span in {a:%Y-%m-%d}..{b:%Y-%m-%d}")
            continue
        taken.append(best)
        # the run's start is UTC midnight; the tower's local day starts 5 h later at BCI
        rows.append((s["name"], best, s.get("role", "cal"), s.get("chain", "cal"), best_score))
    for name, start, role, chain, score in rows:
        print(f'[[windows.list]]\nname = "{name}"\nstart = "{start:%Y-%m-%d}"\nrole = "{role}"\n'
              f'chain = "{chain}"   # coverage score {score:.2f}\n')


# ------------------------------------------------------------------------------------------------
# growth-resp
# ------------------------------------------------------------------------------------------------
def cmd_growth_resp(args):
    import glob
    from netCDF4 import Dataset
    frames = []
    for p in sorted(glob.glob(args.daily)):
        with Dataset(p) as ds:
            if "growth_resp_site" not in ds.variables:
                continue
            when = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:]))
            frames.append(pd.Series(np.asarray(ds["growth_resp_site"][:], float).squeeze(), index=when))
    if not frames:
        raise SystemExit(f"no growth_resp_site in {args.daily}")
    gr = pd.concat(frames).sort_index()
    gr = gr.where(gr.abs() < 1e30) * R.KGC_YR_TO_UMOL_S
    clim = gr.groupby(gr.index.month).mean()
    out = pd.DataFrame({"month": clim.index, "growth_resp_umol": clim.values})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.6g")
    print(out.to_string(index=False))


# ------------------------------------------------------------------------------------------------
# check and smoke
# ------------------------------------------------------------------------------------------------
def setup_model(site, args, work, pool, log, windows=None, params=None, theta=None, states=None):
    obs, fok = observations(site)
    windows = site.windows if windows is None else windows
    record = base_record(site, work, pool, args.meds_main, log)
    params = site.params(record) if params is None else params
    theta = np.array([p.default for p in params]) if theta is None else theta
    if states is None:
        states = run_chains(site, params, theta, work / "chains", pool, args.meds_main, log, windows)
    specs = specs_for(site, windows, obs, fok, log)
    gr = R.load_growth_resp(site.growth_resp_path) if site.growth_resp_path and site.growth_resp_path.exists() else None
    if gr is None:
        log("note: no growth-respiration climatology; night NEE is compared without it")
    model = F.Model(params, windows, specs, states, site.base_main, site.base_pft, site.overrides,
                    args.meds_main, pool, work / "trials", site.tower.utc_offset_h, gr,
                    float(site.fitcfg.get("timeout", 900)), args.keep_netcdf, log=log)
    return model, params, theta, obs, fok


def cmd_check(args):
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "check.log")
    pool = make_pool(args.pool, args.workers, work / "queue")
    cal = [w for w in site.windows if w.role == "cal"]
    model, params, theta, _, _ = setup_model(site, args, work, pool, log, windows=cal)
    t0 = time.time()
    r = model.residuals([theta])[0]
    log(f"base: {model.n_trials} trials, {time.time() - t0:.1f} s wall, median trial "
        f"{np.median(model.seconds):.1f} s; Phi_data = {float(r @ r):.6g} over {len(r)} rows")
    # G2: the same parameters give byte-identical output
    tdir = T.build_trial(site.base_main, site.base_pft, params, theta, cal[0], model.states[cal[0].name],
                         work / "repeat", site.overrides)
    from pool import Task
    pool.run([Task("repeat", T.command(args.meds_main, tdir), str(tdir), str(tdir / "run.log"), 3600)])
    df2 = T.finish(tdir, params, theta, site.tower.utc_offset_h)
    df1 = T.load_series(T.build_trial(site.base_main, site.base_pft, params, theta, cal[0],
                                      model.states[cal[0].name], work / "trials", site.overrides))
    same = all(np.array_equal(df1[v].values, df2[v].values) for v in T.TRIAL_VARIABLES)
    log(f"G2 repeat trial byte-identical: {same}")
    # G1: the stand is identical at the start and the end of a trial
    gdir = T.build_trial(site.base_main, site.base_pft, params, theta, cal[0], model.states[cal[0].name],
                         work / "g1", site.overrides, write_state=True)
    pool.run([Task("g1", T.command(args.meds_main, gdir), str(gdir), str(gdir / "run.log"), 3600)])
    g1 = stand_unchanged(model.states[cal[0].name], sorted((gdir / "out").glob(f"{T.PREFIX}-S-*.nc"))[-1])
    log(f"G1 stand identical at the trial's start and end: {g1['pass']} {g1['differ']}")
    same = same and g1["pass"]
    for name, sc in R.target_scores([model.specs[w.name] for w in cal], r).items():
        log(f"  {name:10s} nrmse {sc['nrmse']:.3f}  n = {sc['n']}")
    pool.close()
    return 0 if same else 1


STAND_VARIABLES = ("pft", "nplant", "dbh", "height", "leaf_area", "leaf_carbon", "fineroot_carbon",
                   "wood_carbon", "overtopping_lai", "patch_area")


def stand_unchanged(state_a, state_b) -> dict:
    from netCDF4 import Dataset
    with T.NC_LOCK, Dataset(state_a) as a, Dataset(state_b) as b:
        differ = [v for v in STAND_VARIABLES if v in a.variables
                  and not np.array_equal(np.asarray(a[v][:]), np.asarray(b[v][:]))]
    return {"pass": not differ, "differ": differ}


def cmd_smoke(args):
    """CTest smoke: a short window restarted from a state, one Jacobian column that must move the
    output, and a repeated trial that must reproduce it byte for byte."""
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "smoke.log")
    pool = make_pool("local", 3)
    w0 = [w for w in site.windows if w.role == "cal"][0]
    w = T.Window("smoke", w0.start, args.days, "cal", w0.chain)
    site.windows = [w]
    site.chains = {w.chain: w.start - dt.timedelta(days=1)}
    model, params, theta, _, _ = setup_model(site, args, work, pool, log, windows=[w])
    key = args.key
    j = next(i for i, p in enumerate(params) if p.name == key)
    prob = F.Problem(model, [j], theta.copy())
    u0 = np.array([params[j].u0])
    r0 = prob.data([u0])[0]
    J, smooth, failed = F.jacobian(prob, u0, r0, log=log)
    moved = bool(np.any(J[:, 0])) and not failed
    log(f"Jacobian column for {key}: |J| = {np.linalg.norm(J[:, 0]):.4g}, smoothness {smooth[0]:.3f}")
    tdir = T.build_trial(site.base_main, site.base_pft, params, theta, w, model.states[w.name],
                         work / "repeat", site.overrides)
    from pool import Task
    pool.run([Task("repeat", T.command(args.meds_main, tdir), str(tdir), str(tdir / "run.log"), 3600)])
    df2 = T.finish(tdir, params, theta, site.tower.utc_offset_h)
    df1 = T.load_series(T.build_trial(site.base_main, site.base_pft, params, theta, w,
                                      model.states[w.name], work / "trials", site.overrides))
    same = all(np.array_equal(df1[v].values, df2[v].values) for v in T.TRIAL_VARIABLES)
    log(f"repeat byte-identical: {same}; Jacobian column moves the output: {moved}")
    pool.close()
    return 0 if (same and moved) else 1


# ------------------------------------------------------------------------------------------------
# fit
# ------------------------------------------------------------------------------------------------
def cmd_fit(args):
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "fit.log")
    pool = make_pool(args.pool, args.workers, work / "queue")
    fc = site.fitcfg
    cal = [w for w in site.windows if w.role == "cal"]
    val = [w for w in site.windows if w.role == "val"]
    t_start = time.time()
    log(f"variant {site.variant}: {len(cal)} calibration and {len(val)} validation windows")
    model, params, theta0, obs, fok = setup_model(site, args, work, pool, log)
    model.windows = cal
    states_default = dict(model.states)
    report = {"variant": site.variant, "site": str(site.path), "windows": {w.name: [str(w.start), w.role]
              for w in site.windows}, "closure": TW.closure_summary(obs)}

    # ----- screening from the first Jacobian at the default ------------------------------------
    everything = F.Problem(model, list(range(len(params))), theta0.copy())
    u_def = everything.u_prior
    r_def = everything.data([u_def])[0]
    if r_def is None:
        raise SystemExit("the default's trials failed; see the trial logs")
    log(f"default: Phi = {everything.cost(u_def, r_def):.6g} over {len(r_def)} rows; screening "
        f"{len(params)} keys ({2 * len(params) * len(cal)} trials)")
    J0, sm0, fl0 = F.jacobian(everything, u_def, r_def, log=log)
    keep, scr = F.screening(everything, J0, r_def, sm0, fl0, [model.specs[w.name] for w in cal],
                            max_free=int(fc.get("max_free", 20)))
    report["screening"] = scr
    if scr["dead"]:
        log(f"WARNING: exact-zero Jacobian columns (a harness bug or a dead key): {scr['dead']}")
    log(f"screening: fit {len(keep)} keys {scr['fitted']}; fixed {len(scr['fixed'])}")
    prob = F.Problem(model, keep, theta0.copy())
    jac_def = (J0[:, keep], sm0[keep], [])

    # ----- the fit from several starts ----------------------------------------------------------
    rng = np.random.default_rng(int(fc.get("seed", 1)))
    starts = [prob.u_prior]
    for _ in range(int(fc.get("starts", 3)) - 1):
        z = np.clip(rng.standard_normal(len(keep)), -1.5, 1.5)
        starts.append(prob.u_prior + z * SIGMA_U)
    results = [None] * len(starts)

    def one(i):
        results[i] = F.lm(prob, starts[i], int(fc.get("max_iter", 15)), float(fc.get("rtol", 1e-3)),
                          log=log, label=f"start {i}", jac0=jac_def if i == 0 else None)
    threads = [threading.Thread(target=one, args=(i,)) for i in range(len(starts))]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    costs = [r["cost"] for r in results if r is not None]
    best = min((r for r in results if r is not None), key=lambda r: r["cost"])
    spread = (max(costs) - min(costs)) / min(costs) if costs else float("nan")
    report["starts"] = [{"u_start": s.tolist(), "cost": r["cost"], "u": r["u"].tolist(),
                         "iterations": len(r["history"]) - 1} for s, r in zip(starts, results)]
    report["start_spread"] = spread
    log(f"starts: Phi = {[round(c, 3) for c in costs]} (spread {100 * spread:.1f} %)")

    # ----- the state refresh (§6.2) --------------------------------------------------------------
    theta_map = prob.theta(best["u"])
    theta_refresh = theta_map.copy()
    n_refresh = int(fc.get("refresh_iter", 5))
    if n_refresh > 0:
        log("refreshing the state chains with the MAP")
        model.states = run_chains(site, params, theta_map, work / "chains", pool, args.meds_main, log)
        best = F.lm(prob, best["u"], n_refresh, float(fc.get("rtol", 1e-3)), log=log, label="refresh")
        theta_map = prob.theta(best["u"])

    # ----- rough keys: held at the default, or set by a 1-D line search (opt-in) -----------------
    #       A rough key's response is a threshold, not a slope (at BCI the xylem vulnerability and the
    #       turgor-loss point): a crude search can move it to a corner of its range, where the
    #       five-year run broke its water budget (G7). So the default holds them.
    rough_rule = fc.get("rough_keys", "default")
    for name in (scr["rough"] if rough_rule == "line_search" else []):
        j = next(i for i, p in enumerate(params) if p.name == name)
        sub = F.Problem(model, [j], theta_map.copy())
        u1, c1 = F.line_search(sub, np.array([params[j].to_u(theta_map[j])]), 0)
        theta_map[j] = float(params[j].to_theta(u1[0]))
        log(f"rough key {name}: line search -> {theta_map[j]:.6g}")
    report["refresh_theta"] = {p.name: float(v) for p, v in zip(params, theta_refresh)}
    report["cost_default"] = everything.cost(u_def, r_def)
    report["scores_cal_default"] = R.target_scores([model.specs[x.name] for x in cal], r_def)
    post_fit(site, model, params, theta0, theta_map, keep, states_default, report, work, log)
    report["trials"] = {"run": model.n_trials, "failed": model.n_failed,
                        "median_s": float(np.median(model.seconds_ok)) if model.seconds_ok else None,
                        "wall_s": time.time() - t_start}
    F.save_json(work / "fit.json", report)
    log(f"done: Phi {report['cost']['default']:.6g} -> {report['cost']['map']:.6g}; "
        f"{model.n_trials} trials ({model.n_failed} failed) in {(time.time() - t_start) / 60:.1f} min")
    pool.close()
    return 0


def post_fit(site, model, params, theta0, theta_map, keep, states_default, report, work, log):
    """Everything after the fit: the covariance and linearity check at the MAP (on the current,
    MAP-refreshed states), the scores, the validation (each set on its own chain's states), the
    gates, and the calibrated configs."""
    cal = [w for w in site.windows if w.role == "cal"]
    val = [w for w in site.windows if w.role == "val"]
    cspecs = [model.specs[x.name] for x in cal]
    prob = F.Problem(model, keep, theta_map.copy())
    u_map = np.array([params[i].to_u(theta_map[i]) for i in keep])
    r_map = prob.data([u_map])[0]
    J_map, _, _ = F.jacobian(prob, u_map, r_map, log=log)
    cov_u, w = F.posterior(prob, J_map, r_map, cspecs)
    cov_t = F.theta_cov(prob, u_map, cov_u)
    report["map"] = {p.name: float(theta_map[i]) for i, p in enumerate(params)}
    report["default"] = {p.name: float(theta0[i]) for i, p in enumerate(params)}
    report["fitted"] = [p.name for p in prob.params]
    report["cov_theta"] = cov_t
    report["corr"] = F.correlation(cov_u)
    report["sigma_ratio"] = {p.name: float(math.sqrt(cov_u[k, k]) / SIGMA_U) for k, p in enumerate(prob.params)}
    report["linearity"] = F.linearity(prob, u_map, cov_u, J_map, log=log)
    report["bounds"] = F.bound_pushers(prob, u_map, J_map, r_map, cspecs)
    report["cost"] = {"default": report.get("cost_default"), "map": prob.cost(u_map, r_map)}
    report["scores_cal"] = {"default": report.get("scores_cal_default"), "map": R.target_scores(cspecs, r_map)}
    if val:
        states_map = dict(model.states)
        model.states = states_default
        rv0 = model.residuals([theta0], windows=val)[0]
        model.states = states_map
        rv1 = model.residuals([theta_map], windows=val)[0]
        vspecs = [model.specs[x.name] for x in val]
        if rv0 is not None and rv1 is not None:
            report["scores_val"] = {"default": R.target_scores(vspecs, rv0), "map": R.target_scores(vspecs, rv1)}
            report["cost_val"] = {"default": float(rv0 @ rv0), "map": float(rv1 @ rv1)}
    gates = {}
    if "scores_val" in report:
        sv = report["scores_val"]
        worse = {k: sv["map"][k]["nrmse"] / sv["default"][k]["nrmse"] - 1.0 for k in sv["map"]}
        gates["G4"] = {"pass": report["cost_val"]["map"] < report["cost_val"]["default"]
                       and all(v <= 0.10 for v in worse.values()), "nrmse_change": worse}
    else:
        gates["G4"] = {"pass": False, "note": "no validation scores"}
    gates["G5"] = {"pass": True, "near_bound": report["bounds"]}
    scr = report.get("screening", {})
    gates["G3"] = {"pass": not scr.get("dead"), "dead": scr.get("dead"), "rough": scr.get("rough")}
    gates["G6"] = {"pass": report.get("start_spread", 1.0) <= 0.05, "spread": report.get("start_spread")}
    report["local_only"] = any(d["local_only"] for d in report["linearity"])
    report["gates"] = gates
    write_calibrated(site, params, theta_map, work)


def cmd_analyze(args):
    """Redo the post-fit analysis of a finished fit from its fit.json: the states are rebuilt (or
    found in the cache) for the default and for the parameters the fit refreshed them with."""
    import json
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    log = log_to(work / "analyze.log")
    pool = make_pool(args.pool, args.workers, work / "queue")
    report = json.loads((work / "fit.json").read_text())
    model, params, theta0, obs, fok = setup_model(site, args, work, pool, log)
    states_default = dict(model.states)
    names = [p.name for p in params]
    keep = [names.index(k) for k in report["fitted"]]
    theta_map = np.array([report["map"][k] for k in names])
    if site.fitcfg.get("rough_keys", "default") != "line_search":
        for name in report["screening"].get("rough", []):
            theta_map[names.index(name)] = theta0[names.index(name)]
            log(f"rough key {name} held at its default ({theta0[names.index(name)]:.6g})")
    if "refresh_theta" in report:
        theta_refresh = np.array([report["refresh_theta"][k] for k in names])
    else:                                    # a fit.json from before refresh_theta was recorded
        best = min(report["starts"], key=lambda s: s["cost"])
        theta_refresh = theta0.copy()
        for k, uk in zip(keep, best["u"]):
            theta_refresh[k] = float(params[k].to_theta(uk))
    if site.fitcfg.get("refresh_iter", 5) > 0:
        model.states = run_chains(site, params, theta_refresh, work / "chains", pool, args.meds_main, log)
    model.windows = [w for w in site.windows if w.role == "cal"]
    if "cost_default" not in report:           # the default, on the default chain's states
        everything = F.Problem(model, list(range(len(params))), theta0.copy())
        model.states, st = states_default, model.states
        r_def = model.residuals([theta0], windows=model.windows)[0]
        model.states = st
        report["cost_default"] = everything.cost(everything.u_prior, r_def)
        report["scores_cal_default"] = R.target_scores([model.specs[x.name] for x in model.windows], r_def)
    post_fit(site, model, params, theta0, theta_map, keep, states_default, report, work, log)
    F.save_json(work / "fit.json", report)
    gates = {k: v["pass"] for k, v in report["gates"].items()}
    log(f"analyzed: Phi {report['cost']['default']:.6g} -> {report['cost']['map']:.6g}; validation "
        f"{report.get('cost_val')}; gates {gates}")
    pool.close()
    return 0


# ------------------------------------------------------------------------------------------------
# the calibrated configs: the base files with the MAP values written into them, comments kept
# ------------------------------------------------------------------------------------------------
def _fmt_num(x) -> str:
    if isinstance(x, str):
        return '"' + x + '"'
    if isinstance(x, bool):
        return "true" if x else "false"
    return f"{x:.6g}"


def set_toml_text(text: str, key: str, value, index=None) -> str:
    """Set `key` (dotted: section.name) in TOML text, keeping every other line and comment. An array
    element is replaced by index (0-based); a missing key is added at its section's end."""
    section, name = key.rsplit(".", 1)
    lines = text.splitlines()
    cur, sec_end, found = "", None, False
    pat = re.compile(rf"^(\s*){re.escape(name)}(\s*=\s*)(\[[^\]]*\]|[^#\s]+)(.*)$")
    for i, line in enumerate(lines):
        m = re.match(r"^\s*\[([^\]]+)\]\s*(#.*)?$", line)
        if m:
            if cur == section:
                sec_end = i
            cur = m.group(1).strip()
            continue
        if cur == section:
            mm = pat.match(line)
            if mm:
                if index is None:
                    new = _fmt_num(value)
                else:
                    vals = [v.strip() for v in mm.group(3).strip("[]").split(",")]
                    vals[index] = _fmt_num(value)
                    new = "[" + ", ".join(vals) + "]"
                lines[i] = f"{mm.group(1)}{name}{mm.group(2)}{new}{mm.group(4)}"
                found = True
                break
    if not found:
        entry = f"{name} = {'[' + _fmt_num(value) + ']' if index is not None else _fmt_num(value)}   # set by calibrate_fast"
        if cur == section and sec_end is None:
            sec_end = len(lines)
        if sec_end is None:
            lines += ["", f"[{section}]", entry]
        else:
            lines.insert(sec_end, entry)
    return "\n".join(lines) + "\n"


def calibrated_header(site: Site, base: Path, changed: list, settings: dict | None = None) -> str:
    """The comment block that opens a calibrated file: where it came from and what the fit set."""
    lines = [f"# CALIBRATED by scripts/calibrate_fast (variant {site.variant or 'none'}) from {base.name}.",
             "# Every line is the base file's, comments included, except these keys the fit set",
             "# (calibrated value, the base file's value):"]
    lines += [f"#   {key:30s} {v:.6g}  ({d:.6g})" for key, v, d in changed] or ["#   (none)"]
    if settings:
        lines.append("# and these settings of the calibration (calibration.toml [overrides], the variant, [calibrated]):")
        lines += [f"#   {key:30s} {json.dumps(v)}" for key, v in settings.items()]
    return "\n".join(lines) + "\n"


def with_header(text: str, header: str) -> str:
    """Put the header after a leading SPDX line, if there is one."""
    first, _, rest = text.partition("\n")
    if first.startswith("# SPDX-License-Identifier"):
        return first + "\n" + header + rest
    return header + text


def write_calibrated(site: Site, params, theta, work: Path):
    main_text = site.main_path.read_text()
    pft_text = site.pft_path.read_text()
    changed = {"pft": [], "main": []}
    for p, v in zip(params, theta):
        if abs(v - p.default) <= 1e-12 * max(1.0, abs(p.default)):
            continue
        if p.file == "pft":
            pft_text = set_toml_text(pft_text, p.key, v, index=p.pft - 1)
        else:
            main_text = set_toml_text(main_text, p.key, v)
        changed["pft" if p.file == "pft" else "main"].append((p.key, float(v), float(p.default)))
    #----- the declaration's overrides (the variant's among them: the calibrated set was fitted
    #      with them), then its [calibrated] keys, e.g. the calibrated PFT file and an output prefix
    settings = dict(site.overrides) | dict(site.decl.get("calibrated", {}))
    for k, v in settings.items():
        main_text = set_toml_text(main_text, k, v)
    (work / "pft_parameters_calibrated.toml").write_text(
        with_header(pft_text, calibrated_header(site, site.pft_path, changed["pft"])))
    (work / "meds_config_calibrated.toml").write_text(
        with_header(main_text, calibrated_header(site, site.main_path, changed["main"], settings)))


def cmd_write_calibrated(args):
    """The calibrated configs from a finished fit's fit.json, into --out."""
    import json
    site = Site(args.site, args.variant)
    rep = json.loads(Path(args.fit).read_text())
    params = load_registry(site.registry_path, site.variant)
    missing = [p.name for p in params if p.name not in rep["map"]]
    if missing:
        raise SystemExit(f"fit.json has no MAP value for {missing}")
    for p in params:
        p.default = rep["default"][p.name]
    theta = np.array([rep["map"][p.name] for p in params])
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    write_calibrated(site, params, theta, out)
    print(f"wrote {out / 'pft_parameters_calibrated.toml'} and {out / 'meds_config_calibrated.toml'}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, needs_exe=True):
        p.add_argument("--site", required=True, help="the site declaration (calibration.toml)")
        p.add_argument("--variant", default=None)
        p.add_argument("--work", required=True, help="the working directory (trials, states, results)")
        if needs_exe:
            p.add_argument("--meds-main", required=True)
            p.add_argument("--pool", choices=("local", "queue"), default="local")
            p.add_argument("--workers", type=int, default=os.cpu_count())
            p.add_argument("--keep-netcdf", action="store_true")
    p = sub.add_parser("select-windows")
    p.add_argument("--site", required=True)
    p.set_defaults(func=cmd_select_windows)
    p = sub.add_parser("growth-resp")
    p.add_argument("--daily", required=True, help="glob of a slow-loop run's daily files")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_growth_resp)
    for name, func in (("check", cmd_check), ("fit", cmd_fit), ("analyze", cmd_analyze)):
        p = sub.add_parser(name)
        common(p)
        p.set_defaults(func=func)
    p = sub.add_parser("smoke")
    common(p)
    p.add_argument("--days", type=int, default=3)
    p.add_argument("--key", default="stomatal_g1")
    p.set_defaults(func=cmd_smoke)
    p = sub.add_parser("write-calibrated")
    p.add_argument("--site", required=True)
    p.add_argument("--variant", default=None)
    p.add_argument("--fit", required=True, help="a finished fit's fit.json")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_write_calibrated)
    p = sub.add_parser("worker")
    p.add_argument("--queue", required=True)
    p.add_argument("--slots", type=int, default=os.cpu_count())
    p.set_defaults(func=lambda a: worker(a.queue, a.slots) or 0)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
