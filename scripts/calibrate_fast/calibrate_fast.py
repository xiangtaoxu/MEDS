#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""calibrate_fast.py -- fit MEDS's fast (sub-daily) parameters to a flux tower, with the stand frozen
at its initial structure (docs/dev_plans/MEDS_FAST_CALIBRATION_PLAN.md and its revision,
MEDS_FAST_CALIBRATION_REVISION_PLAN.md).

The fit is staged along the canopy's causal chain (stages.py): optics with the two-stream alone,
photosynthesis with a canopy of leaf solves, the energy partition with the coupled fast loop, the
water stress with seasonal runs, and a joint polish of every key. Each coupled trial is a frozen
run restarted from a shared state at its window's start. Every stage minimizes the stacked,
filtered, weighted residuals plus Gaussian priors; the Laplace covariance at the MAP gives a rough
uncertainty. Every setting of the site declaration is documented, with its default, in
site_reference.toml.

Commands (each reads the site declaration, e.g. examples/example_flux_tower_bci/calibration.toml):
  select-windows  pick the 10-day windows with the best tower coverage in each declared season
  report          the data report only: every target's rows through every filter, the u* plateau
  check           one trial per window at the start: runs, the parameter record, repeatability
  fit             the data report, the screening, the stages, covariance, validation, gates and
                  the calibrated configs (--stages to run some; --resume to continue from stages.json)
  analyze         redo the post-fit steps from a finished fit's fit.json
  worker          a node's worker for --pool queue (started inside the Slurm allocation)
  smoke           the CTest smoke test: one Jacobian column on a short window, a repeated trial,
                  and gate G8 on the window's canopy drivers (with libmeds)

Trials run through MEDS's Python API (`python -m meds.model`, which needs libmeds.so) unless
--runner names a meds_main executable; the two give the same output, bit for bit. The optics and
photosynthesis stages call libmeds directly (meds.canopy).

Usage:
  calibrate_fast.py fit --site calibration.toml --variant interception_off --work runs/cal \
      --pool local --workers 40
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
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
#----- the meds package of this source tree: meds.config here, meds.model in the trials this starts
sys.path.insert(0, str(ROOT / "python"))
os.environ["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT / "python"),
                                                        os.environ.get("PYTHONPATH")]))

import datarules as DR               # noqa: E402
import fit as F                      # noqa: E402
import obsmodels as OM               # noqa: E402
import priors as PR                  # noqa: E402
import residuals as R                # noqa: E402
import settings as SET               # noqa: E402
import stages as STG                 # noqa: E402
import states as S                   # noqa: E402
import tower as TW                   # noqa: E402
import trials as T                   # noqa: E402
from pool import Task, make_pool, worker   # noqa: E402
from registry import load_registry, resolve_defaults, select   # noqa: E402
from meds.config import RunConfig, load_toml, read_record        # noqa: E402

STAGE_ORDER = ("optics", "photosynthesis", "energy", "water", "polish")
KERNEL_STAGES = ("optics", "photosynthesis")     # stages with a model of their own; their keys can join "energy"


# ------------------------------------------------------------------------------------------------
# The site declaration
# ------------------------------------------------------------------------------------------------
class Site:
    def __init__(self, path, variant=None):
        self.path = Path(path).resolve()
        self.dir = self.path.parent
        raw = load_toml(self.path)
        if "pft" in raw.get("base", {}):
            raise SystemExit("[base].pft is not read: the PFT file is the one the main file names "
                             "in [init].pft_config")
        try:
            self.decl = SET.complete(raw)
        except ValueError as e:
            raise SystemExit(str(e))
        d = self.decl
        b = d["base"]
        #----- the base configuration: the main file and the PFT file it names, inputs made absolute
        self.base = RunConfig.load(self._p(b["main"]))
        self.main_path, self.pft_path = self.base.path, self.base.pft_path
        self.registry_path = self._p(b["registry"])
        self.variant = variant
        self.overrides = dict(d.get("overrides", {}))
        if variant is not None:
            if variant not in d.get("variants", {}):
                raise SystemExit(f"unknown variant '{variant}'; declared: {list(d.get('variants', {}))}")
            self.overrides.update(d["variants"][variant])
        #----- the tower: the facts from its site TOML (the forcing builder's), the choices from [tower]
        tw = d["tower"]
        self.tower = TW.TowerSpec.from_site(
            self._p(tw["site"]), forcing=str(self._p(tw["forcing"])) if tw.get("forcing") else None,
            forcing_qc=tuple(tw["forcing_qc"]), forcing_qc_val=tuple(tw["forcing_qc_val"]),
            forcing_grid=int(tw["forcing_grid"]), daytime_sw=float(tw["daytime_sw"]), closure=dict(d["closure"]))
        #----- the trials write their fast output at the tower's own interval (30 min at BCI)
        dt_fast = seconds(self.base.get("fast.dt_fast"))
        steps = self.tower.step / dt_fast
        if abs(steps - round(steps)) > 1e-9 or round(steps) < 1:
            raise SystemExit(f"the tower's interval ({self.tower.step:g} s) is not a whole number of the model's "
                             f"fast steps (fast.dt_fast = {dt_fast:g} s)")
        self.base.set("output.fast_interval_steps", int(round(steps)))
        #----- the windows: the site's list, else chosen by the rule when the data are read
        #      (load_data). Every window has its own state chain, a frozen run from the initial
        #      stand that starts chain_lead_days before it.
        w = d["windows"]
        self.days = int(w["days"])
        self.skip_hours = int(w["skip_hours"])
        self.chain_lead = int(w["chain_lead_days"])
        self.windows = [T.Window(e["name"], dt.datetime.fromisoformat(e["start"]),
                                 int(e.get("days", self.days)), e["role"], e["name"])
                        for e in w["list"]]
        self.water_windows = [T.Window(e["name"], dt.datetime.fromisoformat(e["start"]),
                                       int(e.get("days", w["seasonal"]["days"])), "water", e["name"])
                              for e in w["seasonal"]["list"]]
        self.derived = {}                     # (file, key, pft) -> value: the Kattge & Knorr shape keys
        self.derived_was = {}                 # ... and the base config's value they replaced
        self.targets_declared = d["targets"]
        self.targets = d["targets"]           # with the u* rules made numbers by load_data
        self.data = None                      # load_data's cache
        self.fitcfg = d["fit"]
        self.stagecfg = d["stages"]
        self.uncertainty = d["uncertainty"]
        self.lat, self.lon = self.tower.lat, self.tower.lon

    def _p(self, rel) -> Path:
        return (self.dir / rel).resolve()

    def menu(self):
        """Every registry key of this variant, with its state."""
        return load_registry(self.registry_path, self.variant)

    def params(self, record=None):
        """The keys this fit moves: the registry's default set, or the site's [fit] list, with the
        site's [priors] over the registry's."""
        try:
            ps = select(self.menu(), self.fitcfg, self.decl.get("priors", {}))
        except ValueError as e:
            raise SystemExit(str(e))
        #----- kappa's sd from the gap between the provider's two partitionings, unless the site sets it
        site_kappa = self.decl.get("priors", {}).get("kappa", {})
        sd = (self.data or {}).get("kappa_sd")
        for p in ps:
            if p.file == "obs" and p.key == "kappa" and sd is not None and not ({"sd", "log_sd"} & set(site_kappa)):
                p.prior = {**p.prior, "sd": sd, "source": "the gap between the provider's night-time and daytime "
                                                          "partitionings (RECO against RECO_DT)"}
                p.prior.pop("log_sd", None)
        climate_priors(self, ps, record)
        return resolve_defaults(ps, self.base, record)

    def setting(self, key, default=None, record=None, file="main", pft=None):
        """A model setting: the base config's, else the base run's parameter record, else default."""
        v = self.base.get(key, file=file, pft=pft)
        if v is None and record is not None:
            hit = record.get((file, key, pft or 0)) or record.get((file, key, 0))
            v = None if hit is None else hit[1]
        return default if v is None else v

    def obs_fixed(self, params) -> dict:
        """The observation keys this fit does not move, at their prior centres (kappa enters the GPP
        residual whether or not it is fitted)."""
        fitted = {p.name for p in params}
        return {p.key: float(p.prior.get("centre", 1.0)) for p in self.menu() if p.file == "obs" and p.name not in fitted}

    def loss(self):
        return (self.fitcfg["loss"], float(self.fitcfg["huber_c"]))


def climate_priors(site, params, record):
    """The priors from the site's climate (priors.py): the Kattge & Knorr shape keys set in the base
    config (so in every trial and chain), and the EEO centres of stomatal_g1 and vcmax25. The values
    go to site.data["climate_priors"]."""
    menu = site.menu()
    fitted = {p.name for p in params}
    need_kk = [p for p in menu if p.fixed_at == "kattge_knorr" and p.name not in fitted]
    eeo_keys = [p for p in params if p.prior.get("centre") == "eeo"]
    if not need_kk and not eeo_keys:
        return
    if not site.tower.forcing:
        raise SystemExit("the climate priors (Kattge & Knorr, EEO) need the forcing file ([tower].forcing), or set "
                         "the keys' prior centres in [priors]")
    clim = PR.growth_climate(site.tower.forcing, site.tower.forcing_grid, site.tower.leaf_on_months,
                             site.tower.daytime_sw)
    out = {"climate": clim}
    value = lambda key, default: site.setting(key, default, record)          # noqa: E731
    if need_kk:
        if value("leaf_physiology.thermal_acclimation", False):
            out["kattge_knorr"] = "the model's own thermal acclimation sets them (leaf_physiology.thermal_acclimation)"
        else:
            kk = PR.kattge_knorr(clim["t_growth_c"])
            for p in need_kk:
                where = (p.file, p.key, p.pft if p.file == "pft" else None)
                site.derived_was.setdefault(where, float(site.setting(p.key, float("nan"), record, p.file, where[2])))
                site.base.set(p.key, float(kk[p.key]), file=p.file, pft=where[2])
                site.derived[where] = float(kk[p.key])
            out["kattge_knorr"] = {p.name: float(kk[p.key]) for p in need_kk}
    if eeo_keys:
        lp = PR.leaf_settings(value)
        ca = PR.co2_ppm(site.base, clim["years"])
        out["co2_ppm"] = ca
        for p in eeo_keys:
            if p.key == "pft.stomatal_g1":
                v = PR.eeo_g1(lp, clim)
            else:
                pft = {k: float(site.setting(f"pft.{k}", None, record, file="pft", pft=p.pft))
                       for k in ("theta_j", "jmax_vcmax_ratio")}
                v = PR.eeo_vcmax25(lp, pft, clim, ca, str(value("leaf_physiology.temp_response_form", "peaked")))
            p.prior = {**p.prior, "centre": float(v), "source": "EEO: " + p.prior.get("source", "")}
            entry = {"eeo": float(v)}
            pt = site.fitcfg.get("plant_type", "")
            meta = p.meta.get(pt)
            if meta:
                sd = meta.get("log_sd") or (meta.get("sd", 0.0) / meta["centre"])
                entry["meta"] = meta
                entry["meta_apart_sd"] = float(abs(math.log(v / meta["centre"])) / sd) if sd else None
            out.setdefault("eeo", {})[p.name] = entry
    if site.data is not None:
        site.data["climate_priors"] = out


def seconds(duration) -> float:
    """A model duration ("900s", "15min", "1h", or a number of seconds) in seconds."""
    m = re.fullmatch(r"\s*([0-9.]+)\s*(s|min|h)?\s*", str(duration))
    if not m:
        raise SystemExit(f"cannot read the duration {duration!r}")
    return float(m.group(1)) * {"s": 1.0, None: 1.0, "min": 60.0, "h": 3600.0}[m.group(2)]


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
    """The tower's records, the observed-forcing masks for calibration and validation windows, and
    the sun's elevation at every record (None without the site's coordinates); with the data rules
    applied once (load_data)."""
    if site.data is None:
        load_data(site)
    return site.data["obs"], site.data["fok"], site.data["elev"]


def load_data(site: Site) -> dict:
    """Read the tower and apply the data rules (datarules.py): each target's u* rule made a number,
    the windows chosen where the site lists none, and the seasonal runs from the water deficit."""
    obs = TW.observations(site.tower)
    fok = {"cal": TW.forcing_observed(site.tower, obs.index),
           "val": TW.forcing_observed(site.tower, obs.index, site.tower.forcing_qc_val)}
    fok["water"] = fok["cal"]
    elev = None
    if site.lat is not None and site.lon is not None:
        elev = pd.Series(TW.solar_elevation(obs.index, float(site.lat), float(site.lon), site.tower.step),
                         index=obs.index)
    tw = site.tower
    try:
        site.targets, ustar = DR.resolve_ustar(site.targets_declared, obs, DR.provider_threshold(tw.provider),
                                               site.decl["ustar"], tw.daytime_sw, tw.utc_offset_h)
    except ValueError as e:
        raise SystemExit(str(e))
    #----- the observation models (obsmodels.py): the closure shares, then each target's sigma
    ccfg = site.decl["closure"]
    test = OM.attribution(obs, ccfg, float(site.decl["ustar"]["min_driver"]["rnet"]), tw.daytime_sw)
    try:
        s_h, s_le, why = OM.closure_shares(ccfg["shares"], test, float(ccfg["rise_min"]))
        obs["h_c"], obs["le_c"] = OM.corrected(obs, s_h)
        sigma = OM.set_sigmas(site.targets, obs, site.decl["sigma"], tw.utc_offset_h, tw.step)
    except ValueError as e:
        raise SystemExit(str(e))
    closure = {**tw.report.get("closure", {}), "attribution": test, "shares": {"h": s_h, "le": s_le}, "reason": why}
    span = TW.forcing_span(tw)
    record_start = span[0] if span else obs.index[0]
    first = pd.Timestamp(record_start).floor("D") + pd.Timedelta(days=site.chain_lead)
    wcfg = site.decl["windows"]
    windows_report = {"source": "the site's list"}
    if not site.windows:
        chosen, windows_report = DR.select_windows(obs, fok, wcfg, tw.leaf_on_months, record_start,
                                                   tw.daytime_sw, tw.utc_offset_h)
        windows_report["source"] = "the rule (datarules.select_windows)"
        site.windows = [T.Window(n, start, site.days, role, n) for n, start, role in chosen]
        if not any(w.role == "cal" for w in site.windows):
            raise SystemExit("no calibration window passes the rule: " + json.dumps(windows_report["slots"]))
    deficit, seasonal_report = None, {"source": "the site's list"}
    if span:
        deficit = DR.water_deficit(DR.forcing_daily(tw.forcing, tw.forcing_grid))
    scfg = wcfg["seasonal"]
    if not site.water_windows and int(scfg["max_runs"]) > 0:
        if deficit is None:
            seasonal_report = {"note": "no forcing file ([tower].forcing): no water-deficit index, no seasonal runs"}
        else:
            scores = DR.day_scores(obs, fok["cal"], tw.daytime_sw, tw.utc_offset_h)
            runs, seasonal_report = DR.seasonal_runs(deficit, scfg, tw.leaf_on_months, first, scores,
                                                     float(wcfg["min_score"]))
            site.water_windows = [T.Window(n, start, days, "water", n) for n, start, days in runs]
    site.data = {"obs": obs, "fok": fok, "elev": elev, "ustar": ustar, "windows": windows_report,
                 "seasonal": seasonal_report, "deficit": deficit, "closure": closure, "sigma": sigma,
                 "kappa_sd": OM.kappa_sd(obs)}
    return site.data


def specs_for(site: Site, windows, obs, fok, elev, targets=None, log=print):
    out = {}
    for w in windows:
        idx = R.window_index(w.start, w.days, site.tower.step)
        e = None if elev is None else elev.reindex(idx).to_numpy()
        tcfg = site.targets if targets is None else {k: v for k, v in site.targets.items() if k in targets}
        out[w.name] = R.build_spec(w.name, idx, obs, fok[w.role], tcfg, site.skip_hours,
                                   site.tower.daytime_sw, e, site.tower.utc_offset_h)
        out[w.name].loss = site.loss()
        if out[w.name].n == 0:
            log(f"WARNING: window {w.name} has no usable records; it contributes nothing")
    return out


def run_chains(site, params, theta, root, pool, runner, log, windows=None):
    """Each window's state: its own frozen chain from the initial stand, chain_lead_days long,
    every chain at once."""
    windows = site.windows if windows is None else windows
    states, errs, threads = {}, [], []

    def one(w):
        try:
            start = w.start - dt.timedelta(days=site.chain_lead)
            states.update(S.run_chain(w.name, start, [w], site.base, params, theta, root,
                                      pool, runner, site.overrides))
        except Exception as e:           # noqa: BLE001 -- reported below
            errs.append(e)
    log(f"state chains: {len(windows)} windows, each from the initial stand {site.chain_lead} days before it")
    for w in windows:
        th = threading.Thread(target=one, args=(w,))
        th.start()
        threads.append(th)
    for th in threads:
        th.join()
    if errs:
        raise errs[0]
    return states


def base_record(site, work, pool, runner, log):
    """Run one short trial at the base configs to read the parameter record: the defaults of keys
    the base TOML leaves to the model."""
    rdir = work / "record"
    rec_csv = rdir / "out" / f"{T.PREFIX}_parameters.csv"
    if not rec_csv.exists():
        cfg = T.with_params(site.base, [], [], site.overrides)
        w = site.windows[0]
        for k, v in {"run.start_time": T.stamp(w.start), "run.end_time": T.stamp(w.start + dt.timedelta(days=1)),
                     "run.slow_on": False, "run.n_threads": 1,
                     "state.write_state": False, "output.enabled": True, "output.dir": str(rdir / "out"),
                     "output.prefix": T.PREFIX, "output.daily.enabled": False, "output.fast.enabled": False,
                     "output.monthly.enabled": False, "output.annual.enabled": False}.items():
            cfg.set(k, v)
        (rdir / "out").mkdir(parents=True, exist_ok=True)
        main = cfg.write(rdir)
        log("reading the parameter record from a one-day base run")
        pool.run([Task("record", T.command(runner, main), str(rdir), str(rdir / "run.log"), 3600)])
    return read_record(rec_csv, {rdir / "main.toml": "main", rdir / "pft.toml": "pft"})


# ------------------------------------------------------------------------------------------------
# select-windows, report
# ------------------------------------------------------------------------------------------------
def cmd_select_windows(args):
    """Print the windows and seasonal runs the rule chooses, as calibration.toml entries (to freeze
    a choice; a site that lists none gets the same choice at every run)."""
    site = Site(args.site)
    site.windows, site.water_windows = [], []
    load_data(site)
    for w in site.windows:
        print(f'[[windows.list]]\nname = "{w.name}"\nstart = "{w.start:%Y-%m-%d}"\nrole = "{w.role}"\n')
    for w in site.water_windows:
        print(f'[[windows.seasonal.list]]\nname = "{w.name}"\nstart = "{w.start:%Y-%m-%d}"\ndays = {w.days}\n')
    print("# " + json.dumps(site.data["windows"]) + "\n# " + json.dumps(site.data["seasonal"]))


def data_report(site, obs, fok, elev, log, specs=None) -> dict:
    """The data report (best-practice plan §2, §6.2): every target's rows over the whole record
    through each step, the u* diagnostics and the rule each target got, the windows and seasonal
    runs, what the calibration windows cover, and the reader's checks."""
    tw = site.tower
    rep = R.record_report(obs, fok["cal"], site.targets, tw.daytime_sw,
                          None if elev is None else elev.reindex(obs.index).to_numpy(), tw.utc_offset_h)
    for t, r in rep.items():
        log(f"data {t:9s} " + " -> ".join(f"{label} {n}" for label, n in r["counts"]))
    for t, u in site.data["ustar"].items():
        d = u["diagnostic"]
        if d is not None:
            b = d.get("bootstrap", {})
            log(f"data {t} u* diagnostic ({d.get('driver', '?')} classes): {d['outcome']}"
                + (f" at {d['threshold']}" if d.get("threshold") is not None else "")
                + (f", bootstrap 5-50-95 % {b['threshold_p05_p50_p95']}, outcomes {b['outcome_share']}"
                   if b.get("threshold_p05_p50_p95") else "") + (f" ({d['note']})" if d.get("note") else ""))
        log(f"data {t} u* filter: {u['ustar_min']} ({u['reason']})")
    for w in site.windows + site.water_windows:
        log(f"data window {w.name:12s} {w.role:5s} {w.start:%Y-%m-%d} + {w.days} d")
    sr = site.data["seasonal"]
    if sr.get("years"):
        log("data water deficit by year (deepest inside the leaf-on months): " + ", ".join(
            f"{y['year']} {y['deficit_mm']:.0f} mm on {y['end']}" for y in sr["years"]))
    if sr.get("note"):
        log(f"data seasonal runs: {sr['note']}")
    #----- what the calibration windows' kept records span of the record's conditions
    cal = [w for w in site.windows if w.role == "cal"]
    specs = specs or specs_for(site, cal + site.water_windows, obs, fok, elev, log=lambda *_: None)
    kept = {k: [] for k in ("sw_in", "tair", "vpd")}
    kept["deficit_mm"] = []
    deficit = site.data["deficit"]
    for w in cal + site.water_windows:
        spec = specs[w.name]
        rows = sorted({int(i) for t in spec.targets for i in t.rows})
        o = obs.reindex(spec.index).iloc[rows]
        for k in ("sw_in", "tair", "vpd"):
            kept[k] += o[k].tolist()
        if deficit is not None:
            kept["deficit_mm"] += deficit.reindex(pd.date_range(w.start, w.end, freq="D")).tolist()
    day = obs["sw_in"] > tw.daytime_sw
    record = {k: obs.loc[day, k].to_numpy() for k in ("sw_in", "tair", "vpd")}
    record["deficit_mm"] = deficit.to_numpy() if deficit is not None else []
    coverage = DR.range_coverage(record, kept)
    log("data the fit's records (calibration windows and seasonal runs) span of the record's daytime range: " + ", ".join(
        f"{k} {v['share']:.0%}" for k, v in coverage.items() if v))
    cl = site.data["closure"]
    if cl.get("f_median") is not None:
        log(f"data closure: daily (H + LE) / (Rnet - G) median {cl['daily_closure_median']:.2f} over {cl['valid_days']} "
            f"days, f median {cl['f_median']:.2f} (10-90 % {cl['f_p10_p90'][0]:.2f}-{cl['f_p10_p90'][1]:.2f})")
    if cl["attribution"].get("rises", {}).get("h") is not None:
        log("data closure attribution (rise from the calmest to the most turbulent third, by VPD class): "
            + "; ".join(f"{k.upper()} " + " ".join(f"{x:+.0%}" for x in v)
                        for k, v in cl["attribution"]["rises_by_vpd_class"].items()))
    log(f"data closure shares: H {cl['shares']['h']}, LE {cl['shares']['le']} ({cl['reason']})")
    for t, e in site.data["sigma"].items():
        log(f"data {t} sigma = {e['sigma_abs']:.3g} + {e['sigma_rel']:.3g} |{e['flux']}| ({e['source']}), "
            f"evaluated at the smoothed {e['flux']}")
    checks = tw.report.get("fluxes", {})
    for failure in checks.get("failures", []):
        log(f"data WARNING the tower's fluxes fail a check: {failure}")
    return {"filters": rep, "ustar": site.data["ustar"], "windows": site.data["windows"],
            "seasonal": sr, "range_coverage": coverage, "closure": cl, "sigma": site.data["sigma"],
            "tower_checks": checks}


def cmd_report(args):
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "report.log")
    load_data(site)
    obs, fok, elev = observations(site)
    F.save_json(work / "data_report.json", data_report(site, obs, fok, elev, log))
    return 0


# ------------------------------------------------------------------------------------------------
# the model, the start, check and smoke
# ------------------------------------------------------------------------------------------------
def setup_model(site, args, work, pool, log, windows=None, params=None, theta=None, states=None):
    obs, fok, elev = observations(site)
    windows = site.windows if windows is None else windows
    record = base_record(site, work, pool, args.runner, log)
    site.data["record"] = record
    params = site.params(record) if params is None else params
    theta = start_theta(params) if theta is None else theta
    if states is None:
        states = run_chains(site, params, theta, work / "chains", pool, args.runner, log, windows)
    specs = specs_for(site, windows, obs, fok, elev, log=log)
    model = F.Model(params, windows, specs, states, site.base, site.overrides,
                    args.runner, pool, work / "trials", site.tower.step,
                    timeout=float(site.fitcfg["timeout"]), keep_netcdf=args.keep_netcdf, log=log,
                    obs_fixed=site.obs_fixed(params))
    return model, params, theta, obs, fok, elev


def area_above_sensor(site, state_file, log) -> dict:
    """The share of the stand's area whose canopy-air top is above the tower's sensor (#350)."""
    rec = site.data.get("record") or {}
    def value(key, default):
        v = site.base.get(key)
        if v is None and ("main", key, 0) in rec:
            v = rec[("main", key, 0)][1]
        return float(default if v is None else v)
    out = DR.area_above_sensor(state_file, value("aerodynamics.canopy_freeboard", 5.0),
                               value("aerodynamics.min_canopy_depth", 5.0), site.tower.sensor_height)
    log(f"data canopy-air tops: {out['share_above']:.0%} of the area above the sensor at {out['sensor_height']:g} m "
        f"(area-weighted mean top {out['top_mean_area_weighted']:.1f} m, range {out['top_range'][0]:.1f}-"
        f"{out['top_range'][1]:.1f} m)")
    return out


def base_theta(params):
    """The base configuration's values."""
    return np.array([p.default for p in params])


def start_theta(params):
    """Where the fit starts: every key at its prior's centre."""
    return np.array([p.centre for p in params])


def cmd_check(args):
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "check.log")
    pool = make_pool(args.pool, args.workers, work / "queue")
    load_data(site)
    cal = [w for w in site.windows if w.role == "cal"]
    model, params, theta, _, _, _ = setup_model(site, args, work, pool, log, windows=cal)
    t0 = time.time()
    r = model.residuals([theta])[0]
    log(f"start: {model.n_trials} trials, {time.time() - t0:.1f} s wall, median trial "
        f"{np.median(model.seconds):.1f} s; Phi_data = {float(r @ r):.6g} over {len(r)} rows")
    # G2: the same parameters give byte-identical output
    tdir = T.build_trial(site.base, params, theta, cal[0], model.states[cal[0].name],
                         work / "repeat", site.overrides)
    pool.run([Task("repeat", T.command(args.runner, tdir / "main.toml"), str(tdir),
                   str(tdir / "run.log"), 3600)])
    df2 = T.finish(tdir, params, theta, site.tower.step)
    df1 = T.load_series(T.build_trial(site.base, params, theta, cal[0],
                                      model.states[cal[0].name], work / "trials", site.overrides))
    same = all(np.array_equal(df1[v].values, df2[v].values) for v in T.TRIAL_VARIABLES)
    log(f"G2 repeat trial byte-identical: {same}")
    # G1: the stand is identical at the start and the end of a trial
    gdir = T.build_trial(site.base, params, theta, cal[0], model.states[cal[0].name],
                         work / "g1", site.overrides, write_state=True)
    pool.run([Task("g1", T.command(args.runner, gdir / "main.toml"), str(gdir),
                   str(gdir / "run.log"), 3600)])
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
    output, a repeated trial that must reproduce it byte for byte, and -- when libmeds.so is
    there -- gate G8, the canopy of leaf solves against the window's own GPP."""
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "smoke.log")
    pool = make_pool("local", 3)
    load_data(site)
    w0 = [w for w in site.windows if w.role == "cal"][0]
    w = T.Window("smoke", w0.start, args.days, "cal", "smoke")
    site.windows, site.water_windows = [w], []
    site.chain_lead = 1
    model, params, theta, _, _, _ = setup_model(site, args, work, pool, log, windows=[w])
    key = args.key
    j = next(i for i, p in enumerate(params) if p.name == key)
    prob = F.Problem(model, [j], theta.copy())
    u0 = np.array([params[j].u0])
    r0 = prob.data([u0])[0]
    J, smooth, failed = F.jacobian(prob, u0, r0, log=log)
    moved = bool(np.all(np.isfinite(J[:, 0])) and np.any(J[:, 0] != 0.0)) and not failed
    log(f"Jacobian column for {key}: |J| = {np.linalg.norm(J[:, 0]):.4g}, smoothness {smooth[0]:.3f}")
    tdir = T.build_trial(site.base, params, theta, w, model.states[w.name],
                         work / "repeat", site.overrides)
    pool.run([Task("repeat", T.command(args.runner, tdir / "main.toml"), str(tdir),
                   str(tdir / "run.log"), 3600)])
    df2 = T.finish(tdir, params, theta, site.tower.step)
    df1 = T.load_series(T.build_trial(site.base, params, theta, w,
                                      model.states[w.name], work / "trials", site.overrides))
    same = all(np.array_equal(df1[v].values, df2[v].values) for v in T.TRIAL_VARIABLES)
    log(f"repeat byte-identical: {same}; Jacobian column moves the output: {moved}")
    g8_ok = True
    if args.g8:
        photo = photo_model(site, model, params, theta, [w], work, log)
        g8 = photo.g8(theta)
        log(f"G8 canopy of leaf solves vs the model's GPP: {g8}")
        g8_ok = g8["pass"]
    pool.close()
    return 0 if (same and moved and g8_ok) else 1


# ------------------------------------------------------------------------------------------------
# the stages
# ------------------------------------------------------------------------------------------------
def driver_data(model, theta, windows):
    """Run (or find) each window's driver trial at theta: {name: drivers}, {name: fast series}."""
    dirs = model.run([theta], windows, kind="driver")[0]
    if dirs is None:
        raise RuntimeError("a driver trial failed; see the trial logs")
    return ({w.name: T.load_drivers(td) for w, td in zip(windows, dirs)},
            {w.name: T.load_series(td) for w, td in zip(windows, dirs)})


def kernel_model(cls, site, model, params, theta, windows, target, work):
    drivers, series = driver_data(model, theta, windows)
    specs = {w.name: model.specs[w.name].subset([target]) for w in windows}
    km = cls(params, windows, specs, drivers, series,
             STG.config_writer(site.base, params, site.overrides, work / "kernel_configs"), model.log)
    km.set_anchor(theta)
    return km


def photo_model(site, model, params, theta, windows, work, log):
    return kernel_model(STG.PhotoModel, site, model, params, theta, windows, "gpp", work)


def free_of(params, stage, run_stages=STAGE_ORDER, skipped=(), seasonal=True):
    """The keys a stage moves:
      - optics, photosynthesis, water: their own;
      - energy (the coupled loop): its own, and those of a kernel stage (optics, photosynthesis)
        that is not run, which then join the coupled fit;
      - polish: every key, except those of a stage skipped for lack of data (its target off) and,
        when the polish has no seasonal runs, the water keys (ten-day windows reward a value that
        dries the soil months later).
    A key the drop-mode screening fixed is in none."""
    def moves(p):
        if stage == "polish":
            return p.stage not in skipped and (seasonal or p.stage != "water")
        if stage == "energy":
            return p.stage == "energy" or (p.stage in KERNEL_STAGES and p.stage not in run_stages)
        return p.stage == stage
    return [i for i, p in enumerate(params) if p.stage != "dropped" and moves(p)]


def fitted_keys(params, run_stages, skipped, seasonal):
    """The keys a whole fit moved: the polish's, or without one, every run stage's together."""
    if "polish" in run_stages:
        return free_of(params, "polish", run_stages, skipped, seasonal)
    keys = {i for s in run_stages for i in free_of(params, s, run_stages, skipped)}
    return sorted(i for i in keys if seasonal or params[i].stage != "water")


def run_stage(stage, site, model, params, theta, cal, water, work, log, report):
    """One stage from theta; returns the new theta (every key; only this stage's moved) and the
    stage's report entry."""
    sc = site.stagecfg.get(stage, {})
    free = free_of(params, stage, report["stages_run"])
    if not free:
        return theta, {"skipped": "no fitted key in this stage"}
    names = [params[i].name for i in free]
    entry = {"keys": names}
    t0 = time.time()
    if stage == "optics":
        if "albedo" not in site.targets or not site.targets["albedo"].get("on", True):
            return theta, {"skipped": "the albedo target is off"}
        km = kernel_model(STG.OpticsModel, site, model, params, theta, cal, "albedo", work)
        prob = F.Problem(km, free, theta.copy())
        out = F.lm(prob, prob.u_of(theta), int(sc["max_iter"]), 1e-4, log=log, label="optics")
        theta = prob.theta(out["u"])
        entry.update(cost=out["cost"], iterations=len(out["history"]) - 1, evaluations=km.n_evals)
    elif stage == "photosynthesis":
        if "gpp" not in site.targets or not site.targets["gpp"].get("on", True):
            return theta, {"skipped": "the GPP target is off"}
        passes = []
        for k in range(int(sc["passes"])):
            km = photo_model(site, model, params, theta, cal, work, log)
            if k == 0:
                entry["G8"] = km.g8(theta)
                log(f"G8 canopy of leaf solves vs the model's GPP (before the anchor): {entry['G8']}")
            prob = F.Problem(km, free, theta.copy())
            out = F.lm(prob, prob.u_of(theta), int(sc["max_iter"]), 1e-4, log=log, label=f"photosynthesis pass {k + 1}")
            new = prob.theta(out["u"])
            move = max(abs(new[i] / theta[i] - 1.0) for i in free)
            passes.append({"theta": {params[i].name: float(new[i]) for i in free}, "cost": out["cost"],
                           "largest_move": move, "evaluations": km.n_evals})
            log(f"photosynthesis pass {k + 1}: " + ", ".join(f"{params[i].name} {new[i]:.4g}" for i in free)
                + f" (largest move {100 * move:.2f} %)")
            theta = new
            if move < float(sc["tol"]):
                break
        entry["passes"] = passes
    elif stage == "energy":
        prob = F.Problem(model, free, theta.copy())
        out = F.lm(prob, prob.u_of(theta), int(sc["max_iter"]), float(sc["rtol"]), log=log, label="energy")
        theta = prob.theta(out["u"])
        entry.update(cost=out["cost"], iterations=len(out["history"]) - 1)
    elif stage == "water":
        if not water:
            return theta, {"skipped": "no seasonal windows ([stages.water].windows)"}
        model.windows = water
        prob = F.Problem(model, free, theta.copy())
        u, c, hist = STG.grid_search(prob, prob.u_of(theta), int(sc["rounds"]), float(sc["grid_sigma"]),
                                     log=log, label="water")
        theta = prob.theta(u)
        model.windows = cal
        entry.update(cost=c, history=hist)
    entry["theta"] = {p.name: float(theta[i]) for i, p in enumerate(params)}
    entry["wall_s"] = time.time() - t0
    log(f"stage {stage} done in {(time.time() - t0) / 60:.1f} min: "
        + ", ".join(f"{params[i].name} {theta[i]:.4g}" for i in free))
    return theta, entry


def ess_weights_at(site, model, theta, specs_list, windows, log, scale_sigma=False):
    """Set each target's effective-sample weight from the residuals at theta (unweighted); with
    scale_sigma, also scale each target's sigma by its misfit (residuals.scale_sigma). Returns the
    weights and the scales."""
    for s in specs_list:
        for t in s.targets:
            t.weight = 1.0
    r = model.residuals([theta], windows)[0]
    if r is None:
        raise RuntimeError("the weights' trials failed")
    w = R.set_weights(specs_list, r)
    log("effective-sample weights n_eff/n: " + ", ".join(
        f"{k} {np.mean([v for kk, v in w.items() if kk.endswith('/' + k)]):.2f}"
        for k in sorted({kk.split('/')[1] for kk in w})))
    scales = None
    if scale_sigma:
        for s in specs_list:                     # the scale is from the plain residuals: weights out
            for t in s.targets:
                t.weight_saved, t.weight = t.weight, 1.0
        scales = R.scale_sigma(specs_list, r, float(site.fitcfg["sigma_scale_max"]))
        for s in specs_list:
            for t in s.targets:
                t.weight = t.weight_saved
        log("sigma scaled by the model's misfit (max(1, sqrt(chi2 per row)), at most "
            f"{site.fitcfg['sigma_scale_max']}): " + ", ".join(f"{k} {v:.2f}" for k, v in scales.items()))
    return w, scales


def cmd_fit(args):
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / "fit.log")
    pool = make_pool(args.pool, args.workers, work / "queue")
    fc = site.fitcfg
    stages = [s.strip() for s in args.stages.split(",")] if args.stages else list(fc["stages"])
    bad = [s for s in stages if s not in STAGE_ORDER]
    if bad:
        raise SystemExit(f"unknown stages {bad}; known: {STAGE_ORDER}")
    load_data(site)
    cal = [w for w in site.windows if w.role == "cal"]
    val = [w for w in site.windows if w.role == "val"]
    water = site.water_windows
    t_start = time.time()
    log(f"variant {site.variant}: {len(cal)} calibration, {len(val)} validation and {len(water)} seasonal "
        f"windows; stages {stages}")
    model, params, theta0, obs, fok, elev = setup_model(site, args, work, pool, log,
                                                        windows=site.windows + water)
    model.windows = cal
    #----- the seasonal runs are scored on the water stage's targets only
    for w in water:
        model.specs[w.name] = model.specs[w.name].subset(site.stagecfg["water"]["targets"])
    theta_base = base_theta(params)
    states_start = dict(model.states)
    report = {"variant": site.variant, "site": str(site.path), "stages_run": stages,
              "windows": {w.name: [str(w.start), w.role, w.days] for w in site.windows + water},
              "menu": {p.name: {"state": p.state, "stage": p.stage, "kind": p.kind, "scope": p.scope,
                                "reason": p.reason} for p in site.menu()},
              "keys": {p.name: {"stage": p.stage, "kind": p.kind, "scope": p.scope, "range": [p.lo, p.hi],
                                "default": p.default, "prior_centre": p.centre, "prior_sigma_u": p.sigma_u,
                                "prior_source": p.prior.get("source", p.source)} for p in params},
              "climate_priors": site.data.get("climate_priors"),
              "derived": [[f, k, i, v] for (f, k, i), v in site.derived.items()]}
    cp = site.data.get("climate_priors") or {}
    if cp.get("kattge_knorr"):
        log(f"Kattge & Knorr at the growth temperature {cp['climate']['t_growth_c']:.1f} C: {cp['kattge_knorr']}")
    for k, e in (cp.get("eeo") or {}).items():
        log(f"EEO prior centre of {k}: {e['eeo']:.4g}" + (f" (the {site.fitcfg.get('plant_type')} meta-analysis: "
            f"{e['meta']['centre']:.4g}, {e['meta_apart_sd']:.1f} sd apart{' -- FLAGGED' if e['meta_apart_sd'] > 2 else ''})"
            if e.get("meta") else ""))
    report["data"] = data_report(site, obs, fok, elev, log, model.specs)
    #----- keys whose process the kept data never sample are fixed (best-practice plan §2.5)
    coverage = DR.process_coverage(model.specs, obs, cal + water, site.tower.daytime_sw)
    fixed = DR.fix_by_coverage(params, coverage, int(fc["min_process_records"]))
    for p in params:
        if p.name in fixed:
            p.stage = "dropped"
    report["process_coverage"], report["fixed_by_coverage"] = coverage, fixed
    log(f"kept records sampling each process: {coverage}")
    for k, why in fixed.items():
        log(f"fixed by coverage: {k}: {why}")
    report["data"]["area_above_sensor"] = area_above_sensor(site, model.states[cal[0].name], log)
    log("fitting " + ", ".join(f"{p.name} [{p.stage}]" for p in params))

    # ----- the weights, at the start --------------------------------------------------------------
    cspecs = [model.specs[w.name] for w in cal]
    wspecs = [model.specs[w.name] for w in water]
    if fc["weights"] == "ess":
        report["weights_start"], _ = ess_weights_at(site, model, theta0, cspecs, cal, log)
    elif fc["weights"] != "none":
        raise SystemExit('[fit].weights must be "ess" or "none"')

    # ----- the screening: the first Jacobian at the start, every fitted key --------------------
    everything = F.Problem(model, list(range(len(params))), theta0.copy())
    u_start = everything.u_prior
    r_start = everything.data([u_start])[0]
    if r_start is None:
        raise SystemExit("the start's trials failed; see the trial logs")
    log(f"start: Phi = {everything.cost(u_start, r_start):.6g} over {len(r_start)} rows; screening "
        f"{len(params)} keys ({2 * len(params) * len(cal)} trials)")
    J0, sm0, fl0 = F.jacobian(everything, u_start, r_start, log=log)
    keep, scr = F.screening(everything, J0, r_start, sm0, fl0, cspecs, max_free=int(fc["max_free"]),
                            mode=fc["screening"], weighted=fc["weights"] == "none")
    report["screening"] = scr
    if scr["dead"]:
        log(f"WARNING: exact-zero Jacobian columns (a harness bug or a dead key): {scr['dead']}")
    if scr["rough"]:
        log(f"WARNING: rough keys (a step response; the gradient cannot steer them): {scr['rough']}")
    log(f"screening ({fc['screening']}): fit {scr['fitted']}; the drop mode would fix {scr['would_drop']}")
    if fc["screening"] == "drop":
        for i in range(len(params)):                     # a dropped key stays at its prior's centre
            if i not in keep:
                params[i].stage = "dropped"
    report["cost_start"] = everything.cost(u_start, r_start)
    report["scores_cal_start"] = R.target_scores(cspecs, r_start)

    # ----- the stages -----------------------------------------------------------------------------
    stages_path = work / "stages.json"
    theta = theta0.copy()
    progress = {"stages": {}}
    if args.resume and stages_path.exists():
        progress = json.loads(stages_path.read_text())
        last = progress.get("theta")
        if last:
            theta = np.array([last.get(p.name, theta[i]) for i, p in enumerate(params)])
            log(f"resuming from stages.json after {list(progress['stages'])}")
    for stage in [s for s in stages if s != "polish"]:
        if args.resume and stage in progress["stages"]:
            continue
        theta, entry = run_stage(stage, site, model, params, theta, cal, water, work, log, report)
        if "skipped" in entry:
            log(f"stage {stage} skipped: {entry['skipped']}")
        progress["stages"][stage] = entry
        progress["theta"] = {p.name: float(theta[i]) for i, p in enumerate(params)}
        F.save_json(stages_path, progress)

    # ----- the polish: every key together, from the stage values --------------------------------
    pc = site.stagecfg["polish"]
    theta_polish = None
    if "polish" in stages:
        theta_polish = {p.name: float(theta[i]) for i, p in enumerate(params)}
        if pc["refresh"]:
            log("refreshing the state chains with the stage values")
            model.states = run_chains(site, params, theta, work / "chains", pool, args.runner, log,
                                      site.windows + water)
        pwin = cal + (water if pc["include_water"] else [])
        pspecs = [model.specs[w.name] for w in pwin]
        if fc["weights"] == "ess":
            report["weights_polish"], report["sigma_scale_refresh"] = ess_weights_at(
                site, model, theta, pspecs, pwin, log, scale_sigma=pc["refresh"])
        model.windows = pwin
        skipped = [s for s, e in progress["stages"].items() if "skipped" in e]
        prob = F.Problem(model, free_of(params, "polish", stages, skipped, seasonal=bool(water and pc["include_water"])),
                         theta.copy())
        out = F.lm(prob, prob.u_of(theta), int(pc["max_iter"]), float(pc["rtol"]), log=log, label="polish")
        theta = prob.theta(out["u"])
        progress["stages"]["polish"] = {"cost": out["cost"], "iterations": len(out["history"]) - 1,
                                        "theta": {p.name: float(theta[i]) for i, p in enumerate(params)}}
        progress["theta"] = progress["stages"]["polish"]["theta"]
        F.save_json(stages_path, progress)
        jac_final = (out["J"], out["r"])
    else:
        pwin = cal
        model.windows = pwin
        jac_final = None
    report["stages"] = progress["stages"]
    report["theta_stages"] = progress.get("theta")
    #----- where the polish started: its weights were set there, and its states (with refresh) ran there
    report["theta_polish_start"] = theta_polish
    report["chains_refreshed"] = bool(theta_polish and pc["refresh"])

    # ----- rough keys: held, or set by a 1-D line search (opt-in) --------------------------------
    if fc["rough_keys"] == "line_search":
        for name in scr["rough"]:
            j = next(i for i, p in enumerate(params) if p.name == name)
            sub = F.Problem(model, [j], theta.copy())
            u1, _ = F.line_search(sub, np.array([params[j].to_u(theta[j])]), 0)
            theta[j] = float(params[j].to_theta(u1[0]))
            log(f"rough key {name}: line search -> {theta[j]:.6g}")
    post_fit(site, model, params, theta_base, theta0, theta, pwin, states_start, report, work, log, jac_final)
    report["trials"] = {"run": model.n_trials, "failed": model.n_failed,
                        "median_s": float(np.median(model.seconds_ok)) if model.seconds_ok else None,
                        "wall_s": time.time() - t_start}
    F.save_json(work / "fit.json", report)
    cv = report.get("cost_val") or {}
    log(f"done: validation Phi {cv.get('default', float('nan')):.6g} (default) -> {cv.get('map', float('nan')):.6g} (MAP); "
        f"gates {({k: v['pass'] for k, v in report['gates'].items()})}; "
        f"{model.n_trials} trials ({model.n_failed} failed) in {(time.time() - t_start) / 60:.1f} min")
    pool.close()
    return 0


def post_fit(site, model, params, theta_base, theta0, theta_map, pwin, states_start, report, work, log,
             jac_final=None):
    """Everything after the stages: the covariance (sigma-scaled) and its intervals, the linearity
    check and the filter sensitivity at the MAP (on the current states), the scores, the validation
    (each set on its own chain's states), the gates, and the calibrated configs."""
    fc = site.fitcfg
    cal = [w for w in site.windows if w.role == "cal"]
    val = [w for w in site.windows if w.role == "val"]
    pspecs = [model.specs[x.name] for x in pwin]
    skipped = [s for s, e in report.get("stages", {}).items() if "skipped" in e]
    seasonal = any(w.role == "water" for w in pwin)
    prob = F.Problem(model, fitted_keys(params, report.get("stages_run", list(STAGE_ORDER)), skipped, seasonal),
                     theta_map.copy())
    u_map = prob.u_of(theta_map)
    model.windows = pwin
    if jac_final is not None:
        J_map, r_map = jac_final
    else:
        r_map = prob.data([u_map])[0]
        J_map, _, _ = F.jacobian(prob, u_map, r_map, log=log)
    chi2 = R.chi2_per_row(pspecs, r_map)
    report["chi2_per_row"] = chi2
    #----- sigma was scaled by the misfit at the refresh; without one, the covariance scales it here
    scales = None if report.get("sigma_scale_refresh") else chi2
    big = {k: round(v, 2) for k, v in chi2.items() if v >= 2.0}
    if big:
        log(f"targets with chi2 per row >= 2 at the MAP (structural misfit): {big}")
    cov_u, _ = F.posterior(prob, J_map, r_map, pspecs, weighted=fc["weights"] == "none", sigma_scale=scales)
    report["map"] = {p.name: float(theta_map[i]) for i, p in enumerate(params)}
    report["default"] = {p.name: float(theta_base[i]) for i, p in enumerate(params)}
    report["start"] = {p.name: float(theta0[i]) for i, p in enumerate(params)}
    report["fitted"] = [p.name for p in prob.params]
    report["cov_theta"] = F.theta_cov(prob, u_map, cov_u)
    report["corr"] = F.correlation(cov_u)
    report["intervals"] = F.intervals(prob, u_map, cov_u)
    report["sigma_ratio"] = {k: v["sigma_ratio"] for k, v in report["intervals"].items()}
    for k, v in report["intervals"].items():
        log(f"  {k:22s} {v['map']:10.4g}  68 % [{v['i68'][0]:.4g}, {v['i68'][1]:.4g}]  95 % "
            f"[{v['i95'][0]:.4g}, {v['i95'][1]:.4g}]  sd ratio {v['sigma_ratio']:.2f}  ({v['stage']})")
    #----- how far the fit stopped from its own optimum: the Gauss-Newton step left, in posterior sd
    rest = F.gn_step(prob, u_map, J_map, r_map) / np.sqrt(np.diag(cov_u))
    report["unconverged_sd"] = {p.name: float(v) for p, v in zip(prob.params, rest)}
    log("the Gauss-Newton step left at the MAP [posterior sd] (0 at a converged fit): " + ", ".join(
        f"{n} {v:+.2f}" for n, v in report["unconverged_sd"].items()))
    report["linearity"] = F.linearity(prob, u_map, cov_u, J_map, int(site.uncertainty["linearity_dirs"]), log=log)
    for d in report["linearity"]:
        c, s = d["curvature_ratio"], d["slope_ratio"]
        log(f"linearity along direction {d['direction']}: "
            + ("a trial failed" if c is None else f"curvature {c:.2f} x the quadratic's, slope left {s:+.2f}")
            + (" (local only)" if d["local_only"] else ""))
    report["bounds"] = F.bound_pushers(prob, u_map, J_map, r_map, pspecs)
    #----- each key's prior z at the MAP; a trait key more than 2 sd from its evidence is a question (§5.1)
    report["prior_z"] = F.prior_z(prob, u_map, J_map, r_map, pspecs)
    for k, v in report["prior_z"].items():
        flag = abs(v["z"]) > 2.0
        log(f"  prior z of {k:22s} {v['z']:+.2f} ({v['kind']}){' pushed by ' + str(v['pushed_by']) if flag else ''}")
    report["cost"] = {"start": report.get("cost_start"), "map": prob.cost(u_map, r_map)}
    report["scores_cal"] = {"start": report.get("scores_cal_start"), "map": R.target_scores(pspecs, r_map)}
    report["filter_sensitivity"] = filter_sensitivity(site, model, prob, u_map, pwin, cov_u, J_map, r_map, log)
    if val:
        states_map = dict(model.states)
        #----- the default is the base configuration, on chains run with it (the start's chains ran
        #      at the prior centres, which differ where a prior is centred off the base value)
        if np.allclose(theta_base, theta0):
            model.states = states_start
        else:
            log("running the state chains with the base configuration for the validation's default")
            model.states = run_chains(site, params, theta_base, work / "chains", model.pool, model.runner,
                                      log, val)
        rv0 = model.residuals([theta_base], windows=val)[0]
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
    gates["G3"] = {"pass": not scr.get("dead"), "dead": scr.get("dead"), "rough": scr.get("rough"),
                   "note": "reported; a requested key is fitted unless [fit].screening = 'drop'"}
    g8 = report.get("stages", {}).get("photosynthesis", {}).get("G8")
    gates["G8"] = g8 if g8 else {"pass": None, "note": "the photosynthesis stage did not run"}
    fs = report["filter_sensitivity"]
    gates["G10"] = ({"pass": all(abs(v["shift_sd"]) < 1.0 for v in fs["keys"].values()), **fs}
                    if fs.get("keys") else {"pass": None, "note": fs.get("note")})
    gates["G11"] = {"pass": True, "near_bound": report["bounds"],
                    "note": "every key near a bound is listed with the target that pushed it and its stage"}
    flagged = {k: v for k, v in report["prior_z"].items() if v["kind"] == "trait" and abs(v["z"]) > 2.0}
    gates["G13"] = {"pass": not flagged, "flagged": flagged,
                    "note": "a trait key more than 2 prior sd from its evidence: diagnose it (plan §5.1) or relabel it effective"}
    report["local_only"] = any(d["local_only"] for d in report["linearity"])
    report["gates"] = gates
    write_calibrated(site, params, theta_map, work)


def filter_sensitivity(site, model, prob, u_map, windows, cov_u, J_map, r_map, log) -> dict:
    """Revision plan §8.5: the MAP's linear response to the [uncertainty].alternative filters, from
    the cached trials of the final point and Jacobian (no new model runs, unless the cache lacks
    them): the Gauss-Newton step on the alternative rows minus the step on the fit's own rows."""
    alt = site.uncertainty.get("alternative") or {}
    if not alt:
        return {"note": "no [uncertainty].alternative"}
    obs, fok, elev = observations(site)
    targets = {k: dict(v) for k, v in site.targets.items()}
    for t, over in alt.items():
        targets.setdefault(t, {}).update(over)
    saved = model.specs
    alt_site_targets, site.targets = site.targets, targets
    try:
        alt_specs = specs_for(site, windows, obs, fok, elev, log=log)
    finally:
        site.targets = alt_site_targets
    #----- the same target weights as the fit's rows (one scalar per window and target)
    for w in windows:
        wmap = {t.name: t.weight for t in saved[w.name].targets}
        for t in alt_specs[w.name].targets:
            t.weight = wmap.get(t.name, 1.0)
    model.specs = {**saved, **alt_specs}
    try:
        r_alt = prob.data([u_map], windows)[0]
        J_alt, _, _ = F.jacobian(prob, u_map, r_alt, log=log)
    finally:
        model.specs = saved
    du = F.shift(prob, u_map, J_map, r_map, J_alt, r_alt)
    sd = np.sqrt(np.diag(cov_u))
    out = {"alternative": alt, "keys": {}}
    for k, p in enumerate(prob.params):
        th0, th1 = float(p.to_theta(u_map[k])), float(p.to_theta(u_map[k] + du[k]))
        out["keys"][p.name] = {"map": th0, "alternative": th1, "shift_sd": float(du[k] / sd[k])}
    log("filter sensitivity " + json.dumps(alt) + ": " + ", ".join(
        f"{n} {v['map']:.4g}->{v['alternative']:.4g} ({v['shift_sd']:+.2f} sd)" for n, v in out["keys"].items()))
    return out


def cmd_analyze(args):
    """Redo the post-fit analysis of a finished fit from its fit.json, on the fit's own states and
    weights: the states are rebuilt (or found in the cache) for the start and, with a refresh, for
    the polish's starting point, and the weights are set at that point."""
    site = Site(args.site, args.variant)
    work = Path(args.work).resolve()
    log = log_to(work / "analyze.log")
    pool = make_pool(args.pool, args.workers, work / "queue")
    report = json.loads((work / "fit.json").read_text())
    load_data(site)
    water = site.water_windows
    model, params, theta0, obs, fok, elev = setup_model(site, args, work, pool, log, windows=site.windows + water)
    states_start = dict(model.states)
    names = [p.name for p in params]
    theta_map = np.array([report["map"][k] for k in names])
    theta_base = base_theta(params)
    pc = site.stagecfg["polish"]
    theta_weights = theta0
    if report.get("theta_polish_start"):
        theta_weights = np.array([report["theta_polish_start"][k] for k in names])
        if report.get("chains_refreshed"):
            model.states = run_chains(site, params, theta_weights, work / "chains", pool, args.runner, log,
                                      site.windows + water)
    cal = [w for w in site.windows if w.role == "cal"]
    pwin = cal + (water if pc["include_water"] and "polish" in report.get("stages_run", []) else [])
    if site.fitcfg["weights"] == "ess":
        _, report["sigma_scale_refresh"] = ess_weights_at(site, model, theta_weights, [model.specs[w.name] for w in pwin],
                                                          pwin, log, scale_sigma=bool(report.get("sigma_scale_refresh")))
    post_fit(site, model, params, theta_base, theta0, theta_map, pwin, states_start, report, work, log)
    F.save_json(work / "fit.json", report)
    gates = {k: v["pass"] for k, v in report["gates"].items()}
    cv = report.get("cost_val") or {}
    log(f"analyzed: validation Phi {cv.get('default', float('nan')):.6g} (default) -> "
        f"{cv.get('map', float('nan')):.6g} (MAP); gates {gates}")
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
    """The comment block that opens a calibrated file: where it came from and what the fit set. An
    effective key (a scheme property, not a measurable trait) is labelled so."""
    lines = [f"# CALIBRATED by scripts/calibrate_fast (variant {site.variant or 'none'}) from {base.name}.",
             "# Every line is the base file's, comments included, except these keys the fit set",
             "# (calibrated value, the base file's value; EFFECTIVE: its value belongs to this model structure;",
             "# KATTGE & KNORR: fixed at the growth temperature's value, not fitted):"]
    lines += [f"#   {key:30s} {v:.6g}  ({d:.6g}){'  ' + label if label else ''}"
              for key, v, d, label in changed] or ["#   (none)"]
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


def write_calibrated(site: Site, params, theta, work: Path, derived=None):
    """The base main and PFT files with the MAP (and the Kattge & Knorr shape keys) written in."""
    main_text = site.main_path.read_text()
    pft_text = site.pft_path.read_text()
    changed = {"pft": [], "main": []}
    for p, v in zip(params, theta):
        if p.file == "obs" or abs(v - p.default) <= 1e-12 * max(1.0, abs(p.default)):
            continue
        if p.file == "pft":
            pft_text = set_toml_text(pft_text, p.key, v, index=p.pft - 1)
        else:
            main_text = set_toml_text(main_text, p.key, v)
        changed["pft" if p.file == "pft" else "main"].append((p.key, float(v), float(p.default),
                                                              "EFFECTIVE" if p.kind == "effective" else ""))
    for (file, key, pft), v in (site.derived if derived is None else derived).items():
        old = site.derived_was.get((file, key, pft), site.setting(key, float("nan"), file=file, pft=pft))
        if file == "pft":
            pft_text = set_toml_text(pft_text, key, v, index=pft - 1)
        else:
            main_text = set_toml_text(main_text, key, v)
        changed[file].append((key, float(v), float(old), "KATTGE & KNORR"))
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
    site = Site(args.site, args.variant)
    rep = json.loads(Path(args.fit).read_text())
    params = [p for p in site.menu() if p.name in rep["map"]]
    for p in params:
        p.default = rep["default"][p.name]
    theta = np.array([rep["map"][p.name] for p in params])
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    write_calibrated(site, params, theta, out, {(f, k, i): v for f, k, i, v in rep.get("derived", [])})
    print(f"wrote {out / 'pft_parameters_calibrated.toml'} and {out / 'meds_config_calibrated.toml'}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--site", required=True, help="the site declaration (calibration.toml)")
        p.add_argument("--variant", default=None)
        p.add_argument("--work", required=True, help="the working directory (trials, states, results)")
        p.add_argument("--runner", default=T.PYTHON_RUNNER,
                       help='how a trial runs: "python" (the Python API, python -m meds.model; the '
                            'default) or the path of a meds_main executable')
        p.add_argument("--pool", choices=("local", "queue"), default="local")
        p.add_argument("--workers", type=int, default=os.cpu_count())
        p.add_argument("--keep-netcdf", action="store_true")
    p = sub.add_parser("select-windows")
    p.add_argument("--site", required=True)
    p.set_defaults(func=cmd_select_windows)
    p = sub.add_parser("report")
    p.add_argument("--site", required=True)
    p.add_argument("--variant", default=None)
    p.add_argument("--work", required=True)
    p.set_defaults(func=cmd_report)
    for name, func in (("check", cmd_check), ("analyze", cmd_analyze)):
        p = sub.add_parser(name)
        common(p)
        p.set_defaults(func=func)
    p = sub.add_parser("fit")
    common(p)
    p.add_argument("--stages", default=None, help="comma-separated stages to run (default: [fit].stages)")
    p.add_argument("--resume", action="store_true", help="continue from <work>/stages.json")
    p.set_defaults(func=cmd_fit)
    p = sub.add_parser("smoke")
    common(p)
    p.add_argument("--days", type=int, default=3)
    p.add_argument("--key", default="stomatal_g1")
    p.add_argument("--g8", action="store_true", help="also check gate G8 (needs libmeds.so)")
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
    #----- a trial runs in its own directory, so an executable's path must not be relative
    if getattr(args, "runner", T.PYTHON_RUNNER) != T.PYTHON_RUNNER:
        args.runner = str(Path(args.runner).resolve())
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
