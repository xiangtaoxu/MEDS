# SPDX-License-Identifier: Apache-2.0
"""Targets and the residual vector (MEDS_FAST_CALIBRATION_PLAN.md §5.2-§5.4).

For each window, the rows a target contributes -- which hours, the observed value and its
observation error sigma -- depend only on the observations, so they are fixed once
(`WindowSpec`) and every trial's residual vector has the same rows in the same order:

    r = (y_model - y_obs) / sigma

The targets, each on the hours the tower and the forcing observed (the forcing mask is applied to
all of them) after the first `skip_hours` of a window:

    albedo     sw_up / sw_in, hours with sw_in > min_sw         sigma = abs
    lw_up      upwelling longwave                               sigma = abs
    rnet       net radiation                                    sigma = abs + rel |obs|
    le, h      closure-corrected (tower.py)                     sigma = abs + rel |obs|
    ef         daily daytime LE / (H + LE), uncorrected         sigma = abs
    gpp        daytime GPP                                      sigma = abs + rel |obs|
    nee_night  night NEE (u* >= ustar_min) + growth respiration sigma = abs
    ustar      friction velocity                                sigma = abs + rel |obs|

A frozen stand has no growth respiration (§3.3), so the model's night NEE is increased by the
stand's monthly growth-respiration climatology from a run with the slow loop on (§5.3).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

TARGETS = ("albedo", "lw_up", "rnet", "le", "h", "ef", "gpp", "nee_night", "ustar")
KGC_YR_TO_UMOL_S = 1000.0 / 12.011 * 1.0e6 / (365.25 * 86400.0)   # kgC m-2 yr-1 -> umol m-2 s-1


@dataclass
class TargetRows:
    name: str
    hours: np.ndarray           # positions in the window's hourly index (or, for ef, day groups)
    obs: np.ndarray
    sigma: np.ndarray
    groups: list = field(default_factory=list)   # ef: per row, the hour positions of that day
    months: np.ndarray | None = None             # nee_night: calendar month of each row


@dataclass
class WindowSpec:
    name: str
    index: pd.DatetimeIndex      # the window's local hours, as the trial writes them
    targets: list                # TargetRows, in TARGETS order (empty targets dropped)

    @property
    def n(self) -> int:
        return sum(len(t.obs) for t in self.targets)

    def slices(self) -> dict:
        out, i = {}, 0
        for t in self.targets:
            out[t.name] = slice(i, i + len(t.obs))
            i += len(t.obs)
        return out


def window_index(start_utc, days: int, utc_offset_h: float) -> pd.DatetimeIndex:
    t0 = pd.Timestamp(start_utc) + pd.Timedelta(hours=utc_offset_h)
    return pd.date_range(t0, periods=24 * days, freq="1h")


def _sigma(cfg: dict, obs: np.ndarray) -> np.ndarray:
    return cfg.get("sigma_abs", cfg.get("sigma", 0.0)) + cfg.get("sigma_rel", 0.0) * np.abs(obs)


def build_spec(name, index, obs: pd.DataFrame, forcing_ok: pd.Series, targets_cfg: dict,
               skip_hours: int = 0, daytime_sw: float = 10.0) -> WindowSpec:
    o = obs.reindex(index)
    ok = np.array(forcing_ok.reindex(index, fill_value=False).to_numpy(), dtype=bool)
    ok[:skip_hours] = False
    sw = o["sw_in"].to_numpy()
    day = sw > daytime_sw
    rows = []

    def add(tname, mask, values):
        mask = mask & ok & np.isfinite(values)
        pos = np.flatnonzero(mask)
        if len(pos):
            v = values[pos]
            rows.append(TargetRows(tname, pos, v, _sigma(targets_cfg[tname], v)))

    for tname in TARGETS:
        if tname not in targets_cfg or not targets_cfg[tname].get("on", True):
            continue
        c = targets_cfg[tname]
        if tname == "albedo":
            with np.errstate(divide="ignore", invalid="ignore"):
                alb = o["sw_up"].to_numpy() / sw
            add(tname, sw > c.get("min_sw", 200.0), alb)
        elif tname == "lw_up":
            add(tname, np.ones(len(index), bool), o["lw_up"].to_numpy())
        elif tname == "rnet":
            add(tname, np.ones(len(index), bool), o["rnet"].to_numpy())
        elif tname in ("le", "h"):
            col = f"{tname}_c" if c.get("closure", True) else tname
            add(tname, np.ones(len(index), bool), o[col].to_numpy())
        elif tname == "gpp":
            add(tname, day, o["gpp"].to_numpy())
        elif tname == "ustar":
            add(tname, np.ones(len(index), bool), o["ustar"].to_numpy())
        elif tname == "nee_night":
            night = sw < c.get("night_sw", 5.0)
            us = o["ustar"].to_numpy()
            add(tname, night & (us >= c.get("ustar_min", 0.2)), o["nee"].to_numpy())
            if rows and rows[-1].name == tname:
                rows[-1].months = index.month.to_numpy()[rows[-1].hours]
        elif tname == "ef":
            h, le = o["h"].to_numpy(), o["le"].to_numpy()
            use = day & ok & np.isfinite(h) & np.isfinite(le)
            groups, vals = [], []
            for d in np.unique(index.normalize()[use]):
                g = np.flatnonzero(use & (index.normalize() == d))
                if len(g) >= c.get("min_hours", 6):
                    sh, sl = h[g].sum(), le[g].sum()
                    if sh + sl > 0:
                        groups.append(g)
                        vals.append(sl / (sh + sl))
            if vals:
                v = np.array(vals)
                rows.append(TargetRows(tname, np.array([g[0] for g in groups]), v,
                                       _sigma(c, v), groups=groups))
    return WindowSpec(name, index, rows)


def model_values(t: TargetRows, df: pd.DataFrame, growth_resp=None) -> np.ndarray:
    """The model side of one target, on its rows."""
    if t.name == "albedo":
        return df["sw_up_fast"].to_numpy()[t.hours] / df["sw_in_fast"].to_numpy()[t.hours]
    if t.name == "lw_up":
        return df["lw_up_fast"].to_numpy()[t.hours]
    if t.name == "rnet":
        return df["rnet_fast"].to_numpy()[t.hours]
    if t.name == "le":
        return df["le_flux_fast"].to_numpy()[t.hours]
    if t.name == "h":
        return df["h_flux_fast"].to_numpy()[t.hours]
    if t.name == "gpp":
        return df["gpp_rate_fast"].to_numpy()[t.hours]
    if t.name == "ustar":
        return df["ustar_fast"].to_numpy()[t.hours]
    if t.name == "nee_night":
        gr = 0.0 if growth_resp is None else np.array([growth_resp[m] for m in t.months])
        return df["nee_fast"].to_numpy()[t.hours] + gr
    if t.name == "ef":
        h, le = df["h_flux_fast"].to_numpy(), df["le_flux_fast"].to_numpy()
        return np.array([le[g].sum() / (h[g].sum() + le[g].sum()) for g in t.groups])
    raise KeyError(t.name)


def residual(spec: WindowSpec, df: pd.DataFrame, growth_resp=None) -> np.ndarray:
    """(model - obs) / sigma for every row of the window, in the spec's fixed order."""
    if not spec.targets:
        return np.zeros(0)
    m = df.reindex(spec.index)
    if m.isna().any().any():
        raise ValueError(f"{spec.name}: the trial's hourly output does not cover the window")
    return np.concatenate([(model_values(t, m, growth_resp) - t.obs) / t.sigma for t in spec.targets])


def lag1(x: np.ndarray) -> float:
    if len(x) < 3:
        return 0.0
    a, b = x[:-1] - x.mean(), x[1:] - x.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def ess_weights(specs: list, r: np.ndarray) -> np.ndarray:
    """Per row, n_eff / n of its target in its window, from the lag-1 autocorrelation of the
    residuals (§5.4): n_eff = n (1 - rho) / (1 + rho), rho clipped to [0, 0.99]."""
    w = np.ones_like(r)
    i = 0
    for spec in specs:
        for t in spec.targets:
            n = len(t.obs)
            rho = min(max(lag1(r[i:i + n]), 0.0), 0.99)
            w[i:i + n] = (1.0 - rho) / (1.0 + rho)
            i += n
    return w


def target_scores(specs: list, r: np.ndarray) -> dict:
    """Normalized RMSE per target, pooled over the windows: sqrt(mean(r^2))."""
    acc = {}
    i = 0
    for spec in specs:
        for t in spec.targets:
            n = len(t.obs)
            s = acc.setdefault(t.name, [0.0, 0])
            s[0] += float((r[i:i + n] ** 2).sum())
            s[1] += n
            i += n
    return {k: {"nrmse": float(np.sqrt(v[0] / v[1])), "n": v[1]} for k, v in acc.items()}


def load_growth_resp(path) -> dict | None:
    """{month: umol m-2 s-1} from a CSV with columns month, growth_resp_umol."""
    if path is None:
        return None
    df = pd.read_csv(path)
    return {int(m): float(v) for m, v in zip(df["month"], df["growth_resp_umol"])}
