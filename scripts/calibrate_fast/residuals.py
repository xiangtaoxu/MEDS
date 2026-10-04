# SPDX-License-Identifier: Apache-2.0
"""Targets, their filters, and the residual vector (MEDS_FAST_CALIBRATION_PLAN.md §5.2-§5.4; the
revision plan MEDS_FAST_CALIBRATION_REVISION_PLAN.md §2, §3, §6).

For each window, the rows a target contributes -- which records, the observed value and its
observation error sigma -- depend only on the observations, so they are fixed once
(`WindowSpec`) and every trial's residual vector has the same rows in the same order:

    r = w_t * rho((y_model - y_obs) / sigma)

with w_t the target's effective-sample weight in that window (1 until the fit sets it) and rho the
identity, or the Huber transform when the fit asks for a robust loss.

The targets, each on the records the tower measured (the site TOML's rule, BCI: FLAG = 1) and the
forcing observed, after the first `skip_hours` of a window. A window's records are the tower's own
interval (30 min at BCI) on UTC, the clock of the model's output:

    albedo     sw_up / sw_in, records with sw_in > min_sw       sigma = abs
    lw_up      upwelling longwave                               sigma = abs
    rnet       net radiation                                    sigma = abs + rel |obs|
    le, h      closure-corrected (tower.py)                     sigma = abs + rel |obs|
    ef         daily daytime LE / (H + LE), uncorrected         sigma = abs
    gpp        daytime GPP                                      sigma = abs + rel |obs|
    nee_night  night NEE + growth respiration                   sigma = abs
    ustar      friction velocity                                sigma = abs + rel |obs|

Then each target's own filters, in this order (FILTERS; every one is a setting of the target,
documented with its default and reason in site_reference.toml):

    ustar_min            u* at or above (not on the u* target: it would bias it); a number, or a
                         table of years -- datarules.resolve_ustar turns "provider" and
                         "diagnostic" into one before the rows are built
    par_min              incident PAR at or above (needs the site TOML's variables.PAR)
    hours                local hours [from, to), on the site TOML's clock
    closure_range        days whose own closure ratio sum(H+LE)/sum(Rnet) is inside [lo, hi]
    min_solar_elevation  the sun at least this high [degrees]
    snow_free_days       no freezing air for this many days before (the tower has no snow sensor)

H, LE and u* are daytime targets: at night they are set mostly by the model's numerical floors and
by stable-air measurement problems. `night = true` keeps their night records too, at the provider's
u* threshold.

A frozen stand has no growth respiration (§3.3), so the model's night NEE is increased by the
stand's monthly growth-respiration climatology from a run with the slow loop on (§5.3).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

import datarules as DR

TARGETS = ("albedo", "lw_up", "rnet", "le", "h", "ef", "gpp", "nee_night", "ustar")
FILTERS = ("ustar_min", "par_min", "hours", "closure_range", "min_solar_elevation", "snow_free_days")
KGC_YR_TO_UMOL_S = 1000.0 / 12.011 * 1.0e6 / (365.25 * 86400.0)   # kgC m-2 yr-1 -> umol m-2 s-1


@dataclass
class TargetRows:
    name: str
    rows: np.ndarray            # positions in the window's index (for ef, each day's first record)
    obs: np.ndarray
    sigma: np.ndarray
    groups: list = field(default_factory=list)   # ef: per row, the record positions of that day
    months: np.ndarray | None = None             # nee_night: calendar month of each row
    weight: float = 1.0                          # the effective-sample weight sqrt(n_eff / n)
    counts: list = field(default_factory=list)   # (step, rows left) through the base mask and filters


@dataclass
class WindowSpec:
    name: str
    index: pd.DatetimeIndex      # the window's records (UTC starts), as the trial writes them
    targets: list                # TargetRows, in TARGETS order (empty targets dropped)
    loss: tuple = ("l2", 2.0)    # ("l2" | "huber", Huber's c)

    @property
    def n(self) -> int:
        return sum(len(t.obs) for t in self.targets)

    def slices(self) -> dict:
        out, i = {}, 0
        for t in self.targets:
            out[t.name] = slice(i, i + len(t.obs))
            i += len(t.obs)
        return out

    def subset(self, names) -> "WindowSpec":
        """The same window with only these targets (the kernel stages fit one target each)."""
        return replace(self, targets=[t for t in self.targets if t.name in names])


def window_index(start_utc, days: int, step: float) -> pd.DatetimeIndex:
    """A window's records: the UTC start of each of the tower's intervals (`step` seconds)."""
    return pd.date_range(pd.Timestamp(start_utc), periods=int(round(days * 86400.0 / step)),
                         freq=f"{int(round(step))}s")


def local_hours(index: pd.DatetimeIndex, utc_offset_h: float) -> np.ndarray:
    """The local clock hour of each record's start."""
    return (index + pd.Timedelta(hours=utc_offset_h)).hour.to_numpy()


def _sigma(cfg: dict, obs: np.ndarray) -> np.ndarray:
    return cfg.get("sigma_abs", cfg.get("sigma", 0.0)) + cfg.get("sigma_rel", 0.0) * np.abs(obs)


def filter_masks(c: dict, o: pd.DataFrame, elev: np.ndarray | None, utc_offset_h: float = 0.0) -> list:
    """The target's filters as (label, mask) in FILTERS order, for the records of `o`."""
    out = []
    if c.get("ustar_min") is not None:
        thr = c["ustar_min"]
        label = f"u* >= {thr}" if not isinstance(thr, dict) else "u* >= the provider's yearly threshold"
        out.append((label, o["ustar"].to_numpy() >= DR.threshold_at(o.index, thr)))
    if c.get("par_min") is not None and float(c["par_min"]) > 0:
        par = o["par"].to_numpy() if "par" in o else np.full(len(o), np.nan)
        out.append((f"PAR >= {c['par_min']}", par >= float(c["par_min"])))
    if c.get("hours") is not None:
        a, b = c["hours"]
        hr = local_hours(o.index, utc_offset_h)
        out.append((f"local hours [{a}, {b})", (hr >= a) & (hr < b)))
    if c.get("closure_range") is not None:
        lo, hi = c["closure_range"]
        cd = o["closure_day"].to_numpy() if "closure_day" in o else np.full(len(o), np.nan)
        out.append((f"day closure in [{lo}, {hi}]", (cd >= lo) & (cd <= hi)))
    if c.get("min_solar_elevation") is not None and float(c["min_solar_elevation"]) > 0:
        if elev is None:
            raise ValueError("min_solar_elevation needs the site's latitude and longitude")
        out.append((f"sun >= {c['min_solar_elevation']} deg", elev >= float(c["min_solar_elevation"])))
    if c.get("snow_free_days") is not None and int(c["snow_free_days"]) > 0:
        if "days_since_frost" not in o:
            raise ValueError("snow_free_days needs the tower's air temperature (variables.Tair)")
        out.append((f"no frost in {c['snow_free_days']} d", o["days_since_frost"].to_numpy() > int(c["snow_free_days"])))
    return out


def base_masks(tname: str, c: dict, o: pd.DataFrame, daytime_sw: float):
    """The rows a target could use before its filters, and their observed values: (steps, values),
    steps a list of (label, mask) of the target's own definition (daytime, night, sw_in > min_sw)."""
    sw = o["sw_in"].to_numpy()
    day = sw > daytime_sw
    every = np.ones(len(o), bool)
    if tname == "albedo":
        with np.errstate(divide="ignore", invalid="ignore"):
            v = o["sw_up"].to_numpy() / sw
        return [(f"sw_in > {c.get('min_sw', 200.0)}", sw > c.get("min_sw", 200.0))], v
    if tname in ("lw_up", "rnet"):
        return [], o[tname].to_numpy()
    if tname in ("le", "h", "ustar"):
        v = o[f"{tname}_c" if tname != "ustar" and c.get("closure", True) else tname].to_numpy()
        if not c.get("night"):
            return [("daytime", day)], v
        thr = DR.threshold_at(o.index, c["night_ustar"])
        return [(f"daytime, or night at u* >= the provider's", day | (o["ustar"].to_numpy() >= thr))], v
    if tname == "gpp":
        return [("daytime", day)], o["gpp"].to_numpy()
    if tname == "nee_night":
        return [(f"sw_in < {c.get('night_sw', 5.0)}", sw < c.get("night_sw", 5.0))], o["nee"].to_numpy()
    raise KeyError(tname)


def build_spec(name, index, obs: pd.DataFrame, forcing_ok: pd.Series, targets_cfg: dict,
               skip_hours: int = 0, daytime_sw: float = 10.0, elev: np.ndarray | None = None,
               utc_offset_h: float = 0.0) -> WindowSpec:
    """The window's rows: each target's base rows (measured, observed forcing, after skip_hours,
    and the target's own definition), then its filters, with the rows left after every step kept
    in TargetRows.counts for the filter report. `elev` is the sun's elevation at the window's
    records (for min_solar_elevation); `utc_offset_h` the site's clock (local hours and days)."""
    o = obs.reindex(index)
    ok = np.array(forcing_ok.reindex(index, fill_value=False).to_numpy(), dtype=bool)
    step = (index[1] - index[0]).total_seconds() if len(index) > 1 else 3600.0
    ok[:int(round(skip_hours * 3600.0 / step))] = False
    sw = o["sw_in"].to_numpy()
    day = sw > daytime_sw
    rows = []
    for tname in TARGETS:
        if tname not in targets_cfg or not targets_cfg[tname].get("on", True):
            continue
        c = targets_cfg[tname]
        if tname == "ustar" and c.get("ustar_min") is not None:
            raise ValueError("[targets.ustar].ustar_min would bias the u* target; filter the other targets")
        if tname == "ef":
            h, le = o["h"].to_numpy(), o["le"].to_numpy()
            use = day & ok & np.isfinite(h) & np.isfinite(le)
            for _, m in filter_masks(c, o, elev, utc_offset_h):
                use &= m
            groups, vals = [], []
            days = (index + pd.Timedelta(hours=utc_offset_h)).normalize()
            for d in np.unique(days[use]):
                g = np.flatnonzero(use & (days == d))
                if len(g) * step / 3600.0 >= c.get("min_hours", 6):
                    sh, sl = h[g].sum(), le[g].sum()
                    if sh + sl > 0:
                        groups.append(g)
                        vals.append(sl / (sh + sl))
            if vals:
                v = np.array(vals)
                rows.append(TargetRows(tname, np.array([g[0] for g in groups]), v, _sigma(c, v),
                                       groups=groups, counts=[("days", len(v))]))
            continue
        steps, values = base_masks(tname, c, o, daytime_sw)
        mask = np.isfinite(values)
        counts = [("measured", int(mask.sum()))]
        mask &= ok
        counts.append(("forcing observed", int(mask.sum())))
        for label, m in steps + filter_masks(c, o, elev, utc_offset_h):
            mask &= m
            counts.append((label, int(mask.sum())))
        pos = np.flatnonzero(mask)
        if len(pos):
            v = values[pos]
            t = TargetRows(tname, pos, v, _sigma(c, v), counts=counts)
            if tname == "nee_night":
                t.months = index.month.to_numpy()[pos]
            rows.append(t)
    return WindowSpec(name, index, rows)


def model_values(t: TargetRows, df: pd.DataFrame, growth_resp=None) -> np.ndarray:
    """The model side of one target, on its rows."""
    if t.name == "albedo":
        return df["sw_up_fast"].to_numpy()[t.rows] / df["sw_in_fast"].to_numpy()[t.rows]
    if t.name == "lw_up":
        return df["lw_up_fast"].to_numpy()[t.rows]
    if t.name == "rnet":
        return df["rnet_fast"].to_numpy()[t.rows]
    if t.name == "le":
        return df["le_flux_fast"].to_numpy()[t.rows]
    if t.name == "h":
        return df["h_flux_fast"].to_numpy()[t.rows]
    if t.name == "gpp":
        return df["gpp_rate_fast"].to_numpy()[t.rows]
    if t.name == "ustar":
        return df["ustar_fast"].to_numpy()[t.rows]
    if t.name == "nee_night":
        gr = 0.0 if growth_resp is None else np.array([growth_resp[m] for m in t.months])
        return df["nee_fast"].to_numpy()[t.rows] + gr
    if t.name == "ef":
        h, le = df["h_flux_fast"].to_numpy(), df["le_flux_fast"].to_numpy()
        return np.array([le[g].sum() / (h[g].sum() + le[g].sum()) for g in t.groups])
    raise KeyError(t.name)


def huber(z: np.ndarray, c: float) -> np.ndarray:
    """The residual whose square is Huber's loss times 2: z inside +-c, sign(z) sqrt(2c|z| - c^2)
    beyond, so least squares on it is the robust fit (each row's pull is capped at c)."""
    a = np.abs(z)
    return np.where(a <= c, z, np.sign(z) * np.sqrt(np.maximum(2.0 * c * a - c * c, 0.0)))


def residual(spec: WindowSpec, df: pd.DataFrame, growth_resp=None, weighted: bool = True) -> np.ndarray:
    """w_t rho((model - obs) / sigma) for every row of the window, in the spec's fixed order."""
    if not spec.targets:
        return np.zeros(0)
    m = df.reindex(spec.index)
    need = [c for c in m.columns]
    if m[need].isna().any().any():
        raise ValueError(f"{spec.name}: the trial's output does not cover the window")
    parts = []
    for t in spec.targets:
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (model_values(t, m, growth_resp) - t.obs) / t.sigma
        if not np.all(np.isfinite(z)):
            raise ValueError(f"{spec.name}: non-finite {t.name} residuals at {int((~np.isfinite(z)).sum())} rows")
        if spec.loss[0] == "huber":
            z = huber(z, float(spec.loss[1]))
        parts.append(z * (t.weight if weighted else 1.0))
    return np.concatenate(parts)


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


def set_weights(specs: list, r: np.ndarray) -> dict:
    """Put each target's effective-sample weight sqrt(n_eff / n) into its rows, from residuals r
    taken WITHOUT weights; returns {window/target: n_eff / n}."""
    w = ess_weights(specs, r)
    out, i = {}, 0
    for spec in specs:
        for t in spec.targets:
            n = len(t.obs)
            t.weight = float(np.sqrt(w[i])) if n else 1.0
            out[f"{spec.name}/{t.name}"] = float(w[i]) if n else 1.0
            i += n
    return out


def unweighted(specs: list, r: np.ndarray) -> np.ndarray:
    """The rows of r with the effective-sample weights taken back out."""
    out, i = r.copy(), 0
    for spec in specs:
        for t in spec.targets:
            n = len(t.obs)
            if t.weight > 0:
                out[i:i + n] = r[i:i + n] / t.weight
            i += n
    return out


def target_scores(specs: list, r: np.ndarray) -> dict:
    """Normalized RMSE per target, pooled over the windows: sqrt(mean(z^2)), z the unweighted
    normalized residual; with each target's share of the fit's objective (weighted)."""
    z = unweighted(specs, r)
    acc, i = {}, 0
    for spec in specs:
        for t in spec.targets:
            n = len(t.obs)
            s = acc.setdefault(t.name, [0.0, 0, 0.0])
            s[0] += float((z[i:i + n] ** 2).sum())
            s[1] += n
            s[2] += float((r[i:i + n] ** 2).sum())
            i += n
    total = sum(v[2] for v in acc.values()) or 1.0
    return {k: {"nrmse": float(np.sqrt(v[0] / v[1])), "n": v[1], "share": float(v[2] / total)}
            for k, v in acc.items()}


def sigma_scales(specs: list, r: np.ndarray) -> dict:
    """Each target's chi^2 per row, mean(z^2) of its unweighted normalized residual: 1 when its sigma
    describes the misfit; 2 or more marks a structural error in that target."""
    return {k: v["nrmse"] ** 2 for k, v in target_scores(specs, r).items()}


def load_growth_resp(path) -> dict | None:
    """{month: umol m-2 s-1} from a CSV with columns month, growth_resp_umol."""
    if path is None:
        return None
    df = pd.read_csv(path)
    return {int(m): float(v) for m, v in zip(df["month"], df["growth_resp_umol"])}


# ------------------------------------------------------------------------------------------------
# the filter report (revision plan §2): the whole record, every target, every filter in turn
# ------------------------------------------------------------------------------------------------
def record_report(obs: pd.DataFrame, forcing_ok: pd.Series, targets_cfg: dict, daytime_sw: float,
                  elev: np.ndarray | None, utc_offset_h: float = 0.0) -> dict:
    """For each target over the whole tower record: the records left after each step, and the
    diurnal mean of the observation (by local hour) before and after the filters."""
    out = {}
    ok = forcing_ok.reindex(obs.index, fill_value=False).to_numpy()
    hr = local_hours(obs.index, utc_offset_h)
    for tname in TARGETS:
        if tname == "ef" or tname not in targets_cfg or not targets_cfg[tname].get("on", True):
            continue
        c = targets_cfg[tname]
        steps, values = base_masks(tname, c, obs, daytime_sw)
        mask = np.isfinite(values)
        counts = [("measured", int(mask.sum()))]
        mask &= ok
        counts.append(("forcing observed", int(mask.sum())))
        for label, m in steps:
            mask &= m
            counts.append((label, int(mask.sum())))
        before = mask.copy()
        for label, m in filter_masks(c, obs, elev, utc_offset_h):
            mask &= m
            counts.append((label, int(mask.sum())))
        diurnal = {}
        for h in range(24):
            b, a = before & (hr == h), mask & (hr == h)
            if b.sum() >= 10:
                diurnal[h] = [float(np.nanmean(values[b])), float(np.nanmean(values[a])) if a.sum() else None,
                              int(b.sum()), int(a.sum())]
        out[tname] = {"counts": counts, "diurnal": diurnal}
    return out
