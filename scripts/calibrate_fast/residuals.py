# SPDX-License-Identifier: Apache-2.0
"""Targets, their filters, and the residual vector (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2.5,
§3).

For each window, the rows a target contributes -- which records, the observed value and its
observation error sigma -- depend only on the observations, so they are fixed once
(`WindowSpec`) and every trial's residual vector has the same rows in the same order:

    r = w_t * rho((y_model - y_obs) / sigma)

with w_t the target's effective-sample weight in that window (1 until the fit sets it) and rho
Huber's transform (or the identity for least squares).

The targets, each on the records the tower measured (the site TOML's rule, BCI: FLAG = 1) and the
forcing observed, after the first `skip_hours` of a window. A window's records are the tower's own
interval (30 min at BCI) on UTC, the clock of the model's output. Each target has an observation
model (obsmodels.py):

    albedo     sw_up / sw_in, records with sw_in > min_sw          identity
    lw_up      upwelling longwave                                  identity
    le, h      daytime                                             closure: the tower's value
                                                                   corrected for the closure gap
    gpp        daytime                                             respiration: the tower's GPP plus
                                                                   (1/kappa - 1) times its respiration
    ustar      daytime friction velocity                           identity, a loose sigma

sigma = sigma_abs + sigma_rel |x|, with x the smoothed observation of the target's error flux (NEE's
for GPP) where obsmodels.set_sigmas set one, else the observation itself.

Then each target's own filters, in this order (FILTERS; every one is a setting of the target,
documented with its default and reason in site_reference.toml):

    ustar_min            u* at or above (not on the u* target: it would bias it); a number, or a
                         table of years -- datarules.resolve_ustar turns "provider" and
                         "diagnostic" into one before the rows are built
    par_min              incident PAR at or above (needs the site TOML's variables.PAR)
    hours                local hours [from, to), on the site TOML's clock
    min_solar_elevation  the sun at least this high [degrees]
    snow_free_days       no freezing air for this many days before (the tower has no snow sensor)

H, LE and u* are daytime targets: at night they are set mostly by the model's numerical floors and
by stable-air measurement problems. `night = true` keeps their night records too, at the provider's
u* threshold.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

import datarules as DR

TARGETS = ("albedo", "lw_up", "le", "h", "gpp", "ustar")
FILTERS = ("ustar_min", "par_min", "hours", "min_solar_elevation", "snow_free_days")
OBS_MODELS = {"albedo": ("identity",), "lw_up": ("identity",), "ustar": ("identity",),
              "le": ("closure", "identity"), "h": ("closure", "identity"), "gpp": ("respiration", "identity")}
MODEL_COLUMN = {"lw_up": "lw_up_fast", "le": "le_flux_fast", "h": "h_flux_fast", "gpp": "gpp_rate_fast",
                "ustar": "ustar_fast"}


@dataclass
class TargetRows:
    name: str
    rows: np.ndarray            # positions in the window's index
    obs: np.ndarray
    sigma: np.ndarray
    reco: np.ndarray | None = None               # gpp's respiration model: the tower's respiration
    weight: float = 1.0                          # the effective-sample weight sqrt(n_eff / n)
    counts: list = field(default_factory=list)   # (step, rows left) through the base mask and filters


@dataclass
class WindowSpec:
    name: str
    index: pd.DatetimeIndex      # the window's records (UTC starts), as the trial writes them
    targets: list                # TargetRows, in TARGETS order (empty targets dropped)
    loss: tuple = ("huber", 2.0)  # ("l2" | "huber", Huber's c)

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
        """The same window with only these targets (the seasonal runs are scored on LE alone)."""
        return replace(self, targets=[t for t in self.targets if t.name in names])


def window_index(start_utc, days: int, step: float) -> pd.DatetimeIndex:
    """A window's records: the UTC start of each of the tower's intervals (`step` seconds)."""
    return pd.date_range(pd.Timestamp(start_utc), periods=int(round(days * 86400.0 / step)),
                         freq=f"{int(round(step))}s")


def local_hours(index: pd.DatetimeIndex, utc_offset_h: float) -> np.ndarray:
    """The local clock hour of each record's start."""
    return (index + pd.Timedelta(hours=utc_offset_h)).hour.to_numpy()


def _sigma(cfg: dict, at: np.ndarray) -> np.ndarray:
    return cfg.get("sigma_abs", cfg.get("sigma", 0.0)) + cfg.get("sigma_rel", 0.0) * np.abs(at)


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
    if c.get("min_solar_elevation") is not None and float(c["min_solar_elevation"]) > 0:
        if elev is None:
            raise ValueError("min_solar_elevation needs the site's latitude and longitude")
        out.append((f"sun >= {c['min_solar_elevation']} deg", elev >= float(c["min_solar_elevation"])))
    if c.get("snow_free_days") is not None and int(c["snow_free_days"]) > 0:
        if "days_since_frost" not in o:
            raise ValueError("snow_free_days needs the tower's air temperature (variables.Tair)")
        out.append((f"no frost in {c['snow_free_days']} d", o["days_since_frost"].to_numpy() > int(c["snow_free_days"])))
    return out


def obs_model(tname: str, c: dict) -> str:
    m = c.get("obs_model", OBS_MODELS[tname][0])
    if m not in OBS_MODELS[tname]:
        raise ValueError(f"[targets.{tname}].obs_model must be one of {OBS_MODELS[tname]}, not {m!r}")
    return m


def base_masks(tname: str, c: dict, o: pd.DataFrame, daytime_sw: float):
    """The rows a target could use before its filters, and their observed values: (steps, values),
    steps a list of (label, mask) of the target's own definition (daytime, sw_in > min_sw)."""
    sw = o["sw_in"].to_numpy()
    day = sw > daytime_sw
    model = obs_model(tname, c)
    if tname == "albedo":
        with np.errstate(divide="ignore", invalid="ignore"):
            v = o["sw_up"].to_numpy() / sw
        return [(f"sw_in > {c.get('min_sw', 200.0)}", sw > c.get("min_sw", 200.0))], v
    if tname == "lw_up":
        return [], o["lw_up"].to_numpy()
    if tname == "gpp":
        v = o["gpp"].to_numpy()
        if model == "respiration":
            v = np.where(np.isfinite(o["reco"].to_numpy()), v, np.nan)
        return [("daytime", day)], v
    v = o[f"{tname}_c" if model == "closure" else tname].to_numpy()
    if not c.get("night"):
        return [("daytime", day)], v
    thr = DR.threshold_at(o.index, c["night_ustar"])
    return [("daytime, or night at u* >= the provider's", day | (o["ustar"].to_numpy() >= thr))], v


def build_spec(name, index, obs: pd.DataFrame, forcing_ok: pd.Series, targets_cfg: dict,
               skip_hours: int = 0, daytime_sw: float = 10.0, elev: np.ndarray | None = None,
               utc_offset_h: float = 0.0) -> WindowSpec:
    """The window's rows: each target's base rows (measured, observed forcing, after skip_hours,
    and the target's own definition), then its filters, with the rows left after every step kept
    in TargetRows.counts for the filter report. `elev` is the sun's elevation at the window's
    records (for min_solar_elevation); `utc_offset_h` the site's clock (local hours)."""
    o = obs.reindex(index)
    ok = np.array(forcing_ok.reindex(index, fill_value=False).to_numpy(), dtype=bool)
    step = (index[1] - index[0]).total_seconds() if len(index) > 1 else 3600.0
    ok[:int(round(skip_hours * 3600.0 / step))] = False
    rows = []
    for tname in TARGETS:
        if tname not in targets_cfg or not targets_cfg[tname].get("on", True):
            continue
        c = targets_cfg[tname]
        if tname == "ustar" and c.get("ustar_min") is not None:
            raise ValueError("[targets.ustar].ustar_min would bias the u* target; filter the other targets")
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
            at = o[c["sigma_at"]].to_numpy()[pos] if c.get("sigma_at") else v
            at = np.where(np.isfinite(at), at, v)
            t = TargetRows(tname, pos, v, _sigma(c, at), counts=counts)
            if tname == "gpp" and obs_model(tname, c) == "respiration":
                t.reco = o["reco"].to_numpy()[pos]
            rows.append(t)
    return WindowSpec(name, index, rows)


def model_values(t: TargetRows, df: pd.DataFrame) -> np.ndarray:
    """The model side of one target, on its rows."""
    if t.name == "albedo":
        return df["sw_up_fast"].to_numpy()[t.rows] / df["sw_in_fast"].to_numpy()[t.rows]
    return df[MODEL_COLUMN[t.name]].to_numpy()[t.rows]


def observed(t: TargetRows, obs_keys: dict | None) -> np.ndarray:
    """The observation side of one target: GPP's respiration model adds (1/kappa - 1) R_tower."""
    if t.reco is None:
        return t.obs
    kappa = float((obs_keys or {}).get("kappa", 1.0))
    return t.obs + (1.0 / kappa - 1.0) * t.reco


def huber(z: np.ndarray, c: float) -> np.ndarray:
    """The residual whose square is Huber's loss times 2: z inside +-c, sign(z) sqrt(2c|z| - c^2)
    beyond, so least squares on it is the robust fit (each row's pull is capped at c)."""
    a = np.abs(z)
    return np.where(a <= c, z, np.sign(z) * np.sqrt(np.maximum(2.0 * c * a - c * c, 0.0)))


def unhuber(h: np.ndarray, c: float) -> np.ndarray:
    """The normalized residual z back from its Huber transform (huber's exact inverse)."""
    a = np.abs(h)
    return np.where(a <= c, h, np.sign(h) * (a * a + c * c) / (2.0 * c))


def residual(spec: WindowSpec, df: pd.DataFrame, obs_keys: dict | None = None, weighted: bool = True) -> np.ndarray:
    """w_t rho((model - obs) / sigma) for every row of the window, in the spec's fixed order."""
    if not spec.targets:
        return np.zeros(0)
    m = df.reindex(spec.index)
    if m.isna().any().any():
        raise ValueError(f"{spec.name}: the trial's output does not cover the window")
    parts = []
    for t in spec.targets:
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (model_values(t, m) - observed(t, obs_keys)) / t.sigma
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
    residuals: n_eff = n (1 - rho) / (1 + rho), rho clipped to [0, 0.99]."""
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


def raw(specs: list, r: np.ndarray) -> np.ndarray:
    """The rows of r as plain normalized residuals (model - obs) / sigma: the effective-sample
    weights taken out and the Huber transform undone."""
    out, i = r.copy(), 0
    for spec in specs:
        for t in spec.targets:
            n = len(t.obs)
            if t.weight > 0:
                out[i:i + n] = r[i:i + n] / t.weight
            if spec.loss[0] == "huber":
                out[i:i + n] = unhuber(out[i:i + n], float(spec.loss[1]))
            i += n
    return out


def target_scores(specs: list, r: np.ndarray) -> dict:
    """Per target, pooled over the windows: the normalized RMSE sqrt(mean(z^2)) of the plain
    residuals (chi^2 per row is its square), and the target's share of the fit's objective."""
    z = raw(specs, r)
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


def chi2_per_row(specs: list, r: np.ndarray) -> dict:
    """Each target's chi^2 per row, mean(z^2) of its plain residuals: 1 when its sigma describes the
    misfit; 2 or more marks a structural error in that target."""
    return {k: v["nrmse"] ** 2 for k, v in target_scores(specs, r).items()}


def scale_sigma(specs: list, r: np.ndarray, cap: float) -> dict:
    """The model-error step (best-practice plan §3.4): multiply each target's sigma by
    max(1, sqrt(chi^2 per row)), at most `cap`, in every window, so a target the model cannot fit
    stops bending the keys. Returns {target: scale}."""
    scales = {k: float(min(max(1.0, np.sqrt(v)), cap)) for k, v in chi2_per_row(specs, r).items()}
    for spec in specs:
        for t in spec.targets:
            t.sigma = t.sigma * scales.get(t.name, 1.0)
    return scales


# ------------------------------------------------------------------------------------------------
# the filter report: the whole record, every target, every filter in turn
# ------------------------------------------------------------------------------------------------
def record_report(obs: pd.DataFrame, forcing_ok: pd.Series, targets_cfg: dict, daytime_sw: float,
                  elev: np.ndarray | None, utc_offset_h: float = 0.0) -> dict:
    """For each target over the whole tower record: the records left after each step, and the
    diurnal mean of the observation (by local hour) before and after the filters."""
    out = {}
    ok = forcing_ok.reindex(obs.index, fill_value=False).to_numpy()
    hr = local_hours(obs.index, utc_offset_h)
    for tname in TARGETS:
        if tname not in targets_cfg or not targets_cfg[tname].get("on", True):
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
