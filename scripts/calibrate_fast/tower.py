# SPDX-License-Identifier: Apache-2.0
"""The observations (MEDS_FAST_CALIBRATION_PLAN.md §5.2): the tower's half hours averaged to local
clock hours, the energy-balance closure correction, and the hours whose forcing was observed.

The tower's half hours are stamped at their start in local time; an hour is used only when both
of its half hours are there. The turbulent fluxes (H, LE, NEE, GPP, u*) are kept only where the
tower's flag says they were measured.

Closure (§5.2.1, FLUXNET2015's correction, Bowen ratio preserved): over a sliding window of
daytime hours, f = sum Rnet / sum (H + LE); daytime H and LE are multiplied by f. Night hours are
not corrected. The ground heat flux is taken as 0 when the file has none.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date

#: the quantities a declaration maps onto its file's column names
QUANTITIES = ("sw_in", "sw_up", "lw_up", "rnet", "le", "h", "nee", "gpp", "ustar")
FLAGGED = ("le", "h", "nee", "gpp", "ustar")


@dataclass
class TowerSpec:
    path: str
    columns: dict                     # quantity -> column name
    time_column: str = "date"
    flag_column: str = "FLAG"
    flag_good: int = 1
    utc_offset_h: float = 0.0
    forcing: str | None = None        # the forcing file, for the observed-hours mask
    forcing_qc: tuple = ("LWdown_qc", "Wind_qc")
    forcing_qc_val: tuple = ("Wind_qc",)   # validation windows: a tower year without longwave still validates
    forcing_grid: int = 1
    closure_days: int = 31            # the sliding window, centred
    daytime_sw: float = 10.0          # [W m-2] daytime = incoming shortwave above this
    extra: dict = field(default_factory=dict)


def load_tower(spec: TowerSpec) -> pd.DataFrame:
    """Hourly local-time observations: one column per quantity, NaN where not usable, plus the
    closure-corrected h_c and le_c and the closure factor."""
    raw = pd.read_csv(spec.path, parse_dates=[spec.time_column], index_col=spec.time_column)
    good = raw[spec.flag_column] == spec.flag_good
    half = pd.DataFrame(index=raw.index)
    for q in QUANTITIES:
        col = spec.columns.get(q)
        if col is None:
            half[q] = np.nan
            continue
        x = pd.to_numeric(raw[col], errors="coerce")
        half[q] = x.where(good) if q in FLAGGED else x
    hourly = half.resample("1h").mean().where(half.resample("1h").count() == 2)
    return add_closure(hourly, spec)


def add_closure(hourly: pd.DataFrame, spec: TowerSpec) -> pd.DataFrame:
    day = hourly["sw_in"] > spec.daytime_sw
    both = day & hourly[["rnet", "h", "le"]].notna().all(axis=1)
    rn = hourly["rnet"].where(both).resample("1D").sum(min_count=1)
    tf = (hourly["h"] + hourly["le"]).where(both).resample("1D").sum(min_count=1)
    rn = rn.rolling(spec.closure_days, center=True, min_periods=spec.closure_days // 3).sum()
    tf = tf.rolling(spec.closure_days, center=True, min_periods=spec.closure_days // 3).sum()
    f_daily = (rn / tf).where(tf > 0)
    f = f_daily.reindex(hourly.index.floor("D")).to_numpy()
    f = np.where(day.to_numpy(), f, 1.0)
    out = hourly.copy()
    out["closure_f"] = f
    out["h_c"] = out["h"] * f
    out["le_c"] = out["le"] * f
    return out


def forcing_observed(spec: TowerSpec, index: pd.DatetimeIndex, qc=None) -> pd.Series:
    """True for local hours whose forcing half hours are both observed (qc 0 in every listed
    variable, `qc` or the spec's calibration list); all True when no forcing file is declared."""
    if not spec.forcing:
        return pd.Series(True, index=index)
    qc = spec.forcing_qc if qc is None else qc
    with Dataset(spec.forcing) as ds:
        t = ds["time"]
        when = pd.to_datetime([d.strftime("%Y-%m-%d %H:%M:%S")
                               for d in num2date(t[:], t.units, only_use_cftime_datetimes=False)])
        ok = np.ones(len(when), dtype=bool)
        for v in qc:
            if v in ds.variables:
                ok &= np.asarray(ds[v][:, spec.forcing_grid - 1]) == 0
    local = when + pd.Timedelta(hours=spec.utc_offset_h)
    half = pd.Series(ok.astype(float), index=local)
    hour = half.resample("1h").agg(["sum", "count"])
    obs = (hour["sum"] == 2) & (hour["count"] == 2)
    return obs.reindex(index, fill_value=False)


def closure_summary(obs: pd.DataFrame) -> dict:
    day = obs["closure_f"] != 1.0
    f = obs.loc[day, "closure_f"].dropna()
    return {"median_f": float(f.median()) if len(f) else None,
            "p10_f": float(f.quantile(0.1)) if len(f) else None,
            "p90_f": float(f.quantile(0.9)) if len(f) else None}
