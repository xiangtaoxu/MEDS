# SPDX-License-Identifier: Apache-2.0
"""The observations (MEDS_FAST_CALIBRATION_PLAN.md §5.2): the tower's half hours averaged to local
clock hours, the energy-balance closure correction, and the hours whose forcing was observed.

The tower's half hours are stamped at their start in local time; an hour is used only when both
of its half hours are there. The turbulent fluxes (H, LE, NEE, GPP, u*) are kept only where the
tower's flag says they were measured.

Closure (§5.2.1, FLUXNET2015's correction, Bowen ratio preserved): over a sliding window of
daytime hours, f = sum Rnet / sum (H + LE); daytime H and LE are multiplied by f. Night hours are
not corrected. The ground heat flux is taken as 0 when the file has none. `closure = "none"` keeps
H and LE as measured. Each day's own closure ratio, sum (H + LE) / sum Rnet over its measured
daytime hours, is kept for the targets' optional `closure_range` filter.

Incident PAR (quantity `par`, optional) is read as is; the GPP target's `par_min` filter uses it.
`solar_elevation` gives the sun's elevation for the albedo target's `min_solar_elevation` filter.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date

#: the quantities a declaration maps onto its file's column names
QUANTITIES = ("sw_in", "sw_up", "lw_up", "rnet", "le", "h", "nee", "gpp", "ustar", "par")
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
    closure: str = "bowen"            # "bowen" (FLUXNET2015's, the Bowen ratio kept) or "none"
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
        if not col:                                  # unmapped, or mapped to "" (the reference's default)
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
    if spec.closure == "none":
        f = np.ones(len(hourly))
    elif spec.closure != "bowen":
        raise ValueError(f'[tower].closure must be "bowen" or "none", not {spec.closure!r}')
    out = hourly.copy()
    out["closure_f"] = f
    out["h_c"] = out["h"] * f
    out["le_c"] = out["le"] * f
    #----- each day's own closure ratio (at least 6 measured daytime hours), for closure_range
    n = both.astype(float).resample("1D").sum()
    tf_d = (hourly["h"] + hourly["le"]).where(both).resample("1D").sum(min_count=1)
    rn_d = hourly["rnet"].where(both).resample("1D").sum(min_count=1)
    ratio = (tf_d / rn_d).where((n >= 6) & (rn_d > 0))
    out["closure_day"] = ratio.reindex(hourly.index.floor("D")).to_numpy()
    return out


def solar_elevation(index_local: pd.DatetimeIndex, lat: float, lon: float, utc_offset_h: float) -> np.ndarray:
    """The sun's elevation [degrees] at the middle of each local hour (NOAA's approximation:
    declination and equation of time from the fractional year; well within a degree)."""
    t = index_local - pd.Timedelta(hours=utc_offset_h) + pd.Timedelta(minutes=30)     # UTC, mid-hour
    doy = t.dayofyear.to_numpy()
    hour = t.hour.to_numpy() + t.minute.to_numpy() / 60.0
    g = 2.0 * np.pi / 365.0 * (doy - 1 + (hour - 12.0) / 24.0)
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    dec = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
           + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = hour * 60.0 + eot + 4.0 * lon                                                 # true solar time [min]
    ha = np.radians(tst / 4.0 - 180.0)
    la = np.radians(lat)
    cosz = np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)
    return np.degrees(np.arcsin(np.clip(cosz, -1.0, 1.0)))


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
