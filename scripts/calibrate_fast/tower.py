# SPDX-License-Identifier: Apache-2.0
"""The observations (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2): the tower as its site TOML declares
it, read by the one reader (scripts/prepare_flux_tower/tower_inputs.py), on its own interval and
the UTC start of each record -- the clock the model's output is on. Facts about the data (the file,
its columns, units, clock, and which values were measured) live in the site TOML; this module holds
the calibration's choices about them ([tower] of calibration.toml): the hours whose forcing was
observed and the daytime threshold.

A flux is used only where the site TOML's rule says it was measured (BCI: FLAG = 1), the radiation
wherever it is present. Each record also carries the closure factor f_d of its day (obsmodels.py:
from the provider's daily sums, gap-filled values included, over days with enough measured records)
and the provider's random uncertainty of a flux where the site declares one (<flux>_randunc).

Incident PAR (quantity `par`, the site's PAR variable) is read as is. `solar_elevation` gives the
sun's elevation for the albedo target's `min_solar_elevation` filter.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "prepare_flux_tower"))
import tower_inputs as ti            # noqa: E402  (the one reader of tower files)
import datarules as DR               # noqa: E402
import obsmodels as OM               # noqa: E402

#: the calibration's quantities, from the site TOML's variables and fluxes
QUANTITIES = {"sw_in": "SWdown", "sw_up": "SW_out", "lw_up": "LW_out", "rnet": "Rnet", "le": "LE", "h": "H",
              "nee": "NEE", "gpp": "GPP", "reco": "RECO", "ustar": "USTAR", "par": "PAR",
              "tair": "Tair", "vpd": "VPD", "rain": "Rainf", "wind": "Wind",
              "gpp_dt": "GPP_DT", "reco_dt": "RECO_DT"}


@dataclass
class TowerSpec:
    site: str                         # the site TOML (the facts)
    forcing: str | None = None        # the forcing file, for the observed-forcing mask
    forcing_qc: tuple = ("LWdown_qc", "Wind_qc")
    forcing_qc_val: tuple = ("Wind_qc",)   # validation windows: a tower year without longwave still validates
    forcing_grid: int = 1
    daytime_sw: float = 10.0          # [W m-2] daytime = incoming shortwave above this
    closure: dict = field(default_factory=lambda: {"window_days": 15, "min_measured": 0.7, "min_days": 5})
    utc_offset_h: float = 0.0         # the site's clock [h ahead of UTC], from the site TOML
    step: float = 3600.0              # the tower's interval [s], from the site TOML
    lat: float | None = None          # the site's coordinates, from the site TOML
    lon: float | None = None
    leaf_on_months: list = field(default_factory=lambda: list(range(1, 13)))
    provider: dict = field(default_factory=dict)  # the site TOML's [provider]
    sensor_height: float = 0.0        # [m] the tower's temperature and humidity sensor
    report: dict = field(default_factory=dict)   # the reader's checks (tower_inputs.read_standard)

    @classmethod
    def from_site(cls, site_toml, **choices) -> "TowerSpec":
        """The spec with the site TOML's clock, interval, coordinates, leaf-on months and provider
        filled in."""
        site = ti.read_site(str(site_toml))
        return cls(site=str(site_toml), utc_offset_h=site.utc_offset, step=site.timestep,
                   lat=site.latitude, lon=site.longitude, leaf_on_months=site.leaf_on_months,
                   provider=dict(site.provider), sensor_height=site.tq_height, **choices)


def observations(spec: TowerSpec) -> pd.DataFrame:
    """The tower's records on the UTC start of each interval: one column per quantity, NaN where not
    measured, the day's closure factor (closure_f), the provider's random uncertainty of a flux
    (<flux>_randunc) and the days since frost. The reader's and the closure's reports go to
    spec.report."""
    table = ti.read_standard(ti.read_site(spec.site))
    spec.report = table.report
    obs = pd.DataFrame(index=table.values.index)
    for q, name in QUANTITIES.items():
        if name in table.values:
            obs[q] = table.values[name].where(table.measured[name])
        else:
            obs[q] = np.nan
    if obs["vpd"].isna().all() and "RH" in table.values:      # VPD from RH where the site gives no VPD
        tc = obs["tair"] - 273.15
        rh = table.values["RH"].where(table.measured["RH"])
        obs["vpd"] = (1.0 - rh) * ti.SATURATION_CURVES["alduchov_eskridge"](tc)
    if obs["tair"].notna().any():
        obs["days_since_frost"] = DR.days_since_frost(obs["tair"], spec.utc_offset_h)
    names = {v: k for k, v in QUANTITIES.items()}
    for flux in table.uncertainty.columns:
        obs[f"{names[flux]}_randunc"] = table.uncertainty[flux]
    obs["closure_f"], spec.report["closure"] = OM.closure_factor(table.values, table.measured,
                                                                 spec.utc_offset_h, spec.closure)
    return obs


def local_days(index: pd.DatetimeIndex, utc_offset_h: float) -> pd.DatetimeIndex:
    """The local calendar day of each UTC record start."""
    return (index + pd.Timedelta(hours=utc_offset_h)).floor("D")


def solar_elevation(index_utc: pd.DatetimeIndex, lat: float, lon: float, step: float) -> np.ndarray:
    """The sun's elevation [degrees] at the middle of each record (NOAA's approximation:
    declination and equation of time from the fractional year; well within a degree)."""
    t = index_utc + pd.Timedelta(seconds=step / 2.0)
    doy = t.dayofyear.to_numpy()
    hour = t.hour.to_numpy() + t.minute.to_numpy() / 60.0 + t.second.to_numpy() / 3600.0
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


def forcing_span(spec: TowerSpec):
    """The forcing file's first and last record (None without a forcing file)."""
    if not spec.forcing:
        return None
    with Dataset(spec.forcing) as ds:
        t = ds["time"]
        first, last = num2date([t[0], t[-1]], t.units, only_use_cftime_datetimes=False)
    return pd.Timestamp(first.strftime("%Y-%m-%d %H:%M:%S")), pd.Timestamp(last.strftime("%Y-%m-%d %H:%M:%S"))


def forcing_observed(spec: TowerSpec, index: pd.DatetimeIndex, qc=None) -> pd.Series:
    """True for tower records whose forcing was observed (qc 0 in every listed variable, `qc` or the
    spec's calibration list): every forcing record inside the tower's interval, or the one forcing
    record containing it when the forcing is coarser. All True when no forcing file is declared."""
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
    fstep = float(np.median(np.diff(when.values).astype("timedelta64[s]").astype(float)))
    if fstep >= spec.step:                       # a coarser (or equal) forcing: the record containing each tower record
        f = pd.Series(ok, index=when)
        return pd.Series(f.reindex(index.floor(f"{int(fstep)}s")).to_numpy() == True, index=index)  # noqa: E712
    bins = pd.Series(ok.astype(float), index=when).groupby(when.floor(f"{int(spec.step)}s")).agg(["sum", "count"])
    need = int(round(spec.step / fstep))
    obs = (bins["sum"] == need) & (bins["count"] == need)
    return obs.reindex(index, fill_value=False)
