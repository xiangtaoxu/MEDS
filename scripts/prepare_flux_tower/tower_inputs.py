# SPDX-License-Identifier: Apache-2.0
"""tower_inputs.py -- read a flux tower's meteorology as its site TOML declares it.

The site TOML DECLARES what the data are: the file and its format, the location, the clock (UTC
offset and which end of each interval a stamp marks), the sensor heights, and for every variable
its column and units. This module reads the file under that declaration and converts every
variable to MEDS units; tower_checks.py then VALIDATES the declaration against the sun and the data.
Nothing is sniffed: a missing declaration is an error naming the key.

Three input formats share one column map and differ only in their defaults:
  csv             any delimited table with one timestamp column (e.g. Barro Colorado Island)
  ameriflux_base  AmeriFlux BASE: '#' header lines, TIMESTAMP_START/END as YYYYMMDDHHMM, -9999
  fluxnet         FLUXNET/ONEFlux: as BASE, with gap-filled *_F columns whose *_QC column says
                  which values the provider filled
"""
import os
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

try:
    import tomllib
except ModuleNotFoundError:          # Python < 3.11
    import tomli as tomllib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import meds_forcing_file as mff  # noqa: E402  (the writer, in scripts/)
import tower_conversions as conv  # noqa: E402  (the model's conversions, this folder)

FORMAT_DEFAULTS = {
    "csv":            dict(timestamp=None, timestamp_format=None, missing=[], comment=None),
    "ameriflux_base": dict(timestamp="TIMESTAMP_START", timestamp_format="%Y%m%d%H%M", missing=[-9999],
                           comment="#"),
    "fluxnet":        dict(timestamp="TIMESTAMP_START", timestamp_format="%Y%m%d%H%M", missing=[-9999],
                           comment=None),
}
TIMESTAMP_STAMP = {"TIMESTAMP_START": "begin", "TIMESTAMP_END": "end"}

# The variables a site TOML may declare: the forcing (the first seven are required, with RH or VPD
# for the humidity) and the ones read for the reports only.
FORCING_VARIABLES = ("Tair", "RH", "VPD", "PSurf", "Rainf", "SWdown", "LWdown", "Wind")
REPORT_VARIABLES = ("PAR",)

# Unit conversions to MEDS units. Rain per interval needs the interval length, applied by the caller.
UNITS = {
    "Tair":   {"degC": lambda x: x + 273.15, "K": lambda x: x},
    "RH":     {"%": lambda x: x / 100.0, "1": lambda x: x},
    "VPD":    {"kPa": lambda x: x * 1000.0, "hPa": lambda x: x * 100.0, "Pa": lambda x: x},
    "PSurf":  {"kPa": lambda x: x * 1000.0, "hPa": lambda x: x * 100.0, "Pa": lambda x: x},
    "Rainf":  {"mm": None, "mm s-1": lambda x: x, "kg m-2 s-1": lambda x: x},
    "SWdown": {"W m-2": lambda x: x},
    "LWdown": {"W m-2": lambda x: x},
    "Wind":   {"m s-1": lambda x: x},
    "PAR":    {"umol m-2 s-1": lambda x: x},
}

# The saturation curves a provider may have made its VPD with [Pa, T in degC]. V3 checks the
# declared one and names the best fit; RH is recovered through the declared one where it is missing.
SATURATION_CURVES = {
    "bolton":            lambda tc: 611.2 * np.exp(17.67 * tc / (tc + 243.5)),
    "alduchov_eskridge": lambda tc: 610.94 * np.exp(17.625 * tc / (tc + 243.04)),
    "tetens":            lambda tc: 610.78 * np.exp(17.27 * tc / (tc + 237.3)),
    "buck":              lambda tc: 611.21 * np.exp((18.678 - tc / 234.5) * tc / (257.14 + tc)),
    "campbell_norman":   lambda tc: 611.0 * np.exp(17.502 * tc / (tc + 240.97)),
}


@dataclass
class Site:
    """A site TOML, read and checked for completeness."""
    toml_path: str
    input_format: str
    input_path: str
    timestamp: str
    timestamp_format: str
    missing: list
    comment: str
    name: str
    latitude: float
    longitude: float
    elevation: float
    utc_offset: float
    stamp: str
    timestep: float
    tq_height: float
    wind_height: float
    pressure_height: float
    variables: dict
    gapfill: dict = field(default_factory=dict)


def _require(table, key, where):
    if key not in table:
        raise SystemExit(f"ERROR: the site TOML has no {where}.{key}; it is required (MEDS never guesses it)")
    return table[key]


def read_site(path):
    """The site TOML at `path`; relative input paths resolve against its folder."""
    with open(path, "rb") as fh:
        t = tomllib.load(fh)
    inp, site, clock, heights = (_require(t, k, "the file") for k in ("input", "site", "clock", "heights"))
    fmt = _require(inp, "format", "input")
    if fmt not in FORMAT_DEFAULTS:
        raise SystemExit(f"ERROR: input.format = {fmt!r}; choose from {list(FORMAT_DEFAULTS)}")
    defaults = FORMAT_DEFAULTS[fmt]
    base = os.path.dirname(os.path.abspath(path))
    input_path = os.path.expanduser(_require(inp, "path", "input"))
    if not os.path.isabs(input_path):
        input_path = os.path.join(base, input_path)
    timestamp = inp.get("timestamp", defaults["timestamp"])
    timestamp_format = inp.get("timestamp_format", defaults["timestamp_format"])
    if timestamp is None or timestamp_format is None:
        raise SystemExit("ERROR: a csv input needs input.timestamp and input.timestamp_format")
    stamp = _require(clock, "stamp", "clock")
    if stamp not in ("begin", "end"):
        raise SystemExit(f"ERROR: clock.stamp = {stamp!r}; it is 'begin' or 'end'")
    if timestamp in TIMESTAMP_STAMP and TIMESTAMP_STAMP[timestamp] != stamp:
        raise SystemExit(f"ERROR: input.timestamp = {timestamp} marks the {TIMESTAMP_STAMP[timestamp]} of each "
                         f"interval, but clock.stamp = {stamp!r}")
    variables = _require(t, "variables", "the file")
    for name in variables:
        if name not in FORCING_VARIABLES + REPORT_VARIABLES:
            raise SystemExit(f"ERROR: variables.{name} is not one of {FORCING_VARIABLES + REPORT_VARIABLES}")
        spec = variables[name]
        _require(spec, "column", f"variables.{name}")
        units = _require(spec, "units", f"variables.{name}")
        if units not in UNITS[name]:
            raise SystemExit(f"ERROR: variables.{name}.units = {units!r}; choose from {list(UNITS[name])}")
    for name in ("Tair", "PSurf", "Rainf", "SWdown", "Wind"):
        _require(variables, name, "variables")
    if "RH" not in variables and "VPD" not in variables:
        raise SystemExit("ERROR: the site TOML declares no humidity: give variables.RH, or variables.VPD with its curve")
    if "VPD" in variables:
        curve = _require(variables["VPD"], "curve", "variables.VPD")
        if curve not in SATURATION_CURVES:
            raise SystemExit(f"ERROR: variables.VPD.curve = {curve!r}; choose from {list(SATURATION_CURVES)}")
    gapfill = t.get("gapfill", {})
    for key in gapfill:
        if key != "short_gap_max":
            raise SystemExit(f"ERROR: gapfill.{key} is not a setting; [gapfill] takes short_gap_max only. The "
                             f"longwave is always the model's synthesis regressed onto the tower, and a long gap "
                             f"in another variable the mean diurnal variation; to fill from ERA5-Land or another "
                             f"source, fill the tower file before the build.")
    return Site(
        toml_path=os.path.abspath(path), input_format=fmt, input_path=input_path, timestamp=timestamp,
        timestamp_format=timestamp_format, missing=list(inp.get("missing", defaults["missing"])),
        comment=inp.get("comment", defaults["comment"]), name=site.get("name", ""),
        latitude=float(_require(site, "latitude", "site")), longitude=float(_require(site, "longitude", "site")),
        elevation=float(_require(site, "elevation", "site")), utc_offset=float(_require(clock, "utc_offset", "clock")),
        stamp=stamp, timestep=float(_require(clock, "timestep", "clock")),
        tq_height=float(_require(heights, "tq_height", "heights")),
        wind_height=float(_require(heights, "wind_height", "heights")),
        pressure_height=float(_require(heights, "pressure_height", "heights")),
        variables=variables, gapfill=gapfill)


@dataclass
class TowerData:
    """The tower's variables in MEDS units on the SOURCE clock, one row per stamp.

    values    DataFrame: Tair [K], RH [1], VPD [Pa], PSurf [Pa, at the barometer], Rainf [kg m-2 s-1],
              SWdown, LWdown [W m-2], Wind [m s-1], PAR [umol m-2 s-1]; NaN where missing
    provider  DataFrame of the same columns: True where the provider filled the value (FLUXNET _QC > 0)
    """
    site: Site
    values: pd.DataFrame
    provider: pd.DataFrame


def read_tower(site):
    """Read the declared columns of the input file and convert them to MEDS units."""
    if not os.path.exists(site.input_path):
        raise SystemExit(f"ERROR: the tower file {site.input_path} does not exist")
    columns = {name: spec["column"] for name, spec in site.variables.items()}
    qc_columns = {}
    if site.input_format == "fluxnet":
        qc_columns = {name: col + "_QC" for name, col in columns.items() if col.endswith("_F")}
    wanted = [site.timestamp] + list(columns.values()) + list(qc_columns.values())
    raw = pd.read_csv(site.input_path, comment=site.comment, na_values=[str(m) for m in site.missing],
                      usecols=lambda c: c in wanted, dtype={site.timestamp: str}, low_memory=False)
    missing_cols = [c for c in wanted if c not in raw.columns]
    if missing_cols:
        raise SystemExit(f"ERROR: {site.input_path} has no column(s) {missing_cols}")
    stamps = pd.to_datetime(raw[site.timestamp].str.strip(), format=site.timestamp_format)
    values = pd.DataFrame(index=pd.DatetimeIndex(stamps, name="stamp"))
    provider = pd.DataFrame(index=values.index)
    for name, col in columns.items():
        x = pd.to_numeric(raw[col], errors="coerce").to_numpy(dtype=float, copy=True)
        for m in site.missing:
            x[x == float(m)] = np.nan
        convert = UNITS[name][site.variables[name]["units"]]
        if convert is None:                              # rain per interval -> a mean rate
            x = x / site.timestep
        else:
            x = convert(x)
        values[name] = x
        flags = np.zeros(len(x), dtype=bool)
        if name in qc_columns:
            q = pd.to_numeric(raw[qc_columns[name]], errors="coerce").to_numpy()
            flags = np.nan_to_num(q, nan=0.0) > 0
        provider[name] = flags
    return TowerData(site=site, values=values, provider=provider)


def interval_bounds(stamps_utc, stamp, timestep):
    """[start, end) of each record's interval, from its UTC stamp and the stamp convention."""
    stamps_utc = np.asarray(stamps_utc, dtype="datetime64[s]")
    dt = np.timedelta64(int(round(timestep)), "s")
    start = stamps_utc if stamp == "begin" else stamps_utc - dt
    return start, start + dt


def to_utc(stamps_local, utc_offset):
    """UTC stamps from stamps on a clock `utc_offset` hours ahead of UTC (Panama: -5)."""
    shift = np.timedelta64(int(round(utc_offset * 3600.0)), "s")
    return np.asarray(stamps_local, dtype="datetime64[s]") - shift


def mean_cosz_of_intervals(stamps_utc, site):
    """The model's window-mean cos z over each record's interval."""
    start, _ = interval_bounds(stamps_utc, site.stamp, site.timestep)
    return conv.window_mean_cosz(start, site.timestep, site.latitude, site.longitude)
