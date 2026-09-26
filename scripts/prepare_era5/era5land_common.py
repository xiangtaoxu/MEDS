# SPDX-License-Identifier: Apache-2.0
"""Shared definitions for the ERA5-Land download and post-processing scripts.

The downloaders only fetch source files, unchanged:
  download_era5land_cds.py   the Copernicus Climate Data Store (GRIB by default, for a box)
  download_era5land_gdex.py  NCAR GDEX dataset d633008 (global 5-day NetCDF files, a shared raw pool)
postprocess_era5land.py turns either source's raw files into uniform NetCDF box files (decode, cut,
merge, trim), keeping ERA5-Land's own variable names, units and time stamps. Conversion to the MEDS
forcing format and the split into regions come after that.
"""
import calendar
import datetime as dt
import json
import math
import os
import time

import numpy as np

# The eight ERA5-Land variables MEDS forcing needs.
#   short name -> (CDS request name, GDEX product, GDEX file code, GRIB step type)
# 'accum' fields are accumulated since 00 UTC in both sources and are de-accumulated downstream.
VARIABLES = {
    "t2m":  ("2m_temperature",                      "instan", "t2m_167",  "instant"),
    "d2m":  ("2m_dewpoint_temperature",             "instan", "d2m_168",  "instant"),
    "sp":   ("surface_pressure",                    "instan", "sp_134",   "instant"),
    "u10":  ("10m_u_component_of_wind",             "instan", "u10_165",  "instant"),
    "v10":  ("10m_v_component_of_wind",             "instan", "v10_166",  "instant"),
    "tp":   ("total_precipitation",                 "accumu", "tp_228",   "accum"),
    "ssrd": ("surface_solar_radiation_downwards",   "accumu", "ssrd_169", "accum"),
    "strd": ("surface_thermal_radiation_downwards", "accumu", "strd_175", "accum"),
}

GRID_STEP = 0.1   # [deg] ERA5-Land native grid spacing


def parse_variables(text):
    """'t2m,tp' -> ['t2m', 'tp']; 'all' -> every variable in catalogue order."""
    if text.strip().lower() == "all":
        return list(VARIABLES)
    names = [v.strip() for v in text.split(",") if v.strip()]
    unknown = [v for v in names if v not in VARIABLES]
    if unknown:
        raise SystemExit(f"unknown variable(s) {unknown}; choose from {list(VARIABLES)} or 'all'")
    return names


def parse_bbox(text):
    """'N,W,S,E' in degrees (longitudes -180..180; W > E means the box crosses the antimeridian)."""
    try:
        north, west, south, east = (float(x) for x in text.split(","))
    except ValueError:
        raise SystemExit(f"--bbox must be N,W,S,E (got {text!r})")
    if not (-90.0 <= south < north <= 90.0):
        raise SystemExit(f"--bbox needs -90 <= S < N <= 90 (got N={north}, S={south})")
    for lon in (west, east):
        if not -180.0 <= lon <= 180.0:
            raise SystemExit(f"--bbox longitudes must be in -180..180 (got {lon})")
    return north, west, south, east


def align_bbox(bbox):
    """Snap the box outward to multiples of the 0.1 deg grid. The CDS returns a shifted grid for
    unaligned edges, and snapping outward never drops a requested cell."""
    north, west, south, east = bbox
    snap_up = lambda x: math.ceil(round(x / GRID_STEP, 6)) * GRID_STEP
    snap_down = lambda x: math.floor(round(x / GRID_STEP, 6)) * GRID_STEP
    return (round(min(snap_up(north), 90.0), 1), round(max(snap_down(west), -180.0), 1),
            round(max(snap_down(south), -90.0), 1), round(min(snap_up(east), 180.0), 1))


def parse_dates(start, end):
    """Inclusive calendar-day range. There is deliberately no default period."""
    try:
        d0, d1 = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    except ValueError as err:
        raise SystemExit(f"--start/--end must be YYYY-MM-DD ({err})")
    if d1 < d0:
        raise SystemExit("--end precedes --start")
    return d0, d1


def interval_stamps(d0, d1):
    """End-stamped hourly intervals covering the days d0..d1: the first stamp is d0 01:00 and the last
    is (d1 + 1 day) 00:00, which closes the last hour of d1."""
    first = dt.datetime.combine(d0, dt.time(1))
    last = dt.datetime.combine(d1 + dt.timedelta(days=1), dt.time(0))
    return first, last


def month_starts(first, last):
    """First day of every month touched by the datetimes first..last."""
    months, y, m = [], first.year, first.month
    while (y, m) <= (last.year, last.month):
        months.append(dt.date(y, m, 1))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return months


# --- NCAR GDEX d633008 file layout ----------------------------------------------------------------
GDEX_BASE_URL = "https://data.gdex.ucar.edu/d633008"
GDEX_BLOCK_START_DAYS = (1, 6, 11, 16, 21, 26)   # one file per variable per 5 days; the last runs to month end


def gdex_blocks(first, last):
    """GDEX 5-day blocks overlapping the end-stamped interval first..last, as (month, block_start,
    block_end); block_start is 01:00 of its first day, block_end 00:00 after its last (files are
    end-stamped, so the stamp closing a month is inside that month's last file)."""
    out = []
    for month in month_starts(first, last):
        ndays = calendar.monthrange(month.year, month.month)[1]
        for i, d in enumerate(GDEX_BLOCK_START_DAYS):
            b0 = dt.datetime(month.year, month.month, d, 1)
            end_day = GDEX_BLOCK_START_DAYS[i + 1] if i + 1 < len(GDEX_BLOCK_START_DAYS) else ndays + 1
            b1 = dt.datetime(month.year, month.month, 1) + dt.timedelta(days=end_day - 1)
            if b0 <= last and b1 >= first:
                out.append((month, b0, b1))
    return out


def gdex_relpath(variable, month, b0, b1):
    """Path of a GDEX file relative to the dataset root; the raw pool mirrors this layout."""
    _, product, code, _ = VARIABLES[variable]
    name = f"e5land.oper.fc.sfc.{product}.{code}.{b0:%Y%m%d%H}-{b1:%Y%m%d%H}.nc"
    return f"e5land.oper.fc.sfc.{product}/{month:%Y%m}/{name}"


def stamps_between(b0, b1):
    return int((b1 - b0).total_seconds() // 3600) + 1


# --- box selection on a source grid -----------------------------------------------------------------
def select_rows(lat, north, south, tol=1e-6):
    """Indices of latitudes inside [south, north], in source order."""
    return np.flatnonzero((lat <= north + tol) & (lat >= south - tol))


def select_cols(lon, west, east, tol=1e-6):
    """Indices of the box's longitude columns, west to east, on a monotonic source axis (any frame:
    -180..180, 0..360, or unwrapped past 360). Works for a global axis whose seam the box crosses and
    for boxes crossing the antimeridian. Returns (indices, longitudes as output)."""
    lon = np.asarray(lon, dtype=float)
    width = (east - west) % 360.0
    if width == 0.0 and east != west:
        width = 360.0
    w = lon[0] + ((west - lon[0]) % 360.0)
    if w > lon[-1] + tol:        # the box begins just west of a regional source axis
        w -= 360.0
    ext = np.concatenate([lon - 360.0, lon, lon + 360.0])
    idx = np.tile(np.arange(len(lon)), 3)
    sel = (ext >= w - tol) & (ext <= w + width + tol)
    cols, vals = idx[sel], ext[sel]
    out = np.where(vals > 180.0, vals - 360.0, np.where(vals <= -180.0, vals + 360.0, vals))
    if np.any(np.diff(out) <= 0):   # crosses the antimeridian: keep a monotonic 0..360 axis instead
        out = vals % 360.0
    return cols, np.round(out, 4)


def use_group_umask():
    """Outputs in a shared group directory stay writable by the group, so collaborators can extend
    each other's downloads and products."""
    os.umask(0o002)


class RunLog:
    """Append-only JSON-lines log of every transfer, kept next to the outputs for provenance and
    throughput records."""

    def __init__(self, out_dir, name):
        self.path = os.path.join(out_dir, name)

    def write(self, **record):
        record["logged_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        with open(self.path, "a") as fh:
            fh.write(json.dumps(record) + "\n")


def human_rate(nbytes, seconds):
    return f"{nbytes / 1e6:.1f} MB in {seconds:.1f} s ({nbytes / 1e6 / max(seconds, 1e-9):.1f} MB/s)"


def monotonic_seconds():
    return time.monotonic()
