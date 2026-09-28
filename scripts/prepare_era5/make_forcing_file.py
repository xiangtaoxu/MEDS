#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""make_forcing_file.py -- write a MEDS single forcing file, the (time, grid) NetCDF that
[forcing].format = "netcdf" reads (src/forcing/meds_met_driver.f90), from either source of ERA5-Land
data this folder produces:

  --data-path DIR   an ED_ERA5land archive (build_era5land_archive.py), for the days --start..--end.
                    The archive already holds hourly rates; this tool picks each site's nearest valid
                    cell within --max-distance-km, the same rule the model's archive reader uses, and
                    takes the cell's orography from the static file.
  --box-dir DIR     the box files of postprocess_era5land.py --split none: a small download from the
                    CDS or GDEX, without an archive. Their fluxes are still accumulated since 00 UTC,
                    and this tool de-accumulates them.

Both inputs give the same file (MEDS_FORCING_DESIGN.md §7.1):
  * dims (time, grid); coordinates time (seconds since the first record), latitude, longitude,
    elevation per grid point;
  * Tair [K], PSurf [Pa], Qair [kg/kg] from the dewpoint by the model's own Bolton (1980) saturation
    form, the 10 m wind vector u10, v10 [m/s] and its speed Wind = sqrt(u10^2 + v10^2) (the model
    floors the speed itself, for every source), and the hour-mean fluxes Rainf
    [kg m-2 s-1], SWdown (total; the model partitions it) and LWdown [W m-2];
  * end-stamped hourly records (avg_convention = "end"), UTC.
No CO2air: CO2 is not meteorology, and the model takes it from [forcing].co2_source (a constant, or
a MEDS CO2 file such as data/co2/), never from this file -- it rejects a file that carries one.
MEDS never gap-fills: a missing value is an error, here as in the model.

De-accumulation of the box files (the 00Z trap, MEDS_FORCING_DESIGN.md §7.3): ERA5-Land accumulates
from 00 UTC and resets daily, and the 00:00 stamp holds the WHOLE previous day. Ordered by valid
time, per-hour(H) = raw(H) - raw(H-1) everywhere except at 01:00 UTC, where per-hour = raw(01:00).
A leading sample that is not 01:00 cannot be differenced and is dropped (a boundary, not a gap).
Negative packing noise is clipped to 0 for precipitation and shortwave, as in the archive.

Usage (Ithaca NY, calendar year 2024, the example_biophysics recycle window):
  # from an archive
  python make_forcing_file.py --data-path $ERA5LAND_ROOT/ED_ERA5land --start 2024-01-01 \\
      --end 2024-12-31 --lat 42.44 --lon -76.50 --out ithaca_forcing.nc
  # without an archive: a small CDS download, decoded into box files, then this tool
  python download_era5land_cds.py --bbox 42.5,-76.6,42.4,-76.4 --start 2024-01-01 --end 2024-12-31 \\
      --variables all --out-dir raw_ithaca
  python postprocess_era5land.py --source cds --raw-dir raw_ithaca --bbox 42.5,-76.6,42.4,-76.4 \\
      --start 2024-01-01 --end 2024-12-31 --variables all --split none --out-dir box_ithaca
  python make_forcing_file.py --box-dir box_ithaca --lat 42.44 --lon -76.50 --out ithaca_forcing.nc
  # several locations -> one grid index each; or every cell of the box (box files only)
  python make_forcing_file.py --data-path ... --start ... --end ... --cells 42.44,-76.50 44.00,-77.00 --out multi.nc
  python make_forcing_file.py --box-dir box_ithaca --all-cells --out multi.nc

Dependencies: numpy + netCDF4.
"""
import argparse
import datetime as dt
import glob
import math
import os

import numpy as np
from netCDF4 import Dataset, num2date

import era5land_common as common

RHO_W = 1000.0          # [kg/m3] water density (1 m of water = 1000 kg/m2)
SEC_PER_HOUR = 3600.0
EARTH_RADIUS_KM = 6371.0

# ERA5-Land NetCDF variable names in the box files (NOT the GRIB shortnames).
NC_NAMES = dict(t2m="t2m", d2m="d2m", sp="sp", u10="u10", v10="v10",
                tp="tp", ssrd="ssrd", strd="strd")


# ---------------------------------------------------------------------------------------------
# Saturation vapour pressure -- Bolton (1980), IDENTICAL to meds_thermo%sat_vapor_pressure so the
# file's Qair reconciles with the Fortran reader (design §4.3, test §9.7). Input K, output Pa.
# ---------------------------------------------------------------------------------------------
def sat_vapor_pressure(t_k):
    tc = t_k - 273.15
    return 611.2 * np.exp(17.67 * tc / (tc + 243.5))


def dewpoint_to_specific_humidity(td_k, p_pa):
    """q [kg/kg] from dewpoint Td [K] and surface pressure P [Pa]. The actual vapour pressure is the
    saturation vapour pressure evaluated AT the dewpoint: e = e_sat(Td)."""
    e = sat_vapor_pressure(td_k)
    return 0.622 * e / (p_pa - 0.378 * e)


# ---------------------------------------------------------------------------------------------
# De-accumulate an ERA5-Land accumulated field to a per-hour amount (interval [H-1, H], stamped at H).
# `hours_utc` is the integer UTC hour of each (ascending) timestamp. Rule (fact sheet C):
#   H==01 UTC -> raw as-is (first step of the daily accumulation, implicit 0 at 00:00)
#   else      -> raw(H) - raw(H-1)   (the 00:00 stamp, = whole prior day, then yields prior day's [23,00])
# The very first sample, if not 01:00 UTC, cannot be differenced -> NaN (drop or ignore).
# GRIB packing noise can make a difference slightly negative; for precipitation and shortwave it is
# clipped to 0, exactly as build_era5land_archive.py does. Only negatives: a positive threshold
# (formerly 1e-5 m) zeroed 0.95% of all rain and 55% of wet hours (MEDS_FORCING_DESIGN.md §7.3).
# ---------------------------------------------------------------------------------------------
def deaccumulate_hourly(accum, hours_utc, clip_negative):
    accum = np.asarray(accum, dtype=float)
    out = np.empty_like(accum)
    out[:] = np.nan
    n = accum.shape[0]
    for i in range(n):
        if hours_utc[i] == 1:
            out[i] = accum[i]
        elif i > 0:
            out[i] = accum[i] - accum[i - 1]
        # else i==0 and hour!=1 -> leave NaN (no previous sample to difference)
    if clip_negative:                                   # keep NaN as NaN
        out[np.isfinite(out) & (out < 0.0)] = 0.0
    return out


# ---------------------------------------------------------------------------------------------
# BOX FILES. Returns times (list[datetime], UTC), 1-D lat/lon arrays, and a dict of (time, lat, lon)
# arrays for each needed variable. Each variable must come from exactly one file, and every file must
# share the time axis and grid. Handles 'valid_time'|'time' and stray singleton dims (e.g.
# 'number'/'expver' the new CDS sometimes adds).
# ---------------------------------------------------------------------------------------------
def read_box_files(paths):
    times = lat = lon = None
    data = {}
    for path in paths:
        ds = Dataset(path)
        tname = "valid_time" if "valid_time" in ds.variables else "time"
        tvar = ds.variables[tname]
        t = num2date(tvar[:], tvar.units,
                     getattr(tvar, "calendar", "standard"),
                     only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        t = list(np.atleast_1d(t))
        latname = "latitude" if "latitude" in ds.variables else "lat"
        lonname = "longitude" if "longitude" in ds.variables else "lon"
        la = np.atleast_1d(ds.variables[latname][:]).astype(float)
        lo = np.atleast_1d(ds.variables[lonname][:]).astype(float)
        if times is None:
            times, lat, lon = t, la, lo
        elif t != times or not (np.array_equal(la, lat) and np.array_equal(lo, lon)):
            raise SystemExit(f"ERROR: {path}: time axis or grid differs from {paths[0]}")
        for k, nc in NC_NAMES.items():
            if nc not in ds.variables:
                continue
            if k in data:
                raise SystemExit(f"ERROR: {nc} appears in more than one input file "
                                 f"(write the box files with postprocess_era5land.py --split none)")
            arr = np.ma.filled(ds.variables[nc][:], np.nan).astype(float)
            # collapse to (time, lat, lon): drop any leading singleton dims (number/expver)
            while arr.ndim > 3:
                arr = arr[0]
            if arr.ndim == 2:                   # (lat, lon): a single time step
                arr = arr[np.newaxis, :, :]
            data[k] = arr
        ds.close()
    missing = [nc for k, nc in NC_NAMES.items() if k not in data]
    if missing:
        raise SystemExit(f"ERROR: variable(s) {missing} not found in the input files")
    return times, lat, lon, data


def nearest_index(lat, lon, tlat, tlon):
    """Nearest (ilat, ilon) grid index to a target lat/lon (planar distance; cells are ~0.1 deg)."""
    ilat = int(np.argmin(np.abs(lat - tlat)))
    ilon = int(np.argmin(np.abs(lon - tlon)))
    return ilat, ilon


def select_box_cells(lat, lon, args):
    """Return a list of (lat, lon, ilat, ilon) selected grid points -> the `grid` dim."""
    cells = []
    if args.all_cells:
        for i, la in enumerate(lat):
            for j, lo in enumerate(lon):
                cells.append((float(la), float(lo), i, j))
    else:
        for tla, tlo in targets(args):
            i, j = nearest_index(lat, lon, tla, tlo)
            cells.append((float(lat[i]), float(lon[j]), i, j))
    return cells


def from_box_files(args):
    """(times, cells, per-cell fields) from box files; fluxes de-accumulated and unit-converted."""
    paths = sorted(glob.glob(os.path.join(args.box_dir, "*.nc")))
    if not paths:
        raise SystemExit(f"ERROR: no NetCDF files in {args.box_dir}")
    times, lat, lon, data = read_box_files(paths)
    hours_utc = np.array([t.hour for t in times], dtype=int)
    cells = select_box_cells(lat, lon, args)
    fields = []
    for cla, clo, i, j in cells:
        t2m = data["t2m"][:, i, j]; d2m = data["d2m"][:, i, j]; sp = data["sp"][:, i, j]
        u10 = data["u10"][:, i, j]; v10 = data["v10"][:, i, j]
        # de-accumulate fluxes (per-hour amount), then convert to rate / mean flux
        tp_hr = deaccumulate_hourly(data["tp"][:, i, j], hours_utc, clip_negative=True)
        ssrd_hr = deaccumulate_hourly(data["ssrd"][:, i, j], hours_utc, clip_negative=True)
        strd_hr = deaccumulate_hourly(data["strd"][:, i, j], hours_utc, clip_negative=False)
        fields.append(dict(
            Tair=t2m, PSurf=sp, Qair=dewpoint_to_specific_humidity(d2m, sp),
            u10=u10, v10=v10, Wind=np.hypot(u10, v10),
            Rainf=tp_hr * RHO_W / SEC_PER_HOUR,       # [m/hr] -> [kg/m2/s]
            SWdown=ssrd_hr / SEC_PER_HOUR,            # [J/m2/hr] -> [W/m2] (total; reader partitions)
            LWdown=strd_hr / SEC_PER_HOUR))           # [J/m2/hr] -> [W/m2]
    grid = [(cla, clo, args.elevation) for cla, clo, _, _ in cells]
    return times, grid, fields, "ERA5-Land hourly (reanalysis-era5-land)"


# ---------------------------------------------------------------------------------------------
# ARCHIVE. The archive's files already hold hourly rates, end-stamped, on the global grid of its
# static file; a site takes the nearest cell with data within max_distance_km, as the model does.
# ---------------------------------------------------------------------------------------------
def great_circle_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    c = np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dl)
    return EARTH_RADIUS_KM * np.arccos(np.clip(c, -1.0, 1.0))


def select_archive_cells(data_path, args):
    """Return a list of (lat, lon, elevation, row, col) for the targets' nearest valid cells."""
    with Dataset(common.static_file(data_path)) as st:
        lat = st["lat"][:].astype(float); lon = st["lon"][:].astype(float)
        cells = []
        for tla, tlo in targets(args):
            # a window of cells around the target is enough: max_distance_km is tens of km
            dlat = args.max_distance_km / 111.0 + 0.2
            rows = np.nonzero(np.abs(lat - tla) <= dlat)[0]
            dlon = dlat / max(math.cos(math.radians(tla)), 0.01)
            dl = (lon - tlo + 180.0) % 360.0 - 180.0
            cols = np.nonzero(np.abs(dl) <= dlon)[0]
            if rows.size == 0 or cols.size == 0:
                raise SystemExit(f"ERROR: ({tla}, {tlo}) is outside the archive grid")
            valid = st["valid"][rows.min():rows.max() + 1, :][:, cols] > 0
            la2, lo2 = np.meshgrid(lat[rows.min():rows.max() + 1], lon[cols], indexing="ij")
            dist = np.where(valid, great_circle_km(tla, tlo, la2, lo2), np.inf)
            k = np.unravel_index(np.argmin(dist), dist.shape)
            if not np.isfinite(dist[k]) or dist[k] > args.max_distance_km:
                raise SystemExit(f"ERROR: no valid ED_ERA5land cell within {args.max_distance_km} km "
                                 f"of ({tla}, {tlo})")
            r, c = rows.min() + k[0], cols[k[1]]
            elev = float(st["elevation"][r, c])
            cells.append((float(lat[r]), float(lon[c]), elev, int(r), int(c)))
            print(f"cell ({lat[r]:.3f}, {lon[c]:.3f}), {dist[k]:.2f} km from ({tla}, {tlo}); "
                  f"orography {elev:.1f} m")
    return cells


def from_archive(args):
    """(times, cells, per-cell fields) from the archive, for [start 01:00, end + 1 day 00:00]."""
    if not (args.start and args.end):
        raise SystemExit("ERROR: --data-path needs --start and --end (YYYY-MM-DD, inclusive)")
    first, last = common.parse_dates(args.start, args.end)
    t0 = dt.datetime(first.year, first.month, first.day, 1)
    t1 = dt.datetime(last.year, last.month, last.day) + dt.timedelta(days=1)
    cells = select_archive_cells(args.data_path, args)
    names = ["Tair", "Tdew", "PSurf", "u10", "v10", "Rainf", "SWdown", "LWdown"]
    series = {n: [[] for _ in cells] for n in names}
    times = []
    y, m = first.year, first.month
    while (y, m) <= (last.year, last.month):          # a day's last record (next 00:00) is in its month's file
        for n in names:
            path = common.archive_file(args.data_path, n, y, m)
            if not os.path.exists(path):
                raise SystemExit(f"ERROR: the archive has no {os.path.basename(path)}")
            with Dataset(path) as ds:
                tv = ds["time"]
                t = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                keep = np.array([t0 <= s <= t1 for s in t])
                if n == names[0]:
                    times += [s for s, k in zip(t, keep) if k]
                for g, (_, _, _, r, c) in enumerate(cells):
                    series[n][g].append(np.ma.filled(ds[n][keep, r, c], np.nan).astype(float))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    fields = []
    for g in range(len(cells)):
        s = {n: np.concatenate(series[n][g]) for n in names}
        fields.append(dict(
            Tair=s["Tair"], PSurf=s["PSurf"], Qair=dewpoint_to_specific_humidity(s["Tdew"], s["PSurf"]),
            u10=s["u10"], v10=s["v10"], Wind=np.hypot(s["u10"], s["v10"]),
            Rainf=s["Rainf"], SWdown=s["SWdown"], LWdown=s["LWdown"]))
    grid = [(la, lo, el) for la, lo, el, _, _ in cells]
    with Dataset(common.archive_file(args.data_path, "Tair", first.year, first.month)) as ds:
        source = f"{getattr(ds, 'product', 'ERA5-Land hourly')}, via the ED_ERA5land archive ({ds.source})"
    return times, grid, fields, source


def targets(args):
    """The requested (lat, lon) locations: --cells, else --lat/--lon."""
    if args.cells:
        return [tuple(float(x) for x in pair.split(",")) for pair in args.cells]
    return [(args.lat, args.lon)]


def write_meds_forcing(path, time_seconds, base_iso, grid_lat, grid_lon, grid_elev, fields, attrs):
    """Write the MEDS multi-grid forcing NetCDF. `fields` maps ALMA var name -> (ntime, ngrid) array,
    plus per-var (units, long_name, cell_methods)."""
    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", None)            # unlimited
    ds.createDimension("grid", len(grid_lat))

    tv = ds.createVariable("time", "f8", ("time",))
    tv.units = f"seconds since {base_iso}"
    tv.calendar = "proleptic_gregorian"
    tv.standard_name = "time"
    tv[:] = time_seconds

    lav = ds.createVariable("latitude", "f8", ("grid",))
    lav.units = "degrees_north"; lav.standard_name = "latitude"; lav[:] = grid_lat
    lov = ds.createVariable("longitude", "f8", ("grid",))
    lov.units = "degrees_east"; lov.standard_name = "longitude"; lov[:] = grid_lon
    elv = ds.createVariable("elevation", "f8", ("grid",))
    elv.units = "m"; elv.long_name = "surface elevation"; elv[:] = grid_elev

    for name, (arr, units, long_name, cell_methods) in fields.items():
        v = ds.createVariable(name, "f4", ("time", "grid"), fill_value=np.float32(1.0e20))
        v.units = units
        v.long_name = long_name
        v.cell_methods = cell_methods
        v[:, :] = arr

    for k, val in attrs.items():
        setattr(ds, k, val)
    ds.close()


def main(argv=None):
    ap = argparse.ArgumentParser(description="An ED_ERA5land archive or ERA5-Land box files -> a MEDS single "
                                             "forcing file ((time, grid) NetCDF).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--data-path", help="an ED_ERA5land archive folder (build_era5land_archive.py)")
    src.add_argument("--box-dir", help="a folder of box files from postprocess_era5land.py --split none")
    ap.add_argument("--start", help="archive only: first day, YYYY-MM-DD (UTC); the file starts at 01:00")
    ap.add_argument("--end", help="archive only: last day, YYYY-MM-DD (UTC, inclusive); the file ends at 00:00 "
                                  "the next day")
    ap.add_argument("--out", required=True, help="output MEDS forcing NetCDF")
    ap.add_argument("--lat", type=float, default=42.44, help="site latitude [deg N] (default Ithaca NY)")
    ap.add_argument("--lon", type=float, default=-76.50, help="site longitude [deg E] (default Ithaca NY)")
    ap.add_argument("--cells", nargs="+", default=None,
                    help='explicit locations "lat,lon" -> one grid index each (multi-grid file)')
    ap.add_argument("--all-cells", action="store_true", help="box files only: every cell of the box as a grid point")
    ap.add_argument("--max-distance-km", type=float, default=15.0,
                    help="archive only: the farthest a site may be from its cell (default 15)")
    ap.add_argument("--elevation", type=float, default=320.0,
                    help="box files only: the elevation [m] to record for each grid point (default 320, "
                         "Ithaca); from the archive it is the cell's orography")
    args = ap.parse_args(argv)
    if args.box_dir and (args.start or args.end):
        ap.error("--start/--end apply to --data-path; box files carry their own period")
    if args.data_path and args.all_cells:
        ap.error("--all-cells applies to --box-dir; name the locations with --cells for an archive")

    if args.data_path:
        times, grid, per_cell, source = from_archive(args)
    else:
        times, grid, per_cell, source = from_box_files(args)
    ng, nt = len(grid), len(times)
    names = ["Tair", "Qair", "PSurf", "u10", "v10", "Wind", "Rainf", "SWdown", "LWdown"]
    arrays = {n: np.column_stack([f[n] for f in per_cell]) for n in names}

    # ---- MEDS never gap-fills (design comment 2). Two steps: -----------------------------------
    # (a) TRIM the leading de-accumulation boundary: the single leading non-01Z sample cannot be
    #     de-accumulated (deaccumulate_hourly -> NaN). That is a structural boundary, not a data gap,
    #     so we DROP those rows (never write them) rather than fill them.
    flux_nan = (np.isnan(arrays["SWdown"]).any(axis=1) | np.isnan(arrays["Rainf"]).any(axis=1)
                | np.isnan(arrays["LWdown"]).any(axis=1))
    start = 0
    if args.box_dir:
        while start < nt and flux_nan[start]:
            start += 1
        if start > 0:
            print(f"note: trimming {start} leading de-accumulation boundary step(s) (not written)")
    times = times[start:]
    arrays = {n: a[start:] for n, a in arrays.items()}
    base_iso = times[0].strftime("%Y-%m-%d %H:%M:%S")
    time_seconds = np.array([(t - times[0]).total_seconds() for t in times], dtype=float)
    # (b) ERROR on any remaining missing value (a real data gap, or a source NaN). MEDS does NOT
    #     fill it -- gap-filling is the user's job, entirely upstream/outside MEDS.
    for name in names:
        arr = arrays[name]
        if np.isnan(arr).any():
            k = np.argwhere(np.isnan(arr))[0]
            raise SystemExit(f"ERROR: {name} has {int(np.isnan(arr).sum())} missing value(s) "
                             f"(first at time index {k[0]}, grid {k[1]}). MEDS does not gap-fill — "
                             f"fix the source upstream (design MEDS_FORCING_DESIGN.md §5.5).")

    meta = {
        "Tair":   ("K",         "air temperature",              "time: point"),
        "Qair":   ("kg kg-1",   "specific humidity",            "time: point"),
        "PSurf":  ("Pa",        "surface pressure",             "time: point"),
        "u10":    ("m s-1",     "eastward wind at 10 m",        "time: point"),
        "v10":    ("m s-1",     "northward wind at 10 m",       "time: point"),
        "Wind":   ("m s-1",     "wind speed at 10 m",           "time: point"),
        "Rainf":  ("kg m-2 s-1","precipitation rate",           "time: mean"),
        "SWdown": ("W m-2",     "downward shortwave (total)",   "time: mean"),
        "LWdown": ("W m-2",     "downward longwave",            "time: mean"),
    }
    fields = {n: (arrays[n],) + meta[n] for n in meta}
    attrs = dict(
        Conventions="MEDS-forcing-1.0",
        title="MEDS meteorological forcing",
        source=source,
        history="make_forcing_file.py",
        timestep_seconds=int(SEC_PER_HOUR),
        avg_convention="end",          # flux vars: mean over the hour ENDING at the stamp
        sw_input_kind="total",         # total SWdown; the Fortran reader partitions (design §5.6)
        time_zone="UTC",
        wind_meas_height_m=10.0,       # ERA5-Land wind is at 10 m (design §5.2/§10)
    )
    write_meds_forcing(args.out, time_seconds, base_iso, [g[0] for g in grid], [g[1] for g in grid],
                       [g[2] for g in grid], fields, attrs)
    print(f"wrote {args.out}: ntime={len(time_seconds)} ({times[0]:%Y-%m-%d %H:%M} .. {times[-1]:%Y-%m-%d %H:%M} UTC), "
          f"ngrid={ng}, cells={[(round(g[0], 3), round(g[1], 3)) for g in grid]}")


if __name__ == "__main__":
    main()
