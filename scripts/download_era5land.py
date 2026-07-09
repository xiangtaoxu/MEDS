#!/usr/bin/env python3
# =========================================================================================
# download_era5land.py -- pull ERA5-Land hourly reanalysis for a grid cell (or box) from the
# Copernicus Climate Data Store (CDS), as the raw input to prep_era5land_forcing.py (which
# converts it into a MEDS multi-grid forcing NetCDF). Stage 1 of 2 (design MEDS_FORCING_DESIGN.md
# section 7).
#
# The eight variables MEDS forcing needs (all ERA5-Land has; CO2 is NOT in ERA5-Land and is
# supplied as a config constant downstream):
#   instantaneous : 2m_temperature, 2m_dewpoint_temperature, surface_pressure,
#                   10m_u_component_of_wind, 10m_v_component_of_wind
#   accumulated   : total_precipitation, surface_solar_radiation_downwards,
#                   surface_thermal_radiation_downwards
#
# Verified against the CDS docs (2026-07-08 research fact sheet):
#   * dataset id            : "reanalysis-era5-land" (hourly, 0.1 deg, UTC)
#   * ~/.cdsapirc (2 lines) : url: https://cds.climate.copernicus.eu/api   (NO /v2)
#                             key: <PERSONAL-ACCESS-TOKEN>                  (single token, no UID)
#   * client                : pip install "cdsapi>=0.7.7"  (Python >= 3.8)
#   * area                  : [North, West, South, East] decimal degrees (W/S may be negative)
#   * ZIP FALLBACK          : a NetCDF request mixing instantaneous + accumulated fields is
#                             returned as a .zip (one .nc per GRIB stepType) even with
#                             download_format="unarchived". This script detects + merges it.
#   * DE-ACCUM trailing day : the 00:00 UTC stamp holds the WHOLE previous day; to de-accumulate
#                             the last requested day's [23Z,00Z] hour, one extra trailing day is
#                             fetched (prep_era5land_forcing.py does the de-accumulation).
#
# You must accept the ERA5-Land licence once on the CDS website before the first API pull.
#
# Usage:
#   python download_era5land.py --start 2020-06-01 --end 2020-08-31 \
#          --lat 42.44 --lon -76.50 --out data/forcing/era5land_ithaca_2020JJA.nc
#   (defaults are the Ithaca NY reference cell)
# =========================================================================================
import argparse
import datetime as dt
import os
import sys
import zipfile

# The eight ERA5-Land variables MEDS needs, grouped by GRIB stepType (the split that triggers
# the CDS zip fallback). Requesting the two groups separately (--split-requests) sidesteps the
# zip entirely; the default single request handles the zip by extract+merge.
INSTANT_VARS = [
    "2m_temperature",
    "2m_dewpoint_temperature",
    "surface_pressure",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
]
ACCUM_VARS = [
    "total_precipitation",
    "surface_solar_radiation_downwards",
    "surface_thermal_radiation_downwards",
]
DATASET = "reanalysis-era5-land"


def daterange_days(start, end):
    """Inclusive list of datetime.date from start to end."""
    n = (end - start).days
    return [start + dt.timedelta(days=i) for i in range(n + 1)]


def ym_day_lists(days):
    """CDS requests are (year, month, day) lists. Collapse a day list into the minimal set of
    full (year, month) requests, each carrying the days that fall in that month. Returns a list
    of dicts {year, month, day:[...]}."""
    from collections import OrderedDict
    buckets = OrderedDict()
    for d in days:
        buckets.setdefault((d.year, d.month), []).append(f"{d.day:02d}")
    return [
        {"year": f"{y:04d}", "month": f"{m:02d}", "day": days_}
        for (y, m), days_ in buckets.items()
    ]


def bbox_from_point(lat, lon, half):
    """area = [North, West, South, East]. A small box (half-degree default) captures the covering
    0.1 deg cell plus neighbours (useful for a future nearest-match / multi-grid file)."""
    return [round(lat + half, 4), round(lon - half, 4),
            round(lat - half, 4), round(lon + half, 4)]


def build_request(variables, ym, area):
    return {
        "variable": variables,
        "year": ym["year"],
        "month": ym["month"],
        "day": ym["day"],
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": area,
        "data_format": "netcdf",
        "download_format": "unarchived",   # may STILL return .zip when stepTypes are mixed
    }


def is_zip(path):
    with open(path, "rb") as fh:
        return fh.read(4)[:2] == b"PK"


def unzip_and_merge(zip_path, out_path):
    """The CDS zip holds one .nc per GRIB stepType (instantaneous, accumulated). Extract and merge
    them on the shared time/lat/lon axes into a single NetCDF at out_path."""
    import xarray as xr
    extract_dir = out_path + ".unzipped"
    os.makedirs(extract_dir, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if m.endswith(".nc")]
        zf.extractall(extract_dir)
    ncs = [os.path.join(extract_dir, m) for m in members]
    dsets = [xr.open_dataset(p) for p in ncs]
    merged = xr.merge(dsets, compat="override", join="outer")
    merged.to_netcdf(out_path)
    for d in dsets:
        d.close()
    return out_path


def retrieve(client, dataset, request, target):
    """One CDS retrieval; returns the (possibly-zip) file path actually written."""
    client.retrieve(dataset, request, target)
    return target


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download ERA5-Land hourly forcing for a grid cell.")
    ap.add_argument("--start", required=True, help="start date YYYY-MM-DD (UTC)")
    ap.add_argument("--end", required=True, help="end date YYYY-MM-DD (UTC, inclusive)")
    ap.add_argument("--lat", type=float, default=42.44, help="site latitude [deg N] (default Ithaca NY)")
    ap.add_argument("--lon", type=float, default=-76.50, help="site longitude [deg E] (default Ithaca NY)")
    ap.add_argument("--half-box", type=float, default=0.10,
                    help="half-width of the lat/lon box [deg] around the point (default 0.10)")
    ap.add_argument("--out", default="era5land_forcing_raw.nc", help="output NetCDF path")
    ap.add_argument("--split-requests", action="store_true",
                    help="issue two single-stepType requests (instant/accum) to avoid the zip fallback")
    args = ap.parse_args(argv)

    try:
        import cdsapi
    except ImportError:
        sys.exit("cdsapi not installed. Run:  pip install 'cdsapi>=0.7.7'  and set up ~/.cdsapirc "
                 "(url: https://cds.climate.copernicus.eu/api ; key: <token>).")

    start = dt.date.fromisoformat(args.start)
    end = dt.date.fromisoformat(args.end)
    if end < start:
        sys.exit("--end precedes --start")
    # Pad by one trailing day so prep can de-accumulate the last day's [23Z,00Z] hour (00Z next day).
    end_padded = end + dt.timedelta(days=1)
    days = daterange_days(start, end_padded)
    area = bbox_from_point(args.lat, args.lon, args.half_box)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    client = cdsapi.Client()

    print(f"ERA5-Land {DATASET}: {start}..{end} (+1 trailing day), area(N,W,S,E)={area}")

    if args.split_requests:
        # Two clean single-stepType requests, merged locally -> no zip.
        import xarray as xr
        parts = []
        for tag, vlist in (("instant", INSTANT_VARS), ("accum", ACCUM_VARS)):
            part_dsets = []
            for ym in ym_day_lists(days):
                tmp = f"{args.out}.{tag}.{ym['year']}{ym['month']}.nc"
                retrieve(client, DATASET, build_request(vlist, ym, area), tmp)
                if is_zip(tmp):
                    tmp = unzip_and_merge(tmp, tmp + ".merged.nc")
                part_dsets.append(xr.open_dataset(tmp))
            parts.append(xr.concat(part_dsets, dim=_time_dim(part_dsets[0])) if len(part_dsets) > 1
                         else part_dsets[0])
        xr.merge(parts, compat="override", join="outer").to_netcdf(args.out)
    else:
        # Single request (all 8 vars). Concatenate per-month pulls; handle the zip fallback.
        import xarray as xr
        month_dsets = []
        for ym in ym_day_lists(days):
            tmp = f"{args.out}.{ym['year']}{ym['month']}.nc"
            retrieve(client, DATASET, build_request(INSTANT_VARS + ACCUM_VARS, ym, area), tmp)
            if is_zip(tmp):
                tmp = unzip_and_merge(tmp, tmp + ".merged.nc")
            month_dsets.append(xr.open_dataset(tmp))
        ds = (xr.concat(month_dsets, dim=_time_dim(month_dsets[0])) if len(month_dsets) > 1
              else month_dsets[0])
        ds.to_netcdf(args.out)

    print(f"wrote {args.out}")


def _time_dim(ds):
    """The CDS NetCDF time coordinate is 'valid_time' (new engine) or 'time' (legacy)."""
    return "valid_time" if "valid_time" in ds.dims else "time"


if __name__ == "__main__":
    main()
