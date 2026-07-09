#!/usr/bin/env python3
# =========================================================================================
# prep_era5land_forcing.py -- convert a raw ERA5-Land hourly NetCDF (from download_era5land.py)
# into a MEDS multi-grid forcing NetCDF the Fortran reader (src/forcing/meds_met_driver.f90)
# consumes. Stage 2 of 2 (design MEDS_FORCING_DESIGN.md sections 5.2, 7).
#
# What it does (all verified against the CDS docs, 2026-07-08 research fact sheet):
#   1. Read the ERA5-Land vars by their NetCDF names: t2m, d2m, sp, u10, v10, tp, ssrd, strd,
#      on (valid_time|time, latitude, longitude), UTC.
#   2. Select one or more grid cells (nearest to requested lat/lon, or the whole box), stack into
#      an unstructured `grid` dimension  ->  MEDS `(time, grid)` layout (multi-grid by construction).
#   3. De-accumulate tp/ssrd/strd. ERA5-Land accumulates from 00 UTC and resets daily; the 00:00
#      stamp is the WHOLE PREVIOUS DAY total. Correct rule (fact sheet C): per-hour(H)=raw(H)-raw(H-1)
#      everywhere EXCEPT at 01:00 UTC where per-hour=raw(01:00) as-is. Clip GRIB-packing negatives.
#   4. Convert units: precip [m/hr] -> [kg/m2/s] (x1000/3600); ssrd/strd [J/m2/hr] -> [W/m2] (/3600).
#   5. Specific humidity from dewpoint + surface pressure, using the SAME saturation form as
#      meds_thermo (Bolton-1980) so the file's Qair matches the reader's dewpoint_to_specific_humidity.
#   6. Wind speed = sqrt(u10^2 + v10^2)  (10 m; height recorded in attrs for a future adjustment).
#   7. CO2: ERA5-Land has none -> a config constant (default 420 umol/mol), attribute-flagged synthetic.
#   8. Write the MEDS forcing NetCDF: dims (time, grid); coords time/latitude/longitude/elevation;
#      total SWdown (reader partitions, §5.6); ALMA-named vars; global provenance attrs.
#
# Dependencies: numpy + netCDF4 (dependency-light; no xarray required).
#
# Usage:
#   python prep_era5land_forcing.py --in era5land_ithaca.nc --out ithaca_forcing.nc \
#          --lat 42.44 --lon -76.50            # single nearest cell (default Ithaca NY)
#   python prep_era5land_forcing.py --in raw.nc --out multi.nc --all-cells   # every cell -> grid
#   python prep_era5land_forcing.py --in raw.nc --out multi.nc \
#          --cells 42.44,-76.50 44.00,-77.00   # explicit list of locations -> grid
# =========================================================================================
import argparse
import sys

import numpy as np
from netCDF4 import Dataset, num2date

RHO_W = 1000.0          # [kg/m3] water density (1 m of water = 1000 kg/m2)
SEC_PER_HOUR = 3600.0
NEG_EPS = 1.0e-5        # [m] clip threshold for GRIB-packing negatives (precip); scaled for fluxes
U_MIN = 0.1             # [m/s] wind floor (M-O similarity stability in the aero kernel)
CO2_DEFAULT = 420.0     # [umol/mol] present-day free-atmosphere CO2 (ERA5-Land has none)

# ERA5-Land NetCDF variable names (NOT the GRIB shortnames).
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
# ---------------------------------------------------------------------------------------------
def deaccumulate_hourly(accum, hours_utc, clip_eps=NEG_EPS):
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
    # Clip GRIB simple-packing negatives to 0 (ECMWF-recommended); keep NaN as NaN.
    neg = np.isfinite(out) & (out <= clip_eps)
    out[neg] = 0.0
    return out


# ---------------------------------------------------------------------------------------------
# Read the raw ERA5-Land file. Returns times (list[datetime], UTC), 1-D lat/lon arrays, and a dict
# of (time, lat, lon) arrays for each needed variable. Handles 'valid_time'|'time' and stray
# singleton dims (e.g. 'number'/'expver' the new CDS sometimes adds).
# ---------------------------------------------------------------------------------------------
def read_era5land(path):
    ds = Dataset(path)
    tname = "valid_time" if "valid_time" in ds.variables else "time"
    tvar = ds.variables[tname]
    times = num2date(tvar[:], tvar.units,
                     getattr(tvar, "calendar", "standard"),
                     only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    times = list(np.atleast_1d(times))
    latname = "latitude" if "latitude" in ds.variables else "lat"
    lonname = "longitude" if "longitude" in ds.variables else "lon"
    lat = np.atleast_1d(ds.variables[latname][:]).astype(float)
    lon = np.atleast_1d(ds.variables[lonname][:]).astype(float)

    def get(v):
        raw = ds.variables[v][:]
        arr = np.ma.filled(raw, np.nan).astype(float)
        # collapse to (time, lat, lon): drop any leading singleton dims (number/expver)
        while arr.ndim > 3:
            arr = arr[0]
        if arr.ndim == 2:                       # (lat, lon): a single time step
            arr = arr[np.newaxis, :, :]
        return arr

    data = {k: get(nc) for k, nc in NC_NAMES.items()}
    ds.close()
    return times, lat, lon, data


def nearest_index(lat, lon, tlat, tlon):
    """Nearest (ilat, ilon) grid index to a target lat/lon (planar distance; cells are ~0.1 deg)."""
    ilat = int(np.argmin(np.abs(lat - tlat)))
    ilon = int(np.argmin(np.abs(lon - tlon)))
    return ilat, ilon


def select_cells(lat, lon, args):
    """Return a list of (name_lat, name_lon, ilat, ilon) selected grid points -> the `grid` dim."""
    cells = []
    if args.all_cells:
        for i, la in enumerate(lat):
            for j, lo in enumerate(lon):
                cells.append((float(la), float(lo), i, j))
    elif args.cells:
        for pair in args.cells:
            tla, tlo = (float(x) for x in pair.split(","))
            i, j = nearest_index(lat, lon, tla, tlo)
            cells.append((float(lat[i]), float(lon[j]), i, j))
    else:
        i, j = nearest_index(lat, lon, args.lat, args.lon)
        cells.append((float(lat[i]), float(lon[j]), i, j))
    return cells


def write_meds_forcing(path, time_seconds, base_iso, grid_lat, grid_lon, grid_elev, fields, attrs):
    """Write the MEDS multi-grid forcing NetCDF. `fields` maps ALMA var name -> (ntime, ngrid) array,
    plus per-var (units, long_name, cell_methods)."""
    ng = len(grid_lat)
    nt = len(time_seconds)
    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", None)            # unlimited
    ds.createDimension("grid", ng)

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
    ap = argparse.ArgumentParser(description="ERA5-Land raw NetCDF -> MEDS multi-grid forcing NetCDF.")
    ap.add_argument("--in", dest="inp", required=True, help="raw ERA5-Land NetCDF (download_era5land.py)")
    ap.add_argument("--out", required=True, help="output MEDS forcing NetCDF")
    ap.add_argument("--lat", type=float, default=42.44, help="site latitude [deg N] (default Ithaca NY)")
    ap.add_argument("--lon", type=float, default=-76.50, help="site longitude [deg E] (default Ithaca NY)")
    ap.add_argument("--cells", nargs="+", default=None,
                    help='explicit locations "lat,lon" -> one grid index each (multi-grid file)')
    ap.add_argument("--all-cells", action="store_true", help="keep every cell of the box as a grid point")
    ap.add_argument("--elevation", type=float, default=320.0, help="site elevation [m] (Ithaca ~320)")
    ap.add_argument("--co2", type=float, default=CO2_DEFAULT, help="constant CO2 [umol/mol] (ERA5-Land has none)")
    ap.add_argument("--source", default="ERA5-Land hourly (reanalysis-era5-land)", help="provenance string")
    args = ap.parse_args(argv)

    times, lat, lon, data = read_era5land(args.inp)
    hours_utc = np.array([t.hour for t in times], dtype=int)
    cells = select_cells(lat, lon, args)
    ng = len(cells)
    nt = len(times)

    # time axis: seconds since the first timestamp
    base = times[0]
    base_iso = base.strftime("%Y-%m-%d %H:%M:%S")
    time_seconds = np.array([(t - base).total_seconds() for t in times], dtype=float)

    # allocate (time, grid)
    Tair = np.empty((nt, ng)); Qair = np.empty((nt, ng)); PSurf = np.empty((nt, ng))
    Wind = np.empty((nt, ng)); Rainf = np.empty((nt, ng)); SWdown = np.empty((nt, ng))
    LWdown = np.empty((nt, ng))
    grid_lat = np.empty(ng); grid_lon = np.empty(ng); grid_elev = np.full(ng, args.elevation)

    for g, (cla, clo, i, j) in enumerate(cells):
        grid_lat[g] = cla; grid_lon[g] = clo
        t2m = data["t2m"][:, i, j]; d2m = data["d2m"][:, i, j]; sp = data["sp"][:, i, j]
        u10 = data["u10"][:, i, j]; v10 = data["v10"][:, i, j]
        Tair[:, g] = t2m
        PSurf[:, g] = sp
        Qair[:, g] = dewpoint_to_specific_humidity(d2m, sp)
        Wind[:, g] = np.maximum(np.hypot(u10, v10), U_MIN)
        # de-accumulate fluxes (per-hour amount), then convert to rate / mean flux
        tp_hr = deaccumulate_hourly(data["tp"][:, i, j], hours_utc, clip_eps=NEG_EPS)
        ssrd_hr = deaccumulate_hourly(data["ssrd"][:, i, j], hours_utc, clip_eps=NEG_EPS * 1.0e6)
        strd_hr = deaccumulate_hourly(data["strd"][:, i, j], hours_utc, clip_eps=-1.0e30)  # LW: don't clip
        Rainf[:, g] = tp_hr * RHO_W / SEC_PER_HOUR       # [m/hr] -> [kg/m2/s]
        SWdown[:, g] = ssrd_hr / SEC_PER_HOUR            # [J/m2/hr] -> [W/m2] (total; reader partitions)
        LWdown[:, g] = strd_hr / SEC_PER_HOUR            # [J/m2/hr] -> [W/m2]

    CO2air = np.full((nt, ng), args.co2, dtype=float)

    # ---- MEDS never gap-fills (design comment 2). Two steps: -----------------------------------
    # (a) TRIM the leading de-accumulation boundary: the single leading non-01Z sample cannot be
    #     de-accumulated (deaccumulate_hourly -> NaN). That is a structural boundary, not a data gap,
    #     so we DROP those rows (never write them) rather than fill them.
    time_vars = [Tair, Qair, PSurf, Wind, Rainf, SWdown, LWdown, CO2air]
    flux_nan = (np.isnan(SWdown).any(axis=1) | np.isnan(Rainf).any(axis=1) | np.isnan(LWdown).any(axis=1))
    start = 0
    while start < nt and flux_nan[start]:
        start += 1
    if start > 0:
        print(f"note: trimming {start} leading de-accumulation boundary step(s) (not written)")
    sl = slice(start, nt)
    time_seconds = time_seconds[sl]
    Tair, Qair, PSurf, Wind, Rainf, SWdown, LWdown, CO2air = (v[sl] for v in time_vars)
    if start > 0:                                    # rebase the time axis to the new first sample
        base = times[start]
        base_iso = base.strftime("%Y-%m-%d %H:%M:%S")
        time_seconds = time_seconds - time_seconds[0]
    # (b) ERROR on any remaining missing value (a real data gap, or a source NaN). MEDS does NOT
    #     fill it -- gap-filling is the user's job, entirely upstream/outside MEDS.
    for name, arr in [("Tair", Tair), ("Qair", Qair), ("PSurf", PSurf), ("Wind", Wind),
                      ("Rainf", Rainf), ("SWdown", SWdown), ("LWdown", LWdown)]:
        if np.isnan(arr).any():
            k = np.argwhere(np.isnan(arr))[0]
            raise SystemExit(f"ERROR: {name} has {int(np.isnan(arr).sum())} missing value(s) "
                             f"(first at time index {k[0]}, grid {k[1]}). MEDS does not gap-fill — "
                             f"fix the source upstream (design MEDS_FORCING_DESIGN.md §5.5).")

    fields = {
        "Tair":   (Tair,   "K",         "air temperature",              "time: point"),
        "Qair":   (Qair,   "kg kg-1",   "specific humidity",            "time: point"),
        "PSurf":  (PSurf,  "Pa",        "surface pressure",             "time: point"),
        "Wind":   (Wind,   "m s-1",     "wind speed",                   "time: point"),
        "Rainf":  (Rainf,  "kg m-2 s-1","precipitation rate",           "time: mean"),
        "SWdown": (SWdown, "W m-2",     "downward shortwave (total)",   "time: mean"),
        "LWdown": (LWdown, "W m-2",     "downward longwave",            "time: mean"),
        "CO2air": (CO2air, "umol mol-1","free-atmosphere CO2 (synthetic constant)", "time: point"),
    }
    attrs = dict(
        Conventions="MEDS-forcing-1.0",
        title="MEDS meteorological forcing",
        source=args.source,
        history="prep_era5land_forcing.py",
        timestep_seconds=int(SEC_PER_HOUR),
        avg_convention="end",          # flux vars: mean over the hour ENDING at the stamp
        sw_input_kind="total",         # total SWdown; the Fortran reader partitions (design §5.6)
        time_zone="UTC",
        co2_const=float(args.co2),
        co2_note="ERA5-Land has no CO2; a single constant was applied",
        wind_meas_height_m=10.0,       # ERA5-Land wind is at 10 m (design §5.2/§10)
    )
    write_meds_forcing(args.out, time_seconds, base_iso, grid_lat, grid_lon, grid_elev, fields, attrs)
    print(f"wrote {args.out}: ntime={len(time_seconds)}, ngrid={ng}, "
          f"cells={[(round(a,3), round(b,3)) for a, b, _, _ in cells]}")


if __name__ == "__main__":
    main()
