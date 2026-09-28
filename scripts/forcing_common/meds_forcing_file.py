# SPDX-License-Identifier: Apache-2.0
"""meds_forcing_file.py -- what every MEDS forcing-preparation tool shares: the writer of an
ED_default forcing file ([forcing].format = "ED_default"), and a Python copy of the model's own
conversions for the tools that need them.

The copies mirror the Fortran term for term, so a value a tool computes is the value the model
would compute from the same inputs:

  sat_vapor_pressure, rh_to_specific_humidity,      meds_therm_lib, meds_forcing_kernels
  dewpoint_to_specific_humidity,
  specific_humidity_to_vapor_pressure
  pressure_at_height                                meds_lapse_rate%lapse_pressure (isothermal form)
  equation_of_time, solar_declination,              meds_forcing_kernels, meds_time
  solar_cosz, window_mean_cosz
  clearness_index, clear_sky_emissivity,            meds_forcing_kernels
  synthesize_lwdown

A forcing file stores the humidity its source measured (RHair, Tdew or Qair) and the model converts
it (docs/science/forcing.md sec. 7), so a tool needs the saturation curve only for work of its own,
such as synthesizing longwave or checking a provider's VPD. The file is always on a UTC clock.

Dependencies: numpy + netCDF4.
"""
import datetime as dt

import numpy as np
from netCDF4 import Dataset

CONVENTIONS = "MEDS-forcing-1.1"

# ---------------------------------------------------------------------------------------------
# Constants, as the model defines them (meds_constants, meds_forcing_kernels).
# ---------------------------------------------------------------------------------------------
GRAV = 9.80665             # [m s-2]
R_DRY = 287.04             # [J kg-1 K-1]
CP_AIR = 1004.6            # [J kg-1 K-1]
STEFAN = 5.670374419e-8    # [W m-2 K-4]
SOLAR_CONSTANT = 1361.0    # [W m-2]
COSZ_MIN = 1.0e-3          # below this cos z it is night
COSZ_BAR_MIN = 1.0e-3      # a window whose mean cos z is at most this is a night window
N_COSZ_SUB = 10            # sub-samples of the window mean
LW_CLOUD_A = 0.22          # [-] the model's default cloud coefficient ([forcing].lw_cloud_a)

# The file's variables and their CF attributes. Every forcing variable is float32 on (time, grid)
# with _FillValue 1e20; states are values at the stamps, fluxes are means over the interval the
# stamp convention names.
VARIABLE_ATTRIBUTES = {
    "Tair":   dict(units="K", long_name="air temperature", standard_name="air_temperature",
                   cell_methods="time: point"),
    "RHair":  dict(units="1", long_name="relative humidity over liquid water",
                   standard_name="relative_humidity", cell_methods="time: point"),
    "Tdew":   dict(units="K", long_name="dew point temperature", standard_name="dew_point_temperature",
                   cell_methods="time: point"),
    "Qair":   dict(units="kg kg-1", long_name="specific humidity", standard_name="specific_humidity",
                   cell_methods="time: point"),
    "PSurf":  dict(units="Pa", long_name="surface pressure", standard_name="surface_air_pressure",
                   cell_methods="time: point"),
    "Wind":   dict(units="m s-1", long_name="wind speed", standard_name="wind_speed",
                   cell_methods="time: point"),
    "u10":    dict(units="m s-1", long_name="eastward wind", standard_name="eastward_wind",
                   cell_methods="time: point"),
    "v10":    dict(units="m s-1", long_name="northward wind", standard_name="northward_wind",
                   cell_methods="time: point"),
    "Rainf":  dict(units="kg m-2 s-1", long_name="total precipitation rate",
                   standard_name="precipitation_flux", cell_methods="time: mean"),
    "SWdown": dict(units="W m-2", long_name="downward shortwave radiation (total)",
                   standard_name="surface_downwelling_shortwave_flux_in_air", cell_methods="time: mean"),
    "LWdown": dict(units="W m-2", long_name="downward longwave radiation",
                   standard_name="surface_downwelling_longwave_flux_in_air", cell_methods="time: mean"),
}
HUMIDITY_VARIABLES = ("RHair", "Tdew", "Qair")


# ---------------------------------------------------------------------------------------------
# Humidity and pressure (meds_therm_lib, meds_forcing_kernels, meds_lapse_rate).
# ---------------------------------------------------------------------------------------------
def sat_vapor_pressure(t_k):
    """Saturation vapour pressure [Pa] over liquid water at t_k [K]: Bolton (1980)."""
    tc = np.asarray(t_k, dtype=float) - 273.15
    return 611.2 * np.exp(17.67 * tc / (tc + 243.5))


def rh_to_specific_humidity(rh, t_k, p_pa):
    """q [kg/kg] from relative humidity (a fraction, clipped to [0, 1]), temperature and pressure."""
    e = np.clip(rh, 0.0, 1.0) * sat_vapor_pressure(t_k)
    return 0.622 * e / (np.asarray(p_pa, dtype=float) - 0.378 * e)


def dewpoint_to_specific_humidity(td_k, p_pa):
    """q [kg/kg] from dewpoint and pressure: the actual vapour pressure is e_sat(Td)."""
    e = sat_vapor_pressure(td_k)
    return 0.622 * e / (np.asarray(p_pa, dtype=float) - 0.378 * e)


def specific_humidity_to_vapor_pressure(q, p_pa):
    """e [Pa] from q and pressure: the exact inverse of the two conversions above."""
    q = np.asarray(q, dtype=float)
    return q * np.asarray(p_pa, dtype=float) / (0.622 + 0.378 * q)


def pressure_at_height(p_pa, t_k, dz):
    """Pressure dz [m] above (dz > 0) or below (dz < 0) a level at p_pa, t_k: the isothermal
    hypsometric form of meds_lapse_rate%lapse_pressure."""
    return np.asarray(p_pa, dtype=float) * np.exp(-GRAV * np.asarray(dz, dtype=float) / (R_DRY * np.asarray(t_k)))


# ---------------------------------------------------------------------------------------------
# Solar geometry on the UTC clock (meds_forcing_kernels, meds_time). `stamps` are numpy
# datetime64 values in UTC; the day of year is the UTC date's, as in the model.
# ---------------------------------------------------------------------------------------------
def _day_of_year_and_seconds(stamps):
    stamps = np.asarray(stamps, dtype="datetime64[s]")
    day = stamps.astype("datetime64[D]")
    year = stamps.astype("datetime64[Y]")
    doy = (day - year.astype("datetime64[D]")).astype(int) + 1
    sec = (stamps - day).astype("timedelta64[s]").astype(float)
    return doy, sec


def equation_of_time(doy):
    """[s] Spencer (1971), as meds_forcing_kernels%equation_of_time."""
    b = 2.0 * np.pi * (np.asarray(doy, dtype=float) - 1.0) / 365.0
    return 229.18 * (0.000075 + 0.001868 * np.cos(b) - 0.032077 * np.sin(b)
                     - 0.014615 * np.cos(2.0 * b) - 0.040849 * np.sin(2.0 * b)) * 60.0


def solar_declination(doy):
    """[rad] Cooper (1969), as meds_time%solar_declination."""
    return 23.45 * np.pi / 180.0 * np.sin(2.0 * np.pi * (284.0 + np.asarray(doy, dtype=float)) / 365.0)


def _cosz(doy, sec_utc, latitude, longitude):
    sec_solar = sec_utc + longitude * 240.0 + equation_of_time(doy)
    lat = np.radians(latitude)
    decl = solar_declination(doy)
    hour_angle = 2.0 * np.pi * (sec_solar / 86400.0 - 0.5)
    return np.maximum(np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(hour_angle), 0.0)


def solar_cosz(stamps, latitude, longitude):
    """cos(solar zenith) at UTC instants, floored at 0 (met_solar_cosz)."""
    doy, sec = _day_of_year_and_seconds(stamps)
    return _cosz(doy, sec, latitude, longitude)


def window_mean_cosz(window_start, window_seconds, latitude, longitude, nsub=N_COSZ_SUB):
    """The mean cos z over [start, start + window) by the model's midpoint rule, night samples 0,
    with every sample on the window start's day of year (cosz_reconstruct_factor)."""
    doy, sec0 = _day_of_year_and_seconds(window_start)
    total = np.zeros(np.shape(sec0))
    for i in range(1, nsub + 1):
        total = total + _cosz(doy, sec0 + (i - 0.5) * window_seconds / nsub, latitude, longitude)
    return total / nsub


def clearness_index(sw_total, cosz):
    """kt = SW / (S0 cos z), clipped to [0, 1]; -1 where the sun is down (the model's sentinel)."""
    sw_total = np.asarray(sw_total, dtype=float)
    cosz = np.asarray(cosz, dtype=float)
    kt = np.clip(np.maximum(sw_total, 0.0) / (SOLAR_CONSTANT * np.maximum(cosz, 1e-30)), 0.0, 1.0)
    return np.where(cosz <= COSZ_MIN, -1.0, kt)


# ---------------------------------------------------------------------------------------------
# Longwave synthesis (meds_forcing_kernels, #182).
# ---------------------------------------------------------------------------------------------
def clear_sky_emissivity(t_k, q, p_pa, form="brutsaert"):
    """Brutsaert (1975) from the vapour pressure of q, or Idso & Jackson (1969); bounded [0.5, 1]."""
    t_k = np.asarray(t_k, dtype=float)
    if form == "idso":
        eps = 1.0 - 0.261 * np.exp(-7.77e-4 * (273.16 - t_k) ** 2)
    else:
        e_hpa = 0.01 * specific_humidity_to_vapor_pressure(np.maximum(q, 0.0), p_pa)
        eps = 1.24 * (np.maximum(e_hpa, 1e-30) / t_k) ** (1.0 / 7.0)
    return np.clip(eps, 0.5, 1.0)


def synthesize_lwdown(t_k, q, p_pa, kt, cloud_a=LW_CLOUD_A, form="brutsaert"):
    """LW = eps_clear sigma T^4 [1 + a (1 - kt)], kt clipped to [0, 1]."""
    eps = clear_sky_emissivity(t_k, q, p_pa, form)
    return eps * STEFAN * np.asarray(t_k, dtype=float) ** 4 * (1.0 + max(cloud_a, 0.0) * (1.0 - np.clip(kt, 0.0, 1.0)))


# ---------------------------------------------------------------------------------------------
# The file.
# ---------------------------------------------------------------------------------------------
QC_FLAG_VALUES = np.array([0, 1, 2, 3, 4], dtype=np.int8)
QC_FLAG_MEANINGS = "observed short_gap_interpolation era5land_regression synthesis_or_mean_diurnal_variation filled_by_provider"


def check_complete(times, arrays):
    """MEDS never gap-fills: a missing value left in the file is an error, named here, not in a run."""
    for name, arr in arrays.items():
        missing = ~np.isfinite(arr)
        if missing.any():
            k, g = np.argwhere(missing)[0]
            raise SystemExit(f"ERROR: {name} has {int(missing.sum())} missing value(s), first at "
                             f"{np.datetime_as_string(times[k], unit='m')} UTC in grid {g}. "
                             f"MEDS does not gap-fill: fill it upstream.")


def write_forcing_file(path, times, grid, arrays, attrs, qc=None, variable_attributes=None,
                       elevation_long_name="elevation of the grid point"):
    """Write an ED_default forcing file.

    times   UTC stamps (numpy datetime64), uniform
    grid    [(latitude, longitude, elevation or None), ...], one per grid point
    arrays  {variable: (ntime, ngrid) float array}, names from VARIABLE_ATTRIBUTES
    attrs   global attributes; Conventions and time_zone = "UTC" are set here
    qc      optional {variable: (ntime, ngrid) int8 array} of QC_FLAG_* codes, written as <var>_qc
    variable_attributes  optional {variable: CF attributes} replacing VARIABLE_ATTRIBUTES' entry
    """
    times = np.asarray(times, dtype="datetime64[s]")
    humidity = [name for name in HUMIDITY_VARIABLES if name in arrays]
    if len(humidity) != 1:
        raise SystemExit(f"ERROR: a forcing file carries exactly one of {HUMIDITY_VARIABLES}, not {humidity}")
    check_complete(times, arrays)
    lat, lon, elevation = zip(*grid)
    t0 = times[0].astype(dt.datetime)
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", None)
        ds.createDimension("grid", len(grid))
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = f"seconds since {t0:%Y-%m-%d %H:%M:%S}"
        tv.calendar = "proleptic_gregorian"
        tv.standard_name = "time"
        tv[:] = (times - times[0]).astype("timedelta64[s]").astype(float)
        v = ds.createVariable("latitude", "f8", ("grid",))
        v.units, v.standard_name, v[:] = "degrees_north", "latitude", lat
        v = ds.createVariable("longitude", "f8", ("grid",))
        v.units, v.standard_name, v[:] = "degrees_east", "longitude", lon
        if None not in elevation:
            v = ds.createVariable("elevation", "f8", ("grid",))
            v.units, v.long_name, v[:] = "m", elevation_long_name, elevation
        for name, values in arrays.items():
            v = ds.createVariable(name, "f4", ("time", "grid"), fill_value=np.float32(1.0e20))
            v.setncatts((variable_attributes or {}).get(name) or VARIABLE_ATTRIBUTES[name])
            v[:, :] = values
        for name, flags in (qc or {}).items():
            v = ds.createVariable(f"{name}_qc", "i1", ("time", "grid"))
            v.long_name = f"how each {name} value was obtained"
            v.flag_values = QC_FLAG_VALUES
            v.flag_meanings = QC_FLAG_MEANINGS
            v[:, :] = flags
        ds.setncatts(dict(attrs, Conventions=CONVENTIONS, time_zone="UTC"))
