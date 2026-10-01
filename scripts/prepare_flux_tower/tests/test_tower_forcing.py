# SPDX-License-Identifier: Apache-2.0
"""Tests of the flux-tower forcing tool (MEDS_FLUX_TOWER_FORCING_PLAN.md sec. 11).

Every input is synthetic, made here from a known sun and a known humidity, so each check has an
answer that does not come from the code under test. Run with:
    python -m pytest scripts/prepare_flux_tower/tests
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HERE, "..", "..", "prepare_era5"))
import make_tower_forcing as mt   # noqa: E402
import tower_gapfill as tg        # noqa: E402
import tower_inputs as ti         # noqa: E402

conv = ti.conv
LAT, LON, OFFSET = 9.15, -79.85, -5.0
DT = 1800.0
DAYS = 45


def synthetic_tower(start="2014-02-01 00:00", days=DAYS):
    """A begin-stamped half-hourly tower on UTC-5: clear-sky shortwave from the model's own sun,
    a diurnal temperature, RH, and the VPD a provider would compute from them with the
    Alduchov-Eskridge curve."""
    local = pd.date_range(start, periods=int(days * 86400 / DT), freq="30min")
    utc = local.values.astype("datetime64[s]") - np.timedelta64(int(OFFSET * 3600), "s")
    cz = conv.window_mean_cosz(utc, DT, LAT, LON)
    hour = local.hour + local.minute / 60.0
    tc = 25.0 + 3.0 * np.sin(2 * np.pi * (hour - 9.0) / 24.0)
    rh = 0.85 - 0.12 * np.sin(2 * np.pi * (hour - 9.0) / 24.0)
    vpd_kpa = (1.0 - rh) * ti.SATURATION_CURVES["alduchov_eskridge"](tc) / 1000.0
    rain = np.where((local.day % 5 == 0) & (local.hour == 15), 6.0, 0.0)      # mm per half hour
    lw = 440.0 + 2.0 * (tc - 25.0)
    return pd.DataFrame({"date": local.strftime("%Y-%m-%d %H:%M"), "tair": tc, "RH": 100.0 * rh,
                         "vpd": vpd_kpa, "p_kpa": 98.8, "PPT": rain, "Rs": 1000.0 * cz, "Rl_dn": lw,
                         "ubar": 2.5})


def write_site(tmp_path, frame, stamp="begin", utc_offset=OFFSET, curve="alduchov_eskridge",
               extra=""):
    csv = tmp_path / "tower.csv"
    frame.to_csv(csv, index=False)
    toml = tmp_path / "site.toml"
    toml.write_text(f'''
[input]
format = "csv"
path = "tower.csv"
timestamp = "date"
timestamp_format = "%Y-%m-%d %H:%M"
missing = ["NaN"]
[site]
name = "synthetic"
latitude = {LAT}
longitude = {LON}
elevation = 100.0
[clock]
utc_offset = {utc_offset}
stamp = "{stamp}"
timestep = 1800
[heights]
tq_height = 30.0
wind_height = 32.0
pressure_height = 30.0
[variables]
Tair = {{ column = "tair", units = "degC" }}
RH = {{ column = "RH", units = "%" }}
VPD = {{ column = "vpd", units = "kPa", curve = "{curve}" }}
PSurf = {{ column = "p_kpa", units = "kPa" }}
Rainf = {{ column = "PPT", units = "mm" }}
SWdown = {{ column = "Rs", units = "W m-2" }}
LWdown = {{ column = "Rl_dn", units = "W m-2" }}
Wind = {{ column = "ubar", units = "m s-1" }}
[gapfill]
{extra}
''')
    return ti.read_site(str(toml))


def build(tmp_path, site, **kw):
    out = str(tmp_path / "forcing.nc")
    report = mt.build(site, out, **kw)
    return Dataset(out), report


# ---------------------------------------------------------------------------------------------
# The declarations are checked against the data.
# ---------------------------------------------------------------------------------------------
def test_true_declaration_passes_the_sun_and_humidity_checks(tmp_path):
    ds, report = build(tmp_path, write_site(tmp_path, synthetic_tower()))
    assert abs(report["V2_sun"]["best_shift_min"]) <= 5.0
    assert report["V2_sun"]["night_shortwave_fraction"] < 1e-3
    assert report["V3_humidity"]["residual_p99_pa"] < 1e-6
    assert report["V3_humidity"]["best_fit"] == "alduchov_eskridge"


def test_wrong_stamp_convention_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="V2"):
        build(tmp_path, write_site(tmp_path, synthetic_tower(), stamp="end"))


def test_local_clock_declared_utc_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="V2"):
        build(tmp_path, write_site(tmp_path, synthetic_tower(), utc_offset=0.0))


def test_wrong_vpd_curve_is_refused_and_names_the_right_one(tmp_path):
    with pytest.raises(SystemExit, match="alduchov_eskridge"):
        build(tmp_path, write_site(tmp_path, synthetic_tower(), curve="bolton"))


def test_a_ragged_axis_is_refused(tmp_path):
    frame = synthetic_tower().drop(index=[100])
    with pytest.raises(SystemExit, match="V1"):
        build(tmp_path, write_site(tmp_path, frame))


def test_percent_units_declared_wrong_are_refused(tmp_path):
    frame = synthetic_tower()
    frame["p_kpa"] = frame["p_kpa"] * 1000.0          # Pa, but the TOML says kPa
    with pytest.raises(SystemExit, match="V4"):
        build(tmp_path, write_site(tmp_path, frame))


# ---------------------------------------------------------------------------------------------
# What the file carries.
# ---------------------------------------------------------------------------------------------
def test_the_file_is_utc_begin_stamped_with_its_heights(tmp_path):
    ds, _ = build(tmp_path, write_site(tmp_path, synthetic_tower()))
    assert ds.time_zone == "UTC" and ds.avg_convention == "begin" and ds.Conventions == "MEDS-forcing-1.1"
    assert ds.tq_height_m == 30.0 and ds.wind_height_m == 32.0 and ds.height_above == "ground"
    assert ds["time"].units == "seconds since 2014-02-01 05:00:00"   # 00:00 UTC-5
    assert "RHair" in ds.variables and "Qair" not in ds.variables and "Tdew" not in ds.variables


def test_relative_humidity_is_the_measured_one(tmp_path):
    frame = synthetic_tower()
    frame["RH"] = 72.0                                 # constant, so re-centring cannot move it
    frame["vpd"] = 0.28 * ti.SATURATION_CURVES["alduchov_eskridge"](frame["tair"].to_numpy()) / 1000.0
    ds, _ = build(tmp_path, write_site(tmp_path, frame))
    assert np.allclose(ds["RHair"][:, 0], 0.72, atol=1e-6)


def test_vpd_that_no_curve_explains_is_refused(tmp_path):
    frame = synthetic_tower()
    frame["RH"] = 72.0                                 # the VPD column no longer describes this air
    with pytest.raises(SystemExit, match="do not describe the same air"):
        build(tmp_path, write_site(tmp_path, frame))


def test_states_are_recentred_to_the_stamps(tmp_path):
    frame = synthetic_tower()
    ds, _ = build(tmp_path, write_site(tmp_path, frame))
    t = ds["Tair"][:, 0].astype(float) - 273.15
    expected = 0.5 * (frame["tair"].to_numpy()[:-1] + frame["tair"].to_numpy()[1:])   # begin: (k-1, k)
    assert np.allclose(t[1:], expected, atol=1e-4)
    assert abs(t[0] - frame["tair"].iloc[0]) < 1e-4                                 # one neighbour


def test_pressure_is_brought_down_to_the_ground(tmp_path):
    ds, _ = build(tmp_path, write_site(tmp_path, synthetic_tower()))
    p = ds["PSurf"][:, 0].astype(float)
    t = ds["Tair"][:, 0].astype(float)
    assert np.allclose(p, 98800.0 * np.exp(conv.GRAV * 30.0 / (conv.R_DRY * t)), rtol=2e-6)


def test_rain_and_shortwave_stay_interval_means(tmp_path):
    frame = synthetic_tower()
    ds, _ = build(tmp_path, write_site(tmp_path, frame))
    assert np.allclose(ds["Rainf"][:, 0], frame["PPT"].to_numpy() / DT, atol=1e-9)
    assert np.allclose(ds["SWdown"][:, 0], frame["Rs"].to_numpy(), atol=1e-3)


# ---------------------------------------------------------------------------------------------
# Gap filling is explicit and flagged.
# ---------------------------------------------------------------------------------------------
def test_short_gaps_are_interpolated_and_flagged(tmp_path):
    frame = synthetic_tower()
    frame.loc[200:202, "tair"] = np.nan                 # 3 records
    ds, _ = build(tmp_path, write_site(tmp_path, frame))
    qc = ds["Tair_qc"][:, 0]
    assert set(qc[200:203]) <= {tg.QC_SHORT} | set(qc[201:204])
    assert (qc[201:203] == tg.QC_SHORT).all()


def test_a_long_longwave_gap_is_filled_by_the_synthesis_regression(tmp_path):
    frame = synthetic_tower()
    frame.loc[500:1300, "Rl_dn"] = np.nan
    ds, report = build(tmp_path, write_site(tmp_path, frame))
    assert (ds["LWdown_qc"][600:1200, 0] == tg.QC_SYNTH_OR_MDV).all()
    lw = ds["LWdown"][:, 0].astype(float)
    truth = 440.0 + 2.0 * (frame["tair"].to_numpy() - 25.0)
    assert np.sqrt(np.mean((lw[600:1200] - truth[600:1200]) ** 2)) < 3.0       # T carries LW here
    assert report["fills"]["LWdown"]["method"] == "synthesis_regression"


def test_a_tower_without_longwave_gets_the_model_synthesis(tmp_path, capsys):
    frame = synthetic_tower().drop(columns=["Rl_dn"])
    site = write_site(tmp_path, frame)
    del site.variables["LWdown"]
    ds, report = build(tmp_path, site)
    assert report["fills"]["LWdown"]["fell_back_unfitted"]
    assert "WARNING" in capsys.readouterr().out
    t = ds["Tair"][:, 0].astype(float)
    # the file's LWdown is re-centred; compare an interior stamp against the synthesis of its two intervals
    p = mt.prepare(site)
    y = p["values"]
    lw = tg.synthesized_longwave(y["Tair"], y["RH"], y["PSurf"], y["SWdown"], p["mean_cosz"])
    assert np.allclose(ds["LWdown"][1:, 0], 0.5 * (lw[:-1] + lw[1:]), atol=1e-3)
    assert (ds["LWdown_qc"][:, 0] == tg.QC_SYNTH_OR_MDV).all()
    assert np.isfinite(t).all()


def test_missing_rain_is_an_error(tmp_path):
    frame = synthetic_tower()
    frame.loc[50, "PPT"] = np.nan
    with pytest.raises(SystemExit, match="does not invent rain"):
        build(tmp_path, write_site(tmp_path, frame))


def test_negative_wind_is_screened_then_filled(tmp_path):
    frame = synthetic_tower()
    frame.loc[300, "ubar"] = -1.0
    ds, report = build(tmp_path, write_site(tmp_path, frame))
    assert report["V4_bounds"]["Wind"]["out_of_bounds"] == 1
    assert ds["Wind"][300, 0] > 0


def test_a_stale_gapfill_setting_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="gapfill.longwave is not a setting"):
        write_site(tmp_path, synthetic_tower(), extra='longwave = "era5"')


# ---------------------------------------------------------------------------------------------
# The other input formats.
# ---------------------------------------------------------------------------------------------
def test_ameriflux_base_reads_timestamp_start_and_missing_values(tmp_path):
    frame = synthetic_tower()
    local = pd.to_datetime(frame["date"])
    base = pd.DataFrame({"TIMESTAMP_START": local.dt.strftime("%Y%m%d%H%M"),
                         "TIMESTAMP_END": (local + pd.Timedelta(minutes=30)).dt.strftime("%Y%m%d%H%M"),
                         "TA_1_1_1": frame["tair"], "RH_1_1_1": frame["RH"], "PA": frame["p_kpa"],
                         "P": frame["PPT"], "SW_IN": frame["Rs"], "LW_IN": frame["Rl_dn"], "WS": frame["ubar"]})
    base.loc[10, "TA_1_1_1"] = -9999
    path = tmp_path / "US-Xyz_HH_201402010000_201403180000.csv"
    with open(path, "w") as fh:
        fh.write("# Site: US-Xyz\n# Version: 1-5\n")
        base.to_csv(fh, index=False)
    toml = tmp_path / "base.toml"
    toml.write_text(f'''
[input]
format = "ameriflux_base"
path = "{path.name}"
[site]
latitude = {LAT}
longitude = {LON}
elevation = 100.0
[clock]
utc_offset = {OFFSET}
stamp = "begin"
timestep = 1800
[heights]
tq_height = 30.0
wind_height = 30.0
pressure_height = 0.0
[variables]
Tair = {{ column = "TA_1_1_1", units = "degC" }}
RH = {{ column = "RH_1_1_1", units = "%" }}
PSurf = {{ column = "PA", units = "kPa" }}
Rainf = {{ column = "P", units = "mm" }}
SWdown = {{ column = "SW_IN", units = "W m-2" }}
LWdown = {{ column = "LW_IN", units = "W m-2" }}
Wind = {{ column = "WS", units = "m s-1" }}
''')
    site = ti.read_site(str(toml))
    ds, report = build(tmp_path, site)
    assert report["V1_axis"]["records"] == len(frame)
    assert ds["Tair_qc"][10, 0] == tg.QC_SHORT or ds["Tair_qc"][11, 0] == tg.QC_SHORT


def test_ameriflux_timestamp_end_with_begin_stamp_is_refused(tmp_path):
    toml = tmp_path / "bad.toml"
    toml.write_text(f'''
[input]
format = "ameriflux_base"
path = "x.csv"
timestamp = "TIMESTAMP_END"
[site]
latitude = 0.0
longitude = 0.0
elevation = 0.0
[clock]
utc_offset = 0.0
stamp = "begin"
timestep = 1800
[heights]
tq_height = 2.0
wind_height = 2.0
pressure_height = 0.0
[variables]
Tair = {{ column = "TA", units = "degC" }}
''')
    with pytest.raises(SystemExit, match="TIMESTAMP_END"):
        ti.read_site(str(toml))


def test_fluxnet_provider_filled_values_are_flagged(tmp_path):
    frame = synthetic_tower()
    local = pd.to_datetime(frame["date"])
    flx = pd.DataFrame({"TIMESTAMP_START": local.dt.strftime("%Y%m%d%H%M"),
                        "TA_F": frame["tair"], "TA_F_QC": 0, "VPD_F": frame["vpd"] * 10.0, "VPD_F_QC": 0,
                        "PA_F": frame["p_kpa"], "PA_F_QC": 0, "P_F": frame["PPT"], "P_F_QC": 0,
                        "SW_IN_F": frame["Rs"], "SW_IN_F_QC": 0, "LW_IN_F": frame["Rl_dn"], "LW_IN_F_QC": 0,
                        "WS_F": frame["ubar"], "WS_F_QC": 0})
    flx.loc[100:140, "TA_F_QC"] = 2                     # the provider filled these
    flx.to_csv(tmp_path / "flx.csv", index=False)
    toml = tmp_path / "flx.toml"
    toml.write_text(f'''
[input]
format = "fluxnet"
path = "flx.csv"
[site]
latitude = {LAT}
longitude = {LON}
elevation = 100.0
[clock]
utc_offset = {OFFSET}
stamp = "begin"
timestep = 1800
[heights]
tq_height = 30.0
wind_height = 30.0
pressure_height = 0.0
[variables]
Tair = {{ column = "TA_F", units = "degC" }}
VPD = {{ column = "VPD_F", units = "hPa", curve = "alduchov_eskridge" }}
PSurf = {{ column = "PA_F", units = "kPa" }}
Rainf = {{ column = "P_F", units = "mm" }}
SWdown = {{ column = "SW_IN_F", units = "W m-2" }}
LWdown = {{ column = "LW_IN_F", units = "W m-2" }}
Wind = {{ column = "WS_F", units = "m s-1" }}
''')
    ds, report = build(tmp_path, ti.read_site(str(toml)))
    qc = ds["Tair_qc"][:, 0]
    assert (qc[101:141] == tg.QC_PROVIDER).all()
    # with no RH column, RH comes from the provider's VPD through its own curve, exactly
    rh = ds["RHair"][:, 0].astype(float)
    rh_true = frame["RH"].to_numpy() / 100.0
    assert np.allclose(rh[1:], 0.5 * (rh_true[:-1] + rh_true[1:]), atol=1e-5)
    assert report["V3_humidity"]["rh_from_vpd_records"] == len(frame)
    assert (ds["RHair_qc"][1:, 0] == tg.QC_FROM_VPD).all()             # flagged, not "observed"


# ---------------------------------------------------------------------------------------------
# The Python copy of the model's math.
# ---------------------------------------------------------------------------------------------
def test_humidity_round_trip_is_exact():
    t = np.array([280.0, 299.0, 305.0])
    p = np.array([95000.0, 98800.0, 101325.0])
    rh = np.array([0.3, 0.9, 1.0])
    q = conv.rh_to_specific_humidity(rh, t, p)
    assert np.allclose(conv.specific_humidity_to_vapor_pressure(q, p) / conv.sat_vapor_pressure(t), rh, atol=1e-13)


def test_solar_noon_follows_the_longitude_and_the_equation_of_time():
    day = np.datetime64("2014-03-21T00:00:00")
    stamps = day + np.arange(0, 86400, 60).astype("timedelta64[s]")
    cz = conv.solar_cosz(stamps, 0.0, -75.0)
    noon = stamps[np.argmax(cz)]
    doy = 80
    expected = 43200.0 + 75.0 * 240.0 - conv.equation_of_time(doy)       # 17:00 UTC minus EoT
    got = (noon - day).astype("timedelta64[s]").astype(float)
    assert abs(got - expected) <= 60.0


def test_window_mean_cosz_matches_a_fine_integral_by_day():
    start = np.array(["2014-06-01T16:00:00"], dtype="datetime64[s]")
    coarse = conv.window_mean_cosz(start, 3600.0, 9.15, -79.85)[0]
    fine = conv.window_mean_cosz(start, 3600.0, 9.15, -79.85, nsub=3600)[0]
    assert abs(coarse - fine) < 1e-3


def test_make_forcing_file_writes_the_dewpoint(tmp_path):
    import make_forcing_file as mf
    times = [np.datetime64("2024-01-01T01:00:00") + np.timedelta64(3600 * k, "s") for k in range(5)]
    times = [t.astype("datetime64[s]").astype(object) for t in times]
    series = {name: np.full(5, value) for name, value in
              dict(Tair=290.0, Tdew=285.0, PSurf=99000.0, u10=3.0, v10=-4.0, Rainf=0.0, SWdown=0.0,
                   LWdown=300.0).items()}
    arrays = mf.file_arrays([series])
    path = str(tmp_path / "era5.nc")
    mf.write_meds_forcing(path, times, [(42.4, -76.5, 367.5)], arrays,
                          dict(avg_convention="end", sw_input_kind="total", timestep_seconds=3600,
                               wind_meas_height_m=10.0, tq_height_m=2.0, height_above="zero_plane"))
    with Dataset(path) as ds:
        assert "Tdew" in ds.variables and "Qair" not in ds.variables
        assert ds.time_zone == "UTC" and ds.Conventions == "MEDS-forcing-1.1"
        assert np.allclose(ds["Wind"][:, 0], 5.0)
        assert ds["Tdew"].height == "2 m"
