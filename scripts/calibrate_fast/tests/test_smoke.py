# SPDX-License-Identifier: Apache-2.0
"""The calibration tool's smoke test with MEDS itself (MEDS_FAST_CALIBRATION_PLAN.md §7 P3): a
self-contained site -- the demography example's census, a synthetic forcing file and tower -- through
`calibrate_fast.py smoke`: a state chain, a trial restarted from it with its parameter record
checked, one Jacobian column that must move the output, and a repeated trial that must reproduce it
byte for byte. It runs once with meds_main (MEDS_MAIN names it) and once through the Python API
(MEDS_LIB names libmeds.so); CTest sets both, and a run whose model is missing is skipped. The
Python-API run also checks gate G8: the canopy of leaf solves (meds.canopy) on the window's own
hourly drivers against the model's GPP (revision plan §7.2)."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "python"))
from meds.config import RunConfig  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
import meds_forcing_file as mff  # noqa: E402

MEDS_MAIN = os.environ.get("MEDS_MAIN")
MEDS_LIB = os.environ.get("MEDS_LIB")
RUNNERS = [pytest.param(MEDS_MAIN, id="meds_main",
                        marks=pytest.mark.skipif(not MEDS_MAIN or not Path(MEDS_MAIN).exists(),
                                                 reason="MEDS_MAIN is not set to a meds_main executable")),
           pytest.param("python", id="python_api",
                        marks=pytest.mark.skipif(not MEDS_LIB or not Path(MEDS_LIB).exists(),
                                                 reason="MEDS_LIB is not set to a built libmeds.so"))]


def synthetic_tower(path, start="2001-06-01", days=6):
    idx = pd.date_range(start, periods=48 * days, freq="30min")
    h = idx.hour + idx.minute / 60.0
    sw = np.clip(800 * np.sin(np.pi * (h - 6) / 12), 0, None)
    df = pd.DataFrame({"date": idx, "Rs": sw, "Rs_dn": 0.13 * sw, "Rl_dn": 450.0, "Rnet": 0.7 * sw - 40,
                       "LE": np.where(sw > 0, 0.45 * sw, 5.0), "H": np.where(sw > 0, 0.2 * sw, -15.0),
                       "NEE": np.where(sw > 0, -12.0, 6.0), "gpp": np.where(sw > 0, 18.0, 0.0),
                       "ustar": 0.45, "FLAG": 1})
    df.to_csv(path, index=False)


def synthetic_forcing(path, start="2001-05-31", days=7):
    """An hourly ED_default file with a clear diurnal cycle (the fast output needs a forcing file)."""
    t = np.arange(np.datetime64(start), np.datetime64(start) + np.timedelta64(days, "D"),
                  np.timedelta64(1, "h"))
    hour = (t - t.astype("datetime64[D]")).astype("timedelta64[h]").astype(int).astype(float)
    sw = np.clip(900 * np.sin(np.pi * (hour - 6) / 12), 0, None)[:, None]
    tair = (293.0 + 5.0 * np.sin(np.pi * (hour - 9) / 12))[:, None]
    n = len(t)
    arrays = {"Tair": tair, "Qair": np.full((n, 1), 0.011), "PSurf": np.full((n, 1), 1.0e5),
              "Wind": np.full((n, 1), 2.5), "Rainf": np.zeros((n, 1)), "SWdown": sw,
              "LWdown": np.full((n, 1), 380.0)}
    attrs = {"title": "smoke-test forcing", "source": "synthetic", "history": "test_smoke.py",
             "timestep_seconds": 3600, "avg_convention": "begin", "sw_input_kind": "total",
             "wind_meas_height_m": 10.0}
    mff.write_forcing_file(path, t, [(42.4, -76.5, 300.0)], arrays, attrs)


@pytest.mark.parametrize("runner", RUNNERS)
def test_smoke(tmp_path, runner):
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    main = cfg.main
    main["fast"]["fast_biophysics_on"] = True
    synthetic_forcing(tmp_path / "forcing.nc")
    main["site"] = {"latitude": 42.4, "longitude": -76.5, "elevation": 300.0, "apply_elevation_lapse": False}
    main["forcing"] = {"forcing_on": True, "format": "ED_default", "path": str(tmp_path / "forcing.nc"),
                       "grid_index": 1, "grid_match": "explicit", "tq_height": 10.0, "wind_height": 10.0,
                       "height_above": "ground", "wind_exposure": "local", "timestep": "3600s",
                       "avg_convention": "begin", "sw_partition": "clearidx", "lwdown_source": "file",
                       "co2_source": "file",
                       "co2_file": str(ROOT / "data/co2/co2_cmip7_global_annual_1000-2022.txt"),
                       "recycle": False, "start_clamp": "error"}
    main["init"].update(init_mode=1, census_file=str(ROOT / "examples/example_demography/census_example.csv"))
    cfg.write(tmp_path)
    (tmp_path / "registry.toml").write_text(
        '[stomatal_g1]\nfile = "pft"\nkey = "pft.stomatal_g1"\npft = 3\nrange = [1.0, 12.0]\n'
        'transform = "log"\n')
    synthetic_tower(tmp_path / "tower.csv")
    decl = f"""
[base]
main = "main.toml"
registry = "registry.toml"
[overrides]
"trait_dynamics.trait_plasticity_on" = true
[tower]
path = "tower.csv"
utc_offset_h = 0.0
[tower.columns]
sw_in = "Rs"
sw_up = "Rs_dn"
lw_up = "Rl_dn"
rnet = "Rnet"
le = "LE"
h = "H"
nee = "NEE"
gpp = "gpp"
ustar = "ustar"
# the synthetic sun is not the site's (its noon is 12 UTC at 76.5 W): only the three targets the smoke checks
[targets.albedo]
on = false
[targets.lw_up]
on = false
[targets.rnet]
on = false
[targets.ef]
on = false
[targets.ustar]
on = false
[targets.le]
sigma_abs = 10.0
sigma_rel = 0.15
[targets.h]
sigma_abs = 10.0
sigma_rel = 0.15
[targets.gpp]
sigma_abs = 1.5
sigma_rel = 0.15
[windows]
days = 3
[windows.chains]
cal = "2001-06-01"
[[windows.list]]
name = "w"
start = "2001-06-02"
role = "cal"
chain = "cal"
"""
    (tmp_path / "calibration.toml").write_text(decl)
    g8 = ["--g8"] if runner == "python" else []
    res = subprocess.run([sys.executable, str(HERE / "calibrate_fast.py"), "smoke", "--site",
                          str(tmp_path / "calibration.toml"), "--work", str(tmp_path / "work"),
                          "--runner", runner, "--days", "2", "--key", "stomatal_g1"] + g8,
                         capture_output=True, text=True, timeout=1200)
    assert res.returncode == 0, res.stdout[-4000:] + res.stderr[-4000:]
    assert "repeat byte-identical: True; Jacobian column moves the output: True" in res.stdout
    if g8:
        assert "G8 canopy of leaf solves vs the model's GPP: {'pass': True" in res.stdout, res.stdout[-3000:]
