# SPDX-License-Identifier: Apache-2.0
"""meds.model runs a config as meds_main does.

- A run through the Python API (`python -m meds.model`) writes the same output files as meds_main,
  bit for bit, and the same parameter record.
- Two runs in one Python process each write their own parameter record.

Each run is two days from the demography example's census on the fast loop's constant reference
climate. MEDS_MAIN names the executable for the comparison (CTest sets it); without it that test is
skipped.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from meds.config import RunConfig, read_record
from meds.model import COMPLETED, run

ROOT = Path(__file__).resolve().parents[2]
MEDS_MAIN = os.environ.get("MEDS_MAIN")


def two_days(run_dir: Path) -> Path:
    """A two-day census run with hourly and daily output, written into run_dir."""
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    for key, value in {"run.start_time": "2001-06-01 00:00:00", "run.end_time": "2001-06-03 00:00:00",
                       "run.n_threads": 1, "fast.fast_biophysics_on": True, "forcing.forcing_on": False,
                       "init.init_mode": 1,
                       "init.census_file": str(ROOT / "data/census_example.csv"),
                       "state.write_state": False, "output.enabled": True, "output.dir": str(run_dir / "out"),
                       "output.prefix": "m", "output.fast.enabled": True, "output.daily.enabled": True,
                       "output.monthly.enabled": False, "output.annual.enabled": False}.items():
        cfg.set(key, value)
    (run_dir / "out").mkdir(parents=True, exist_ok=True)
    return cfg.write(run_dir)


def outputs(run_dir: Path) -> dict:
    """{file name: {variable: values}} for every netCDF file the run wrote."""
    out = {}
    for path in sorted((run_dir / "out").glob("*.nc")):
        with Dataset(path) as ds:
            out[path.name] = {v: np.array(ds[v][:]) for v in ds.variables}
    return out


def record(run_dir: Path) -> dict:
    """The run's parameter record, with the run's own directory written as <run>."""
    rec = read_record(run_dir / "out" / "m_parameters.csv",
                      {run_dir / "main.toml": "main", run_dir / "pft.toml": "pft"})
    return {k: (p, v.replace(str(run_dir), "<run>") if isinstance(v, str) else v)
            for k, (p, v) in rec.items()}


@pytest.mark.skipif(not MEDS_MAIN or not Path(MEDS_MAIN).exists(),
                    reason="MEDS_MAIN is not set to a meds_main executable")
def test_python_run_matches_meds_main(tmp_path):
    logs = {}
    for name, argv in (("exe", [MEDS_MAIN]), ("py", [sys.executable, "-m", "meds.model"])):
        main = two_days(tmp_path / name)
        res = subprocess.run(argv + [str(main)], cwd=tmp_path / name, capture_output=True, text=True,
                             timeout=600)
        assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
        logs[name] = res.stdout
    for log in logs.values():
        assert log.rstrip().splitlines()[-1] == COMPLETED
    a, b = outputs(tmp_path / "exe"), outputs(tmp_path / "py")
    assert a.keys() == b.keys() and a
    differ = [(f, v) for f in a for v in a[f]
              if a[f][v].shape != b[f][v].shape or not np.array_equal(a[f][v], b[f][v], equal_nan=True)]
    assert not differ, f"{len(differ)} variables differ between meds_main and meds.model: {differ[:8]}"
    assert record(tmp_path / "exe") == record(tmp_path / "py")


def test_each_run_in_a_process_keeps_its_own_record(tmp_path):
    for name in ("first", "second"):
        run(two_days(tmp_path / name), verbose=False)
    sources = {k[0] for k in read_record(tmp_path / "second" / "out" / "m_parameters.csv")}
    assert not [s for s in sources if str(tmp_path / "first") in s], sources
    assert record(tmp_path / "first").keys() == record(tmp_path / "second").keys()
