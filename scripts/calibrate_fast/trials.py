# SPDX-License-Identifier: Apache-2.0
"""One trial: a frozen MEDS run over one window, restarted from that window's shared state with
a candidate parameter set (MEDS_FAST_CALIBRATION_PLAN.md §5.1, §7 P1 "trials").

A trial is a directory holding its own main and PFT TOML, made from the base configuration with
meds.config by setting keys. Its name is a hash of the two, so a repeated candidate reuses the
finished run. A trial runs through the Python API (`python -m meds.model`) or the meds_main
executable, whichever the fit's runner names; both write the same files. After the run the trial's
parameter record (<prefix>_parameters.csv) must show every key the trial set, marked as set in the
file, with the value written: a key that is missing or defaulted means a misspelling or a key the
model does not read, and the trial fails rather than silently running the default.

A trial writes its fast output at the tower's own interval (the site sets [output].fast_interval_steps
on the base configuration), stamped by each record's UTC start, the index the observations are on.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import re
import shutil
import sys
import threading
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from meds.config import RunConfig, read_record

#: the fast variables a trial writes: every target's model side (residuals.py)
TRIAL_VARIABLES = ("sw_in_fast", "sw_up_fast", "lw_up_fast", "rnet_fast", "le_flux_fast",
                   "h_flux_fast", "gpp_rate_fast", "nee_fast", "ustar_fast")
PREFIX = "t"
#: the runner that runs a config through the Python API instead of an executable
PYTHON_RUNNER = "python"
#: the line meds_main and meds.model print when a run ends well (meds.model.COMPLETED)
COMPLETED = "OK: simulation completed"
#: netCDF4 and HDF5 are not thread-safe, and the fit's starts run in threads: every read takes this
NC_LOCK = threading.RLock()


class TrialError(RuntimeError):
    pass


@dataclass
class Window:
    name: str
    start: dt.datetime          # UTC, the run's start_time
    days: int
    role: str                   # "cal" or "val"
    chain: str

    @property
    def end(self) -> dt.datetime:
        return self.start + dt.timedelta(days=self.days)


def stamp(t: dt.datetime) -> str:
    return t.strftime("%Y-%m-%d %H:%M:%S")


def digest(*objs) -> str:
    """A content hash that ignores key order: the same settings give the same digest however the
    tables were assembled."""
    text = json.dumps(objs, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha1(text.encode()).hexdigest()


def set_param(cfg: RunConfig, p, value: float) -> None:
    """Set one registry parameter in a run's configuration."""
    cfg.set(p.key, float(value), file=p.file, pft=p.pft if p.file == "pft" else None)


def with_params(base: RunConfig, params, theta, overrides: dict | None = None) -> RunConfig:
    """A copy of the base configuration with a parameter set and the calibration's overrides."""
    cfg = base.copy()
    for p, v in zip(params, theta):
        if p.file != "obs":                    # an observation key enters the residuals, not the model
            set_param(cfg, p, v)
    for k, v in (overrides or {}).items():
        cfg.set(k, v)
    return cfg


def build_trial(base: RunConfig, params, theta, window: Window, state_file: str, root: Path,
                overrides: dict | None = None, write_state: bool = False) -> Path:
    """Write a trial directory (or find the finished one) and return its path."""
    cfg = with_params(base, params, theta, overrides)
    run = {"run.start_time": stamp(window.start), "run.end_time": stamp(window.end),
           "run.slow_on": False, "run.n_threads": 1,
           "init.init_mode": 2, "init.restart_file": str(state_file), "init.reacclimate_traits": True,
           "state.write_state": write_state, "state.output_prefix": PREFIX,
           "state.interval_years": 1000, "output.enabled": True, "output.prefix": PREFIX,
           "output.fast.enabled": True, "output.fast.file_chunk": "year",
           "output.daily.enabled": False, "output.monthly.enabled": False,
           "output.annual.enabled": False}
    for k, v in run.items():
        cfg.set(k, v)
    tag = digest(cfg.main, cfg.pft)
    tdir = root / f"{window.name}-{tag[:16]}"
    if (tdir / "series.npz").exists():
        return tdir
    (tdir / "out").mkdir(parents=True, exist_ok=True)
    (tdir / "output_variables.toml").write_text(
        "[variables]\n" + "".join(f'{v} = "F"\n' for v in TRIAL_VARIABLES))
    for k, v in {"output.dir": str(tdir / "out"), "state.output_dir": str(tdir / "out"),
                 "output.io_config": str(tdir / "output_variables.toml")}.items():
        cfg.set(k, v)
    cfg.write(tdir)
    return tdir


def command(runner: str, config: Path) -> list[str]:
    """The command that runs a config: through the Python API when the runner is "python", else
    the meds_main executable the runner names."""
    if runner == PYTHON_RUNNER:
        return [sys.executable, "-m", "meds.model", str(config)]
    return [runner, str(config)]


# ----- the parameter record ---------------------------------------------------------------------
def check_record(tdir: Path, params, theta) -> None:
    rec = read_record(tdir / "out" / f"{PREFIX}_parameters.csv",
                      {tdir / "main.toml": "main", tdir / "pft.toml": "pft"})
    bad = []
    for p, v in zip(params, theta):
        if p.file == "obs":
            continue
        idx = p.pft if p.file == "pft" else 0
        hit = rec.get((p.file, p.key, idx))
        if hit is None:
            bad.append(f"{p.name} ({p.key}): not read by the model")
        elif not hit[0]:
            bad.append(f"{p.name} ({p.key}): read as a default, not from the trial's file")
        elif not math.isclose(hit[1], float(v), rel_tol=1e-14, abs_tol=0.0):
            bad.append(f"{p.name} ({p.key}): the model read {hit[1]!r}, the trial wrote {float(v)!r}")
    if bad:
        raise TrialError(f"{tdir.name}: parameter record check failed:\n  " + "\n  ".join(bad))


# ----- the output --------------------------------------------------------------------------------
def read_series(out_dir: Path, step: float) -> pd.DataFrame:
    """The trial's fast records on the UTC start of each record; they must be `step` seconds apart,
    the tower's interval, or the pairing with the observations would be wrong."""
    frames = []
    for path in sorted(out_dir.glob(f"{PREFIX}-F-*.nc")):
        with NC_LOCK, Dataset(path) as ds:
            cols = {v: np.asarray(ds[v][:], dtype=float).squeeze() for v in TRIAL_VARIABLES}
            when = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:],
                                       hour=ds["hour"][:], minute=ds["minute"][:]))
        frames.append(pd.DataFrame(cols, index=when))
    if not frames:
        raise TrialError(f"no fast output in {out_dir}")
    df = pd.concat(frames).sort_index()
    gaps = np.unique(np.diff(df.index.values).astype("timedelta64[s]").astype(float))
    if len(df) > 1 and not np.allclose(gaps, step):
        raise TrialError(f"{out_dir}: the fast output is {gaps} s apart, not the tower's {step:g} s")
    return df.where(df.abs() < 1e30)


def finish(tdir: Path, params, theta, step: float, keep_netcdf: bool = False) -> pd.DataFrame:
    """Check a completed trial and cache its fast series (series.npz); raise TrialError if it
    failed. `step` is the tower's interval."""
    log = (tdir / "run.log").read_text(errors="replace") if (tdir / "run.log").exists() else ""
    if COMPLETED not in log:
        tail = "\n".join(log.splitlines()[-15:])
        raise TrialError(f"{tdir.name}: the run did not complete\n{tail}")
    #----- a set that breaks conservation is not a good run, however well it fits
    for which, fails in re.findall(r"budget\[(whole_\w+)\].*fails = (\d+)/", log):
        if int(fails) > 0:
            raise TrialError(f"{tdir.name}: the {which} budget breached tolerance {fails} times")
    check_record(tdir, params, theta)
    df = read_series(tdir / "out", step)
    if df.isna().any().any():
        raise TrialError(f"{tdir.name}: missing or non-finite values in the fast output")
    np.savez(tdir / "series.npz", index=df.index.values.astype("datetime64[s]").astype(np.int64),
             **{v: df[v].values for v in TRIAL_VARIABLES})
    if not keep_netcdf:
        for f in (tdir / "out").glob(f"{PREFIX}-F-*.nc"):
            f.unlink()
    return df


def load_series(tdir: Path) -> pd.DataFrame:
    z = np.load(tdir / "series.npz")
    idx = pd.to_datetime(z["index"], unit="s")
    return pd.DataFrame({v: z[v] for v in TRIAL_VARIABLES}, index=idx)


def clean(tdir: Path) -> None:
    shutil.rmtree(tdir, ignore_errors=True)
