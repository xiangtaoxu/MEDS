# SPDX-License-Identifier: Apache-2.0
"""One trial: a frozen MEDS run over one window, restarted from that window's shared state with
a candidate parameter set (MEDS_FAST_CALIBRATION_PLAN.md §5.1, §7 P1 "trials").

A trial is a directory holding its own main and PFT TOML, made by parsing the base configs and
setting keys. Its name is a hash of the two files, so a repeated candidate reuses the finished
run. After the run the trial's parameter record (<prefix>_parameters.csv, written by meds_main)
must show every key the trial set, marked as set in the file, with the value written: a key that
is missing or defaulted means a misspelling or a key the model does not read, and the trial fails
rather than silently running the default.
"""
from __future__ import annotations

import csv
import re
import datetime as dt
import math
import os
import shutil
import threading
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

import tomlio

#: the hourly variables a trial writes: every target's model side (residuals.py)
TRIAL_VARIABLES = ("sw_in_fast", "sw_up_fast", "lw_up_fast", "rnet_fast", "le_flux_fast",
                   "h_flux_fast", "gpp_rate_fast", "nee_fast", "ustar_fast")
PREFIX = "t"
#: netCDF4 and HDF5 are not thread-safe, and the fit's starts run in threads: every read takes this
NC_LOCK = threading.RLock()
#: main-TOML keys holding paths, resolved against the base config's directory
PATH_KEYS = ("init.census_file", "init.pft_config", "forcing.path", "forcing.co2_file",
             "output.io_config")


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


def absolutize(main: dict, base_dir: Path) -> dict:
    """Resolve the base config's relative paths against its own directory."""
    for key in PATH_KEYS:
        v = tomlio.deep_get(main, key)
        if isinstance(v, str) and v not in ("", "none") and not os.path.isabs(v):
            tomlio.deep_set(main, key, str((base_dir / v).resolve()))
    return main


def set_param(main: dict, pft: dict, p, value: float, npft: int) -> None:
    """Deep-set one registry parameter in the parsed configs."""
    if p.file == "main":
        tomlio.deep_set(main, p.key, float(value))
        return
    arr = tomlio.deep_get(pft, p.key)
    if arr is None:
        if npft != 1:
            raise TrialError(f"{p.name}: '{p.key}' is not in the base PFT file, and with {npft} PFTs "
                             "the other PFTs' values are unknown -- add it to the base PFT file")
        arr = [float(value)]
    else:
        arr = [float(x) for x in arr]
        arr[p.pft - 1] = float(value)
    tomlio.deep_set(pft, p.key, arr)


def build_trial(base_main: dict, base_pft: dict, params, theta, window: Window, state_file: str,
                root: Path, overrides: dict | None = None, write_state: bool = False) -> Path:
    """Write a trial directory (or find the finished one) and return its path."""
    main, pft = tomlio.clone(base_main), tomlio.clone(base_pft)
    npft = len(tomlio.deep_get(pft, "pft.vcmax25", [0]))
    for p, v in zip(params, theta):
        set_param(main, pft, p, v, npft)
    for k, v in (overrides or {}).items():
        tomlio.deep_set(main, k, v)
    run = {"run.start_time": stamp(window.start), "run.end_time": stamp(window.end),
           "run.slow_on": False, "run.n_threads": 1,
           "init.init_mode": 2, "init.restart_file": str(state_file), "init.reacclimate_traits": True,
           "state.write_state": write_state, "state.output_prefix": PREFIX,
           "state.interval_years": 1000, "output.enabled": True, "output.prefix": PREFIX,
           "output.fast.enabled": True, "output.fast.file_chunk": "year",
           "output.daily.enabled": False, "output.monthly.enabled": False,
           "output.annual.enabled": False}
    for k, v in run.items():
        tomlio.deep_set(main, k, v)
    pft_text = tomlio.dumps(pft)
    tdir = root / f"{window.name}-{tomlio.digest(main, pft)[:16]}"
    if (tdir / "series.npz").exists():
        return tdir
    tdir.mkdir(parents=True, exist_ok=True)
    (tdir / "pft.toml").write_text(pft_text)
    (tdir / "output_variables.toml").write_text(
        "[variables]\n" + "".join(f'{v} = "F"\n' for v in TRIAL_VARIABLES))
    for k, v in {"init.pft_config": str(tdir / "pft.toml"), "output.dir": str(tdir / "out"),
                 "state.output_dir": str(tdir / "out"),
                 "output.io_config": str(tdir / "output_variables.toml")}.items():
        tomlio.deep_set(main, k, v)
    tomlio.write(tdir / "main.toml", main)
    (tdir / "out").mkdir(exist_ok=True)
    return tdir


def command(exe: str, tdir: Path) -> list[str]:
    return [exe, str(tdir / "main.toml")]


# ----- the parameter record ---------------------------------------------------------------------
def read_record(path: Path, main_path: Path, pft_path: Path) -> dict:
    """{(file, key, index): (present, value)} with file = "main" | "pft" | the source path."""
    out = {}
    names = {str(main_path): "main", str(pft_path): "pft"}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            src = names.get(row["source"], row["source"])
            try:
                val = float(row["value"])
            except ValueError:
                val = row["value"]
            out[(src, row["key"], int(row["index"]))] = (row["present"] == "true", val)
    return out


def check_record(tdir: Path, params, theta) -> None:
    rec = read_record(tdir / "out" / f"{PREFIX}_parameters.csv", tdir / "main.toml", tdir / "pft.toml")
    bad = []
    for p, v in zip(params, theta):
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
def read_hourly(out_dir: Path, utc_offset_h: float) -> pd.DataFrame:
    """The trial's hourly records on LOCAL time (the start of each hour)."""
    frames = []
    for path in sorted(out_dir.glob(f"{PREFIX}-F-*.nc")):
        with NC_LOCK, Dataset(path) as ds:
            cols = {v: np.asarray(ds[v][:], dtype=float).squeeze() for v in TRIAL_VARIABLES}
            when = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:],
                                       hour=ds["hour"][:], minute=ds["minute"][:]))
        frames.append(pd.DataFrame(cols, index=when + pd.Timedelta(hours=utc_offset_h)))
    if not frames:
        raise TrialError(f"no hourly output in {out_dir}")
    df = pd.concat(frames).sort_index()
    return df.where(df.abs() < 1e30)


def finish(tdir: Path, params, theta, utc_offset_h: float, keep_netcdf: bool = False) -> pd.DataFrame:
    """Check a completed trial and cache its hourly series (series.npz); raise TrialError if it
    failed."""
    log = (tdir / "run.log").read_text(errors="replace") if (tdir / "run.log").exists() else ""
    if "OK: simulation completed" not in log:
        tail = "\n".join(log.splitlines()[-15:])
        raise TrialError(f"{tdir.name}: the run did not complete\n{tail}")
    #----- a set that breaks conservation is not a good run, however well it fits
    for which, fails in re.findall(r"budget\[(whole_\w+)\].*fails = (\d+)/", log):
        if int(fails) > 0:
            raise TrialError(f"{tdir.name}: the {which} budget breached tolerance {fails} times")
    check_record(tdir, params, theta)
    df = read_hourly(tdir / "out", utc_offset_h)
    if df.isna().any().any():
        raise TrialError(f"{tdir.name}: missing or non-finite values in the hourly output")
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
