# SPDX-License-Identifier: Apache-2.0
"""The shared states (MEDS_FAST_CALIBRATION_PLAN.md §5.1, D7).

A chain is a frozen run from the census at the chain's start date that stops at each of its
windows' start dates in turn, writing a state there and restarting from it to the next. Every trial
of a window restarts from that window's state, so the soil water and temperature a window starts
with are those of a spun-up run, not the census's initialization. The chain runs with one
parameter set (the default, or the MAP when the fit refreshes it, §6.2).
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import tomlio
from trials import stamp, set_param, TrialError

PREFIX = "c"


def chain_segments(chain_start: dt.datetime, starts: list) -> list[tuple]:
    """(from, to) spans of a chain passing each window start, in time order."""
    out, t0 = [], chain_start
    for t in sorted(set(starts)):
        if t < chain_start:
            raise ValueError(f"window start {t} is before its chain's start {chain_start}")
        if t > t0:
            out.append((t0, t))
            t0 = t
    return out


def run_chain(name, chain_start, windows, base_main, base_pft, params, theta, root: Path, pool, exe,
              overrides=None, timeout=7200) -> dict:
    """Run one chain; return {window name: state file}. A chain already run with the same
    parameters is reused."""
    starts = [w.start for w in windows]
    key = tomlio.digest(name, str(chain_start), [str(s) for s in sorted(set(starts))],
                        [(p.name, float(v)) for p, v in zip(params, theta)], overrides or {},
                        base_main, base_pft)
    cdir = root / f"chain-{name}-{key[:12]}"
    cdir.mkdir(parents=True, exist_ok=True)
    npft = len(tomlio.deep_get(base_pft, "pft.vcmax25", [0]))
    state_at = {}
    prev_state = None
    for i, (t0, t1) in enumerate(chain_segments(chain_start, starts)):
        sdir = cdir / f"seg{i:02d}"
        state = sdir / "out" / f"{PREFIX}-S-{t1.strftime('%Y%m%d%H%M%S')}.nc"
        if not state.exists():
            main, pft = tomlio.clone(base_main), tomlio.clone(base_pft)
            for p, v in zip(params, theta):
                set_param(main, pft, p, v, npft)
            for k, v in (overrides or {}).items():
                tomlio.deep_set(main, k, v)
            sdir.mkdir(parents=True, exist_ok=True)
            (sdir / "out").mkdir(exist_ok=True)
            (sdir / "pft.toml").write_text(tomlio.dumps(pft))
            settings = {"run.start_time": stamp(t0), "run.end_time": stamp(t1), "run.slow_on": False,
                        "run.n_threads": 1, "init.init_mode": 1 if prev_state is None else 2,
                        "init.restart_file": "none" if prev_state is None else str(prev_state),
                        "init.pft_config": str(sdir / "pft.toml"), "state.write_state": True,
                        "state.output_dir": str(sdir / "out"), "state.output_prefix": PREFIX,
                        "state.interval_years": 1000, "output.enabled": False}
            for k, v in settings.items():
                tomlio.deep_set(main, k, v)
            tomlio.write(sdir / "main.toml", main)
            from pool import Task
            res = pool.run([Task(f"{name}-{i}", [exe, str(sdir / "main.toml")], str(sdir),
                                 str(sdir / "run.log"), timeout)])
            status = next(iter(res.values()))[0]
            if status != "ok" or not state.exists():
                tail = "\n".join((sdir / "run.log").read_text(errors="replace").splitlines()[-15:])
                raise TrialError(f"chain {name} segment {i} ({t0} -> {t1}) failed: {status}\n{tail}")
        prev_state = state
        for w in windows:
            if w.start == t1:
                state_at[w.name] = state
    missing = [w.name for w in windows if w.name not in state_at]
    if missing:
        raise TrialError(f"chain {name}: no state for windows {missing} (a window at the chain start?)")
    return state_at
