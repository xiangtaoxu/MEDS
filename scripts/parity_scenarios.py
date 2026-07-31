#!/usr/bin/env python3
"""Emit the Phase-B scenario base configs for the integrator-parity study.

docs/dev_plans/docs/dev_plans/archive/MEDS_INTEGRATOR_PARITY.md [RETIRED], Phase B.  Four cells, a 2x2 over SEASON x STAND
STRUCTURE, because those are the two axes that independently drive integrator stress:

  season          the stiffness regime.  Winter is cold + leaf-off, and it is where RK45's explicit
                  march goes stiff and hybrid-rescues to split; summer is high LAI with strong
                  transpiration and a fast (tau ~ 130 s) CAS.
  stand structure the leaf<->CAS coupling strength.  A near-bare patch has almost none (and almost
                  no GPP, so score those cells on soil temperature / freeze front / theta, NOT on
                  carbon); the spun-up stand has the dense-canopy coupling that meds_fast_split's
                  own comment calls "marginally unstable".

Each scenario is ONE base config.  The scheme axis is then applied by numerics_sweep.py, which
already knows how to sweep schemes, build a reference and score the daily stream -- this script
deliberately does not duplicate that.

DEMOGRAPHY IS FROZEN (`[run].slow_on = false`) in all four.  The target of Phase B is biophysics
fidelity; leaving the slow tier on would let a month of growth, phenology and the monthly
recruit/fuse cadence back in and reintroduce exactly the confound the phase exists to remove.

STATE RESTART, not census, for the spun-up cells.  A census carries no fast-loop state, so all
three integrators would spend the month relaxing a cold soil column and the initialization
transient would dominate a one-month window.  The known cost is that the 50-yr spin-up itself ran
on `split`, so a state restart starts every scheme from split's own converged fast state -- a
shared starting point, so the comparison is fair, but split starts at home.  Say so in any writeup.

Usage
-----
  scripts/parity_scenarios.py --out runs/parity_phaseB
  scripts/numerics_sweep.py --base runs/parity_phaseB/b3_stand_winter.toml \
      --out runs/parity_phaseB/sweep_b3 --parity --schemes split ark rk45 --dt 1800
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from numerics_sweep import deep_set, dumps, tomllib  # noqa: E402  (same-dir helper module)

REPO = Path(__file__).resolve().parent.parent
EXAMPLE = REPO / "examples" / "example_biophysics"

# The four cells.  Winter runs January and summer runs July of model year 2074, for BOTH stand
# types, so the forcing is identical across the stand axis and only the vegetation/soil state
# differs.  2074 is 50 recycle wraps past the forcing file's own year -- safe since PR #69 made the
# recycle window declared and validated (it reproduces the exact-forcing diel cycle to 0.03 W/m2).
SCENARIOS = {
    "b1_bare_winter":  {"start": "2074-01-01", "end": "2074-02-01", "init_mode": 0, "restart": None},
    "b2_bare_summer":  {"start": "2074-07-01", "end": "2074-08-01", "init_mode": 0, "restart": None},
    "b3_stand_winter": {"start": "2074-01-01", "end": "2074-02-01", "init_mode": 2,
                        "restart": "spinup-S-20740101000000.nc"},
    "b4_stand_summer": {"start": "2074-07-01", "end": "2074-08-01", "init_mode": 2,
                        "restart": "spinup-S-20740701000000.nc"},
}


def build(base: dict, name: str, spec: dict, out_dir: Path) -> dict:
    cfg = {k: (dict(v) if isinstance(v, dict) else v) for k, v in base.items()}
    cfg = {k: ({kk: vv for kk, vv in v.items()} if isinstance(v, dict) else v) for k, v in cfg.items()}

    deep_set(cfg, "run.start_time", spec["start"])
    deep_set(cfg, "run.end_time", spec["end"])
    deep_set(cfg, "run.slow_on", False)          # biophysics fidelity: freeze demography

    deep_set(cfg, "init.init_mode", spec["init_mode"])
    if spec["restart"] is not None:
        deep_set(cfg, "init.restart_file", str(EXAMPLE / spec["restart"]))

    # Absolute paths: numerics_sweep writes generated configs into its own --out tree, so anything
    # relative in the base config would resolve against the wrong directory and fail obscurely.
    deep_set(cfg, "forcing.path", str(REPO / "data" / "forcing" / "ithaca_forcing.nc"))
    deep_set(cfg, "init.pft_config", str(EXAMPLE / "pft_parameters.toml"))
    deep_set(cfg, "output.io_config", str(EXAMPLE / "output_variables.toml"))

    deep_set(cfg, "fast.dt_fast", "1800s")
    # Every diagnostic Phase B scores, plus the Phase-A health counters.  numerics_sweep sets these
    # too, but the base config must carry them so a scenario is runnable on its own.
    deep_set(cfg, "output.enabled", True)
    deep_set(cfg, "output.energy_fluxes", True)
    deep_set(cfg, "output.water_fluxes", True)
    deep_set(cfg, "output.numerics", True)
    deep_set(cfg, "output.daily.enabled", True)
    deep_set(cfg, "output.monthly.enabled", False)
    deep_set(cfg, "output.annual.enabled", False)
    # HOURLY sub-daily output: 2 * dt_fast(1800 s) = 3600 s.  A one-month daily series is only ~30
    # points and averages away the diel structure that separates the schemes in the first place.
    deep_set(cfg, "output.fast.enabled", True)
    deep_set(cfg, "output.fast_interval_steps", 2)
    deep_set(cfg, "output.fast.file_chunk", "month")
    deep_set(cfg, "output.dir", str(out_dir / name))
    deep_set(cfg, "output.prefix", "s")
    deep_set(cfg, "io.output_dir", str(out_dir / name))
    deep_set(cfg, "io.output_prefix", name)
    deep_set(cfg, "io.write_state", False)
    return cfg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=REPO / "runs" / "parity_phaseB")
    ap.add_argument("--base", type=Path, default=EXAMPLE / "meds_config_july.toml",
                    help="config to derive from (its forcing/soil/PFT blocks are inherited verbatim)")
    args = ap.parse_args(argv)

    base = tomllib.load(open(args.base, "rb"))
    args.out.mkdir(parents=True, exist_ok=True)
    for name, spec in SCENARIOS.items():
        cfg = build(base, name, spec, args.out)
        path = args.out / f"{name}.toml"
        path.write_text(dumps(cfg))
        print(f"  {path}   {spec['start']}..{spec['end']}  init_mode={spec['init_mode']}")
    print(f"\n{len(SCENARIOS)} scenario configs in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
