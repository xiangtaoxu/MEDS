#!/usr/bin/env python3
"""Score integrator FIDELITY on the sub-daily (FAST) stream of a numerics_sweep run.

docs/dev_plans/docs/dev_plans/archive/MEDS_INTEGRATOR_PARITY.md [RETIRED], Phase B.

Reads the hourly FAST-tier output of every scheme cell in a sweep directory and reports:

  PAIRWISE   how far the schemes are from EACH OTHER at production dt_fast.  This is the primary
             Phase-B number.  All three share the same Category-0 coefficient freeze at a common
             dt_fast, so a pairwise difference is integrator disagreement with the dominant error
             source held common -- which is exactly the question "do they agree?".
  VS REF     how far each is from the refined reference cell.  Context, not the headline: it mixes
             time-stepping error with the freeze error that refining dt_fast also removes.

Why the FAST tier and not the daily stream numerics_sweep itself scores: the sub-daily
temperatures (cas_temp_site, soil_temp_top_site, leaf_temp_cohort_fast) exist ONLY here, and a
one-month daily series is ~30 points that average away the diel structure the schemes differ on.

Cells are aligned on their (year, month, day, hour) stamps rather than by index, so a reference at
a different dt_fast is compared only where the two axes genuinely coincide.

Usage
-----
  scripts/parity_fidelity.py runs/parity_phaseB/sweep_ifx_b4_stand_summer
  scripts/parity_fidelity.py runs/parity_phaseB/sweep_ifx_*        # one block per sweep
"""

from __future__ import annotations

import argparse
import glob
import math
import sys
from pathlib import Path

try:
    import netCDF4
except ModuleNotFoundError:
    netCDF4 = None

# (variable, label, unit).  Chosen to span the coupled system: the two prognostic reservoirs whose
# disagreement propagates (CAS, soil surface), the two fluxes that carry it (LE, H), and the carbon
# response that any biophysical difference ultimately shows up in (GPP).
VARS = [
    ("cas_temp_site",      "CAS temp",   "K"),
    ("soil_temp_top_site", "soil-top T", "K"),
    ("le_flux_fast",       "LE",         "W/m2"),
    ("h_flux_fast",        "H",          "W/m2"),
    ("gpp_rate_fast",      "GPP",        "umol/m2/s"),
]


def read_fast(cell_dir: Path) -> dict:
    """Concatenate the FAST-tier files in time order -> {var: [values]} plus a 'stamp' key."""
    out: dict[str, list] = {}
    files = sorted(glob.glob(str(cell_dir / "s-F-*.nc")))
    for f in files:
        with netCDF4.Dataset(f) as ds:
            n = len(ds.dimensions["time"])
            for v, _, _ in VARS:
                if v in ds.variables:
                    out.setdefault(v, []).extend(float(x) for x in ds.variables[v][:])
            # An (y,m,d,h) stamp per record, so cells on different dt_fast align by TIME rather than
            # by position -- an index-based comparison would silently pair unrelated hours.
            cal = {k: [int(x) for x in ds.variables[k][:]] for k in ("year", "month", "day", "hour")
                   if k in ds.variables}
            if cal:
                out.setdefault("stamp", []).extend(zip(*(cal[k] for k in ("year", "month", "day", "hour"))))
            else:
                out.setdefault("stamp", []).extend([(0, 0, 0, i) for i in range(n)])
    return out


def paired(a: dict, b: dict, var: str):
    """Values of `var` from a and b at their COMMON stamps, in stamp order."""
    if var not in a or var not in b:
        return [], []
    ia = {s: i for i, s in enumerate(a["stamp"])}
    ib = {s: i for i, s in enumerate(b["stamp"])}
    common = sorted(set(ia) & set(ib))
    return [a[var][ia[s]] for s in common], [b[var][ib[s]] for s in common]


def stats(x, y):
    if not x:
        return None
    d = [u - v for u, v in zip(x, y)]
    rmse = math.sqrt(sum(t * t for t in d) / len(d))
    return rmse, max(abs(t) for t in d), len(d)


def report(sweep: Path):
    cells = {}
    ref = None
    for d in sorted(p for p in sweep.iterdir() if p.is_dir()):
        data = read_fast(d)
        if not data:
            continue
        if d.name.startswith("ref_"):
            ref = data
        else:
            cells[d.name.split("_")[0]] = data       # "split_1800_full_i" -> "split"
    if not cells:
        print(f"  (no FAST output under {sweep})")
        return

    print(f"\n===== {sweep.name} =====")
    order = [s for s in ("split", "ark", "rk45") if s in cells]
    pairs = [(order[i], order[j]) for i in range(len(order)) for j in range(i + 1, len(order))]

    head = f"  {'variable':12s} {'unit':10s}"
    for a, b in pairs:
        head += f" {a[:4]+'~'+b[:4]:>16s}"
    if ref:
        head += "".join(f" {s+' vs ref':>16s}" for s in order)
    print(head)

    for var, label, unit in VARS:
        line = f"  {label:12s} {unit:10s}"
        any_data = False
        for a, b in pairs:
            st = stats(*paired(cells[a], cells[b], var))
            line += f" {(f'{st[0]:.4g}' if st else '-'):>16s}"
            any_data = any_data or bool(st)
        if ref:
            for s in order:
                st = stats(*paired(cells[s], ref, var))
                line += f" {(f'{st[0]:.4g}' if st else '-'):>16s}"
        if any_data:
            print(line)

    # n is reported once: it is the same for every variable of a given pair, and printing it per
    # row would crowd out the numbers that matter.
    st = stats(*paired(cells[order[0]], cells[order[-1]], VARS[0][0]))
    if st:
        print(f"  (RMSE over {st[2]} paired hourly records; pairwise = same freeze cadence, "
              f"so integrator-only)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweeps", nargs="+", type=Path)
    args = ap.parse_args(argv)
    if netCDF4 is None:
        print("netCDF4 not available", file=sys.stderr)
        return 1
    for s in args.sweeps:
        report(s)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
