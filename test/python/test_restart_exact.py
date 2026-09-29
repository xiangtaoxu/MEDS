# SPDX-License-Identifier: Apache-2.0
"""A restart continues exactly as the run that wrote the state would have, and a restart that
re-acclimates its traits ([init].reacclimate_traits) gets the traits a census start gives.

Each case runs meds_main from the demography example's census on the fast loop's constant
reference climate (no forcing file) and compares the final state files variable by variable,
bit for bit. MEDS_MAIN names the executable; CTest sets it.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

try:
    import tomllib
except ModuleNotFoundError:  # py3.10 and older
    import tomli as tomllib
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[2]
MEDS_MAIN = os.environ.get("MEDS_MAIN", str(ROOT / "build-ifx" / "meds_main"))


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in v) + "]"
    raise TypeError(type(v).__name__)


def dumps(d, prefix=""):
    out = []
    scalars = {k: v for k, v in d.items() if not isinstance(v, dict)}
    if prefix:
        out.append(f"[{prefix}]")
    out += [f"{k} = {_fmt(v)}" for k, v in scalars.items()] + [""]
    for k, v in d.items():
        if isinstance(v, dict):
            out.append(dumps(v, f"{prefix}.{k}" if prefix else k))
    return "\n".join(out)


def base_config(run_dir, start, end, *, slow_on, restart_file=None, reacclimate=False,
                plasticity=True, canopy_water=False):
    cfg = tomllib.loads((ROOT / "meds_config_main.toml").read_text())
    cfg["run"].update(start_time=start, end_time=end, slow_on=slow_on, n_threads=1)
    cfg["fast"]["fast_biophysics_on"] = True
    cfg["fast"]["canopy_water_on"] = canopy_water
    cfg.setdefault("forcing", {})["forcing_on"] = False
    cfg["init"].update(init_mode=2 if restart_file else 1, restart_file=restart_file or "none",
                       census_file=str(ROOT / "examples/example_demography/census_example.csv"),
                       pft_config=str(run_dir / "pft.toml"))
    if reacclimate:
        cfg["init"]["reacclimate_traits"] = True
    cfg["state"].update(output_dir=str(run_dir / "out"), output_prefix="s", write_state=True)
    cfg["output"]["enabled"] = False
    cfg.setdefault("trait_dynamics", {})["trait_plasticity_on"] = plasticity
    return cfg


def run(run_dir, cfg, pft_edit=None):
    run_dir.mkdir(parents=True, exist_ok=True)
    pft = (ROOT / "meds_config_pft.toml").read_text()
    if pft_edit:
        pft = pft_edit(pft)
    (run_dir / "pft.toml").write_text(pft)
    (run_dir / "main.toml").write_text(dumps(cfg))
    res = subprocess.run([MEDS_MAIN, str(run_dir / "main.toml")], cwd=run_dir, capture_output=True,
                         text=True, timeout=600)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    states = sorted((run_dir / "out").glob("s-S-*.nc"))
    assert states, f"no state file in {run_dir / 'out'}"
    return states[-1]


def read_state(path):
    with Dataset(path) as ds:
        return {v: np.array(ds[v][:]) for v in ds.variables}, ds.getncattr("restructure_pending")


def assert_same_state(a, b):
    va, pa = read_state(a)
    vb, pb = read_state(b)
    assert pa == pb
    assert va.keys() == vb.keys()
    differ = [k for k in va if not np.array_equal(va[k], vb[k], equal_nan=True)]
    assert not differ, f"{len(differ)} variables differ between {a.name} and {b.name}: {differ}"


def vcmax_times(factor):
    """A PFT-file edit scaling every PFT's vcmax25."""
    def edit(text):
        out = []
        for line in text.splitlines():
            key = line.split("=")[0].strip()
            if key == "vcmax25" and "[" in line:
                head, rest = line.split("[", 1)
                vals, tail = rest.split("]", 1)
                new = ", ".join(repr(float(x) * factor) for x in vals.split(","))
                line = f"{head}[{new}]{tail}"
            out.append(line)
        return "\n".join(out) + "\n"
    return edit


@pytest.mark.parametrize("slow_on,canopy_water", [(False, True), (True, False)])
def test_split_matches_unsplit(tmp_path, slow_on, canopy_water):
    whole = run(tmp_path / "whole", base_config(tmp_path / "whole", "2001-06-01", "2001-06-03",
                                                  slow_on=slow_on, canopy_water=canopy_water))
    first = run(tmp_path / "first", base_config(tmp_path / "first", "2001-06-01", "2001-06-02",
                                                  slow_on=slow_on, canopy_water=canopy_water))
    second = run(tmp_path / "second", base_config(tmp_path / "second", "2001-06-02", "2001-06-03",
                                                    slow_on=slow_on, canopy_water=canopy_water,
                                                    restart_file=str(first)))
    assert whole.name == second.name
    assert_same_state(whole, second)


@pytest.mark.parametrize("plasticity", [True, False])
def test_reacclimate_matches_census_start(tmp_path, plasticity):
    """A frozen stand: the state written by a census start with vcmax25 x1, restarted with x1.3 and
    re-acclimated, carries the traits a census start with x1.3 gives, and keeps its leaf area."""
    day = ("2001-06-01", "2001-06-02")
    written = run(tmp_path / "a", base_config(tmp_path / "a", *day, slow_on=False, plasticity=plasticity))
    census_b = run(tmp_path / "b", base_config(tmp_path / "b", *day, slow_on=False, plasticity=plasticity),
                   pft_edit=vcmax_times(1.3))
    restart_b = run(tmp_path / "r", base_config(tmp_path / "r", "2001-06-02", "2001-06-03", slow_on=False,
                                                  plasticity=plasticity, restart_file=str(written),
                                                  reacclimate=True),
                    pft_edit=vcmax_times(1.3))
    a, _ = read_state(written)
    b, _ = read_state(census_b)
    r, _ = read_state(restart_b)
    for trait in ("sla", "vcmax25", "rd25", "llspan"):
        assert np.array_equal(r[trait], b[trait]), trait
    assert not np.array_equal(r["vcmax25"], a["vcmax25"]), "the changed vcmax25 did not reach the cohorts"
    assert np.array_equal(r["leaf_area"], a["leaf_area"]), "re-acclimation moved the leaf area"
