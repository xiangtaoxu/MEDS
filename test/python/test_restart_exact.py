# SPDX-License-Identifier: Apache-2.0
"""A restart continues exactly as the run that wrote the state would have, and a restart that
re-acclimates its traits ([init].reacclimate_traits) gets the traits a census start gives.

Each case runs meds_main from the demography example's census on the fast loop's constant
reference climate (no forcing file) and compares the final state files variable by variable,
bit for bit. MEDS_MAIN names the executable; CTest sets it.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from meds.config import RunConfig  # noqa: E402

MEDS_MAIN = os.environ.get("MEDS_MAIN", str(ROOT / "build-ifx" / "meds_main"))


def base_config(run_dir, start, end, *, slow_on, restart_file=None, reacclimate=False,
                plasticity=True, canopy_water=False):
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    for key, value in {"run.start_time": start, "run.end_time": end, "run.slow_on": slow_on,
                       "run.n_threads": 1, "fast.fast_biophysics_on": True,
                       "fast.canopy_water_on": canopy_water, "forcing.forcing_on": False,
                       "init.init_mode": 2 if restart_file else 1,
                       "init.restart_file": restart_file or "none",
                       "init.census_file": str(ROOT / "data/census_example.csv"),
                       "state.output_dir": str(run_dir / "out"), "state.output_prefix": "s",
                       "state.write_state": True, "output.enabled": False,
                       "trait_dynamics.trait_plasticity_on": plasticity}.items():
        cfg.set(key, value)
    if reacclimate:
        cfg.set("init.reacclimate_traits", True)
    return cfg


def run(run_dir, cfg, pft_edit=None):
    if pft_edit:
        pft_edit(cfg)
    res = subprocess.run([MEDS_MAIN, str(cfg.write(run_dir))], cwd=run_dir, capture_output=True,
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
    """A change to the PFT file: every PFT's vcmax25 times `factor`."""
    def edit(cfg):
        cfg.set("pft.vcmax25", [v * factor for v in cfg.get("pft.vcmax25", file="pft")], file="pft")
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
    # a restart whose PFT file gives the traits the state already has changes nothing
    same = run(tmp_path / "s", base_config(tmp_path / "s", "2001-06-02", "2001-06-03", slow_on=False,
                                             plasticity=plasticity, restart_file=str(written), reacclimate=True))
    s, _ = read_state(same)
    for v in ("sla", "vcmax25", "rd25", "llspan", "leaf_area", "leaf_carbon", "nonstructural_carbon"):
        assert np.array_equal(s[v], a[v]), v
