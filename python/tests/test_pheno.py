"""Smoke tests for the meds.pheno Python API (needs libmeds_plant_c built; see python/README.md).

Skips itself cleanly if the shared library hasn't been built, so `pytest` never hard-fails on a
machine that only has the Python sources.
"""
import pytest

import meds.pheno as pheno


def _lib_or_skip():
    try:
        pheno.self_test()
    except FileNotFoundError as exc:
        pytest.skip(f"libmeds_plant_c not built: {exc}")


def _drive(params, days, **env):
    """Run one strategy for `days` steps under a constant environment; return the final Out."""
    ph = pheno.Phenology(params)
    out = None
    for doy in range(1, days + 1):
        out = ph.step(doy=doy, **env)
    return out


def test_self_test_passes():
    _lib_or_skip()


def test_evergreen_holds_flush_no_shed():
    _lib_or_skip()
    p = pheno.temperate_evergreen()
    out = _drive(p, 30, temp_day=298.15)
    assert abs(out.leaf_flush_rate - p.k_flush_max) < 1e-9   # permissive flush at the max rate
    assert out.leaf_shed_rate == 0.0                          # no active shed cue -> never sheds


def test_deciduous_sheds_in_cold_flushes_in_warm():
    _lib_or_skip()
    p = pheno.temperate_deciduous()
    cold = _drive(p, 60, temp_day=265.0, soil_temp=265.0, daylength=9.0)
    assert cold.leaf_shed_rate > 0.0                          # autumn/winter cold-drop raises shed
    # a warm, long-day year drives the GDD flush up and the shed back to ~0
    warm = _drive(p, 200, temp_day=295.0, soil_temp=295.0, daylength=14.0)
    assert warm.leaf_flush_rate > 0.5 * p.k_flush_max
    assert warm.leaf_shed_rate < 0.1 * p.k_shed_max


def test_drought_deciduous_sheds_when_dry():
    _lib_or_skip()
    p = pheno.drought_deciduous()                             # tlp = -1.5 MPa
    wet = _drive(p, 30, dmax_leaf_psi=-0.4)
    dry = _drive(p, 30, dmax_leaf_psi=-3.0)
    assert wet.leaf_shed_rate < 0.05 * p.k_shed_max           # facultatively evergreen when watered
    assert dry.leaf_shed_rate > 0.5 * p.k_shed_max            # sheds under sustained drought
    assert dry.leaf_flush_rate < wet.leaf_flush_rate


def test_light_exchanging_flush_high_shed_tracks_light():
    _lib_or_skip()
    p = pheno.light_exchanging()
    dim = _drive(p, 40, rad=80.0)
    bright = _drive(p, 40, rad=520.0)
    # flush is permissive (== k_flush_max) regardless of light ...
    assert abs(dim.leaf_flush_rate - p.k_flush_max) < 1e-9
    assert abs(bright.leaf_flush_rate - p.k_flush_max) < 1e-9
    # ... while the active shed RISES with radiation (both rates > 0 under high light).
    assert bright.leaf_shed_rate > dim.leaf_shed_rate
    assert bright.leaf_shed_rate > 0.0 and bright.leaf_flush_rate > 0.0


def test_integrate_lai_bare_to_full_and_snap():
    _lib_or_skip()
    # A pure flush at k_flush_max fills a bare canopy toward full in ~1/k days.
    lai = 0.0
    for _ in range(30):
        lai = pheno.integrate_lai(lai, 1.0 / 15.0, 0.0)
    assert lai > 0.9
    # A pure shed empties it and snaps to exactly bare.
    for _ in range(40):
        lai = pheno.integrate_lai(lai, 0.0, 1.0 / 20.0)
    assert lai == 0.0
