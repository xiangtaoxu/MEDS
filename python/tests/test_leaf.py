"""Smoke tests for the meds.leaf Python API (needs libmeds_leaf_c built; see python/README.md).

Skips itself cleanly if the shared library hasn't been built, so `pytest` never hard-fails on a
machine that only has the Python sources.
"""
import pytest

import meds.leaf as leaf


def _lib_or_skip():
    try:
        return leaf.self_test()
    except FileNotFoundError as exc:
        pytest.skip(f"libmeds_leaf_c not built: {exc}")


def test_self_test_passes():
    flux = _lib_or_skip()
    assert flux.converged
    assert flux.a_net > 0.0
    assert 0.0 < flux.ci < flux.cs


def test_tref_identity():
    _lib_or_skip()
    # peaked/arrhenius return k25 exactly at the 25 degC reference.
    assert abs(leaf.peaked(60.0, 65330.0, 200000.0, 650.0, 298.15) - 60.0) < 1e-9
    assert abs(leaf.arrhenius(0.9, 46390.0, 298.15) - 0.9) < 1e-9


def test_gs_falls_with_vpd():
    _lib_or_skip()
    p = leaf.c3_params()
    dry = leaf.gas_exchange(par=1500.0, leaf_temp=298.15, vpd=2500.0, ca=400.0, params=p)
    wet = leaf.gas_exchange(par=1500.0, leaf_temp=298.15, vpd=500.0, ca=400.0, params=p)
    assert dry.gs < wet.gs


def test_c4_less_ci_sensitive_than_c3():
    _lib_or_skip()
    # A C4 leaf should out-assimilate a C3 leaf of the same Vcmax at low Ci / high light.
    env = dict(par=1800.0, leaf_temp=303.15, vpd=1200.0, ca=400.0)
    c3 = leaf.gas_exchange(**env, params=leaf.c3_params(vcmax25=40.0))
    c4 = leaf.gas_exchange(**env, params=leaf.c4_params(vcmax25=40.0))
    assert c3.converged and c4.converged
    assert c4.a_net > c3.a_net
