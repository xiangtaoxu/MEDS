# SPDX-License-Identifier: Apache-2.0
"""meds.canopy — the canopy's kernels over a frozen stand, parameterized by a run's own config.

The leaf solve over many leaves and the two-stream radiation over a stand, with every parameter
built from a run's main TOML (and the PFT file it names) by the model's own code, so a caller
outside the model computes what the fast loop computes:

    from meds.canopy import Canopy
    c = Canopy("run/main.toml")
    vc, rd = c.plastic_traits(pft, lai_above)          # each cohort's plastic Vcmax25, Rd25
    a = c.leaf(pft, vc, rd, par=..., leaf_temp=..., vpd=..., ca=..., pressure=...,
               psi_leaf=..., gb=..., psi=...)           # dict of arrays: A_gross, A_net, gs, ci
    r = c.radiation(cosz, vis_beam, vis_diff, nir_beam, nir_diff,
                    pft, height, lai, wai)              # one patch, many hours

The fast calibration's optics and photosynthesis stages use it (scripts/calibrate_fast). The leaf
drivers are those the model reports per cohort as the `gx_*_cohort_fast` outputs; the stand's
geometry is its restart state plus the `wai_cohort` output.

The library holds ONE configuration at a time. Each Canopy re-opens its own before a call when
another has been opened since, so several can be used in turn (not from threads). A configuration
with thermal acclimation on is refused: its leaf table depends on model state.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import _ffi

__all__ = ["Canopy"]

_OPEN = {"path": None}


class Canopy:
    def __init__(self, config):
        self.path = str(Path(config).resolve())
        self._open()

    def _open(self):
        if _OPEN["path"] == self.path:
            return
        status = _ffi.open_config(self.path)
        if status == -1:
            _OPEN["path"] = None
            raise ValueError(f"{self.path}: thermal acclimation is on; its leaf table depends on model "
                             "state, so meds.canopy cannot reproduce it")
        if status != 0:
            _OPEN["path"] = None
            raise RuntimeError(f"{self.path}: meds_canopy_open returned {status}")
        _OPEN["path"] = self.path

    @property
    def n_pft(self) -> int:
        self._open()
        return _ffi.n_pft()

    def plastic_traits(self, pft, lai_above):
        """Each cohort's Vcmax25 and Rd25 from its PFT (1-based) and the leaf area above it, as a
        restart re-acclimates them; the PFT's own values when trait plasticity is off."""
        self._open()
        return _ffi.plastic(pft, lai_above)

    def leaf(self, pft, vcmax25, rd25, *, par, leaf_temp, vpd, ca, pressure, psi_leaf=0.0, gb=0.0,
             psi=0.0) -> dict:
        """The coupled leaf solve for every leaf (arrays broadcast): par [umol/m2 leaf/s, incident-
        equivalent], leaf_temp [K], vpd [Pa], ca [umol/mol], pressure [Pa], psi_leaf and psi [MPa],
        gb [mol/m2/s]. Returns A_gross, A_net [umol/m2 leaf/s], gs [mol/m2/s] and ci [umol/mol]."""
        self._open()
        return _ffi.leaf(pft, vcmax25, rd25, par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi)

    def radiation(self, cosz, vis_beam, vis_diff, nir_beam, nir_diff, pft, height, lai, wai) -> dict:
        """The two-stream for one patch's cohorts (any order) over the given hours of shortwave
        [W/m2]: up_vis, up_nir (at the canopy top), leaf_vis, wood_vis (absorbed), ground_vis (down at
        the ground), each per hour, and cohort_leaf_vis (hours x cohorts). The ground is the
        configuration's bare soil."""
        self._open()
        return _ffi.radiation(cosz, vis_beam, vis_diff, nir_beam, nir_diff, np.asarray(pft), height, lai, wai)
