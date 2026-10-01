# SPDX-License-Identifier: Apache-2.0
"""The parameter registry: which keys the fit may move, where each lives, its range, and how it is
transformed (MEDS_FAST_CALIBRATION_PLAN.md §4, §6.1).

A parameter bounded to [a, b] is fitted as u = logit((g(theta) - g(a)) / (g(b) - g(a))), with g the
identity for a "linear" key and log for a "log" key (a positive scale). The fit is unconstrained in
u and every theta stays inside its range. The prior is Gaussian in u, centred on the default, with
SIGMA_U chosen so that a default at the centre of its range has the central 95 % of the range
inside +-2 sigma.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from meds.config import load_toml

SIGMA_U = math.log(0.975 / 0.025) / 2.0     # 1.832: logit(0.975) / 2
H_U = 0.04                                   # the central-difference step in u (~1 % of the range)
TRANSFORMS = ("linear", "log")
FILES = ("pft", "main")


@dataclass
class Param:
    name: str                  # registry id
    file: str                  # "pft" (the PFT TOML) or "main"
    key: str                   # dotted TOML key, e.g. "pft.stomatal_g1" or "aerodynamics.z0m_ratio"
    lo: float
    hi: float
    transform: str = "linear"
    pft: int | None = None     # 1-based element of a per-PFT array (file = "pft")
    default: float | None = None
    group: str = ""
    source: str = ""
    variants: tuple = ()       # empty: every variant; else only these
    notes: str = ""
    extra: dict = field(default_factory=dict)

    # ----- the transform ------------------------------------------------------------------------
    def _g(self, x):
        return np.log(x) if self.transform == "log" else np.asarray(x, dtype=float)

    def _ginv(self, y):
        return np.exp(y) if self.transform == "log" else y

    def to_u(self, theta):
        ga, gb = self._g(self.lo), self._g(self.hi)
        p = (self._g(theta) - ga) / (gb - ga)
        if np.any(p <= 0.0) or np.any(p >= 1.0):
            raise ValueError(f"{self.name}: {theta} is not strictly inside [{self.lo}, {self.hi}]")
        return np.log(p / (1.0 - p))

    def to_theta(self, u):
        ga, gb = self._g(self.lo), self._g(self.hi)
        p = 1.0 / (1.0 + np.exp(-np.asarray(u, dtype=float)))
        return self._ginv(ga + p * (gb - ga))

    def dtheta_du(self, u):
        """d theta / d u at u (the Jacobian of the transform, for covariances in theta)."""
        ga, gb = self._g(self.lo), self._g(self.hi)
        p = 1.0 / (1.0 + np.exp(-np.asarray(u, dtype=float)))
        dg = p * (1.0 - p) * (gb - ga)
        return dg * self.to_theta(u) if self.transform == "log" else dg

    @property
    def u0(self) -> float:
        return float(self.to_u(self.default))


def load_registry(path, variant: str | None = None) -> list[Param]:
    """Read a registry TOML: one table per parameter. Entries restricted to other variants are
    dropped."""
    raw = load_toml(path)
    out = []
    for name, e in raw.items():
        if not isinstance(e, dict):
            continue
        variants = tuple(e.get("variants", ()))
        if variants and variant not in variants:
            continue
        lo, hi = (float(x) for x in e["range"])
        p = Param(name=name, file=e["file"], key=e["key"], lo=lo, hi=hi,
                  transform=e.get("transform", "linear"), pft=e.get("pft"),
                  default=None if e.get("default") is None else float(e["default"]),
                  group=e.get("group", ""), source=e.get("source", ""), variants=variants,
                  notes=e.get("notes", ""))
        if p.file not in FILES:
            raise ValueError(f"{name}: file must be one of {FILES}")
        if p.transform not in TRANSFORMS:
            raise ValueError(f"{name}: transform must be one of {TRANSFORMS}")
        if not lo < hi:
            raise ValueError(f"{name}: range must be increasing")
        if p.transform == "log" and lo <= 0.0:
            raise ValueError(f"{name}: a log transform needs a positive range")
        if p.file == "pft" and p.pft is None:
            p.pft = 1
        out.append(p)
    return out


def resolve_defaults(params: list[Param], base, record: dict | None = None):
    """Give every parameter its default: the registry's, else the base configuration's (a
    meds.config.RunConfig), else the value the base run's parameter record says the model read (a
    key the base TOML leaves to its compiled default). A key the model never read is an error: a
    misspelling, or a dead key."""
    for p in params:
        if p.default is not None:
            continue
        v = base.get(p.key, file=p.file, pft=p.pft if p.file == "pft" else None)
        if v is None and record is not None:
            v = record_value(record, p)
        if v is None:
            raise KeyError(f"{p.name}: '{p.key}' is neither in the base {p.file} TOML nor in the "
                           "model's parameter record -- misspelt, or not read by this model")
        p.default = float(v)
        if not p.lo < p.default < p.hi:
            raise ValueError(f"{p.name}: default {p.default} is not strictly inside [{p.lo}, {p.hi}]")
    return params


def record_value(record: dict, p: Param):
    """The value a parameter record (meds.config.read_record) holds for p, or None."""
    idx = p.pft if (p.file == "pft" and p.pft is not None) else 0
    hit = record.get((p.file, p.key, idx))
    if hit is None and p.file == "pft":
        hit = record.get((p.file, p.key, 0))
    return None if hit is None else hit[1]


def prior_u(params: list[Param]) -> tuple[np.ndarray, np.ndarray]:
    return np.array([p.u0 for p in params]), np.full(len(params), SIGMA_U)
