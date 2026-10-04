# SPDX-License-Identifier: Apache-2.0
"""The parameter registry: which keys the fit may move, where each lives, its range, its prior,
its stage, and how it is transformed (MEDS_FAST_CALIBRATION_PLAN.md §4, §6.1; the revision plan
MEDS_FAST_CALIBRATION_REVISION_PLAN.md §4).

A parameter bounded to [a, b] is fitted as u = logit((g(theta) - g(a)) / (g(b) - g(a))), with g the
identity for a "linear" key and log for a "log" key (a positive scale). The fit is unconstrained in
u and every theta stays inside its range.

The prior is Gaussian in u. Its centre is the registry's `prior.centre`, else the key's default
(the base configuration's value). Its width is set by `prior.sd` (in theta's units) or
`prior.log_sd` (in log theta, for a log key), mapped into u at the centre. Without either it is
SIGMA_U, chosen so that a default at the centre of its range has the central 95 % of the range
inside +-2 sigma.

Each key has a state, the registry's recommendation:
  fit       in the default fitted set
  optional  available, fitted only when the site declaration asks for it
  fixed     not tunable in the fast fit (`reason` says why); a site may still ask for it
and a stage, the part of the staged fit that sets it (stages.py). A key that acts only through
one process (the wet canopy, night, snow, drought) names it in `process`: the fit fixes it when the
kept data sample that process in too few records (datarules.process_coverage).
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
STATES = ("fit", "optional", "fixed")
STAGES = ("optics", "photosynthesis", "energy", "water")
PROCESSES = ("", "wet_canopy", "night", "snow", "drought")
PRIOR_KEYS = ("centre", "sd", "log_sd", "source")


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
    state: str = "fit"
    stage: str = "energy"
    reason: str = ""           # why a fixed key is fixed
    process: str = ""          # the one process it acts through, if any (fixed without coverage)
    prior: dict = field(default_factory=dict)   # centre, sd or log_sd, source
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

    # ----- the prior ----------------------------------------------------------------------------
    @property
    def centre(self) -> float:
        """The prior's centre in theta: prior.centre, else the default."""
        c = self.prior.get("centre")
        return float(self.default if c is None else c)

    @property
    def u0(self) -> float:
        """The prior's centre in u (also where the fit starts)."""
        return float(self.to_u(self.centre))

    @property
    def sigma_u(self) -> float:
        """The prior's sd in u: prior.sd (theta units) or prior.log_sd (log theta) mapped into u at
        the centre, else SIGMA_U."""
        u0 = self.u0
        if self.prior.get("sd") is not None:
            return float(self.prior["sd"]) / float(self.dtheta_du(u0))
        if self.prior.get("log_sd") is not None:
            return float(self.prior["log_sd"]) * self.centre / float(self.dtheta_du(u0))
        return SIGMA_U


def _check(p: Param):
    if p.file not in FILES:
        raise ValueError(f"{p.name}: file must be one of {FILES}")
    if p.transform not in TRANSFORMS:
        raise ValueError(f"{p.name}: transform must be one of {TRANSFORMS}")
    if p.state not in STATES:
        raise ValueError(f"{p.name}: state must be one of {STATES}")
    if p.stage not in STAGES:
        raise ValueError(f"{p.name}: stage must be one of {STAGES}")
    if p.process not in PROCESSES:
        raise ValueError(f"{p.name}: process must be one of {PROCESSES[1:]}")
    if not p.lo < p.hi:
        raise ValueError(f"{p.name}: range must be increasing")
    if p.transform == "log" and p.lo <= 0.0:
        raise ValueError(f"{p.name}: a log transform needs a positive range")
    bad = set(p.prior) - set(PRIOR_KEYS)
    if bad:
        raise ValueError(f"{p.name}: unknown prior settings {sorted(bad)}; known: {PRIOR_KEYS}")
    if p.prior.get("sd") is not None and p.prior.get("log_sd") is not None:
        raise ValueError(f"{p.name}: give the prior's sd or its log_sd, not both")
    if p.prior.get("log_sd") is not None and p.transform != "log":
        raise ValueError(f"{p.name}: prior.log_sd needs a log transform")


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
                  notes=e.get("notes", ""), state=e.get("state", "fit"), stage=e.get("stage", "energy"),
                  reason=e.get("reason", ""), process=e.get("process", ""), prior=dict(e.get("prior", {})))
        _check(p)
        if p.file == "pft" and p.pft is None:
            p.pft = 1
        out.append(p)
    return out


def select(params: list[Param], fitcfg: dict, priors: dict) -> list[Param]:
    """The keys a fit moves, from the registry's menu and the site declaration's [fit] and
    [priors] tables:

      [fit].keys       the exact list (any state, a fixed key included)
      [fit].add        keys added to the registry's default set (state "fit")
      [fit].remove     keys removed from it
      [priors.<key>]   a key's prior and range: centre, sd or log_sd, source, range

    An unknown name is an error, so a misspelling cannot silently drop a key."""
    names = {p.name for p in params}
    for k in ("keys", "add", "remove"):
        unknown = [n for n in fitcfg.get(k, []) if n not in names]
        if unknown:
            raise ValueError(f"[fit].{k}: {unknown} are not in the registry (or not in this variant)")
    unknown = [n for n in priors if n not in names]
    if unknown:
        raise ValueError(f"[priors]: {unknown} are not in the registry (or not in this variant)")
    if "keys" in fitcfg and ("add" in fitcfg or "remove" in fitcfg):
        raise ValueError("[fit]: give `keys` (the exact list) or `add`/`remove`, not both")
    if "keys" in fitcfg:
        chosen = set(fitcfg["keys"])
    else:
        chosen = {p.name for p in params if p.state == "fit"} | set(fitcfg.get("add", []))
        chosen -= set(fitcfg.get("remove", []))
    out = []
    for p in params:
        if p.name not in chosen:
            continue
        over = dict(priors.get(p.name, {}))
        if "range" in over:
            p.lo, p.hi = (float(x) for x in over.pop("range"))
        p.prior = {**p.prior, **over}
        if "sd" in over:
            p.prior.pop("log_sd", None)
        if "log_sd" in over:
            p.prior.pop("sd", None)
        _check(p)
        out.append(p)
    return out


def resolve_defaults(params: list[Param], base, record: dict | None = None):
    """Give every parameter its default: the registry's, else the base configuration's (a
    meds.config.RunConfig), else the value the base run's parameter record says the model read (a
    key the base TOML leaves to its compiled default). A key the model never read is an error: a
    misspelling, or a dead key. The prior's centre must lie strictly inside the range; the default
    need not (a base value at a bound is fine: the fit starts from the prior's centre)."""
    for p in params:
        if p.default is None:
            v = base.get(p.key, file=p.file, pft=p.pft if p.file == "pft" else None)
            if v is None and record is not None:
                v = record_value(record, p)
            if v is None:
                raise KeyError(f"{p.name}: '{p.key}' is neither in the base {p.file} TOML nor in the "
                               "model's parameter record -- misspelt, or not read by this model")
            p.default = float(v)
        if not p.lo < p.centre < p.hi:
            raise ValueError(f"{p.name}: the prior's centre {p.centre} is not strictly inside "
                             f"[{p.lo}, {p.hi}]")
    return params


def record_value(record: dict, p: Param):
    """The value a parameter record (meds.config.read_record) holds for p, or None."""
    idx = p.pft if (p.file == "pft" and p.pft is not None) else 0
    hit = record.get((p.file, p.key, idx))
    if hit is None and p.file == "pft":
        hit = record.get((p.file, p.key, 0))
    return None if hit is None else hit[1]


def prior_u(params: list[Param]) -> tuple[np.ndarray, np.ndarray]:
    return np.array([p.u0 for p in params]), np.array([p.sigma_u for p in params])


def interval(p: Param, u: float, sd_u: float, z: float) -> tuple[float, float]:
    """A +-z sd interval in u mapped to theta: asymmetric in theta, and inside the range."""
    a, b = float(p.to_theta(u - z * sd_u)), float(p.to_theta(u + z * sd_u))
    return (min(a, b), max(a, b))
