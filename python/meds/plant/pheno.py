"""meds.plant.pheno — leaf-phenology SIGNAL kernel (part of the plant-ecophysiology package).

A Pythonic front end to the MEDS Fortran phenology kernel (meds_phenology.f90, exposed through the
same libmeds_plant_c as meds.plant's gas exchange -- one C-API for the whole plant module). Given daily
environmental cues + per-PFT traits it returns two RELATIVE rate tendencies -- leaf_flush_rate and
leaf_shed_rate [1/day] -- from two governor accumulators it advances in place.

    import meds.plant.pheno as pheno

    ph = pheno.Phenology(pheno.temperate_deciduous())     # a stateful driver (holds params + memory)
    out = ph.step(temp_day=290.0, soil_temp=290.0, daylength=13.0, doy=150)
    print(out.leaf_flush_rate, out.leaf_shed_rate)        # [1/day] flush / shed TENDENCIES

    # the kernel is signal-only (no leaf mass). Track relative LAI + realized litter yourself:
    lai, litter = pheno.leaf_step(lai, out.leaf_flush_rate, out.leaf_shed_rate, baseline_turnover=8e-4)

`leaf_step` is a compact relative-unit analogue of the Fortran carbon leaf update
(meds_plant_carbon_dynamics): the flush TENDENCY fills toward full, the shed TENDENCY + a small
baseline turnover remove leaves, and it reports the REALIZED litter (leaf actually shed that step),
which is not the same as the shed tendency -- e.g. a bare deciduous canopy has a high winter shed
tendency but zero realized litter. Requires the compiled libmeds_plant_c (see meds.plant._ffi).
"""
from __future__ import annotations

import ctypes
from ctypes import c_double, c_int, byref, POINTER
from dataclasses import dataclass, asdict
from enum import IntFlag

from ._ffi import _lib   # the shared libmeds_plant_c loader (also used by leaf gas exchange)

__all__ = [
    "Cue", "Params", "State", "Out", "Phenology", "step", "leaf_step", "integrate_lai",
    "temperate_deciduous", "temperate_evergreen", "drought_deciduous", "light_exchanging",
    "self_test",
]


#===========================================================================================#
#  ctypes mirrors — field order MUST match the bind(c) types in src/plant/meds_plant_capi.f90.#
#===========================================================================================#
_ENV_FIELDS = [
    ("temp_day", c_double), ("soil_temp", c_double), ("avail_water", c_double),
    ("dmax_leaf_psi", c_double), ("rad", c_double), ("daylength", c_double),
    ("doy", c_int), ("hemis_north", c_int),
]
_PARAM_FIELDS = [
    ("flush_cue_mask", c_int), ("shed_cue_mask", c_int),
    ("cue_sharpness", c_double), ("k_flush_max", c_double), ("k_shed_max", c_double),
    ("tau_flush", c_double), ("tau_shed", c_double),
    ("gdd_base_temp", c_double), ("chill_base_temp", c_double),
    ("phen_a", c_double), ("phen_b", c_double), ("phen_c", c_double),
    ("cold_drop_daylength", c_double), ("cold_drop_soiltemp1", c_double),
    ("cold_drop_soiltemp2", c_double),
    ("water_use_potential", c_int),
    ("water_off_threshold", c_double), ("water_on_threshold", c_double),
    ("water_window", c_double), ("water_width", c_double),
    ("leaf_psi_tlp", c_double), ("low_psi_threshold", c_double), ("high_psi_threshold", c_double),
    ("photo_crit", c_double), ("photo_slope", c_double),
    ("light_on_threshold", c_double), ("light_width", c_double), ("light_window", c_double),
]
_STATE_FIELDS = [
    ("flush_drive", c_double), ("shed_drive", c_double), ("gdd", c_double), ("chill", c_double),
    ("water_avg", c_double), ("low_psi_days", c_double), ("high_psi_days", c_double),
    ("light_avg", c_double),
]
_OUT_FIELDS = [("leaf_flush_rate", c_double), ("leaf_shed_rate", c_double), ("cue_limiting", c_int)]
_INT_NAMES = {"doy", "hemis_north", "flush_cue_mask", "shed_cue_mask", "water_use_potential"}


class _EnvC(ctypes.Structure):
    _fields_ = _ENV_FIELDS


class _ParamsC(ctypes.Structure):
    _fields_ = _PARAM_FIELDS


class _StateC(ctypes.Structure):
    _fields_ = _STATE_FIELDS


class _OutC(ctypes.Structure):
    _fields_ = _OUT_FIELDS


_ENV_NAMES = tuple(n for n, _ in _ENV_FIELDS)
_PARAM_NAMES = tuple(n for n, _ in _PARAM_FIELDS)
_STATE_NAMES = tuple(n for n, _ in _STATE_FIELDS)
_PHENO_BOUND = False


def _pheno_lib():
    """The shared libmeds_plant_c (loaded by meds.plant._ffi), with meds_phenology_step bound once."""
    global _PHENO_BOUND
    lib = _lib()
    if not _PHENO_BOUND:
        lib.meds_phenology_step.restype = None
        lib.meds_phenology_step.argtypes = [POINTER(_EnvC), POINTER(_ParamsC), c_double,
                                            POINTER(_StateC), POINTER(_OutC)]
        _PHENO_BOUND = True
    return lib


def _make(struct_cls, names, src):
    return struct_cls(**{n: (int(src[n]) if n in _INT_NAMES else float(src[n])) for n in names})


class Cue(IntFlag):
    """Cue-enable bits (mirror meds_plant_types). Flush and shed sides are selected independently."""
    NONE = 0
    TEMP = 1     # temperature: GDD flush + autumn cold-drop shed
    WATER = 2    # soil-water running mean       (driver not wired in the standalone model yet)
    HYDRO = 4    # daily-max leaf water potential (dmax_leaf_psi)
    PHOTO = 8    # photoperiod (gates the temperature flush)
    LIGHT = 16   # radiation: active shed rises with running-mean light


@dataclass
class Params:
    """Per-PFT phenology traits. Rates are RELATIVE, per day; the two masks pick which cues drive
    each side. k_flush_max = 1/15 d => a bare canopy fills in ~15 days at full drive."""
    flush_cue_mask: int = Cue.NONE
    shed_cue_mask: int = Cue.NONE
    cue_sharpness: float = 2.0
    k_flush_max: float = 1.0 / 15.0
    k_shed_max: float = 1.0 / 20.0
    tau_flush: float = 5.0
    tau_shed: float = 5.0
    gdd_base_temp: float = 278.15
    chill_base_temp: float = 278.15
    phen_a: float = -68.0
    phen_b: float = 638.0
    phen_c: float = -0.01
    cold_drop_daylength: float = 10.9
    cold_drop_soiltemp1: float = 284.3
    cold_drop_soiltemp2: float = 275.15
    water_use_potential: bool = False
    water_off_threshold: float = 0.2
    water_on_threshold: float = 0.5
    water_window: float = 10.0
    water_width: float = 0.1
    leaf_psi_tlp: float = -2.0
    low_psi_threshold: float = 10.0
    high_psi_threshold: float = 10.0
    photo_crit: float = 11.0
    photo_slope: float = 2.0
    light_on_threshold: float = 200.0
    light_width: float = 50.0
    light_window: float = 10.0


@dataclass
class State:
    """The prognostic phenological memory (two governor drives + cue sub-accumulators). Born
    flushing (flush_drive=1) with no active shed — the evergreen fixed point."""
    flush_drive: float = 1.0
    shed_drive: float = 0.0
    gdd: float = 0.0
    chill: float = 0.0
    water_avg: float = 0.0
    low_psi_days: float = 0.0
    high_psi_days: float = 0.0
    light_avg: float = 0.0


@dataclass(frozen=True)
class Out:
    """One step's phenology signal (TENDENCIES; the realized fluxes come from leaf_step)."""
    leaf_flush_rate: float    # [1/day] relative flush tendency (0 => dormant)
    leaf_shed_rate: float     # [1/day] relative active-shed tendency (0 => no active shed)
    cue_limiting: int         # the strongest active shed cue (a Cue bit; diagnostic)


_ENV_DEFAULTS = dict(temp_day=298.15, soil_temp=298.15, avail_water=0.5, dmax_leaf_psi=0.0,
                     rad=400.0, daylength=12.0, doy=1, hemis_north=True)


def step(env, params, state, dt=1.0):
    """Low-level one-step call. `env` is a dict (missing keys default), `params`/`state` are Params/
    State (or dicts). Returns (Out, new_State); `state` is NOT mutated -- feed new_State back."""
    env_full = {**_ENV_DEFAULTS, **dict(env)}
    p = asdict(params) if isinstance(params, Params) else dict(params)
    s = asdict(state) if isinstance(state, State) else dict(state)
    env_c = _make(_EnvC, _ENV_NAMES, env_full)
    params_c = _make(_ParamsC, _PARAM_NAMES, p)
    state_c = _make(_StateC, _STATE_NAMES, s)
    out_c = _OutC()
    _pheno_lib().meds_phenology_step(byref(env_c), byref(params_c), float(dt),
                                     byref(state_c), byref(out_c))
    out = Out(leaf_flush_rate=out_c.leaf_flush_rate, leaf_shed_rate=out_c.leaf_shed_rate,
              cue_limiting=int(out_c.cue_limiting))
    new_state = State(**{n: getattr(state_c, n) for n in _STATE_NAMES})
    return out, new_state


class Phenology:
    """A stateful phenology driver: holds the per-PFT params + the advancing memory."""

    def __init__(self, params: Params, state: State | None = None):
        self.params = params
        self.state = state if state is not None else State()

    def step(self, dt: float = 1.0, **env) -> Out:
        out, self.state = step(env, self.params, self.state, dt=dt)
        return out


def leaf_step(elongf, leaf_flush_rate, leaf_shed_rate, dt=1.0,
              baseline_turnover=0.0, elongf_min=0.02):
    """Advance relative LAI (canopy fullness, elongf in [0,1]) ONE step; return (new_elongf, litter).

    A compact relative-unit analogue of the Fortran carbon leaf update (get_plant_flux_slow):
      * flush FILLS toward full at leaf_flush_rate (linear, capped by the deficit 1-elongf),
      * the ACTIVE shed removes leaf_shed_rate per day (linear toward bare), NON-replaceable,
      * `baseline_turnover` [1/day] is a small proportional, REPLACEABLE background loss,
      * a net decline that crosses below elongf_min SNAPS to bare (ED2's fully-abscised state).
    `litter` is the REALIZED relative leaf carbon shed to litter this step (= flush-adjusted leaf
    removed), which differs from the shed TENDENCY: a bare canopy sheds nothing however high the
    tendency, and a full evergreen canopy litters via baseline turnover with zero shed tendency.
    """
    flush_gain = min(leaf_flush_rate * dt, max(0.0, 1.0 - elongf))
    loss = min(leaf_shed_rate * dt + baseline_turnover * elongf * dt, elongf)
    e = elongf + flush_gain - loss
    if leaf_shed_rate > 0.0 and loss > flush_gain and e < elongf_min:
        e = 0.0
    e = min(1.0, max(0.0, e))
    litter = max(0.0, (elongf + flush_gain) - e)      # leaf present (incl. this step's flush) that left
    return e, litter


def integrate_lai(elongf, leaf_flush_rate, leaf_shed_rate, dt=1.0,
                  baseline_turnover=0.0, elongf_min=0.02):
    """Relative LAI after one step (canopy fullness in [0,1]); the scalar half of `leaf_step`."""
    e, _ = leaf_step(elongf, leaf_flush_rate, leaf_shed_rate, dt, baseline_turnover, elongf_min)
    return e


#===========================================================================================#
#  The four target phenological strategies (design §1a / §6.1), as ready-to-use Params.       #
#===========================================================================================#
def _preset(defaults, overrides) -> Params:
    return Params(**{**defaults, **overrides})


def temperate_deciduous(**overrides) -> Params:
    """Cold-deciduous: flush on spring GDD, shed on the autumn cold-drop (both masks = TEMP)."""
    return _preset(dict(flush_cue_mask=Cue.TEMP, shed_cue_mask=Cue.TEMP,
                        k_flush_max=1.0 / 15.0, k_shed_max=1.0 / 18.0), overrides)


def temperate_evergreen(**overrides) -> Params:
    """Evergreen: permissive flush (no onset cue), no active shed — the canopy is held full year
    round, thinned only by the small baseline turnover the caller applies in leaf_step. Pass
    flush_cue_mask=Cue.TEMP for a mild seasonal flush that still never actively sheds."""
    return _preset(dict(flush_cue_mask=Cue.NONE, shed_cue_mask=Cue.NONE, k_flush_max=1.0 / 15.0),
                   overrides)


def drought_deciduous(**overrides) -> Params:
    """Facultative drought-deciduous (tropical): flush + shed keyed on the daily-max leaf water
    potential vs the turgor-loss point (both masks = HYDRO). Evergreen when never droughted."""
    return _preset(dict(flush_cue_mask=Cue.HYDRO, shed_cue_mask=Cue.HYDRO,
                        leaf_psi_tlp=-1.5, low_psi_threshold=10.0, high_psi_threshold=10.0,
                        k_flush_max=1.0 / 15.0, k_shed_max=1.0 / 20.0), overrides)


def light_exchanging(**overrides) -> Params:
    """Light-driven leaf-exchanging (tropical evergreen): permissive high flush + an active shed that
    RISES with running-mean radiation, so the canopy stays ~full while turning leaves over fast."""
    return _preset(dict(flush_cue_mask=Cue.NONE, shed_cue_mask=Cue.LIGHT,
                        light_on_threshold=280.0, light_width=70.0, light_window=10.0,
                        k_flush_max=1.0 / 12.0, k_shed_max=1.0 / 25.0), overrides)


def self_test() -> None:
    """Smoke test: evergreen holds flush=k_flush_max / shed=0; a cold spell drives a deciduous shed
    up; a full canopy under a hard shed goes bare and its litter matches the loss."""
    ev = Phenology(temperate_evergreen())
    o = ev.step(temp_day=298.15, doy=180)
    assert abs(o.leaf_flush_rate - ev.params.k_flush_max) < 1e-9, "evergreen flush should be k_flush_max"
    assert o.leaf_shed_rate == 0.0, "evergreen should have no active shed"

    dec = Phenology(temperate_deciduous())
    for doy in range(1, 60):
        out = dec.step(temp_day=265.0, soil_temp=265.0, daylength=9.0, doy=doy)
    assert out.leaf_shed_rate > 0.0, "cold winter should raise the deciduous shed rate"

    lai, total_litter = 1.0, 0.0
    for _ in range(40):
        lai, lit = leaf_step(lai, 0.0, dec.params.k_shed_max)
        total_litter += lit
    assert lai < 0.05, f"canopy should shed to ~bare, got {lai}"
    assert abs(total_litter - 1.0) < 0.05, f"a full canopy shed to bare should litter ~1.0, got {total_litter}"
    print("meds.plant.pheno.self_test: OK")
