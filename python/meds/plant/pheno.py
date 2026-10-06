# SPDX-License-Identifier: Apache-2.0
"""meds.plant.pheno — the leaf-phenology kernel (part of the plant-ecophysiology package).

A Pythonic front end to the MEDS Fortran phenology kernel (meds_phenology.f90, exposed through the
same libmeds.so as meds.plant's gas exchange -- one C-API for the whole model). Each day it takes
the cues (air temperature, day length, PAR, predawn leaf water potential) and per-PFT traits, and
advances two smoothed tendencies in [0,1]: leaf_flush_tendency and leaf_shed_tendency.

    import meds.plant.pheno as pheno

    ph = pheno.Phenology(pheno.temperate_deciduous())     # params + cue memory + leaf cover
    day = ph.step(temp_day=290.0, daylength=pheno.daylength(42.5, 150), doy=150)
    print(day.leaf_cover, day.senescence)                 # canopy fullness, leaf lost today

`Phenology.step` applies the coupled model's own leaf rule (`leaf_step`, the Fortran
leaf_turnover_step through the C-API) with carbon never limiting the flush: senescence at
shed_rate_max * leaf_shed_tendency down to min_leaf_cover, background turnover while the canopy
flushes, and a flush of up to flush_rate_max * leaf_flush_tendency per day toward full cover.
Requires the compiled libmeds.so (see meds._libmeds).
"""
from __future__ import annotations

import ctypes
from ctypes import c_double, c_int, byref, POINTER
from dataclasses import dataclass, asdict, fields
from enum import IntFlag

from ._ffi import _lib   # the shared libmeds.so handle (also used by leaf gas exchange)

__all__ = [
    "Cue", "Params", "State", "Out", "Day", "Phenology", "step", "leaf_step", "daylength",
    "temperate_deciduous", "boreal_evergreen", "drought_deciduous", "light_exchanging",
    "self_test",
]


#===========================================================================================#
#  ctypes mirrors — field order MUST match the bind(c) types in                             #
#  src/c_api/meds_c_api_phenology.f90.                                                      #
#===========================================================================================#
_ENV_FIELDS = [
    ("temp_day", c_double), ("daylength", c_double), ("par", c_double),
    ("predawn_leaf_psi", c_double), ("doy", c_int), ("hemis_north", c_int),
]
_PARAM_FIELDS = [
    ("flush_cue_mask", c_int), ("shed_cue_mask", c_int),
    ("flush_cue_timescale", c_double), ("shed_cue_timescale", c_double),
    ("flush_rate_max", c_double), ("shed_rate_max", c_double),
    ("flush_base_temp", c_double), ("flush_degree_days", c_double), ("flush_temp_sharpness", c_double),
    ("shed_base_temp", c_double), ("shed_degree_days", c_double), ("shed_temp_sharpness", c_double),
    ("flush_daylength_threshold", c_double), ("flush_daylength_sharpness", c_double),
    ("shed_daylength_threshold", c_double), ("shed_daylength_sharpness", c_double),
    ("flush_par_threshold", c_double), ("flush_par_sharpness", c_double),
    ("shed_par_threshold", c_double), ("shed_par_sharpness", c_double), ("par_window", c_double),
    ("leaf_psi_tlp", c_double), ("flush_water_sum", c_double), ("flush_water_sharpness", c_double),
    ("shed_water_sum", c_double), ("shed_water_sharpness", c_double),
]
_STATE_FIELDS = [
    ("leaf_flush_tendency", c_double), ("leaf_shed_tendency", c_double),
    ("growing_degree_days", c_double), ("cold_degree_days", c_double),
    ("wet_psi_sum", c_double), ("dry_psi_sum", c_double), ("par_mean", c_double),
]
_OUT_FIELDS = [("leaf_flush_potential", c_double), ("leaf_shed_potential", c_double)]


class _EnvC(ctypes.Structure):
    _fields_ = _ENV_FIELDS


class _ParamsC(ctypes.Structure):
    _fields_ = _PARAM_FIELDS


class _StateC(ctypes.Structure):
    _fields_ = _STATE_FIELDS


class _OutC(ctypes.Structure):
    _fields_ = _OUT_FIELDS


_BOUND = False


def _pheno_lib():
    """The shared libmeds.so (via meds.plant._ffi), with the phenology entries bound once."""
    global _BOUND
    lib = _lib()
    if not _BOUND:
        lib.meds_phenology_step.restype = None
        lib.meds_phenology_step.argtypes = [POINTER(_EnvC), POINTER(_ParamsC), c_double,
                                            POINTER(_StateC), POINTER(_OutC)]
        lib.meds_leaf_turnover_step.restype = None
        lib.meds_leaf_turnover_step.argtypes = [c_double] * 10 + [POINTER(c_double)] * 3
        lib.meds_daylength.restype = c_double
        lib.meds_daylength.argtypes = [c_double, c_int]
        _BOUND = True
    return lib


def _make(struct_cls, src):
    return struct_cls(**{n: (int(src[n]) if t is c_int else float(src[n]))
                         for n, t in struct_cls._fields_})


class Cue(IntFlag):
    """Cue-enable bits (mirror meds_phenology_types). flush_cue_mask and shed_cue_mask pick the
    cues of each side: the flush signal is the PRODUCT of its cues' switches, the shed signal the
    larger of the seasonal trigger (TEMP x DAYLENGTH x PAR) and the water trigger."""
    NONE = 0
    TEMP = 1        # warmth sum from midwinter (flush), cold sum from midsummer (shed)
    DAYLENGTH = 2   # photoperiod [h]
    WATER = 4       # predawn leaf water potential summed above / below the turgor-loss point
    PAR = 8         # running-mean PAR at the cohort's top [umol/m2/s]


@dataclass
class Params:
    """Per-PFT phenology traits. Every switch is sigma(s (x - x*)): a centre x* and a signed
    sharpness s; its 12-88 % transition spans x* +- 2/|s|. The defaults (no cues: always flushing,
    never senescing) are those of a PFT without a [phenology] section.

    The last three fields belong to the carbon layer (leaf_step), not the cue kernel: the
    background turnover while flushing, the leaf cover senescence stops at (an evergreen floor),
    and the cover below which a dormant canopy goes bare."""
    flush_cue_mask: int = Cue.NONE
    shed_cue_mask: int = Cue.NONE
    flush_cue_timescale: float = 5.0          # [day]
    shed_cue_timescale: float = 5.0           # [day]
    flush_rate_max: float = 0.06667           # [1/day] leaf growth at full tendency
    shed_rate_max: float = 0.05               # [1/day] senescence at full tendency
    flush_base_temp: float = 278.15           # [K]
    flush_degree_days: float = 100.0          # [K day] warmth requirement
    flush_temp_sharpness: float = 0.04        # [1/(K day)]
    shed_base_temp: float = 290.15            # [K]
    shed_degree_days: float = 50.0            # [K day] cold requirement
    shed_temp_sharpness: float = 0.1          # [1/(K day)]
    flush_daylength_threshold: float = 12.0   # [h]
    flush_daylength_sharpness: float = 1.0    # [1/h]; > 0: long days permit flushing
    shed_daylength_threshold: float = 11.0    # [h]
    shed_daylength_sharpness: float = -1.0    # [1/h]; < 0: short days trigger senescence
    flush_par_threshold: float = 300.0        # [umol/m2/s]
    flush_par_sharpness: float = 0.02         # [m2 s/umol]; > 0: bright light permits flushing
    shed_par_threshold: float = 300.0         # [umol/m2/s]
    shed_par_sharpness: float = 0.02          # > 0: bright light triggers senescence; < 0: dim does
    par_window: float = 10.0                  # [day] PAR running mean
    leaf_psi_tlp: float = -2.0                # [MPa] turgor-loss point
    flush_water_sum: float = 10.0             # [MPa day]
    flush_water_sharpness: float = 0.5        # [1/(MPa day)]
    shed_water_sum: float = 10.0              # [MPa day]
    shed_water_sharpness: float = 0.5         # [1/(MPa day)]
    leaf_turnover_rate: float = 0.0           # [1/yr] background loss while flushing
    min_leaf_cover: float = 0.0               # [-] senescence stops here
    bare_leaf_cover: float = 0.02             # [-] a dormant canopy below this goes bare


@dataclass
class State:
    """The prognostic memory: two smoothed tendencies + the cue accumulators. Born flushing and
    not senescing."""
    leaf_flush_tendency: float = 1.0
    leaf_shed_tendency: float = 0.0
    growing_degree_days: float = 0.0
    cold_degree_days: float = 0.0
    wet_psi_sum: float = 0.0
    dry_psi_sum: float = 0.0
    par_mean: float = 0.0


@dataclass(frozen=True)
class Out:
    """One step's potential relative rates (rate_max * tendency)."""
    leaf_flush_potential: float    # [1/day]
    leaf_shed_potential: float     # [1/day]


@dataclass(frozen=True)
class Day:
    """One day of a Phenology run: the tendencies and what the leaf rule made of them."""
    leaf_flush_tendency: float     # [-]
    leaf_shed_tendency: float      # [-]
    leaf_cover: float              # [-] leaf / full leaf, after today's loss and flush
    senescence: float              # [-] leaf cover lost to senescence today
    background: float              # [-] leaf cover lost to background turnover today


_ENV_DEFAULTS = dict(temp_day=298.15, daylength=12.0, par=800.0, predawn_leaf_psi=0.0,
                     doy=1, hemis_north=True)


def step(env, params, state, dt=1.0):
    """Low-level one-step call. `env` is a dict (missing keys default), `params`/`state` are Params/
    State (or dicts). Returns (Out, new_State); `state` is NOT mutated -- feed new_State back."""
    p = asdict(params) if isinstance(params, Params) else dict(params)
    s = asdict(state) if isinstance(state, State) else dict(state)
    env_c = _make(_EnvC, {**_ENV_DEFAULTS, **dict(env)})
    params_c = _make(_ParamsC, p)
    state_c = _make(_StateC, s)
    out_c = _OutC()
    _pheno_lib().meds_phenology_step(byref(env_c), byref(params_c), float(dt),
                                     byref(state_c), byref(out_c))
    out = Out(leaf_flush_potential=out_c.leaf_flush_potential,
              leaf_shed_potential=out_c.leaf_shed_potential)
    return out, State(**{n: getattr(state_c, n) for n, _ in _STATE_FIELDS})


def leaf_step(leaf_cover, leaf_flush_tendency, leaf_shed_tendency, params, dt=1.0):
    """Advance leaf cover (leaf / full leaf) ONE step with the coupled model's leaf rule.

    Returns (new_leaf_cover, senescence, background). The loss terms come from the Fortran
    leaf_turnover_step (the carbon layer calls the same routine); the flush then refills the
    post-loss canopy toward full by up to flush_rate_max * leaf_flush_tendency * dt, with carbon
    never limiting. Litter is senescence * (1 - resorption) + background.
    """
    sen, bg, cap = c_double(), c_double(), c_double()
    _pheno_lib().meds_leaf_turnover_step(
        float(leaf_cover), 1.0, float(leaf_flush_tendency), float(leaf_shed_tendency),
        float(params.flush_rate_max), float(params.shed_rate_max), float(params.leaf_turnover_rate),
        float(params.min_leaf_cover), float(params.bare_leaf_cover), float(dt),
        byref(sen), byref(bg), byref(cap))
    post = leaf_cover - sen.value - bg.value
    return post + min(max(0.0, 1.0 - post), cap.value), sen.value, bg.value


def daylength(lat_deg, doy):
    """Day length [h] at latitude `lat_deg` on day of year `doy` (meds_time::daylength)."""
    return _pheno_lib().meds_daylength(float(lat_deg), int(doy))


class Phenology:
    """A stateful driver: per-PFT params, the cue memory, and the leaf cover.

    The C structs are built once and kept, so a long daily run costs two library calls a day.
    `params` is fixed at construction; make a new Phenology to change it."""

    def __init__(self, params: Params, state: State | None = None, leaf_cover: float = 1.0):
        self.params = params
        self.leaf_cover = leaf_cover
        self._lib = _pheno_lib()
        self._params_c = _make(_ParamsC, asdict(params))
        self._state_c = _make(_StateC, asdict(state if state is not None else State()))
        self._env_c, self._out_c = _EnvC(), _OutC()
        self._sen, self._bg, self._cap = c_double(), c_double(), c_double()

    @property
    def state(self) -> State:
        """The cue memory now (a copy)."""
        return State(**{n: getattr(self._state_c, n) for n, _ in _STATE_FIELDS})

    def step(self, dt: float = 1.0, **env) -> Day:
        """Advance one step: the cues (env keys as in pheno_env_t; missing keys take their
        defaults), then the leaf rule."""
        for name, value in {**_ENV_DEFAULTS, **env}.items():
            setattr(self._env_c, name, value)
        st, p = self._state_c, self.params
        self._lib.meds_phenology_step(byref(self._env_c), byref(self._params_c), float(dt),
                                      byref(st), byref(self._out_c))
        self._lib.meds_leaf_turnover_step(
            self.leaf_cover, 1.0, st.leaf_flush_tendency, st.leaf_shed_tendency,
            p.flush_rate_max, p.shed_rate_max, p.leaf_turnover_rate, p.min_leaf_cover,
            p.bare_leaf_cover, float(dt), byref(self._sen), byref(self._bg), byref(self._cap))
        sen, bg = self._sen.value, self._bg.value
        post = self.leaf_cover - sen - bg
        self.leaf_cover = post + min(max(0.0, 1.0 - post), self._cap.value)
        return Day(st.leaf_flush_tendency, st.leaf_shed_tendency, self.leaf_cover, sen, bg)


#===========================================================================================#
#  Leaf habits as parameter sets, from examples/example02_canopy_phenology: Harvard Forest     #
#  (fitted), Hyytiala (fitted), the BCI light exchanger (fitted) and drought-deciduous species. #
#===========================================================================================#
def _preset(defaults, overrides) -> Params:
    names = {f.name for f in fields(Params)}
    unknown = set(overrides) - names
    if unknown:
        raise TypeError(f"unknown phenology parameter(s): {sorted(unknown)}")
    return Params(**{**defaults, **overrides})


def temperate_deciduous(**overrides) -> Params:
    """Cold-deciduous (Harvard Forest): flush on warmth while days are longer than a photoperiod
    threshold, senesce on cold once days shorten; no leaf-cover floor, so the canopy goes bare."""
    return _preset(dict(flush_cue_mask=Cue.TEMP | Cue.DAYLENGTH, shed_cue_mask=Cue.TEMP | Cue.DAYLENGTH,
                        flush_rate_max=0.03386, shed_rate_max=0.3329,
                        leaf_turnover_rate=0.0003072, flush_degree_days=93.98,
                        shed_base_temp=290.1, shed_degree_days=46.18,
                        flush_daylength_threshold=13.75, flush_daylength_sharpness=7.99,
                        shed_daylength_threshold=9.127), overrides)


def boreal_evergreen(**overrides) -> Params:
    """Evergreen conifer (Scots pine, Hyytiala): flush on warmth and bright PAR, senesce as PAR
    dims in autumn; senescence stops at min_leaf_cover -- the needles kept -- so the canopy never
    goes bare."""
    return _preset(dict(flush_cue_mask=Cue.TEMP | Cue.PAR, shed_cue_mask=Cue.TEMP | Cue.PAR,
                        flush_degree_days=88.0, shed_base_temp=299.6, shed_degree_days=455.8,
                        flush_rate_max=0.2866, shed_rate_max=0.03123, leaf_turnover_rate=0.2742,
                        min_leaf_cover=0.895, flush_par_threshold=306.5,
                        flush_par_sharpness=0.3694, shed_par_threshold=124.2,
                        shed_par_sharpness=-0.01328, par_window=16.09), overrides)


def drought_deciduous(**overrides) -> Params:
    """Drought-deciduous tropical tree (BCI example): flush and senesce on the predawn water
    potential, with a high threshold and no PAR response, so it is leafless in the dry season.
    Its water threshold is on example02's soil-water surrogate scale."""
    return _preset(dict(flush_cue_mask=Cue.WATER, shed_cue_mask=Cue.PAR | Cue.WATER,
                        leaf_psi_tlp=-0.33, flush_water_sum=3.0, flush_water_sharpness=6.667,
                        shed_water_sum=1.0, shed_water_sharpness=20.0, flush_cue_timescale=2.0,
                        flush_rate_max=0.06667, shed_rate_max=0.1, leaf_turnover_rate=0.0,
                        min_leaf_cover=0.0, shed_par_threshold=5000.0), overrides)


def light_exchanging(**overrides) -> Params:
    """Light-driven leaf exchanger (BCI example): a water threshold it never reaches, and senescence
    on bright PAR while the canopy refills, so leaves turn over in the dry season and the canopy
    stays full."""
    return _preset(dict(flush_cue_mask=Cue.WATER, shed_cue_mask=Cue.PAR | Cue.WATER,
                        leaf_psi_tlp=-1.5, flush_water_sum=3.0, flush_water_sharpness=6.667,
                        shed_water_sum=1.0, shed_water_sharpness=20.0, flush_rate_max=0.1699,
                        shed_rate_max=0.002864, leaf_turnover_rate=0.6122,
                        min_leaf_cover=0.8349, shed_par_threshold=445.3,
                        shed_par_sharpness=0.4996, par_window=3.684), overrides)


def self_test() -> None:
    """Smoke test: no cues hold the canopy full; a deciduous canopy senesces to ~bare in a cold,
    short-day autumn; an evergreen floor stops the same senescence at min_leaf_cover."""
    none = Phenology(Params())
    for doy in range(1, 30):
        day = none.step(temp_day=298.15, doy=doy)
    assert abs(day.leaf_cover - 1.0) < 1e-12 and day.senescence == 0.0, "no cues should hold full"

    autumn = dict(temp_day=275.0, daylength=8.5, par=60.0)          # cold, short and dim
    dec = Phenology(temperate_deciduous())
    for doy in range(250, 330):
        day = dec.step(doy=doy, **autumn)
    assert day.leaf_cover < 0.01, f"a deciduous canopy should be ~bare, got {day.leaf_cover}"

    evg = Phenology(boreal_evergreen())
    for doy in range(250, 330):
        day = evg.step(doy=doy, **autumn)
    floor = evg.params.min_leaf_cover
    assert abs(day.leaf_cover - floor) < 1e-6, f"an evergreen should stop at {floor}, got {day.leaf_cover}"
    print("meds.plant.pheno.self_test: OK")
