"""meds.pheno — leaf-phenology SIGNAL kernel.

A Pythonic front end to the MEDS Fortran phenology kernel (meds_phenology.f90). Given daily
environmental cues + per-PFT traits it returns two RELATIVE rate tendencies — leaf_flush_rate and
leaf_shed_rate [1/day] — from two governor accumulators it advances in place. The compiled Fortran
does the arithmetic; this module exposes it with dataclasses + an IntFlag cue mask, so callers never
touch ctypes.

    import meds.pheno as pheno

    ph = pheno.Phenology(pheno.temperate_deciduous())     # a stateful driver (holds params + memory)
    out = ph.step(temp_day=290.0, soil_temp=290.0, daylength=13.0, doy=150)
    print(out.leaf_flush_rate, out.leaf_shed_rate)        # [1/day]

    # the module does NOT track leaf area; integrate relative LAI (canopy fullness in [0,1]) yourself:
    lai = pheno.integrate_lai(lai, out.leaf_flush_rate, out.leaf_shed_rate, baseline_turnover=8e-4)

The kernel is signal-only: it never touches carbon or leaf mass. The actual leaf/storage carbon
update lives in the Fortran meds_plant_carbon_dynamics; `integrate_lai` is a compact relative-LAI
analogue of it for demos (see examples/example_phenology). Requires the compiled shared library
(see meds.leaf._ffi for the one-time CMake build).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from enum import IntFlag

from . import _ffi

__all__ = [
    "Cue", "Params", "State", "Out", "Phenology", "step", "integrate_lai",
    "temperate_deciduous", "temperate_evergreen", "drought_deciduous", "light_exchanging",
    "self_test",
]


class Cue(IntFlag):
    """Cue-enable bits (mirror meds_plant_types). OR them for a mask; flush and shed sides are
    selected independently (min over flush cues, max over shed cues)."""
    NONE = 0
    TEMP = 1     # temperature: GDD flush + autumn cold-drop shed
    WATER = 2    # soil-water running mean            (driver not wired in the standalone model yet)
    HYDRO = 4    # daily-max leaf water potential (dmax_leaf_psi)
    PHOTO = 8    # photoperiod (gates the temperature flush)
    LIGHT = 16   # radiation: active shed rises with running-mean light


#----- Per-PFT parameters. Defaults mirror pheno_params_t (meds_plant_types). ---------------#
@dataclass
class Params:
    """Per-PFT phenology traits. Rates are RELATIVE, per day; the two masks pick which cues drive
    each side. k_flush_max = 1/15 d => a bare canopy fills in ~15 days at full drive."""
    flush_cue_mask: int = Cue.NONE
    shed_cue_mask: int = Cue.NONE
    cue_sharpness: float = 2.0
    k_flush_max: float = 1.0 / 15.0     # [1/day] max relative flush rate (~full in 15 d)
    k_shed_max: float = 1.0 / 20.0      # [1/day] max relative active-shed rate (~bare in 20 d)
    tau_flush: float = 5.0              # [day] flush governor low-pass timescale
    tau_shed: float = 5.0               # [day] shed governor low-pass timescale
    # thermal (CUE_TEMP): GDD flush a+b*exp(c*chill) + autumn cold drop
    gdd_base_temp: float = 278.15
    chill_base_temp: float = 278.15
    phen_a: float = -68.0
    phen_b: float = 638.0
    phen_c: float = -0.01
    cold_drop_daylength: float = 10.9
    cold_drop_soiltemp1: float = 284.3
    cold_drop_soiltemp2: float = 275.15
    # water (CUE_WATER)
    water_use_potential: bool = False
    water_off_threshold: float = 0.2
    water_on_threshold: float = 0.5
    water_window: float = 10.0
    water_width: float = 0.1
    # hydraulic (CUE_HYDRO): dmax_leaf_psi vs turgor-loss point
    leaf_psi_tlp: float = -2.0
    low_psi_threshold: float = 10.0
    high_psi_threshold: float = 10.0
    # photoperiod (CUE_PHOTO)
    photo_crit: float = 11.0
    photo_slope: float = 2.0
    # light (CUE_LIGHT): active shed rises with running-mean radiation
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
    """One step's phenology signal."""
    leaf_flush_rate: float    # [1/day] relative flush tendency (0 => dormant)
    leaf_shed_rate: float     # [1/day] relative active-shed tendency (0 => no active shed)
    cue_limiting: int         # the strongest active shed cue (a Cue bit; diagnostic)


#----- Daily environment defaults (a warm, well-lit, well-watered, long day). ----------------#
_ENV_DEFAULTS = dict(temp_day=298.15, soil_temp=298.15, avail_water=0.5, dmax_leaf_psi=0.0,
                     rad=400.0, daylength=12.0, doy=1, hemis_north=True)


def step(env, params, state, dt=1.0):
    """Low-level one-step call. `env` is a dict (missing keys default), `params`/`state` are Params/
    State (or dicts). Returns (Out, new_State). `state` is NOT mutated; feed new_State back."""
    env_full = {**_ENV_DEFAULTS, **dict(env)}
    p = asdict(params) if isinstance(params, Params) else dict(params)
    s = asdict(state) if isinstance(state, State) else dict(state)
    out, new_s = _ffi.step(env_full, p, s, dt=dt)
    return Out(**out), State(**new_s)


class Phenology:
    """A stateful phenology driver: holds the per-PFT params + the advancing memory. Call `.step`
    each day with the environment; it advances the memory and returns the two rates."""

    def __init__(self, params: Params, state: State | None = None):
        self.params = params
        self.state = state if state is not None else State()

    def step(self, dt: float = 1.0, **env) -> Out:
        out, self.state = step(env, self.params, self.state, dt=dt)
        return out


def integrate_lai(elongf, leaf_flush_rate, leaf_shed_rate, dt=1.0,
                  baseline_turnover=0.0, elongf_min=0.02):
    """Advance relative LAI (canopy fullness, elongf in [0,1]) ONE step from the two phenology rates.

    A compact relative-unit analogue of the Fortran carbon leaf update (get_plant_flux_slow):
      * flush FILLS toward full at leaf_flush_rate (linear, capped by the deficit 1-elongf),
      * the ACTIVE shed removes leaf_shed_rate per day (linear toward bare), NON-replaceable,
      * `baseline_turnover` [1/day] is a small proportional, REPLACEABLE background loss,
      * when a net decline crosses below elongf_min it SNAPS to bare (ED2's fully-abscised state).
    Returns the new elongf in [0,1]. LAI tracking is the caller's job — the kernel is signal-only.
    """
    flush_gain = min(leaf_flush_rate * dt, max(0.0, 1.0 - elongf))
    loss = min(leaf_shed_rate * dt + baseline_turnover * elongf * dt, elongf)
    e = elongf + flush_gain - loss
    if leaf_shed_rate > 0.0 and loss > flush_gain and e < elongf_min:
        e = 0.0
    return min(1.0, max(0.0, e))


#===========================================================================================#
#  The four target phenological strategies (design §1a / §6.1), as ready-to-use Params.       #
#===========================================================================================#
def _preset(defaults, overrides) -> Params:
    """Build a Params from a preset's fixed fields, letting caller overrides win."""
    return Params(**{**defaults, **overrides})


def temperate_deciduous(**overrides) -> Params:
    """Cold-deciduous: flush on spring GDD, shed on the autumn cold-drop (both masks = TEMP)."""
    return _preset(dict(flush_cue_mask=Cue.TEMP, shed_cue_mask=Cue.TEMP,
                        k_flush_max=1.0 / 15.0, k_shed_max=1.0 / 18.0), overrides)


def temperate_evergreen(**overrides) -> Params:
    """Evergreen: permissive flush (no onset cue), no active shed — the canopy is held full year
    round, thinned only by the small baseline turnover the caller applies in integrate_lai. Pass
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
    """Smoke test: evergreen holds flush=k_flush_max / shed=0; a cold winter drives a deciduous
    cohort's shed drive up and its flush drive down."""
    ev = Phenology(temperate_evergreen())
    o = ev.step(temp_day=298.15, doy=180)
    assert abs(o.leaf_flush_rate - ev.params.k_flush_max) < 1e-9, "evergreen flush should be k_flush_max"
    assert o.leaf_shed_rate == 0.0, "evergreen should have no active shed"

    dec = Phenology(temperate_deciduous())
    for doy in range(1, 60):                         # deep-winter cold spell
        out = dec.step(temp_day=265.0, soil_temp=265.0, daylength=9.0, doy=doy)
    assert out.leaf_shed_rate > 0.0, "cold winter should raise the deciduous shed rate"

    lai = 1.0
    for _ in range(40):                              # a full canopy under a hard shed goes bare
        lai = integrate_lai(lai, 0.0, dec.params.k_shed_max)
    assert lai < 0.05, f"canopy should shed to ~bare, got {lai}"
    print("meds.pheno.self_test: OK")
