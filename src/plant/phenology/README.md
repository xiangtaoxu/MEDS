# plant/phenology

Per-individual **leaf phenology** — a stateless-in-leaf-mass directional signal generator. Given
plant traits and environmental conditions it predicts the phenological **status** (the direction of
leaf-display change), and nothing else. Links `meds_shared` only, so it compiles and tests standalone
(`cmake --build … --target meds_phenology`), exactly like `src/plant/leaf/` and `src/plant/hydraulics/`.

Design: `archive/MEDS_PHENOLOGY_DESIGN.md`.

## Public seam

`meds_plant_phenology%update_phenology(env, params, dt, state, out)` advances a cohort's cue
accumulators over one daily step `dt` [day] and returns the phenological status. Given the daily
environmental drivers (`pheno_env_t`: temperature, soil temperature, available water / soil potential,
leaf water potential, daylength, day-of-year), the per-PFT traits (`pheno_params_t`), and the
prognostic cue memory (`pheno_state_t`, passed by argument), it returns a `pheno_out_t`:

- `phenology_status ∈ { PHEN_ON, PHEN_OFF, PHEN_DORMANT }` — ON = seek leaf growth, OFF = drop leaves,
  DORMANT = hold the current state (a neutral deadband).
- `cue_limiting` — which `CUE_*` gave the lowest favorability (the governing cue).

It does **not** compute leaf growth, leaf drop, a target leaf level, or carbon: the displayed leaf
fraction *emerges downstream* from this status plus leaf growth/shed rates (a future leaf-dynamics
module).

## Cues (per-PFT enable mask — the only "strategy")

`cue_mask` is an OR of `CUE_TEMP` (GDD/CDD cold-deciduous), `CUE_WATER` (soil-water drought),
`CUE_HYDRO` (leaf-ψ hydraulic), and `CUE_PHOTO` (photoperiod, a gate on the temperature cue). Evergreen
is `CUE_NONE` (perpetually ON). Every ED2 habit is a mask + parameter special-case. Each active cue
produces a favorability in [0,1]; the **most-limiting** cue governs (`Φ = min`), banded into the
tri-state via `phen_on_threshold` / `phen_off_threshold`.

## Modules

- `meds_pheno_types`      — the interface types + cue-mask bits + status codes.
- `meds_pheno_engine`     — the `pure`/`elemental` cue kernels (thermal GDD/chill + cold-drop, water,
  hydraulic, photoperiod), the accumulator update, the most-limiting combination + banding, and the
  `daylength(lat, doy)` helper (ED2 polar-branch bug fixed).
- `meds_plant_phenology`  — the sealed public seam.

## Status

Implemented: the standalone directional kernel + `test/test_plant_phenology.f90` (evergreen, cold
GDD/CDD annual cycle, drought, hydraulic, deadband hysteresis, multi-cue limiting factor, degenerate
drivers, daylength polar fix). Validated on ifx and nvfortran multicore.

Follow-ups (see the design doc): the downstream **leaf-dynamics module** that turns the status into
leaf growth/drop and the emergent leaf fraction; the cohort-SoA `pheno_state_t` columns + demography
wiring; the `[phenology]` config block + TOML traits + the config-driven seam; a `daylength` +
latitude helper in `meds_time`; richer generic multi-driver co-limitation.
