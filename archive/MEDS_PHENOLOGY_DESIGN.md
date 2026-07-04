# MEDS Phenology Module — Design & Implementation Plan

## Context

MEDS needs a **leaf phenology** module. Scoped tightly (this design): phenology is a **pure directional
signal generator** — given plant traits and environmental conditions, it predicts the **phenological
status**, a discrete direction of change:

- **`PHEN_ON`** — conditions favor leaves; actively seeking leaf growth.
- **`PHEN_OFF`** — conditions unfavorable; actively dropping leaves.
- **`PHEN_DORMANT`** — neutral; hold the current leaf state, seek no change.

It does **not** compute leaf growth, leaf drop, or a target leaf level. The actual leaf-display fraction
(`target_frac`) **emerges downstream** from these states plus the leaf-growth / leaf-shed rates, in a
**separate leaf-dynamics module** (§8). The `src/plant/phenology/` directory is currently an empty
placeholder.

**Goal (from the request):** cover ED2's main pheno-schemes — **GDD/CDD (cold-deciduous)** and
**water-potential-driven (drought / hydraulic)** phenology — without ED2's tiered `iphen_scheme` global
switch. One **generic engine where PFTs differ only in parameters**; the only "strategy" is a per-PFT
**cue-enable mask**. (Integrating several simultaneously-active drivers *generically* — richer than the
limiting-factor default used here — is deferred to future development.)

### The framing — phenology is a directional controller; a downstream module integrates

MEDS's cohort layer runs on a **rate provider → integrator** seam (`vital_rates` → the demography engine).
Leaf phenology is a **controller**, one level down and split into two responsibilities:

- **Phenology (this module) = the direction.** Environment + traits + phenological *memory* → grow / drop /
  hold. Pure signal; no leaf mass in, no rates out. It owns only its cue *accumulators* (degree-day sums,
  chilling count, running-mean water, dry/wet-day counters) — *phenological* memory, not *leaf-mass* state.
- **Leaf-dynamics (a separate future module) = the process.** Given the direction, the current leaf state,
  and available carbon, it applies the state-dependent **leaf growth** (ON), **leaf drop** (OFF), or
  background-turnover **hold** (DORMANT). The displayed leaf fraction `target_frac` is the *integral* of
  these — emergent and path-dependent, never set directly.

**Why directional, not a level setpoint?** Phenology is fundamentally about *transitions* (leaf-out /
leaf-fall events). A directional controller models exactly that, and it makes **hysteresis intrinsic**:
DORMANT is a **deadband** between the ON and OFF thresholds, so a canopy that just dropped must climb back
above the ON threshold to re-flush (no separate anti-flicker timers, unlike ED2/FATES/CLM). And it maps
directly onto ED2's prognostic `phen_status` (1 flushing / −1 dropping / {0,−2} holding).

### Scope — deliberately narrow (leaf-display *direction* only)
- **Leaf only.** Root/wood phenology dropped (a future tissue-phenology concern).
- **Direction only.** No leaf growth/drop rates, no carbon, no leaf-mass state, no target level — those are
  the downstream leaf-dynamics module (§8).
- **Quantity only.** Leaf *quality* plasticity (Vcmax / SLA / leaf-lifespan — ED2's Kim-2012 "light
  phenology") is a **separate future trait-plasticity module**; in ED2 light phenology is almost entirely
  trait plasticity and the leaf *display* for those PFTs is water-driven, so it leaves this module cleanly.
- **Cues here:** temperature (GDD/CDD), soil-water drought, leaf-ψ hydraulic, photoperiod.
- **Multi-driver combination:** the most-limiting active cue governs (§4). Richer co-limitation of
  simultaneously-active drivers is future development.

### What the reference survey found
- ED2's real branch key is the per-PFT integer `phenology(ipft)`; every habit drives one scalar
  `elongf ∈ [0,1]` and a prognostic status machine. The instant-drought habit is legacy — dropped. Each
  ED2 PFT is governed by essentially one cue, so a limiting-factor combination reproduces it.
- FATES & CLM dispatch on a small per-PFT flag set (evergreen/cold/drought), converge on **chilling-
  adaptive GDD thresholds** and prognostic status flags.

**Design decisions:**
1. **First cut = standalone stateless kernel + CTest only** — links `meds_shared` only, drivers passed by
   argument, *not* wired into demography (MEDS has no met forcing, soil water, `psi_leaf`, or latitude
   yet). Mirrors how `src/plant/leaf/` and `src/plant/hydraulics/` shipped.
2. **Home = sibling `src/plant/phenology/`** (rescoped leaf-only; parallels `hydraulics/`).
3. **Output = a directional tri-state** `phenology_status ∈ {ON, OFF, DORMANT}` (+ `cue_limiting`). No
   `target_frac` output — it emerges downstream. Internally a continuous favorability `Φ` is computed and
   **banded** into the state (§4).
4. **Multi-cue combination = the most-limiting cue** (`Φ = min` of active-cue favorabilities) — a
   deliberately simple, ED2-faithful default; generic multi-driver co-limitation is deferred to future
   development.
5. **Cues = TEMP / WATER / HYDRO / PHOTO** (light deferred to the trait-plasticity module; water and
   hydraulic kept as separate cues).
6. **No `opts` type** (folded into per-PFT params); **no packed `N_PHENO` array** (a small named state
   record); no root phenology; no leaf-quality outputs; no rates; no target level.
7. **Accept event-based deciduousness** — stable *partial* canopies under steady moderate drought (ED2's
   continuous `elongf = paw_avg`) are not reproduced; partial canopies emerge dynamically/path-dependently
   under fluctuating forcing.
8. **Forward target = optimality overlay reserved** — sign of marginal carbon gain → grow/hold/drop is the
   native tri-state (Caldararu 2014 / P-model).

Template/style: `archive/MEDS_HYDRAULICS_DESIGN.md`. Kinds `wp/ik`, `_wp` literals, `implicit none`,
`pure`/`elemental` kernels, `error stop`, ≤132 cols; **potentials in MPa**.

---

## 1. Module layout — `src/plant/phenology/`

One static lib `libmeds_phenology.a` (`PUBLIC meds_shared`), one public seam. **Three modules.**

| File | Role | Analogue in `hydraulics/` |
|------|------|---------------------------|
| `meds_pheno_types.f90` | `pheno_env_t` (drivers), `pheno_params_t` (flat per-PFT: cue mask + selectors + cue params + on/off thresholds), `pheno_state_t` (the accumulator record), `pheno_out_t` (`phenology_status` + `cue_limiting`); cue-mask bits, status codes | `meds_hydro_types.f90` |
| `meds_pheno_engine.f90` | the `pure`/`elemental` cue kernels (thermal GDD/chill + cold-drop, water, hydraulic, photoperiod, `logistic`, `running_mean`, `daylength`) **and** `update_phenology`'s core: accumulate → per-cue favorability → most-limiting `Φ` → band into ON/OFF/DORMANT | `meds_hydro_conductance.f90` |
| `meds_plant_phenology.f90` | **THE sealed seam** `update_phenology(env, params, dt, state, out)`; re-exports the public types + constants | `meds_plant_hydraulics.f90` |
| `test/test_plant_phenology.f90` | CTest target (links `meds_phenology` only) | `test/test_plant_hydraulics.f90` |
| `src/plant/phenology/README.md` + `archive/MEDS_PHENOLOGY_DESIGN.md` | seam doc + this design | hydraulics README + design |

CMake (mirror the hydraulics block, `CMakeLists.txt` ~L95–104 and the test block ~L222–226):
```cmake
file(GLOB PHENO_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/phenology/*.f90)
add_library(meds_phenology STATIC ${PHENO_SOURCES})
target_link_libraries(meds_phenology PUBLIC meds_shared)
meds_fortran_flags(meds_phenology)
add_executable(test_plant_phenology test/test_plant_phenology.f90)
target_link_libraries(test_plant_phenology PRIVATE meds_phenology)
meds_fortran_flags(test_plant_phenology)
add_test(NAME plant_phenology COMMAND test_plant_phenology)
```

---

## 2. The stateless contract (the seam)

```fortran
subroutine update_phenology(env, params, dt, state, out)
   type(pheno_env_t),    intent(in)    :: env      ! today's raw environmental drivers (read-only)
   type(pheno_params_t), intent(in)    :: params   ! flat per-PFT traits + selectors (self-contained)
   real(wp),             intent(in)    :: dt       ! [day] phenology step (normally 1)
   type(pheno_state_t),  intent(inout) :: state    ! prognostic cue accumulators — cohort-owned, passed in
   type(pheno_out_t),    intent(out)   :: out      ! phenology_status {ON,OFF,DORMANT}, cue_limiting
end subroutine
```

* `state` is the **only mutable thing** — a small named record of cue accumulators (phenological memory),
  owned by the cohort (the FATES `*Mem`/compute split). Disjoint per cohort ⇒ reentrant. No `save`, no
  module vars, no allocatables.
* **No leaf-mass state, no rates, no target level.** The current leaf amount is *not* an input; growth/drop
  are *not* outputs. The downstream leaf-dynamics module (§8) consumes `out%phenology_status` with the leaf
  state and carbon.
* `pheno_params_t` is flattened from `cfg%pft(ipft)` by the (future) config-driven overload, exactly like
  `leaf_photo_params_t`; the engine never reaches into `meds_config`.

---

## 3. Types

### `pheno_env_t` — raw daily environmental drivers (read-only)
| Field | Units | Notes |
|-------|-------|-------|
| `temp_day` | K | daily-mean air/canopy temperature (thermal sums) |
| `soil_temp` | K | shallow-layer soil temperature (cold-drop trigger) |
| `avail_water` | – or MPa | available-water fraction ∈[0,1] **or** soil water potential (sign-selected) — CUE_WATER |
| `psi_leaf` | MPa (<0) | leaf/xylem water potential — CUE_HYDRO (bridges to `src/plant/hydraulics`) |
| `daylength` | h | photoperiod — CUE_PHOTO; caller-supplied (pure `daylength(lat,doy)` helper provided for tests) |
| `doy`, `hemis_north` | –, logical | day-of-year + hemisphere for season gating of thermal sums |

### `pheno_state_t` — cue accumulators (the prognostic memory; named record, no packed array)
| Field | Meaning | Units | Cue |
|-------|---------|-------|-----|
| `gdd` | growing-degree-day sum | K·day | TEMP |
| `chill` | chilling-day count | day | TEMP |
| `water_avg` | running-mean available water (≈10-day) | – or MPa | WATER |
| `low_psi_days` | consecutive days `psi_leaf < ψ_tlp` | day | HYDRO |
| `high_psi_days` | consecutive days `psi_leaf ≥ ½ψ_tlp` | day | HYDRO |

(In the future SoA these become five named cohort columns, matching the `growth_accum`/`growth_count`
idiom — not a packed block. `phenology_status` is derived each step, not stored.)

### `pheno_params_t` — flat per-PFT traits + selectors (self-contained)
| Group | Fields | Notes |
|-------|--------|-------|
| **selector** | `cue_mask` (int, OR of `CUE_*`), `cue_sharpness`, `phen_on_threshold` (θ_on), `phen_off_threshold` (θ_off) | mask = the only "strategy"; θ_on > θ_off define the DORMANT deadband |
| **thermal** (CUE_TEMP) | `gdd_base_temp`, `chill_base_temp`, `phen_a`, `phen_b`, `phen_c`, `cold_drop_daylength`, `cold_drop_soiltemp1`, `cold_drop_soiltemp2` | GDD thresh `= a+b·exp(c·chill)` (chilling-adaptive); autumn drop by daylength+soil-temp (White 1997 / Botta 2000) |
| **water** (CUE_WATER) | `water_use_potential` (sign flag), `water_off_threshold`, `water_on_threshold`, `water_window` | FATES sign-trick (moisture vs soil-ψ); ~10-day mean |
| **hydraulic** (CUE_HYDRO) | `leaf_psi_tlp` [MPa], `low_psi_threshold`, `high_psi_threshold` [day] | Xu 2016 dry/wet consecutive-day counters vs ψ_tlp |
| **photoperiod** (CUE_PHOTO) | `photo_crit` [h], `photo_slope` | latitude-gated optional cue (extratropics) |

(No rate scales, no target — `k_grow`/`k_shed`/`turnover_base` live in the downstream leaf-dynamics
module.) `derive_pheno_params`: unit conversions (ψ_tlp, daylength) — pure trait-table transform.

### `pheno_out_t` — the phenological status
| Field | Units | Notes |
|-------|-------|-------|
| `phenology_status` | – | the direction of change: `PHEN_ON` \| `PHEN_OFF` \| `PHEN_DORMANT` |
| `cue_limiting` | flag | which cue gave the lowest favorability (the governing cue, `argmin_i f_i`) |

Cue-mask bits: `CUE_NONE=0`, `CUE_TEMP=1`, `CUE_WATER=2`, `CUE_HYDRO=4`, `CUE_PHOTO=8`. Status codes:
`PHEN_ON=1`, `PHEN_DORMANT=0`, `PHEN_OFF=-1`.

---

## 4. The generic engine (`update_phenology`)

Daily update, two steps:

**(1) Accumulate** into `state` from `env`:
- Thermal (CUE_TEMP): hemisphere/season-gated `gdd += (temp_day − gdd_base_temp)·dt` on warm growing-season
  days; `chill += dt` on cold dropping-season days; seasonal resets (ED2 `update_thermal_sums`).
- Water (CUE_WATER): `water_avg` exponential running mean (weight `dt/water_window`).
- Hydraulic (CUE_HYDRO): `low_psi_days += dt` while `psi_leaf < ψ_tlp` else 0; `high_psi_days += dt` while
  `psi_leaf ≥ ½ψ_tlp` else 0.

**(2) Per-cue favorability → most-limiting `Φ` → band into the tri-state.** Each active cue → `f_i ∈ [0,1]`
via `logistic(cue_sharpness·(x−x*))`:
- **CUE_TEMP:** `f = logistic(gdd − gdd_thresh(chill)) · [1 − cold_drop(daylength, soil_temp)]`.
- **CUE_WATER:** `f = ramp(water_avg between off/on thresholds)` (moisture or soil-ψ by sign).
- **CUE_HYDRO:** `f` rises with `high_psi_days` (vs `high_psi_threshold`), falls with `low_psi_days`
  (vs `low_psi_threshold`).
- **CUE_PHOTO:** `f_photo = logistic(photo_slope·(daylength − photo_crit))` — **multiplies** the temp gate.

`Φ = min_i f_i` over the active cues (the most-limiting cue governs; `cue_limiting = argmin_i f_i`), then
**band** (memoryless; the deadband is the hysteresis, so no latch/state needed):
```
if      (Φ > phen_on_threshold)   phenology_status = PHEN_ON        ! favorable — seek growth
else if (Φ < phen_off_threshold)  phenology_status = PHEN_OFF       ! unfavorable — drop
else                              phenology_status = PHEN_DORMANT   ! neutral deadband — hold
```
Evergreen (`CUE_NONE`) ⇒ `Φ ≡ 1` ⇒ perpetually `PHEN_ON` (the leaf-dynamics module then holds a full canopy
via background-turnover replacement — ON means "growth favored," not "not yet full"). Everything downstream
— realizing ON/OFF/DORMANT against actual leaf mass and carbon — is the leaf-dynamics module (§8), **not**
here. (Generic co-limitation of simultaneously-active cues, beyond this limiting-factor `min`, is future
development.)

---

## 5. Coverage — every retained ED2 habit is a mask + parameter special-case (no code branches)

| ED2 habit | `cue_mask` | Phenological status over the year |
|-----------|------------|-----------------------------------|
| evergreen | `CUE_NONE` | perpetually `ON` (canopy held full downstream) |
| cold-GDD (Botta/White) | `CUE_TEMP` (+`CUE_PHOTO` opt.) | winter `OFF`; spring GDD thresh `a+b·exp(c·chill)` → `ON`; summer `ON` (held full downstream); autumn daylength+soil-temp → `OFF` |
| drought (10-day) | `CUE_WATER` | wet → `ON`; dry → `OFF`; marginal deadband → `DORMANT` (hold) |
| hydraulic (Xu) | `CUE_HYDRO` | sustained ψ_leaf≥½ψ_tlp → `ON`; sustained ψ_leaf<ψ_tlp → `OFF` |
| light (deferred) | — | trait plasticity → separate leaf-quality module |

(Multiple cues on one PFT already work via the limiting-factor `min`; a richer generic multi-driver
integration is future development — §8.)

---

## 6. Numerics / GPU (nvfortran gate)

Kernels `pure`/`elemental`, arithmetic + intrinsics only; flat params; fixed-size named state ⇒ no
allocatables/runtime shapes. Branch-light: cues are smooth logistics; only the season-gate, the counter
resets, and the final 3-band `if` are data-dependent (cheap, warp-friendly). **FPE-safe in Debug**
(`-fpe0`/`-Ktrap=fp`): `logistic` uses `safe_exp` (`meds_constants`), running-mean weights bounded, no
`0/0`, no `(neg)**real`; use `if` not `merge` for the banding. A green **ifx** run is not sufficient —
build **nvfortran multicore** on the new module (CLAUDE.md rule).

---

## 7. Validation — `test/test_plant_phenology.f90`

Tests exercise the **status** directly (no leaf-mass integration — that is the leaf-dynamics module's
test), driving the kernel over synthetic forcing series.

1. **Evergreen** — `CUE_NONE`: perpetually `PHEN_ON`, invariant under any driver sequence.
2. **Cold-deciduous GDD/CDD** — seasonal temperature: `PHEN_ON` from spring (once `gdd ≥ a+b·exp(c·chill)`)
   through summer, `PHEN_OFF` in autumn (short daylength + cold soil) and winter, `PHEN_DORMANT` in the
   shoulder transitions; more chilling lowers the heat requirement; reproduces ED2 leaf-on/off timing in the
   sharp-slope (`cue_sharpness→∞`) limit.
3. **Drought (10-day)** — dry-down/wet-up: `OFF` under sustained low water, `ON` on recovery, `DORMANT` in
   the deadband between; moisture and soil-ψ modes (sign trick).
4. **Hydraulic (Xu)** — `OFF` after `psi_leaf < ψ_tlp` for ≥`low_psi_threshold` days; `ON` after
   `psi_leaf ≥ ½ψ_tlp` for ≥`high_psi_threshold` days.
5. **Deadband hysteresis** — a driver oscillating across θ_off then up toward θ_on shows the OFF→DORMANT→ON
   path with **no chatter** (DORMANT holds between); confirms the intrinsic anti-flicker.
6. **Multi-cue limiting factor** — cold + drought both enabled: the more restrictive cue governs the status
   (`Φ = min`), and `cue_limiting` reports which; a mid-summer drought forces `OFF` even while the thermal
   cue is favorable.
7. **Range / FPE / degenerate** — status ∈ {ON,OFF,DORMANT} always; no NaN/trap under degenerate drivers
   (constant, empty mask, extreme temperatures/potentials).
8. **nvfortran multicore build green** (portability gate).

---

## 8. Reserved follow-ups (documented; **not** in this first cut)

- **Leaf-dynamics module (THE consumer of the status):** the state- and carbon-dependent half. Given
  `phenology_status`, the current leaf state, and available carbon, it applies **leaf growth** (ON, rate
  `k_grow`, capped at the allometric maximum), **leaf drop** (OFF, rate `k_shed`), or a background-turnover
  **hold** (DORMANT). The displayed leaf fraction `target_frac` (≡ `elongf`) is the *integral* of these —
  emergent, path-dependent. Emits `leaf_drop` → litter and `leaf_pheno_frac·leaf_area` → the LAI
  competition sweep in `vital_rates.f90`.
- **Generic multi-driver integration:** richer combination of simultaneously-active cues than the
  limiting-factor `min` (e.g. smooth co-limitation / weighted integration), so several drivers jointly
  shape the transition timing. Deferred to future development.
- **Cohort-SoA accumulators + demography wiring:** add the five `pheno_state_t` fields as named cohort
  columns; touch the **7 lockstep routines** in `meds_demography_types.f90` (`cohort_alloc`, `site_free`,
  `cohort_ensure_capacity`, `move_alloc_block`, `cohort_reorder`, `cohort_compact`, `copy_cohort_slot`);
  init at every creation site (grep `assign_cohort_id`); fusion/fission blend accumulators by `nplant`.
- **PFT traits + config:** add the phenology fields to `pft_table_t` (`alloc_pft_table` +
  `derive_pheno_params` + `validate_config`, incl. θ_on > θ_off), `[pft]` rows in `meds_config_pft.toml`,
  and a `[phenology]` block in `meds_config_main.toml` (mirroring `[leaf_physiology]`); add the
  config-driven overload `update_phenology(env, cfg, ipft, …)`.
- **`meds_time` daylength + latitude:** add `daylength(lat, doy)` (fix ED2's polar-branch `≤ 1` → `≤ −1`
  bug) + a site latitude field (MEDS has none) so CUE_PHOTO can be driven internally.
- **Stepper coupling:** call `update_phenology` **daily, before the leaf-dynamics update and `vital_rates`**,
  in `meds_stepper` — needs a met forcing source, soil water, and `psi_leaf` (deferred with the fast loop).
- **Leaf-quality (trait-plasticity) module:** Vcmax/SLA/leaf-lifespan acclimation (ED2 Kim-2012 "light
  phenology") as a sibling module that **consumes the leaf-dynamics turnover rate** as the acclimation
  timescale — the leaf-quality twin of leaf-quantity dynamics.
- **Optimality overlay (forward direction):** replace the fixed-threshold banding with a marginal-carbon-
  gain signal from `meds_leaf_physiology` + `meds_plant_hydraulics` — the **sign** of marginal gain is the
  native tri-state (grow / hold / drop). Same seam and output; only how the direction is decided changes
  (Caldararu 2014 / P-model lineage).
- **C-API:** `meds_pheno_capi.f90` + `-DMEDS_BUILD_PYLIB` → `libmeds_pheno_c` for a `meds.pheno` Python
  package (as `src/plant/leaf/meds_leaf_capi.f90` does).

---

## 9. Open questions / decisions during implementation

1. **ED2 defects to fix (recommended):** the stale `phenology_status` doc (code uses −2, not 2) and the
   `daylength` polar-branch bug (`arg ≤ 1` should be `≤ −1`).
2. **Default `cue_sharpness`** — visibly smooth by default; the sharp limit is exercised only in the
   ED2-comparison test.
3. **Deadband width** — default `phen_on_threshold`/`phen_off_threshold` (e.g. 0.6 / 0.4?) balancing
   responsiveness against flicker, per cue-family.
4. **Example PFT set** — seed the three MEDS PFTs (pioneer/mid/climax) with an evergreen + cold-deciduous +
   drought-deciduous demonstrator set (values from ED2 `init_pft_phen_params` + FATES defaults) for tests.

---

## 10. Verification

- **Build (netCDF-free, strict):** `cmake -S . -B build-debug -DCMAKE_Fortran_COMPILER=ifx
  -DCMAKE_BUILD_TYPE=Debug -DMEDS_ENABLE_IO=OFF && cmake --build build-debug --target test_plant_phenology`,
  then `ctest --test-dir build-debug -R plant_phenology --output-on-failure`.
- **Standalone lib compiles:** `cmake --build build-debug --target meds_phenology`.
- **Portability gate (required):** rebuild the module + test under **nvfortran multicore**
  (`-DCMAKE_Fortran_COMPILER=nvfortran -DMEDS_GPU=multicore`) — a green ifx suite is not sufficient.
- **Science checks:** the 8 CTest cases in §7.

---

## 11. References
- **ED2** `dynamics/phenology_driv.f90`, `phenology_aux.f90`, `memory/{phenology_coms,pft_coms}.f90`,
  `init/ed_params.f90` (`init_phen_coms`, `init_pft_phen_params`) — the algorithmic reference.
- **FATES** `biogeochem/EDPhysiologyMod.F90` + tech note — chilling-adaptive GDD, the soil-water sign trick,
  prognostic status flags.
- **CLM5** `CNPhenologyMod` tech note — seasonal (GDD + `crit_dayl`) vs stress (soil-ψ ±0.6/−0.8 MPa)
  deciduous, onset/offset flags.
- **Chuine (2000)** unified budburst; **Delpierre (2009)** autumn (CDD×photoperiod); **turgor-loss-point**
  shedding (ψ_tlp) — the smooth-cue formulations.
- **Caldararu (2014)** / P-model — the reserved optimality overlay.
- Template: `archive/MEDS_HYDRAULICS_DESIGN.md`.
```
