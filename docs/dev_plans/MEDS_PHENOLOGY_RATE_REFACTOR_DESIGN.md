# MEDS Phenology — Rate-Based Refactor Design & Plan

**Status:** **P0–P2 IMPLEMENTED** (branch `feature/phenology-rate-refactor`); P3–P4 deferred. Supersedes
parts of `MEDS_PHENOLOGY_DESIGN.md` (see §1). Done: the signal-only two-rate kernel with the full cue set
(two masks + `CUE_LIGHT`, all four target patterns unit-tested), the drive state + lockstep, the config
(two masks + rate params; WATER/HYDRO/LIGHT still rejected by validate until their drivers are wired), the
driver fold into `meds_vegetation_dynamics`, and the carbon single-authority (two loss channels, linear
active shed + snap, flush cap, phenology-off bit-identical). ifx 32/32; nvfortran-mc phenology/carbon green
(`plant_hydraulics` is a pre-existing `-Ktrap=fp` flake). **Deferred:** P3 = thread the WATER/HYDRO/LIGHT
drivers (soil water, `dmax_leaf_psi`, `rad_avg`) from the fast loop into `advance_leaf_phenology` + lift the
validate rejection; P4 = `retained_carbon_fraction`.
**Scope:** replace the phenology module's *directional tri-state* output with **two relative rate
tendencies** (`leaf_flush_rate`, `leaf_shed_rate`, both **per day**), driven by **two prognostic
accumulators**, generalizing every ED2 phenological habit through per-PFT parameters only. The phenology
module is a **pure signal generator** — it touches no carbon and no leaf/storage state; the *actual*
leaf/storage carbon update lives entirely in `plant_carbon_dynamics`. Relocate four utility helpers to
`shared/util`.

This plan was produced by a multi-agent design pass (4 specialists → synthesis → 3 adversarial
reviewers); the reviewers' two blockers + majors are folded in (flagged ⚑). It was then revised per two
rounds of user direction (locked decisions below).

### Decisions locked (user, 2026-07-19)
1. **Phenology is signal-only — no carbon, no state mutation.** The module outputs *only* the two relative
   rate tendencies. It never reads or writes `leaf_carbon`/`nonstructural`/leaf area. **All** leaf/storage
   carbon update (flush fill, active shed, retained→storage, litter outflow) happens in
   `plant_carbon_dynamics` (§5).
2. **No `elongf` in the module.** `elongf` (leaf_carbon ÷ full-canopy leaf carbon) is **not** an input and
   **not** stored by phenology; it is a **diagnostic** computed outside the module (in `plant_carbon_dynamics`
   for its pool bounds, and in IO for output). The kernel is environment-cues + traits only.
3. **Two relative rates, per day; absolute flux set by `plant_carbon_dynamics`.** `leaf_flush_rate` and
   `leaf_shed_rate` are **relative** tendencies `[1/day]` (`= k_·_max · drive`). `plant_carbon_dynamics`
   turns each into the actual carbon flux (× a reference pool × dt, clamped to the pool). `leaf_flush_rate
   = 0` ⇔ dormancy (no active flush); `leaf_shed_rate = 0` ⇔ no active shedding (baseline turnover, a
   *separate* concern in `plant_carbon_dynamics`, still runs).
4. **Controller shape = pure two-rate.** No `elongf_target` output. Deciduous-under-stress is handled by
   the **min-flush / max-shed** cue combine (a drought-deciduous PFT carries `CUE_WATER` in its flush cues,
   so `min` drives flush→0 under drought) — no explicit "flush-off-while-shedding" coupling is needed (§4).
5. **Active shed is LINEAR toward bare, snapped at `leaf_carbon_min`.** In `plant_carbon_dynamics`, the
   active shed carbon flux is `leaf_shed_rate · leaf_carbon_full · dt_day` (referenced to the **full-canopy**
   allometric leaf carbon `leaf_carbon_full = size2leaf_carbon(dbh,h,sla)`), a **constant absolute decline**
   — a true linear full→bare traversal in exactly `1/leaf_shed_rate` days, **no exponential tail** — clamped
   to the pool and snapped to true bare below `leaf_carbon_min`. Flush fills symmetrically:
   `leaf_flush_rate · leaf_carbon_full · dt_day` toward the target (§5).
6. **Litter = untracked outflow (as today).** The refactor guarantees only **leaf+storage** closure;
   shed/turnover carbon leaves untracked exactly as baseline turnover does now. A diagnostic litter
   accumulator + the soil-carbon wiring are reserved follow-ups, out of scope (§5).
7. **Must reproduce four target patterns** (acceptance criteria, §1a). These require **flush and shed to
   be driven by potentially *different* cues**, so `cue_mask` is **split into `flush_cue_mask` +
   `shed_cue_mask`** (each side takes `min`/`max` over its own mask), and a **`CUE_LIGHT`** cue is added
   whose shed signal **rises with light** (running-mean radiation, ED2 `rad_avg`) to drive light-driven
   leaf exchange while flush stays high.
8. **Fold `meds_phenology_driver.f90` into `meds_vegetation_dynamics.f90`.** The phenology driver is a thin
   slow-loop marshalling seam whose two rates `carbon_growth` (already in `meds_vegetation_dynamics`)
   consumes in the *same* step; ED2 likewise calls `phenology_driver` inside the vegetation-dynamics slow
   loop. The phenology advance moves in as the **first step** of `vegetation_dynamics`, the file is deleted,
   and `meds_stepper` drops its separate call. This keeps `meds_vegetation_dynamics` as "THE single place
   plant and demography meet" (§7.2).

---

## 1. Motivation and what changes

### The current design (superseded)
`meds_phenology` (`src/plant/meds_phenology.f90`) folds each active cue into one favorability
`f_i ∈ [0,1]`, takes `phi = min_i f_i`, and **bands** `phi` into `PHEN_ON / PHEN_OFF / PHEN_DORMANT`
(the deadband is the hysteresis). The **only** consumer is `plant_carbon_allocation`'s
`can_flush = (phenology_status == PHEN_ON)` boolean, which gates leaf/fine-root growth **from storage**
(priority P2). **There is no active leaf shedding anywhere in MEDS today** — only baseline turnover
(`loss%leaf = rates%leaf·leaf_carbon·dt` in `get_plant_flux_slow`) removes leaves, so deciduousness is
effectively absent. `MEDS_PHENOLOGY_DESIGN.md` §8 deliberately **deferred** the leaf-dynamics rates
(`k_grow`/`k_shed`), leaf-shed→litter, and retained→storage.

### The pivot
Turn phenology into a **pure rate-tendency generator**: env + traits → two relative rates
(`leaf_flush_rate`, `leaf_shed_rate` `[1/day]`). The tri-state and `PHEN_*` constants leave the module.
The **actual** leaf/storage carbon update — deferred by `MEDS_PHENOLOGY_DESIGN.md` §8 — is folded into
`plant_carbon_dynamics` (not a new module, not the driver): it multiplies the relative rates by the
allometric leaf pool and `dt`, clamps to the pool, splits shed carbon to storage/litter, and reconciles
with baseline turnover.

**Why relative rates are the right generalization.** ED2 schemes 1–4 set the leaf-display fraction
`elongf` **instantaneously** to a conditions-derived target; ED2 scheme 5 (plant-hydraulic deciduous)
**moves** `elongf` toward its target at per-day `leaf_grow_rate`/`leaf_shed_rate`. Both are the same
object — a rate-limited relaxation of leaf display — and the instantaneous case is the rate→∞ limit. MEDS
makes **everything** rate-limited (no instantaneous LAI jumps). This is also required by the sub-step fast
loop: a discontinuous `leaf_carbon` jump injects an un-budgeted shock into the leaf/CAS energy balance,
two-stream RT, and soil-water root sink (all of which close to machine precision per sub-step). ED2's
sharp-drop habits are recovered as the large-`k` limit for regression comparison.

### Supersession map (vs `MEDS_PHENOLOGY_DESIGN.md`)
| Old section | Fate |
|---|---|
| §3 tri-state output `phenology_status` + `cue_limiting` | → `pheno_out_t{leaf_flush_rate, leaf_shed_rate, cue_limiting}` (§3) |
| §4 banding into ON/OFF/DORMANT; §6 `PHEN_*` semantics | **obsolete** — replaced by the two-governor low-pass + rate map (§4) |
| §8 bullet 1 (deferred leaf-dynamics `k_grow`/`k_shed`, shed→litter, retained→storage) | the **rate tendencies** become the phenology output (§4); the **carbon application** lands in `plant_carbon_dynamics` (§5) |
| §8 cohort-SoA accumulators, PFT traits/config, `meds_time` daylength | now **required work** (§6, §7) |
| Kim-2012 leaf **quality** plasticity (SLA/Vcmax/llspan); prescribed phenology (`iphen_scheme=1`) | **out of scope** — quality is a sibling module consuming `leaf_shed_rate` as its acclimation timescale; prescribed phenology is a `src/forcing` prescribed-`elongf` overlay, not a cue |

### §1a. The four target patterns (acceptance criteria)
The refactor must reproduce these four canopy behaviors from per-PFT params alone (decision 7). The key is
that **flush and shed can be driven by different cues** (hence the split `flush_cue_mask` / `shed_cue_mask`,
§3) — evergreen flushes on temperature but never actively sheds; a leaf-exchanger flushes constantly but
sheds on light. Exact param settings are in §6.1; the mechanics are §3–§5.

| # | Pattern | leaf_flush_rate | leaf_shed_rate (active) | canopy result |
|---|---|---|---|---|
| 1 | **Temperate evergreen** | high & ~constant (full in ~15 d), *or* seasonal on temperature | 0 (baseline turnover only), *or* mild | stays full year-round; winter retention comes from cold-suppressed **baseline** turnover pairing with (possibly paused) flush |
| 2 | **Temperate deciduous** | temperature (GDD): spring pulse, **0 in dormancy** | temperature (cold-drop): autumn pulse; **0 or very high in dormancy** (irrelevant once bare) | full in growing season, bare in winter |
| 3 | **Facultative drought-deciduous** | hydraulic: rises with sustained `dmax_leaf_psi ≥ ½ψ_tlp` | hydraulic: rises with sustained `dmax_leaf_psi < ψ_tlp` | full when watered (facultatively evergreen), sheds under drought, reflushes on rewet — ED2 scheme 5 |
| 4 | **Light-driven leaf-exchanging** | high & ~constant (full in ~15 d), *or* matched to shed | **rises with light** (`CUE_LIGHT`, running-mean radiation) | stays ~full while turning leaves over fast — flush keeps up with the light-driven shed |

Design consequences: pattern 1 needs flush-from-TEMP **without** shed-from-TEMP ⇒ separate masks; pattern 4
needs an **active light shed with sustained high flush** ⇒ `CUE_LIGHT` on the shed side + permissive/high
flush; pattern 3 needs the **daily-max** leaf psi `dmax_leaf_psi` (not instantaneous), a fast-loop reduction.
Pattern 1's winter canopy hold is not a phenology trick — it is the existing evergreen cold-suppression of
baseline turnover (`f_t` in `tissue_turnover_rates`, §5.2) going to ~0 in the cold, so the canopy persists
even while seasonal flush pauses.

---

## 2. The two responsibilities, in one picture

```
                    [ meds_phenology — pure SIGNAL kernel, src/plant ]
env cues (temp/water/dmax_psi/rad/daylength/doy) ─┐   NO leaf_carbon, NO storage, NO elongf in/out
                                          ▼
  (1) accumulate cue memory (gdd, chill, water_avg, dry/wet-day counters, light_avg)
  (2) per-cue  s_flush ∈[0,1] (onset)         s_shed ∈[0,1] (stress/light)
      combine: s_flush = MIN over flush_cue_mask   s_shed = MAX over shed_cue_mask   (SEPARATE masks)
  (3) low-pass the two GOVERNORS (the prognostic pheno_state):
        flush_drive ← flush_drive + w_f·(s_flush − flush_drive)
        shed_drive  ← shed_drive  + w_s·(s_shed  − shed_drive)
  (4) emit RELATIVE rates [1/day]:
        leaf_flush_rate = k_flush_max · flush_drive        (0 ⇒ dormant)
        leaf_shed_rate  = k_shed_max  · shed_drive         (0 ⇒ no active shed)
                                          │
                                          ▼   two relative rates [1/day] + cue_limiting diagnostic
   ─────────────────────────────────────────────────────────────────────────────────────────────
                    [ plant_carbon_dynamics — the ACTUAL leaf/storage update, src/plant, §5 ]
   receives the two rates + leaf_carbon_full (=size2leaf_carbon) + dt_day + current pools:
     • flush fill   :  min(target − leaf_carbon,  leaf_flush_rate · leaf_carbon_full · dt_day)
     • baseline turn:  rates%leaf(tissue_temp) · leaf_carbon · dt_day     (tissue_turnover_rates; replaceable)
     • active shed  :  min(leaf_shed_rate · leaf_carbon_full · dt_day,  leaf_carbon)   (NON-replaceable,
                       linear, snapped to bare < leaf_carbon_min); retained→storage, rest→litter (untracked)
     • updates leaf_carbon / nonstructural  (elongf = leaf_carbon/leaf_carbon_full is a diagnostic here)
```

The two boxes are the two responsibilities the user separated: **phenology decides *how fast* (relative);
`plant_carbon_dynamics` decides *how much* (absolute) and moves the carbon.**

---

## 3. The phenology contract (types, units, boundary semantics)

The kernel is environment-cues + traits only. It emits the two **relative rates** `[1/day]` and carries
the two governor accumulators as its prognostic memory. `elongf` and every carbon quantity are absent.

Cue-mask bits (unchanged plus `CUE_LIGHT`): `CUE_NONE=0`, `CUE_TEMP=1`, `CUE_WATER=2`, `CUE_HYDRO=4`,
`CUE_PHOTO=8`, **`CUE_LIGHT=16`** (new — leaf-exchange driver). `PHEN_*` are removed.

```fortran
type :: pheno_env_t                 ! raw daily environmental drivers only (read-only) — NO leaf status
   real(wp)    :: temp_day       = 0.0_wp   ! [K]  daily-mean air/canopy temp (thermal sums)
   real(wp)    :: soil_temp      = 0.0_wp   ! [K]  shallow soil temp (autumn cold-drop)
   real(wp)    :: avail_water    = 0.0_wp   ! [-]|[MPa] soil water (CUE_WATER)
   real(wp)    :: dmax_leaf_psi  = 0.0_wp   ! [MPa,<=0] DAILY-MAX leaf psi (CUE_HYDRO; ED2 dmax_leaf_psi,
                                            !           NOT instantaneous — a fast-loop daily reduction)
   real(wp)    :: rad            = 0.0_wp   ! [W/m2] daily-mean radiation (CUE_LIGHT; SW or PAR)
   real(wp)    :: daylength      = 12.0_wp  ! [h]  photoperiod (CUE_PHOTO / autumn short-day)
   integer(ik) :: doy            = 1_ik
   logical     :: hemis_north    = .true.
end type

type :: pheno_state_t               ! the TWO governors (prognostic) + cue feedstock accumulators
   real(wp) :: flush_drive   = 1.0_wp   ! [-] smoothed flush permission in [0,1]; 0 => dormant
   real(wp) :: shed_drive    = 0.0_wp   ! [-] smoothed active-shed pressure in [0,1]; 0 => no active shed
   real(wp) :: gdd           = 0.0_wp   ! [K day]   (CUE_TEMP)
   real(wp) :: chill         = 0.0_wp   ! [day]     (CUE_TEMP)
   real(wp) :: water_avg     = 0.0_wp   ! [-]|[MPa] running mean (CUE_WATER)
   real(wp) :: low_psi_days  = 0.0_wp   ! [day] consecutive dry, dmax_leaf_psi<ψ_tlp        (CUE_HYDRO)
   real(wp) :: high_psi_days = 0.0_wp   ! [day] consecutive wet, dmax_leaf_psi>=½ψ_tlp       (CUE_HYDRO)
   real(wp) :: light_avg     = 0.0_wp   ! [W/m2] running-mean radiation (CUE_LIGHT; ED2 rad_avg)
end type

type :: pheno_out_t                 ! RELATIVE rate tendencies [1/day] + a diagnostic
   real(wp)    :: leaf_flush_rate = 0.0_wp   ! [1/day] relative flush tendency; 0 == dormancy
   real(wp)    :: leaf_shed_rate  = 0.0_wp   ! [1/day] relative ACTIVE-shed tendency; 0 == no active shed
   integer(ik) :: cue_limiting    = CUE_NONE ! diagnostic: strongest active shed cue (argmax)
end type
```

`pheno_params_t`: **split** the single `cue_mask` into **`flush_cue_mask`** and **`shed_cue_mask`** (each
side takes its combine over its own mask — §4, the mechanism for the four target patterns); keep
`cue_sharpness` and all thermal/water/hydro/photo cue fields; **add** `k_flush_max [1/day]`,
`k_shed_max [1/day]` (per-PFT max relative rates), `tau_flush [day]`, `tau_shed [day]`, `water_width`, and
the `CUE_LIGHT` fields `light_on_threshold [W/m2]`, `light_width [W/m2]`, `light_window [day]`; **remove**
`phen_on_threshold`, `phen_off_threshold` (banding is gone). There is **no** baseline-turnover param here —
baseline leaf turnover stays a *carbon* trait (`leaf_turnover_rate`) applied by `tissue_turnover_rates` in
`plant_carbon_dynamics`. **Leaf-exchange flush option** (pattern 4): either a fixed high `k_flush_max`
(refills as fast as the light shed removes ⇒ canopy stays full — the recommended default), or an optional
`flush_matches_shed` flag that sets `leaf_flush_rate = leaf_shed_rate` for an exactly balanced exchange.

### Units — relative, per day
Both output rates and the `k_·_max` params are **relative** and in **`[1/day]`** (per locked decision 3).
A rate of `1/N` day⁻¹ means the flux would traverse a full canopy in `N` days (the user's "N days
bare↔full"). The kernel's internal `dt` is in days and `tau_·` are in days — all phenology arithmetic is
day-based. `plant_carbon_dynamics` multiplies by `dt_day` (= `dt_slow / day_sec`) when it forms the
absolute carbon flux; the pre-existing *carbon* turnover path keeps its own `[1/yr]`×`dt_yr` convention
(the two never mix — see §5.2).

### Boundary semantics (the two values the user specified)
- `leaf_flush_rate == 0`  ⇔  `flush_drive == 0`  ⇔  **dormancy** (no active flushing; canopy is not built).
- `leaf_shed_rate == 0`   ⇔  `shed_drive == 0`   ⇔  **no active shedding**. Baseline leaf turnover still
  removes leaves in `plant_carbon_dynamics` (a separate, temperature-live, replaceable channel) — phenology
  simply commands no *extra* senescence. Evergreen sits here permanently and is bit-identical to today.

### The rate map lives in the kernel
`pheno_out_t` carries the assembled relative rates (`k_·_max · drive`), computed inside `phenology_kernel`
from the governor drives and the per-PFT `k_·_max`. There is no separate baseline term in the phenology
output (locked decision 1) — the earlier "drives-only + downstream floor" split is dropped because the
floor is no longer phenology's concern.

---

## 4. The kernel algorithm (two governors → two relative rates)

Signature **unchanged** (drop-in): `pure subroutine phenology_kernel(env, params, dt, state, out)`,
`dt` in days. Four steps; environment-cues only.

**(1) Accumulate cue memory — as today's `accumulate()` plus `light_avg`:** season-gated GDD/chill (reset
outside their half-years), exponential `water_avg` running mean, consecutive `low/high_psi_days` counters
against `dmax_leaf_psi` (reset to 0 on a single opposite day), and the exponential `light_avg` running mean
(weight `dt/light_window`, ED2 `rad_avg`). The consecutive-day counters are **not** merged into signed
counters — the reset robustness ED2 scheme 5 relies on would be lost.

**(2) Per-cue `(s_flush, s_shed)`, each ∈ [0,1]** — each cue supplies a flush signal and a shed signal;
step (3) selects which side each cue feeds via the two masks:
- **CUE_TEMP:** `s_flush = logistic(k·(gdd − (phen_a + phen_b·exp(phen_c·chill)))/gdd_width)`;
  `s_shed = max(g_day·g_st1, g_st2)`, `g_day = logistic(k·(cold_drop_daylength − daylength)/daylen_width)`,
  `g_st1/2 = logistic(k·(cold_drop_soiltemp1/2 − soil_temp)/soiltemp_width)`. (Old `favor_temp`
  `f = f_gdd·(1 − drop)`, un-multiplied: `flush = f_gdd`, `shed = drop`.)
- **CUE_WATER:** `s_flush = logistic(k·(water_avg − water_on_threshold)/max(water_width,tiny))`;
  `s_shed = logistic(k·(water_off_threshold − water_avg)/max(water_width,tiny))`. `on > off` gives a
  hold band `(off, on)` where both ≈ 0 — hysteresis, no latch.
- **CUE_HYDRO:** `s_flush = clamp01(high_psi_days/max(high_psi_threshold,tiny))`;
  `s_shed = clamp01(low_psi_days/max(low_psi_threshold,tiny))` — ED2 scheme-5's consecutive-day model over
  the **daily-max** leaf psi, smoothed.
- **CUE_PHOTO:** `f_photo = logistic(photo_slope·(daylength − photo_crit))` **multiplies** the flush side
  (spring gate); short days also enter cold `s_shed` via `g_day`. (A modifier, not a standalone side.)
- **CUE_LIGHT (new):** `s_shed = logistic(k·(light_avg − light_on_threshold)/max(light_width,tiny))` — shed
  **rises with light** (leaf-exchange); `s_flush = 1` (non-limiting, so light never *blocks* flush). This is
  the leaf-DISPLAY analogue of ED2's Kim-2012 light phenology (whose leaf-QUALITY half stays out of scope).

**(3) Combine over the two masks (asymmetric):**
`s_flush = MIN over the cues in flush_cue_mask` (**conjunctive** — build canopy only when *every* flush cue
is clear; empty mask ⇒ `s_flush = 1`, permissive);
`s_shed = MAX over the cues in shed_cue_mask` (**disjunctive** — *any* shed cue commands senescence; empty
mask ⇒ `s_shed = 0`, no active shed). `cue_limiting = argmax` shed cue. **Splitting the mask is what
expresses the four patterns** (§1a, §6.1): a cue can drive one side without the other — evergreen puts TEMP
in `flush_cue_mask` but leaves `shed_cue_mask` empty (seasonal flush, no active shed); a leaf-exchanger
leaves `flush_cue_mask` empty (permissive high flush) and puts `CUE_LIGHT` in `shed_cue_mask`. **No explicit
mutual-exclusion term is needed** (locked decision 4): a PFT that must not flush under a given stress simply
lists that cue in `flush_cue_mask` too, so the `min` drives `s_flush → 0` there (a drought-deciduous PFT has
`CUE_WATER`/`CUE_HYDRO` in both masks ⇒ drought kills flush while raising shed).

**(4) Integrate the two governors (low-pass) and emit relative rates:**
```
w_f = dt/max(tau_flush, dt);  flush_drive = clamp01(flush_drive + w_f*(s_flush - flush_drive))
w_s = dt/max(tau_shed , dt);  shed_drive  = clamp01(shed_drive  + w_s*(s_shed  - shed_drive))
out%leaf_flush_rate = params%k_flush_max * flush_drive     ! [1/day]
out%leaf_shed_rate  = params%k_shed_max  * shed_drive       ! [1/day]
```
The guarded weight `dt/max(tau,dt)` caps at 1 (no overshoot when `tau < dt`) and is FPE-safe. `tau ≈ dt`
nearly disables the second filter. The governors are the "two internal pheno_state (one for flushing, one
for shedding)" of the original request.

**The four patterns in the kernel** (params in §6.1): **(1) evergreen** — `flush_cue_mask={TEMP}` or `{}`
(⇒ `flush_drive→` seasonal or 1), `shed_cue_mask={}` (⇒ `shed_drive=0`, `leaf_shed_rate=0`). **(2) temperate
deciduous** — both masks `{TEMP(+PHOTO)}`: spring `flush_drive→1`, autumn `shed_drive→1`. **(3) facultative
drought-deciduous** — both masks `{HYDRO}`: `high_psi_days` drives flush, `low_psi_days` drives shed (ED2
scheme 5). **(4) light-driven leaf-exchanging** — `flush_cue_mask={}` (permissive ⇒ `leaf_flush_rate =
k_flush_max`, high), `shed_cue_mask={LIGHT}` (⇒ `leaf_shed_rate` rises with `light_avg`); the sustained high
flush refills what the light shed removes, so the canopy stays ~full while turning over — the interplay
completes in `plant_carbon_dynamics` (§5.3), at a steady leaf pool.

**GPU/pure compliance:** `pure`, scalar-only, arithmetic + intrinsics, fixed-size (no array temps ⇒
issue-#7 N/A). `safe_exp` inside `logistic`; every divide guarded by `max(width,tiny)` **not** `merge`
(`merge` evaluates both arms ⇒ divide-by-zero on the dead arm under `-fpe0`). Integer/logical `iand` cue
tests and counter resets stay `if`. Phenology runs in the daily host driver loop (no `target` region) but
stays offload-clean; the nvfortran-multicore gate applies to `meds_phenology` + `meds_numerics` +
`meds_time` + the test.

---

## 5. `plant_carbon_dynamics` — the actual leaf/storage update (single authority)

All leaf/storage carbon movement lives here (locked decision 1). The phenology rates are pure relative
tendencies; this section turns them into carbon and mutates the pools (via the MEDS pattern: the kernel
returns per-pool NPP, the core engine applies it). Authored from the verified seam
(`get_plant_flux_slow` in `meds_plant_interface`; `tissue_turnover_rates` + `plant_carbon_allocation` in
`meds_plant_carbon_dynamics`).

### 5.1 Two loss channels — baseline (replaceable) vs active shed (non-replaceable) (⚑ blocker 1)
A single lumped shed flux would be wrong: `plant_carbon_allocation` P1 funds `loss%leaf` via
`fund(loss%leaf, allow_store=.true.)` **unconditionally**, so a large shed would be immediately refilled
from NPP/storage — a deciduous canopy would drop-then-replace forever, draining storage, never shedding
(evergreen/leaf-exchange mask this because continuous replacement is *correct* for them). The two-channel
split fixes it — and here it is **natural**, because the two rates come from different owners:

| Channel | Flux | P1/P2 refill? | Owner |
|---|---|---|---|
| baseline turnover | `rates%leaf(tissue_temp) · leaf_carbon · dt_yr` | **yes** (replaceable) | `tissue_turnover_rates` — a *carbon* trait, temperature-live, **unchanged** ⇒ evergreen bit-identical |
| active shed | `min(leaf_shed_rate · leaf_carbon_full · dt_day, leaf_carbon)` | **no** — never refilled | phenology `leaf_shed_rate` |

Because the active channel is non-replaceable and (§5.3) linear-in-`leaf_carbon_full`, it drives the
canopy down even if baseline replacement is still running; the flush cap (§5.3) governs only growth toward
the target, so it cannot re-grow the shed carbon within the step.

### 5.2 Baseline turnover stays temperature-live and in the carbon layer (⚑ major)
`tissue_turnover_rates` keeps computing `rates%leaf = leaf_turnover_rate · f_t`,
`f_t = 1/(1 + safe_exp(evg_slope·(evg_ref_temp − tissue_temp)))` — a function of the **live** `tissue_temp`.
Phenology never sees it (locked decision 1). This removes the whole "temperature-live floor baked into a
static param" hazard the reviewers flagged: the floor was never in the phenology output to begin with.
(Note the unit boundary: baseline turnover keeps `[1/yr]`×`dt_yr`; the phenology-driven flush/shed use
`[1/day]`×`dt_day`. They are separate terms and never mix.)

### 5.3 Flush fill + active shed (linear, clamped, snap-to-bare) (locked decision 5, ⚑ major)
Inside `plant_carbon_dynamics` (the flush cap moves here from the driver, per locked decision 1):
```
leaf_carbon_full = size2leaf_carbon(dbh, h, sla)        ! full-canopy allometric leaf C (= leaf_target)
! flush: cap the leaf growth demand at the relative flush rate (linear toward full)
demand%leaf   = min(leaf_carbon_full - leaf_carbon,  leaf_flush_rate * leaf_carbon_full * dt_day)   ! >=0
demand%fineroot = min(root_to_leaf_ratio*leaf_carbon_full - fineroot_carbon,
                      leaf_flush_rate * root_to_leaf_ratio*leaf_carbon_full * dt_day)               ! >=0
can_flush     = (leaf_flush_rate > 0.0_wp)              ! replaces the PHEN_ON boolean
! active shed: linear constant-absolute decline, clamped to the pool, snapped to bare
shed_linear   = leaf_shed_rate * leaf_carbon_full * dt_day
loss%leaf_shed = min(shed_linear, leaf_carbon - loss%leaf_baseline)                                ! >=0
if (leaf_carbon - loss%leaf_baseline - loss%leaf_shed < leaf_carbon_min .and. leaf_shed_rate > 0)  &
     loss%leaf_shed = leaf_carbon - loss%leaf_baseline  ! snap remaining canopy to true bare
```
The linear (constant-absolute) form empties a full canopy in exactly `1/leaf_shed_rate` days with **no
exponential tail**; the pool clamp protects the large-`k` (near-instantaneous ED2) limit from driving
`leaf_carbon`/LAI negative; the `leaf_carbon_min` snap gives ED2's fully-abscised `−2` (no residual winter
canopy for RT/energy/root sink). The **baseline** channel stays proportional (`rates%leaf·leaf_carbon·dt_yr`)
— do **not** convert it to the linear form (that would perturb the evergreen regression). `elongf`, if
needed for a bound or a diagnostic, is computed **here** as `leaf_carbon/leaf_carbon_full` — never in
phenology (locked decision 2).

### 5.4 Retained fraction — correct closure (⚑ major; P4)
When `retained_carbon_fraction` (a carbon trait) is introduced, split the **shed** carbon (not baseline):
```
npp%leaf          = a_leaf - (loss%leaf_baseline + loss%leaf_shed)          ! FULL removal
npp%nonstructural = a_store - draw + retained_carbon_fraction * loss%leaf_shed
litter (untracked)= (1 - retained_carbon_fraction) * loss%leaf_shed + loss%leaf_baseline
```
Assert per call the testable half: `leaf_delta + storage_delta + wood + repro + litter == net_carbon`
(leaf+storage+wood+repro closes once the untracked outflow is accounted). The synthesizer's "credit
retained to storage, only litter as negative `npp%leaf`" **creates** carbon — avoid it.

### 5.5 Litter — untracked outflow, as today (⚑ blocker 2; locked decision 6)
**There is no vegetation→litter sink in MEDS today:** `compute_slow_derivatives`
(`meds_vegetation_dynamics.f90:141`) does `lc_new = max(leaf_carbon + npp%leaf, 0)` and the lost carbon
**vanishes**; `meds_soil_biogeochem`'s `build_litter_input`/`soil_carbon_step` exist but are **never
called** by any driver. Baseline turnover carbon is *already* unconserved today. So the refactor promises
only *"leaf + storage closes; shed/turnover carbon leaves untracked, exactly as today."* A diagnostic
per-patch litter accumulator (`kgC/m²`) and the full `build_litter_input` wiring are **reserved
follow-ups**, out of scope (they belong with the demography→litter→Rh driver seam `meds_soil_biogeochem`
still lacks).

### 5.6 Phenology-off fallback (behavior-preserving)
When `cfg%phenology_on == .false.`, `plant_carbon_dynamics` uses `leaf_flush_rate` = large (uncapped
fill), `leaf_shed_rate` = 0 (no active shed; baseline turnover unchanged), reproducing today's behavior
**bit-identically**. Use `if/else`, **not** `merge` (⚑ `merge` evaluates both arms; a rate-compute arm
may divide and would trip `-fpe0` when phenology is off).

### 5.7 Fine-root / root phenology
Phenology stays **strictly leaf** (locked decision 1). ED2 ties fine-root senescence to leaf via
`root_phen_factor` (scheme 5); if adopted, `plant_carbon_dynamics` (not phenology) applies it to the
fine-root loss. Reserve a `root_phen_factor` param.

---

## 6. Coverage

### 6.1 The four target patterns — parameterization (acceptance criteria, §1a)
`k_turn = leaf_turnover_rate` is the *carbon-layer* baseline (cold-suppressed via `f_t`, §5.2), **not** a
phenology output. `{}` = empty mask (permissive flush / no active shed).

| # | Pattern | `flush_cue_mask` | `shed_cue_mask` | key params | resulting canopy |
|---|---|---|---|---|---|
| 1 | temperate evergreen | `{}` (high flush) or `{TEMP}` (seasonal) | `{}` | `k_flush_max≈1/15 d⁻¹`; `k_turn` low + cold-suppressed | full year-round; winter hold from cold-suppressed `k_turn` |
| 2 | temperate deciduous | `{TEMP, PHOTO}` | `{TEMP, PHOTO}` | `phen_a/b/c`, `cold_drop_*`; `k_flush_max≈1/15`, `k_shed_max≈1/20 d⁻¹` | full in season, bare in winter |
| 3 | facultative drought-decid | `{HYDRO}` | `{HYDRO}` | `leaf_psi_tlp`, `low/high_psi_threshold`; `k_flush_max≈leaf_grow_rate`, `k_shed_max≈leaf_shed_rate` | evergreen when watered, sheds on drought (ED2 sched 5) |
| 4 | light-driven leaf-exchanging | `{}` (high flush) | `{LIGHT}` | `light_on_threshold`, `light_width`, `light_window`; `k_flush_max` high (≈1/15) ≥ `k_shed_max` | ~full while turning over fast under high light |

Pattern 4's flush must keep up with the light shed: either a fixed high `k_flush_max` (default) or the
optional `flush_matches_shed` flag (`leaf_flush_rate = leaf_shed_rate`, exactly balanced exchange).

### 6.2 ED2 habit → MEDS mapping (still covered)
ED2's two orthogonal selectors — the global `iphen_scheme` (driver path) and the per-PFT `phenology(ipft)`
(habit) — collapse into the per-PFT masks + rate params. `flush`/`shed` are the two **relative** outputs.

| ED2 habit `phenology(ipft)` | flush_cue_mask / shed_cue_mask | leaf_flush_rate `[1/day]` | leaf_shed_rate `[1/day]` |
|---|---|---|---|
| 0 evergreen | `{}` (or `{TEMP}`) / `{}` | `k_flush_max` while pool < full | 0 (baseline `k_turn` only) |
| 1 drought-decid (instant paw) | `{WATER}` / `{WATER}` (window≈dt) | `k_flush_max` on wet recovery | large `k_shed_max` on dry (instant = large-`k` + snap) |
| 2 cold-decid (Botta/White) | `{TEMP,PHOTO}` / `{TEMP,PHOTO}` | spring pulse when `gdd`>`a+b·exp(c·chill)`; 0 winter | `k_shed_max` in autumn cold-drop; 0 in season |
| 3 light-decid + drought | `{WATER}` / `{WATER}` | leaf-**display** half = habit 4 | = habit 4 |
| 4 drought-decid (10-day paw) | `{WATER}` / `{WATER}` (window=10) | ramps with running-mean recovery | `k_shed_max·shed_drive` tracks dry-down |
| 5 hydraulic decid (**the model**) | `{HYDRO}` / `{HYDRO}` | `k_flush_max ≈ leaf_grow_rate` (direct) | `k_shed_max ≈ leaf_shed_rate` (direct) |
| light-exchanging (Kim-2012 display half) | `{}` / `{LIGHT}` | high `k_flush_max` (permissive) | rises with `light_avg` — **pattern 4** |
| `iphen_scheme=1` prescribed | — | `src/forcing` prescribed-`elongf` overlay, **not** a cue |

Seed magnitudes: `k_flush_max` from ED2 `leaf_grow_rate` (~1/15–1/24 d⁻¹); `k_shed_max` from ED2
`leaf_shed_rate` (~1/10–1/20 d⁻¹); `k_turn = leaf_turnover_rate ≈ 1/leaf_lifespan` (carbon layer);
`CUE_LIGHT` thresholds from ED2 `radto_min/max` / `rad_avg`. Habit-3's leaf-**quality** half
(`update_turnover` SLA/Vcmax/llspan), prescribed phenology, and the N-cycle split leave cleanly.

---

## 7. Migration — state, driver, config, tests

### 7.1 Cohort state (`src/core/meds_core_state_types.f90`)
**Store the two governor accumulators** `pheno_flush_drive`, `pheno_shed_drive` (the user's "two internal
pheno_state"); the relative rates are re-derived where consumed (`rate = k_·_max · drive`, a trivial pure
helper in the plant lib) — nothing new is stored for the rates, and `elongf` is stored nowhere (it is a
diagnostic, locked decision 2). Net P1 change: **remove** `phenology_status(:)` (int) +
`PHENOLOGY_STATUS_INIT`; **add** `pheno_flush_drive(:)`, `pheno_shed_drive(:)` (real); **keep** `pheno_gdd`,
`pheno_chill`. The four cue columns `water_avg`/`low_psi_days`/`high_psi_days`/`light_avg` are added in
**P3** when the WATER/HYDRO/LIGHT drivers (soil water, daily-max leaf psi, radiation) are threaded.

Thread the two new columns through the **seven lockstep routines** (line map verified by read):

| Routine | What to touch |
|---|---|
| `cohort_alloc` | `allocate(...)` (~L275) **and** zero-init (~L279): `flush_drive=1.0`, `shed_drive=0.0` |
| `site_free` | dealloc list (~L254) |
| `cohort_ensure_capacity` | `1:m` copy block (~L355–357) |
| `move_alloc_block` | `move_alloc` pair (~L398–400) |
| `cohort_reorder` | permute (`perm`) (~L482–484); `cohort_compact` rides this |
| `copy_cohort_slot` | slot copy (~L540–542) |
| `init_cohort` | birth init (~L605–607): `flush_drive=1.0`, `shed_drive=0.0` (replaces `phenology_status = PHENOLOGY_STATUS_INIT`) |

Creation sites also stamping the new fields: `init_cohort`, `add_cohort`/`init_bare_ground`
(`src/init/meds_init.f90`), `apply_recruitment` + `split_cohorts` (`meds_core_cohort_fusefiss.f90`),
`apply_patch_disturbance` (`meds_core_patch_fusefiss.f90`). **Fusion:** `fuse_2_cohorts`
(`meds_core_cohort_fusefiss.f90` ~L163–220) is a **by-hand** merge that today does **not** blend
`pheno_gdd`/`chill`/`phenology_status` (implicit survivor-keep). Decide explicitly: nplant-weight the two
drives (add to the accumulator block ~L189–192) — recommended — vs survivor-keep (a no-op, but must be
stated). See §9.

### 7.2 Driver — fold `meds_phenology_driver` into `meds_vegetation_dynamics` (locked decision 8)
**Delete `src/driver/meds_phenology_driver.f90`** and move its two routines into
`meds_vegetation_dynamics` as the first step of the slow loop. Both files already live in `src/driver/`
(`libmeds_aux`) and `meds_vegetation_dynamics` already imports `meds_core_interface` (`site_t`),
`meds_plant_interface` (pheno types + `update_phenology`), `meds_config`, and `meds_allometry` — so the
fold adds **no** library edge and **no** CMake change (the `src/driver/*.f90` glob simply drops the deleted
file). Concretely:

- **`advance_leaf_phenology(site, cfg, doy)`** (the renamed `leaf_phenology`) + **`flatten_pheno_params`**
  become module procedures of `meds_vegetation_dynamics`. `advance_leaf_phenology` builds `pheno_env_t` from
  the site daily-mean temp + latitude/`doy` in P1 (and `dmax_leaf_psi` + running-mean `rad` in P3) **only**
  (no `size2leaf_carbon`, no `elongf` — locked decision 2), packs/unpacks
  `pheno_flush_drive`/`pheno_shed_drive`, and calls `update_phenology`. `flatten_pheno_params` gains
  `flush_cue_mask`, `shed_cue_mask`, `k_flush_max`, `k_shed_max`, `tau_flush`, `tau_shed`, `water_width`,
  and the `CUE_LIGHT` params; drop the `phen_on/off_threshold` packing. It touches **no** leaf/storage carbon.
- **`vegetation_dynamics` gains an (optional) `doy` argument** and calls `advance_leaf_phenology` as **step
  0** when `cfg%phenology_on` (before `carbon_growth`, so the drives are current when the carbon seam reads
  them). Ordering is **preserved**: the fast loop still runs first (it fills the daily-mean-temp accumulator),
  then `vegetation_dynamics` advances phenology, then `carbon_growth`. `doy` is optional so test callers that
  don't enable phenology are unchanged; `phenology_on` with absent `doy` is an `error stop`.
- **`meds_stepper.advance_one_step`** drops the `use meds_phenology_driver` import and the standalone
  `call leaf_phenology(...)` block; it computes `day_of_year(step_start)` (it already imports `meds_time`)
  and passes it as `vegetation_dynamics(..., doy=…)` when `phenology_on` (keeping the "phenology_on needs
  step_start" guard). This reinforces the module's charter: `carbon_growth` calls `get_plant_flux_slow` and
  now `advance_leaf_phenology` calls `update_phenology` — **both** plant-kernel calls sit in the one
  plant↔demography meeting-point module.
- **Keep `advance_leaf_phenology` public** so `test_phenology_driver` can drive it directly (without running
  the full demographic step). The relative rates still flow to `plant_carbon_dynamics` via the carbon seam
  (`carbon_growth` → `get_plant_flux_slow`), which reads the stored drives and derives the rates.

### 7.3 Config (`meds_pft_params.f90`, `meds_config.f90`, `meds_config_io.f90`)
A **compile dependency** across four sites:
1. `pft_table_t` field decls — **replace** `pheno_cue_mask` with `pheno_flush_cue_mask` +
   `pheno_shed_cue_mask`; add `pheno_k_flush_max`, `pheno_k_shed_max`, `pheno_tau_flush`, `pheno_tau_shed`,
   `pheno_water_width`, the `CUE_LIGHT` fields `pheno_light_on_threshold`/`pheno_light_width`/
   `pheno_light_window`, the optional `pheno_flush_matches_shed` flag (and `pheno_retained_carbon_fraction`
   in P4); drop `pheno_on_threshold`, `pheno_off_threshold`.
2. `alloc_pft_table` — the `allocate(...)` list **and** the default-install assignments.
3. `load_phenology_pft` — the `req_pa` calls (two masks now).
4. `validate_config` — **remove** the `off<on` check (`meds_config.f90` ~L344–346); **add**
   `k_flush_max>0`, `k_shed_max>=0` (per day; and `retained∈[0,1]` in P4); **replace** the single
   `cue_mask` range check with range checks on **both** masks; **keep** the **WATER/HYDRO/LIGHT rejection**
   for any cue whose fast-loop driver is not yet threaded — in P0–P2 only `TEMP`/`PHOTO` are permitted in
   either mask; P3 lifts the rejection as each driver (soil water, `dmax_leaf_psi`, radiation) lands.

### 7.4 Carbon-seam types (`meds_plant_types.f90`, `meds_plant_interface.f90`,
`meds_plant_carbon_dynamics.f90`)
Rewrite the PHENOLOGY section of `meds_plant_types` per §3; remove `PHEN_OFF/PHEN_DORMANT/PHEN_ON` (`CUE_*`
stay). `carbon_env_t`: replace `phenology_status(int)` with `leaf_flush_rate`, `leaf_shed_rate` (real,
`[1/day]`). `plant_carbon_allocation` gains the flush cap + active-shed split (§5.3) — its integer
`phenology_status` arg becomes the logical `can_flush` **plus** the new rate/`leaf_carbon_full`/`dt_day`
arguments (a public-seam break, §7.5). `meds_plant_interface`: drop the `PHEN_*` and `daylength`
re-exports (daylength now from `meds_time`); keep `pheno_*` types + `CUE_*`.

### 7.5 Tests to migrate (⚑ four files)
Removing `phenology_status`/`PHEN_*` and changing `plant_carbon_allocation`'s signature is a **public-seam
break** (`plant_carbon_allocation` re-exported at `meds_plant_interface.f90:59`). All four
phenology-touching test files must be rewritten or the P1/P2 "full suite green" gates fail:

| File | Change | Lands with |
|---|---|---|
| `test/test_plant_phenology.f90` | tri-state → relative-rate semantics (§8-below) | P0 |
| `test/test_phenology_driver.f90` | call `advance_leaf_phenology` from `meds_vegetation_dynamics` (folded, §7.2); assert `pheno_flush_drive`/`pheno_shed_drive`; drop `PHEN_*` | P1 |
| `test/test_carbon_growth.f90` | set `env%leaf_flush_rate`/`leaf_shed_rate` instead of `phenology_status` | P2 |
| `test/test_plant_carbon_dynamics.f90` | new `plant_carbon_allocation` signature (`can_flush` + rates); flush-cap + active-shed cases | P2 |

No IO/output module consumes `phenology_status` (grep-confirmed), so IO is unaffected.

---

## 8. Shared-util migration

Four helpers leave `meds_phenology` for the DAG root so respiration/hydrology/energy/forcing stop
re-implementing them. **DAG-legal** (verified): `shared/util` may not depend on `plant`; every move is a
`plant → shared` or intra-`shared` edge.

**Into `src/shared/util/meds_numerics.f90`** (add `use meds_constants, only: safe_exp` — safe:
`meds_constants` uses `meds_kinds` only, and `meds_time` already carries the identical edge; update the
module's "uses `meds_kinds` only" header note): `logistic(z) = 1/(1+safe_exp(-z))` and
`clamp01(x) = min(1,max(0,x))`, both `elemental pure`; add a general `clamp(x,lo,hi)` companion
(`adaptive_step_update` clamps inline today).

**Into `src/shared/util/meds_time.f90`** (has `day_of_year` + `solar_cosz` + `pi`; no new `use`):
`daylength(lat_deg, doy)` [h] beside `solar_cosz`, and `doy_effective(doy, hemis_north)` beside
`day_of_year`, both `elemental`.

**Declination reconciliation — deferred (recommended).** `meds_time.solar_cosz` uses Cooper-1969
`23.45·sin(2π(284+doy)/365)`; `meds_phenology.daylength` uses White-1997 `−23.44·cos(2π(doy+9)/365)`. Same
curve to within ~1.25-day phase, but switching `daylength` onto Cooper **changes** its outputs and forces
re-baselining the daylength/cold-drop-timing golden values — a *science* choice, not a mechanical move.
**Keep White-1997 verbatim** (polar-day fix `arg<=-1 ⇒ 24h` preserved); add a header comment that
`meds_time` now hosts two standard declination approximations; **defer** unification (a single private
`solar_declination(doy)` both call) as a follow-up (§9). Keeps the migration bit-behavior-identical.

**Edits:** `meds_phenology` deletes the four helper bodies, `use meds_constants, only: safe_exp`, adds
`use meds_time, only: doy_effective` + `use meds_numerics, only: logistic, clamp01`, drops
`public :: logistic, daylength`. `meds_plant_interface` drops the `daylength` re-export.
`meds_vegetation_dynamics` (the folded phenology advance, §7.2) + `test/test_plant_phenology` import
`daylength` from `meds_time` (an astronomical function stops masquerading as a plant symbol — a DAG
improvement). **CMake: none** — both destinations are
existing files under the `src/shared/*.f90` glob; no new file, no link edge (`meds_plant` already links
`meds_shared` PUBLIC). Cross-domain dedup (`logistic` in `meds_plant_carbon_dynamics:60`; `clamp01` across
forcing/biophysics/biogeochemistry) is a **separate** P5 sweep; do **not** touch `meds_hydro_curve`'s
`1/(1+r**kexp)` or `meds_temp_response` denominators (they look similar, are not `logistic`/`clamp01`).

### Test plan — the four target patterns are the acceptance tests
Rewrite `test/test_plant_phenology.f90` to **relative-rate** semantics, gating on ifx (`-stand f18 -check
all -fpe0`) **and** nvfortran multicore. The four §1a patterns are the headline cases:
(1) **temperate evergreen** (`flush={}`/`shed={}`, or `flush={TEMP}`): `leaf_flush_rate = k_flush_max`
(or seasonal), `leaf_shed_rate == 0` always. (2) **temperate deciduous** (`flush={TEMP,PHOTO}` /
`shed={TEMP,PHOTO}`): spring flush pulse once `gdd ≥ a+b·exp(c·chill)`, autumn `leaf_shed_rate = k_shed_max`
cold-drop, both 0 in dormancy; more chilling lowers the heat requirement; ED2-timing check in the
sharp-slope limit. (3) **facultative drought-deciduous** (`flush={HYDRO}` / `shed={HYDRO}`): reproduces ED2
scheme-5 `leaf_grow_rate`/`leaf_shed_rate` over a `dmax_leaf_psi` series; stays evergreen (both rates flat)
when `dmax_leaf_psi` never crosses `ψ_tlp`. (4) **light-driven leaf-exchanging** (`flush={}` /
`shed={LIGHT}`): `leaf_flush_rate = k_flush_max` (high, ~constant) while `leaf_shed_rate` **rises with
`light_avg`** — assert both > 0 under high light, and (in the P2 carbon test) that the canopy stays ~full
(leaf_carbon near target) with nonzero turnover. Plus: (5) **hysteresis** — a driver oscillating across the
water/light hold-band shows no flicker (`tau` + hold-band). (6) **boundary values** —
`leaf_flush_rate==0 ⇔ flush_drive==0`; `leaf_shed_rate==0 ⇔ shed_drive==0`. (7) **range/FPE** — rates finite
and non-negative under degenerate drivers; no NaN/trap under `-fpe0`. The **carbon** tests (flush cap,
linear active shed, `leaf_carbon_min` snap, `k_shed_max·dt_day>1` pool-clamp, bare→full→bare leaf+storage
conservation, pattern-4 full-canopy-under-shed, phenology-off bit-identical) live with P2 in
`test_plant_carbon_dynamics`/`test_carbon_growth`.

---

## 9. Phasing

Each phase is independently buildable + testable and gates on **ifx** (`-stand f18 -check all -fpe0`)
**and** **nvfortran multicore** (a green ifx run is not sufficient — the whole-program-optimizer trap).

- **P0 — Signal kernel + shared-util relocation** (standalone, no demography wiring). Rewrite the PHENOLOGY
  section of `meds_plant_types` (env = cues only, no `elongf`; state = two drives + `gdd`/`chill`; params
  rate fields, drop banding thresholds; out = two relative rates). Rewrite `phenology_kernel` (per-cue
  `(s_flush,s_shed)`, min/max combine, low-pass governors, `rate = k·drive` — pure/FPE/GPU-clean). Relocate
  `logistic`/`clamp01` → `meds_numerics`, `daylength`/`doy_effective` → `meds_time` (behavior-identical;
  keep White-1997). Rewrite `test_plant_phenology` to relative-rate semantics. **Gate:** standalone
  `meds_phenology` target builds; both compilers green.
- **P1 — State + driver fold + config.** Replace `phenology_status(int)` with
  `pheno_flush_drive`/`pheno_shed_drive(real)` across the 7 lockstep + 5 creation sites; decide
  `fuse_2_cohorts` blend policy. **Fold `meds_phenology_driver.f90` into `meds_vegetation_dynamics`** (§7.2):
  move `advance_leaf_phenology` + `flatten_pheno_params` in, add the optional `doy` arg to
  `vegetation_dynamics` (phenology as step 0), update `meds_stepper` to drop the separate call and pass `doy`,
  delete the file. `advance_leaf_phenology` builds a cues-only env (no `elongf`) and packs/unpacks the drives.
  Add the PFT keys `[1/day]` (two masks), drop thresholds, update `alloc_pft_table` + `validate_config`.
  Rewrite `test_phenology_driver` against the folded routine. Still TEMP/PHOTO only. **Gate:** full suite
  green; fusion-conservation test on the drives.
- **P2 — `plant_carbon_dynamics` applies the rates (single authority).** `carbon_env_t`: two rates replace
  `phenology_status`. Move the flush cap into `plant_carbon_allocation` (linear toward full); add the active
  shed channel (linear-in-`leaf_carbon_full`, non-replaceable, pool-clamped, `leaf_carbon_min` snap);
  baseline turnover stays from `tissue_turnover_rates(tissue_temp)`. `can_flush = (leaf_flush_rate>0)`;
  `if/else` phenology-off fallback (bit-identical). `elongf` is diagnosed here where needed. Rewrite
  `test_carbon_growth` + `test_plant_carbon_dynamics`. **Gate:** bare→full→bare **leaf+storage(+wood+repro)**
  conservation with litter as the computed untracked outflow; phenology-off bit-identical regression;
  linear reach-bare (`leaf_carbon_min` snap, exact `1/leaf_shed_rate`-day traversal); `k_shed_max·dt_day>1`
  pool-clamp test.
- **P3 — WATER + HYDRO + LIGHT cues wired (completes the four patterns).** Add
  `water_avg`/`low_psi_days`/`high_psi_days`/`light_avg` cohort columns (7 lockstep + creation sites again);
  thread the soil-water running mean, the **daily-max** leaf psi `dmax_leaf_psi` (a fast-loop daily
  reduction), and the running-mean radiation `rad_avg` from the fast loop into the daily phenology
  accumulators; move `soil_temp` off the daily-air-temp v1 proxy; lift the `validate_config`
  WATER/HYDRO/LIGHT rejection as each driver lands. Enables **pattern 3** (facultative drought-deciduous,
  `CUE_HYDRO`) and **pattern 4** (light-driven leaf-exchanging, `CUE_LIGHT`). **Gate:** all four §1a patterns
  reproduced — drought/hydraulic timing vs ED2 scheme 4/5; leaf-exchanger holds a full canopy under rising
  `light_avg` with nonzero turnover.
- **P4 — Retained-fraction (litter stays untracked).** New `retained_carbon_fraction` *carbon* trait; split
  active-shed carbon retained→`nonstructural` / rest→untracked litter outflow with the **full-removal**
  closure (§5.4). Optional fine-root `root_phen_factor` coupling (in `plant_carbon_dynamics`; phenology stays
  leaf-only). **Gate:** `leaf+storage+wood+repro+outflow == net_carbon` closes to <1% per call. (Diagnostic
  litter accumulator + `build_litter_input` wiring are reserved follow-ups, out of scope — §5.5.)
- **P5 — Cross-domain dedup + declination unification (follow-up, separate PR).** `logistic`/`clamp01`
  call-site sweep; optional single private `solar_declination` shared by `solar_cosz` + `daylength` with a
  deliberate daylength-golden re-baseline.

---

## 10. Open questions

1. **Fusion blend.** nplant-weight the two drives (recommended; and, for consistency, `gdd`/`chill` — a
   minor change from today's survivor-keep), or blend *only* the new drives and leave `gdd`/`chill`
   survivor-kept?
2. **Evergreen path.** Keep the `flush={} , shed={}` (empty-mask) `cycle` in `advance_leaf_phenology`
   (correctness relies silently on the init `flush_drive=1`/`shed_drive=0` seed at every creation site), or
   drop the cycle so the kernel actively maintains the fixed point?
2b. **Pattern-4 flush.** Default to a fixed high `k_flush_max` (refills faster than the light shed removes ⇒
   canopy stays ~full — recommended), or wire the optional `flush_matches_shed` flag
   (`leaf_flush_rate = leaf_shed_rate`, exactly balanced exchange)? And should `CUE_LIGHT` use a running-mean
   radiation (ED2 `rad_avg`, `light_window ~10 d`) or instantaneous daily radiation?
3. **`tau_flush`/`tau_shed` defaults.** A fixed few-day anti-flicker timescale (keeps the dominant
   timescale = `1/k_flush_max`), or derive `tau` from the rate so one param controls both smoothness and
   traverse speed?
4. **Fine-root coupling.** Phenology strictly leaf-only + `plant_carbon_dynamics` applies `root_phen_factor`
   (recommended), or also emit a `fineroot_shed_rate` from the module? ED2 scheme 5 couples them.
5. **Seed parameters.** Which `k_flush_max`/`k_shed_max` `[1/day]` and `k_turn` `[1/yr]` values (ED2
   `init_pft_phen_params` vs FATES) for the four exemplar habits (cold-decid, drought-decid, leaf-exchanger,
   evergreen) across the three MEDS demonstrator PFTs?
6. **Declination.** Defer unification (recommended, no golden re-baseline), or unify `solar_cosz` +
   `daylength` onto one `solar_declination(doy)` now?

---

## 11. References
- **ED2** `dynamics/phenology_driv.f90` (per-PFT `phenology(ipft)` habit switch; scheme 5 = the per-day
  `leaf_grow_rate`/`leaf_shed_rate` model), `phenology_aux.f90` (`update_thermal_sums`, `daylength`,
  `update_turnover`, `pheninit_balive_bstorage`), `memory/phenology_coms.f90`
  (`retained_carbon_fraction`, `root_phen_factor`, `elongf_min/flush`, Botta/White params).
- **MEDS** `MEDS_PHENOLOGY_DESIGN.md` (the superseded tri-state design),
  `MEDS_PLANT_CARBON_DYNAMICS_DESIGN.md` (PARTEH allocation ladder), `CLAUDE.md` (DAG wall, lockstep
  discipline, issue-#7 nvfortran trap).
- **FATES** `EDPhysiologyMod.F90`; **CLM5** `CNPhenologyMod`; **Chuine (2000)** budburst; **Delpierre
  (2009)** autumn CDD×photoperiod; **Xu et al. (2016)** turgor-loss-point shedding.
