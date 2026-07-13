# MEDS FAST (sub-daily) output tier — implementation plan

**Status:** IMPLEMENTED 2026-07-13 (P0+P1+P2 combined). ifx Debug 33/33 + nvfortran multicore 33/33;
FAST netCDF output cross-checked ifx-vs-nvfortran to roundoff (worst 3e-11) incl. per-cohort/soil slabs;
diurnal cycle verified on 2024 forcing (SW 0→744, Rnet/LE/H/Tcas). Adversarial review NOT yet run (pending
user decision). One deviation from the plan below: **h_flux is surfaced as an out-arg of `column_fast_step`
(next to `le_flux`), not computed in the fast loop** — the post-call `aero`/`aenv` read miscompiled to 0.0
on nvfortran (issue #7 class); computing it in-call (where `le_flux` is) fixed it.
**Goal:** emit an opt-in, per-`dt_fast` netCDF stream so a user can analyze **diurnal cycles** of
site-level canopy/soil/energy/water state directly from a MEDS run (CAS temperature, GPP, ET, soil
temperature, and — from P1 — the sensible/net-radiation/incident-SW/ustar fluxes).
**Scope:** MEDS codebase only. Site-scalar (and, P1, soil-column) variables only — **no** per-cohort /
per-patch sub-daily output (see §8, out of scope).

---

## 1. Motivation & what "FAST output" means

The diagnostic-aggregation subsystem (`archive/MEDS_IO_DESIGN.md`) already emits DAILY / MONTHLY / ANNUAL
tiers. The **FAST tier** (`FREQ_FAST`, tier index 1) is modeled everywhere in the data structures and
config but is **never produced** — it was deferred at P0 of the IO work. Diurnal-cycle analysis (the
shape of the day: dawn/midday/dusk temperature, the GPP and ET curve, soil-surface thermal wave) needs
output at the fast-loop cadence (`dt_fast`, ~30 min), which no current tier provides. A bespoke opt-in
CSV probe (`[fast].fast_probe` → `write_fast_probe`) was added as a stopgap during the integrator
evaluation; this plan delivers the **proper netCDF FAST tier** through the existing aggregation pipeline.

## 2. What already exists (do NOT rebuild)

The workflow read the whole subsystem. The FAST tier is far more provisioned than "deferred" implies:

- **Tier model** — `FREQ_FAST = 1`, `N_FREQ = 4`, `freq_tier_index(FREQ_FAST)=1`, `freq_letter='F'`
  (`meds_output_config.f90:30,79,94`). `freq_on` default `[.false.,.true.,.true.,.true.]` — FAST is the
  only tier off by default.
- **Config surface already parsed** — `[output.fast].enabled` → `freq_on(1)` (default false),
  `[output.fast].file_chunk` → `file_chunk(1)` (default `'day'`), and `[output].fast_interval_steps`
  (default 4) → `cfg%output%fast_interval_steps` → **`mgr%fast_interval_steps`**
  (`meds_config_io.f90:355,364,368`; `meds_output_registry.f90:274`). It is threaded end-to-end but
  **never consumed** — nothing closes the FAST window on it yet.
- **Storage provisioned** — `manager_alloc_buffers` loops `do t=1,N_FREQ` and allocates `buf(:,1)`,
  `pending(1)`, `stream(1)` whenever any registered var carries `FREQ_FAST`
  (`meds_output_registry.f90:279+`). No new allocation is needed for FAST.
- **Reduction engine is tier-agnostic** — `integrate_scalar` / `integrate_slab` (fold),
  `normalize_scalar` / `normalize_slab` (close), `reset_buffer`/`seed_of` (re-seed), and `close_tier`
  (`meds_output_integrate.f90:139,155,196,114,358`) all work on any tier index; `close_tier(mgr,1)`
  already yields `FREQ_FAST`. **Reuse them verbatim.**
- **Serializer is already sub-daily-capable** — `output_serialize_pending` drains `pending(t)` for
  `do t=1,N_FREQ` (`meds_output_manager.f90:27`), so a staged FAST record is written with no change;
  the netCDF `time` axis is `time_to_decimal_year` which **already encodes hour/min/sec**; `FC_DAY` is
  the finest file chunk (`<prefix>-F-YYYYMMDD.nc`, ~`n_fast_per_slow` records/day on the unlimited time
  dimension); `t_open` (period start) + `cell_methods='time: mean'` are already stamped.

### The two real gaps

1. **No fast-cadence tick.** `output_integrate` (the per-step fold) is `do t = 2_ik, N_FREQ`
   (`meds_output_integrate.f90:340`) — it hard-skips FAST — and it is called **once per slow step**
   from `meds_main.f90:178`, *after* `advance_one_step` has already collapsed the entire `dt_fast`
   loop into slow accumulators. Sub-daily resolution exists **only inside `run_fast_biophysics`'s
   `isub` loop**, so the FAST tier must be fed from there.
2. **Fluxes are the wrong quantity.** The existing `gpp_site`/`et_site` are `AGG_SUM` vars whose
   reductions read **slow-step accumulators** (`cohort%gpp_accum`, `site%et_accum`) that are reset once
   per slow step and grow monotonically across sub-steps — they cannot be sampled at fast cadence
   (they would never reset per fast window). FAST needs **instantaneous rate twins** (W/m², µmol/m²/s).
   State reservoirs (CAS temp, soil temp) are fine as-is: they are `AGG_TMEAN` and instantaneously
   valid at any sub-step.

Plus one shutdown bug to fix: `output_manager_close`'s partial-flush loop is `do t = 2_ik, N_FREQ`
(`meds_output_manager.f90:44`) — it would drop the final partial FAST window at run end.

## 3. Core strategy — buffer, do not restructure

The fast loop is **patch-outer / sub-step-inner** (`meds_fast_loop.f90:213` / `:267`). At any instant
inside the inner loop only the *current* patch has reached sub-step `isub`; earlier patches already
finished all sub-steps and later patches have not started. **No point in the loop nest holds a
site-wide state at a single sub-step**, so a site FAST record cannot be ticked inside the inner loop.

Two ways out were weighed:

- **Strategy A — assembly buffer (CHOSEN).** Keep the loop order. At each `(ip, isub)`, area-weight the
  live per-sub-step site quantities into a manager-owned, netCDF-free buffer indexed by sub-step. After
  the patch loop, that buffer holds the complete site-wide value at each of the `n_fast_per_slow`
  sub-steps. `meds_main` then replays the buffer through the unchanged reduction engine.
- **Strategy B — restructure to sub-step-outer / patch-inner (REJECTED).** Would require keeping every
  patch's `bio`/`coh`/`forc`/`aero` working bundle (incl. `leaf_temp`/`psi` + CAS/soil reservoirs)
  alive simultaneously — a full rewrite of gather/seed/write-back and ~`n_patch`× peak memory, for no
  benefit the buffer doesn't already give.

Strategy A reuses the proven area-weighted `et_accum` pattern (`meds_fast_loop.f90:287`), only resolved
onto a sub-step axis instead of a single slow-step total.

## 4. Architecture & data flow (P0)

```
run_fast_biophysics(site, ctx, cfg, met_drv, step_start, [mgr])       [meds_aux, netCDF-FREE]
  before patch loop:  if FAST active -> zero mgr%fast(1:n_fast_per_slow)
  do ip:  ... do isub:
      column_fast_step(... bio ..., le_flux)            (evolves bio in place)
      assemble a fast_sample_t s from LIVE state:
          s%cas_temp      = bio%cas%can_temp
          s%soil_temp_top = bio%soil_e%soil_temp(1)
          gpp_patch       = sum_j gpp_coh(j)*coh%nplant(j)   (bind to named var: nvfortran #7)
          s%gpp_rate      = gpp_patch
          s%le_flux       = le_flux
      area-weight into the sub-step slot (areas sum to 1 -> site mean):
          mgr%fast(isub) += site%patch%area(ip) * s ;  mgr%fast_time(isub) = t_sub
  after patch loop:   mgr%fast_ready = .true.

meds_main (per slow step, AFTER advance_one_step returns)              [drives netCDF]
  if FAST enabled .and. mgr%fast_ready:
     do isub = 1..n_fast_per_slow:
        call output_integrate_fast(mgr, mgr%fast(isub), mgr%fast_time(isub), dt_fast)   [core, netCDF-free]
        fast_step_total += 1
        if mod(fast_step_total, mgr%fast_interval_steps) == 0:
           call close_tier(mgr, 1)                 (normalize buf(:,1) -> pending(1))
           call output_serialize_pending(mgr)      (drain pending(1) -> <prefix>-F-YYYYMMDD.nc)   [serializer]
     mgr%fast_ready = .false.
  ... existing DAILY/MONTHLY/ANNUAL tick (meds_main:178) UNCHANGED ...
```

**DAG hygiene (hard constraint).** `run_fast_biophysics` lives in `meds_aux`, which links the
netCDF-free `meds_output_core` but **not** the serializer. It may write `mgr%fast(:)` (netCDF-free
buffers on `output_manager_t`) but **must not** call `output_serialize_pending` / `stream_write_record`
— those stay in `meds_main`. The `output_manager_t` type itself is in `meds_output_types.f90`
(netCDF-free core), so threading `mgr` in as an `optional intent(inout)` argument is DAG-legal.

**Why the `fast_sample_t` array (not a bare `[var,substep]` matrix):** the fast loop fills a fixed
physical sample (known fields) and knows nothing about the registry; the registry→sample mapping happens
in `output_integrate_fast` via `extract_fast_scalar(source_id, sample)`. This keeps the fast driver
decoupled from the variable table and mirrors how `extract_scalar_source` already works for the slow
tiers.

**Single-slot `pending(1)`:** interleaving `close_tier` + `output_serialize_pending` inside the replay
loop drains each FAST window before the next overwrites `pending(1)` — no pending ring needed.

## 5. Variable set

| var | dim | agg | group | source | how |
|---|---|---|---|---|---|
| `cas_temp_site` | scalar | TMEAN | energy | `SRC_S_CAS_TEMP` | **reuse** — OR `FREQ_FAST` into its mask |
| `soil_temp_top_site` | scalar | TMEAN | energy | `SRC_S_SOIL_TEMP_TOP` | **reuse** — OR `FREQ_FAST` into its mask |
| `gpp_rate_fast` | scalar | TMEAN | carbon | `SRC_F_GPP_RATE` (new) | **new FAST-only** var, µmol m⁻² s⁻¹ |
| `le_flux_fast` | scalar | TMEAN | water | `SRC_F_LE` (new) | **new FAST-only** var, W m⁻² |

The two existing state vars simply gain the `FREQ_FAST` bit (they are already `AGG_TMEAN` reservoirs,
valid at any sub-step). The two flux vars are **new FAST-only** (`streams = FREQ_FAST`) instantaneous
rate twins — deliberately *not* the `AGG_SUM` accumulator vars `gpp_site`/`et_site`, which stay
slow-tier period totals. Extend this set in P1 (§7).

## 6. Aggregation semantics & timestamping

- **Operator:** `AGG_TMEAN` for every FAST var — the dt-weighted mean over the `fast_interval_steps *
  dt_fast` window (with uniform `dt_fast` this equals the simple mean of the window's sub-steps). Use
  `AGG_LAST` only if a strict instantaneous snapshot is ever wanted; TMEAN is the natural diurnal read
  and keeps CF `cell_methods='time: mean'` consistent with the coarse tiers.
- **Window timestamp:** `t_open` = window **start** (already stamped by `close_tier`); the CF window is
  `[t_open, t_open + fast_interval_steps*dt_fast)`.
- **Sub-step timestamp:** `t_sub = step_start + (isub-0.5)*dt_fast` (interval **midpoint**, matching the
  probe) recorded per sub-step for the fold.
- **Cadence knob:** `fast_interval_steps` (already on `mgr`) sets records-per-day = `n_fast_per_slow /
  fast_interval_steps × steps_per_day`. **Recommend `fast_interval_steps = 1` (or 2)** in the run
  config for diurnal work — the finest, cleanest curve (48/day at `dt_fast=1800 s`). It **must divide**
  `n_fast_per_slow` for a stable local-hour axis (validate at config load; §9).

## 7. Phased implementation

### P0 — MVP: working site-scalar FAST tier (4 variables)

Deliverable: enabling `[output.fast]` writes `<prefix>-F-YYYYMMDD.nc` with `n_fast_per_slow /
fast_interval_steps` records/day of area-weighted site **CAS temperature, top-soil temperature, GPP
rate, latent-heat flux**, with DAILY/MONTHLY/ANNUAL output **bit-identical**.

1. **`meds_output_types.f90`** — add a netCDF-free `fast_sample_t` derived type (`cas_temp`,
   `soil_temp_top`, `gpp_rate`, `le_flux`, all `real(wp)`). Add to `output_manager_t` (after line 163):
   `type(fast_sample_t), allocatable :: fast(:)`, `type(meds_time_t), allocatable :: fast_time(:)`,
   `integer(ik) :: n_fast_sub = 0`, `logical :: fast_ready = .false.`. Provide an elemental/`+` helper
   or explicit field accumulation for the area-weighting.
2. **`meds_output_integrate.f90`** — add ids `SRC_F_GPP_RATE`, `SRC_F_LE` near the `SRC_S_*` block
   (~line 70). Add `pure function extract_fast_scalar(source_id, s)` mapping `SRC_S_CAS_TEMP →
   s%cas_temp`, `SRC_S_SOIL_TEMP_TOP → s%soil_temp_top`, `SRC_F_GPP_RATE → s%gpp_rate`, `SRC_F_LE →
   s%le_flux`. Add `subroutine output_integrate_fast(mgr, sample, now, dt)`: if `reg%nidx(1)==0` return;
   if `.not. has_data(1)` set `t_open(1)=now`; for each `j=1..nidx(1)`, `k=idx_freq(j,1)`, call
   `integrate_scalar(buf(k,1), extract_fast_scalar(reg%var(k)%source_id, sample), dt)`; set
   `has_data(1)=.true.`. All FAST vars are `DIM_SCALAR` → no slab path in P0. **Leave the `do
   t=2,N_FREQ` fold (line 340) and every coarse close trigger untouched.**
3. **`meds_output_registry.f90`** — define `FAST_ONLY = FREQ_FAST` (near the `DAY_MON_YR` mask,
   ~line 58). OR `FREQ_FAST` into the `streams` arg of the existing `cas_temp_site` and
   `soil_temp_top_site` rows. Add two `add_variable` rows: `gpp_rate_fast` (`SRC_F_GPP_RATE`,
   `DIM_SCALAR`, `AGG_TMEAN`, `GRP_CARBON`, `'umol m-2 s-1'`, `FAST_ONLY`) and `le_flux_fast`
   (`SRC_F_LE`, `DIM_SCALAR`, `AGG_TMEAN`, `GRP_WATER`, `'W m-2'`, `FAST_ONLY`). No
   `extract_scalar_source` cases needed — FAST-only vars are never sourced by the slow tick.
4. **`meds_fast_loop.f90`** — add `optional, intent(inout) :: mgr` (type `output_manager_t`) to
   `run_fast_biophysics`. Before the patch loop (~line 211), if `present(mgr) .and. mgr%enabled .and.
   mgr%reg%nidx(1)>0 .and. do_forcing`: lazily allocate `mgr%fast`/`mgr%fast_time` to
   `cfg%n_fast_per_slow`, set `n_fast_sub`, zero the samples. Inside the `isub` loop at the
   `write_fast_probe`/`et_accum` site (~line 287): assemble the sample from live state, `gpp_patch =
   sum(gpp_coh(1:ncoh)*coh%nplant(1:ncoh))` **bound to a named var** (nvfortran #7), and accumulate
   `mgr%fast(isub) += site%patch%area(ip) * sample`; set `mgr%fast_time(isub)=t_sub`. After the patch
   loop set `mgr%fast_ready=.true.`.
5. **`meds_stepper.f90`** — add `optional, intent(inout) :: mgr` to `advance_one_step` (line 33) and
   forward it to `run_fast_biophysics` at both the forcing-on and forcing-off call sites (~line 50/52).
6. **`meds_main.f90`** — declare `integer(ik) :: fast_step_total = 0`. Pass `mgr` into both
   `advance_one_step` calls (166/169) when `cfg%output%enabled .and. cfg%output%freq_on(1) .and.
   cfg%fast_biophysics_on`. After `advance_one_step` returns and before the existing tick (line 178),
   add the replay/drain block from §4 (loop `isub`; `output_integrate_fast`; on `fast_interval_steps`
   boundary `close_tier(mgr,1)` + `output_serialize_pending(mgr)`). Keep line 178/179 unchanged.
7. **`meds_output_manager.f90`** — change the partial-flush loop from `do t = 2_ik, N_FREQ` (line 44) to
   `do t = 1_ik, N_FREQ` so the final partial FAST window is closed + drained at run end (already
   guarded by `if (mgr%has_data(t))`). **Land this in P0** — the cross-slow-step partial window relies
   on it.
8. **Gate:** the whole FAST path is off unless `freq_on(1) .and. fast_biophysics_on .and. forcing_on`
   (the diurnal signal only exists when live forcing drives the sub-step loop). Config with
   `freq_on(1)=.true.` but forcing off must not error — gate cleanly, optionally warn.

### P1 — energy-balance fluxes, soil profile, sub-daily axis, sync policy

Deliverable: add the diurnal energy/water fluxes and the soil profile at FAST cadence, plus a
human-readable time axis and a FAST-appropriate sync policy.

- **`meds_column_dynamics.f90`** — surface new `optional intent(out)` args from `column_fast_step`:
  sensible heat `H` (`rho*ustar*temp1*(can_temp-theta_atm)`), net radiation `Rn` (`coh_rnet` + net
  ground), incident SW (`ctx_now%rad_sw_top`), `ustar` (`aero%ustar`). Only `le_flux` is returned today
  (~line 285).
- **`meds_fast_loop.f90`** — extend `fast_sample_t` + the fill with `h_flux`, `rnet`, `sw_in`, `ustar`
  (area-weighted). Add the site `DIM_SOIL` columns (`bio%soil_e%soil_temp(:)`, `bio%soil_w%theta(:)`) —
  state reservoirs, instantaneously valid, site-level per layer, so no cohort/patch guard trips.
- **`meds_output_integrate.f90`** — add `SRC_F_H`, `SRC_F_RNET`, `SRC_F_SW_IN`, `SRC_F_USTAR` cases to
  `extract_fast_scalar`; add a `DIM_SOIL` slab replay path in `output_integrate_fast` using
  `integrate_slab`/`normalize_slab` (fixed-slot within the ≤1-day window). Widen the manager sample to a
  slab variant sized from `cfg`.
- **`meds_output_registry.f90`** — register the new flux vars (`AGG_TMEAN`, `GRP_ENERGY` for H/Rn/SW/
  ustar, `GRP_WATER` for LE) and the `DIM_SOIL` FAST vars (`soil_temp_site`, `soil_water_site` gain
  `FREQ_FAST`), each with its `add_variable` row + `extract_fast_scalar` case (two-edit coupling).
- **`meds_output_stream.f90`** — optional human-readable sub-daily companion coords: define
  `v_hour`/`v_minute`/`v_second` (or a single `seconds_into_day`) in `stream_open_file` gated on
  `tier==1`, mirroring `v_year`/`v_month`/`v_day` (168–171); write them from `pr%t_open` in
  `write_one_record` (mirroring 268–270); add fields to `stream_file_t`. Purely additive — the decimal
  `time` axis already disambiguates records, so this is convenience for group-by-hour analysis.
- **`meds_output_config.f90` / `meds_config_io.f90`** — add an optional per-tier sync policy (or default
  the FAST tier to `SYNC_NEVER` with an end-of-file sync) so ~`n`/day records do not each `nc_sync`
  (`meds_output_stream.f90:55`).
- **Retire/flag the probe** — `write_fast_probe` is now superseded for production diurnal analysis. Keep
  it behind its own `[fast].fast_probe` flag as an independent integrator/numerics debug tool (cheap,
  self-contained); add a one-line comment noting the netCDF FAST tier supersedes it. Do not delete.

### P2 — per-cohort / per-patch sub-daily slabs (deferred; out of scope for diurnal deliverable)

Only if later required. A FAST cohort/patch var forces `FC_DAY` (one file per sim-day) so the live
count is stable within a file (no intra-day fuse/fission) — sidestepping the serializer's `FC_MONTH`
cap (`meds_output_stream.f90:46-47`) and count-grew `error stop` (`:262-265`) that a >1-day file would
trip. Fixed-slot slabs, **no** `global_id` keying (a sub-daily window is ≤1 day). `enforce_annual_guard`
already permits `FREQ_FAST` on every dim, and the per-variable `meds_io_config.toml` `'F'` override is
already parseable (`parse_stream_mask`). Multi-day cohort/patch FAST files (needing id-keyed slabs) stay
explicitly unimplemented.

## 8. Constraints, risks & gotchas

- **nvfortran array-temp trap (CLAUDE.md #7):** `gpp_patch = sum(gpp_coh*nplant)` and any soil-column
  fill must bind array results to a named var before passing into a call. A green ifx run is **not**
  sufficient — build nvfortran multicore on every touched module and diff FAST output ifx vs nvfortran.
- **`fast_interval_steps` must divide `n_fast_per_slow`** — else records/day and the local-hour axis
  drift. Validate at config load and document. Size `mgr%fast` from `cfg%n_fast_per_slow`, **never a
  hard-coded 48** (`dt_slow` need not be one day).
- **Partial window across slow steps:** if `fast_interval_steps` does not divide `n_fast_per_slow`, a
  partial FAST window carries across the slow-step boundary (correct — accumulation continues), and the
  final partial relies on the `output_manager_close` `t=1` fix (P0 step 7).
- **DAG hygiene:** a single stray `output_serialize_pending`/`stream_write_record` call from `meds_aux`
  breaks the acyclic library DAG. The fast loop writes `mgr%fast(:)` only.
- **Bit-identical coarse tiers:** any accidental edit to `output_integrate`'s `do t=2,N_FREQ` fold or
  the DAILY/MONTHLY/ANNUAL close triggers breaks the golden anchors. Keep the FAST driver strictly
  separate; add a byte-match regression guard (§9).
- **FAST = site-scalar + soil-column only:** the serializer's count-grew guard + `FC_MONTH` cap forbid
  cohort/patch slabs at sub-daily cadence. This matches the diurnal use case exactly — treat it as a
  hard registry-build constraint, not a limitation to lift.
- **UTC day-boundary:** `FC_DAY` file roll uses the model calendar `t_open`; a UTC run clock splits the
  diurnal file mid-local-night. This is an upstream calendar convention, not a serializer bug — document
  it (users analyzing local diurnal cycles should be aware).
- **Forcing-off:** FAST produces nothing without live forcing (no diurnal signal, no `t_sub`); the gate
  handles it.

## 9. Testing plan

- **Unit** (`test/`): `output_integrate_fast` folds a hand-built `fast_sample_t` sequence into `buf(:,1)`
  and `close_tier(mgr,1)` normalizes to the expected TMEAN; `extract_fast_scalar` maps every source id.
- **Integration:** enable `[output.fast]` on a 2-day forced Ithaca run; assert `-F-YYYYMMDD.nc` exists
  with `n_fast_per_slow/fast_interval_steps` records/day, physical diurnal shape (GPP 0 at night, peak
  midday; CAS temp/LE tracking SW), and area-weighted site means. Cross-check against the CSV probe on
  the same run (values must agree within rounding).
- **Regression (golden):** with `[output.fast].enabled=false`, DAILY/MONTHLY/ANNUAL netCDF output is
  byte-identical to pre-change (coarse tiers untouched).
- **Portability:** ifx Debug (`-check all`) **and** nvfortran multicore both build and pass the suite;
  FAST output diffed ifx vs nvfortran (guards the array-temp trap).

## 10. Out of scope

Per-cohort/per-patch sub-daily output (P2, deferred); FAST feeding DAILY/MONTHLY/ANNUAL (tier chaining
stays deferred — each coarse tier integrates raw slow state independently); variance/second-moment on
FAST; an `FC_HOUR` file chunk (`FC_DAY` + unlimited time axis is the intended granularity); a separate
sub-daily time coordinate (`time_to_decimal_year` already resolves it); removing `write_fast_probe`
(kept as a debug switch); any forcing/LWdown work (orthogonal). **Adversarial review** is intentionally
excluded here — do it as a later pass after this design is edited.
