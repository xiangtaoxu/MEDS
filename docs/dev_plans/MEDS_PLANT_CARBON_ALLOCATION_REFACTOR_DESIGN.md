# MEDS plant carbon-allocation refactor — DESIGN + IMPLEMENTED (2026-07-20)

Status: **IMPLEMENTED** (+ round-2 refinements). Supersedes `MEDS_PLANT_CARBON_DYNAMICS_DESIGN.md`.
Validation (§7): ifx Debug 32/32 + nvfortran multicore 32/32; a 40-yr bare-ground spin-up completes
clean (75 cohorts, LAI 12.5, AGB 8.87 kgC/m², area conserved, no NaNs). Not bit-identical, as designed.

**Round-2 refinements (2026-07-20):** (a) `growth_respiration` relocated from `meds_plant_respiration`
into `meds_plant_carbon_allocation`, taking the *realized* growth `npp_growth` (not the pre-allocation
balance) — the seam for a future mechanistic construction cost; the kernel's `(1+g)` charge calls it.
(b) `plant_carbon_allocation` is now **growth-only**: turnover/shed is applied UPSTREAM by the driver's
`update_biomass_turnover` (shed-first), so growth demands see the post-shed pool and evergreen replacement
is exact within the step; the kernel no longer takes shed args, and `flush_growth_cap` was inlined.
(c) **Storage funds leaf/fine-root growth even when net < 0** (spring leaf-out from reserves): the
maintenance debt is paid from storage first, then P1 growth draws remaining reserves — fixes a leaf-out
deadlock latent under the fast loop (leafless ⇒ GPP≈0 ⇒ net<0 ⇒ no growth). (d) Closure identity sign
corrected: `Σ(growth pools + npp_store) − deficit = (gpp − resp_maint) − growth_resp`.
Companion science docs: **new** `docs/science/plant_carbon_allocation.md` (§8); **edit**
`docs/science/plant_phenology.md` to add the baseline-turnover floor (§4c/§8).

## 1. Motivation

`src/plant/meds_plant_carbon_dynamics.f90` today mixes three concerns: (a) the pure PARTEH-H1
distribution of the daily carbon budget, (b) tissue-turnover *rate* computation (with evergreen
cold-suppression), and (c) the active-shed rate→amount conversion with snap-to-bare. After the
phenology rate-refactor, (b) and (c) are really *phenology* — they decide how much leaf/root display
a plant carries — while only (a) is genuinely "carbon allocation." This refactor separates the three,
renames the module for what it does, and folds the whole daily NPP calculation (GPP − respiration −
allocation) into one **elemental** kernel.

Accepted decisions (review 2026-07-20, incl. answers to the v1 open questions):

1. **Rename** file/module to `meds_plant_carbon_allocation`; `plant_carbon_allocation` is the master
   entry, at `dt_slow`. It now **owns the full NPP calculation** — it takes raw cohort scalars
   (GPP, maintenance resp, growth-resp factor, current storage, shed amounts, capped demands) and
   returns per-pool NPP + growth respiration. It is **`elemental pure`** over the cohort SoA. The
   `carbon_gain_t/loss_t/demand_t/npp_t/env_t` structs are **retired**.
2. **Move tissue turnover into `meds_phenology`.** Baseline turnover is a degenerate phenology: a
   tissue with no active shed scheme sheds at its baseline rate. `leaf_shed_rate` carries a
   **baseline-turnover floor** via `max(active, baseline)`; the separate `baseline_loss` channel and
   the replaceable/non-replaceable split are **removed** (replaceability is emergent, §3). Root and
   (future) wood get the same treatment; disturbance/damage losses are a later module.
3. **Growth respiration is charged on realized growth**, inside the kernel — on carbon allocated to
   growth pools (leaf/fine-root/wood[/repro]), **excluding storage**. Resolves the current
   inconsistency where growth resp was charged on the whole pre-allocation balance (§3).
4. **Litter stays in the driver** (`meds_vegetation_dynamics`): the kernel returns `npp` + `growth_resp`
   only; the driver derives `leaf_litter = leaf_shed + fineroot_shed` for the biogeochem seam and routes
   `growth_resp` to the autotrophic-respiration diagnostics.
5. **No hard-coded tunables:** evergreen constants, the bare-snap fraction, and the phenology cue
   widths (`gdd/daylen/soiltemp`) all move to **per-PFT** parameters.
6. **`leaf_shed_amount` lives in the carbon-allocation file** (a public `elemental` helper) even though
   it is, strictly, a phenology→carbon *bridge* (rate → mass), not allocation proper — kept here for
   cohesion; noted as such in code + science doc.

## 2. Target architecture (three layers)

```
meds_phenology            (plant lib)  SIGNAL + RATE authority
   emits per-tissue rates [1/day]:  leaf_flush_rate, leaf_shed_rate, fineroot_shed_rate
   leaf_shed_rate     = max( k_shed_max·shed_drive ,  baseline_leaf_turnover·cold_factor )
   fineroot_shed_rate = baseline_fineroot_turnover·cold_factor         (no active root phenology yet)
   owns: turnover baseline, evergreen cold-suppression, all cue widths (now PFT params)

meds_plant_carbon_allocation  (plant lib)  KERNELS (all elemental pure)
   plant_carbon_allocation(gpp, resp_maint, growth_resp_frac, storage,          & ! MASTER
                           leaf_shed, fineroot_shed,                            &
                           leaf_demand, fineroot_demand, storage_demand, repro_frac, &
                           npp_leaf, npp_fineroot, npp_wood, npp_store, npp_repro,    &
                           growth_resp, deficit, starving)
   leaf_shed_amount(shed_rate, flush_rate, leaf_carbon, leaf_carbon_full, bare_snap_frac, dt_day)
   flush_growth_cap(flush_rate, leaf_carbon_full, dt_day)              ! max leaf growth this step
   fill_carbon_demand(...)                                            ! private; charges (1+g) for growth

src/driver/ orchestrator   (meds_aux) ACCUMULATION + CALL SEQUENCING
   per cohort: gather GPP/resp (fast-loop accumulators or stub); phenology rates→amounts via the
   plant helpers (leaf_shed_amount, flush_growth_cap, ×root_to_leaf); allometric demands; ONE
   elemental call into plant_carbon_allocation over the SoA; derive litter; route growth_resp.
```

The plant seam `get_plant_flux_slow` in `meds_plant_interface` is **deleted** (its body split between
the driver orchestrator and the elemental kernel). `update_phenology` and the extended
`pheno_drives_to_rates` remain the phenology seams.

## 3. Science of the kernel

### 3a. Emergent replaceability (the load-bearing idea)

Today leaf loss is two channels: *baseline* (replaceable — P1 refunds it) and *active shed*
(non-replaceable). The refactor collapses them into one `leaf_shed_rate` and lets replaceability
emerge from whether **flush** is simultaneously filling the deficit:

- **Evergreen / growing season:** `shed = baseline`, `flush_rate > 0`. Shed opens a leaf deficit; the
  flush-capped growth step refills it (NPP, then storage). Canopy held ~full → baseline turnover is
  *replaced*, with no special branch.
- **Deciduous dormancy:** `shed = active` (large), `flush_rate → 0` → flush cap → 0 → leaf demand
  capped to 0 → no refill → canopy drives to bare → shed is *non-replaceable*.

So the ladder loses the old P1 baseline-replacement step and `can_flush` disappears (the flush cap on
`leaf_demand/fineroot_demand` already encodes the gate; storage is simply allowed at the growth step).

**Snap-to-bare re-gates:** today it fires on `shed_rate > 0`, which a baseline floor makes always true
(would snap a slowly-turning-over evergreen). New rule (in `leaf_shed_amount`): snap only when the
canopy is **dormant** (`flush_rate ≤ ε`) *and* remaining leaf would fall below
`bare_snap_frac · leaf_carbon_full`.

### 3b. Growth respiration on realized growth

`net = gpp − resp_maint`. When `net ≥ 0` the kernel distributes it down the priority ladder; building
one unit of growth tissue consumes `(1+g)` carbon (`g = growth_resp_frac`), storage refill consumes
`1` (no construction cost). This charges growth resp on exactly the realized growth and resolves the
circularity non-iteratively.

**Priority ladder** (`net ≥ 0`):
1. grow leaf + fine-root toward the **flush-capped** deficit — storage allowed, cost `(1+g)`
2. refill storage toward target — NPP only, cost `1`
3. reproduction — fraction of residual, NPP only, cost `(1+g)` *(open Q, §9.1)*
4. wood — residual sink, cost `(1+g)`; leftover → storage

`growth_resp = g·(leaf + fineroot + wood [+ repro] built)`. When `net < 0`: pay the maintenance debt
from storage; no growth/growth-resp; set `deficit`/`starving` if storage is short (unchanged).

**Exact carbon closure** (both branches, with `deficit` on the negative branch):

```
Σnpp = gpp − resp_maint − growth_resp − leaf_shed − fineroot_shed
```

i.e. `plant-pool change = GPP − (maintenance + growth resp) − litter`. This is now self-consistent and
subsumes both the old header-comment bug (missing shed sink) and the growth-resp inconsistency.

**Not bit-identical.** Removing growth resp from the storage-bound carbon frees marginally more carbon,
and leaf replacement moves from an unconditional P1 to the flush-capped growth step. Spin-up equilibrium
pools shift slightly and must be validated, not asserted equal (§7).

## 4. Interface changes

### 4a. Types (`meds_plant_types`)

- **Retire** `carbon_gain_t`, `carbon_loss_t`, `carbon_demand_t`, `carbon_npp_t`, `carbon_env_t` — the
  elemental kernel is all scalars, and the driver writes results straight into the `carbon_flux_block`
  SoA / `npp_repro` array.
- **Move** `turnover_env_t/params_t/rates_t` into the phenology section, folded into the pheno types:
  - `pheno_params_t`: add `leaf_turnover_rate`, `fineroot_turnover_rate` [1/yr], `evergreen`,
    `evg_ref_temp`, `evg_slope`, `gdd_width`, `daylen_width`, `soiltemp_width`.
  - `pheno_env_t`: add `tissue_temp` (evergreen cold-suppression driver).
  - `pheno_out_t`: add `fineroot_shed_rate`.

### 4b. `meds_plant_carbon_allocation` (renamed module)

Public surface: `plant_carbon_allocation` (master), `leaf_shed_amount`, `flush_growth_cap` — all
`elemental pure`; private `fill_carbon_demand` (renamed from `fund`, now `(1+g)`-aware, `elemental
pure`). No module-constant parameters; no turnover; no `can_flush`; no `WOOD_DEMAND_BIG`.

> Elemental note: every helper the kernel calls must itself be `elemental pure` (Fortran forbids an
> elemental procedure calling a non-elemental one). All are scalar arithmetic — no `exp`, no arrays —
> so this holds. The wide (~19) argument list is accepted; the single driver call site uses **keyword
> arguments** to keep it non-positional. Build the nvfortran multicore back end (issue #7).

### 4c. `meds_phenology`

`pheno_drives_to_rates` gains the turnover baseline + tissue temperature and emits three rates:

```fortran
! cold  = evergreen ? 1/(1+exp(evg_slope·(evg_ref_temp − tissue_temp))) : 1
! leaf_flush_rate     = k_flush_max · flush_drive
! leaf_shed_rate      = max( k_shed_max · shed_drive ,  (leaf_turnover_rate/yr_day) · cold )
! fineroot_shed_rate  = (fineroot_turnover_rate/yr_day) · cold
```

`phenology_kernel` reads `gdd_width/daylen_width/soiltemp_width` from `params` (were module constants).

### 4d. Driver orchestrator (`src/driver/meds_vegetation_dynamics`)

`carbon_growth` rewritten: per cohort gather `gpp`/`resp_maint` (accumulators or stub), compute
allometric deficits, cap leaf/root by `flush_growth_cap(...)` (`×root_to_leaf` for root), convert shed
rates via `leaf_shed_amount(...)`, then **one elemental** `plant_carbon_allocation` call writing the
`carbon_flux_block` arrays; set `leaf_litter(j) = leaf_shed + fineroot_shed`; accumulate `growth_resp`
into the resp diagnostic. Phenology-OFF path: `flush_rate = PHENOLOGY_OFF_FLUSH`,
`shed_rate = baseline turnover` (OFF still applies leaf lifespan, now via the shed channel — a small,
intended change). `flatten_pheno_params` extended with the new fields.

## 5. Parameter migration (decision 5)

New per-PFT fields + TOML keys + CSV columns + `build_test_config` defaults + `flatten_pheno_params`
wiring (all REQUIRED by the presence-map — land these first, in every shipped TOML and the test config):

| current constant | file | new PFT field | default |
|---|---|---|---|
| `evg_ref_temp` 278.15 | carbon | `pheno_evg_ref_temp` | 278.15 |
| `evg_slope` 0.4 | carbon | `pheno_evg_slope` | 0.4 |
| `ELONGF_MIN` 0.02 | carbon | `bare_snap_frac` (renamed) | 0.02 |
| `gdd_width` 50 | phenology | `pheno_gdd_width` | 50 |
| `daylen_width` 1.0 | phenology | `pheno_daylen_width` | 1.0 |
| `soiltemp_width` 2.0 | phenology | `pheno_soiltemp_width` | 2.0 |

`leaf_turnover_rate`, `fineroot_turnover_rate`, `evergreen` already exist as PFT fields — **re-routed**
from the interface's turnover flattening into `flatten_pheno_params` (no new keys).

## 6. File-by-file change list

- `src/plant/meds_plant_carbon_dynamics.f90` → **rename** `meds_plant_carbon_allocation.f90`; elemental
  kernel + `leaf_shed_amount` + `flush_growth_cap` + private `fill_carbon_demand`. (CMake globs
  `src/plant/*.f90` — rename is transparent.)
- `src/plant/meds_phenology.f90` — baseline floor + third rate + widths from params.
- `src/plant/meds_plant_types.f90` — retire the five carbon structs; extend the three pheno types (§4a).
- `src/plant/meds_plant_interface.f90` — delete `get_plant_flux_slow` + carbon/turnover re-exports;
  keep/adjust the phenology seams.
- `src/plant/meds_plant_capi.f90` + `examples/example_phenology/` — add `fineroot_shed_rate` (and any
  new params the example exercises) to the phenology C struct.
- `src/driver/meds_vegetation_dynamics.f90` — rewrite `carbon_growth`; extend `flatten_pheno_params`;
  derive litter; route `growth_resp`.
- `src/shared/config/meds_pft_params.f90` — new fields + allocations + defaults (§5).
- `src/io/meds_config_io.f90` — `load_phenology_pft` new `req_pa` keys; CSV header/row columns.
- `test/test_plant_carbon_dynamics.f90` → **rename** `test_plant_carbon_allocation.f90`; rebuild around
  the elemental kernel — closure identity (§3b), emergent replaceability, snap-gating, growth-resp-on-
  realized-growth, negative-NPP debt.
- `test/test_carbon_growth.f90`, `test/test_plant_phenology.f90` — new rate set + widths + floor.
- `test/meds_test_support.f90` (`build_test_config`) — supply all new params.
- Example / `meds_config` PFT TOMLs — add the six new keys.
- `docs/science/plant_carbon_allocation.md` — **new** (§8).
- `docs/science/plant_phenology.md` — **edit**: document `leaf_shed_rate = max(active, baseline)` floor
  + `fineroot_shed_rate`.
- `docs/dev_plans/MEDS_PLANT_CARBON_DYNAMICS_DESIGN.md` — mark superseded.

## 7. Risks & validation — OUTCOME

- **Not bit-identical** (§3), as designed. A 40-yr bare-ground spin-up (ifx Debug) reaches a sane canopy
  (93 cohorts, LAI 12.5, AGB 8.85 kgC/m², N=0.94, Dmean 7.7 cm; area conserved; no NaNs). The 250-yr
  shipped spin-up `error stop`s at the FINAL state write on `io_cohort_max` (>2048 cohorts): the
  uncalibrated **stub GPP** (`gpp_ref·leaf_area`, no light-limitation) accumulates biomass unboundedly
  over centuries — orthogonal to this refactor (repro now pays growth resp, which *lowers* recruitment).
  Bump `io_cohort_max` or add light-limited stub GPP separately if a full 250-yr state write is needed.
- **Test coverage:** `test_plant_carbon_allocation` (closure identity, growth-resp-on-realized-growth,
  storage exemption, wood residual, flush cap, shed decay + dormancy snap, no-snap-while-flushing,
  negative-NPP starvation, turnover floor `max(active, baseline)`); `test_carbon_growth` (closure with
  shed, negative branch, carbon-mode wood→dbh flip); `test_plant_phenology` unchanged and green.
- **Snap-to-bare regression** (subtlest point): test that a small evergreen at low leaf with
  `flush_rate>0` does **not** snap, while a dormant deciduous (`flush_rate≈0`) does.
- **Growth-resp accounting:** confirm total autotrophic respiration (maintenance + the new growth_resp)
  still closes the site carbon budget and that the CO2/NEE diagnostics receive `growth_resp`.
- **Presence-map breakage:** every new REQUIRED key must land in all TOMLs and `build_test_config` in
  the same commit or the loader `error stop`s. Land params first.
- **nvfortran / elemental:** all kernel-called helpers elemental+scalar; no array-valued function result
  fed into a call (issue #7). Build the nvfortran multicore back end, not just ifx.
- **Fine-root coupling:** root growth still tracks leaf via `root_to_leaf_ratio` under the flush cap;
  the root shed floor (= leaf turnover for now) must not strip roots when leaves are bare.

## 8. Science docs

### 8a. New — `docs/science/plant_carbon_allocation.md`

Mirror `plant_phenology.md` style (GitHub-protected math: ```math``` fences for display, `` $`…`$ ``
for conflicting inline). Sections:
1. **The daily carbon budget** — GPP → maintenance resp → growth resp → net; units (kgC plant⁻¹ / `dt_slow`).
2. **PARTEH-H1 allocation** — allometric targets, deficits, the priority ladder (§3), wood as residual.
3. **Growth respiration on realized growth** — the `(1+g)` construction charge, storage exemption, the
   non-iterative circularity resolution.
4. **Leaf display as emergent replaceability** — one flush rate + one shed rate; evergreen maintenance
   vs deciduous drop; flush cap; snap-to-bare on dormancy. The `leaf_shed_amount` rate→mass bridge.
5. **Carbon closure & litter** — the exact `Σnpp = gpp − resp_maint − growth_resp − shed` identity;
   litter to biogeochemistry; the negative-NPP debt/starvation branch.
6. **Parameters** — PFT traits consumed, cross-referencing the phenology doc.

### 8b. Edit — `docs/science/plant_phenology.md`

Add the **baseline-turnover floor**: `leaf_shed_rate = max(k_shed_max·shed_drive, baseline·cold)` —
turnover as degenerate phenology — and the new `fineroot_shed_rate`; note the evergreen cold-suppression
now lives here.

## 9. Resolved decisions (were open questions)

1. **Reproduction incurs growth respiration** — repro carbon is new biomass, charged `(1+g)` like
   leaf/fine-root/wood.
2. **Elemental signature width accepted** (~19 scalar args) with a **keyword call site**; to be
   revisited if the driver-folder review finds it unwieldy.
3. **Snap-to-bare dormancy threshold `ε` is a shared module constant** in the carbon-allocation file
   (e.g. `1e-6 day⁻¹`), not per-PFT.

Not-bit-identical is accepted; validation is "document the equilibrium delta" (§7).
