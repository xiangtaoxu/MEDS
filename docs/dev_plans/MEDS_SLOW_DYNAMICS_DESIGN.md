# MEDS slow-loop dynamics — design plan

**Status:** design-only (no code changed). **Scope:** the slow (daily → annual) dynamics tier. Two parts:

- **Part I — vegetation-driver cleanup** (`src/driver/meds_vegetation_dynamics.f90`): unconditional
  phenology + a vanilla-evergreen equivalence, four routine renames, and the compute/apply split (§1–§6).
- **Part II — biogeochemistry dynamics**: wire the CENTURY soil-carbon matrix into the slow loop as a new
  `meds_biogeochem_dynamics` domain driver (§7–§10).

**Out of scope** (explicitly): any `*_time_derivs` rename (that suffix is reserved for the fast tier's pure
re-entrant ARK RHS — the slow loop is orchestration + discrete events + advance-and-commit; see
`MEDS_DRIVER_REORG_DESIGN.md` §4.1). *(A thin `meds_slow_dynamics` coordinator, earlier deferred as
premature, is now **in scope for Part II** — see §10a — because Part II adds the second slow domain that
motivates it.)*

Part I is a **behavior-preserving-where-possible cleanup**; Part II is a **new feature**, not a reorg.
Part I and Part II are independent — Part I should land first, but neither blocks the other. Everything is
grounded in a direct read of the code; load-bearing facts are cited by `file:line`.

---

# Part I — Vegetation-driver cleanup

Three author requests:
1. Make **phenology unconditional** and express "no phenology" as a **vanilla-evergreen** parameterization.
2. **Rename** four routines to drop the legacy `carbon_` prefix.
3. **Revisit** the `compute_slow_derivatives` (driver) / `update_cohort_states` (core) split.

---

## 1. Decision summary

| # | Request | Decision |
|---|---|---|
| A | Always run phenology; "no phenology" = a vanilla evergreen | ✅ **Adopt the realistic 15-day flush as the new default** (option A2). A vanilla evergreen already reduces to the old phenology-off behavior *except* for the flush rate; we take the rate-limited flush deliberately — a documented science change, not a bit-identical no-op |
| B | Run model / validator | ✅ **Fast biophysics is always on**; phenology always runs (its daily-temperature source is guaranteed). Remove the `phenology_on` gate + the validator precondition. The **slow tier can be frozen** (static vegetation + soil); a **slow-only / empirical-demography** run is the Python C-API path |
| C | `advance_trait_dynamics` rename | ✅ **`advance_plant_traits`** (keeps the `advance_*` mutator verb + the `advance_leaf_phenology` twin). Public with 2 external callers — update them |
| D | `carbon_growth` rename | ✅ **`compute_carbon_allocation`** (`compute_*` producer verb) |
| E | `carbon_rates` rename | ✅ **`compute_vital_rates`** (`compute_*` producer verb) |
| F | `compute_slow_derivatives` rename | ✅ **`update_cohort_derivatives`** — verb-correct: it *does* update the cohort `cohort_deriv_block` (+ the growth ring buffer), `intent(inout)` |
| G | `compute_slow_derivatives` / `update_cohort_states` split | ✅ **Keep as-is** — load-bearing (mode-agnostic applier + offload split); only its documented rationale needs correcting (§5) |

The revised renames follow the codebase's own verb convention: **`compute_*` = producers** (`intent(in)`,
out-only: `compute_carbon_allocation`, `compute_vital_rates`), **`update_*` / `advance_*` = mutators**
(`intent(inout)`: `update_cohort_derivatives` fills the deriv block; `advance_plant_traits` advances the
trait state). Two nuances remain (§4): `advance_plant_traits` is public with external callers, and the
empirical mode-pair twin `compute_empirical_derivatives`.

---

## 2. Request A — unconditional phenology, "no phenology" = a vanilla evergreen

### 2.1 Why a vanilla evergreen is (almost) the old phenology-off path

Today the ON/OFF split lives entirely in `cohort_carbon_demand`
([meds_vegetation_dynamics.f90:344-354](../../src/driver/meds_vegetation_dynamics.f90#L344-L354)):

```
ON : call pheno_drives_to_rates(...)         ! flush_rate = k_flush_max * flush_drive
OFF: call turnover_shed_rates(...)           ! the baseline turnover floor
     flush_rate = PHENOLOGY_OFF_FLUSH (1e6)  ! instant fill
```

Decompose the two branches (`pheno_drives_to_rates`,
[meds_phenology.f90:214-229](../../src/plant/meds_phenology.f90#L214-L229); `turnover_shed_rates`,
[:191-203](../../src/plant/meds_phenology.f90#L191-L203)):

- **fine-root shed:** ON `= root_base`; OFF `= root_base`. **Identical, unconditionally.**
- **leaf shed:** ON `= max(k_shed_max·shed_drive, leaf_base)`; OFF `= leaf_base`. **Equal iff
  `shed_drive = 0`** (no active shed) — i.e. an evergreen.
- **`cold` suppression** (`1/(1+exp(evg_slope·(evg_ref_temp − tissue_temp)))`): both branches pass the
  same `STUB_TISSUE_TEMP = 298.15 K` ([:41](../../src/driver/meds_vegetation_dynamics.f90#L41)) →
  **identical**, and *temperature-independent*.
- **flush:** ON `= k_flush_max·flush_drive`; OFF `= 1e6` (instant). **This is the only real difference.**

The prognostic drives already **initialize to the evergreen fixed point** —
`PHENO_FLUSH_INIT = 1.0`, `PHENO_SHED_INIT = 0.0`
([meds_core_state_types.f90:41-42](../../src/core/meds_core_state_types.f90#L41-L42), commented "born
flushing (evergreen fixed point)") — and the default PFT cue masks are both `0` ("permissive flush / no
active shed = evergreen", [meds_pft_params.f90:179-193](../../src/shared/config/meds_pft_params.f90#L179-L193)).
So a vanilla evergreen already sheds at exactly the old turnover floor; only its **leaf flush** differs —
the old OFF path refilled the leaf deficit *instantly* (`flush = 1e6`), the ON evergreen refills over
~15 days (`pheno_k_flush_max = 0.06667`, [:195](../../src/shared/config/meds_pft_params.f90#L195)).

### 2.2 Decision — adopt the realistic 15-day flush (A2)

Make phenology unconditional and **keep the default `pheno_k_flush_max = 0.06667`**: the vanilla evergreen
fills its canopy over ~15 days instead of instantly. This is a **deliberate, documented science change**
from the old instant-fill, not a bit-identical migration — and it is the natural default now that
phenology always runs (a canopy does not refill in a single slow step).

**Consequence:** the slow-loop regression baseline **changes** at this step (the old no-phenology numbers
no longer match). Re-baseline the slow-loop smoke, document the one-time shift, and gate on the fast-on
carbon smoke (both compilers). There is **no** bit-identical requirement here.

(If the exact old instant-fill behavior is ever wanted, it is recoverable per-PFT with
`pheno_k_flush_max = 1e6` — a "phenology-off-equivalent" evergreen preset — but that is not the default.)

---

## 3. Request B — the run model (fast always on) + removing the phenology gate

**New run-model policy (this plan adopts it): the standalone Fortran MEDS always runs the fast
biophysics loop.** The fast loop is the always-on floor; the slow loop is the optional layer on top — the
*inverse* of today's opt-in-fast default. Because the fast loop supplies the sub-daily GPP/energy/water
**and** the daily-mean temperature phenology needs:

- **Phenology always runs.** Delete the `cfg%phenology_on` gate
  ([meds_vegetation_dynamics.f90:73-77](../../src/driver/meds_vegetation_dynamics.f90#L73-L77)) — and the
  whole per-PFT validator question dissolves. With fast always on, `forcing_on` (hence the daily
  temperature + `step_start`) is always present, so a deciduous PFT is never at risk of the
  frozen-init-drives trap. The old validator (`phenology_on ⟹ fast+forcing`,
  [meds_config.f90:349-353](../../src/shared/config/meds_config.f90#L349-L353)) is **removed, not relaxed**.
- **`doy` is always available.** Thread `step_start` (→ `doy`) **unconditionally** from `meds_main`
  (`prev`, the sim clock, is always known — it is just not passed on the fast-off branch today,
  [meds_main.f90:191-194](../../src/driver/meds_main.f90#L191-L194)). The two `doy`-absent error stops
  ([meds_stepper.f90:62-63](../../src/driver/meds_stepper.f90#L62-L63);
  [meds_vegetation_dynamics.f90:74-75](../../src/driver/meds_vegetation_dynamics.f90#L74-L75)) become dead
  and are removed.
- **Turning the SLOW tier off** (frozen vegetation + soil) is the supported way to run "biophysics only":
  a config toggle skips `vegetation_dynamics` (and the future biogeochem step, Part II), holding the
  demographic + soil-carbon state static while the fast loop runs.
- **Slow-only / empirical-demography runs move to Python.** A pure-demography run (no fast biophysics,
  demographic rates supplied externally) is the C-API path: `Site.apply_rates` (empirical — Python
  supplies growth/mortality/recruitment) or `Site.advance_slow` (the carbon engine), exercised by
  `examples/example_demography/{empirical_spinup,run_carbon}.py`. This is why the empirical
  `compute_empirical_derivatives` path **stays** (§4) — it *is* the slow-only-via-Python mechanism.

This policy is recorded in `CLAUDE.md` and the project memory (updated alongside this plan).

---

## 4. Requests C–F — the renames

Final scheme (the author's revised choices — all verb-correct against the routine's `intent`):

| Routine | `site` intent | Role | **New name** | Verb fit |
|---|---|---|---|---|
| `carbon_growth` ([:249](../../src/driver/meds_vegetation_dynamics.f90#L249)) | `intent(in)` | producer (out `npp`, `npp_repro`) | **`compute_carbon_allocation`** | ✅ `compute_` = producer |
| `carbon_rates` ([:193](../../src/driver/meds_vegetation_dynamics.f90#L193)) | `intent(in)` | producer (out `mortality`, `recruitment`) | **`compute_vital_rates`** | ✅ `compute_` = producer |
| `compute_slow_derivatives` ([:154](../../src/driver/meds_vegetation_dynamics.f90#L154)) | `intent(inout)` | fills `site%deriv` + advances the growth ring buffer | **`update_cohort_derivatives`** | ✅ `update_` = mutator — it *does* update the `cohort_deriv_block` |
| `advance_trait_dynamics` ([:465](../../src/driver/meds_vegetation_dynamics.f90#L465)) | `intent(inout)` | mutator (acclimates sla/vcmax25/rd25/llspan) | **`advance_plant_traits`** | ✅ `advance_` = step-a-prognostic-state; keeps the `advance_leaf_phenology` twin |

`update_cohort_derivatives` → `update_cohort_states` now reads as a clean **compute-then-apply pair** (the
driver updates the derivatives; the core applies them to the states).

**Two things to handle:**

1. **`advance_plant_traits` is public** ([public list :37](../../src/driver/meds_vegetation_dynamics.f90#L37))
   with two external callers (`meds_main.f90`, `test_carbon_growth.f90`) — update the public list + both
   call sites. (The other three are private single-file swaps; no collisions — none of the four names
   exist in `src/`/`test/`.) Minor echo to note: the trait-plasticity *kernel* module is
   `meds_plant_trait_dynamics`; `advance_plant_traits` is the driver routine that calls it — different
   namespaces, acceptable.
2. **The empirical twin — dissolve it (capi-only, author).** The carbon `update_cohort_derivatives` has a
   structural twin in the capi, `compute_empirical_derivatives`
   ([meds_demography_capi.f90:169](../../src/capi/meds_demography_capi.f90#L169)), whose *only* caller is
   `meds_apply_rates` ([:136](../../src/capi/meds_demography_capi.f90#L136)). The empirical path stays live
   (it is the Request-B Python slow-only mechanism, `Site.apply_rates`, exercised by
   `examples/example_demography/empirical_spinup.py`), but the author's decision is to **dissolve
   `compute_empirical_derivatives` and inline its body directly into `meds_apply_rates`** — the empirical
   rate-application is thin capi glue (a forward-allometry loop over `fill_cohort_deriv`), not a first-class
   peer deserving its own named routine, and there should be **no empirical twin** shadowing the carbon
   side. So `update_cohort_derivatives` stands alone. This change is **confined to
   `meds_demography_capi.f90`** — it does *not* touch `meds_vegetation_dynamics.f90` or the carbon split
   (§5); `fill_cohort_deriv` (core) remains the shared builder both paths call.

---

## 5. Request G — the `compute_slow_derivatives` / `update_cohort_states` split

*(This section uses `compute_slow_derivatives` by its current name; §4 renames it to
`update_cohort_derivatives`.)*

**Keep the split. The organization is near-optimal; only its rationale needs correcting.** (The author's
"dissolve … and inline to `meds_apply_rates`" comment targets the *empirical* twin
`compute_empirical_derivatives` in the capi — **not** this carbon routine, and not
`meds_vegetation_dynamics.f90` at all; see §4.2. The carbon compute→apply split, and its
`update_cohort_derivatives` rename, stand.)

The author's recollection — the split exists "to keep a stateless plant module" — **misattributes a real
invariant to the wrong split.** Neither computer touches a plant kernel: `compute_slow_derivatives` calls
only `carbon_to_structure` (shared allometry) + `fill_cohort_deriv` (core); the empirical twin calls only
inlined allometry + `fill_cohort_deriv`. Plant statelessness is upheld by a *different* split —
`carbon_growth`/`carbon_rates` do the site-level sweep and call the per-individual kernels, so the plant
library stays state-free.

The split's **true** two-fold reason is stated verbatim in the core headers:
1. **Mode-agnostic applier.** `update_cohort_states` "is MODE-AGNOSTIC: whether the tendencies came from
   the carbon flip or the empirical growth law is the driver's concern, invisible here"
   ([meds_core_state_update.f90:6-8](../../src/core/meds_core_state_update.f90#L6-L8)); `fill_cohort_deriv`
   is "the shared per-cohort TENDENCY BUILDER used by BOTH the carbon (driver) and empirical (capi)
   computers" ([:12-15](../../src/core/meds_core_state_update.f90#L12-L15)). Two mode-specific producers
   feed one law-free `cohort_deriv_block`; the core applies it. This *is* the mechanism/policy wall.
2. **Offload discipline.** The applier is a bare-array, arithmetic-only OpenMP-`target` kernel
   ([:100-141](../../src/core/meds_core_state_update.f90#L100-L141)); the tendency *computation* is
   host-only because it is branchy (the `wood_carbon→dbh` flip). The split cleanly separates the
   offloadable applier from the un-offloadable computer.

**The ring-buffer side effect is real but well-placed.** `compute_slow_derivatives` is `intent(inout)`
because, via `fill_cohort_deriv`, it advances the prognostic growth moving-average ring buffer
([:63-71](../../src/core/meds_core_state_update.f90#L63-L71)). Ordering is load-bearing: `carbon_rates`
reads `growth_avg` for Camac mortality **before** `compute_slow_derivatives` refreshes it
([call-site comment :104](../../src/driver/meds_vegetation_dynamics.f90#L104)). This side effect is
honestly declared and delegated to a single shared core routine — the right granularity.

**Rejected alternatives:** *relocate into core* (would put carbon policy into the mode-agnostic engine and
drag empirical policy in too — breaks CORE-applies/DRIVER-computes); *merge the two computers* (needs a
carbon-vs-empirical branch inside one routine, undoing the deliberate edge-split; they live in different
libraries with different build gates); *de-side-effect the ring buffer* (net-negative: a second cohort
loop + second core call in both modes, separating the realized-rate sample from the new-state it derives
from, for only nominal "purity" the current honest `intent(inout)` already communicates).

**Actionable item = documentary only:** correct any comment that implies this split is about plant
statelessness (it is a mode-agnostic-applier + offload split), and let this finding kill the `get_` rename
in §4 (the routine is not a pure getter).

---

## 6. Phased implementation sequence (Part I)

Standing discipline (every phase): build **both** ifx Debug and nvfortran multicore; run the full `ctest`
suite; and a **fast-on carbon smoke** run (now the primary slow-loop regression, since fast is the
run-model floor). The fast-loop golden anchor (`tc_split(54)=292.450065`) is untouched — none of this
touches the fast loop.

**Phase 0 — baseline capture.** Record, before any change: the `ctest` set on both compilers, and a
fast-on carbon-mode smoke summary (AGB/LAI/cohorts + a deterministic slow-loop scalar) — the reference for
the Phase 1 rename bit-identity check and the Phase 2 re-baseline.

**Phase 1 — renames (pure, behavior-preserving).** `carbon_growth → compute_carbon_allocation`;
`carbon_rates → compute_vital_rates`; `compute_slow_derivatives → update_cohort_derivatives`;
`advance_trait_dynamics → advance_plant_traits` (+ its public list & the 2 external callers); and
**dissolve `compute_empirical_derivatives`** by inlining it into `meds_apply_rates` (capi-only, §4.2).
Refresh the `CLAUDE.md` references to the renamed symbols. *Verify:* both compilers, identical `ctest`,
slow-loop scalar **bit-identical** (renames + inline are behavior-preserving), fast-on carbon smoke
byte-identical, and the Python `empirical_spinup.py` still runs (the capi inline).

**Phase 2 — run model: fast always on + unconditional phenology (a behavior change).** (a) Thread
`step_start`/`doy` unconditionally (`meds_main → advance_one_step → vegetation_dynamics`); remove the two
`doy`-absent error stops. (b) Remove the `cfg%phenology_on` gate
([:73-77](../../src/driver/meds_vegetation_dynamics.f90#L73-L77)) and the now-dead OFF branch in
`cohort_carbon_demand` ([:349-354](../../src/driver/meds_vegetation_dynamics.f90#L349-L354)); retire
`phenology_on`, `PHENOLOGY_OFF_FLUSH`, and the validator precondition
([meds_config.f90:349-353](../../src/shared/config/meds_config.f90#L349-L353)). (c) Make fast the run
default and add a **master `slow_on` switch** (`slow_on=.false.` skips the whole slow tier — vegetation +
biogeochem — holding demographic + soil state static; one flag for now, not per-domain). Default
`pheno_k_flush_max` stays `0.06667` (the 15-day flush, §2 / A2). *Verify:* the fast-on carbon smoke runs;
**the slow-loop baseline changes here** (the A2 re-baseline — document the one-time shift, this is *not*
bit-identical); both compilers; confirm the frozen-slow path holds vegetation + soil static.
*Risk:* the one behavior-contract change — the re-baseline is expected, not a regression.

**Phase 3 — documentary.** Correct the split rationale comments per §5 (mode-agnostic applier + offload,
not plant-statelessness).

---

# Part II — Biogeochemistry dynamics (soil-carbon wiring)

## 7. Decision + the reserved seam

**Decision: SOUND — but a new FEATURE, not a reorg.** Wiring the CENTURY soil-carbon matrix into the slow
loop is desirable, the insertion point is pre-reserved, and the kernels are complete — but
`soil_carbon_step` has **zero driver callers today** and needs four pieces of new plumbing plus a strict
double-counting contract.

**Target module / file: `meds_biogeochem_dynamics.f90`** (module `meds_biogeochem_dynamics`, author's
choice) — a new slow-tier **domain driver**, peer of `meds_vegetation_dynamics` under a thin
`meds_slow_dynamics` coordinator (§10a), following the `*_dynamics` convention (NOT `*_time_derivs`). The
`biogeochem` abbreviation matches the existing biogeochemistry modules (`meds_soil_biogeochem`,
`meds_biogeochem_types`). `MEDS_DRIVER_REORG_DESIGN.md` §8 had a placeholder `meds_biogeochemistry_dynamics`
— now reconciled to `meds_biogeochem_dynamics` in both docs.

**The reserved seam.** The slow loop already marks where soil carbon joins: the comment "the soil-carbon
step will join here" sits at the `update_patch_states` call
([meds_vegetation_dynamics.f90:114-115](../../src/driver/meds_vegetation_dynamics.f90#L114-L115)), and
`update_patch_states` (core) is the documented per-patch apply-seam — mirroring the cohort path exactly
(driver computes the tendency; the mode-agnostic core engine applies — the §5 wall).

**What already exists (kernels, complete + reachable).** `meds_soil_biogeochem` exports the full CENTURY
machinery — `soil_carbon_step`, `heterotrophic_respiration_matrix`, `build_litter_input`,
`solve_soil_carbon_steady_state` — and `meds_aux` already links `meds_biogeochemistry`
([CMakeLists.txt:190](../../CMakeLists.txt#L190)). `soil_carbon_step(pools intent(inout), u, lignin_in,
xi_int, opts, rh_today, audit)` is a once-a-day **advance-and-commit** step: it consumes the fast loop's
day-accumulated environmental factor and commits the 7-pool state
([meds_soil_biogeochem.f90:279+](../../src/biogeochemistry/meds_soil_biogeochem.f90#L279)). *Naming
(author):* **keep the name `xi_int`** (it is the matrix-model term `ξ`, integrated over the day), but
**annotate it** at the declaration — a comment stating `xi_int = ∫_day ξ dt`, the day-integral of the
decomposition environmental scalar `ξ` (temperature × moisture limitation on decay), units [day] — so the
jargon is documented rather than renamed. The pool type `soil_carbon_t` is a 7-pool vector whose index-2
field keeps the name `fast_soil_carbon`
([meds_biogeochem_types.f90:107-123](../../src/biogeochemistry/meds_biogeochem_types.f90#L107-L123)) — the
pool the fast loop reads for its Rh, **held frozen across the day** and updated only by this end-of-day step
(§9).

This is the CLAUDE.md / README-labeled "P3" biogeochemistry item.

## 8. The four plumbing pieces (the real work)

1. **Prognostic state.** Add a per-patch `soil_carbon_t` to the patch state + the lockstep
   grow/copy/sort/fuse reorder. **DAG wrinkle:** `soil_carbon_t` lives in `meds_biogeochem_types`, but
   `meds_core` links `meds_shared` **only** ([CMakeLists.txt:119](../../CMakeLists.txt#L119) — biogeochem
   is a shared-only leaf). Putting the type directly on `patch_block` would create a forbidden
   `core → biogeochem_types` edge. **Decision (author):** relocate `soil_carbon_t` into `shared/state`
   (mirroring the `meds_column_state_types` precedent), so `patch_block` carries it with no
   `core → biogeochem` back-edge.
2. **Config + provenance + restart.** A `[soil_carbon]` TOML block, ed_params-style provenance for the
   decomposition parameters, and netCDF restart of the pools.
3. **The demography→litter seam.** `compute_carbon_allocation` (the renamed `carbon_growth`) already
   computes this step's `leaf_shed_c` / `fineroot_shed_c` litter amounts
   ([meds_vegetation_dynamics.f90:299-300](../../src/driver/meds_vegetation_dynamics.f90#L299-L300)) but
   currently **discards** them — route them (plus wood/mortality litter) into `build_litter_input` to form
   `u` / `lignin_in`.
4. **The `xi_int` accumulator** (annotated, §7). Thread the fast loop's day-integral of the environmental
   decomposition scalar (`∫_day ξ_j dt`, per pool) out of the fast loop into the daily slow step —
   accumulated over the **same frozen pools** the fast Rh respires (§9).

## 9. The double-counting contract (the dominant hazard)

The fast loop is the **sole** Rh emitter to the canopy-air CO₂ / NEE today:
`heterotrophic_respiration_flux(ccfg%fast_soil_carbon, ...)`
([meds_fast_ark.f90:841](../../src/driver/meds_fast_ark.f90#L841)), feeding `budg%nee_last`. When the slow
`soil_carbon_step` also debits the pools, the CO₂ must be emitted **exactly once**.

**Design (author) — the fast loop respires the FROZEN daily pool; the pool is updated only at day-end.**
The soil-carbon pools are held constant across the day's fast sub-steps (like the other frozen pre-pass
quantities); the fast loop accumulates its Rh over that frozen state, and the daily `soil_carbon_step`
commits the pool change **once, at day-end**, from the **same** frozen pools + the same accumulated
`xi_int`. So the day's total fast Rh **equals** the slow pool debit **by construction** — the
`rh_seam_gap` correction goes to ~0 on its own, and `rh_today` (= `litter_in − ΔC_pool`) stays
**audit-only** (never a second CAS/NEE source). This is the clean version of the reconciliation, not a
correction term bolted on afterward.

- To make that identity hold across **all** pools (not just the lumped metabolic pool the fast loop reads
  today via the single scalar `ccfg%fast_soil_carbon`), the fast Rh moves to the **matrix form**
  (`heterotrophic_respiration_matrix`) over the frozen 7-pool vector with the same K/ξ the slow matrix
  uses. The prescribed constant `fast_soil_carbon = 5.0 kgC/m²`
  ([meds_fast_types.f90:68](../../src/driver/meds_fast_types.f90#L68)) is retired in favor of the frozen
  prognostic pool.
- The reconciliation fields (`rh_fast_accum`, `rh_seam_gap`) are **stubbed today**
  ([meds_soil_biogeochem.f90:344-345](../../src/biogeochemistry/meds_soil_biogeochem.f90#L344-L345)); with
  the frozen-pool design they close to ~0 automatically, so `rh_seam_gap` becomes a cheap **assertion
  guard** rather than a live correction.

**None of this is exercised by ordinary `ctest`** (the default is fast-off) — a **fast-on carbon smoke run
is mandatory** to prove the closed whole-system CO₂/carbon budget.

## 10. DAG + phasing (Part II)

**DAG-clean.** `libmeds_biogeochemistry` links `meds_shared` only and is already linked by `meds_aux`; a
`meds_biogeochem_dynamics` driver above the vegetation driver + soil-biogeochem kernels is acyclic.
The single real edge concern is the `soil_carbon_t`-on-`patch_block` wrinkle (§8.1, resolved by the
`shared/state` relocation).

### 10a. Slow-tier orchestration — a thin `meds_slow_dynamics` coordinator (author: yes)

Today `vegetation_dynamics` calls `update_patch_states` internally (patch aging), and the reserved
soil-carbon seam sits at that same call ([:114-115](../../src/driver/meds_vegetation_dynamics.f90#L114-L115)).
Once biogeochem is a **second** per-patch slow domain, introduce a **thin `meds_slow_dynamics` coordinator**
that `advance_one_step` (`meds_stepper`) fans out to: it sequences the domain drivers
(`vegetation_dynamics` → `meds_biogeochem_dynamics`) and owns the **shared per-patch slow-state
application** hoisted out of `vegetation_dynamics`. So the two slow domains are **peers** rather than
biogeochem nesting inside the vegetation driver. This is the trigger the header deferred `meds_slow_dynamics`
to; it stays *thin* (a coordinator, not a numerical stepper) and preserves the "core applies, drivers
compute" wall (`update_patch_states` remains the shared core applier, now fed by both domains). **Do the
hoist + introduce `meds_slow_dynamics` as part of B2** (when biogeochem joins the patch seam), not before —
it is churn with no benefit while vegetation is the only slow domain.

Part II is **independent of and later than Part I** (it touches neither phenology nor the renames).
Suggested sequence — a multi-part feature PR, gated per project discipline (both compilers + a fast-on
carbon smoke at each step):

- **B0 — state + spin-up** (no fast coupling yet): add the per-patch `soil_carbon_t` state (in
  `shared/state`, §8.1) + the `[soil_carbon]` config. Support **both** initialization paths (author): a
  `solve_soil_carbon_steady_state` cold-start **and** a netCDF pool restart; verify pools persist +
  round-trip. No Rh double-count risk yet (the slow step is not emitting).
- **B1 — litter seam:** route `compute_carbon_allocation`'s `leaf_shed_c`/`fineroot_shed_c` (+ wood/mortality
  litter) through `build_litter_input` into the pools; verify litter mass balance (litter out of
  vegetation = litter into soil pools).
- **B2 — daily `soil_carbon_step`** at the hoisted per-patch slow seam (§10a), consuming `xi_int`,
  `rh_today` audit-only. Move the fast Rh to the **frozen-daily-pool** matrix form so the day's fast Rh =
  the slow pool debit **by construction** (§9); `rh_seam_gap` becomes the confirming assertion. Do the
  `update_patch_states` hoist here. **This is the double-count gate** — the fast-on carbon smoke must show a
  closed whole-system carbon budget.
- **B3 — diagnostics/output:** soil-carbon pools + Rh into the output subsystem.

---

## 11. Resolved decisions (all open questions answered)

All author decisions, folded into the plan:

- **Phenology (§2):** unconditional; realistic **15-day flush** the default (A2).
- **Run model (§3):** fast biophysics **always on**; phenology gate + validator removed; a **master
  `slow_on` switch** freezes the whole slow tier (vegetation + biogeochem) — one flag, not per-domain;
  slow-only / empirical runs = Python C-API.
- **Renames (§4):** `carbon_growth → compute_carbon_allocation`, `carbon_rates → compute_vital_rates`,
  `compute_slow_derivatives → update_cohort_derivatives`, `advance_trait_dynamics → advance_plant_traits`.
- **Carbon split (§5):** **kept** (mode-agnostic applier + offload). **Empirical twin (§4.2):** **dissolve**
  `compute_empirical_derivatives`, inline into `meds_apply_rates` (capi-only; no twin).
- **`soil_carbon_t` → `shared/state` (§8.1).** **`xi_int` name kept + annotated** (§7).
- **Fast Rh = frozen daily pool**, committed end-of-day → budget closes by construction (§9).
- **Both** spin-up paths — steady-state cold-start + netCDF restart (§10 B0).
- **Slow-tier orchestration (§10a):** a **thin `meds_slow_dynamics` coordinator** (introduced in B2) owns
  the domain sequence + the hoisted `update_patch_states`, making vegetation and biogeochem peers.
- **Biogeochem module = `meds_biogeochem_dynamics`** (§7); reconciled in `MEDS_DRIVER_REORG_DESIGN.md`.

No open questions remain — the plan is ready to implement (Part I first).
