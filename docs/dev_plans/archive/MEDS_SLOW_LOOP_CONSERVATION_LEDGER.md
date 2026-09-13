# MEDS slow-loop conservation ledger — the measurement record

> # 🗃️ ARCHIVED — 2026-09-13. The working record behind PRs #132–#137.
>
> **This was §10.2 of `MEDS_CODE_STRUCTURE_DESIGN.md`**, 549 lines of it — the largest single
> section of that plan, and self-contained. It is the item-by-item account of building the slow-loop
> conservation ledger and of every leak the ledger then found: growth respiration that was 23 % of
> growth carbon and never exhaled, tissue thermal mass appearing rather than being exchanged,
> time-averaged diagnostics that did not close, mortality valued on what was requested rather than
> on what was removed, and reproduction carbon that simply vanished.
>
> **Its sub-numbering is cited from source comments** (`§10.2.11`, `§10.2.12`, `§10.2.13`,
> `§10.2.2 items 1–3`), so the original numbering is kept verbatim.
>
> **Shipped as PRs #132–#137** (2026-09-09/10). The lesson worth carrying forward is in
> `CLAUDE.md`: a closed budget proves bookkeeping, not plausibility — these leaks all sat behind
> ledgers that balanced, because the ledger declared the same quantity it consumed.
>
> **Live description:** `docs/science/` has no slow-loop conservation page; the ledger itself is
> `src/slow_dynamics/driver/meds_slow_ledger.f90`, whose header states the invariants.

### 10.2 Slow-loop conservation (review item 1B #4–#10) → **after migration step 6, in `slow_dynamics/`**

Physics, so outside this plan's "no numerics" scope — but the *placement* is this plan's business, and
the order matters: these land after step 6 so they are written once into their final home instead of
being moved a week later. **Surveyed against the code 2026-09-09**; what that survey found is recorded
below, because several of the review's original items turned out to be larger, smaller or differently
shaped than the review's one-line description of them.

#### 10.2.0 Why none of this is currently visible

`budget_t` lives in `util/meds_budget_check`, **not** `base/` (an earlier draft of this section said
`base/`; the ledger reuses it from `util/`). It accumulates **per-fast-step flux residuals**, which
`fast_dynamics` merges up into the run-level `run_energy_budget` / `run_water_budget`. Nothing anywhere
compares a **store** across the slow step. A discontinuous jump between the end of one fast window and
the start of the next is therefore invisible **by construction**, not by oversight — which is why every
item below has survived a full review and a conservation-focused PR.

Two consequences for the design:

- The slow ledger is a **snapshot** ledger — store before, declared boundary terms, store after — not a
  flux accumulator. It is the fast loop's ledger turned inside out.
- It must be **site-level and area-weighted** (`Σ_p area_p × store_p`). Patch identity does not survive
  the step: `apply_patch_disturbance` creates a patch, `fuse_2_patches` destroys one, and
  `terminate_patches` renormalizes every area. A per-patch ledger has nothing to compare against.

#### 10.2.1 The stores

`meds_fast_ark`'s `whole_water` / `whole_energy` check already names the canonical list; the slow ledger
spans the same stores so the two tiers cannot disagree about what exists, and adds carbon, which the fast
whole-column ledger does not track.

| currency | stores (each area-weighted to the site) |
|---|---|
| water  | soil column; pond `w_surface`; snow `swe`; CAS vapour `cas_mass_capacity·shv`; tissue water (`leaf_water_mass`+`wood_water_mass`, ×nplant); interception films (already ground-referenced — **not** ×nplant) |
| energy | soil energy; pond enthalpy; snow enthalpy; CAS `cas_mass_capacity·enth`; tissue heat |
| carbon | live pools (leaf/fineroot/wood/nonstructural, ×nplant); the 7 CENTURY pools; CAS CO2 (`cas_molar_capacity·can_co2`); `recruit_pool` as carbon-in-transit |

#### 10.2.2 Carbon terms, largest first

1. **Growth respiration is destroyed, not exhaled** (this is 1B #5, and it is much the largest item in
   this section). `plant_carbon_allocation` charges `growth_resp` against the plant's carbon and its own
   header states the identity it must satisfy; nothing consumes it — `growth_resp` reaches only the
   `CS_GROWTH_RESP` diagnostic. Meanwhile `nee_biotic = ra_leaf + ra_stem + ra_root + rh − gpp` carries
   **maintenance** respiration only. With `growth_resp_factor = 0.3` that is `0.3/1.3` = **23 % of every
   unit of carbon entering growth**, debited from the plant and never reaching the atmosphere. Same
   defect class as the deleted soil-carbon fallback (#128), opposite sign: the CAS runs too low, so
   photosynthesis is under-fertilized and the site reports as a larger sink than its own carbon flow
   implies.
2. **The starvation `deficit` is the mirror image.** When storage cannot cover maintenance, the fast loop
   has *already* respired that carbon into the CAS, but no pool is debited. Carbon from nothing, per
   cohort, on every step a cohort starves.
3. **Recruitment endows ~5–6× what it debits** (1B #4). The debit is `min_cohort_carbon` =
   `dbh_to_agb(...)`, **AGB only**; `init_cohort` → `set_cohort_size` endows `wood_carbon = agb /
   aboveground_frac` **plus** leaf, fineroot and storage. Measured on the shipped 3-PFT table
   (`min_cohort_height` 2 m): total/agb = 6.43, 5.52, 5.18. `repro_carbon_efficiency = 1e-3` then
   over-corrects, so **net ~99.4 % of reproduction carbon vanishes**. Both halves need declaring: the
   establishment loss is physically real but it is necromass, not nothing, and the endowment/debit ratio
   should be 1.
4. **Silent floors.** `max(pool, 0)` appears in both `update_cohort_derivatives` and
   `update_cohort_states_kernel`; `nplant` is floored at `negligible_nplant`. Each is creation. The
   nplant floor also **double-counts**: `accumulate_mortality_litter` values litter on the *unfloored*
   `died_nplant`, so litter is produced for individuals the floor then refuses to kill.
5. **Litter is gated on `soil_carbon_on`** in four places (turnover, continuous mortality, cull,
   disturbance kill). With soil carbon off — the default — all necromass vanishes. Defensible as an
   *export*, but after #128 it is declared as one rather than left implicit.
6. **Operator-split offset.** Mortality litter is valued on pre-growth pools while the nplant decrement
   lands on post-growth pools; the error is `died_nplant × npp`, one-signed.
7. **External seed rain** (`seed_rain_recruits`, fully endowed recruits) is a legitimate boundary
   **import**. Declared, not eliminated.
8. **Reproduction is aliased, not leaked** — recorded here so it is not mistaken for a leak later.
   `npp_repro` is debited every slow step; `apply_recruitment` runs only under `is_new_month` and treats
   that one day's rate as the whole month. Unbiased for a steady rate, but reproduction tracks NPP
   seasonally, so it is a 12-sample estimator of a daily flux. A correctness question for the demography
   cadence, not a ledger term.

#### 10.2.3 Water terms

- `shed_turnover_water` → `patch%shed_water_rate` → the fast loop's `precip_ground` is **the model of
  what every other seam here should look like**: the mass leaves one store, a named variable carries it,
  and another store receives it. Already closed; the ledger just declares it.
- `reconcile_tissue_water_capacity` **already returns `seeded` and `discarded`**, and its header already
  says they are booked nowhere and become ledger terms when this lands. Free.
- **Cull and disturbance-kill water is discarded** (1B #7). `terminate_cohorts` routes carbon to litter
  and drops `leaf_water_mass`, `wood_water_mass` and both films; `apply_patch_disturbance` does the same
  for the killed canopy. A pure sink.

#### 10.2.4 Energy terms — and two that are bugs, not bookkeeping

- **`blend_cas` is wrong for patches of different height** (1B #9, and it needs fixing whether or not the
  ledger lands). `can_enthalpy`, `can_shv` and `can_co2` are *specific*; the extensive content is
  `area × depth × ρ × value`. The function area-weights the intensive values while *separately*
  area-weighting `can_depth`, so fusion drops the covariance term `w₁w₂(d₁−d₂)(h₁−h₂)`. Fusing a tall
  patch with a gap corrupts canopy-air energy, humidity **and** CO2 at once.
- **`cas_set_depth`'s open-volume term is computed nowhere.** The routine already carries the `de_open`
  hook and a header explaining that the jumps that matter are disturbance and fusion — a 20 m canopy
  becoming a 1 m gap in one step. `refresh_canopy_depth` calls it without `rho_air`/`de_open`, so the
  term is never formed.
- **Fusion weights both tissue temperatures by leaf area.** `fuse_cohort_fast_state` blends `leaf_temp`
  *and* `wood_temp` on leaf area. For leaves that roughly tracks heat capacity; for **wood** it does not
  — wood heat capacity follows sapwood+fine-root carbon and the water in it, which does not scale with
  leaf area. And `leaf_water_mass` is nplant-weighted while its temperature is leaf-area-weighted, so
  their product — the actual energy — is conserved by neither weighting.
- Shed water leaves tissue at tissue temperature and re-enters the ground at `t_film_valuation`
  (deliberate, and documented as such — still an energy term).
- Cull and disturbance-kill tissue heat, alongside the water above.

#### 10.2.5 Tolerances

`fuse_2_cohorts` asserts its carbon pools at `conservation_tol` (`1e-3` in the shipped configs). That
operation is exact algebra and should hold at round-off — `split_cohorts` already asserts its film water
at `1e-12`. A 0.1 % window on every fusion is eleven orders too loose. Per 1B #10: fuse/split assert at
~1e-12 relative, `conservation_tol` is **reserved for the site ledger**, where the declared terms
genuinely carry approximation, and `size_tol` is deleted.

#### 10.2.6 Implementation order

**The skeleton lands before any fix.** Snapshot the stores, declare the terms that are already honest
(shed water, reconcile seed/discard, seed rain, litter export), and route everything else into an
explicit per-currency `unattributed` bucket. That makes each gap *measurable* before deciding what it
deserves — the lesson from #128 and #129, where measurement inverted the expected ranking twice. The
prediction on record is that growth respiration dominates the carbon bucket by an order of magnitude and
that `blend_cas` dominates energy only in runs with active patch fusion; the numbers, not this
paragraph, decide the order of the fixes that follow.

Then, in their final homes:

- `slow_dynamics/driver/meds_slow_ledger` — `site_ledger_t`, the snapshot, and the daily/annual assert.
- `slow_dynamics/plant/`: recruitment carbon debit vs `init_cohort` endowment (1B #4); `growth_resp` and
  `deficit` routed into the CAS carbon balance (1B #5).
- `slow_dynamics/demography/`: mortality, cull and disturbance-kill hand tissue water, film water and
  tissue heat to a declared sink (1B #7); `blend_cas` depth blending and the `terminate_patches`
  survivor rescale fixed (1B #9); the fuse/split tolerances of 10.2.5 (1B #10).
- `slow_dynamics/soil/`: CENTURY transfer-matrix column conservation asserted; `rh_seam_gap` asserted in
  production; `audit%litter_in` includes cull/disturbance litter (1B #8, #6).
- Not structure at all, listed so it is not lost: `docs/ed2_comparison.md` soil-temperature statements
  predate the PR #119 root-heat-sink fix and must be re-run.

#### 10.2.7 What the skeleton measured (2026-09-09)

The skeleton shipped and ran. Numbers are **site totals per m², cumulative over a 1096-day (3-year)
Ithaca run** from a small establishing stand, ifx, byte-identical to `main` on all 75 outputs. Read
`residual` against `declared`, which is the gross boundary flux the phase did account for.

| phase | carbon [kgC] | declared | water [kg] | energy [J] | marks |
|---|---|---|---|---|---|
| allocate         | –          | 9.1e-4 | **+1.3e-12** | **+1.3e-6** | 1096 |
| grow + mortality | **−3.0e-3** (−3.4e-4 with soil C on) | 2.5e-3 (5.2e-3) | **−4.0e-4** | **−272** | 1096 |
| recruit          | **+7.7e-3** | 0 | 1e-13 | **+1.08e5** | 36 |
| cohort fuse/fiss | 3e-18 | 0 | 1e-13 | **−6.8e4** | 36 |
| disturbance      | 2e-18 | 0 | 5e-13 | **+3345** | 3 |
| patch fuse/term  | 2e-17 | 0 | 1e-13 | **−3345** | 3 |
| canopy depth     | **−9e-18** | 3.6e-3 | **+7e-13** | **+3.6e-7** | 1096 |
| soil carbon      | **+9.0e-19** | 3.5e-3 | – | – | 1096 |

**Three phases are verified closed**, which is what makes the rest trustworthy: the canopy-depth
open-volume declaration closes to 1e-6 J against **5.77e6 J** declared — eleven orders — so the
entrainment term is both large and now fully accounted; the CENTURY step closes to 9e-19 against
3.5e-3 declared once `rh_today` is declared; and the turnover-shed seam closes to 1e-12.

**The prediction in 10.2.6 was wrong, and this is the record of it.** It said growth respiration
would dominate the carbon bucket by an order of magnitude. It does not. **Recruitment creates
+7.7e-3 kgC/m², twenty-three times the growth phase's −3.4e-4**, and it does so in 36 events rather
than 1096. Two honest caveats before that is read as settled: this is a 3-year *establishing* stand,
where recruitment is proportionally far larger than in a mature forest, and the growth-phase figure
is a *net* of terms with opposite signs (growth respiration destroys, the starvation deficit and the
floors create). Both need a mature run and a decomposition before the fix order is set. What is not
in doubt is that the ranking was not the one predicted.

**With soil carbon off — the default — the necromass export is 2.7e-3 kgC/m²**, the difference
between the two grow-phase figures, and the largest single carbon term in a stock configuration.
That is item 5 of 10.2.2, now measured. With soil carbon on, the growth phase loses **6.4 % of the
carbon handed to it by the fast tier** (−3.4e-4 against 5.2e-3 declared).

**Two findings the survey did not predict:**

- **Recruits are born with tissue heat from nowhere: +1.08e5 J over 36 events.** `init_cohort` sets
  `leaf_temp`/`wood_temp` to `LEAF_TEMP_INIT` and `set_cohort_size` immediately gives the cohort a
  heat capacity, so a recruitment event creates sensible heat in proportion to the biomass it also
  creates. This is the energy twin of 10.2.2 item 3 and was not on the list.
- **Disturbance (+3345 J) and patch fusion (−3345 J) are exactly equal and opposite**, to every
  digit printed. The gap patch's canopy air is created by `blend_cas` and reabsorbed by it; the sign
  symmetry says the two are the same arithmetic run forwards and backwards, which is a useful
  constraint on any fix to `blend_cas`.

Cohort fusion's **−6.8e4 J** is 10.2.4's leaf-area-weighted temperature blend, measured; mortality's
**−4.0e-4 kg** is 1B #7's discarded tissue water, measured.

**One correction to this section's own design.** `shed_water_rate` is *not* a store, though 10.2.3
implies it is. It is a handoff whose lifetime is the FAST window, not the slow step: by the time the
next slow step opens, that water is already in the soil and the rate variable still holds it, so
carrying it as a store double-counts at every open — a steady one-signed ~8e-7 kg/step phantom leak
in the allocate phase, which is exactly how it was found. It is declared instead, the mirror of the
fast→slow carbon handover.

**Deferred from the skeleton, on purpose.** `slow_site_store` duplicates four store loops from
`meds_column_state_ops` (soil water, soil energy, plant water, canopy film) and the tissue
heat-capacity construction from `meds_fast_frozen`, because those live in `meds_fast` and `meds_slow`
must not link it. Moving them down to `state/column`, so both tiers value the stores with one piece
of code rather than two that agree today, is the natural next commit and belongs with the energy fix.

#### 10.2.8 The mature stand overturns 10.2.7's ranking (2026-09-09)

10.2.7's headline — that recruitment dominates the carbon bucket — was measured on a 3-year
*establishing* stand and 10.2.7 said in as many words that the stand flattered recruitment. It did.
Re-run from the `runs/ithaca_ark30` 2074 spin-up restart (114 cohorts, LAI 5.6, AGB 16.7 kgC/m²,
mean dbh 37 cm, 1 PFT, soil carbon on), same 3-year span, same ledger:

| | establishing | **mature** |
|---|---|---|
| grow+mortality carbon | −3.02e-3 kgC/m² | **−2.5102 kgC/m²** |
| …as % of declared | 119 % | **25.3 %** |
| recruit carbon | +7.67e-3 | +7.61e-3 |
| **grow : recruit** | **0.39 : 1** | **330 : 1** |

The establishing stand nets **0.85 gC/m²/yr**; the mature stand **837 gC/m²/yr lost in the growth
phase against ~2925 gC/m²/yr GPP**. Seed rain is a fixed 0.01 plant/m²/yr, so on an unproductive
stand it swamps everything and on a real one it is noise. **10.2.6's original prediction was right:
growth respiration leads, and the measured 25.3 % sits right beside the `g/(1+g)` = 23.1 % the
construction cost implies.** The lesson is not about growth respiration; it is that a conservation
ranking measured on a stand that is not growing measures the stand, not the model.

The soil-carbon phase closes to **−1.1e-13 against 3.41 kgC/m² declared** on a productive forest,
and the canopy-depth phase to 1e-6 J against 4.42e6 J. Both hold at the mature scale.

#### 10.2.9 Birth and death become paired transfers (2026-09-09)

Implementing the author's framing: a recruit **draws** what it arrives with from outside the model,
and a death **hands on** what it carried. Measured on the mature stand, before → after:

| term | before | after | |
|---|---|---|---|
| grow+mortality water | −1.2325 kg | **+6.0e-12** | closed |
| disturbance water | −0.8852 kg | **−3.8e-6** | closed |
| recruit energy | +9.54e4 J | **−7.7e-7** | closed |
| disturbance energy | −3.498e6 J | **+9.9e3** | 99.7 % |
| recruit carbon | +7.61e-3 | **+9.60e-4** | 87 % |
| cohort fuse energy | −8.38e4 J | −8.30e4 J | unchanged (weighting, not a transfer) |
| grow+mortality energy | +5.18e6 J | **+1.00e7 J** | *exposed*, see below |

**Why birth draws externally.** Recruitment stands in for everything between a seed and a 2 m
sapling — germination and the seedling's own photosynthesis, transpiration and energy balance —
and the model tracks no cohort below `min_cohort_height`. What a recruit arrives with was fixed and
absorbed by a size class that is not represented, so it is genuinely external. Drawing it from the
free atmosphere rather than the patch's canopy air is deliberate: crediting the CAS with seedling
uptake while representing none of the seedling's respiration, transpiration or shading would add one
term of a missing process and call it an improvement.

Only the part the model did **not** already pay for is declared. `recruit_pool` carries reproduction
carbon at `carbon_min` per plant and that much *is* debited from the parents, so the draw is
`endowment − carbon_min`, plus the baseline seed rain at its true entry point (the monthly pool
credit — which arrives whether or not anything is born that month, a distinction worth 2/3 of the
term). The residual 9.6e-4 that remains is the productivity-driven reproduction carbon: debited from
parents in the growth phase, re-created here at a quantity the `repro_carbon_efficiency / carbon_min`
conversion does not preserve. That is 10.2.2 item 3, deliberately left visible rather than declared
away.

**Birth temperature was the real energy bug.** `init_cohort` stamped `LEAF_TEMP_INIT = 288.15 K`, one
global constant, so a sapling appearing in an Ithaca January was born ~20 K warmer than the air it
stood in. Recruits now start at their patch's canopy-air temperature (guarded: an unstepped CAS falls
back to the constant). Declaring the old term would have been declaring an artifact.

**Death hands its water to the ground down the channel turnover shedding already uses** — one
verified path rather than a second mechanism — for all three death paths (background mortality, the
cull, the disturbance kill). The tissue **heat** leaves the thermal system with the necromass and is
reported, because the CENTURY pools it becomes carry no temperature; a litter thermal store is where
it would belong if one existed.

**What this exposed.** The growth phase's energy residual *rose*, from +5.18e6 to +1.00e7 J, because
mortality's heat was partly cancelling it. Growing biomass raises the tissue heat capacity at
constant temperature, so `cap × T` rises with no flux: **the model creates thermal mass out of
carbon**. That is the birth-side twin of the death-side term just fixed, it is now the largest energy
item in the ledger, and it was invisible until the offsetting term was removed.

**Still open after this**, in measured order: growth respiration and the growth phase's carbon
(−2.51 kgC/m², 25 % of the handover); the growth-side thermal mass (+1.00e7 J); cohort fusion's
leaf-area-weighted `wood_temp` blend (−8.3e4 J) and `blend_cas`; the item-3 reproduction conversion
(+9.6e-4).

#### 10.2.10 The allocator's outputs all get a destination (2026-09-09)

`plant_carbon_allocation` produces seven outputs. Four are pools the driver commits. **Three left the
plant and went nowhere**, and they are one sentence, not three problems: the allocator's outputs had
no destinations. All three are fixed together.

- **Growth respiration → the canopy air.** Charged against the plant, reaching only the
  `CS_GROWTH_RESP` diagnostic, while `nee_biotic` carried maintenance respiration alone. It is now
  handed down as a frozen daily rate on a new per-patch `slow_co2_rate`, added to `nee_biotic` by the
  prepass — the carbon twin of `shed_water_rate`, built to the same read-only, seeded-once discipline.
- **The starvation `deficit` → the same channel, opposite sign.** The fast loop has *already* exhaled
  the full maintenance respiration; when storage could not fund it, `deficit` is the part no pool paid
  for. Real maintenance respiration is substrate-limited, so the honest reading is that the fast loop
  **over-reported**, and an operator-split model corrects an over-report on the next step rather than
  rewriting the last one. Netting it into one signed channel keeps one sign convention.
- **The unestablished seed fraction → litter.** Reproduction carbon is debited in full from the parent
  and only `repro_carbon_efficiency` of it establishes. The remaining ~99.9 % is dead seed and dead
  seedling — **necromass, not nothing** (author's decision). It enters `necromass_to_litter` as
  `storage_c`, which pools with the canopy and splits on `f_labile_leaf`: seed tissue is labile and
  canopy-derived, which is what that argument means.

**Result on the mature stand.** The growth phase's carbon residual falls from **−2.5102 to −7.337e-4
kgC/m²**, a factor of **3420**, and the worst single step from 5.06e-3 to 1.81e-6. What remains is
0.006 % of the declared flux: the pool and `nplant` floors (item 4) and the pre/post-growth offset in
the mortality valuation (item 6). The soil-carbon phase still closes at 7e-15, now against 4.50
declared rather than 3.41 — the seed litter is real carbon arriving in real pools.

**The NEE shift is the missing respiration, to within 0.4 %.** Site NEE moves from −5.988 to −4.616
µmol/m²/s, **+1.372**. The growth respiration the ledger said was absent is ~1.55 kgC/m² over 3 years
= 0.517 kgC/m²/yr = **1.366 µmol/m²/s**. Two independent routes to the same number: the ledger's
residual before the fix, and the flux difference after it.

| | main | branch | |
|---|---|---|---|
| `nee_site` | −5.988 | **−4.616** | +22.9 % — the site is a much smaller sink |
| `soilc_total_site` | 1.775 | 2.181 | +22.8 % — seed necromass accumulates |
| `rh_site` | 0.01091 | 0.01555 | +42.6 % — and decomposes |
| `gpp_site` | 0.243774 | 0.243946 | **+0.070 % — the CO2 fertilization feedback** |
| `agb_site` / `lai_site` / `nplant_site` | — | — | +0.005…0.007 % over 3 yr |

The GPP rise is small but it is the feedback predicted when the soil-carbon fallback was deleted in
#128, running the other way: putting CO2 back into the canopy air raises photosynthesis. Structure
barely moves, which is right — the carbon *through* the plant is unchanged; only its fate after
leaving is.

**A test the ledgers cannot replace.** Neither ledger can catch a break in this channel: the slow
ledger declares the **handoff**, so it closes whether or not the fast loop ever picks the rate up, and
the fast CAS ledger closes around whatever `nee_biotic` it is handed. Only a test binds the two ends.
`test_slow_ledger` gains three differential assertions — growth respiration reaches the channel and a
zero construction cost empties it, the unestablished seed fraction reaches litter, and a **starving**
stand owes a **negative** flux. All three mutation-tested.

Two flaws in the first version of that test, worth recording because both made it pass while asserting
nothing: it called `vegetation_dynamics` twice on one site, and since that call *commits* growth the
two variants compared different forests; and its fixture cohort was below `min_reproduction_height`,
so no seed carbon existed to lose. Each variant now runs on a fresh stand, sized above the threshold
and given a leaf lifespan long enough that turnover does not consume the supply before reproduction is
reached.

**Still open**, in measured order: the growth-side thermal mass (+1.00e7 J — biomass growth raises the
tissue heat capacity at constant temperature, the birth-side twin of the death-side term 10.2.9
fixed); cohort fusion's leaf-area-weighted `wood_temp` blend (−8.3e4 J) and `blend_cas`; item 3's
remaining half, the recruit pool as a carbon quantity, together with the monthly-sampling aliasing
(+9.6e-4); the pool and `nplant` floors (−7.3e-4).

#### 10.2.11 Tissue thermal mass (2026-09-09)

A cohort's heat capacity is a function of its biomass and its density, so growing or dying changes
`cap × T` **with no flux at all**: the model makes thermal mass out of carbon, and unmakes it. Both
directions were undeclared and they partly cancelled, which is why the growth side stayed hidden
until 10.2.9 declared the death side and the residual *rose* from +5.18e6 to +1.00e7 J.

**Declared as one exchange, not two.** The whole change across the growth commit is one mechanism,
and splitting it into a growth part and a mortality part needs a cross-term convention the physics
does not supply (the `hcap_min` floor is not linear in density). Both `shed_mortality_water` and
`update_cohort_states` leave the temperatures alone, so the change across them is *purely* thermal
mass. 10.2.9's separate mortality-heat term is therefore withdrawn — subsumed, not repealed — and its
water term stays, because water going to the ground is a transfer and this is not.

**Why it is an exchange and not a leak.** New tissue is assembled from CO2 and water at the plant's
own temperature and arrives carrying the sensible heat of that mass. The model tracks no thermal
content for CO2 — the canopy air's capacity is dry air alone — so that heat genuinely crosses the
boundary of the modelled thermal system. A fuller treatment would give CO2 a heat capacity in the CAS
and carry enthalpy on root water uptake; both are far larger than this, and naming the approximation
here is better than burying it.

**Result: +1.00e7 → +2.06e-6 J against 8.00e6 declared.** Twelve orders. The growth phase now closes
on all three currencies to round-off except carbon's −7.3e-4 (the pool and `nplant` floors).

**Byte-identical**, because nothing here changes the model — `slow_tissue_heat` is a pure read, and
the only other change is deleting an output that is now computed elsewhere. Mutation-tested: with the
declaration removed the phase carries +5.20e6 J, which is the *net* of the two directions and exactly
the figure 10.2.7 reported before the death side was declared.

**Also folded in:** `slow_site_store` carried its own copy of the tissue heat-capacity formula from
before the shared one existed. It now calls `cohort_tissue_heat_capacity`, so the ledger, the
demography operators and the slow driver cannot drift on what a cohort's thermal mass is. Two of the
duplications 10.2.7 flagged are gone; the four store loops in `meds_column_state_ops` remain.

**Still open**, in measured order — and the character of the list has changed. What is left is no
longer *missing transfers* but **wrong averages**: cohort fusion's leaf-area-weighted `wood_temp`
blend (−8.3e4 J), `blend_cas`'s depth-blended intensive quantities (disturbance +9.9e3 J and patch
fusion −1.0e4 J, still equal and opposite), item 3's remaining half (+9.6e-4), and the pool and
`nplant` floors (−7.3e-4). Those need a fix to the weighting, not a new declaration.

#### 10.2.12 The averages (2026-09-09)

Three wrong averages, and the thermal-mass declaration that 10.2.11 established extended to the
structural phases. **Energy and water now close everywhere in the ledger; only two carbon terms
remain.**

**`blend_cas` weights on AIR MASS, not ground area.** Enthalpy, specific humidity and CO2 mixing
ratio are per kg of air, and the air mass is `area × depth`, so the conserving weight is `area ×
depth`. Weighting on area alone dropped the covariance `w₁w₂(d₁−d₂)(v₁−v₂)`, corrupting canopy-air
energy, humidity **and** CO2 together by the same relative error. It was exact only when the two
depths matched — precisely the case that needs no blend. `rho` cancels (site-uniform), and
`can_depth` stays area-weighted because that is what a depth is.

**Cohort fusion weights the tissue temperatures on HEAT CAPACITY.** They were leaf-area weighted,
which is roughly right for leaves — leaf capacity tracks leaf carbon, and leaf area tracks that
through `sla` — and not even approximately right for **wood**, whose capacity follows wood carbon
and the sapwood ring. Capacity weighting conserves `cap × T` exactly for the additive part, since
every carbon pool and both tissue waters are nplant-weighted and so add across the merge.

The diagnostic twins keep their leaf-area weights, and 10.2.4's claim that the two tables share a
weight is withdrawn: a per-leaf-area *diagnostic* really is leaf-area weighted; a prognostic
*temperature* is per unit heat capacity. Treating them as one kind is what made the wood wrong.

**`terminate_patches` MERGES the sliver instead of deleting it.** It used to drop any patch under
`min_patch_area` and renormalize the survivors' areas back to 1, which silently redistributed the
whole site: every conserved quantity changed by `(a/(1−a))·Σ_kept a_i X_i − a·X_dropped`, zero only
if the doomed patch held the survivors' mean — and a doomed patch is atypical by construction,
usually a fresh gap. It is now fused into the largest survivor with `fuse_2_patches`, which
conserves exactly and reuses the operator patch fusion already depends on. The renormalisation
becomes a round-off correction rather than a redistribution.

**And the structural phases now bracket their thermal mass**, as the growth commit does. Fusion and
fission re-derive the sapwood ring from the merged or perturbed diameter, and `sapwood_fraction` is
**nonlinear in dbh**, so the merged capacity is not the sum of the two even when every carbon pool
is. That is a thermal-mass change of exactly the kind growth makes. It is declared *only after* the
weighting was fixed — otherwise the declaration would have been hiding the wrong average rather than
accounting for what the right one leaves behind. 10.2.9's separate cull and disturbance-kill heat
reports are withdrawn into it, as 10.2.11's mortality report already was; their **water** reports
stay, because water reaching the ground is a transfer.

| phase, energy [J] | before | after |
|---|---|---|
| cohort fuse/fiss | −8.30e4 | **−1.30e-6** |
| disturbance | +1.01e4 | **+1.49e-7** |
| patch fuse/term | −1.00e4 | **+3.43e-7** |

Water at those phases closed too: patch fusion −2.2e-5 → 3.4e-13, disturbance −3.8e-6 → 1.1e-13,
and `patch fuse/term` carbon 8.3e-8 → −1.5e-11.

**What is left in the entire ledger is two carbon terms**, both known and both named here already:
the recruit pool's productivity-driven credit (+9.61e-4, item 3's remaining half, needing the pool
to become a carbon quantity and the monthly sampling to stop aliasing) and the growth phase's pool
and `nplant` floors (−7.34e-4, items 4 and 6). Energy and water close on every phase.

**A bug this nearly shipped with.** The first version of the sliver merge rebuilt the CSR map once
after the loop. `fuse_2_patches` reads the *receptor's* CSR slice to rescale its cohorts, and
`patch_fuse_pass` — the existing caller — rebuilds after **every** fusion for exactly that reason.
With two slivers the second merge would have read a stale slice. The test now uses two.

**And a blind fixture.** `test_patch`'s CAS-fusion assertion had both patches at the default 20 m
depth, where the area-weighted and mass-weighted answers agree exactly. It passed against the wrong
code and would have passed against the right one. It now uses 30 m and 10 m, checks the
mass-weighted value, and separately asserts the extensive content survives — the property, not the
formula. `test_fusion_cohort` was rewritten the same way: it asserted the old leaf-area formula, and
now asserts that leaf tissue energy is conserved **exactly** and wood to within the sapwood
re-derivation.

#### 10.2.13 Mortality is valued on what the applier removed (2026-09-09)

The live carbon store is `nplant × pool`, and its change over the growth commit decomposes exactly:

    n₁p₁ − n₀p₀  =  n₀(p₁ − p₀)  +  (n₁ − n₀)p₁

— **growth at the old density, mortality at the new pools.** Every declaration was on the other
side of both terms: `accumulate_mortality_litter` valued the litter on the pre-growth pools `p₀`,
and the fast→slow handover was taken after the commit, at `n₁`. Both now match the decomposition:
the mortality routines run **after** `update_cohort_states` on the committed pools, and the handover
is captured **before** it.

The density drop is taken as `n₀ − n₁` rather than recomputed from the Camac hazard, which fixes
item 4's other half for free: when the `negligible_nplant` floor stops a cohort dying, `n₀ − n₁` is
smaller than the hazard implies, so litter is no longer credited for individuals still standing.

**What this leaves is one term, not several.** The growth phase's carbon residual is now
**−9.170e-4** and the recruit phase's **+9.607e-4** — the same reproduction carbon, debited from
parents in one phase and re-created in the other, with nothing connecting them. They sum to
**+4.4e-5**, which is the seed rain, the pool floors and the monthly-sampling asymmetry. Energy and
water close on every phase; the *entire* remaining ledger is item 3.

**Two mistakes in getting here, both recorded because both were invisible to everything else.**

The first estimate had item 6's sign backwards. Re-basing the litter on `p₁` makes the declared
export *larger* (the pools grew), which pushes the residual positive; re-basing the handover on `n₀`
makes the declared import larger, which pushes it negative — and the handover term dominates. The
measured residual moved from −7.34e-4 to −9.17e-4, in the direction the arithmetic says once both
changes are counted rather than one.

The second was a real bug that survived a green suite: moving `shed_mortality_water` after the
commit left it calling `cohort_tissue_water`, which carries the cohort's *current* nplant — so the
loss was scaled by `n₁/n₀` instead of being `(n₀−n₁) × per-plant`. The growth phase's water went
from 1e-12 to **1.22e-4** and nothing but the ledger noticed. `test_slow_ledger` now asserts the
identity directly — every kg that leaves tissue arrives in the patch shed channel — and both that
bug and the no-routing case fail it.

#### 10.2.14 Reproduction carbon becomes a flow, and §10.2 closes (2026-09-09)

The last item, and it needed **no new state**. `recruit_pool` was credited inside
`apply_recruitment`, monthly, from whatever recruitment rate the driver had computed on that one
day, scaled up to stand for the month. Crediting it **every step instead** fixes both halves at
once:

- the **cadence**: a 12-point sample of a quantity the model computes 365 times a year, whose value
  depended on which days happened to be month boundaries, becomes the exact integral;
- the **carbon link**: `recruitment × dt_yr` is `n·npp_repro·efficiency / carbon_min`, so the pool
  valued at `carbon_min` is exactly the establishing share of the reproduction carbon the parents
  were debited for *in that same step*. Debit and credit land in the same phase, and the ledger sees
  a transfer rather than carbon vanishing in one phase and appearing in another.

The seed-rain declaration moves with it, from monthly to daily.

**And the gap inherits the seed bank.** `apply_patch_disturbance` zeroed `recruit_pool` on the new
gap while inheriting `soil_carbon` and `xi_accum` from the same donors. A treefall gap does not
sterilise the ground it opens: the carry-forward pool sits in the soil with the litter and the
CENTURY carbon. Zeroing it destroyed the pool's carbon on the disturbed fraction — the last
non-round-off term in the ledger, **−2.96e-6 kgC/m²** over three events.

| carbon [kgC/m²] | before 10.2.13 | after 10.2.13 | **now** |
|---|---|---|---|
| grow + mortality | −7.34e-4 | −9.170e-4 | **−1.98e-13** |
| recruit | +9.61e-4 | +9.607e-4 | **+1.30e-14** |
| disturbance | −4.02e-6 | −2.96e-6 | **−1.13e-11** |

**§10.2 is closed.** Every phase closes on every currency to round-off; the largest residual
anywhere in the ledger is 1.5e-11 against declared fluxes of order 1 to 12. The report now ends in
a **verdict**: each phase and currency is judged against the flux it declared plus a per-mark
absolute floor, and the run says so in one line. Reported, not fatal — the same choice the fast
loop's whole-column ledgers make, because a run that breaches this is telling you something and
stopping it half-way tells you less than finishing it.

**One correction to what 10.2.13 said was coming.** It called the monthly sampling a case where
"353 of 365 days' reproduction carbon is never sampled". That overstates it: the estimator is a
12-node rectangle rule targeting the annual total, not a wholesale loss, so the error was quadrature
rather than a leak — and the conservation break was the separate fact that the debit and the credit
were computed from different quantities. Both are fixed by the same one-line move, but they were
two faults, not one.

**The verdict's own tolerance had to be corrected on first use**, which is worth recording because
it is the same class of mistake the ledger keeps finding. The first version used a per-mark
absolute floor plus a relative test on the declared flux — and immediately flagged two phases.
It was right to: `disturbance` and `patch fuse/term` declare *no* carbon, so their whole allowance
was `3 × 1e-12`, while the structural operators permute, merge and renormalise a ~25 kgC/m² store
and cost ~1.3e-11 in arithmetic doing it. Round-off scales with the **store**, not with the number
of checks. The tolerance now carries a store-proportional round-off allowance of 1e-11 per mark.

That is *not* the store-relative tolerance `budget_check`'s header warns against: 1e-11 is a
round-off bound five orders tighter than the 1e-6 that let a sustained 1 W/m² leak hide, and the
smallest **real** term this ledger ever caught — the disturbance seed bank, 2.96e-6 — is still four
orders above it.

**Verification.** 42/42 on ifx and nvfortran. Not byte-identical, and cannot be — mortality water now
reaches the soil and recruits are born at a different temperature. Site integrals move by ≤4e-8
relative (`veg_carbon_site`, `gpp_site`, `nee_site`, `agb_site`, `nplant_site` all unchanged to six
digits); soil carbon −0.004 %, Rh +0.04 %. The per-cohort `dmax_psi_leaf` distribution moves by a
median of 6e-7 with a handful of cohorts — the recruits whose birth temperature changed — moving up
to 0.068 across a range of 0.35, which is the intended change and confined to them.

