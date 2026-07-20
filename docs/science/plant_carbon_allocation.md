# Plant carbon allocation

MEDS turns each cohort's daily carbon budget into tissue growth with a single **stateless, elemental**
kernel, `plant_carbon_allocation` (`meds_plant_carbon_allocation`). It is the mechanistic replacement for
the phenomenological growth engine, covering the science of ED2's `growth_balive` + `structural_growth`
**unified into one daily step**, and it follows FATES PARTEH Hypothesis-1 ("Allometrically Guided, Carbon
Only"): every pool has an allometric target, the daily net gain fills toward the targets in **priority
order**, and the residual advances stature (wood). All quantities are carbon $`[\mathrm{kgC\,plant^{-1}}]`$
over one slow step; every biomass↔carbon conversion is folded into the PFT traits once at initialization,
so the kernel never converts. See `docs/dev_plans/MEDS_PLANT_CARBON_ALLOCATION_REFACTOR_DESIGN.md`.

## 1. The daily carbon budget

The step's gross primary production $`G`$ and maintenance respiration $`R_m`$ (both accumulated by the fast
biophysics loop, or a stub when it is off) give the **net carbon** available:

```math
C_{\mathrm{net}} = G - R_m \qquad(1)
```

which may be negative. Growth (construction) respiration is **not** subtracted here — it is charged inside
the ladder, on realized growth only (§3). The kernel also receives the current storage (nonstructural
carbon) $`S`$, this step's shed carbon $`\ell_{\mathrm{leaf}},\ell_{\mathrm{root}}`$ (the phenology→carbon
bridge, see the phenology doc §5), and the four allometric demands.

## 2. PARTEH-H1 allocation ladder

Targets come from allometry: leaf $`L^{*}=`$ `size2leaf_carbon`, fine root $`q\,L^{*}`$ (with $q$ the
root:leaf ratio), storage $`c_s\,L^{*}`$. A **demand** is the deficit `target − pool`; the leaf and
fine-root demands arrive already **flush-capped** ($`\Delta_{\mathrm{fl}} = r_{\mathrm{fl}} L^{*} dt`$),
so a dormant canopy ($`r_{\mathrm{fl}}=0`$) presents zero growth demand. When $`C_{\mathrm{net}}\ge0`$ the
kernel funds, in order:

1. **leaf + fine-root growth** toward the flush-capped demand — funded from NPP, then storage;
2. **storage refill** toward $`c_s L^{*}`$ — NPP only, no construction cost;
3. **reproduction** — a fraction $`f_r`$ of the post-storage residual (zero below the maturity height);
4. **wood** — the residual sink: everything left becomes structural growth.

Wood needs no explicit demand (it simply absorbs the remainder), and any cohort below the maturity height
sets $`f_r=0`$.

When $`C_{\mathrm{net}}<0`$ (e.g. a leafless canopy at bud-break, where GPP≈0 but stem/root maintenance
still runs), the maintenance debt is paid **from storage first**, and then leaf/fine-root growth may still
draw the **remaining** reserves — so spring leaf-out is storage-funded, not deadlocked. Only a plant whose
storage cannot even cover maintenance is flagged `starving`, with the shortfall reported as `deficit` for
the stateful updater to resolve by destroying tissue (this pure kernel never mutates a pool). The full
priority order is therefore **maintenance debt → leaf/fine-root growth (NPP then storage) → storage refill
(NPP only) → reproduction → wood**.

## 3. Growth respiration on realized growth

Building one unit of **growth tissue** (leaf, fine root, wood, or reproduction) consumes $`(1+g)`$ carbon,
where $g$ = `growth_resp_factor`: the fraction $g$ is respired as **growth (construction) respiration**.
Storage refill is a 1:1 carbon transfer (nonstructural sugar has no construction cost), so it is exempt.
Charging $`(1+g)`$ *inside* the funding step resolves — exactly and without iteration — the circularity
that growth respiration reduces the carbon available for growth, which changes the growth. If the pools
built are $`a_{\mathrm{leaf}},a_{\mathrm{root}},a_{\mathrm{wood}},a_{\mathrm{repro}}`$, then

```math
R_g = g\,(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}}) \qquad(2)
```

This corrects the pre-refactor engine, which charged growth respiration on the whole pre-allocation balance
(including the part bound for storage).

## 4. Leaf display as emergent replaceability

There is one **shed rate** and one **flush rate** per tissue (phenology doc §1′), not a
replaceable/non-replaceable split. Whether a shed is replaced is **emergent**: the shed opens a leaf
deficit, and the flush-capped growth step (P1) refills it *iff* flushing is active. An evergreen holds a
full canopy because its baseline-floor shed is continuously refilled; a deciduous canopy in dormancy
($`r_{\mathrm{fl}}\to0`$ ⇒ flush cap $`\to0`$) is not refilled and drives to bare. The shed carbon decays
the current pool and snaps to bare only during dormancy (phenology doc eq 8), so a partially-built canopy
is never over-shed.

## 5. Carbon closure and litter

The kernel is **growth-only** — it returns the per-pool growth $`a_{\mathrm{leaf}},a_{\mathrm{root}},
a_{\mathrm{wood}},a_{\mathrm{repro}}\ge0`$ and the net storage change $`\mathrm{npp}_{\mathrm{store}}=`$
refill − drawdown. The **turnover/shed is applied upstream** by the driver
(`update_biomass_turnover`): the pools handed to the kernel are already net of this step's shed, and the
driver forms the net leaf/root change $`\mathrm{npp}_{\mathrm{leaf}}=a_{\mathrm{leaf}}-\ell_{\mathrm{leaf}}`$
(and likewise for root). The **growth-side** budget the kernel closes on every call is

```math
\big(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}} + \mathrm{npp}_{\mathrm{store}}\big) \;-\; \mathrm{deficit} \;=\; (G - R_m) \;-\; R_g \qquad(3)
```

with `deficit` the unpaid maintenance on the starving branch (it *adds back* — carbon the plant owes but
has not yet removed from any pool). Adding the driver's shed, the full **plant-pool change = GPP −
(maintenance + growth respiration) − litter**, where the litter $`\ell_{\mathrm{leaf}}+\ell_{\mathrm{root}}`$
feeds the (deferred) demography→litter→$`R_h`$ biogeochemistry seam. Because growth respiration is charged
on realized growth and the shed decays the current pool, this path **closes in carbon but is not
bit-identical** to the pre-refactor engine.

## Parameters (per-PFT traits consumed)

| Symbol | Config key | Meaning |
|---|---|---|
| $g$ | `growth_resp_factor` | construction-cost fraction charged on growth |
| $`c_s`$ | `storage_cushion` | storage target as a multiple of the leaf target |
| $q$ | `root_to_leaf_ratio` | fine-root : leaf target ratio |
| $`f_r`$ | `reproduction_investment_fraction` | fraction of the residual → reproduction (above maturity) |
| $`L^{*}`$ inputs | `sla`, `hgt_max`, … | full-canopy leaf carbon via `size2leaf_carbon` |

The turnover / flush / shed **rates** and the snap-to-bare fraction $`e_{\min}`$ are phenology traits
(see `plant_phenology.md` §"Parameters"); the allocation kernel consumes their *carbon amounts*.

## Interface with other modules

`meds_plant_carbon_allocation` is a **pure kernel library**: it `use`s only `meds_kinds` — no `site_t`,
no config, no PFT table. Everything is passed as plain scalars, so the **driver**
(`meds_vegetation_dynamics.carbon_growth`) is the single place that assembles the inputs and disposes of
the outputs, calling the kernel once per cohort as an `elemental` sweep over the cohort Structure-of-Arrays.
This is the *only* plant-flux call on the slow carbon path.

**Inputs the driver gathers per cohort:**

| input | source |
|---|---|
| `gpp`, `resp_maint` | fast-loop accumulators on the cohort SoA (`gpp_accum`, `*_resp_accum`) when `fast_biophysics_on`; otherwise the `gpp_ref·leaf_area` stub |
| `leaf_demand`, `fineroot_demand` | allometric deficit (`meds_allometry.size2leaf_carbon`, using the **plastic** `cohort%sla`) **flush-capped** by the phenology flush rate |
| `storage_demand`, `storage`, `repro_frac` | cohort SoA (`nonstructural_carbon`) + PFT traits (`storage_cushion`, `reproduction_investment_fraction`) |
| `growth_resp_frac` | PFT trait `growth_resp_factor` |
| this step's **shed** (litter) | `update_biomass_turnover` converts the phenology shed rate (`meds_phenology.turnover_shed_rates` / `pheno_drives_to_rates`, floored at the baseline leaf turnover `1/llspan`) into a carbon amount and pre-subtracts it from the pools |

**Outputs the driver disposes of:**

| output | destination |
|---|---|
| per-pool growth `growth_{leaf,fineroot,wood,repro}`, `npp_store` | the driver forms the net per-pool change `growth − shed`, assembles the `carbon_flux_block`, and hands it to `compute_slow_derivatives` → the `wood_carbon → dbh` flip (`meds_allometry.carbon_to_structure`) → the core engine's `update_cohort_states` |
| `growth_resp` | autotrophic-respiration accounting (with maintenance resp) |
| `deficit`, `starving` | flags for the stateful updater to resolve by destroying tissue (not yet acted on) |
| leaf + fine-root **litter** (`shed`) | kept in the driver for the (deferred) demography→litter→$`R_h`$ biogeochemistry seam |

`growth_respiration` is **re-exported** through `meds_plant_interface`; the whole kernel is orthogonal to
the core engine (`demography ⊥ plant`), which only ever *applies* the tendency arrays the driver backs out.

## Code map

| Concept | Routine |
|---|---|
| daily GROWTH allocation (master, elemental) | `meds_plant_carbon_allocation`: `plant_carbon_allocation` (+ private `fill_carbon_demand`, the $`(1+g)`$ funder) |
| growth (construction) respiration | `meds_plant_carbon_allocation`: `growth_respiration(npp_growth, g)` |
| baseline turnover as a shed rate | `meds_phenology`: `turnover_shed_rates`, `pheno_drives_to_rates` |
| turnover → shed carbon amount + snap-to-bare | `meds_vegetation_dynamics`: `update_biomass_turnover` (+ private `leaf_shed_amount`) |
| per-cohort orchestration (slow loop) | `meds_vegetation_dynamics`: `carbon_growth` (rates → shed-first → flush-capped demands → one elemental call → net npp + litter) |
| geometry flip (wood_carbon → dbh) | `meds_vegetation_dynamics`: `compute_slow_derivatives`; `meds_allometry`: `carbon_to_structure` |
