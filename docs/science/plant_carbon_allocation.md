# Plant carbon allocation

MEDS turns each cohort's daily carbon budget into tissue growth with a single **stateless, elemental**
kernel, `plant_carbon_allocation` (`meds_plant_carbon_allocation`). It covers the science of ED2's
`growth_balive` + `structural_growth` **unified into one daily step**, and it follows FATES PARTEH
Hypothesis-1 ("Allometrically Guided, Carbon Only"): every pool has an allometric target, the daily net
gain fills toward the targets in **priority order**, and the residual advances stature (wood). All
quantities are carbon $`[\mathrm{kgC\,plant^{-1}}]`$ over one slow step; every biomass↔carbon conversion is
folded into the PFT traits once at initialization, so the kernel never converts.

## 1. The daily carbon budget

The step's gross primary production $`G`$ and maintenance respiration $`R_m`$ (both accumulated by the fast
biophysics loop, or a stub when it is off) give the **net carbon** available:

```math
C_{\mathrm{net}} = G - R_m \qquad(1)
```

which may be negative. Growth (construction) respiration is **not** subtracted here — it is charged inside
the ladder, on realized growth only (§3). The kernel also receives the current storage (nonstructural
carbon) $`C_{store}`$, three demands (leaf and fine root, each already net of this step's shed carbon,
which the driver subtracts beforehand: `plant_phenology.md` §5 and §5 below; and the deficit toward the
storage target) and the reproduction fraction `repro_frac`; wood takes the residual.

### Storage maintenance (charged first)

The non-structural pool is charged a **fractional turnover** each step, before the ladder below runs:

```math
M_s = C_{store} \cdot \min\!\bigl(1,\ k_s\,\Delta t\bigr) \qquad(2)
```

with $`k_s`$ = `storage_turnover_rate` [yr⁻¹] per PFT. It is ED2's `growth_balive.f90` form, and like
ED2 it carries **no temperature dependence** — ED2 sets its `maintenance_temp_dep` to 1 here and
leaves the temperature form commented out as "experimental and arbitrary", so adopting one would go
beyond the reference rather than follow it.

Without it, storage would be the one live carbon pool that costs nothing to hold: a cohort could
carry an arbitrarily large reserve for free.

**Default $`k_s = 0`$** (no storage maintenance; whether to turn it on by default is open, see
`docs/ROADMAP.md`). ED2's own values are temperate broadleaf **0.6243**, temperate grass and conifer
**0**, tropical non-grass **1/6**, tropical grass **1/3** — so zero is a legitimate ED2 setting, not an
absence of physics. The charge
is substantial where it is turned on: at ED2's temperate-broadleaf rate an Ithaca run loses 35 % of
GPP and 43 % of AGB over five years from cold start, because the drain compounds through stand
development.

$`M_s`$ is **netted into the step's storage tendency**, not written to the pool directly. That is the
standing rule (the driver computes, the engine applies) and it is also what keeps the carbon ledger
closed: the phase that declares the CO₂ efflux is the phase where the pool has to drop. The allocator
still sees the post-maintenance reserve, so maintenance is paid before growth is funded from it —
ED2's ordering.

### Leaf and fine-root resorption

A fraction $`f_r`$ = `retained_carbon_fraction` of **all** leaf and fine-root loss returns to the
non-structural pool instead of entering litter. With $B$ the background turnover and $S$ the
senescence of `leaf_turnover_step` (`plant_phenology.md` §5, eq. (9)), and $R$ the fine-root turnover:

```math
\Delta C_{leaf} = -(B + S), \qquad
\Delta C_{root} = -R, \qquad
\Delta C_{store} \mathrel{+}= f_r (B + S + R), \qquad
\text{litter} = (1-f_r)(B + S + R) \qquad(3)
```

**Each pool loses the full loss either way.** Crediting storage while removing only the litter
share would *create* carbon. The split is exact by construction.

Background turnover is resorbed like senescence: `leaf_turnover_rate_toc` and `fineroot_turnover_rate`
are tissue loss rates, not litterfall, so the litter is the unresorbed share of the loss.

**Default $`f_r = 0.3`$.** Most measured resorption is of N and P rather than C, so the carbon
fraction is uncertain and is a per-PFT key; examples 03–05 set it to 0, as calibrated. Resorption moves
carbon from the litter path to the plant reserve, and the reserve then funds growth, so a larger
$`f_r`$ raises storage, GPP and LAI and lowers soil carbon.

## 2. PARTEH-H1 allocation ladder

Targets come from allometry: leaf $`L^{*}=`$ `size2leaf_carbon`, fine root $`q\,L^{*}`$ (with $q$ the
root:leaf ratio), storage $`c_s\,(1+q)\,L^{*}`$, the carbon of $`c_s`$ full flushes. A **demand** is the deficit `target − pool`; the leaf and
fine-root demands arrive already **flush-capped** ($`\Delta_{\mathrm{fl}} = r_{\mathrm{fl}} L^{*} dt`$),
so a dormant canopy ($`r_{\mathrm{fl}}=0`$) presents zero growth demand. When $`C_{\mathrm{net}}\ge0`$ the
kernel funds, in order:

1. **leaf + fine-root growth** toward the flush-capped demand — funded from NPP, then storage;
2. **storage refill** toward $`c_s (1+q) L^{*}`$ — NPP only, no construction cost;
3. **reproduction** — a fraction $`f_{\mathrm{repro}}`$ of the post-storage residual (zero below the maturity height);
4. **wood** — the residual sink: everything left becomes structural growth.

Wood needs no explicit demand (it simply absorbs the remainder), and any cohort below the maturity height
sets $`f_{\mathrm{repro}}=0`$. The maturity height is a share, `reproduction_height_fraction` (default 0.33), of the **local maximum
height** $`H_{\max}`$, the height of a 100 cm tree on the run's height curve, capped at `hgt_max`
(`meds_allometry::local_hmax`). Where the local allometry makes trees short, they mature short: this is
the relative size at onset of maturity of Thomas (1996), whose 37 Malaysian species matured at 0.20–0.75
of their asymptotic height, 0.42 on average. Under Chave et al. (2014)'s curve the stress index E scales
every height by the same factor, so a given share of $`H_{\max}`$ falls at the same diameter at every site.

### Sink limitation of wood growth

Carbon supply is not the only limit on growth: cambial cell division and expansion run at a bounded
rate however much carbon the leaves fix (sink limitation; Körner 2015). Wood growth is therefore capped
by `growth_sink_limitation` (`meds_sink_limitation`), which returns the most wood the step can build.
It has two forms, each a per-PFT trait that is off at zero; with both set the tighter applies:

```math
a_{\mathrm{wood}}^{\max} = \min\!\Big(r_{\max}\,C_{\mathrm{wood}}\,\Delta t,\;\;
W\big(D + \max(g_{\max}D^{c}\,T - \Delta D_T,\,0)\big) - W(D)\Big) \qquad(4)
```

The first is **relative**: $`r_{\max}`$ is the PFT's `max_relative_growth_rate` [yr⁻¹] and
$`C_{\mathrm{wood}}`$ the cohort's wood carbon at the step's start. The second is **absolute**:
$`g_{\max}`$ is the PFT's `max_absolute_growth_rate` [cm yr⁻¹], the most a 1 cm stem's diameter can grow,
$c$ its `max_absolute_growth_exponent` (default 0) with $D$ in cm, and $`W(D)`$ the wood carbon of a stem of
diameter $D$ on the model's allometry (`size2wood_carbon`). The absolute cap is a **rate over the growth
window** $T$ (`growth_memory_days`, 365 days in the shipped config, the window mortality averages growth over), not
over each step: $`\Delta D_T`$ is the diameter the stem has already grown in the window, read from the
cohort's ring buffer of realized growth, and $T$ includes this step (for a cohort younger than the window,
its age). A step the limit binds therefore brings the window's growth to exactly $`g_{\max}D^{c}T`$, and a
stem that grew nothing through a dry season may make up the year's diameter in the wet season. A per-step
cap, $`g_{\max}D^{c}\Delta t`$, would instead throw away every day without spare carbon: at Palo Verde it
holds saplings that had carbon to spare all wet season to about half of the cap's annual rate. The cambium lays down a bounded width of wood, so the absolute form
scales with the stem's surface rather than its mass: per unit wood it is far more generous to a sapling
than to a large tree. At BCI the upper quantiles of census diameter growth at a given light rise about
as $`D^{0.5}`$ (example 05), between a size-independent width ($`c=0`$) and the relative form ($`c\approx1`$). When the residual would build more than that, the kernel builds
$`a_{\mathrm{wood}}^{\max}`$ and turns the rest of the residual into **root exudate**,
$`E = C_{\mathrm{resid}}/(1+g) - a_{\mathrm{wood}}^{\max}`$. Making and exporting exudate costs the same
$`(1+g)`$ per unit as tissue, so growth respiration (eq. (5)) is the same whether or not the limit binds;
the limit only changes where the carbon ends up. The driver puts $E$ in the patch's below-ground labile
litter (`labile_soil`), where soil respiration returns it to the air. A
cohort below the limit is unaffected, so the limit only bites where carbon is plentiful (open-grown
trees), and $`r_{\max}=0`$ (the default) means no limit at all. The form is deliberately simple; a
mechanistic sink (temperature- or water-limited cambial activity) replaces the body of
`growth_sink_limitation` without touching the allocator.

When $`C_{\mathrm{net}}<0`$ (e.g. a leafless canopy at bud-break, where GPP≈0 but stem/root maintenance
still runs), the maintenance debt is paid **from storage first**, and then leaf/fine-root growth may still
draw the **remaining** reserves — so spring leaf-out is storage-funded, not deadlocked. Only a plant whose
storage cannot even cover maintenance is flagged `starving`, with the shortfall reported as `deficit`: the
maintenance respiration no pool paid for. The driver nets it out of the next day's CO₂ release
(`plant_respiration.md` §6); no tissue is destroyed (see `docs/ROADMAP.md`), and this pure kernel never
mutates a pool. The full
priority order is therefore **maintenance debt → leaf/fine-root growth (NPP then storage) → storage refill
(NPP only) → reproduction → wood**.

### A cohort with no leaves and no storage

When $`C_{\mathrm{net}} < 0`$ the maintenance debt is paid from storage first; what storage cannot
cover is the **deficit**. No pool pays it and no tissue is destroyed: the fast loop has already
respired it, and the driver corrects that on the next step (`plant_respiration.md` §6).

One case is different: a cohort whose deficit is larger than **all the carbon in its leaves**. Its
storage is spent, and burning every leaf it has could not pay one step of its upkeep, so nothing is
left to pay for a flush and it cannot recover. It dies. Its mortality for the step
(`decline_to_floor_mortality`) takes its density to the density floor, `[demography].negligible_nplant`,
down the path every death takes: its carbon goes to litter and its water to the ground. The next
monthly cull removes the cohort, so it is gone within a month.

The test needs no threshold: its scale is the cohort's own maintenance respiration. A canopy holds far
more leaf carbon than one step's upkeep, and a cohort with storage pays it, so a leafless cohort living
on its reserves through a dry season or a winter is not affected. Every other cohort dies at the
growth-driven hazard; there is no separate starvation hazard. This case is interim, until a growth- or
storage-based mortality replaces it (`docs/ROADMAP.md`).

## 3. Growth respiration on realized growth

Building one unit of **growth tissue** (leaf, fine root, wood, or reproduction), or exporting one unit of
root exudate (§2), consumes $`(1+g)`$ carbon,
where $g$ = `growth_resp_factor`: the fraction $g$ is respired as **growth (construction) respiration**.
Storage refill is a 1:1 carbon transfer (nonstructural sugar has no construction cost), so it is exempt.
Charging $`(1+g)`$ *inside* the funding step resolves — exactly and without iteration — the circularity
that growth respiration reduces the carbon available for growth, which changes the growth. If the pools
built are $`a_{\mathrm{leaf}},a_{\mathrm{root}},a_{\mathrm{wood}},a_{\mathrm{repro}}`$ and the exudate is $E$, then

```math
R_g = g\,(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}} + E) \qquad(5)
```

A charge on the whole pre-allocation balance would also bill the carbon bound for storage, which builds
nothing (`plant_respiration.md` §5).

## 4. Leaf display as emergent replaceability

The leaf loss (background turnover plus senescence) and the flush cap come from one routine,
`leaf_turnover_step` (`plant_phenology.md` §5), not a replaceable/non-replaceable split. Whether a loss is
replaced is **emergent**: the loss opens a leaf deficit, and the flush-capped growth step (step 1 of the
ladder, §2) refills
it *iff* flushing is active. A flushing canopy holds full cover because its turnover is continuously
refilled; a deciduous canopy in dormancy (flush tendency $`\to0`$ ⇒ flush cap $`\to0`$) is not refilled
and senesces to near bare; an evergreen's senescence stops at its leaf-cover floor. The losses decay
the current pool, so a partially-built canopy is never over-shed.

## 5. Carbon closure and litter

The kernel is **growth-only** — it returns the per-pool growth $`a_{\mathrm{leaf}},a_{\mathrm{root}},
a_{\mathrm{wood}},a_{\mathrm{repro}}\ge0`$ and the net storage change $`\mathrm{npp}_{\mathrm{store}}=`$
refill − drawdown. The **turnover/shed is applied upstream** by the driver
(`cohort_carbon_demand`): the pools handed to the kernel are already net of this step's shed, and the
driver forms the net leaf/root change $`\mathrm{npp}_{\mathrm{leaf}}=a_{\mathrm{leaf}}-\ell_{\mathrm{leaf}}`$
(and likewise for root). The **growth-side** budget the kernel closes on every call is

```math
\big(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}} + \mathrm{npp}_{\mathrm{store}} + E\big) \;-\; \mathrm{deficit} \;=\; (G - R_m) \;-\; R_g \qquad(6)
```

with `deficit` the unpaid maintenance on the starving branch (it *adds back* — carbon the plant owes but
has not yet removed from any pool) and $E$ the root exudate of the sink limit (zero without one). Adding the
driver's terms closes the whole plant: the **plant-pool change = GPP − (maintenance + storage maintenance +
growth respiration) − litter − exudate − reproduction**, where the litter is the unresorbed share of the
shed, $`(1-f_r)(\ell_{\mathrm{leaf}}+\ell_{\mathrm{root}})`$ (eq. (3)). The driver puts it in the patch's
litter, from which the soil-carbon pools and $`R_h`$ follow (`soil_carbon.md`).

## Parameters (per-PFT traits consumed)

| Symbol | Config key | Meaning |
|---|---|---|
| $g$ | `growth_resp_factor` | construction-cost fraction charged on growth |
| $`c_s`$ | `storage_cushion` | storage target as a multiple of the leaf + fine-root target |
| $q$ | `root_to_leaf_ratio` | fine-root : leaf target ratio |
| $`f_{\mathrm{repro}}`$ | `reproduction_investment_fraction` | fraction of the residual → reproduction (above maturity) |
| — | `reproduction_height_fraction` (scalar, optional, default 0.33) | maturity height as a share of the local maximum height |
| $`k_s`$ | `storage_turnover_rate` (optional, default 0 = off) | storage maintenance, eq. (2) [yr⁻¹] |
| $`f_r`$ | `retained_carbon_fraction` (optional, default 0.3) | share of leaf and fine-root loss resorbed to storage, eq. (3) |
| $`r_{\max}`$ | `max_relative_growth_rate` (optional, default 0 = off) | sink limit on wood growth per unit wood carbon [yr⁻¹]; the residual above it is exuded |
| $`g_{\max}`$ | `max_absolute_growth_rate` (optional, default 0 = off) | sink limit on the diameter growth of a 1 cm stem [cm yr⁻¹], a budget over the growth window; the residual above it is exuded |
| $T$ | `growth_memory_days` (`[demography]`) | the window the absolute limit (and mortality's growth average) runs over |
| $c$ | `max_absolute_growth_exponent` (optional, default 0) | how that limit scales with diameter, $`g_{\max}D^{c}`$ |
| $`L^{*}`$ inputs | `sla`, `hgt_max`, … | full-canopy leaf carbon via `size2leaf_carbon` |

The flush and senescence rates, `min_leaf_cover` and `bare_leaf_cover` are phenology traits (see
`plant_phenology.md` §"Parameters"); the allocation kernel consumes their *carbon amounts*.

## Interface with other modules

`meds_plant_carbon_allocation` is a **pure kernel library**: it `use`s only `meds_kinds` — no `polygon_t`,
no config, no PFT table. Everything is passed as plain scalars, so the **driver**
(`meds_vegetation_dynamics.compute_carbon_allocation`) is the single place that assembles the inputs and disposes of
the outputs, calling the kernel once per cohort as an `elemental` sweep over the cohort Structure-of-Arrays.
This is the *only* plant-flux call on the slow carbon path.

**Inputs the driver gathers per cohort:**

| input | source |
|---|---|
| `gpp`, `resp_maint` | fast-loop accumulators on the cohort SoA (`gpp_accum`, `*_resp_accum`) when `fast_biophysics_on`; otherwise the `gpp_ref·leaf_area` stub |
| `leaf_demand`, `fineroot_demand` | allometric deficit (`meds_allometry.size2leaf_carbon`, using the **plastic** `cohort%sla`) **flush-capped** by the phenology flush rate |
| `storage_demand`, `storage`, `repro_frac` | cohort SoA (`nonstructural_carbon`) + PFT traits (`storage_cushion`, `reproduction_investment_fraction`) |
| `growth_resp_frac` | PFT trait `growth_resp_factor` |
| this step's **shed** (litter) | `cohort_carbon_demand` calls `meds_phenology.leaf_turnover_step` (background turnover at the cohort's `leaf_turnover_rate` while flushing and no leaves drop, plus senescence down to `min_leaf_cover`) and turns fine roots over on the same `turnover_gate`, for the carbon amounts it pre-subtracts from the pools |

**Outputs the driver disposes of:**

| output | destination |
|---|---|
| per-pool growth `growth_{leaf,fineroot,wood,repro}`, `npp_store` | the driver forms the net per-pool change `growth − shed`, assembles the `pool_change_block` (`meds_vegetation_dynamics`; the net pool change, not the gross NPP the `CO_NPP_*` diagnostics report), and hands it to `update_cohort_derivatives` → the `wood_carbon → dbh` flip (`meds_allometry.carbon_to_structure`) → the core engine's `update_cohort_states` |
| `growth_resp` | autotrophic-respiration accounting (with maintenance resp) |
| `deficit`, `starving` | `deficit` is netted, with the opposite sign, into the patch's `slow_co2_rate` (`plant_respiration.md` §6); no tissue is destroyed. A deficit larger than the cohort's leaf carbon means it has no leaves and no storage, and it dies (§2) |
| leaf + fine-root **litter** (`shed`) | the resorbed share to storage, the rest to the patch's litter (`add_necromass`); the exudate to `labile_soil` |

`growth_respiration` lives in `meds_plant_carbon_allocation` beside the allocation kernel; the whole kernel is orthogonal to
the core engine (`demography ⊥ plant`), which only ever *applies* the tendency arrays the driver backs out.

## Code map

| Concept | Routine |
|---|---|
| daily GROWTH allocation (master, elemental) | `meds_plant_carbon_allocation`: `plant_carbon_allocation` (+ private `fill_carbon_demand`, the $`(1+g)`$ funder) |
| growth (construction) respiration | `meds_plant_carbon_allocation`: `growth_respiration(npp_growth, g)` |
| sink limit on wood growth | `meds_sink_limitation`: `growth_sink_limitation` (relative and absolute limits, the tighter applies) |
| leaf loss (turnover + senescence) + flush cap | `meds_phenology`: `leaf_turnover_step` |
| turnover-first carbon demand per cohort | `meds_vegetation_dynamics`: `cohort_carbon_demand` |
| per-cohort orchestration (slow loop) | `meds_vegetation_dynamics`: `compute_carbon_allocation` (rates → shed-first → flush-capped demands → one elemental call → net npp + litter) |
| geometry flip (wood_carbon → dbh) | `meds_vegetation_dynamics`: `update_cohort_derivatives`; `meds_allometry`: `carbon_to_structure` |
