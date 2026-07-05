# MEDS Plant Carbon Dynamics — Design

`meds_plant_carbon_dynamics.f90` — the mechanistic replacement for the phenomenological growth
engine. A per-cohort **carbon budget + allocation** module that turns assimilation/respiration into
tissue growth. Covers the science of ED2's `growth_balive.f90` + `structural_growth.f90`, but
**unified into one daily allocation** (no separate monthly structural step — see §0). Follows FATES
**PARTEH Hypothesis-1** ("Allometrically Guided, Carbon Only"): every pool has an allometric target,
allocation fills toward targets by priority, and the residual advances stature. NOT CLM's fixed-ratio
partitioning (CLM has no prognostic size axis and is N-gated).

Part III of the plant-ecophysiology roadmap (after leaf gas exchange, hydraulics, phenology,
respiration; see `MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md`). Reuses `growth_respiration`
(`meds_plant_respiration`) and actuates the ON/OFF/DORMANT signal from `meds_phenology`.

---

## 0. Locked scope (user decisions 2026-07-04, refined)

1. **All state is CARBON.** Pool names are carbon-explicit (`leaf_carbon`, `wood_carbon`,
   `nonstructural_carbon`, …); every biomass↔carbon conversion (C2B, wood density → carbon density,
   SLA in carbon units) is done ONCE at parameter initialization. The kernels never convert.
2. **Carbon-prognostic, wood-anchored.** Four prognostic carbon pools:
   `leaf_carbon, fineroot_carbon, wood_carbon, nonstructural_carbon`. `wood_carbon` (total woody =
   sapwood + heartwood) is the size anchor: `dbh = wood_to_dbh(wood_carbon)`, derived (flips today's
   size-prognostic design). `sapwood_carbon` and `heartwood_carbon` are **diagnostic** (§2).
3. **Sapwood by Huber value, not the qsw pipe model.** `sapwood_carbon` is slaved to leaf area via
   `huber_value`; `heartwood_carbon = wood_carbon − sapwood_carbon`. One trait feeds both the carbon
   budget and the hydraulics boundary condition (`sap_area`).
4. **Woody growth is DAILY**, expressed as `wood_carbon_demand` and allocated as the lowest-priority
   residual sink — no hardcoded monthly cadence, no `cb` monthly buffers. If no carbon remains after
   higher priorities, `npp_wood = 0`. (This is the intended MEDS end-state, and more FATES-faithful.)
5. **No `cb` carbon-balance buffer.** ED2's relative-carbon-balance is empirical/non-measurable.
   Carbon starvation propagates through the physical chain instead: negative NPP → storage drawdown →
   wood growth → 0 → dbh growth-rate → 0 → the existing Camac low-growth mortality hazard rises.
6. **Reproduction + mortality stay phenomenological this round** (`meds_vital_rates`) — no
   double-count. Mechanistic seed allocation is a later, optional PR.
7. **No nutrients** (validated: PARTEH-H1 is pure carbon).
8. **This PR = the two PURE functions only.** The stateful state application (apply npp to the SoA
   pools, re-derive `dbh` from `wood_carbon`) lands as a demography DATA seam next PR; the plant kernel
   is invoked one layer up, in the **driver** (`demography` never links `plant`) — see §5–§6.

### The causal model

```
TODAY (size-prognostic):   dbh --allometry--> {height, agb, leaf_area}
NEW  (carbon-prognostic):  NPP --allocation--> {leaf,fineroot,wood,nonstructural}_carbon
                           wood_carbon --wood_to_dbh--> dbh --> {height, agb, leaf_area}
                           sapwood_carbon = huber(leaf_area,height)   ; heartwood = wood - sapwood
```

---

## 1. Prognostic pools & diagnostics (SoA — next PR, stated here for design closure)

**Prognostic carbon pools** (add to `cohort_block`, `real(wp) :: (:)`, [kgC/plant]):

| Field | Role |
|---|---|
| `leaf_carbon` | living leaf; own turnover; flushed from storage on phenology ON |
| `fineroot_carbon` | living fine root; target `= root_to_leaf_ratio · leaf_target` |
| `wood_carbon` | **total woody structural carbon — the size anchor**; `dbh = wood_to_dbh(wood_carbon)` |
| `nonstructural_carbon` | storage / TNC; target `= storage_cushion · leaf_target` |

**Diagnostic** (cached, re-derived in `set_cohort_size`):
```
dbh             = wood_to_dbh(wood_carbon, wood_carbon_density, hgt_max)
height          = dbh_to_height(dbh, hgt_max)
leaf_area       = leaf_carbon * sla
sapwood_carbon  = huber_value * leaf_area * height * wood_carbon_density   ! Huber
heartwood_carbon= wood_carbon - sapwood_carbon
sap_area        = huber_value * leaf_area        ! hydraulics BC (k_plant / sapflow)
basal_area      = pio4 * dbh^2
agb             = leaf_carbon + aboveground_frac * wood_carbon
```
No `cb` fields (dropped). Every new field threads through the lockstep machinery
(`cohort_alloc`/`site_free`/`cohort_ensure_capacity`/`move_alloc_block`/`cohort_reorder`/
`cohort_compact`/`copy_cohort_slot`/`set_cohort_size`) — the "one place knows every field" invariant.

---

## 2. Allometry additions — `src/allometry/meds_allometry.f90` (next PR)

All `elemental pure`, in `x = dbh²·height`, coefficients installed by `set_allometry`.

**Leaf target (un-fold SLA, behavior-preserving).** Add `sla` as a real trait; set
`leaf_c_b1 = lai_b1 / sla` at load, so `size2leaf_carbon(dbh,h) = leaf_c_b1 · x^lai_b2`. LAI is
unchanged: `nplant · leaf_carbon · sla = nplant · lai_b1 · x^lai_b2`.

**Wood target + analytic inverse.** `wood_carbon` is the prognostic; define a single power law so the
inverse is analytic (like today's `agb_to_dbh`):
```
size2wood_carbon(dbh,h,rho_c) = wood_c_b1 * rho_c^agb_c2 * x^wood_c_b2
wood_to_dbh(wood_carbon, rho_c, hgt_max)      ! analytic inverse + capped-height branch
```
`wood_c_b1, wood_c_b2` **calibrated once at load** so the diagnostic `agb = leaf_carbon +
aboveground_frac·wood_carbon` reproduces today's Chave `dbh_to_agb` curve within tolerance (keeps the
AGB currency for validation and spin-up continuity).

**Sapwood (Huber, diagnostic):** `sapwood_carbon = huber_value · (leaf_carbon·sla) · height ·
wood_carbon_density`. No new allometry function — it's a product of existing quantities.

---

## 3. New PFT traits — `src/shared/meds_pft_params.f90` + `meds_config_pft.toml` (next PR)

Carbon-explicit, conversions folded in at init. Each needs a TOML key + `write_pft_params_csv` column.

| Trait | Units | Meaning |
|---|---|---|
| `sla` | m²/kgC | specific leaf area (un-folded from `lai_b1`) |
| `root_to_leaf_ratio` | – | fine-root : leaf target ratio (ED2 `q`) |
| `huber_value` | m² sap / m² leaf | sapwood-area : leaf-area (ties sapwood + hydraulics) |
| `aboveground_frac` | – | aboveground fraction of woody carbon (ED2 `agf_bs`) |
| `storage_cushion` | – | `nonstructural_carbon` target as a multiple of leaf target |
| `growth_resp_factor` | – | construction cost (already used by `growth_respiration`) |
| `leaf_turnover_rate` | 1/yr | baseline leaf turnover (constant now; §4a extensible) |
| `fineroot_turnover_rate` | 1/yr | baseline fine-root turnover |
| `wood_carbon_density` | kgC/m³ | for the Huber sapwood-carbon conversion |
| `evergreen` | flag | selects the evergreen maintenance temperature factor |

---

## 4. The module — `src/plant/meds_plant_carbon_dynamics.f90`

Links `meds_shared` only; types in `meds_plant_types`; seam re-exported via `meds_plant_interface`
(same pattern as respiration/hydraulics/phenology). Two `pure` public procedures. Stateless — nothing
persists between calls; the caller owns the pools.

### 4a. `tissue_turnover_rates(env, params, rates)`

Returns per-tissue turnover **rates** [1/yr]. Trivial now (constant from traits), but its own function
so it can later depend on light / tropical leaf phenology (the light-driven leaf-turnover you flagged)
without touching allocation.
```
rates%leaf     = params%leaf_turnover_rate      ! * (future light/phenology modifier)
rates%fineroot = params%fineroot_turnover_rate
```
Evergreen maintenance temperature factor `f_T = 1/(1+exp(0.4·(278.15 − tissue_temp)))` applied here
(or passed through) when `params%evergreen`.

### 4b. `plant_carbon_allocation(gains, losses, demands, flags, npp)`

The allocation core: **gains, losses, demands, flags → net npp per pool.** PARTEH priority ladder,
daily, wood as residual sink. Pure arithmetic (allometry/targets computed upstream and passed as
`demands`), so it stays GPU/SIMD-friendly.

**Inputs (all per plant, [kgC/plant/day] unless noted):**
- `gains`: `net_carbon` = NPP = GPP − R_leaf − R_wood − R_root − R_growth (assembled upstream;
  `R_growth` via `growth_respiration`), and the current `nonstructural_carbon` available for drawdown.
- `losses`: turnover fluxes `turnover_leaf = rates%leaf·leaf_carbon·dt`, `turnover_fineroot = …`
  (from §4a; also the litter flux the caller routes to biogeochem).
- `demands`: allometric deficits (target − current), computed upstream from §2 targets scaled by
  phenology `leaf_on` and drought `elongf`: `leaf_deficit`, `fineroot_deficit`, `storage_deficit`,
  and `wood_demand` (default = residual/unbounded; optionally a max structural-growth rate cap).
- `flags`: `phenology_status` (PHEN_ON/OFF/DORMANT).

**Algorithm (priority ladder):**
```
available = net_carbon + nonstructural_carbon
1. REPLACE TURNOVER   toward leaf/fineroot (highest priority; keeps foliage/roots whole)
2. FILL LEAF & FINEROOT allometric deficits  (flush from storage if net_carbon short & PHEN_ON)
3. FILL NONSTRUCTURAL toward storage target   (storage_deficit)
4. GROW WOOD          = min(remaining_available, wood_demand)   ! residual sink; 0 if nothing left
   remainder (if any) -> nonstructural_carbon
NEGATIVE net_carbon (starvation): draw from nonstructural_carbon first (PHEN_ON/DORMANT) or shed
   leaf/fineroot first (PHEN_OFF); wood is never consumed; set flags%starving if storage exhausted.
```
**Outputs — net carbon flux per pool [kgC/plant/day], signed:**
`npp_leaf, npp_fineroot, npp_wood, npp_nonstructural`
where `npp_i = allocated_i − turnover_i` (net pool change). Carbon closes:
`Σ npp_i = net_carbon − (turnover_leaf + turnover_fineroot)`, and the turnover fluxes are the litter
input to biogeochem. `npp_wood ≥ 0` and `= 0` on days with no surplus (satisfies decision #4).

### Types (add to `meds_plant_types`, carbon-explicit)
```fortran
type :: carbon_env_t          ! per-plant state + drivers
   real(wp) :: net_carbon                 ! NPP after growth resp [kgC/plant/day]
   real(wp) :: leaf_carbon, fineroot_carbon, wood_carbon, nonstructural_carbon
   real(wp) :: dbh, height, leaf_area
   real(wp) :: leaf_on, elongf, tissue_temp
   integer(ik) :: phenology_status
end type
type :: carbon_params_t       ! per-PFT (flattened by the seam)
   real(wp) :: sla, root_to_leaf_ratio, huber_value, aboveground_frac, storage_cushion
   real(wp) :: leaf_turnover_rate, fineroot_turnover_rate, wood_carbon_density, growth_resp_factor
   logical  :: evergreen
end type
type :: turnover_rates_t
   real(wp) :: leaf, fineroot        ! [1/yr]
end type
type :: carbon_npp_t          ! outputs
   real(wp) :: npp_leaf, npp_fineroot, npp_wood, npp_nonstructural   ! net, signed
   real(wp) :: turnover_leaf, turnover_fineroot                       ! litter -> biogeochem
   logical  :: starving
end type
```

---

## 5. State application — `apply_carbon_npp`, `src/demography/` (NEXT PR)

Not in this PR. A **data-array seam** — the mirror of `update_demography(growth[], …)`. It takes the
per-cohort `npp` arrays as PLAIN DATA and applies them to the SoA daily: `pool += npp_pool`, re-derive
`dbh = wood_to_dbh(wood_carbon)` + `set_cohort_size`, route the turnover fluxes to litter, and feed the
dbh increment into the existing `growth_avg` → Camac mortality. Because woody growth is daily (§0.4)
there is no monthly cadence — this replaces `growth_step` one-for-one, and hosts the OpenMP-`target`
kernel.

**Crucially it does NOT call the plant kernel.** `demography` links `state + allometry` only (no plant
edge, verified in `CMakeLists.txt`), so it never *computes* `npp` — it only *applies* it. This preserves
demography's self-containment (standalone spin-up, Python-callable, the "two halves" data seam). The
`npp` is computed one layer up, in the driver (§6).

---

## 6. Orchestration & seam — the driver layer (NEXT PR)

The plant kernel and the demography engine meet in the **driver/stepper layer** (`src/driver`, compiled
into `meds_aux`) — the only layer above BOTH libraries. `demography ⊥ plant` is preserved: neither calls
the other; the driver weaves them. Demography must never gain a `→ plant` edge (it would break the
standalone spin-up, the data-array seam, and the Python-callability of the engine).

```
  meds_aux (driver)  ── links ──►  meds_plant  +  meds_demography  (+ meds_allometry, meds_biophys)
     1. demands   = carbon targets − current pools          (allometry; a demography/allometry helper)
     2. net_carbon = GPP − R_leaf − R_wood − R_root − R_growth   (R_growth = growth_respiration(…))
     3. npp[]     = plant carbon seam(env, cfg, ipft, demands)   (meds_plant_interface)
     4. apply_carbon_npp(site, npp[], dt)                        (meds_demography — DATA in, SoA mutated)
```

A thin coupling module (e.g. `src/driver/meds_carbon_coupling.f90`, part of `meds_aux`) owns steps 1–4;
`meds_aux`'s link line gains `meds_plant` (it links only `meds_demography` today). The DAG stays acyclic:
`aux → {plant, demography, allometry}`, with `plant ⊥ demography`.

**The `meds_plant_interface` seam (step 3):** add a cfg-flattening `carbon_allocation(env, cfg, ipft,
npp)` wrapper mirroring `leaf_gas_exchange` (which flattens `cfg%pft` into `leaf_photo_params_t`). Its
CALLER is the driver, never demography. Note: `meds_plant_interface` links `shared` only, so it flattens
`cfg%pft` traits but CANNOT compute the allometric `demands` (no allometry link) — those come from the
driver side (step 1), which is exactly why `plant_carbon_allocation` takes `demand` as an INPUT (§4b).

Driver-assembled inputs:
- `net_carbon`: `GPP − R_leaf − R_wood − R_root`, then `R_growth = growth_respiration(that, factor)`,
  `net_carbon = that − R_growth`. `growth_respiration` stays the single owner of the `Rg` formula.
- GPP is stubbed (constant/prescribed) until canopy RT → per-cohort absorbed PAR + leaf energy balance
  land, so the allocation engine is testable ahead of the biophysics chain. **Caveat:** with stubbed
  GPP the carbon path does not reproduce light-competition growth spread — do not benchmark against a
  spun-up phenomenological run until GPP is real.
- `elongf` (drought down-regulation of targets) comes from `meds_plant_hydraulics` (`psi_leaf`) /
  water stress; `leaf_on` from the deciduous fraction; `phenology_status` from `meds_phenology`.

---

## 7. Fusion/fission conservation — `meds_demography_structure.f90` (next PR)

Today conserves AGB only. With pools: conserve every prognostic pool
(`leaf_carbon, fineroot_carbon, wood_carbon, nonstructural_carbon`) as the nplant-weighted mean,
re-anchor `dbh = wood_to_dbh(wood_carbon)`, extend the 1% guard per pool. Sapwood/heartwood need no
conservation (diagnostic).

---

## 8. Mortality (unchanged this round)

Camac growth-driven hazard, keyed on `growth_avg` of the dbh rate — which now comes from daily wood
growth → daily dbh increment. No `cb`. Carbon starvation manifests as reduced dbh growth → higher
low-growth hazard. Light competition acts through GPP once RT is wired (transiently absent under stub).

---

## 9. Invariants & double-count guards

- **Single source of truth:** `wood_carbon` prognostic; `dbh`, `agb`, `sapwood_carbon`,
  `heartwood_carbon`, `leaf_area` all diagnostic. Never integrate `dbh`.
- **LAI one way:** `leaf_area = leaf_carbon · sla`; retire the `lai_b1` LAI path once `leaf_carbon`
  exists.
- **Reproduction/mortality:** phenomenological only this round; kernel reserves no seed carbon and no
  `cb`.
- **Growth respiration:** kernel-adjacent caller **calls** `growth_respiration`; single owner.
- **Carbon closure assert:** `Σ npp_i == net_carbon − Σ turnover_i` within tolerance, per cohort/step.
- **All-carbon:** no C2B or density conversion inside the kernels (done at param-init).

---

## 10. Testing (this PR)

Unit-test the two pure functions like the other plant kernels:
(a) on-allometry cohort + constant GPP → transfers → 0 as pools hit targets, `npp_wood` steady > 0;
(b) carbon closes (`Σ npp = net_carbon − turnover`); (c) negative-GPP episode → storage drawdown then
`starving`, `npp_wood == 0`; (d) phenology OFF→ON edge → leaf flush funded from storage; (e) Huber
round-trip `sapwood_carbon`/`sap_area` sane. Build the nvfortran multicore back end too (a green ifx
run is not sufficient — CLAUDE.md portability trap).

---

## 11. Phased PRs

1. **THIS PR — `meds_plant_carbon_dynamics.f90`:** the two pure functions (`tissue_turnover_rates`,
   `plant_carbon_allocation`) + types in `meds_plant_types` + seam in `meds_plant_interface` + unit
   tests. No SoA change, no allometry change yet (targets/demands supplied by the test harness).
2. **Allometry + traits (behavior-preserving):** un-fold `sla`; add `size2leaf_carbon`,
   `size2wood_carbon`, `wood_to_dbh`; add the PFT traits; calibrate `wood_c_b1/b2`. Assert `agb`
   reproduced.
3. **Pools in state:** add the 4 pools to the SoA + all lockstep routines; initialize on-allometry
   from `dbh`; extend fusion/fission + IO with per-pool conservation.
4. **State application + orchestration:** a demography DATA seam `apply_carbon_npp(site, npp[], dt)`
   that mutates the SoA (no plant edge) + a driver-layer coupling module
   (`src/driver/meds_carbon_coupling.f90`, links plant + demography) that computes demands + net_carbon,
   calls the plant `carbon_allocation` seam, and hands `npp` to demography; add `meds_plant` to
   `meds_aux`'s link line; drive with stub GPP; validate. Demography never links plant.
5. **(Later, optional) mechanistic reproduction:** `nonstructural_carbon → seed`, disabling the
   phenomenological reproduction array.

Deferred (documented): monotonic-heartwood prognostic (deciduous refinement); bark + fire thermal
protection; AG/BG split; trait/turnover plasticity beyond constant rates; nutrient (N/P) limitation.
