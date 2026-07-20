# MEDS plant trait dynamics — DESIGN (review draft, not implemented)

Status: **design-only**. Proposes `src/plant/meds_plant_trait_dynamics.f90` — the home for processes
that CHANGE a plant's traits within its lifetime. First cut = **light-driven plasticity** of the four
leaf traits SLA, Vcmax, Rd, and leaf lifespan, porting ED2's `update_cohort_plastic_trait`
(`../ED2/ED/src/utils/update_derived_utils.f90`). Thermal **acclimation** is a **deferred** sibling (§5).

## 1. Motivation & scope

`meds_pft_params` holds the **static** per-PFT trait table. But leaf traits are not static — they
**acclimate to the light environment**: a shaded leaf is thinner (higher SLA), lower-capacity (lower
Vcmax/Rd), and longer-lived (larger leaf lifespan). MEDS has no home for these *dynamic* trait
processes today. This module is that home, kept separate from the static table (hence NOT
`meds_plant_traits`, which would read as the table).

Scope of the first cut = light plasticity of **SLA, Vcmax25, Rd25, and leaf lifespan (`llspan`),
implemented together** (they co-vary and partly offset — doing one alone biases the response, §3d). All
four become **cohort-level dynamic traits** (§3c). Thermal acclimation is deferred (§5).

## 2. The ED2 process (what we port)

`update_cohort_plastic_trait` scales top-of-canopy (`_toc`) trait values by the **cumulative LAI above**
the cohort (`max_cum_lai`), the light-competition proxy:

```
new_sla    = sla_toc    * exp(kplastic_sla    * cum_lai_above)      (shaded => higher SLA)
new_vcmax  = vcmax_toc  * exp(kplastic_vm0    * cum_lai_above)      (kplastic_vm0 < 0 => shaded => lower Vcmax)
new_rd, new_llspan  scale analogously; exponents clamped to [lnexp_min, lnexp_max].
```

Update is either instantaneous (snap to target) or gradual (a fractional step toward target). ED2 then
re-derives geometry from `bleaf` because its size anchor is leaf biomass — MEDS does **not** (§3d).

## 3. MEDS design

Follows the established split — **stateless kernel in `plant`, orchestration in the driver, state in
`core`** — exactly like phenology and carbon allocation.

### 3a. Kernel (`meds_plant_trait_dynamics`, `plant` library)
Pure/elemental, scalar-in/scalar-out (GPU/SoA-safe):

```fortran
elemental pure subroutine light_plastic_traits(cum_lai_above,                                  &
        sla_toc, vcmax25_toc, rd25_toc, llspan_toc,                                            &
        kplastic_sla, kplastic_vm0, kplastic_rd, kplastic_llspan, lnexp_min, lnexp_max,        &
        sla_target, vcmax25_target, rd25_target, llspan_target)   ! the light-acclimated targets
```
returns the acclimated TARGET traits. A sibling `relax_trait(current, target, tau, dt)` gives the
**gradual** mode (a rate toward the target, governor-style like phenology) so the kernel stays free of an
`is_instant` switch; `tau -> 0` recovers the instantaneous snap.

### 3b. Light input — reuse `overtopping_lai`
MEDS already maintains the cumulative-LAI-above diagnostic per cohort (`update_overtopping_lai`,
`meds_core`) — the direct analogue of ED2's `max_cum_lai`. The driver passes `cohort%overtopping_lai(j)`
into the kernel; no new light machinery is needed.

### 3c. Cohort-level trait STATE + the `p_` naming convention
These four traits become genuine per-cohort STATE, so the cohort block gains **`sla`, `vcmax25`, `rd25`,
`llspan`** — with **no `p_` prefix**. Convention (new, to be documented in `meds_core_state_types`):

- **`p_<trait>`** = a per-cohort cache of a *static* PFT constant, copied at birth and never varying
  (`p_wood_density`, `p_hgt_max`, `p_aboveground_frac`, `p_root_to_leaf_ratio`).
- **bare `<trait>`** = a *dynamic* cohort trait that plasticity/acclimation mutates.

So the existing `p_sla` is **renamed `sla`** (it is now plastic), and `vcmax25/rd25/llspan` join it
bare. `leaf_gas_exchange`'s flattening reads `cohort%vcmax25(j)` / `cohort%rd25(j)` instead of the PFT
table. Each new field must be threaded through the ONE centralized lockstep machinery in
`meds_core_state_types` (`copy_cohort_slot` / `cohort_reorder` / `set_cohort_size` / capacity) and,
crucially, **cohort FUSION** (§3d) — this SoA + fusion work, not the kernel, is the bulk of the effort.

### 3d. Leaf-AREA-conserving allocation (the key modeling choice)
MEDS is **already leaf-area conserving**: the leaf carbon target is
`size2leaf_carbon(dbh, h, sla) = dbh_to_leaf_area(dbh, h) / sla`. The allometry fixes the leaf **AREA**
from size (crown/light-capture); SLA sets how much **carbon** that area costs. So once `sla` is the
cohort's plastic value:

- a shaded cohort (higher SLA) reaches the **same allometric leaf area with LESS leaf carbon** — thinner,
  cheaper leaves, i.e. **smaller leaf biomass**, exactly the intended behavior. The freed carbon flows to
  other pools through the normal allocation ladder.
- the leaf-area diagnostic stays `leaf_area = leaf_carbon · sla`; at target it returns the allometric
  area. No `bl2dbh`/`bl2h` re-derivation (ED2's leaf-biomass anchoring) — MEDS anchors size on
  `wood_carbon`, so SLA plasticity touches leaf display/cost only, never dbh/height. The one feedback is
  SLA → leaf_area → `overtopping_lai` → next step's plasticity (a benign slow loop).

**Fusion** keeps the existing conserved invariants — total leaf/wood **carbon** (AGB) + plant number —
and **leaf-area-weights the four intensive traits** (`weight = cohort leaf_area = leaf_carbon·sla`), the
physically correct combination for per-leaf-area properties. Leaf area itself stays a derived quantity
(`Σ leaf_carbon·sla`), not an independently conserved one. **This weighting is the load-bearing step**: a
new per-cohort trait silently left unaveraged on fusion is precisely the bug class the centralized
machinery exists to prevent.

### 3e. `llspan` drives the baseline leaf turnover
Leaf lifespan and leaf turnover are the same fact twice, so **baseline leaf turnover rate = `1 / llspan`
[1/yr]** rather than an independent PFT number. Concretely: replace the static PFT `leaf_turnover_rate`
with **`leaf_lifespan_toc` [yr]** (top-of-canopy); initialize `cohort%llspan` from it; and have the
phenology turnover floor consume `1/cohort%llspan(j)`. The seam already exists — `turnover_shed_rates` /
`pheno_drives_to_rates` (`meds_phenology`) take the leaf turnover rate as a scalar argument, and the
driver currently passes `pft%leaf_turnover_rate`; it will pass `1/cohort%llspan(j)` instead. So a shaded,
longer-lived canopy automatically turns its leaves over more slowly. Fine-root turnover stays a PFT rate
(no root plasticity in this cut).

### 3f. Orchestration & cadence (driver)
A slow-loop step in `meds_vegetation_dynamics`, sequenced next to `advance_leaf_phenology`, AFTER growth +
the `overtopping_lai` refresh so the light environment is current. It flattens the per-PFT plasticity
coefficients + `_toc` values, calls the kernel per cohort for the acclimated targets, and `relax_trait`s
the cohort's live traits toward them.

## 4. Parameters (new per-PFT, `[trait_dynamics]` block)
`kplastic_sla`, `kplastic_vm0`, `kplastic_rd`, `kplastic_llspan` (light-response slopes), `lnexp_min/max`
(exponent clamps), `trait_relax_tau` (gradual-mode timescale); and **`leaf_lifespan_toc` [yr]** replaces
`leaf_turnover_rate` (§3e). Defaults reproduce ED2, and **plasticity OFF** (all `kplastic_* = 0`) must be
a supported, static-equivalent configuration (opt-in, like phenology).

## 5. Thermal acclimation — DEFERRED
Out of scope for this cut. Noted for the future: the reserved `t_acclim` running-means on
`wood_env_t`/`root_env_t` ("RESERVED; unused v1") are temperature-acclimation state, and this module is
their eventual home (a running-mean tissue temperature shifting the peaked-Arrhenius reference). Not
implemented now.

## 6. Risks & phasing
Implemented as **one increment** (all four plastic traits together, per §1), with plasticity **OFF by
default** so the state/fusion machinery can land bit-identical to the static path before any coefficient
is turned on.

- **Load-bearing risk: fusion weighting + leaf-area bookkeeping** (§3c/§3d). Add all four cohort fields to
  the lockstep machinery AND the fusion leaf-area weighting in the same change; a green suite with
  `kplastic_* = 0` proves the plumbing before behavior turns on.
- **Validation:** OFF ⇒ bit-identical to the static path; a shaded vs sunlit cohort diverges in the ED2
  direction (higher SLA + smaller leaf biomass, lower Vcmax, longer llspan ⇒ slower turnover); fusion of
  two cohorts with different traits yields the leaf-area-weighted mean and conserves leaf/wood carbon +
  plant number. Build nvfortran multicore too (new SoA fields on the offloaded appliers).

## 7. Open questions
1. **Gradual vs instant default:** recommend gradual (`relax_trait`, governor-style) for numerical
   smoothness; ED2's instantaneous snap is `tau -> 0`.
2. **`leaf_lifespan_toc` vs keeping `leaf_turnover_rate`:** §3e ties turnover to lifespan (recommended,
   removes a redundant parameter). Confirm no other consumer needs the standalone rate.
3. **Rd plasticity coupling:** Rd is derived from Vcmax via `rd_vcmax_ratio` in the static path — decide
   whether plastic Rd follows plastic Vcmax through that ratio (simplest) or gets its own `kplastic_rd`.
