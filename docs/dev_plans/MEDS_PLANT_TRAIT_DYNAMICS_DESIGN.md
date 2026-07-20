# MEDS plant trait dynamics — DESIGN + IMPLEMENTED (2026-07-20)

Status: **IMPLEMENTED**. `meds_plant_trait_dynamics` kernel (`light_plastic_traits` + replacement-weighted
`relax_trait`); the four traits `sla/vcmax25/rd25/llspan` are now per-cohort SoA state (dropped `p_`
prefix), seeded from PFT top-of-canopy values, threaded through the lockstep machinery, and **leaf-area-
weighted on fusion**. `leaf_lifespan_toc` replaced `leaf_turnover_rate` (baseline leaf turnover =
`1/cohort%llspan`); `kplastic_*` derived per-PFT in `derive_pft_rates` (ED2 scheme-2). The slow-loop
driver step `advance_trait_dynamics` (gated on `[trait_dynamics].trait_plasticity_on`, default off) runs
before `carbon_growth`, which consumes `cohort%sla` (leaf-area-conserving target) and `1/cohort%llspan`;
`leaf_gas_exchange` takes per-cohort `vcmax25/rd25` overrides. Validation: ifx Debug + nvfortran multicore
**33/33**, OFF path bit-identical; a plasticity-ON 40-yr spin-up is stable (area conserved, no NaNs) with a
denser shade-acclimated canopy. Thermal acclimation still deferred (§5).

Original design follows.

Status (original): **design-only**. Proposes `src/plant/meds_plant_trait_dynamics.f90` — the home for processes
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

The update is **gradual and leaf-replacement-limited**: a leaf's trait is fixed once it flushes, so the
cohort-mean trait relaxes toward the target only as fast as leaves are **replaced** — weighted by the
leaf turnover, i.e. the current leaf lifespan. ED2 then re-derives geometry from `bleaf` because its size
anchor is leaf biomass — MEDS does **not** (§3d).

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
returns the acclimated TARGET traits (`trait_toc * exp(clamp(kplastic·cum_lai_above))`). A sibling
**replacement-weighted** relaxer then moves each live trait toward its target at the leaf-turnover rate:

```fortran
elemental pure subroutine relax_trait(current, target, llspan, dt, updated)
   ! f = 1 - exp(-dt / llspan)   ! the fraction of leaves REPLACED this step (uses the cohort's CURRENT
   ! updated = current + f*(target - current)   ! plastic llspan) -- ED2's assumption: only new leaves
end subroutine                                  ! carry the new trait, so plasticity is turnover-limited
```

Weighting by the **current `llspan`** (not a free `tau`) is the mechanism ED2 uses and comment-driven
here: a long-lived (deeply shaded) canopy acclimates slowly because it replaces leaves slowly; `llspan`
being itself plastic makes the relaxation self-consistent. There is no `is_instant`/`tau` knob.

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
the cohort's live traits toward them using the cohort's **current `llspan`** as the replacement rate — so
the trait update and the leaf turnover it rides on advance at the same cadence.

## 4. Parameters (new per-PFT, `[trait_dynamics]` block)
Four **separate** light-response slopes — `kplastic_sla`, `kplastic_vm0`, `kplastic_rd`, `kplastic_llspan`
— plus the exponent clamps `lnexp_min/max`. **No relaxation-timescale parameter**: the update rate is the
leaf turnover `1/llspan` (§3a). And **`leaf_lifespan_toc` [yr]** replaces `leaf_turnover_rate` (§3e). The
relaxation is turnover-limited, so there is no separate acclimation timescale to set.

**Rd gets its own `kplastic_rd`** (ED2 keeps `kplastic_rd0` a distinct slot, defaulting it to
`kplastic_vm0` — matching that here).

**Defaults reproduce ED2's `trait_plasticity_scheme = 2` (Lloyd et al. 2010, canopy trait gradients),
derived per-PFT in `derive_pft_rates`** (the same pattern as the wood-density → mortality derivation), so
the PFT file carries base traits + `leaf_lifespan_toc` and the slopes are derived unless overridden:

| slope | ED2 default (derived from base traits) |
|---|---|
| `kplastic_vm0` | `-exp(-2.788 + 0.01439·Vcmax25)`, clamped to `[lnexp_min, lnexp_max]` (negative ⇒ Vcmax ↓ in shade) |
| `kplastic_rd`  | defaults to `kplastic_vm0` (own slot; ED2 `kplastic_rd0`) |
| `kplastic_sla` | from the SLA canopy-gradient exponent `eplastic_sla ≈ -1.18` (positive ⇒ SLA ↑ in shade) |
| `kplastic_llspan` | `0.2126 - 0.062·ln(12 / leaf_turnover)` — i.e. from leaf lifespan |

(ED2 exponents: `eplastic_vm0 ≈ -1.10`, `eplastic_sla ≈ -1.18`; the BCI-tropical `scheme = 3` fits are a
documented alternative.) **Plasticity OFF** (all `kplastic_* = 0`) must be a supported, static-equivalent
configuration (opt-in, like phenology).

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

## 7. Resolved decisions (were open questions)
1. **Update = gradual + leaf-replacement-weighted** by the current `llspan` (§3a) — ED2's assumption that
   only newly-flushed leaves carry the new trait. No `is_instant`/`tau` knob.
2. **`leaf_lifespan_toc` replaces `leaf_turnover_rate`** (§3e); baseline leaf turnover = `1/llspan`.
3. **Rd has its own `kplastic_rd`** (§4), defaulting to `kplastic_vm0` as in ED2.
4. **Defaults are ED2's** (`trait_plasticity_scheme = 2`, Lloyd et al. 2010), derived per-PFT in
   `derive_pft_rates` (§4).

Remaining to confirm at implementation: nothing outside `turnover_shed_rates` still needs the standalone
`leaf_turnover_rate` once it is derived from `1/llspan`.
