# MEDS plant trait dynamics — DESIGN (review draft, not implemented)

Status: **design-only**. Proposes `src/plant/meds_plant_trait_dynamics.f90` — the home for processes
that CHANGE a plant's traits within its lifetime: light-driven **plasticity** and (later) thermal
**acclimation**. First cut reproduces ED2's `update_cohort_plastic_trait`
(`../ED2/ED/src/utils/update_derived_utils.f90`).

## 1. Motivation & scope

`meds_pft_params` holds the **static** per-PFT trait table. But several traits are not static — SLA,
Vcmax (Vm0), leaf dark respiration Rd, and leaf lifespan **acclimate to the light environment** (a
shaded leaf is thinner / higher-SLA, lower-Vcmax), and photosynthetic capacities acclimate to
temperature. MEDS has no home for these *dynamic* trait processes today. This module is that home,
kept separate from the static table (hence NOT `meds_plant_traits`, which would read as the table).

Scope of the first cut = **light plasticity** (ED2 `update_cohort_plastic_trait`): given a cohort's
light environment, adjust SLA / Vcmax / Rd / leaf-lifespan away from their top-of-canopy PFT values.
Thermal acclimation (the reserved `t_acclim` running-means, §5) is a planned sibling in the same module.

## 2. The ED2 process (what we port)

`update_cohort_plastic_trait` scales top-of-canopy traits by the **cumulative LAI above** the cohort
(`max_cum_lai`), the light-competition proxy:

```
new_sla   = sla_toc  * exp(kplastic_sla * max_cum_lai)          (shaded => higher SLA)
new_vm0   = vm0_toc  * exp(kplastic_vm0 * max_cum_lai)          (kplastic_vm0 < 0 => shaded => lower Vcmax)
new_rd, new_llspan  scale analogously; exponents clamped to [lnexp_min, lnexp_max].
```

Two update modes: `is_instant` (snap to the target) or gradual (`trait_change_frac`, a fractional step
toward the target). ED2 then re-derives geometry from `bleaf` because its size anchor is leaf-biomass.

## 3. MEDS design

Follows the established split — **stateless kernel in `plant`, orchestration in the driver, state in
`core`** — exactly like phenology and carbon allocation.

### 3a. Kernel (`meds_plant_trait_dynamics`, `plant` library)
Pure/elemental, scalar-in/scalar-out (GPU/SoA-safe, issue-#7 N/A):

```fortran
elemental pure subroutine light_plastic_traits(cum_lai_above, sla_toc, vcmax25_toc, rd25_toc, &
        llspan_toc, kplastic_sla, kplastic_vm0, kplastic_rd, kplastic_llspan,                   &
        lnexp_min, lnexp_max, sla, vcmax25, rd25, llspan)   ! -> the light-acclimated traits
```
It takes the top-of-canopy PFT values + plasticity coefficients and returns the acclimated traits. A
sibling `relax_trait(current, target, tau, dt)` gives the **gradual** mode (a rate toward the target,
matching the phenology governor style) so the kernel need not know `is_instant`.

### 3b. Light input — reuse `overtopping_lai`
MEDS already maintains the cumulative-LAI-above diagnostic per cohort (`update_overtopping_lai`,
`meds_core`) — the direct analogue of ED2's `max_cum_lai`. No new light machinery is needed; the driver
passes `cohort%overtopping_lai(j)` into the kernel.

### 3c. Per-cohort trait STATE (the real work)
Plastic traits must be **per-cohort state**, not PFT-static:
- **SLA is already per-cohort** — `cohort%p_sla` exists — so plastic SLA has its substrate today.
- **Vcmax / Rd / leaf-lifespan need NEW per-cohort fields** on `cohort_block` (`p_vcmax25`, `p_rd25`,
  `p_llspan`). Adding them means updating the ONE centralized lockstep machinery in
  `meds_core_state_types` (`copy_cohort_slot` / `cohort_reorder` / `set_cohort_size` / capacity) and,
  crucially, **cohort FUSION weighting** (leaf-area-weighted, like the other cached per-cohort params) —
  this SoA + fusion work, not the kernel, is the bulk of the effort.
- `leaf_gas_exchange` currently flattens `cfg%pft%vcmax25(ipft)`; once Vcmax is plastic it must read
  `cohort%p_vcmax25(j)` instead — a small change at the flattening seam.

### 3d. A MEDS simplification (do NOT port ED2's bleaf coupling)
ED2 re-derives `dbh/hite` from `bleaf` when SLA changes. MEDS anchors size on **`wood_carbon`**, so SLA
plasticity changes **leaf display** (`leaf_area = leaf_carbon·sla` → LAI → light competition) but **not**
the dbh/height size anchor. That decoupling is cleaner and avoids ED2's `bl2dbh`/`bl2h` circularity —
the only feedback is SLA → leaf_area → overtopping_lai → next step's plasticity (a benign slow loop).

### 3e. Orchestration & cadence (driver)
A slow-loop step in `meds_vegetation_dynamics`, sequenced next to `advance_leaf_phenology` (traits
acclimate on the same daily/monthly cadence, AFTER growth + the overtopping-LAI refresh so the light
environment is current). It flattens the per-PFT plasticity coefficients, calls the kernel per cohort,
and writes the acclimated traits back to the cohort SoA (gradual relax by default).

## 4. Parameters (new per-PFT, `[trait_dynamics]` block)
`kplastic_sla`, `kplastic_vm0`, `kplastic_rd`, `kplastic_llspan` (light-response slopes), `lnexp_min/max`
(clamps), and `trait_relax_tau` (gradual-mode timescale). Defaults reproduce ED2; **plasticity OFF**
(all `kplastic_* = 0`) must be a supported, bit-identical-to-static configuration (opt-in, like phenology).

## 5. Thermal acclimation (planned sibling, not first cut)
The reserved `t_acclim` running-means on `wood_env_t`/`root_env_t` ("RESERVED; unused v1") are temperature
acclimation state. This module is their home: a `thermal_acclimation` kernel advancing a running-mean
tissue temperature and shifting the peaked-Arrhenius reference (ED2/Kumarathunge form), consolidating all
trait-change processes in one place.

## 6. Risks & phasing
- **P0 (kernel + SLA only):** `light_plastic_traits` + `relax_trait`, wired to the existing `p_sla`;
  no new SoA fields. Validates the light→trait→display feedback with minimal state churn.
- **P1 (Vcmax/Rd/llspan):** add the per-cohort fields + lockstep + **fusion weighting** (the load-bearing
  step — a new per-cohort field silently unaveraged on fusion is the classic bug this repo guards against),
  repoint the `leaf_gas_exchange` flattening at `cohort%p_vcmax25`.
- **P2:** thermal acclimation (§5).
- **Validation:** plasticity-OFF must be bit-identical to the static path; a shaded vs sunlit cohort must
  diverge in the ED2 direction; conserve nothing special (traits are intensive) but confirm fusion of two
  cohorts with different plastic traits gives the leaf-area-weighted mean. Build nvfortran multicore too.

## 7. Open questions
1. **Name:** `meds_plant_trait_dynamics` (chosen) vs `meds_plant_acclimation`. Trait_dynamics is the
   broad home for plasticity + acclimation; acclimation alone is narrower.
2. **Gradual vs instant** default: recommend gradual (`relax_trait`, governor-style) for numerical
   smoothness; ED2's `is_instant` is a config option.
3. **Vcmax plasticity in P0 or defer to P1?** P0-SLA-only keeps the first PR free of SoA/fusion changes;
   Vcmax is where the state work lands.
