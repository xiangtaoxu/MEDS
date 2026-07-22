# Plant trait dynamics

Leaf traits are not static PFT constants — they **acclimate to the light environment** within a plant's
lifetime. MEDS models this with a small **stateless, elemental** kernel, `meds_plant_trait_dynamics`, that
ports ED2's `update_cohort_plastic_trait` (Lloyd et al. 2010 canopy trait gradients). The first cut covers
**light-driven plasticity** of four cohort-level leaf traits — specific leaf area (SLA), photosynthetic
capacity ($`V_{\mathrm{cmax25}}`$), leaf dark respiration ($`R_{\mathrm{d25}}`$), and leaf lifespan
($`\ell`$). Ma et al. (2025) show that constraining exactly this plasticity with canopy observations
markedly improves predicted tropical-forest demography, structure, and biomass. Thermal acclimation is a
deferred sibling. See `docs/dev_plans/MEDS_PLANT_TRAIT_DYNAMICS_DESIGN.md`.

The kernel is two pure functions (targets, then a relaxation step); the per-cohort trait **state**, the
orchestration, and every carbon consequence live in the driver (`meds_vegetation_dynamics`) and the core
Structure-of-Arrays — the same kernel/driver/state split as phenology and carbon allocation.

## 1. The canopy-light trait gradient

A shaded leaf is thinner (higher SLA), lower-capacity (lower $`V_{\mathrm{cmax}}`$/$`R_\mathrm{d}`$), and
longer-lived. Each trait acclimates toward a **target** that scales its top-of-canopy value ($`\cdot_{\mathrm{toc}}`$)
by an exponential gradient in the **cumulative LAI above** the cohort, $`L_{\uparrow}`$ (the light-competition
proxy MEDS already tracks as `overtopping_lai`):

```math
X^{\mathrm{tgt}} = X_{\mathrm{toc}}\,\exp\!\big(k_X\, L_{\uparrow}\big), \qquad X \in \{\mathrm{SLA}, V_{\mathrm{cmax25}}, R_{\mathrm{d25}}, \ell\} \qquad(1)
```

with the exponent clamped to $`[-30, 30]`$ (an overflow guard, not science) and $`L_{\uparrow}`$ floored at
0 (a top-of-canopy cohort sees no shade, so its targets equal the `_toc` values). The **sign** of each
slope $`k_X`$ sets the direction: $`k_{\mathrm{SLA}}>0`$ and $`k_{\ell}\ge0`$ (rise in shade),
$`k_{V}, k_{R}<0`$ (fall in shade).

## 2. Replacement-weighted update

Plasticity is **leaf-turnover-limited**: a leaf's trait is fixed once it flushes, so the cohort-mean trait
can only move toward its target as fast as leaves are **replaced**. Each slow step the live trait relaxes by
the fraction of the canopy replaced, set by the cohort's *current* leaf lifespan $`\ell`$:

```math
X \leftarrow X + f\,\big(X^{\mathrm{tgt}} - X\big), \qquad f = 1 - e^{-\Delta t/\ell} \qquad(2)
```

A long-lived (deeply shaded) canopy therefore acclimates slowly, and — because $`\ell`$ is itself plastic —
the relaxation is self-consistent. There is no free acclimation-timescale parameter and no `is_instant`
switch; instead a single `instant` flag jumps straight to the target ($`f=1`$), used **once** on a census
restart (§5), where the cohorts sit in an established stand but carry no trait history.

## 3. The four slopes (derived per PFT)

The slopes $`k_X`$ reproduce ED2's `trait_plasticity_scheme = 2` (Lloyd et al. 2010), **derived per PFT** in
`derive_pft_rates` from the base traits — the same pattern MEDS uses for the wood-density→mortality
derivation — so the PFT file carries only base traits and `leaf_lifespan_toc` unless a slope is overridden:

```math
k_{V} = -\exp\!\big(\mathrm{clamp}(-2.788 + 0.01439\,V_{\mathrm{cmax25}})\big), \qquad k_{R} = k_{V} \qquad(3)
```

```math
k_{\mathrm{SLA}} = -k_{V}\,\tfrac{1.18}{1.10}, \qquad k_{\ell} = \max\!\Big(0,\; 0.2126 - 0.062\,\ln(12\,\ell_{\mathrm{toc}})\Big) \qquad(4)
```

$`R_\mathrm{d}`$ has its own slope (defaulting to $`k_V`$, as ED2's `kplastic_rd0`). The **cap at 0** on
$`k_{\ell}`$ is a MEDS choice: ED2's fitted lifespan relationship crosses zero near $`\ell_{\mathrm{toc}}\approx2.6`$ yr,
which would make *long-lived* PFTs shorten their leaves in shade — unphysical — so those PFTs instead keep
one canopy-invariant lifespan, while short-lived PFTs still lengthen.

## 4. Carbon consequences

Two of the four traits feed back into the carbon budget, so plasticity is not a passive diagnostic:

- **Leaf lifespan → turnover.** Baseline leaf turnover is not an independent parameter; it is
  $`1/\ell`$ [yr⁻¹]. The driver feeds $`1/\mathrm{cohort\%llspan}`$ to the phenology turnover floor
  (`turnover_shed_rates`), so a shaded, longer-lived canopy sheds and rebuilds its leaves more slowly.

- **SLA is leaf-area-conserving.** The leaf carbon target is $`L^{*} = A_{\mathrm{leaf}}(\mathrm{dbh},h)/\mathrm{SLA}`$
  (`size2leaf_carbon`): allometry fixes the leaf **area** (crown light-capture), and SLA sets how much
  **carbon** that area costs. When acclimation raises SLA, a cohort reaches the same allometric area with
  *less* carbon — the leaf pool now **overshoots** the lower target, so `advance_trait_dynamics` resorbs the
  excess leaf carbon to **storage** (leaf area stays at allometry; carbon is conserved). The reverse
  (SLA falling on un-shading) opens a deficit that normal flush-limited growth fills over time, not an
  instantaneous storage draw — a deliberate asymmetry matching the growth path.

Plastic $`V_{\mathrm{cmax25}}`$ and $`R_{\mathrm{d25}}`$ enter the fast loop through leaf gas exchange (§5):
$`J_{\mathrm{max}}`$ and TPU scale with the plastic $`V_{\mathrm{cmax}}`$ via their fixed ratios, while
$`R_\mathrm{d}`$ carries its own value.

## 5. Interface with other modules

`meds_plant_trait_dynamics` `use`s only `meds_kinds` + `meds_constants` — a pure kernel. State, cadence, and
consumers live outside it.

- **State (`meds_core_state_types`).** The four traits are per-cohort SoA fields **without** the `p_`
  prefix — the convention that marks them *dynamic* (`p_<trait>` = a cached static PFT constant). They are
  seeded from the PFT top-of-canopy values at birth, threaded through the one centralized lockstep
  machinery, and **leaf-area-weighted on cohort fusion** (intensive per-leaf-area properties).

- **Orchestration (`meds_vegetation_dynamics.advance_trait_dynamics`).** A slow-loop step, opt-in via
  `[trait_dynamics].trait_plasticity_on` (default off ⇒ traits stay at top-of-canopy, **bit-identical** to
  the static path). It runs **after** the competition sweep refreshes `overtopping_lai` and **before**
  `carbon_growth`, so this step's leaf demand and turnover see the updated traits. Per cohort it calls
  `light_plastic_traits` (eq 1) then `update_plastic_trait` (eq 2), then resorbs any SLA overshoot (§4).

- **Consumers.** `carbon_growth` uses `cohort%sla` for the leaf-area-conserving target; the phenology
  turnover floor uses `1/cohort%llspan`; and `meds_plant_interface.leaf_gas_exchange` takes per-cohort
  `vcmax25`/`rd25` overrides (threaded through the fast-loop column buffer in `meds_fast_split`).

- **Persistence (`meds_io`).** The four trait states are written to and read from the state checkpoint
  (tolerant of pre-feature files, which restart at top-of-canopy), so a **state** restart recovers the
  acclimated traits exactly. A **census** restart instead acclimates instantaneously in `meds_main`
  (`advance_trait_dynamics(..., instantaneous=.true.)`), because census cohorts carry no history.

## Parameters

| Symbol | Config key | Meaning |
|---|---|---|
| $`k_{\mathrm{SLA}}, k_{V}, k_{R}, k_{\ell}`$ | `kplastic_sla`, `kplastic_vm0`, `kplastic_rd`, `kplastic_llspan` | light-response slopes — **derived** per PFT (§3), overridable |
| $`\ell_{\mathrm{toc}}`$ | `leaf_lifespan_toc` | top-of-canopy leaf lifespan [yr]; baseline turnover = $`1/\ell`$ (replaced `leaf_turnover_rate`) |
| $`\mathrm{SLA}_{\mathrm{toc}}, V_{\mathrm{toc}}, R_{\mathrm{toc}}`$ | `sla`, `vcmax25`, `rd25` | top-of-canopy trait values (base PFT traits) |
| — | `[trait_dynamics].trait_plasticity_on` | opt-in master switch (default off) |

## Code map

| Concept | Routine |
|---|---|
| light target traits (eq 1) | `meds_plant_trait_dynamics`: `light_plastic_traits` (+ private `light_gradient`) |
| replacement-weighted / instant update (eq 2) | `meds_plant_trait_dynamics`: `update_plastic_trait` |
| derived slopes (eqs 3–4) | `meds_pft_params`: `derive_pft_rates` |
| per-cohort trait state + fusion weighting | `meds_core_state_types` (SoA `sla/vcmax25/rd25/llspan`); `meds_core_cohort_fusefiss`: `fuse_2_cohorts` |
| slow-loop orchestration + SLA-overshoot resorption | `meds_vegetation_dynamics`: `advance_trait_dynamics` |
| carbon / gas-exchange consumers | `carbon_growth` (`cohort%sla`, `1/cohort%llspan`); `meds_plant_interface`: `leaf_gas_exchange` (`vcmax25`/`rd25`) |
| state persistence + census instant | `meds_io`: `io_write_state`/`io_read_state`; `meds_main` (census restart) |

## References
- **Ma, Y., Moorcroft, P. R., Wright, S. J., Rogers, A., Lamour, J., Davidson, K. J., Serbin, S. P.,
  Detto, M., & Xu, X. (2025).** Constraining light-driven plasticity in leaf traits with observations
  improves the prediction of tropical forest demography, structure, and biomass dynamics. *Journal of
  Geophysical Research: Biogeosciences*, 130(6), e2025JG008814. https://doi.org/10.1029/2025JG008814 —
  the motivating study: how light-driven leaf-trait plasticity, constrained by canopy observations,
  reshapes predicted forest demography, structure, and biomass.
- **Lloyd et al. (2010)**, *Biogeosciences* 7:1833 — within-canopy gradients of foliar traits (the
  plasticity relationships and `trait_plasticity_scheme = 2` defaults).
- ED2 `ED/src/utils/update_derived_utils.f90` (`update_cohort_plastic_trait`), `ED/src/init/ed_params.f90`
  (the `kplastic_*` derivations).
