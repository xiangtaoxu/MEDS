# Plant trait dynamics

Leaf traits are not static PFT constants — they **acclimate to the light environment** within a plant's
lifetime. MEDS models this with a small **stateless, elemental** kernel, `meds_plant_trait_dynamics`, that
ports ED2's `update_cohort_plastic_trait` (Lloyd et al. 2010 canopy trait gradients). It covers
**light-driven plasticity** of four cohort-level leaf traits — specific leaf area (SLA), photosynthetic
capacity ($`V_{\mathrm{cmax25}}`$), leaf dark respiration ($`R_{\mathrm{d25}}`$), and the leaf turnover
rate $`\gamma`$ [yr⁻¹], the share of the leaf carbon replaced a year (`plant_phenology.md`'s
$`k_{\mathrm{turn}}`$). Ma et al. (2025) show that constraining exactly this plasticity with canopy observations
markedly improves predicted tropical-forest demography, structure, and biomass. Thermal acclimation of
photosynthetic capacity is a separate, site-level mechanism (`leaf_gas_exchange.md` §"Thermal
acclimation"); per-cohort acclimation is open (see `docs/ROADMAP.md`).

The kernel is two pure functions (targets, then a relaxation step); the per-cohort trait **state**, the
orchestration, and every carbon consequence live in the driver (`meds_vegetation_dynamics`) and the core
Structure-of-Arrays — the same kernel/driver/state split as phenology and carbon allocation.

## 1. The canopy-light trait gradient

A shaded leaf is thinner (higher SLA), lower-capacity (lower $`V_{\mathrm{cmax}}`$/$`R_\mathrm{d}`$), and
longer-lived, so a shaded canopy turns its leaves over more slowly. Each trait acclimates toward a
**target** that scales its top-of-canopy value ($`\cdot_{\mathrm{toc}}`$) by an exponential gradient in
the **cumulative LAI above** the cohort, $`L_{\uparrow}`$ (the light-competition proxy MEDS already tracks
as `overtopping_lai`):

```math
X^{\mathrm{tgt}} = X_{\mathrm{toc}}\,\exp\!\big(k_X\, L_{\uparrow}\big), \qquad X \in \{\mathrm{SLA}, V_{\mathrm{cmax25}}, R_{\mathrm{d25}}, \gamma\} \qquad(1)
```

with the exponent clamped to $`[-30, 30]`$ (an overflow guard, not science) and $`L_{\uparrow}`$ floored at
0 (a top-of-canopy cohort sees no shade, so its targets equal the `_toc` values). The **sign** of each
slope $`k_X`$ sets the direction: $`k_{\mathrm{SLA}}>0`$ (rises in shade), $`k_{V}, k_{R}<0`$ and
$`k_{\gamma}\le0`$ (fall in shade).

For the turnover rate, eq. (1) reads $`\gamma = \gamma_{\mathrm{toc}}\,e^{-|k_\gamma| L_{\uparrow}}`$. It is
algebraically the reciprocal of the form ED2 and Lloyd et al. (2010) give for the leaf lifespan
$`\ell = 1/\gamma`$, which lengthens in shade as $`\ell = \ell_{\mathrm{toc}}\,e^{|k_\gamma| L_{\uparrow}}`$: the
target is the same; only the quantity carried is the rate.

## 2. Replacement-weighted update

Plasticity is **leaf-turnover-limited**: a leaf's trait is fixed once it flushes, so the cohort-mean trait
can only move toward its target as fast as leaves are **replaced**. Each slow step the live trait relaxes by
the fraction of the canopy replaced, set by the cohort's *current* leaf turnover rate $`\gamma`$:

```math
X \leftarrow X + f\,\big(X^{\mathrm{tgt}} - X\big), \qquad f = 1 - e^{-\gamma\,\Delta t} \qquad(2)
```

A slow-turnover (deeply shaded) canopy therefore acclimates slowly, and — because $`\gamma`$ is itself
plastic — the relaxation is self-consistent. The update is the leaf-carbon mean of the leaves kept and the
leaves replaced; for $`\gamma`$ that mean is the canopy's turnover, the leaf carbon the two sets of leaves
lose a year divided by their carbon. There is no free acclimation-timescale parameter; a single `instant`
flag jumps straight to the target ($`f=1`$), used **once** on a census restart (§5), where the cohorts sit
in an established stand but carry no trait history.

## 3. The four slopes (per PFT)

The slopes $`k_X`$ are PFT traits, `[pft].kplastic_vcmax`, `kplastic_rd`, `kplastic_sla` and
`kplastic_leaf_turnover`. The shipped PFTs take ED2's `trait_plasticity_scheme = 2` fits (Lloyd et al.
2010) to their base traits, which `scripts/pft_traits` computes outside the model:

```math
k_{V} = -\exp\!\big(\mathrm{clamp}(-2.788 + 0.01439\,V_{\mathrm{cmax25}})\big), \qquad k_{R} = k_{V} \qquad(3)
```

```math
k_{\mathrm{SLA}} = -k_{V}\,\tfrac{1.18}{1.10}, \qquad k_{\gamma} = \min\!\Big(0,\; 0.062\,\ln(12/\gamma_{\mathrm{toc}}) - 0.2126\Big) \qquad(4)
```

$`R_\mathrm{d}`$ has its own slope ($`k_V`$ in the fit, as ED2's `kplastic_rd0`). $`k_\gamma`$ is minus ED2's
slope for the leaf lifespan, $`0.2126 - 0.062\ln(12\,\ell_{\mathrm{toc}})`$ with $`12\,\ell_{\mathrm{toc}} = 12/\gamma_{\mathrm{toc}}`$
the lifespan in months. The **cap at 0** is a MEDS choice: ED2's fit crosses zero near
$`\gamma_{\mathrm{toc}}\approx0.39`$ yr⁻¹ (a 2.6-yr lifespan), which would make *slow-turnover* PFTs turn their
leaves over faster in shade — unphysical — so those PFTs instead keep one canopy-invariant rate, while
fast-turnover PFTs still slow down.

## 4. Carbon consequences

Two of the four traits feed back into the carbon budget, so plasticity is not a passive diagnostic:

- **The turnover rate is the turnover.** Baseline leaf turnover is the cohort's $`\gamma`$ [yr⁻¹]. The
  driver feeds `cohort%leaf_turnover_rate` to the background turnover of `leaf_turnover_step`, so a
  shaded canopy sheds and rebuilds its leaves more slowly. With
  `[trait_dynamics].fineroot_turnover_plastic` the fine roots turn over more slowly in the same
  proportion, `fineroot_turnover_rate` $`\times\,\gamma/\gamma_{\mathrm{toc}}`$.

- **SLA is leaf-area-conserving.** The leaf carbon target is $`L^{*} = A_{\mathrm{leaf}}(\mathrm{dbh},h)/\mathrm{SLA}`$
  (`size_to_leaf_carbon`): allometry fixes the leaf **area** (crown light-capture), and SLA sets how much
  **carbon** that area costs. When acclimation raises SLA, a cohort reaches the same allometric area with
  *less* carbon — the leaf pool now **overshoots** the lower target, so `advance_plant_traits` resorbs the
  excess leaf carbon to **storage** (leaf area stays at allometry; carbon is conserved). The reverse
  (SLA falling on un-shading) opens a deficit that normal flush-limited growth fills over time, not an
  instantaneous storage draw — a deliberate asymmetry matching the growth path.

Plastic $`V_{\mathrm{cmax25}}`$ and $`R_{\mathrm{d25}}`$ enter the fast loop through leaf gas exchange (§5):
$`J_{\mathrm{max}}`$ and TPU scale with the plastic $`V_{\mathrm{cmax}}`$ via their fixed ratios, while
$`R_\mathrm{d}`$ carries its own value.

## 5. Interface with other modules

`meds_plant_trait_dynamics` `use`s only `meds_kinds` + `meds_constants` — a pure kernel. State, cadence, and
consumers live outside it.

- **State (`meds_cohort_state_vars`).** The four traits are per-cohort SoA fields, section 6 of
  `cohort_block` — *dynamic*, one value per cohort, where a static PFT trait is read from the block's
  PFT table (`pft_table%x(pft(i))`) and never copied onto a cohort. They are
  seeded from the PFT top-of-canopy values at birth, threaded through the one centralized lockstep
  machinery, and blended on cohort fusion on the quantity each is per: $`V_{\mathrm{cmax25}}`$ and
  $`R_{\mathrm{d25}}`$ are per leaf area, so leaf area weights them; SLA is leaf area per leaf carbon, so
  leaf carbon weights it and the fused cohort shows the leaf area its two parents showed; $`\gamma`$ is
  the share of the leaf carbon lost a year, so leaf carbon weights it too and the fused cohort turns
  over the leaf carbon its two parents did, which a leaf-area mean does not when their SLAs differ.

- **Orchestration (`meds_vegetation_dynamics.advance_plant_traits`).** A slow-loop step, opt-in via
  `[trait_dynamics].trait_plasticity_on` (default off ⇒ every cohort keeps its PFT's
  top-of-canopy traits). It runs **after** the competition sweep refreshes `overtopping_lai` and **before**
  `compute_carbon_allocation`, so this step's leaf demand and turnover see the updated traits. Per cohort it calls
  `light_plastic_traits` (eq. (1)) then `update_plastic_trait` (eq. (2)), then resorbs any SLA overshoot (§4).

- **Consumers.** `compute_carbon_allocation` uses `cohort%sla` for the leaf-area-conserving target; the background
  leaf turnover uses `cohort%leaf_turnover_rate`; and the fast loop's `leaf_gas_exchange_batch` (called from
  `meds_fast_prepass.canopy_leaf_gas_exchange`) takes per-cohort `vcmax25`/`rd25` (threaded through the
  fast-loop column buffer in `meds_fast_types`), with `leaf_params_at_capacity` scaling
  $`J_{\mathrm{max25}}`$ and TPU with the plastic $`V_{\mathrm{cmax25}}`$. The single-leaf wrapper
  `meds_fast_config.leaf_gas_exchange` takes the same two overrides.

- **Persistence (`meds_init_restart`).** The four trait states are written to and read from the state checkpoint,
  so a **state** restart recovers the acclimated traits exactly; `[init].reacclimate_traits` instead re-derives them from this
  run's PFT file. A **census** start acclimates instantaneously when the run opens (`region_open`,
  `advance_plant_traits(..., instantaneous=.true.)`), because census cohorts carry no history.

## Parameters

| Symbol | Config key | Meaning |
|---|---|---|
| $`k_{\mathrm{SLA}}, k_{V}, k_{R}, k_{\gamma}`$ | `[pft]` `kplastic_sla`, `kplastic_vcmax`, `kplastic_rd`, `kplastic_leaf_turnover` | light-response slopes per PFT [1/(m² m⁻²)] (§3; `scripts/pft_traits` for the fits) |
| $`\gamma_{\mathrm{toc}}`$ | `leaf_turnover_rate_toc` | top-of-canopy baseline leaf turnover rate [yr⁻¹], > 0 |
| $`\mathrm{SLA}_{\mathrm{toc}}, V_{\mathrm{toc}}, R_{\mathrm{toc}}`$ | `sla`, `vcmax25`, `rd25` | top-of-canopy trait values (base PFT traits) |
| — | `[trait_dynamics].trait_plasticity_on` | opt-in master switch (default off) |
| — | `[trait_dynamics].fineroot_turnover_plastic` | fine-root turnover scaled by $`\gamma/\gamma_{\mathrm{toc}}`$ (default off) |

## Code map

| Concept | Routine |
|---|---|
| light target traits (eq. (1)) | `meds_plant_trait_dynamics`: `light_plastic_traits` (+ private `light_gradient`) |
| replacement-weighted / instant update (eq. (2)) | `meds_plant_trait_dynamics`: `update_plastic_trait` |
| the slopes' fits (eqs. (3)–(4)) | `scripts/pft_traits/pft_traits.py`: `plasticity` |
| per-cohort trait state + fusion weighting | `meds_cohort_state_vars` (SoA `sla/vcmax25/rd25/leaf_turnover_rate`); `meds_demography_cohort_fusefiss`: `fuse_2_cohorts` |
| slow-loop orchestration + SLA-overshoot resorption | `meds_vegetation_dynamics`: `advance_plant_traits` |
| carbon / gas-exchange consumers | `compute_carbon_allocation` (`cohort%sla`, `cohort%leaf_turnover_rate`); `meds_leaf_gas_exchange`: `leaf_gas_exchange_batch`, `leaf_params_at_capacity` (`vcmax25`/`rd25`) |
| state persistence + census instant | `meds_init_restart`: `restart_file_write`/`read_restart_file`; `meds_init`: `init_polygon` (census start); `meds_vegetation_dynamics`: `reacclimate_plant_traits` |

## References
- **Ma, Y., Moorcroft, P. R., Wright, S. J., Rogers, A., Lamour, J., Davidson, K. J., Serbin, S. P.,
  Detto, M., & Xu, X. (2025).** Constraining light-driven plasticity in leaf traits with observations
  improves the prediction of tropical forest demography, structure, and biomass dynamics. *Journal of
  Geophysical Research: Biogeosciences*, 130(6), e2025JG008814. https://doi.org/10.1029/2025JG008814 —
  the motivating study: how light-driven leaf-trait plasticity, constrained by canopy observations,
  reshapes predicted forest demography, structure, and biomass.
- **Lloyd et al. (2010)**, *Biogeosciences* 7:1833 — within-canopy gradients of foliar traits (the
  plasticity relationships and `trait_plasticity_scheme = 2` defaults).
- ED2 `ED/src/utils/update_derived_utils.f90` (`update_cohort_plastic_trait`), `ED/src/io/init/ed_params.f90`
  (the `kplastic_*` derivations).
