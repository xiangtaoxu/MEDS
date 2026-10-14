# MEDS phenology: leaf habit from phenology parameters alone

**Status:** written 2026-10-05 against `beta` at `413b587`. **Revised after review 1 the same day:
§8 supersedes §2.2, §2.3 and the open questions of §6. §9–§10 (the generic model) were implemented
on 2026-10-06 on `dev/example02-canopy-phenology`; §11 records what the implementation found.** The
fits of §8–§9 are offline fits of a Python emulation, kept in `~/claude_workspace/meds_pheno_proto/`;
`examples/example02_canopy_phenology` refits the compiled kernel.

**Goal:** a PFT's leaf habit (deciduous, evergreen, semi-deciduous) emerges from its phenology
parameters, with no `evergreen` flag, and leaf lifespan becomes an output rather than an input.

## 1. Decisions (owner, 2026-10-05)

| # | Decision |
|---|---|
| D1 | Evergreen is not a parameter; the habit emerges from the phenology parameters. |
| D2 | The model's leaf parameter is a background **`leaf_turnover_rate`**; leaf lifespan is emergent. |
| D3 | A new per-PFT **`senescence_fraction`**: the largest fraction of the canopy one shed episode removes. |
| D4 | The evergreen **cold suppression** of turnover is dropped: it is phenomenological, and D3 does its job. |
| D5 | Phenology drivers are climate or MEDS-simulated plant variables; no soil data in the examples. |

## 2. The scheme

### 2.1 What is unchanged: cues, masks, governors

Each day, for each cohort, the kernel (`meds_phenology.f90`):

1. **Accumulates the cue memory:** growing-degree days above `gdd_base_temp` (Jan–Aug, northern
   equivalent), chilling days below `chill_base_temp` (Nov–Jun), running means of root-zone water
   and shortwave, and counts of consecutive days with predawn leaf ψ below `leaf_psi_tlp` or above
   half of it.
2. **Forms per-cue signals in [0, 1]:**

   | Cue | Flush signal | Shed signal |
   |---|---|---|
   | `TEMP` | degree days past `phen_a + phen_b·exp(phen_c·chill)` (Botta et al. 2000) | days < `cold_drop_daylength` with soil T < `cold_drop_soiltemp1`, or soil T < `cold_drop_soiltemp2` (White et al. 1997) |
   | `WATER` | running-mean water above `water_on_threshold` | below `water_off_threshold` |
   | `HYDRO` | consecutive days with ψ ≥ ½ `leaf_psi_tlp` (Xu et al. 2016) | consecutive days with ψ < `leaf_psi_tlp` |
   | `PHOTO` | day length above `photo_crit` (gates `TEMP`) | — |
   | `LIGHT` | — | running-mean shortwave above `light_on_threshold` (Kim et al. 2012) |

3. **Combines them:** flush = MIN over `flush_cue_mask` (every cue must allow it; empty = 1),
   shed = MAX over `shed_cue_mask` (any cue can force it; empty = 0).
4. **Smooths** each into a governor drive over `tau_flush` / `tau_shed` days and returns
   `r_fl = k_flush_max · flush_drive` and `r_sh = k_shed_max · shed_drive` [day⁻¹].

### 2.2 New: the senescence fraction (D3)

A shed episode is a run of days with `shed_drive > ε`. The kernel carries one new state,
`senescence_exposure` E, the cumulative relative shed of the current episode, and one new
parameter, `senescence_fraction` f ∈ (0, 1]:

```
if shed_drive < ε:  E = 0                                   ! the episode is over
r_sh = min(r_sh, max(0, -ln(1 - f) - E) / dt)               ! f = 1: no limit (bit-identical)
E    = E + r_sh · dt
```

Because the carbon layer sheds a fraction of the current canopy, a canopy that is not regrowing,
stepped daily, keeps about `exp(-E)` of what it had at the episode's start, so the limit `-ln(1 - f)`
stops the episode once about a fraction f is gone. If the flush refills the canopy during the
episode, the episode turns over `-ln(1 - f)` canopies, not f; with `dt_slow` > 1 day the explicit
shed removes more than `1 - exp(-E)`. Review 1 (§8) found this placement unworkable in the coupled
model.

**Basis.** Senescence is programmed for a cohort of leaves, not a rate applied forever: a conifer
holds several needle age classes and sheds mostly the oldest each autumn (Scots pine in Finland:
August–October), a deciduous tree senesces its whole canopy, and a semi-deciduous tropical tree
drops part of its crown in drought. f is the share of the crown that is senescence-competent in one
episode. The cues set **when** leaves senesce, `k_shed_max` **how fast**, and f **how much**.

### 2.3 Changed: the carbon layer (D2, D4)

| | Now | Proposed |
|---|---|---|
| Background turnover | `1/llspan` [yr⁻¹], `llspan` from `leaf_lifespan_toc` | `leaf_turnover_rate` [yr⁻¹] from `leaf_turnover_rate_toc` |
| Cold suppression | × `1/(1+exp(0.4·(5 °C − T)))` for `evergreen = 1`, evaluated at a 25 °C stub (factor 0.9997) | removed, with `evergreen`, `evg_ref_temp`, `evg_slope` |
| Leaf shed rate | `max(r_sh, turnover/365.2425)` × current leaf carbon | unchanged in form |
| Fine-root shed | `fineroot_turnover_rate`, cold-suppressed for evergreens | `fineroot_turnover_rate` |
| Dormant snap, flush cap, resorption (`retained_carbon_fraction`) | — | unchanged |

`leaf_turnover_rate` now means **non-seasonal** background loss (damage, herbivory, shading inside
the crown); for temperate trees it is small (fits: 0.03 yr⁻¹ deciduous; the boreal pine fit used 0,
and the one fit at 0.02 yr⁻¹ is not valid: it shed only 8 % of the canopy a year).

### 2.4 Changed: trait plasticity (D2)

| | Now | Proposed |
|---|---|---|
| Plastic trait | lifespan `ℓ = ℓ_toc · exp(k_ℓ · LAI_above)` | turnover `γ = γ_toc · exp(−k_γ · LAI_above)` |
| Slope | `k_ℓ = max(0, 0.2126 − 0.062 ln(12 ℓ_toc))` (Lloyd et al. 2010, ED2 scheme 2) | `k_γ = max(0, 0.2126 − 0.062 ln(12 / γ_toc))`, the same relation |
| Acclimation per step | `1 − exp(−Δt / ℓ)` | the fraction of the canopy actually regrown this step (leaf growth / leaf carbon) |

With realized replacement, a deciduous canopy takes its acclimated traits as it flushes, and an
evergreen acclimates at the rate it actually replaces leaves.

### 2.5 New diagnostic

`leaf_lifespan_cohort` [yr]: the realized mean leaf residence time, a year's running mean of leaf
carbon divided by a year's running mean of leaf shed. It replaces the `llspan_cohort` output.

## 3. Leaf habits as parameter sets

| Habit | Flush mask | Shed mask | `senescence_fraction` | `leaf_turnover_rate` |
|---|---|---|---|---|
| Temperate deciduous | `TEMP` (+ `PHOTO`) | `TEMP` | 1 | ~0.03 |
| Boreal/temperate conifer | `TEMP` | `TEMP` | ~1 / needle age classes (0.2–0.35) | ~0–0.02 |
| Aseasonal evergreen | — | — | — | 1 / leaf lifespan |
| Tropical drought-deciduous | `HYDRO` | `HYDRO` | 1 | site-specific |
| Semi-deciduous | `HYDRO` | `HYDRO` | 0.3–0.7 | site-specific |
| Tropical leaf exchanger | — | `LIGHT` | 1 | site-specific |

## 4. Evidence

**Harvard Forest, deciduous** (MODIS LAI, HF003 leaf fall, HF069 litter baskets, HF001 air T;
2003–2023). Unaffected by the proposal (f = 1, turnover 0.03). RMSE 0.29 → 0.19 (MODIS relative
LAI), 0.27 → 0.065 (leaf fall), 0.14 → 0.073 (baskets) against the preset.

**Hyytiälä FI-Hyy, Scots pine** (ICOS L2 needle litter, 74 collections of ~30 days, 2018–2024;
ICOS daily air T). Interval needle-fall shares, fits held to ~30 % of the canopy a year:

| Design | Flag | r | Jan–Apr | Sep | Oct | Nov |
|---|---|---|---|---|---|---|
| Observed | — | — | 0.06–0.21 | 3.70 | 3.52 | 0.27 |
| Current kernel | — | 0.26 | 1.08–1.21 | 1.35 | 1.38 | 1.32 |
| Shed × cold factor | yes | 0.69 | 0.33–0.96 | 3.47 | 2.85 | 1.38 |
| **f = 0.23, no turnover, no cold suppression** | **no** | **0.79** | 0.00 | 5.75 | 3.35 | 0.00 |
| f = 0.18 + cold-suppressed turnover 0.19 | no | 0.81 | 0.04–0.26 | 4.84 | 2.72 | 0.17 |

In the no-cold-suppression fit the shed becomes possible from 11 August (days < 16.1 h, air
< 17 °C) at 0.005 day⁻¹. With f = 0.23 the canopy averages about 0.9 and loses 0.25 a year: a
realized needle life of about 3.6 years, not 1/f. Neither design reproduces a June litter peak (0.94), which may be bud scales,
pollen cones or old needles dropped during shoot growth.

## 5. What moves numbers

- **D4** is nearly a no-op today but not bit-identical: the cold factor is evaluated at the 25 °C
  stub (0.99966), so leaf and fine-root turnover of every `evergreen = 1` PFT rise by 3.4e-4
  (relative) in the default, demography and BCI configs.
- **D3** changes nothing until a PFT sets f < 1.
- **D2** moves the plasticity runs (acclimation at realized replacement) and the census tool
  (`make_census.py`), whose steady-state leaf litter `leaf_carbon / lifespan` must become the
  realized annual leaf loss (a whole canopy for a deciduous PFT). Converting configs with
  `leaf_turnover_rate_toc = 1 / leaf_lifespan_toc` is otherwise exact.

## 6. Open questions for review

1. **Shed rate: max or sum?** The leaf shed is `max(active, background)`. With background now
   meaning non-seasonal damage loss, the two are independent processes and might add.
2. **Episode boundary.** E resets when `shed_drive < ε`. With `tau_shed` = 5 d the drive takes
   about a month to decay below 1e-3 after its signal stops. A warm spell in winter that turns the
   `TEMP` shed off could end an episode early and allow a second, spring senescence. Alternatives:
   reset only after the signal has been off for N days, or at a seasonal boundary.
3. **Light exchanger.** The `LIGHT` cue acts only on the shed side; in tropical forests new-leaf
   flushing tracks light (BCI litter peaks a month before the model's). Not addressed here.
4. **Cold-drop temperature terms.** At both sites the soil-temperature threshold ran to its bound
   and day length alone set the autumn shed. Whether to simplify is not addressed here.
5. **Exactness of `-ln(1 - f)`** when flush and shed overlap (light exchanger, drought episodes
   with rain).

## 7. Implementation

1. **Kernel and types:** `senescence_fraction` (params), `senescence_exposure` (state), the limit
   in `phenology_kernel`; the C-API structs and `meds.plant.pheno` (`Params`, `State`).
2. **Carbon layer:** `pheno_drives_to_rates` / `turnover_shed_rates` without `evergreen` and the
   cold factor; `cohort_carbon_demand` without `STUB_TISSUE_TEMP`.
3. **Traits:** `leaf_turnover_rate` cohort state and plasticity (`meds_plant_trait_dynamics`,
   `advance_plant_traits`), realized-replacement acclimation.
4. **Config:** `leaf_turnover_rate_toc`, `[phenology] senescence_fraction`; `leaf_lifespan_toc`,
   `evergreen`, `evg_ref_temp`, `evg_slope` refused with a message naming the replacement; the
   reference `meds_config_pft.toml` and every shipped config converted.
5. **I/O:** `leaf_lifespan_cohort` (realized) and `leaf_turnover_rate_cohort` replace `llspan`.
6. **Census tool:** realized annual leaf loss for the steady-state litter.
7. **Tests:** kernel f = 1 bit-identical; f < 1 stops the episode at `-ln(1 - f)` and resets;
   plasticity under realized replacement; config refusals; the Python mirror.
8. **Docs and CHANGELOG:** `plant_phenology.md` (§1′, §4, §5), `plant_traits.md` (§2, §4),
   before/after numbers for every shipped example.

## 8. Review 1 (2026-10-05) and the revised design (v2)

An independent review read this plan against the code. Its main findings, each checked:

1. **The kernel limit never reaches the carbon.** `advance_leaf_phenology` writes back only the
   governor drives (`meds_vegetation_dynamics.f90:1094-1105`); `cohort_carbon_demand` recomputes
   the shed rate from `shed_drive` (`:971`). A cap inside `phenology_kernel` acts only on the
   C-API/Python path the prototypes used.
2. **The episode reset fails with realistic soil temperature.** The logistic cue never reaches 0,
   so `shed_drive < ε` is a knife edge. Confirmed offline: the §2.2 design driven by ICOS soil
   temperature (instead of air) sheds **nothing** after its first episode.
3. **`-ln(1-f)` is not what the carbon layer removes** (refill during the episode, explicit Euler
   for `dt_slow` > 1 d, kernel rate vs realized shed). A budget in realized carbon is exact.
4. **The evidence for dropping the temperature response was thin:** the no-background fit puts no
   litter in Nov–Jul, where ~22 % of the observed annual needle fall lies (mostly May–Aug).
5. **`leaf_turnover_rate` had two meanings** (damage loss in §2.3, 1/lifespan in §3), the Lloyd
   slope on `ln(12/γ)` is undefined at γ → 0, and resorption excludes baseline turnover.
6. **Realized-replacement acclimation:** traits are updated before this step's allocation
   (`:103-109`), carbon-starved cohorts would never acclimate, and only SLA is fixed at leaf
   expansion (Vcmax/Rd acclimate within leaves over weeks).
7. **Plumbing** the plan omitted: cohort-state allocate/permute/fuse/split, restart I/O with
   `llspan` conversion, C-API structs, required `[phenology]` keys, other consumers of `llspan`.

### 8.1 The revised design

| | v2 |
|---|---|
| Senescence budget | Cohort state **B** ∈ [0, f], in the carbon layer: realized active shed draws it down (`active_shed_c / leaf_target`), realized leaf growth refills it, capped at f = `senescence_fraction`. f = 1 bypasses it (bit-identical). Exact at any `dt_slow`; no episode boundary. |
| `TEMP` shed window | Only while days shorten (`doy_effective` ≥ 172), beside the kernel's existing GDD/chilling calendar gates. Autumn senescence answers to shortening days; without the window the unspent budget is shed in early spring. |
| Background loss | `leaf_turnover_rate · flush_drive · leaf`: old leaves are replaced while new ones grow. No temperature function, no `evergreen` flag. An empty flush mask (drive ≡ 1) keeps today's constant turnover. |
| Combination | Active shed and background loss **add** (independent processes). |

### 8.2 Evidence (FI-Hyy, soil temperature drives the cold-drop cue)

| Design | r | RMSE | Jan–Apr | May–Aug | Sep | Oct | Nov | Needle life |
|---|---|---|---|---|---|---|---|---|
| Observed | — | — | 0.06–0.21 | 0.20–1.06 | 3.70 | 3.52 | 0.27 | — |
| §2.2 episode budget | fails: no shedding after the first episode | | | | | | | |
| Carbon budget, all-year `TEMP` shed, cold-suppressed background | 0.70 | 1.39 | 0.01–2.28 (April peak) | 0.25–0.67 | 3.78 | 2.64 | 0.06 | 2.0 yr |
| Carbon budget + autumn window | 0.78 | 1.31 | 0.00 | 0.00–2.18 | 5.45 | 3.30 | 0.00 | 3.7 yr |
| + flat background (0.19 yr⁻¹) | 0.75 | 1.28 | 0.45–0.48 | 0.55–0.72 | 1.90 | 3.62 | 0.47 | 2.9 yr |
| + cold-suppressed background (0.25 yr⁻¹) | 0.80 | 1.18 | 0.04–0.27 | 0.52–0.87 | 4.19 | 3.20 | 0.18 | 3.0 yr |
| **+ flush-coupled background (0.27 yr⁻¹)** | **0.80** | **1.18** | 0.00–0.07 | 0.50–1.12 | 4.41 | 3.11 | 0.00 | 3.5 yr |

The flush-coupled background fits as well as the temperature response it replaces and reproduces
the June litter (0.91 against 0.94). Fitted: the shed can start from 2 September (day length < 13.9 h),
`k_shed_max` 0.005 d⁻¹, f = 0.12, `leaf_turnover_rate` 0.27 yr⁻¹. All rows are offline fits of
monthly needle-fall shares with annual loss held near 30 %; none has run in the coupled model.

### 8.3 Open for the owner

1. **Resorption of background loss.** If background loss now includes replacement senescence of
   old leaves, should `retained_carbon_fraction` apply to it?
2. **Plasticity.** Keep the Lloyd slope (on realized lifespan) or make `k_γ` an explicit parameter
   (default 0); acclimate SLA at realized replacement and Vcmax/Rd with a time constant?
3. **Realized lifespan** as an output computed from a cohort leaf-shed diagnostic, rather than two
   new running means of state.
4. **Census tool:** steady-state leaf litter from realized annual loss, net of resorption, taken
   from a spin-up or declared, since it cannot be computed offline for `HYDRO` PFTs.
5. **Sum** moves the light-exchanger and the shipped deciduous configs (bit-identical where the
   shed mask is empty).

### 8.4 Leaf cover instead of a budget state (owner, 2026-10-05)

The owner proposed tracking leaf cover (leaf carbon / full-canopy leaf carbon) instead of a separate
budget: senescence stops when leaf cover reaches a floor, `1 − senescence_fraction`. No new cohort
state is needed. Offline at FI-Hyy (soil T, flush-coupled background) it fits exactly as well as
the budget state: with the autumn window r 0.80, RMSE 1.18, floor 0.85, needle life 2.9 yr; without
the window r 0.71 and the April peak returns (warm April days start a small flush that short days
and cold soil then shed), so the window stays. Naming of states and parameters is under discussion.

### 8.5 Owner naming decisions (2026-10-05)

States: `leaf_flush_tendency`, `leaf_shed_tendency` (smoothed, 0–1), `leaf_cover` (leaf carbon / full
leaf carbon), `leaf_flush_rate` and `leaf_shed_rate` (realized carbon fluxes; resorption applies to
the shed); cue memory `growing_degree_days`, `chilling_days`, `dry_days`, `wet_days`,
`soil_water_mean`, `shortwave_mean`. Parameters: `flush_cue_timescale`, `shed_cue_timescale`,
`flush_rate_max`, `shed_rate_max`, `leaf_turnover_rate_toc` (top-of-canopy, subject to plasticity),
`min_leaf_cover` proposed for the senescence floor.

### 8.6 A symmetric temperature cue (tested 2026-10-05)

Today the `TEMP` flush integrates warmth (degree days vs a chilling-dependent threshold) while the
`TEMP` shed is instantaneous (short days with cool soil, or cold soil). The symmetric alternative
integrates on both sides, air temperature only: warmth `W = Σ max(0, T − T_flush)` and cold
`C = Σ max(0, T_shed − T)` from midsummer (Delpierre et al. 2009 style), both reset at midwinter;
flush signal = σ(W vs threshold) × (1 − σ(C vs `C*`)), shed signal = σ(C vs `C*`).

| Site | Structure | Fit | Autumn date, between years |
|---|---|---|---|
| Harvard | current | RMSE 0.186 MODIS / 0.065 leaf fall / 0.073 baskets | r 0.51, RMSE 2.7 d |
| Harvard | symmetric | RMSE 0.194 / 0.085 / 0.122 | r 0.59, RMSE 3.8 d |
| Hyytiälä | current + autumn window (soil T) | r 0.80, RMSE 1.18 | — |
| Hyytiälä | symmetric (air T) | r 0.78, RMSE 1.22 | — |

Comparable, not better. It drops five parameters for two (`cold_drop_daylength`,
`cold_drop_soiltemp1/2`, two widths → `T_shed`, `C*`) and soil temperature as a phenology driver.
Fitted `T_shed` is 16.9 °C at Harvard and 23.6 °C at Hyytiälä, where the cold sum acts almost as a
clock from midsummer.

### 8.7 Chilling and light control in the symmetric cue (tested 2026-10-05)

Offline fits of the symmetric temperature cue (§8.6), best of three optimizer seeds each
(`pheno_variants.py`). Light variants weight the sums: warmth counts under long/bright conditions,
cold under short/dark ones, through a logistic in day length D or 10-day mean shortwave R.

| Variant | Extra params | Harvard leaf fall / baskets / MODIS | Harvard spring, autumn timing | Harvard AIC | Hyytiälä r | Hyytiälä AIC |
|---|---|---|---|---|---|---|
| with chilling | b | 0.087 / 0.114 / 0.195 | 7.2 d, 3.8 d | −7258 | 0.790 | 43.4 |
| no chilling | — | 0.087 / 0.115 / 0.197 | 8.5 d, 3.8 d | −7223 | 0.788 | 41.8 |
| P1: one day length D₀ | D₀ | **0.072 / 0.101 / 0.188** | **5.8 d, 2.3 d** | **−7502** | 0.819 | 32.9 |
| P2: separate day lengths | D_f, D_s | 0.071 / 0.102 / 0.188 | 5.7 d, 2.4 d | −7506 | 0.822 | 34.4 |
| R1: one radiation level R₀ | R₀ | 0.079 / 0.100 / 0.199 | 10.3 d, 3.0 d | −7255 | **0.877** | **7.3** |
| H1: D flush, R shed | D_f, R_s | 0.077 / 0.102 / 0.192 | 7.1 d, 2.9 d | −7389 | 0.877 | 9.3 |
| H2: R flush, D shed | R_f, D_s | 0.071 / 0.098 / 0.191 | 7.7 d, 2.3 d | −7456 | 0.822 | 34.5 |

Chilling changes neither site. One shared light threshold gives most of the gain; a second adds
nothing (Harvard's P2 day lengths came out equal, 12.4 and 12.3 h). The gain is on the shed side at
both sites, through day length at Harvard (D₀ = 12.3 h, about the equinox) and through radiation at
Hyytiälä (R₀ = 131 W m⁻², reached in late summer). `T_shed` ran to its 27 °C bound in several fits, so the cold
sum behaves as a light-started clock accelerated by cold. AIC is indicative only (residuals are
autocorrelated); Hyytiälä has six years.

### 8.8 What ends the flush (tested 2026-10-05)

The warmth sum only grows from midwinter, so once past W* the flush signal would stay on into
autumn and refill a shedding canopy; today's kernel stops it with a hard-coded reset of the warmth
sum on 31 August. The symmetric cue ends it with `(1 − σ(C − C*))`, i.e. **senescence ends the
flush**. The alternative that keeps the flush free of the cold sum, `flush = σ(W − W*) · g(L)`
(long days end the flush, no new parameter), was fitted at both sites:

| | Harvard leaf fall / baskets / MODIS | spring, autumn timing | AIC | Hyytiälä r | AIC |
|---|---|---|---|---|---|
| senescence ends the flush (P1) | 0.072 / 0.101 / 0.188 | 5.8 d, 2.3 d | −7502 | **0.819** | **32.9** |
| long days end the flush (P1g) | 0.073 / **0.089** / **0.184** | **5.3 d**, 2.7 d | **−7584** | 0.797 | 40.9 |

Harvard prefers the day-length end; Hyytiälä prefers the senescence end, because one shared L₀
cannot place both the end of the pine's flush (and its flush-coupled summer litter) and the start
of cold counting: with the flush gate the fit moved L₀ to 12.2 h and shifted the needle-fall pulse
from September to October.

### 8.9 Turning the flush down: dormancy literature and tests (2026-10-05)

Literature (notes and Crossref-checked references in
`~/claude_workspace/meds_pheno_proto/dormancy_literature_review.md`): growth cessation and bud set
come weeks before leaf senescence and have their own cue (Fracheboud et al. 2009); the cue is a
critical day length that lengthens with latitude of origin in most temperate/boreal trees (Böhlenius
et al. 2006; Heide 1974), low temperature in Rosaceae (Heide & Prestrud 2005), and completion of
temperature-driven preformed shoot growth in pines (Ekberg et al. 1979; Schiestl-Aalto et al. 2015).
Budburst follows an alternating chilling–forcing rule, heat requirement `a + b·exp(c·chill)`
(Murray et al. 1989), which matters under warm winters (Fu et al. 2015) but cannot be constrained
at cold-winter sites, so its parameters should come from experiments (Chuine et al. 2016).

Tested ways to end the flush (best of three seeds; shed side as in §8.7 P1):

| Flush ends by | Extra | Harvard MODIS / fall / baskets | AIC | spring, autumn | Hyytiälä r | June share (obs 0.94) |
|---|---|---|---|---|---|---|
| senescence (P1) | — | 0.188 / 0.072 / 0.101 | −7502 | 5.8, 2.3 d | 0.819 | 0.75 |
| shared day length (P1g) | — | 0.184 / 0.073 / 0.089 | −7584 | 5.3, 2.7 d | 0.797 | 1.26 |
| own day length L_gc (P1c) | L_gc | 0.171 / 0.064 / 0.077 | −7926 | 11.8, 2.4 d | 0.836 | 1.42 |
| + chilling-released dormancy (P1d) | L_gc, chill | 0.171 / 0.064 / 0.075 | −7927 | 11.0, 2.4 d | 0.835 | 1.38 |
| **warmth window ΔW (P1w)** | ΔW | **0.154 / 0.062 / 0.079** | **−8305** | 7.3, 2.2 d | 0.818 | **1.16** |
| warmth window + Botta chilling fixed (P1wc) | ΔW | 0.153 / 0.062 / 0.075 | −8334 | 6.8, 2.3 d | 0.818 | 1.08 |

The warmth window, `flush = σ(W − W*) · (1 − σ(W − W* − ΔW))` (growth stops once a heat sum ΔW
past budburst is complete), fits Harvard best and Hyytiälä as well as P1, with the needle-fall
shape right: the flush-coupled background makes the June litter (old needles shed during shoot
growth). Fitted ΔW: 274 K·d at Harvard, 174 K·d at Hyytiälä (flush window late May–late June).
Chilling at Botta's fixed values costs nothing at either site. P1c's day-length cessation fitted
16.4 h at Harvard (longer than any day there: the flush stops at midsummer) and 11.3 h at Hyytiälä
(late September, too late for pine bud set). Chilling never limits at either site (P1d fitted
14 and 1.4 chill days).

### 8.10 Separate photoperiod gates (owner's preference, tested 2026-10-05)

P1c (§8.9) was not a fair test of a photoperiod gate: its growth-cessation gate was a hard on/off
switch applied only after midsummer, and with background loss tied to the flush the gate also set
how long summer litter continued. At Harvard it fitted 16.4 h, beyond the longest day, so the flush
stopped at midsummer and the optimizer slowed the flush to its bound (spring error 11.8 d); at
Hyytiälä it fitted 11.3 h (late September), the long flush made background loss carry the litter
and senescence took only 6 % (September 1.73 against 3.70).

G2 uses smooth, symmetric gates on the signals, sums unweighted:
`flush = σ(W − W*) · σ(D − L_flush)`, `shed = σ(C − C*) · σ(L_shed − D)`.

| | Harvard MODIS / fall / baskets | AIC | spring, autumn | Hyytiälä r | AIC |
|---|---|---|---|---|---|
| P1 | 0.188 / 0.072 / 0.101 | −7502 | 5.8, 2.3 d | 0.819 | 32.9 |
| P1w (warmth window) | 0.154 / 0.062 / 0.079 | **−8305** | 7.3, 2.2 d | 0.818 | 35.6 |
| **G2** | 0.184 / 0.061 / 0.079 | −7693 | 6.1, 2.8 d | **0.850** | **22.5** |
| G2, Harvard gates within 9–15.2 h | 0.185 / 0.061 / 0.079 | −7674 | 5.9, 2.8 d | — | — |
| G2 with constant background | 0.184 / 0.061 / 0.079 | −7693 | 5.9, 2.8 d | 0.823 | 33.7 |

At Hyytiälä the gates are physiologically plausible: flush from 20 May to 24 July (17.65 h, the
pine's shoot-growth period), senescence from 20 August (15.16 h); September 3.77 and October 3.61
against 3.70 and 3.52. At Harvard the data do not identify thresholds: the flush gate sits at the
longest day (15.2 h) and the shed gate at 9.2 h (late November), so both act as graded day-length
multipliers and the temperature sums set the timing. Background loss tied to the flush matters at
Hyytiälä (constant background: r 0.823, flat winter litter).

### 8.11 G2 adopted; the flush gate applies all year (owner, 2026-10-05)

G2 as tested (§8.10) already gates the flush all year. Against the same gate applied only after
midsummer (G2s), and a Harvard fit with the flush threshold held to ≤ 14 h (G2h):

| | Harvard MODIS / fall / baskets | AIC | spring r, error | Hyytiälä r | AIC |
|---|---|---|---|---|---|
| **G2, gate all year** | 0.184 / 0.061 / 0.079 | −7693 | 0.09, 6.1 d | **0.850** | **22.5** |
| G2s, gate after midsummer | 0.174 / 0.062 / 0.074 | −7887 | 0.34, 11.9 d | 0.812 | 37.8 |
| G2h, flush gate ≤ 14 h | 0.185 / 0.061 / 0.080 | −7672 | 0.29, 8.0 d | — | — |

- Hyytiälä: spring flushing is limited by warmth (W* 177 K·d met about 28 May; the gate opens
  20 May), the gate ends the flush on 24 July and allows senescence from 20 August.
- Harvard does not identify the flush threshold: at 16.9 h (G2) the warmth requirement is tiny
  (23 K·d, met in March) and the gate (never above 0.10) sets spring timing, identical every year
  (r 0.09); held to ≤ 14 h the fit takes 10.35 h, the gate stops binding, and warmth sets spring
  (92 K·d, about 24 April; r 0.29) at the same overall fit. G2s repeats P1c's midsummer stop.
- A year-round gate ties the spring opening and the late-summer closing symmetrically about the
  solstice (14 h at 42.5° N: open 4 May, close 9 August). Where the data cannot pin it, the flush
  threshold should come from the literature for the PFT and latitude.

## 9. Current recommendation for temperate forests (2026-10-05)

Supersedes the temperature cue of §2.1 and the senescence limit of §2.2/§8.1. Offline evidence
only (Harvard Forest, Hyytiälä); none of it has run in the coupled model.

### 9.1 Cue equations (daily, per cohort)

```
W  += max(0, T_air − T_flush)                    from midwinter; reset at midwinter
C  += max(0, T_shed − T_air)                     from midsummer; reset at midwinter
W* = flush_degree_days                           a constant: no chilling (owner, 2026-10-05)

flush signal = σ(s_W (W − W*)) · σ(s_Lf (D − L_flush))   warmth reached AND long days
shed signal  = σ(s_C (C − C*)) · σ(s_Ls (D − L_shed))    cold reached AND short days (s_Ls < 0)

leaf_flush_tendency += (flush signal − tendency) · dt / flush_cue_timescale
leaf_shed_tendency  += (shed signal  − tendency) · dt / shed_cue_timescale
```

σ(z) = 1/(1 + e^−z). Each switch has one **integrated sharpness** s (units 1/[x]; its sign gives the
direction; the 0.12–0.88 transition spans x* ± 2/|s|), no separate width and no global
`cue_sharpness` (owner, 2026-10-05). Values used in the fits: s_W = 0.04 (K·d)⁻¹, s_C = 0.1 (K·d)⁻¹,
s_Lf = 1 h⁻¹, s_Ls = −1 h⁻¹ (the kernel's `photo_slope` default is 2 h⁻¹). D is day length [h];
midwinter/midsummer are hemisphere-aware
(`doy_effective`). No soil temperature; air temperature only.

### 9.2 Carbon layer

```
senescence  = shed_rate_max · leaf_shed_tendency · leaf     stops when leaf_cover reaches min_leaf_cover
background  = leaf_turnover_rate · leaf_flush_tendency · leaf
leaf_shed_rate  = senescence + background                    (realized; to litter, net of resorption)
leaf_flush_rate = min(full − post-shed leaf, flush_rate_max · leaf_flush_tendency · full · dt)
                                                              (realized; carbon permitting)
```

A dormant canopy (flush tendency ~0) below `bare_leaf_cover` snaps to bare.

### 9.3 Parameters

| Parameter | Meaning | Harvard (deciduous) | Hyytiälä (pine) | Source |
|---|---|---|---|---|
| `flush_base_temp` | warmth base | 5 °C | 5 °C | fixed (literature) |
| `flush_degree_days` (W*) | warmth requirement for flushing | 92 K·d* | 177 K·d | fit |
| `critical_daylength_flush` | days must be longer to flush | not identified (10–17 h)* | 17.65 h (20 May–24 Jul) | literature where data cannot pin it |
| `shed_base_temp` | cold base | 17.1 °C | 6.7 °C | fit |
| `shed_degree_days` | cold requirement after midsummer | 46–48 K·d | 32 K·d | fit |
| `critical_daylength_shed` | days must be shorter to shed | 9.2–9.8 h (graded) | 15.16 h (from 20 Aug) | fit / literature |
| `flush_rate_max` | max flush [d⁻¹] | 0.02–0.30* | 0.18 | fit |
| `shed_rate_max` | max senescence [d⁻¹] | 0.32 | 0.13 | fit |
| `flush_cue_timescale`, `shed_cue_timescale` | smoothing | 5 d | 5 d | fixed |
| `leaf_turnover_rate_toc` | background loss while flushing [yr⁻¹] | ~0 | 0.42 | fit |
| `min_leaf_cover` | senescence floor | 0 | 0.82 | fit / observable |

\* Harvard trades the flush day length against the warmth requirement and flush rate: 16.9 h
gives W* = 23 K·d (spring set by day length, r 0.09); ≤ 14 h gives 10.35 h, W* = 92 K·d (spring
set by warmth, r 0.29); same overall fit. The sharpnesses s are global constants per cue type (§9.5).

### 9.4 Performance and what remains uncertain

Harvard: RMSE 0.184–0.185 MODIS, 0.061 leaf fall, 0.079–0.080 baskets; autumn timing 2.8 d
(r 0.46–0.48), spring 6.1–8.0 d. Hyytiälä: needle-fall shares r 0.85, September 3.77 and October
3.61 against 3.70 and 3.52. Uncertain: the Harvard flush day length (set it from the literature);
year-to-year autumn skill (r 0.46 against 0.64 for P1); two sites; resorption of background loss
(§8.3). Chilling is left out by the owner's decision: neither site constrains it (§8.7), but the
literature argues it matters where winters are warm enough to leave chilling short (Fu et al. 2015;
Chuine et al. 2016, §8.9) — revisit if MEDS is applied there.

### 9.5 Switch widths: fixed, not fitted per site (tested 2026-10-05)

Each switch σ(k(x − x*)/w) depends only on k/w, so `cue_sharpness` and a width are one parameter:
each switch has a threshold and one width (today: warmth 25 K·d = 50/2, cold 10 K·d, day-length
gates 1 h). Fitting all four widths per site (G2w):

| | params | fit | AIC | seeds agree |
|---|---|---|---|---|
| Harvard G2 | 8 | 0.184 / 0.061 / 0.079 | −7693 | yes |
| Harvard G2w | 12 | 0.187 / 0.063 / 0.070 | −7625 | no (warmth width 45–135 K·d, cold 4.5–48 K·d, W* 26–138 K·d at losses 0.0437–0.0449) |
| Hyytiälä G2 | 9 | r 0.850 | 22.5 | — |
| Hyytiälä G2w | 13 | r 0.845 | 32.5 | no (best seed does not reach the fixed-width optimum) |

Neither site constrains the widths. **Decided (owner): one integrated sharpness per switch,
σ(s(x − x*)), replacing both the width and `cue_sharpness`.** Recommendation: a global constant set
from physical reasoning (stand-level spread of budburst and senescence, a day-length gate sharp
enough to act as a threshold — the kernel's `photo_slope` 2 h⁻¹ is half the prototypes' width), or
estimated jointly across many sites; drop `cue_sharpness`. The 5-day cue smoothing also spreads
each transition and partly overlaps the widths.

## 10. The water cue (owner proposal, 2026-10-05)

Accumulate predawn leaf water potential against the turgor-loss point, the water analogue of the
warmth and cold sums:

```
M_wet += max(0, ψ_pd − ψ_tlp)        flush signal = σ(s_wet (M_wet − M_wet*))
M_dry += max(0, ψ_tlp − ψ_pd)        shed signal  = σ(s_dry (M_dry − M_dry*))
```

Water has no calendar, so each sum needs a reset rule. Each resets when the opposite switch
crosses 0.5 upward (M_dry when the wet switch does, M_wet when the dry switch does): the crossing
event, not the level, or a long wet season would hold the wet switch high and wipe every drought.
A brief rain then adds a little wet credit without wiping the drought, unlike the old
consecutive-day counters, which reset on a single opposite day. ψ_tlp can be fixed at the PFT's measured
turgor-loss point, leaving the two centres to fit.

**Not fitted offline.** The predawn ψ from the Palo Verde test run (§8, `pv_daily.csv`) has no
shed → ψ feedback: the placeholder BCI PFT has no `[phenology]` section, so its LAI stays 5.2–5.4
all year while ψ falls to −3.1 MPa in March. A drought-deciduous canopy that sheds transpires
less and its predawn ψ recovers toward the soil's, so fitting the cue to that driver would be
fitting to the wrong water potential. The water cue has to be tested in the coupled model, with the
cue on and a dry-forest stand (the Guanacaste example), after implementation.

## 11. Implementation (2026-10-06, `dev/example02-canopy-phenology`)

The generic model of §9–§10 replaces the old kernel, its config keys and its Python front end. What
the implementation added or found:

1. **The compiled kernel reproduces the emulation.** With the §9 values (G2h, G2) the example scores
   exactly what the Python prototype did: Harvard MODIS / leaf-fall / basket RMSE 0.185 / 0.061 /
   0.080; Hyytiälä share RMSE 1.031, r 0.850, annual needle fall 0.379.
2. **The water sums reset on the crossing.** §10 as first written ("resets when the opposite switch
   passes 0.5") was implemented as a level test, which deadlocks: a wet season holds the wet switch
   on and wipes the dry sum every day. The kernel resets a sum when the other switch crosses 0.5
   upward; `test_plant_phenology` covers the wet-season and brief-rain cases.
3. **A gradual flush gate churns leaves in autumn.** With the G2h values the Harvard flush tendency
   is still 0.64 in October (gate σ(D − 10.35), D ≈ 11 h): the canopy refills while it senesces and
   drops 1.89 canopies of leaves a year (1.13 in October alone). The basket target, normalised
   within each year, cannot see it. Adding the annual leaf fall to the loss (one canopy, weighted as
   Hyytiälä's 0.30):
   - with the gate sharpness fixed at 1 h⁻¹ the fit pushes the flush day length to its bound
     (14.1 h) and W* to 6 K·d: litter 1.04, but spring timing collapses (r 0.05, RMSE 18.7 d);
   - with the sharpness fitted it runs to the 8 h⁻¹ bound (threshold 13.75 h): litter 1.01, spring
     RMSE 5.5 d (r 0.34), autumn 2.9 d (r 0.44), loss 0.0436, below G2h's 0.0446 without the term.
     Hyytiälä fitted the same way picks 6.6 h⁻¹ at 13.1 h (share r 0.853, needle fall 0.304).
   The example fits `flush_light_sharpness` at both sites; the kernel default stays 1.
4. **The dormant snap needs the sharp gate.** With a gradual gate the flush tendency never falls
   below the 10⁻⁶ day⁻¹ dormancy rate, so a deciduous canopy sits at a small residual cover instead
   of going bare; with the step gate it goes bare.
5. **Deferred:** renaming `leaf_lifespan_toc` to `leaf_turnover_rate_toc` (D2); the water cue on
   the model's own predawn leaf potential (Guanacaste, with the shed-to-water feedback); scoring a
   coupled run's leaf cycle. The biophysics example's PFT file takes the Harvard values, but its
   committed outputs predate the scheme. (The per-cohort light and the tropical examples: §12.)

## 12. Day length or PAR, and the tropical examples (owner, 2026-10-06)

The owner asked whether one light variable could serve all habits (PAR in the temperate fits, or day
length in the tropics) and, if not, to separate the photoperiod and PAR controls; and for a tropical
example at BCI with a soil-moisture surrogate for predawn water potential.

**Test** (`~/claude_workspace/meds_pheno_proto/light_test.py`, a Python emulator of the kernel checked
against the compiled one to 1e-16, so each side can read either variable). Same loss and data as the
example; the PAR fits carry the averaging window as one more parameter.

| Light cue, flush / senescence | Harvard | Hyytiälä | BCI exchanger, senescence only | BCI exchanger, both sides |
|---|---|---|---|---|
| day length / day length | **0.0429** | 1.03 (example fit) | 0.138 | 0.101 |
| PAR / PAR | 0.0524 | **0.699** | **0.067** | **0.065** |
| day length / PAR | 0.0491 | 0.750 | — | — |
| PAR / day length | 0.0461 | 0.810 | — | — |

- Harvard: a PAR flush gate nearly doubles the spring timing error (RMSE 5.5 → 9.6 days, r 0.31 →
  0.05); PAR on the senescence side costs 0.6–0.9 days in autumn.
- Hyytiälä: PAR raises the needle-fall correlation from 0.85 to 0.91.
- BCI (GLiMP 2003–2019, Lutz shortwave × 2.115 as PAR, canopy floor ≥ 0.8): day length reaches 0.051
  without the floor only by stripping an evergreen canopy to 27 %; with it, PAR wins and follows the
  year-to-year dry-season litter (r 0.59 against 0.08).
- No mixed pair beats the best pure one at any site.

**Decision: separate cues.** The light cue becomes two: DAYLENGTH (bit 2) and PAR (bit 8; WATER stays
4). PAR is the running mean of the daily PAR reaching the cohort's top: the two-stream solver now
returns the beam plus downward diffuse at each layer's top (`incid_top`), the fast loop sums it per
cohort (`par_accum`, the gpp_accum pattern, no reduction), and the driver divides by the step. So the
cue carries the canopy light gradient; the site shortwave reduction for phenology is gone. Example02
moves Hyytiälä to PAR (refit: r 0.90, 0.30 canopy a year).

**BCI example.** Drivers from the BCI tower (July 2012 – August 2017): PAR = 2.11 × shortwave;
predawn potential surrogate ψ = −exp(−3.74) (SWC/0.829)^−2.58 MPa, a Campbell curve fitted to Kupers
et al.'s 1,020 paired samples (r −0.42; spatial scatter dominates) with the tower/plot moisture ratio
on their four dates as the bulk density. The surrogate spans −0.05 to −0.9 MPa, so the water
thresholds are on its scale, not a leaf TLP. Two species, the same masks (flush WATER; senescence
PAR + WATER):

- light exchanger, fitted to GLiMP (share RMSE 0.19, r 0.91, one canopy a year, after Leigh 1999):
  water threshold −1.5 MPa (never reached), senescence when the 3.7-day PAR mean passes 445; cover
  stays 1, half the year's leaf loss in January–April;
- drought-deciduous, set by hand: threshold −0.33 MPa, senescence after 1 MPa day of drought
  (sharpness 20), reflush after 3 MPa day of wetting, flush smoothing 2 days, no PAR response;
  leafless 29–100 days each dry season, 1.16 canopies a year. Sharpness 5 left a 0.7 % senescence
  floor that, refilled by the flush, cost 0.25 canopy a year; 5-day flush smoothing overlapped the
  drop by another 0.1.

## 13. One light cue: the hours of light (owner, 2026-10-06)

**Proposal (owner).** One light variable for every PFT: the hours a day the PAR exceeds a per-PFT
`par_min`, like a degree-day sum for light. A low `par_min` makes it the photoperiod; a high one, the
bright hours. Temperate PFTs flush as the hours rise and shed as they fall; a tropical light
exchanger sheds (and flushes) as they rise in the bright dry season.

**Tests (emulator; `~/claude_workspace/meds_pheno_proto/parhours_test.py`).** On tower PAR, with
`par_min` fitted, the hours of light matched the better of the two §12 cues at each site: Harvard
0.046 (`par_min` 33; day length 0.0445, PAR 0.052), Hyytiälä 0.59 (280; PAR 0.70), BCI 0.067 (700;
PAR 0.065). On ERA5-Land PAR through MEDS's reader, thresholds fitted on tower PAR do not carry over
(ERA5's PAR is 9–16 % brighter, mostly MEDS's SW→PAR split, #369); refitted on the forcing the model
runs with, Harvard matched the tower fit (0.046) and the Hyytiälä fit agreed across drivers once the
optimizer converged. A free `par_min` at BCI drifted to the photoperiod, a calendar; tropical fits
need `par_min` bounded to bright hours. The forcing reader dropped the shortwave of a sunrise or
sunset hour whose midpoint was dark, which biased low-`par_min` hours by up to −0.65 h a day (#371,
fixed in #372).

**Decision.** Day length and PAR (§12) are replaced by one cue, LIGHT (bit 2): `par_min`,
`flush_light_hours`, `flush_light_sharpness`, `shed_light_hours`, `shed_light_sharpness`,
`light_window`. The fast loop counts the hours per cohort (`light_hours_accum`: each fast step whose
PAR at the cohort's top, in photons, exceeds `par_min` adds its length), the driver divides by the
step, and the kernel keeps a running mean (`light_hours_mean`). A cohort with no light memory (a
cold start, a recruit) starts the mean from its first day, so its first days do not read as
darkness; fusion keeps the survivor's memory, as for the other cues.

**Example02 (owner).** Four sites, one model, all on ERA5-Land PAR-hours through MEDS's reader:
Harvard Forest and Hyytiälä (TEMP + LIGHT), the BCI light exchanger (WATER + LIGHT; BCI no longer
carries a drought-deciduous species) and Palo Verde, Costa Rica (drought-deciduous; WATER + LIGHT),
driven by a canopy predawn leaf ψ from a MEDS run of the BCI 2010 census at Palo Verde on ERA5-Land
and scored against the Xu et al. (2016) litter traps (g per 0.25 m² trap) and MODIS LAI, 2009–2013.
Fits use bounded sharpness and windows and four optimizer seeds, each with a local polish (the
example README); leave-one-year-out checks against a climatology were run in the prototypes only,
and the owner accepted fits of the mean seasonal cycle.

**Result (2026-10-06).** Best losses, with the four seeds' spread: Harvard 0.04392 (to 0.04393),
Hyytiälä 0.00397 (to 0.00429), BCI 0.00056 (to 0.00075), Palo Verde 0.0234 (to 0.0369). `par_min`
came out at 2 (Harvard, held to 1–100), 99 (Hyytiälä, 50–600), 1081 (BCI, 100–1200) and 5 (Palo
Verde, free): the deciduous canopies count the photoperiod, the evergreen ones the bright hours. At
Palo Verde, under the same objective, water alone reaches 0.064 (leaves fall in December–January,
before the stand-in ψ drops below −3 MPa in mid-February to mid-March), the photoperiod alone 0.025
and both 0.023, on the monthly leaf litter.

## References

Botta et al. (2000) *Glob. Change Biol.*; White et al. (1997) *Glob. Biogeochem. Cycles*; Xu et
al. (2016) *New Phytol.* 212:80; Kim et al. (2012) *Glob. Change Biol.*; Lloyd et al. (2010)
*Biogeosciences* 7:1833; Finnish ICP Forests Level II litterfall (Luke); ICOS ETC L2 FI-Hyy.
