# Plant respiration

Autotrophic respiration returns a large fraction of gross photosynthesis to the atmosphere, and MEDS charges
it in **three separate places on two timescales**: leaf dark respiration inside the sub-daily gas-exchange
solver, stem and fine-root **maintenance** respiration in `meds_plant_respiration` (also sub-daily), and
**growth (construction)** respiration once a day inside the carbon allocator. This page collects the three
limbs, states the accounting rule that keeps them from being double-counted, and gives the equations for the
two non-leaf maintenance terms — ED2 forms (Chambers et al. 2004 for stem), re-referenced to the 25 °C datum.

## 1. The autotrophic budget in one place

```math
R_a \;=\; \underbrace{R_d}_{\text{leaf, sub-daily}} \;+\; \underbrace{R_{\mathrm{stem}} + R_{\mathrm{root}}}_{\text{maintenance, sub-daily}} \;+\; \underbrace{R_g}_{\text{growth, daily}} \qquad(1)
```

- **Leaf dark respiration** $R_d$ is part of the photosynthesis solution, not a separate flux: the leaf kernel
  scales a per-PFT (and plastic) $`R_{d25}`$ to leaf temperature and returns $`A = A_g - R_d`$. It is derived
  in `docs/science/leaf_gas_exchange.md` and only *placed* here.
- **Stem and fine-root maintenance** are computed every `dt_fast` from the two kernels in §2 and §3.
- **Growth respiration** $R_g$ is charged once per slow step by the allocator, on realized growth (§5).

The fast loop integrates the three sub-daily limbs onto the cohort as `leaf_resp_accum`,
`stem_resp_accum` and `root_resp_accum` $`[\mathrm{kgC\,plant^{-1}}]`$ (converting with `umol_2_kgC` at
each sub-step), and the daily step sums exactly those three into the allocator's `resp_maint`, so

```math
\mathrm{NPP} \;=\; G \;-\; \big(R_d + R_{\mathrm{stem}} + R_{\mathrm{root}}\big) \qquad(2)
```

**NPP as reported is GPP net of MAINTENANCE respiration only.** Growth respiration is deliberately not
subtracted here: it is charged *inside* the allocator, on the growth actually built, so subtracting it
again would double-count it. That one sentence is the whole accounting rule for the autotrophic budget.

## 2. Stem maintenance respiration

Woody tissue respires per unit of **stem surface area** (the ED2 / Chambers et al. 2004 form), not per
unit sapwood or structural carbon — so a cohort's stem burden follows its geometry (diameter, height,
branch area) rather than its wood pool:

```math
R_{\mathrm{stem}} \;=\; r^{\mathrm{stem}}_{25}\,10^{\,\sigma\,\mathrm{dbh}}\;\; s(T_w)\;\; A_{\mathrm{stem}},
\qquad
A_{\mathrm{stem}} \;=\; \frac{\pi\,(10^{-2}\,\mathrm{dbh})\,h \;+\; \pi\,\mathrm{WAI}/n}{f_{\mathrm{ag}}} \qquad(3)
```

with $`r^{\mathrm{stem}}_{25}`$ = `stem_resp_factor25` the per-PFT baseline
$`[\mathrm{\mu mol\,CO_2\,m^{-2}\,stem\,s^{-1}}]`$ at 25 °C, $h$ the cohort height [m], $n$ = `nplant`
$`[\mathrm{plant\,m^{-2}}]`$, and WAI the **per-ground** wood area index (hence the $`/n`$ returning the
branch term to a per-plant area). The result is per plant
$`[\mathrm{\mu mol\,CO_2\,plant^{-1}\,s^{-1}}]`$; the caller multiplies by `nplant` for the patch total.
Three details are load-bearing:

- **Temperature driver.** $`T_w`$ is the cohort's own **wood temperature** (`patch_biophys_t%wood_temp`),
  the prognostic store solved by the vegetation energy balance — not canopy air, not leaf temperature. Its
  larger heat capacity damps the diurnal cycle (`docs/science/vegetation_energy_dynamics.md`).
- **Aboveground fraction.** The cylinder-plus-branch area is the *aboveground* surface; dividing by
  $`f_{\mathrm{ag}}`$ = `aboveground_frac` scales it to the whole woody plant, so belowground structure
  respires in proportion. It is a per-PFT trait gathered per cohort — the same number demography and
  cohort fusion read — never a run-uniform constant.
- **The woody gate.** `is_woody = .false.` (a grass) returns exactly zero, whatever the baseline.

The size scaler $\sigma$ = `stem_resp_size_scaler` ships at **0** (flat in diameter) — a module default
on `wood_params_t`, not a TOML key. ED2's positive DBH exponent (~0.0041 cm⁻¹) is noisy and pulls against
the shrinking-live-fraction argument for large stems, so MEDS leaves the seam present and switched off.

## 3. Fine-root maintenance respiration

Fine roots respire **per unit fine-root carbon**, the ED2 per-`broot` form:

```math
R_{\mathrm{root}} \;=\; r^{\mathrm{root}}_{25}\;\; s(T_s)\;\; B_{\mathrm{root}} \qquad(4)
```

with $`r^{\mathrm{root}}_{25}`$ = `root_resp_factor25` $`[\mathrm{\mu mol\,CO_2\,kgC^{-1}\,s^{-1}}]`$ at
25 °C and $`B_{\mathrm{root}}`$ the cohort's fine-root carbon pool `fineroot_carbon`
$`[\mathrm{kgC\,plant^{-1}}]`$. Unlike the stem limb, this one tracks a **pool** directly: a cohort that
sheds fine roots stops paying for them the same day.

**Temperature driver.** $`T_s`$ is `soil_temp_root`, the **root-fraction-weighted mean soil temperature**
over the active layers, formed once per patch per sub-step by `root_zone_environment` from the soil
column's `root_frac` profile — the same scalar the heterotrophic-respiration kernel is handed. The
kernel itself sees one temperature, so a later per-layer sum can be added without touching the seam.
A cohort with no fine roots respires nothing.

## 4. The shared temperature response

Both maintenance limbs — and leaf $R_d$ — scale from a 25 °C baseline through the **same two functions**
in `meds_temp_response`:

```math
s_{\mathrm{arr}}(T) = \exp\!\Big[\frac{E_a}{R\,T_{25}}\Big(1 - \frac{T_{25}}{T}\Big)\Big],
\qquad
s(T) = s_{\mathrm{arr}}(T)\,\frac{f_H(T_{25})}{f_H(T)}, \quad f_H(T) = 1 + \exp\!\Big[\frac{\Delta S\,T - H_d}{R\,T}\Big] \qquad(5)
```

The plain Arrhenius form rises without bound; the **peaked** form multiplies it by a high-temperature
deactivation envelope (Medlyn et al. 2002), normalised so that $`s(T_{25}) = 1`$ — a thermal optimum with a
roll-off above it. $`T_{25}`$ = `t_ref_photo` = 298.15 K is the model's single reference, so
`*_resp_factor25` **is** the rate at 25 °C and the kernels perform no reference conversion; a parameter
quoted at another datum (ED2 references stem and root respiration at 15 °C) is converted **once, where the
PFT file is written**, never at runtime.

Leaf $R_d$ reaches these functions through the dispatcher `temp_response(form, …)`, so its form follows the
run's configured selection; stem and root call `peaked_arrhenius_scale` directly and are always peaked.
Both routes end in the same code, and that is the point: `meds_temp_response` sits in the shared library at
the root of the dependency DAG precisely so leaf, stem and root cannot drift into three near-copies of an
Arrhenius — a change to the deactivation normalisation or to the clamped `safe_exp` reaches every tissue at
once. $`E_a, H_d, \Delta S`$ live separately on `wood_params_t` and `root_params_t`, defaulted to the leaf
$R_d$ values (46390 / 200000 / 490), so the three tissues share one thermal shape but may diverge.

## 5. Growth respiration on realized growth

Building one unit of growth tissue costs $`(1+g)`$ carbon, of which the fraction $g$ =
`growth_resp_factor` is respired. The allocator applies that cost **inside** each funding step, per pool,
so the charge is

```math
R_g \;=\; g\,\big(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}}\big) \qquad(6)
```

over the growth $`a_p`$ the cohort actually built this step. **Storage refill is exempt** — moving carbon
into a nonstructural pool is a transfer of sugar, not construction.

Charging on *realized* growth rather than on *available* carbon is the correct choice for three reasons.
It resolves, exactly and without iteration, the circularity that construction cost reduces the carbon
available for growth, which changes the growth. It gets the degenerate cases right: a dormant or
carbon-starved cohort that funds no growth pays no construction cost, whereas a charge on the
pre-allocation balance would bill it for tissue it never built. And it keeps carbon bound for storage out
of the charge, which a balance-based form cannot. This is a real term in the ecosystem carbon balance,
not a rounding detail: in the reference run it was roughly **23 %** of all carbon entering growth.

## 6. Where the respired carbon goes

All four limbs reach the canopy air space as part of the biotic CO₂ source driving the prognostic
`can_co2` twin, assembled in the fast pre-pass as

```math
\mathrm{NEE}_{\mathrm{biotic}} \;=\; R_d + R_{\mathrm{stem}} + R_{\mathrm{root}} + R_h \;-\; G \;+\; \dot{c}_{\mathrm{slow}} \qquad(7)
```

$`[\mathrm{\mu mol\,m^{-2}\,ground\,s^{-1}}]`$, positive toward the canopy air. The three sub-daily limbs
are exhaled in the very sub-step that computes them. Growth respiration cannot be — it is known only once a
day — so the allocator commits the per-patch total as `slow_co2_rate`, a frozen rate the *next* fast window
emits, exactly the handoff the shed-tissue water uses. That rate carries an opposite-signed starvation
correction (`growth_resp - deficit`): when storage could not fund the maintenance the fast loop had already
reported, the over-report is netted back here rather than by rewriting the previous step.

$`R_h`$ is the heterotrophic limb from the CENTURY pools — see `docs/science/soil_carbon.md`; the
canopy-air budget these fluxes enter is in `docs/science/canopy_air_space_biophysics.md`.

## 7. What is not here

- **No thermal acclimation** — no running-mean-temperature reference shift for any tissue: the 25 °C
  baselines are fixed traits, and the kernels carry no acclimation-temperature argument.
- **No storage or reproductive-tissue maintenance** — only leaf, stem and fine root carry a maintenance term;
  nonstructural carbon and seed carbon respire nothing.
- **No vertically resolved root respiration** — one root-fraction-weighted bulk soil temperature drives the
  whole fine-root pool, so warm surface roots cannot respire differently from cold deep ones.
- **No substrate or moisture limitation of maintenance respiration** — it runs at the temperature-set rate
  whatever the carbon status; §6's starvation correction is after-the-fact bookkeeping, not a rate law.
- **No size scaling of the stem baseline** in shipped runs ($\sigma = 0$, §2). All tracked in `docs/ROADMAP.md`.

## Parameters (per-PFT traits consumed)

| Symbol | Config key | Meaning |
|---|---|---|
| — | `is_woody` | woody gate; `.false.` ⇒ stem respiration is identically 0 |
| $`r^{\mathrm{stem}}_{25}`$ | `stem_resp_factor25` | stem baseline $`[\mathrm{\mu mol\,CO_2\,m^{-2}\,stem\,s^{-1}}]`$ at 25 °C (shipped 0.06) |
| $`r^{\mathrm{root}}_{25}`$ | `root_resp_factor25` | fine-root baseline $`[\mathrm{\mu mol\,CO_2\,kgC^{-1}\,s^{-1}}]`$ at 25 °C (shipped 0.30) |
| $`f_{\mathrm{ag}}`$ | `aboveground_frac` | aboveground fraction of woody carbon; scales the stem area (eq 3) |
| $g$ | `growth_resp_factor` | construction-cost fraction charged on realized growth (shipped 0.3) |
| $`R_{d25}`$ | `rd25` (plastic per cohort) | leaf dark-respiration capacity — see `plant_traits.md`, `leaf_gas_exchange.md` |
| $\sigma$ | `stem_resp_size_scaler` | DBH size effect; a `wood_params_t` default of 0, not a TOML key |
| $`E_a, H_d, \Delta S`$ | `wood_params_t` / `root_params_t` | peaked-Arrhenius parameters, defaulted to the leaf $R_d$ values |

## Where the code is

| Concept | Routine |
|---|---|
| stem maintenance respiration (eq 3) | `meds_plant_respiration`: `stem_maintenance_respiration` (`elemental pure`: scalar = one cohort, array = one patch) |
| fine-root maintenance respiration (eq 4) | `meds_plant_respiration`: `fine_root_maintenance_respiration` |
| shared temperature response (eq 5) | `meds_temp_response`: `arrhenius_scale`, `peaked_arrhenius_scale`, `temp_response` |
| leaf dark respiration $R_d$ | `meds_leaf_gas_exchange` (inside the Ci solve; see `leaf_gas_exchange.md`) |
| per-sub-step call site + patch totals | `meds_fast_prepass`: `canopy_maintenance_respiration`, `root_zone_environment` (the $`T_s`$ weighting), `column_prepass` (eq 7) |
| sub-daily → daily integration | `meds_fast_dynamics` → `leaf_resp_accum` / `stem_resp_accum` / `root_resp_accum` on the cohort SoA |
| daily maintenance sum handed to the allocator | `meds_vegetation_dynamics`: `cohort_carbon_demand` (`resp_maint`) |
| growth respiration (eq 6) | `meds_plant_carbon_allocation`: `growth_respiration`, charged inside `fill_carbon_demand` / `plant_carbon_allocation` |
| growth respiration → canopy air | `meds_vegetation_dynamics`: `compute_carbon_allocation` → `patch%slow_co2_rate` → `meds_fast_prepass` |
| per-PFT traits + TOML loading | `meds_pft_params` (`pft_table_t`), `meds_config_io` (`[pft]` block, every key required) |
| tests | `test/test_plant_respiration.f90` (grass gate, 25 °C identity, T-response, WAI branch term, `aboveground_frac` inverse scaling, size scaler, root pool scaling, per-PFT trait differentiation); growth respiration in `test/test_plant_carbon_allocation.f90` |
