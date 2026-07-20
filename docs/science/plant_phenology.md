# Leaf phenology

MEDS models leaf phenology as a **pure signal generator**: from daily environmental cues and per-PFT
traits it emits two **relative rate tendencies** — a leaf-**flush** rate and a leaf-**shed** rate, both
in $`\mathrm{day^{-1}}`$ — and nothing else. It touches no carbon and no leaf mass; the actual change in
leaf area is applied downstream by the carbon layer (`meds_plant_carbon_dynamics`, §5). The kernel's only
prognostic state is a small **phenological memory** — two smoothed governor drives plus the cue
accumulators — carried per cohort and advanced once per day.

The one generic engine covers every strategy through **two per-PFT cue masks** (which cues drive
flushing, which drive shedding); there is no `iphen_scheme`-style global switch. This generalizes ED2's
per-PFT `phenology(ipft)` habits — in particular the plant-hydraulic scheme 5, which already advances the
leaf-display fraction at per-day rates — to all habits (Botta 2000 cold-deciduous, drought-deciduous,
light-driven leaf exchange). See `docs/dev_plans/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md`.

## 1. The contract: two governors → two rates

Let $f,s\in[0,1]$ be the **flush** and **shed** governor drives (`flush_drive`, `shed_drive`) — the
prognostic memory. Each day the kernel forms a per-cue flush signal $`s^{\mathrm{fl}}_c\in[0,1]`$ and shed
signal $`s^{\mathrm{sh}}_c\in[0,1]`$ for each active cue $c$ (§2), combines them over the two masks
$`M_{\mathrm{fl}},M_{\mathrm{sh}}`$ (§3), low-passes the governors, and returns the two rates:

```math
r_{\mathrm{fl}} = k^{\max}_{\mathrm{fl}}\, f, \qquad r_{\mathrm{sh}} = k^{\max}_{\mathrm{sh}}\, s \qquad(1)
```

with $`k^{\max}_{\mathrm{fl}},k^{\max}_{\mathrm{sh}}`$ the per-PFT maximum relative rates
$`[\mathrm{day^{-1}}]`$. A rate $`1/N\ \mathrm{day^{-1}}`$ traverses a full canopy in $N$ days, so
$`k^{\max}_{\mathrm{fl}}=1/15`$ means "a bare canopy fills in ~15 days at full drive." The two boundary
values are the semantic anchors: $`r_{\mathrm{fl}}=0`$ (⇔ $f=0$) is **dormancy** — no flushing; and
$`r_{\mathrm{sh}}=0`$ (⇔ $s=0$) is **no active shedding** — the canopy still loses leaves at the separate
baseline turnover rate (§5), but phenology commands no extra senescence.

## 2. Per-cue signals

Each cue maps its (accumulated) driver to a flush and/or shed signal through the logistic
$`\sigma(z)=1/(1+e^{-z})`$; the shared dimensionless sharpness $k$ (`cue_sharpness`) sets how abruptly the
signal switches (large $k$ → an ED2-like hard threshold). The transition widths $`w_G,w_D,w_T,w_W,w_R`$
normalize each driver so $k$ is dimensionless.

**Temperature (`CUE_TEMP`).** Flushing follows a chilling-adaptive growing-degree-day sum $G$ against the
Botta (2000) threshold; shedding is the White (1997) autumn cold-drop. With $D$ the day length, $T_s$ the
shallow soil temperature, and $`\mathrm{chill}`$ the chilling-day count:

```math
s^{\mathrm{fl}}_{\mathrm{TEMP}}=\sigma\!\Big(\tfrac{k}{w_G}\big(G-G^{*}\big)\Big),\quad
G^{*}=a+b\,e^{\,c\,\mathrm{chill}};\qquad
s^{\mathrm{sh}}_{\mathrm{TEMP}}=\max\!\Big(g_D\,g_{T1},\ g_{T2}\Big) \qquad(2)
```

with $`g_D=\sigma(\tfrac{k}{w_D}(D_{\mathrm{drop}}-D))`$ (short day), $`g_{T1}=\sigma(\tfrac{k}{w_T}(T_1-T_s))`$
(cool soil) and $`g_{T2}=\sigma(\tfrac{k}{w_T}(T_2-T_s))`$ (very cold soil, unconditional). The sums are
**season-gated** (ED2 `update_thermal_sums`): in the growing half of the (hemisphere-adjusted) year, warm
days add to the GDD sum $`G \mathrel{+}= \max(0,\,T-T_{\mathrm{base}})\,dt`$ and $G$ resets outside it;
cold days in the cool half add to $`\mathrm{chill}`$. More chilling lowers $`G^{*}`$, so a colder winter
brings an earlier spring flush.

**Soil water (`CUE_WATER`).** An exponential running mean $`\bar W`$ (window $`W_{\mathrm{win}}`$) of the
available water (moisture fraction, or soil potential if `water_use_potential`) ramps between an off and an
on threshold; the gap $`(W_{\mathrm{off}},W_{\mathrm{on}})`$ is a **hold band** (hysteresis, no latch):

```math
s^{\mathrm{fl}}_{\mathrm{WATER}}=\sigma\!\Big(\tfrac{k}{w_W}\big(\bar W-W_{\mathrm{on}}\big)\Big),\qquad
s^{\mathrm{sh}}_{\mathrm{WATER}}=\sigma\!\Big(\tfrac{k}{w_W}\big(W_{\mathrm{off}}-\bar W\big)\Big) \qquad(3)
```

**Plant-hydraulic (`CUE_HYDRO`).** Consecutive-day counters on the **daily-maximum** leaf water potential
$`\psi^{\max}_{\mathrm{leaf}}`$ (ED2 `dmax_leaf_psi`), against the turgor-loss point $`\psi_{\mathrm{tlp}}`$:
a dry day ($`\psi^{\max}_{\mathrm{leaf}}<\psi_{\mathrm{tlp}}`$) increments $`n_{\mathrm{lo}}`$ and resets it
otherwise; a wet day ($`\psi^{\max}_{\mathrm{leaf}}\ge\tfrac12\psi_{\mathrm{tlp}}`$) increments
$`n_{\mathrm{hi}}`$. Their thresholds give

```math
s^{\mathrm{fl}}_{\mathrm{HYDRO}}=\mathrm{clamp}_{01}\!\big(n_{\mathrm{hi}}/n^{*}_{\mathrm{hi}}\big),\qquad
s^{\mathrm{sh}}_{\mathrm{HYDRO}}=\mathrm{clamp}_{01}\!\big(n_{\mathrm{lo}}/n^{*}_{\mathrm{lo}}\big) \qquad(4)
```

This is ED2 scheme 5 (Xu 2016), smoothed. A PFT that never droughts (its $`\psi^{\max}_{\mathrm{leaf}}`$
stays above $`\psi_{\mathrm{tlp}}`$) never sheds — it is **facultatively evergreen**.

**Photoperiod (`CUE_PHOTO`).** A day-length gate $`\sigma\!\big(m_p(D-D_c)\big)`$ (critical day length
$`D_c`$, slope $`m_p`$) **multiplies** the temperature flush signal (an extratropical spring gate).

**Light (`CUE_LIGHT`).** An exponential running mean $`\bar R`$ (window; ED2 `rad_avg`) of the incident
radiation drives an **active shed that rises with light**, while contributing $1$ (non-limiting) to the
flush side:

```math
s^{\mathrm{sh}}_{\mathrm{LIGHT}}=\sigma\!\Big(\tfrac{k}{w_R}\big(\bar R-R_{\mathrm{on}}\big)\Big),\qquad
s^{\mathrm{fl}}_{\mathrm{LIGHT}}=1 \qquad(5)
```

This is the leaf-*display* half of ED2's Kim (2012) light phenology (the leaf-*quality* half — SLA / Vcmax
/ lifespan plasticity — is a separate, deferred concern).

## 3. Combining the cues and advancing the governors

The two masks select which cues drive each side. Flushing is **conjunctive** (build the canopy only when
*every* flush cue is clear ⇒ a `min`); shedding is **disjunctive** (any stress commands senescence ⇒ a
`max`). Empty masks are the identities — permissive flush, no active shed:

```math
S_{\mathrm{fl}}=\min_{c\in M_{\mathrm{fl}}} s^{\mathrm{fl}}_c\ \ (\equiv 1\ \text{if empty}),\qquad
S_{\mathrm{sh}}=\max_{c\in M_{\mathrm{sh}}} s^{\mathrm{sh}}_c\ \ (\equiv 0\ \text{if empty}) \qquad(6)
```

Splitting the mask is what makes the four strategies (§4) fall out of parameters alone — a cue can drive
one side without the other. The governors are then a bounded low-pass (the anti-chatter memory; timescales
$`\tau_{\mathrm{fl}},\tau_{\mathrm{sh}}`$), with the guarded weight $`w=dt/\max(\tau,dt)\le1`$ so there is
no overshoot when $`\tau<dt`$:

```math
f \leftarrow \mathrm{clamp}_{01}\!\Big(f+\tfrac{dt}{\max(\tau_{\mathrm{fl}},dt)}\big(S_{\mathrm{fl}}-f\big)\Big),
\qquad
s \leftarrow \mathrm{clamp}_{01}\!\Big(s+\tfrac{dt}{\max(\tau_{\mathrm{sh}},dt)}\big(S_{\mathrm{sh}}-s\big)\Big) \qquad(7)
```

The kernel is `pure`, scalar, arithmetic-only and fixed-size (GPU/SIMD-friendly); `logistic` uses the
FPE-safe `safe_exp`, and every width is floored so no divide traps under `-fpe0`.

## 4. The four strategies (masks + parameters only)

$`M_{\mathrm{fl}}/M_{\mathrm{sh}}`$ are the flush/shed masks; $`k_{\mathrm{turn}}`$ is the baseline leaf
turnover (a *carbon-layer* trait, §5), not a phenology output.

| strategy | $`M_{\mathrm{fl}}`$ | $`M_{\mathrm{sh}}`$ | driver | canopy |
|---|---|---|---|---|
| **temperate evergreen** | `{}` (or `TEMP`) | `{}` | temperature | held full; loses leaves only at $`k_{\mathrm{turn}}`$ |
| **temperate deciduous** | `TEMP` (+`PHOTO`) | `TEMP` | temperature | spring flush, autumn cold-drop shed, bare in winter |
| **tropical drought-deciduous** | `HYDRO` | `HYDRO` | daily-max leaf ψ | full when watered, sheds under drought, reflushes on rewet |
| **tropical light-exchanging** | `{}` | `LIGHT` | radiation | stays ~full while shedding+reflushing fast under high light |

The `min`/`max` asymmetry keeps flush and shed decoupled, so the light-exchanger holds a full canopy
(permissive flush) while its shed rate swings with radiation, and a drought-deciduous PFT carrying `HYDRO`
in *both* masks stops flushing exactly when it starts shedding.

## 5. From rates to leaf area (the carbon layer)

Phenology emits only rates; `meds_plant_carbon_dynamics` turns them into the actual leaf-carbon change,
keeping the leaf-display fraction $`e=L/L_{\mathrm{full}}\in[0,1]`$ a **diagnostic** ($L$ = current leaf
carbon, $`L_{\mathrm{full}}`$ = the full-canopy allometric leaf carbon `size2leaf_carbon`). Leaf loss has
**two channels** with different owners:

- **baseline turnover** — replaceable, temperature-live (`tissue_turnover_rates`; evergreen cold-suppressed):
  $`\ell_{\mathrm{base}} = k_{\mathrm{turn}}(T)\,L\,dt`$;
- **active shed** — the phenology contribution, NON-replaceable, LINEAR in the full canopy so the full→bare
  traversal time is deterministic, clamped to the pool and snapped to bare below $`e_{\min}`$:

```math
\ell_{\mathrm{sh}} = \min\!\big(r_{\mathrm{sh}}\,L_{\mathrm{full}}\,dt,\ L-\ell_{\mathrm{base}}\big),
\qquad L \to 0 \ \text{if the step would leave } e<e_{\min} \qquad(8)
```

Flushing caps the leaf-growth demand at the flush rate (linear toward the target),
$`\Delta_{\mathrm{fl}}=\min\!\big(L_{\mathrm{full}}-L,\ r_{\mathrm{fl}}\,L_{\mathrm{full}}\,dt\big)`$, and is
funded from NPP then storage; the active shed is never refilled. The **realized leaf litter**,
$`\ell_{\mathrm{base}}+(1-\text{retained})\,\ell_{\mathrm{sh}}`$, is therefore *not* the shed tendency: a
bare deciduous canopy has a high winter shed tendency but zero litter, and a full evergreen litters via
baseline turnover with zero shed tendency (see `examples/example_phenology`). With `phenology_on=.false.`
the flush cap is non-binding and the active shed is zero, so the model is bit-identical to the
no-phenology path.

## Parameters (config names, `[phenology]` per-PFT block)

| Symbol | Config key | Meaning |
|---|---|---|
| $`M_{\mathrm{fl}}`$ | `flush_cue_mask` | OR of `CUE_*` bits driving flushing |
| $`M_{\mathrm{sh}}`$ | `shed_cue_mask` | OR of `CUE_*` bits driving shedding |
| $k$ | `cue_sharpness` | shared logistic slope [–] |
| $`k^{\max}_{\mathrm{fl}}`$, $`k^{\max}_{\mathrm{sh}}`$ | `k_flush_max`, `k_shed_max` | max relative flush / shed rate [day⁻¹] |
| $`\tau_{\mathrm{fl}}`$, $`\tau_{\mathrm{sh}}`$ | `tau_flush`, `tau_shed` | governor low-pass timescales [day] |
| $`T_{\mathrm{base}}`$ | `gdd_base_temp`, `chill_base_temp` | GDD / chilling base temperature [K] |
| $a,b,c$ | `phen_a`, `phen_b`, `phen_c` | chilling-adaptive GDD threshold $`a+b\,e^{c\,\mathrm{chill}}`$ (Botta 2000) |
| $`D_{\mathrm{drop}},T_1,T_2`$ | `cold_drop_daylength`, `cold_drop_soiltemp1/2` | autumn cold-drop thresholds (White 1997) |
| $`W_{\mathrm{off}},W_{\mathrm{on}},W_{\mathrm{win}}`$ | `water_off_threshold`, `water_on_threshold`, `water_window`, `water_width` | soil-water ramp + hold band |
| $`\psi_{\mathrm{tlp}}, n^{*}_{\mathrm{lo}}, n^{*}_{\mathrm{hi}}`$ | `leaf_psi_tlp`, `low_psi_threshold`, `high_psi_threshold` | turgor-loss point + dry/wet-day thresholds (Xu 2016) |
| $`D_c, m_p`$ | `photo_crit`, `photo_slope` | photoperiod gate |
| $`R_{\mathrm{on}}, w_R`$ | `light_on_threshold`, `light_width`, `light_window` | light-driven shed onset + running mean |
| $`k_{\mathrm{turn}}`$ | `leaf_turnover_rate` | baseline (replaceable) leaf turnover [yr⁻¹] — a *carbon* trait, §5 |

The `CUE_TEMP`/`CUE_PHOTO` cues are wired into the standalone model today; the `CUE_WATER`/`CUE_HYDRO`/
`CUE_LIGHT` drivers (soil water, daily-max leaf ψ, radiation) are threaded from the fast loop in a later
phase, so the standalone runs currently accept only `TEMP`/`PHOTO` masks (the Python
`meds.plant.pheno` example drives all five directly).

## References
- **Botta et al. (2000)**, *Glob. Change Biol.* — chilling-adaptive GDD budburst.
- **White et al. (1997)**, *Glob. Biogeochem. Cycles* — autumn day-length / soil-temperature offset.
- **Xu et al. (2016)**, *New Phytologist* — turgor-loss-point (plant-hydraulic) drought deciduousness.
- **Kim et al. (2012)**, *Glob. Change Biol.* — light phenology (the quality half is out of scope).
- ED2 `ED/src/dynamics/phenology_driv.f90`, `phenology_aux.f90`;
  `docs/dev_plans/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md`.

## Code map

| Concept | Routine |
|---|---|
| cue signals + governors → two rates | `meds_phenology`: `phenology_kernel` (`accumulate` + per-cue signals + min/max combine + low-pass), `pheno_drives_to_rates` |
| types (env / params / state / out) + cue bits | `meds_plant_types` (§ PHENOLOGY: `pheno_env_t`, `pheno_params_t`, `pheno_state_t`, `pheno_out_t`, `CUE_*`) |
| helpers | `meds_numerics`: `logistic`, `clamp01`; `meds_time`: `daylength`, `doy_effective` |
| per-cohort advance (slow loop) | `meds_vegetation_dynamics`: `advance_leaf_phenology` (folds the ED2 phenology driver), `flatten_pheno_params` |
| rates → leaf carbon (two channels + litter) | `meds_plant_carbon_dynamics`: `tissue_turnover_rates`, `active_leaf_shed`, `plant_carbon_allocation`; the seam `meds_plant_interface%get_plant_flux_slow` |
| Python front end + demo | `meds.plant.pheno` (`python/meds/plant/pheno.py`); `examples/example_phenology/` |
