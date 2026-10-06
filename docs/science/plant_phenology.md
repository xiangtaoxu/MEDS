# Leaf phenology

MEDS models leaf phenology in two layers. The **cue kernel** (`meds_phenology`) turns daily
environmental cues and per-PFT traits into two smoothed **tendencies** in $[0,1]$: a leaf **flush**
tendency $f$ and a leaf **shed** (senescence) tendency $s$. The **carbon layer** turns the tendencies
into leaf growth and leaf loss (`leaf_turnover_step`, §5). The kernel's state is the two tendencies
plus a small **cue memory** (temperature sums, water-potential sums, a shortwave running mean),
carried per cohort and advanced once per day.

One kernel covers every leaf habit. Two per-PFT **cue masks** pick which cues drive flushing and which
drive senescence, and the habit follows from the parameters: there is no evergreen flag and no global
phenology scheme switch. An evergreen is a PFT whose senescence stops at a **leaf-cover floor**
(`min_leaf_cover`), not a different model (§6).

## 1. The contract: cues → two tendencies → two potential rates

Each day the kernel (i) accumulates the cue memory (§2), (ii) passes each cue through a switch, (iii)
combines the switches into a flush signal $`S_{\mathrm{fl}}`$ and a shed signal $`S_{\mathrm{sh}}`$ (§3),
and (iv) smooths each signal into its tendency (§4). The potential relative rates are

```math
r_{\mathrm{fl}} = k_{\mathrm{fl}}\, f, \qquad r_{\mathrm{sh}} = k_{\mathrm{sh}}\, s \qquad(1)
```

with $`k_{\mathrm{fl}}`$ = `flush_rate_max` and $`k_{\mathrm{sh}}`$ = `shed_rate_max` [day⁻¹]. A rate of
$`1/N`$ day⁻¹ builds (or senesces) a full canopy in about $N$ days.

Every switch has the same form, a logistic with a **centre** $`x^{*}`$ and a **signed sharpness** $s_x$:

```math
\sigma\big(s_x\,(x - x^{*})\big), \qquad \sigma(z) = \frac{1}{1+e^{-z}} \qquad(2)
```

The switch goes from 12 % to 88 % over $`x^{*}\pm 2/|s_x|`$. Positive $s_x$ switches *on* as $x$ rises,
negative $s_x$ as it falls, so one form covers "warm enough", "short enough days" and "bright enough".
There is no separate width parameter: the sharpness is the only scale.

## 2. The cues

**Temperature (`CUE_TEMP`).** Two degree-day sums of the daily-mean air temperature $T$:

```math
W \mathrel{+}= \max(0,\ T - T_{\mathrm{fl}})\,dt, \qquad
C \mathrel{+}= \max(0,\ T_{\mathrm{sh}} - T)\,dt \ \ \text{(from midsummer on)} \qquad(3)
```

The warmth sum $W$ (`growing_degree_days`, base `flush_base_temp`) counts all year; the cold sum $C$
(`cold_degree_days`, base `shed_base_temp`) counts only once days shorten, from the summer solstice.
Both restart at midwinter. In the southern hemisphere the calendar is shifted by half a year
(`doy_effective`), so the same parameters describe the same seasons. There is no chilling
requirement: the photoperiod gate (below) does the work chilling was meant to do, holding the flush
until days are long enough, and adding chilling did not improve the fits (§7).

**Light (`CUE_LIGHT`).** The light variable $\Lambda$ is either the **day length** $D$ [h]
(`light_variable = 1`, from latitude and day of year) or a **running mean of incident shortwave**
$`\bar R`$ [W m⁻²] (`light_variable = 2`; ED2 `rad_avg`) over `light_window` days:

```math
\bar R \mathrel{+}= w\,(R - \bar R), \qquad w = \min\!\big(1,\ dt/\max(\texttt{light\_window},\,dt)\big) \qquad(4)
```

Day length is the extratropical cue (photoperiodic budburst and growth cessation); shortwave is the
tropical one (leaf exchange in the bright dry season, Kim et al. 2012).

**Water (`CUE_WATER`).** Two sums of the predawn leaf water potential $`\psi_{\mathrm{pd}}`$ against the
turgor-loss point $`\psi_{\mathrm{tlp}}`$:

```math
M_{\mathrm{wet}} \mathrel{+}= \max(0,\ \psi_{\mathrm{pd}} - \psi_{\mathrm{tlp}})\,dt, \qquad
M_{\mathrm{dry}} \mathrel{+}= \max(0,\ \psi_{\mathrm{tlp}} - \psi_{\mathrm{pd}})\,dt \qquad(5)
```

Water has no calendar, so each sum restarts when the *other* side's switch crosses 0.5 upward: the dry
sum when the wet switch turns on, the wet sum when the dry switch does. The crossing, not the level:
a long wet season would otherwise hold the wet switch on and wipe every drought. A brief rain adds
some wet credit without erasing the drought sum, unlike a consecutive-day counter. $`\psi_{\mathrm{tlp}}`$
is not a phenology key: it is derived from the same pressure–volume curve the leaf stress arrestor
uses (`pv_psi_tlp`), so the cue and the arrestor share one threshold.

## 3. Combining the switches

```math
S_{\mathrm{fl}} = \underbrace{\sigma\big(s_W (W - W^{*})\big)}_{\text{TEMP}}
\cdot \underbrace{\sigma\big(s_{\Lambda,\mathrm{fl}} (\Lambda - \Lambda^{*}_{\mathrm{fl}})\big)}_{\text{LIGHT}}
\cdot \underbrace{\sigma\big(s_{\mathrm{wet}} (M_{\mathrm{wet}} - M^{*}_{\mathrm{wet}})\big)}_{\text{WATER}}
\qquad(6)
```

```math
S_{\mathrm{sh}} = \max\Big(
\underbrace{\sigma\big(s_C (C - C^{*})\big)}_{\text{TEMP}}
\cdot \underbrace{\sigma\big(s_{\Lambda,\mathrm{sh}} (\Lambda - \Lambda^{*}_{\mathrm{sh}})\big)}_{\text{LIGHT}},\
\underbrace{\sigma\big(s_{\mathrm{dry}} (M_{\mathrm{dry}} - M^{*}_{\mathrm{dry}})\big)}_{\text{WATER}}
\Big) \qquad(7)
```

Only the cues in each side's mask enter; an absent factor is 1 in a product and absent from the max.
An empty flush mask always flushes ($`S_{\mathrm{fl}}=1`$); an empty shed mask never senesces
($`S_{\mathrm{sh}}=0`$).

- **Flushing needs every flush cue** (a product): warm enough *and* long enough days *and* wet enough.
  The day-length gate applies all year, so it both holds back budburst in a warm spell in late winter
  and turns the flush off in late summer, when the warmth sum is still high. For the second it has to
  be sharp, a photoperiod threshold (§7).
- **Senescence has two triggers** (a max): the **seasonal** one, cold *and* short days (Delpierre et
  al. 2009), and the **water** one, sustained drought on its own. A drought need not wait for autumn.

For the shed side a negative $`s_{\Lambda,\mathrm{sh}}`$ means *short* days (or dim light) trigger
senescence; a positive one means *bright* light does, as in tropical leaf exchange.

## 4. Smoothing into tendencies

```math
f \leftarrow \mathrm{clamp}_{01}\!\big(f + w_{\mathrm{fl}}\,(S_{\mathrm{fl}} - f)\big), \quad
s \leftarrow \mathrm{clamp}_{01}\!\big(s + w_{\mathrm{sh}}\,(S_{\mathrm{sh}} - s)\big), \quad
w = \frac{dt}{\max(\tau, dt)} \qquad(8)
```

with $`\tau`$ = `flush_cue_timescale`, `shed_cue_timescale` [day]. The guarded weight cannot overshoot
when $`\tau < dt`$. A cohort is born flushing and not senescing ($`f=1,\ s=0`$).

The kernel is `pure`, scalar and arithmetic-only; `logistic` uses the FPE-safe `safe_exp`.

## 5. From tendencies to leaf area: `leaf_turnover_step`

The carbon layer calls `leaf_turnover_step` once per cohort and step, *before* growth is allocated
(turnover first). With $L$ the leaf carbon, $`L_{\mathrm{full}}`$ the full-canopy (allometric) leaf
carbon and $`c = L/L_{\mathrm{full}}`$ the **leaf cover**:

```math
B = \min\!\Big(\tfrac{k_{\mathrm{turn}}}{365.2425}\, f\, L\, dt,\ L\Big), \qquad
S = \min\!\Big(k_{\mathrm{sh}}\, s\, L\, dt,\ \max\big(0,\ L - B - c_{\min} L_{\mathrm{full}}\big)\Big), \qquad
\Delta_{\mathrm{fl}} = k_{\mathrm{fl}}\, f\, L_{\mathrm{full}}\, dt \qquad(9)
```

- **Background turnover** $B$: old leaves replaced while new ones grow, at $`k_{\mathrm{turn}}`$ [yr⁻¹]
  (the reciprocal of the cohort's leaf lifespan, `leaf_lifespan_toc` with light plasticity). It is
  scaled by the flush tendency, so a dormant canopy does not turn over.
- **Senescence** $S$ at the shed potential, which stops where the leaf cover reaches
  $`c_{\min}`$ = `min_leaf_cover`. The two losses add.
- **Flush cap** $`\Delta_{\mathrm{fl}}`$: the leaf growth demand is $`\min(L_{\mathrm{full}} - L_{\mathrm{post}},\ \Delta_{\mathrm{fl}})`$
  against the post-loss pool $`L_{\mathrm{post}} = L - B - S`$, funded from NPP and then storage by the
  allocation ladder (`plant_carbon_allocation.md`). Fine roots follow at `root_to_leaf_ratio`.
- **Snap to bare**: a dormant canopy (flush rate $`k_{\mathrm{fl}} f \le 10^{-6}`$ day⁻¹) whose
  senescence would leave it below $`c_{\mathrm{bare}}`$ = `bare_leaf_cover` goes bare, when its floor
  lies below that cover. A sharp flush day-length gate (§7) turns the flush off in autumn and the
  snap acts; with a gradual one the flush tendency can stay above that rate, and a deciduous canopy
  then keeps a small residual cover, about $`k_{\mathrm{fl}} f / (k_{\mathrm{sh}} s)`$, that it keeps
  refilling and shedding.

**Resorption** (`retained_carbon_fraction`) returns a share of the senescence loss $S$ to storage; the
background turnover $B$ is not resorbed, because the turnover rate is calibrated against observed
litterfall, which already has resorption in it. The litter is $`B + (1-f_r)\,S`$. Fine roots turn over
at `fineroot_turnover_rate` with no phenology.

The same routine is exposed through the C-API (`meds_leaf_turnover_step`), so the Python front end
(`meds.plant.pheno.leaf_step`) applies exactly the carbon layer's rule.

## 6. Leaf habits as parameter sets

| habit | flush mask | shed mask | `min_leaf_cover` | what happens |
|---|---|---|---|---|
| no cues (default) | – | – | – | always flushing; loses leaves only to background turnover |
| temperate deciduous | TEMP + LIGHT | TEMP + LIGHT | 0 | flush on warmth once days are long; senesce to ~bare on cold once days are short |
| evergreen conifer | TEMP + LIGHT | TEMP + LIGHT | the share kept through winter (e.g. 0.8) | the same cues; autumn needle fall stops at the floor |
| drought deciduous | WATER | WATER | 0 | full while watered, senesces in a drought, flushes again on rewetting |
| light leaf exchange | – | LIGHT (shortwave, $`s>0`$) | – | always flushing; senescence rises with light, so leaves turn over while the canopy stays full |

Evergreen and deciduous differ in one number. The emergent leaf lifespan follows from the rates and
the floor; it is not a parameter of the kernel.

## 7. Evidence and status

The temperate scheme was fitted offline, before implementation, with a Python emulation of these
equations (`docs/dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md` §8–§9), and
`examples/example02_canopy_phenology` refits it with the compiled kernel:

- **Harvard Forest** (temperate deciduous): MODIS LAI, the HF003 leaf-fall observations and the HF069
  litter baskets.
- **Hyytiälä** (Scots pine): the ICOS monthly needle litter.

Three results shaped the design. Removing the chilling requirement changed neither fit, while light
control improved both, mostly on the senescence side. The flush day-length gate has to be sharp: a
gradual one (1 per hour) is still partly open in October, so the canopy refills while it senesces
and a deciduous stand drops about 1.9 canopies of leaves a year, which no timing observation shows.
With the annual leaf fall held to one canopy at Harvard and 0.3 at Hyytiälä, both fits make the
gate a step (sharpness 8.0 and 6.6 per hour, thresholds 13.8 and 13.1 h), and Harvard's spring
timing improves (RMSE 5.5 against 8.0 days). The other sharpnesses are not identifiable from these data (different values fit equally well), so they keep
fixed defaults (0.04 per K day for warmth, 0.1 per K day for cold, −1 per hour for the senescence
day length) and the centres are fitted.

What is **not** validated: no coupled MEDS run's leaf-area cycle has yet been scored against
observations, and the water and shortwave cues are covered by unit and driver tests only. The tropical
habits (drought-deciduous Guanacaste, light-exchanging BCI) are the next examples.

## Parameters (`[phenology]` per-PFT block)

| Symbol | Config key | Meaning |
|---|---|---|
| | `flush_cue_mask`, `shed_cue_mask` | sum of cue bits: TEMP 1, LIGHT 2, WATER 4 |
| $`\tau_{\mathrm{fl}},\tau_{\mathrm{sh}}`$ | `flush_cue_timescale`, `shed_cue_timescale` | smoothing of the tendencies [day] |
| $`k_{\mathrm{fl}},k_{\mathrm{sh}}`$ | `flush_rate_max`, `shed_rate_max` | relative rates at full tendency [day⁻¹] |
| $`c_{\min}`$ | `min_leaf_cover` | leaf cover senescence stops at [–] |
| $`c_{\mathrm{bare}}`$ | `bare_leaf_cover` | a dormant canopy below this goes bare [–] |
| $`T_{\mathrm{fl}}, W^{*}, s_W`$ | `flush_base_temp`, `flush_degree_days`, `flush_temp_sharpness` | warmth sum base [K], centre [K day], sharpness [1/(K day)] |
| $`T_{\mathrm{sh}}, C^{*}, s_C`$ | `shed_base_temp`, `shed_degree_days`, `shed_temp_sharpness` | cold sum base [K], centre [K day], sharpness [1/(K day)] |
| | `light_variable` | 1 day length [h], 2 running-mean shortwave [W m⁻²] |
| $`\Lambda^{*}_{\mathrm{fl}}, s_{\Lambda,\mathrm{fl}}`$ | `flush_light_threshold`, `flush_light_sharpness` | flush light gate (> 0: more light permits) |
| $`\Lambda^{*}_{\mathrm{sh}}, s_{\Lambda,\mathrm{sh}}`$ | `shed_light_threshold`, `shed_light_sharpness` | shed light trigger (< 0: short days; > 0: bright light) |
| | `light_window` | shortwave running mean [day] (optional) |
| $`M^{*}_{\mathrm{wet}}, s_{\mathrm{wet}}`$ | `flush_water_sum`, `flush_water_sharpness` | wet sum that permits flushing [MPa day] (optional) |
| $`M^{*}_{\mathrm{dry}}, s_{\mathrm{dry}}`$ | `shed_water_sum`, `shed_water_sharpness` | dry sum that triggers senescence [MPa day] (optional) |
| $`\psi_{\mathrm{tlp}}`$ | *(derived)* | turgor-loss point from `leaf_pi0` and `leaf_elastic_mod` |
| $`k_{\mathrm{turn}}`$ | `leaf_lifespan_toc` ([pft]) | background turnover $`1/\ell`$ [yr⁻¹] |

The `[phenology]` section is gated by its presence: absent, every PFT keeps the no-cue defaults;
present, every key above is required except the optional four water keys and `light_window`. A sharpness
of 0 is rejected (the switch would be a constant 0.5). The keys of the previous scheme (`k_flush_max`,
`tau_flush`, `cue_sharpness`, `phen_a/b/c`, the cold-drop and chilling keys, `pft.evergreen`, …) are
refused by name with a pointer to their replacement.

### Drivers

| cue | driver | source |
|---|---|---|
| TEMP | daily-mean air temperature | fast-loop `pheno_tair` reduction |
| LIGHT, day length | day length | `meds_time::daylength(lat, doy)` |
| LIGHT, shortwave | daily-mean incident shortwave at the canopy top, area-weighted over patches | fast-loop `pheno_rad` reduction |
| WATER | the cohort's predawn (daily-maximum) leaf water potential | `dmax_psi_leaf`, the same value that drives the stomatal stress limb |

No cue reads soil data. All the cue sums and the shortwave mean are per-cohort state and are written to
restart files.

## References

- **Botta et al. (2000)**, *Glob. Change Biol.* 6:709 — growing-degree-day budburst.
- **Delpierre et al. (2009)**, *Agric. For. Meteorol.* 149:938 — autumn senescence from a cold-degree-day
  sum under shortening days.
- **White et al. (1997)**, *Glob. Biogeochem. Cycles* 11:217 — day-length and temperature leaf offset.
- **Xu et al. (2016)**, *New Phytol.* 212:80 — drought deciduousness keyed on the turgor-loss point.
- **Kim et al. (2012)**, *Glob. Change Biol.* 18:1322 — light-driven leaf phenology in the tropics.
- ED2 `ED/src/dynamics/phenology_driv.f90`, `phenology_aux.f90`; the design record,
  `docs/dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md`.

## Code map

| Concept | Routine |
|---|---|
| cues → tendencies → potential rates | `meds_phenology`: `phenology_kernel` (+ private `accumulate`, the water switches) |
| tendencies → leaf loss + flush cap | `meds_phenology`: `leaf_turnover_step` |
| types (env / params / state / out) + cue bits | `meds_phenology_types` |
| helpers | `meds_numerics`: `logistic`, `clamp01`; `meds_time`: `daylength`, `doy_effective` |
| per-cohort advance (slow loop) | `meds_vegetation_dynamics`: `advance_leaf_phenology`, `flatten_pheno_params` |
| turnover-first carbon demand | `meds_vegetation_dynamics`: `cohort_carbon_demand` (called from `compute_carbon_allocation`) |
| C-API | `meds_c_api_phenology`: `meds_phenology_step`, `meds_leaf_turnover_step`, `meds_daylength` |
| Python front end + example | `meds.plant.pheno`; `examples/example02_canopy_phenology/` |
