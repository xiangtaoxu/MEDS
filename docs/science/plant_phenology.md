# Leaf phenology

MEDS models leaf phenology in two layers. The **cue kernel** (`meds_phenology`) turns daily
environmental cues and per-PFT traits into two smoothed **tendencies** in $[0,1]$: a leaf **flush**
tendency $f$ and a leaf **shed** (senescence) tendency $s$. The **carbon layer** turns the tendencies
into leaf growth and leaf loss (`leaf_turnover_step`, §5). The kernel's state is the two tendencies
plus a small **cue memory** (temperature sums, water-potential sums, a running mean of the hours of light),
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
requirement: the light gate (below; at a low `par_min`, the photoperiod) does the work chilling was
meant to do, holding the flush until days are long enough, and adding chilling did not improve the
fits (§7).

**Light (`CUE_LIGHT`).** The hours of light $H$ [h day⁻¹]: the time each day the PAR reaching the
cohort's top exceeds the PFT's `par_min` [µmol m⁻² s⁻¹], as a running mean $`\bar H`$ over
`light_window` days:

```math
\bar H \mathrel{+}= w\,(H - \bar H), \qquad w = \min\!\big(1,\ dt/\max(\texttt{light\_window},\,dt)\big) \qquad(4)
```

In the coupled model the fast loop counts $H$: each fast step adds its length when the visible light,
beam plus downward diffuse, reaching the top of the cohort's layer in the two-stream canopy radiation
exceeds `par_min` in photons. A cohort under others counts fewer hours than one in the sun. A cohort
with no light memory yet (a cold start, a recruit) starts its mean from its first day's hours, so its
first days do not read as darkness.

One cue serves every PFT, and `par_min` sets what it measures. With a low `par_min` (a few
µmol m⁻² s⁻¹), $H$ is the photoperiod, give or take twilight and the heaviest overcast: a calendar,
the extratropical cue of photoperiodic budburst and growth cessation. With a high `par_min` (hundreds
of µmol m⁻² s⁻¹) it counts the bright hours, with their clouds, their year-to-year swings and the
cohort's place in the canopy: the cue of tropical leaf exchange in the bright dry season (Kim et al.
2012). At three sites one such variable, with a fitted threshold, did as well as the better of day
length and running-mean PAR (§7).

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
\cdot \underbrace{\sigma\big(s_{H,\mathrm{fl}} (\bar H - H^{*}_{\mathrm{fl}})\big)}_{\text{LIGHT}}
\cdot \underbrace{\sigma\big(s_{\mathrm{wet}} (M_{\mathrm{wet}} - M^{*}_{\mathrm{wet}})\big)}_{\text{WATER}}
\qquad(6)
```

```math
S_{\mathrm{sh}} = \max\Big(
\underbrace{\sigma\big(s_C (C - C^{*})\big)}_{\text{TEMP}}
\cdot \underbrace{\sigma\big(s_{H,\mathrm{sh}} (\bar H - H^{*}_{\mathrm{sh}})\big)}_{\text{LIGHT}},\
\underbrace{\sigma\big(s_{\mathrm{dry}} (M_{\mathrm{dry}} - M^{*}_{\mathrm{dry}})\big)}_{\text{WATER}}
\Big) \qquad(7)
```

Only the cues in each side's mask enter; an absent factor is 1 in a product and absent from the max.
An empty flush mask always flushes ($`S_{\mathrm{fl}}=1`$); an empty shed mask never senesces
($`S_{\mathrm{sh}}=0`$).

- **Flushing needs every flush cue** (a product): warm enough *and* enough hours of light *and* wet
  enough.
  The light gate applies all year, so at a low `par_min` it both holds back budburst in a warm spell
  in late winter and turns the flush off in late summer, when the warmth sum is still high. For the
  second it has to be sharp, a photoperiod threshold (§7).
- **Senescence has two triggers** (a max): the **seasonal** one, cold *and* the light condition
  (Delpierre et al. 2009), and the **water** one, sustained drought on its own. A drought need not
  wait for autumn.

On the shed side a negative light sharpness means few hours of light trigger senescence (short days
at a low `par_min`, dim days at a high one); a positive one means many bright hours do, as in
tropical leaf exchange.

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
  lies below that cover. A sharp flush light gate (§7) turns the flush off in autumn and the
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
| temperate deciduous (Harvard) | TEMP + LIGHT | TEMP + LIGHT | 0 | low `par_min`: flush on warmth while days are long; senesce to bare on cold once days shorten |
| evergreen conifer (Hyytiälä) | TEMP + LIGHT | TEMP + LIGHT | the share kept through winter (0.9) | high `par_min`: flush on warmth and many bright hours; autumn needle fall stops at the floor |
| drought deciduous (Palo Verde) | WATER + LIGHT | WATER + LIGHT | the share kept through the dry season (0.27) | low `par_min`: senesces as the days shorten and after a long drought; flushes once the days lengthen and the leaves are wet |
| light leaf exchange (BCI) | WATER | LIGHT + WATER | 0.95 | high `par_min`: a water threshold it never reaches; senescence on many bright hours while the canopy refills, so it stays full |

Evergreen and deciduous differ in one number. The emergent leaf lifespan follows from the rates and
the floor; it is not a parameter of the kernel.

## 7. Evidence and status

`examples/example02_canopy_phenology` fits the compiled kernel at four sites, every one with its
light from ERA5-Land through MEDS's forcing reader:

- **Harvard Forest** (temperate deciduous; `par_min` 2, the photoperiod): MODIS LAI, the HF003
  leaf-fall observations and the HF069 litter baskets.
- **Hyytiälä** (Scots pine; `par_min` 99): the ICOS monthly needle litter.
- **Barro Colorado Island** (a light leaf exchanger; `par_min` 1081): GLiMP litter traps. The water
  driver is a surrogate, the tower's soil water content through a retention curve fitted to the
  paired soil samples of Kupers et al. (2019).
- **Palo Verde** (a drought-deciduous dry forest; `par_min` 5, the photoperiod): the monthly leaf litter
  of Xu et al. (2016) and MODIS LAI. The water driver is hypothetical, the canopy predawn leaf ψ of a
  MEDS run at the site.

Five results shaped the design (`docs/dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md` §8–§12):

- **No chilling.** Removing the chilling requirement changed neither temperate fit.
- **One light cue, the hours of light above a fitted `par_min`.** Day length and running-mean PAR
  each fit some sites and fail others: Harvard is best on day length (loss 0.0445 against 0.052 on
  PAR), Hyytiälä on PAR (0.70 against 1.03) and the BCI leaf exchanger on PAR (0.065 against 0.101).
  The hours of light with a fitted `par_min` came within 4 % of the better of the two at every site
  (Harvard 0.046 at 33 µmol m⁻² s⁻¹, Hyytiälä 0.59 at 280, BCI 0.067 at 700; tower PAR), so one
  variable serves all three. Refitted on ERA5-Land, the deciduous canopies count the photoperiod
  (at Palo Verde with `par_min` free) and the evergreen ones the bright hours (within bounds that
  hold them to bright light).
- **Light makes the drought-deciduous canopy.** At Palo Verde water alone fits poorly (loss 0.064):
  the leaves fall in December–January, before the water potential drops. The photoperiod alone
  reaches 0.025 and both cues 0.023, as
  photoperiodic leaf fall and spring flushing are known in tropical dry forests (Borchert & Rivera
  2001; Rivera et al. 2002).
- **The flush gate must be sharp where it ends the flush.** A gradual day-length gate (1 per hour) is
  still partly open in October, so the canopy refills while it senesces and a deciduous stand drops
  about 1.9 canopies of leaves a year, which no timing observation shows. Held to one canopy a year,
  the Harvard fit makes the gate a step (8 per hour). The same holds for the water switches: a
  sharpness of 5 per MPa day leaves a 0.7 % senescence floor that, refilled by the flush, sheds a
  quarter canopy a year.
- **Other sharpnesses are not identifiable** from these data (different values fit equally well), so
  they keep fixed defaults (0.04 per K day for warmth, 0.1 per K day for cold) and the centres are
  fitted.

What is **not** validated: no coupled MEDS run's leaf-area cycle has yet been scored against
observations; the water cue has been exercised only on a soil surrogate (BCI) and on a MEDS run's
leaf potential without the shed-to-water feedback (Palo Verde), not on a run's own predawn leaf
potential with it. The fits capture the mean seasonal cycle; four to five scored years per site do
not constrain its year-to-year variation.

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
| | `par_min` | the PAR at the cohort's top that counts as an hour of light [µmol m⁻² s⁻¹] |
| $`H^{*}_{\mathrm{fl}}, s_{H,\mathrm{fl}}`$ | `flush_light_hours`, `flush_light_sharpness` | flush light gate [h day⁻¹], [h⁻¹] (> 0: long or bright days permit) |
| $`H^{*}_{\mathrm{sh}}, s_{H,\mathrm{sh}}`$ | `shed_light_hours`, `shed_light_sharpness` | shed light trigger (< 0: short or dim days; > 0: bright days) |
| | `light_window` | running mean of the hours of light [day] |
| $`M^{*}_{\mathrm{wet}}, s_{\mathrm{wet}}`$ | `flush_water_sum`, `flush_water_sharpness` | wet sum that permits flushing [MPa day] (optional) |
| $`M^{*}_{\mathrm{dry}}, s_{\mathrm{dry}}`$ | `shed_water_sum`, `shed_water_sharpness` | dry sum that triggers senescence [MPa day] (optional) |
| $`\psi_{\mathrm{tlp}}`$ | *(derived)* | turgor-loss point from `leaf_pi0` and `leaf_elastic_mod` |
| $`k_{\mathrm{turn}}`$ | `leaf_lifespan_toc` ([pft]) | background turnover $`1/\ell`$ [yr⁻¹] |

The `[phenology]` section is gated by its presence: absent, every PFT keeps the no-cue defaults;
present, every key above is required except the optional water keys. A sharpness
of 0 is rejected (the switch would be a constant 0.5). The keys of the previous scheme (`k_flush_max`,
`tau_flush`, `cue_sharpness`, `phen_a/b/c`, the cold-drop and chilling keys, `pft.evergreen`, …) are
refused by name with a pointer to their replacement.

### Drivers

| cue | driver | source |
|---|---|---|
| TEMP | daily-mean air temperature | fast-loop `pheno_tair` reduction |
| LIGHT | the hours a day the visible light (beam + downward diffuse) reaching the top of the cohort's layer exceeds `par_min` | two-stream `incid_top`, counted per fast step in `light_hours_accum` |
| WATER | the cohort's predawn (daily-maximum) leaf water potential | `dmax_psi_leaf`, the same value that drives the stomatal stress limb |

No cue reads soil data. All the cue sums and the light mean are per-cohort state and are written to
restart files.

## References

- **Botta et al. (2000)**, *Glob. Change Biol.* 6:709 — growing-degree-day budburst.
- **Delpierre et al. (2009)**, *Agric. For. Meteorol.* 149:938 — autumn senescence from a cold-degree-day
  sum under shortening days.
- **White et al. (1997)**, *Glob. Biogeochem. Cycles* 11:217 — day-length and temperature leaf offset.
- **Xu et al. (2016)**, *New Phytol.* 212:80 — drought deciduousness keyed on the turgor-loss point.
- **Borchert & Rivera (2001)**, *Tree Physiol.* 21:213 — photoperiodic control of leaf shedding and
  dormancy in tropical trees.
- **Rivera et al. (2002)**, *Trees* 16:445 — increasing day length induces flushing of tropical dry
  forest trees before the rains.
- **Kim et al. (2012)**, *Glob. Change Biol.* 18:1322 — light-driven leaf phenology in the tropics.
- ED2 `ED/src/dynamics/phenology_driv.f90`, `phenology_aux.f90`; the design record,
  `docs/dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md`.

## Code map

| Concept | Routine |
|---|---|
| cues → tendencies → potential rates | `meds_phenology`: `phenology_kernel` (+ private `accumulate`, the water switches) |
| tendencies → leaf loss + flush cap | `meds_phenology`: `leaf_turnover_step` |
| types (env / params / state / out) + cue bits | `meds_phenology_types` |
| helpers | `meds_numerics`: `logistic`, `clamp01`; `meds_time`: `doy_effective` |
| per-cohort advance (slow loop) | `meds_vegetation_dynamics`: `advance_leaf_phenology`, `flatten_pheno_params` |
| hours of light at each cohort's top | `meds_canopy_radiation`: `solve_band` (`incident`), `rad_flux_t%incid_top`; `meds_fast_dynamics`: `light_hours_accum` |
| turnover-first carbon demand | `meds_vegetation_dynamics`: `cohort_carbon_demand` (called from `compute_carbon_allocation`) |
| C-API | `meds_c_api_phenology`: `meds_phenology_step`, `meds_leaf_turnover_step`, `meds_daylength` |
| Python front end + example | `meds.plant.pheno`; `examples/example02_canopy_phenology/` |
