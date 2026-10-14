# MEDS fast calibration, revision 2: user-defined filters and keys, and a staged fit

**Status:** approved by the owner on 2026-10-02 (rev 3, after two rounds of decisions, §0).
R1–R5 implemented on 2026-10-02 (#345); R6 has run once at BCI, and its findings await the owner
(§13). It builds on `MEDS_FAST_CALIBRATION_PLAN.md` (§1–14, the method and the fits up to the
shipped 11-key set B) and on the BCI GPP diagnosis of 2026-10-01/02 (§1). Section numbers below are
this document's own.

**Scope:** a fast, rough, strictly site-based calibration against one flux tower:
- **No remote sensing.** A regional calibration with remote-sensing data is a separate, later effort.
- **Rough uncertainty only.** A robust analysis is left to dedicated tools (PEcAn) later.

**The changes:**
1. **Every data filter is a user setting** with a documented, recommended default.
2. **Which keys are tuned, and how strongly the data may move them, is a user setting** too: the
   key list and each key's prior.
3. **Default priors come from PFT-level meta-analyses,** not from the code's defaults.
4. **The fit runs in stages along the causal chain,** each stage with the cheapest model that can
   see its keys.
5. **Uncertainty is a rough estimate** from the fit's own Jacobian, at almost no extra cost.

## 0. Owner decisions (2026-10-02)

| # | decision |
|---|---|
| 1 | No observation model for the tower's respiration bias (soil-carbon dynamics are not developed enough). Filter low-quality data and give GPP a larger σ (§6). |
| 2 | Every filter is user-definable, with a recommended default and its explanation (§3). |
| 3 | Night NEE leaves the fast fit; `rd_vcmax_ratio` is fixed. |
| 4 | Cohort-level sunlit/shaded leaves: filed as #343 to revisit; no implementation now (§5). |
| 5 | `theta_j` becomes an optional key, range 0.7–0.9. The user defines or controls which keys are tuned, including through their priors (§4). |
| 6 | "Screening" explained (§4.3); it reports and no longer drops keys by default. |
| 7 | Priors: a general default from PFT-level meta-analyses (§4.4). The `vcmax25` prior comes from Kattge et al. 2009, the easier source. Emergent-trait priors are left to PEcAn. |
| 8 | No GOSIF or other remote-sensing target. |
| 9 | Rough uncertainty at low cost, not a robust analysis (§8). |
| 10 | The plan is its own file. |
| 11 | The water-stress stage is in the first implementation (§7.4). |
| 12 | Dim-light GPP: a larger σ rather than dropping hours, in the simple σ_abs + σ_rel form: 2.5 + 0.15·GPP (§6). |
| 13 | The GPP u* default is 0.4 (§3.3). |
| 14 | The plan is approved and goes in as its own file. |

## 1. What the BCI diagnosis showed

Ten-day wet-season window (2015-09-11), onset build (#341), default parameters unless stated.
"Ratio" is model/tower daytime GPP.

| finding | number | what it means for the calibration |
|---|---|---|
| default GPP is high | ratio 1.51; morning (6–10 h) 1.88, midday 1.36 | the fit lowers `vcmax25` to its floor (25) to cut GPP |
| accounting, PAR timing and amount, light budget | all exact or within 2 % | not the cause |
| leaf quantum yield | the formula is CLM5's; ~15 % above leaf measurements at matched CO2 and temperature (Ehleringer & Björkman 1977) | not the main cause |
| Γ* | Bernacchi 2001; would need ×1.5–2.3 to explain the gap | not the cause (pressure scaling is #342, negligible at BCI) |
| each cohort's leaves get the cohort's mean light | sunlit/shaded split: ratio 1.51 → 1.34 | the largest single structural effect; deferred to #343 (§5) |
| leaf clumping 0.8 | floor PAR 2.4 % under closed canopy against 1–2 % measured in tropical understories (Chazdon & Fetcher 1984); 1.0 gives 1.1 %, ratio 1.40 | clumping is informed by the light profile |
| θ_J 0.9 (CLM5: 0.7) | ratio 1.42 at 0.7 | trades with `vcmax25`; optional key (§4) |
| sunlit/shaded + clumping 1.0 + θ_J 0.7 | ratio 1.16; midday 1.01; dawn and dusk 1.3–1.4 | the left-over excess sits where tower GPP is smallest |
| Vcmax plasticity | 45·e^(−0.118 L), matching the Panama gradient (top 40–46, understory ~27, as quoted in the diagnosis) | not a fit target |
| tower GPP = daytime NEE + one respiration value per day | that respiration averages 4.0 µmol m⁻² s⁻¹, about what the soil chambers measure for the soil alone; the provider notes it "appeared underestimated", since above-ground respiration "can contribute up to 40-50%" of the total (corrected 2026-10-03: this read "40–50 % below soil chambers") | tower GPP is low by a roughly constant amount, which weighs most in dim light: a larger σ (§6) |
| low turbulence in the morning | 7 h GPP/PAR +15 % (u* ≥ 0.3), +26 % (u* ≥ 0.5) over FLAG alone; no change from 9 h | a u* filter on GPP (§3) |
| FLAG | removes every rain half-hour; gap-filled hours share the measured light response | keep FLAG = 1 |
| energy-balance closure | daily median 0.72; 11 % of days within 0.8–1.2 | no closure-based day mask |
| fit B | `vcmax25` and `rd_vcmax_ratio` at their floors | both pulled by the tower's low respiration |

## 2. Principle: the user decides, the tool recommends and reports

**Where the settings live:**
- **The site declaration** (`calibration.toml` in the example) holds every choice that changes what
  is fitted: the filters (§3), the targets and their σ (§6), the keys and their priors (§4), and the
  stages (§7).
- **A commented reference file,** `scripts/calibrate_fast/site_reference.toml`, documents every
  setting: its default, the recommended range, and why. It follows `meds_config_main.toml`'s style.
  - A key left out of the site file takes the reference default.
  - An unknown key stops the run, as N-10 does for the model config.

**What the tool reports at the start of every fit:**
- per target, the hours kept by each filter in turn;
- the target's diurnal mean before and after filtering;
- the keys being fitted, with their ranges and priors.

The same is written into the fit report.

## 3. Filtering the tower data

### 3.1 Settings and defaults

Filters are set per target. The tool refuses a filter that would bias its own target: a u* filter on
the u* target.

| setting | default | recommended | why | BCI cost |
|---|---|---|---|---|
| `[tower].flag_column`, `flag_good` | `FLAG`, 1 | measured only | gap-filled hours are a neural network's prediction | ~40 % of daytime hours |
| `[tower].forcing_qc` | `LWdown_qc`, `Wind_qc` | observed forcing only | a trial driven by gap-filled forcing scores the gap-filling | site-dependent |
| `[tower].daytime_sw` | 10 W m⁻² | 10 | defines "daytime" for every daytime target | — |
| `[tower].closure` | `"bowen"` | `"bowen"` or `"none"` | FLUXNET2015's correction keeps the Bowen ratio | — |
| `[targets.gpp].ustar_min` | **0.4 m s⁻¹** | 0.3–0.5; each fit reports the site's plateau value (§3.3) | low turbulence in the morning under-reads uptake; the effect stops by 9 h | −17 % of daytime hours at 0.3, −29 % at 0.4, −41 % at 0.5 |
| `[targets.gpp].par_min` | 0 | 0: the σ of §6 down-weights dim hours instead | dropping hours is the blunter alternative | −6.7 % at 200 |
| `[targets.gpp].hours` | all daylight | all daylight | a local-hour window (e.g. 9–17) as a blunter alternative to u* | |
| `[targets.le/h].ustar_min` | 0 | 0 | available; not needed at BCI | |
| `[targets.<t>].closure_range` | off | off | a day mask on daily (H+LE)/Rnet; kept for sites with good closure | at BCI it would keep 11 % of days |
| `[targets.nee_night].on` | **false** | false (decision 3) | night fluxes under-read in low turbulence; frozen soil carbon | — |
| `[targets.albedo].min_sw`, `min_solar_elevation` | 200 W m⁻², 0° | 200, 0–30° | low-sun albedo is the least reliable | small |
| rain | (FLAG) | nothing more | at BCI, FLAG already removes every rain half-hour | 0 |

### 3.2 Not adopted, and why

- **A global u\* mask:** it would delete the low-u\* hours from the u\* target.
- **Closure day masks at BCI:** 11 % of days pass.
- **Gap-filled hours as targets:** they are model output, even though they show the same light
  response.

### 3.3 The GPP u* threshold

**The data.** BCI, whole record, FLAG = 1. The tower's GPP/PAR (ratio of sums) by u\* class:

| u\* class [m s⁻¹] | 7 h | 8 h | 11–13 h (control) |
|---|---|---|---|
| 0.15–0.20 | 0.0212 | 0.0192 | — |
| 0.20–0.30 | 0.0197–0.0205 | 0.0173–0.0183 | 0.0195–0.0202 |
| 0.30–0.40 | 0.0234–0.0253 | 0.0179–0.0188 | 0.0181 |
| 0.40–0.50 | 0.0253 | 0.0193 | 0.0163 |
| 0.50–0.60 | 0.0318 | 0.0223 | 0.0160 |
| 0.60–0.80 | 0.0295 | 0.0215 | 0.0150 |
| > 0.80 | 0.0300 | 0.0214 | 0.0129 |

- **Mornings:** at 7 h the efficiency is flat below u\* 0.3, rises through 0.3–0.5, and levels off
  near 0.030 from 0.5.
- **Plateau test:** the test used for night-time u\* thresholds (Reichstein et al. 2005; Papale et
  al. 2006: the lowest class within 95 % or 99 % of the mean of the classes above) gives 0.5 for 7 h.
- **Midday control:** efficiency falls with u\*, because windy hours are sunnier and nearer light
  saturation. The morning rise is therefore not a confounded light effect; if anything it is
  understated.

**The literature:**
- u\* filtering is standard for night fluxes (Reichstein 2005; Papale 2006), not for daytime ones.
- The morning deficit is the low-turbulence transition after sunrise.
- The provider applied u\* > 0.4 to its own respiration estimate for this tower.

**Recommendation: 0.4.**
- It is the provider's own threshold for this tower.
- It removes about 60 % of the 7 h deficit: 0.0225 with FLAG alone, 0.0269 at 0.4, against a
  plateau near 0.030.
- It costs 29 % of daytime hours, against 41 % at 0.5.
- What remains at 7 h (about 1 µmol m⁻² s⁻¹ at GPP 8) is a quarter of that hour's σ under §6.

0.5 is the data's own threshold and the choice for a stricter fit. Every fit prints this table and
the plateau-test value for its tower, so the user can see what that site suggests.

## 4. Keys and priors

### 4.1 How the user controls what is tuned

There are two controls.

**1. The key list.** The registry is the menu. Every key there has a range, a transform, a default
prior and a recommended state:
- **fit:** in the default set;
- **optional:** available, off unless the user turns it on;
- **fixed:** not tunable in the fast fit, with the reason.

The site file can set the list explicitly (`[fit].keys = [...]`) or edit the default
(`[fit].add`, `[fit].remove`).

**2. The prior of each key** (`[priors.<key>]`: centre, sd, range, source) sets how far the data may
move it. A narrow prior holds a key near its centre but still lets strong evidence nudge it; a wide
prior lets the data decide. A key fitted with a very narrow prior is close to fixed, and its
posterior-to-prior σ ratio near 1 shows that the tower taught nothing.

**Precedence:** site file > registry default. The fit report lists every key with the source of its
range and prior.

### 4.2 The menu (defaults for a tropical evergreen PFT)

| stage | key | state | range | default prior: centre (± sd = range/4 unless given) | source |
|---|---|---|---|---|---|
| 1 optics | `leaf_reflect_nir` | fit | 0.30–0.55 | 0.45 | CLM5 broadleaf evergreen tropical |
| | `leaf_transmit_nir` | **fixed** | 0.15–0.40 | 0.25 | collinear with the reflectance (−0.996): the albedo sees their sum |
| | `leaf_clumping` | fit | 0.50–1.00 | 0.80 | current default; BCI's floor light points to 0.9–1.0 (§1) |
| | `leaf_angle_mean` | fit | 30–65° | 45° | current default (spherical is 57.3°) |
| 2 photosynthesis | `vcmax25` (top of canopy) | fit | 25–80 (log) | the PFT mean, sd = its across-species spread | Kattge et al. 2009, read from the paper's PFT table in R2; PEcAn's trait meta-analysis (LeBauer et al. 2013) can replace it later |
| | `theta_j` | **optional** | 0.70–0.90 | 0.80 ± 0.05 | decision 5; CLM5 uses 0.7, MEDS 0.9 |
| | `jmax_vcmax_ratio` | optional | 1.4–2.2 | 1.67 ± 0.15 | Medlyn et al. 2002; Panama leaves 1.57–1.78 |
| | `phi_psii` | fixed | — | 0.85 | Björkman & Demmig 1987; CLM5. Without a respiration observation model (decision 1), fitting it would absorb the tower's bias |
| 3 energy and aerodynamics | `stomatal_g1` | fit | **1.5–8.0** (log) | 3.77, log-sd 0.35 | Lin et al. 2015, tropical rainforest trees (167 species); the species spread reaches ~8, and set B sat at the old 6.0 ceiling |
| | `stomatal_g0` | fit | 0.001–0.05 (log) | 0.01 | current |
| | `z0m_ratio` | fit | 0.05–0.20 | 0.13 | current |
| | `leaf_width` | fit | 0.02–0.15 m (log) | 0.04 | current |
| | `dsl_dmax` | fit | 0.005–0.05 m (log) | 0.015 | current |
| | `dewmx` (interception on) | fit | 0.05–0.3 (log) | 0.1 | current |
| | `intercept_k` | fixed | 0.3–0.8 | 0.5 | rough response, σ ratio 0.89 |
| 4 water stress | `wstress_sref_stomata` | fit | 0.5–5.0 | 2.0 | Sabot et al. 2022 |
| | `stomata_psi_onset` (#341) | fit | 2ψ_tlp–0 | ψ_tlp/2 | #341 |
| | `root_beta` | optional | e-folding 0.25–2 m | current | plan §4.4 |
| fixed | `rd_vcmax_ratio` | fixed | — | 0.015 | decision 3 |
| fixed | Γ*, Kc, Ko and their Ea; `leaf_absorptance` (cancels); the plasticity slope; `ds_vcmax`; soil and structure | fixed | | | plan §4.6 |

The default fitted set is 3 optics + 1 photosynthesis + 5 energy + 2 water-stress keys: 11 keys, 12
with interception. The optional keys add 2 to stage 2 (`theta_j`, `jmax_vcmax_ratio`) and 1 to
stage 4 (`root_beta`).

### 4.3 What "screening" means

**What it is:** before a fit, the tool runs the model once with each key nudged up and once nudged
down, at its prior centre. That is the first Jacobian: 2 runs per key per window, run in parallel.
From it, for every key, it reports four things:
- **sensitivity:** how much each target moves per unit of the key's prior width;
- **information:** how much the tower could shrink the key's prior. A posterior-to-prior σ ratio of
  0.9 or more means the tower can't tell;
- **collinearity:** keys the targets cannot tell apart, like the two NIR optics;
- **smoothness:** whether the up and down responses agree, or the key acts like a switch.

**What the user controls:** today the screening also drops keys automatically. In this revision it
**only reports, by default** (`[fit].screening = "report"`). A key the user asked for is fitted
unless the user sets `screening = "drop"`. To "screen" a key means to put it through this test
first and decide from the report.

### 4.4 Default priors from meta-analyses (decision 7)

Each key's default prior is a distribution from a PFT-level synthesis (§4.2):
- normal or log-normal in the transformed space;
- centred on the PFT mean, with the across-species spread as the width (range/4 where no synthesis
  exists);
- recorded with its `source`, so a user can see and replace it.

The format is one PEcAn's trait meta-analysis can supply later (LeBauer et al. 2013). Priors on
emergent leaf traits (Asat, ci/ca), which would constrain combinations of keys, are left to that
PEcAn work (owner, round 2).

## 5. Sunlit and shaded leaves: deferred to #343

**Owner, round 2:** filed as #343 to revisit; no implementation now. The issue holds:
- the comparison with ED2 (cohort-mean light, as MEDS), FATES (sunlit and shaded leaves per leaf
  layer) and CLM5 (two big leaves);
- the prototype's BCI numbers: GPP −11 %, run time +3–8 %;
- the pros and cons.

**For this plan:** the calibration runs on the mean-light canopy. Part of the mean-light bias is
therefore absorbed by `vcmax25`, `leaf_clumping` and, if turned on, `theta_j`, and their calibrated
values belong to that setting. The fit report says so.

## 6. The cost function

Weighted least squares, as today, with three changes:

```
Φ(u) = Σ_t w_t Σ_{i∈t} ((y_model,i(u) − y_obs,i) / σ_t,i)²  +  Σ_k ((u_k − μ_k) / σ_k)²
```

1. **GPP σ = 2.5 + 0.15·GPP µmol m⁻² s⁻¹** (decisions 1 and 12; was 1.5 + 0.15·GPP).
   - **Why:** the tower's single daily respiration value is likely too low: the provider notes it "appeared underestimated" against soil chambers, since above-ground respiration "can contribute up to 40-50%" of the total (corrected 2026-10-03: this read "is reported 40–50 % low"). That is about
     2 µmol m⁻² s⁻¹ at **every** daytime hour, not only in dim light, on top of the random error.
   - **Why this form:** independent errors add in quadrature, √((1.5 + 0.15·GPP)² + 2²). The linear
     form that stands in for it keeps σ_rel and raises σ_abs to 2.5. It equals the quadrature sum at
     GPP 0 and is 12 % wider at GPP 20 (5.5 against 4.9).
     - Adding the terms linearly (3.5 + 0.15·GPP) would count the systematic error twice.
     - A step in σ below a PAR threshold leaves the 200–500 class, the 7–9 h mornings, as exposed as
       before.
   - **Where a constant 2 µmol bias weighs in the fit** (χ² per hour it would add; BCI, FLAG = 1,
     u\* ≥ 0.4; share of rows in brackets):

     | σ form | PAR 0–200 (12 %) | 200–500 (15 %) | 500–1000 (22 %) | > 1000 (51 %) |
     |---|---|---|---|---|
     | 1.5 + 0.15·GPP (before) | 1.04 | 0.50 | 0.27 | 0.19 |
     | 3.0 below PAR 200, else 1.5 (+0.15·GPP) | 0.33 | 0.50 | 0.27 | 0.19 |
     | 3.5 + 0.15·GPP | 0.25 | 0.17 | 0.11 | 0.09 |
     | √((1.5 + 0.15·GPP)² + 2²) | 0.49 | 0.32 | 0.21 | 0.16 |
     | **2.5 + 0.15·GPP (adopted)** | **0.45** | **0.27** | **0.17** | **0.13** |

   - **Settings:** the existing `sigma_abs` and `sigma_rel`; no new σ setting. The reference file
     records 2.5 as the recommendation for a tower whose GPP is built from a daily respiration, with
     this reason.
2. **Effective-sample weights in the fit,** not only in the covariance.
   - w_t = n_eff,t / n_t, from the lag-1 autocorrelation of each target's residuals at the default,
     refreshed once at the stage-5 start.
   - Hourly targets would otherwise outweigh daily ones by row count alone. The fit report shows
     each target's share of Φ.
3. **σ scale report** (no extra runs): each target's χ²/n_eff at the MAP. A value of 2 or more marks
   structural error in that target. It feeds the uncertainty (§8).

Targets stay as in plan §5.2, without night NEE. A robust (Huber) loss is available as an option,
default off.

## 7. A staged, process-guided fit

Within the fast loop, the coupling is mostly one-directional:

```
optics → light per leaf → photosynthesis (A, ci) → stomata → energy balance (LE, H, leaf T)
                               ↑                                                 |
                               └─────────── leaf temperature, VPD (weak) ───────┘
```

Each stage fits its keys with the upstream results fixed and the downstream drivers taken from the
latest full run. A short joint polish then settles the weak feedback and gives the covariance. A
user can run any stage alone (`--stage`).

### 7.1 Stage 1, optics: the radiation solver alone

- **What runs:** the two-stream only, over the frozen stand (cohort LAI, WAI and heights per patch),
  with the forcing's sun angle and beam/diffuse split. No physiology, no energy balance.
- **Needs:** a Python binding, `meds.canopy.radiation`, like `meds.plant.leaf`.
- **Data:** albedo over the whole record (~44,000 hours, not 8 windows); understory PAR where a site
  measures it.
- **Cost:** seconds per evaluation.

### 7.2 Stage 2, photosynthesis: a canopy of leaf kernels

With the leaf drivers fixed (absorbed PAR per leaf, leaf temperature, leaf VPD, canopy-air CO2,
predawn ψ, boundary-layer conductance), canopy GPP is a sum of leaf-kernel calls:

```
GPP_site,h = Σ_patches area_p Σ_cohorts LAI_c · A_gross(PAR_c,h, T_c,h, VPD_c,h, CO2_h, ψ_c; θ)
```

- **Vcmax per cohort:** the top-of-canopy value × e^(−k·LAI above), recomputed inside the sum.
- **Fitted:** `vcmax25`, plus `theta_j` and `jmax_vcmax_ratio` if the user turns them on.
- **Data:** all 13,940 measured daytime hours (fewer after the u\* filter), not 8 windows.
- **This is the answer to "Vcmax from a day or two":** yes, without running the coupled model. With
  a batch kernel and leaf-area-weighted cohort sampling (as the BCI diagnosis did), one evaluation
  of the whole record takes about a second, so there is no reason to limit it to a day or two.
- **Needs:**
  - the per-cohort leaf drivers as a supported hourly output (today a scratch hijack only; the
    `cdiag` fields are already computed);
  - a batch leaf call in Python (`leaf_gas_exchange_batch` exists in Fortran; one Python call per
    leaf costs ~190 µs).
- **Outer passes:** the feedback gs → LE → leaf temperature → A is second order. Each pass reruns
  the frozen windows once with the new leaf set, refreshes the drivers and refits; stop when
  `vcmax25` moves less than 1 %. Expect 2–3 passes.
- **Gate G8:** with the full run's own drivers and parameters, the canopy of kernels reproduces that
  run's hourly GPP within 1 % (median). The kernel already matches the model's leaves within 4.3 %
  with CO2 fixed at 400 ppm and no water stress; the real CO2 and stress should close most of that.

### 7.3 Stage 3, energy partition and aerodynamics: the coupled fast loop

- **What runs:** the existing method (Levenberg–Marquardt, central-difference Jacobian, worker pool,
  10-day frozen windows), on the stage-3 keys with stages 1–2 fixed.
- **Cost:** 5–6 keys give (2k+1)·W ≈ 88–104 trials per iteration, against 184 for set B, from a
  better start.
- **To test:** whether 5-day windows hold the diurnal energy partition as well. That would halve
  the trial time.

### 7.4 Stage 4, water stress (in the first implementation)

- **What runs:** frozen dry-season runs, January to May (~50 s each), chained from the wet-season
  state.
- **Targets:** dry-season LE, GPP and evaporative fraction.
- **Keys:** `wstress_sref_stomata` and `stomata_psi_onset`, plus `root_beta` if turned on.
- **Method:** derivative-free (a coarse grid, then a local quadratic). The responses can be rough at
  thresholds.
- **Why separate:** these keys act through soil-water drawdown that a 10-day window does not see
  (plan §13.3, §14.1).

### 7.5 Stage 5, joint polish

- **What runs:** all fitted keys, Levenberg–Marquardt from the stage values, 2–3 iterations on the
  10-day windows (plus stage 4's runs as residual blocks if used).
- **Ends with:** the final Jacobian, which gives the uncertainty (§8). Then validation and the
  five-year slow-on run (G7).

## 8. Rough uncertainty at low cost (decision 9)

1. **Laplace covariance from the final Jacobian:** Σ_post ≈ (JᵀWJ + Σ_prior⁻¹)⁻¹, no extra runs.
   Each target's block is scaled by its σ factor (§6), so a target the model cannot fit does not
   make the parameters look better known than they are.
2. **Intervals in physical units.** The Gaussian is in the transformed space (logit or log-logit),
   so intervals mapped back through the transform are asymmetric and stay inside the bounds.
3. **Per key:** MAP, 68 % and 95 % intervals, prior and posterior σ and their ratio, and the
   correlation matrix. A key at a bound is flagged with the target that pushed it (G11).
4. **Linearity check** (existing): ±1σ along the 3 leading directions, 6 runs per window, about
   50 trials. It marks the covariance "local only" where the quadratic fails.
5. **Filter sensitivity at zero model cost.** The MAP's response to a different filter, such as GPP
   u\* ≥ 0.4 instead of 0.3, is estimated from the same Jacobian by one linear update,
   Δu ≈ (JᵀWJ + Σ_prior⁻¹)⁻¹ JᵀW Δr (the rows added or removed). A full refit is needed only where
   the predicted shift exceeds the posterior σ (G10).
6. **Not in scope:** MCMC, ensemble smoothers, emulators. That is the later PEcAn work. Stage 2 is
   cheap enough for a quick MCMC on `vcmax25` alone, which can be an option for checking the Laplace
   interval.

## 9. Phases

| phase | work | main files | size |
|---|---|---|---|
| R1 | filters as per-target settings (§3); `site_reference.toml`; filter report; night NEE off; GPP σ; effective-sample weights in the fit | `scripts/calibrate_fast/tower.py`, `residuals.py`, `fit.py`, `calibrate_fast.py` | small |
| R2 | key list and priors as user settings (§4.1); registry states (fit / optional / fixed); meta-analysis priors entered with sources; screening report-only by default | `scripts/calibrate_fast/registry.py`, `parameters.toml` | small |
| R3 | model and API: per-cohort leaf drivers as hourly outputs; `meds.canopy.radiation`; a batch `meds.plant.leaf` call; tests | `src/fast_dynamics/`, `src/io/meds_output_registry.f90`, `python/meds/` | medium |
| R4 | stages 1–5 (§7) with `--stage`, and stage 2's outer passes | `scripts/calibrate_fast/` | medium |
| R5 | rough uncertainty (§8): σ-scaled Laplace, transformed intervals, filter sensitivity | `scripts/calibrate_fast/fit.py` | small |
| R6 | the BCI refit on the merged #340/#341 base; README and figures. Only when the owner asks | `examples/example04_column_biophysics/` | small |

Each phase comes with tests:
- the synthetic-model checks in `scripts/calibrate_fast/tests`;
- G8 as a ctest;
- the binding tests beside `python_api`;
- a test that every site setting is documented in `site_reference.toml`.

## 10. Gates

G1–G7 stay (plan §8). Changed or new:

| gate | requirement |
|---|---|
| G3 | (changed) every fitted key has a non-zero, smooth Jacobian column; a user-requested key that fails is reported, not dropped, unless `screening = "drop"` |
| G8 | the canopy of kernels reproduces the full run's hourly GPP within 1 % (median), same parameters and drivers |
| G10 | each key's predicted MAP shift under the alternative filter (§8.5) is below its posterior σ, or a full refit with that filter is reported |
| G11 | a key at a bound is reported with the target that pushed it there and its stage |

## 11. Cost

Estimated from the measured trial times: 6.6 s on an idle core and ~17 s on a loaded one for a
10-day window.

| stage | evaluations | core-hours |
|---|---|---|
| 1 optics | ~10³ two-stream sweeps of the record, seconds each | < 0.5 |
| 2 canopy of kernels, 3 passes | 3 × 8 frozen windows + ~10³ kernel sums | ~1 |
| 3 energy and aerodynamics | ~5 × 100 trials | ~2.5 |
| 4 water stress | ~30 seasonal runs | ~0.5 |
| 5 polish + linearity check | ~3 × 200 + 50 trials | ~3 |
| **total** | | **~7–8, against 66 for set B and 151 for the 28-key fit A** |

## 12. Next

The owner's decisions on §13.3: the optics stage's keys, the water stage's targets, the polish's
cost, and whether the BCI example ships the staged fit's set. Then the `vcmax25` prior (Kattge et
al. 2009), and the model bug of §13.3 item 5.

## 13. Implementation (2026-10-02)

R1–R5 are implemented (#345), and R6 has run at BCI. The example's shipped calibrated set is not
replaced: the findings in §13.3 need the owner's decisions first.

### 13.1 Where the implementation differs from this plan

1. **Stages 1 and 2 score the calibration windows' hours, not the whole record.** The kernels need
   each hour's per-cohort drivers from a full run. The windows' driver trials provide them from the
   same state chains; the whole record would need a five-year run with hourly per-cohort output.
   After the filters, stage 2 has 745 GPP hours in the 8 windows.
2. **The `vcmax25` prior is not entered.** Kattge et al. (2009)'s PFT table could not be read here.
   Until it is, the prior is the range (its ±2 sd band). For reference, CLM4.5 uses 55 for tropical
   broadleaf evergreen trees, set above Kattge's value on purpose (Bonan et al. 2012).
3. **`--stages a,b` (a list) and `--resume`**, in place of `--stage`.
4. **The u\* plateau test runs within four PAR classes.** Pooled over PAR, the calm hours are the
   dim ones, whose GPP/PAR is high, and the test returned 0. Within PAR classes BCI's threshold is
   0.5 (0.35 in the dimmest class), as §3.3 found hour by hour.
5. **The polish runs up to 10 iterations** (`rtol` ends it sooner), not 2–3. At BCI, 3 iterations
   left it still falling 4 % per iteration, a Gauss–Newton step of up to 19 posterior sd short of
   its optimum; 9 more iterations converged it (step left ≤ 1.2 sd).
6. **The uncertainty diagnostics** (§8):
   - The filter sensitivity subtracts the fit's own Gauss–Newton step at the MAP; §8.5's formula
     assumes a converged fit.
   - Every fit reports that step, per key in posterior sd, as its convergence.
   - The linearity check judges the mean of its two sides (the curvature) and reports half their
     difference (the slope left).
7. **The model's existing outputs do not move:** with the new outputs off, a 60-day BCI run with
   FAST output is bit-identical to the parent commit's.

### 13.2 R6: the BCI run

Setup:
- variant `interception_on`, on the onset build (#341);
- 8 calibration and 8 validation windows, plus the 2016 and 2017 dry seasons (120 days each);
- 2 nodes.

It ran twice. The first run stopped the polish at the plan's 3 iterations; a `--resume` continued
it to convergence. Fit D is the last joint fit (13 keys, the v0.3.2 build, before the onset),
shipped with #340.

**Keys** (start = the prior centres):

| key | stage | start | after its stage | MAP | 95 % interval | sd ratio | fit D |
|---|---|---|---|---|---|---|---|
| `leaf_reflect_nir` | optics | 0.45 | 0.332 | 0.308 | 0.304 to 0.317 | 0.20 | 0.349 |
| `leaf_clumping` | optics | 0.8 | 0.518 | 0.997 | 0.923 to 1 | 0.93 | 1 |
| `leaf_angle_mean` | optics | 45 | 63.6 | 62.5 | 59.2 to 64 | 0.27 | 59.5 |
| `vcmax25` | photosynthesis | 45 | 25.1 | 25.1 | 25 to 25.5 | 0.55 | 25 |
| `stomatal_g1` | energy | 3.77 | 3.77 | 5.06 | 4.32 to 5.77 | 0.27 | 5.97 |
| `stomatal_g0` | energy | 0.01 | 0.0201 | 0.0105 | 0.00877 to 0.0125 | 0.05 | 0.0336 |
| `z0m_ratio` | energy | 0.13 | 0.0522 | 0.0522 | 0.0504 to 0.0624 | 0.51 | 0.0561 |
| `leaf_width` | energy | 0.04 | 0.087 | 0.132 | 0.0742 to 0.147 | 0.58 | 0.0474 |
| `dsl_dmax` | energy | 0.015 | 0.0067 | 0.0406 | 0.0138 to 0.0491 | 0.71 | 0.0479 |
| `dewmx` | energy | 0.1 | 0.0687 | 0.0536 | 0.0504 to 0.0826 | 0.63 | 0.0754 |
| `wstress_sref_stomata` | water | 2 | 4.55 | 0.514 | 0.5 to 0.898 | 0.94 | — |
| `stomata_psi_onset` | water | -0.857 | -0.0719 | -0.114 | -1.05 to -0.00913 | 0.71 | — |

**Fit, on the calibration windows** (normalized RMSE per target, start → MAP; fit D's σ for GPP
was 1.5 + 0.15 GPP, so its GPP column is not comparable):

| target | staged fit | fit D |
|---|---|---|
| albedo | 5.99 → 1.29 | 5.99 → 1.39 |
| upwelling longwave | 2.19 → 1.81 | 2.02 → 1.87 |
| net radiation | 1.53 → 1.12 | 1.56 → 1.09 |
| LE | 1.70 → 1.52 | 1.74 → 1.76 |
| H | 2.28 → 2.50 | 2.43 → 2.01 |
| evaporative fraction | 3.10 → 2.39 | 3.75 → 2.90 |
| GPP | 2.21 → 0.98 | 2.61 → 1.43 |
| u\* | 1.99 → 0.66 | 1.99 → 0.65 |

**Validation:** Φ 59,067 at the default → 29,687 at the MAP. G4 passes: every target is better
except H (+6.1 %) and the evaporative fraction (+7.3 %), both under the 10 % limit.

**Gates:**
- **G3, G5, G8, G11 pass.** G8 is the kernel against the model's GPP: median 0.43 % over 745
  hours, mean bias +2.0 %.
- **G10 fails.** With GPP u\* ≥ 0.3 in place of 0.4, the linear update moves `stomatal_g0`
  +20 sd, `leaf_clumping` +19 sd (toward its bound) and `stomatal_g1` +4 sd. A shift this large is
  outside the linear regime (below), so it needs the full refit with 0.3 that §8.5 calls for.
- **The covariance is "local only".** Along the three leading directions, the objective's curvature
  is 3–9× the Gauss–Newton quadratic's. The reported intervals are therefore wider than the curvature
  along those directions implies.

**Five years** (the converged MAP, `cal` run of the example):

| five years | default (onset build) | fit D (the #340 build, before the onset) | staged fit | tower |
|---|---|---|---|---|
| GPP mean (RMSE), µmol m⁻² s⁻¹ | 11.61 (7.24) | 7.49 (3.56) | 7.52 (3.35) | 7.46 |
| LE mean (RMSE), W m⁻² | 71.79 (38.91) | 81.86 (43.33) | 71.37 (38.95) | 75.47 |
| H mean (RMSE), W m⁻² | 66.79 (50.41) | 67.77 (52.27) | 69.84 (52.43) | 32.43 |
| net radiation mean (RMSE), W m⁻² | 121.28 (44.68) | 135.52 (22.71) | 133.49 (24.60) | 136.32 |
| April 2014 GPP | 9.73 | 5.11 | 6.62 | 6.93 |
| April 2016 GPP | 7.29 | 3.66 | 4.98 | 6.26 |
| April 2017 GPP | 10.20 | 5.26 | 6.96 | 7.03 |
| albedo (days with sw_in > 200) | 0.179 | 0.121 | 0.123 | 0.129 |
| whole-site budget breaches (energy, water) | 0, 0 | 0, 0 | 0, 0 | |

**Cost:**
- 4,804 trials over the two runs, against fit D's 10,224.
- 31.8 core-hours of trials: 4,264 ten-day trials at 15 s, and 858 seasonal ones at 57 s. Fit D's
  were about 43.
- Wall time: 71 min on 2 nodes, against 36. The chains, the refresh and the seasonal runs are
  serial, and each polish iteration waits for its 120-day trials.
- §11's estimate of 7–8 core-hours was wrong. The stages cost under 10 minutes in all; the cost
  is the polish, which moves all 12 keys over 10 windows for 12 iterations.

### 13.3 Findings for the owner

1. **The albedo alone cannot place clumping or leaf angle.**
   - Stage 1 sent clumping to its floor (0.52) and the leaf angle to 63.6°; the polish then moved
     clumping to its ceiling (0.997), where fit D also put it.
   - At the stage-1 result, clumping's Jacobian column correlates 0.986 with the NIR reflectance's.
   - The model's albedo at the start is 0.186 against the tower's 0.128 on the same hours. The NIR
     reflectance alone, at its floor (0.30), still leaves the model too bright.

   Options:
   - fit only `leaf_reflect_nir` in stage 1, and move clumping and angle to the energy stage;
   - or give clumping a prior (e.g. 0.9 ± 0.05, from the understory light of §1).
2. **Every key that lowers GPP ends at a bound.** `vcmax25` sits at its floor (25.1; fit D 25.0).
   The water stage first pushed the stress to its strongest (`wstress_sref_stomata` 4.55,
   `stomata_psi_onset` −0.07 MPa); the polish then moved the sensitivity to its floor (0.51), with
   the onset at −0.11 MPa. Seven of the 12 keys end at or within a few percent of a bound. The GPP
   gap of §1 (model 1.5× the tower at the defaults) is still the dominant signal; the larger σ and
   the u\* filter did not remove it.
3. **The water stage scores GPP.** Its keys then answer to GPP's level as much as to the drought.
   Option: score the seasonal runs on LE and the evaporative fraction only.
4. **The polish is the cost.** Options:
   - polish without the seasonal runs (`include_water = false`), with the water keys held at their
     stage values;
   - or polish only the keys whose stage the coupling affects.
5. **A model bug, found by three failed trials.** In the 2017 dry season, near the MAP
   (`stomata_psi_onset` −0.09 to −0.15, `dsl_dmax` ~0.04, interception on), the fast loop blew up
   at a mid-afternoon drop in shortwave (2017-04-17 18 UTC): H rose to 631 W m⁻² under 12 W m⁻² of
   incoming shortwave. Then u\* (until the run's end) and GPP (30 hours) were NaN. The run still
   ended "no NaNs", and only the whole-site budget check (a NaN cumulative residual) caught it.
   Not yet filed.

### 13.4 After R6: the owner's decisions and two quick fits (2026-10-02)

**Decisions.**
1. **Observations.**
   - No closure correction: LE and H are used as measured. At BCI, H + LE = 0.72 Rnet, and the gap
     behaves like missing sensible heat. As u\* rises from 0.3 to 1.0 (10–15 h, within a VPD class),
     H/Rnet climbs from 0.21 to 0.38, while LE/Rnet stays at 0.42–0.48.
   - LE is kept at u\* ≥ 0.4. H is kept at u\* ≥ 0.6 and 9–16 h only, with σ 10 W m⁻² + 30 %.
   - Net radiation is off: the four-component sum repeats the albedo and the upwelling longwave.
     The evaporative fraction is off: it is biased high where closure is poor. Soil water is not
     used.
2. **Stages.** The optics and photosynthesis keys are fitted in the coupled `energy` stage: the
   default is `energy`, `water`, `polish` (finding 1). The water stage scores LE and GPP.
3. **Keys.**
   - Fixed:
     - `leaf_clumping` at 0.80 (the base value; FATES uses 0.85);
     - `leaf_width`;
     - `dsl_dmax`;
     - the film capacities. The tower does not see the wet canopy: FLAG gap-fills every raining
       hour.
   - Fitted:
     - `theta_j`, prior 0.80 ± 0.05;
     - `jmax_vcmax_ratio`, prior 1.70 ± 0.15;
     - `ds_vcmax`, prior 641 ± 5;
     - `ds_jmax`, prior 640 ± 4.

     The last three priors are Kattge & Knorr (2007) acclimated at BCI's 25.5 °C. That matches
     Slot & Winter (2017): four Panama species have Vcmax optima at 32.9–39.7 °C. The model's
     default `ds_vcmax` (650) peaks at 31.7 °C.
   - Optional: `ea_vcmax`, `ea_jmax`.

**Two quick fits.** These are the `energy` stage alone, over the eight calibration windows: about
2,100 trials and 24 min on 2 nodes each.
- Fit 2: the combined stage, before the key decisions.
- Fit 3: with `leaf_clumping` and `leaf_width` fixed and the photosynthesis keys fitted.
  `dsl_dmax` was still fitted.

| key | fit 2 | fit 3 | bound |
|---|---|---|---|
| `vcmax25` | 25.1 | 25.3 | 25 |
| `theta_j` | — | 0.73 | 0.70 |
| `jmax_vcmax_ratio` | — | 1.45 | 1.4 |
| `ds_vcmax` | — | 639.4 (prior 641) | |
| `ds_jmax` | — | 640.3 (prior 640) | |
| `stomatal_g1` | 3.46 | 4.11 | |
| `leaf_clumping` | 0.51 | 0.80 (fixed) | 0.50 |
| `leaf_width` | 0.143 | 0.04 (fixed) | 0.15 |

| model GPP / tower GPP (calibration hours) | all hours | 6–10 h | 10–14 h |
|---|---|---|---|
| fit 3 at the MAP | 1.02 | 0.99 | 0.99 |
| fit 3 with `vcmax25` = 45 | 1.44 | 1.35 | 1.44 |

**Findings.**
- `vcmax25` stays at its floor. `theta_j` and Jmax/Vcmax go to their GPP-lowering floors too. Doing
  so they remove the morning excess, but what remains is a flat offset: at a physiological
  `vcmax25` of 45 the model is 44 % above the tower at every hour.
- The temperature terms stay at their priors, and they cannot absorb the offset:
  - the GPP hours' air temperature spans only 26.0–29.1 °C (10th–90th percentile);
  - within each light class, model/tower GPP falls by 3.8 % per °C, so the data want a steeper
    response, not a weaker one.
- Fit 3 fails G4 on H (+18 % on validation). Fit 2 passes it, with H +3 %.
- Splitting cohorts into thinner layers (`cohort_lai_cap` 0.1) lowers GPP by about 2 %. That test
  exposed the leaf-PAR floor `max(LAI, 0.1)` (#346), fixed in #347 (+0.33 % GPP at BCI).
- Still open for the owner:
  - how to treat `vcmax25`: a literature prior, or an effective value that absorbs the offset;
  - sunlit and shaded leaves (#343);
  - finding 5.


## Sources

- BCI flux data and processing notes: Detto, Dryad, doi:10.5061/dryad.3tx95x6j5
- ED2 source: github.com/EDmodel/ED2, `ED/src/dynamics/photosyn_driv.f90`
- FATES source: github.com/NGEET/fates, `radiation/FatesNormanRadMod.F90`, `biogeophys/FatesPlantRespPhotosynthMod.F90`
- Björkman & Demmig (1987), Planta 170, 489–504 (PSII photon yield)
- Chazdon & Fetcher (1984), Journal of Ecology 72, 553–564 (understory light)
- Ehleringer & Björkman (1977), Plant Physiology 59, 86–90 (C3 quantum yield)
- Bonan et al. (2012), Journal of Geophysical Research 117, G02026 (CLM4.5's tropical Vcmax25)
- Kattge et al. (2009), Global Change Biology 15, 976–991 (Vcmax by PFT)
- Kattge & Knorr (2007), Plant, Cell & Environment 30, 1176–1190 (thermal acclimation of Vcmax and Jmax)
- LeBauer et al. (2013), Ecological Monographs 83, 133–154 (PEcAn trait meta-analysis)
- Lin et al. (2015), Nature Climate Change 5, 459–464 (g1 by PFT; Fig. 2e)
- Papale et al. (2006), Biogeosciences 3, 571–583 (u* filtering and flux processing)
- Reichstein et al. (2005), Global Change Biology 11, 1424–1439 (u* threshold, flux partitioning)
- Medlyn et al. (2002), Plant, Cell & Environment 25, 1167–1179 (Jmax/Vcmax)
- Sabot et al. (2022), Plant, Cell & Environment (stomatal water stress)
- Slot & Winter (2017), Plant, Cell & Environment 40, 3055–3068 (temperature response of Panama canopy leaves)
- The BCI diagnosis (2026-10-01/02): its numbers are reproduced in §1, §3.3 and §6; its scripts
  are not part of the repository.
- Deferred issues: #342 (Γ* pressure scaling), #343 (cohort-level sunlit/shaded leaves)
