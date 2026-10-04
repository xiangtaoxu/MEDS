# MEDS fast calibration: a cross-site protocol

**Status:** revision 4, 2026-10-03. The owner and Claude reviewed it step by step, and every
decision below is the owner's (§0). Not yet implemented; the phases are in §9.
- **What it replaces:** where the two plans differ from it, `MEDS_FAST_CALIBRATION_PLAN.md` (the
  method, the first fits) and `MEDS_FAST_CALIBRATION_REVISION_PLAN.md` (the staged fit). Their
  history stays in them.
- **Shipped beside the review:** PR #351 (`phi_psii` 0.74, θ_J 0.7; merged into beta) and issue
  #350 (the forcing's reference height).

**Scope:**
- the fast parameters, a frozen stand, flux towers one at a time, a rough uncertainty;
- no remote-sensing target;
- a full Bayesian analysis stays with PEcAn.

**Design rule:** rules that work at any tower, not BCI's numbers.
- **Site-specific code** lives only in the adapter that reads a tower's file (§2.1).
- **A tower's numbers** (its u\* threshold, σ, closure, windows, priors) come from its own data and
  climate by the rules below. The site file can override any of them.
- **Code stays simple.** No fallback methods, and one fit (§6).

## 0. Owner decisions (2026-10-03)

| step | decision | § |
|---|---|---|
| layers | one site TOML per tower (the forcing builder's), one reader (`tower_inputs.py`); facts in the site file, choices in `calibration.toml` | 2.1 |
| metadata | checks in the adapter; the forcing build stops on a mismatch, the calibration reports | 2.1 |
| hours used | measured hours only; observed forcing only in calibration windows | 2.2 |
| u\* | per target, by its own diagnostic; for CO₂ the provider's threshold by default, the daytime plateau as the alternative | 2.3 |
| transitions | cut hours only for a measurement or definition reason; no transition cut for H and LE; every tower flux treated as turbulent (storage ignored) | 2.4 |
| heights | fluxes compared as they are; the forcing's neutral move to each patch's canopy-air top is #350, no model change now | 2.4 |
| closure | an observation model with the shortfall measured over whole days; shares from the attribution test; Bowen as the alternative | 3.3 |
| targets | albedo, upwelling longwave, H, LE, GPP, u\*; H, LE and u\* daytime only; keys fixed when the kept data never sample their process; no optional targets for now; NEE and respiration later | 2.5 |
| interval | the tower's native interval (30 min at BCI) | 2.6 |
| random error | provider's uncertainty, then paired days, then defaults; σ at a smoothed observation; Huber loss | 3.1 |
| autocorrelation | effective-sample weights as now | 3.2 |
| GPP | observation model: NEE as measured, the tower's respiration × κ; κ prior centre 1 at other sites, 0.65 ± 0.10 at BCI | 3.3 |
| model error | σ scaled once by model misfit, never below the measurement σ, at most 3× | 3.4 |
| priors | by plant type; EEO centres with a wide range (log-sd 0.5); bounds kept apart from priors | 3.5 |
| leaf light use | `phi_psii` 0.74 (the low-light electron yield), θ_J 0.7 (mean-light cohorts), both fixed (#351) | 3.6 |
| chains | frozen (`slow_on = false`) until the slow processes are split; the site declares its leaf-on months | 6.1 |
| windows | 8 ten-day windows spread over the year; validation in other years; seasonal runs from a water-deficit index | 6.2 |
| keys | the default set of §4.3, with kinds and scope tags | 4 |
| the fit | one joint fit of all keys; the separate water stage, the polish and the kernel stages go | 6.3 |
| uncertainty | Laplace, declared alternatives (refit above 1 sd), structural variants | 7.1 |
| acceptance | gates, validation, full record, diagnosis of conflicts, labels | 7.2 |

## 1. What the fits taught

| # | lesson | evidence | rule |
|---|---|---|---|
| L1 | Calibration finds model bugs before it finds parameters. | H reported ~105 W m⁻² low (#328); the stomatal latch (#332 → #335); gravity head read as drought (#341); leaf PAR / max(LAI, 0.1) (#346 → #347); undeclared water floors (#333 → #336). | Fix first; a key at a bound is a question, not an answer (§5.1). |
| L2 | A rough key is a switch in the model. | `leaf_pi0`, `wood_psi50` under the latch; their line-searched values broke the five-year water budget (659 breaches). | Make the model continuous. |
| L3 | Filter choices move the answer more than noise does. | GPP u\* 0.4 → 0.3 moved `stomatal_g0` +20 posterior sd. | Declared alternatives (§7.1). |
| L4 | Find where the closure gap goes before correcting it. | At BCI the gap behaves like missing H: H/Rnet rises with u\*, LE/Rnet does not. | The attribution test (§3.3). |
| L5 | Some targets repeat others. | Rnet = the four components; EF repeats H's bias. | Redundant targets stay off. |
| L6 | Tower GPP carries a shared error that a wider σ does not handle. | One respiration value per day, about the soil chambers' alone; `vcmax25` still at its floor with σ_abs 2.5. | The GPP observation model with κ (§3.3). |
| L7 | Calm air makes morning GPP read low. | 7 h GPP/PAR +26 % at u\* ≥ 0.5; the plateau test works only within PAR classes. | The u\* diagnostic (§2.3). |
| L8 | Keys whose process the kept data never sample can't be fitted. | FLAG removes every rain half-hour: film capacities unconstrained. | Fix from coverage (§2.5). |
| L9 | Fewer keys, better long runs. | 11 keys lost 4 % on validation against 28, at half the cost, and did better over five years. | Parsimony; the full record decides. |
| L10 | Ten days can't set keys that act through soil drying. | Frozen windows reward values that dry the soil months later. | Seasonal runs in the fit (§6.3). |
| L11 | A stage works only if its target separates its keys. | The optics stage sent clumping to its floor, the polish to its ceiling. | No kernel stages (§6.3). |
| L12 | A stage scored on a biased target takes on the bias. | The water stage scored on GPP set the strongest stress. | Seasonal runs scored on LE only. |
| L13 | Shape keys absorb a level error. | θ_J, Jmax/Vcmax and `vcmax25` all at their floors in fit 3. | Shape keys fixed (§4.3). |
| L14 | The Laplace covariance is only local. | Curvature 3–9× the quadratic's. | Three uncertainties side by side (§7.1). |
| L15 | The harness must prove each trial. | glibc math in Python, refused keys, stale states, idle-node timeouts. | Parameter record, bit-identity, digests (in place). |
| L16 | A trial can fail without saying so. | 2017 dry season: H 631 W m⁻² under 12 W m⁻² SW, then NaN; "no NaNs" printed. | G12; file the bug. |
| L17 | The 120-day runs inside the polish were the cost. | 31.8 core-hours, mostly the polish. | Fewer keys, cheaper iterations (§6.3). |
| L18 | A paraphrased data note gets copied as fact. | "40–50 % below soil chambers" in four files; the provider wrote something else (§3.3). | Quote the provider in the site file; correct the four places (P0). |

## 2. The tower data

### 2.1 One site declaration, one reader

- **The site TOML:** the forcing builder's (`scripts/prepare_flux_tower`), declaring the input format, clock, heights and variables. It gains:
  - `[fluxes]`: SW_out, LW_out, Rnet, LE, H, NEE, GPP, RECO, u\*, PAR, each with column, units and its own "measured" rule (BCI: `FLAG = 1`; FLUXNET: `_QC = 0`);
  - `[provider]`: how GPP was made, with the provider's words quoted; the respiration column (BCI: gpp + NEE; FLUXNET: `RECO_NT_VUT_REF`); the published u\* threshold; any random-uncertainty columns; what the provider's quality flag screened on.
- **The reader:** `tower_inputs.py` (csv, AmeriFlux BASE, FLUXNET/ONEFlux). Its output is one standard table: UTC, start-stamped, SI units, a measured mask per variable. `calibrate_fast`'s own reader goes.
- **What goes where:** facts about the data in the site TOML; choices about the fit in `calibration.toml`, whose `[tower]` block becomes a pointer to the site TOML.
- **Metadata checks:** in the adapter, beside V1–V5, extended to the fluxes:
  - Rnet against the sum of its components;
  - the sign of LW_out − σT⁴, which catches a swapped column;
  - SW_out and LW_out bounds.

  On a mismatch the forcing build stops; the calibration reports.

### 2.2 Which hours

- **Measured hours only,** by each variable's mask. Gap-filled values are another model's prediction ("a model–model comparison", Fer et al. 2018).
- **Observed forcing only** in calibration windows. Validation may relax this.
- **Daytime:** SW_in above 10 W m⁻².

### 2.3 u\*: a diagnostic per target

The usual u\* filtering is for night-time NEE. REddyProc filters night only by default, and ONEFlux applies its thresholds to day and night NEE. For other targets and hours, the tool decides from each target's own data:
1. **The diagnostic:** within classes of the target's driver (PAR for GPP or NEE; Rnet for LE and H), the flux-to-driver ratio across u\* classes. The plateau test finds the lowest class within 95 % of the mean of the classes above (Papale et al. 2006), with a bootstrap for its spread (Barr et al. 2013).
2. **The outcome:**
   - **a plateau:** filter;
   - **a flat ratio:** no daytime filter;
   - **still rising at the top classes:** no filter; it is the target's observation model's job (H at BCI, §3.3).
3. **For CO₂ the provider's threshold is the default** (BCI 0.4; the ONEFlux convention); the daytime plateau (BCI 0.5) is the declared alternative.
4. **Fallback** when a diagnostic has too few hours: the night threshold, on CO₂ only.
5. **Always:**
   - the tower's measured u\*, never the model's;
   - never on the u\* target itself;
   - per-year thresholds where the provider gives them.

At BCI: GPP filtered at 0.4; LE not filtered (LE/Rnet flat); H not filtered (handled by §3.3).

### 2.4 Transitions, storage and heights

- **Hours are cut only for a measurement reason or a definition mismatch,** never because the model fits badly there. That misfit is reported (residual by hour), as model structure.
- **Every tower flux is treated as a turbulent flux.** Storage is ignored, being small at long time scales. So:
  - no transition cut for H and LE (BCI's 9–16 h window for H goes);
  - no CO₂ sunrise cut.
- **Heights:**
  - **The fluxes are comparable as they are:** in a constant-flux layer, a patch's flux at its canopy-air top equals its contribution at the sensor.
  - **The air that drives them is not:** MEDS moves the forcing to each patch's canopy-air top along a neutral profile. At BCI (patch tops 15–51 m, area-weighted mean 38.3 m, sensor 41 m) this gives midday H errors of −84 to +100 W m⁻² per patch, which net to +1.5 to +7.5 W m⁻² for the site. Where the sensor sits far from the patch tops (a grassland tower below a 5 m minimum canopy-air depth) the error has one sign.
  - **Filed as #350;** no model change now. The data report prints the share of area whose canopy-air top is above the sensor.

### 2.5 Targets

| target | observation model | hours |
|---|---|---|
| albedo (SW_out / SW_in) | identity | SW_in > 200 W m⁻², the sun at least 20° high, no snow |
| upwelling longwave | identity | day and night |
| H, LE | closure model (§3.3) | daytime |
| GPP | respiration model (§3.3) | daytime, CO₂ u\* filter |
| u\* | identity, loose σ (roughness sublayer) | daytime |

- **Not targets:**
  - Rnet: an input to the closure model only;
  - the evaporative fraction: repeats H and LE;
  - night NEE and respiration: a later version, once the soil carbon is developed;
  - soil water: a relative index;
  - optional targets (reflected PAR, canopy temperature): not for now.
- **Night H, LE and u\*** are set mostly by numerical floors (`ustmin`, `canopy_freeboard`) and stable-air measurement problems, so they're off by default, with night as an option at the provider's threshold.
- **Keys fixed from coverage:** the report counts, for each process (wet canopy, drought, night, snow), the kept hours that sample it. A key whose process has almost none is fixed, with the count as its reason. At BCI the film capacities are fixed: the flag removes every rain half-hour.

### 2.6 Time interval

The fit runs at the tower's native interval: 30 minutes at BCI and FLUXNET's half-hourly sites, hourly at "HR" sites.
- **The model's output** is written at the same interval (`fast_interval_steps = 2` at a 15-minute step), and the tool checks the match.
- **The old hourly pairing goes:** it kept an hour only when both of its half-hours were measured, and dropped 5.8 % of BCI's measured half-hours (12.3 % of the daytime ones).

## 3. Error assumptions

### 3.1 Random error

- **Form:** σ = σ_abs + σ_rel·|flux| per target (Richardson et al. 2006).
- **Source, in order:**
  1. the provider's per-half-hour random uncertainty (FLUXNET2015 `*_RANDUNC`);
  2. the tool's paired-day estimate (Hollinger & Richardson 2005): the same half-hour on neighbouring days with similar radiation, temperature, VPD and wind, binned by flux size;
  3. defaults.

  Albedo and upwelling longwave keep instrument-level defaults (0.01; 5 W m⁻²).
- **σ is evaluated at a smoothed observation:** the mean of the measured values at the same time of day within ±7 days.
  - **Why:** weighting each row by σ(its own observation) gives randomly low values more weight.
  - **Its size:** a simulation put the fitted GPP 6–14 % low (7 % at GPP 10 with BCI's old 2.5 + 0.15·GPP).
- **Loss:** Huber (c = 2), since half-hourly errors are heavy-tailed (Lasslop et al. 2008; PEcAn uses a Laplace likelihood). Least squares stays an option. The χ² report uses raw residuals.

### 3.2 Autocorrelation

- **The weights:** n_eff/n for each target and window, from the lag-1 autocorrelation of its residuals.
- **When they are set:** at the start and at the mid-fit refresh (§6.3), fixed in between. That avoids rewarding a worse fit during the search.

### 3.3 Systematic errors: observation models

Each target's observation model says how the tower's number relates to the true flux, with each known bias as its own term and prior.

**GPP (the owner's proposal):**

```
NEE_obs,h   = NEE_true,h + ε_h            identity; ε is NEE's random error
R_tower,d   = κ · R_true,d                one multiplicative error on the tower's respiration
GPP_tower,h = R_tower,d − NEE_obs,h       the provider's definition

r_h = [ GPP_model,h − GPP_tower,h − (1/κ − 1) · R_tower,d ] / σ_NEE,h
```

- **Closed with the tower's own respiration,** not MEDS's, whose soil carbon is frozen and which has no growth respiration. This is the same as fitting daytime NEE.
- **κ is one global observation key.** Its gradient is known exactly. It is fitted with a prior and never alongside free shape keys.
- **The prior:**
  - **at other sites:** centre 1, with sd from the gap between the provider's two partitionings (FLUXNET `RECO_NT` against `RECO_DT`), else 0.15. The night-time partitioning is the default respiration, the daytime one the alternative.
  - **at BCI: 0.65 ± 0.10, bounded 0.4–1.0.** The tower's respiration (4.09 µmol m⁻² s⁻¹) is about the soil chambers' alone (1,613 gC m⁻² yr⁻¹ ≈ 4.3; Rubio & Detto 2017). The provider's note says it "appeared underestimated … considering that RECO includes also above ground respiration which can contribute up to 40-50% of total respiration". The correction runs 1.7 µmol m⁻² s⁻¹ in the dry season and 2.7 in May.

**H and LE (closure):**

```
f_d       = median over ±15 days of daily (Rnet − G) / (H + LE)      days with ≥ 70 % measured hours
H_true,h  = H_obs,h  + s_H  · (f_d − 1) · (H_obs,h + LE_obs,h)
LE_true,h = LE_obs,h + s_LE · (f_d − 1) · (H_obs,h + LE_obs,h)       s_H + s_LE = 1
```

- **Whole days:** over a day, ground heat and storage roughly net out, so no estimate of them is needed. G = 0 where not measured. The provider's gap-filled values fill the daily sums for this estimate only.
- **The shares come from the attribution test:** within VPD classes, which of H/Rnet and LE/Rnet rises with u\*. One rising takes the share, both rising means Bowen, neither means as measured. **Bowen is the declared alternative** (the FLUXNET convention).
- **The assumption:** the tower misses the same fraction of turbulence at every hour.
- **At BCI:** s_H = 1; f ≈ 1.33 (daily closure median 0.75); midday H target 287 W m⁻² against 160 as measured; LE unchanged.
- **This replaces** the `closure` flags, H's 30 % σ and BCI's H filters.

### 3.4 Model error

- **When:** once, at the mid-fit refresh (§6.3).
- **What:** each target's σ is multiplied by max(1, √(χ²/n)), never below the measurement σ and at most 3×, and the scale is reported. Targets the model structurally cannot fit then stop bending the keys.

PEcAn samples each target's variance; ORCHIDEE sets it to the prior model's misfit × 30 (§11).

### 3.5 Priors

- **Order of preference:**
  1. site leaf measurements;
  2. **eco-evolutionary optimality (EEO) values from the site's climate** (forcing only, never the tower's fluxes);
  3. a meta-analysis by plant type, reported beside EEO and flagged when more than 2 sd apart;
  4. another model's default;
  5. the range alone.
- **EEO centres with a wide range** (log-sd 0.5), so one rule covers plant types:
  - **`stomatal_g1`:** the least-cost value. Medlyn's g1 equals ξ = √(β(K + Γ\*)/1.6η\*), with β = 146 (Prentice et al. 2014; Stocker et al. 2020). **BCI: 2.83 kPa^0.5.**
  - **`vcmax25`:** the coordination condition (Wang et al. 2017; Smith et al. 2019), solved with **MEDS's own leaf equations**. That is, the vcmax25 at which MEDS's Rubisco-limited and light-limited rates are equal at the site's growing-season daytime climate, with the configured `phi_psii` and θ_J, converted to 25 °C by MEDS's temperature response. **BCI: 54.**
- **Jmax/Vcmax, `ds_vcmax`, `ds_jmax`:** fixed at Kattge & Knorr (2007) for the site's growth temperature, computed from the forcing, not hard-coded.
- **Plant types:** the PFT file names its class, and the registry holds a prior per class.
- **Bounds apart from priors:** the hard bound is what is physically possible, the prior what the evidence says. Each key's prior z is reported; |z| > 2 is flagged (G13).

### 3.6 Leaf light use (#351)

- **`phi_psii` 0.85 → 0.74:** the electron yield of linear transport in low light. J's initial slope is 0.5·phi_psii per absorbed photon.
  - **The problem with 0.85:** with it and 4ci + 8Γ\*, a leaf fixed CO₂ at the O₂-evolution yield (0.106 per absorbed photon; Björkman & Demmig 1987).
  - **What 0.74 gives:** measured CO₂ fixation is 0.093 (Long et al. 1993), which 0.74 gives. In normal air, 0.064 → 0.056 against the measured 0.052 ± 0.003 (Skillman 2008).
  - **The 4/8 form stays.**
- **θ_J 0.85/0.90 → 0.7** for the reference C3 types: the effective curvature of a mean-light cohort (CLM5, ED2). To revisit with sunlit/shaded leaves (#343, comment posted).
- **Both fixed:** the tower's dim-light GPP can't separate them from κ. It gives 0.046 per absorbed photon as reported, and 0.060–0.080 corrected with κ 0.75–0.55.
- **The examples** keep 0.85 and 0.90 until they are revisited (P7).

## 4. Keys

### 4.1 Seven questions per key

A key is fitted only if every answer is yes.

| # | question | if no |
|---|---|---|
| T1 | a property of the plant or surface, not of the numerics? | never calibrate |
| T2 | acts within the fitted runs (the windows, or the seasonal runs)? | fixed |
| T3 | a kept target sees it (posterior/prior σ ratio < 0.9)? | fixed at the prior |
| T4 | separable from the other keys (\|correlation\| < 0.95)? | fix one of the pair |
| T5 | a smooth response? | fix the model; until then fixed |
| T6 | the target that sees it measures that process without its own bias? | fixed, or a tight prior |
| T7 | the data move it beyond its prior? | fixed at the prior |

### 4.2 Kinds and scope

- **Kinds:**
  - trait: a measurable property, with a prior from evidence;
  - effective: a scheme property, whose value belongs to this model structure;
  - numerical: never calibrated;
  - observation: a term of a target's observation model, never written to a MEDS config.
- **Scope tags:**
  - plant type (`vcmax25`, g1), shared by all towers of that type;
  - site (`z0m_ratio`, clumping);
  - observation (κ).

  A later multi-tower fit becomes a configuration change.

### 4.3 The default set

| kind | fitted | fixed |
|---|---|---|
| optics | `leaf_reflect_nir`, `leaf_angle_mean` | NIR transmittance (moves with the reflectance), visible optics (unseen), `leaf_clumping` (the albedo can't separate it) |
| photosynthesis | `vcmax25` | `phi_psii`, θ_J (#351); Jmax/Vcmax, `ds_vcmax`, `ds_jmax` (Kattge & Knorr at growth temperature); `rd_vcmax_ratio` |
| stomata | `stomatal_g1`; `stomatal_g0` if the screening finds the daytime data inform it | |
| aerodynamics | `z0m_ratio` (effective, flagged) | `d_ratio`, `leaf_width`, `dsl_dmax`, film capacities |
| water stress | `wstress_sref_stomata`, `stomata_psi_onset` (only with a drawdown in the record) | `root_beta` (optional) |
| observation | κ | closure shares (§3.3) |

## 5. Model structure

### 5.1 The calibration's signals as a model audit

| signal | first action |
|---|---|
| a key at a bound, or \|prior z\| > 2 | find the target that pushes it; look at that target's residual by hour, light and season |
| a rough key | find the threshold in the model; make it continuous |
| a budget breach or a NaN | file and fix; never fit around it |
| a target's χ²/n ≥ 2 | report it as structure, after the data checks |
| validation windows fit, the full record doesn't | fewer keys |

### 5.2 Known structural errors

| error | status |
|---|---|
| mean-light cohorts (sunlit/shaded: −11 % GPP at BCI) | #343, deferred; θ_J set for it (§3.6) |
| understory too bright (clumping) | `leaf_clumping` fixed; set from understory light when revisited |
| frozen LAI | until the slow processes are split (§6.1) |
| neutral move of the forcing to patch tops | #350 |
| single-slab canopy air; stiff sub-canopy conductance | #269, #265 |
| Γ\* in Pa; `veg_coupling_floor` in ground units | #342, #349 |
| the 2017 dry-season blow-up | to file (P0) |

### 5.3 Structural variants

The same keys, priors and filters under the structures declared in `[variants]`: interception on and off by default, and sunlit/shaded leaves once #343 exists. The spread of each key's MAP across them is reported (§7.1).

### 5.4 The BCI GPP gap, as it now stands (rough)

| step | model/tower daytime GPP |
|---|---|
| the old defaults | 1.51 (measured) |
| θ_J 0.7 | 1.42 (measured) |
| `phi_psii` 0.74 | ~1.33 (estimated) |
| κ = 0.65 raises the tower | ~1.16 (estimated) |
| sunlit/shaded leaves (not in) | ~1.03 (estimated) |

So the fit should land `vcmax25` near its EEO value (54), not at its floor. The first BCI run will tell.

## 6. The fit

### 6.1 State chains

- **Frozen** (`slow_on = false`), from the initial stand: a census or a spun-up state, declared in the site file. They write a state at each window's start, and every trial restarts from it.
- **Phenology** runs only inside the whole slow tier. Chains stay frozen until the slow processes are split, and the site declares its leaf-on months (all months for an evergreen site).

### 6.2 Windows

- **8 calibration windows of 10 days, spread evenly over the year:** one per 1.5-month slot inside the leaf-on months, each the slot's best-covered window (measured hours, observed forcing).
- **Validation windows:** the same slots in other years, where the record has at least 2 years.
- **The report shows what the windows cover:** the share of the record's radiation, temperature, VPD and soil-water ranges they span.
- **Seasonal runs (soil drying):**
  - **Where:** up to 2 runs of 120 days, each ending at a year's deepest cumulative water deficit (rain minus Priestley–Taylor potential evaporation, from the forcing), inside the leaf-on months.
  - **Skipped** when the deepest deficit is under about 100 mm; the water keys are then fixed with that reason.
- **Nothing is hand-picked.**

### 6.3 One joint fit

1. **Screening:** the first gradient matrix, at the prior centres, with the triage report (§4.1).
2. **Levenberg–Marquardt on every fitted key at once.**
   - **Rows:** the 8 windows (all targets) plus the seasonal runs (LE only).
   - **Start:** the prior centres, one start.
   - **Iterations:** one-sided differences, with the gradient matrix reused between iterations (a Broyden update) and recomputed in full every third iteration or after a rejected step.
   - **Failures:** a failed trial rejects its step; a NaN fails the trial (G12).
3. **Once, after the first convergence:**
   - re-run the chains with the current values;
   - refresh the effective-sample weights;
   - scale σ by model misfit (§3.4);
   - continue to convergence. At most 10 iterations each time; stop when the cost falls by under 0.1 %.
4. **A final full central-difference gradient matrix** for the uncertainty. Every key, water keys included, gets an interval and its correlations.

**What goes, to keep the code simple:**
- the separate water stage and its grid search;
- the polish;
- the optics and photosynthesis kernel stages and their gate (G8);
- multiple starts as a default.

## 7. Uncertainty and acceptance

### 7.1 Three uncertainties side by side

1. **Laplace, from the final gradient matrix,** with the σ scales applied.
   - **Intervals:** mapped back through the transform, so they are asymmetric and inside the bounds.
   - **Linearity check:** along the 3 leading directions; where it fails, the covariance is marked "local only".
2. **Declared alternatives, the same at every site:**
   - GPP u\*: provider against daytime plateau;
   - closure: attribution shares against Bowen;
   - partitioning: night-time against daytime, where both exist.

   Each alternative's shift is first estimated from the final gradient matrix. Above 1 posterior sd, the fit is rerun from the MAP with that alternative, and both MAPs are reported.
3. **Structural variants** (§5.3).

### 7.2 Gates and acceptance

| gate | requirement |
|---|---|
| G1 | the stand is identical at the start and end of every trial |
| G2 | the same parameters give byte-identical output |
| G3 | every fitted key has a non-zero, smooth gradient column (reported) |
| G4 | on validation windows the MAP beats the default; no target more than 10 % worse |
| G5 | keys near a bound are reported with the target that pushed them |
| G7 | the full record with the slow tier on: closed budgets; dry-season GPP and LE no worse than the default's |
| G10 | every declared alternative's shift is under 1 posterior sd, or its refit is reported |
| G12 | no scored output has a NaN |
| G13 | every trait key with \|prior z\| > 2 is diagnosed (§5.1) or relabelled effective |

G6 (multiple starts) is only on request; G8 is removed with the kernel stages; G11 is merged into G5.

**A calibrated set ships when** G1–G13 pass, any conflict has been diagnosed, the effective keys are labelled in the calibrated files, and the report carries the three uncertainties.

**The report also gives:**
- per target: χ²/n, the σ scale, and the model/tower ratio by hour and by light class;
- per key: kind, scope, prior source and z, posterior/prior σ ratio;
- κ, with the respiration and full-record GPP it implies;
- the data report (§2): rows through each filter, the u\* diagnostics, the closure attribution, the window coverage, and the share of area above the sensor.

## 8. For a new tower

1. Write the site TOML (format, clock, heights, variables, `[fluxes]`, `[provider]` with quotes, leaf-on months, the initial stand) and build the forcing.
2. `report`: read the metadata checks, u\* diagnostics, closure attribution, σ estimates and window coverage. Override a rule only with a reason.
3. `check`: trial integrity (G1, G2, the parameter record).
4. `fit`: screening, then the joint fit.
5. Read the report. Any key at a bound or with \|z\| > 2 goes through §5.1.
6. Run the full record with the slow tier on; ship with the labels.

## 9. Work

| phase | work | size |
|---|---|---|
| P0 | correct the misread respiration note (revision plan §1, §6; `site_reference.toml`; BCI `calibration.toml`); file the 2017 dry-season blow-up; G12 | small |
| P1 | the adapter: one site TOML through `tower_inputs.py` (`[fluxes]`, `[provider]`, measured masks, metadata checks); the native interval; `calibrate_fast`'s reader removed; BCI's declarations moved into `bci_site.toml` | medium |
| P2 | data rules: the u\* diagnostic with bootstrap; daytime-only turbulent targets; albedo rules; coverage-based fixing; window selection and coverage report; the water-deficit index and seasonal runs | medium |
| P3 | errors: `obs_model` per target (identity, closure, respiration with κ); σ sources and smoothed evaluation; Huber by default; σ scaling at the refresh | medium |
| P4 | keys and priors: kinds, scope tags, plant-type priors, bounds apart, prior z and G13; EEO g1 and vcmax25 (MEDS's leaf); Kattge & Knorr from growth temperature | small–medium |
| P5 | the joint fit with one refresh; one-sided differences and gradient reuse; the final central matrix; **removal** of the water stage, polish, kernel stages, G8 and default multi-start | medium; removes code |
| P6 | uncertainty: declared alternatives with the refit rule, structural variants, the report | small |
| P7 | BCI: the examples with #351's values, the refit, the full-record run, the README | small; Slurm |

P0 first; P1–P3 before P5; P7 last. Each phase adds or changes settings with documented defaults, and the tests change with it.

## 10. Later

- **Splitting the slow processes,** so phenology can run in chains (owner).
- **Sunlit/shaded leaves (#343),** with θ_J back toward a leaf value; and clumping from the understory light.
- **#350:** the reference height for the forcing.
- **NEE and respiration as targets,** once MEDS's soil carbon and growth respiration are developed. GPP's observation model then closes on MEDS's respiration.
- **Optional targets:** reflected PAR (visible optics), canopy temperature.
- **Multi-tower fits** by scope tag (Kuppel et al. 2012); emulator posteriors (Fer et al. 2018).

## 11. How PEcAn and ORCHIDEE handle the same problems

Read from the full text of Fer et al. (2018), Kuppel et al. (2012), Lasslop et al. (2008) and Pastorello et al. (2020); from the abstracts for MacBean et al. (2022) and Ingwersen et al. (2015).

| | PEcAn (Fer et al. 2018) | ORCHIDEE (Kuppel et al. 2012) | MEDS (this protocol) |
|---|---|---|---|
| targets | NEE, LE, half-hourly | NEE, LE, daily means | albedo, LW_up, H, LE, GPP, u\*; native interval |
| filters | u\* ≥ 0.40; no gap-filled hours | daily means of gap-filled data | measured hours; u\* per target |
| random error | heteroscedastic Laplace | — | σ_abs + σ_rel·\|x\| from provider/paired days; Huber |
| error size | variance sampled with the parameters | the prior model's misfit × 30 | measurement σ, scaled once by misfit |
| autocorrelation | N/N_eff from a time-series model of the data | the ×30 | N/N_eff from residuals |
| known bias | multiplicative k on soil-chamber respiration (lognormal prior) | KsoilC on initial soil carbon | κ on the tower's respiration; closure shares |
| structure | inside the variance; across models | multi-site, multi-stream | variants; prior z; labels |
| priors | trait meta-analysis | expert ranges | EEO, then meta-analysis |

- **Fer et al. (2018) on their bias term:** "While the introduction of the bias term makes it impossible for these data to constrain the magnitude of soil carbon fluxes, it does provide information on the shape of the functional response". κ makes the same trade.
- **Lasslop et al. (2008):** systematic errors bias the parameters, and "the real value of the parameter is not within the uncertainty range"; they "need to be addressed more thoroughly … since otherwise uncertainties will be vastly underestimated".
- **Ingwersen et al. (2015):** the closure methods' spread (up to 110 W m⁻² for LE) is an uncertainty band; "a single post-closing method might result in severe misinterpretations".
- **ONEFlux (Pastorello et al. 2020):**
  - NEE = FC + SC, with SC from a profile or one point;
  - u\* thresholds applied "to daytime and nighttime data", together with dropping the first half-hour after a calm period;
  - H and LE are turbulent only, with the closure factor computed away from sunrise and sunset.

## Sources

- Barr et al. (2013), Agricultural and Forest Meteorology 171–172, 31–45 (u\* threshold uncertainty)
- Björkman & Demmig (1987), Planta 170, 489–504 (photon yield of O₂ evolution; Fv/Fm)
- Brynjarsdóttir & O'Hagan (2014), Inverse Problems 30, 114007 (model discrepancy)
- Detto, BCI eddy covariance flux data 2012–2017, Dryad doi:10.5061/dryad.3tx95x6j5 (the provider's note)
- Fer et al. (2018), Biogeosciences 15, 5801–5830 (PEcAn emulator calibration)
- Hollinger & Richardson (2005), Tree Physiology 25, 873–885 (paired-day random error)
- Ingwersen et al. (2015), Biogeosciences 12, 2311 (closure uncertainty band)
- Kattge & Knorr (2007), Plant, Cell & Environment 30, 1176–1190 (thermal acclimation)
- Kattge et al. (2009), Global Change Biology 15, 976–991 (Vcmax by plant type)
- Keenan et al. (2011), Oecologia 167, 587–597 (model–data fusion assuming certain data)
- Kuppel et al. (2012), Biogeosciences 9, 3757–3776 (ORCHIDEE multi-site)
- Lasslop et al. (2008), Biogeosciences 5, 1311 (observation errors and parameter estimates)
- Lasslop et al. (2010), Global Change Biology 16, 187–208 (daytime partitioning)
- Lin et al. (2015), Nature Climate Change 5, 459–464 (g1 by plant type)
- Long, Postl & Bolhár-Nordenkampf (1993), Planta 189, 226–234 (maximum CO₂ quantum yield)
- MacBean et al. (2022), Global Biogeochemical Cycles 36 (15 years of ORCHIDEE assimilation)
- Malhi, Doughty & Galbraith (2011), Phil. Trans. R. Soc. B 366, 3225–3245 (cited by the provider)
- Ögren & Evans (1993), Planta 189, 182–190 (curvature of leaf light responses)
- Papale et al. (2006), Biogeosciences 3, 571–583 (u\* filtering)
- Pastorello et al. (2020), Scientific Data 7, 225 (FLUXNET2015, ONEFlux)
- Prentice et al. (2014), Ecology Letters 17, 82–91 (least-cost stomatal behaviour)
- Richardson et al. (2006), Agricultural and Forest Meteorology 136, 1–18 (flux random error)
- Rubio & Detto (2017), Ecology and Evolution 7 (BCI soil respiration)
- Skillman (2008), Journal of Experimental Botany 59, 1647–1661 (quantum yields)
- Smith et al. (2019), Ecology Letters 22, 506–517 (photosynthetic capacity optimized to environment)
- Stocker et al. (2020), Geoscientific Model Development 13, 1545–1581 (P-model v1.0)
- Wang et al. (2017), Nature Plants 3, 734–741 (coordination hypothesis)
- The earlier plans, their PRs (#329, #330, #340, #344, #345, #348) and the BCI diagnosis of 2026-10-01/02
