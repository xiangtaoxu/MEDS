# Changelog

All notable changes to MEDS. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
MEDS does not yet promise semantic versioning, because v0.x is explicitly pre-benchmark.

**This file is the only place change history lives.** Source comments state present-tense
rationale; what changed and when goes here; what is deferred goes in
[`docs/ROADMAP.md`](docs/ROADMAP.md). See `CLAUDE.md` for the rule.

Entries cite the pull request that shipped them. Anything that moves a number states the
before and after.

---

## [Unreleased]

### Added
- **The targets' observation models** (`scripts/calibrate_fast/obsmodels.py`, best-practice plan P3):
  - **The closure model for H and LE.** The daily closure factor over ±15 days, with shares from the
    attribution test and Bowen as the alternative.
  - **The respiration model for GPP.** κ is an observation key: never written to a config, and
    never fitted beside a shape key.
  - **σ from the provider's uncertainty or the paired days**, evaluated at a smoothed observation.
  - **Huber loss by default**, and σ scaled by the model's misfit at the refresh (max(1, √χ²/n),
    at most 3×).
  - **`[fluxes]` may declare G, GPP_DT and RECO_DT.** κ's prior sd is then the gap between the
    provider's two partitionings.

  At BCI:
  - daily closure median 0.75 (f = 1.33); the gap is H's;
  - κ prior 0.65 ± 0.10;
  - σ from the paired days: LE 10.5 + 0.29|LE|, H 7.5 + 0.14|H|, NEE 2.1 + 0.17|NEE|. These
    replace H's 30 % σ and GPP's 2.5 + 0.15 GPP.

  Net radiation, the evaporative fraction and night NEE are no longer targets, and the
  `growth-resp` command and `[base].growth_resp` go (#359).
- **The calibration's data rules** (`scripts/calibrate_fast/datarules.py`, best-practice plan P2):
  - **u\* per target.** A diagnostic on each turbulent target's own data, with a bootstrap over days.
    LE and H take its outcome (`ustar_min = "diagnostic"`). GPP takes the provider's threshold
    (`"provider"`, from the site TOML, one value or one per year).
  - **Windows by rule.** One per 1.5-month slot in the leaf-on months, each the slot's best covered,
    with validation windows in other years. Every window gets its own 180-day chain from the initial
    stand.
  - **Seasonal runs from a water-deficit index.** Rain minus Priestley–Taylor evaporation from the
    forcing; the runs are scored on LE only.
  - **Keys fixed from process coverage**, and the albedo's sun and frost rules.
  - **The report** gives the share of the record's conditions the fit spans, and the share of area
    whose canopy-air top is above the sensor.

  At BCI:
  - LE is flat in u\* and H still rising, so neither is filtered (they were at 0.4 and 0.6, with
    H also cut to 9–16 h).
  - The calibration windows are now picked by rule (they were picked by hand).
  - The seasonal runs are again 2016 and 2017, chosen by rule.
  - Daytime LE rows go from 17,737 (u\* ≥ 0.4, day and night) to 13,781, and H rows from 4,739 to
    13,781 (#358).
- **The site TOML's `[fluxes]` and `[provider]`** (`scripts/prepare_flux_tower`): each flux the
  calibration fits (SW_out, LW_out, Rnet, LE, H, NEE, GPP, RECO, USTAR) with its column, units and the
  rule that says when it was measured (`measured = { column, equals }`; without one, a FLUXNET flux's
  `_QC` = 0, any other flux wherever present), RECO optionally as a sum of declared fluxes, and what
  the provider did. `tower_inputs.read_standard()` is the one table both tools read: UTC interval
  starts, MEDS units, a measured mask for every column (#357).
- **Flux metadata checks F1–F3** (`tower_checks.check_fluxes`): net radiation against its four
  components (median residual within 25 W m-2), the upwelling longwave above the downwelling at night
  (a swapped pair of columns), and the reflected shortwave and upwelling longwave within physical
  bounds. They stop the forcing build and are reported by the calibration (`tower_checks` in the data
  report). BCI passes with its swapped longwave labels declared: F1 median residual 0, night LW_out −
  LW_in +32 W m-2 (#357).
- **`[leaf_physiology].medlyn_vpd_min`** (0.05 kPa, optional): the least leaf-to-air VPD the
  Medlyn stomatal model uses, in the equation's own units (g1 is in kPa^0.5). It was a bare 50 Pa in
  the code; CLM5 uses the same 50 Pa. At BCI (one year) 1.8 % of the lit leaf-hours sit below it,
  with 0.6 % of the gross assimilation; 0.02 or 0.1 kPa moves the year's GPP by at most 0.01 % and LE
  by at most 0.05 % (#347).
- **ctest `clamps_named`** (`scripts/lint/check_clamps.py`): a `max`/`min` against a real literal
  other than 0 or ±1 fails the suite unless the line says why (`! clamp-ok: <reason>`). A guard uses
  `tiny_num` or `safe_exp`; a physical threshold is a named setting with its source (CLAUDE.md, "No
  bare thresholds") (#347).
- **`meds.config`**, the Python API's view of a run's configuration: it reads the main TOML and the
  PFT (plant trait) TOML that `[init].pft_config` names, sets keys (one PFT's element of a trait
  array included), writes the pair for a run, and reads back the run's parameter record. It needs no
  compiled library (#340).
- **`meds.model.run(config)` and `python -m meds.model CONFIG`** run a config to its end through the
  Python API, in place of `meds_main CONFIG`: the same output files and the same closing line
  (#340).
- **ctest `python_api`** (with `MEDS_BUILD_PYLIB`): the package's tests, among them a run through
  `meds.model` compared with the same run of `meds_main`, bit for bit, and two runs in one process
  (#340).
- **`meds.canopy`**, the canopy's fast pieces on their own through the Python API: the two-stream
  radiation over a stand of cohorts, the light-plastic Vcmax25 and Rd25 per cohort, and a batch of
  leaf solves with every driver given. A new C API module (`meds_c_api_canopy.f90`) opens a config
  once and serves all three. `meds.plant.leaf.gas_exchange_batch` solves many leaves in one call
  (#345).
- **The leaf solve's per-cohort drivers as hourly outputs** (FAST tier, off unless listed):
  `gx_par_cohort_fast`, `gx_leaf_temp_cohort_fast`, `gx_vpd_cohort_fast`, `gx_ca_cohort_fast`,
  `gx_pressure_cohort_fast`, `gx_psi_leaf_cohort_fast`, `gx_psi_predawn_cohort_fast`,
  `gx_gb_cohort_fast`, `gx_agross_cohort_fast` and `abs_par_cohort_fast`, with the cohort's wood
  area `wai_cohort`. With them a canopy of leaf solves reproduces a run's hourly GPP: at BCI to
  0.43 % (median hour, 745 hours) from the hourly means of the four 15-minute steps. With the new
  outputs off, a run is bit-identical to one before (#345).
- **calibrate_fast: a staged, user-controlled fit** (`MEDS_FAST_CALIBRATION_REVISION_PLAN.md`,
  #345). The fit runs optics, then photosynthesis, energy, water stress and a joint polish:
  - optics uses the two-stream alone against the albedo;
  - photosynthesis uses a canopy of leaf solves against GPP (gate G8 checks it against the model);
  - water stress uses frozen 120-day dry-season runs and a grid search.

  It also adds:
  - `site_reference.toml`: every site setting with its default and the reason, checked against
    each site file;
  - per-target data filters (`ustar_min`, `par_min`, `hours`, `closure_range`,
    `min_solar_elevation`);
  - a `report` command: each target's rows through each filter, and the morning u* plateau test;
  - a rough uncertainty: σ-scaled Laplace intervals in physical units, and the MAP's shift under an
    alternative filter from the final Jacobian (gate G10).

  At BCI it ran 4,804 trials (31.8 core-hours), against 10,224 (about 43) for the joint fit shipped
  with #340. Over the five years its set gives GPP 7.52 µmol m⁻² s⁻¹ (tower 7.46; the shipped set's
  7.49), and April GPP of 6.6, 5.0 and 7.0 in 2014, 2016 and 2017 (tower 6.9, 6.3 and 7.0; shipped
  set 5.1, 3.7 and 5.3). The example still ships the joint fit's set; the decisions the run raises
  are in the plan's §13.

### Changed
- **`calibrate_fast` reads the tower through its site TOML** (plan P1). `[tower].site` names it, and
  the facts it holds leave `calibration.toml`: `[tower].path`, `time_column`, `flag_column`,
  `flag_good`, `utc_offset_h` and `[tower.columns]` are gone (a site that still sets them is refused).
  BCI's declarations moved into `bci_site.toml`. The observations stay on the tower's own interval:
  the tool sets `[output].fast_interval_steps` from the tower's interval and `fast.dt_fast` (2 at
  BCI) and a trial whose output is at another spacing fails; the hourly pairing, which needed both
  half hours measured, is gone. Local hours and days come from the site's clock. At BCI the rows left
  after every filter go from hours to half hours: GPP 4990 → 10272, LE 8498 → 17737, H 2253 → 4739,
  albedo 5448 → 10830, upwelling longwave 17078 → 34180 (#357).
- **`[leaf_physiology].phi_psii` 0.85 → 0.74 in the reference config, and its meaning.** It is the
  electron yield of linear transport in low light (J's initial slope is 0.5·phi_psii per absorbed
  photon), not Fv/Fm. With 0.85 and Aj's 4ci + 8Γ\*, a leaf fixed CO2 at most at 0.106 per absorbed
  photon, the O2-evolution yield (Björkman & Demmig 1987); measured CO2 fixation is 0.093 (Long et
  al. 1993), which 0.74 gives. In normal air (30 °C, ci 300 ppm) the leaf's yield goes from 0.064 to
  0.056, against the measured 0.052 ± 0.003 (Skillman 2008). The examples keep their own values. The
  Python leaf binding's default follows the reference config (#351).
- **`theta_j` 0.85/0.90 → 0.70 for the reference C3 PFTs.** Every leaf of a cohort gets the cohort's
  mean light, and a curve averaged over the light the leaves really see bends more than one leaf's,
  so a mean-light cohort's effective curvature sits below a leaf's (0.6–0.95 measured, Ögren & Evans
  1993). 0.7 is CLM5's and ED2's value. At the bend (absorbed light equal to the capacity), J goes from
  0.76 Jmax (0.9) to 0.65 Jmax. To revisit with sunlit/shaded leaves (#343). The calibration registry
  holds it fixed at the PFT file's value; the Python leaf binding's default follows. The examples keep
  their own values (#351).
- **The canopy films' water capacity is a plant trait** (breaking). `[soil].dewmx`, CLM's
  interception capacity, held only the leaf and wood films. It is replaced by the optional PFT
  traits `leaf_surf_water_max` [kg m⁻² leaf] and `wood_surf_water_max` [kg m⁻² wood], both 0.1 by
  default (the old `dewmx`). A config that sets `soil.dewmx` stops with a message naming them. The
  calibration registry's `dewmx` is `leaf_surf_water_max`, fixed. Output is bit-identical with
  interception off; with it on, the means agree to round-off (LE −0.0007 W m⁻²) (#347).
- **The stomatal water stress begins at an onset** (`pft.stomata_psi_onset`, new and optional).
  β_stomata is 1 while the predawn leaf potential stays above the onset and
  exp(sref · (ψ − onset)) below it; before, it fell from any negative potential. Without the key the
  onset is half the PFT's leaf turgor-loss point (−0.857 MPa at the reference leaf traits), the
  recommended value. Sabot et al. (2022, Eq. 5) apply no stress while the soil is at field
  capacity, and MEDS's predawn leaf potential carries a tree's gravity head even in wet soil: at
  BCI the top of the canopy sat at β = 0.51 in the wet season (predawn −0.34 MPa at 35 m), the
  middle at 0.67, the understory at 0.76; all three are now 1. On the BCI five-year run with the
  default parameters: LE 56.4 → 71.8 W m⁻² (tower 75.5; RMSE 50.1 → 38.9), H 77.3 → 66.8 (tower
  32.4), GPP 10.70 → 11.61 µmol m⁻² s⁻¹ (tower 7.46), April 2016 GPP 9.7 → 7.3 (tower 6.3); the
  budgets still close. The calibrated set shipped with the example was fitted without the onset
  (#341).
- **The fast calibration runs through the Python API** (`scripts/calibrate_fast`). Every trial,
  state chain and base record is built with `meds.config` and run by `python -m meds.model`;
  `--runner <meds_main>` runs the executable instead, and replaces `--meds-main`. The tool's own
  TOML module (`tomlio.py`) is gone, and the site declaration's `[base].pft` is refused: the PFT
  file is the one the main file names. `test/python/test_restart_exact.py` builds its configs with
  `meds.config` too, in place of its own TOML writer (#340).
- **The fast-calibration registry holds the 11 keys a ten-day window can set** (13 with canopy
  interception), down from 28 (plan §14). Left at their defaults: the eight keys the tower cannot
  inform, the three that act through soil water a ten-day window does not draw down
  (`wstress_sref_stomata`, `root_beta`, `leaf_pi0`), `d_ratio` (it repeats `z0m_ratio`), the
  numerical `ustmin` and `canopy_freeboard`, and the leaf biochemistry `jmax_vcmax_ratio`, `theta_j`
  and `ds_vcmax` (#340).
- **The BCI example's calibration is refitted** on v0.3.2 with the 11 keys, through the Python API.
  Against the 28-key registry refitted on v0.3.2: half the trials (9,824 against 18,504), 66
  core-hours against 151, and a validation objective 4 % higher (30,969 against 29,862), the loss in
  GPP and the evaporative fraction. Over the five years the new set does as well or better: GPP 7.56
  against the tower's 7.46 (the 28-key set's 6.77; v0.3.1's shipped set 6.83), NEE −3.63 (−3.31;
  −3.14), LE 80.5 (81.9), and April 2016 GPP 4.1 (3.9) against the tower's 6.3. The interception-on
  set now passes the five-year water budget (v0.3.1's had 53 breaches) (#340).

- **The fast-calibration registry is a menu the site chooses from** (#345). Each key is `fit`
  (the default set), `optional` or `fixed` (with the reason), has a stage and a prior, and a site
  edits the set (`[fit].keys`, `add`, `remove`) and any prior (`[priors.<key>]`). The default set
  changes:
  - in: `wstress_sref_stomata` and `stomata_psi_onset`, fitted on the dry-season runs;
  - out: `leaf_transmit_nir` (collinear with `leaf_reflect_nir`), `intercept_k` (rough) and
    `rd_vcmax_ratio` (fixed at 0.015 with night NEE out);
  - optional: `theta_j` (0.7–0.9) and `jmax_vcmax_ratio`.

  Default priors come from syntheses where there is one: `stomatal_g1` 3.77 kPa^0.5 (Lin et al.
  2015, tropical rainforest trees) in place of the base value with the range as its ±2 sd. The
  fit starts at the priors' centres.
- **calibrate_fast fits the optics and photosynthesis keys in the coupled stage, against the
  turbulent fluxes as measured** (#348, plan §13.4):
  - The default stages are `energy`, `water` and `polish`. A kernel stage (`optics`,
    `photosynthesis`) that is not listed has its keys fitted in `energy`.
  - No closure correction (`[tower].closure = "none"`, was `"bowen"`). At BCI the gap behaves like
    missing sensible heat.
  - LE is kept at u\* ≥ 0.4 m s⁻¹. H is kept at u\* ≥ 0.6 m s⁻¹ and 9–16 h only, with
    σ 10 W m⁻² + 30 % (was 15 %).
  - Net radiation and the evaporative fraction are off. The water stage scores LE and GPP.
  - Every target accepts every filter.
  - Registry:
    - fixed: `leaf_clumping` (0.80), `leaf_width` and `dsl_dmax`;
    - fitted: `theta_j`, `jmax_vcmax_ratio` (prior 1.70 ± 0.15), `ds_vcmax` (641 ± 5) and the
      new `ds_jmax` (640 ± 4);
    - new and optional: `ea_vcmax` and `ea_jmax`.

    The temperature priors are Kattge & Knorr (2007) acclimated at BCI's 25.5 °C. Slot & Winter
    (2017) put four Panama species' Vcmax optima at 32.9–39.7 °C; the default `ds_vcmax`, 650,
    peaks at 31.7 °C.
  - Two quick BCI fits of the coupled stage leave `vcmax25` at its floor (25.1, 25.3). `theta_j` and
    Jmax/Vcmax go to their floors too. At `vcmax25` 45 the model's GPP stays 1.44 times the tower's
    at every hour.
- **calibrate_fast's default targets and weights** (#345):
  - GPP's σ is 2.5 + 0.15 GPP (was 1.5 + 0.15 GPP), and its hours need u* ≥ 0.4 m s⁻¹ (was no
    filter): the tower's GPP carries its one-per-day respiration's error at every daytime hour, and
    calm mornings under-read uptake.
  - Night NEE is off (was on), with `rd_vcmax_ratio` fixed.
  - The screening reports by default and fixes nothing (was: fix the keys the tower cannot inform);
    `[fit].screening = "drop"` restores it.
  - The effective-sample-size weights apply in the fit itself (`[fit].weights = "ess"`), not only
    in the covariance.
  - A target the model cannot fit (χ² per row above 1) has its σ scaled up in the covariance.
  - A trial with non-finite residuals fails instead of entering the Jacobian.

### Fixed
- **The fast loop blew up when heavy rain hit dry soil** (#352). The ARK integrator takes each
  step's soil water from a separate soil-water solve, and that solve went wrong at rain onset.
  - **The cause:** after each linear solve, the solver recomputed the flows between layers with the
    conductivity at the NEW water potentials. A top layer dried near residual holds almost no water
    per metre of potential, so 30 mm/h of rain lifts its potential by hundreds of metres to
    saturation within one sub-step. The saturated conductivity times that gradient drained the layer
    *below residual* while it rained. Step-doubling could not see it, because the full step and the
    two half-steps overshot alike. From there the potential clamped at −3e15 m, the solve ran out of
    sub-steps, and the "flows" reached metres per second, clipped and floored back into ~1e7 kg m⁻² of
    water. The soil-heat step turned that into a top-soil temperature of −1e10 K, then NaN. RK45 never
    used these flows, so it survived.
  - **The fix:** the flows keep the conductivities and bottom-flux slope the solve used, at the
    solve's new potentials, so each layer changes by exactly the water the solve moved. Only the root
    uptake is re-read at the new potentials: it stays within the plant's request, which is all the
    plant side is credited (re-reading it linearly instead leaked water).
  - **BCI 2017 (trial 7fe525… of the cross-site calibration):** the run that stopped at the
    2017-04-17 storm now completes, with no energy or water budget breach and none under
    `[energy].debug_error`. In normal weather the change is small. On another trial over Jan–Apr
    2017, hourly LE moves by 0.33 W m⁻² RMS (mean 89), H by 0.09 W m⁻² and GPP by 0.0008
    µmol m⁻² s⁻¹. That is 25–300 times smaller than the usual gap between ARK and RK45.
  - **Test:** `test_column_hydrology` replays the failing solve's inputs from that patch (it failed
    with flows of 5.8 m s⁻¹ and 2.5e7 kg m⁻² created by the floor).
- **ARK's stages read the top soil too warm during rain** (found with #352). The stages carried the
  enthalpy of the water moving between layers, but held the water itself at its start-of-step amount.
  So the next stage read that enthalpy as heat in a dry layer: tens of kelvin late in a heavy-rain step
  (up to ~50 K at BCI), which the ground skin then saw. The stages now move soil water at the soil
  solve's steady rate, ending the step at the same committed amount as before. Over Jan–Aug 2017 at
  BCI, the ARK − RK45 gap in rain hours shrinks from 0.52 to 0.23 K for top-soil temperature and from
  6.1 to 5.5 W m⁻² for LE; dry hours are unchanged.
- **The wood heated itself whenever it refilled with water** (#355). Root uptake brings water into the
  wood and sapflow takes it out to the leaves, each carrying its liquid enthalpy, counted from the
  liquid datum (~4.4e5 J/kg). The wood's heat store is its heat capacity times its temperature, with no
  term for water mass coming or going. So whenever uptake exceeded sapflow, the arriving water's whole
  enthalpy was read as heat; the leaves had the same error on a smaller scale.
  - **At the BCI storm front of 2017-04-17,** wet soil refilled the dry-season wood with ~7 kg m⁻² in
    a quarter hour. The wood went 13 K above the canopy air (313 K in 295 K air) and gave the heat
    back as H ≈ +600 W m⁻² and LE ≈ +420 W m⁻² at zero net radiation, where the tower measured
    H ≈ −17 W m⁻². Every night it added ~10 W m⁻² of heat as the wood refilled.
  - **The fix** values each tissue's own water at the tissue's start-of-step temperature, the way the
    canopy film is valued at the liquid enthalpy its water arrived with. The tissue's temperature now
    sees only how far the arriving water's temperature is from its own, and both energy ledgers book
    the water stores' enthalpy change. The leaf still pays the full vapour enthalpy of what it
    transpires; its water store's outflow is counted at the step's transpiration demand.
  - **BCI 2017:** at the storm front H falls from 581 to 118 W m⁻² at 18 UTC and from 363 to 35 at
    19 UTC, and ARK and RK45 now agree. Over Jan–Jul, mean night-time H (01–09 UTC) goes from +2.5 to
    −0.9 W m⁻² (tower −23), mean H from 72.1 to 70.5 W m⁻²; GPP is unchanged.
  - **Test:** `test_column_ark` refills dried wood from moist soil at one temperature and checks the
    wood gains no heat, and that what the soil gives up equals what the tissues and their water stores
    gain.
- **Rain on bare ground was valued at the canopy-air temperature** (#355). It falls through the air
  above the canopy, so it now takes the air temperature at the canopy-air top (the forcing moved there),
  the same value rain landing on a snowpack and sub-threshold snowfall already used. No recorded reason
  for the old choice was found. At the BCI storm front the canopy air was 2.3 K warmer than that air;
  the change lowers H by 2 W m⁻² and LE by 3 W m⁻² at 18 UTC, and Jan–Jul means by under 0.1 W m⁻².
  With one rule for all precipitation, the snowpack and the bare ground now share one function for it
  (`precip_enthalpy` in `meds_therm_lib`) instead of writing the same expression twice (BCI output
  bit-identical; snow sites change at round-off).
- **The rain/snow split used a different temperature from the precipitation's enthalpy.** The split
  was made once at ingest, at the forcing's height (ERA5-Land's 2 m air, a clearing's), while the
  enthalpy uses each patch's canopy-air top. The forcing record now carries only the total
  precipitation, and it is split once, where each patch's forcing is filled (`fill_forcing`), at the
  patch's canopy-air top; the ingest split is gone. At freezing, a 30 m canopy top 0.12 K cooler than
  ERA5-Land's 2 m air turns 6 % of the rain to snow; BCI, never near freezing, is bit-identical.
  - **Outputs:** `rainf_fast`, `snowfall_fast` and `snowfall_site` report the split the patches
    received, area-weighted, instead of a split at the forcing's height that the model never used;
    the new `precip_fast` gives the total at the sub-daily tier, beside `precip_site`.
- **A soil-water solve that failed was used without a word.** Its `converged` flag went unread; on
  the #352 storm day 26 failed solves went into the state. A failed solve now counts as a failed
  check of the soil column: the end-of-run report prints a warning with the count, and
  `[energy].debug_error` stops the run.
- **The adaptive soil-water sub-stepping ignored whether its iteration converged.** Under
  `[soil] linearize = "picard"`, a sub-step whose iteration did not converge was accepted when its
  error estimate passed. It is now rejected and retried at a quarter of the size, like the fixed-count
  path already did. The default single linear solve always counts as converged, so default runs are
  unchanged by this.
- **A NaN in the fast loop went unreported, and the run still ended "no NaNs"** (#352).
  - **The old guard:** four cohort-structure fields (diameter, density, AGB, wood carbon), and only on
    year boundaries. Meanwhile the fast loop's integrators commit a non-finite state when even their
    smallest step fails.
  - **The new guard runs every step,** on:
    - the stand;
    - every cohort's leaf and wood temperature;
    - every patch's canopy air (enthalpy, humidity, CO2);
    - the run's energy and water ledgers, whose running residual turns non-finite the step a NaN
      enters them.

    It names the step it fired in.
  - **The #352 case** (BCI, a storm front on 2017-04-17) now stops in the step ending 2017-04-18
    00:00 instead of finishing with NaN fluxes and NaN budgets. `meds.model` raises the same.
- **The BCI provider's note on the tower's respiration was misquoted** in `site_reference.toml`
  (`[targets.gpp]`), the BCI `calibration.toml` and the calibration revision plan (§1, §6) as "40–50 %
  below soil chambers". The note says the respiration "appeared underestimated" against soil chambers,
  "considering that RECO includes also above ground respiration which can contribute up to 40-50% of
  total respiration". The tower's mean (4.09 µmol m⁻² s⁻¹) is about what the chambers measure for
  the soil alone (≈ 4.3; Rubio & Detto 2017). No setting changes.
- **A cohort with less than 0.1 m² m⁻² of leaf got too little light per leaf** (since v0.1.0). Its
  leaf PAR was its absorbed PAR divided by max(LAI, 0.1), so a cohort of LAI 0.01 saw a tenth of its
  light. It now divides by the cohort's own LAI. A stand of seedlings could not grow:
  - The biophysics example's spin-up from bare ground (Ithaca) grew nothing after #321 halved the
    leaf-area scale: under the floor, that halved a seedling's light per leaf as well as its leaf.
    At the end of its 50 years: LAI 0.000 → 3.76, AGB 0.001 → 9.42 kgC m⁻², cohorts 2 → 28 (before
    #321: LAI 4.1, AGB 9.6). After five years: LAI 0.000 → 0.048, cohorts 2 → 8.
  - BCI, five years from the 2010 census, default parameters: GPP 11.24 → 11.28 µmol m⁻² s⁻¹
    (+0.33 %), LE 68.56 → 68.68 W m⁻², H 66.79 → 66.73 W m⁻², LAI at the end 5.72 → 5.76.

  With it (#347):
  - the within-canopy wind treated a crown smaller than 1 % of the ground as 1 %; it now uses the
    crown's own area;
  - the forcing's wind had a 0.1 m s⁻¹ floor of its own under the 0.65 m s⁻¹
    (`aerodynamics.ubmin`) that every use applies; only `ubmin` is left.
- **A run through the Python API did not match `meds_main`** in an Intel build. Inside Python,
  `libmeds.so`'s calls to `exp`, `sin`, `pow` and the other math functions reached glibc's versions,
  which Python had loaded first, instead of Intel's. A ten-day BCI trial differed by up to 7e-4 W
  m⁻² in LE and 2e-5 µmol m⁻² s⁻¹ in GPP; the two-day test case differs in 38 output variables.
  `libmeds.so` now links Intel's runtime into itself and keeps its names inside (`-static-intel`,
  `--exclude-libs,ALL`), and both runs are bit-identical. The library grows from 2.5 to 4.0 MB
  (#340).
- **A second run in one process wrote the first run's settings into its parameter record**, and the
  record would have stopped growing after 4,096 rows. Loading a config now starts a new record
  (#340).
- **24 optional per-PFT keys could not be set since v0.3.2** (N-10): the per-PFT plant-hydraulics
  overrides (`pft.leaf_pi0`, `pft.wood_psi50`, `pft.k_plant_max` and ten more),
  `pft.storage_turnover_rate`, `pft.retained_carbon_fraction`, and the nine optional WATER, HYDRO
  and LIGHT phenology-cue keys. The loader reads them, but `meds_config_pft.toml` did not list them,
  so a config that set one stopped at the start. They are listed now, and the check that the
  references match the loader (`config_keys_listed`) reads the `opt_*` readers too; it had matched
  only `toml_*` and `req_*` calls. The calibration registry's `pft.leaf_pi0` is one of them, so the
  BCI calibration as shipped could not start on v0.3.2 (#340).

## [0.3.2] — 2026-10-01

An **efficiency and consolidation** release. The fast loop no longer slows down past four
threads, a region's polygons run on threads, and much of the code that listed the same thing in
two or three places now lists it once. The sub-daily output gains a patch axis and a skin
temperature, and a config key MEDS does not read now stops the run instead of being ignored.
- **Speed.** The BCI example's five years take 6 min 43 s on one thread and 1 min 25 s on 16; at
  4 threads, 2 min 23 s against v0.3.1's 5 min 26 s (#325). A region's 100 cells for a year take
  81 s on 40 threads, against 1,012 s serial in v0.3.1 (#183, #310).
- **Consolidation.** The fast-loop state's fields, the ARK storage, the tissue water curves, the
  forcing echo, the region and site serializers, the output set-up, and the fast tier's sample are
  each listed once (#146, #188, #195, #310, #311, #312). The met reader is split into its two
  sources (#311).
- **Output.** A patch axis on the FAST tier (#270) and a skin temperature (#275). Leaf and wood
  temperatures are within-step means, and a cohort recruited in the slow step, which has no
  samples yet, no longer reads 0 K (N-11).

**Upgrading from v0.3.1:**
- **Unknown config keys stop the run** (N-10). Every key MEDS reads is listed in
  `meds_config_main.toml` or `meds_config_pft.toml`. The error names every other key a config
  holds, with what replaced a retired one or the key a misspelling most likely meant: delete or
  rename each. Retired here are the `[output]` aliases `carbon_fluxes`, `water_fluxes` and
  `energy_fluxes` (use `carbon`, `water`, `energy`), `[fast].ark_niter` (use `ark_coupled`),
  `output.strict_caps`, and the PFT growth curve keys.
- **Removed outputs.** The four `*_var_site` variances are gone (#275); an `[output].io_config`
  that lists one is refused.
- **Outputs that move.** Leaf and wood temperatures change by up to several kelvin (the end-of-step
  sample is now a within-step mean). Cohort means read fill, not 0, on bare ground. Elsewhere, the
  thread and code-generation changes move results at rounding level.

### Changed

- **A config key MEDS does not read stops the run** (plan item N-10). Unknown keys used to be
  ignored in silence. `[forcing] dt_forcing = 7200` (the key is `timestep`) ran at the file's own
  spacing, and every shipped example, and the reference itself, carried keys nothing read.
  - **The rule.** Every key MEDS reads is listed in `meds_config_main.toml` or
    `meds_config_pft.toml`, set or commented out at its default. A key a config holds that its
    reference does not list is an error. The report names every such key at once, together with
    any missing required ones, before the run starts.
  - **The suggestions.** A retired key is named with what replaced it. A key that moved section is
    named with its new home (`state.cohort_max` → `output.cohort_max`). A misspelling is named with
    the key it most likely meant.
  - **The list cannot drift.** The build reads the list from the two references, and a new test
    (`config_keys_listed`) holds the references equal to the keys the loader reads, in both
    directions.
  - **The reference is complete.** `meds_config_main.toml` gains the 86 keys it lacked, at their
    defaults: a `[soil_carbon]` block, `[trait_dynamics]`, the `[fast]` solver settings and process
    mask, the thermal-acclimation keys, the soil optics, `[run].slow_on`, the `[init]` soil seeds,
    `[output].fast_interval_steps` and two longwave-synthesis keys.
  - **The reference was wrong in places.**
    - It documented `[output.fast] interval_steps`, which nothing read; the key is
      `[output].fast_interval_steps`.
    - Its reserved `fast.soil_water_coupling` did nothing.
    - The PFT reference listed four growth keys nothing read.
  - **Retired, with a message.**
    - The `[output]` aliases `carbon_fluxes`, `water_fluxes` and `energy_fluxes` (use `carbon`,
      `water`, `energy`).
    - `[fast].ark_niter` (use `ark_coupled`).
    - `[carbon].growth_source` and `[phenology].phenology_on`.
    - The PFT growth curve `growth_dbh_slope`, `growth_dbh_cap`, `growth_dbh_max` and
      `growth_lai_slope`.
  - **One table.** The refusals that were scattered through the loader are now entries in the
    retired-key table (`meds_config_keys`), with their messages.
  - **Migrating.** Delete or rename each key the error names. The shipped examples are cleaned.

- **The met reader is split into its sources** (#311 F6). `meds_met_driver` was 1,415 lines, with
  three backends behind ten branches. It is now four modules:
  - `meds_met_file_source`: a MEDS forcing file. It holds the open and one record's values as
    stored.
  - `meds_met_archive_source`: the ED_ERA5land archive. It holds the open, the month prefetch and
    one record's values.
  - `meds_met_source_common`: what the two share. That is the status codes, the time axis with its
    recycle window, and the check of a file against its config.
  - `meds_met_driver`: what is the same for every backend. That is opening, the cursor, the stepping,
    and the one ingest (`read_record`) that checks, converts, lapses and partitions a record.

  The routines moved unchanged, and every r1 output is bit-identical. A forcing that `met_open`
  rejects now prints its reason and then stops with one fixed message, because ifx garbles a stop
  code that is not a constant.

- **The FAST tier reads the patch block's own row** (#270). The fast sample was a list of its own
  (13 fields, 13 source ids, a 13-line fold) beside the patch block's table (`PD_*`) of the same
  quantities.
  - **One row.** `patch_diag_row` now fills that row once per patch and sub-step. The patch
    diagnostics accumulate it, the FAST tier stages it per patch, and the FAST site variables read
    its area-weighted sum (`SRC_F_PD0 + PD_*`).
  - **New rows.** The block gains NPP, Reco and the two-band upwelling shortwave, which the FAST tier
    needs. Its four water-flux rows (root uptake, infiltration, drainage, runoff) were never written
    or registered, and are gone.
  - **The residuals are rates in the row.** They become rates there (the step's residual over dt),
    so the whole row is dt-weighted alike.
  - **Staging follows the stand.** It grows with the live patch and cohort counts. The per-cohort
    slabs were `output.cohort_max` long, about 4.7 MB a polygon at the defaults.
  - **The caps are checked in the output layer.** The FAST tier checks both caps there, by the
    coarse tiers' rule: only when a cohort or patch variable is live. A run with the FAST tier on
    used to stop past `cohort_max` even when it wrote nothing per cohort.
  - **Rounding.** On the r1 cases every value is bit-identical except `et_rate_site`, at 4e-16.
    Before, its expression was one line that the compiler could reassociate.
  - **The listing.** `meds_io_config.toml` is regenerated. That also brings in the five forcing
    echoes that O6 moved to the `forcing` group.

- **The forcing echo is listed once** (#312 O6). The forcing each fast sub-step used was listed three
  times over: 15 fields of the fast sample, their copies in the fast loop, and 15 source ids with
  their cases in the output layer, beside the polygon block's own table (`PY_*`).
  - **One routine.** `forcing_echo` fills the polygon block's table once per sub-step; the coarse
    tiers accumulate it, and the FAST tier stages it as it is (`fast_forcing`). The FAST tier's
    forcing variables read their row of that table (`SRC_F_PY0 + PY_*`).
  - **One more row.** The table gains the liquid rain (`PY_RAINF`), which `rainf_fast` needs.
  - **Exact echoes.** `sw_in_fast`, `air_temp_fast` and `atm_co2_fast` were area-summed over the
    patches although they are the same everywhere. They are now the sample itself (a 1e-16 change).
  - **Group.** They, `sw_in_site` and `precip_site` join the `forcing` group. A config with
    `[output].forcing = false` that wants them lists them in its `io_config`, as both examples do.
  - Every other r1 value is bit-identical.

- **Forcing files are written in one contiguous block per variable, and read about 3 s faster**
  (MEDS_EFFICIENCY_SWEEP_PLAN.md, item N-1). The shared writer made the time dimension unlimited
  and took netCDF's default chunking: one record per chunk. The model reads a site's whole series at
  open, chunk by chunk, so the BCI file (90,528 records) took about 3.3 s of every run's start-up,
  most of a 10-day calibration trial's 8.5 s.
  - **The writer** gives the time dimension a fixed length, so each variable is stored contiguously.
    Rewritten this way, the BCI file shrinks from 63 MB to 4.2 MB, a 2-day run starts 3.1–3.4 s
    sooner, and every output variable is bit-identical.
  - **`met_open` prints a note** when a forcing file stores one time record per chunk, saying why
    it is slow and giving the `nccopy -c time/8760,grid/1` line that rewrites it with the same
    values. Files written before this change can be rewritten that way, or rebuilt.
  - **One script, not a folder.** `scripts/forcing_common/meds_forcing_file.py` becomes
    `scripts/meds_forcing_file.py`, the writer only, shared by `prepare_era5/make_forcing_file.py`
    and the flux-tower tool. The tower tool's Python copies of the model's conversions (humidity,
    pressure, solar geometry, the longwave synthesis) move to
    `scripts/prepare_flux_tower/tower_conversions.py`.
  - **Tests:** the flux-tower tool's 25 tests and the calibration smoke test pass unchanged.

- **The plant-hydraulics work counters include the corrector's solves** (plan item N-2). ARK
  re-solves the plant hydraulics on every step attempt, rejected ones included (the transpiration
  corrector in `advance_water_mass_full`), but that solve's sub-step and non-convergence counts were
  thrown away, so `work_hydro_nsub_site` and `work_nonconv_site` counted only the pre-pass solve. They
  now count both. On the established Ithaca stand in July `work_hydro_nsub_site` goes from 679 to
  2,460 sub-steps a day: about 2.4 step attempts per `dt_fast`, each re-solving every cohort.
  `work_hydro_thrash_site` still tests the pre-pass solve, and the physics is unchanged.

- **Leaf water potential converges at `dt_fast` = 900 s, and the documentation says so** (#162, plan
  item N-3). Four places still quoted the figure measured the day before the transpiration corrector
  (#91): "daytime mean −0.23 MPa at 12.5 s against −1.19 MPa at 900 s, not converged".
  - **Re-measured** on the established Ithaca stand over July at 12.5, 75 and 900 s: daily-mean leaf
    water potential at 900 s matches 12.5 s to within 0.001 MPa on every day after the first, each
    cohort's daily maximum to 0.008 MPa, and July GPP and ET to 0.1%. The first day after a restart
    carries a start-up transient (−1.83 against −0.28 MPa).
  - **Corrected:** `numerical_scheme.md` §5a, §6 and §7, the `dt_fast > 900 s` warning and the
    comments in `meds_config`, `examples/example_biophysics/meds_config_july.toml`, and ROADMAP §4.
    The capacity-limb table and the RK45 warning keep their figures, now labelled as measured before
    the corrector: that limb's table has not been re-measured, and RK45 has no corrector.
  - **ROADMAP §4–§5** mark #158, #159 and #167 closed, to revisit with the numerical scheme, and #104
    as planned (`MEDS_EFFICIENCY_SWEEP_PLAN.md` Phase 3).
  - **Removed:** `TISSUE_STORE_SCALE`, a source switch fixed at 1 whose banner still explained why the
    tissue heat store was "not on yet". The store is simply always on; outputs are bit-identical.

- **The fast loop's parameters are no longer copied into the frozen record every step** (#188).
  `column_frozen_t` carried `column_params_t`, a copy of the column's soil, thermal and hydraulics
  parameters taken from `column_config_t` once per `dt_fast` and patch: about 3.5 KB of fixed data
  plus a heap copy of the per-PFT hydraulics table (n_pft × 4.2 KB). The routines that read it
  (`column_be_stage`, `advance_water_mass_full`, `column_derivs`, both marches and both steps, the
  RK4 oracle) now take `col_config` beside the frozen record, the two clamps take it in its place,
  and `column_params_t` is gone. The frozen record now holds only what freezes each step.
  Bit-identical on every r1 case.

- **The fast-loop state's fields are listed once** (#146). Each of the seven state combinators
  (`state_init`, `state_axpy`, `state_accum`, `state_extrap`, `state_sub`, `state_err_diff`,
  `zero_like`) and the step-size error norm named every field of `column_state_t` itself, so adding a
  field meant finding and editing eight routines by hand. Now `state_to_array`, `array_to_state` and
  `tend_to_array` list the fields once, `state_entry_rules` says how each field enters the error norm
  and the embedded error, and the combinators work on the flat array. The combinators fill their
  output in place instead of allocating a new state.
  - **Outputs move at rounding level** (the sums are evaluated in a different order): at most 2e-9
    relative on the r1 cases' fluxes.
  - **Test:** `state_combinators` checks the round trip through the flat array and the rule counts.

- **The ARK march reuses its storage from step to step** (#195). Every ARK step allocated its three
  stage states, their surface tendencies and the march's trial states afresh, and each stage then
  copied its tendencies to the caller. Each thread now keeps one `ark_workspace_t` (in the patch
  loop's per-thread pool), and the stages, `column_be_stage` and `surface_derivs` fill it in place;
  its arrays are allocated again only when a thread moves to a patch with a different number of
  cohorts.
  - **Allocations:** 287 per patch-step before, 104 after (Ithaca, 14 cohorts, serial build; 304
    before #146). The allocator's share of a serial BCI run falls from 14.9% to 7.7% of CPU time.
  - **Timings**, 60 days of the BCI example, the second of two runs on an idle `R128C40` node: 1
    thread 25.9 → 23.7 s, 4 threads 14.6 → 13.6 s; 8 and 16 threads unchanged (11.5 s, mostly the
    serial start-up).
  - **Not reused:** the frozen record (about 50 of the remaining 104 allocations, 4% of the serial
    fast loop) is still built afresh each step. It relies on a new record's defaults for its
    scalars, and reusing it would need a reset that lists them all.
  - Outputs are bit-identical on every r1 case.

- **The canopy radiation lists the cohorts bottom-up without sorting them** (plan item N-4).
  `apply_rt_forcing` ordered the cohorts by height with a selection sort every sub-step (cost ∝
  cohorts²), although the cohort block is kept tallest first; the canopy aerodynamics already took
  the reverse in O(n) and kept the sort as a fallback. Both now call one routine,
  `ascending_order` (`meds_numerics`): the reverse when the heights do not increase, the sort
  otherwise, with the same order on ties. Bit-identical on every r1 case; no measurable change at
  BCI. **Test:** `numerics` checks both paths and the tie rule.

- **A tissue's pressure-volume traits travel as one record, and the storage curves have their own
  module** (plan Phase 3, step 2). About a dozen calls passed `water_content`, `capacitance` and
  `psi_from_water_content` the same four loose traits (π₀, ε, apoplastic fraction, saturated water)
  plus the biomass.
  - **The record.** `water_curve_t` holds the four traits, and `hydro_params_t` carries one per
    tissue (`leaf_curve`, `wood_curve`) instead of eight loose fields. The curve routines and
    `clamp_water_to_capacity` take the record and the biomass.
  - **The module split.** `meds_hydr_lib` held four jobs. The storage curves, tissue and soil
    (with the soil conductivity that follows from the retention curve), move to
    `shared/functions/meds_water_retention.f90`. The vulnerability curve, the Kirchhoff flux, its
    lookup table and the root profile stay in `meds_hydr_lib`.
  - **Removed:** `pv_water_cap_from_traits`, which nothing called.
  - Bit-identical on every r1 case. The science pages' code maps now name the new module, and a
    stale row that still pointed at `gauss_legendre_7` and `bisect_root` (deleted in #325) is fixed.

- **A region's polygons run side by side on threads** (#183 R3, #310 R4 and R5; plan Phase 5).
  `[run].n_threads` now means polygon threads in a region run, where it had to be 1. Each month's
  polygons run in one OpenMP loop, costliest (by last month's time) first, each with a
  single-threaded patch loop; the forcing load before the month and the output after it stay on one
  thread. A polygon's patch-loop thread count is set by `polygon_prepare` (the fast context's new
  `patch_threads`), not read from the config, so a site run keeps threading its patches.
  - **Timings**, the 100 cells of a 1° box around Ithaca, one year, each run alone on an idle
    `R128C40` node: 1 thread 785 s (v0.3.1: 1,012 s), 10 threads 131 s, 20 threads 92 s, 40 threads
    81 s. The one detail polygon costs about four times an ordinary one (29 s a year, for its hourly
    site files and per-cohort diagnostics), so on 40 threads it sets each month's pace: without it
    the run takes 53 s, and with no output at all 51 s. The region's own files cost about 2 s.
  - **The same results at any thread count:** all 403 files of that run are identical at 1, 10, 20
    and 40 threads, and `test_region` now runs its region on four threads against the site runs.
  - **A failed polygon no longer stops the month** (R4). It is reported and stops; the others
    finish the month, its output is written, and the region moves on. Before, the month was left
    unfinished, and a caller that went on stepped the earlier polygons through it again.
  - **A month's steps are listed once** (R5), instead of once for a forcing check and again for the
    polygons. The region's write-only step counters are gone (R12).

- **One output set-up, and one call into the stepper** (#310 R2, R6, R9, R11). Bit-identical on
  every r1 case.
  - **R2.** The six-call sequence that lays out an output file set was written three times (the
    site run, the region's files, a detail polygon's files). It is now `open_output_files`, and
    each polygon joins a set with `attach_output`. Both sit beside `polygon_prepare`, with
    `apply_io_overrides` and `ensure_output_dir`, which the region no longer borrows from the site
    driver.
  - **R11.** `activate_site_diag` adds to what a site already accumulates, so a detail polygon keeps
    the blocks both of its file sets need, not only the last one's.
  - **R6.** `advance_one_step` hands every optional argument on instead of branching on which are
    present (six calls down to three), and the polygon's step makes one call instead of two. The
    branch that dropped the latitude is gone. `advance_slow_dynamics` takes the step's start and
    reads its day of year itself.
  - **R9.** The polygon carries its location; the forcing cursor and the leaf phenology read it
    from there, not from the cursor.

- **One serializer for site and region files** (#312 O7, O8). The region writer was a copy of the
  site writer, about 130 lines: the dimensions, the calendar and axis coordinates, the variable
  definitions and the record writer. One `write_record` now writes every file set. A site's (or a
  detail polygon's) set has one polygon's buffers and no polygon axis; a region's has every
  polygon's, packed into one hyperslab per variable.
  - **The entry points.** The manager's four entry points become two, `output_serialize` and
    `output_close`, which take the file set whole instead of nine of its fields. A site run and a
    detail polygon hold their buffers as an array of one, as a region already did.
  - **What a reader sees is unchanged.** A site file now writes the fill value into the slab
    entries past a record's live length, where it used to leave them unwritten; netCDF returned
    the fill value for them either way.
  - Every r1 case and `test_region`'s region-against-site comparison are bit-identical.

- **Small consolidations and checks** (#312 O10, O12; #311 F9, F11, F12):
  - **O10.** `validate_config` refuses an `[output].fast_interval_steps` that does not divide the
    fast steps in a slow step, so no fast-tier window straddles two slow steps (the stand can be
    restructured between them). Every shipped config uses 4 or 24 of 96.
  - **O12.** `apply_patch_disturbance` no longer grows the patch diagnostic block a second time:
    `patch_ensure_capacity` already grows it with the patch arrays.
  - **F9.** Deleted `wind_log_profile`, which nothing has called since #305, with its three test
    checks, and `forcing_config_t%rad_sw_ground_const`, which nothing read.
  - **F11.** The ED_ERA5land reader refuses a month file whose `_FillValue` is a number: it
    recognises a missing value only as NaN, so such a number would have been read as data. The
    archive stores NaN. **Test:** `met_era5land` writes a month declaring `_FillValue = -9999`.
  - **F12.** `docs/science/forcing.md` says that a run cannot start in the archive's first month:
    its first step reads the 00:00 record, which lives in the previous month's file.

- **The forcing reader's names and literals, each in one place** (#311 F7, F8). Bit-identical on
  every r1 case.
  - **F7.** The archive's hourly spacing is one constant, `ARCHIVE_DT_SEC`, where it was a literal
    in five places and a sixth encoding (24 records a day). The two hand-written "seconds since"
    parsers are one, `time_units_base` in `meds_time`.
  - **F8.** The 15 fields a MEDS forcing file may carry are listed once (`MEDS_FIELD`, with named
    indices). `met_open` finds each field's column once, where the reader used to match names as
    strings for every field of every record, and one `series_value` replaces `read_scalar` and
    `read_scalar_default`. When the forcing bracket slides by one record, as it does once an hour,
    its old upper record becomes the new lower one instead of being read again.

- **`output.strict_caps` is gone, and a run stops at the step it outgrows a cap** (#312 O9).
  `strict_caps` promised a warn-and-truncate mode that was never built, so it parsed and did
  nothing; a config that sets it is now refused, naming what to do. Nothing checked
  `cohort_max` or `patch_max` until a record was written, up to a month late, while the fast tier's
  per-cohort slabs, `cohort_max` long, were written by site slot every sub-step. A run now stops at
  the step its live cohort or patch count exceeds the cap of an axis some live variable uses (or,
  for the fast tier, the cohort cap), with a message naming the cap. Removed from the shipped
  configs; bit-identical on every r1 case. **Test:** `test_region` checks the refusal.

### Added

- **A skin temperature** (#275). `skin_temp_site` is the skin temperature as land models define it
  (CLM's `TSKIN`). It also comes per patch (`skin_temp_patch`), sub-daily (`skin_temp_fast`) and both
  (`skin_temp_patch_fast`).
  - **What it is.** The black-body temperature of the longwave leaving the canopy top, what an
    infrared thermometer or a satellite land-surface temperature sees.
  - **Why it is a row.** It is a fourth root, so it is formed per patch and sub-step, in the patch
    block, like the VPD.
  - **What it is not.** MEDS has no separate ground skin, although #275 assumed one: snow-free, the
    ground surface is the top soil layer (`soil_temp_top_site`). The patch row that holds it, which
    was labelled "ground/skin temperature", is now `PD_SOIL_TEMP_TOP`.

- **The FAST tier has a patch axis** (#270). A site mean over a closed canopy and a gap can describe
  neither: in the biophysics example the gap's surface soil ran 13 K above the air at midday and the
  closed patch's 0.7 K below it.
  - **Fifteen variables.** Each FAST quantity of the patch block has a per-patch twin, named after
    its coarse patch variable plus `_fast`: `cas_temp_patch_fast`, `soil_temp_top_patch_fast`,
    `le_patch_fast`, `h_patch_fast`, `rnet_patch_fast` and `nee_patch_fast`.
  - A quantity with no coarse patch variable keeps its site stem: `gpp_rate_`, `sw_up_`, `lw_up_`,
    `ustar_`, `npp_rate_`, `reco_` and `cas_co2_patch_fast`.
  - The soil columns are `soil_temp_layer_patch_fast` and `soil_water_layer_patch_fast`.
  - They follow `[output].axes_patch` (on) and `axes_soil_patch` (off), as the coarse patch
    variables do; a region writes none of them. The site `*_fast` value is their area-weighted sum.

### Removed

- **The integrator-study scripts** `scripts/numerics_sweep.py`, `scripts/parity_scenarios.py` and
  `scripts/parity_fidelity.py`. They served the integrator selection and parity studies, which are
  finished: ARK is the production integrator, and the parity plan is retired.
  - `parity_scenarios.py` no longer ran. It wrote the `[io]` block refused since v0.3.0, and
    started from restarts spun up on the retired `split` scheme.
  - `parity_fidelity.py` only scored `numerics_sweep.py` output. Nothing else (no test, example
    or script) used either.
  - They stay in git history; `scripts/calibrate_fast/` covers driving trial configs.

- **The four variance outputs and their operator** (#275). `cas_temp_var_site`,
  `soil_temp_top_var_site`, `cas_vpd_var_site` and `leaf_temp_var_site` are gone, and so is
  `AGG_VARIANCE` (`cell_methods = "time: variance"`).
  - **Why.** Nothing used them. They squared one end-of-step sample per slow step, which is the
    day-to-day spread of the state at one hour, not the within-step spread their names suggested.
  - **Migrating.** An `[output].io_config` that still lists one of them is refused, as for any
    name the registry does not have. Delete the line.
  - **Gone with them.** The six end-of-step accessors only they read: canopy-air temperature,
    humidity, CO2 and VPD, soil-top temperature and surface water. Their time means were already
    rows of the patch block (#264), so the registry test that kept time means off those accessors
    goes too.
  - **Adding one back** needs a sum of squares inside the step, and first a choice of definition
    (`docs/ROADMAP.md`).

### Fixed

- **Leaf and wood temperatures are within-step means** (plan item N-11). `leaf_temp_site`,
  `wood_temp_site`, `leaf_temp_cohort` and `wood_temp_cohort` averaged the end-of-step temperature:
  one sample per slow step, taken at the boundary's local hour. That is the bias #264 removed from
  the canopy air and the soil.
  - **The new source.** They now read the cohort block's rows (`CD_LEAF_TEMP`, `CD_WOOD_TEMP`),
    sampled every sub-step. Water output already keeps that block on, so default runs do no extra
    work.
  - **Slow-only runs.** As for the other fast-loop variables, a run without the fast loop reads
    these as fill.
  - **A cohort with no samples is left out.** That is a cohort recruited in the slow step, after the
    fast loop. It is fill on the cohort axis and out of every mean and sum. It used to read 0, so a
    recruit pulled every leaf-area-weighted cohort mean towards 0 (0 K for a temperature).
  - **A mean over nothing is fill.** A site mean with no weight (bare ground, say) reported 0, and
    `leaf_temp_site` read 0 K on bare ground. It is now fill, as the empty-set rule says; sums are
    unchanged.

- **A soil-by-patch variable is handled like the other patch variables** (plan item N-9). Three rules
  named the cohort and patch axes and missed the soil-by-patch profiles:
  - **The record's patch count.** A tier whose only patch output was a `*_layer_patch` variable
    wrote one patch, with every value as fill.
  - **The file cap.** Such a tier's files were not capped at a month.
  - **The annual guard.** An `[output].io_config` could put one on the annual stream.
  - **One rule now.** All three use one test for these axes (`ragged_dim`).

- **More than four threads no longer slow the fast loop** (#325). Two fast-loop routines handed one of
  their contained functions to another routine: `flux_potential` passed its Kirchhoff integrand to the
  quadrature, and `solve_leaf_gas_exchange` passed its Ci residuals to the bisection. Under ifx each
  such call allocates a lock-guarded record, and at 8 threads those allocations took 83% of all CPU
  time (232 of 280 CPU-seconds on 60 days of the BCI example).
  - **What changed.** The quadrature is a written-out 7-point sum (`kirchhoff_integral`). The leaf
    solver's residuals are module functions that take the leaf's problem as one explicit record
    (`ci_problem_t`), with a bisection written for it; one field (`gs_rule`) replaces the two flags
    that chose between the model, cuticular and pinned conductance passes. `bisect_root`,
    `gauss_legendre_7` and `phi_inverse` are deleted: nothing else used them.
  - **Timings**, the BCI example, each run alone on an idle `R128C40` node, the second of two runs:

    | `n_threads` | 60 days, v0.3.1 | 60 days, now | five years, now |
    |---|---|---|---|
    | 1 | 29.3 s | 25.5 s | 6 min 43 s |
    | 4 | 21.9 s | 15.1 s | 2 min 23 s (v0.3.1: 5 min 26 s) |
    | 8 | 40.9 s | 12.3 s | 1 min 36 s |
    | 16 | 54.4 s | 12.1 s | 1 min 25 s |

    Sixty days is mostly the serial start-up (about 8 s), and a site run gains little past its
    patch count (BCI keeps 14–25 patches).
  - **Outputs move at rounding level**, from the changed code generation: on the r1 cases the ARK
    fluxes differ by at most 2e-9 relative over a month (hourly H by 1e-8 W m⁻²); the stage-clamp
    counter `work_clamp_stage_site` flips by up to 14% because it counts threshold crossings. RK45
    amplifies the same first-day difference (9e-11 W m⁻²) to a few W m⁻² of hourly H within a week.
  - **Guard:** a new ctest, `no_procedure_arguments`, fails if any `procedure(...)` declaration
    appears under `src/`. `plant_hydraulics` now checks the quadrature against a composite-Simpson
    reference and the closed form, instead of through `phi_inverse`.
- **A discarded RK45 step no longer stops a debug run on its budget check** (#189, item 3). Under
  `[energy].debug_error`, `column_fast_step_rk45` stopped the run on a whole-column ledger breach
  before `column_fast_step` decided whether to keep the step or redo it on ARK, so a step about to be
  thrown away could end the run. The step now records its checks without stopping, and the
  dispatcher stops (`rk45_ledgers_stop`) only for a step it keeps. The failure counts were already
  right, because a rescued step's budget is rolled back. Outputs are unchanged. **Test:** `numerics`
  checks the new `last_check_closed` on a breached and a closed check.
- **The daily tissue-water reconcile uses each PFT's own curve** (plan item N-7, found during Phase
  3). `reconcile_tissue_water_capacity`, which seeds an empty leaf or wood store at the starting
  potential and caps a store above saturation, read the shared `[hydraulics]` traits, while the
  fast loop reads the same water back on the PFT's own curve (`[pft]` overrides, #179). A PFT with
  its own saturated water or pressure-volume traits was therefore seeded off its curve and capped at
  the wrong ceiling. It now takes the per-PFT hydraulics table. No shipped configuration overrides
  these traits per PFT, so every r1 case is bit-identical. **Test:** `pft_optics_config` seeds and
  caps a PFT whose saturated leaf water differs from the shared value.
- **A slow-only run no longer reports stale values for fast-loop diagnostics** (found while
  merging the serializers). With `fast_biophysics_on = false` the per-cohort and per-patch
  diagnostic blocks hold no entries, so their readers left the caller's array unset, and the
  site-level reduction then read whatever the previous variable had put there. On r1's
  `demography_30yr`, `soil_temp_top_site` read the canopy-air depth, a "5 K" soil temperature, and
  `w_surface_site` read 0 only by chance. A slot the block holds no entry for now reads 0, as one
  with no weight already did. #299 remains: those variables should read as missing in a slow-only
  run, and the slow operators' own rows should report their values.
- **A slow-only run reports its slow rates, and reads fast-loop variables as missing** (#299). With
  `fast_biophysics_on = false`:
  - **Fast-loop variables.** 39 of them (the fluxes, canopy air, ground, forcing echo) read 0 on
    r1's `demography_30yr`, a value nothing computed. They now read as `_FillValue`. The
    per-cohort fast block and the polygon block stay off in such a run; the patch block reports
    only the slow operators' rows (`patch_diag_slow_row`). A variable the run does not simulate is
    caught once, before extraction, and the accumulator skips a missing sample, so a window closes
    as exactly the fill value. Averaging the fill over patches used to give a number just off it.
  - **Slow rows.** The slow step already weighted the patch block by `dt_slow`, but only the fast
    loop set the block's patch count, so the slow rows were unreadable until a disturbance set it.
    The step now sets it too: `nplant_recruit_site` reports its first year (0.03 plant/m²/yr
    instead of 0). The disturbed area (0.0139/yr), background and disturbance mortality carbon
    were already reported.
  - **Litterfall** stays 0 in a run with soil carbon off, which does not accumulate litter.
  - **Test:** `slow_diag_units` checks the slow-only weight, patch count and recruitment rate.

## [0.3.1] — 2026-09-29

A **flux-tower** release. MEDS now runs from a tower's own meteorology, starts from a forest census,
and calibrates its fast parameters against the tower's fluxes, with Barro Colorado Island as the
worked example throughout.
- **Forcing from towers.** `scripts/prepare_flux_tower/` declares, validates and builds a forcing
  file from tower data, and the reader is UTC-only (#320). The longwave gap fill is the model's own
  synthesis regressed onto the tower, and BCI's swapped longwave columns are read the right way
  round (#326).
- **A census start.** The BCI example starts from the 2010 census of the 50-ha plot, restructured
  before the first step (#323), with a ceiling on the patch-fusion tolerance (#322).
- **Calibration.** `scripts/calibrate_fast/` fits the sub-daily parameters to a tower with the stand
  held fixed, and the BCI example ships a calibration scored on windows the fit never saw (#330). It
  needed exact restarts, trait re-acclimation at restart, a parameter record and hourly upwelling
  radiation (#329).

**Physics fixes, several of which move results:**
- The reported sensible heat was about 100 W m⁻² too low (#328).
- The fast loop's emissivity, wood temperature, hydraulics keys, turgor-loss point and ground
  optics are now per PFT or live (#327).
- Allometry works in carbon throughout (#321).
- Soil-water faces take ED2's geometric rule, so a dried surface re-wets (#320).
- Stomata close linearly from the turgor-loss point to twice it, instead of shutting in one step
  (#335).

**Water budgets.** The calibration's parameter sets exposed three water-budget faults, all fixed:
- the per-face soil check subtracted the requested root uptake instead of the uptake the solver
  removed (#334);
- the ARK water ledger did not declare the water the tissue-water floor creates;
- interception discarded film water above capacity (both #336).

**Portability and output:** nvfortran works again (#317), slab output is sized after the
`io_config` overrides (#318), and OpenMP is compiled in by default (#324).

**Upgrading.** Forcing configs change in four ways (#320, `MEDS_FLUX_TOWER_FORCING_PLAN.md` §5):
- rename `[forcing].format = "netcdf"` to `"ED_default"` and `"era5land"` to `"ED_ERA5land"`;
- delete `[site].utc_offset` and `[site].apply_solar_longitude`, and build every forcing file on a
  UTC clock with `time_zone = "UTC"` (`make_forcing_file.py` always has);
- with `[site].apply_elevation_lapse = false`, delete `lapse_rate_tair` and `grid_elevation`;
- a file that states `wind_meas_height_m` (every `make_forcing_file.py` file does, as 10 m) must agree
  with `[forcing].wind_height`.

Each old form stops at startup with a message naming the fix.

The root profile moves to `[hydraulics]` (#327). `[soil_column].root_beta`, an exponential decay per
metre, is refused; `[hydraulics].root_beta = exp(−b · root_depth)` gives the same decay b, and a config
that sets neither keeps the default profile. A `[hydraulics].root_beta` copied from the old
`meds_config_main.toml` (0.96) used to be ignored and now takes effect: it puts 19% of the roots in
the top 0.37 m instead of 53%. Delete it to keep the default.

### Changed

- **Stomata close gradually at low leaf water potential, not in one step** (#335, #332).
  `[leaf_physiology].low_water_potential_control = "linear_decline"` is optional, and the only option.
  - **What it does.** The conductance the stomatal model calculates, g0 included, is multiplied by a
    factor that falls linearly from 1 at the leaf's turgor-loss point, ψ_tlp, to 0 at twice it. The
    factor is set from the previous day's predawn leaf potential.
  - **The solve stays coupled.** Leuning and Medlyn apply the factor inside the Ci solve. Katul
    re-solves with gs pinned at the factor times its optimum.
  - **A fully closed leaf** exchanges no CO₂ or water by day (net assimilation 0, its respiration
    refixed), and respires at night.
  - **The carbon consequence.** Rd is unchanged and still charged in full. As gs goes to 0 the
    coupled solve drives net assimilation to 0, so gross assimilation, which the canopy counts as
    GPP, tends to Rd. A fully closed leaf is therefore carbon-neutral by day (GPP = Rd, respiration
    Rd), where the former shutdown gave it GPP 0 and a loss of Rd. At night it loses Rd, as before.
    A tower's GPP, partitioned from NEE, cannot see refixed CO₂.
  - **It replaces a hard shutdown at 2·ψ_tlp** (`ARREST_GS_CLAMP`), which no config could change. That
    shut a cohort completely below the threshold and left it untouched above. The step made the
    fluxes jump as a parameter moved the threshold or the predawn potential across it, and no
    gradient-based calibration could see past it (the BCI calibration, #330).
  - **What changes.** A run whose cohorts stay above ψ_tlp is bit-identical. Between ψ_tlp and 2·ψ_tlp
    the conductance now falls. Below 2·ψ_tlp a leaf's daytime net assimilation is 0, where it was −Rd.
    A thermodynamic limit on transpiration (#96) would make this control matter less.
  - **On the BCI example, five years:**
    - The default run changes only in the fifth digit (GPP 10.7023 → 10.7022 µmol m⁻² s⁻¹).
    - With the calibrated set (#330), closing from ψ_tlp keeps the plants from drying out. The
      lowest predawn potential goes from −21.6 to −6.8 MPa, and cohort-days below 2·ψ_tlp from
      2.1 % to 0.12 %.
    - April GPP in 2014, 2016 and 2017 rises from 4.3, 2.9 and 4.8 to 5.2, 3.9 and 5.3 (the tower:
      6.9, 6.3 and 7.0).
    - Both runs close their budgets.
  - **Tests:** `leaf_physiology` checks:
    - the factor itself: 1 at ψ_tlp, linear, 0 at 2·ψ_tlp;
    - half way through the band, gs is half the Medlyn conductance of the solved leaf;
    - for Leuning, Medlyn and Katul, gs is continuous, never rises as ψ falls, and stays consistent
      with A through diffusion;
    - past 2·ψ_tlp, a leaf exchanges nothing by day and respires at night.

    The old shutdown fails it. `test_region` refuses any other value of the key.

- **The BCI example's leaf traits follow the canopy's light gradient**
  (#330; `[trait_dynamics].trait_plasticity_on = true`; `MEDS_FAST_CALIBRATION_PLAN.md` D6). Each cohort's
  Vcmax25, Rd25, SLA and leaf lifespan are its PFT's top-of-canopy values scaled by the leaf area
  above it.
  - **The stand.** Over the five tower years LAI now holds at 5.6 where it fell to 4.8, and AGB
    rises from 16.1 to 18.1 kgC m⁻² where it reached 17.6.
  - **The fluxes.** GPP is 10.70 µmol m⁻² s⁻¹ where it was 11.11, and NEE is −4.17 where it was
    −4.07. LE, H and net radiation move by 0.1 W m⁻² or less.
  - The README and `evaluation.png` are regenerated.

- **The flux-tower tool fills the longwave by the synthesis regression only; its ERA5-Land fill is
  removed** (#326). Filling from ERA5-Land or another source is the user's to do in the tower file before
  the build. `make_tower_forcing.py` loses `--lw-fill`, `--states-fill` and `--era5-file`, and a
  site TOML's `gapfill.longwave`, `gapfill.states` or `gapfill.era5_file` stops the build with that
  message. A rain gap now stops the build instead of taking ERA5-Land's rain, the V5 report no
  longer gives the barometer height ERA5-Land implied, and qc code 2 is unused.
  `compare_longwave_fill.py` scores the synthesis regression beside a monthly climatology and the
  synthesis as MEDS computes it. Scored against the ED_ERA5land archive before its removal, the
  ERA5-Land fill tied with the synthesis at BCI (RMSE 14.8 against 14.5 W m⁻² on 7,200 hidden half
  hours). The BCI example builds one forcing file, `data/bci_forcing.nc`. Tests: two ERA5-Land
  cases removed, and one added for a refused `gapfill` key.

- **OpenMP is compiled in by default: `MEDS_OPENMP` is now `ON`** (#324). The thread count stays a
  run-time setting, `[run].n_threads`, default 1, so a default build runs serially until a config asks for
  threads. `-DMEDS_OPENMP=OFF` builds serial, and a compiler with no Fortran OpenMP now falls back to
  serial with a CMake warning where it used to stop the configure. The Python wheel stays serial
  (`python/pyproject.toml`). An existing build directory keeps its cached value. In a default build
  `test_fast_loop`'s 4-thread check now runs threaded; in a serial build the directives are comments
  and it passes trivially.
  - **Numbers move at rounding level.** An OpenMP build gives the same bytes at every thread count,
    but not the bytes of a serial build: compiling the patch loop as a parallel region changes its
    rounding. On the BCI census example the two builds part in the 13th significant digit after 62
    hours, and single hourly fluxes differ by up to 2% of their largest value after five years. The
    tower statistics and the stand agree to every printed digit, and the test suite passes
    unchanged.
  - Measured on the five-year BCI census example, each run alone on a 40-core node: serial 7 min
    11 s; OpenMP at 1 thread 7 min 17 s, 4 threads 5 min 44 s, 8 threads 13 min 15 s, 16 threads
    20 min 10 s. More than four threads are slower because ifx's bound-procedure-value allocations in
    `flux_potential` and `solve_leaf_gas_exchange` serialize the threads (`docs/building.md`,
    "Parallel builds").
  - Tests: ifx Release, ifx Debug, gfortran Release and a serial ifx build, 60/60 each; the ifx
    OpenMP suite also under an 8 MB stack.

- **The BCI flux-tower example starts from the 2010 census of the BCI 50-ha plot, with no spin-up**
  (#323; `MEDS_BCI_CENSUS_INIT_PLAN.md`). `meds_config_spinup.toml` and its 50-year run are gone.
  `bci_census.toml` declares the census; `run_example.py` builds the census file with
  `scripts/prepare_census` and runs the five tower years from it, with the soil at 298.65 K and soil
  carbon in steady state with the stand's litter. MEDS fuses the 1,250 quadrat patches and 84,937
  rows to 25 patches and 419 cohorts before the first step. `plot_evaluation.py` now shows the mean
  diurnal and seasonal cycles of GPP, NEE, latent and sensible heat and net radiation against the
  tower, and writes their statistics. The README and `evaluation.png` are regenerated from this run.

- **A census stand is restructured before the first step** (#323; `MEDS_BCI_CENSUS_INIT_PLAN.md` §5.4).
  After `init_from_census`, the driver applies the slow step's own monthly cohort block (fusion, cull,
  fission, sort) and yearly patch block (fusion, cull, cohort fusion), without recruitment or
  disturbance and under the same switches. A census used to run its first month with every row a
  cohort and its first year with every cell a patch. The run log prints the patch and cohort counts
  before and after. Test: `init_census` fuses 20 identical cells to one patch with one cohort per
  size, conserving the site's stems and biomass, and does nothing with the switches off.

- **The census reader matches columns by name and reads `patch_area`**
  (#323; `MEDS_BCI_CENSUS_INIT_PLAN.md` §5.3). `init_from_census` takes its columns from the header line
  in any order: `patch_id`, `dbh`, `pft` and `nplant` required, `patch_area`, `site_id`,
  `cohort_id` and `height` optional. With `patch_area` the patches take their areas normalized to
  the site, where every census patch used to get an equal share. A header-less file of seven
  numbers per row still reads positionally, so existing census files load unchanged. An unknown or
  repeated column, a missing required one, and a patch whose rows disagree on `patch_area` stop
  the run with the name. Tests: `init_census` (reordered columns and areas),
  `init_census_refuses_area`, `init_census_refuses_column`.

- **The initial soil state is configurable** (#323; `MEDS_BCI_CENSUS_INIT_PLAN.md` §5.1).
  `[init].soil_temp` [K] and `[init].soil_theta` [m³ m⁻³] set every soil layer of every patch at the
  start of a run that does not restore the soil from a state file. They default to the constants the
  fast context carried, 288 K and 0.30, so no existing config changes. The loader refuses a
  temperature outside 233–333 K and a water content outside (`theta_res`, `theta_sat`]. Test:
  `init_soil_state`.

- **Patch fusion's tolerance has a ceiling, `[demography].patch_light_tol_max`, default 0.15** (#322).
  The light-profile tolerance steps geometrically from `patch_light_tol` to the ceiling over
  `n_patch_fusion_iter` passes, as the cohort tolerance does, and goes no further. Patches more
  different than the ceiling stay apart even when the count is still above `max_patch`, which is
  now a target. Before, the tolerance grew by a fixed 1.5 per pass with no bound, from 0.10 to 0.76
  over six passes, and the last pass's largest-difference limit, 1.14, could never bind. So
  `max_patch` acted as a hard limit, paid for with heterogeneity.
  - The key is optional. Absent, it is 0.15, or `patch_light_tol` if that is larger. A ceiling of
    0.759375 reproduces the old schedule exactly: a 50-year biophysics spin-up is identical to the
    build before this change.
  - **No shipped example moves.** Their patch counts stay within `max_patch`, so fusion never gets
    past the first pass, which uses `patch_light_tol` under both schedules. The biophysics spin-up
    (50 years, with the allometry before #321) and the demography golden are identical under both.
  - **Where it binds:** a stand more heterogeneous than `max_patch` patches can hold. Simulated on
    the Barro Colorado Island 2010 census with one patch per 20 m quadrat, six passes leave 24
    patches, where the old schedule went down to 9.
  - A run whose restructuring left more than `max_patch` patches, or more than `max_cohort`
    cohorts in a patch, says so at the end, for example
    `NOTE: patch fusion left 24 patches on ... (max_patch = 12)`.
  - `test_patch` gains the case: two unlike patches under `max_patch = 1` stay two under the
    ceiling and fuse under the old schedule. It fails with the ceiling removed.
  - Documented in `docs/configuration.md`, "Cohort and patch fusion".

- **The demography example's golden is recaptured for #321** (#322). #321 moved
  `test/golden/empirical_spinup_golden.csv` without recapturing it, so `empirical_spinup.py`
  reported a maximum relative error of 4.9e-1 in `total_agb` and 1.7e-1 in `total_nplant`; it
  reports 0 again. At year 40 the stand's AGB goes from 10.45 to 8.48 kgC m⁻² and its LAI from
  6.62 to 5.97. The example README records the table.

- **The default biomass law is Chave et al. (2014), and both it and the leaf-area scale are in
  carbon** (#321). In `meds_config_pft.toml`, every example's PFT file and the `meds_allometry`
  initializers:
  - `[allometry].agb_c1` goes from 0.06080334 to 0.03365 and `agb_c2` from 1.0044785 to 0.976:
    Chave's eq. 4, AGB = 0.0673 (ρD²H)^0.976 kg dry mass, divided by `C2B = 2`. The MEDS form
    `agb_c1·ρ^agb_c2·(D²H)^agb_c2` is Chave's exactly. The old values were ED2's `IALLOM = 3` refit
    of Chave (`c14f15_bs_tf`), in dry mass but read as carbon.
  - `lai_b1` goes from 0.46769540 to 0.23384770, with `lai_b2` unchanged. It is ED2's BAAD leaf fit
    `c14f15_bl_xx` as ED2 applies it (`size2bl` divides by `C2B`); MEDS used it undivided, so leaf area
    was twice ED2's.
  - **Per tree**, at a given diameter and the default heights: a 100 cm, 42 m tree goes from 16,261
    to 6,315 kgC of AGB and from 1,887 to 944 m² of leaf; a 10 cm tree from 43 to 20 kgC. Wood,
    leaf, fine-root and storage carbon follow. The recruit unit carbon falls by about half too, so a
    given reproduction flux makes about twice the recruits.
  - **Per stand**, on the Barro Colorado Island 2010 census (207,259 trees) the stand's AGB goes
    from 39.8 to 16.1 kgC m⁻² against the census's own 15.1, and its LAI from 11.2 to 5.6.
  - **Upgrading:** a PFT file copied from an earlier `meds_config_pft.toml` keeps the old values;
    copy the three new ones into it. A file with its own fitted values needs nothing.
  - **Every example's figures and quoted numbers predate this change** and are regenerated with the
    next release.
  - The new test `allometry_defaults` checks that the shipped `[allometry]` block equals the
    `meds_allometry` initializers, that the biomass law is Chave's in carbon, and that the leaf
    scale is ED2's over `C2B`. It fails on the old values.

- **Forcing files carry the humidity their source measured, and MEDS converts it** (#320,
  `MEDS_FLUX_TOWER_FORCING_PLAN.md` D2). An `ED_default` file carries exactly one of `RHair`
  (a fraction; a flux tower), `Tdew` (a reanalysis) or `Qair`, and the reader turns it into specific
  humidity at each stamp with the model's own Bolton curve, as it always did for the ED_ERA5land
  archive's dewpoint. A file with none, with two, or with `RHair` above 1.5 (a percentage) is
  refused. The point is the saturation curve: a tower's `vpd` column, or a `q` made from it
  offline, carries the provider's curve, and Barro Colorado Island's is Alduchov–Eskridge, whose
  saturation pressure is 5.7 Pa below Bolton's at 25 °C. With `RHair` in the file, the model's
  relative humidity at the forcing temperature is the tower's to 4e-15, and a saturated record reads
  back as VPD = 0. Qair files still load unchanged; `make_forcing_file.py` now writes `Tdew`.
- **Every forcing clock is UTC** (#320, D1). `[site].utc_offset` and `apply_solar_longitude` are refused,
  and an `ED_default` file whose `time_zone` attribute is missing or not `"UTC"` stops at open,
  because a local-time file read as UTC keeps its daily totals and moves its sun. Solar time is the
  UTC clock plus the longitude and the equation of time.
- **The two file formats are named `"ED_default"` and `"ED_ERA5land"`** (#320, D3). The old names stop
  with the name that replaced them.
- **A file's stated heights are checked against `[forcing]`** (#320). `tq_height_m`, `wind_height_m`,
  `wind_meas_height_m` (within 0.01 m) and `height_above`, when present, must match
  `tq_height`, `wind_height` and `height_above`, because a disagreement moves every sample to the
  canopy-air top from the wrong height.
- **The terrain-lapse keys are read only with the lapse on** (#320). `[site].lapse_rate_tair` and
  `grid_elevation` are required with `apply_elevation_lapse = true` and refused with it off; they
  used to be required either way and did nothing.
- **`specific_humidity_to_vpd` is the exact inverse of the forcing conversions** (#320). It used the
  molar-mass ratio 0.621987 where every forward conversion and `sat_specific_humidity` use 0.622,
  so a humidity round trip was off by 2e-5 in relative humidity (0.06 Pa of vapour pressure at
  3 kPa). It is used only by output diagnostics (`cas_vpd_site` and `cas_vpd_var_site`), which move by
  that much; nothing in the model state changes.

- **Soil-water faces take ED2's geometric rule** (#320, `MEDS_FLUX_TOWER_FORCING_PLAN.md` §13).
  - **Between layers:** conductivity is ln K interpolated linearly between the two nodes to the face,
    the thickness-weighted geometric mean of `rk4_derivs`. It replaces the upstream pick.
  - **At the surface:** the infiltration capacity is the geometric mean of K_sat and the top layer's
    K, times the gradient to the top node, where it was the top layer's own K. That rule sealed a
    dried surface: near residual water content K is 3×10⁻¹⁵ of K_sat for the default loam, so after
    a dry season the pond overflowed and the soil never re-wet.
  - **The aquifer bottom boundary** stays upstream-weighted.
  - **Barro Colorado Island** is where it showed. Before, the top layer sat at θ = 0.081 through
    the 1963 wet season, 3-year ET was 16 % of rain, and the 50-year spin-up from bare ground ended
    at LAI 0.009. Now the top layer re-wets to θ = 0.25–0.33, ET is 26 % of rain, and the stand
    reaches LAI 4.8 and AGB 15.3 kgC m⁻².
  - **Ithaca** is barely moved: its 50-year spin-up ends at AGB 9.956 against 9.936 kgC m⁻², with
    the same stem density and LAI.
  - **Capillary rise into a dry profile is slower.** In `test_column_hydrology`'s aquifer case, the
    residual bottom flux passes 10⁻⁵ kg m⁻² s⁻¹ at about 1,000 h instead of about 450 h, so that
    test's relaxation window is now 1,200 h.
  - **New tests:** the face flux against the known log-linear answer, and a dried top layer
    re-wetting under 12 h of 2 mm/h rain. Both are mutation-checked: with the old surface rule, 79 %
    of the rain runs off.

### Added

- **`scripts/calibrate_fast`: calibration of the fast parameters against a flux tower**
  (#330; `MEDS_FAST_CALIBRATION_PLAN.md` P1, P3). It fits the sub-daily parameters (radiation,
  photosynthesis and stomata, aerodynamics, water stress, respiration) to a tower's albedo,
  upwelling longwave, net radiation, LE and H corrected for closure with the Bowen ratio kept,
  daily evaporative fraction, daytime GPP, night NEE and u\*, with the stand held fixed. Every
  trial is a 10-day `slow_on = false` restart with `reacclimate_traits`, from the state that a chain
  of frozen runs wrote at its window's start.
  - **The method.** Levenberg–Marquardt with Gaussian priors on logit-transformed parameters, a
    central-difference Jacobian whose trials all run at once, three damping values tried at once,
    and three starts. A screening step keeps the keys the tower can inform. The covariance is a
    Laplace approximation weighted by each target's effective sample size, and a linearity check
    tests it.
  - **Every trial proves what it ran.** Its parameter record must list every key the trial set,
    and a whole-site budget breach fails the trial.
  - **Running it.** Trials run on a local pool or, on a cluster, on a directory queue served by one
    worker per node inside one allocation. The commands are `select-windows`, `growth-resp`, `check`
    (gates G1 and G2), `fit`, `analyze` and `write-calibrated`.
  - **Tests:** 17 unit tests, and a smoke test through `meds_main` on the demography census with a
    synthetic tower (ctest `calibrate_fast`).
- **The BCI example has a calibration** (#330; `MEDS_FAST_CALIBRATION_PLAN.md` P2). `calibration.toml`
  sets it up, and `calibration/` ships the fit and the calibrated configs. The fit used 8 ten-day
  windows from 2015–17 and was scored on 8 it never saw, with interception off and on.
  - **What the shipped set does.** It is interception off, with 20 fitted keys. It lowers the
    validation objective from 68,018 to 29,790 and every target's error, by 7 % (upwelling
    longwave) to 70 % (albedo).
  - **Over the five tower years:** GPP 6.75 against the tower's 7.46 µmol m⁻² s⁻¹ (default 10.70),
    LE 82.4 against 75.5 W m⁻² (56.4), net radiation 135.3 against 136.3 (120.6), albedo 0.17
    against 0.13 (0.26), and night u\* 0.50 against 0.41 m s⁻¹ (0.86). H is still 36 W m⁻² high.
  - **A structural limit.** Eight keys end at a bound of their range, among them `vcmax25`,
    `stomatal_g1` and `z0m_ratio`, so part of the misfit is not in the parameters.
  - **Gates.** G1–G7 pass. The five-year run closes its energy and water budgets and its slow ledger.
    The interception-on fit scores the same but fails G7's water budget (#333).
  - **Known limits.** The late dry season is too stressed: April GPP is 3.9 against the tower's 6.3
    in 2016 (2.9 before the stomata closed gradually, #335). The fit held the two hydraulic keys at
    their defaults because the former whole-day shutdown made them rough.
  - **Running it.** `run_example.py` runs the calibrated five years beside the default, and
    `evaluation.png` and `calibration.png` draw both. `run_example.py --calibrate` redoes the fit,
    about 100 core-hours.

- **A restart can take this run's leaf traits: `[init].reacclimate_traits`** (#329; default false,
  restart only). The plastic traits (`sla`, `vcmax25`, `rd25`, leaf lifespan) are then set from this
  run's PFT file as a census start sets them: acclimated to each cohort's LAI above it, as the state
  holds it, with plasticity on, and the PFT's top-of-canopy values with it off. Leaf area stays as
  read, and leaf carbon scales by the SLA's change, storage taking the difference, so unchanged traits
  change nothing. Without it a restart
  keeps the state's traits, so a changed `vcmax25` never reached the cohorts. This is what lets a
  calibration trial restart from a shared state (`MEDS_FAST_CALIBRATION_PLAN.md` P0b). Tests:
  `restart_exact` checks that a restart with `vcmax25` × 1.3 carries the traits a census start with
  × 1.3 gives, with and without plasticity, and keeps the leaf area.
- **The parameter record, `<prefix>_parameters.csv`** (#329), beside the diagnostic output and
  beside the state: one row per key the loader read from any file, `source,key,index,present,value`,
  with `present` saying whether it was set in the file or defaulted and `value` the value used, to
  17 digits. A key nothing reads is absent, which is how a misspelt key in an optional block, until
  now silently ignored, can be caught. The BCI example's record has 375 rows. Tests:
  `test_biophysics_opts_config` checks a set key, a defaulted key and a misspelt one.
- **Hourly `sw_up_fast` and `lw_up_fast`** (#329): the shortwave (VIS + NIR) and the longwave,
  emission included, leaving the canopy top, beside `rnet_fast`. Until now they were daily only.

- **`[soil]` sets the bare ground's optics** (#327): `ground_albedo_vis` (0.15), `ground_albedo_nir`
  (0.30) and `ground_emissivity` (0.95), the values the code had fixed. The canopy radiation solver
  reads them, and snow still covers them by its fraction. The albedos must lie in [0, 1) and the
  emissivity in (0, 1]. Tests: `test_biophysics_opts_config` reads them and `test_pft_optics_config`
  checks that they reach the fast loop.

- **`scripts/prepare_census/make_census.py`: a ForestGEO tree table to a MEDS census**
  (#323; `MEDS_BCI_CENSUS_INIT_PLAN.md` §6). It only maps trees to patches: one patch per square plot cell
  at its true area, one row per distinct (cell, diameter), `nplant` the count over the area. It
  keeps live trees with a diameter of at least the declared minimum, and counts every exclusion. Given
  a PFT file and an earlier census, its summary also carries the stand's steady-state litter input for
  `[soil_carbon].spinup_steady`. On the Barro Colorado Island 2010 census at 20 m it writes 84,937
  rows in 1,250 patches, conserving all 207,259 stems. Tests: `prepare_census` (pytest, run by CTest
  when the Python it finds has numpy, pandas and pytest). `environment.yml` gains pandas and pytest,
  which this tool and the flux-tower tool need.

- **`scripts/prepare_flux_tower/`: MEDS forcing from flux-tower data** (#320; AmeriFlux BASE,
  FLUXNET/ONEFlux or any CSV). The tool works from a site TOML that declares the file, location,
  clock, stamp convention, sensor heights and every column's units.
  - **Checks (V1–V5).** It validates each declaration and stops on a disagreement: a uniform axis
    (V1); the clock against the model's own sun, within 10 min (V2); a provider VPD against RH
    under the declared saturation curve, naming the curve that fits (V3); physical bounds (V4). V5
    is a JSON report.
  - **What it writes.** An `ED_default` file on a UTC clock with the measured `RHair`, pressure
    brought down to the ground, states re-centred to the stamps, and the tower's heights stated.
  - **Gap filling** is explicit, with a `<Var>_qc` flag on every value: short gaps interpolated,
    long ones by the mean diurnal variation. Longwave is filled from the model's synthesis,
    regressed onto the tower in its clear-sky and cloud parts. (An ERA5-Land fill added here was
    removed before release, #326; see Changed.)
  - **`compare_longwave_fill.py`** scores the fills on held-out observations. `tests/` has 24
    pytest cases on synthetic towers, run by CTest as `prepare_flux_tower` when the Python it finds
    has the dependencies; two mutations of the tool (the UTC sign, the re-centring) fail 15 and 2
    of them.
- **`scripts/forcing_common/meds_forcing_file.py`** (#320), the one writer of an `ED_default` file and the
  Python copy of the model's conversions (humidity, hypsometric pressure, solar geometry, window-mean
  cos z, clearness index, longwave synthesis). `make_forcing_file.py` now writes through it and
  stores ERA5-Land's dewpoint as `Tdew`, with `tq_height_m = 2` and `height_above = "zero_plane"`.
- **`examples/example_flux_tower_bci/`** (#320): the worked example at Barro Colorado Island (AmeriFlux
  PA-Bar, a 41 m tower).
  - **Data.** It downloads the CC0 data from Zenodo and checks the md5; nothing is committed.
  - **The build.** It passes V2 5 min from the declared UTC−5 begin-stamped clock, and V3 to
    0.000 Pa under the Alduchov–Eskridge curve.
  - **Longwave.** The observed longwave is 61 % missing. The synthesis regression fills it with
    RMSE 9.1 W m⁻² on hidden records, where the model's `lwdown_source = "synthesize"` would be
    36.3 (bias −22.5). Those numbers were measured on the upwelling column, which the file labels
    as downwelling; the correction is under Fixed (#326).
  - **The model stages.** A 50-year spin-up and a five-year evaluation compare MEDS with the tower in
    local time. They needed the soil-water fix above: before it, rain did not infiltrate a dried
    top layer and no stand grew.
- **`test_met_tower`** (#320), the flux-tower contract of an `ED_default` file: the three humidity forms
  and their rejections, the UTC requirement, stated heights, rain and shortwave from the interval
  containing the instant on end- and begin-stamped files, and the tower round trip (relative
  humidity, VPD = 0 at saturation, the move from a 41 m tower to a canopy-air top). Each fix above
  was mutation-checked: restoring the old rain read fails the three end-stamped rain checks, and
  restoring the old ratio fails the humidity round trip by 1.9e-5.
- **`test_region` refusals** (#320) for the old format names, `utc_offset`, `apply_solar_longitude` and the
  lapse keys with the lapse off, each run through `meds_main`.

### Fixed

- **Canopy interception discarded film water above capacity** (#336, #333). A cohort can start a step
  holding more film water than its capacity, `dewmx·(LAI + WAI)`, when it loses leaf area under a full
  film (the daily slow step sheds leaves). `intercept_canopy_layer` clipped the film to capacity, and
  the clipped water reached neither the ground nor any flux.
  - **The fix.** The excess now drips to the layer below with the rest of the drip.
  - **Where it showed.** With interception on, in five BCI years under the calibrated set, 25 `whole_water
    (ark)` breaches in heavy rain lost 1e-6 to 4e-5 kg m⁻² each. At every one, rain minus the film's
    interception minus the throughfall equalled the residual. With the fix: 0 breaches.
  - **What does not change.** Runs with interception off (the default) are unaffected, and a film
    within capacity behaves as before.
  - **Test:** `column_hydrology` starts a film 0.05 kg m⁻² above capacity. The storage is capped, and
    throughfall plus the storage change equals the rain. The old kernel fails it.
- **The ARK water ledger did not declare the water the tissue-water floor creates** (#336, #333). When one
  step's transpiration debit would take a cohort's leaf or wood water below zero,
  `advance_water_mass_full` floors the store and creates water (#148). The mass was reported
  (`work_clamp_mass_site`) but not entered in `whole_water (ark)`, so every firing breached the ledger
  by exactly the water made.
  - **The fix.** The ledger now takes it as an input, as it already takes the soil's θ_res floor. The
    adaptive march used to add the floor's mass on every attempted sub-step, rejected ones included,
    and now counts accepted sub-steps only, so the declared mass is what the committed state received.
    `work_clamp_mass_site` therefore now reports committed water only.
  - **Where it showed.** In five BCI years under the calibrated set with `dsl_dmax` = 0.015, all 91
    breaches fell on the three days the floor fired (April 2016). The daily residual equalled the
    floor's mass. With the fix: 0 breaches, worst residual 5e-13 kg m⁻².
  - **What is left.** Only the ledger changes: the floor still creates water where the plant
    hydraulics collapse (#104).
  - **Test:** `test_column_derivs` checks that on a floored step the plant water changes by its fluxes
    plus the reported floor mass, exactly.

- **The soil column's per-face check reported a wilting-limited root sink as a face error** (#334).
  `advance_soil_water_column` checks each layer's change of water against its face fluxes and its
  root sink. It subtracted the plant's requested uptake, `forcing%root_uptake(k)·dt`. In a layer
  drier than `psi_open` the solver removes less than that, because the wilting ramp cuts the sink,
  and the check reported the difference as a face inconsistency. With `[energy].debug_error` that
  stopped the run.
  - **The fix.** `soil_water_advance` now returns the sink it removed from each layer, summed over
    the accepted sub-steps, and the check subtracts that.
  - **Only the diagnostic changes.** A 31-day run of the BCI census example matches the previous
    build on every output variable.
  - **Under the BCI calibration's parameter sets.** Over five years, `faces[soil_layer_mass]` had
    read 1.1–2.6 kg m⁻², against 0.0011 with the defaults. For the shipped calibrated set it now
    reads 4.6e-13. A `debug_error` run of the set with the most water-budget breaches now gets past
    October 2012 and stops at the first real breach, in March 2016.
  - **Test:** `column_hydrology` puts a top layer at θ = 0.10, where the ramp passes 74 % of the
    demand, and the face residual reads 1.9e-14. The previous check fails it with 2.4e-2.

- **The slow ledger missed the tissue heat that trait plasticity moves** (#330). With plasticity on,
  `advance_plant_traits` moves leaf carbon above the new SLA's allometric target into storage, and
  the leaves' heat content changes with it. The allocate phase then left an undeclared energy
  residual: −5,903 J m⁻² over the BCI example's five years. `vegetation_dynamics` now declares the
  change, and the ledger closes. The fluxes are unchanged.

- **A restart did not continue the run that wrote the state** (#329). The state file kept each
  cohort's dbh but not its carbon pools or geometry, and the reader rebuilt them on the allometry.
  Cohort fusion keeps the pools and leaves a fused cohort below the allometry for its dbh, so the
  restart moved the leaf area (BCI: LAI 5.6380 written, 5.6389 read) and, with plasticity on, reset
  the fine-root carbon, which sets root respiration and rhizosphere conductance, by up to 2× in the
  most shaded cohorts. The file also lacked the LAI above each cohort (so the traits' light
  environment restarted as an open sky), the canopy's interception film, the growth-rate buffer
  behind the mortality predictor, and the slow loop's CO₂ hand-off to the fast NEE. All are stored
  now, optional on read, so an older state file still restarts as before.
  - Tests: `restart_exact` checks that a run split at a restart matches the unsplit run on every
    state variable, bit for bit, with the slow loop off (interception on) and on. The previous
    binary fails it in 13 and 28 variables.
  - The BCI example, which starts from a census, is unchanged: a 31-day run matches `beta` on every
    output variable.

- **The reported sensible heat was about 100 W m⁻² too low, and every surface-layer solve too
  stable** (#328). The reference air's potential temperature was referenced to the ground,
  `T + (g/cp)·zref`, and the canopy air's was not, while the canopy air's energy budget exchanged heat
  on the two actual temperatures. So the reported H, `g_ah·cp·(T_cas − θ_atm)`, sat
  `g_ah·g·zref` below the flux the budget booked, and the Monin–Obukhov solve saw `(g/cp)·zref`
  (0.37 K at a 38 m canopy-air top) of stable stratification that was not there. Both potential
  temperatures are now referenced to `zref`, the canopy-air top the forcing is moved to, so each is
  its actual temperature there (`set_aero_env_atm`; `canopy_aerodynamics.md` §2). The budgets closed
  before and still do. Every H MEDS reported before this fix is low by `g_ah·g·zref`, with zref the
  fixed reference height (30 m by default) before 774d040 moved the forcing to each patch's
  canopy-air top, and that top since.
  - On the five-year BCI census example:

    | | before | after | tower |
    |---|---|---|---|
    | sensible heat, mean [W m⁻²] | −39.1 | 77.3 | 32.4 |
    | sensible heat, night (19–05 h) | −95.5 | −0.1 | −23.6 |
    | sensible heat, midday (10–14 h) | 89.7 | 244.2 | 163.0 |
    | sensible heat, r hourly / r seasonal | 0.86 / −0.59 | 0.93 / 0.95 | |
    | Rnet − H − LE, all hours, five-year mean | +104.6 | −11.4 | |
    | friction velocity, night [m s⁻¹] | 0.78 | 0.87 | 0.41 |

    GPP, NEE, latent heat and net radiation move by less than 1%. The miss the example shows is now
    the partition, not the sign: a midday Bowen ratio of 1.5 against the tower's 0.69, and the
    README says so.
  - Tests: `test_column_dynamics` RUN 5b checks that canopy air at the reference air's temperature
    and humidity is a neutral surface layer with no temperature scale, and that a dark isothermal
    column reports |H| under 5% of `g_ah·g·zref` under both integrators (−0.4 against 26 W m⁻²).
    With the old line restored all four checks fail, at H = −16.9 W m⁻².

- **The tissue energy balance gave every leaf and all wood an emissivity of 0.95** (#327). The
  radiation solver takes each PFT's `leaf_emissivity` (0.97 by default) and `wood_emissivity` (0.90) and
  emits at the canopy-air temperature. The energy balance then adds the change in emission for the
  tissue's departure from that temperature, 4εσT³ per unit leaf or wood area, and took ε from a fixed
  `leaf_emiss = 0.95` for both. It now takes the cohort's PFT values, so the two halves of the longwave
  agree, and `veg_thermal_params_t` loses `leaf_emiss`. Tests: `test_column_ark` checks that the frozen
  tissue carries each cohort's own emissivities.
  - Together with the next three fixes, on the five-year BCI census example: sensible heat moves
    from −39.01 to −39.08 W m⁻², and its midday mean from 90.3 to 89.7. GPP, NEE, latent heat, net
    radiation, the albedo and the stand agree to the digits the README prints. At BCI only this fix,
    the wood temperature and the root profile apply: it has one PFT with no pressure–volume traits of
    its own, the default ground optics, and whole-plant conductance.
- **The wood's boundary layer took its free convection from the leaf temperature** (#327).
  `boundary_gbh_mos` adds free convection driven by an element's difference from the canopy-air
  temperature, and `aero_bottom_to_top` passed the leaf's temperature for the wood too. The wood now
  gets its own. Tests: `test_column_dynamics` checks that the wood conductance ignores the leaf
  temperature and rises with the wood's.
- **Four `[hydraulics]` settings did nothing** (#327). `root_beta` and `root_depth` were read and
  never used, because the root profile came from `[soil_column].root_beta`. `wood_kmax` and
  `vessel_curl` feed the segment conductance, which no key could select.
  - The root profile is now ED2's β^(d/D) from `[hydraulics]`, integrated over each layer and
    renormalized over the column; a layer below `root_depth` holds no roots (Upgrading above;
    `plant_hydraulics.md` §4). The default β = e⁻⁴ with D = 2 m is the old e^(−2d) decay, integrated
    over each layer instead of sampled at its centre. On the default ten-layer grid that moves 0.4% of
    the roots, and no layer by more than 0.15 percentage points. ED2's own default, β = 0.001, would
    put 72% of them in the top 0.37 m.
  - `[hydraulics].conductance = "whole_plant"` (default) or `"segment"` selects the conductance form.
  - `validate_config` refuses `root_beta` outside (0, 1), `root_depth` ≤ 0 and, in segment mode, a
    non-positive `wood_kmax` or `vessel_curl`, shared or per PFT.
  - Tests: `test_plant_hydraulics` checks the segment conductance; `test_soil_column_config` checks
    that the profile sums to one, that the default is the integrated e^(−2d), and that a shallow
    `root_depth` leaves the deeper layers empty; `test_region` refuses `[soil_column].root_beta` and an
    unknown `conductance`.
- **Stomatal closure and the drought-phenology cue used one turgor-loss point for every PFT**
  (#327). Both took ψ_tlp from the shared `[hydraulics]` `leaf_pi0` and `leaf_elastic_mod`, while the
  plant-water solver used a PFT's own values where the PFT file sets them. `pft_leaf_psi_tlp` now
  gives both the PFT's own curve by the solver's rule. With the default traits ψ_tlp stays −1.71 MPa.
  Tests: `test_pft_optics_config` checks a PFT with its own `leaf_pi0` against one without.
- **Stale config text** (#327). The PFT files said a mean leaf inclination of 45° is spherical; a
  spherical distribution has mean 57.3° and standard deviation 21.6°, and the comment in
  `meds_config_pft.toml` and the three example PFT files now says so. Three example configs set
  `[phenology].phenology_on`, which nothing has read since phenology became unconditional; it is
  removed, with the opt-in comment in `meds_config_main.toml`.

- **The BCI flux-tower example forced MEDS with the canopy's upwelling longwave** (#326). `BCI_v5.1.csv`
  labels its downwelling longwave `Rl_up` and its upwelling `Rl_dn`, and `bci_site.toml` took the
  label at its word. Three independent checks agree: the provider's `Rnet` equals
  Rs − Rs_dn + Rl_up − Rl_dn (RMS residual 1.7 W m⁻² after the bounds screen, 87.6 as labelled); at
  night `Rl_dn` averages 1.022 times the air's blackbody emission σT⁴ and `Rl_up` 0.948; and
  ERA5-Land's daily longwave follows `Rl_up` (r 0.86) and not `Rl_dn` (0.12). `LWdown` now reads
  `Rl_up`, and both the site TOML and the example README say why. The observed longwave the
  forcing carries averages 429 W m⁻² instead of 466, and the synthesis regression is re-fitted to
  it. The five-year census run against the tower, before and after:

  | | before | after | tower |
  |---|---|---|---|
  | net radiation, mean [W m⁻²] | 153.0 | 120.4 | 136.3 |
  | net radiation, night (20–05 h local) | +6.4 | −21.8 | −33.5 |
  | sensible heat, mean | −14.0 | −39.0 | 32.4 |
  | latent heat, mean | 59.0 | 56.5 | 75.5 |
  | GPP, mean [µmol m⁻² s⁻¹] | 10.86 | 11.07 | 7.46 |

  The mislabelled forcing had hidden an albedo error: the model's canopy reflects 0.26 of the
  shortwave against the tower's 0.13, and net radiation now falls short by day. The story the
  example told about BCI's longwave, that it fell with daytime cloudiness so the synthesis's cloud
  term had the wrong sign, came from the same mislabelling and is gone: the fitted cloud coefficient
  is +0.10.

- **Rain arrived one forcing record late on end-stamped files, including every ED_ERA5land run**
  (#320).
  `met_instant` held `rec_prev%rainf` over each interval whatever the stamp convention, while it
  took shortwave from the record whose interval contains the instant (`rec_next` on an `"end"`
  file). On ERA5-Land, whose records are means over the hour ending at the stamp, each hour's rain
  fell in the following hour; totals were unchanged, and a begin-stamped file was read correctly.
  Rain now comes from the same record as shortwave. No test caught it because every fixture's rain
  was zero.
- **The longwave synthesis held the wrong interval's clearness on begin-stamped files** (#320). With
  `lwdown_source = "synthesize"`, `met_advance` remembered the clearness of `rec_next` whatever the
  stamp convention. On a `"begin"` file that is the next interval, dark after the last daylight
  one, so `kt_last_day` was 0 every night and every night took the full cloud term, about +40 W m⁻²
  of longwave at a humid tropical site. Rain, shortwave and the remembered clearness now take their
  record from one function, `interval_mean_record`.

- **Under nvfortran, v0.3.0 could not open a site run, and four tests failed** (#317). v0.3.0
  was verified on ifx and gfortran only. Under nvfortran 25.11, with `MEDS_GPU=multicore` and
  `MEDS_GPU=gpu` alike, 49 of 53 tests passed, and `meds_main` stopped in `driver_open` on any
  config, because code new in v0.3.0 uses four constructs that nvfortran gets wrong. Each
  reproduces in a small standalone program, and each is invisible to ifx. They are now in the
  portability traps of `docs/building.md` and `CLAUDE.md`.
  - `driver_open` reset the run's output buffers with `output_buffers_t()`. nvfortran compiles that
    constructor to an ALLOCATE of a garbage size, for example 72,340,172,838,076,672 bytes, because
    the type has fixed-size array components whose own type has allocatable components. Every
    site run, through `meds_main` or the C API, stopped there. The reset now copies a
    default-initialised local.
  - The CO₂ file reader matched the `timestep` unit with `findloc(UNIT_NAME, trim(unit))`. nvfortran
    returns 0 for a value shorter than the names, so every CO₂ file was refused on its `timestep`
    line. The reader now passes the untrimmed token.
  - `era5land_select_box` built the columns of a box that crosses 180° as
    `[pack([(i, …)], m1), pack([(i, …)], m2)]`. nvfortran returns wrong elements, and the vector
    subscript that follows segfaulted, so a region across 180° could not open. The column indices
    are now a named array.
  - The NaN scan of a loaded ERA5-Land month called `findloc` on a LOGICAL array. The nvfortran
    runtime does not implement that and aborts, so a NaN in the archive stopped the run with
    "FINDLOC: unimplemented for data type" instead of a message naming the variable, hour and
    cell. The scan now searches an integer mask.
  - `test_met_driver` no longer segfaults when a CO₂ file fails to open. It segfaulted because it
    read the unloaded series, and the crash discarded the buffered output that said why the open
    failed. The test now reports FAIL lines instead. With the `trim()` restored, the test reports
    `15 FAILED` and exits instead of crashing.
  - **Checked:** 53/53 on ifx Release and Debug, nvfortran multicore and nvfortran gpu. With
    `OMP_TARGET_OFFLOAD=MANDATORY`, the gpu suite launches the offloaded cohort-update kernel 201
    times. ifx output is unchanged, compared on a 13-month coupled run from bare ground: all
    50,330 values of the variables whose output repeats between two runs of the unmodified
    binary are identical, bit for bit. Four variables in that run do not repeat between two runs
    of the unmodified binary (`soil_temp_site_fast`, `soil_water_site_fast`, `area_patch`,
    `lai_patch`), so they cannot be compared this way. That is a separate defect.

- **Slab output variables switched on by an `[output].io_config` were written from unwritten
  memory, and they corrupted their neighbours** (#318). `manager_setup` sized the shared
  pending-record slab (`max_slab`) from the variables live at that point, and the `io_config`
  overrides ran after it. A config that switches tiers, groups or axes off and picks its variables
  through `io_config` (the file's stated purpose) therefore got a slab too short for them. On the
  example spin-up config with output on, every tier is off, so `max_slab` was 1 while the soil
  slabs need 20 rows and the patch slabs 6.
  - At every period close, `normalize_slab` wrote each slab past its own column of the scratch
    record, into the columns of the variables registered after it. A long enough slab near the
    end of the registry, such as a per-cohort fast variable in a stand with many cohorts, would
    write past the end of the array. The ifx Debug build stops at the first such write:
    "Subscript #1 of the array OUT has value 2 which is greater than the upper bound of 1".
  - The v0.3.0 output queue then kept one row per variable, and the writer read the full slab
    length from it. So `soil_temp_site_fast` layer *i* held `soil_water_site_fast` layer *i* − 1
    and then unwritten memory, up to 1e270 K. `area_patch` and `lai_patch` held denormals, and
    `cohort_count` and `cohort_offset` held the fill value or another patch's value in place of
    their own. The first four differed between two runs of the same binary on the same config.
  - `max_slab` is now computed in `manager_finalize`, the one step every caller (site, region, and
    a region's detail polygons) runs after the overrides and before any buffer is sized.
    `normalize_slab` stops the run if a slab is longer than its record, instead of writing past
    it.
  - The fast tier's soil slabs also folded all `n_soil_layer_max` layers, so the inactive tail read
    back as 0 K and 0 m³/m³. It now folds the active layers only, and the tail is the fill value,
    as on the coarse tiers since #246.
  - `test_output_integrate` runs the driver's order (setup, an override that switches a soil slab
    on, finalize, allocate) and checks every layer in the queued record that the writer reads.
    `test_fast_loop` checks that the fast tier folds only the active layers. Both tests fail with
    their fix reverted, and with the `max_slab` fix reverted the new guard stops the run.
  - **Checked** on a 13-month coupled run from bare ground: the example spin-up config with output
    on and `output_variables.toml`, 411 netCDF files.
    - Two runs are now byte-identical in every file. Before, 409 of 412 files differed.
    - The soil output is physical, 268–325 K and 0.23–0.41 m³/m³, and patch areas sum to 1 in
      every record.
    - Only the six variables named above differ from v0.3.0.
    - 53/53 on ifx Release and Debug.
  - **Not affected:** the example's July stage. With the daily and fast tiers on, its setup
    already sizes the slab to `cohort_max`, and it runs without error even with the fix reverted.
    Its fast-tier soil tail was the 0 described above.

## [0.3.0] — 2026-09-28

A **regional simulation** release. `[run].mode = "region"` runs every ED_ERA5land cell of a
latitude–longitude box as its own polygon, in one process (#289): one forcing reader loads each
month for every cell at once, each polygon takes exactly the step a site run at its cell takes, and
the output is region files with a `polygon` dimension, plus full single-site files for the cells
named in `detail_polygons`. A 1° box around Ithaca, 100 polygons, runs a year in 16 minutes on one
core. The regions are fed by the new global **ED_ERA5land archive**, hourly ERA5-Land at 0.1° in one
file per variable per month, which `scripts/prepare_era5/` downloads and builds (#279–#281, #288)
and which the reader also serves to single sites (`[forcing].format = "era5land"`, #282).

Around it, the forcing now reaches each patch at the top of its own canopy air space (#305), CO₂ is
prescribed as a constant or a time series (#184, #301), the forcing a run used is written to its
output (#293), and every record is dated by the period it covers (#294, #296, #297). The
pre-release review (#314) closed the whole-column energy ledger (#290), made restarts exact (#298),
made the Debug suite and the gfortran builds pass, and removed the `[io]` block (#309).

**Upgrading.** Configs need the new `[forcing]` height keys, and `[site].reference_height` and its
three siblings are refused (#305). A forcing file that carries `CO2air` (#301), a config with an
`[io]` block (#309) and an `[output].io_config` that lists `ground_temp_site` (#275) are refused
too, and `make_forcing_file.py` needs a location (#313). Every entry below that moves a number
states its before and after.

### Added

- **Per-patch forcing output** (#305). Each patch now has its own forcing (see Changed), so the patch
  tier writes it: `wind_cas_top_patch` and `air_temp_cas_top_patch`, with the inputs that produced
  them — `cas_depth_patch` (the canopy-air depth the fast loop used), `rough_patch` and
  `displace_patch`. They go to the daily and monthly tiers, in the `forcing` group, beside the
  polygon echo.

- **Prescribed CO₂** (#184, #301). `[forcing].co2_source` sets the free-atmosphere CO₂, the same way
  for every backend and every polygon:
  - `"const"`, the default, holds `co2_const`, so existing configs run unchanged;
  - `"file"` reads `co2_file`, a MEDS CO₂ file (new module `meds_co2_series`).

  The file format (format 1, `docs/science/forcing.md` §12) is plain text:
  - a `timestep <n> year|month|day|hour|minute` line and a `units umol/mol` line;
  - then one `<period start> <value>` row per period, each value the period's mean.

  Each value sits at its period's middle, and the CO₂ is linear between middles on **model** time,
  so it keeps rising while recycled met repeats. Rows must be consecutive, and a run the file does
  not cover stops at startup.

  The repository now ships `data/co2/co2_cmip7_global_annual_1000-2022.txt`: CMIP7 input4MIPs
  `CR-CMIP-1-0-0` global annual means for 1000–2022 (CC BY 4.0). Its header carries the format, so it
  doubles as a template. `scripts/prepare_co2/make_co2_file.py` builds it from the ESGF files.

  Checked on the r1 regression cases:
  - every case at `co2_const = 420` is bitwise unchanged against `beta`;
  - a two-year run from 1990, with recycled 2024 met and the shipped file, writes an `atm_co2_site`
    equal to the series at the sub-step samples (to 2e-13 µmol/mol). The air temperature repeats
    exactly from 1990 to 1991, while the CO₂ rises 1.1–1.3 µmol/mol.

- **`docs/science/order_of_processes.md`, the order of processes** (#300). It covers:
  - a run, and a region's month loop;
  - one slow step (`polygon_step`), and the fast loop's per-step setup and per-sub-step sequence;
  - the slow dynamics;
  - the calendar boundary's restructuring and what it means for records and checkpoints;
  - output and I/O, with a code map.

  No other page gave the sequence: each science page covers one process.

- **The forcing the run used is written to the output** (#293; §6.7 of `MEDS_FORCING_DESIGN.md`), as a new
  `forcing` group (`[output].forcing`, on by default). It is recorded after the reader's shortwave
  partition, rain/snow split and optional corrections:
  - polygon means at the daily, monthly and yearly tiers: air temperature, specific humidity,
    surface pressure, wind speed, downward longwave, the four shortwave streams, snowfall, CO₂, the
    solar zenith cosine and air density;
  - the same fields plus liquid rain at the fast tier, one value per record from the sub-step
    samples themselves.

  The means are accumulated once per polygon, in a new polygon diagnostic block (`site%diag`, fields
  `PY_*`, beside the patch block's `PD_*` and the cohort blocks' `CD_*`/`CS_*`), since the forcing is
  the same across the polygon. `sw_in_site` and `precip_site` moved there from the patch block. That
  fixes them on a year-boundary step: a patch created by the year's disturbance had a cleared
  diagnostic slot and counted as zero in the area-weighted mean. In one case the daily record of
  31 December read 272.31 K instead of 276.15 K, and the January means were low by a 1/31 share of
  that. Every other existing variable is unchanged (six regression cases, 1,267 files).
  `meds_io_config.toml` is regenerated, and the variable inventory in `docs/science/diagnostics.md`
  is recounted (252 variables in this release).

- **`scripts/prepare_era5/make_forcing_file.py`** (#291) writes the single forcing file
  (`[forcing].format = "netcdf"`) from either input:
  - **an ED_ERA5land archive** (`--data-path`, `--start`, `--end`): it picks each site's nearest
    valid cell within `--max-distance-km`, as the model's reader does, and records the cell's
    orography. A file cut this way reproduces a run on the archive itself to about 1e-7 (July 2024,
    Ithaca), with identical precipitation;
  - **box files** from `postprocess_era5land.py --split none` (`--box-dir`), for a small download
    without an archive.

  One behaviour change for box files: only negative packing noise is clipped, the rule the archive
  already follows (`MEDS_FORCING_DESIGN.md` §7.3). The old threshold zeroed every hour with up to
  0.01 mm of rain; at Ithaca in 2024 that was 0.38% of the year's precipitation. Otherwise the
  output equals `prep_era5land_forcing.py`'s, variable for variable.

- **Region runs** (#289; R2 of `MEDS_POLYGON_RUNTIME_PLAN.md`). `[run].mode = "region"` simulates every
  selected ED_ERA5land cell of a `[region].box_nwse` as its own polygon, in one process:
  - one forcing reader serves all polygons, loading a month for every cell at once;
  - each polygon runs the site run's own step, so it computes exactly what a site run at its cell
    computes;
  - the output is region files with a `polygon` dimension (site totals, per-PFT, per-size-class and
    per-soil-layer variables), plus ordinary single-site files for `detail_polygons`;
  - cells below `land_fraction_min` (default 0.5) are skipped.

  A new CTest, `region`, checks a 3-polygon region on a synthetic archive against 3 site runs, bit
  for bit, and checks the region-mode config rules. Regions start from bare ground and write no
  checkpoints until region restarts exist (R4). A 1° box around Ithaca (100 polygons) runs 2016
  in 16 minutes on one core, at 0.81 s per polygon-month: the same cost as the physics of a site run.
  Memory is about 1.3 MB per polygon. See `docs/configuration.md`, "Regional runs".

- **ERA5-Land forcing tools in `scripts/prepare_era5/`** (#279). Downloading and post-processing are
  separate tools:
  - `download_era5land_gdex.py` fetches NSF NCAR GDEX d633008's global 5-day files unchanged into a
    raw pool, in parallel (capped at GDEX's per-user limit of 10 streams), verified and resumable;
  - `download_era5land_cds.py` fetches the Copernicus CDS, in GRIB by default: one variable ×
    12 months per request, which is half the requests NetCDF needs;
  - `postprocess_era5land.py` turns either source's raw files into NetCDF box files.

  They run in their own `meds-era5` environment (`scripts/prepare_era5/environment.yml`). On a New
  York State box, the CDS and GDEX outputs agree to within 0.00024 K. They replace
  `scripts/download_era5land.py` (see Removed).
- **Forcing-data design and a polygon runtime plan** (#279).
  - `MEDS_FORCING_DESIGN.md` gains Part II (§11–§19): a global per-variable monthly `ED_ERA5land_`
    archive, and a reader upgrade that reads a site or a box from it one month at a time and
    converts dewpoint to specific humidity inside the model. Its header now marks what is already
    implemented.
  - `MEDS_POLYGON_RUNTIME_PLAN.md` designs regional runs as an OpenMP loop over polygons, without
    MPI.
- **The global monthly `ED_ERA5land` archive builder** (#280).
  - `build_era5land_static.py` writes the static file: a valid-data mask defined from the data,
    ERA5-Land orography and land fraction.
  - `build_era5land_archive.py` writes one global file per variable per month,
    `ED_ERA5land_<Var>_<YYYYMM>.nc`, flat in the archive folder:
    - `Tair`, `Tdew`, `PSurf`, `u10`, `v10` as delivered; `Rainf`, `SWdown`, `LWdown`
      de-accumulated;
    - month-long chunks, and quantized;
    - it checks every hour against the static mask and plausibility bounds, records each file in a
      manifest with its checksum, and deletes raw files once verified.
  - It reads GDEX raw files or global CDS GRIB (`--source cds`). The CDS path indexes the GRIB
    headers and decodes only the valid cells, one field at a time.
  - The July 2022 pilot (GDEX) built in 457 s on 8 cores (16.0 GB). Its New York values match the
    independent box output to 0.0039 K, the quantization.
  - June 2022 was built from CDS in 411 s on 8 cores; its `Tair` and `Rainf` are bit-identical to a
    GDEX build of the same month.
- **`download_era5land_cds.py --bbox global`** (#280) requests the native global grid, and `--parallel`
  (default 3) keeps several requests in the CDS queue at once (#280).
- **The forcing reader reads the ED_ERA5land archive** (#282; `[forcing].format = "era5land"`,
  `MEDS_FORCING_DESIGN.md` §15). A run names the archive folder in `data_path`, and the reader:
  - binds the site to its cell, or to the nearest valid cell within `max_distance_km` when the site
    falls on a no-data cell, and takes the cell's elevation from the archive's static file;
  - reads one month at a time, one chunk column per variable, over the recycle window or the run
    period, and checks at open that every file the run needs exists;
  - converts dewpoint to specific humidity with the model's own saturation curve, and forms the
    wind speed from the stored components.
  - Keys that do not apply to the archive (`path`, `grid_index`, `grid_match`,
    `[site].grid_elevation`) are rejected rather than ignored. Path keys hold up to 1024
    characters.
  - Box selection, including boxes across 180°, is what `[region].box_nwse` selects (see Region
    runs).
  - With the archive, every output file carries a `forcing_qair` global attribute naming the
    humidity formula, because a run's humidity then depends on the model version.
- **The forcing record carries the wind vector** (#282; `wind_u`, `wind_v`, `has_wind_vector`) beside
  the speed when the source supplies components: the archive always, a MEDS forcing file when it
  carries `u10` and `v10`. The components interpolate linearly and take the same height correction
  as the speed, so the direction is preserved.

### Changed

- **The forcing is moved to each patch's canopy-air top** (#305; `docs/science/forcing.md` §8). All
  vertical corrections now live in `src/forcing/meds_lapse_rate.f90`.
  - **Each sample goes to the top of each patch's canopy air space,** per patch and sub-step, instead
    of being applied at one fixed reference height. That top grows with the stand.
    - Potential temperature and specific humidity are conserved.
    - The wind follows the patch's own log profile, the one its aerodynamics starts from.
    - ERA5's 10 m wind, an open-terrain diagnostic, is first returned to its 40 m blending height.
    - Pressure stays at the ground.
  - **The forcing's own heights are declared,** in new required `[forcing]` keys: `tq_height`,
    `wind_height`, `height_above` (`"zero_plane"` or `"ground"`) and `wind_exposure`
    (`"open_terrain"` with `wind_exposure_z0` and `wind_blending_height`, or `"local"`).
    `[site].reference_height`, `wind_meas_height`, `apply_wind_profile` and `wind_roughness_z0` are
    removed and rejected. Nothing has to clear the canopy any more.
  - **The terrain lapse is broadened:**
    - specific humidity at constant relative humidity (it was held fixed);
    - file longwave by the clear-sky ε·T⁴ ratio;
    - `[site].lapse_rate_tair` accepts twelve monthly rates.
  - **Before and after**, on the r1 cases at Ithaca:
    - **Established stand, annual mean:** friction velocity 0.281 → 0.375 m/s (+33%), sensible heat
      2.93 → 0.62 W/m², latent heat +0.36 W/m², GPP +0.6%.
    - **Regrowth, July:** friction velocity +21%, sensible heat 4.5 → 11.8 W/m².
    - **With the terrain lapse on** (the site is 47.5 m below its ERA5-Land cell): +0.31 K air
      temperature, +2.2 W/m² longwave and +0.56 kPa pressure.
    - **Unchanged:** water still closes to 10⁻¹¹ kg/m². The energy ledger leak of #290, fixed later
      in this release, was about the same size (−0.31 → −0.38 W/m² in the year run).
    - **Slow-only runs are bitwise unchanged.**
  - **Configs migrated:** the examples, `meds_config_main.toml` (with ERA5-Land's values and the
    Kunkel 1989 monthly rates) and `meds_io_config.toml`.

- **CO₂ no longer comes from the met file** (#184, #301).
  - A MEDS forcing file that carries `CO2air` is rejected at open. Before, its value overrode
    `co2_const` without a message, and it repeated with recycled met. Drop it with
    `ncks -x -v CO2air in.nc out.nc`, or rebuild the file.
  - `make_forcing_file.py` writes no `CO2air`, which was a constant (`--co2`, default 420), and no
    `co2_const` or `co2_note` attribute. The `--co2` option is gone.
  - The `const` backend uses `co2_const`. It returned 420 whatever the config said.
  - The test forcing files lose `CO2air`. `test_fast_loop` sets `co2_const = 415`, the value its file
    carried, so its numbers are unchanged.
  - In `docs/science/forcing.md` the longwave section is renumbered §11 (it was a second §8), ahead of
    the new §12.

  Before and after: the r1 `bare_july_single` case, whose file carried `CO2air = 420` with
  `co2_const = 420`, is bitwise unchanged once the variable is dropped.

- **The stand's calendar restructuring runs between two slow steps, after the output has read the one
  that ends on the boundary** (#297). Monthly recruitment, cohort fusion, culling and splitting, and
  yearly patch disturbance and fusion, used to run inside the step that ends on the boundary, after
  its fast loop and before its patch ageing, canopy-depth refresh, soil biogeochemistry and output
  tick. Now `polygon_step` sets `restructure_pending` at the end of that step, and the next step
  restructures (`advance_boundary`) before its fast loop. `restructure_stand` is split out of
  `vegetation_dynamics` for this.
  - **Every record now describes one stand, and a boundary's events belong to the period it opens.**
    The #296 exception is gone: monthly cohort and patch records are exactly the mean of their daily
    records, and the writer's queued-record axis sizing is reverted. On `est_year`:
    - the year's disturbance moves from the December 2074 record (#296) to January 2075
      (`disturb_area_site` 0.1638);
    - the biomass it removes first shows in the 1 January daily state (10.0037 → 9.8645 kgC/m²);
    - `gpp_site` no longer misses the killed canopy's last day. It read 1.39 % below the patch-sourced
      `gpp_rate_site` on the disturbance day; the ratio is now constant to 1.3e-4 on every day.
  - **The model moves slightly.** The fast loops see the same stands as before. What moves is the
    boundary step's patch ageing, canopy-depth refresh and soil biogeochemistry, which now run on the
    stand before its restructuring. On `est_year`, monthly values move by at most 4.5e-4 in `agb_site`,
    2.2e-3 in `nplant_site` and 4.0e-4 in `soilc_total_site`, and fluxes by 1e-9. In a bare-ground
    run the first recruits appear one period later in the records: `demography_30yr`'s year-2000 AGB
    (1.3e-6) now belongs to 2001.
  - **Restarts.** A checkpoint on a boundary holds the stand before the restructuring, and records
    what is owed as the global attribute `restructure_pending` (`none` | `month` | `year`). A resumed
    run performs it first. A state file without the attribute, which includes every file written
    before this change, owes none. Resuming `ckpt_2yr` from its 1 January 2025 checkpoint reproduces
    the continuous run's January disturbance exactly. The differences that remain are fast-tier
    fluxes on the first resumed day (up to 0.03 W/m²) and the integrator's work counters. Both
    predate this change and are a subset of #296's.
  - **Diagnostic blocks.**
    - They are zeroed after the output tick (`reset_step_diagnostics`) instead of at the start of
      the fast loop, so a boundary's events reach the next record.
    - They are zeroed at allocation.
    - `cohort_diag_grow` and `patch_diag_grow` keep every existing slot, not just the first `n`. A
      resumed run's restructuring writes before any fast loop has set `n`, and the grow used to drop
      those rows.
  - **Interfaces.**
    - `advance_one_step`, `advance_slow_dynamics` and `vegetation_dynamics` no longer take
      `is_new_month` / `is_new_year`; `advance_boundary`, `advance_boundary_dynamics` and
      `restructure_stand` do.
    - The C API's `meds_advance_slow` keeps its signature and restructures after the step.
    - `output_integrate` takes the step's start only.

- **The fast loop reads the forcing record directly** (#292; §6.2 and Q6 of `MEDS_FORCING_DESIGN.md`).
  `fill_forcing` and `fill_aenv` take the sub-step's `met_forcing_t`: the reader's sample, or the
  context's reference climate without a forcing source. The `apply_met_to_ctx` shim, and the
  per-thread copies of the whole fast context it wrote into, are gone. Outputs are identical in the
  six regression cases and with 1 or 4 patch threads.

- **`MEDS_POLYGON_RUNTIME_PLAN.md` revised** (#283). A run's polygons form a *region*, always contiguous;
  scattered site networks run as separate processes (a job array, or one allocation filled with GNU
  parallel) rather than in one process. The plan adds the output-performance analysis, two blockers
  found in the code (one pending record per output frequency; site location in the shared config),
  and a step-by-step plan for R0 (measurements), R1 (the compute/I-O split) and R2 (the serial region
  container). R0 is measured: an established Ithaca stand costs about 1.1 s and 0.7 MB per
  polygon-month and polygon, 18.6% of CPU goes to allocation, and a site run spends 24% of its time
  in archive reads.
- **No netCDF inside a time step** (#285; R1 of `MEDS_POLYGON_RUNTIME_PLAN.md`). A step now computes only;
  file work happens in an I/O phase when a calendar month closes and at the end of the run:
  - closed output records wait in per-tier queues (stored compactly) and are written then;
  - the yearly checkpoint moved into that phase;
  - forcing is in memory before each step: the ED_ERA5land reader loads the month a step reads plus
    the record before it, reading each month once per pass, and a MEDS forcing file's records for
    the run are read at open.
  Outputs are identical. Two visible changes: `[output].sync_every` now takes effect at month
  boundaries, so a crash loses at most the current month's output; and `format = "era5land"` needs
  daily steps from midnight (`dt_slow = "1d"`, `start_time` at 00:00:00).
- **The forcing reader is split into a shared source and per-polygon cursors** (#286; R2 of
  `MEDS_POLYGON_RUNTIME_PLAN.md`). `met_source_t` holds the file, the time axis, the cells and the
  loaded month; `met_cursor_t` holds one polygon's cell, location and bracketing records, so one
  archive read serves every polygon of a region. `met_open` takes an optional cell list, and the
  new `met_cursor_init` places a cursor on a cell. Leaf phenology takes the polygon's latitude for
  day length and hemisphere, defaulting to `[site]`. Built with `-fp-model consistent`, outputs are
  identical to before in six regression cases. With the default ifx flags, the compiler optimizes the
  day-length call differently, so runs with a growing stand differ at round-off from the first autumn
  on: at most 5e-12 kg C in cohort AGB after one year.
- **The output manager is split into a file set and per-polygon buffers** (#287; R2 of
  `MEDS_POLYGON_RUNTIME_PLAN.md`). `output_files_t` (registry, file settings, streams) is one set of
  output files, shared by every polygon writing into it. `output_buffers_t` holds one polygon's
  reductions, records and fast-tier staging for one file set. The fast loop and the stepper now see
  only the buffers. Outputs are identical.
- **Faster ED_ERA5land reads** (#284). The reader reads the archive as float32, as it is stored, and only
  the cells a run needs within each 16 × 16 chunk, instead of whole chunks converted to double. A
  site year from a spun-up stand runs about 10% faster; outputs are identical.
- **The GDEX downloader retries transient network failures** (#288; `era5land_common.http_download`,
  `http_head_ok`): a reset or timed-out connection, a short transfer, or HTTP 5xx or 429, with waits
  from 15 s to 8 min. A single refused connection used to abort a whole month. That happened to
  every task in two rounds of the archive build, because each failed task started the next straight
  away while GDEX was refusing connections. A missing file (HTTP 404) still fails at once.
- **`build_era5land_archive.py --work-dir`** (#281) writes each output on another disk, such as a compute
  node's local drive, then copies the finished file into the archive in one sequential pass that
  also computes its checksum. Writing HDF5 chunks directly over a network filesystem made builds
  3–4 times slower once several ran at once.

- **`ground_temp_site` is removed, and the four variances say what they measure** (#275,
  #314). `ground_temp_site` read the top soil layer under the name "ground (skin) temperature";
  `soil_temp_top_site` is the same series, so a config or script that used it should switch to that
  name. An `[output].io_config` that still lists it is refused at startup.
  - `cas_temp_var_site`, `soil_temp_top_var_site`, `cas_vpd_var_site` and `leaf_temp_var_site`
    square their partner's end-of-step state, once per slow step. With daily steps that is the
    day-to-day spread at one hour, not the diurnal cycle. Their `long_name`s now say "variance of
    end-of-step samples of …", and `docs/science/diagnostics.md` explains the difference.
  - A skin temperature and within-step variances stay on #275.
  - The registry holds 252 variables; `meds_io_config.toml` is regenerated.

- **The FAST tier stamps each record by its period's start** (review O2, #314), as the daily,
  monthly and annual tiers do since #294. It stamped the instant the period's first forcing sample
  was taken, `forcing_sample_frac · dt_fast` in (7 min 30 s at the defaults), under a `long_name`
  that said "period start". The `time`, `minute` and `second` of every `-F-` record move back by
  that much; no value changes.

- **`met_instant` has one exit** (review F2, #314). The constant backend no longer repeats the
  longwave synthesis and air density in an early return; only the file interpolation is skipped.
  Unchanged bit for bit (`test_fast_loop` block 8 runs the constant backend).

- **`test_disturbance` covers a disturbance at a calendar boundary** (review O5, #314). There the
  step's diagnostics were already read and reset, so every weight is 0 and the patch slots hold only
  that boundary's events; the gap keeps its share of them. The comments on `patch_diag_inherit` and
  its caller describe both that case and a disturbance inside a step (the C API).

- **`make_forcing_file.py` writes the wind vector** (review P4, #314). A file it cut from the
  archive carried only the speed, and floored it at 0.1 m/s as the reader already does for every
  source. It now writes `u10` and `v10` beside the unfloored `Wind`, so the reader takes the speed
  from the vector, as it does with the archive. Files written before still load as they did.

- **The ERA5-Land tools define what they share once** (#313, #314), in `era5land_common.py`:
  - the de-accumulation (`deaccumulate`), which the archive builder and `make_forcing_file.py
    --box-dir` both call;
  - the archive's variable table, whose factors, clip flags, units and CF attributes
    `make_forcing_file.py` now uses too;
  - the epoch, the month of an end-stamped hour, the 0…360 → ±180 reorder, the GRIB grid decoding,
    the CDS file names and `PROCESSING_VERSION`.

  The long functions (`build_one`, the builders' and `make_forcing_file.py`'s `main`) are split into
  single steps. No output value changes: the builder rebuilds the published July 2024 `Rainf` and
  `LWdown` bit for bit (1.6 × 10⁹ cell-hours each, attributes and storage too) and a band of `Tair`,
  and every other script writes the data it wrote before. What does change:
  - a forcing file describes each variable as the archive does: `long_name`, `standard_name` and
    `height` are the archive's (`Qair` and `Wind` have their own), and the units are unchanged;
  - box files spell their time units as the archive does (`seconds since 1970-01-01 00:00:00`);
  - the GRIB decoder refuses a scan order other than rows of west-to-east points, which its reshape
    to (`Nj`, `Ni`) assumes, and the GDEX reader checks every raw file's latitudes and longitudes,
    not only the last file's longitudes;
  - `build_era5land_static.py` takes the mask from the file's one (`valid_time`, `latitude`,
    `longitude`) field, not its first non-coordinate variable.

- **`make_forcing_file.py` needs a location** (review P9, #313, #314). `--lat` and `--lon`
  defaulted to Ithaca, so a command that named no site wrote an Ithaca file, and `--cells` or
  `--all-cells` overrode them without a word. The location is now required, in exactly one form:
  `--lat` with `--lon`, `--cells`, or `--all-cells` (box files). `--elevation` (default 320 m, box
  files) is gone: `elevation(grid)` records the source cell's orography, which the archive supplies
  and box files do not, so a file made from box files has none; the model reads
  `[site].grid_elevation`, never this variable. Box files must have an hourly time axis without
  gaps, which the de-accumulation assumes, and a box cell without data stops with an error naming
  it (it used to trim every record away and crash).

- **A failed download stops only itself** (review P7, #313, #314). A GDEX transfer that failed after
  its retries, or any failed CDS request, ended the run with a traceback, and the transfers still
  queued ran on unreported and unverified. Each failure is now reported and logged, the other
  downloads finish, the exit status counts the failures, and a rerun fetches only what is missing.
  The CDS downloader checks its files on the main thread, as the GDEX one does, since HDF5 is not
  thread-safe; the CDS client retries transient HTTP failures itself.

- **The documentation describes v0.3.0** (review §5, §10, #314).
  - Six finished design plans move to `docs/dev_plans/archive/` with tombstones (biogeochemistry,
    snow, the GPU evaluation, the veg-energy plan, the 2026-09-13 docs review, the v0.2 release
    plan); the structure design becomes a Reference document, and the other headers say what is
    open.
  - `src/README.md` counts the tree as built (35.5 k lines, 91 modules, 21 libraries) and names the
    new modules; the rules files and source comments state the current order of the output tick and
    the restructuring, in the present tense.
  - `docs/science/forcing.md` describes the source/cursor reader, `docs/science/order_of_processes.md`
    names `driver_open`, `docs/ed2_comparison.md` is pinned to v0.3.0, and `README.md` lists regional
    runs, the archive and prescribed CO₂.

### Removed

- **The `[io]` config block** (#309, #314), renamed `[state]` in v0.2.0 (#173) and loaded since then
  behind a deprecation warning. A config that still carries it is refused at load, with a message
  naming the keys that replaced it (`io.state_interval_years` is `state.interval_years`); before, it
  loaded. `test_region` checks the refusal.
- **`scripts/prep_era5land_forcing.py`** (#291), replaced by `scripts/prepare_era5/make_forcing_file.py`
  (F5 of `MEDS_FORCING_DESIGN.md`). Its box-file code moved into the new tool. The config comment,
  the forcing README, the science doc, the ED2 comparison and the `example_biophysics` instructions
  show the new commands.
- **`scripts/download_era5land.py`** (#280), replaced by `scripts/prepare_era5/download_era5land_cds.py`
  and `postprocess_era5land.py`, whose box files `make_forcing_file.py --box-dir` reads. The forcing
  README, science doc, config comment and the `example_biophysics` instructions show the new
  commands.

### Fixed

- **Daily, monthly and annual records were dated one slow step late** (#294, #296). The tick ran after
  a step from `prev` to `now`, closed the periods `now` had left, and opened the next window at
  `now`. So a step's fluxes, accumulated over `[prev, now)`, landed in the period that starts at
  `now`. Every daily record held the day before its stamp: a July run wrote records stamped 2 July to
  1 August, and the record stamped *d* matched the hourly mean of day *d − 1* to 1e-13, while against
  the hourly records of its own stamped date it was off by up to 4 K and 202 W/m². Each monthly
  and annual record likewise covered its last day of the previous period, not its own last day.
  - **The fix.** The tick now folds the step into the window that holds `prev`, then closes each
    period `now` has left. A record stamped *d* holds the steps that start in *d*, and a state is the
    value at the end of each of the period's steps. On the six regression cases:
    - each daily record is the old one stamped a day earlier, bitwise;
    - it equals the mean of the fast records of its own date to 1.4e-15;
    - each monthly and annual site value is the aggregate of the new daily or monthly records of its
      period to 1.8e-15;
    - the fast tier is unchanged;
    - no record is stamped at or after the run's end, so the stray file past the end is gone. For
      example, a July run wrote a `D-…08` file before and now doesn't.
  - **Monthly means move by up to** 0.64 K in `air_temp_site` (Nov 2074: 279.80 → 279.16 K), 1.8 W/m²
    in `le_site` (May 2075: 58.35 → 60.13) and 5.5 W/m² in `sw_in_site` (May 2075: 226.9 → 232.3), on
    the one-year established stand.
  - The restructuring at a boundary now runs after the output has read that boundary's step (#297,
    under Changed), so every record holds one cohort/patch slot set and needs no exception.
  - `test_output_roundtrip` checks the daily stamps, and `test_output_integrate` the boundary step.
    The convention is written up in `docs/science/diagnostics.md` §4.

- **Patch-sourced diagnostics read low by the disturbed fraction on the year-boundary step** (#295).
  `apply_patch_disturbance` carves one gap from every donor at a uniform `frac = 1 - exp(-rate dt)`
  and shrinks the donors to `(1 - frac)` of their area, but it *cleared* the gap's patch-diagnostic
  slot, so the gap read 0 for the whole step. Every site mean built as Σ area·value from `PD_*`
  fields — the fluxes, the canopy air, the ground and litter diagnostics — was therefore low by
  exactly `frac` on that step. The gap now inherits the area-weighted donor slot
  (`patch_diag_inherit`), as its canopy air, soil, snow and carbon reservoirs already do, and the
  `/(1 - frac)` that `PD_MORT_C_DISTURB` carried to compensate is removed. On the six regression
  cases only records stamped 1 January change, and no state variable does: the daily record rises by
  1/(1 - frac) = 1.014098 (e.g. `cas_temp_site` 271.75 → 275.58 K on the est_year case), the
  January monthly means by 1.0001–1.0015 (`le_site` 1.000138, `gpp_rate_site` 1.001485,
  `cas_temp_site` 1.000457), `disturb_area_site` by 1.014098, and `mort_carbon_disturb_site` is
  unchanged. `test_disturbance` checks that the gap and the site mean of a two-donor fixture read the
  donors' area-weighted value.

- **ED_ERA5land recycling stopped at the first seam unless the window started at 01:00** (review
  F1, #314). When the seam fell inside a daily step, the one-month buffer held the window's last
  month, so the window's first record had not been read, nor, for a window starting at 02:00–23:00,
  the rest of its first day: `met_driver: a forcing record was not prefetched before the step
  (internal error)`. A window starting at 00:00 on 1 January stopped at 23:00 on the last day of the
  first cycle.
  - The source now keeps the window's first day (its first record through the next midnight, at
    most 24 records per cell) from open. A month whose previous record is the window's first takes
    it from there instead of reading the month before again.
  - Every config in the tree starts its window at 01:00, where the seam opens a step; the six
    regression cases and the regional smoke run are unchanged, bit for bit.
  - A region loads each month's forcing once, before the month, so `validate_config` now requires a
    region's window to start at 00:00 or 01:00 on the 1st of a month. Other anchors used to stop
    inside the first month they affected.
  - `test_met_era5land` walks windows starting at 2021-01-02 00:00, 2021-01-01 06:00 and 2021-02-01
    00:00 across their seams; the synthetic archive gains January 2022. The unfixed reader stops on
    the first of them.

- **gfortran builds again** (review 6.1, #314). Six lines past column 132, a hard error in
  gfortran 11 and 15 (ifx only warns), had stopped every gfortran build since v0.2.2. They are
  wrapped.

- **Under gfortran, a region run segfaulted at the first step after a month's output** (review 6.2,
  #314). The region handed its polygons' output buffers to the serializer as the component
  section `reg%poly(:)%out_bufs`. `output_buffers_t` has allocatable components, and gfortran passes
  such a section through a temporary whose copy-out leaves the buffers' allocations dangling.
  - The buffers now sit beside the polygons as a contiguous `meds_region_t%out_bufs(:)` (a site run:
    `meds_run_t%out_bufs`), and `polygon_step` takes the polygon's buffers as an argument.
  - ifx output is unchanged, bit for bit. gfortran Release and Debug pass `region`.
  - The construct joins the compiler traps in `CLAUDE.md` and `docs/building.md`.

- **A region month in which one polygon failed stopped the I/O phase and lost the month for every
  polygon** (review O1, #314). The serializer required every polygon to hold the same number of
  closed records, but a polygon that fails has closed fewer, and the polygons after it none.
  - It now writes as many records as the longest queue holds, with the fill value where a polygon
    holds none, and takes each record's calendar from the first polygon that holds it.
  - Nothing changes when every polygon completes the month.
  - `test_region` makes one polygon's soil carbon impossible in mid-month and checks the daily file.

- **A reused run inherited the previous run's owed restructuring and output files** (review R1,
  #314). A run that ends on the 1st leaves its boundary's restructuring pending, and the C API
  hands its freed slot to the next run, whose first step then restructured a stand that owed
  nothing. A next run with `[output].enabled = false` also still saw the previous run's files as
  enabled and ticked their buffers. `driver_open` now clears both; a restart still sets the
  restructuring flags from its state file. `test_region` reuses one run across both cases.

- **`met_open` left a rejected MEDS forcing file's handle behind** (review F5, #314). When the
  file's record spacing or attributes contradicted `[forcing]`, the file was closed but its handle
  kept, so a caller's `met_close` closed it again and stopped. The rejection now goes through
  `met_close`, as the other three do.

- **A region reported an unknown `detail_polygons` id only after building every polygon**, and left
  the forcing source open (review R8, #314). The ids are now checked against the box's cells
  before anything is opened.

- **The Debug suite passes: three tests failed under `-check all -fpe0` on every branch** (#308,
  #314), and the first failure was a real configuration gap.
  - **Campbell retention could not be configured.** `validate_config` required
    `[soil_column].curve_par_a > 0` for both retention families, but Campbell's is the air-entry
    suction ψ_sat, which is negative; a positive one made the retention curve raise a negative base
    to a fractional power, so field capacity and wilting point came out NaN. Each family is now
    validated on its own terms (van Genuchten α > 0 and n > 1; Campbell ψ_sat < 0 and b > 0), and
    `test_soil_column_config` builds a Campbell column with ψ_sat = −0.26 m, b = 5.65.
  - **The soil-carbon NaN guard raised an FPE on the NaN it exists to catch:** the NaN test shared
    an `.or.` with ordered comparisons, and ifx at `-O0` signals even on `x /= x`. It now tests
    `ieee_is_nan` first, alone. No result changes.
  - **`test_column_ark` asserted on two variables it never set** (`tw_tiny`, `tw_diag`); the
    diagnostic-wood reference that check was meant to compare against no longer exists, so it is
    gone. The two real assertions stay.
  - ifx and gfortran Debug now pass 53/53, as do both Release builds.

- **The whole-column energy ledger failed at nearly every fast step, by −0.3 to −4 W/m²** (#290,
  #314). Under `[energy].bottom_bc = "dirichlet"`, the shipped examples' choice since v0.2.1
  (#267), the soil solve conducts heat to `deep_temp` across the column's bottom face, but every
  ledger booked that face as the Neumann flux, which is 0. The model was right; the books were not.
  - The BE stage now books the face its implicit solve committed, in the soil and whole-column
    ledgers alike, and RK45 b-weights the face its explicit tendency applied (`soil_energy_time_deriv`
    returns it). Under the adiabatic base the face is 0 and nothing changes.
  - The energy budget now closes at every check of the regression cases:
    - `est_july`: cumulative residual −1.060e7 J/m² (mean −3.96 W/m², 5952 of 5952 checks failing)
      → −6.0e−7 J/m² (none failing);
    - `bare_july_archive` and `bare_july_single`: −1.065e7 J/m² (−3.98 W/m², 17856/17856) →
      at most 1.2e−5 J/m² in size (none);
    - `est_year`: −1.186e7 J/m² (−0.38 W/m², 69978/70080) → 1.7e−5 J/m² (none);
    - `ckpt_2yr`: −1.418e8 J/m² (−2.25 W/m², 210638/210816) → −4.1e−5 J/m² (none);
    - the regional smoke run: 17856 failing checks per polygon → 0.
  - `resid_energy_site` is the only output that moves, by up to 9.6 W/m², to about 1e−12. No state
    variable changes.

- **A slow-only run reported 0 for every patch-sourced rate** (#299, #314). The patch block is
  normalised by the weight the fast loop accumulates, and with `fast_biophysics_on = false`
  nothing added one. The slow step now weights the block with its own `dt_slow` when the fast loop
  is off.
  - On `demography_30yr` the annual `disturb_area_site` is 0.0139 /yr (was 0),
    `nplant_recruit_site` 0.030 plant/m²/yr (was 0), and `mort_carbon_background_site` averages
    0.016 kgC/m²/yr (was 0). Its litter rows stay 0 because that case runs without soil carbon,
    which forms no litter.
  - `test_mortality_pathways` no longer stands in the fast loop's weight for its slow-only check.
  - The other half of #299, fast-only rows that read 0 instead of `_FillValue` in a slow-only run,
    stays open.

- **Every run's first day ran its canopy air at the 20 m type default** (#306, #314). Only the slow
  step set the canopy-air depth, so a run from bare ground, a census or a restart used the default
  until its first slow step, and since #305 the forcing is moved to the top of that depth.
  `polygon_prepare` now sets it from the stand, as the slow step does: the tallest cohort plus the
  freeboard, and at least the floor.
  - `bare_july_single`, day 1: `cas_depth_patch` 20 → 5 m (the bare-ground floor),
    `wind_cas_top_patch` 3.57 → 2.51 m/s, `h_site` 4.22 → 8.27 W/m², `le_site` 19.25 → 19.89 W/m².
    Every later day moves a little; the July means by at most 1.1 % (`h_site`).
  - `est_july`, restarted from a spun-up stand: day-1 depths 20 → 23.23 and 25.51 m; the July means
    by at most 0.8 % (`h_site`), and the integrator's rejection count by 1.3 %.
  - Slow-only runs do not change: their slow step set the depth before the first output.
  - `test_region` checks that a bare-ground stand opens at the floor.

- **A restart was not exact** (#298, #314). The state file did not carry each patch's
  `adapt_dt_last`, the step the fast integrators' adaptive controller last accepted and starts the
  next step from, so a resumed run cold-started the controller and took different sub-steps. It is
  now written, and read when present; an older state file cold-starts the controller as before.
  - `ckpt_2yr` resumed from its 2025-01-01 checkpoint now reproduces the continuous run bit for bit
    in every tier: 546 fast, 18 daily and 18 monthly files.
  - Before, every fast file differed from the first resumed day (`h_flux_fast` by up to 1.8e−5 W/m²,
    `cas_temp_fast` by 3 µK), and so did every monthly file: the integrator's work counters, and
    `global_cohort_id` by one cohort from January 2026.

## [0.2.2] — 2026-09-25

An **open-source and layout** release. MEDS is now licensed under the Apache License 2.0, and the
source tree is laid out so that each folder's purpose is clearer from its name. **Nothing in the
model changes.** The biophysics example's output is byte-identical to v0.2.1, and `libmeds.so`
exports the same C interface.

It is also the first release made through the **`beta` integration branch**. Pull requests now
target `beta`, which collects merged work until a release merges it into `main`; see
[`CONTRIBUTING.md`](CONTRIBUTING.md).

### Added

- **MEDS is open source under the Apache License 2.0** (#277). Before this, the repository had no
  license at all, so nobody could legally reuse the code, while `python/pyproject.toml` declared
  MIT. MEDS now ships:
  - `LICENSE`;
  - a `NOTICE` holding the copyright line, *The MEDS Authors*, and the attribution to ED2, whose
    CC BY 4.0 terms cover the portions MEDS adapts from it;
  - an `AUTHORS` file listing the copyright holders;
  - a `CONTRIBUTING.md`.

  The license also covers v0.1.0 through v0.2.1. Every source file now starts with an
  `SPDX-License-Identifier: Apache-2.0` line. The wheel declares `License-Expression: Apache-2.0`
  and carries `LICENSE` and `NOTICE`, which raises the build floor to scikit-build-core 0.11, the
  first release that reads a license expression.

### Changed

- **The foundation library is one folder again: `src/shared/`** (#276). `base/`, `functions/` and
  `util/` together are `libmeds_shared`, and they were the only library whose sources spanned
  three top-level folders. They now sit under `src/shared/{base,functions,util}`, while `config/`
  and `state/` stay top-level. `src/README.md` states the admission test: a module belongs in
  `shared/` only if it uses nothing outside it. That test keeps the regrouped folder from becoming
  the catch-all that the #125 restructure dissolved. No module, library target or Python name
  changes.
- **`src/capi/` is now `src/c_api/`, and its four shims are `meds_c_api_{leaf,phenology,demography,run}`**
  (#276). Their ctest targets are now `test_c_api_*` (select them with `ctest -R c_api`). The
  `bind(c)` names are unchanged, so `libmeds.so` exports the same C symbols and the Python package
  is unaffected. The ctypes-mirror banner in `python/meds/plant/pheno.py` now points at the
  phenology shim; it named a file an earlier restructure had removed.

## [0.2.1] — 2026-09-15

A diagnostics-and-boundaries release. Nothing here changes the demographic core; what it changes is
**what the model reports and what its lower boundaries do**, and two of the four move numbers a
reader would otherwise have trusted.

The thread that produced it is worth stating, because it is a methodological caution as much as a
set of fixes. A question about why the biophysics example's soil ran warmer than the air turned out
to have a **site-mean artefact** as its main answer: that stand carries a disturbance gap, patch LAI
spanned 0.64 to 5.36, and the gap supplied ~70 % of the apparent anomaly from 21 % of the area. Per
patch, a closed canopy was already behaving correctly — 2 cm soil sitting at the daily-mean air
temperature and running *below* air at midday. **In a demography model a site mean of a nonlinear
surface quantity is not a coarse answer, it is a wrong one**, and the example now plots one patch
rather than the stand. Chasing that question nevertheless turned up three real defects, all fixed
here, and one output gap that hid it (#270, open).

### Changed

- **The shipped examples now use the absorbing thermal bottom boundary** (#267). #145 built
  `[energy].bottom_bc = "dirichlet"`, `deep_temp`, the derived `deep_depth = 3.12 m` optimum and
  `test_soil_annual_damping`, and then nothing switched to it: the default is still `geothermal`,
  which holds the bottom heat flux at zero and makes the base an adiabatic wall that **reflects**
  the annual temperature wave. Measured in `examples/example_biophysics`, switching it moves the
  July soil temperature at the 1.73 m base node from **20.5 °C to 16.2 °C** and its annual range
  from **3.1 K to 1.2 K**; the analytic damping for that depth is 0.42 of the surface amplitude,
  against 0.80 for the adiabatic base and 0.46 for the Dirichlet one. The surface moves only
  0.29 K, so this matters for deep soil temperature — soil-carbon Q10, root-zone temperature, the
  timing of spring thaw — not for the canopy energy balance. `deep_temp` is a SITE constant with no
  defensible global value, so the default is unchanged and `meds_config_main.toml` now documents
  all three keys instead of omitting them.
- **`examples/example_biophysics` regenerated** against that boundary and against the #266 ground-
  evaporation fix. The stand ends at 14 cohorts / 2 patches, peak LAI 4.16, AGB 9.6 kgC m⁻², mean
  dbh 24.7 cm, 16.2 kgC m⁻² soil carbon (was 20 / 3, 4.06, 9.0, 23.5, 15.4). July gross uptake
  442.6 gC m⁻² against 408.5, net 206.6 against 177.7.
- **The example README is trimmed to what a reader needs.** The `dt_fast` 150 s-vs-900 s
  convergence table and the Python-driver-vs-executable agreement study were development evidence,
  not reader-facing; the substance of the latter is already recorded here under #139. The
  reader-facing rules both carried — `dt_fast` is not the output cadence, and long runs compare
  through site aggregates rather than cohort by cohort — are kept. Also corrected: the opening
  sentence said the energy balance is solved every 30 minutes; `dt_fast` is 900 s, so it is 15.

### Fixed

- **Ground evaporation used the wrong air-filled porosity, which inverted its moisture response**
  (#266). `ground_evaporation` referenced CLM5's air-filled pore space to the **bulk** top-layer
  moisture, `phi_air = phi - theta1`. CLM5 eq 5.80 references it to `theta_air`, the **air-dry**
  water content of eq 5.78 — a texture constant obtained by inverting the retention curve at
  `psi = -1e4 m`, not a state. With `phi_air` tied to `theta1` the tortuosity falls as
  `(phi - theta1)^(10/3)`, faster than the dry-surface-layer thickness falls linearly, so `r_soil`
  **rises** with wetness and **a wetter soil evaporates less**: for the shipped loam the ground
  latent flux had a *minimum* at `theta = 0.31`, and the approach to `theta_init` was a 60x cliff
  rather than a limit. At the biophysics example's actual state — `theta` 0.246, `psi` −0.0106 MPa,
  soil-surface relative humidity 0.99992, i.e. water not limiting — `r_soil` was **9167 s/m**
  against an aerodynamic 1389 s/m, pinning ground LE at 3.0 W m⁻² (3.6 % of site ET). Corrected,
  the response is monotone, joins the DSL-free limit continuously, and the closed-canopy floor's
  daily-mean warm bias drops from +1.42 K to +0.99 K. The DSL thickness denominator is corrected to
  `(theta_init - theta_air)` (eq 5.77) at the same time; the tortuosity keeps the Millington–Quirk
  exponent rather than CLM5's `phi_air^2 (phi_air/phi)^(3/B1)`, because `B1` is Clapp–Hornberger and
  the shipped curve is van Genuchten — with `phi_air` constant the two differ only by a constant
  factor, so the response *shape* is unaffected. **Note that simply disabling the DSL also removes
  the bias but desiccates the top layer to a daily-mean −17.8 MPa**: the DSL does real work and was
  referenced to the wrong porosity. `test_column_hydrology` now asserts the physics rather than the
  formula — evaporation monotone in soil moisture, and continuous into the DSL-free limit — and
  both assertions fail on the previous code. Nothing in the suite had constrained this: the only
  soil-evaporation assertions were snow-cover area weighting, which is why a 9x error in ground
  latent flux sat behind a green suite.
- **`examples/example_biophysics` plots the closed-canopy patch, not a site mean across a gap.**
  The stand's patch LAI spans 0.64 to 5.36 and the gap — 21 % of the area, passing ~70 % of incident
  shortwave straight to the ground — behaves like bare soil and runs more than 10 K above air at
  midday, supplying ~70 % of the site-mean soil warm anomaly on its own. The figure now selects the
  highest-LAI patch, and the soil series is relabelled **"Surface soil layer"**: it is the top soil
  layer (node ~1.8 cm, 0–4 cm thick), not a skin or litter temperature, and MEDS has no surface
  organic horizon to confuse it with. Per-patch sub-daily data comes from the opt-in
  `[fast].fast_probe`; the carbon and soil figures remain site means (#270).
- **Sub-daily state variables declared `AGG_TMEAN` were once-per-window snapshots, not time
  means** (#264). A family of `FLD_P_*` output sources read `site` at the output tick — once per
  `dt_slow` — and were then folded by `AGG_TMEAN` as though they had been time-integrated. For a
  quantity with a diurnal cycle that makes the "mean" one instantaneous sample per window, taken at
  whatever local time the boundary falls on, **so the bias is a function of the site's longitude**.
  Measured over one July at Ithaca: `soil_temp_top_site` read **27.79 °C** against the identical
  field accumulated properly at **26.39 °C**, and `cas_temp_site` **23.71 °C** against **22.80 °C**.
  `cas_temp`, `cas_shv`, `cas_co2`, `cas_vpd`, `soil_temp_top` and `w_surface` (site and patch
  axes, nine registry entries) now read the dt-weighted `PD_*` accumulators. After the fix
  `soil_temp_top_site` agrees with the sub-step reference to **0.005 K** and `cas_temp_site` to
  0.03 K. Stocks that are correct as instants — area, age, the soil-carbon pools, SWE, snow depth —
  are unchanged. `PD_CAS_VPD` and `PD_W_SURFACE` are new accumulator rows; VPD is nonlinear in
  temperature, so it is formed per sub-step rather than derived from the averaged twins.
  `test_output_registry` now asserts over the **whole registry** that no `AGG_TMEAN` variable is
  sourced from a once-per-tick state read, so a new variable cannot reintroduce it.
- `air_vpd` and `specific_humidity_to_vpd` moved from `src/io/meds_diagnostic_kernels` to
  `meds_therm_lib`, where the thermodynamics belongs: the fast loop now needs the same formula the
  output layer uses, and a kernel may not depend on the output layer. `meds_diagnostic_kernels`
  re-exports them, so no caller changes. The io layer's private `PRSS_REF` duplicate of `p_std` is
  retired in the same change.

- **Exactly-tied cohorts are now co-dominant in the two-stream, not stacked** (#207).
  `update_overtopping_lai` already treated equal-height cohorts as co-dominant, but
  `canopy_radiation` builds one discrete RT layer per cohort and its input is merely sorted, so
  tied cohorts **stacked** and whichever sorted first was placed above. The tie is exact and
  routine, not hypothetical: `apply_recruitment` uses **one scalar** `recruit_dbh` for every PFT,
  height is re-derived from it, and the sort is stable — so at every recruitment event the winner
  was the **lower PFT index**. Nothing about the ecology chose that. Measured on two optically
  identical tied cohorts, the lower slot absorbed **12.6 % less shortwave at zero overstory,
  rising to 15.1 % under LAI 6** — worst exactly where a whole-canopy albedo or GPP check is least
  likely to notice it, because the total is roughly conserved and only the split between PFTs
  moves. `canopy_radiation` now takes `height` and merges exactly-tied cohorts into one RT layer,
  weighting the blended optics by area index — the same weight `blend_cohort_optics` already uses
  to mix leaf against wood inside a cohort — and scattering the layer's absorbed flux back by
  absorptivity-weighted area, which is identically the quantity `leaf_frac` is built from. A
  one-member layer takes its cohort's values verbatim, so a canopy with no ties is bit-identical.
  After the fix the positional penalty is **0.0 % at every overstory LAI**.
  `test_canopy_radiation` asserts both that optically identical twins absorb equally and that
  permuting the tied pair's PFT order moves neither one; all five assertions fail on the previous
  code. **#1 stays open** — breaking the tie *ecologically* is its standing question, and this
  change does not answer it.

## [0.2.0] — 2026-09-14

**Read this before comparing a v0.2.0 run against a v0.1.0 one.** The release moved real numbers,
and the largest change is that **phenology now runs at all**. See
[`docs/ed2_comparison.md` §0](docs/ed2_comparison.md) for the before/after table.

### Known limitations, stated

Three things this release does **not** fix, carried here because they would otherwise be found by
surprise:

- **Leaf water potential does not converge at the shipped `dt_fast`** (#162, open). Daytime mean
  −0.23 MPa at 12.5 s against −1.19 MPa at 900 s: the plant water-mass update is an explicit step
  with frozen sapflow and root uptake, so the per-step excursion grows with the step. The error is
  **inherited from the canopy air and amplified about 4×**, and the residual relocates to
  `psi_wood` through the frozen uptake seam. Every other state and flux converges. Any study keyed
  to leaf water potential — hydraulic stress, potential-driven mortality — should run at ≤ 150 s
  regardless of what the carbon budget looks like.
- **Phenology is selectable and self-consistent. It is not validated.** No MEDS leaf-area cycle has
  been scored against an observation, at any site, under any strategy. The thresholds are literature
  values for the biome each strategy describes, not site calibrations — selecting `CUE_LIGHT` with
  its default 200 W/m² onset at Ithaca strips the canopy every summer, because temperate summer
  insolation sits above a threshold chosen for a tropical dry season. And it did not run at all
  before this release (#245), so **no MEDS result published before v0.2.0 had a leaf-area cycle**.
- **MEDS has still never been benchmarked.** No EDTS-equivalent regression suite, no site compared
  flux-for-flux, no output scored against observations. What is verified is internal: 49 CTest
  targets on two compilers, per-step conservation ledgers that close to machine precision, a
  per-layer face-closure check, and output byte-identical at any thread count.

Four features that move numbers substantially ship **off**, because the rebaseline window closed
before they landed: Kattge–Knorr thermal acclimation, storage-pool maintenance respiration, leaf
resorption on shed, and the non-stomatal water-stress limb. Per-PFT hydraulic traits are selectable
but uncalibrated, so MEDS ships hydraulically identical PFTs.

### Added

- **An evaluation notebook and a PFT / size-class plotter** (#175). The PFT and DBH-class output
  axes shipped in v0.1 with no reference consumer, and the v0.1 IO plan's worked evaluation against
  the Ithaca test bed was never written. Both now exist in [`post_proc/`](post_proc/).

  `plot_pft_size.py` draws composition by PFT and structure by size class — six panels — and prints
  the closure identities those axes promise (`Σ_pft agb_pft == agb_site`,
  `Σ_class agb_size == agb_site`, and the same for stem density) **before** drawing anything. All
  four hold to roundoff (worst 6.5e-16). A figure built on axes that do not partition the stand is a
  picture of a bug.

  `evaluate_ithaca.ipynb` is a worked pass over one run, using each tier for what that tier can
  answer: diurnal energy and carbon from FAST, seasonal cycles from DAILY, demographic trajectories
  and the conservation identities from MONTHLY/ANNUAL. It is explicitly **not** a benchmark — MEDS
  has never been scored against flux-tower or inventory data — but it checks what is knowable
  without observations, and that turned out to be plenty.

  **Its first run found three defects, all now filed:**

  - **#245 — the phenology cue never reaches the model.** Every PFT has been running evergreen
    whatever it declares, because `load_phenology_pft` skips its whole block unless a key the
    shipped `meds_config_pft.toml` does not document (`flush_cue_mask`) is present, which bypasses
    the presence map and leaves eleven required keys silently missing. The Ithaca reference stand,
    declared cold-deciduous, holds LAI 5.28–5.66 through every January of a 50-year run. Supplying
    the real key names produces a textbook cycle (LAI 0.34 → 5.28 → 2.89), so the phenology works —
    it has never been switched on.
  - **#247 — an out-of-bounds write in the output path.** `cohort_diag_reorder` permutes the
    diagnostic rows but never updates the block's `n`, so after a cull the block's count exceeds the
    live cohort count while `extract_variable` sizes its scratch buffer from the live one. Memory
    corruption in Release. Invisible until now because the reference stand only ever grows.
  - **#246 — the soil axis is padded to `n_soil_layer_max` with `0` / `NaN`** instead of
    `_FillValue`, so reducing over it gives a 136 K soil column, and `soil_psi_site` divides by zero
    on the padding.

  It also documents two things that are easy to get wrong reading these files by hand: turning a
  rate into an amount needs the width of **its own** window (pairing it with the previous window's
  width turned one disturbance event into a 10.8 kgC/m² phantom), and **`litter_*_site` is not all
  the carbon entering the soil** — the cull and disturbance pathways bypass that accumulator, so the
  soil-carbon budget closes from file to **0.28%** with the `mort_carbon_*_site` variables of #169
  and misses by **29%** without them.

- **A per-layer face-closure check, for the one defect class both whole-column ledgers are blind to**
  (#189, item 1 of 3). The fast loop's energy and water ledgers close to machine precision against
  the column boundary — and a purely *vertical* error survives that untouched, because enthalpy put
  in the wrong **layer** still sums correctly. This repository has paid for that three times (PRs
  #77, #85, #86): a scheme advecting soil enthalpy on a mass flux borrowed from the frozen scratch
  solve while committing θ from its own stages, so heat moves for water that never did. What found
  the last instance was an implausible temperature, which is not a detector.

  `faces[soil_layer_mass]` is: per interior layer, |the mass the soil-energy equation was charged
  for − the mass the committed θ actually moved|. It is measured where each scheme commits — on ARK
  from the scratch solve whose θ it commits verbatim, on RK45 at its own b-weighted commit, before
  the clip/floor guards edit θ. Reported next to the two budgets and asserted in the suite.

  **Measured, and mutation-verified.** Clean, it is machine-zero: **2.29e-13 kg/m²** worst over
  420,480 checks (ARK, one Ithaca year) and **2.29e-13** over 207,360 (RK45). Reintroducing the #85
  defect — handing the energy forcing `w_flux_frozen` instead of the stage's own faces:

  | | clean | with the defect |
  |---|---|---|
  | `whole_energy` | 0 fails / 103,680 | **0 fails / 103,680** |
  | `whole_water` | 0 fails / 103,680 | **0 fails / 103,680** |
  | `faces[soil_layer_mass]` | 2.29e-13 | **1.15e-01 kg/m²** |

  Both ledgers are *exactly* as clean with the defect present as without it. The new check moves
  twelve orders of magnitude. And this is at **Ithaca**, where the original temperature symptom was
  invisible (the historical note records 276.49 → 276.50 K) — so the check is strictly more
  sensitive than the accident that first caught the bug.

  Also surfaces `face_mass_resid`, which `advance_soil_water_column` had been computing and
  publishing to a field nothing read.

  **Items 2 and 3 of #189 are not in this change** and the issue stays open for them: per-cohort
  tissue residuals, and moving the RK45 ledger assertion after the rail decision. On the second,
  the concern is confirmed real — `budget_check(..., halt_budgets)` runs inside
  `column_fast_step_rk45`, *before* `column_fast_step` decides whether to roll the step back and
  rescue it on ARK, so under `[energy].debug_error` a step that was about to be discarded can abort
  the run on a ledger breach it would never have committed. The `n_fail` counting is already safe
  (the rollback restores the whole budget); only the hard stop is exposed.

- **The shipped `meds_io_config.toml` example is regenerated, and a test now keeps it that way**
  (#241). That file is generated by `meds_main --dump-io-config` but is also *tracked*, because it
  is the shipped example of the per-variable override surface that
  [`docs/science/diagnostics.md`](docs/science/diagnostics.md) points users at. Nobody regenerated
  it for five pull requests, so it had fallen **twenty variables** behind the registry — the four
  variance variables (#174), `storage_resp_site` (#177), `disturb_area_site` (#170), the five
  top-of-canopy radiative fluxes (#171), the whole five-variable fast-tier CO₂ family, and #169's
  four — plus two stale long names.

  It degrades *quietly*, which is why it drifted: an unknown name in the file is a hard error, but a
  **missing** one is not an error at all — that variable just keeps its registry-default streams. A
  user copying the example as a starting point and switching everything off would silently fail to
  switch off the variables the file forgot.

  New `io_config_example` test asserts every registered variable appears in the tracked file, and
  names the ones that do not. It checks name *presence*, not a byte diff: regeneration rewrites
  every comment and mask line, so a byte comparison would fail on cosmetic churn and end up
  suppressed.

- **Mortality carbon split by pathway** (#169). The output carried mortality *rates* and the total
  litter flux, but nothing separated the three ways a MEDS plant can die, so a stand thinning
  continuously and a stand being knocked over looked alike. Three variables now carry it:
  `mort_carbon_background_site` (the continuous hazard), `mort_carbon_cull_site` (cohorts dropping
  below the tracked size floor), `mort_carbon_disturb_site` (canopy killed when treefall opens a
  gap), plus `mort_carbon_background_patch`. Each is the whole individual — leaf, fine root, wood
  and storage — because that is what a death removes.

  **Emitted whatever `[soil_carbon].soil_carbon_on` says.** Two of the three litter sites sit behind
  that switch, and how much biomass died is a demographic question that does not stop being asked
  when the soil pools are off. The write is placed on the demography side of each guard.

  Measured on a spun-up Ithaca stand (AGB 17.5 kgC/m², `veg_carbon_site` 26.2 kgC/m²), annual means:
  **background 0.463, disturbance 0.348, cull 0.000 kgC/m²/yr** — so treefall is **42.8%** of the
  stand's mortality carbon.

  That headline decomposes rather than standing on its own. The disturbance term is the 1.39%/yr
  treefall hazard applied to 25.07 kgC/m², which is **95.8% of the stand's plant carbon** — i.e. to
  essentially the whole canopy, since only small recruits sit below the 10 m
  `disturbance_survive_height`. And the background term is an effective **1.771%/yr** of the live
  pool against `agb_mort_site`'s **1.758%/yr**, two numbers computed on entirely different code
  paths (the `PD_` patch block and the `CS_` cohort block) agreeing to 0.75%.

  **The cull reading a hard zero is physical, not a silent zero** — the failure mode #170 was. At
  the shipped `negligible_nplant = 1e-8` no cohort ever reaches the floor, because cohort fusion
  consolidates small cohorts long before they decay that far. Confirmed by raising the floor to
  5e-3 in a scratch run, where the variable reports and the stand duly collapses (AGB 16.7 → 1.3).

  The disturbance term is valued on the donor patch's *surviving* ground, so it carries a
  `1/(1-frac)` factor: the carbon died on the `frac` of the donor that became the gap, and the site
  aggregation weights by area at read time, after the shrink. Without it the site total is one-signed
  low by 1.4%/yr. Writing it on the gap instead is not an option — that slot is cleared, and carries
  no weight, for the remainder of the step.

- **Top-of-canopy radiative fluxes on the diagnostic path** (#171). The two-stream forms the upward
  flux every step and the output discarded it, so a run could not be compared against a radiometer or
  a satellite product without re-deriving it offline. Five variables now carry it: `sw_in_vis_site`,
  `sw_in_nir_site`, `sw_up_vis_site`, `sw_up_nir_site`, `lw_up_site` (surface emission included).

  **Fluxes, not a time-averaged albedo — deliberately.** A period-mean albedo is the mean of a
  *ratio*, which is not the ratio of the means, and at night the shortwave ratio is 0/0. The albedo a
  reader wants, and what a satellite product *is*, is `Σ up / Σ down` over whatever period they
  choose. Emitting a mean albedo would have produced a number that looks right and is wrong.

  Measured on a spun-up Ithaca stand, annual sums: **VIS 0.083, NIR 0.320, broadband 0.208**. The
  VIS/NIR contrast is the vegetation red-edge that every vegetation index is built on. Winter months
  rise to 0.23 / 0.45 as the snow albedo ramps in, and `lw_up_site` runs 305–468 W/m², consistent
  with σT⁴ across the annual surface-temperature range — so the numbers are checkable against
  something external, which is the point of the issue.

  No second RT solve: `canopy_radiation` already forms `albedo(b)`, so the upwelling is that ratio
  times the incident it was formed from.

- **Variance output** (#174). `AGG_MEANSQ` had existed since the IO design with no consumer:
  `normalize_scalar` computed a variance into an optional argument that nothing passed, and
  `normalize_slab` had no variance path at all. It is now `AGG_VARIANCE`, emitting the dt-weighted
  $\langle x^2\rangle - \langle x\rangle^2$ directly, with the dead `out2`/`has2` plumbing removed.

  **Registered as ordinary variables sharing their partner's source id**, not as a companion slot
  bolted to the mean. Two consequences, both deliberate: the existing per-variable buffer / normalize
  / serialize path carries them with **no new machinery** — no change to the pending record, the
  serializer, or `close_tier` — and each variance is **independently switchable** through the
  `[variables]` override, exactly like every other output. A bolted-on slot would have been neither.

  Four ship, monthly and annual only, off unless requested: `cas_temp_var_site`,
  `leaf_temp_var_site`, `soil_temp_top_var_site`, `cas_vpd_var_site`. Units are the partner's
  **squared**; take the square root for a standard deviation. Emitting the standard deviation
  directly would have lost the additivity that lets a variance combine across periods.

  Measured on a spun-up Ithaca stand: canopy-air standard deviation **2.3 K in July against 5.4 K in
  December**, and VPD 341 Pa in July against 70 Pa in January. A monthly mean cannot tell you that,
  and it is most of what a sub-daily process actually sees.

  The variance is floored at zero — the two moments accumulate independently, so round-off can put
  the difference a hair below zero for a near-constant series, and a negative variance in an output
  file is worse than a zero. `test_output_integrate` covers the dt-weighting (18.75, not the
  equal-weight 25) and the constant-series case.

### Changed

- **`[io]` is now `[state]`** (#173), with a deprecation path. The block was named for a legacy
  diagnostic writer retired at v0.1; what remained was the restart stream, so `io` named the one
  output path it did *not* cover — every diagnostic goes through `[output]`.

  | old | new |
  |---|---|
  | `[io].output_dir` | `[state].output_dir` |
  | `[io].output_prefix` | `[state].output_prefix` |
  | `[io].write_state` | `[state].write_state` |
  | `[io].state_interval_years` | `[state].interval_years` |

  The last one loses its `state_` prefix, which was stuttering once the block itself is called
  `state`.

  **`[io]` still loads** in v0.2.x, printing **one** deprecation warning for the block rather than
  one per key, and will be removed in a later release. A config that silently stopped being read
  would fall back to whatever the missing-key report produced rather than to the values the user
  wrote, which is why this is a shim and not a straight edit. Verified byte-identical: a run under
  each spelling produced 4 identical netCDF files.

  `test_capi_run` now derives a legacy-spelling config **from the shipped one** and asserts the two
  parse to the same values — so the deprecation path cannot rot into a block that reads nothing, and
  the test cannot drift from what the shipped config actually contains.

  Internally `cfg%io_*` became `cfg%state_*`. One incidental fix: `write_derived_config` in
  `test_capi_run` emitted `[io] write_state = false` and then appended the shipped config, so once
  the shipped file moved to `[state]` the new key would have won and silently flipped it back on.

### Added

- **Per-PFT hydraulic traits** (#179). All thirteen — the leaf and wood pressure–volume curves, the
  xylem vulnerability, and the conductance parameterization — can now be given per-PFT arrays in the
  `[pft]` table under the **same key names** the `[hydraulics]` block uses. Until now every PFT shared
  one parameter set, so wood density was the only axis on which PFTs could differ hydraulically, in a
  model whose point is that plant strategies differ.

  **All optional**: an absent key falls back to the `[hydraulics]` scalar, so a config can make *one*
  trait per-PFT without restating the other twelve. A run that sets none is **byte-identical** — 13
  netCDF files compared against a binary built from the parent commit.

  **Each table entry carries its own Kirchhoff lookup**, rebuilt from that PFT's `wood_kexp`. Sharing
  one across PFTs would have given every PFT the first one's vulnerability *shape* while every budget
  still closed — the defect class this repository keeps finding.

  **`solve_plant_water` itself is untouched.** The batch selects the cohort's table entry and hands
  the solver the same `hydro_params_t` it always took, so the numerically delicate part of the
  hydraulics did not change at all; only the parameter selection moved.

  **The first cut segfaulted five tests**, because the fixtures build `col_config` by hand and did not
  populate the new table — the fifth instance of the fixture-does-not-mirror-the-driver trap this
  release. Rather than patch five fixtures, `apply_hydraulics_config` now *is* the table builder and
  there is deliberately **no** PFT-uniform companion field: a fixture that forgets it fails to
  **compile**. That converts a runtime segfault into a build error, which is the only durable answer
  to a trap that has recurred five times.

- **Longwave synthesis for a forcing source without `LWdown`** (#182). `lwdown_source =
  "synthesize"` was a declared selector that parsed and then took longwave from the file anyway;
  PR #149 made `validate_config` reject it rather than let it lie. It is now implemented and the
  rejection is gone:

  ```
  LW↓ = ε_clear · σ · Tₐ⁴ · [1 + a·(1−kt)]
  ```

  `lw_clear_form` selects Brutsaert (1975), a power law in screen-level vapour pressure, or Idso &
  Jackson (1969), temperature alone — the fallback when a source's humidity is not trustworthy. `kt`
  is the **same** clearness index the shortwave partition uses, exposed rather than recomputed, so
  there is one definition of "how clear is the sky" in the forcing path. `lw_cloud_a` defaults to
  0.22.

  **The cloud term is not cosmetic.** Against the Ithaca ERA5-Land file, clear-sky Brutsaert alone
  underestimates `strd` by a mean of **−29.9 W/m²** (RMSE 43.6) on a file mean of 312.3 — clear-sky
  formulations miss the cloud enhancement and a temperate site is cloudy most of the time.

  **At night `kt` is undefined**, and `clearness_index` returns a **negative sentinel** rather than
  zero: zero would read as fully overcast and give every night the maximum cloud correction. The
  driver holds the last *daytime* index through the dark. That scalar is updated in `met_advance`,
  not `met_instant` — `met_instant` is `intent(in)` and runs per sub-step, while the advance sits
  outside the patch loop, so the state cannot become a data race when the patch axis is threaded.

  **It is a fallback, not a substitute**: driving Ithaca from synthesis instead of the file's `strd`
  leaves the soil surface **1.37 K cooler** in the annual mean. Use the file's longwave when the file
  has it.

  `LWdown` is now read *optionally* when synthesizing — demanding a variable the feature exists to
  do without would have defeated it. The two selector codes live in `meds_forcing_config`, beside
  `LW_FILE`/`LW_SYNTHESIZE`, not with the kernel that evaluates them: `meds_forcing` links
  `meds_config`, so putting them in the kernel would have pointed the config layer at the forcing
  layer and closed a dependency cycle.

- **Leaf resorption on shed** (#151), per-PFT `retained_carbon_fraction`. A fraction of the **active**
  (senescence) leaf shed returns to the non-structural pool instead of entering litter; until now
  every gram of shed leaf carbon became litter, which overstates litter input and understates the
  plant's retained reserve.

  **The leaf pool loses the full shed either way** — the design note warns that crediting storage
  while removing only the litter share *creates* carbon, so `npp%leaf` keeps the complete removal and
  the split happens between the two destinations:

  ```
  leaf   -= S_base + S_act          storage += f · S_act          litter = S_base + (1−f)·S_act
  ```

  **Only the active excess is resorbed, not the baseline.** The shed rate is
  `max(k_shed·drive, k_turn)` — a *max*, not a sum — so the active excess is exactly
  `shed_rate − base_rate`, which decomposes the max correctly in both regimes. Baseline turnover is
  excluded because `leaf_turnover_rate` is calibrated against observed **litterfall**, which already
  has resorption in it; resorbing it again would double-count. `pheno_drives_to_rates` now reports the
  baseline share so the carbon layer can make that split.

  **Default 0**, reproducing the earlier behaviour. Most measured resorption is of N and P rather than
  C, so carbon fractions are modest — 0.1–0.2 is defensible. On a **temperate-deciduous** stand at
  `f = 0.35`, Ithaca, 5 years:

  | | f = 0 | f = 0.35 | |
  |---|---|---|---|
  | storage pool | 0.000424 | 0.000688 | **+62.2 %** |
  | soil carbon | 0.001516 | 0.001449 | **−4.4 %** |
  | NPP to storage | 0.000375 | 0.000286 | −23.9 % |
  | GPP / LAI | | | +22.0 % / +12.2 % |

  The storage/soil-carbon pair is the direct signature — carbon moved from the litter path to the
  plant reserve — and less NPP then has to be diverted to refill storage. On the **evergreen** shipped
  config the effect is *exactly zero*, because that stand never actively sheds; that is the intended
  behaviour and `test_plant_phenology` asserts it.

  The slow-loop carbon ledger closes at **−4.9E-17** against a declared 7.41E-02 with `f = 0.35`.

- **Storage-pool maintenance respiration** (#177), per-PFT `storage_turnover_rate` [yr⁻¹]. The
  non-structural pool was the one live carbon store that cost nothing to hold — a cohort could carry
  an arbitrarily large reserve for free. It is now charged a fractional turnover each step, ED2's
  `growth_balive.f90` form:

  ```
  M_s = C_s · min(1, storage_turnover_rate · dt)
  ```

  **No temperature dependence**, because ED2 has none here: it sets `maintenance_temp_dep = 1.0` for
  storage and leaves the temperature form commented out as "experimental and arbitrary". Adopting one
  would go beyond the reference rather than follow it.

  **Default 0**, which reproduces the earlier behaviour — verified not by a byte compare (a new
  output variable changes the file layout) but by comparing **1677 variable instances** across 13
  monthly files against a binary built from the parent commit: worst absolute difference **exactly
  0.000e+00**, with `storage_resp_site` the only addition. ED2's own values are temperate broadleaf
  0.6243, temperate grass and conifer 0, tropical non-grass 1/6, tropical grass 1/3 — so zero is a
  legitimate ED2 setting, not an absence of physics. Where it is turned on the charge is large: at
  ED2's temperate-broadleaf rate an Ithaca run loses **35 % of GPP and 43 % of AGB** over five years,
  because the drain compounds through stand development. Whether it should be on by default is a
  v0.3.0 rebaseline question, recorded on `docs/ROADMAP.md`.

  **The ledger caught the first implementation.** Decrementing `nonstructural_carbon` in place — the
  obvious reading of ED2, which does exactly that — left the slow-loop carbon ledger with
  **−2.4339E-03 kgC** undeclared in the `allocate` phase and **+2.4339E-03** over-declared in
  `grow+mortality`, equal and opposite. The phase that *declares* the CO₂ efflux has to be the phase
  where the pool drops. Routing the charge through `npp%nonstructural` as a **tendency** — the
  standing "driver computes, engine applies" rule — closes it: the residual is now **−6.0E-17**
  against a declared 5.13E-02. The allocator still sees the post-maintenance reserve, so ED2's
  maintenance-before-growth ordering is preserved even though the pool moves one phase later.

  Reported as `storage_resp_site` [kgC/m²/yr], and it rides `co2_owed` — the same per-patch channel
  growth respiration already uses to reach the fast loop's `nee_biotic`.

### Changed

- **Fine-root maintenance respiration is summed over soil layers, not taken at a mean temperature**
  (#178). The model resolves a soil temperature per layer; the root response collapsed it to a
  root-weighted mean *before* the temperature function, throwing that resolution away. It now
  evaluates `Σ_k root_frac_k · f(T_k)` instead of `f(Σ_k root_frac_k · T_k)`.

  This is the Jensen argument #145 makes for Rh, with one extra turn worth stating: the **peaked**
  Arrhenius form is convex below its optimum and **concave** near it, so the error changes sign with
  season instead of biasing one way. Measured at Ithaca on the default 2 m column with `root_beta =
  2`:

  | | Dec | Mar | Jun | Sep | annual mean |
  |---|---|---|---|---|---|
  | layered vs mean-temperature | **+4.8 %** | +1.8 % | **−8.3 %** | −0.3 % | **−0.17 %** |

  So it is a **seasonal** correction to root respiration, not an annual-budget one — reporting only
  the annual figure would understate it by a factor of thirty.

  The weighted scale is patch-uniform, so it is evaluated once per patch and broadcast over the
  cohort array; `fine_root_maintenance_respiration` accordingly takes a temperature *scale* rather
  than a temperature, which keeps it `elemental` over cohorts. `test_plant_respiration` asserts that
  a **uniform** profile reproduces the old answer exactly (so the change is inert without a
  gradient) and that a **split** profile with the identical weighted mean does not.

### Added

- **Thermal acclimation of photosynthetic capacity** (#176), opt-in via
  `[leaf_physiology].thermal_acclimation`. Kattge & Knorr (2007): the peaked form's entropy term and
  the capacity ratio follow a running-mean growth temperature $T_g$ (°C),

  ```
  dS_v = 668.39 − 1.07·Tg     dS_j = 659.70 − 0.75·Tg     Jmax25/Vcmax25 = 2.59 − 0.035·Tg
  ```

  so a warm-grown and a cold-grown stand of the same PFT no longer share one temperature response.
  Measured across a 20 K range of growth temperature, the Vcmax optimum moves **28.2 °C → 35.0 °C**.

  **Measured in the coupled model, and decomposed** — because the headline number is misleading on
  its own. Ithaca, 5 years, growth temperature ≈ 10 °C:

  | | off | dS only | dS + ratio | dS | ratio | total |
  |---|---|---|---|---|---|---|
  | GPP | 0.004286 | 0.004823 | 0.006118 | +12.5 % | +26.8 % | **+42.8 %** |
  | NPP | 0.003675 | 0.004154 | 0.005332 | +13.0 % | +28.4 % | +45.1 % |
  | LAI | 0.048615 | 0.053505 | 0.062691 | +10.1 % | +17.2 % | +29.0 % |

  **Two thirds of the effect is the capacity ratio, not the optimum shift.** At Ithaca's ~10 °C
  growth temperature Kattge & Knorr put Jmax25/Vcmax25 at **2.24** against the PFT file's fixed
  **1.7** — a 32 % higher Jmax. The fixed value is a single global number being applied at a cold
  site; acclimation replaces it with a climate-dependent one. That is the same class of finding as
  #118's Vcmax review: a preset that was never climate-specific.

  Three deliberate choices, each recorded in `docs/science/leaf_gas_exchange.md`:

  - **The ratio acclimates too.** Shifting dS alone would acclimate the *shape* of each response
    while pinning the two branches in fixed proportion, which is not what the study measured.
  - **$T_g$ is the growth temperature, not the leaf temperature** — the fit is calibrated against the
    mean **air** temperature of the preceding weeks, so MEDS tracks an exponential running mean of
    the daily mean (window `acclim_window_days`, default 30 d).
  - **Site-level, not per cohort**, for the same reason. Per-cohort acclimation would need a
    relation calibrated on leaf rather than air temperature; it is on the ROADMAP, together with
    acclimation of *respiration*, which this does not touch.

  Applied by refreshing the per-PFT leaf table once per slow step from the **config** reference
  coefficients, which are never overwritten — so no kernel signature changes, the law lives in one
  place, and the refresh is idempotent. The running mean is prognostic and slow, so it is written to
  the state file; rebuilding a month of memory on every restart would make the optimum jump.
  `thermal_acclimation` **requires** `temp_response_form = "peaked"` — the Arrhenius form has no dS
  to shift, and the loader refuses rather than letting the flag be silently inert. Off by default,
  so the shipped path is unchanged.

- **All five phenology cues are wired; the drought-deciduous and light-exchanging strategies now run
  from configuration alone** (#150). `validate_config` used to reject the `WATER(2)`, `HYDRO(4)` and
  `LIGHT(16)` cue bits in either mask, because the kernel computed all five cues while the driver fed
  three of them zeros. Both halves are closed.

  **Four cue accumulators became per-cohort state.** `water_avg`, `low_psi_days`, `high_psi_days` and
  `light_avg` lived only inside `advance_leaf_phenology`, where `state = pheno_state_t()` re-zeroed
  them every slow step. A 10-day running mean reset daily is just its own instantaneous input, and a
  *consecutive*-dry-day counter reset daily never exceeds one — so those cues could not have worked
  even with their drivers present. They now carry the lockstep reorder, every creation site, and a
  **declared** fusion policy (survivor-keeps, stated explicitly rather than left implicit).

  **The drivers**, each a fast-loop daily reduction unless noted:

  | cue | driver |
  |---|---|
  | `CUE_WATER` | root-weighted fraction of extractable water, `Σ f_root,k · clamp01[(θ−θ_wp)/(θ_fc−θ_wp)]` |
  | `CUE_HYDRO` | the cohort's published predawn leaf potential (`dmax_psi_leaf`, already reduced for #95) vs a turgor-loss point **derived from the same pressure–volume curve** the leaf stress arrestor uses |
  | `CUE_LIGHT` | daily-mean incident shortwave at the canopy top (ED2 `rad_avg`) |
  | cold-drop | **top-layer soil temperature**, replacing the air-temperature proxy |

  The soil-temperature swap is a **number-mover for the already-shipped temperate-deciduous
  strategy**, not only an unlock: soil lags and damps air, so an air proxy crosses the 284.3 K and
  275.15 K cold-drop thresholds earlier in autumn than the soil does.

  **Nine `[phenology]` cue parameters became configurable**, through a new optional per-PFT array
  loader (`opt_pa`). Five had no table entry at all and four had a table entry but no loader, so the
  WATER/HYDRO/LIGHT cues could not have been tuned even once their drivers landed. They are optional
  — absent keys take the table defaults — because those cues are opt-in and a temperature-strategy
  config should not have to supply five numbers it never reads. A wrong-length array is still an
  error.

  **Acceptance, measured.** Both strategies run from a TOML config on the Ithaca driver, and
  reproduce the design's patterns 3 and 4 end to end: drought-deciduous holds a full canopy when
  watered, sheds under sustained drought (shed drive 0.996) and **reflushes on rewet** (0.996) — the
  reflush being the part that needs the persisted counters; light-exchanging sheds with light (0.9999
  bright vs 0.004 dim) while its flush stays permissive (1.000).

  **Selectable is not validated, and the release notes say so.** The kernel is unit-tested for all
  four strategies, the four drivers each have a hand-computed unit test, and both new strategies are
  asserted end to end — but **no MEDS run's leaf-area cycle has ever been scored against a phenology
  observation, under any strategy**, including the two that have shipped since v0.1.0. The thresholds
  are literature values for the biome each strategy describes, not site calibrations: selecting
  `CUE_LIGHT` with its default 200 W/m² onset at Ithaca strips the canopy every summer, because
  temperate summer insolation sits above a threshold chosen for a tropical dry season.

### Fixed

- **Soil-dimensioned output wrote its inactive padding as data** (#246). Every soil variable was
  emitted over all `n_soil_layer_max` (20) slots of what is normally a 10-layer column, because
  `layer_source_field` set `nlayer = n_soil_layer_max` against a comment that already said "the
  ACTIVE layer count". The `DIM_SOIL` branch marks every slot it is handed as valid, so the padding
  was written as data rather than left for netCDF to fill.

  Two consequences. Reducing over the soil axis gave nonsense — `soil_temp_site` averaged over its
  own soil dimension read **136 K**, ten real layers averaged with ten zeros. And
  `soil_matric_potential` was evaluated on padding whose `theta_sat` is 0: a divide by zero that a
  Debug build (`-fpe0`) **aborts** on and a Release build turns into NaN in the file. That abort was
  only reachable with `[output].water_fluxes = true`, which the shipped Ithaca config has off — so
  the one configuration that crashes was not one the suite or the test bed exercised.

  `nlayer` is now `dp%n_soil`, and `soil_z` is written over the same extent. The padding is
  therefore never written and comes back as `_FillValue`, masked by any CF-aware reader — which is
  exactly how the cohort and patch axes have always handled their own unused capacity. Measured
  after: `soil_temp_site` over the unmasked layers reads **273.8 K**, no unmasked non-finite value
  anywhere, and the Debug build completes with `soil_psi_site` enabled.

  `test_output_roundtrip` now runs a 6-layer column against the 20-slot ceiling and asserts the tail
  is fill for both the variable and the coordinate. Its source array is 290 K in *every* slot
  including the padding, so a writer that emitted the ceiling would read 290 and pass for the wrong
  reason; it has to read `_FillValue`.

- **Phenology was silently disabled in every run** (#245). **This is the largest behavioural change
  in v0.2.0: any PFT with a `[phenology]` block now actually has a leaf-area cycle.**

  `load_phenology_pft` gated the whole section on one key inside it —
  `toml_has(t, 'phenology.flush_cue_mask')` — and that is a key the shipped `meds_config_pft.toml`
  never documented. It documented `cue_mask`, a name no reader consumes. So a config that wrote a
  full, deliberate `[phenology]` block was **skipped in silence**: the presence map never saw
  twenty-three required keys go missing, and every PFT fell back to `CUE_NONE` / `CUE_NONE` —
  `flush = 1`, `shed = 0`, the evergreen fixed point — whatever leaf habit it declared.

  The Ithaca reference stand is declared cold-deciduous, has sane thresholds, and reaches 270.5 K
  soil temperature in winter. It held **LAI 5.28–5.66 through every January of a 50-year run**, with
  `flush = 1.0000, shed = 0.0000` in every month. **No MEDS run has ever had a leaf-area cycle.**

  Three changes:

  - **The gate is the section, not a key inside it.** New `toml_has_section` asks whether the author
    meant to configure phenology at all, rather than whether they spelled one particular key the way
    the loader wanted. A `[phenology]` block missing keys is now a hard error naming all of them.
  - **`cue_mask` is rejected by name**, with a migration message giving both replacements and the cue
    bits. It was retired when the flush and shed cues became independently selectable; ignoring it
    silently changed a PFT's leaf habit, which is the one outcome worse than stopping.
  - **The shipped `meds_config_pft.toml` block is corrected and completed** (all twenty-three keys,
    and a note that `evergreen = [0]` alone does not give a PFT a season), as is
    `examples/example_biophysics/pft_parameters.toml`, which used the stale name live along with two
    more keys — `on_threshold` / `off_threshold` — left over from the retired status/deadband design.

  With the real key names, the Ithaca stand runs as it always should have: **LAI 0.001 (Feb) → 4.72
  (Sep) → 0.70 (Dec)**, shed governor 0.97 in January and 0.00 through summer. The phenology kernel
  was never broken — it was never switched on.

  This is the **second** occurrence of the absent-key trap that `.claude/rules/config.md` already
  describes, in a different block. The rule text described it exactly.

  **Migration:** a config carrying `cue_mask` now stops with instructions. A config with no
  `[phenology]` section is unchanged (evergreen defaults), which is the documented fallback. A config
  that had a complete block under the real key names was already working and is unaffected.

- **An out-of-bounds write in the output path: the diagnostic blocks kept a stale count after a
  cull** (#247). `cohort_diag_reorder` and `patch_diag_reorder` permuted the diagnostic rows but
  never updated the block's own `n`, while `cohort_reorder` set `cohort%n = m` right after calling
  them. So after any operator that **shrinks** the array, the block claimed more live slots than the
  array had — and `extract_variable` sizes its scratch buffer from the live count while
  `cohort_diag_value` filled `1..d%n` into it. A write past the end: silent in Release, `subscript
  97 > upper bound 96` under Debug, a SIGSEGV some steps later.

  Invisible until now because the reference stand only ever **grows** — Ithaca goes 114 → 125
  cohorts over three years, so the block's count was never the larger one. It becomes routine the
  moment phenology works (#245): shedding drives cohorts below the tracking floor,
  `terminate_cohorts` culls them, and a two-year run dies in its first December.

  Fixed at both ends. The count now rides the lockstep like every other field, and
  `cohort_diag_value` / `patch_diag_value` bound their write by `size(x)` as well as by the block —
  the second is what turns the next count mismatch into a short read instead of memory corruption.
  `test_containers` asserts both after its cull, and fails on the unfixed code.

- **Slow per-patch diagnostics emitted a per-SECOND rate under a per-YEAR label** (#239). **Six
  shipped variables read a factor `yr_sec = 3.1557e7` too small**: `litter_leaf_site`,
  `litter_leaf_patch`, `litter_fineroot_site`, `litter_struct_site`, `nplant_recruit_site` and
  `disturb_area_site`. Anyone who plotted litterfall or recruitment from a MEDS run before this
  should re-read those series.

  The patch diagnostic block is `(value, weight)` and the reader returns `value/weight`. The weight
  is `Σ dt` in **seconds** — `accumulate_patch_diag` adds `dt_fast` every fast step — so every row
  must contribute *(its rate, in the units the registry declares)* × *(dt in seconds)*. The fast rows
  always did, with `flux * dt_fast`. The five slow rows, added later, contributed a bare per-step
  **amount** (and, for recruitment, `rate * dt_years`) — dimensionally a per-second rate.

  On a spun-up Ithaca stand (114 cohorts, AGB 16.7 kgC/m², LAI 5.56), annual means before → after:
  `litter_leaf_site` **9.80e-09 → 0.309 kgC/m²/yr**, `litter_fineroot_site` **1.21e-08 → 0.380**,
  `litter_struct_site` **2.07e-08 → 0.654**, `nplant_recruit_site` **1.11e-09 → 0.0349 plant/m²/yr**.
  The corrected leaf number is the independent check: LAI 5.56 at SLA ≈ 20 m²/kgC is ≈ 0.28 kgC/m² of
  leaf carbon turned over annually.

  **Why nothing caught it.** A diagnostic is downstream of every conservation ledger, and the slow
  ledger's own litter term reads `patch%litter_in` directly rather than this slot, so the two never
  had to agree. The one test that did touch a slow slot — `test_disturbance`, added with #170 —
  asserted the **raw accumulator** rather than the number `patch_diag_value` emits, so it certified
  the contribution and was blind to the convention it was written against. It now asserts through the
  reader, and a new `slow_diag_units` test asserts, for all five slow slots, the identity the
  per-year label *means*: `reported rate × elapsed years == the amount that actually flowed`, against
  oracles outside the diagnostic path (`patch%litter_in`, the recruit pool, the configured hazard).

- **The phenology memory now survives a restart** (#150). No phenology state was written to the state
  file at all — not the two governors, not the GDD and chilling sums. A restart resurrected every
  cohort at its **birth** values (`flush_drive = 1`, `shed_drive = 0`, `gdd = chill = 0`, the
  evergreen fixed point), so a temperate-deciduous stand restarted in January came back with flushing
  permitted and no chilling accumulated — it leafed out in midwinter and then had to rebuild weeks of
  thermal memory. This affected the two strategies that have shipped **since v0.1.0**, not only the
  two #150 adds. None of it is derivable from the instantaneous state: these are time integrals over
  the preceding weeks. All eight columns are now written, and optional on read, so an older state
  file still restarts exactly as it did.

### Examples

- **Every example regenerated against v0.2.0**, and two of them were broken.
  `example_leaf_gas_exchange` raised `TypeError` on every run: it passed `theta=p.theta_j` to
  `assimilation_demand_c3`, which #118 renamed when it split the C3 co-limitation curvatures out —
  literally the confusion #118 was about, sitting in the example. `example_demography` is the entry
  below. `example_phenology` now genuinely exercises all four strategies; two could not be selected
  before #150.

  **`example_biophysics` is re-spun and genuinely cold-deciduous for the first time** (#245). It
  ends at 20 cohorts / 3 patches, peak LAI 4.06, AGB 9.00 kgC/m², mean dbh 23.5 cm, 15.4 kgC/m² soil
  carbon — against 128 / 12 / 5.32 / 15.75 / 35.2 / 24.2 when it ran evergreen. Establishment is
  about twice as slow: a deciduous seedling rebuilds its canopy from storage every spring and spends
  its early growing seasons near break-even, so little happens until year 15. Its trajectory figure
  now plots the **annual maximum** LAI — it sampled at the year boundary, which is 1 January, and for
  a deciduous stand that plotted a bare canopy topping out at 0.28 for a stand whose July LAI is 4.1.

  **Its July stage drops from `dt_fast = 150 s` to the 900 s production default.** The short step was
  justified on the argument that a diel figure needs sub-daily fidelity and a long step smears the
  energy partitioning. Measured over that exact July, it does not: canopy air 22.69 vs 22.67 °C mean
  and 9.34 vs 9.37 K diel amplitude, leaf 23.32 vs 23.31 and 12.92 vs 12.91, soil surface 23.42 vs
  23.39, leaf−air +4.04 K by day either way, `dmax_psi_leaf` −0.2951 vs −0.2954 MPa. Every number the
  figures report agrees to 0.03 K, for 14.3 s of wall time against 4.0 s. `dt_fast` is an **accuracy**
  parameter and not the output cadence — the hourly records come from `fast_interval_steps` — and the
  config and README now say so, along with what 900 s *is* wrong for (the sub-daily `psi_leaf`
  excursion, #162, for which MEDS emits no diagnostic).

- **The demography example is slow-scale demography only** (#260), which is what it always claimed
  to be: cohort and patch dynamics, fusion and fission, growth, mortality, recruitment and treefall.
  **No carbon dynamics and no soil carbon.**

  It had stopped being that without anyone noticing. Its vital-rate laws are LAI-driven and
  empirical, and the reorg moved them out of Fortran into `empirical_laws.py`; the Fortran model has
  only the *carbon* path. So the example's config, brought up to schema in September and never
  re-run, was silently pointing `meds_main` at **a different model** — stub GPP
  (`gross_gpp = gpp_ref * leaf_area`) with no light competition, hence no negative feedback on leaf
  area. It diverged to LAI 501 and AGB 907 kgC/m² and died in `cohort_reorder`. The committed output
  beside it was older still, from the deleted *Fortran* empirical model.

  The example's driver is now the Python one that actually implements its laws, and it writes its own
  output: new `meds_site_get_patch_real` / `meds_site_get_patch_int` in the demography C-API expose
  the patch areas, ages and the cohort→patch CSR map, and `_write_nc.py` writes the ragged
  cohort/patch netCDF `post_proc/` already reads.

  **The stand now equilibrates and is physically sensible**: 250 years settles at ~0.97 stems m⁻²,
  **AGB 16.5 kgC m⁻², LAI 7.4**, 266–355 cohorts over 12 patches, with a textbook inverse-J size
  distribution — 0.43 stems m⁻² below 1 cm DBH down to 0.0008 above 50 cm, and **76% of the biomass
  in stems over 20 cm**. The previous committed output had AGB 121 kgC m⁻², about four times the
  densest forest on Earth. The run takes **22 seconds**, against 2 min 18 s for the carbon run that
  crashed, and `empirical_spinup.py` still reproduces its golden exactly.

  The stale `example_output_pft_parameters.csv` is deleted: it described the deleted Fortran model
  (`growth_lai_slope`, `mort_gamma/alpha/beta`, no carbon traits), and a provenance record for a
  model that does not exist is worse than none.

### Documentation

- **The leaf water-stress divergence from ED2 is now recorded as a decision, not an open question**
  (#47), in `docs/science/leaf_gas_exchange.md` §4.3. MEDS keeps the Sabot et al. (2022) two-limb
  scheme; ED2's `farq_katul.f90` uses a Manzoni-style single `stoma_beta` plus a turgor-loss
  capacity term. The three divergences and the reason each one stands: the capacity limb's shape
  (linear $\psi_{leaf}$ ramp vs ED2's 6th-power turgor-loss) and its targets ($V_{cmax}$, $J_{max}$
  and TPU vs ED2's $J_{max}$ and $\alpha$ only, which acts only on the light-limited branch); the
  stomatal limb's scope (Leuning/Medlyn $g_1$ **and** Katul $\lambda$ vs Katul only, which would
  leave the two explicit schemes with no drought response at all); and the parameter split, where
  MEDS's $s_{ref}\,e$ is exactly ED2's `stoma_beta`, so a calibrated value transfers as
  `stoma_beta = -sref_stomata * lambda_psi_exp`.

  **The `psi_soil` wiring defect the issue also filed was already closed by #95** and needed no work
  here: the field is `psi`, the driver passes `dmax_psi_leaf` on both call paths, and unset cohorts
  are seeded from the surface-layer soil potential. Verified rather than assumed.

### Fixed

- **Saturation over a frozen surface now uses the ice curve** (#89). `sat_vapor_pressure`,
  `sat_specific_humidity` and both their temperature derivatives take an optional `fliq` (liquid
  fraction). Absent means pure liquid, which is exactly what every caller did before, so omitting it
  is **bit-identical**. Supplied, they blend:

  ```
  e_sat = fliq * 611.2 exp(17.67 Tc/(Tc+243.5))  +  (1-fliq) * 611.2 exp(21.87 Tc/(Tc+265.5))
  ```

  Both branches carry the same 611.2 Pa constant, so they cross **exactly** at `Tc = 0` and the
  blend is continuous in temperature *and* in `fliq` — a pack that freezes or melts slides between
  the curves instead of stepping, which matters because a step here lands in the right-hand side an
  adaptive controller integrates.

  **What was wrong.** `snow_surface_fluxes` drove sublimation with `sat_specific_humidity` on the
  liquid curve, while the enthalpy side of the same routine was already ice-aware — removing
  `enthalpy_vapor` from an ice-referenced layer debits sublimation (vaporization + fusion)
  automatically. So the model treated the snow surface as ice for *energy* and as liquid for *vapour
  pressure*. Over ice `e_sat` is **10 % lower at −10 °C, 22 % at −20 °C and 34 % at −30 °C**, so the
  driving gradient was overstated by those factors. The same applied to ground evaporation from
  frozen soil, which now uses the top layer's `soil_fliq`.

  **Measured** (Ithaca, 5 years, `dt_fast = 900 s`), over the 18 months with snow on the ground:

  | | liquid curve | ice branch | change |
  |---|---|---|---|
  | latent heat flux | 0.96083 | 0.91988 | **−4.26 %** |
  | snow water equivalent | 4.85938 | 4.87776 | **+0.38 %** |
  | sensible heat flux | −2.69802 | −2.69890 | −0.03 % |
  | soil surface temperature | 276.595 | 276.596 | +0.000 % |

  Peak SWE rises 0.18–0.26 % in each of the five winters; the February shoulder gains 2.1 %. The
  flux responds by less than `e_sat` does because what drives it is the *gradient*
  `q_sat(T_s) − q_CAS`, and the winter canopy air is often near saturation — and because the Ithaca
  pack spends much of its life near 0 °C, where the two curves nearly coincide. At a cold
  continental site holding −20 °C for months the correction is a much larger fraction.

  **Not applied to dewpoint conversion or diagnostic VPD.** Dewpoint is *defined* over liquid, so an
  ice branch there would mis-convert the forcing. Also not applied to the canopy, which has no ice
  state at all — intercepted snowfall is held as liquid with no fusion debit, which is a separate
  known gap.

  The ED2 proposal this issue points at (EDmodel/ED2#442) concerns the liquid formula. Measured
  against Murphy & Koop (2005), MEDS's existing Bolton form is within **0.2 %** from −30 to +40 °C,
  so the liquid curve was never the problem and is left alone. The ice form added here is within
  0.1 % at −10 °C and 0.9 % at −30 °C.

### Changed

- **One solar declination for the whole model** (#152). `solar_cosz` used Cooper (1969),
  `daylength` used White (1997); both now call `solar_declination(doy)`, which is the Cooper form.
  The two were the same function offset by 1.25 days — $-\cos x = \sin(x-\pi/2)$ puts White's
  ascending zero crossing at doy 82.25 against Cooper's 81 — and Cooper's is the closer to the true
  vernal equinox near doy 79–80. It is also the form the radiation path already used every
  `dt_fast`, so `daylength` moved rather than `solar_cosz`.

  At Ithaca (42.44 °N) daylength changes by at most **3.7 minutes**, at the equinoxes, and by
  essentially nothing at the solstices; the 10.5 h autumn phenology cue fires **2 days earlier**
  (doy 297 → 295). `solar_cosz` and therefore the radiation are unchanged.

  `test_time` now asserts the shared source by **inverting each consumer back to a declination**
  rather than re-implementing the formula, so a future divergence is caught even if the formula
  itself changes. Verified to fail on the pre-fix code.

### Added

- **C3 gets its own co-limitation curvatures** (#118): `theta_cj_c3` (default 0.98, the
  $A_c$/$A_j$ transition) and `theta_ip_c3` (0.95, the transition with $A_p$). Both are new
  **required** keys in the `[pft]` table. C3 previously passed `theta_j` — the curvature of the
  electron-transport hyperbola — into both of its smoothings, while C4 already had a dedicated pair.
  A co-limitation curvature and a light-saturation curvature share units and a functional form and
  nothing else; CLM and FATES put the C3 $A_c$/$A_j$ curvature near 0.98 against the 0.7–0.9 that
  fits the $J$ hyperbola.

  **At the leaf**, PFT-1 kinetics, $C_i = 280$ µmol mol⁻¹, saturating light
  ($A_c = 14.40$, $A_j = 17.53$, $A_p = 45.0$, so $\min = 14.40$):

  | curvatures | after $A_c$/$A_j$ | after $A_p$ | total vs min() |
  |---|---|---|---|
  | `theta_j` = 0.85 (before) | 11.31 (−21.4 %) | 10.80 (−4.5 %) | **−25.0 %** |
  | 0.98 / 0.95 (now) | 13.49 (−6.3 %) | 13.22 (−2.0 %) | **−8.2 %** |

  i.e. **+22.4 % gross assimilation** at that point. The shortfall was nearly independent of
  $V_{cmax}$ (30, 31, 31, 32 % at $V_{cmax,25}$ = 60, 90, 120, 150), so it was a systematic offset,
  not a regime effect — no measured $V_{cmax}$ reproduced a measured rate.

  **In the coupled model**, Ithaca, one year from a common spun-up stand with **demography frozen**,
  so stand structure cannot diverge and the difference is physiology:

  | | `theta_j` | `theta_*_c3` | change |
  |---|---|---|---|
  | GPP | 0.07633 | 0.09039 | **+18.4 %** |
  | NPP | 0.06496 | 0.07869 | **+21.1 %** |
  | NEE | −1.4888 | −1.8274 | −22.8 % |
  | LAI (allocation still runs) | 0.9796 | 1.0128 | +3.4 % |

  With demography **live**, 10 years from cold start, the same change reads GPP +66 % and AGB +89 %
  — that is the *compounded* figure, because a persistent rate increase accelerates stand
  development, and it is not an equilibrium sensitivity. The controlled number above is the one to
  quote.

  **The $V_{cmax,25}$ presets were re-examined and deliberately left alone.** The concern was that
  they might have been tuned against the over-smoothing, in which case the compensation would have
  been hidden inside a parameter named for a different process. They were not: `vcmax25 = [60, 45,
  40]` entered in the commit that first added the leaf module and has never been revised, the values
  sit mid-range for their PFT descriptions, and the only other GPP-facing knob (`gpp_ref`) is the
  stub used when the fast loop is off. So this is a correction, not the unwinding of a calibration,
  and the presets should not be lowered to absorb it.

  `leaf_photo_params_t%theta_cj`/`theta_ic` are renamed `theta_cj_c4`/`theta_ic_c4` so the pathway is
  visible at the point of use — which is the whole defect — and the C API mirror, `meds/plant/_ffi.py`
  and `LeafParams` follow. `meds_assimilation_demand_c3` now takes two curvatures instead of one.

- **A Dirichlet bottom thermal boundary for the soil column** (#145), selected by
  `[energy].bottom_bc = "dirichlet"` with `deep_temp` [K] and `deep_depth` [m]. The bottom node
  conducts to a plane held at `deep_temp`, a distance `deep_depth - |z_node(n)|` below it. The
  default stays `geothermal`, and the Neumann path is bit-identical to before.

  **The premise, re-measured against an exact oracle.** A homogeneous column at constant `theta` has
  constant `kappa` and `C`, so the analytic semi-infinite solution `exp(-z/d)` with
  `d = sqrt(2*alpha/omega)` holds exactly — an oracle independent of the model's own machinery. The
  harness is validated by a 12 m column with the adiabatic base, which reproduces that profile to
  0.1 % over the top three damping depths. At the default 2 m column, amplitude relative to the
  surface at the bottom node (−1.73 m):

  | bottom BC | amplitude ratio | vs analytic 0.421 | RMS error over the profile |
  |---|---|---|---|
  | `geothermal` (adiabatic) | 0.764 | **+82 %** | 0.145 |
  | `dirichlet`, `deep_depth = 3.12` | 0.412 | **−2 %** | 0.015 |

  **The anchor depth is derived, not fitted.** A resistive termination reflects least when its
  impedance `kappa/l` matches the magnitude of the half-space impedance `kappa*(1+i)/d`, i.e. when
  `l = d/sqrt(2)`. For the default column that puts the anchor at `1.727 + 1.973/sqrt(2) = 3.12` m; a
  measured sweep puts the optimum at 3.0–3.1 m, so the derivation is right to within one sweep step.

  **In the coupled model** (Ithaca, 10 years from cold start, `dt_fast = 900 s`, `deep_temp = 283.24 K`
  = the forcing's mean annual air temperature; final-year monthly means):

  | | `geothermal` | `dirichlet` | change |
  |---|---|---|---|
  | annual swing at −1.73 m | 26.25 K | 15.39 K | **−41 %** |
  | annual mean at −1.73 m | 290.47 K | 286.41 K | **−4.06 K** |
  | Rh, annual mean | 0.00240 | 0.00230 | **−4.2 %** |
  | GPP / NPP / AGB / LAI | — | — | +0.6 % / +0.8 % / +1.0 % / +0.7 % |

  The adiabatic profile *flattens* below −0.45 m (swing 30.0, 28.6, 27.9, 27.1, 26.3 K) — the
  reflection signature — while the anchored one keeps decaying (29.0, 26.7, 23.7, 20.2, 15.4 K). The
  Rh change is seasonal in shape, not a level shift: the August peak falls 0.0057 → 0.0051 while
  April rises 0.0016 → 0.0018, which is the Jensen argument in #145 made visible. The soil-carbon
  difference (+3.6 %) is a 10-year transient-accumulation difference, not an equilibrium one — the
  pool is still climbing at the end of the run.

  **What it does not fix, stated plainly.** A purely resistive termination cannot reflect less than
  0.41 in amplitude at any `l`: matching a complex impedance with a real one leaves the phase wrong by
  45°. Closing the rest needs heat *capacity* below the column — passive deep thermal layers under the
  hydrologically active one — which is on `docs/ROADMAP.md` as the #145 follow-up. The anchor buys
  about a factor of ten in this metric, not exactness.

  `deep_temp` is **required** when the Dirichlet BC is selected. An error in it is a steady flux
  `kappa/l * error` into the column base, so a silent default would reintroduce a mean-annual
  deep-soil bias — the very thing this boundary condition removes.

### Fixed

- **The PFT-parameter CSV dump lost its last two columns and printed a bit pattern** (#118). The
  `write_pft_params_csv` format had 44 item slots against a 46-item output list. A Fortran format
  shorter than its list does not fail — it reverts to the last repeat group and keeps going — so an
  integer edit descriptor silently received a real (`fineroot_turnover_rate` printed as
  `4605380978949069210`) and `f_labile_stem` / `struct_lignin_frac` vanished. Introduced by adding
  the two curvature columns in this same release and caught by a new column-count assertion in
  `test_pft_optics_config`, which is verified to fail on the broken format.

- **The forcing file and the config are now checked against each other** (#185). Five items, decided
  individually:
  - **`dt_forcing` is validated against the file's actual record spacing**, and the axis against
    itself for uniformity. This was the real hazard and it is *not* the one the issue named: the
    value was read straight from the config and **never compared with the file at all**, while
    placing interval midpoints, disaggregating shortwave and bracketing the recycle seam. A config
    saying 3600 s against a half-hourly file mis-timed all three, plausibly. Checking the spacing is
    strictly stronger than checking the `timestep_seconds` attribute, which stays provenance.
  - **`avg_convention` and `sw_input_kind` are read and validated** when present, skipped when
    absent so older files still load. `sw_input_kind = "total"` against `sw_partition =
    "passthrough"` used to crash on a missing variable; it is now a clear rejection.
  - **`avg_convention = "instant"` and `"center"` are rejected.** Both parsed and then ran the
    end-of-interval path anyway, because only `METAVG_BEGIN` has a branch — so selecting either got
    a scheme the user did not ask for, silently. Same precedent as `lwdown_source`.
  - **`SWPART_SIB` deleted.** A reserved code with no implementation and no TOML spelling: it could
    never be selected, and `partition_shortwave`'s default branch would have routed it to Erbs.
  - **`elevation(grid)` stays unread**, deliberately: site elevation is a `[site]` property and the
    variable is provenance about the source grid.

  Reported through `met_open`'s existing status channel rather than `error stop`, which is what
  makes each rejection testable — the module's own comment says that is what the channel is for.
  The synthetic test fixture now writes the two global attributes the prep script writes; without
  that the new checks would have been skipped in the test while firing in production.
- **The tissue-water floor now reports the water it creates** (#148). `advance_water_mass_full`
  floors `leaf_water_mass` and `wood_water_mass` at a tiny positive value so the linear mass Euler
  step cannot go negative — and creates water doing so. The code's own comment called the case
  "unobserved in this pass's test scenarios", which was **a belief, not a measurement**: the floor
  fires per cohort per tissue, while the whole-column water ledger sums leaf + wood over all
  cohorts, so water created in one cohort's wood is indistinguishable from a redistribution between
  cohorts. The budget closed to ~4×10⁻¹² kg m⁻² with the floor entirely unmonitored. It is now
  reported through the existing `budget%clamp_mass` / `clamp_commit_n` commit-clamp channel, which
  already reduces to `site%work_clamp_mass` and already has an output variable — so it surfaces
  end-to-end with no new reporting surface, which is what the issue asked for. Measured on a forced
  fixture: 2 activations creating 2.92×10⁻² kg m⁻², and **zero on an ordinary step**. The claim in
  `ark2_column_step` that its commit counter "stays 0 by construction" was corrected — that claim
  was the whole of the defect.
- **`PD_DISTURB_AREA` is written and emitted** (#170). The slot was declared in the patch
  diagnostic block with no writer and no registry row, so the disturbed-area flux read as a **silent
  zero** rather than a missing variable — the harder failure to notice, and one no conservation
  check can see, because zero disturbed area is a perfectly conservative answer. Written on the
  donor patches inside `apply_patch_disturbance` *before* the gap is appended, which is where the
  diag slots still line up with the donors; writing after the append is the out-of-bounds trap this
  file has already paid for once. New output variable `disturb_area_site`.
- **A run selecting `time_integrator = "rk45"` above `dt_fast` = 300 s is now warned** (#160). The
  transpiration corrector that cut a ~1 MPa `psi_leaf` error by 314× (PR #91) lives in
  `advance_water_mass_full`, which ARK calls and RK45 does not, so an RK45 production run silently
  carried an error the default path does not. A warning rather than an error, because RK45 is the
  accuracy baseline and is meant for a fine step.
- **Γ\* now responds to `o2_mol_frac`, so the O₂ knob propagates to both places oxygen enters the
  C3 demand** (#117). The compensation point is set by Rubisco's CO₂/O₂ specificity and is
  proportional to the O₂ partial pressure, but it was computed with no O₂ dependence at all — so
  raising O₂ correctly inhibited carboxylation through `Kc(1 + O/Ko)` while leaving the entire
  photorespiratory penalty on `Aj` untouched. Scaled by `o2_mol_frac / 0.209`, the O₂ the shipped
  `gstar25` was measured at (Bernacchi et al. 2001). **The factor is exactly 1 at the default, so
  every shipped configuration is bit-identical** (10.18931 µmol m⁻² s⁻¹ before and after on the test
  fixture). Away from it, measured on a midday tropical leaf: halving O₂ raises A_net by 22.0 %
  against 8.2 % before, and 35 % O₂ cuts it by 23.7 % against 10.5 %.
- **Two roadmap items described code that no longer exists**, found by re-measuring every filed
  premise before scheduling it (#168, #166). `bsap` stopped being a placeholder in PR #125 —
  `set_cohort_wood_geometry` derives it from ED2's real `b1SA`/`b2SA` sapwood-area allometry. The
  old `0.10 * wood_carbon` placeholder made the wood thermal time constant **6.5–10× too short**
  across the whole size range (`f_sap` runs 1.00 at dbh ≤ 19 cm to 0.655 at 117 cm, against 0.10).
  `veg_energy_step_implicit` was deleted in PR #120, and `veg_energy_diagnostic` does not exist
  either; `veg_energy_balance` is the single closure. Both entries corrected in `docs/ROADMAP.md`
  and their four stale source/doc references removed.
- **`scripts/numerics_sweep.py` could not run a cross-scheme comparison.** The `--parity` preset
  pinned three config keys that no longer exist (`fast.integration_scheme`,
  `fast.leaf_energy_model`, `fast.wood_energy_model` — #199 named two), and the `SCHEMES` table
  still offered `split` and `picard`, which are now a hard error. The preset is removed rather than
  repointed: every difference it pinned has since been closed by making the schemes agree.
- **The C-API demography shim inlined the allometry** instead of calling `meds_allometry`, so a
  coefficient change updated the model and not the shim (#200). It now calls `dbh_to_height`,
  `dbh_to_agb`, `dbh_to_leaf_area`, `size2leaf_carbon` and `size2wood_carbon`.
  `examples/example_demography/empirical_laws.py` reads the four allometry coefficients from the
  `[allometry]` block of its shipped PFT config instead of hard-coding them; the values are
  unchanged, so the example's behaviour is unchanged.

### Documentation

- **`fuse_cohort_fast_state` says what it is: the one declaration of the per-cohort fast-state
  policy, centralised but not enforced** (#190, deferred). Four of the five benefits that justified
  deleting `column_cohort_t` have already landed piecemeal — `column_cohort_init` gives the test
  fixtures allometric consistency (the `bwood`-on-uninitialized-memory bug is fixed), the three
  hard-coded canopy constants are PFT parameters, the derived geometry is on the cohort block, and
  `reconcile_tissue_water_capacity` took the seed and clamp out of the fast gather. The fusion and
  scaling policy is centralised too, in `fuse_cohort_fast_state` and `scale_cohort_ground_fields`.
  What is left is **completeness you cannot forget** — a table the blend *iterates* cannot omit a
  field that a hand-written routine can — and that is #146's silent-omission class. The two are now
  **paired for v0.3.0**: one packed, policy-carrying layout should serve the fast state vector and
  the cohort slice together, because separately each is a large refactor buying a fraction of one
  property.
- **The RK45 stiff rescue keeps its whole-step rollback** (#161), measured rather than optimised.
  E5 proposed snapshotting at the point of failure to avoid redoing accepted sub-steps. Over a full
  simulated year at Ithaca on `rk45` at `dt_fast` = 900 s the rescue fires **zero times**
  (`work_rk45_rescue_site` = 0, against 70 577 integrator sub-steps on the same run — the counter is
  live and the zero is real), so there is no work to save on the reference workload. It is also not
  free to build: retaining part of RK45's boundary-flux accumulation while ARK finishes the interval
  means one `dt_fast`'s ledger summing two schemes' contributions, and it is well-defined only for
  the `stiff_bail` trigger — `rk45_state_railed` tests the *final* state, so on that path no
  sub-step is identified to resume from. Recorded at the site.
- **The FAST output tier's staging path is not a duplicate switchboard, and is kept** (#172). It was
  slated for deletion on the belief that it duplicated the general extraction path. It does not:
  the tier is *already* on the general machinery — it shares the registry, the buffers, `close_tier`
  and the serializer — and the only bespoke part is the extraction *source*, which is forced by when
  the values exist. Sub-daily quantities are sampled **during** the fast loop; by the time `main`
  folds them, `site` holds the post-fast-loop, end-of-slow-step snapshot, so resolving them against
  live site state would silently emit end-of-day values on a sub-daily axis. The two functions
  cannot overlap either: `SRC_S_*` occupy 4000–4999 and `SRC_F_*` 5000–5999, disjoint by
  construction. Deleting the staging would have deleted sub-daily sampling. Recorded at the site.
- **The frozen-seam contract is written down** (#201):
  [`docs/dev_plans/MEDS_FROZEN_SEAM_CONTRACT.md`](docs/dev_plans/MEDS_FROZEN_SEAM_CONTRACT.md). The
  Λ = F·dt/S criterion and its case split (linear-in-store is scale-free and sound; a prescribed
  flux against a prognostic store is unsound, with the +30.7 µmol m⁻² s⁻¹ NEE scar to prove it), the
  four seams classified against it, why Λ is meaningless for a *rate* seam and what replaces it
  (debit-before-credit), and the arbitration rule — scale all demands by `min(1, S/D)` — for a store
  with several consumers, which is order-independent where per-process clamping is not.
- **The ED2 two-stream defects found during the port are folded into
  [`docs/ed2_comparison.md`](docs/ed2_comparison.md) §5a** (#6), with the MEDS ↔ ED2 RT structure
  mapping. All six are reported upstream. Two of them — the stale-PFT diffuse index and the missing
  clumping factor in the longwave split — are *impossible by construction* in MEDS, and that is why
  its RT is shaped the way it is.
- **`test/` stays flat, deliberately** (#193, structure-plan decision #13), recorded in
  `src/README.md` with the reason: the discipline that matters is the link line, not the directory.
- **The year-rollover `rh_seam_gap` residual is attributed** (#192). 8.370×10⁻⁴ kgC m⁻² at a year
  boundary is the annual patch cadence, not a leak: the fast loop accumulates against one patch
  composition and the daily step debits another, blended through a matrix nonlinear in the lignin
  fraction. Documented at the field, with the instruction not to widen a tolerance to absorb it.
- **The build files no longer describe GPU offload as the parallel path** (#194). Measured, the
  offload build runs **1.4× slower** than the CPU (49.8 s against 36.2 s), one kernel sits at 0.4 %
  occupancy, and the device treated as 20 slow cores is 26× slower than 4 CPU cores. Patch-axis CPU
  threading is the parallel path. `MEDS_GPU=gpu` is kept as a reproducible experiment.
  `CMakeLists.txt`, `docs/building.md`, `src/README.md` and `docs/ed2_comparison.md` corrected.
- **A v0.2.0 release plan** in [`docs/dev_plans/archive/MEDS_V02_RELEASE_PLAN.md`](docs/dev_plans/archive/MEDS_V02_RELEASE_PLAN.md):
  all 66 open issues triaged into six phases plus release, 48 in scope and 18 deferred to v0.3+,
  with twelve decisions recorded.
- **The documentation was reorganized against the restructured source tree.** Thirty design
  plans moved to `docs/dev_plans/archive/` with tombstones; eleven stay live or as reference.
  The README dropped from 310 lines to a reader's entry point, with building, configuration
  and post-processing moved to their own pages. A new [`src/README.md`](src/README.md)
  documents the source layout, the placement rules and the library graph.
  `CLAUDE.md` split into a short always-loaded file plus path-scoped rules. This file and
  `docs/ROADMAP.md` were created. Three missing science pages were written:
  `docs/science/soil_carbon.md`, `forcing.md`, `plant_respiration.md`.

### Changed

- **The soil-energy adaptive-substep surface is deleted** (#163), after measuring rather than
  assuming. `soil_energy_step_implicit` hard-codes `flux%nsub = 1` and backward Euler is
  unconditionally stable, so an inner substepper could only ever buy *accuracy* — and the accuracy
  of the integrated state is already owned by the outer adaptive march through `GRP_SE` (soil
  internal energy, which is what the state vector carries). A second controller on a quantity the
  first one already controls is not worth wiring.
  The dead surface was wider than filed: besides `substep`, `h_init` and `max_substep`, `[energy].rtol`
  fed only `GRP_SOIL_T` — **a tolerance group with no member in `state_wrms_grouped`**, which sums
  `GRP_SE` and `GRP_THETA` and never `GRP_SOIL_T`. The group, its `ATOL_SOIL_T_DEF` default and the
  `ENERGY_SOLVER_BE` / `ENERGY_SUBSTEP_ADAPTIVE` single-valued enums went too; `N_TOL_GROUP` is 7.
  `[energy].atol` and `debug_error` stay — both are live, `atol` as the budget-closure threshold for
  the debug halt, not a step tolerance. `bottom_bc` stays because #145 wires the Dirichlet thermal
  anchor onto it in this same release.
  One coupling was preserved deliberately: `atol_scale` used to reach `[energy].atol` by
  round-tripping config → `tols(GRP_SOIL_T)` → config, so deleting the dead group would have
  silently dropped it. It is now applied directly, and a test asserts it.
  **Verified byte-identical**: all 14 output files of a 1-year Ithaca ARK run — restart state and
  every monthly diagnostic — are unchanged.
- **The two unreachable heterotrophic-respiration kernels are deleted** (#153).
  `heterotrophic_respiration_damm` (Davidson 2012) and `heterotrophic_respiration_flux` (Q10 / ED2
  capped exponential) were tested but **could not be selected**: there was no `hr_model` TOML key
  anywhere and `co2_opts_t` was never carried by `meds_config`, so the selector could not be set
  and the kernels could not be called. A tested-but-unreachable kernel is the worst of both worlds —
  maintenance and review weight for nothing, while reading to a newcomer as an available option.
  There is **one** production Rh authority: the CENTURY matrix.
  The dead surface was wider than the issue recorded: the `HR_*` selector codes, `co2_opts_t`,
  `damm_params_t`, the shared `water_modifier` helper, three DAMM-only constants in
  `meds_constants`, and unused `co2_opts_t` imports in `meds_fast_types` and `meds_fast_prepass`.
  Net **−237 lines**. The implementations are preserved on branch `archive/damm-hr` (at `2fcb647`).
  Five subtests in `test_column_co2` went with the kernels they exercised; in
  `test_soil_biogeochem` only part (a) of the fast/slow seam test went — part (b), the diurnal
  accumulation and the Jensen counter-check, builds its own factors inline and is untouched.
- **The IMEX-Euler oracle tier is retired** (#198). It was not an oracle: it returned no reference
  trajectory, and `imex_euler_column_step` was a two-line wrapper around `column_be_stage` plus
  `advance_water_mass_full` — the production scheme's own kernels — so its independence was in the
  tableau, not in the machinery it was supposed to check. All five of its consumers kept their
  coverage: four now call a test-local `be_euler_step` (the same two-line composition, living in the
  code that degrades the scheme), and `test_adaptive_march` was **ported to `adaptive_ark_march`**,
  so it now exercises the production controller instead of one that existed only to serve the tier
  (8 sub-steps at `rtol` 1e-3 against 25 at 1e-6). `meds_fast_rk4_oracle` holds one oracle and its
  name is accurate again.
- **The RK4 oracle's independence is now structural.** With the tier gone, the module no longer
  imports `meds_fast_be_stage` at all — no `column_be_stage`, no `newton_surface_solve`, no
  `advance_water_mass_full`. It sees the pure right-hand side and the state algebra and nothing
  else, which is what makes agreement between it and an implicit scheme rule out a shared-bug false
  pass. Before, the module imported the BE machinery for the tier's benefit while the oracle itself
  never touched it.
- **One implementation of each test assertion helper** (#191). Twenty-three local copies across
  nineteen files, consolidated behind generic interfaces so all ~1000 call sites compile unchanged;
  net −403 lines. The two families (fatal condition-first, accumulating name-first) are kept
  deliberately — they differ in failure behaviour, not just signature. A new `meds_test_assert`
  holds the assertions and depends on nothing but a kind, because seventeen tests deliberately link
  one narrow library each. Found and closed one coverage hole: `test_biogeochem_dynamics`' local
  helper error-stopped while the shared one accumulates, which turned four assertions into no-ops
  until a verdict call was added.

### Fixed

- **`soil_carbon_on` now actually defaults to on.** PR #144 flipped the in-type default to
  `.true.` but the TOML loader still passed `.false.` as its absent-key default, and
  `toml_logical` returns the supplied default when the key is absent. Every config that
  omitted the key, including the shipped `meds_config_main.toml`, therefore ran with soil
  carbon off. Off is not a coarser soil-carbon model, it is no soil carbon: litter is
  discarded and `rh = 0`.
- **`forcing.lwdown_source = "synthesize"` is rejected rather than silently ignored.** The
  value parsed to `LW_SYNTHESIZE` but the met driver always read LWdown from the file, so it
  selected the file path under a name promising Brutsaert/Idso synthesis. `validate_config`
  now stops on it until the synthesis exists.

### Changed

- Test coverage for the three untested fast-loop state combinators, with two silent
  exclusions named (#147).
- `soil_carbon_on` default flipped to `.true.` in the type (#144). Measured over a year at
  Ithaca, off reports annual-mean NEE at −4.657 against −2.464 µmol m⁻² s⁻¹, an 89 % stronger
  apparent sink, with Rh identically zero against 0.833 kgC m⁻² yr⁻¹, for no wall-clock saving.
- `[energy].phase_change` deleted; ice-aware soil thermal properties are unconditional (#143).
  The flag only ever gated ice-aware `κ_sat(f_liq)` and `C_eff(f_liq)`; the difference was
  ≤ 0.61 K and 35 692 substeps either way.

### Added

- **The full coupled model is drivable from Python** (`meds.model.Run`): open, step, finalize
  through a C-API shim (#139). The Python and executable paths differ only by libm-versus-libimf
  interposition: worst relative difference ~1×10⁻¹² over a simulated day.
- **One `libmeds.so`** built by scikit-build-core, with a mandatory ctest target per C-API shim
  so an ABI change is a build failure in a default build rather than a silent break in an
  optional one (#138).
- **A slow-loop conservation ledger**, with birth and death as paired transfers (#132).

### Fixed (slow-loop conservation, #132–#137)

The ledger closed a series of real carbon, energy and water leaks that a green test suite and
the per-store fast budgets had both missed:

- **Growth respiration was 23 % of growth carbon that was never exhaled.** The allocator's
  outputs now all have a destination (#133).
- Tissue thermal mass is an exchange, not an appearance (#134).
- The time-averaged diagnostics now close for energy and water everywhere (#135).
- Mortality is valued on what the applier actually removed, not on what was requested (#136).
- Reproduction carbon became a flow rather than a disappearance, which closed the ledger
  (#137). Side effect: the recruit pool now accrues every slow step instead of as a monthly
  lump, so the first cohorts appear about a month later and the early trajectory is offset.
  The offset decays as the stand fills: 70 % in AGB at year 2, 1.6 % at year 10, 0.1 % at
  year 40.
- **Soil respiration reached the atmosphere 964× too small** — a CENTURY Rh unit error, hidden
  because no shipped config turned soil carbon on. The seam check that would have caught it was
  computed but never reported (#140).
- An out-of-bounds litter read: a per-patch quantity indexed outside the patch block when
  disturbance added a patch mid-step. The ledger declared the same garbage it consumed, so
  conservation balanced on it (#139).

### Changed (source-tree restructure, #125–#127, #141)

The tree is now **timescale-first for processes** over a **state layer in two halves**. Moving a
file changes only CMake wiring, because Fortran `use` is by module name, so every step was
verified byte-identical on both back ends.

- Steps 1–6: state as a layer (`state/column`, `state/site`), processes timescale-first
  (`fast_dynamics/`, `slow_dynamics/`) (#125). `src/shared/` dissolved into `base/`,
  `functions/`, `config/`, `state/`.
- Steps 9–10: every façade module dissolved; the fast-loop argument vocabulary renamed so
  names tell the truth — `wcap`/`ccap` → `cas_mass_capacity`/`cas_molar_capacity`,
  `gah`/`gaw`/`gac` → `g_atm_heat`/`g_atm_vapour`/`g_atm_co2`,
  `uext_to_temp`/`temp_to_uext` → `internal_energy_to_temp`/`temp_to_internal_energy`, and the
  three different things called `hydro` split into `soil_water_opts` versus
  `hydraulics_params`/`hydraulics_opts` (#126). netCDF registry strings and TOML keys were
  deliberately not renamed, so output files and configs are unchanged.
- Steps 0 and 7: the fast-loop state vector (#127).
- The library name `demography` was renamed to `core` in July and back to `demography` here:
  `core` stopped carrying information once its state half became `state/site`.
- Structure plan §15 phased remainder, Phase 1 executed (#141).

### Changed (configuration, #129, #130)

- **The soil column comes from config, not from literals in the driver** (#129): the new
  `[soil_column]` block carries layer count, depth, grid growth, hydraulic texture, retention
  family, root profile and the three thermal properties, validated at load. `depth` is the knob
  for the known too-shallow-column defect (2.0 m against a ~2.5 m annual damping depth).
- **The fast driver's hard-coded parameters found their homes** (#130). `agf_bs` duplicated the
  per-PFT `aboveground_frac` with a hard-coded 0.7, so a run whose PFTs differed in allocation
  used their values everywhere except stem respiration (#128). Canopy optics — leaf and wood
  reflectance and transmittance per band, clumping, leaf-angle mean and standard deviation —
  became `[pft]` traits, so two PFTs can finally differ in how they intercept light (#131).
  Longwave is configured as emissivity, with reflectance derived as `1 − emissivity` and
  transmittance zero, because a leaf is opaque at thermal wavelengths.

### Fixed (2026-09 review, #119–#124)

- Conservation fixes: a closed whole-column energy ledger, RK45 and threading defects,
  film-water and seam conservation (#119). A 50-year spin-up's energy imbalance went from
  +6.75 MJ m⁻² to 3.6×10⁻¹¹ W m⁻².
- Dead code removed, shared column-state algebra extracted, `t_ground` made explicit, one CAS
  box kernel instead of two (#120).
- The fast-loop argument vocabulary renamed; over-long lines wrapped (#122).
- `column_prepass` split into the five processes it fused; a per-PFT leaf photosynthesis table;
  scalar signatures; shared assemblers; ledger-consistent latent and sensible heat (#123).
- The frozen work record decomposed by physical content; `integrator_opts_t`;
  `apply_process_mask` (#124).

---

## [0.1.0] — 2026-08-02

First tagged release. **Unbenchmarked**: no EDTS-equivalent regression suite has been run, no
site compared flux-for-flux, no output scored against observations. What is verified is
internal — the test suite on two compilers, per-step conservation ledgers, and thread-invariant
output.

### Added

- **The v0.1 diagnostic output layer** (#111): per-variable, per-timescale output control.
  ~208 registered variables across 8 groups and 7 axes (cohort, patch, site, soil, PFT,
  DBH size class, and the 2-D soil×patch slab), each switchable individually per timescale
  from TOML. Extensive quantities carry their own aggregation weight, so one registry line
  emits a cohort field's patch, site, PFT and size-class rollups. `meds_main --dump-io-config`
  lists every variable. Verified byte-identical at 1 versus 4 threads across all 75 files of a
  3-year run.
- **The ED2-to-MEDS comparison** for people who already run ED2 (#113, #115),
  [`docs/ed2_comparison.md`](docs/ed2_comparison.md).
- **Patch-axis threading** (#109): byte-identical output at 1, 2, 4 and 8 threads. 2.03× at
  4 threads against a measured 3.03× hardware ceiling. Three silent compiler traps were found
  and worked around: Intel's `-auto-scalar` default placing local arrays in static storage,
  nvfortran rejecting `BLOCK` inside a parallel region, and ifx building `private` copies of a
  derived type through a compiler-generated static mold.

### Changed

- **`dt_fast` became an accuracy parameter, not a stability one** (#90). The period-2
  canopy-air oscillation was traced to one frozen coefficient — the canopy-air-to-atmosphere
  conductance — and re-solving it at every integrator stage removed the oscillation while
  *reducing* integrator work. The production default went from 150 s to **900 s**. The
  non-stomatal water-stress limb was gated off by default in the same change.
- The legacy `[io]` diagnostic writer was retired at v0.1: an annual-cadence, instantaneous,
  hard-coded 21-variable schema that duplicated `[output]` and collided with it on the `-D-`
  filename prefix. `meds_io` now carries the state (restart) stream only.

### Fixed

- **The transpiration seam in the plant water update** (#91): a pure flux inconsistency,
  exactly `dt·(transp_pp − transp_bw)`, corrected by a transpiration corrector. Improved
  `psi_leaf` convergence 314×. ARK only; RK45 does not carry the corrector.
- Stomatal water-stress closure (#98, issue #95): `beta_stomata` was identically 1 until fixed.
  Restart persistence for the daily-maximum leaf water potential.
- Pathological hydraulics sub-stepping is now detected (#105): a collapsed, floored wood store
  burned 13× the wall clock, silently (issue #104).
- The adaptive warm start got per-patch storage (#108, issue #106) — it had been loop-carried,
  so patch 1 cold-started and patches 2..N inherited a neighbour's step size.

### Documentation

- **GPU offload evaluated against measurement and found not viable as scoped** (#110). The GPU
  build ran 1.4× *slower* than the CPU (49.8 s against 36.2 s), with one kernel at 0.4 %
  occupancy; the GPU treated as 20 slow cores was 26× slower than 4 CPU cores. The
  recommendation is to thread and vectorise the cohort axis on the CPU instead.

---

## Pre-0.1 development — 2026-06-23 to 2026-08-02

MEDS began on 2026-06-23. The sections below group the first six weeks by subsystem rather than
by date, because the work proceeded as a dozen parallel subsystem builds.

### Demographic core

- Source tree organized by process domain; per-PFT `hgt_max` (#3).
- The engine became **carbon-driven**: carbon pools on the cohort structure-of-arrays (#17),
  carbon-driven growth wired in (#19), the empirical vital-rate laws moved out to the Python
  example and the engine reduced to law-free apply-primitives (#16, #18, #43).
- The demography engine reorganized into `core` with a tendency-seam growth/mortality
  interface (#44, #45).
- **Run-model decision** (#63): the fast biophysics loop is always on, and phenology is
  unconditional. A slow-only, empirical-demography run is the Python C-API path.

### Fast-loop biophysics

- **Canopy radiative transfer** (#5): a faithful modernized ED2 two-stream (`icanrad = 2`) with
  a unified multi-band solver, SCOPE/4SAIL leaf-angle scattering over a Beta leaf-angle
  distribution, and an in-house block-tridiagonal solve. ED2 bugs found during the port are
  recorded in issue #6.
- **Plant hydraulics** (#9, #49): a stateless matrix-exponential network solver, shared
  constitutive curves (pressure-volume, Kirchhoff conductance, xylem vulnerability) extracted
  into `meds_hydr_lib`, and multi-layer root water uptake.
- **Leaf gas exchange** (#2, #46): FvCB C3 and Collatz C4 demand, Leuning / Medlyn / Katul
  stomatal models, a bracketed C_i solver, and Sabot two-limb water stress.
- **Column soil hydrology** (#22, #23): implicit backward-Euler Thomas Richards with Celia
  modified-Picard linearization, upstream-weighted conductivity, adaptive step-doubling,
  infiltration and ponding, and a free-drain / bedrock / aquifer bottom boundary.
- **Energy balance** (#24, #25): four stateless per-store thermal kernels carrying prognostic
  **internal energy rather than temperature**, so the freeze/thaw plateau is a read-off of the
  enthalpy inverter rather than a special case.
- **Canopy air space CO₂** (#26, #27): `can_co2` as the third prognostic CAS twin, with a DAMM
  heterotrophic-respiration option.
- **Fast-loop coupling capstone** (#28): the sub-daily loop coupled and owning its own state.
- **Prognostic leaf and wood energy** (#41) with a separate wood temperature; FAST diurnal
  diagnostics in the same change.
- **Snow and temporary surface water** (#42): stateless kernels and a conserving fast-loop
  coupling that closes whole-column mass and energy budgets to machine precision through
  accumulation, sublimation, melt into infiltration, and the snow-albedo ramp.

### Numerics

- **IMEX-ARK** (2026-07-11, e77f7f0, no PR number): an L-stable ESDIRK2 on the ARS(2,2,2)
  tableau with an arrowhead Newton surface solve. The config string stays `"ark"`; the explicit
  part is empty, so despite the name it is a diagonally implicit scheme.
- **Error-control infrastructure** (#65) and tolerance unification, a process mask, a sweep
  harness (#66).
- **ED2-faithful RK45** (#67, #68): an adaptive Cash-Karp 5(4) march, internal water carried as
  **mass** rather than potential, and a mass-conserving demographic seam.
- **Integrator parity** (#77, #80–#87): a long sequence making the schemes solve the same
  physics. Along the way: the reference used to score them was mis-timed; ARK and RK45 dropped
  snowfall all winter; `veg_coupling_floor` destroyed energy; the split scheme's soil-energy
  budget was a tautology; the pond held no enthalpy so water crossing into it shed its heat
  while the books still closed.
- **The operator-split integrator was retired** (#88). It converged to a different limit
  (~0.45 K in canopy-air temperature) that refinement never removed and nobody ever attributed,
  and it could not carry the coupled tissue heat store. `ark` became the default and `rk45` the
  accuracy baseline, with the RK45 stiff rescue redoing the step on `ark`. The tissue heat
  store was turned on in the same change, integrated by an exact exponential with two weights
  (endpoint and step-average) because the tissue ODE is linear under the frozen coefficients.

### Slow loop

- **Phenology**: a stateless signal module (#10), wired into the run loop (#39), then refactored
  to a **rate-based** signal-only kernel emitting two relative tendencies instead of a
  directional tri-state (#51).
- **Plant carbon**: pure carbon-dynamics kernels (#14), carbon-pool allometry and PFT traits
  (#15), then an **elemental growth-allocation kernel** following FATES PARTEH Hypothesis-1
  (#52), with growth respiration charged inside the kernel on realized growth.
- **Trait plasticity** (#54): light-driven acclimation of specific leaf area, V_cmax, dark
  respiration and leaf lifespan.
- **Non-leaf maintenance respiration** (#13): stem and fine-root, an ED2 Chambers-2004 port.
- **Slow soil-carbon biogeochemistry** (#35): ED2's CENTURY decomposition reorganized as the
  carbon matrix ODE `dX/dt = B·I + A·ξ·K·X`, with a 7-pool state, a lignin sub-tracer, an exact
  augmented matrix exponential for accelerated steps, and a SASU steady-state solve. Wired into
  the slow loop (#64): the fast loop's heterotrophic respiration respires the same frozen pool
  the daily step debits, so the day's total fast Rh equals the daily debit by construction.

### Forcing

- **Meteorological forcing** (#36): a single-site NetCDF reader over a multi-grid `(time, grid)`
  file, `pure`/`elemental` disaggregation kernels including an interval-mean-conserving
  shortwave reconstruction, the canopy-RT join, net longwave, Weiss-Norman band-specific
  shortwave, multi-year calendar recycling, nearest-grid matching, and wind-height and
  elevation lapse. MEDS never gap-fills: a missing or NaN required value is a hard error.
- **The recycle window became declared rather than inferred** (#69). The previous classifier
  accepted only Jan-1 00:00 files and silently fell back to an absolute-seconds span wrap
  otherwise. Real ERA5-Land records are stamped at the *end* of each interval, so their first
  record is 01:00:00 and they always took that fallback: the Ithaca file spans 366 d 22 h, and
  a 29-year run ended up reading late May at a ~10 h offset. Nothing caught it, because the
  cosz reconstruction is mean-conserving, so daily-mean shortwave stayed correct and the slow
  demography looked healthy.

### Output and I/O

- **Diagnostic aggregation and output subsystem** (#38, #40): the registry, the per-tier
  temporal integrators, and the netCDF serializer, written through the netCDF **C** library via
  `iso_c_binding` so the output layer builds under ifx and nvfortran (netCDF-Fortran's module
  format is gfortran-only).
- A **FAST (sub-daily) output tier** (2026-07-13, 01347bb, no PR number) for diurnal-cycle
  analysis.

### Code review and bug fixes

- **Adversarial code review, 2026-07-06**, over ~10.4k lines and 47 modules, targeting
  physical-process bugs, numerical defects and organization. All sections addressed across
  #29–#34: 8 critical/high physical-process bugs, then medium, then low-priority guards, then
  performance and solver issues, then organization.

### Tooling and infrastructure

- CMake with automatic Fortran module-dependency resolution — the deliberate fix for ED2's
  "run `make` six times" hack and per-platform `include.mk` files.
- **nvfortran portability trap documented** (#8, issue #7): never pass an array-valued function
  result straight into a call. nvfortran's whole-program optimizer miscompiles the temporary
  descriptor — silently wrong values at `-O2`, segfault at `-O0` — while
  `ifx -stand f18 -check all` tolerates it. A green ifx run is not sufficient.
- 3D visualization of forest structure (#4); the biophysics example (#70).
- Design plans relocated to `docs/dev_plans/` (#48); GitHub math rendering fixed in the science
  pages (#50).

---

[Unreleased]: https://github.com/xiangtaoxu/MEDS/compare/v0.3.2...beta
[0.3.2]: https://github.com/xiangtaoxu/MEDS/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/xiangtaoxu/MEDS/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/xiangtaoxu/MEDS/compare/v0.2.2...v0.3.0
[0.2.2]: https://github.com/xiangtaoxu/MEDS/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/xiangtaoxu/MEDS/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/xiangtaoxu/MEDS/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/xiangtaoxu/MEDS/releases/tag/v0.1.0
