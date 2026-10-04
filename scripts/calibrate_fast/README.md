# calibrate_fast — fast-parameter calibration against a flux tower

Fits MEDS's sub-daily parameters (radiation, photosynthesis and stomata, aerodynamics, water stress)
to a flux tower's fluxes, with the stand frozen at its initial structure. It also gives a rough
uncertainty for each parameter, the starting point for a longer-term calibration. The protocol, and
every decision in it, is
[`docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md):
rules that work at any tower, with the tower's numbers coming from its own data and climate.

Barro Colorado Island is the worked example: [`examples/example_flux_tower_bci/calibration.toml`](../../examples/example_flux_tower_bci/calibration.toml).

**The user decides, the tool recommends and reports.** Every setting of a site's `calibration.toml`
is listed in [`site_reference.toml`](site_reference.toml) with its default and the reason for it.
That file is also the schema: a setting it does not list stops the tool, so a misspelt key cannot
silently fall back to the default.

## How it works

### The tower

- **The site TOML holds the facts, `calibration.toml` the choices.** `[tower].site` names the
  tower's site TOML, the file the forcing build reads
  ([`../prepare_flux_tower`](../prepare_flux_tower/README.md)). It declares:
  - the file, the clock, the columns and their units;
  - in `[fluxes]`, each target's column with the rule that says when it was measured (BCI:
    `FLAG = 1`);
  - in `[provider]`, what the provider did;
  - the leaf-on months.

  Both tools read it through one reader, `tower_inputs.read_standard`.
- **The tower's own interval, on UTC starts.** The observations stay on the tower's interval (30 min
  at BCI). The tool sets `[output].fast_interval_steps` so each trial writes one record per tower
  interval, and a trial whose output is at another spacing fails.
- **The metadata checks** (F1–F3): net radiation against its four components, the upwelling longwave
  above the downwelling at night (a swapped pair of columns), and the reflected shortwave and
  upwelling longwave within bounds. The forcing build stops on a failure; `report` and `fit` list
  it under `tower_checks`.

### Targets, filters and error models

- **Targets:** the upwelling longwave, H, LE, GPP and u\*, on the records the tower measured and the
  forcing observed.
  - **The albedo is off for now** (`[targets.albedo].on`). At realistic leaf optics the model's canopy
    reflects too much NIR, so the albedo could be met only by moving the leaf's NIR reflectance past
    its evidence (at BCI to 0.32, prior z −3.1). The NIR reflectance is fixed with it.
  - **Not targets:** net radiation (an input to the closure model), the evaporative fraction (it
    repeats H and LE) and night NEE.
  - **H, LE and u\* are daytime targets.** At night they are set mostly by the model's numerical
    floors and stable-air measurement problems; `night = true` keeps night records at the provider's
    u\* threshold.
  - **No hours are cut for transitions.** Canopy storage is ignored: every tower flux is treated as
    turbulent.
- **u\* per target, by its own diagnostic** (`datarules.py`, `[ustar]`).
  - **The test:** within classes of the target's driver (PAR for GPP, Rnet for LE and H), the
    flux-to-driver ratio across u\* classes. The plateau test (Papale et al. 2006) finds the lowest
    class within 95 % of the mean of the classes above, and a bootstrap over days gives its spread
    (Barr et al. 2013).
  - **The outcome:** a plateau is filtered at its threshold. A flat ratio, or one still rising at
    the top classes, is not filtered: a rising one is the observation model's to handle.
  - **GPP** takes the provider's CO₂ threshold (`ustar_min = "provider"`, one value or one per year).
  - **At BCI:** LE is flat and H rising, so neither is filtered; GPP takes the provider's 0.4. The
    all-day GPP diagnostic finds a plateau at 0.33 (bootstrap 0.15–0.35).
- **The other filters** (`par_min`, `hours`, `min_solar_elevation`, `snow_free_days`) are per-target
  settings. The albedo's, for when it is on, are the plan's rule: SW_in > 200 W m⁻², the sun at least
  20° high, a week without frost.
- **Each target has an observation model** (`obsmodels.py`): how the tower's number relates to the
  true flux, each known bias as its own term.
  - **H and LE: closure.** The tower misses the same fraction of the turbulent flux at every hour.
    - **The fraction:** f_d = median over ±15 days of the daily (Rnet − G)/(H + LE), from days with at
      least 70 % of their records measured; over a day, storage nets out.
    - **The shares come from the attribution test:** within VPD classes, which of H/Rnet and LE/Rnet
      rises from the calmest to the most turbulent third of the records.
    - **At BCI:** H/Rnet rises 44–64 % in every class and LE/Rnet does not, so H takes the whole gap.
      With f_d = 1.33 (daily closure 0.75), H rises to close the day and LE stays as measured.
    - **Bowen** (`[closure].shares = "bowen"`) is the declared alternative.
  - **GPP: respiration.** The provider made GPP from its own respiration, R_tower = κ R_true. The
    residual is [GPP_model − GPP_tower − (1/κ − 1) R_tower]/σ_NEE.
    - **κ is an observation key.** It is fitted with a prior, never written to a MEDS config, and its
      gradient costs no trials.
    - **At BCI** the prior is 0.65 ± 0.10: the tower's respiration is about the soil chambers' alone.
    - **κ is never fitted beside a free shape key** of the leaf's light response.
- **Random error,** σ = `sigma_abs` + `sigma_rel` |x|.
  - **Its source,** in order: the provider's per-record uncertainty, the paired days (Hollinger &
    Richardson 2005), then the defaults.
  - **It is evaluated at a smoothed observation:** the mean at the same time of day within ±7 days.
    GPP's σ is NEE's.
  - **At BCI** the paired days give LE 10.5 + 0.29|LE|, H 7.5 + 0.14|H| and NEE 2.1 + 0.17|NEE|.
- **Huber loss by default** (c = 2): half-hourly errors are heavy-tailed. χ² is reported on the
  plain residuals.

### Windows and states

- **Trials.** A trial is a frozen run (`[run].slow_on = false`) restarted from its window's state, with
  `[init].reacclimate_traits = true` so a changed trait reaches the cohorts. Its main and PFT files
  are the base configuration with keys set, through `meds.config`. Its directory is named by a hash
  of them, so a repeated candidate reuses the finished run.
  - **Trials run through the Python API** (`python -m meds.model`, which needs `libmeds.so`), or
    through `meds_main` with `--runner`. The two give the same output, bit for bit.
  - **Every trial proves what it ran.** The model's parameter record must list every key the trial
    set, marked as set in the file, with the value written.
  - **A trial fails** when it breaches a whole-site budget, ends without `OK: simulation completed`,
    has a NaN in its output (G12), or runs past 3× the median trial time. A failed trial rejects its
    step; it does not stop the fit.
- **Calibration windows, by rule** (`[windows]`; `select-windows` prints the choice).
  - There is one ten-day window per 1.5-month slot of the year inside the leaf-on months, each the
    slot's best covered: the share of daytime records with H, LE and GPP measured, times the share
    with observed forcing.
  - Validation windows are in the same slots of other years.
  - A site may list its own windows.
- **Seasonal runs, by rule** (`[windows.seasonal]`).
  - There are up to two 120-day runs, each ending at a year's deepest cumulative water deficit (rain
    minus Priestley–Taylor evaporation, from the forcing) inside the leaf-on months.
  - They are scored on LE only.
  - A year under 100 mm gives none; without any, the drought keys are fixed.
  - At BCI the rule picks the 2016 and 2017 dry seasons (1,331 and 698 mm).
- **States.** Every window has its own chain: a frozen run from the initial stand (the base config's
  census or state) that starts `[windows].chain_lead_days` (180) before it. The chains run at once.

### Keys and priors

- **The registry** (`parameters.toml`) is a menu.
  - **States:** `fit` (the default set), `optional` (fitted when the site asks) and `fixed` (with the
    reason).
  - **Kinds:** trait (measurable, prior from evidence), effective (a scheme property, labelled in the
    calibrated files), numerical (never calibrated) and observation (κ).
  - **Scopes:** plant type, site, observation.
  - The site chooses the keys in `[fit]` (`keys`, or `add`/`remove`) and sets any prior or range in
    `[priors.<key>]`.
- **The range is what is physically possible, the prior what the evidence says.** Every fitted key
  has its own sd. Its prior z at the MAP is reported, and gate G13 asks for a diagnosis of any trait
  key beyond 2 sd.
- **The priors from the site's climate** (`priors.py`; the forcing only, never the tower's fluxes):
  - **`stomatal_g1`:** the least-cost value, ξ = √(β (K + Γ\*)/(1.6 η\*)) with β = 146, at the
    growing-season daytime climate with MEDS's Rubisco kinetics; log-sd 0.5. At BCI: 2.80. Lin et
    al.'s 3.77 for tropical rainforest is 0.9 sd away.
  - **`vcmax25`:** the coordination of the Rubisco- and light-limited rates, solved with MEDS's own
    leaf (`meds.plant.leaf`, so libmeds); log-sd 0.5. At BCI: 41 with phi_psii 0.74 and θ_J 0.7.
  - **Jmax/Vcmax, `ds_vcmax`, `ds_jmax`:** fixed at Kattge & Knorr for the growth temperature, set in
    every run and written into the calibrated files. At BCI (25.5 °C): 1.70, 641.1 and 640.6.
- **Priors by plant type** (`[fit].plant_type`, the registry's `meta`) are the prior of a key without
  an EEO centre. Beside an EEO centre, the type's prior is reported and flagged when more than 2 sd
  apart. Site leaf data in `[priors]` come first.
- **Keys fixed from coverage.** A key that acts through one process (the registry's `process`: wet
  canopy, night, snow, drought) is fixed when the kept records sample that process fewer than
  `[fit].min_process_records` (48) times. At BCI the flag removes every rain half hour.

### The joint fit

One fit of every key at once (best-practice plan §6.3):

1. **Screening.** The first gradient matrix, by central differences at the prior centres. With
   `[fit].screening = "triage"`, a key is fixed at its prior's centre when it is:
   - dead (a zero column);
   - rough (the two one-sided slopes disagree 3×: make the model continuous);
   - uninformed (posterior/prior sd ratio ≥ 0.9);
   - or the less informed of a pair correlated beyond 0.95.
2. **Levenberg–Marquardt on every fitted key at once.**
   - **Rows:** the calibration windows (every target) and the seasonal runs (LE).
   - **Start:** the prior centres, one start.
   - **Iterations:** one-sided differences. Broyden's update carries the gradient matrix between
     full recomputations, every third accepted step and after a rejected one.
   - **Stopping:** at most 10 iterations, or when the cost falls by under 0.1 %.
3. **Once, after the first convergence, a refresh:**
   - re-run the chains at the current values;
   - refresh the effective-sample weights;
   - scale each target's σ by the model's misfit, max(1, √(χ²/n)) at most 3×;
   - then continue to convergence.
4. **A final central-difference gradient matrix** for the uncertainty.

The separate water stage, the polish, the kernel stages (and their gate G8) and multiple starts are
gone.

### Three uncertainties side by side

1. **Laplace, from the final gradient matrix,** with the σ scales applied.
   - **Intervals** (68 % and 95 %), mapped back through the transform: asymmetric, and inside the
     bounds.
   - **The linearity check** along the 3 leading directions: the objective at ±1 sd against the
     quadratic. Outside 0.5–2× the quadratic's curvature, the covariance is marked "local only".
   - **Convergence:** the Gauss–Newton step left at the MAP, per key in posterior sd.
2. **The declared alternatives,** the same at every site (`[uncertainty].alternatives`):
   - **GPP u\*:** the provider's threshold against the daytime plateau;
   - **closure:** the attribution shares against Bowen;
   - **partitioning:** night-time against daytime, where both exist.

   Each alternative's shift of the MAP is first estimated from the final gradient matrix, from the
   cached trials, so it costs no runs. Above 1 posterior sd (`refit_sd`), the fit is rerun from the
   MAP under that alternative and both MAPs are reported (gate G10).
3. **Structural variants.** The same keys, priors and filters under each structure in
   `[variants]` (interception on and off). `calibrate_fast.py variants` puts the variants' fits side
   by side: each key's MAP, and their spread in posterior sd.

**The report** (`report.md` beside `fit.json`) gives:
- the gates;
- per key: kind, scope, MAP and interval, prior centre, source and z, and posterior/prior σ ratio;
- per target: χ²/n, the σ scale, and the model/tower ratio overall, by local hour and by light
  class;
- κ, with the respiration and GPP it implies;
- the three uncertainties.

### Validation and gates

The default and the MAP are scored on validation windows the fit never saw, each from its own
chain's states.

| gate | what |
|---|---|
| G1, G2 | `check`: the stand is identical at a trial's start and end; the same parameters give byte-identical output |
| G3 | every fitted key has a non-zero, smooth gradient column |
| G4 | the MAP beats the default on validation, no target more than 10 % worse |
| G5 | keys near a bound, each with the target that pushed it |
| G7 | the full record with the slow tier on (run the calibrated configs) |
| G10 | every declared alternative's shift is under 1 posterior sd, or its refit is reported |
| G12 | no scored output has a NaN (a trial with one fails) |
| G13 | every trait key with \|prior z\| > 2 is diagnosed or relabelled effective |

## Commands

```
calibrate_fast.py select-windows --site calibration.toml                 # print the rule's windows and seasonal runs
calibrate_fast.py report --site calibration.toml --work runs/report      # the data report only
calibrate_fast.py check --site calibration.toml --variant interception_on --work runs/check \
    --workers 8                                                         # runs, record, G1, G2
calibrate_fast.py fit   --site calibration.toml --variant interception_on --work runs/on --workers 40
calibrate_fast.py analyze --site calibration.toml --variant interception_on --work runs/on \
    --workers 40                                                        # redo the post-fit steps
calibrate_fast.py write-calibrated --site calibration.toml --variant interception_on \
    --fit runs/on/fit.json --out calibration                            # the configs with the MAP
calibrate_fast.py variants runs/off/fit.json runs/on/fit.json --out runs/variants.json   # the structures side by side
calibrate_fast.py smoke --site ... --work ...                           # the ctest smoke test
```

- **Runners:** add `--runner ../../build-ifx/meds_main` to `check`, `fit` or `analyze` to run the trials
  with the executable. The EEO vcmax25 prior needs `libmeds.so` in any case.
- **Re-running:** a re-run of `fit` replays its points from the trial cache. `analyze` reruns
  everything after the fit from the MAP in `fit.json`.
- **On a Slurm cluster,** `--pool queue` uses a directory queue under `--work`. Start one
  `calibrate_fast.py worker --queue <work>/queue --slots <cores>` per node inside the same
  allocation; the driver hands trials to them.

```
srun --ntasks-per-node=1 --cpus-per-task=40 python calibrate_fast.py worker --queue runs/on/queue --slots 40 &
python calibrate_fast.py fit ... --pool queue
touch runs/on/queue/STOP      # the driver writes it too, when the fit ends
```

Ask Slurm for the memory with `--mem=0` (the whole node). A trial takes up to ~1 GB, and a cluster
whose default is 1 GB per job kills the workers.

## What `fit` writes (in `--work`)

| file | what |
|---|---|
| `fit.json` | everything below |
| `report.md` | the results for a reader: gates, keys, targets, κ, the three uncertainties |
| `pft_parameters_calibrated.toml`, `meds_config_calibrated.toml` | the base configs with the MAP and the Kattge & Knorr values written in, comments kept, effective keys labelled |
| `fit.log` | the run's log |
| `trials/`, `chains/` | every trial and state (cached: a re-run reuses them) |

`fit.json` holds:
- the data report, the climate priors, the coverage, and the screening;
- each Levenberg–Marquardt phase's history, and the refresh's weights and σ scales;
- the MAP, start and default values, the intervals, the covariance and correlations;
- each key's posterior-to-prior σ ratio and prior z;
- the linearity check, the declared alternatives, the model/tower ratios, κ's implications, and
  the scores on the calibration and validation windows;
- the gates, and the trial count and timing.

## Files

| file | role |
|---|---|
| `calibrate_fast.py` | the commands, the joint fit, the post-fit analysis, the calibrated configs |
| `site_reference.toml` | every site setting, its default and why (the schema) |
| `settings.py` | reads a site declaration against `site_reference.toml` |
| `parameters.toml` | the registry: each key's state, kind, scope, file, TOML key, range, transform, prior and source |
| `registry.py` | the transforms, the priors, the key selection |
| `priors.py` | the growth climate, the EEO centres (least-cost g1, coordination vcmax25 with MEDS's leaf), Kattge & Knorr |
| `datarules.py` | the u\* diagnostic, the process coverage, the window rule, the water-deficit index and seasonal runs |
| `obsmodels.py` | the closure model and its attribution test, the respiration model's κ prior, the random error at a smoothed observation |
| `tower.py` | the tower's records (through `tower_inputs`), the observed-forcing mask, the sun's elevation |
| `residuals.py` | the targets, the filters, the residual vector, the weights, the σ scale, the filter report |
| `fit.py` | Levenberg–Marquardt with Broyden updates, the gradient matrices, the triage, the covariance and intervals, the linearity check, prior z, the filter shift |
| `trials.py` | the trial writer and runner, the parameter-record check, the fast-series reader |
| `states.py` | the state chains |
| `pool.py` | the local pool and the directory queue |
| `tests/` | unit tests, a smoke test and an end-to-end fit that run the model (ctest `calibrate_fast`) |
