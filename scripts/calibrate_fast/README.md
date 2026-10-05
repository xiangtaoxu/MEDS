# calibrate_fast — fast-parameter calibration against a flux tower

Fits MEDS's sub-daily parameters (radiation, photosynthesis and stomata, aerodynamics, water stress)
to a flux tower's fluxes, with the stand frozen at its initial structure. It also gives a rough
uncertainty for each parameter, the starting point for a longer-term calibration. The protocol, and
every decision in it, is
[`docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md):
rules that work at any tower, with the tower's numbers coming from its own data and climate.

Barro Colorado Island is the worked example:
[`examples/example_flux_tower_bci/calibration.toml`](../../examples/example_flux_tower_bci/calibration.toml),
with what the rules chose there and the results in the example's README.

**The user decides, the tool recommends and reports.** Every setting of a site's `calibration.toml`
is listed in [`calibration_reference.toml`](calibration_reference.toml) with its default and the
reason for it. That file is also the schema: a setting it does not list stops the tool, so a
misspelt key cannot silently fall back to the default. A calibration file holds only the choices
that differ from the defaults.

**When something is wrong, the tool stops and says why.** It stops when:
- a trial fails;
- a window has no usable record;
- a declared forcing qc variable is missing;
- the tower's flux checks fail;
- a target's observation model lacks its data.

It does not work around a problem silently.

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
- **The forcing is the base configuration's own** (`[forcing].path` and `grid_index` of the base
  main file), so the windows, the water-deficit index and the climate priors read the file the
  runs read.
- **The tower's own interval, on UTC starts.** The observations stay on the tower's interval (30 min
  at BCI). The tool sets `[output].fast_interval_steps` so each trial writes one record per tower
  interval, and a trial whose output is at another spacing fails.
- **The flux checks** (F1–F3) are the same as in the forcing build, and they stop both tools:
  - net radiation against its four components;
  - the upwelling longwave above the downwelling at night (a swapped pair of columns);
  - the reflected shortwave and upwelling longwave within bounds.

### Targets, filters and observation models

- **Targets:** the upwelling longwave, H, LE, GPP and u\*, on the records the tower measured and the
  forcing observed.
  - **The albedo is off for now** (`[targets.albedo].on`): at realistic leaf optics the model's
    canopy reflects too much NIR, so the albedo could be met only by moving the leaf's NIR
    reflectance past its evidence. The NIR reflectance is fixed with it.
  - **Not targets:** net radiation (an input to the closure model), the evaporative fraction (it
    repeats H and LE) and night NEE.
  - **H, LE, GPP and u\* are daytime targets.** At night they are set mostly by the model's
    numerical floors and by stable-air measurement problems.
  - **No hours are cut for transitions.** Canopy storage is ignored: every tower flux is treated as
    turbulent.
- **u\* per target, by its own diagnostic** (`data_rules.py`, `[ustar]`).
  - **The test:** within classes of the target's driver (PAR for GPP, Rnet for LE and H), the
    flux-to-driver ratio across u\* classes. The plateau test (Papale et al. 2006) finds the lowest
    class within `plateau_fraction` of the mean of the classes above, and a bootstrap over days gives
    its spread (Barr et al. 2013).
  - **The outcome:** a plateau is filtered at its threshold. A flat ratio, one still rising at the
    top classes, or too few records for the test leaves the target unfiltered; a rising ratio is the
    observation model's to handle.
  - **GPP** takes the provider's CO₂ threshold (`ustar_min = "provider"`, one value or one per
    year). Its own diagnostic's plateau is the declared alternative.
- **The other filters** (`par_min`, `hours`, `min_solar_elevation`, `snow_free_days`) are settings of
  any target. The albedo's, for when it is on, are the plan's rule: SW_in > 200 W m⁻², the sun at
  least 20° high, a week without frost.
- **Each target has an observation model** (`observation_models.py`): how the tower's number relates
  to the true flux, each known bias as its own term.
  - **H and LE: closure.** The tower misses the same fraction of the turbulent flux at every hour.
    - **The fraction:** f_d = median over ±15 days of the daily (Rnet − G)/(H + LE), from days with
      at least 70 % of their records measured; over a day, storage nets out.
    - **The shares come from the attribution test:** within VPD classes, which of H/Rnet and LE/Rnet
      rises from the calmest to the most turbulent third of the records. One rising takes the whole
      gap; both rising means Bowen; neither means as measured.
    - **Bowen** (`[closure].shares = "bowen"`) is the declared alternative; `"none"` keeps H and LE
      as measured.
  - **GPP: respiration.** The provider made GPP from its own respiration, R_tower = κ R_true. The
    residual is [GPP_model − GPP_tower − (1/κ − 1) R_tower]/σ_NEE.
    - **κ is an observation key.** It is fitted with a prior, never written to a MEDS config, and its
      gradient costs no trials. Without κ in the fit, κ = 1 and the tower's GPP is taken as given.
    - **Without a site prior,** κ is centred at 1, with its sd from the gap between the provider's
      two partitionings where both exist.
    - **κ is never fitted beside a free shape key** of the leaf's light response.
- **Random error,** σ = `sigma_abs` + `sigma_rel` |x|.
  - **Its source,** the first with enough records: the provider's per-record uncertainty, the
    paired days (Hollinger & Richardson 2005), then the target's defaults.
  - **It is evaluated at a smoothed observation:** the mean at the same time of day within ±7 days.
    GPP's σ is NEE's.
- **Huber's loss** (`huber_c` = 2): half-hourly errors are heavy-tailed. χ² is reported on the plain
  residuals.

### Windows, seasonal runs and chains

- **Trials.** A trial is a frozen run (`[run].slow_on = false`) restarted from its window's state,
  with `[init].reacclimate_traits = true` so a changed trait reaches the cohorts. Its main and PFT
  files are the base configuration with keys set, through `meds.config`. Its directory is named by a
  hash of them, so a repeated set of values reuses the finished run.
  - **How trials run:** through the Python API (`python -m meds.model`, which needs `libmeds.so`),
    or through `meds_main` with `--runner`. The two give the same output, bit for bit.
  - **When a trial passes:** it ends with `OK`, no whole-site budget breaches its tolerance, the
    model's parameter record lists every key the trial set (read from the trial's file, with the
    value written), and the output has no missing or non-finite value.
  - **A trial that fails** (or times out: `[fit].timeout_per_day`) **stops the fit with its log.**
    A model that cannot run a set of values inside the keys' ranges has a bug to fix.
- **Calibration and validation windows, by rule** (`[windows]`; `report` prints the choice).
  - There is one ten-day calibration window per 1.5-month slot of the year inside the leaf-on
    months, each the slot's best covered: the share of daytime records with H, LE and GPP measured,
    times the share with observed forcing (`[tower].forcing_qc`).
  - Validation windows are in the same slots of other years, under the same rule.
  - A calibration may list its own windows.
- **Seasonal runs, by rule** (`[seasonal_runs]`).
  - There are up to two 120-day runs, each ending at a year's deepest cumulative water deficit (rain
    minus Priestley–Taylor evaporation, from the forcing) inside the leaf-on months.
  - They are scored on LE only.
  - A year under 100 mm gives none; without any, the drought keys are fixed.
- **States.** Every window has its own chain: a frozen run from the initial stand (the base
  config's census or state) that starts `[windows].chain_lead_days` (180) before it. The chains
  run at once.

### Keys and priors

- **The registry** (`parameters.toml`) lists the keys a fit may move.
  - **States:** `fit` (the default set), `optional` (fitted when the calibration asks: a tower rarely
    informs it) and `fixed` (with the reason).
  - **Kinds:** trait (measurable, prior from evidence), effective (a scheme property, labelled in
    the calibrated files), numerical (never calibrated) and observation (κ).
  - **Scopes:** plant type, site, observation.
  - The calibration chooses the keys in `[fit]` (`keys`, or `add`/`remove`) and sets any prior or
    range in `[priors.<key>]`.
- **The range is what is physically possible, the prior what the evidence says.** Every fitted key
  has its own sd, and its prior z at the MAP is reported.
- **The priors from the site's climate** (`priors.py`; the forcing only, never the tower's fluxes):
  - **`stomatal_g1`:** the least-cost value, ξ = √(β (K + Γ\*)/(1.6 η\*)) with β = 146, at the
    growing-season daytime climate with MEDS's Rubisco kinetics; log-sd 0.5.
  - **`vcmax25`:** the coordination of the Rubisco- and light-limited rates, solved with MEDS's own
    leaf (`meds.plant.leaf`, so libmeds); log-sd 0.5.
  - **Jmax/Vcmax, `ds_vcmax`, `ds_jmax`:** fixed at Kattge & Knorr for the growth temperature, set in
    every run and written into the calibrated files.
- **Priors by plant type** (`[fit].plant_type`, the registry's `meta`) are the prior of a key without
  an EEO centre. Beside an EEO centre, the type's prior is reported and flagged when more than 2 sd
  apart. Site leaf data in `[priors]` come first.
- **Keys fixed from coverage.** A key that acts through one process (the registry's `process`: wet
  canopy, night, snow, drought) is fixed when the kept records sample that process fewer than
  `[fit].min_process_records` (48) times.

### The joint fit

One fit of every key at once (best-practice plan §6.3):

1. **Screening.** The first gradient matrix, by central differences at the prior centres. A key is
   fixed at its prior's centre when it is:
   - dead (a zero column);
   - rough (the two one-sided slopes disagree 3×: make the model continuous);
   - uninformed (posterior/prior sd ratio ≥ `uninformed_sd_ratio`, 0.9);
   - or the less informed of a pair correlated beyond `max_correlation` (0.95).
2. **Levenberg–Marquardt on every remaining key at once.**
   - **Rows:** the calibration windows (every target) and the seasonal runs (LE), each target's rows
     weighted by its effective sample size.
   - **Start:** the prior centres.
   - **Iterations:** one-sided differences. Broyden's update carries the gradient matrix between
     full recomputations, every third accepted step and after a rejected one.
   - **Stopping:** at most 10 iterations, or when the cost falls by under 0.1 %.
3. **Once, after the first convergence, a refresh:**
   - re-run the chains at the current values;
   - refresh the effective-sample weights;
   - scale each target's σ by the model's misfit, max(1, √(χ²/n)) at most 3×;
   - then continue to convergence.
4. **A final central-difference gradient matrix** for the uncertainty.

### Three uncertainties side by side

1. **Laplace, from the final gradient matrix** (`uncertainty.py`).
   - **Intervals** (68 % and 95 %), mapped back through the transform: asymmetric, and inside the
     bounds.
   - **The linearity check** along the 3 leading directions: the cost at ±1 sd against the
     quadratic. Outside 0.5–2× the quadratic's curvature, the covariance is marked "local only".
   - **Convergence:** the Gauss–Newton step left at the MAP, per key in posterior sd.
   - **The key table:** each key's prior z, whether it is near a bound, and the target that pushes
     it there.
2. **The declared alternatives,** the same at every site (`[uncertainty].alternatives`):
   - **GPP u\*:** the daytime plateau of GPP's own diagnostic, against the provider's threshold;
   - **closure:** Bowen, against the attribution test's shares;
   - **partitioning:** the provider's daytime partitioning, against its night-time one, where both
     exist.

   Each alternative's shift of the MAP is first estimated from the final gradient matrix, from the
   cached trials, so it costs no runs. Above 1 posterior sd (`refit_beyond_sd`), the fit is rerun
   from the MAP under that alternative and both MAPs are reported. An alternative the fit already
   used is skipped.
3. **Structural variants.** The same keys, priors and filters under each structure in
   `[variants]` (interception on and off). `calibrate_fast.py variants` puts the variants' fits side
   by side: each key's MAP, and their spread in posterior sd.

### Validation, and what to check

- **Validation.** The default and the MAP are scored on validation windows the fit never saw, each
  from its own chain's states. It **passes** when the calibrated cost is below the default's and no
  target's RMSE is more than 10 % worse.
- **The key table flags two things to look at:**
  - **a trait key more than 2 prior sd from its evidence:** diagnose it (plan §5.1), or relabel it
    effective;
  - **a key near a bound,** with the target that pushes it there.
- **The full record is the next check:** run the calibrated configs with the slow tier on. The
  budgets must close, and the dry-season GPP and LE must be no worse than the default's.
- **The tool and the model's restart are tested by `check`** (CTest runs it):
  - a repeated trial must match byte for byte;
  - the stand must be unchanged at a trial's end;
  - one gradient column must move the output.

## Commands

```
calibrate_fast.py report --config calibration.toml --work runs/report       # the data report; prints the windows
calibrate_fast.py check  --config calibration.toml --work runs/check --workers 8   # the tool and the model's restart
calibrate_fast.py fit    --config calibration.toml --variant interception_off --work runs/off --workers 40
calibrate_fast.py variants runs/off/fit.json runs/on/fit.json --out runs/variants.json   # the structures side by side
```

- **Runners:** add `--runner ../../build-ifx/meds_main` to `check` or `fit` to run the trials with
  the executable. The EEO vcmax25 prior needs `libmeds.so` in any case.
- **Re-running:** a re-run of `fit` replays every trial it has run from the cache in `--work`. It
  costs only the work that is new, such as the post-fit steps after a change in the tool.
- **On a Slurm cluster,** `--queue` hands the trials to a directory queue under `--work`. Start one
  `calibrate_fast.py worker --queue <work>/queue --slots <cores>` per node inside the same
  allocation; the driver hands trials to them.

```
srun --ntasks-per-node=1 --cpus-per-task=40 python calibrate_fast.py worker --queue runs/on/queue --slots 40 &
python calibrate_fast.py fit ... --queue
touch runs/on/queue/STOP      # the driver writes it too, when the fit ends
```

Ask Slurm for the memory with `--mem=0` (the whole node). A trial takes up to ~1 GB, and a cluster
whose default is 1 GB per job kills the workers.

## What `fit` writes (in `--work`)

| file | what |
|---|---|
| `fit.json` | everything below |
| `report.md` | the results for a reader: the validation, the keys with their flags, the targets, κ, the uncertainty |
| `pft_parameters_calibrated.toml`, `meds_config_calibrated.toml` | the base configs with the MAP and the Kattge & Knorr values written in, comments kept, effective keys labelled |
| `fit.log` | the run's log |
| `trials/`, `chains/` | every trial and state (cached: a re-run reuses them) |

`fit.json` holds:
- the data report, the climate priors, the coverage, and the screening;
- each Levenberg–Marquardt phase's history, and the refresh's weights and σ scales;
- the MAP, start and default values, the intervals, the covariance and correlations;
- the key table (prior z, bounds, which target pushes each key) and each key's posterior-to-prior σ
  ratio;
- the linearity check, the declared alternatives, the model/tower ratios, κ's implications;
- the scores on the calibration and validation windows, and the validation's verdict;
- the trial count and timing.

## Files

| file | role |
|---|---|
| `calibrate_fast.py` | the commands and the fit's sequence |
| `calibration.py` | a calibration: its settings, the base configuration, the tower's facts, the data with the rules applied |
| `calibration_reference.toml` | every calibration setting, its default and why (the schema) |
| `settings.py` | reads a calibration's settings against `calibration_reference.toml` |
| `parameters.toml` | the registry: each key's state, kind, scope, file, TOML key, range, transform, prior and source |
| `parameters.py` | the transforms, the priors, the key selection |
| `priors.py` | the growth climate, the EEO centres (least-cost g1, coordination vcmax25 with MEDS's leaf), Kattge & Knorr |
| `tower.py` | the tower's records (through `tower_inputs`), the observed-forcing mask, the sun's elevation |
| `data_rules.py` | the u\* diagnostic, the process coverage, the window rule, the water-deficit index and seasonal runs |
| `observation_models.py` | the closure model and its attribution test, the random error at a smoothed observation, κ's sd |
| `targets.py` | the targets, their filters, the residuals, the weights, the σ scale, the filter report |
| `trials.py` | the trial writer, its checks, and the trial runner |
| `chains.py` | the state chains |
| `workers.py` | the local workers and the directory queue |
| `fit.py` | Levenberg–Marquardt with Broyden updates, the gradient matrices, the triage, the covariance and intervals |
| `uncertainty.py` | the linearity check, the key table, the declared alternatives |
| `report.py` | the data report, the model/tower ratios, κ, the validation's verdict, `report.md` |
| `calibrated_files.py` | the calibrated configs |
| `tests/` | unit tests, and the `check` and an end-to-end fit with MEDS (ctest `calibrate_fast`) |
