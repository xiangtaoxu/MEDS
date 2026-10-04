# calibrate_fast — fast-parameter calibration against a flux tower

Fits MEDS's sub-daily parameters (radiation, photosynthesis and stomata, aerodynamics, water stress)
to a flux tower's fluxes, with the stand frozen at its initial structure. It also gives a rough
uncertainty for each parameter, the starting point for a longer-term calibration. The design, and
every decision in it, is in two plans:
- [`docs/dev_plans/MEDS_FAST_CALIBRATION_PLAN.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_PLAN.md): the method;
- [`MEDS_FAST_CALIBRATION_REVISION_PLAN.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_REVISION_PLAN.md):
  user-defined filters, keys and priors, the staged fit and the rough uncertainty.

Barro Colorado Island is the worked example: [`examples/example_flux_tower_bci/calibration.toml`](../../examples/example_flux_tower_bci/calibration.toml).

**The user decides, the tool recommends and reports.** Every setting of a site's `calibration.toml`
is listed in [`site_reference.toml`](site_reference.toml) with its default and the reason for it.
That file is also the schema: a setting it does not list stops the tool, so a misspelt key cannot
silently fall back to the default.

## How it works

### Trials and states

- **Trials.** A trial is a 10-day frozen run (`[run].slow_on = false`) restarted from a shared state at
  its window's start, with `[init].reacclimate_traits = true` so a changed trait reaches the cohorts.
  Its main and PFT files are the base configuration with keys set, through `meds.config`, and its
  directory is named by a hash of them, so a repeated candidate reuses the finished run.
- **Trials run through the Python API.** Each trial is `python -m meds.model main.toml` in its own
  process, which needs `libmeds.so` (build with `-DMEDS_BUILD_PYLIB=ON`; `MEDS_LIB` names it if it is
  not in a build directory of this tree). `--runner <path to meds_main>` runs the executable instead.
  The two give the same output, bit for bit (ctest `python_api`).
- **Every trial proves what it ran.** After the run, the trial's parameter record
  (`<prefix>_parameters.csv`, written by the model) must list every key the trial set, marked as set
  in the file, with the value written. A missing or defaulted key stops the fit instead of silently
  running the default. A trial that breaches a whole-site budget (a `budget[whole_*]` line with
  `fails` above 0) or does not end with `OK: simulation completed` fails.
- **A failed trial rejects its step; it does not stop the fit.** A trial also fails when it runs
  past 3 times the median of the trials that succeeded, or `[fit].timeout` seconds if that is
  longer. So one pathological candidate cannot hold up an iteration.
- **States.** Every window has its own chain: a frozen run from the initial stand (the base
  config's census or state) that starts `[windows].chain_lead_days` (180) before the window and
  writes a state at its start. The window's soil water and temperature are then the model's own.
  The chains run at once. Before the polish they are re-run with the stage values
  (`[stages.polish].refresh`).
- **Windows, by rule** (`[windows]`; `select-windows` prints the choice). One ten-day calibration
  window per 1.5-month slot of the year inside the site TOML's leaf-on months, each the slot's best
  covered: the share of daytime records with H, LE and GPP measured, times the share with observed
  forcing. Validation windows are in the same slots of other years. A site may list its own.
- **Seasonal runs, by rule** (`[windows.seasonal]`). Up to two 120-day runs, each ending at a year's
  deepest cumulative water deficit (rain minus Priestley–Taylor evaporation, from the forcing)
  inside the leaf-on months, the deepest covered years first. They are scored on LE only. A year
  under 100 mm gives none; without any, the drought keys are fixed. At BCI the rule picks the 2016
  and 2017 dry seasons (1,331 and 698 mm).

### The tower

- **The site TOML holds the facts, `calibration.toml` the choices.** `[tower].site` names the
  tower's site TOML, the file the forcing build reads
  ([`../prepare_flux_tower`](../prepare_flux_tower/README.md)). It declares the file, the clock, the
  columns and their units, and in `[fluxes]` each target's column with the rule that says when it
  was measured (BCI: `FLAG = 1`). Both tools read it through one reader, `tower_inputs.read_standard`.
- **The tower's own interval, on UTC starts.** The observations stay on the tower's interval (30 min
  at BCI). The tool sets `[output].fast_interval_steps` so each trial writes one record per tower
  interval (2 × 900 s), and stops if the interval is not a whole number of `fast.dt_fast`. A trial
  whose output is at another spacing fails. The `hours` filter reads local hours, and the closure
  and EF local days, from the site's clock.
- **The metadata checks** (F1–F3): net radiation against its four components, the upwelling longwave
  above the downwelling at night (a swapped pair of columns), and the reflected shortwave and
  upwelling longwave within bounds. The forcing build stops on a failure; `report` and `fit` list
  it under `tower_checks`.

### Targets and filters

- **Targets:** the albedo, the upwelling longwave, H, LE, GPP and u\*. Net radiation (an input to the
  closure model), the evaporative fraction (it repeats H and LE) and night NEE are not targets.
  - Each target uses only the records the tower measured and the forcing observed.
- **Each target has an observation model** (`obsmodels.py`): how the tower's number relates to the
  true flux, each known bias as its own term.
  - **H and LE: closure.** The tower misses the same fraction of the turbulent flux at every hour.
    That fraction comes from whole days, f_d = median over ±15 days of the daily (Rnet − G)/(H + LE),
    using days with at least 70 % of records measured. The provider's gap-filled values fill those
    days' sums.
  - **The closure shares come from the attribution test:** within VPD classes, which of H/Rnet and
    LE/Rnet rises from the calmest to the most turbulent third of the records. At BCI, H/Rnet rises
    44–64 % in every class and LE/Rnet does not (−4 %, median), so H takes the whole gap. With
    f_d = 1.33 (daily closure 0.75) H rises to close the day and LE stays as measured.
  - **Bowen** (`[closure].shares = "bowen"`) is the declared alternative.
  - **GPP: respiration.** The provider made GPP from its own respiration, R_tower = κ R_true. The
    residual is [GPP_model − GPP_tower − (1/κ − 1) R_tower]/σ_NEE, and κ is an observation key: it
    is fitted with a prior, never written to a MEDS config, and its gradient costs no trials.
    - At BCI the prior is 0.65 ± 0.10: the tower's respiration is about the soil chambers' alone.
    - κ is never fitted beside a free shape key of the leaf's light response.
- **Random error,** σ = `sigma_abs` + `sigma_rel` |x|. The source is, in order: the provider's
  per-record uncertainty, the paired days (Hollinger & Richardson 2005), then the defaults.
  - σ is evaluated at a smoothed observation: the mean at the same time of day within ±7 days. GPP's
    σ is NEE's.
  - At BCI the paired days give LE 10.5 + 0.29|LE|, H 7.5 + 0.14|H| and NEE 2.1 + 0.17|NEE|.
- **Huber loss by default** (c = 2): half-hourly errors are heavy-tailed. χ² is reported on the
  plain residuals.
- **The model-error step:** at the refresh, each target's σ is multiplied by
  max(1, √(χ²/n)), at most 3×, so a target the model structurally cannot fit stops bending the keys.
  - H, LE and u\* are daytime targets: at night they are set mostly by the model's numerical floors
    and stable-air measurement problems (`night = true` keeps night records at the provider's u\*
    threshold). No hours are cut for transitions: canopy storage is ignored, every tower flux is
    treated as turbulent.
- **u\* per target, by its own diagnostic** (`datarules.py`, `[ustar]`). Within classes of the
  target's driver (PAR for GPP, Rnet for LE and H), the flux-to-driver ratio across u\* classes. The
  plateau test (Papale et al. 2006) finds the lowest class within 95 % of the mean of the classes
  above, and a bootstrap over days gives its spread (Barr et al. 2013).
  - A plateau is filtered at its threshold. A flat ratio or one still rising at the top classes
    is not filtered: a rising one is the observation model's to handle.
  - GPP takes the provider's CO₂ threshold (`ustar_min = "provider"`, one value or one per year).
  - At BCI: LE flat, H rising, so neither is filtered. GPP takes the provider's 0.4. The all-day
    diagnostic finds a GPP plateau at 0.33 (bootstrap 0.15–0.35); the morning hours alone gave 0.5.
- **Filters are per-target settings:** `ustar_min`, `par_min`, `hours`, `closure_range`,
  `min_solar_elevation` and `snow_free_days`, each off unless set. The albedo's defaults are the
  plan's rule: SW_in > 200 W m⁻², the sun at least 20° high, a week without frost.
- **Keys fixed from coverage.** A key that acts through one process (the registry's `process`: wet
  canopy, night, snow, drought) is fixed when the kept records sample that process fewer than
  `[fit].min_process_records` (48) times. At BCI the flag removes every rain half hour.
- **The data report** (`report`, and the start of every `fit`) gives:
  - every target's rows through each filter step;
  - the u\* diagnostics, and the rule each target got;
  - the windows, the years' water deficits and the seasonal runs;
  - the share of the record's daytime radiation, temperature, VPD and deficit range that the fit's
    records span;
  - at the fit, the share of the stand's area whose canopy-air top is above the sensor (#350).

### Keys and priors

- **Kinds and scopes** (the registry's `kind`, `scope`):
  - **kinds:** trait (measurable, prior from evidence), effective (a scheme property, labelled in
    the calibrated files), numerical (never calibrated), observation (κ);
  - **scopes:** plant type, site, observation.
- **The range is what is physically possible, the prior what the evidence says.** Every fitted
  key has its own sd. Its prior z at the MAP is reported, and gate G13 asks for a diagnosis of any
  trait key beyond 2 sd.
- **The priors from the site's climate** (`priors.py`; the forcing only, never the tower's fluxes):
  - **`stomatal_g1`:** the least-cost value, ξ = √(β (K + Γ\*)/(1.6 η\*)) with β = 146, at the
    growing-season daytime climate with MEDS's Rubisco kinetics. Log-sd 0.5. At BCI: 2.80
    (Lin et al.'s 3.77 for tropical rainforest is 0.9 sd away).
  - **`vcmax25`:** the coordination of the Rubisco- and light-limited rates, solved with MEDS's own
    leaf (`meds.plant.leaf`, so libmeds). Log-sd 0.5. At BCI: 41 with phi_psii 0.74 and θ_J 0.7; 94
    with the example's 0.85 and 0.9.
  - **Jmax/Vcmax, `ds_vcmax`, `ds_jmax`:** fixed at Kattge & Knorr for the growth temperature, set
    in every run and written into the calibrated files. At BCI (25.5 °C): 1.70, 641.1 and 640.6.
- **Priors by plant type** (`[fit].plant_type`, the registry's `meta`): the prior of a key without an
  EEO centre. Beside an EEO centre the type's prior is reported, and flagged when more than 2 sd
  apart. Site leaf data in `[priors]` come first.

- **The registry** (`parameters.toml`) is a menu. Each key has a state:
  - `fit`: in the default set;
  - `optional`: fitted when the site asks;
  - `fixed`: left at the base value, with the reason.

  Each key also has a stage, a range, a transform (logit on the range, log first for a positive
  scale) and a prior.
- **The site chooses the keys** in `[fit]`: `keys` (the exact list) or `add` / `remove` (edits of the
  default set). It sets any key's prior or range in `[priors.<key>]`.
- **Default priors** come from PFT-level syntheses where there is one, with the source given. For
  example, `stomatal_g1` is 3.77 kPa^0.5 with log-sd 0.35 (Lin et al. 2015, tropical rainforest
  trees). A key without a prior sd has the range as its ±2 sd band. The fit starts at the priors'
  centres.

### The staged fit

Within the fast loop the coupling is mostly one-directional: optics → light per leaf →
photosynthesis → stomata → energy balance → leaf temperature (a weak feedback). Each stage fits its
own keys with the cheapest model that can see them, keeping the upstream stages' values
(`stages.py`, `[fit].stages`):

| stage | model | target | method |
|---|---|---|---|
| `optics` | the two-stream alone (`meds.canopy`), over each window's stand and shortwave | albedo | Levenberg–Marquardt; seconds per evaluation |
| `photosynthesis` | a canopy of leaf solves (`meds.canopy`), over each window's per-cohort drivers | GPP | Levenberg–Marquardt; outer passes re-run the windows and refresh the drivers |
| `energy` | the coupled fast loop on the 10-day windows | every target | Levenberg–Marquardt |
| `water` | frozen seasonal runs (`[stages.water].windows`, e.g. 120 dry-season days) | LE, GPP | a grid in u, then a local quadratic (the responses can be rough) |
| `polish` | the coupled loop, every fitted key, from the stage values | every target (plus the seasonal runs) | Levenberg–Marquardt, up to `[stages.polish].max_iter` (10) iterations |

- **The default runs `energy`, `water` and `polish`.** A kernel stage (`optics`, `photosynthesis`)
  that is not listed has its keys fitted in `energy`, against every target. At BCI the albedo alone
  could not place the clumping: the `optics` stage sent it to its floor, and the joint fit then moved
  it to its ceiling.

- **The kernels are anchored** to the full run they take their drivers from. Their value for a
  parameter set is the kernel's, times, per record, the model's over the kernel's at the drivers'
  parameters. At the anchor they reproduce the model exactly; elsewhere the kernel carries the
  response.
- **Gate G8** checks the photosynthesis kernel before the anchor. With the full run's own
  parameters, the canopy of leaf solves must match that run's GPP within 1 % (median). A record
  averages the 15-minute steps in it (two at BCI), and the kernel solves the record's mean drivers
  once.
- **`["polish"]` alone** is the original joint fit of every key at once.
- **The stage values** are saved in `stages.json` after each stage. `--resume` continues from it,
  and `--stages` runs some stages only.

### Screening, weights and the cost

- **Screening** is the first Jacobian, at the start. It gives each key's sensitivity per target, its
  posterior-to-prior σ ratio, the collinear pairs, dead columns (a harness bug) and rough keys.
  - `[fit].screening = "report"` (the default) only reports it: every requested key is fitted.
  - `"drop"` fixes the keys the tower cannot inform, the dead ones and the rough ones.
- **Rough keys.** A rough key is one whose response over the finite-difference step is not smooth,
  so a Jacobian cannot steer it. With `"drop"` it stays at its default (`[fit].rough_keys = "default"`).
  `rough_keys = "line_search"` sets it by a 1-D search at the MAP instead; at BCI that search broke
  the five-year water budget (gate G7, #333).
- **Weights.** With `[fit].weights = "ess"` each target's rows count by its effective sample size,
  from the lag-1 autocorrelation of its residuals. The weights are set at the start, and again
  before the polish.
- **The cost** is the stacked weighted residuals `(model − obs)/σ` plus the Gaussian priors in u.
  `[fit].loss = "huber"` caps each row's pull at `huber_c` σ.

### Rough uncertainty

All of it comes from the final Jacobian, with almost no extra runs (revision plan §8):
- **The Laplace covariance at the MAP.** A target whose χ² per row is above 1 has its σ multiplied
  by √(χ² per row), so a target the model cannot fit does not make the keys look better known.
- **Intervals** (68 % and 95 %) mapped back through the transform: asymmetric, and inside the bounds.
- **The linearity check** along the leading directions: the objective at ±1 sd against the quadratic.
  The mean of the two sides is the actual curvature, and it marks the covariance "local only" where
  that is outside 0.5–2× the quadratic's. Half their difference is the slope left, 0 at a converged fit.
- **Convergence:** the Gauss–Newton step left at the MAP, per key in posterior sd. A fit stopped by
  `max_iter` while still descending shows it here.
- **The filter sensitivity.** The MAP's shift under `[uncertainty].alternative` (e.g. GPP u\* ≥ 0.3)
  is the Gauss–Newton step of the alternative rows minus that of the fit's own rows, both from the
  MAP (gate G10).

### Validation and gates

The default and the MAP are scored on windows the fit never saw, each from its own chain's states.
The default is the base configuration, on chains run with it.

| gate | what |
|---|---|
| G3 | screening: no dead Jacobian column (reported) |
| G4 | the MAP beats the default on validation, no target more than 10 % worse |
| G5, G11 | keys near a bound, each with the target that pushed it and its stage |
| G8 | the photosynthesis kernel against the model's GPP |
| G10 | the filter sensitivity: every key's predicted shift below 1 posterior sd |

## Commands

```
calibrate_fast.py select-windows --site calibration.toml                 # print the rule's windows and seasonal runs
calibrate_fast.py report --site calibration.toml --work runs/report      # the data report only
calibrate_fast.py check --site calibration.toml --variant interception_on --work runs/check \
    --workers 8                                                         # runs, record, G1, G2
calibrate_fast.py fit   --site calibration.toml --variant interception_on --work runs/on --workers 40
calibrate_fast.py fit   ... --stages optics,photosynthesis              # some stages only
calibrate_fast.py fit   ... --resume                                    # continue from stages.json
calibrate_fast.py analyze --site calibration.toml --variant interception_on --work runs/on \
    --workers 40                                                        # redo the post-fit steps
calibrate_fast.py write-calibrated --site calibration.toml --variant interception_on \
    --fit runs/on/fit.json --out calibration                            # the configs with the MAP
calibrate_fast.py smoke --site ... --work ... --g8                      # the ctest smoke test
```

Add `--runner ../../build-ifx/meds_main` to `check`, `fit` or `analyze` to run the trials with the
executable. The optics and photosynthesis stages need `libmeds.so` in any case.

`analyze` reruns everything after the stages from the MAP in `fit.json`, reusing the cached trials:
the covariance, the intervals, the linearity check, the filter sensitivity, the validation and the
gates. Use it after changing the validation windows, a target's σ or `[uncertainty]`.

`--pool local` runs `--workers` trials at once on this machine. On a Slurm cluster, `--pool queue`
uses a directory queue under `--work`. Start one `calibrate_fast.py worker --queue <work>/queue
--slots <cores>` per node inside the same allocation, and the driver hands trials to them. So one
allocation holds the workers for the whole fit, and no trial waits in the scheduler:

```
srun --ntasks-per-node=1 --cpus-per-task=40 python calibrate_fast.py worker --queue runs/on/queue --slots 40 &
python calibrate_fast.py fit ... --pool queue
touch runs/on/queue/STOP      # the driver writes it too, when the fit ends
```

The kernel stages run in the driver's process, on the first node.

Ask Slurm for the memory with `--mem=0` (the whole node). A trial takes up to ~1 GB, and a cluster
whose default is 1 GB per job kills the workers. A node running a trial on every core runs each
trial about 2.5× slower than an idle node does: 16 s against 6.4 s at BCI.

## What `fit` writes (in `--work`)

| file | what |
|---|---|
| `fit.json` | the data report, the screening, each stage's result, the MAP, start and default values, the intervals, the covariance and correlations in θ, each key's posterior-to-prior σ ratio, the σ scales, the linearity check, the filter sensitivity, the scores per target on the calibration and validation windows, the gates, and the trial count and timing |
| `stages.json` | the values after each stage (for `--resume`) |
| `pft_parameters_calibrated.toml`, `meds_config_calibrated.toml` | the base configs with the MAP written in, comments kept |
| `fit.log` | the run's log |
| `trials/`, `chains/`, `kernel_configs/` | every trial, state and kernel configuration (cached: a re-run reuses them) |

## Files

| file | role |
|---|---|
| `calibrate_fast.py` | the commands, the stages' driver, the post-fit analysis |
| `site_reference.toml` | every site setting, its default and why (the schema) |
| `settings.py` | reads a site declaration against `site_reference.toml` |
| `parameters.toml` | the registry: each key's state, stage, file, TOML key, range, transform, prior and source |
| `registry.py` | the transforms, the priors, the key selection |
| `stages.py` | the optics and photosynthesis kernel models, gate G8, the water stage's grid search |
| `trials.py` | the trial writer and runner, the parameter-record check, the fast-series reader (on the tower's interval), the driver trials |
| `states.py` | the state chains |
| `tower.py` | the tower's records (through `tower_inputs`), the closure correction, the observed-forcing mask, the sun's elevation |
| `datarules.py` | the u\* diagnostic, the process coverage, the window rule, the water-deficit index and seasonal runs |
| `priors.py` | the growth climate, the EEO centres (least-cost g1, coordination vcmax25 with MEDS's leaf), Kattge & Knorr |
| `obsmodels.py` | the closure model and its attribution test, the respiration model's κ prior, the random error (provider, paired days) at a smoothed observation |
| `residuals.py` | the targets, the filters, the residual vector, the weights, the filter report |
| `fit.py` | Levenberg–Marquardt, the Jacobian, screening, the covariance and intervals, the linearity check, the filter shift |
| `pool.py` | the local pool and the directory queue |
| `tests/` | unit tests, and a smoke test that runs the model and checks G8 (ctest `calibrate_fast`) |
