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
- **States.** A chain is a frozen run from the census that stops at each window's start in turn and
  writes a state there. The windows' soil water and temperature are then a spun-up run's. Before the
  polish the chains are re-run with the stage values (`[stages.polish].refresh`).

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

- **Targets.** On by default: the albedo, upwelling longwave, LE, H, daytime GPP and u\*.
  - LE and H are as measured (`[tower].closure = "none"`). At BCI the closure gap behaves like
    missing sensible heat, so the Bowen-ratio correction (`"bowen"`) would raise LE in a way the
    data do not show.
  - Off by default, each with its reason in `site_reference.toml`:
    - net radiation: from a four-component radiometer it repeats the albedo and the upwelling
      longwave;
    - the daily evaporative fraction: biased high where closure is poor;
    - night NEE.
  - Each target uses only the records the tower measured and the forcing observed. Its observation
    error is σ = `sigma_abs` + `sigma_rel` |obs|.
- **Filters are per-target settings:** `ustar_min`, `par_min`, `hours`, `closure_range` and
  `min_solar_elevation`, each off unless set. The recommended values, with the reasons, are in
  `site_reference.toml`. Every target takes every filter. The defaults:
  - GPP: u\* ≥ 0.4 m s⁻¹ and σ = 2.5 + 0.15 GPP. The tower's GPP is built from one respiration
    value per day, so it carries a systematic error at every daytime hour, on top of the random one.
  - LE: u\* ≥ 0.4 m s⁻¹.
  - H: u\* ≥ 0.6 m s⁻¹, 9 to 16 h, and σ = 10 W m⁻² + 30 %. H is under-measured in low
    turbulence, and canopy heat storage holds it back in the morning.
- **The data report** (`report`, and the start of every `fit`) gives:
  - every target's rows through each filter step;
  - the GPP diurnal mean before and after the filters;
  - the u\* plateau test on the morning GPP/PAR (Reichstein et al. 2005, within PAR classes).

### Keys and priors

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
calibrate_fast.py select-windows --site calibration.toml                 # the best-covered windows
calibrate_fast.py growth-resp --daily "output/eval-D-*.nc" --out calibration/growth_resp_monthly.csv
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
| `residuals.py` | the targets, the filters, the residual vector, the weights, the data report, the u\* plateau |
| `fit.py` | Levenberg–Marquardt, the Jacobian, screening, the covariance and intervals, the linearity check, the filter shift |
| `pool.py` | the local pool and the directory queue |
| `tests/` | unit tests, and a smoke test that runs the model and checks G8 (ctest `calibrate_fast`) |
