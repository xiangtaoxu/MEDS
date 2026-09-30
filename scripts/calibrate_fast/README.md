# calibrate_fast — fast-parameter calibration against a flux tower

Fits MEDS's sub-daily parameters (radiation, photosynthesis and stomata, aerodynamics, water stress,
respiration) to a flux tower's fluxes, with the stand frozen at its initial structure. It also yields
the parameters' covariance, the starting point for a longer-term calibration. The design, and every
decision in it, is [`docs/dev_plans/MEDS_FAST_CALIBRATION_PLAN.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_PLAN.md).
Barro Colorado Island is the worked example: [`examples/example_flux_tower_bci/calibration.toml`](../../examples/example_flux_tower_bci/calibration.toml).

## How it works

- **Trials.** A trial is a 10-day frozen run (`[run].slow_on = false`) restarted from a shared state at
  its window's start, with `[init].reacclimate_traits = true` so a changed trait reaches the cohorts.
  Its configs are made by parsing the base TOMLs and setting keys, and its directory is named by a
  hash of them, so a repeated candidate reuses the finished run.
- **Every trial proves what it ran.** After the run, the trial's parameter record
  (`<prefix>_parameters.csv`, written by `meds_main`) must list every key the trial set, marked as set
  in the file, with the value written. A missing or defaulted key stops the fit instead of silently
  running the default. A trial that breaches a whole-site budget (a `budget[whole_*]` line with
  `fails` above 0) or does not end with `OK: simulation completed` fails.
- **A failed trial rejects its step; it does not stop the fit.** A trial also fails when it runs
  past 3 times the median of the trials that succeeded, or `[fit].timeout` seconds if that is
  longer. So one pathological candidate cannot hold up an iteration.
- **States.** A chain is a frozen run from the census that stops at each window's start in turn and
  writes a state there. The windows' soil water and temperature are then a spun-up run's. The
  chains are re-run with the MAP set once the fit has one.
- **Targets.** The targets are albedo, upwelling longwave, net radiation, LE, H (both closure-corrected
  with the Bowen ratio kept), the daily evaporative fraction, daytime GPP, night NEE plus the stand's
  growth respiration, and u\*. Each uses only the hours the tower measured and the forcing observed,
  and each has an observation error σ (`tower.py`, `residuals.py`).
- **The fit.** Levenberg–Marquardt on the stacked residuals `(model − obs)/σ` of every target in
  every window, plus Gaussian priors in a transformed space (logit on the range, log first for a
  positive scale). The Jacobian comes from central differences, and all 2k·W trials of a Jacobian run
  at once. Each iteration tries three damping values at once. There are three starts, and a failed
  or timed-out trial rejects its step. `fit.py` has the details.
- **Screening.** The first Jacobian, at the default, gives each key's sensitivity per target, its
  posterior-to-prior σ ratio, the collinear pairs, dead columns (a harness bug) and rough keys. Keys
  the tower cannot inform are fixed at their default.
- **Rough keys stay at their default** (`[fit].rough_keys = "default"`). A rough key is one whose
  response over the finite-difference step is not smooth, so a Jacobian cannot steer it. At BCI
  these were `wood_psi50` and `leaf_pi0`. `rough_keys = "line_search"` sets each by a 1-D search at
  the MAP instead. At BCI that search found values that fit the 10-day windows better but broke the
  five-year run's water budget (gate G7, #333), which is why it is not the default. Both BCI keys are
  rough because they move the threshold of MEDS's whole-day stomatal latch (#332).
- **The covariance.** At the MAP, the Laplace covariance weights each target by its effective sample
  size, from the lag-1 autocorrelation of its residuals. A linearity check along the three leading
  directions compares the actual change in the objective at ±1σ with the quadratic prediction.
- **Validation.** The default and the MAP are scored on windows the fit never saw, each from its own
  chain's states. The gates of plan §8 follow.

## Commands

```
calibrate_fast.py select-windows --site calibration.toml                 # the best-covered windows
calibrate_fast.py growth-resp --daily "output/eval-D-*.nc" --out calibration/growth_resp_monthly.csv
calibrate_fast.py check --site calibration.toml --variant interception_off --work runs/check \
    --meds-main ../../build-ifx/meds_main --workers 8                   # runs, record, G1, G2
calibrate_fast.py fit   --site calibration.toml --variant interception_off --work runs/off \
    --meds-main ../../build-ifx/meds_main --workers 40
calibrate_fast.py analyze --site calibration.toml --variant interception_off --work runs/off \
    --meds-main ../../build-ifx/meds_main --workers 40                  # redo the post-fit steps
calibrate_fast.py write-calibrated --site calibration.toml --variant interception_off \
    --fit runs/off/fit.json --out calibration                           # the configs with the MAP
```

`analyze` reruns everything after the iterations from the MAP in `fit.json`, reusing the cached
trials: the rough keys, the covariance, the linearity check, the validation and the gates. Use it
after changing `[fit].rough_keys`, the validation windows or a target's σ.

`--pool local` runs `--workers` trials at once on this machine. On a Slurm cluster, `--pool queue`
uses a directory queue under `--work`. Start one `calibrate_fast.py worker --queue <work>/queue
--slots <cores>` per node inside the same allocation, and the driver hands trials to them. So one
allocation holds the workers for the whole fit, and no trial waits in the scheduler:

```
srun --ntasks-per-node=1 --cpus-per-task=40 python calibrate_fast.py worker --queue runs/off/queue --slots 40 &
python calibrate_fast.py fit ... --pool queue
touch runs/off/queue/STOP      # the driver writes it too, when the fit ends
```

Ask Slurm for the memory with `--mem=0` (the whole node). A trial takes up to ~1 GB, and a cluster
whose default is 1 GB per job kills the workers. A node running a trial on every core runs each
trial about 2× slower than an idle node does: 18.6 s against 8.5 s at BCI.

## What `fit` writes (in `--work`)

| file | what |
|---|---|
| `fit.json` | the screening, every start's path, the MAP and default values, the covariance and correlations in θ, each key's posterior-to-prior σ ratio, the linearity check, the scores per target on the calibration and validation windows, the gates, and the trial count and timing |
| `pft_parameters_calibrated.toml`, `meds_config_calibrated.toml` | the base configs with the MAP written in, comments kept |
| `fit.log` | the run's log |
| `trials/`, `chains/` | every trial and state (cached: a re-run reuses them) |

## Files

| file | role |
|---|---|
| `calibrate_fast.py` | the commands |
| `parameters.toml` | the registry: each key's file, TOML key, range, transform and source |
| `registry.py` | the transforms and the prior |
| `trials.py` | the trial writer and runner, the parameter-record check, the hourly reader |
| `states.py` | the state chains |
| `tower.py` | the tower's hours, the closure correction, the observed-forcing mask |
| `residuals.py` | the targets and the residual vector, the effective-sample-size weights |
| `fit.py` | Levenberg–Marquardt, the Jacobian, screening, the covariance, the linearity check |
| `pool.py` | the local pool and the directory queue |
| `tomlio.py` | TOML in and out |
| `tests/` | unit tests, and a smoke test that runs `meds_main` (ctest `calibrate_fast`) |
