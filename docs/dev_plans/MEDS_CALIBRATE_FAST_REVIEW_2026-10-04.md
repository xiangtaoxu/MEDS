# calibrate_fast: a clean-up review (2026-10-04)

**Status:** executed the same day, on branch `cleanup/calibrate-fast`. The owner took the recommended choice in each of D1–D4 (§6).

The review covers `scripts/calibrate_fast` at beta `d8a82b5`, after the best-practice protocol was merged. That is 12 modules (3,972 lines), the two TOMLs, the README and the tests, with the BCI example's `calibration.toml` and `run_example.py`. It answers the owner's four questions:

1. Which gates before and after the calibration are unnecessary, and which were set up for this site? Would it be simpler to let the tool fail with an error?
2. Which code sections are unnecessary? Would an error message be simpler?
3. Are the names readable for an ecologist who knows flux towers and forests but not necessarily Python?
4. What else can be consolidated, and which special cases can go?

Decisions D1–D4 (§6) are the owner's. Everything else is a recommendation that follows from the code principles: readable and modular, plain-language comments, fewer special cases, no fallbacks.

## Summary

- **Gates:** none of the nine stays a gate.
  - **Two are real checks.** G4 (validation) stays as the report's one pass/fail line. G13 (prior z) becomes a mark in the keys table.
  - **Five are pass/fail labels that cannot fail, or that report nothing:** G3, G5, G7, G10, G12.
  - **G1 and G2 test the tool, not the fit.** They move into the CTest check.
- **Let it fail:**
  - **A failed trial stops the fit (D1)**, instead of being absorbed by retries and rejected steps.
  - **Five places where data are dropped silently become errors.** For example, a misspelt forcing qc name currently switches the observed-forcing mask off.
- **Remove:**
  - four commands: `analyze`, `write-calibrated`, `smoke` (merged into `check`), and `select-windows` (folded into `report`);
  - eight options the protocol has already decided;
  - `[tower].forcing` and `forcing_grid`, which the base config already holds;
  - one BCI-only diagnostic;
  - nine unused helpers and five dead branches.
- **Rename:**
  - `calibration.toml` stops being "the site declaration", because "site" already names the tower's site TOML;
  - module aliases (`DR`, `F`, `OM` …) become full names;
  - the seasonal runs get one name;
  - about a dozen settings and code terms get plain names.
- **Consolidate:**
  - one lookup for "the base config's value, else the parameter record's";
  - one helper for each target's rows;
  - one "which target pushes this key" function;
  - the tower spec folded into the calibration object;
  - no swapping of the tool's state in a `try/finally`.
- **Size:** roughly 500 of the 3,972 lines go (about 13 %, a rough estimate), and about 60 lines of `site_reference.toml`. BCI's fit is unchanged, except under D3.

## 1. Gates and checks

### Before the fit

| check | where | now | recommendation |
|---|---|---|---|
| unknown or missing settings | `settings.complete` | stops the tool | **keep**: it is what makes a misspelt key an error |
| `[base].pft` refused | `Site.__init__` | stops the tool | **remove**: the schema already refuses any key not in `site_reference.toml` |
| unknown `--variant` | `Site.__init__` | stops | keep |
| tower interval a whole number of fast steps | `Site.__init__` | stops | keep |
| tower flux checks F1–F3 (Rnet's parts, longwave at night, bounds) | `tower_inputs.read_standard(strict=False)` | a `WARNING` line and `tower_checks` in the report | **read strictly**: the forcing build already stops on the same site TOML, so the calibration should behave the same way |
| registry entries | `registry._check` | stops | keep, trimmed with the fields removed (§3) |
| a prior centre outside the range | `registry.resolve_defaults` | stops | keep |
| κ beside a shape key | `registry.select` | stops | keep |
| no calibration window passes the rule | `load_data` | stops | keep |
| G1, G2 (stand unchanged; byte-identical repeat) | the `check` command | a separate command the user must run | **move into the CTest check** (§3): they test the tool and the model's restart, not a site's fit |

### During the fit

The trial checks in `trials.finish` stay. A trial must:
- end with `OK`;
- show no whole-site budget breach;
- show in its parameter record every key the trial set;
- have no NaN in its output (G12);
- write its output at the tower's interval.

These checks found #352 and #363.

The timeout is now three times the median trial time, with a floor at `[fit].timeout`. It becomes simply `[fit].timeout` once a failure stops the fit (D1).

### After the fit

| gate | what the code does | recommendation |
|---|---|---|
| G3 | passes when no fitted key is dead or rough. In triage mode the triage has already fixed those keys, yet the gate still fails whenever a dead key existed (an inconsistency) | **remove**: the triage's report lists what it fixed and why |
| G4 | MAP beats the default on validation; no target more than 10 % worse | **keep as the report's one verdict** ("validation: passed / failed, H +4 %") |
| G5 | `pass = True` always; it lists the keys near a bound | **a "near a bound, pushed by H" mark in the keys table** |
| G7 | `pass = None` always; a reminder to run the five years | **remove from the tool**; a step in the README and the example |
| G10 | passes when every alternative moved less than 1 sd or was refitted. The tool refits whenever one moves more, so it always passes | **remove**: the alternatives section reports each shift and refit |
| G12 | `pass = True` always; failed trials never score | **remove**: `trials.finish` enforces it |
| G13 | flags trait keys with \|prior z\| > 2 | **a mark in the keys table** ("beyond 2 sd: diagnose (§5.1) or relabel effective") |

The report's "Gates" section becomes "Validation": the default against the MAP per target, and the verdict. README §"Validation and gates" and plan §7.2 follow.

### Set up for this site

- **`[tower].forcing_qc_val`** exists because BCI's tower measured longwave only from 2015. A second qc list lets the validation windows sit in 2012–2014, where the forcing's longwave is synthesized. That is a special case, and its validation partly scores the longwave fill (D3).
- **The registry's `fixed` states that record BCI's screening:**
  - VIS optics: "the tower cannot inform it";
  - NIR transmittance and `d_ratio`: collinear;
  - `leaf_width` and `dsl_dmax`: the BCI fits drove them to a bound;
  - `intercept_k`: rough.

  At another tower the triage and the coverage rule decide these per site (D2).
- **BCI's numbers in the generic files.** `site_reference.toml`, the tool's README and several docstrings quote BCI results ("H/Rnet rises 44–64 %", "κ 0.65", "LE 10.5 + 0.29|LE|"). The rules belong there, and BCI's outcome belongs in the example's README.
- **`area_above_sensor`** (#350) is a diagnostic of BCI's forcing height against its canopy-air tops. It changes nothing in the fit. **Remove it from the tool**; #350 has the analysis.

## 2. Let it fail

### A failed trial (D1)

Today a failed trial is absorbed:
- `Model.run` and `Model.residuals` return `None` for it;
- `Problem.data` passes the `None` on;
- `jacobian` retries a column at h/2, and `jacobian_one_sided` retries it backward;
- `lm` gives the candidate an infinite cost;
- the triage calls a key with a failed column "rough";
- `linearity` and `alternatives` each have a "a trial failed" branch;
- `n_failed` is counted and reported.

At BCI v2 no trial failed, and the failures in earlier fits were model bugs (#352, #363), each better seen at once. **Recommendation: a failed trial stops the fit.** The error names the trial directory, its parameter values and the last lines of its log. The `None` plumbing and both retry loops go (about 60 lines).

What it costs: a 25-minute fit stops at a candidate the model cannot run. The steps are capped at 2 in u, so that candidate is within the keys' ranges. A model that fails there has a bug to fix.

### Silent drops that should be errors

| where | now | recommendation |
|---|---|---|
| `obsmodels.closure_factor` without Rnet, H or LE | f_d is NaN, so every H and LE row drops and the targets vanish | error: "the closure model needs Rnet, H and LE in the site TOML's [fluxes]" |
| the respiration model without RECO | every GPP row drops | error: "the GPP observation model needs RECO" |
| `tower.forcing_observed`, a listed qc variable missing from the forcing file | skipped (`if v in ds.variables`): a misspelt qc name switches the mask off | error naming the variable |
| `specs_for`, a window with no usable records | `WARNING` | error: the window rule, or the site's list, chose a window with no data |
| `datarules.resolve_ustar`, GPP's diagnostic with too few records | GPP alone falls back to the provider's threshold | the same as the other targets (no filter, reported), or an error. The default rule is `"provider"` anyway |

### Fallbacks that stay

These are rule outcomes, not errors, and each is reported:
- σ's source, in order (the provider, the paired days, the defaults);
- the attribution test with too few records, giving "as measured";
- the u\* diagnostic's outcomes;
- no year deep enough for a seasonal run, so the drought keys are fixed.

## 3. Unnecessary code

### Never used

| what | where |
|---|---|
| `to_datetime` | `datarules.py` |
| `local_days` | `tower.py` |
| `prior_u` | `registry.py` |
| `clean` | `trials.py` |
| `WindowSpec.slices` | `residuals.py` |
| `Window.chain` (always the window's name) | `trials.py` |
| `Param.group`, `notes`, `extra` (read from the registry, never used) | `registry.py` |
| `specs_for(targets=…)` (never passed) | `calibrate_fast.py` |
| `NC_LOCK`: no thread reads netCDF since the multiple starts went | `trials.py` |

### Branches never taken

| branch | why |
|---|---|
| `elev is None` (3 places), `if site.lat is not None` | `tower_inputs` requires the site's latitude and longitude |
| `days_since_frost` without Tair; "snow_free_days needs Tair" | Tair is a required variable of the site TOML |
| no forcing file: `forcing_span` None, `forcing_observed` all True, no deficit, no seasonal runs, no climate priors | the forcing is the base config's `[forcing].path` and is always there (below) |
| `fit.row_scales`: the ESS and σ-scale factors | never taken with `weights = "ess"` and `refresh = true`, so the covariance is plainly (J'J + P)⁻¹ |
| `states.chain_segments` and the loop over a chain's windows | every chain now has one window |

### Documented but missing

`site_reference.toml` documents a `closure_range` filter that no code implements. Remove the lines.

### Commands

| command | recommendation |
|---|---|
| `analyze` | **remove.** Re-running `fit` replays every trial from the cache, and `analyze` repeats `fit`'s setup (40 lines) |
| `write-calibrated` | **remove.** `fit` already writes both calibrated configs into `--work`; `run_example.py` copies them as it copies `fit.json` |
| `check` and `smoke` | **merge into one `check`**, which CTest runs. It runs a chain, a trial with its record checked, one gradient column that must move the output, a byte-identical repeat (G2) and an unchanged stand (G1) |
| `select-windows` | **fold into `report`**, which already shows the windows and can print them as `[[windows.list]]` entries |
| `report`, `fit`, `variants`, `worker` | keep |

### Options the protocol has already decided

Each is a second code path, and the BCI fit uses one side only.

| option | recommendation |
|---|---|
| `[fit].weights` ("ess" / "none") | ESS always |
| `[fit].refresh` (true / false) | always |
| `[fit].screening` ("triage" / "report") | triage always |
| `[fit].loss` ("huber" / "l2") | Huber always; `huber_c` stays |
| `[targets.*].obs_model` | remove. For H and LE, `[closure].shares = "none"` already means "as measured". For GPP, taking κ out of the fit (`[fit].remove = ["kappa"]`) leaves κ = 1, which is "as given" |
| `[targets.*].sigma_source` | remove: the protocol's order always |
| `[targets.*].night` and its `night_ustar` | remove: H, LE and u\* are daytime targets (owner decision, 2026-10-03) |
| `sigma` as an alias of `sigma_abs` | one name, `sigma_abs` |
| `[tower].forcing`, `[tower].forcing_grid` | read the base config's `[forcing].path` and `grid_index`, so a calibration cannot point at another forcing than its runs use |

### The declared alternatives' reverse cases

`alternative_data` also handles two cases BCI never reaches:
- a fit that already used Bowen, whose alternative is the attribution's shares;
- a fit that used GPP's diagnostic, whose alternative is the provider's threshold.

With the defaults, each alternative is one direction. **Skip the alternative when the fit already used it** (about 15 lines).

### Kept, though BCI never exercises them

These are protocol rules for FLUXNET towers, and they are unit-tested:
- κ's sd from the provider's two partitionings (RECO against RECO_DT);
- the partitioning alternative;
- σ from the provider's random uncertainty.

## 4. Names

### "Site" means two files

| file | now called | proposed |
|---|---|---|
| `calibration.toml` | "the site declaration", `--site`, class `Site`, `site_reference.toml` | "the calibration settings": `--config`, class `Calibration`, `calibration_reference.toml` |
| `bci_site.toml` | "the site TOML", `[tower].site` | unchanged |

### Module aliases

`import datarules as DR`, `fit as F`, `obsmodels as OM`, `priors as PR`, `residuals as R`, `settings as SET`, `states as S`, `tower as TW` and `trials as T` make a reader look up a letter at every call (`DR.select_windows`, `F.lm`, `T.build_trial`). **Import them by their names.**

### Files

| now | proposed | why |
|---|---|---|
| `datarules.py` | `data_rules.py` | readable |
| `obsmodels.py` | `observation_models.py` | readable |
| `residuals.py` | `targets.py` | it defines the targets, their filters and their residuals |
| `states.py` | `chains.py` | it runs the state chains |
| `pool.py` | `workers.py` | the local and Slurm workers |
| `registry.py` | `parameters.py` | it reads `parameters.toml` |
| `fit.py`, `trials.py`, `tower.py`, `priors.py`, `settings.py` | unchanged | `fit.Model`, which runs trials, moves to `trials.py` as `TrialRunner` ("model" means MEDS) |

### Code terms

`theta`, `u` and `Phi` stay in `fit.py`'s mathematics, with the module docstring saying what they are. Elsewhere:

| now | proposed |
|---|---|
| `theta0`, `theta_base`, `theta_map` | `start`, `default`, `map` (the report's own words) |
| `WindowSpec`, `spec`, `fspecs` | `WindowRows`, `rows`, `fit_rows` |
| `fok` | `forcing_observed` |
| `elev` | `sun_elevation` |
| `obs_keys`, `obs_fixed` | `observation_keys` |
| `derived`, `derived_was` | `kattge_knorr`, `kattge_knorr_replaced` |
| `menu()` | `registry()` |
| `Problem` | `FreeKeys` |
| `water_windows`, role `"water"`, `[windows.seasonal]`, `dry2016` | one name, **seasonal runs**: `seasonal_runs` and role `"seasonal"` |
| "Jacobian" and "gradient matrix" (both used) | "gradient matrix" in prose, `J` in the mathematics |
| `cand`, `J0`, `sm0`, `fl0`, `rv0`, `rv1`, `kr`, `cv`, `ccfg`, `scfg`, `wcfg`, `fc`, `uc`, `tw` | spelled out |

### Settings

| now | proposed |
|---|---|
| `[ustar].crit` | `plateau_fraction` |
| `[ustar].edges` | `classes` |
| `[closure].rise_min` | `min_rise` |
| `[windows].slots` | `per_year` |
| `[windows].min_score` | `min_coverage` |
| `[fit].rtol` | `min_cost_drop` |
| `[fit].jacobian_refresh` | `recompute_every` |
| `[fit].informed_ratio` | `uninformed_sd_ratio` |
| `[fit].max_corr` | `max_correlation` |
| `[uncertainty].refit_sd` | `refit_beyond_sd` |

Only BCI's `calibration.toml` uses the tool, so now is the cheap time to rename (D4).

## 5. Consolidation

- **One lookup:** "the base config's value, else the base run's parameter record's". It is now written four times: `Site.setting`, the `value` helper in `area_above_sensor`, the `value` lambda in `climate_priors`, and `registry.record_value`.
- **One helper for each target's rows** in the stacked residuals. The same loop is written in `fit.triage`, `fit.bound_pushers` and `fit.prior_z`.
- **One "which target pushes this key" function.** `bound_pushers` and `prior_z` both compute each target's share of a key's gradient. It feeds the keys table's prior z, bound and push columns.
- **`TowerSpec` folds into the calibration object.** It copies `tower_inputs.Site`'s fields (clock, interval, coordinates, leaf-on months, provider, sensor height), and its defaults repeat `site_reference.toml`'s.
- **No hidden state:**
  - `load_data` rewrites `site.targets`;
  - `site.data` holds the observations, the masks, the record and the σ scales;
  - `alternative_specs` swaps `site.targets` and `site.data["obs"]` inside a `try/finally`.

  Instead, `load_data` returns the data, and `specs_for` takes the targets and observations it builds from.
- **`post_fit`'s 12 positional arguments** become one small object for the fit's state.
- **The data report is written three times:** log lines, `data_report.json` and `fit.json`. Write it once into `fit.json` and a "Data" section of `report.md`; the log keeps one line per step.
- **Stale text:**
  - references to the superseded `MEDS_FAST_CALIBRATION_PLAN.md` and `…_REVISION_PLAN.md` in five modules and three TOMLs;
  - "the fit's parallel starts" (`pool.py`, `trials.py`);
  - "the hydraulic cost cliff" (`fit.Model._timeout`);
  - `states.py`'s chain through several windows;
  - the README's "the separate water stage … are gone" (history, for the CHANGELOG).

## 6. Decisions

- **D1. A failed trial stops the fit** (§2). Recommended: yes.
- **D2. The registry's states that record BCI's screening** (VIS optics, NIR transmittance, `d_ratio`, `leaf_width`, `dsl_dmax`, `intercept_k`). Recommended: make them `optional`, with reasons stated as rules rather than BCI results. BCI's fit does not change, since they stay out of the default set; a site that adds them gets the triage. The other choice is to keep them fixed as now.
- **D3. One forcing-qc list for every window** (drop `forcing_qc_val`). Recommended: yes. The validation then scores only observed forcing. BCI's validation windows move into 2015–2017, the years with observed longwave, which takes a refit of both variants (~25 min each).
- **D4. The renames' scope** (§4). Recommended: all of them (files, flags, settings, code terms), with BCI's `calibration.toml` updated in the same change.

## 7. Order of work

1. Remove the dead code and the gates, and turn the silent drops into errors (no behaviour change at BCI).
2. Remove the options and commands, read the forcing from the base config, and merge `check` with `smoke`.
3. Make the D1–D3 changes.
4. Rename (D4), then consolidate (§5).
5. Update the README, `site_reference.toml` (rules only) and plan §7.2, and move BCI's numbers to the example's README.
6. Run the tests (unit tests, the CTest check, the end-to-end fit) and refit BCI where D3 asks for it. A re-run of the shipped fit must reproduce its MAP from the trial cache.
