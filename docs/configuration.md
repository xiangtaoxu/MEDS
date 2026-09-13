# Configuring a MEDS run

A run is described by **two [TOML](https://toml.io) files**, and the first is the only command-line
argument:

```bash
./build/meds_main path/to/run_main.toml
```

With no argument, `meds_main` reads `./meds_config_main.toml`. Relative paths inside the config
resolve against the directory you launch from, not the config's own location.

## The parameter philosophy

**MEDS hard-codes no model parameters.** The source defines only true *constants* — numerical,
geometric and calendar. Every parameter is user-mutable in TOML and **required**: a missing
parameter, or a missing file, is a **hard error**. The reader builds a presence map while loading
and aborts listing every absent key, rather than falling back to a default nobody chose.

This is deliberate and it is not free — a MEDS config is long. The reason is that a silent default
is indistinguishable from a considered value once a run is a year old, and the model has already
been bitten by exactly that: a flag whose loader default disagreed with its type default ran an
entire subsystem off for every config that omitted the key.

Two exceptions, both explicit:

- **Derived quantities** are computed from the primary ones — the step in years, the height class
  edges, the mortality-hazard coefficients from wood density. Set `[options].override_derived =
  true` and supply a `[derived]` block to pin them instead.
- **Defaulted blocks.** A few late-added blocks (`[soil_column]`, `[snow]`, `[soil_carbon]`) read
  each key with a fallback to its in-type default, so the feature needs no TOML edits. Where this
  applies, the block's comment says so.

Every run writes `<output_dir>/<prefix>_pft_parameters.csv` — one row per PFT with every per-PFT
trait actually used. It is a provenance record: what the run really did, not what you meant.

## The two files

### The main file

All non-PFT settings. Named on the command line; it names the PFT file via `[init].pft_config`.

| Block | What it sets |
|---|---|
| `[run]` | The slow timestep, and the run span as **calendar dates** (`start_time`, `end_time`, leap-year-aware Gregorian, `"YYYY-MM-DD[ HH:MM:SS]"`). Thread count. |
| `[fast]` | The sub-daily loop: `dt_fast`, the integrator, tolerances, error control. |
| `[init]` | How the run starts, and the path to the PFT file. |
| `[demography]` | Cohort and patch fusion/fission, the cadence switches. |
| `[disturbance]`, `[recruitment]` | Treefall rates; the seed-rain and reproduction settings. |
| `[soil_column]` | The **physical ground**: layer count, depth, grid growth, hydraulic texture, retention family, root profile, thermal properties. |
| `[soil]`, `[energy]`, `[snow]`, `[aerodynamics]` | The **solvers over it**: selectors and tolerances. |
| `[soil_carbon]` | The CENTURY decomposition: selectors, rate parameters, cold-start spin-up. |
| `[hydraulics]` | Plant water transport. |
| `[forcing]`, `[site]` | The meteorological driver, and where the site is. |
| `[output]` | Which diagnostics are written, on which axes, at which timescales. |
| `[io]` | Output directory and prefix, and restart checkpointing. |
| `[options]` | `override_derived` and other run switches. |

**`[soil_column]` is the ground; `[soil]` is the solver over it.** The split matters: a layer count
above the compile-time ceiling or a saturated water content below the residual produces a silently
wrong column rather than a crash, so the physical block is validated at load. `depth` is the knob
for the known too-shallow-column defect — the default 2.0 m sits against a roughly 2.5 m annual
thermal damping depth, and no depth in 2–3 m is converged (issue #145).

### The PFT file

Everything PFT-specific: the `[pft]` trait table, the `[camac]` mortality-hazard derivation
coefficients, and the `[allometry]` coefficients. **The number of PFTs is the length of
`[pft].wood_density`** — every other trait array must match it.

Wood density is the primary PFT axis. The Camac-2018 mortality coefficients are derived from it as
power laws, so low-density PFTs get a higher growth-independent hazard and a steeper low-growth
penalty. Growth, competition and reproduction parameters are per-PFT but ship uniform.

## How a run starts

`[init].init_mode` selects one of three, and the files for the unselected modes are ignored rather
than removed:

| Mode | Start from | Needs |
|---|---|---|
| `0` | **near-bare ground** (the default) | nothing |
| `1` | a **cohort census** | `[init].census_file` — a CSV with one row per cohort: `site_id, patch_id, cohort_id, dbh, height, pft, nplant`. `dbh` drives the allometry. |
| `2` | a **state checkpoint** | `[init].restart_file` — a `<prefix>-S-*.nc` written by a previous run. Continues the exact instantaneous state. |

A census is how you start from a field inventory; see
[`examples/example_demography/census_example.csv`](../examples/example_demography/census_example.csv)
and `init_from_census` in [`../src/init/meds_init.f90`](../src/init/meds_init.f90). Unusable input
falls back to near-bare ground with a warning.

## Output

Two streams, both under `[io]` with the stem `<output_dir>/<output_prefix>`:

- **Diagnostic timeseries** — the `[output]` subsystem. Around 208 variables across 8 groups and
  7 axes, each switchable individually per timescale (sub-daily, daily, monthly, annual). Run
  `meds_main --dump-io-config` to generate a file listing every available variable name; the
  override mechanism always worked, what was missing was any way to learn what exists.
  See [`science/diagnostics.md`](science/diagnostics.md).
- **State checkpoints** — `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, enabled with `[io].write_state`. The
  instantaneous prognostic state only, no diagnostics, written every `state_interval_years` and at
  run end. **The timestamp is the simulated date**, so pointing `[init].restart_file` at one
  resumes from exactly that date.

A checkpoint is raw prognostic state at an instant and a diagnostic record is a time average; the
two streams are deliberately separate because conflating them is how a restart silently starts from
a mean.

## The forcing file

`[forcing].forcing_on` turns the meteorological driver on; `[site]` says where the site is.
Without it the fast loop runs against a constant reference climate, which is useful for tests and
useless for science.

The forcing file format, the ERA5-Land preparation recipe, and the recycling rules are documented
in [`science/forcing.md`](science/forcing.md). Two things to know before writing a config:

- **The recycle window is declared, never inferred.** If `recycle = true`, then `recycle_start` and
  `recycle_end` are required, and the span must be an exact whole number of calendar years. A
  window of any other length drifts both hour-of-day and day-of-year on every wrap while the daily
  mean stays correct, so nothing downstream complains.
- **MEDS never gap-fills.** A missing or NaN required value is a hard error, not an interpolation.

## Worked examples

The example configs are the fastest way in, and each is a complete pair:

- [`examples/example_demography/`](../examples/example_demography/) — a 250-year demographic
  spin-up from near-bare ground.
- [`examples/example_biophysics/`](../examples/example_biophysics/) — a 50-year spin-up at Ithaca
  with real forcing, then one July restarted at high output resolution.
