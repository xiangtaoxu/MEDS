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
| `[run]` | The slow timestep, and the run span as **calendar dates** (`start_time`, `end_time`, leap-year-aware Gregorian, `"YYYY-MM-DD[ HH:MM:SS]"`). Thread count. The run `mode`: one site, or a region. |
| `[fast]` | The sub-daily loop: `dt_fast`, the integrator, tolerances, error control. |
| `[init]` | How the run starts, and the path to the PFT file. |
| `[demography]` | Cohort and patch fusion/fission, the cadence switches. |
| `[disturbance]`, `[recruitment]` | Treefall rates; the seed-rain and reproduction settings. |
| `[soil_column]` | The **physical ground**: layer count, depth, grid growth, hydraulic texture, retention family, root profile, thermal properties. |
| `[soil]`, `[energy]`, `[snow]`, `[aerodynamics]` | The **solvers over it**: selectors and tolerances. |
| `[soil_carbon]` | The CENTURY decomposition: selectors, rate parameters, cold-start spin-up. |
| `[hydraulics]` | Plant water transport. |
| `[forcing]`, `[site]` | The meteorological driver, and where the site is. |
| `[region]` | For `mode = "region"` only: the box of forcing cells and which of them to simulate. |
| `[output]` | Which diagnostics are written, on which axes, at which timescales. |
| `[state]` | Restart checkpointing: output directory, prefix, interval. *(Renamed from `[io]`; the old name still loads with a warning.)* |
| `[options]` | `override_derived` and other run switches. |

**`[soil_column]` is the ground; `[soil]` is the solver over it.** The split matters: a layer count
above the compile-time ceiling or a saturated water content below the residual produces a silently
wrong column rather than a crash, so the physical block is validated at load.

### The bottom thermal boundary — `[energy].bottom_bc`

The default 2.0 m column sits against an annual thermal damping depth of about 2 m, so the annual
temperature wave has barely attenuated by the time it reaches the base. `bottom_bc` decides what
happens to it there.

| key | values | meaning |
|---|---|---|
| `bottom_bc` | `geothermal` (default) \| `dirichlet` | prescribed bottom flux (held at zero — an adiabatic wall that **reflects** the annual wave) vs conduction to a fixed deep temperature |
| `deep_temp` | [K] | the anchor temperature. **Required** when `bottom_bc = "dirichlet"`; there is no default |
| `deep_depth` | [m], default 3.12 | depth of the anchor plane below the surface; must lie below the bottom soil node |

`deep_temp` is the mean annual soil temperature below the damping depth — close to the mean annual
air temperature of the forcing that drives the run, and a site property like latitude. It is required
rather than defaulted because an error in it is a steady flux into the column base, and so a
mean-annual bias in deep-soil temperature; a silent default would put back, somewhere else, the bias
this boundary condition exists to remove.

`deep_depth` is a physical choice, not a numerical one, and the default is derived: a resistive
termination reflects least when the anchor sits `d/sqrt(2)` below the bottom node, `d` being the
annual damping depth. For the default column that is 1.727 + 1.973/√2 = 3.12 m. **Recompute it if you
change `[soil_column].depth` or the thermal texture** — the derivation, the measured sweep and what
the anchor does and does not fix are in [`docs/science/soil_biophysics.md`](science/soil_biophysics.md).

Deepening the column with `[soil_column].depth` remains available and is the other way out, at the
cost of layers; it also moves the root profile and the drainage, so it is not a thermal-only change.

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

Two streams, both with the stem `<output_dir>/<output_prefix>` from `[state]`:

- **Diagnostic timeseries** — the `[output]` subsystem. Around 208 variables across 8 groups and
  7 axes, each switchable individually per timescale (sub-daily, daily, monthly, annual). Run
  `meds_main --dump-io-config` to generate a file listing every available variable name; the
  override mechanism always worked, what was missing was any way to learn what exists.
  See [`science/diagnostics.md`](science/diagnostics.md).
- **State checkpoints** — `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, enabled with `[state].write_state`. The
  instantaneous prognostic state only, no diagnostics, written every `[state].interval_years` and at
  run end. **The timestamp is the simulated date**, so pointing `[init].restart_file` at one
  resumes from exactly that date.

A checkpoint is raw prognostic state at an instant and a diagnostic record is a time average; the
two streams are deliberately separate because conflating them is how a restart silently starts from
a mean.

> **`[io]` was renamed to `[state]`.** The block was named for a legacy diagnostic writer that was
> retired at v0.1; what remained was the restart stream, so `io` named the one output path it did
> *not* cover. The old spelling still loads in v0.2.x and prints one deprecation warning naming the
> keys; `io.state_interval_years` becomes `state.interval_years` (the `state_` prefix was stuttering
> once the block itself was called `state`). It will be removed in a later release — a 0.x minor is
> the cheapest moment a rename like this will ever have.

## The forcing file

`[forcing].forcing_on` turns the meteorological driver on; `[site]` says where the site is.
Without it the fast loop runs against a constant reference climate, which is useful for tests and
useless for science.

`[forcing].format` picks the source:

- **`"era5land"`** reads the global ED_ERA5land archive: one file per variable per month, built
  with the tools in `scripts/prepare_era5/`. Give it `data_path` (the archive folder) and
  `max_distance_km`. The reader finds the site's cell itself and takes that cell's elevation from
  the archive, so `path`, `grid_index`, `grid_match` and `[site].grid_elevation` do not apply and
  are rejected if present. Every month file the run needs must exist when it starts, and the run
  must step daily from midnight (`dt_slow = "1d"`, `start_time` at 00:00:00): the reader loads a
  month at a time before each step.
- **`"netcdf"`** reads one MEDS forcing file, named by `path`.

The file formats, the ERA5-Land preparation recipe, and the recycling rules are documented in
[`science/forcing.md`](science/forcing.md). Two things to know before writing a config:

- **The recycle window is declared, never inferred.** If `recycle = true`, then `recycle_start` and
  `recycle_end` are required, and the span must be an exact whole number of calendar years. A
  window of any other length drifts both hour-of-day and day-of-year on every wrap while the daily
  mean stays correct, so nothing downstream complains.
- **MEDS never gap-fills.** A missing or NaN required value is a hard error, not an interpolation.

**The forcing is moved to each patch's canopy-air top.** There is no fixed reference height:
every sample is taken from the forcing's own heights to the top of each patch's canopy air space, which
grows with the stand, and the aerodynamics runs from there ([`science/forcing.md`](science/forcing.md) §8).
- **Declare the forcing's own heights in `[forcing]`:** `tq_height` and `wind_height`; `height_above`
  (`"zero_plane"` for a reanalysis, `"ground"` for a flux tower); and `wind_exposure` (`"open_terrain"`
  for ERA5's 10 m wind, with `wind_exposure_z0 = 0.03` and `wind_blending_height = 40`, or `"local"`).
  `meds_config_main.toml` shows ERA5-Land's values.
- **The terrain lapse** (`[site].apply_elevation_lapse`) moves temperature, pressure, humidity (at constant
  relative humidity) and file longwave from the forcing cell's elevation to the site's.
  `[site].lapse_rate_tair` takes one rate or twelve monthly rates.
- The old `[site].reference_height`, `wind_meas_height`, `apply_wind_profile` and `wind_roughness_z0`
  are rejected, with a message naming these keys.

**CO₂ is set in `[forcing]`, never by the met file.**
- **`co2_source = "const"`**, the default, holds `co2_const` for the whole run.
- **`co2_source = "file"`** reads `co2_file`, a MEDS CO₂ file looked up on model time, so the CO₂
  keeps rising while recycled met repeats.
  - The repository ships `data/co2/co2_cmip7_global_annual_1000-2022.txt`, CMIP7 global annual
    means for 1000–2022.
  - The path is relative to where MEDS runs, like `path`.
  - A run the file does not cover is refused at startup.
- **Rejected combinations.** The other mode's key is rejected, and so is a met file that carries
  `CO2air`.

The format (a `timestep` line, a `units` line, then one `<period start> <value>` row per period) and
the shipped series are in [`science/forcing.md`](science/forcing.md) §12.

## Regional runs

`[run].mode = "region"` simulates every selected cell of a box of the ED_ERA5land archive as its own
polygon, in one process. Each polygon is an independent stand at its cell's centre, at the cell's
orography, in UTC, and it computes exactly what a site run at that point would: the region steps
every polygon with the site run's own step. Scattered sites are not a region; run them as separate
site runs (a job array, or one allocation filled with GNU parallel).

```toml
[run]
mode = "region"

[region]
box_nwse          = [42.95, -76.95, 41.95, -75.95]   # [N, W, S, E]; may cross 0 or 180 degrees
land_fraction_min = 0.5                              # skip lakes and fractional coastal cells
detail_polygons   = [1714634]                        # optional: these also write single-site files
```

- **Selection.** Every cell of the box that the archive has data for and whose static land fraction
  is at least `land_fraction_min` (default 0.5) becomes a polygon, in row-major order. A polygon's
  id is its cell's index on the global grid, `row * 3600 + col` at 0.1 degrees, so ids agree across
  regions.
- **Output.** One file per tier per time chunk for the whole region, with a `polygon` dimension and
  the coordinates `polygon_id`, `lat`, `lon`, `row` and `col`. Region files hold the variables that
  have the same shape everywhere: site totals, per PFT, per size class and per soil layer. Cohort
  and patch variables, and the sub-daily tier, are written only for `detail_polygons`, as ordinary
  single-site files named `<prefix>-p<polygon id>-...`.
- **Rules.** A region needs `[forcing].format = "era5land"` with forcing and the fast loop on. It
  takes its locations from its cells, so `[site].latitude`, `longitude`, `utc_offset`, `elevation`
  and `[forcing].max_distance_km` are refused; the rest of `[site]` (reference heights, profile and
  lapse switches) still applies. Until restarts of regions exist, a region starts from bare ground
  (`init_mode = 0`) and writes no checkpoints (`[state].write_state = false`), and it runs one patch
  thread without the fast probe.
- **Cost.** Work and memory grow with the polygon count. Output is written between months, so a
  crash loses at most the current month. See `docs/dev_plans/MEDS_POLYGON_RUNTIME_PLAN.md` for the
  measured cost per polygon-month and the roadmap to threads, restarts and tiles.

## Worked examples

The example configs are the fastest way in, and each is a complete pair:

- [`examples/example_demography/`](../examples/example_demography/) — a 250-year demographic
  spin-up from near-bare ground.
- [`examples/example_biophysics/`](../examples/example_biophysics/) — a 50-year spin-up at Ithaca
  with real forcing, then one July restarted at high output resolution.
